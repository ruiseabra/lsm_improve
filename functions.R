fes.tides <- function(loc, t_range, t_res) {
	Sys.setenv(HDF5_DISABLE_VERSION_CHECK = "2")
	# prepare
	ORIGIN  <- ymd("1950-01-01")
	T_RANGE <- (t_range + c(0, t_res)) %>%
		julian(origin = ORIGIN) %>%
		as.numeric %>%
		formatC(format = "f")
	T_RES <- t_res / 60
	CALL  <- str_c("fes_slev", loc$lat, loc$lon, T_RANGE[1], T_RANGE[2], T_RES, sep = " ")
	# run fes_slev
	tides <- system(CALL, intern = TRUE)[-(1:2)] %>%
		str_split(",")
	# extract timestamps
	times <- map_chr(tides, 1) %>%
		as.numeric %>%
		"*"(., (24 * 3600)) %>%
		round
	times <- (times - (times %% 60)) %>%
		as.POSIXct(origin = ORIGIN)
	times <- times - (as.numeric(times) %% 60)
	steady <- diff(times) %>%
		unique %>%
		length %>%
		"=="(., 1)
	if (!steady) stop("the period of 'times' is irregular")
	# extract tide elevation
	tides <- map_chr(tides, 2) %>% as.numeric
	# combine
	tides <- tibble(time = times, tide = tides)
	# filter to ensure that the data return does not exceed the t_range supplied
	tides <- filter(tides, time %within% interval(t_range[1], t_range[2]))
	# return
	tides
}

# collect robolimpet data
get.loggers <- function(T0, T1, DT, ENV) {
	fn1 <- str_c(ENV, "robolimpet/logger_info.csv")
	fn2 <- str_c(ENV, "robolimpet/logger_info_mod.csv")
	if (file.exists(fn2)) {
		dat <- suppressMessages(read_csv(fn2))
	}else{
		dat <- suppressMessages(read_csv(fn1))
		dat$lvl <- str_sub(dat$micro, 1, 1)
		dat$lat <- dat$lat_original
		dat$lon <- dat$lon_original

		### a manual adjustment must be done (once) so that the heights of the pixels
		#     for the robolimpets, extracted from the topo raster, match the real heights
		#     recorded with GPS in the field
		### this is run once, stored in the csv file and later the data is imported and the
		#     adjusted xx and yy are used instead of the original xTRUE and yTRUE
		BRG <- c("black", "red", "grey")
		for (i in 1:nrow(dat)){
			# i <- 1
			h <- dat$h[i]
			x <- dat$lon[i]
			y <- dat$lat[i]
			RNG1 <- 0.05
			RNG2 <- c(-0.25,0.25)
			ext <- extent(c(x + 20 * res(fake_colors)[1] * c(-1,1), y + 10 * res(fake_colors)[1] * c(-1,1)))

			m <- c(-10, h - RNG1, 0,
						 h - RNG1, h + RNG1, 1,
						 h + RNG1, 10, 2)
			rclmat <- matrix(m, ncol = 3, byrow = TRUE)
			fake_colors <- reclassify(r$REFdem, rclmat)
			fake_colors <- crop(fake_colors, ext)

			loc <- list(x = x, y = y)
			repeat {
				old_loc  <- loc
				hnow     <- round(extract(r$REFdem,   cbind(loc$x, loc$y)), 3)
				slopenow <- round(extract(r$REFslope, cbind(loc$x, loc$y)), 1)
				plot(fake_colors, col = BRG, xlim = x + RNG2, ylim = y + RNG2, axes = FALSE, main = paste(dat$micro[i], "\nheight now:", hnow, " - real height:", h, "  - diff to real:", round(h - hnow, 3), " - slope:", slopenow), col.main = ifelse(abs(h - hnow) < RNG1, "green", "black"))
				points(x,     y,     col = "yellow", cex = 3, pch = 10, lwd = 3)
				points(loc$x, loc$y, col = "green",  cex = 3, pch = 10, lwd = 3)
				loc <- locator(1)
				if (loc$y > par("usr")[4]) {
					dat$lon[i] <- old_loc$x
					dat$lat[i] <- old_loc$y
					break
				}
			}
		}
		lonlat <- cbind(dat$lon, dat$lat)
		dat$h     <- round(extract(r$REFdem,   lonlat), 3)
		dat$slope <- round(extract(r$REFslope, lonlat), 3)
		dat$svf   <- round(extract(r$REFsvf,   lonlat), 3)
		dat$ind   <- cellFromXY(r$REFdem, lonlat)

		dat <- dplyr::select(dat, micro, lvl, lat, lon, ind, h, slope, svf)
		write_csv(dat, path = fn2)
	}

	paths <- str_c(str_c(ENV, "robolimpet/"), dat$micro, ".txt")
	if (any(!file.exists(paths))) stop("some logger files are missing")
	dat$log <- map(paths, ~get.logger(T0, T1, DT, .x))
	dat
}

get.logger <- function(T0, T1, DT, path) {
	x <- suppressMessages(read_csv(path, col_names = c("time", "temp"), skip = 2))
	x <- xts(x$temp, x$time)
	# trim xout data to match the T_RANGE used in the LSMs
	xout <- seq(max(first(time(x)), ymd_h(T0)), min(last(time(x)), ymd_h(T1)), DT)
	if (length(xout) != nrow(x)) {
		yout <- approx(time(x), x, xout, method = "linear", rule = 1)$y
		x <- xts(yout, xout)
	}
	colnames(x) <- gsub(".txt", "", basename(path))
	x
}

# extract water temperature from loggers
extract.water <- function(loc, T0, T1, dt, ref) {
	tid <- fes.tides(loc, c(ymd_h(T0, T1)), dt)
	tid$hi <- rollapply(tid$tide, 5, function(x) which.max(x) == 3, fill = FALSE)
	tid <- filter(tid, hi)$time

	wat <- ref[time(ref) %in% tid, ] %>%
		apply(1, mean) %>%
		xts(tid)
	colnames(wat) <- "sst"

	wat
}

list.lsm.run <- function(dat) {
	t  <- list()
	pb <- txtProgressBar(1, nrow(dat), style = 3)
	for (i in 1:nrow(dat)) {
		W <- w
		# use modified radiation
		W$sw <- dat$sw[[i]]
		W$lw <- dat$lw[[i]]

		# transform tide heights into tide in or tide out
		# 1 = underwater, 0 = out-of-water
		W$tide <- ifelse(W$tide > dat$h[i], 1, 0)

		# run lsm
		t[[i]] <- run.lsm_r_v11(W)

		setTxtProgressBar(pb, i)
	}
	close(pb)
	t
}


run.lsm_r_v11 <- function(W) {
	# set up matrix to store soil layer temperatures at each time step
	t <- matrix(NA, nrow = nrow(W), NSOIL + 1)

	### in loop
	for (i in 1:NRUN) {
		# i <- 1
		# grab line [i] of the forcing data tibble
		read.env(i, W)

		# calculate land-surface physics
		housekeeping()
		sflx()

		# store layer temperatures
		t[i,] <- c(TSKIN, STC)
	}

	colnames(t) <- c("tskin", str_c("l", 1:(ncol(t) - 1)))
	xts(t - 273.15, TIMESTAMPS)
}

load.rasters <- function(ENV, aggregFACT = NULL) {
	# lisbon projection (see http://www.spatialreference.org/ref/epsg/)
	lisbon <- "+proj=tmerc +lat_0=39.66666666666666 +lon_0=1 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +pm=lisbon +units=m +no_defs"

	# load dem
	dem <- raster(str_c(ENV, "3d/dem.asc"))
	if (!is.null(aggregFACT)) dem <- aggregate(dem, fact = aggregFACT, fun = mean)

	# load svf
	# svf fields calculated by importing the 3D model to GRASS and using the function r.skyview
	svf <- raster(str_c(ENV, "3d/svf.asc"))
	if (!is.null(aggregFACT)) svf <- aggregate(svf, fact = aggregFACT, fun = mean)
	# match NA values in dem
	svf[is.na(dem[])] <- NA

	# project raster with lisbon projection
	EXTdem <- extent(dem)

	dem <- projectRaster(dem, crs = projection(lisbon))
	svf <- projectRaster(svf, crs = projection(lisbon))

	# vectors normal to each grid cell
	normal <- cgrad(dem)

	# slope & aspect
	slope  <- normal %>%  slope(degrees = TRUE) %>% raster(crs = projection(lisbon))
	aspect <- normal %>% aspect(degrees = TRUE) %>% raster(crs = projection(lisbon))
	extent(slope) <- extent(aspect) <- extent(dem)
	slope[is.na(dem[])] <- aspect[is.na(dem[])] <- NA

	REFdem   <- dem
	REFsvf   <- svf
	REFslope <- slope
	extent(REFdem) <- extent(REFsvf) <- extent(REFslope) <- EXTdem

	list(proj = lisbon, dem = dem, REFdem = REFdem, svf = svf, REFsvf = REFsvf, normal = normal, slope = slope, REFslope = REFslope, aspect = aspect)
}

compute.rad <- function() {
	# lisbon projection (see http://www.spatialreference.org/ref/epsg/)
	lisbon <- "+proj=tmerc +lat_0=39.66666666666666 +lon_0=1 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +pm=lisbon +units=m +no_defs"

	# compute sun vector
	sunV <- sunvector(JD(TIMESTAMPS), LOC$lat, LOC$lon, timezone = 0)

	SWrad <- LWrad <- xts(matrix(NA, length(TIMESTAMPS), nrow(dat)), TIMESTAMPS)
	colnames(SWrad) <- colnames(LWrad) <- dat$micro

	pb <- txtProgressBar(1, NRUN, style = 3)
	for (h in 1:NRUN) {
		# h <- 35; TIMESTAMPS[h]

		Tnow    <- TIMESTAMPS[h]
		sunVnow <- sunV[h,, drop = FALSE]
		azimuth <- sunpos(sunVnow)[,"azimuth"]
		zenith  <- sunpos(sunVnow)[,"zenith"]

		## project light on the rock
		# light intensity
		light <- hillshading(r$normal, sunVnow) %>% raster(crs = lisbon)
		extent(light) <- extent(r$dem)
		# cast shadows
		shadow <- doshade(r$dem, sunVnow)
		# remove shaded areas
		light.minus.shade <- light * shadow

		## load current atmospheric flux
		# shortwave (UV and Visible, 280-3000 nm) global radiation at this location
		sw <- as.numeric(w$sw[Tnow])
		# longwave (IR, >3000 nm) global radiation at this location
		lw <- as.numeric(w$lw[Tnow])

		## determine the ratio of diffuse to total solar radiation
		#   as described in Tian et al 2001, Agricultural and Forest Metereology 109:67-74
		# Iqbal 1983 pg 61, yields the instantaneous radiation for 1 hour centered around the Zenith angle for the midpoint of the hour
		# Solar constant of 1367 W/m2, London & Frohlich (1982), cited by Pons & Ninyerola (2008) Int J Clim 28:1821-1834  #1367 W/m2 also in Iqbal 1983
		Isc <- 1367
		# Iqbal Eqn 1.2.2
		# E = eccentricity of the earth's orbit
		angle <- 2 * pi * (yday(Tnow) - 1) / 365
		E <- 1.000110 +
			0.034221 * cos(angle) +
			0.001280 * sin(angle) +
			0.000719 * cos(2 * angle) +
			0.000077 * sin(2 * angle)
		# cosine of incidence of extraterrestrial radiation
		I0 <- Isc * E * cos(radians(zenith))
		I0 <- replace(I0, I0 < 0.0001, 0)

		# shortwave radiation cannot exceed radiation at the top of the atmosphere
		if (sw > I0) sw <- I0

		## get diffuse radiation
		# piecewise function following the ER method described in Bindi et al (1992) Clim Res 2:47-54.
		#   they cite Erbs, D. G., Klein, S. A., Duffie, J. A. (1982).
		#   estimation of the diffuse radiation fraction for hourly, daily and monthly average global radiation
		#   Solar Energy 28: 293-302
		# % of diffuse
		if (sw == 0) {
			trans    <- 0
			diffFrac <- 1
		}else{
			trans <- sw / I0
			if (trans <= 0.22) diffFrac <- 1 - (0.09 * trans)
			if (trans > 0.22 & trans <= 0.80) {
				diffFrac <- 0.9511 -
					0.1604 * (trans)   +
					4.388  * (trans^2) -
					16.638  * (trans^3) +
					12.336  * (trans^4)
			}
			if(trans > 0.8) diffFrac <- 0.165
		}

		# decomposing total radiation into direct (beam) and diffuse
		sw.indirect <- sw * diffFrac
		sw.direct   <- sw * (1 - diffFrac)

		## projecting radiation into surfaces
		# direct radiation
		sw.direct.ground <- sw.direct * light.minus.shade

		# indirect
		# Klutcher 1979
		F0 <- 1 - (diffFrac)^2 # he uses ^2 here
		F1 <- 1 + F0 * sin(radians(r$slope / 2))^3 # horizon brightning
		F2 <- 1 + F0 * light^2 * sin(radians(zenith))^3 # circumsolar brightning, note the + here
		sw.indirect.ground <- sw.indirect * r$svf * F1 * F2 # note the + here

		# longwave is considered isotropic, and it is shaded by topography
		# proportional to sky view factor
		lw.ground <- lw * r$svf

		# reflected
		# Temps and Coulson, described in p 6 from
		#   Measuring and Modeling Solar Irradiance on Vertical Surfaces, Maxwell, Stoffel and Bird
		# cast shadows of reflected light
		shadow.reflection <- doshade(r$dem, (sunVnow * c(0.5, 0.5, -1)))

		M5 <- (1 - cos(radians(r$slope / 2))^2)
		surf.azimuth <- abs(azimuth - r$aspect) %% 360
		change <- surf.azimuth[] > 180
		surf.azimuth[change] <- 360 - surf.azimuth[change]
		# surface receives more reflected if sun is low and if it is facing the sun
		M6 <- (1 + sin(radians(zenith / 2))^2) * abs(cos(radians(surf.azimuth)))
		sw.reflected <- sw * ALBEDO * M5 * M6
		sw.reflected <- sw.reflected * shadow.reflection

		# total
		sw <- sw.direct.ground + sw.indirect.ground + sw.reflected
		lw <- lw.ground
		SWrad[h,] <- sw[dat$ind]
		LWrad[h,] <- lw[dat$ind]

		setTxtProgressBar(pb, h)
	}
	close(pb)
	list(sw = SWrad, lw = LWrad)
}
