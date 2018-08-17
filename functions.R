# fes tides ####
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

# robolimpet data ####
get.ref <- function(T0, T1, DT, paths) {
	refs <- list()
	for (p in 1:length(paths)) {
		x <- suppressMessages(read_csv(paths[p], col_names = c("time", "temp"), skip = 2))
		x <- xts(x$temp, x$time)
		# trim xout data to match the T_RANGE used in the LSMs
		xout <- seq(max(first(time(x)), ymd_h(T0)), min(last(time(x)), ymd_h(T1)), DT)
		if (length(xout) != nrow(x)) {
			yout <- approx(time(x), x, xout, method = "linear", rule = 1)$y
			x <- xts(yout, xout)
		}
		colnames(x) <- gsub(".txt", "", basename(paths[p]))
		refs[[p]] <- x
	}
	do.call(cbind, refs)
}

# water temperature ####
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

# lsm ####
run.lsm_r_v11 <- function() {
	# set up matrix to store soil layer temperatures at each time step
	t <- matrix(NA, nrow = nrow(w), NSOIL + 1)

	### in loop
	pb <- txtProgressBar(1, NRUN, style = 3)
	for (i in 1:NRUN) {
		# i <- 1
		# grab line [i] of the forcing data tibble
		read.env(i, w)

		# calculate land-surface physics
		housekeeping()
		sflx()

		# store layer temperatures
		t[i,] <- c(TSKIN, STC)
		setTxtProgressBar(pb, i)
	}
	close(pb)

	colnames(t) <- c("tskin", str_c("l", 1:(ncol(t) - 1)))
	xts(t - 273.15, TIMESTAMPS)
}
