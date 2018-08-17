## original file: export_radiations_to_model_v5

OBJECTS <- ls()

nbytes <- 2

# projections
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
# plate caree projection http://spatialreference.org/ref/epsg/wgs-84-plate-carree/
plate.carree <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "
# lisbon projection (see http://www.spatialreference.org/ref/epsg/)
lisbon <- "+proj=tmerc +lat_0=39.66666666666666 +lon_0=1 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +pm=lisbon +units=m +no_defs"
# etrs89
etrs89 <- "+init=epsg:4258"

# get xyz data
topo <- raster(paste0(IO_folder, "3d_model/squared/3d_model_squared_pixel.asc"))
topo <- aggregate(topo, fact=aggregFACT_3dmodel, fun=mean)

# get sky view factor data
# svf fields calculated by importing the 3D model to GRASS and using the function r.skyview
svf <- raster(paste0(IO_folder, "3d_model/squared/Skyview_factor.asc"))
svf <- aggregate(svf, fact=aggregFACT_3dmodel, fun=mean)
# correct NA values
data2change <- values(svf)
data2change[is.na(values(topo))] <- NA
svf[] <- data2change

# central coordinates
central.long <- mean(topo@extent[1:2])
central.lat  <- mean(topo@extent[3:4])

# project raster with lisbon projection
projected.topo <- projectRaster(topo, crs=projection(lisbon))
projected.svf  <- projectRaster(svf,  crs=projection(lisbon))

# slope
slope.raster <- slope(cgrad(projected.topo), degrees=T)
slope.raster <- raster(slope.raster, crs=projection(lisbon))
extent(slope.raster) <- extent(projected.topo)

# timestamps
datejulian.UTC <- JD(matrix(times1, nrow=1, byrow=T))
nsteps         <- length(datejulian.UTC)

# sun vector
sv.UTC <- sunvector(datejulian.UTC, central.lat, central.long, location$tz)

# vectors normal to each grid cell
normal.vectors <- cgrad(projected.topo)


# find the indexes of pixels to save
# this has to be done using total.rad.sw because lw has more non-NA values
#   (it has not been processed using 'this.light', which due to the calculations of the projection of light results in the trimming of some of the border pixels)
# a <- total.rad.lw; a[!is.na(total.rad.sw[])] <- NA; plot(a, main="these pixels have values\nin 'total.rad.lw' but not in 'total.rad.sw'")
ind <- raster(hillshading(normal.vectors, sv.UTC[1,]), crs=lisbon)
extent(ind) <- extent(projected.topo)
ind <- which(!is.na(ind[]))

#robolimpet data
rCoords <- read.csv(paste0(IO_folder, "logger_info.csv"))
rCoords$xTRUE <- rCoords$x
rCoords$yTRUE <- rCoords$y
rCoords$heightTRUE <- rCoords$height

# a manual adjustment must be done (once) so that the heights of the pixels for the robolimpets, extracted from the topo raster, match the real heights recorded with GPS in the field
# this is run once, stored in the csv file and later the data is imported and the adjusted xx and yy are used instead of the original xTRUE and yTRUE
if(F)
{
  for(r in 1:nrow(rCoords))
  {
    #r <- 1
    h <- rCoords$height[r]
    x <- rCoords$xTRUE[r]
    y <- rCoords$yTRUE[r]
    RNG1 <- 0.05
    RNG2 <- c(-0.25,0.25)

    m <- c(-10, h - RNG1, 0,
       h - RNG1, h + RNG1, 1,
       h + RNG1, 10, 2)
    rclmat <- matrix(m, ncol=3, byrow=T)
    fake_colors <- reclassify(projected.topo, rclmat)
    #plot(projected.topo, xlim=x+RNG2, ylim=y+RNG2)
    hnow <- round(extract(projected.topo, cbind(x,y)), 3)
    plot(fake_colors, col=c("black", "red", "grey"), xlim=x+RNG2, ylim=y+RNG2, axes=F, main=paste("have:", hnow, " - want:", h))

    loc <- list(x=x,y=y)
    repeat
    {
      old_loc  <- loc
      hnow     <- round(extract(projected.topo, cbind(loc$x,loc$y)), 3)
      slopenow <- round(extract(slope.raster, cbind(loc$x,loc$y)), 1)
      plot(fake_colors, col=c("black", "red", "grey"), xlim=x+RNG2, ylim=y+RNG2, axes=F, main=paste(rCoords$micro1[r], "\nheight now:", hnow, " - real height:", h, "  - diff to real:", round(h - hnow, 3), " - slope:", slopenow), col.main=ifelse(abs(h - hnow) < RNG1, "green", "black"))
      points(x, y, col="yellow", cex=3, pch=10, lwd=3)
      points(loc$x, loc$y, col="green", cex=3, pch=10, lwd=3)
      loc <- locator(1)
      if(loc$y > par("usr")[4])
      {
        rCoords$x[r] <- old_loc$x
        rCoords$y[r] <- old_loc$y
        rCoords$height[r] <- round(extract(projected.topo, cbind(old_loc$x, old_loc$y)), 3)
        rCoords$slope[r]  <- round(extract(slope.raster,   cbind(old_loc$x, old_loc$y)), 3)
        break
      }
    }
  }
  write.csv(rCoords, file=paste0(IO_folder, "logger_info.csv"), row.names=F)
}
#rPoints <- cbind(rCoords$lon, rCoords$lat)
#rownames(rPoints) <- rCoords$micro1
#rPoints <- SpatialPoints(rPoints)
#projection(rPoints) <- etrs89
#robolimpet.proj     <- spTransform(rPoints, CRS(lisbon))

# compute adjusted shortwave and longwave fluxes for each pixel of the 3d model
swaveFiles <- paste0(OUTPUTS$rad_time_slices, "total.rad.sw.", format(times1, "%Y%m%d%H%M"), ".bin")
lwaveFiles <- paste0(OUTPUTS$rad_time_slices, "total.rad.lw.", format(times1, "%Y%m%d%H%M"), ".bin")

t0 <- Sys.time()
cat(paste("\n\n----", t0, "\n     computing shortwave and longwave flux\n     for", nsteps, "time steps\n\n"))
# parallel processing
foreach(h=1:nsteps, .inorder=F, .packages=package.list) %dopar%
{
  if(!file.exists(check)) stop('process halted by user')

  #h <- 1
  #times1[h]
  # projection of light on the rock
  this.date  <- times1[h]
  this.sv    <- sv.UTC[h,] #sun vector
  this.light <- raster(hillshading(normal.vectors, this.sv), crs=lisbon) #light intensity
  extent(this.light) <- extent(projected.topo)
  shadow <- doshade(projected.topo, this.sv) #cast shadows
  light.minus.shade <- round(this.light * shadow, 4) #remove shade areas

  # shortwave (UV and Visible, 280-3000 nm) global radiation at this location
  sw <- as.numeric(forcing$swave[this.date])
  # longwave (IR, >3000 nm) global radiation at this location
  lw <- as.numeric(forcing$lwave[this.date])

  # determining the ratio of diffuse to total solar radiation described in Tian et al 2001, Agricultural and Forest Metereology 109:67-74
  # using package sirad
  #require(sirad)
  #ET.rad<-extrat(dn,radians(central.lat))
  # conversion from megajoules/m2 to kw.h/m2
  #c(ET.rad$ExtraTerrestrialSolarRadiationHourly*1000000/3600)
  ## Iqbal 1983 pg 61, yields the instantaneous radiation for 1 hour centered around the Zenith angle for the midpoint of the hour
  # Solar constant of 1367 W/m2, London & Frohlich (1982), cited by Pons & Ninyerola (2008) Int J Clim 28:1821-1834  #1367 W/m2 also in Iqbal 1983
  Isc <- 1367
  # E=eccentricity of the earth's orbit : Iqbal Eqn 1.2.2
  # Iqbal Eqn 1.2.2
  day_angle <- 2 * pi * (daydoy(this.date) - 1) / 365
  # E=eccentricity of the earth's orbit : Iqbal Eqn 1.2.1
  E  <- 1.000110 + 0.034221 * cos(day_angle) + 0.001280 * sin(day_angle) + 0.000719 * cos(2 * day_angle) + 0.000077 * sin(2 * day_angle)
  # calculate cosine of incidence of extraterrestrial radiation
  I0 <- Isc * E * cos(radians(sunpos(sv.UTC[h,,drop=F])[,'zenith']))
  I0 <- replace(I0, I0 < 0.0001, 0)

  ## shortwave radiation cannot exceed radiation at the top of the atmosphere
  if(sw > I0) sw <- I0

  # piecewise function to get diffuse radiation following the ER method described in Bindi et al (1992) Clim Res 2:47-54. They cite Erbs, D. G., Klein, S. A., Duffie, J. A. (1982). Estimation of the diffuse radiation fraction for hourly, daily and monthly average global radiation. Solar Energy 28: 293-302
  # % of diffuse
  if(sw == 0)
  {
    trans <- 0
    diffFrac    <- 1
  }else{
    trans <- sw / I0
    if(trans <= 0.22) diffFrac <- 1 - (0.09 * trans)
    if(trans > 0.22 & trans <= 0.80)
    {
      diffFrac <- 0.9511 -
        0.1604 * (trans)   +
        4.388  * (trans^2) -
        16.638 * (trans^3) +
        12.336 * (trans^4)
    }
    if(trans > 0.8) diffFrac <- 0.165
  }

  # decomposing total radiation into direct (beam) and diffuse
  sw.indirect <- sw * diffFrac
  sw.direct   <- sw * (1 - diffFrac)

  # projecting radiation into surfaces
  # direct radiation
  sw.direct.ground <- round(sw.direct * light.minus.shade, 2)

  # indirect
  sw.IDB <- sw.indirect * projected.svf #proportional to sky view factor
  F0     <- 1 - (diffFrac)^2
  F1     <- 1 + (F0 * sin(radians(slope.raster/2))^3)
  sun.zenith <- sunpos(rbind(this.sv))[, "zenith"]
  sun.zenith <- ifelse(sun.zenith > 90, 90, sun.zenith)
  F2     <- 1 + (F0 * (this.light)^2) * (sin(radians(sun.zenith))^3)
  sw.indirect.ground   <- round(sw.IDB * F1 * F2, 2)

  # longwave is considered isotropic, and it is shaded by topography
  lw.IDB <- lw * projected.svf #proportional to sky view factor

  # but, in addition to blocking background longwave radiation, surrounding rocks also emit their own longwave
  # the code below computes the emited value, but we disregarded this because the values were negligenciable (at 20 C it's about 6 W/m2)
  # # Stefan-Boltzmann law, https://en.wikipedia.org/wiki/Stefan-Boltzmann_law
  #sigma <- 5.670373 * 10^-8 # Stefan-Boltzmann constant
  #E <- 0.75 # emissivity of the rock
  #temp <- 20
  #Emited_LW <- E * sigma * temp^4 * 1000 # in W/m2

  lw.ground <- round(lw.IDB, 2)

  # reflected

  #temps and coulson anisotropic
  surf.angle <- 0 #general surface azimuth angle relative to the sky dome (zero on horizontal ground)
  M5 <- (1 - cos(radians(slope.raster/2))^2)
  M6 <- (1 + sin(radians(sun.zenith/2))^2 * (cos(radians(surf.angle))))
  sw.reflected <- round(sw * ground_albedo * M5 * M6 * projected.svf, 2)

  #total
  total.rad.sw  <- sw.direct.ground + sw.indirect.ground + sw.reflected
  total.rad.lw  <- lw.ground
  if(h == 1) save(total.rad.sw, file="tmp.RData")

  #plot(total.rad.swlw,col=grey(seq(0,1,0.01)),main=paste0("total radiation at ", substr(times1[h],1,16)))

  #export radiation data to binary files
  #swave
  file.con <- file(swaveFiles[h], "wb")
  total.rad.sw <- round(total.rad.sw[ind] * 10)
  writeBin(as.integer(total.rad.sw), file.con, size=nbytes)
  close(file.con)
  #lwave
  file.con <- file(lwaveFiles[h], "wb")
  total.rad.lw <- round(total.rad.lw[ind] * 10)
  writeBin(as.integer(total.rad.lw), file.con, size=nbytes)
  close(file.con)
}

#export topography for biophysical model, tiff format
fname <- paste0(IO_folder, "topo.tif")
writeGDAL(as(projected.topo,"SpatialPixelsDataFrame"), fname, mvFlag=-9999)

#export topography for biophysical model, now in binary
#CAUTION##############
#to reconstruct heights, values must be divided by 100!!!!!
nbytes   <- 2
file.con <- file(paste0(IO_folder, "topo.bin"), "wb")
writeBin(as.integer(round(na.omit(as.vector(projected.topo))*100)), file.con, size=nbytes)
close(file.con)

## the radiation calculations seem to be trimming the area of the projected.topo by a few pixels (most likely the calculation of the slope requires pixels surrounding the target pixels, and therefore the pixels at the edge of the topo file do not receive a slope value, and are thus ommited from the radiation calculations)
## in order to make the topo and radiation files match in extend, a trimmed version of the projected.topo is saved as well as the original one
x  <- projected.topo
x[!is.na(x[])] <- 1
x[is.na(x[])]  <- 0
load("tmp.RData")
y  <- total.rad.sw
unlink("tmp.RData")
y[!is.na(y[])] <- 1
y[is.na(y[])]  <- 0
xy <- x-y
topo.trimmed <- projected.topo
topo.trimmed[xy[] == 1] <- NA

#export topography for biophysical model, tiff format
fname <- paste0(IO_folder, "topo_trimmed.tif")
writeGDAL(as(topo.trimmed, "SpatialPixelsDataFrame"), fname, mvFlag=-9999)

#export topography for biophysical model, now in binary
#CAUTION##############
#to reconstruct heights, values must be divided by 100!!!!!
nbytes   <- 2
file.con <- file(paste0(IO_folder, "topo_trimmed.bin"), "wb")
writeBin(as.integer(round(na.omit(as.vector(topo.trimmed))*100)), file.con, size=nbytes)
close(file.con)

## identify the indexes of robolimpets in "topo_trimmed.bin"
topo_seq  <- topo.trimmed
topo_seq[!is.na(topo_seq[])] <- seq(sum(!is.na(topo_seq[])))
refFile1  <- paste0(IO_folder, "logger_info.csv")
refFile2  <- paste0(IO_folder_OUTPUTS, "logger_info.csv")
REF       <- read.csv(refFile1)
ind       <- extract(topo_seq, rCoords[,c("x","y")])
REF$ind2  <- ind
REF$slope <- round(extract(slope.raster,  rCoords[,c("x","y")]), 3)
REF$svf   <- round(extract(projected.svf, rCoords[,c("x","y")]), 3)
write.csv(REF, file=refFile1, row.names=F)
write.csv(REF, file=refFile2, row.names=F)
if(F)
{
  t <- topo_seq
  t[] <- NA
  t[topo_seq[] %in% ind] <- 1
  plot(t)
  points(robolimpet.proj)
}

elapsed <- round(difftime(Sys.time(), t0, units="hours"), 2)
cat(paste("\n\n----", Sys.time(), "\n     done computing shortwave and longwave flux\n     elapsed time:", elapsed, "hours\n"))
units(elapsed) <- "mins"
cat(paste("                  ", elapsed, "mins\n\n"))

OBJECTS <- setdiff(ls(), OBJECTS)
rm(list=OBJECTS)
invisible(gc())
