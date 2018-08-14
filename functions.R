get.VARLIST <- function() {
  x <- vector(length = length(VARLIST))
  for (i in 1:length(VARLIST)) {
    xx <- first(get(VARLIST[i]))
    x[i] <- if (length(xx) > 1) xx[5] else xx
  }
  x
}

lsm.r <- function(t0, t1, dt, loc, debug = FALSE) {
  print("running lsm --> r version")
  LS <- ls(envir = .GlobalEnv)
  
  T0  <<- ymd_h(t0)
  T1  <<- ymd_h(t1)
  DT  <<- dt
  LOC <<- loc
  
  # load parameters
  for (f in dir(pattern = "lsm.PARAMS.")) source(f)

  # load functions
  for (f in dir(pattern = "lsm.FUNS.")) source(f)

  # prepare forcing data
  w <<- forcing.data()
  
  # update water temperature from logger data
  w$sst <- approx(wat$time, wat$temp, w$time, method = "linear", rule = 2)$y %>% "+"(273.15) %>% round(2)
  w <<- w
  
  # set up matrix to store soil layer temperatures at each time step
  t <- matrix(NA, nrow = nrow(w), NSOIL + 1)
  
  # prepare debug tibble
  if (debug) {
    DEBUG <- matrix(NA, ncol = length(VARLIST), nrow = NRUN) %>% as_tibble
    colnames(DEBUG) <- VARLIST
  }
    
  ### in loop
  pb <- txtProgressBar(1, NRUN, style = 3)
  for (i in 1:NRUN) {
    # i <- 1
    # grab line [i] of the forcing data tibble
    read.env(i)
    
    # calculate land-surface physics
    housekeeping()
    sflx()
    
    # store layer temperatures
    t[i,] <- c(TSKIN, STC)
    
    # dump states of main variables to debug tibble
    if (debug) DEBUG[i,] <- get.VARLIST()
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  t <- cbind(TIMESTAMPS, t - 273.15) %>% as_tibble
  colnames(t) <- c("time", str_c("l", 0:(ncol(t) - 2)))
  t$time <- as.POSIXct(t$time, origin = origin)

  if (debug) write_csv(DEBUG, path = "DEBUG.csv")
  rm(list = setdiff(ls(envir = .GlobalEnv), c(LS, "t")), envir = .GlobalEnv)
  t
}

lsm.f <- function(layers = 1, t0, t1, dt, loc) {
  print("running lsm --> fortran version")
  LS <- ls(envir = .GlobalEnv)
  
  t_range <<- ymd_h(t0, t1)
  t_res   <<- dt
  loc     <<- loc
  
  source("fortran_functions.R")
  timestamps <<- seq.POSIXt(t_range[1], t_range[2], by = t_res)
  
  # set up temporary folder for running the lsm model
  tmp <<- "tmp/"
  unlink(tmp, recursive = TRUE)
  dir.create(tmp)
  
  forcing <<- str_c(tmp, "forcing.in")
  forcing.file(loc, t_range, t_res, forcing)
  
  t <- list()
  for (l in layers) {
    lsm <- compile.lsm(tmp, layer = l)
    out <- system(str_c("cd ", tmp, "; ./", lsm), intern = TRUE, ignore.stderr = TRUE)
    t[[str_c("l", l)]] <- as.numeric(out)
  }
  
  t <- do.call(cbind, t)
  t <- tibble(time = timestamps) %>% cbind(t - 273.15) %>% as_tibble
  
  rm(list = setdiff(ls(envir = .GlobalEnv), c(LS, "t")), envir = .GlobalEnv)
  t
}

lsm.f.debug <- function() {
  print("running lsm --> fortran version")

  t_range <<- ymd_h(t0, t1)
  t_res   <<- dt
  loc     <<- loc
  
  source("fortran_functions.R")
  timestamps <<- seq.POSIXt(t_range[1], t_range[2], by = t_res)
  
  # set up temporary folder for running the lsm model
  tmp <<- "tmp/"
  unlink(tmp, recursive = TRUE)
  dir.create(tmp)
  
  forcing <<- str_c(tmp, "forcing.in")
  forcing.file(loc, t_range, t_res, forcing)
  
  lsm <- compile.lsm(tmp, layer = 1)
  out <- system(str_c("cd ", tmp, "; ./", lsm), intern = TRUE, ignore.stderr = TRUE)
  t <- as.numeric(out)
  
  print(head(t, 20))
}
