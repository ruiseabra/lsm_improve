Sys.setenv(TZ = "UTC")
options(digits = 10)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xts))
source("TEST_functions.R")

t0  <- "2017-10-10 00"
t1  <- "2018-03-01 00"
dt  <- 1800
loc <- list(lon = -8.876, lat = 41.839, height = 30)
H <- tibble(
  l  = c("m", "t", "s"),
  nm = set_names(c("1_M", "2_T", "3_S"), l),
  h  = c(30, 50, 70))
hh <- 1

temps <- tibble(
  nm  = H$nm,
  l   = H$l,
  h   = H$h,
  t.r = as.list(rep(NA, nrow(H))),
  t.f = t.r,
  ref = t.r)

# load robolimpet data
ref <- tibble(
  path = dir("TEST_robolimpet_data/", full.names = TRUE),
  l    = basename(path) %>% str_sub(1,1),
  nm   = H$nm[match(l, H$l)])
ref$t <- as.list(rep(NA, nrow(ref)))

for (i in 1:nrow(ref)) {
  x <- suppressMessages(read_csv(ref$path[i], col_names = c("time", "temp"), skip = 2))
  # trim xout data to match the T_RANGE used in the LSMs
  xout <- seq(max(first(x$time), ymd_h(t0)), min(last(x$time), ymd_h(t1)), dt)
  if (length(xout) != nrow(x)) {
    yout <- approx(x$time, x$temp, xout, method = "linear", rule = 1)$y
    x <- tibble(time = xout, temp = yout)
  }
  ref$t[[i]] <- x
}
for (i in 1:nrow(temps)) {
  x  <- filter(ref, nm == temps$nm[i])
  xx <- x$t[[1]]
  if (nrow(x) > 1) {
    for (ii in 2:nrow(x)) {
      xx[[str_c("ref_", ii)]] <- x$t[[ii]]$temp
    }
  }
  colnames(xx) <- c("time", str_c("ref_", 1:nrow(x)))
  temps$ref[[i]] <- xx
}

# get water temperature from loggers
source("lsm.FUNS.env_tide.R")
tid <- fes.tides(loc, c(ymd_h(t0, t1)), dt)
tid$hi <- rollapply(tid$tide, 5, function(x) which.max(x) == 3, fill = FALSE)
tid <- filter(tid, hi)$time

WAT <- filter(ref, l %in% c("l", "m"))
wat <- list()
for (i in 1:nrow(WAT)) wat[[i]] <- filter(WAT$t[[i]], time %in% tid)$temp %>% as.numeric
wat <- do.call(cbind, wat)
wat <- tibble(time = tid, temp = apply(wat, 1, mean))

# reduce the number of heights during first stages of testing
# to make script faster
ref   <-   ref[  ref$l %in% H$l[hh], ]
temps <- temps[temps$l %in% H$l[hh], ]

# run LSMs
for (i in 1:nrow(temps)) {
  loc$height <- temps$h[i]
  temps$t.r[[i]] <- lsm.r(t0, t1, dt, loc, TRUE)
}
DEBUG.r <- read_csv("DEBUG.csv")

for (i in 1:nrow(temps)) {
  loc$height <- temps$h[i]
  temps$t.f[[i]] <- lsm.f(c(1), t0, t1, dt, loc) # lsm.f(c(1,3,5), t0, t1, dt, loc)
}

# merge data
t <- list()
for (i in 1:nrow(temps)) {
  t[[i]] <- tibble(
    nm   = temps$nm[i],
    time = temps$ref[[i]]$time, 
    ref  = temps$ref[[i]]$ref_1, 
    lsmR = temps$t.r[[i]]$l1, 
    lsmF = temps$t.f[[i]]$l1)
}
t <- do.call(rbind, t)
REF_mean <- mean(t$ref) - 3
tr <- select(t, -lsmF) %>% cbind("R", .)
tf <- select(t, -lsmR) %>% cbind("fortran", .)
td <- select(t, -lsmR, -lsmF) %>% add_column(lsm = t$lsmR - t$lsmF + REF_mean) %>% cbind("R - fortran", .)
colnames(tr) <- colnames(tf) <- colnames(td) <- c("type", "nm", "time", "ref", "lsm")
t <- rbind(tr, tf, td) %>% as_tibble
t$lsm <- round(t$lsm, 2)

# plot all data
tt <- filter(t, time > ymd_h("2017-11-01 00") & time < ymd_h("2017-11-20 00"))
# ggplot(tt) +
#   geom_line(aes(time, ref), col = "darkgrey", size = 1) +
#   geom_line(aes(time, lsm), col = "blue") +
#   facet_grid(nm ~ type) +
#   xlab("") + ylab("") +
#   theme(axis.title.x = element_blank(),
#         axis.text.x  = element_blank(),
#         axis.ticks.x = element_blank())

bias_all  <- mean((td$lsm - REF_mean)) %>% round(3)
bias_here <- mean((filter(tt, type == "R - fortran")$lsm - REF_mean)) %>% round(3)
ggplot(tt) +
  geom_hline(aes(yintercept = REF_mean)) +
  geom_line(aes(time, lsm, color = type)) +
  xlab("") + ylab("") + ggtitle(str_c("bias all: ", bias_all, " / bias here: ", bias_here))
  
