Sys.setenv(TZ = "UTC")
options(digits = 10)
pkgs <- c("tidyverse", "lubridate", "stringr", "xts", "dygraphs", "scales", "shiny", "raster", "insol")
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))

# load general functions ####
source("functions.R")

# MAKE saved files "remember" the aggregation factor, so that both raster loading and more importantly adjusted latlon for loggers only has to be done once

# ADD evaluating stats to the dygraph to assess log vs lsm performance

fn <- "lsm.RData"
if (file.exists(fn)) {
	load(fn)
}else{
	ENV <- str_c(dirname(getwd()), "/io/")

	### verificar times locais e utc (quando 1 e qd outro)
	# set main run params ####
	T0  <- "2017-10-10 00"
	T1  <- "2018-03-01 00"
	T1  <- "2017-12-01 00"
	DT  <- 1800
	LOC <- list(lon = -8.876, lat = 41.839)
	aggregFACT <- 8


	# load lsm functions ####
	LSM <- str_c(dirname(getwd()), "/lsm_r_v11/")
	# load parameters
	for (f in dir(LSM, pattern = "PARAMS.", full.names = TRUE)) source(f)
	# load functions
	for (f in dir(LSM, pattern = "FUNS.", full.names = TRUE)) source(f)

	# load 3d model rasters ####
	r <- load.rasters(ENV, aggregFACT = aggregFACT)

	# logger data ####
	dat <- get.loggers(T0, T1, DT, ENV)

	# water temperature ####
	# ... from loggers
	ref <- filter(dat, lvl == "m")$log %>% do.call(merge, .)
	wat <- extract.water(LOC, T0, T1, DT, ref)

	# load weather data ####
	w <- forcing.data(ENV)
	# replace sst from sat by sst from logger data
	w$sst <- approx(time(wat), as.numeric(wat$sst), time(w), method = "linear", rule = 2)$y %>%
		"+"(273.15) %>%
		round(2)

	# adjust sw & lw fields ####
	# raw sw and lw data are modified according to radiation model
	rad <- compute.rad()

	dat$sw <- dat$lw <- dat$log
	for (i in 1:nrow(dat)) {
		dat$sw[[i]] <- rad$sw[,i]
		dat$lw[[i]] <- rad$lw[,i]
	}

	# run lsm ####
	dat$lsm <- list.lsm.run(dat)

	# save environment ####
	save.image(file = "lsm.RData")
}

# prepare data ####
# ... for shiny
# tskin = 0 (aka 'l0'), l1 = 1, ..., l15 = 15
layer <- 1

lsm <- map(dat$lsm, ~.x[,layer + 1]) %>% do.call(merge, .)
colnames(lsm) <- str_c(dat$micro, "_lsm")

log <- do.call(merge, dat$log)
colnames(log) <- str_c(dat$micro, "_log")

out  <- cbind(lsm, log)
fixed.ylim <- range(out)
W <- apply(w, 2, function(x) rescale(x, fixed.ylim * c(1, 0.5)))
out  <- cbind(out, W)
colW <- colnames(w)

# visualize ####
ui <- fluidPage(
	fluidRow(
		column(2,
					 checkboxInput(inputId  = "ylim", label = strong("fixed ylim"), value = TRUE),
					 checkboxInput(inputId  = "hi",   label = strong("highlight"),  value = FALSE),
					 hr(),
					 column(6,
					 			 radioButtons(inputId  = "logger", label = strong("logger"),
					 			 						 choices  = dat$micro, selected = dat$micro[1])),
					 column(6,
					 			 radioButtons(inputId  = "col3", label = strong("env"),
					 			 						 choices = colW, selected = colW[1]))),
		mainPanel(
			fluidRow(
				column(10, dygraphOutput("plot", height = "700px")),
				column(2,  textOutput("legendDivID")))))
)

server <- function(input, output) {
	# Subset data
	subset.data <- reactive({
		col <- c(str_c(input$logger, "_lsm"), str_c(input$logger, "_log"))
		out[, c(col, input$col3)]
	})

	output$plot <- renderDygraph({
		df <- subset.data()
		d <- dygraph(df) %>%
			dyRangeSelector(dateWindow = c(ymd_h(T0), ymd_h(T0) + (6 * 24 * 3600))) %>%
			dyLegend(
				show = "always",
				hideOnMouseOut = FALSE,
				labelsSeparateLines = TRUE,
				labelsDiv = "legendDivID") %>%
			dyOptions(
				retainDateWindow = TRUE,
				colors = c("red", "green", "blue"))
				if (input$ylim) d <- d %>% dyAxis("y", valueRange = fixed.ylim)
		if (input$hi) d <- d %>% dyHighlight(highlightSeriesOpts = list(strokeWidth = 2))
		d
	})
}

shinyApp(ui = ui, server = server)
