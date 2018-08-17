Sys.setenv(TZ = "UTC")
options(digits = 10)
pkgs <- c("tidyverse", "lubridate", "stringr", "xts", "dygraphs", "scales", "shiny")
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))

T0  <- "2017-10-10 00"
T1  <- "2018-03-01 00"
# T1  <- "2017-12-01 00" # must fix how NRUN is computed
DT  <- 1800
LOC <- list(lon = -8.876, lat = 41.839, height = -20)

source("functions.R")

LSM <- str_c(dirname(getwd()), "/lsm_r_v11/")
# load parameters
for (f in dir(LSM, pattern = "PARAMS.", full.names = TRUE)) source(f)
# load functions
for (f in dir(LSM, pattern = "FUNS.", full.names = TRUE)) source(f)

# robolimpet data ####
ENV   <- str_c(dirname(getwd()), "/io/")
paths <- dir(str_c(ENV, "robolimpet/"), full.names = TRUE)
ref <- get.ref(T0, T1, DT, paths)

# water temperature ####
# ... from loggers
wat <- extract.water(LOC, T0, T1, DT, ref[, grepl("msu", colnames(ref))])

# load weather data ####
w <- forcing.data(ENV)
# replace sst from sat by sst from logger data
w$sst <- approx(time(wat), as.numeric(wat$sst), time(w), method = "linear", rule = 2)$y %>%
	"+"(273.15) %>%
	round(2)
# adjust raw shortwave data according to radiation model


# run lsm ####
t <- run.lsm_r_v11()

# prepare data ####
# ... for shiny
out  <- cbind(t, ref)
fixed.ylim <- range(out)
W <- apply(w, 2, function(x) rescale(x, fixed.ylim))
out  <- cbind(out, W)
col1 <- colnames(t)
col2 <- colnames(ref)
col3 <- colnames(w)

# visualize ####
ui <- fluidPage(
	fluidRow(
		column(2,
					 checkboxInput(inputId  = "ylim", label = strong("fixed ylim"), value = FALSE),
					 hr(),
					 column(6,
					 			 checkboxGroupInput(inputId  = "col1", label = strong("lsm temp"),
					 			 									 choices  = col1, selected = col1[c(1,2,4,6)])),
					 column(6,
					 			 checkboxGroupInput(inputId  = "col2", label = strong("ref temp"),
					 			 									 choices  = col2, selected = col2),
					 			 hr(),
					 			 radioButtons(inputId  = "col3", label = strong("env"),
					 			 						 choices = col3, selected = col3[1]))),
		mainPanel(
			fluidRow(
				column(9, dygraphOutput("plot", height = "700px")),
				column(3,  textOutput("legendDivID")))))
)

server <- function(input, output) {
	# Subset data
	subset.data <- reactive({
		out[, c(input$col1, input$col2, input$col3)]
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
				colors = c(heat.colors(length(input$col1)), rep("darkgrey", length(input$col2)), "blue"))
		if (input$ylim) d <- d %>% dyAxis("y", valueRange = fixed.ylim)
		d
	})
}

shinyApp(ui = ui, server = server)
