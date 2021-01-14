# Sam Lambrick, 2019-20
#
# Based on work by Matthew Bergin in the SMF group at the Cavendish Laboratory.
# doi:10.17863/CAM.37853
# doi:10.1016/j.ultramic.2019.112833
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(magrittr)
library(tidyverse)
library(plotly)

# Physical variables used here:
# f      - working distance
# beta   - Virtual source size, assuming small skimmers
# d      - pinhole diameter
# lambda - wavelength of helium


# Define physical constants and other fixed parameters, SI units used throughout
lambda <- 5.63e-11     # Room temperature helium wavelength
B_peak <- 1e20         # Optomised peak brightness (not calculated yet)
gamma <- pi^2*B_peak/4 # Pinhole flux

# Define mathematical functions

# Standard deviation for a pinhole
theta_p <- function(d, beta, f) {
    pinhole <- d/(2*sqrt(3))       # Geometric pattern from pinhole approximated as Gaussian
    source <- beta*f/sqrt(3)       # Broadening due to finite source size
    diffraction <- 0.42*lambda*f/d # Pinhole diffraction, Airy function
    return(sqrt(pinhole^2 + source^2 + diffraction^2))
}

# Optimal pinhole diameter
d_o <- function(sig) sqrt(6)*sig

# Optimal virtual source size
beta_o <- function(sig, f) {
    a <- 0.42*f*lambda
    return((sqrt(3)/f)*sqrt(sig^2/2 - a^2/(d_o(sig)^2)))
}

# Optimal flux
F_o <- function(beta, d) gamma*d^2*beta^2

# Virtual source size from skimmer-pinhole distance
beta_calc <- function(skim, dist) atan(skim/dist)

#produce_plots <- function(xs, ys, )

working_fixed_plot <- function(source_range, pinhole_range, working_dist, skim_x, 
                               flux_base, flux_scaling, axis_scale) {
    # Values for the source distance
    s_range <- source_range
    source_dist <- seq(s_range[1], s_range[2], length.out = 50)*1e-2
    N_s <- length(source_dist)
  
    # Values for the pinhole diameter
    pin_range <- pinhole_range
    pinhole_d <- seq(pin_range[1], pin_range[2], length.out = 50)*1e-6
    N_p <- length(pinhole_d)
  
    # Create data
    pin_df <- tibble(pinhole = rep(pinhole_d, each = N_s),
                   source_dists = rep(source_dist, times = N_p))
  
    pin_df %<>% mutate(
        betas = beta_calc(skim_x, source_dists),
        standard_deviation = theta_p(pinhole, betas, working_dist*1e-3),
        flux = F_o(betas, pinhole),
        FWHM = standard_deviation*2*sqrt(2*log(2))
    )
  
    # Interactive plot of the resolution using plotly
    wide_df <- pin_df %>% select(pinhole, source_dists, FWHM) %>% 
        spread(pinhole, FWHM) %>% mutate(source_dists=NULL)
    FWHM_mat <- as.matrix(wide_df)
  
  
    wide_df2 <- pin_df %>% select(pinhole, source_dists, flux) %>% 
        spread(pinhole, flux) %>% mutate(source_dists=NULL)
    flux_mat <- as.matrix(wide_df2)
  
    if (flux_scaling == "Log")
        flux_mat <- log10(flux_mat)
  
    # Create plot of the resolution
    if (axis_scale == "Linear")
        ax <- list(title = "Pinhole diameter/&mu;m", typ2 = "linear")
    else
        ax <- list(title = "Pinhole diameter/&mu;m", type = "log")
    ay <- list(title = "Source distance/cm")
  
    p1 <- plot_ly(
        x = pinhole_d*1e6,
        y = source_dist*1e2,
        z = FWHM_mat*1e6,
        type = "contour"
    ) %>% colorbar(title = "FWHM/&mu;m") %>%
        layout(xaxis = ax, yaxis = ay)
  
    p2 <- plot_ly(
        x = pinhole_d*1e6,
        y = source_dist*1e2,
        z = flux_mat/flux_base,
        type = "contour"
    ) %>% colorbar(title = "Relative flux") %>%
        layout(xaxis = ax, yaxis = ay)
  
  return(list("plot1" = p1, "plot2" = p2))
}

pinhole_fixed_plot <- function(source_range, pinhole_r, working_dist, skim_x,
                               flux_base, flux_scaling, axis_scale) {
  # Values for the source distance
  s_range <- source_range
  source_dist <- seq(s_range[1], s_range[2], length.out = 50)*1e-2
  N_s <- length(source_dist)
  
  # Values for the working distance
  working_d <- seq(working_dist[1], working_dist[2], length.out = 50)*1e-3
  N_w <- length(working_d)
  
  # Create data
  pin_df <- tibble(working_dists = rep(working_d, each = N_s),
                   source_dists = rep(source_dist, times = N_w))
  
  pin_df %<>% mutate(
    betas = beta_calc(skim_x, source_dists),
    standard_deviation = theta_p(pinhole_r*1e-6, betas, working_dists),
    flux = F_o(betas, pinhole_r*1e-6),
    FWHM = standard_deviation*2*sqrt(2*log(2))
  )
  
  # Interactive plot of the resolution using plotly
  wide_df <- pin_df %>% select(working_dists, source_dists, FWHM) %>% 
    spread(working_dists, FWHM) %>% mutate(source_dists=NULL)
  FWHM_mat <- as.matrix(wide_df)
  
  
  wide_df2 <- pin_df %>% select(working_dists, source_dists, flux) %>% 
    spread(working_dists, flux) %>% mutate(source_dists=NULL)
  flux_mat <- as.matrix(wide_df2)
  
  if (flux_scaling == "Log")
    flux_mat <- log10(flux_mat)
  
  # Create plot of the resolution
  if (axis_scale == "Linear")
    ax <- list(title = "Working distance/mm", typ2 = "linear")
  else
    ax <- list(title = "Working distance/mm", type = "log")
  ay <- list(title = "Source distance/cm")
  
  p1 <- plot_ly(
    x = working_d*1e3,
    y = source_dist*1e2,
    z = FWHM_mat*1e6,
    type = "contour"
  ) %>% colorbar(title = "FWHM/&mu;m") %>%
    layout(xaxis = ax, yaxis = ay)
  
  p2 <- plot_ly(
    x = working_d*1e3,
    y = source_dist*1e2,
    z = flux_mat/flux_base,
    type = "contour"
  ) %>% colorbar(title = "Relative flux") %>%
    layout(xaxis = ax, yaxis = ay)
  
  return(list("plot1" = p1, "plot2" = p2))
}

source_fixed_plot <- function(source_dist, pinhole_r, working_dist, skim_x,
                               flux_base, flux_scaling, axis_scale) {
  # Values for the pinhole diameter
  pin_range <- pinhole_r
  pinhole_d <- seq(pin_range[1], pin_range[2], length.out = 50)*1e-6
  N_p <- length(pinhole_d)
  
  # Values for the working distance
  working_d <- seq(working_dist[1], working_dist[2], length.out = 50)*1e-3
  N_w <- length(working_d)
  
  # Create data
  pin_df <- tibble(working_dists = rep(working_d, each = N_p),
                   pinhole = rep(pinhole_d, times = N_w))
  
  pin_df %<>% mutate(
    betas = beta_calc(skim_x, source_dist*1e-2),
    standard_deviation = theta_p(pinhole, betas, working_dists),
    flux = F_o(betas, pinhole),
    FWHM = standard_deviation*2*sqrt(2*log(2))
  )
  
  # Interactive plot of the resolution using plotly
  wide_df <- pin_df %>% select(working_dists, pinhole, FWHM) %>% 
    spread(working_dists, FWHM) %>% mutate(pinhole=NULL)
  FWHM_mat <- as.matrix(wide_df)
  
  
  wide_df2 <- pin_df %>% select(working_dists, pinhole, flux) %>% 
    spread(working_dists, flux) %>% mutate(pinhole=NULL)
  flux_mat <- as.matrix(wide_df2)
  
  if (flux_scaling == "Log")
    flux_mat <- log10(flux_mat)
  
  # Create plot of the resolution
  if (axis_scale == "Linear")
    ax <- list(title = "Working distance/mm", typ2 = "linear")
  else
    ax <- list(title = "Working distance/mm", type = "log")
  ay <- list(title = "Pinhole diameter/&mu;m")
  
  p1 <- plot_ly(
    x = working_d*1e3,
    y = pinhole_d*1e6,
    z = FWHM_mat*1e6,
    type = "contour"
  ) %>% colorbar(title = "FWHM/&mu;m") %>%
    layout(xaxis = ax, yaxis = ay)
  
  p2 <- plot_ly(
    x = working_d*1e3,
    y = pinhole_d*1e6,
    z = flux_mat/flux_base,
    type = "contour"
  ) %>% colorbar(title = "Relative flux") %>%
    layout(xaxis = ax, yaxis = ay)
  
  return(list("plot1" = p1, "plot2" = p2))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
    # Application title
    titlePanel("SHeM Constrained Optomisation"),
    
    # Sidebar with a slider input for the 
    sidebarLayout(
    sidebarPanel(
    selectInput("fixed_param",
                        "Parameter to fix",
                        choices = c("Working distance", "Pinhole size", "Source distance")),
    
    conditionalPanel(
        "input.fixed_param == 'Working distance'",
        sliderInput("working_dist",
                    "Working distance (mm):",
                    min = 0.1,
                    max = 5,
                    value = 3)
    ),
    conditionalPanel(
        "input.fixed_param != 'Working distance'",
        sliderInput("working_dist2",
                    "Working distance range (mm):",
                    min = 0.1,
                    max = 5,
                    value = c(0.5,3))
    ),
    
    conditionalPanel(
        "input.fixed_param == 'Pinhole size'",
        sliderInput("pinhole_range",
                    "Pinhole diameter (um):",
                    min = 0.01,
                    max = 10,
                    value = 1)
    ),
    conditionalPanel(
        "input.fixed_param != 'Pinhole size'",
        sliderInput("pinhole_range2",
                    "Pinhole range (um):",
                     min = 0.01,
                     max = 10,
                     value = c(0.2, 2))
    ),
    
    conditionalPanel(
        "input.fixed_param == 'Source distance'",
        sliderInput("source_range",
                    "Source distance (cm):",
                    min = 5,
                    max = 100,
                    value = 25)
    ),
    conditionalPanel(
        "input.fixed_param != 'Source distance'",
        sliderInput("source_range2",
                    "Source distance range (cm):",
                    min = 5,
                    max = 100,
                    value = c(10, 30))
    ),
    
    sliderInput("x_skimmer",
                "Skimmer diameter (um):",
                min = 10,
                max = 500,
                value = 100),
    
    radioButtons("axis_scale",
                 "Pinhole axis scaling:",
                 choices = c("Linear", "Log"),
                 selected = "Linear",
                 inline = TRUE),
    radioButtons("flux_scaling",
                 "Relative flux scaling:",
                 choices = c("Linear", "Log"),
                 selected = "Linear",
                 inline = TRUE),
    radioButtons("dist_type",
                 "Actual (beam) working distance or perpendicular distance:",
                 choices = c("Beam", "Perpendicular"),
                 selected = "Beam",
                 inline = TRUE),
    ),
    
    mainPanel(
    plotlyOutput("resolution_plot"),
    verbatimTextOutput("event")
    )
    )
)

# Define server logic required to draw a heat map
server <- function(input, output, session) {


   
    output$resolution_plot <- renderPlotly({
      
        # Value for the skimmer radius
        skim_x <- input$x_skimmer*1e-6
        
        # Normalisation value for the flux
        flux_base <- F_o(beta_calc(skim_x, 0.23), 0.38e-6)
        
        if (input$dist_type == "Beam") 
            factor <- 1
        else
            factpr <- sqrt(2)
        # Which parameter is fixed
        if (input$fixed_param == "Working distance") {
            ps <- working_fixed_plot(input$source_range2, input$pinhole_range2, input$working_dist, 
                                     skim_x, flux_base, input$flux_scaling, input$axis_scale)
        } else if (input$fixed_param == "Pinhole size") {
          ps <- pinhole_fixed_plot(input$source_range2, input$pinhole_range, input$working_dist2, 
                                   skim_x, flux_base, input$flux_scaling, input$axis_scale)
        } else { # This case is the source distance
          ps <- source_fixed_plot(input$source_range, input$pinhole_range2, input$working_dist2, 
                                   skim_x, flux_base, input$flux_scaling, input$axis_scale)
        }

        # Draw two plots
        p <- subplot(ps[["plot1"]], ps[["plot2"]], shareY = TRUE, titleX = TRUE)
        # Hack the damn thing to have sub plot titles
        p %>% layout(annotations = list(
          list(x = 0.2, 
               y = 1.08, 
               text = "Resolution", 
               showarrow = FALSE, 
               xref='paper', 
               yref='paper', 
               font = list(size = 16)
          ),
          list(x = 0.8, 
               y = 1.08, 
               text = "Relative flux", 
               showarrow = FALSE, 
               xref='paper', 
               yref='paper', 
               font = list(size = 16)
          ))
        )
    })
   
    output$event <- renderPrint({
        d <- event_data("plotly_click")
        if (is.null(d)) "Click on a point!" else d
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

