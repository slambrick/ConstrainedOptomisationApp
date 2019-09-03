# Sam Lambrick, 2019
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

# Define physical constants and other fixed parameters, SI units used throughout
lambda <- 5.63e-11     # Room temperature helium wavelength
B_peak <- 1e20         # Optomised peak brightness (not calculated yet)
gamma <- pi^2*B_peak/4 # Pinhole flux
skim <- 100e-6         # Skimmer diameter

# Define mathematical functions

# Standard deviation for a pinhole
theta_p <- function(d, beta, f) {
    pinhole <- d/(2*sqrt(3))
    source <- beta*f/sqrt(3)
    diffraction <- 0.42*lambda*f/d
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
F_o <- function(sig, f) 3*gamma*(3*sig^4/f^2 - (0.42*lambda)^2)

# Virtual source size from skimmer-pinhole distance
beta_calc <- function(dist) atan(skim/dist)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("SHeM Constrained Optomisation"),
   
   # Sidebar with a slider input for the 
   sidebarLayout(
      sidebarPanel(
         sliderInput("working_dist",
                     "Working distance (mm):",
                     min = 0.5,
                     max = 5,
                     value = 3)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("resolution_plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$resolution_plot <- renderPlot({
        # Values for the source distance
        source_dist <- seq(10, 30, by = 0.2)*1e-2
        N_s <- length(source_dist)
        
        # Values for the pinhole diameter
        pinhole_d <- seq(0.1, 1, by = 0.02)*1e-6
        N_p <- length(pinhole_d)
        
        # Create data
        pin_df <- tibble(pinhole = rep(pinhole_d, each = N_s),
                         source_dists = rep(source_dist, times = N_p))
        
        pin_df %<>% mutate(
            betas = beta_calc(source_dists),
            standard_deviation = theta_p(pinhole, betas, input$working_dist*1e-3),
            flux = F_o(standard_deviation, input$working_dist*1e-3),
            FWHM = standard_deviation*2*sqrt(2*log(2))
        )
        
        # Plot of the resolution
        pin_df %>% ggplot(aes(x = pinhole*1e6, y = source_dists*1e2, z = FWHM*1e6)) +
            geom_raster(aes(fill = FWHM*1e6), interpolate = TRUE) +
            geom_contour(colour = "white", linemitre = 20) +
            scale_fill_continuous(type = "viridis") + 
            labs(x = expression(paste("Pinhole diameter/", mu, "m", sep="")),
                 y = "Source distance/cm", 
                 fill = expression(paste("FWHM/", mu, "m", sep=""))) +
            theme_classic() +
            theme(legend.key.height = unit(2, "cm"))
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

