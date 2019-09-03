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
library(plotly)

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
F_o <- function(beta, d) gamma*d^2*beta^2

# Virtual source size from skimmer-pinhole distance
beta_calc <- function(dist) atan(skim/dist)

# Normalisation value for the flux
flux_base <- F_o(beta_calc(0.23), 0.38e-6)

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
                        value = 3),
            sliderInput("pinhole_range",
                        "Pinhole range (um):",
                        min = 0.01,
                        max = 10,
                        value = c(0.1, 1)),
            sliderInput("source_range",
                        "Source distance range (cm):",
                        min = 5,
                        max = 100,
                        value = c(10, 30)),
            radioButtons("axis_scale",
                         "Pinhole axis scaling:",
                         choices = c("Linear", "Log"),
                         selected = "Linear",
                         inline = TRUE),
            radioButtons("flux_scaling",
                         "Relative flux scaling:",
                         choices = c("Linear", "Log"),
                         selected = "Linear",
                         inline = TRUE)
        ),
      
        # Show a plot of the generated distribution
        mainPanel(
            plotlyOutput("resolution_plot"),
            verbatimTextOutput("event")
        )
    )
)

# Define server logic required to draw a heat map
server <- function(input, output) {


   
    output$resolution_plot <- renderPlotly({
        # Values for the source distance
        s_range <- input$source_range
        source_dist <- seq(s_range[1], s_range[2], length.out = 100)*1e-2
        N_s <- length(source_dist)
        
        # Values for the pinhole diameter
        pin_range <- input$pinhole_range
        pinhole_d <- seq(pin_range[1], pin_range[2], length.out = 100)*1e-6
        N_p <- length(pinhole_d)
        
        # Create data
        pin_df <- tibble(pinhole = rep(pinhole_d, each = N_s),
                         source_dists = rep(source_dist, times = N_p))
        
        pin_df %<>% mutate(
            betas = beta_calc(source_dists),
            standard_deviation = theta_p(pinhole, betas, input$working_dist*1e-3),
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
        
        if (input$flux_scaling == "Log")
            flux_mat <- log10(flux_mat)
        
        # Create plot of the resolution
        if (input$axis_scale == "Linear")
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
        
        # Draw two plots
        p <- subplot(p1, p2, shareY = TRUE, titleX = TRUE)
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
            )
        ))
    })
   
    output$event <- renderPrint({
        d <- event_data("plotly_click")
        if (is.null(d)) "Click on a point!" else d
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

