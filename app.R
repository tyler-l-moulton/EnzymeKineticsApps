#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Enzyme Kinetic Demo"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput(inputId = "km",
                        "Km Value",
                        min = 1e-4,
                        max = 1,
                        value = .5,
                        step = 1e-4),
            sliderInput(inputId = "Vmax",
                        label = "Vmax Value",
                        min = 1e-6,
                        max = 1e-3,
                        value = 3e-5,
                        step = 1e-6),
            
            #michaelis menten xy limits
            numericInput(inputId = "mm.Xmax",
                      label = "Michaelis-Menten X max",
                      value = 5,
                      min = 1e-20,
                      max = 1e20,
                      step = 1e-2),
            
            numericInput(inputId = "mm.Ymax",
                         label = "Michaelis-Menten Y max",
                         value = 5e-5,
                         min = 1e-20,
                         max = 1e20,
                         step = 1e-6),
            
            #lineweaver burk xy limits
            numericInput(inputId = "lb.Xmax",
                         label = "Lineweaver-Burk X max",
                         value = 10,
                         min = 1e-20,
                         max = 1e20,
                         step = 1e-2),
            
            numericInput(inputId = "lb.Ymax",
                         label = "Lineweaver-Burk Y max",
                         value = 4e5,
                         min = 1e-20,
                         max = 1e20,
                         step = 100),
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("mm.plot"),
           plotOutput("lb.plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$mm.plot <- renderPlot({


        curve(expr = (input$Vmax * x) / (input$km + x),
              lwd =2 , col = "blue",
              xlab = "Substrate Concentration [S]", 
              ylab = "Initial Reaction Velocity (v)",
              main = "Michaelis-Menten Plot",
              xlim = c(0, input$mm.Xmax),
              ylim = c(0, input$mm.Ymax))
        
        #Horizontal line
        lines(x = c(-4,input$km), 
              y = c(input$Vmax/2,input$Vmax/2),
              lty = 2)
        
        text(x = input$km/2, y = 0.6*input$Vmax,
             labels = expression('V'['max']*'/2'))
        
        #Vertical line
        lines(x = c(input$km, input$km),
              y = c(-1, input$Vmax/2), 
              col = "red", lty = 2)
        text(x = input$km * 1.4, y = 0.25*input$Vmax,
             labels = expression('K'['m']), col = "red")

    })
    
    output$lb.plot <- renderPlot({
        # generate bins based on input$bins from ui.R
        lb.Xmin <- -input$lb.Xmax*0.3
        lb.Ymin <- -input$lb.Ymax*0.3
        
        
        curve(expr = (input$km*x +1)*(1/input$Vmax) ,
              lwd =2 , col = "blue",
              main = "Lineweaver-Burk Plot",
              xlab = "1 / [S]", 
              ylab = "1 / v",
              xlim = c(lb.Xmin, input$lb.Xmax),
              ylim = c(lb.Ymin, input$lb.Ymax))
        abline(0,0)
        abline(v=0)
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
