#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyalert)
#####
data <- mtcars
###
#Function to convert [Maltose] into OD reading
DNSA.maltose <- function(x,sd.reader = 0.03){
  fas <- numeric(length = length(x))
  
  for(i in 1:length(x)){
    if(is.na(x[i]) == T){
      absorbance <- NA
    } else if(x[i] < 125){
      
      absorbance <- 0.065
      
    }else if(x[i] > 6000){
      
      absorbance <- 3.011 + (1.875e-05)*x[i]  
      
    }else if(x[i] <= 6000 & x[i] >= 2500){
      
      absorbance <- log10(x[i]*.3395 - 748.75)
      
    }else{
      
      absorbance <- x[i]*(8.298e-04)-(3.794e-02)
      
    }
    er.adj <- rnorm(n = 1, mean = 0, sd = sd.reader)
    fin.abs <- absorbance + er.adj
    fas[i] <- fin.abs
    
  }
  
  return(fas)
  
}

#Enzyme efficiency curve for pH optimization based on
#Gangadharan et al 2009 Appl Biochem Biotech
#Returns % efficiency of amylase
ph.adj <- function(x){
  mod.coefs <- c(  2.403317e+03,
                   -2.943847e+03,
                   1.546476e+03,
                   -4.552539e+02,
                   8.221168e+01,
                   -9.334084e+00,
                   6.512833e-01,
                   -2.555783e-02,
                   4.322819e-04
  )
  exponents <- c(0:8)
  kme <-numeric(length(x))
  for(i in 1:length(x)){
    km.eff <- sum(mod.coefs * x[i]^exponents)  
    kme[i]<-km.eff
  }
  return(kme)
}

#liberation of reducing ends by alpha amylase from starch concentration.
#Values based on Heitman et al 1997 Enzyme Microb. Technol.
red.liber <- function(st.conc, Km = 2, Vm = 3e3, time =3, pH =7){
  #units#
  #st.conc = ug/ml
  #Km = mg/ml
  #Vm = umol/(min*mg)
  #time = min
  #values based on Heitman et al 1997 Enzyme Microb. Technol.
  st.conc.gl <- st.conc*.001  #this puts conc in mg/ml
  #Adjusting Km and Vm for pH
  Km.adj <- Km / ph.adj(pH)
  Vm.adj <- Vm *sqrt(ph.adj(pH))
  
  v = (Vm * st.conc.gl) / (Km.adj + st.conc.gl)
  
  lib.val = time*v*342.3e-3
  
  return(lib.val)
  
  #km values calculated from these data will be in ug/ml, not the input value of mg/ml
}

plate.read <- function(substrate.data, substrate, pH,
                       pH.data = NULL, polym.fac = 15,
                       sd.reader = 0.01, sd.experimenter.pct = 0.05,
                       amylase = FALSE, Km = 4, Vm = 3e3, time = 3){
  
  # polym.fac is a proxy for how long the starch chain is. It is
  # used to determine how many reducing ends there are per mass of
  # starch relative to the number of reducing ends per mass maltose. 
  # i.e. if polym.fac = 15, 15ug of starch will cause the same DNSA 
  # color change as 1ug of maltose.
  
  d.plate.long <- reshape(data = substrate.data,
                          idvar = "Row",
                          varying = list(2:13),
                          v.names = "Conc",
                          timevar = "Column",
                          direction = "long")
  
  
  d.plate.long$Row.Num <- as.numeric(d.plate.long$Row)*-1
  
  ## add experimenter error to concentrations##
  ## specified as sd.experimenter.pct ==> standard deviation of % error!##
  d.plate.long$Conc.Exp <- numeric(length(d.plate.long$Conc))
  for(i in 1:length(d.plate.long[,1])){
    d.plate.long$Conc.Exp[i] <- rnorm(n=1,
                                      mean = d.plate.long$Conc[i],
                                      sd = d.plate.long$Conc[i]*sd.experimenter.pct
    )
  }
  
  if(substrate == "maltose"){
    

    d.plate.long$Reading <- DNSA.maltose(x = d.plate.long$Conc.Exp,
                                         sd.reader = sd.reader)
    
  }else if(substrate == "starch"){

    if(amylase == T){
      #If amylase is present in reaction
      if(is.null(pH.data) == F){
        #If pH values are given in a 96 well plate format
        d.pH.long <- reshape(data = pH.data,
                             idvar = "Row",
                             varying = list(2:13),
                             v.names = "pH",
                             timevar = "Column",
                             direction = "long") 
        
        d.plate.long$pH <- d.pH.long$pH
        
        d.plate.long$Maltose.Conc <- red.liber(st.conc = d.plate.long$Conc.Exp,
                                               Km = Km, Vm = Vm, time = time,
                                               pH = d.plate.long$pH)
        
        
        
      }else{
        
        # pH is constant value for all reactions
        d.plate.long$Maltose.Conc <- red.liber(st.conc = d.plate.long$Conc.Exp,
                                               Km = Km, Vm = Vm, time = time,
                                               pH = pH)
      }
      #Amount of reducing ends expressed in terms of equivalent maltose conc
      # IF AMYLASE IS ADDED
      d.plate.long$Red.Conc <- (d.plate.long$Conc.Exp/15) + d.plate.long$Maltose.Conc
      
    }else{
      #Amount of reducing ends expressed in terms of equivalent maltose conc
      # IF NO AMYLASE
      d.plate.long$Red.Conc <- (d.plate.long$Conc.Exp/15)
      
    }
    # Plate reader values
    d.plate.long$Reading <- DNSA.maltose(x = d.plate.long$Red.Conc,
                                         sd.reader = sd.reader)
  }
  
  
  yrPal <- colorRampPalette(c("yellow","darkred"))
  
  est.col.vec <- c(d.plate.long$Reading,0,3.3)
  est.col.vec.val <- yrPal(1000)[as.numeric(cut(est.col.vec, breaks = 1000))]
  
  d.plate.long$colval <- est.col.vec.val[c(1:96)]

  
  d.plate.long.sel <- d.plate.long[,c("Row", "Column","Reading")]
  
  d.plate.wide <- reshape(data = d.plate.long.sel,
                          v.names = "Reading",
                          idvar = "Row",
                          timevar = "Column",
                          direction = "wide")
  
  names(d.plate.wide) <- c("Row",as.character(c(1:12)))
  row.names(d.plate.wide) <- NULL
  dfs <- list(d.plate.long, d.plate.wide)
  return(dfs)
  
}

template.96 <- as.data.frame(matrix(nrow = 8, ncol =13))
names(template.96) <- c("Row",as.character(c(1:12))) 
template.96$Row <- c("A","B","C","D","E","F","G","H")


# Define UI for application that draws a histogram
ui <- fluidPage(
   useShinyalert(),
   # Application title
   titlePanel("DNSA Reaction Plate Reader Simulation"),
   
   # Sidebar with a slider input for number of bins 
   verticalLayout(
      
        HTML("<p>This tool only accepts .csv files in a specific format. Please download a 96-well plate template spreadsheet below.</p>"),

        
        downloadButton(outputId = "template",
                       label = "Download 96-well plate template"),
        
        br(),
        hr(),
        br(),
        
        HTML("<p>Next, upload a filled in 96-well plate template. Each cell should specify the solute concentration in the corresponding well. Not all cells need be filled. All concentrations are in ug/ml (= mg/l). All solute concentrations must be < 50,000 ug/ml.</p></p>"),

        
        fileInput(inputId = "conc.data",
                     label = "Choose [Solute] .csv file",
                     multiple = FALSE,
                     accept = ".csv"),
         
         br(),
         HTML("<p>Specify solute type</p>"),
         selectInput(inputId = "solute",
                     label = "Solute",
                     choices = c("maltose","starch")),
         
         br(),
         HTML("<p>Will you be digesting with amylase?</p>"),
        selectInput(inputId = "amylaseTF",
                    label = "Amylase?",
                    choices = c("Yes","No")),
        
         
         HTML("<p>If you are digesting the solute with amylase you must specify the pH. You may specify pH by uploading a file with pH values in 96 well plate or set constant pH using slider. Be sure to select appropriate option. If you select -Variable- you must upload a 96 well plate with pH values specifying pH values for each well that contains solute.</p>"),

                br(),
        
        selectInput(inputId = "ph.sel",
                    label = "pH specification method",
                    choices = c("Constant","Variable")),
        
        sliderInput(inputId = "ph.const",label = "Constant pH level",
                      min = 4, max = 10, value=7, step = 0.1),
        
        fileInput(inputId = "PH.DATA",
                  label = "Choose pH .csv file",
                  multiple = FALSE,
                  accept = ".csv"),
        hr(),
        actionButton(inputId = "run.sim",
                     label = "Run Simulation"),
        hr(),
        
        textInput(inputId = "plate.name",
                  label = "Plate Reading File Name--must be filled for download of OD readings to work!"),
        downloadButton(outputId = "plate.reading",
                       label = "Download OD Readings"),
        
  
     
      # Show a plate
      mainPanel(
         plotOutput("Plate")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
 
  vals <- reactiveValues(data.out = NULL,
                         the.data = NULL,
                         the.ph.data = NULL)
  
   output$template <- downloadHandler(
     
     filename = function(){
        "template96well.csv"
      },
     content = function(file){
       write.csv(x = template.96,
                 na = "",
                 file = file,
                 row.names = FALSE)
     }
   )

   observeEvent(eventExpr = input$conc.data,{
     req(input$conc.data)
     vals$the.data <- read.csv(file = input$conc.data$datapath, header = TRUE)
    
   })
   
   observeEvent(eventExpr = input$PH.DATA,{
     req(input$PH.DATA)
     vals$the.ph.data <- read.csv(file = input$PH.DATA$datapath, header = TRUE)
    
   })
  
  observeEvent(eventExpr = input$run.sim, {
      if(is.null(vals$the.data) == TRUE){
        shinyalert(title = "No Solute Data!",
                   text = "You forgot the reaction substrate...",
                   type = "error")

      }else{
        
        if(any(vals$the.data[,c(2:13)] < 0, na.rm=TRUE) == TRUE |
           any(vals$the.data[,c(2:13)] > 5e4, na.rm=TRUE) == TRUE){
          
          shinyalert(title = "Solute Concentration Error",
                     text = "Solute concentration must be non-negative and <50,000",
                     type = "error")
          
        }else{
          
          
          if(input$ph.sel == "Variable"){
            if(is.null(vals$the.ph.data) == TRUE){
              shinyalert(title = "No pH Data!",
                         text = "If 'variable' pH is specified, upload pH data",
                         type = "error")
              
            }else{
              
              phd <- vals$the.ph.data
              if(any(phd[,c(2:13)] < 4, na.rm=TRUE) == TRUE |
                 any(phd[,c(2:13)] > 10, na.rm=TRUE) == TRUE){
                phd <- NULL
                shinyalert(title = "pH out of bounds",
                           text = "pH must be 4 < pH <10 ",
                           type = "error")
              }
              
            }
            
          }else if(input$ph.sel == "Constant"){
            phd <- NULL
          }
          if(input$amylaseTF == "Yes"){
            amTF <- TRUE
          }else if(input$amylaseTF == "No"){
            amTF <- FALSE
          }
          dd <- plate.read(substrate.data = vals$the.data,
                           substrate = input$solute, 
                           pH=input$ph.const,
                           pH.data = phd, polym.fac = 15,
                           sd.reader = 0.01, sd.experimenter.pct = 0.05,
                           amylase = amTF, Km = 4, Vm = 3e3, time = 3)
          
          output$Plate<- renderPlot({plot(data = dd[[1]],
                                          Row.Num ~ Column,
                                          pch = 21,
                                          cex = 4,
                                          bg = dd[[1]]$colval,
                                          ylim = c(-8.5,-.5),
                                          yaxt = "n", ylab="", xaxt = "n", xlab = "")
            axis(2, at = c(-8:-1),
                 labels = c("H","G","F","E","D","C","B","A"), las = 2)
            axis(1, at = c(1:12),
                 labels = c(1:12), las = 1)
          })
          
          vals$data.out <- dd 
          
          
        } 

        
      }

  
   })
  
   # observeEvent(eventExpr = output$plate.reading,{
   #   if(input$plate.name == ""){
   #     
   #     shinyalert(title = "Error",
   #                text = "Provide name for downloaded data",
   #                type ="error")
   #     
   #   }else if(is.null(vals$data.out) == TRUE){
   #     
   #     shinyalert(title = "Error",
   #                text = "Run simulation first!",
   #                type = "error")
   #     
   #   }else{
        output$plate.reading <- downloadHandler(
         filename = function(){
           paste(input$plate.name, ".csv", sep = "")
         
       },
       content = function(file){
         
           write.csv(x = vals$data.out[[2]],
                     na = "",
                     file = file,
                     row.names = FALSE)
       })
   #   }
   # })
  
   
   
}

# Run the application 
shinyApp(ui = ui, server = server)

