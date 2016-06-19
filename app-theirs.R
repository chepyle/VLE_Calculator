#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(reshape2)
library(ggplot2)

pkgLoad <- function(x)

{

  if (!require(x,character.only = TRUE))

  {

    install.packages(x,dep=TRUE,repo="http://cran.rstudio.com")

    if(!require(x,character.only = TRUE)) stop("Package not found")

  }

}
pkgLoad('rbokkeh')

source('UNIFAC_calc.R')

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("VLE-UNIFAC"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("P",
                     "Pressure(bar):",
                     min = 0.001,
                     max = 1,
                     value = 1,
                     step= 0.1),
         selectInput('mol1', 'Species 1', Pvap_coef$Name[!is.na(Pvap_coef$Equation.Form)],selected= 'water'),
         selectInput('mol2', 'Species 2', Pvap_coef$Name[!is.na(Pvap_coef$Equation.Form)],selected = 'methanol'),
         downloadButton('downloadData','Download')
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(type = "tabs", 
                    tabPanel("Txy", plotOutput("distPlot")), 
                    tabPanel("xy", plotOutput("xyplot")),
                    tabPanel("Txy_interactive", rbokehOutput("txyBokeh"))
        )
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {

  vle_data<-reactive({

    return(Txy(input$P,c(input$mol1,input$mol2)))
  })
  
   output$distPlot <- renderPlot({
      txy<-melt(data.frame(X=vle_data()$x1,
                           Y=vle_data()$y1,Temperature=vle_data()$Temperature),
                id.vars = 'Temperature')
      p<-ggplot(txy,aes(x=value,y=Temperature,color=variable))+geom_line()+geom_point()
      return(p)
   })
   
   output$xyplot <- renderPlot({
     xy<-melt(vle_data(),id.vars=c('x1','y1'),value.name=('Temperature'))
     p<-ggplot(xy,aes(x=x1,y=y1,color=Temperature))+geom_line()+geom_point()+geom_abline(slope=1,intercept=0)
     return(p)
   })
   
   output$txyBokeh<-renderRbokeh({
     tools <- c("pan", "wheel_zoom", "box_zoom", "box_select", "resize", "reset")
     txy<-signif(vle_data(),digits=4)
     
     figure(tools=tools, title = paste0("Txy Diagram of ",input$mol1," & ", input$mol2)) %>%
       ly_points(x1,Temperature,data=txy,line_width=2,color='blue',
                 hover=c(x1,y1,Temperature)) %>%
       ly_lines(x1,Temperature,data=txy,color='blue') %>%
       ly_points(y1,Temperature,data=txy,line_width=2,color='red',
                 hover=c(x1,y1,Temperature)) %>%
       ly_lines(y1,Temperature,data=txy,color='red')
   })
   
   output$downloadData <- downloadHandler(
     filename = function() { 
       paste('Txy_',input$P,'bar_',input$mol1,'_',input$mol2, '.csv', sep='') 
     },
     content = function(file) {
       write.csv(vle_data(), file)
     }
   )
})

# Run the application 
shinyApp(ui = ui, server = server)

