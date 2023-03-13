library(VGAM)
library(stats4)
library(splines)
library(gmp)
library(hms)

shinyUI(fluidPage(
  titlePanel("Required sample size to estimate the location parameter in Laplace distribution"),
  
  sidebarLayout(
    sidebarPanel(
      h1("Input Data"),
      
      sliderInput(inputId="f",
                  label="Precision (f):",
                  min=0.1,
                  max=0.3,
                  value=0.2,
                  step=0.05
      ),
      helpText("Precision (f) is in [0.1, 0.3]."),
      
      sliderInput(inputId="c",
                  label="Confidence Level (c) :",
                  min=0.85,
                  max=0.99,
                  value=0.95,
                  step=0.05
      ),
      helpText("Confidence level (c) is in [0.85, 0.99]."),
      
      
      submitButton("Update View")
    ),
    
    mainPanel(
      h1("Summary"),
      tableOutput("values")
    )
  )
))
