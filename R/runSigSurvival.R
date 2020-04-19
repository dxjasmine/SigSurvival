library(shiny)
library(shinydashboard)

mydb <- read.csv("my datasets",header = T,sep = "\t",stringsAsFactors = F)

data_info <- as.matrix(mydb)
discard <- c("platform","Title","Source")
data_info <- data_info[,!colnames(data_info) %in% discard]

if(interactive()){
  shinyApp(
    ui <- dashboardPage(
      
      dashboardHeader(title = "SigSurvival: validation"),
      dashboardSidebar(width = "500px",
                       br(),
                       br(),
                       # Input: gene signature ----
                       
                       textInput(inputId = "gene_id",
                                 label = "Gene signature:",
                                 value = "ATP1B1,TRIM14,FAM64A"),
                       # Input: Selector for choosing tissue ----
                       radioButtons(inputId = "Tissue",
                                    label = "Choose tissue type:",
                                    choices = c("Lung(NSCLC)","Other" ),inline = T),
                       radioButtons(inputId = "endp",
                                    label = "Choose end point:",
                                    choices = c("Overall Survival", "Disease Free Survival", "Recurrance"),inline = T),
                       # Input: Selector for choosing dataset ----
                       checkboxGroupInput(inputId = "dataset",
                                          label = "Choose dataset:",
                                          choices = data_info[,2],inline = T),
                       
                       
                       
                       tableOutput("data1")),
      dashboardBody()
    ),
    
    server = function(input, output, session){
      output$data1 <- renderTable({
        data_info
      })
      
      
      
    }
  )
}
