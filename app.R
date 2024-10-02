rm(list = ls())

library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinydashboard)

library(tidyverse)
library(igraph)
library(ggnetwork)
library(ggsci)

library(this.path)
setwd(this.dir())

source('plot.R')

# about
about = "DDRnet utilises the xgboost model to explore genetic interactions with 
and amongst DDR genes, offering a complementary approach to CRISPR screen 
analyses. The xgboost model identifies interactions with a target gene, such as 
ATM, by predicting its z-scores based on those of other genes across all screens
within the DDRcs. Through a training process, the xgboost model autonomously 
constructs an ensemble of regression trees with a select few important genes, 
which are considered to interact with the target gene (e.g. ATM). DDRnet applies
this procedure to all genes in the human genome to infer genetic interactions 
with DDR genes.<br><br>

The DDRnet web app offers two modes of network visualization: (1) gene query and
(2) cluster. In gene query mode, users select a list of genes of interest and 
query all the genes interacting with these input genes. Users have the 
flexibility to adjust the network confidence by setting frequency and importance
threshold. Higher frequency and importance yield more stringent analyses but 
result in fewer interacting genes. Additionally, users can choose from various 
display options to customize the network visualisation. In cluster mode, users 
can explore subsets of Pendragon genes, which have been pre-clustered offline 
(frequency \u2265 4 and importance \u2265 0.08). Adjusting the frequency and importance 
threshold results in different display of gene interactions within the clusters.
<br><br>

Since we filtered out some genes to maintain high consistency of the data across
multiple CRISPR screens, DDRnet does not encompass a complete genetic network of
all gene pairs. If you cannot locate your genes of interest in the networks or 
if they only interact with few others, it is likely that these genes were not 
included in our analysis. The list of genes utilised in the DDRnet can be found 
in the “Genes” tab.
"

load('DDRnetData.RData')

# Define UI for app that draws a gene interaction network ----
ui <- navbarPage(
  theme = shinytheme("sandstone"),
  title = "DDRnet",
  header = tagList(
    useShinydashboard()
  ), 
  
  ## something about the app
  tabPanel(title = 'Home', 
           fluidRow(
             column(width = 2), 
             column(width = 8, 
                    HTML("<h1 style='text-align: center;'>DDRnet: A xgboost Method to Infer Genetic Interactions on DDR CRISPR Screen Data</h1>"),
                    
                    br(),
                    br(),
                    
                    HTML(about)), 
             column(width = 2)
           )
  ), 
  
  ## plot network and cluster of genes ----
  tabPanel(title = 'Queries', 
           sidebarLayout(
             
             ### Sidebar panel for inputs ----
             sidebarPanel(
               width = 2,
               
               #### general plot options ----
               tags$div(style = 'text-align: center',
                        tags$h3('Plot options')),
               
               ##### Input: Slider for frequency, display as confidence ----
               sliderInput(inputId = "freq",
                           label = "Frequency",
                           min = 2,
                           max = 5,
                           value = 5),
               
               ##### Input: Slider for importance score, display as importance ----
               sliderInput(inputId = "impt",
                           label = "Importance",
                           min = 0,
                           max = 0.2,
                           step = 0.001,
                           value = 0),
               
               ##### Input: Slider for recurrence, display as Expansion level ----
               sliderInput(inputId = "recr",
                           label = "Expansion level",
                           min = 1,
                           max = 3,
                           step = 1,
                           value = 1),
               
               ##### Input: checkbox group for plot options, display as plot options ----
               checkboxGroupInput(inputId = 'plotOpts',
                                  label = 'Display options',
                                  choices = list('Only in GO:DDR' = 'GO:DDR',
                                                 'Only in Pendragon' = 'Pendragon')),
               
               ##### Input: selectised inputs for queried genes, display as genes ----
               # selectizeInput(inputId = 'genes',
               #                label = 'Queried genes',
               #                selected = NULL,
               #                multiple = TRUE,
               #                options = NULL,
               #                choices = NULL),
               
               selectizeInput(inputId = 'genes',
                              label = 'Queried genes',
                              selected = NULL,
                              multiple = TRUE,
                              options = NULL,
                              choices = unique(c(netTab$gene_1,
                                                 netTab$gene_2))),
               
               ##### Input: action button to query gene interactions, display as plot ----
               tags$div(style = 'text-align: center',
                        actionButton(inputId = 'plotQueries',
                                     label = 'Plot'))
             ),
             
             #### Main panel for displaying outputs ----
             mainPanel(
               width = 10,
               style = "height: calc(100vh - 80px);", # adjust height automatically
               
               ##### Output: network ----
               plotOutput(outputId = "plot",
                          height = '100%')
             )
           )
  ), 
  
  ## plot pre-clusterd Pendragon genes
  tabPanel(title = 'Clusters', 
           sidebarLayout(
             
             ### Sidebar panel for inputs ----
             sidebarPanel(
               width = 2,
               
               #### general plot options ----
               tags$div(style = 'text-align: center',
                        tags$h3('Plot options')),
               
               ##### Input: Slider for frequency, display as confidence ----
               sliderInput(inputId = "freqClust",
                           label = "Frequency",
                           min = 4,
                           max = 5, 
                           step = 1,
                           value = 4),
               
               ##### Input: Slider for importance score, display as importance ----
               sliderInput(inputId = "imptClust",
                           label = "Importance",
                           min = 0.08,
                           max = 0.2,
                           step = 0.001,
                           value = 0.08),
               
               ##### Input: Slider for recurrence, display as Expansion level ----
               # sliderInput(inputId = "recrClust",
               #             label = "Expansion level",
               #             min = 1,
               #             max = 3,
               #             step = 1,
               #             value = 1),
               
               ##### Input: checkbox group for plot options, display as plot options ----
               checkboxGroupInput(inputId = 'plotOptsClust',
                                  label = 'Display options',
                                  choices = list('Only in GO:DDR' = 'GO:DDR')),
               
               ##### Input: select input for cluster plot, display as cluster number ----
               selectInput(inputId = 'cluster',
                           label = 'Cluster ID',
                           choices = resultsFreq4$queriedSet),
               
               ##### Input: action button to gene clusters, display as plot ----
               tags$div(style = 'text-align: center',
                        actionButton(inputId = 'plotCluster',
                                     label = 'Plot')
               )
             ),
             
             #### Main panel for displaying outputs ----
             mainPanel(
               width = 10,
               style = "height: calc(100vh - 80px);", # adjust height automatically
               
               ##### Output: network ----
               plotOutput(outputId = "plotClust",
                          height = '100%')
             )
           )
  ),
  
  ## list of genes used to build the xgboost models ----
  tabPanel(title = 'Genes', 
           
           ### list of Pendragon genes ----
           tags$div(style = 'text-align: center',
                    tags$h3('Pendragon genes')), 
           
           textOutput(outputId = 'Pendragon'),
           
           br(),
           br(),
           
           ### list of other genes ----
           tags$div(style = 'text-align: center',
                    tags$h3('Other genes')), 
           
           textOutput(outputId = 'otherGenes')
  )
)

# Define server for app that draws a gene interaction network ----
server <- function(input, output, session) {
  # Reactively generate plot based on selected button
  observeEvent(input$plotQueries,
               {
                 isolatedInput = isolate(input)
                 
                 net = getNet(netTab, geneList, isolatedInput)
                 
                 # plot network
                 output$plot = renderPlot({
                   plotNet(net, isolatedInput)
                 })
               })
  
  observeEvent(input$plotCluster,
               {
                 isolatedInput = isolate(input)
                 
                 net = getClust(clustFreq4List, isolatedInput)
                 
                 # plot network
                 output$plotClust = renderPlot({
                   plotClust(net, isolatedInput)
                 })
               })
  
  # print out list of genes included in DDRnet
  output$Pendragon = renderText({
    paste(geneList[[11]][geneList[[11]] %in% geneXgboost], 
          sep = ', ')
  })
  
  output$otherGenes = renderText({
    paste(geneXgboost[!(geneXgboost %in% geneList[[11]])], 
          sep = ', ')
  })
}


shinyApp(ui = ui, server = server)
