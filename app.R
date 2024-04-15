rm(list = ls())

library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinydashboard)

library(tidyverse)
library(igraph)
library(graphlayouts)
library(ggnetwork)
library(ggsci)

# setwd('~/OneDrive - CRUK Cambridge Institute/DDRcs_exploratoryAnalysis/src/03_permImptPendragonDDRcsZscore/')
# source('xgbCvFeatureImportance.R')
# source('xgbProcessNet.R')
source('plot.R')

# Load clusters on network ----
load('netClusterAllShinyApp.RData')

netFreq = list()
netFreq[[3]] = netFreq3
netFreq[[4]] = netFreq4
netFreq[[5]] = netFreq5

netAll = list()
netAll[[3]] = netAllFreq3
netAll[[4]] = netAllFreq4
netAll[[5]] = netAllFreq5

clustFreq = list()
clustFreq[[3]] = clustFreq3
clustFreq[[4]] = clustFreq4
clustFreq[[5]] = clustFreq5

igraphFreq = list()
igraphFreq[[3]] = igraphFreq3
igraphFreq[[4]] = igraphFreq4
igraphFreq[[5]] = igraphFreq5

clustStat = list()
clustStat[[3]] = clustStatFreq3
clustStat[[4]] = clustStatFreq4
clustStat[[5]] = clustStatFreq5

# about
about = "DDRnet utilises xgboost model to explore genetic interactions with and amongst DDR genes, 
offering a complementary approach to CRISPR screen analyses. The xgboost model identifies interactions 
with a target gene, such as ATM, by predicting its z-scores based on those of other genes across all 
screens within the DDRcs. Through a training process, the xgboost model autonomously constructs an 
ensemble of regression trees with a select few important genes, which are considered to interact with 
the target gene (e.g. ATM). DDRnet applies this procedure to all genes in the Pendragon to infer genetic 
interactions with DDR genes.<br><br>The DDRnet web app offers two modes of network visualization: 
(1) gene query and (2) cluster. In gene query mode, users select a list of genes of interest and query
all the genes interacting with these input genes. Users have the flexibility to adjust the network 
confidence by setting frequency and importance threshold. Higher frequency and importance yield more
stringent analyses but result in fewer interacting genes. Additionally, users can choose from various
display options to customize the network visualisation. In cluster mode, users can explore subsets of
Pendragon genes, which have been pre-clustered offline. Adjusting the frequency threshold results in 
display of different gene clusters. However, changing the importance threshold does not affect the cluster
display.<br><br>Since we only focused on finding interactions with DDR genes in the Pendragon, DDRnet 
does not encompass a complete genetic network of all gene pairs. If you cannot locate your genes of 
interest in the networks or if they only interact with few others, it is likely that these genes were
not included in our analysis. The list of genes utilised in the DDRnet can be found in the “Genes” tab."

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
  tabPanel(title = 'Network', 
           sidebarLayout(
             
             ### Sidebar panel for inputs ----
             sidebarPanel(
               width = 2,
               
               #### general plot options ----
               tags$div(style = 'text-align: center',
                        tags$h3('General plot options')),
               
               ##### Input: Slider for frequency, display as confidence ----
               sliderInput(inputId = "freq",
                           label = "Frequency",
                           min = 3,
                           max = 5,
                           value = 5),
               
               ##### Input: checkbox group for plot options, display as plot options ----
               checkboxGroupInput(inputId = 'plotOpts',
                                  label = 'Plot options',
                                  choices = list('Only in GO:DDR' = 'GO:DDR',
                                                 'Only in Pendragon' = 'Pendragon',
                                                 'Only with direct edges' = 'Direct'),
                                  selected = "Direct"),
               
               #### options to plot interaction between queried genes ----
               br(),
               br(),
               br(),
               tags$div(style = 'text-align: center',
                        tags$h3('Plot queried genes')),
               
               ##### Input: Slider for importance score, display as importance ----
               sliderInput(inputId = "impt",
                           label = "Importance",
                           min = 0,
                           max = 1,
                           step = 0.001,
                           value = 0),
               
               ##### Input: Slider for recurrence, display as Expansion level ----
               sliderInput(inputId = "recr",
                           label = "Expansion level",
                           min = 1,
                           max = 3,
                           step = 1,
                           value = 1),
               
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
                              choices = unique(c(cv.perm.impt.network$gene_1,
                                                 cv.perm.impt.network$gene_2))),
               
               ##### Input: action button to query gene interactions, display as plot ----
               tags$div(style = 'text-align: center',
                        actionButton(inputId = 'plotQueries',
                                     label = 'Plot')),
               
               #### options to plot interaction between queried genes ----
               br(),
               br(),
               br(),
               br(),
               tags$div(style = 'text-align: center',
                        tags$h3('Plot clusters')),
               
               ##### Input: select input for cluster plot, display as cluster number ----
               selectInput(inputId = 'cluster',
                           label = 'Cluster ID',
                           choices = NULL),
               
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
               plotOutput(outputId = "plot",
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
  ## select queried genes ----
  # updateSelectInput(session,
  #                   inputId = 'genes',
  #                   choices = unique(c(cv.perm.impt.network$gene_1,
  #                                      cv.perm.impt.network$gene_2)))

  ## render selectInput to select cluster ----
  observeEvent(input$freq, {
    choices <- switch(input$freq,
                      "1" = "",
                      "2" = "",
                      "3" = clustStatFreq3 %>% arrange(desc(DdrProp)) %>% pull(cluster),
                      "4" = clustStatFreq4 %>% arrange(desc(DdrProp)) %>% pull(cluster),
                      "5" = clustStatFreq5 %>% arrange(desc(DdrProp)) %>% pull(cluster))

    updateSelectInput(session, "cluster",
                      choices = choices,
                      selected = NULL)
  })
  
  # Reactively generate plot based on selected button
  observeEvent(input$plotQueries,
               {
                 isolatedInput = isolate(input)

                 net = getNet(igraphFreq, netAll, isolatedInput)

                 # plot network
                 output$plot = renderPlot({
                   plotNet(net, isolatedInput)
                 })
                })

  observeEvent(input$plotCluster,
               {
                 isolatedInput = isolate(input)

                 clust = getClust(igraphFreq, isolatedInput)

                 # plot network
                 output$plot = renderPlot({
                   plotClust(clust, isolatedInput)
                 })
                })
  
  # print out list of genes included in DDRnet
  output$Pendragon = renderText({
    paste(pendragon %>% 
            filter(Original_symbol %in% geneXgboost) %>% 
            pull(Original_symbol), 
          sep = ', ')
  })
  
  output$otherGenes = renderText({
    paste(geneXgboost[!(geneXgboost %in% pendragon$Original_symbol)], sep = ', ')
  })
}


shinyApp(ui = ui, server = server)
