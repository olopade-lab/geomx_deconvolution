fluidPage(useShinyjs(),
          theme = shinytheme("lumen"),
          titlePanel("GeoMx Cell-type Deconvolution"),
          sidebarLayout(
            sidebarPanel(
              selectInput(inputId = "slide",
                          label = strong("Slide Name"),
                          choices = unique(ROI$SlideName),
                          selected = ROI$SlideName[1]),
              checkboxInput(inputId = "immune",
                            label = strong("Immune only"),
                            value = TRUE),
              conditionalPanel(condition = "input.immune == false",
                               selectInput(inputId = "tumor_only_ROI",
                                           label = strong("Tumor-only ROI"),
                                           choices = c(),
                                           multiple = TRUE)),
              conditionalPanel(condition = "input.tabs == 'Clustering'",
                               numericInput(inputId = "clusters",
                                            label = strong("Number of Clusters"),
                                            value = 2,
                                            min = 2,
                                            step = 1)),
              actionButton("run", "Run analysis")),
            mainPanel(
              tabsetPanel(id = "tabs",
                          type = "tabs",
                          tabPanel("Deconvolution", fluidRow(column(6, plotOutput("fluorescent_deconv", width = "100%")),
                                   column(6, plotOutput("outline_deconv", width = "100%")))),
                          tabPanel("Clustering", fluidRow(column(6, plotOutput("fluorescent_cluster", width = "100%")),
                                                          column(6, plotOutput("outline_cluster", width = "100%"))),
                                   fluidRow(column(12, plotOutput("heatmap", width = "100%")))),
                          tabPanel("Results Table", column(12, DT::dataTableOutput("table")))
                          )
              )
            )
          )

