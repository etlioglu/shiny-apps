library(data.table)
library(shiny)
library(geneSynonym)
library(ggplot2)

source(
  "https://raw.githubusercontent.com/etlioglu/bioinformatics/main/templates/customize-ggplot.r"
)

exp_meta_df <-
  fread("tcga-gtex-normalized-ovc-colorectal-gbm-bc-protein-coding-only.csv")

genes_present <- colnames(exp_meta_df)[-1]

check_gene <- function(x) {
  if (!(x %in% genes_present)) {
    gene_syns <- unlist(humanSyno(x))
    for (syn in gene_syns) {
      if (syn %in% genes_present) {
        x <- syn
        return(x)
      } else {
        x <- NULL
        return(x)
      }
    }
  } else {
    return(x)
  }
}

ui <- fluidPage(
  titlePanel(
    "Gene expession visualiazion with the TCGA-TARGET-GTEx transcriptomics dataset"
  ),
  
  sidebarLayout(
    sidebarPanel(
      textInput(
        "gene_symbol",
        "Gene symbol",
        value = "ACE",
        width = NULL,
        placeholder = NULL
      ),
      
      # download button
      downloadButton("download_data", "Download"),
      
      width = 2
    ),
    
    mainPanel(
      textOutput("gene"),
      
      tabsetPanel(
        tabPanel("Four cancers", plotOutput("gene_exp_plot"))
      )
    )
  )
)


server <- function(input, output) {
  # save check_gene(input$gene_symbol) in a variable in order not to call this multiple times
  
  # output$gene <- check_gene(reactive(input$gene_symbol))
  
  output$gene_exp_plot <- renderPlot({
    ggplot(
      exp_meta_df,
      aes_string(
        x = "TCGA_GTEX_main_category",
        y = check_gene(input$gene_symbol),
        color = "TCGA_GTEX_main_category"
      )
    ) +
      geom_jitter(size = 0.5) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ),
        legend.position = "none"
      ) +
      coord_flip()
  })
  
  output$to_download <-
    
    output$download_data <- downloadHandler(
      filename = function() {
        paste0(check_gene(input$gene_symbol), ".csv")
      },
      content = function(file) {
        write.table(
          exp_meta_df[, c("sample",
                          "TCGA_GTEX_main_category",
                          check_gene(input$gene_symbol))],
          file = file,
          quote = FALSE,
          sep = ",",
          row.names = FALSE
        )
      }
    )
  
}

# Run
shinyApp(ui = ui, server = server)
