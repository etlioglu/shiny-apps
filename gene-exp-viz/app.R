library(shiny)
library(geneSynonym)
library(ggplot2)

source(
  "https://raw.githubusercontent.com/etlioglu/bioinformatics/main/templates/customize-ggplot.r"
)

exp_meta_df <-
  readRDS("deseq2-normalized-expression-for-ovc-colon-gbm-and-all-gtex.rds")
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

filter_data_generate_plot <- function(data, gene, cancer) {
  if (cancer == "ovarian") {
    data <-
      data[!(
        data$TCGA_GTEX_main_category %in% c("Colorectal cancer", "Glioblastoma Multiforme")
      ), ]
  } else if (cancer == "colorectal") {
    data <-
      data[!(
        data$TCGA_GTEX_main_category %in% c(
          "Ovarian Serous Cystadenocarcinoma",
          "Glioblastoma Multiforme"
        )
      ), ]
  } else if (cancer == "gbm") {
    data <-
      data[!(
        data$TCGA_GTEX_main_category %in% c("Ovarian Serous Cystadenocarcinoma", "Colorectal cancer")
      ), ]
  }
  
  ggplot(
    data,
    aes_string(x = "TCGA_GTEX_main_category",
               y = gene,
               color = "TCGA_GTEX_main_category")
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
      tabsetPanel(
        type = "tabs",
        tabPanel("All cancers", plotOutput("gene_exp_plot_all")),
        tabPanel("Ovarian Cancer", plotOutput("gene_exp_plot_ovarian")),
        tabPanel("Colorectal Cancer", plotOutput("gene_exp_plot_colorectal")),
        tabPanel("GBM", plotOutput("gene_exp_plot_gbm"))
      )
    )
  )
)


server <- function(input, output) {

  output$gene_exp_plot_all <- renderPlot({
    filter_data_generate_plot(
      data = exp_meta_df,
      gene = check_gene(input$gene_symbol),
      cancer = "all"
    )
  })
  
  output$gene_exp_plot_ovarian <- renderPlot({
    filter_data_generate_plot(
      data = exp_meta_df,
      gene = check_gene(input$gene_symbol),
      cancer = "ovarian"
    )
  })
  
  output$gene_exp_plot_colorectal <- renderPlot({
    filter_data_generate_plot(
      data = exp_meta_df,
      gene = check_gene(input$gene_symbol),
      cancer = "colorectal"
    )
  })
  
  output$gene_exp_plot_gbm <- renderPlot({
    filter_data_generate_plot(
      data = exp_meta_df,
      gene = check_gene(input$gene_symbol),
      cancer = "gbm"
    )
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
