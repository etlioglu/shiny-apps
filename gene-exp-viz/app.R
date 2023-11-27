library(shiny)
library(geneSynonym)
library(ggplot2)

exp_meta_df <- readRDS("downsampled-deseq2-normalized-expression-for-tcgaovc-and-all-gtex.rds")
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
      width = 2
    ),
    
    mainPanel(textOutput("gene"),
              plotOutput("gene_exp_plot"),
              width = 11)
  )
)


server <- function(input, output) {
  
  output$gene <- renderText({
    check_gene(input$gene_symbol)
  })
  
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
}

# Run the application
shinyApp(ui = ui, server = server)
