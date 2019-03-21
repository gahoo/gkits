#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyFiles)
library(pamr)
library(DT)
library(magrittr)
library(org.Hs.eg.db)
library(pheatmap)

get_pamr_result_df <- function(x){
  df <- data.frame(
    threshold = format(round(x$threshold, 3)),
    nonzero = format(trunc(x$size)), 
    errors = trunc(x$error * nrow(x$yhat))
  )
  
  if (!is.na(x$pvalue.survival[1])) {
    df$pvalue <- pvalue = round(x$pvalue.survival, 6)
  }
  df
}

pamr.plotfdr2 <- function (fdrfit, call.win.metafile = FALSE) 
{
  om = fdrfit$results[, "Number of significant genes"] > 0
  na.min = function(x) {
    min(x[!is.na(x)])
  }
  na.max = function(x) {
    max(x[!is.na(x)])
  }
  plot(fdrfit$results[om, "Number of significant genes"], fdrfit$results[om, 
                                                                         "Median FDR"], log = "x", xlab = "Number of genes called significant", 
       ylab = "False discovery rate (median and 90th percentile)", 
       type = "b", ylim = c(na.min(fdrfit$results[om, "Median FDR"]), 
                            na.max(fdrfit$results[om, "90th percentile of FDR"])))
  x = fdrfit$results[om, "Number of significant genes"]
  xlim <- range(x)
  barw <- abs((log(x))) * 1.2
  upper = fdrfit$results[om, "90th percentile of FDR"]
  lower = fdrfit$results[om, "Median FDR"]
  segments(x, upper, x, lower, lty = 2)
  segments(x - barw, upper, x + barw, upper, lty = 2)
  axis(3, at = fdrfit$results[om, "Number of significant genes"], 
       labels = round(fdrfit$results[om, "Threshold"], 2))
  mtext("Threshold", 3, 2, cex = 1)
  if (call.win.metafile) {
    savePlot("", type = "wmf")
    dev.off()
  }
  return()
}

add_external_links <- function(base, query){
  na_idx <- is.na(query)
  links <- sprintf('<a href="%s%s" target="_blank">%s</a>',
          base, query, query)
  links[na_idx] <- NA
  links
}

add_links <- function(df){
  has_ensembl <- sum(grepl('ENS', df$id)) > 0
  if(has_ensembl){
    df %>%
      dplyr::mutate(
        id =  add_external_links('http://asia.ensembl.org/Homo_sapiens/Gene/Summary?g=', id),
        symbol = add_external_links('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', symbol) )
  }else{
    df
  }
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Pam: prediction analysis for microarrays"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        shinyFilesButton('sampleAnnotFile', label='Select sample annotation file', title='Please select a file', multiple=FALSE),
        textOutput('sampleAnnotFilePath'),
        uiOutput('sampleColumn'),
        uiOutput('classColumn'),
        shinyFilesButton('matrixFile', label='Select matrix file for pamr', title='Please select a file', multiple=FALSE),
        textOutput('matrixFilePath'),
        shinyFilesButton('matrixAnnoFile', label='Select feature annotation file', title='Please select a file', multiple=FALSE),
        textOutput('matrixAnnoFilePath'),
        numericInput('threshold', 'threshold', min=0, value=0, step=0.1),
        shinySaveButton('save', 'Save Plots and Genes', 'Save as ...')
      ),
      # Show a plot of the generated distribution
      mainPanel(
         tabsetPanel(
           tabPanel(
             'samples',
             dataTableOutput('annot_preview')
             ),
           tabPanel(
             'overview',
             htmlOutput('pamr')
           ),
           tabPanel(
             'results',
             h2('cross-validate the nearest shrunken centroid classifier'),
             dataTableOutput('results')
           ),
           tabPanel(
             'cv',
             h2('cross-validated error curves'),
             plotOutput('plotcv', height = "1000px")
           ),
           tabPanel(
             'confusion',
             h2('true versus predicted values'),
             htmlOutput('confusion')
           ),
           tabPanel(
             'cvprob',
             h2('cross-validated sample probabilities'),
             plotOutput('plotcvprob', height = "600px")
           ),
           tabPanel(
             'centroids',
             h2('shrunken class centroids'),
             plotOutput('plotcen', height = "600px")
           ),
           tabPanel(
             'geneplot',
             h2('the genes that surive the thresholding from the nearest shrunken centroid classifier'),
             plotOutput('geneplot', height = "1000px")
           ),
           tabPanel(
             'fdr',
             h2('estimate false discovery rates'),
             plotOutput('fdr', height = "800px")
           ),
           tabPanel(
             'genes',
             h2('the genes that survive the thresholding'),
             dataTableOutput('genelist')
           ),
           tabPanel(
             'heatmap',
             h2('heatmap for the genes that survive the thresholding'),
             plotOutput('heatmap')
           )
         )
      )
   )
)

# Define server logic required to draw a histogram
server <- function(session, input, output) {
  sampleAnnotFilePath <- reactive({
    res <- parseFilePaths(c(root='/'), input$sampleAnnotFile)$datapath
    ifelse(length(res) > 0, res, '')
  })
  
  shinyFileChoose(input,
                  'sampleAnnotFile',
                  root=c(root='/'),
                  defaultPath = '/data2/lijiaping/Methylome/lung_adenocarcinoma/meth/',
                  filetypes=c('', 'txt', 'csv', 'tsv'))
  output$sampleAnnotFilePath <- renderText({sampleAnnotFilePath()})
  
  sampleAnnotFileContent <- reactive({
    if(sampleAnnotFilePath() != ''){
      read.csv(sampleAnnotFilePath())
    }else{
      NULL
    }
  })
  
  output$annot_preview <- renderDT({
    datatable(sampleAnnotFileContent())
  })
  
  output$sampleColumn <- renderUI({
    choices <- colnames(sampleAnnotFileContent())
    selectInput('sampleColumn', 'Choose sample column', choices = choices)
  })
  
  output$classColumn <- renderUI({
    choices <- colnames(sampleAnnotFileContent())
    selectInput('classColumn', 'Choose class column', choices = choices)
  })
  
  matrixAnnoFilePath <- reactive({
    res <- parseFilePaths(c(root='/'), input$matrixAnnoFile)$datapath
    ifelse(length(res) > 0, res, '')
  })
  
  shinyFileChoose(input, 
                  'matrixAnnoFile',
                  root=c(root='/'),
                  defaultPath = '/data2/lijiaping/Methylome/lung_adenocarcinoma/meth/',
                  filetypes=c('', 'txt', 'csv', 'tsv'))
  output$matrixAnnoFilePath <- renderText({matrixAnnoFilePath()})
  
  matrixAnno <- reactive({
    if(matrixAnnoFilePath() != ""){
      read.delim(matrixAnnoFilePath())
    }else{
      data.frame()
    }
  })
  
  matrixFilePath <- reactive({
    res <- parseFilePaths(c(root='/'), input$matrixFile)$datapath
    ifelse(length(res) > 0, res, '')
  })
  
  shinyFileChoose(input, 
                  'matrixFile',
                  root=c(root='/'),
                  defaultPath = '/data2/lijiaping/Methylome/lung_adenocarcinoma/meth/',
                  filetypes=c('', 'txt', 'csv', 'tsv'))
  output$matrixFilePath <- renderText({matrixFilePath()})
  
  matrix_ncol <- reactive({
    ncol(read.delim(matrixFilePath(), nrows = 1, header=F))
  })
  
  matrixPamr <- reactive({
    if(matrixFilePath() == '') {
      return(NULL)
    }
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Loading Data", value = 0.1)
    
    pamr.from.excel(matrixFilePath(), matrix_ncol()) %>%
      pamr.knnimpute() -> m_data
    idx <- sampleAnnotFileContent()[[input$sampleColumn]] %in% m_data$y
    is_the_same <- identical(as.character(sampleAnnotFileContent()[[input$sampleColumn]][idx]), m_data$y)
    if(!is_the_same){
      stop('Samples not match in annotation file.')
    }
    m_data$sample_labels <- m_data$y
    m_data$y <- sampleAnnotFileContent()[[input$classColumn]][idx]
    m_data
  })
  
  output$pamr <- renderUI({
    txt <- capture.output(str(matrixPamr()))
    paste(txt, collapse = '</br>') %>%
      HTML()
  })
  
  train <- reactive({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Traning Data", value = 0.2)
    
    pamr.train(matrixPamr())
  })
  results<- reactive({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "cross validating", value = 0.3)
    
    pamr.cv(train(), matrixPamr())
  })
  
  output$results <- renderDT({
    results() %>%
      get_pamr_result_df() ->
      res
    
    idx <- head(which(res$errors == min(res$errors)), n=1)
    val <- res$threshold[idx]
    updateNumericInput(session, "threshold", value = val)
    
    datatable(res)
  })
  
  output$plotcv <- renderPlot({
    pamr.plotcv(results())
  })
  
  output$confusion <- renderUI({
    txt <- capture.output(pamr.confusion(results(), threshold=input$threshold))
    paste(txt, collapse = '</br>') %>%
      HTML()
    tt <- pamr.confusion(results(), threshold=input$threshold, extra = F)
    tt1 <- tt
    diag(tt1) <- 0
    overall_err_rate <- round(sum(tt1)/sum(tt))
    list(
      datatable(tt, rownames = F),
      div('Overall error rate:', overall_err_rate)
    )
  })
  
  output$plotcvprob <- renderPlot({
    pamr.plotcvprob(results(), matrixPamr(), threshold=input$threshold)
  })
  
  output$plotcen <- renderPlot({
    pamr.plotcen(train(), matrixPamr(), threshold=input$threshold)
  })
  
  output$geneplot <- renderPlot({
    pamr.geneplot(train(), matrixPamr(), threshold=input$threshold)
  })
  
  output$fdr <- renderPlot({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    progress$set(message = "caculating fdr", value = 0.9)
    fdr.obj <- pamr.fdr(train(), matrixPamr())
    progress$set(message = "caculating fdr", value = 0.95)
    pamr.plotfdr2(fdr.obj)
  })
  
  survived <- reactive({ 
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "get survived genes", value = 0.5)
    
    pamr.listgenes(train(), matrixPamr(), threshold=input$threshold) %>%
      as.data.frame(stringsAsFactors=F) -> survived

    has_ensembl <- sum(grepl('ENS', survived$id)) > 0
    has_anno_file <- nrow(matrixAnno()) > 0
    if(has_anno_file){
      matrixAnno()[survived$id,] %>%
        tibble::rownames_to_column('id') %>%
        dplyr::mutate(Length = End - Start) %>%
        dplyr::right_join(survived, by='id')
    }else if(has_ensembl){
      progress$set(message = "annotate survived genes", value = 0.8)
      anno <- select(org.Hs.eg.db,
                     keys=as.character(survived$id),
                     columns=c("SYMBOL","GENENAME"),
                     keytype="ENSEMBL") %>%
        dplyr::rename(id=ENSEMBL,
                      symbol=SYMBOL) %>%
      dplyr::right_join(survived, by='id')
    }else{
      survived
    }
    })
  
  genelist <- reactive({
    survived() %>%
      add_links() %>%
      datatable(escape=F)
  })
  
  output$genelist <- renderDataTable({
    genelist()
  })
  
  survivedMatrixPamr <- reactive({
    idx <- which(matrixPamr()$geneid %in% survived()$id)
    m <- matrixPamr()$x[idx, ]
    rownames(m) <- survived()$id
    colnames(m) <- matrixPamr()$sample_labels
    m
  })
  
  output$heatmap <- renderPlot({
    anno_df <- sampleAnnotFileContent()[input$classColumn]
    row.names(anno_df) <- sampleAnnotFileContent()[[input$sampleColumn]]
    pheatmap(survivedMatrixPamr(), annotation_col = anno_df, show_colnames = F)
  })
  
  shinyFileSave(input, 'save', 
                root=c(root='/'),
                defaultPath = '/data2/lijiaping/Methylome/lung_adenocarcinoma/meth/')
  
  observeEvent(input$save, {
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    save <- parseSavePath(root=c(root='/'), input$save)
    prefix <- paste(
      save$datapath,
      input$threshold,
      sep='-'
    )
    
    if(length(save$name) != 0){
      progress$set(message = "get survived genes", value = 0.1)
      write.csv(survived(), file=paste(prefix, 'csv', sep='.'))
      DT::saveWidget(genelist(), paste(prefix, 'html', sep='.'))
      
      pdf(paste(prefix, 'pdf', sep='.'), height = 16, width = 16)
      progress$set(message = "plot cross-validated error curves", value = 0.2)
      pamr.plotcv(results())
      progress$set(message = "cross-validated sample probabilities", value = 0.3)
      pamr.plotcvprob(results(), matrixPamr(), threshold=input$threshold)
      progress$set(message = "shrunken class centroids", value = 0.4)
      pamr.plotcen(train(), matrixPamr(), threshold=input$threshold)
      progress$set(message = "plot survived genes", value = 0.5)
      pamr.geneplot(train(), matrixPamr(), threshold=input$threshold)
      progress$set(message = "caculating fdr", value = 0.6)
      fdr.obj<- pamr.fdr(train(), matrixPamr())
      progress$set(message = "caculating fdr", value = 0.7)
      pamr.plotfdr2(fdr.obj)
      progress$set(message = "ploting heatmap", value = 0.8)
      anno_df <- sampleAnnotFileContent()[input$classColumn]
      row.names(anno_df) <- sampleAnnotFileContent()[[input$sampleColumn]]
      pheatmap(survivedMatrixPamr(), annotation_col = anno_df, show_colnames = F)
      dev.off()
    }

  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

