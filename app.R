# Courtesy by Dr. Michal Kouril for Data science Course, just use it for learning
# In RStudio, block comment, Cmd+Shift+C

# 1. fluidPage is a prepackaged layout, others like fixedPage(rarely used), navbarPage,dashboard ...
# 2. sidebarlayout is a prepackaged layout, it contains sidebarPanel() and mainPanel()
# 3. in each panel, there are dozens of Input function and output function, function has InputID as first argument, label as second argument
# inputID will be used for tracking in server function, label is what it will manifest itself on webpage, other
# specific function argument follows this two mandatory arguments
# 4. you can insert html5 tags using tags$p, tags$code, tags$a(href="","Rstudio"), some html5 tags can be directly
# deployed, like hr() as below. 
# 5. ui: will run each time you start R session;
#    server: will run each time you start the browser, click run app;
#    reactive function(render-esque): will run each reaction;
#   So we should put least code on reactive function, and most code on ui portion for the sake of being parsimonious
# 6. server:
#   6.1 render function, a dozens of them. They wrap the reactive value(input$L1000, or create a reactive value
#           as the values as below). When the reactive value change, the object that this render function created
#           wil respond. Isolate could stop them from responding. Refer to the official tuturial example
#   6.2 observeEvent() and oberve could be used to associate button input, if command could be used to
#           associate other input
#   6.3 Cache the data and reuse them, data <- reactive({rnorm()/input$num}), then using data() to replace input$num
#   6.4 session as below
# 7. Publish to shinyapps.io or shiny server
# 8. other layout idea: fluidRow(column(4,offset=8),Plotoutput), example refer to official tutorial.
# 9. Add CSS to the code, example refer to official tutorial
# 10. img and .css file should be stored in a subdirectory named www. www folder and app.R are in the same folder, in this way,
#.  you could just specify the img name, like rt.img instead specifying the whole path.
          



library(shiny)
library(DT)
library(httr)

options(shiny.maxRequestSize=70*1024^2)

ui <- fluidPage(

  # Application title
  titlePanel("File upload"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose TSV File', accept=c('text/tsv','.tsv')),   
      selectInput("variable", "Grouping Variable:", choices=c()),
      selectizeInput('group1', "Group1", choices = NULL, multiple = TRUE),
      selectizeInput('group2', "Group2", choices = NULL, multiple = TRUE),
      selectInput("difffunction", "Differential function:", choices=c("t-test")),
      checkboxInput("L1000",label="limit Genes to L1000", value = FALSE),
      selectInput("limit","Limit signature",choices = c("All","Top 100")),
      actionButton("compute", "Compute and Submit Signature")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      dataTableOutput("signature_data"),
      hr(),
      dataTableOutput("correlated_data"),
      hr(),
      dataTableOutput("sample_data")
    )
  )
)


server <- function(input, output, session) {
  
  values <- reactiveValues(data=NULL)
  
  observe({
    # handle file upload 
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    isolate({ 
      file <- (inFile$datapath)

      values$header <- scan(file, nlines = 1, sep="\t", what = character())
      values$data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
      names(values$data) <- values$header 
      values$header2 <- data.frame(scan(file, skip = 1, nlines = 1, sep="\t", what = character()))
    })
  })
  
  # show sample metadata
  output$sample_data <- renderDataTable({
    head(values$data[,1:10],n=50)
  }, caption = "50 genes, 50 samples")
  
  observe({
     # fill in variable selection from the input file
     updateSelectInput(session, "variable", choices=as.character(values$header2[1,]))
  })
  
  # handle UI group values update based on selected variable
  observe({
    if (input$variable !="") {
      
      updateSelectizeInput(session, 'group1', choices=unique(values$header2[-1,values$header2[1,]==input$variable]), server = TRUE)
      updateSelectizeInput(session, 'group2', choices=unique(values$header2[-1,values$header2[1,]==input$variable]), server = TRUE)
    }
  })
  
  # Create signature, upload to API, display results
  observeEvent(input$compute, {
    withProgress(message = 'Creating signature', value = 0, {
      
      #
      # Filter into two groups
      #
      group1 <- names(values$data)[values$header2==input$group1]
      group2 <- names(values$data)[values$header2==input$group2]

      #
      # ensure the values are numeric
      #
      values$values <- values$data[complete.cases(values$data), -1]
      values$values[] <- lapply(values$values, function(x) { as.numeric(as.character(x)) })
      
      #
      # select differential function
      #
      incProgress(1/3, detail = paste0("Running ",input$difffunction))
      if (input$difffunction=="t-test") {
        diff_result <- as.data.frame(apply(values$values, 1, 
                                           function(x) t.test(unlist(x[group1], use.names = FALSE),
                                                              unlist(x[group2], use.names = FALSE))$p.value))
        
        # FIXME: fill in
        # Add t statistic as the measure of differential expression (Value_LogDiffExp)
        diff_result$Value_LogDiffExp <- apply(values$values, 1, 
                                              function(x) t.test(unlist(x[group1], use.names = FALSE),
                                                                 unlist(x[group2], use.names = FALSE))$statistic)
      }
      
      #
      # format signature output
      #
      output_id_column_name <- paste0(values$header[1],"_GeneSymbol")
      diff_result <- data.frame(values$data[values$header[1]], diff_result[,1:2])
      colnames(diff_result) <- c(output_id_column_name, "Significance_pvalue", "Value_LogDiffExp")
      
      # FIXME: fill in
      # choose only L1000 genes
      # l1000genes <- ...
      # diff_result <- ...
      if(input$L1000){
        genes <- read.csv2("L1000.txt",sep="\t")
        L1000 <- genes[genes$pr_is_lm=='1',]$pr_gene_symbol
        diff_result <- diff_result[diff_result[,1] %in% L1000] 
      }
      
      
      # FIXME: fill in
      # choose top 100 most differentially expressed genes
      # diff_result <- ...
      if(input$limit == "Top 100"){
        diff_result <- head(diff_result[order(-abs(diff_result$Value_LogDiffExp)),],100)
      }
      
      
      #
      # show signature in a table
      #
      output$signature_data <- DT::renderDataTable({
        diff_result
      }, caption = "Signature to submit to iLincs")
      
      incProgress(1/3, detail = paste("Submitting the signature to iLincs"))
      
      #
      # create temporary csv file to submit into API
      #
      ftemp <- tempfile(pattern = "file", fileext=".csv", tmpdir = tempdir())
      write.csv(diff_result, ftemp, row.names = F, quote=F)
      cat(ftemp)
      
      r <- POST("http://www.ilincs.org/api/SignatureMeta/uploadAndAnalyze?lib=LIB_5", body = list(file = upload_file(ftemp)))
      
      l <- lapply(content(r)$status$concordanceTable, function(x) unlist(x))
      ilincs_result <- data.frame(t(sapply(l,c)))
      
      #
      # show correlation results
      #
      output$correlated_data <- DT::renderDataTable({
        datatable( ilincs_result, rownames = TRUE, caption = "Correlated signatures")
      })
    })
  })
}

shinyApp(ui = ui, server = server)

