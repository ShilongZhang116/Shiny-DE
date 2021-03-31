options(shiny.maxRequestSize = 30*1024^2)
# packages used in this app
library(shiny)
library(shinydashboard)
library(DiffBind)
library(bslib)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(org.Hs.eg.db)
# packages used in this app
library(DT)
source("volcano_plot.R")

annoData <<- toGRanges(EnsDb.Hsapiens.v75, feature="gene")


# Build the head
header <- dashboardHeader(
  title = "ATAC-Diff"
)
###########################################################
# Build the sidebar
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Instruction", tabName="instruction", 
             icon=icon("glyphicon glyphicon-book", lib="glyphicon")),
    menuItem("Upload Your Data",tabName="dataUpload",selected=TRUE,
             icon=icon("glyphicon glyphicon-open", lib="glyphicon")),
    menuItemOutput("dataReport"),
    menuItemOutput("filter"),
    menuItemOutput("DE_Analysis"),
    menuItemOutput("TF_Enrichment")
    
  )
)
###########################################################
# Build the body
body <- dashboardBody(
  
  tabItems(
    tabItem(tabName = "Instruction",
            h5("Get a quick start of ATAC-Diff")
      
    ),
    
    tabItem(tabName = "dataUpload",
            # Buile the Page of Data Upload
            fluidRow(
              box(
                title = "Upload your data", status = "primary", solidHeader = TRUE,
                collapsible = TRUE,width=12,
                selectInput("data_sel",label="Select the data",choices=c("Example data" ,"Your data")),
                fluidRow(
                  column(6,fileInput("userSample", "Please upload your dba object after contrast",multiple = FALSE)),
                  column(6,fileInput("peakSample", "Please upload your Peak Calling",multiple = TRUE,accept=".bed")),
                  column(6,selectInput("factor","Single or Multiple Factor analysis",
                                       choices=c("Single Factor","Multiple Factor"),selected = 1))
                ),
                actionButton("start", "Do Analysis!"),
                
              ),
              
              uiOutput("multipleFactorPanel"),
              
              box(
                title="Check your Sample Sheet",status="primary",solidHeader=TRUE,
                collapsible = TRUE, width=12,
                tableOutput("sampleSheet"),
              )
              
            )),
    
######################################    
    # Data Report
    tabItem(tabName="dataReport",
            tabPanel("Report",DT::dataTableOutput("dbaReport")),
            ),
    
#######################################    
    tabItem(tabName="DE_Analysis",
            tabsetPanel(
              tabPanel("Venn Plot",
                       fluidRow(
                         box(title="Description of Venn Plot",status="primary",solidHeader = TRUE,
                             collapsible = TRUE,width=12,
                             p("Venn diagrams illustrate overlaps between different sets of peaks. 
                             For example, amongst the differentially bound sites, we can see 
                             the differences between the 'Gain' sites (those that 
                             increase binding enrichment in the Resistant condition) 
                             and the 'Loss' sites (those with lower enrichment) as follows:")),
                         
                         box(title="Venn Plot",status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=8,
                             plotOutput("dba_venn"))
                       )),
              
              tabPanel("PCA",
                       fluidRow(
                         box(title="Description of PCA Plot",status = "primary", solidHeader = TRUE,
                             collapsible = TRUE, width=12,
                             p("While the correlation heatmaps already seen are 
                             good for showing clustering, plots based on principal 
                             components analysis can be used to give a deeper insight 
                             into how samples are associated.")),
                         
                         box(title="PCA Plot",status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=8,
                             plotOutput("dba_PCA"))
                       )),
              
              tabPanel("Volcano",
                       fluidRow(
                         box(title="Description of Volcano Plot",status = "primary",solidHeader = TRUE,
                             collapsible = TRUE,width=12,
                             p("Similar to MA plots, Volcano plots also highlight 
                             significantly differentially bound sites and show their 
                             fold changes. Here, however, the confidence statistic
                             (FDR or p-value) is shown on a negative log scale, 
                             helping visualize the relationship between the magnitude 
                             of fold changes and the confidence that sites are 
                               differentially bound.")),
                         
                         box(title="Volcano Plot",status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=8,
                             plotOutput("dba_Volcano",click = "plot_click")),
                         box(title="Information",statue="primary",solidHeader=TRUE,
                             collapsible=TRUE,width=4,
                             verbatimTextOutput("info"))
                       )),
              
              tabPanel("Heatmap",
                       fluidRow(
                         box(title="Description of Heatmap",status = "primary",solidHeader = TRUE,
                             collapsible = TRUE,width=12,
                             p("Another way to view the patterns of binding affinity 
                             directly in the differentially bound sites is via a 
                             binding affinity heatmap, showing the read scores for 
                             some or all of the binding sites.")),
                         
                         box(title="Heatmap",status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=8,
                             plotOutput("dba_heatmap"))
                       ))
              
              
              
              )),
            
#######################################    
    tabItem(tabName="TF_Enrichment",
            tabsetPanel(
              tabPanel(title="Venn Plot",
                       fluidRow(
                         box(title="Venn Plot", status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=8,
                             plotOutput("ChIPpeakAnno_Venn")),
                         tabBox(title ="Option", width=4,id = "overlapOption", height = "250px",
                                tabPanel("Venn",
                                         uiOutput("vennPlotOption"),
                                         htmlOutput("vennPlotOption_Message"),
                                         uiOutput("vennPlotOption_SelectColor"),
                                         actionButton("vennPlot",tags$b("Plot"),class="btn-primary",style="color:white")),
                                tabPanel("Overlaps",
                                         uiOutput("overlapOption"),
                                         sliderInput("DPATSS_radius", # DPATSS is Distribution of Peaks Around Transcript Start Sites
                                                     p("The radius of the longest distance to feature site"),
                                                     min=2000,max=8000,value=5000,step=500),
                                         sliderInput("DPATSS_nbins",
                                                     p("The number of bins"),
                                                     min=10,max=90,value=50),
                                         actionButton("Distribution_of_Peaks_Around_Transcript_Start_Sites",
                                                      "Plot",class="btn-primary",style="color:white")
                                         
                                         )),
                         
                         box(title="Binding Site Distribution Relative to Features",
                             status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=8,
                             p("You can plot the distribution of the distance of 
                             overlapped peaks to the nearest feature such as the 
                               transcription start sites (TSS)."),
                             plotOutput("bindingSiteDistribution"))
                         
                        
                        )
                       ),
              
              tabPanel("Genomic Element Distribution",
                       fluidRow(
                         box(title="Genomic Element Distribution of Duplicates", status = "primary", solidHeader = TRUE,
                             collapsible = TRUE, width=10,
                             plotOutput("ChIPpeakAnno_DuplicateDistribution")),
                         box(title="Option",status="primary",solidHeader = TRUE,
                             collapsible = TRUE,width=2,
                             uiOutput("duplicateDistribution_Option"),
                             actionButton("duplicateDistributionPlot","Plot",class="btn-primary",style="color:white")),
                         box(title="Genomic Element Distribution of Overlaps", status = "primary", solidHeader = TRUE,
                             collapsible = TRUE, width=6,
                             plotOutput("ChIPpeakAnno_OverlapDistribution"))
                       )),
              
              tabPanel("Annotate Peaks",
                       fluidRow(
                         box(title="Annotate the Peaks to the Promoter Regions", status = "primary", solidHeader = TRUE,
                             collapsible = TRUE, width=8,
                             verbatimTextOutput("ChIPpeakAnno_PeakAnnotation_Table")),
                         box(title="Pie Plot",status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,width=4,
                             plotOutput("ChIPpeakAnno_PeakAnnotation_Pie"))
                       ))
              
            )
            )
  )
)

ui <- dashboardPage(
  header,
  sidebar,
  body
)


###########################################################
# Service Section

server <- function(input, output) {
###########################################################
  #Upload user's data
  output$userSample <- renderUI({
    if(input$data_sel == "Your data"){
      fileInput("userSample", "Please upload your dba object after contrast",multiple = FALSE)
    }
  })
  
  # Side bar menu output: Data Report
  output$dataReport <- renderMenu({
    if(input$start!=0){
      menuItem("Data Report",tabName="dataReport",icon=icon("glyphicon glyphicon-list-alt",lib="glyphicon"))
    }
  })
  
  # Side bar menu output: DE analysis
  output$DE_Analysis <- renderMenu({
    if(input$start!=0){
      menuItem("DE Analysis",tabName="DE_Analysis",icon=icon("glyphicon glyphicon-equalizer",lib="glyphicon"))
    }
  })
  
  # Side bar menu output: TF Enrichment
  output$TF_Enrichment <- renderMenu({
    if(input$start!=0){
      menuItem("TF Enrichment",tabName="TF_Enrichment",icon=icon("glyphicon glyphicon-stats",lib="glyphicon"))
    }
  })
  
  # Side bar menu output: Control
  output$filter <- renderMenu({
    if(input$start!=0){
      menuItem("Filter",tabName="Filter",icon=icon("glyphicon glyphicon-filter",lib="glyphicon"),
               sliderInput("fdr",label="FDR",min=1e-10, max=0.05, value=0.005),
               sliderInput("fold_change",label="log2FoldChange",min=0.5,max=5,value=0.5),
               selectInput("plot_method",label="Analysis Method",choices=c("DESeq2"=DBA_DESEQ2,"EdgeR"=DBA_EDGER)))
    }
  })
  # Multiple Factor Choice Panel
  output$multipleFactorPanel <- renderUI({
    if(input$factor == "Multiple Factor"){
    box(title="Select the Multiple Factors",
        status = "primary", solidHeader = TRUE,
        collapsible = TRUE,width=12,
        checkboxGroupInput("multipleFactorSelectBox","Choose the factors",
                           choices = c("Tissue","Factor","Condition","Treatment","Replicate","Caller"),
                           selected = 2))
        
        
    }
  })
  
  # Input the dba object after conting
  sampleInput <- reactive({
    if(input$data_sel == "Example data"){
      sample <- load("Data/dba_object")
      tamoxifen_sample <- get(sample)
      return(tamoxifen_sample)
    }else if(input$data_sel == "Your data"){
      req(input$userSample)
      userSample <- load(input$userSample$datapath)
      userSample <- get(userSample)
      return(userSample)
    }
  })
  
  # Dynamic ui:contrast model
  output$contrastModel <- renderUI({
    sampleSheet <- sampleInput()
    sampleSheet <- sampleSheet$samples
    condition <- unique(sampleSheet$Condition)
    selectInput("contrastModel",label="Contrast Model",
                choices=condition)
  })
  
  # Do Diffbind analysis
  dbaAnalysis <- reactive({
    samples <- sampleInput()
    dbaObject <- dba(samples)
    
    # Normalizing the data
    dbaObject <- dba.normalize(dbaObject, method = DBA_DESEQ2,
                               normalize = DBA_NORM_DEFAULT, library = DBA_LIBSIZE_DEFAULT, 
                               background = FALSE, bRetrieve=FALSE)
    
    if(input$factor == "Single Factor"){
    # Establishing a model design and contrast
    # Performing the differential analysis
      dbaObject <- dba.contrast(dbaObject,reorderMeta=list(Condition="Resistant"))
      dbaObject <- dba.analyze(dbaObject,bBlacklist=FALSE,bGreylist=FALSE)
    }else if(input$factor == "Multiple Factor"){
      multipleFactorDesign <- vector()
      for(i in input$multipleFactorSelectBox){multipleFactorDesign <- paste(multipleFactorDesign,i,sep=" + ")}
      multipleFactorDesign <- paste("~",multipleFactorDesign,sep="")
      dbaObject <- dba.contrast(dbaObject,design=multipleFactorDesign)
      dbaObject <- dba.analyze(dbaObject,method=DBA_ALL_METHODS,bBlacklist=FALSE,bGreylist=FALSE)
    }
    
    
    
    return(dbaObject)
    
  })
  
  output$sampleSheet <- renderTable({
    sampleSheet_list <- sampleInput()
    sampleSheet <- sampleSheet_list$samples
    sampleSheet <- cbind(sampleSheet[,1:6],sampleSheet[c("ControlID")],sampleSheet[,10:11])
    sampleSheet
  })
  
  output$dbaReport_heatmap <- renderPlot({
    samples <- sampleInput()
    dbaObject <- dba(samples)
    plot(dbaObject)
  })
  
  # Retrieving the differentially bound sites
  dbareport <- reactive({
    withProgress(message = 'Generating the report', value = 0.5,{
    dbaObject <- dbaAnalysis()
    #str(dbaObject)
    tamoxifen.DB <- dba.report(dbaObject, th=input$fdr, bUsePval=FALSE, 
                               fold=input$fold_change,
                               precision=3:5,
                               file="tamoxifen_report.csv",DataType=DBA_DATA_FRAME)
    
    })
    tamoxifen.DB
  })
  
  # Report 
  output$dbaReport = DT::renderDataTable({
    dbareport()
  })
  
  # Venn diagrams
  output$dba_venn <- renderPlot({
    dbaObject <- dbaAnalysis()
    if(input$factor == "Single Factor"){
      dba.plotVenn(dbaObject, contrast=1, bDB=TRUE,bGain=TRUE, bLoss=TRUE, bAll=FALSE)
    }else if(input$factor == "Multiple Factor"){
      dba.plotVenn(dbaObject,contrast=1,method=DBA_ALL_METHODS,bDB=TRUE)
    }
    
    
  })
  
  # PCA plots
  output$dba_PCA <- renderPlot({
    dbaObject <- dbaAnalysis()
    dba.plotPCA(dbaObject,DBA_TISSUE,label=DBA_CONDITION)
    # dba.plotPCA(tamoxifen, contrast=1, label=DBA_TISSUE)
  })
  
  
  # Volcano plots
  output$dba_Volcano <- renderPlot({
    dbaObject <- dbaAnalysis()
    tamoxifen.DB <-  dbareport()
    dba.plotVolcano(dbaObject)
    #DBAplotVolcano(dbaObject)
  })
  
  output$info <- renderPrint({
    tamoxifen.DB <-  dbareport()
    nearPoints(volcano_plot_matrix, input$plot_click, threshold = 10, maxpoints = 1,addDist = TRUE)
  })
  
  # Boxplots
  output$dba_BoxPlot <- renderPlot({
    dbaObject <- dbaAnalysis()
    pvals <- dba.plotBox(dbaObject)
    pvals
  })
  
  #Less Fold2change
  output$dba_LessFoldChange <- renderTable({
    dbaObject <- dbaAnalysis()
    tamoxifen.DB <- dba.report(dbaObject)
    sum(tamoxifen.DB$Fold<input$fold_change)
  })
  
  #More Fold2change
  output$dba_MoreFoldChange <- renderTable({
    dbaObject <- dbaAnalysis()
    tamoxifen.DB <- dba.report(dbaObject)
    sum(tamoxifen.DB$Fold>input$fdr)
  })
  
  output$foldChange <- renderText({
    print(input$fold_change)
  })
  
  # Heatmaps
  output$dba_heatmap <- renderPlot({
    dbaObject <- dbaAnalysis()
    # corvals <- dba.plotHeatmap(dbaObject)
    hmap <- colorRampPalette(c("red", "black", "blue"))(n = 13)
    readscores <- dba.plotHeatmap(dbaObject, contrast=1, correlations=FALSE,
                                  scale="row", colScheme = hmap)
  })

 
###########################################################
  # ChIPpeakAnno section
  # peakData is used to read data and return a summaty list including peakfile 
  # annotion list
  peakData <- reactive({
    
    if(input$data_sel == "Example data"){
      
      peakdata <- read.csv("Data/peak_path.csv")
      peakDataFile <- list()
      ChIPpeakAnno_Grange <- list()
      ChIPpeakAnno_Replicate <-list()
      peakName <- vector()
      
      for(i in 1:length(peakdata$datapath)){
        name <- paste("peak",i,sep="")
        peakDataFile <- c(peakDataFile,name)
        ChIPpeakAnno_Grange <- c(ChIPpeakAnno_Grange,name)
        peakName <- append(peakName,peakDataFile[[i]])
        peakDataFile[[i]] <- file(peakdata$datapath[i])
        ChIPpeakAnno_Grange[[i]] <- toGRanges(peakDataFile[[i]], format="BED", header=FALSE) 
        ChIPpeakAnno_Replicate <- c(ChIPpeakAnno_Replicate,ChIPpeakAnno_Grange[[i]])
      }
      names(ChIPpeakAnno_Grange) <- peakdata$name
      names(ChIPpeakAnno_Replicate) <- peakdata$name
      exampleSummary <<- list(peakDataFile = peakDataFile,
                              ChIPpeakAnno_Grange = ChIPpeakAnno_Grange,
                              ChIPpeakAnno_Replicate = ChIPpeakAnno_Replicate,
                              peakName = peakdata$name)
      return(exampleSummary)
      
    }else if(input$data_sel == "Your data"){
      # User upload the data
      peakDataFile <- list()
      ChIPpeakAnno_Grange <- list()
      ChIPpeakAnno_Replicate <-list()
      peakName <- vector()
      
      for(i in 1:length(input$peakSample$datapath))
        {
        name <- paste("peak",i,sep="")
        peakDataFile <- c(peakDataFile,name)
        ChIPpeakAnno_Grange <- c(ChIPpeakAnno_Grange,name)
        peakName <- append(peakName,peakDataFile[[i]])
        peakDataFile[[i]] <- file(input$peakSample$datapath[i])
        ChIPpeakAnno_Grange[[i]] <- toGRanges(peakDataFile[[i]], format="BED", header=FALSE) 
        ChIPpeakAnno_Replicate <- c(rep,ChIPpeakAnno_Grange[[i]])
      }
      
      names(ChIPpeakAnno_Grange) <- peakdata$name
      names(ChIPpeakAnno_Replicate) <- peakdata$name
      userSummary <<- list(peakDataFile = peakDataFile,
                      ChIPpeakAnno_Grange = ChIPpeakAnno_Grange,
                      ChIPpeakAnno_Replicate = ChIPpeakAnno_Replicate,
                      peakName = peakdata$name)
      return(userSummary)
      
    }
  })
  
  # Return the overlap
  generateOverlapObject <- reactive({
    peakSummary <- peakData()
    overlaps <- findOverlapsOfPeaks(peakSummary$ChIPpeakAnno_Grange[input$vennPlotOption])
    overlaps <-addMetadata(overlaps, colNames="score", FUN=mean)
    return(overlaps)
  })
  
  # Venn Option
  #########################################################
  output$vennPlotOption <- renderUI({
    peakSummary <- peakData()
    peakName <- peakSummary$peakName
    
    checkboxGroupInput("vennPlotOption", 
                       h4("Choose no more than 5"),
                       choices = peakName,
                       selected = 2)
  })
  
  vennPlotOption_Return <- reactive({
    vennPlotOption_UserInput <- input$vennPlotOption
    if(length(vennPlotOption_UserInput)<2){
      message <- "Please choose more samples."
      userInput <- NULL
      resultReturn <- list(message=message,userInput=userInput)
      return(resultReturn)
    }else if(length(vennPlotOption_UserInput)>=2 & length(vennPlotOption_UserInput)<=5){
      message <- "The number of samples is alright."
      userInput <- vennPlotOption_UserInput
      resultReturn <- list(message=message,userInput=userInput)
      return(resultReturn)
    }else if(length(vennPlotOption_UserInput)>5){
      message <- "There are too many samples"
      userInput <- NULL
      resultReturn <- list(message=message,userInput=userInput)
      return(resultReturn)
    }
  })
  
  output$vennPlotOption_Message <- renderText({
    vennPlotOption_UserInput <- vennPlotOption_Return()
    vennPlotOption_UserInput_Message <- vennPlotOption_UserInput[["message"]]
    
    if(is.null(vennPlotOption_UserInput[["userInput"]])){
      paste("<p style='color:red'><b>",vennPlotOption_UserInput_Message,"</b></p>")
      
    }else{
      paste("<p style='color:blue'><b>",vennPlotOption_UserInput_Message,"</b></p>")
    }
    
  })
  
  output$vennPlotOption_SelectColor <- renderUI({
    selectInput("selectColor", h4("Select Color"), 
                choices = list("Color 1" = "color1", 
                               "Color 2" = "color2",
                               "Color 3" = "color3"), 
                selected = 1)
  })
  
  vennPlotOption_Color <- reactive({
    vennPlotOption_UserInput <- input$vennPlotOption
    vennPlotOption_UserInput_Length <- length(vennPlotOption_UserInput)
    if(vennPlotOption_UserInput_Length == 2){
      switch(input$selectColor,
             color1 = c("#8dd3c7","#bebada"),
             color2 = c("#66c2a5","#fc8d62"),
             color3 = c("#a6cee3","#b2df8a"))
    }
    else if(vennPlotOption_UserInput_Length == 3){
      switch(input$selectColor,
             color1 = c("#8dd3c7","#fdb462","#bebada"),
             color2 = c("#66c2a5","#fc8d62","#8da0cb"),
             color3 = c("#a6cee3","#1f78b4","#b2df8a"))
    }
    else if(vennPlotOption_UserInput_Length == 4){
      switch(input$selectColor,
             color1 = c("#8dd3c7","#fdb462","#bebada","#fb8072"),
             color2 = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3"),
             color3 = c("#a6cee3","#1f78b4","#b2df8a","#33a02c"))
    }
    else if(vennPlotOption_UserInput_Length == 5){
      switch(input$selectColor,
             color1 = c("#8dd3c7","#fdb462","#bebada","#fb8072","#80b1d3"),
             color2 = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854"),
             color3 = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99"))
    }
    
  })
  
  # Draw Venn Plot
  output$ChIPpeakAnno_Venn <- renderPlot({
    input$vennPlot
    if (input$vennPlot == 0)return()
    isolate({
      
    peakSummary <- peakData()
    vennPlotOption_UserInput <- vennPlotOption_Return()
    
    withProgress(message = 'Making plot', value = 0.3,{
      # Judge whether user choose right number
      if(is.null(vennPlotOption_UserInput[["userInput"]])){
        return()
      }else{
        overlaps <- generateOverlapObject()
        makeVennDiagram(overlaps, NameOfPeaks = input$vennPlotOption,
                        fill=vennPlotOption_Color(), 
                        col=vennPlotOption_Color(), 
                        cat.col=vennPlotOption_Color(),
                        ) 
      }
      })
    })
    
  })
  
  # Distribution Plot:naccording to genome
  ######################################################### 
  output$ChIPpeakAnno_DuplicateDistribution <- renderPlot({
    
    input$duplicateDistributionPlot
    if (input$duplicateDistributionPlot == 0)return()
    isolate({
      
      peakSummary <- peakData()
      
      withProgress(message = 'Making plot', value = 0.3,{
        # Judge whether user choose right number
        if(length(input$duplicateDistribution_Option)<2){
          return()
        }else{
          peaks <- GRangesList(peakSummary$ChIPpeakAnno_Replicate[input$duplicateDistribution_Option])
          genomicElementDistribution(peaks, 
                                     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                     promoterRegion=c(upstream=2000, downstream=500),
                                     geneDownstream=c(upstream=0, downstream=2000))
        }
      })
    })
    
  })
  
  output$duplicateDistribution_Option <- renderUI({
    peakSummary <- peakData()
    peakName <- peakSummary$peakName
    
    checkboxGroupInput("duplicateDistribution_Option", 
                       h4("Choose no less than 2"),
                       choices = peakName,
                       selected = 2)
    
  })
  
  # Distribution Plot: according to overlaps
  #########################################################
  # Overlap Option: This option will return the combination of different overlaps
  output$overlapOption <- renderUI({
    temp <- vector()
    overlaps <- generateOverlapObject()
    peakNameOverlap <- names(overlaps$peaklist)
    for(i in 1:length(peakNameOverlap)){
      if(grepl("///",peakNameOverlap[i])){
        temp <- c(temp,peakNameOverlap[i])
      }
    }
    peakNameOverlap<-vector()
    peakNameOverlap <- temp
    selectInput("overlapSelection", h4("Select overlap sample"), 
                choices = peakNameOverlap, selected = 1 )
  })
  
  output$bindingSiteDistribution <- renderPlot({
    input$Distribution_of_Peaks_Around_Transcript_Start_Sites
    if (input$Distribution_of_Peaks_Around_Transcript_Start_Sites == 0)return()
    isolate({
      withProgress(message = 'Making plot', value = 0.3,{

    overlaps <- generateOverlapObject()
    overlaps_plot <- overlaps$peaklist[[input$overlapSelection]]
    if(is.null(overlaps_plot)){
      return()
    }else{
    binOverFeature(overlaps_plot, annotationData=annoData,
                   radius=input$DPATSS_radius, nbins=input$DPATSS_nbins, FUN=length, errFun=0,
                   xlab="distance from TSS (bp)", ylab="count", 
                   main="Distribution of aggregated peak numbers around TSS")
    }
      })
    })
  })
  
  output$bindingSiteDistribution_Message <- renderText({
    overlaps <- generateOverlapObject()
    overlaps_plot <- overlaps$peaklist[[input$overlapSelection]]
    if(is.null(overlaps_plot)){
      p("Please use 'Option' to plot Binding Site Distribution plot")
    }
    
  })
  
  # Peak distribution over different genomic features. 
  #########################################################
  output$ChIPpeakAnno_OverlapDistribution <- renderPlot({
    overlaps <- generateOverlapObject()
    overlaps_plot <- overlaps$peaklist[[input$overlapSelection]]
    ## check the genomic element distribution for the overlaps
    ## the genomic element distribution will indicates the 
    ## the best methods for annotation.
    ## The percentages in the legend show the percentage of peaks in 
    ## each category.
    genomicElementDistribution(overlaps_plot, 
                                      TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      promoterRegion=c(upstream=2000, downstream=500),
                                      geneDownstream=c(upstream=0, downstream=2000),
                                      promoterLevel=list(
                                        # from 5' -> 3', fixed precedence 3' -> 5'
                                        breaks = c(-2000, -1000, -500, 0, 500),
                                        labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", 
                                                   "upstream <500b", "TSS - 500b"),
                                        colors = c("#FFE5CC", "#FFCA99", 
                                                   "#FFAD65", "#FF8E32")))
  })
  
  # Peak Annotation
  #########################################################
  ChIPpeakAnno_PeakAnnotation <- reactive({
    overlaps <- generateOverlapObject()
    overlaps_plot <- overlaps$peaklist[[input$overlapSelection]]
    overlaps.anno <- annotatePeakInBatch(overlaps_plot, 
                                         AnnotationData=annoData, 
                                         output="nearestBiDirectionalPromoters",
                                         bindingRegion=c(-2000, 500))
    overlaps.anno <- addGeneIDs(overlaps.anno,
                                "org.Hs.eg.db",
                                IDs2Add = "entrez_id")
    return(overlaps.anno)
  })
  
  output$ChIPpeakAnno_PeakAnnotation_Table <- renderPrint({
    overlaps.anno<- ChIPpeakAnno_PeakAnnotation()
    head(overlaps.anno)
  })
  
  output$ChIPpeakAnno_PeakAnnotation_Pie <- renderPlot({
    overlaps.anno<- ChIPpeakAnno_PeakAnnotation()
    pie1(table(overlaps.anno$insideFeature))
  })
  
}


shinyApp(ui = ui, server = server)
