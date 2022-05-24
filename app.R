####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","gtsummary","tidyr",
              "dplyr","DT","ggplot2","ggpubr","tibble","survival","pheatmap",
              "readr","shinycssloaders","survminer","gridExtra")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("GSVA","clusterProfiler")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))


####----Input----####

ProjectName <- "Pan ICI Checkpoint Atlas"

ExpressionMatrix_file <- "Pan_ICI_iAtlas_ExpressionMatrix.zip"

MetaData_file <- "Pan_ICI_iAtlas_MetaData.txt"

MetaParam_File <- "Pan_ICI_iAtlas_MetaData_Params.txt"

GeneSet_File <- "GeneSet_List_HS.RData"





####----Read In Files----####

##--Meta--##
meta <- as.data.frame(read_delim(MetaData_file,delim = '\t', col_names = T))

MetaParam <- as.data.frame(read_delim(MetaParam_File,delim = '\t',col_names = F))
## Subset meta columns by category
metacol_samplenames <- MetaParam[which(MetaParam[,2] == "SampleName"),1]
metacol_feature <- MetaParam[which(MetaParam[,2] == "Feature"),1]
metacol_description <- MetaParam[which(MetaParam[,2] == "Description"),1]
if (("SampleType" %in% MetaParam[,2]) == TRUE) {
  metacol_sampletype <- MetaParam[which(MetaParam[,2] == "SampleType"),1]
}
if (("SampleType" %in% MetaParam[,2]) == FALSE) {
  metacol_sampletype <- NULL
}
metacol_survtime <- MetaParam[which(MetaParam[,2] == "SurvivalTime"),1]
metacol_survid <- MetaParam[which(MetaParam[,2] == "SurvivalID"),1]

colnames(meta)[which(colnames(meta) == metacol_samplenames)] <- "SampleName"
meta <- meta %>%
  relocate(SampleName)

##--Expression--##
expr <- as.data.frame(read_delim(ExpressionMatrix_file,delim = '\t', col_names = T))
rownames(expr) <- expr[,1]
expr <- expr[,-1]

##--Gene Set--##
# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Read Gene set File provides GMT file
ext <- tools::file_ext(GeneSet_File)
if (ext == "gmt") {
  gmt <- read.gmt(GeneSet_File)
  colnames(gmt) <- c("term","gene")
  gs <- list()
  for (i in unique(gmt[,1])){
    gs[[i]] <- gmt[gmt[,1] == i,]$gene
  }
}
# If user provides RData list file
if (ext == "RData") {
  gs <- loadRData(GeneSet_File)
}
# If user provides tab-delim two-col file
if  (ext == "tsv" | ext == "txt"){
  gmt <- as.data.frame(read_delim(GeneSet_File, delim = '\t'))
  colnames(gmt) <- c("term","gene")
  gs <- list()
  for (i in unique(gmt[,1])){
    gs[[i]] <- gmt[gmt[,1] == i,]$gene
  }
}




# Gene - "Gene Set"
exprGenes <- rownames(expr)
GeneGS_table <- data.frame(Genes = exprGenes)
# Gene Set Table
if (length(gs) == 94973) {
  GeneSetTable_File <- "GeneSet_CatTable.tsv"
  
}
if (exists("GeneSetTable_File") == TRUE) {
  if (file.exists(GeneSetTable_File) == TRUE) {
    
    GeneSetTable <- as.data.frame(read_delim(GeneSetTable_File, delim = '\t', col_names = T))
    gsTab = TRUE
    
  }
  if (file.exists(GeneSetTable_File) == FALSE) {
    
    GeneSets <- names(gs)
    GeneSetTable <- data.frame(GeneSets)
    gsTab = FALSE
    
  }
}
if (exists("GeneSetTable_File") == FALSE) {
  
  GeneSets <- names(gs)
  GeneSetTable <- data.frame(GeneSets)
  gsTab = FALSE
  
}


####----Functions----####

## Quartile Conversion
quartile_conversion = function(mat) {
  new_mat = mat;
  new_mat[mat <=  quantile(mat)[2]] = "Q1_Low";
  new_mat[mat > quantile(mat)[2] & mat <= quantile(mat)[3]] = "Q2_MedLow";
  new_mat[mat > quantile(mat)[3] & mat <= quantile(mat)[4]] = "Q3_MedHigh";
  new_mat[mat > quantile(mat)[4]] = "Q4_High";
  return (new_mat)
}

## High-Low
highlow = function(mat) {
  new_mat = mat;
  new_mat[mat > quantile(mat)[3]] = "High_AboveMedian";
  new_mat[mat <= quantile(mat)[3]] = "Low_BelowMedian";
  return (new_mat)
}

quantile_conversion = function(mat,cutoff) {
  new_mat = mat;
  new_mat[mat >= quantile(mat,1-cutoff)] = "High_AboveCutoff";
  new_mat[mat <= quantile(mat,cutoff)] = "Low_BelowCutoff";
  new_mat[mat > quantile(mat,cutoff) & mat < quantile(mat,1-cutoff)] = "BetweenCutoff";
  return (new_mat)
}



####----UI----####

ui <-
  navbarPage(paste("{ ",ProjectName," Survival Analysis }",sep = ""),
             
             ####----Overall Survival Tab----####
             
             tabPanel("Survival Analysis",
                      fluidPage(
                        sidebarLayout(
                          sidebarPanel(
                            tabsetPanel(
                              id = "survside",
                              tabPanel("Sample Parameters",
                                       p(),
                                       uiOutput("rendSampleTypeSelection"),
                                       uiOutput("rendFeatureSelection"),
                                       uiOutput("rendSubFeatureSelection"),
                                       uiOutput("rendScoreMethodBox"),
                                       tabsetPanel(
                                         id = "GeneSetTabs",
                                         tabPanel("Gene Sets",
                                                  uiOutput("rendGeneSetTable"),
                                                  value = 1
                                         ),
                                         tabPanel("Single Genes",
                                                  uiOutput("rendGeneGeneSetTable"),
                                                  value = 2
                                         ),
                                         tabPanel("User Gene Set",
                                                  fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData")),
                                                  uiOutput("rendUserGeneSetTable"),
                                                  value = 3
                                         )
                                       )
                              ),
                              tabPanel("Survival Parameters",
                                       p(),
                                       h3("Survival Data Type"),
                                       column(6,
                                              uiOutput("rendSurvivalType_time")
                                       ),
                                       column(6,
                                              uiOutput("rendSurvivalType_id")
                                       ),
                                       h3("Quantile Survival Plot Parameter"),
                                       numericInput("QuantPercent","High Risk Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                       h3("Survival Box Plot & Heatmap Parameters"),
                                       fluidRow(
                                         column(6,
                                                numericInput("cutoffTime1","High-Risk Survival Time Cutoff:", value = 364, min = 0, step = 1),
                                                selectInput("survStatus1","Survival Status Below Cutoff:", choices = c("1","0","Both"), selected = "1")
                                         ),
                                         column(6,
                                                numericInput("cutoffTime0","Low-Risk Survival Time Cutoff:", value = 365, min = 0, step = 1),
                                                selectInput("survStatus0","Survival Status Above Cutoff:", choices = c("1","0","Both"), selected = "0")
                                         )
                                       )
                              ),
                              tabPanel("Figure Parameters",
                                       p(),
                                       h3("Boxplot Parameters"),
                                       fluidRow(
                                         column(6,
                                                numericInput("boxplotFont","Boxplot Font Size:", value = 15, step = 1)),
                                         column(6,
                                                numericInput("boxplotDot", "Boxplot Dot Size:", value = 0.75, step = 0.25))
                                       ),
                                       uiOutput("rendBoxoptselec"),
                                       h3("Heatmap Parameters"),
                                       selectInput("ClusterMethod", "Select Cluster Method",
                                                   choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                       fluidRow(
                                         column(6,
                                                numericInput("heatmapFontR", "Heatmap Row Font Size:", value = 9, step = 1)),
                                         column(6,
                                                numericInput("heatmapFontC", "Heatmap Column Font Size:", value = 10, step = 1))
                                       )
                              ),
                              tabPanel("Meta Data",
                                       p(),
                                       uiOutput("rendMetaTableCols"),
                                       uiOutput("rendMetaTable"),
                                       fluidRow(
                                         column(6,
                                                uiOutput("DnldMetaButon")
                                         ),
                                         column(6,
                                                uiOutput("DnldExprButon"))
                                       )
                              )
                            )
                          ),
                          mainPanel(
                            tabsetPanel(
                              id = "surv",
                              tabPanel("Survival Plot",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput("Splot", width = "100%", height = "500px")), type = 6),
                                       fluidRow(
                                         column(3,
                                                offset = 1,
                                                checkboxInput("ShowQaurtHR","Display Quartile Hazard Ratio Table")),
                                         column(4,
                                                uiOutput("rendQuartHRtab"))
                                       ),
                                       fluidRow(
                                         tags$head(
                                           tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
                                         hr()
                                       ),
                                       withSpinner(jqui_resizable(plotOutput("SplotBIN", width = "100%", height = "500px")), type = 6),
                                       fluidRow(
                                         column(3,
                                                offset = 1,
                                                checkboxInput("ShowBINHR","Display Binary Hazard Ratio Table")),
                                         column(4,
                                                uiOutput("rendBINHRtab"))
                                       ),
                                       fluidRow(
                                         tags$head(
                                           tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
                                         hr()
                                       ),
                                       withSpinner(jqui_resizable(plotOutput("SquantPlot", width = "100%", height = "500px")), type = 6),
                                       fluidRow(
                                         column(3,
                                                offset = 1,
                                                checkboxInput("ShowQuantHR","Display Quantile Hazard Ratio Table")),
                                         column(4,
                                                uiOutput("rendQuantHRtab"))
                                       ),
                                       fluidRow(
                                         tags$head(
                                           tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
                                         hr()
                                       ),
                                       value = 1),
                              tabPanel("Survival Boxplot",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput("Sboxplot", width = "100%", height = "500px")), type = 6),
                                       div(DT::dataTableOutput("SboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                       p(),
                                       downloadButton("dnldSBoxplotTab","Download Table"),
                                       value = 2),
                              tabPanel("Survival Heatmap",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
                                       downloadButton("dnldSheatmapexpr","Download Expression Matrix From Heatmap"),
                                       value = 3),
                              tabPanel("Feature Boxplot",
                                       p(),
                                       selectInput("BoxplotFeature","Select Feature to View in Boxplot",
                                                   choices = metacol_feature),
                                       withSpinner(jqui_resizable(plotOutput("Featureboxplot", width = "100%", height = "500px")), type = 6),
                                       div(DT::dataTableOutput("FeatureboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                       p(),
                                       downloadButton("dnldFeatureboxplotTab","Download Table"),
                                       value = 4),
                              tabPanel("Feature Heatmap",
                                       p(),
                                       selectInput("HeatmapFeature","Select Feature to View in Heatmap",
                                                   choices = metacol_feature),
                                       withSpinner(jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
                                       downloadButton("dnldFheatmapexpr","Download Expression Matrix From Heatmap"),
                                       value = 5)
                            )
                          )
                        )
                      )
             )
  )



####----Server----####


server <- function(input, output, session) {
  
  ## Select survival type selection based on meta columns (ex. OS vs EFS)
  output$rendSurvivalType_time <- renderUI({
    
    ## Only show if more than one option
    if (length(metacol_survtime > 1)) {
      
      selectInput("SurvivalType_time","Select Survival Time Data:", choices = metacol_survtime)
      
    }
    
  })
  
  ## Select survival type selection based on meta columns (ex. OS vs EFS)
  output$rendSurvivalType_id <- renderUI({
    
    ## Only show if more than one option
    if (length(metacol_survid > 1)) {
      
      selectInput("SurvivalType_id","Select Survival ID Data:", choices = metacol_survid)
      
    }
    
  })
  
  ## Select sample type to subset samples by - only render if more than one sample type
  output$rendSampleTypeSelection <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      #meta[,metacol_sampletype] = meta[,metacol_sampletype] %>% replace_na("NA")
      SampleTypeChoices <- unique(meta[,metacol_sampletype])
      SampleTypeChoices <- c(SampleTypeChoices,"All_Sample_Types")
      selectInput("SampleTypeSelection","Select Sample Type:", choices = SampleTypeChoices)
      
    }
    
  })
  
  ## Select primary feature to look at - All not working yet
  output$rendFeatureSelection <- renderUI({
    
    FeatureChoices <- c(metacol_feature,"All_Features")
    selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices)
    
  })
  
  ## Select primary features condition to look at - All not working yet
  output$rendSubFeatureSelection <- renderUI({
    
    req(input$FeatureSelection)
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleType <- input$SampleTypeSelection
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleType <- "All_Sample_Types"
    }
    #SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    
    if (SampleType == "All_Sample_Types") {
      meta <- meta
    }
    if (SampleType != "All_Sample_Types") {
      meta <- meta[which(meta[,metacol_sampletype] == SampleType),]
    }
    
    if (Feature != "All_Features") {
      
      SubFeatureChoices <- unique(meta[,Feature])
      # Sort options, will put 1,TRUE,yes before 0,FASLE,no, so the 'positive' value is the initial selected - puts NA last
      SubFeatureChoices <- sort(SubFeatureChoices, decreasing = T, na.last = T)
      selectInput("subFeatureSelection","Feature Condition:", choices = SubFeatureChoices)
      
    }
    
  })
  
  ## Select ssGSEA function scoring method
  output$rendScoreMethodBox <- renderUI({
    
    selectInput("ScoreMethod","Select Scoring Method",choices = c("ssgsea","gsva","zscore","plage"))
    
  })
  
  ## Select stat compare method for boxplots
  output$rendBoxoptselec <- renderUI({
    
    # if on tab of primary feature boxplot
    if (input$surv != 4) {
      
      selectInput("boxoptselec","Boxplot Stat Compare Method:", choices = c("wilcox.test","t.test","none"))
      
    }
    # if on tab of any feature boxplot
    else if (input$surv == 4) {
      
      selectInput("boxoptselec","Boxplot Stat Compare Method:", choices = c("none","wilcox.test","t.test","kruskal.test","anova"))
      
    }
    
  })
  
  ## Select column names of meta to view in sample box
  output$rendMetaTableCols <- renderUI({
    
    meta <- ssGSEAmeta()
    MetaColChoices <- colnames(meta)[c(2:ncol(meta))]
    selectInput("MetaTableCols","Select Meta Descriptors to View:", choices = MetaColChoices, selected = "", multiple = T)
    
  })
  
  ## View of meta table with choices selected
  output$rendMetaTable <- renderUI({
    
    div(DT::dataTableOutput("MetaTable"), style = "font-size:10px; height:450px; overflow-X: scroll; overflow-Y: scroll")
    
    
  })
  
  ## View Gene Set table
  output$rendGeneSetTable <- renderUI({
    
    div(DT::dataTableOutput("GeneSetTable"), style = "font-size:10px; height:450px; overflow-Y: scroll")
    
    
  })
  
  ## View Gene - "Gene Set" table
  output$rendGeneGeneSetTable <- renderUI({
    
    div(DT::dataTableOutput("geneGeneSetTable"), style = "font-size:10px; height:450px; overflow-Y: scroll")
    
  })
  
  ## View User Gene Set table
  output$rendUserGeneSetTable <- renderUI({
    
    req(input$userGeneSet)
    div(DT::dataTableOutput("userGeneSetTable"), style = "font-size:10px; height:450px; overflow-Y: scroll")
    
  })
  
  ## Render Meta Table
  output$MetaTable <- renderDataTable({
    
    # User selections
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    SubFeature <- input$subFeatureSelection
    meta <- ssGSEAmeta()
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    # Initial columns selected is NULL
    if (is.null(input$MetaTableCols) == TRUE) {
      
      userMetaCols <- NULL
      
    }
    # When meta tab is selected
    else if (is.null(input$MetaTableCols) == FALSE) {
      
      # If columns selected, add to list
      if (input$MetaTableCols != "") {
        userMetaCols <- input$MetaTableCols
      }
      else if (input$MetaTableCols == "") {
        userMetaCols <- NULL
      }
      
    }
    
    metaCols <- colnames(meta)[1] #select sample name column automatically
    userMetaCols <- userMetaCols[userMetaCols != Feature] #remove condition column from user selection because it is automatically added
    metaCols <- c(metaCols,surv_time_col,surv_id_col,Feature,userMetaCols) #combine column names selected
    meta_sub <- meta[,metaCols]
    DT::datatable(meta_sub,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 fixedColumns = list(leftColumns = 1)),
                  rownames = F)
    
  })
  
  ## Render Gene Set Selection Table
  output$GeneSetTable <- renderDataTable({
    
    DT::datatable(GeneSetTable,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  ## Render Gene - "Gene Set" selection table
  output$geneGeneSetTable <- renderDataTable({
    
    DT::datatable(GeneGS_table,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  ## Render User Gene Set Table - Backend
  userGeneSetTable_backend <- reactive({
    
    gs.u <- input$userGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
    
    # If user provides GMT file
    if (ext == "gmt") {
      gmt <- read.gmt(gs.u$datapath)
      uGS_table <- as.data.frame(unique(gmt[,1]))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    
    # If user provides RData list file
    else if (ext == "RData") {
      gs_u <- loadRData(gs.u$datapath)
      uGS_table <- as.data.frame(names(gs_u))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    
    # If user provides tab-delim two-col file
    else {
      gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
      uGS_table <- as.data.frame(unique(gmt[,1]))
      colnames(uGS_table)[1] <- "GeneSet"
    }
    uGS_table
    
  })
  
  ## Render User Gene Set Table
  output$userGeneSetTable <- renderDataTable({
    
    req(input$userGeneSet)
    uGS_table <- userGeneSetTable_backend()
    DT::datatable(uGS_table,
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 10,
                                 scrollX = T),
                  selection = list(mode = 'single', selected = 1),
                  rownames = F)
    
  })
  
  ## User gene set list reactive
  user_gs <- reactive({
    
    gs.u <- input$userGeneSet
    ext <- tools::file_ext(gs.u$datapath)
    req(gs.u)
    validate(need(ext == c("gmt","tsv","txt", "RData"), "Please upload .gmt, .tsv, .txt, or .RData file"))
    
    # If user provides GMT file
    if (ext == "gmt") {
      gmt <- read.gmt(gs.u$datapath)
      colnames(gmt) <- c("term","gene")
      gs_u <- list()
      for (i in unique(gmt[,1])){
        gs_u[[i]] <- gmt[gmt[,1] == i,]$gene
      }
      gs_u
    }
    
    # If user provides RData list file
    else if (ext == "RData") {
      gs_u <- loadRData(gs.u$datapath)
    }
    
    # If user provides tab-delim two-col file
    else {
      gmt <- as.data.frame(read_delim(gs.u$datapath, delim = '\t'))
      colnames(gmt) <- c("term","gene")
      gs_u <- list()
      for (i in unique(gmt[,1])){
        gs_u[[i]] <- gmt[gmt[,1] == i,]$gene
      }
      gs_u
    }
    
    gs_u
    
  })
  
  ## Reactive to represent the chosen gene set
  gs_react <- reactive({
    
    if (input$GeneSetTabs == 1) {
      if (gsTab == TRUE) {
        geneset <- gs[(GeneSetTable[input$GeneSetTable_rows_selected,3])]
      }
      if (gsTab == FALSE) {
        geneset <- gs[(GeneSetTable[input$GeneSetTable_rows_selected,1])]
      }
    }
    if (input$GeneSetTabs == 2) {
      gene <- GeneGS_table[input$geneGeneSetTable_rows_selected,1]
      geneset <- list(gene = gene)
      names(geneset)[1] <- gene
    }
    if (input$GeneSetTabs == 3) {
      req(input$userGeneSet)
      gs_u <- user_gs()
      uGS_table <- userGeneSetTable_backend()
      geneset <- gs_u[(uGS_table[input$userGeneSetTable_rows_selected,1])]
    }
    geneset
    
  })
  
  ## Download button for subset meta
  output$DnldMetaButon <- renderUI({
    
    downloadButton("dnldMeta", "Download Meta Subset")
    
  })
  
  ## Download handler for meta
  output$dnldMeta <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_meta.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_meta.txt",sep = "")
      }
    },
    content = function(file) {
      meta <- ssGSEAmeta()
      write_delim(meta,file,delim = '\t')
    }
  )
  
  ## Download button for subset expression
  output$DnldExprButon <- renderUI({
    
    downloadButton("dnldExpr", "Download Expression Subset")
    
  })
  
  ## Download handler for expression
  output$dnldExpr <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      # Reformat to make sure genes show in file
      expr$Gene <- rownames(expr)
      expr <- expr %>%
        relocate(Gene)
      write_delim(expr,file,delim = '\t')
    }
  )
  
  ## Download handler for expression
  output$dnldSheatmapexpr <- downloadHandler(
    filename = function() {
      # Variables
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"SurvivalCutoff_expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"SurvivalCutoff_expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      GSgenes <- geneset[[geneset_name]]
      # Include only genes from gene set
      expr <- expr[which(rownames(expr) %in% GSgenes),]
      # Reformat to make sure genes show in file
      expr$Gene <- rownames(expr)
      expr <- expr %>%
        relocate(Gene)
      write_delim(expr,file,delim = '\t')
    }
  )
  
  ## Download handler for expression
  output$dnldFheatmapexpr <- downloadHandler(
    filename = function() {
      # Variables
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      Feature2 <- input$HeatmapFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      # Make file name
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"Featuring_",Feature2,"_expr.txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"Featuring_",Feature2,"_expr.txt",sep = "")
      }
    },
    content = function(file) {
      expr <- exprSub()
      geneset <- gs_react()                 #Chosen Gene Set
      geneset_name <- names(geneset)        #Name of chosen gene set
      GSgenes <- geneset[[geneset_name]]
      # Include only genes from gene set
      expr <- expr[which(rownames(expr) %in% GSgenes),]
      # Reformat to make sure genes show in file
      expr$Gene <- rownames(expr)
      expr <- expr %>%
        relocate(Gene)
      write_delim(expr,file,delim = '\t')
    }
  )
  
  
  
  ## Meta subset reactive - "All Primary features" not working yet
  metaSub <- reactive({
    
    req(input$FeatureSelection)
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleType <- input$SampleTypeSelection
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleType <- "All_Sample_Types"
    }
    Feature <- input$FeatureSelection
    SubFeature <- input$subFeatureSelection
    
    if (SampleType == "All_Sample_Types") {
      meta <- meta
    }
    if (SampleType != "All_Sample_Types") {
      meta <- meta[which(meta[,metacol_sampletype] == SampleType),]
    }
    
    if (Feature != "All_Features") {
      meta <- meta[which(meta[,Feature] == SubFeature),]
    }
    if (Feature == "All_Features") {
      meta <- meta
    }
    meta
    
  })
  
  ## Expression subset reactive
  exprSub <- reactive({
    
    meta <- metaSub()
    samples <- meta[,1]
    expr <- expr[,which(colnames(expr) %in% samples), drop = F] #Subset by sample names in subset meta table
    expr
    
  })
  
  ## Perform ssGSEA and functions on new meta table
  ssGSEAmeta <- reactive({
    
    meta <- metaSub()                     #Subset meta
    expr <- exprSub()                     #Subset expression
    geneset <- gs_react()                 #Chosen Gene Set
    geneset_name <- names(geneset)        #Name of chosen gene set
    scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
    quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    
    ## Remove rows with NA in survival column
    meta <- meta[!is.na(meta[,surv_time_col]),]
    meta <- meta[!is.na(meta[,surv_id_col]),]
    
    ## Re-subset expression matrix
    samples <- meta[,1]
    expr_mat <- as.matrix(expr[,samples])
    
    #expr_mat <- as.matrix(expr)
    
    ## Perform ssGSEA with gs and new subset data
    ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod)
    
    ## Transform
    ssGSEA <- as.data.frame(t(ssGSEA))
    ssGSEA$SampleName <- rownames(ssGSEA)
    
    ## Perform further functions
    ssGSEA$VAR_Q <- quartile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
    ssGSEA$Quartile <- paste("", ssGSEA$VAR_Q, sep="")
    ssGSEA$BIN <- highlow(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
    ssGSEA$Quantile <- quantile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff)
    
    ## Merge with meta
    meta_ssGSEA <- merge(meta,ssGSEA, by = "SampleName", all = T)
    meta_ssGSEA
    
  })
  
  ## Survival Plot - Quartile
  output$Splot <- renderPlot({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    expr <- exprSub()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Remove rows with NA in survival column
    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
    
    ## Subset columns needed for plot and rename for surv function
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"Quartile")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ Quartile, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients\n", geneset_name," (",scoreMethod," score in quartiles)", sep = ""),
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer", pval=TRUE, ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain")
    )
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
  })
  
  output$rendQuartHRtab <- renderUI({
    
    if (input$ShowQaurtHR == T) {
      div(withSpinner(tableOutput("SQuartileHRtab"), type = 7, size = 0.5), style = "font-size:10px")
    }    
    
  })
  output$rendBINHRtab <- renderUI({
    
    if (input$ShowBINHR == T) {
      div(withSpinner(tableOutput("SBinaryHRtab"), type = 7, size = 0.5), style = "font-size:10px")
    }    
    
  })
  output$rendQuantHRtab <- renderUI({
    
    if (input$ShowQuantHR == T) {
      div(withSpinner(tableOutput("SQuantileHRtab"), type = 7, size = 0.5), style = "font-size:10px")
    }    
    
  })
  
  output$SQuartileHRtab <- renderTable({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    expr <- exprSub()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Remove rows with NA in survival column
    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
    
    ## Subset columns needed for plot and rename for surv function
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"Quartile")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    tab <- coxph(Surv(OS.time,OS.ID) ~ Quartile, data = meta_ssgsea_sdf) %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  ## Survival Plot - Binary
  output$SplotBIN <- renderPlot({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    expr <- exprSub()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Remove rows with NA in survival column
    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"BIN")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ BIN, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients\n",
                                       geneset_name," (",scoreMethod, " Median Cutoff)", sep = ""),
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer", pval=TRUE, ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain")
    )
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
  })
  
  output$SBinaryHRtab <- renderTable({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    expr <- exprSub()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Remove rows with NA in survival column
    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
    
    ## Subset columns needed for plot and rename for surv function
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"BIN")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    tab <- coxph(Surv(OS.time,OS.ID) ~ BIN, data = meta_ssgsea_sdf) %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  ## Survival Plot - Quantile
  output$SquantPlot <- renderPlot({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    quantCutoffOG <- input$QuantPercent
    quantCutoff <- input$QuantPercent/100
    labelQuantCutoff <- paste(quantCutoffOG,"%",sep = "")
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    expr <- exprSub()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Remove rows with NA in survival column
    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
    
    ## Subset columns needed for plot
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"Quantile")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Remove between cutoff samples
    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$Quantile != "BetweenCutoff"),]
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ Quantile, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    ## Generate plot
    ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                         title = paste("Survival curves of ", Feature,SampleTypeLab,"\n",
                                       geneset_name," (Top (Bottom) ",labelQuantCutoff," Patients Split based on ",scoreMethod," score)", sep = ""),
                         xscale = c("d_y"),
                         break.time.by=365.25,
                         xlab = "Years", 
                         ylab = paste(SurvDateType,"Survival Probability"),
                         submain = "Based on Kaplan-Meier estimates",
                         caption = "created with survminer", pval=TRUE, ggtheme = theme_bw(),
                         font.title = c(16, "bold"),
                         font.submain = c(12, "italic"),
                         font.caption = c(12, "plain"),
                         font.x = c(14, "plain"),
                         font.y = c(14, "plain"),
                         font.tickslab = c(12, "plain")
    )
    
    ggsurv$table <- ggsurv$table + theme_cleantable()
    ggsurv
    
  })
  
  output$SQuantileHRtab <- renderTable({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    expr <- exprSub()
    meta_ssgsea <- ssGSEAmeta()
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Remove rows with NA in survival column
    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
    
    ## Subset columns needed for plot and rename for surv function
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"Quantile")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Remove between cutoff samples
    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$Quantile != "BetweenCutoff"),]
    
    ## Survival Function
    tab <- coxph(Surv(OS.time,OS.ID) ~ Quantile, data = meta_ssgsea_sdf) %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  ## Boxplot reactive to determine hig/low risk samples based on user input
  SboxplotReact <- reactive({
    
    ## Assing variables
    cutoff_1 <- input$cutoffTime1                   #High-risk cutoff time
    cutoff_0 <- input$cutoffTime0                   #Low-risk cutoff time
    OS_choice_1 <- input$survStatus1                #High-risk survival status
    OS_choice_0 <- input$survStatus0                #Low-risk survival status
    ssGSEA_meta <- ssGSEAmeta()                     #Meta with ssGSEA scores and function calculations
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    
    # If both survival choices are numeric (1 or 0)
    if (is.na(as.numeric(OS_choice_1)) == F & is.na(as.numeric(OS_choice_0)) == F) {
      ssGSEA_meta <- ssGSEA_meta %>%
        mutate(SurvivalCutoff = case_when(
          ssGSEA_meta[,surv_time_col] <= cutoff_1 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_1) ~ "High-Risk [Below Survival Time Cutoff]",
          ssGSEA_meta[,surv_time_col] >= cutoff_0 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_0) ~ "Low-Risk [Above Survival Time Cutoff]"
        ))
    }
    # If either choice selects "Both" as survival outcome
    if (is.na(as.numeric(OS_choice_1)) == T | is.na(as.numeric(OS_choice_0)) == T) {
      # if both & (1 or 0)
      if (is.na(as.numeric(OS_choice_1)) == T & is.na(as.numeric(OS_choice_0)) == F) {
        ssGSEA_meta <- ssGSEA_meta %>%
          mutate(SurvivalCutoff = case_when(
            ssGSEA_meta[,surv_time_col] <= cutoff_1 ~ "High-Risk [Below Survival Time Cutoff]",
            ssGSEA_meta[,surv_time_col] >= cutoff_0 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_0) ~ "Low-Risk [Above Survival Time Cutoff]"
          ))
      }
      # if (1 or 0) and both
      else if (is.na(as.numeric(OS_choice_1)) == F & is.na(as.numeric(OS_choice_0)) == T) {
        ssGSEA_meta <- ssGSEA_meta %>%
          mutate(SurvivalCutoff = case_when(
            ssGSEA_meta[,surv_time_col] <= cutoff_1 & ssGSEA_meta[,surv_id_col] == as.numeric(OS_choice_1) ~ "High-Risk [Below Survival Time Cutoff]",
            ssGSEA_meta[,surv_time_col] >= cutoff_0 ~ "Low-Risk [Above Survival Time Cutoff]"
          ))
      }
      # if both and both
      else if (is.na(as.numeric(OS_choice_1)) == T & is.na(as.numeric(OS_choice_0)) == T) {
        ssGSEA_meta <- ssGSEA_meta %>%
          mutate(SurvivalCutoff = case_when(
            ssGSEA_meta[,surv_time_col] <= cutoff_1 ~ "High-Risk [Below Survival Time Cutoff]",
            ssGSEA_meta[,surv_time_col] >= cutoff_0 ~ "Low-Risk [Above Survival Time Cutoff]"
          ))
      }
    }
    
    ssGSEA_meta <- ssGSEA_meta[which(is.na(ssGSEA_meta$SurvivalCutoff) == F),]
    ssGSEA_meta$SurvivalCutoff <- factor(ssGSEA_meta$SurvivalCutoff,
                                         levels = c("High-Risk [Below Survival Time Cutoff]","Low-Risk [Above Survival Time Cutoff]"))
    ssGSEA_meta
    
  })
  
  output$Sboxplot <- renderPlot({
    
    ## Assign Variables
    ssGSEA_meta <- SboxplotReact()
    geneset <- gs_react()
    GeneSet <- names(geneset)
    Feature <- input$FeatureSelection
    dot <- input$boxplotDot
    font <- input$boxplotFont
    SampleType <- input$SampleTypeSelection
    scoreMethod <- input$ScoreMethod
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    
    ## Determine type of survival data - OS/EFS/PFS?
    SurvDateType <- sub("\\..*","",surv_time_col)
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    ggplot(ssGSEA_meta, aes(factor(SurvivalCutoff), ssGSEA_meta[,GeneSet], fill = SurvivalCutoff)) +
      geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
      geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
      labs(x = "Group", y = paste(GeneSet," ",scoreMethod, " Score", sep = ""),
           title = paste(GeneSet," Gene Set ",scoreMethod," Score: ",Feature,SampleTypeLab,"Patients",sep = "")) +
      theme_bw() +
      stat_compare_means(method = input$boxoptselec) +
      theme(text = element_text(size = font))
    
    
  })
  
  output$SboxplotTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    ssGSEA_meta <- SboxplotReact()
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    boxTab <- ssGSEA_meta[,c("SampleName",surv_time_col,surv_id_col,GeneSet,"SurvivalCutoff")]
    DT::datatable(boxTab,
                  options = list(paging = F),
                  rownames = F)
    
  })
  
  output$dnldSBoxplotTab <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      geneset <- gs_react()
      GeneSet <- names(geneset)
      method <- input$ScoreMethod
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",GeneSet,"_",method,".txt",sep = "")
    },
    content = function(file) {
      ssGSEA_meta <- SboxplotReact()
      geneset <- gs_react()
      GeneSet <- names(geneset)
      if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
        surv_time_col <- metacol_survtime[1]
        surv_id_col <- metacol_survid[1]
      }
      if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
        surv_time_col <- input$SurvivalType_time
        surv_id_col <- input$SurvivalType_id
      }
      boxTab <- as.data.frame(ssGSEA_meta[,c("SampleName",surv_time_col,surv_id_col,GeneSet,"SurvivalCutoff")])
      write_delim(boxTab,file,delim = '\t')
    }
  )
  
  output$Sheatmap <- renderPlot({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    heatgenes <- geneset[[GeneSet]]
    meta <- SboxplotReact()
    expr_start <- exprSub()
    expr <- expr_start[which(rownames(expr_start) %in% heatgenes),meta$SampleName]
    clmethod <- input$ClusterMethod
    rowfont <- input$heatmapFontR
    colfont <- input$heatmapFontC
    
    dataset <- expr
    dataset <- log2(dataset + 1)
    zdataset <- apply(dataset, 1, scale)
    zdataset <- apply(zdataset, 1, rev)
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    meta2 <- meta[order(meta[,2]),]
    meta3 <- meta2[order(meta2$SurvivalCutoff),]
    samporder <- meta3$SampleName
    dataset2 <- dataset[,samporder]
    type <- meta3$SurvivalCutoff
    meta4 <- as.data.frame(type)
    rownames(meta4) <- meta3[,1]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    pheatmap(dataset2,
             cluster_col = F,
             cluster_row = T,
             fontsize_row = rowfont,
             fontsize_col = colfont,
             show_rownames = T ,
             show_colnames = T,
             annotation_col = meta4,
             clustering_method = clmethod,
             color=hmcols,
             border_color = NA)
    
    
  })
  
  ## Feature Boxplot
  output$Featureboxplot <- renderPlot({
    
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    dot <- input$boxplotDot
    font <- input$boxplotFont
    StatMethod <- input$boxoptselec
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$BoxplotFeature
    meta_ssGSEA <- ssGSEAmeta()
    boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
      geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
      geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
      labs(x = "Group", y = paste(GeneSet," ",scoreMethod," Score", sep = ""),
           title = paste(GeneSet," Gene Set ",scoreMethod," Score: ","\n",
                         Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
           fill = FeatureSelec) +
      theme_bw() +
      stat_compare_means(method = StatMethod) +
      theme(text = element_text(size = font))
    
  })
  
  ## Feature bloxplot table
  output$FeatureboxplotTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$BoxplotFeature
    meta_ssGSEA <- ssGSEAmeta()
    boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    DT::datatable(boxTab,
                  options = list(paging = F),
                  rownames = F)
    
  })
  
  output$dnldFeatureboxplotTab <- downloadHandler(
    filename = function() {
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      FeatureSelec <- input$BoxplotFeature
      geneset <- gs_react()
      GeneSet <- names(geneset)
      method <- input$ScoreMethod
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste("_",SampleType,"_",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- "_"
      }
      
      paste(gsub(" ","",ProjectName),SampleTypeLab,Feature,"_",FeatureSelec,"_",GeneSet,"_",method,".txt",sep = "")
    },
    content = function(file) {
      meta_ssGSEA <- ssGSEAmeta()
      FeatureSelec <- input$BoxplotFeature
      geneset <- gs_react()
      GeneSet <- names(geneset)
      boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
      write_delim(boxTab,file,delim = '\t')
    }
  )
  
  output$FeatureHeatmap <- renderPlot({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$HeatmapFeature
    meta_ssGSEA <- ssGSEAmeta()
    meta <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    
    heatgenes <- geneset[[GeneSet]]
    expr_start <- exprSub()
    expr <- expr_start[which(rownames(expr_start) %in% heatgenes),meta$SampleName, drop = F]
    clmethod <- input$ClusterMethod
    rowfont <- input$heatmapFontR
    colfont <- input$heatmapFontC
    
    dataset <- expr
    dataset <- as.data.frame(log2(dataset + 1))
    zdataset <- as.data.frame(apply(dataset, 1, scale))
    zdataset <- as.data.frame(apply(zdataset, 1, rev))
    colnames(zdataset) <- names(dataset)
    dataset <- as.matrix(zdataset)
    dataset[is.na(dataset)] <- 0
    dataset = dataset[apply(dataset[,-1], 1, function(x) !all(x==0)),]
    minimum = -5;
    maximum = 5;
    if (abs(min(dataset)) > abs(max(dataset))) {
      dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
    } else {
      dataset[dataset > abs(min(dataset))] = abs(min(dataset))
    }
    #meta2 <- meta[order(meta[,2]),]
    meta3 <- meta[order(meta[,FeatureSelec]),]
    meta3[,FeatureSelec] <- as.character(meta3[,FeatureSelec])
    meta3[,FeatureSelec] = meta3[,FeatureSelec] %>% replace_na("NA")
    samporder <- meta3$SampleName
    dataset2 <- dataset[,samporder]
    type <- meta3[,FeatureSelec]
    meta4 <- as.data.frame(type)
    rownames(meta4) <- meta3[,1]
    bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
    hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    pheatmap(dataset2,
             cluster_col = F,
             cluster_row = T,
             fontsize_row = rowfont,
             fontsize_col = colfont,
             show_rownames = T ,
             show_colnames = T,
             annotation_col = meta4,
             clustering_method = clmethod,
             color=hmcols,
             border_color = NA)
    
  })
  
  
  
  
}



# Run the application 
shinyApp(ui = ui, server = server)



