


####----User Input----####

ProjectName <- "Pan ICI Checkpoint Atlas"

ExpressionMatrix_file <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip"

MetaData_file <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt"

MetaParam_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt"


##--Advanced Set-Up--##

## Pre-Selected Inputs
# An option from the meta, All, or NULL
PreSelect_SamplyType <- NULL
PreSelect_Feature <- "All"
# An option from the meta or NULL
PreSelect_SubFeature <- NULL
PreSelect_SecondaryFeature <- NULL

# DO NOT CHANGE when using provided gene set data - only adjust file path if needed
GeneSet_File <- "GeneSet_Data/GeneSet_List.RData"
GeneSetTable_File <- "GeneSet_Data/GeneSet_Table.zip"




####----Install and load packages----####

## Check if Immune deconvolution package is installed
immudecon <- "immunedeconv"
immudecon_check <- immudecon %in% rownames(installed.packages())
if (immudecon_check == TRUE) {
  library(immunedeconv)
}
packages <- c("shiny","shinythemes","shinyjqui","gtsummary","tidyr","RColorBrewer",
              "dplyr","DT","ggplot2","ggpubr","tibble","survival","pheatmap","stringr",
              "plotly","readr","shinycssloaders","survminer","gridExtra","viridis")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
#bioconductor packages
bioCpacks <- c("GSVA")
installed_packages_BIOC <- bioCpacks %in% rownames(installed.packages())
if (any(installed_packages_BIOC == FALSE)) {
  BiocManager::install(bioCpacks[!installed_packages_BIOC], ask = F)
}
invisible(lapply(bioCpacks, library, character.only = TRUE))




####----Read In Files----####

##--Meta--##
meta <- as.data.frame(read_delim(MetaData_file,delim = '\t', col_names = T))

##--Meta Param--##
MetaParam <- as.data.frame(read_delim(MetaParam_File,delim = '\t',col_names = F))
MetaParam[,2] <- gsub(" ","",MetaParam[,2])
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
## Get surv time and ID columns - move OS to the front of the list
metacol_survtime <- MetaParam[which(MetaParam[,2] == "SurvivalTime"),1]
os_time <- grep("os",metacol_survtime, ignore.case = T, value = T)
metacol_survtime <- c(os_time,metacol_survtime[metacol_survtime != os_time])

metacol_survid <- MetaParam[which(MetaParam[,2] == "SurvivalID"),1]
os_id <- grep("os",metacol_survid, ignore.case = T, value = T)
metacol_survid <- c(os_id,metacol_survid[metacol_survid != os_id])

## Rename and move Sample Name column to the front
colnames(meta)[which(colnames(meta) == metacol_samplenames)] <- "SampleName"
meta <- meta %>%
  relocate(SampleName)
## Replace any special characters to make uniform with expression
meta[,1] <- gsub("[[:punct:]]","_",meta[,1])


##--Expression--##
expr <- as.data.frame(read_delim(ExpressionMatrix_file,delim = '\t', col_names = T))
rownames(expr) <- expr[,1]
expr <- expr[,-1]
## Replace any special characters to make uniform with expression
colnames(expr) <- gsub("[[:punct:]]","_",colnames(expr))



##--Gene Set--##
# R Data list load function for naming
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
gs <- loadRData(GeneSet_File)
# Gene - "Gene Set"
exprGenes <- rownames(expr)
GeneGS_table <- data.frame(Genes = exprGenes)
# Gene Set Table
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



## Pre-Selected Inputs
# An option from the meta, All, or NULL
if (is.null(PreSelect_SamplyType) == FALSE) {
  if (grepl("all",PreSelect_SamplyType, ignore.case = T) == TRUE) {
    PreSelect_SamplyType <- "All_Sample_Types"
  }
}
if (is.null(PreSelect_Feature) == FALSE) {
  if (grepl("all",PreSelect_Feature, ignore.case = T) == TRUE) {
    PreSelect_Feature <- "All_Features"
  }
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

quantile_conversion2 = function(mat,cutoff) {
  new_mat = mat;
  new_mat[mat > quantile(mat,cutoff)] = "High_AboveCutoff";
  new_mat[mat <= quantile(mat,cutoff)] = "Low_BelowCutoff";
  return (new_mat)
}


##--Perform Immune Deconvolution--##

if (immudecon_check == TRUE) {
  
  estimate_decon <- as.data.frame(deconvolute(expr, "estimate"))
  mcp_counter_decon <- as.data.frame(deconvolute(expr, "mcp_counter"))
  
  rownames(estimate_decon) <- estimate_decon[,1]
  rownames(mcp_counter_decon) <- mcp_counter_decon[,1]
  
  estimate_decon <- estimate_decon[,-1]
  mcp_counter_decon <- mcp_counter_decon[,-1]
  
  estimate_decon <- as.data.frame(t(estimate_decon))
  mcp_counter_decon <- as.data.frame(t(mcp_counter_decon))
  
  colnames(estimate_decon) <- paste(gsub(" ","_",colnames(estimate_decon)),"ESTIMATE",sep = "_")
  colnames(mcp_counter_decon) <- paste(gsub(" ","_",colnames(mcp_counter_decon)),"MCPcount",sep = "_")
  
  decon_score_cols <- c(colnames(estimate_decon),colnames(mcp_counter_decon))
  
  ## Make Score label rows to add to gene set table
  est_decon_gstab <- data.frame(GeneSet_Category = "Immune Deconvolution Gene Sets",
                                GeneSet_Sub_Category = "ESTIMATE Deconvolution Method",
                                GeneSet_Name = colnames(estimate_decon))
  mcp_decon_gstab <- data.frame(GeneSet_Category = "Immune Deconvolution Gene Sets",
                                GeneSet_Sub_Category = "MCP Counter Deconvolution Method",
                                GeneSet_Name = colnames(mcp_counter_decon))
  if (gsTab == TRUE) {
    GeneSetTable <- rbind(GeneSetTable,est_decon_gstab,mcp_decon_gstab)
  }
  if (gsTab == FALSE) {
    col_name <- colnames(GeneSetTable)[1]
    est_decon_gstab2 <- est_decon_gstab[,3, drop = F]
    colnames(est_decon_gstab2)[1] <- col_name
    mcp_decon_gstab2 <- mcp_decon_gstab[,3, drop = F]
    colnames(mcp_decon_gstab2)[1] <- col_name
    GeneSetTable <- rbind(GeneSetTable,est_decon_gstab2,mcp_decon_gstab2)
  }
  
  estimate_decon$SampleName <- rownames(estimate_decon)
  estimate_decon <- estimate_decon %>%
    relocate(SampleName)
  mcp_counter_decon$SampleName <- rownames(mcp_counter_decon)
  mcp_counter_decon <- mcp_counter_decon %>%
    relocate(SampleName)
  
  meta <- merge(meta,estimate_decon,by = "SampleName",all.x = T)
  meta <- merge(meta,mcp_counter_decon,by = "SampleName",all.x = T)
  
  decon_bin_cols <- c()
  for (i in decon_score_cols) {
    new_colname <- paste(i,"_HiLoScore",sep = "")
    meta[,paste(i,"_HiLoScore",sep = "")] <- highlow(meta[,which(colnames(meta) == i)])
    decon_bin_cols <- c(decon_bin_cols,new_colname)
  }
  
  metacol_feature <- c(metacol_feature,decon_bin_cols,decon_score_cols)
  
}


####----UI----####

ui <-
  navbarPage(paste("{ ",ProjectName," Survival Analysis }",sep = ""),
             
             ####----Overall Survival Tab----####
             
             tabPanel("Survival Analysis",
                      fluidPage(
                        sidebarLayout(
                          
                          ####----Sidebar Panel----####
                          
                          sidebarPanel(
                            tabsetPanel(
                              id = "survside",
                              
                              ##--Sample Parameters--##
                              
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
                                                  radioButtons("RawOrSS","Survival Analysis By:",
                                                               choices = c("Raw Gene Expression","Rank Normalized"),
                                                               selected = "Raw Gene Expression", inline = T),
                                                  uiOutput("rendGeneGeneSetTable"),
                                                  value = 2
                                         ),
                                         tabPanel("User Gene Set",
                                                  fileInput("userGeneSet","Gene Set Upload", accept = c(".gmt",".tsv",".txt",".RData")),
                                                  uiOutput("rendUserGeneSetTable"),
                                                  value = 3
                                         )
                                       ),
                                       uiOutput("rendViewGeneSetGenes"),
                                       uiOutput("rendGenesInGeneSetTab")
                              ),
                              
                              ##--Survival Parameters--##
                              
                              tabPanel("Survival Parameters",
                                       p(),
                                       h3("Survival Data Type"),
                                       column(6,
                                              uiOutput("rendSurvivalType_time")
                                       ),
                                       column(6,
                                              uiOutput("rendSurvivalType_id")
                                       ),
                                       hr(),
                                       h3("Risk Stratification Plot Parameters"),
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
                              
                              ##--Figure Parameters--##
                              
                              tabPanel("Figure Parameters",
                                       p(),
                                       h3("Survival Plot Parameters"),
                                       fluidRow(
                                         column(6,
                                                uiOutput("rendSurvXaxis")
                                         ),
                                         column(6,
                                                checkboxInput("ShowPval","Show P.Value",value = T)
                                         )
                                       ),
                                       hr(),
                                       h3("Boxplot Parameters"),
                                       fluidRow(
                                         column(6,
                                                numericInput("boxplotFont","Boxplot Font Size:", value = 15, step = 1),
                                                uiOutput("rendBoxoptselec")
                                         ),
                                         column(6,
                                                numericInput("boxplotDot", "Boxplot Dot Size:", value = 0.75, step = 0.25),
                                                selectInput("boxplotTextAngle","X-Axis Text Orientation",
                                                            choices = c("Horizontal (0 degrees)" = "0","Angled (45 degrees)" = "45","Vertical (90 degrees)" = "90","Stagger"))
                                         )
                                       ),
                                       hr(),
                                       h3("Heatmap Parameters"),
                                       selectInput("ClusterMethod", "Select Cluster Method",
                                                   choices = c("complete", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")),
                                       fluidRow(
                                         column(6,
                                                numericInput("heatmapFontR", "Heatmap Row Font Size:", value = 9, step = 1)
                                         ),
                                         column(6,
                                                numericInput("heatmapFontC", "Heatmap Column Font Size:", value = 10, step = 1)
                                         )
                                       ),
                                       selectInput("ColorPaletteHeat", "Select Color Palette:",
                                                   choices = c("Red/Blue" = "original",
                                                               "OmniBlueRed" = "OmniBlueRed",
                                                               "LightBlue/BlackRed" = "LightBlueBlackRed",
                                                               "Green/Black/Red" = "GreenBlackRed",
                                                               "Yellow/Green/Blue" = "YlGnBu","Inferno" = "Inferno",
                                                               "Viridis" = "Viridis","Plasma" = "Plasma",
                                                               "Reds" = "OrRd","Blues" = "PuBu","Greens" = "Greens")
                                       ),
                                       hr(),
                                       h3("Forest Plot Parameters"),
                                       numericInput("ForestFontSize","Font Size",value = 1),
                                       hr(),
                                       h3("Linearity Plot Parameters"),
                                       fluidRow(
                                         column(4,
                                                numericInput("linAxisFont","X/Y Axis Font Size",
                                                             value = 14, step = 1)
                                         ),
                                         column(4,
                                                numericInput("linTickFont","Axis Tick Font Size",
                                                             value = 10, step = 1)
                                         ),
                                         column(4,
                                                numericInput("linMainFont","Title Font Size",
                                                             value = 16, step = 1)
                                         )
                                       )
                              )
                            )
                          ),
                          
                          ####----Main Panel----####
                          
                          mainPanel(
                            tabsetPanel(
                              id = "surv",
                              
                              ####----Survival Analysis Tab----####
                              
                              tabPanel("Pathway Level Survival Analysis",
                                       p(),
                                       withSpinner(jqui_resizable(plotOutput("Splot", width = "100%", height = "500px")), type = 6),
                                       fluidRow(
                                         downloadButton("dnldSplot_SVG","Download as SVG"),
                                         downloadButton("dnldSplot_PDF","Download as PDF")
                                       ),
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
                                         downloadButton("dnldSplotBIN_SVG","Download as SVG"),
                                         downloadButton("dnldSplotBIN_PDF","Download as PDF")
                                       ),
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
                                       numericInput("QuantPercent","High/Low Risk Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                       fluidRow(
                                         downloadButton("dnldSquantPlot_SVG","Download as SVG"),
                                         downloadButton("dnldSquantPlot_PDF","Download as PDF")
                                       ),
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
                                       withSpinner(jqui_resizable(plotOutput("SquantPlot2", width = "100%", height = "500px")), type = 6),
                                       numericInput("QuantPercent2","Above/Below Risk Quantile Cutoff (%)", value = 25, min = 0, max = 100),
                                       fluidRow(
                                         downloadButton("dnldSquantPlot2_SVG","Download as SVG"),
                                         downloadButton("dnldSquantPlot2_PDF","Download as PDF")
                                       ),
                                       fluidRow(
                                         column(3,
                                                offset = 1,
                                                checkboxInput("ShowQuantHR2","Display Quantile Hazard Ratio Table")),
                                         column(4,
                                                uiOutput("rendQuantHRtab2"))
                                       ),
                                       fluidRow(
                                         tags$head(
                                           tags$style(HTML("hr {border-top: 1px solid #000000;}"))),
                                         hr()
                                       ),
                                       value = 1),
                              
                              ####----Univariate Survival----####
                              
                              tabPanel("Univariate Survival Analysis",
                                       p(),
                                       fluidRow(
                                         column(4,
                                                uiOutput("rendSurvivalFeatureSingle"),
                                                ## Allows all select inputs to be wide enough to read the contents
                                                tags$head(
                                                  tags$style(HTML('
                                                                  .selectize-input {
                                                                      white-space: nowrap;
                                                                  }
                                                                  .selectize-dropdown {
                                                                      width: 500px !important;
                                                                  }'
                                                  )
                                                  )
                                                ),
                                                fluidRow(
                                                  column(6,
                                                         checkboxInput("UniVarContCheck","Continuous Feature",value = F)
                                                  ),
                                                  column(6,
                                                         checkboxInput("UniVarNAcheck","Remove NA/Unknown",value = T)
                                                  )
                                                ),
                                                uiOutput("rendSurvFeatVariableUni")
                                         ),
                                         column(8,
                                                verbatimTextOutput("UnivarSummExpl")
                                         )
                                       ),
                                       tabsetPanel(
                                         id = "UniVarPlots",
                                         
                                         ##--Survival Plot--##
                                         
                                         tabPanel("Survival Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("featSplot", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldfeatSplot_SVG","Download as SVG"),
                                                    downloadButton("dnldfeatSplot_PDF","Download as PDF")
                                                  )
                                         ),
                                         
                                         ##--Coxh Tables--##
                                         
                                         tabPanel("Coxh Table",
                                                  p(),
                                                  fluidRow(
                                                    column(6,
                                                           div(withSpinner(tableOutput("SSingleFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("UnivarSummary")
                                                    )
                                                  )
                                         ),
                                         
                                         ##--Forest Plot--##
                                         
                                         tabPanel("Forest Plot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("SinglevarForestPlot", width = "100%", height = "800px")), type = 6)
                                         ),
                                         
                                         ##--Linearity Check--##
                                         
                                         tabPanel("Linearity Check",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           selectInput("ResidualTypeUni","Select Residual Type",
                                                                       choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                    ),
                                                    column(3,
                                                           selectInput("linPredict1", "X-axis Scale:",
                                                                       choices = c("linear.predictions","observation.id","time"))
                                                    )
                                                  ),
                                                  uiOutput("timewarnmessage1"),
                                                  withSpinner(jqui_resizable(plotOutput("UnivarLinearityPlot", width = "100%", height = "500px")), type = 6)
                                         )
                                       ),
                                       value = 6),
                              
                              ####----Multivariate Survival----####
                              
                              tabPanel("Multivariate Coxh Analysis",
                                       tabsetPanel(
                                         id = "multivariate",
                                         
                                         ####----Bivariate Additive Survival----####
                                         
                                         tabPanel("Bivariate Additive Survival Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           uiOutput("rendSurvivalFeatureBi1"),
                                                           fluidRow(
                                                             column(6,
                                                                    checkboxInput("BiVarAddContCheck1","Continuous Feature",value = F)
                                                             ),
                                                             column(6,
                                                                    checkboxInput("BiVarAddNAcheck1","Remove NA/Unknown",value = T)
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi1")
                                                    ),
                                                    column(3,
                                                           uiOutput("rendSurvivalFeatureBi2"),
                                                           fluidRow(
                                                             column(6,
                                                                    checkboxInput("BiVarAddContCheck2","Continuous Feature",value = F)
                                                             ),
                                                             column(6,
                                                                    checkboxInput("BiVarAddNAcheck2","Remove NA/Unknown",value = T)
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi2")
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("BivarAddSummExpl")
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "BiVarPlots",
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Cox HR Table",
                                                             p(),
                                                             fluidRow(
                                                               column(6,
                                                                      div(withSpinner(tableOutput("BiFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      verbatimTextOutput("bivarSummary"),
                                                                      fluidRow(
                                                                        column(6,
                                                                               verbatimTextOutput("bivarAnova1")
                                                                        ),
                                                                        column(6,
                                                                               verbatimTextOutput("bivarAnova2")
                                                                        )
                                                                      )
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("BivarForestPlot", width = "100%", height = "800px")), type = 6)
                                                    ),
                                                    
                                                    ##--Linearity Check--##
                                                    
                                                    tabPanel("Linearity Check",
                                                             p(),
                                                             fluidRow(
                                                               column(3,
                                                                      selectInput("ResidualTypeBi","Select Residual Type",
                                                                                  choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                               ),
                                                               column(3,
                                                                      selectInput("linPredict2", "X-axis Scale:",
                                                                                  choices = c("linear.predictions","observation.id","time"))
                                                               )
                                                             ),
                                                             uiOutput("timewarnmessage2"),
                                                             withSpinner(jqui_resizable(plotOutput("BivarLinearityPlot", width = "100%", height = "500px")), type = 6)
                                                    )
                                                  )
                                         ),
                                         
                                         ####----Bivariate Interaction Survival----####
                                         
                                         tabPanel("Bivariate Interaction Survival Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(3,
                                                           uiOutput("rendSurvivalFeatureBi1Inter"),
                                                           fluidRow(
                                                             column(6,
                                                                    checkboxInput("BiVarIntContCheck1","Continuous Feature",value = F)
                                                             ),
                                                             column(6,
                                                                    checkboxInput("BiVarIntNAcheck1","Remove NA/Unknown",value = T)
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi1Inter")
                                                    ),
                                                    column(3,
                                                           uiOutput("rendSurvivalFeatureBi2Inter"),
                                                           fluidRow(
                                                             column(6,
                                                                    checkboxInput("BiVarIntContCheck2","Continuous Feature",value = F)
                                                             ),
                                                             column(6,
                                                                    checkboxInput("BiVarIntNAcheck2","Remove NA/Unknown",value = T)
                                                             )
                                                           ),
                                                           uiOutput("rendSurvFeatVariableBi2Inter")
                                                    ),
                                                    column(6,
                                                           verbatimTextOutput("BivarIntSummExpl")
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "BiVarInterTabs",
                                                    
                                                    ##--Survival Plot--##
                                                    
                                                    tabPanel("Survival Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("featSplotBi", width = "100%", height = "500px")), type = 6),
                                                             fluidRow(
                                                               downloadButton("dnldfeatSplotBi_SVG","Download as SVG"),
                                                               downloadButton("dnldfeatSplotBi_PDF","Download as PDF")
                                                             )
                                                    ),
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Cox HR Table",
                                                             p(),
                                                             fluidRow(
                                                               column(6,
                                                                      div(withSpinner(tableOutput("BiFeatureHRtabInter"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      verbatimTextOutput("bivarSummaryInter"),
                                                                      verbatimTextOutput("bivarAnovaInter1")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("BivarForestPlotInter", width = "100%", height = "800px")), type = 6)
                                                    ),
                                                    
                                                    ##--Linearity Check--##
                                                    
                                                    tabPanel("Linearity Check",
                                                             p(),
                                                             fluidRow(
                                                               column(3,
                                                                      selectInput("ResidualTypeInter","Select Residual Type",
                                                                                  choices = c("deviance", "martingale", "score", "schoenfeld", "dfbeta", "dfbetas", "scaledsch", "partial"))
                                                               ),
                                                               column(3,
                                                                      selectInput("linPredict3", "X-axis Scale:",
                                                                                  choices = c("linear.predictions","observation.id","time"))
                                                               )
                                                             ),
                                                             uiOutput("timewarnmessage3"),
                                                             withSpinner(jqui_resizable(plotOutput("BivarLinearityPlotInter", width = "100%", height = "500px")), type = 6))
                                                    
                                                  )
                                         ),
                                         
                                         ####----Multivariate Survival----####
                                         
                                         tabPanel("Multivariate Coxh Analysis",
                                                  p(),
                                                  fluidRow(
                                                    column(4,
                                                           uiOutput("rendSurvivalFeature"),
                                                           checkboxInput("MultiVarNAcheck","Remove NA/Unknown",value = T)
                                                    )
                                                  ),
                                                  tabsetPanel(
                                                    id = "multivartabstwo",
                                                    
                                                    ##--Coxh Tables--##
                                                    
                                                    tabPanel("Coxh Tables",
                                                             fluidRow(
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Categorical)"),
                                                                      verbatimTextOutput("multivarSummary"),
                                                                      div(withSpinner(tableOutput("SFeatureHRtab"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               ),
                                                               column(6,
                                                                      h4("Coxh Hazard Ratio (Continuous)"),
                                                                      verbatimTextOutput("multivarSummaryCont"),
                                                                      div(withSpinner(tableOutput("SFeatureHRtabCont"), type = 7, size = 0.5), style = "font-size:12px; width:500px; overflow-X: scroll")
                                                               )
                                                             )
                                                    ),
                                                    
                                                    ##--Forest Plot--##
                                                    
                                                    tabPanel("Forest Plot",
                                                             p(),
                                                             withSpinner(jqui_resizable(plotOutput("MultivarForestPlot", width = "100%", height = "800px")), type = 6)
                                                    )
                                                  )
                                         )
                                       ),
                                       value = 7),
                              
                              ####----Data Exploration----####
                              
                              tabPanel("Data Exploration",
                                       tabsetPanel(
                                         id = "DataExploration",
                                         tabPanel("Download Survival Data",
                                                  p(),
                                                  uiOutput("rendMetaTableCols"),
                                                  uiOutput("rendMetaTable"),
                                                  fluidRow(
                                                    column(6,
                                                           uiOutput("DnldMetaButon")
                                                    ),
                                                    column(6,
                                                           uiOutput("DnldExprButon")
                                                    )
                                                  ),
                                                  value = 5),
                                         tabPanel("Score Density",
                                                  p(),
                                                  fluidRow(
                                                    numericInput("densityPercent","User Defined Percentile (Red)",value = 15, width = "200px"),
                                                    checkboxInput("QuartileLinesCheck","Show Quartile Lines (Blue)",value = T)
                                                  ),
                                                  withSpinner(jqui_resizable(plotOutput("ssgseaDensity", width = "100%", height = "500px")), type = 6),
                                                  fluidRow(
                                                    downloadButton("dnldssgseaDensity_SVG","Download as SVG"),
                                                    downloadButton("dnldssgseaDensity_PDF","Download as PDF")
                                                  ),
                                                  div(DT::dataTableOutput("ssgseaDensityTable"), style = "font-size:12px"),
                                                  downloadButton("dnldssgseaDensityTable","Download Table"),
                                                  value = 6),
                                         tabPanel("Risk Straification Boxplot",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("Sboxplot", width = "100%", height = "500px")), type = 6),
                                                  div(DT::dataTableOutput("SboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                  p(),
                                                  downloadButton("dnldSBoxplotTab","Download Table"),
                                                  value = 1),
                                         tabPanel("Risk Straification Heatmap",
                                                  p(),
                                                  withSpinner(jqui_resizable(plotOutput("Sheatmap", width = "100%", height = "2000px")), type = 6),
                                                  downloadButton("dnldSheatmapexpr","Download Expression Matrix From Heatmap"),
                                                  value = 2),
                                         tabPanel("Feature Boxplot",
                                                  p(),
                                                  uiOutput("rendBoxplotFeature"),
                                                  withSpinner(jqui_resizable(plotOutput("Featureboxplot", width = "100%", height = "500px")), type = 6),
                                                  div(DT::dataTableOutput("FeatureboxplotTable"), style = "font-size:12px; height:450px; overflow-Y: scroll"),
                                                  p(),
                                                  downloadButton("dnldFeatureboxplotTab","Download Table"),
                                                  value = 3),
                                         tabPanel("Feature Heatmap",
                                                  p(),
                                                  uiOutput("rendHeatmapFeature"),
                                                  withSpinner(jqui_resizable(plotOutput("FeatureHeatmap", width = "100%", height = "2000px")), type = 6),
                                                  downloadButton("dnldFheatmapexpr","Download Expression Matrix From Heatmap"),
                                                  value = 4)
                                       )
                              )
                            )
                          )
                        )
                      )
             )
  )



####----Server----####


server <- function(input, output, session) {
  
  
  ####----Render UI----####
  
  ## Select sample type to subset samples by - only render if more than one sample type
  output$rendSampleTypeSelection <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      SampleTypeChoices <- unique(meta[,metacol_sampletype])
      SampleTypeChoices <- c(SampleTypeChoices,"All_Sample_Types")
      selectInput("SampleTypeSelection",paste("Select Sample Type (",metacol_sampletype,"):",sep = ""),
                  choices = SampleTypeChoices, selected = PreSelect_SamplyType)
      
    }
    
  })
  
  ## Select primary feature to look at - All not working yet
  output$rendFeatureSelection <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection == "All_Sample_Types") {
        
        FeatureChoices <- c(metacol_sampletype,metacol_feature,"All_Features")
        selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices, selected = PreSelect_Feature)
        
      }
      else if (input$SampleTypeSelection != "All_Sample_Types") {
        
        FeatureChoices <- c(metacol_feature,"All_Features")
        selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices, selected = PreSelect_Feature)
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      FeatureChoices <- c(metacol_feature,"All_Features")
      selectInput("FeatureSelection","Select Feature:", choices = FeatureChoices, selected = PreSelect_Feature)
      
    }
    
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
      selectInput("subFeatureSelection","Feature Condition:", choices = SubFeatureChoices, selected = PreSelect_SubFeature)
      
    }
    
  })
  
  output$rendSurvivalFeatureSingle <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          metacol_feature <- metacol_feature[-which(metacol_feature == input$FeatureSelection)]
        }
        metacol_feature <- c(metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
        selectInput("SingleSurvivalFeature","Select Feature:",
                    choices = metacol_feature, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        
        SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature)
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        SurvFeatChoices2 <- c(SurvFeatChoices2,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
        selectInput("SingleSurvivalFeature","Select Feature:",
                    choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature)
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        SurvFeatChoices2 <- c(SurvFeatChoices2,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
      }
      selectInput("SingleSurvivalFeature","Select Feature:",
                  choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
    }
    
  })
  
  output$rendSurvivalFeatureBi1 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1","Select Feature 1:",
                    choices = SurvFeatChoices, selected = "HiLoPScore")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1","Select Feature 1:",
                    choices = SurvFeatChoices2, selected = "HiLoPScore")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi1","Select Feature 1:",
                  choices = SurvFeatChoices2, selected = "HiLoPScore")
    }
    
  })
  
  output$rendSurvivalFeatureBi2 <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2","Select Feature 2:",
                    choices = SurvFeatChoices, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2","Select Feature 2:",
                    choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi2","Select Feature 2:",
                  choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
    }
    
  })
  
  output$rendSurvivalFeatureBi1Inter <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                    choices = SurvFeatChoices, selected = "HiLoPScore")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                    choices = SurvFeatChoices2, selected = "HiLoPScore")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi1Inter","Select Feature 1:",
                  choices = SurvFeatChoices2, selected = "HiLoPScore")
    }
    
  })
  
  output$rendSurvivalFeatureBi2Inter <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                    choices = SurvFeatChoices, selected = PreSelect_SecondaryFeature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                    choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeatureBi2Inter","Select Feature 2:",
                  choices = SurvFeatChoices2, selected = PreSelect_SecondaryFeature)
    }
    
  })
  
  output$rendSurvivalFeature <- renderUI({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SurvFeatChoices <- c(metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    SurvFeatChoices2 <- c(metacol_sampletype,metacol_feature,geneset_name,"QuartilePScore","HiLoPScore","QuantCutoffPScore")
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices <- SurvFeatChoices[-which(SurvFeatChoices == input$FeatureSelection)]
        }
        selectInput("SurvivalFeature","Select Feature(s):",
                    choices = SurvFeatChoices, multiple = T, selected = "HiLoPScore")
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        if (input$FeatureSelection != "All_Features") {
          SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
        }
        selectInput("SurvivalFeature","Select Feature(s):",
                    choices = SurvFeatChoices2, multiple = T, selected = "HiLoPScore")
        
      }
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      if (input$FeatureSelection != "All_Features") {
        SurvFeatChoices2 <- SurvFeatChoices2[-which(SurvFeatChoices2 == input$FeatureSelection)]
      }
      selectInput("SurvivalFeature","Select Feature(s):",
                  choices = SurvFeatChoices2, multiple = T, selected = "HiLoPScore")
    }
    
  })
  
  output$rendSurvXaxis <- renderUI({
    
    meta_ssgsea <- ssGSEAmeta()
    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
      surv_time_col <- metacol_survtime[1]
      surv_id_col <- metacol_survid[1]
    }
    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
      surv_time_col <- input$SurvivalType_time
      surv_id_col <- input$SurvivalType_id
    }
    max_time <- ceiling(max(meta_ssgsea[,surv_time_col])/365.25)
    numericInput("SurvXaxis","X-Axis Limit (years)", value = max_time)
    
    
  })
  
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
  
  
  
  output$rendBoxplotFeature <- renderUI({
    
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("BoxplotFeature","Select Feature:",
                    choices = metacol_feature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("BoxplotFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature))
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("BoxplotFeature","Select Feature:",
                  choices = metacol_feature)
      
    }
    
  })
  
  output$rendHeatmapFeature <- renderUI({
    
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      
      if (input$SampleTypeSelection != "All_Sample_Types") {
        
        selectInput("HeatmapFeature","Select Feature:",
                    choices = metacol_feature)
        
      }
      else if (input$SampleTypeSelection == "All_Sample_Types") {
        
        selectInput("HeatmapFeature","Select Feature:",
                    choices = c(metacol_sampletype,metacol_feature))
        
      }
      
    }
    else if (length(unique(meta[,metacol_sampletype])) <= 1) {
      
      selectInput("HeatmapFeature","Select Feature:",
                  choices = metacol_feature)
      
    }
    
  })
  
  
  
  output$rendSurvFeatVariableUni <- renderUI({
    
    if (input$UniVarContCheck == FALSE) {
      
      Feature <- input$SingleSurvivalFeature
      metaSub <- ssGSEAmeta()
      Var_choices <- unique(metaSub[,Feature])
      if (input$UniVarNAcheck == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableUni","Select Coxh Feature Reference:",
                  choices = Var_choices)
      
    }
    
  })
  
  output$rendSurvFeatVariableBi1 <- renderUI({
    
    if (input$BiVarAddContCheck1 == FALSE) {
      
      Feature <- input$SurvivalFeatureBi1
      metaSub <- ssGSEAmeta()
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi1","Select Coxh Feature 1 Reference:",
                  choices = Var_choices)
      
    }
    
  })
  
  output$rendSurvFeatVariableBi2 <- renderUI({
    
    if (input$BiVarAddContCheck2 == FALSE) {
      
      Feature <- input$SurvivalFeatureBi2
      metaSub <- ssGSEAmeta()
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi2","Select Coxh Feature 2 Reference:",
                  choices = Var_choices)
      
    }
    
  })
  
  output$rendSurvFeatVariableBi1Inter <- renderUI({
    
    if (input$BiVarIntContCheck1 == FALSE) {
      
      Feature <- input$SurvivalFeatureBi1Inter
      metaSub <- ssGSEAmeta()
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi1Inter","Select Coxh Feature 1 Reference:",
                  choices = Var_choices)
      
    }
    
  })
  
  output$rendSurvFeatVariableBi2Inter <- renderUI({
    
    if (input$BiVarIntContCheck2 == FALSE) {
      
      Feature <- input$SurvivalFeatureBi2Inter
      metaSub <- ssGSEAmeta()
      Var_choices <- unique(metaSub[,Feature])
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        Var_choices <- Var_choices[which(is.na(Var_choices) == FALSE)]
        Var_choices <- Var_choices[grep("unknown",Var_choices,ignore.case = T, invert = T)]
        
      }
      Var_choices <- sort(Var_choices, decreasing = T, na.last = T)
      selectInput("SurvFeatVariableBi2Inter","Select Coxh Feature 2 Reference:",
                  choices = Var_choices)
      
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
    
    div(DT::dataTableOutput("MetaTable"), style = "font-size:12px")
    
    
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
  
  output$rendGenesInGeneSetTab <- renderUI({
    
    if (input$GeneSetTabs == 1) {
      
      if (input$ViewGeneSetGenes == TRUE) {
        
        div(DT::dataTableOutput("GenesInGeneSetTab"), style = "font-size:10px; height:450px; overflow-Y: scroll")
        
      }
      
    }
    
    else if (input$GeneSetTabs == 3) {
      
      req(input$userGeneSet)
      if (input$ViewGeneSetGenes == TRUE) {
        
        div(DT::dataTableOutput("GenesInGeneSetTab"), style = "font-size:10px; height:450px; overflow-Y: scroll")
        
      }
      
    }
    
  })
  
  output$rendViewGeneSetGenes <- renderUI({
    
    if (input$GeneSetTabs == 1) {
      
      checkboxInput("ViewGeneSetGenes","View Genes in Selected Gene Set", value = F)
      
    }
    
    else if (input$GeneSetTabs == 3) {
      
      req(input$userGeneSet)
      checkboxInput("ViewGeneSetGenes","View Genes in Selected Gene Set", value = F)
      
    }
    
  })
  
  ## Download button for subset meta
  output$DnldMetaButon <- renderUI({
    
    downloadButton("dnldMeta", "Download Meta Subset")
    
  })
  
  ## Download button for subset expression
  output$DnldExprButon <- renderUI({
    
    downloadButton("dnldExpr", "Download Expression Subset")
    
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
  
  
  ####----Data Tables----####
  
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
    if (Feature == "All_Features") {
      #userMetaCols <- userMetaCols[userMetaCols != Feature] #remove condition column from user selection because it is automatically added
      metaCols <- c(metaCols,surv_time_col,surv_id_col,userMetaCols) #combine column names selected
    }
    else if (Feature != "All_Features") {
      userMetaCols <- userMetaCols[userMetaCols != Feature] #remove condition column from user selection because it is automatically added
      metaCols <- c(metaCols,surv_time_col,surv_id_col,Feature,userMetaCols) #combine column names selected
    }
    meta_sub <- meta[,metaCols]
    DT::datatable(meta_sub,
                  extensions = "FixedColumns",
                  options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                 pageLength = 20,
                                 fixedColumns = list(leftColumns = 1),
                                 scrollX = T),
                  rownames = F)
    
  })
  
  output$GenesInGeneSetTab <- renderDataTable({
    
    if (input$GeneSetTabs == 1) {
      
      if (input$ViewGeneSetGenes == TRUE) {
        
        geneset <- gs_react()
        gs_df <- as.data.frame(geneset)
        DT::datatable(gs_df, options = list(paging = F), rownames = F)
        
      }
      
    }
    
    else if (input$GeneSetTabs == 3) {
      
      if (input$ViewGeneSetGenes == TRUE) {
        
        req(input$userGeneSet)
        geneset <- gs_react()
        gs_df <- as.data.frame(geneset)
        DT::datatable(gs_df, options = list(paging = F), rownames = F)
        
      }
      
    }
    
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
  
  output$ssgseaDensityTable <- renderDataTable({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    table <- ssgsea_meta[,c("SampleName",GeneSet)]
    DT::datatable(table,
                  options = list(scrollY = T),
                  rownames = F)
    
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
  
  ####----Reactives----####
  
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
        geneset_name <- GeneSetTable[input$GeneSetTable_rows_selected,3]
        if (geneset_name %in% decon_score_cols) {
          geneset <- list(meta[,geneset_name])
          names(geneset) <- geneset_name
        }
        else {
          geneset <- gs[geneset_name]
        }
      }
      if (gsTab == FALSE) {
        geneset_name <- GeneSetTable[input$GeneSetTable_rows_selected,1]
        if (geneset_name %in% decon_score_cols) {
          geneset <- list(meta[,geneset_name])
          names(geneset) <- geneset_name
        }
        else {
          geneset <- gs[geneset_name]
        }
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
    quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
    quantCutoff2 <- input$QuantPercent2/100 #Quantile cutoff given by user
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
    expr_sub <- expr[,colnames(expr) %in% samples]
    expr_mat <- as.matrix(expr_sub)
    rownames(expr_mat) <- rownames(expr_sub)
    colnames(expr_mat) <- colnames(expr_sub)
    
    if (geneset_name %in% decon_score_cols) {
      
      ssGSEA <- meta[,c("SampleName",geneset_name)]
      
    }
    else if ((geneset_name %in% decon_score_cols) == FALSE) {
      
      if (input$GeneSetTabs == 1 | input$GeneSetTabs == 3) {
        
        ## Perform ssGSEA with gs and new subset data
        scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
        ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod)
        ## Transform
        ssGSEA <- as.data.frame(t(ssGSEA))
        ssGSEA$SampleName <- rownames(ssGSEA)
        
      }
      else if (input$GeneSetTabs == 2) {
        
        if (input$RawOrSS == "Raw Gene Expression") {
          
          expr_sub2 <- expr_sub[geneset_name,]
          expr_sub3 <- as.data.frame(t(expr_sub2))
          expr_sub3$SampleName <- rownames(expr_sub3)
          ssGSEA <- expr_sub3
          
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          
          ## Perform ssGSEA with gs and new subset data
          scoreMethod <- input$ScoreMethod      #ssGSEA score method chosen
          ssGSEA <- gsva(expr_mat,geneset,method = scoreMethod)
          ## Transform
          ssGSEA <- as.data.frame(t(ssGSEA))
          ssGSEA$SampleName <- rownames(ssGSEA)
          
        }
        
      }
      
    }
    
    ## Perform further functions
    ssGSEA$VAR_Q <- quartile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
    ssGSEA$QuartilePScore <- paste("", ssGSEA$VAR_Q, sep="")
    ssGSEA$HiLoPScore <- highlow(ssGSEA[, which(colnames(ssGSEA) == geneset_name)])
    ssGSEA$QuantCutoffPScore <- quantile_conversion(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff)
    ssGSEA$AboveBelowCutoffPScore <- quantile_conversion2(ssGSEA[, which(colnames(ssGSEA) == geneset_name)], quantCutoff2)
    
    ## Merge with meta
    if (geneset_name %in% decon_score_cols) {
      ssGSEA <- ssGSEA[,-2] #remove score column so on merge the column does not duplicate
    }
    meta_ssGSEA <- merge(meta,ssGSEA, by = "SampleName", all = T)
    meta_ssGSEA
    
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
  
  
  ####----Survival Plots----####
  
  ## Survival Plot - Quartile
  Splot_react <- reactive({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    xaxlim <- input$SurvXaxis * 365.25
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"QuartilePScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ QuartilePScore, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    if (input$GeneSetTabs == 2) {
      if (input$RawOrSS == "Raw Gene Expression") {
        scoreMethodLab <- "Raw Gene Expression"
      }
      else if (input$RawOrSS == "Rank Normalized") {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs != 2) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    if (is.null(input$SurvXaxis) == TRUE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients\n", geneset_name," (",scoreMethodLab," in quartiles)", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain")
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    else if (is.null(input$SurvXaxis) == FALSE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients\n", geneset_name," (",scoreMethodLab," score in quartiles)", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain"),
                           xlim = c(0,xaxlim)
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
  })
  
  output$Splot <- renderPlot({
    plot <- Splot_react()
    plot
  })
  
  ## Survival Plot - Binary
  SplotBIN_react <- reactive({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    xaxlim <- input$SurvXaxis * 365.25
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"HiLoPScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ HiLoPScore, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    if (input$GeneSetTabs == 2) {
      if (input$RawOrSS == "Raw Gene Expression") {
        scoreMethodLab <- "Raw Gene Expression"
      }
      else if (input$RawOrSS == "Rank Normalized") {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs != 2) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    if (is.null(input$SurvXaxis) == TRUE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients\n",
                                         geneset_name," (",scoreMethodLab, " Median Cutoff)", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain")
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    else if (is.null(input$SurvXaxis) == FALSE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients\n",
                                         geneset_name," (",scoreMethodLab, " Median Cutoff)", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain"),
                           xlim = c(0,xaxlim)
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    
    
  })
  
  output$SplotBIN <- renderPlot({
    
    plot <- SplotBIN_react()
    plot
    
  })
  
  ## Survival Plot - Quantile
  SquantPlot_react <- reactive({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    xaxlim <- input$SurvXaxis * 365.25
    quantCutoffOG <- input$QuantPercent
    show_pval <- input$ShowPval
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"QuantCutoffPScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Remove between cutoff samples
    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$QuantCutoffPScore != "BetweenCutoff"),]
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ QuantCutoffPScore, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    if (input$GeneSetTabs == 2) {
      if (input$RawOrSS == "Raw Gene Expression") {
        scoreMethodLab <- "Raw Gene Expression"
      }
      else if (input$RawOrSS == "Rank Normalized") {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs != 2) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    if (is.null(input$SurvXaxis) == TRUE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ", Feature,SampleTypeLab,"\n",
                                         geneset_name," (Top (Bottom) ",labelQuantCutoff," Patients Split based on ",scoreMethodLab," )", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain")
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    else if (is.null(input$SurvXaxis) == FALSE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ", Feature,SampleTypeLab,"\n",
                                         geneset_name," (Top (Bottom) ",labelQuantCutoff," Patients Split based on ",scoreMethodLab," score)", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain"),
                           xlim = c(0,xaxlim)
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    
    
  })
  
  output$SquantPlot <- renderPlot({
    plot <- SquantPlot_react()
    plot
  })
  
  ## Survival Plot - Quantile
  SquantPlot2_react <- reactive({
    
    ## Assign variables
    geneset <- gs_react()
    geneset_name <- names(geneset)
    SampleType <- input$SampleTypeSelection
    Feature <- input$FeatureSelection
    scoreMethod <- input$ScoreMethod
    show_pval <- input$ShowPval
    xaxlim <- input$SurvXaxis * 365.25
    quantCutoffOG <- input$QuantPercent2
    quantCutoff <- input$QuantPercent2/100
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"AboveBelowCutoffPScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Remove between cutoff samples
    #meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$QuantCutoffPScore != "BetweenCutoff"),]
    
    ## Survival Function
    fit <- survfit(Surv(OS.time,OS.ID) ~ AboveBelowCutoffPScore, data = meta_ssgsea_sdf, type="kaplan-meier")
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    if (input$GeneSetTabs == 2) {
      if (input$RawOrSS == "Raw Gene Expression") {
        scoreMethodLab <- "Raw Gene Expression"
      }
      else if (input$RawOrSS == "Rank Normalized") {
        scoreMethodLab <- paste(scoreMethod, " score", sep = "")
      }
    }
    else if (input$GeneSetTabs != 2) {
      scoreMethodLab <- paste(scoreMethod, " score", sep = "")
    }
    
    if (is.null(input$SurvXaxis) == TRUE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ", Feature,SampleTypeLab,"\n",
                                         geneset_name," (Above (Below) ",labelQuantCutoff," Patients Split based on ",scoreMethodLab," )", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain")
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    else if (is.null(input$SurvXaxis) == FALSE) {
      
      ## Generate plot
      ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                           title = paste("Survival curves of ", Feature,SampleTypeLab,"\n",
                                         geneset_name," (Above (Below) ",labelQuantCutoff," Patients Split based on ",scoreMethodLab," score)", sep = ""),
                           xscale = c("d_y"),
                           break.time.by=365.25,
                           xlab = "Years", 
                           ylab = paste(SurvDateType,"Survival Probability"),
                           submain = "Based on Kaplan-Meier estimates",
                           caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                           font.title = c(16, "bold"),
                           font.submain = c(12, "italic"),
                           font.caption = c(12, "plain"),
                           font.x = c(14, "plain"),
                           font.y = c(14, "plain"),
                           font.tickslab = c(12, "plain"),
                           xlim = c(0,xaxlim)
      )
      
      ggsurv$table <- ggsurv$table + theme_cleantable()
      ggsurv
      
    }
    
    
    
  })
  
  output$SquantPlot2 <- renderPlot({
    plot <- SquantPlot2_react()
    plot
  })
  
  ## Survival Plot - SINGLE FEATURE
  featSplot_react <- reactive({
    
    if (length(input$SingleSurvivalFeature > 0)) {
      
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SingleSurvivalFeature
      scoreMethod <- input$ScoreMethod
      show_pval <- input$ShowPval
      xaxlim <- input$SurvXaxis * 365.25
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$UniVarNAcheck == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      #select_cols <- c("SampleName",surv_time_col,surv_id_col,"Prim.Tumor.Grade")
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$UniVarNAcheck == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature <- gsub("[[:punct:]]","_",Feature)
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
      
      form <- paste("Surv(OS.time,OS.ID) ~ ",Feature,sep = "")
      form2 <- as.formula(form)
      fit <- eval(substitute(survfit(form2,data = meta_ssgsea_sdf, type="kaplan-meier")))
      
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- " "
      }
      
      if (is.null(input$SurvXaxis) == TRUE) {
        
        ## Generate plot
        ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                             title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients", sep = ""),
                             xscale = c("d_y"),
                             break.time.by=365.25,
                             xlab = "Years", 
                             ylab = paste(SurvDateType,"Survival Probability"),
                             submain = "Based on Kaplan-Meier estimates",
                             caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                             font.title = c(16, "bold"),
                             font.submain = c(12, "italic"),
                             font.caption = c(12, "plain"),
                             font.x = c(14, "plain"),
                             font.y = c(14, "plain"),
                             font.tickslab = c(12, "plain")
        )
        
        ggsurv$table <- ggsurv$table + theme_cleantable()
        ggsurv
        
      }
      
      else if (is.null(input$SurvXaxis) == FALSE) {
        
        ## Generate plot
        ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                             title = paste("Survival curves of ",Feature,SampleTypeLab,"Patients", sep = ""),
                             xscale = c("d_y"),
                             break.time.by=365.25,
                             xlab = "Years", 
                             ylab = paste(SurvDateType,"Survival Probability"),
                             submain = "Based on Kaplan-Meier estimates",
                             caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                             font.title = c(16, "bold"),
                             font.submain = c(12, "italic"),
                             font.caption = c(12, "plain"),
                             font.x = c(14, "plain"),
                             font.y = c(14, "plain"),
                             font.tickslab = c(12, "plain"),
                             xlim = c(0,xaxlim)
        )
        
        ggsurv$table <- ggsurv$table + theme_cleantable()
        ggsurv
        
      }
      
      
      
    }
    
  })
  
  output$featSplot <- renderPlot({
    plot <- featSplot_react()
    plot
  })
  
  ## Survival Plot - TWO FEATURE
  featSplotBi_react <- reactive({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter> 0)) {
      
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      scoreMethod <- input$ScoreMethod
      show_pval <- input$ShowPval
      xaxlim <- input$SurvXaxis * 365.25
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarIntNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarIntNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
      colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
      
      form <- paste("Surv(OS.time,OS.ID) ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")
      form2 <- as.formula(form)
      fit <- eval(substitute(survfit(form2,data = meta_ssgsea_sdf, type="kaplan-meier")))
      
      ## Adjust 'Sample Type' for label 
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
      }
      if (length(unique(meta[,metacol_sampletype])) <= 1) {
        SampleTypeLab <- " "
      }
      
      if (is.null(input$SurvXaxis) == TRUE) {
        
        ## Generate plot
        ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                             title = paste("Survival curves of ",Feature1," and\n",Feature2," in",SampleTypeLab,"Patients", sep = ""),
                             xscale = c("d_y"),
                             break.time.by=365.25,
                             xlab = "Years", 
                             ylab = paste(SurvDateType,"Survival Probability"),
                             submain = "Based on Kaplan-Meier estimates",
                             caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                             font.title = c(16, "bold"),
                             font.submain = c(12, "italic"),
                             font.caption = c(12, "plain"),
                             font.x = c(14, "plain"),
                             font.y = c(14, "plain"),
                             font.tickslab = c(12, "plain")
        )
        
        ggsurv$table <- ggsurv$table + theme_cleantable()
        ggsurv
        
      }
      
      else if (is.null(input$SurvXaxis) == FALSE) {
        
        ## Generate plot
        ggsurv <- ggsurvplot(fit, data = meta_ssgsea_sdf, risk.table = TRUE,
                             title = paste("Survival curves of ",Feature1," and ",Feature2," in",SampleTypeLab,"Patients", sep = ""),
                             xscale = c("d_y"),
                             break.time.by=365.25,
                             xlab = "Years", 
                             ylab = paste(SurvDateType,"Survival Probability"),
                             submain = "Based on Kaplan-Meier estimates",
                             caption = "created with survminer", pval=show_pval, ggtheme = theme_bw(),
                             font.title = c(16, "bold"),
                             font.submain = c(12, "italic"),
                             font.caption = c(12, "plain"),
                             font.x = c(14, "plain"),
                             font.y = c(14, "plain"),
                             font.tickslab = c(12, "plain"),
                             xlim = c(0,xaxlim)
        )
        
        ggsurv$table <- ggsurv$table + theme_cleantable()
        ggsurv
        
      }
      
      
      
    }
    
  })
  
  output$featSplotBi <- renderPlot({
    plot <- featSplotBi_react()
    plot
  })
  
  ####----HR Tables----####
  
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"QuartilePScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    tab <- coxph(Surv(OS.time,OS.ID) ~ QuartilePScore, data = meta_ssgsea_sdf) %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"HiLoPScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Survival Function
    tab <- coxph(Surv(OS.time,OS.ID) ~ HiLoPScore, data = meta_ssgsea_sdf) %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
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
    meta_ssgsea_sdf <- meta_ssgsea[,c("SampleName",surv_time_col,surv_id_col,"QuantCutoffPScore")]
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_time_col)] <- "OS.time"
    colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == surv_id_col)] <- "OS.ID"
    
    ## Remove between cutoff samples
    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(meta_ssgsea_sdf$QuantCutoffPScore != "BetweenCutoff"),]
    
    ## Survival Function
    tab <- coxph(Surv(OS.time,OS.ID) ~ QuantCutoffPScore, data = meta_ssgsea_sdf) %>% 
      gtsummary::tbl_regression(exp = TRUE) %>%
      as_gt()
    
    tab_df <- as.data.frame(tab)
    
    tab_df <- tab_df %>%
      select(label,estimate,ci,p.value)
    colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
    
    tab_df
    
  })
  
  output$SSingleFeatureHRtab <- renderTable({
    
    if (input$UniVarContCheck == TRUE) {
      
      if (length(input$SingleSurvivalFeature > 0)) {
        ## Assign variables
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SampleType <- input$SampleTypeSelection
        Feature <- input$SingleSurvivalFeature
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
        
        quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
        if (input$UniVarNAcheck == TRUE) {
          
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
          ## Re-Perform Stat functions
          meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
          meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
          meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
          meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
          
        }
        
        ## Subset columns needed for plot and rename for surv function
        select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
        meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
        
        #if (input$UniVarNAcheck == TRUE) {
        #  
        #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature]) == FALSE),]
        #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature],ignore.case = T, invert = T),]
        #  
        #}
        
        colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
        Feature <- gsub("[[:punct:]]","_",Feature)
        
        
        ## Survival Function
        tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature,sep = "")),
                     data = meta_ssgsea_sdf) %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        
        tab_df <- as.data.frame(tab)
        tab_df[is.na(tab_df)] <- ""
        tab_df <- tab_df %>%
          select(label,n_obs,estimate,std.error,ci,p.value)
        colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
        
        tab_df
      }
      
    }
    
    else if (input$UniVarContCheck == FALSE) {
      
      if (length(input$SingleSurvivalFeature > 0)) {
        ## Assign variables
        geneset <- gs_react()
        geneset_name <- names(geneset)
        SampleType <- input$SampleTypeSelection
        Feature <- input$SingleSurvivalFeature
        ref_Feature <- input$SurvFeatVariableUni
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
        
        quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
        if (input$UniVarNAcheck == TRUE) {
          
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
          ## Re-Perform Stat functions
          meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
          meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
          meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
          meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
          
        }
        
        ## Subset columns needed for plot and rename for surv function
        select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
        meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
        
        #if (input$UniVarNAcheck == TRUE) {
        #  
        #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature]) == FALSE),]
        #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature],ignore.case = T, invert = T),]
        #  
        #}
        
        meta_ssgsea_sdf[,Feature] <- as.factor(meta_ssgsea_sdf[,Feature])
        meta_ssgsea_sdf[,Feature] <- relevel(meta_ssgsea_sdf[,Feature], ref = ref_Feature)
        colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
        Feature <- gsub("[[:punct:]]","_",Feature)
        
        
        ## Survival Function
        tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature,sep = "")),
                     data = meta_ssgsea_sdf) %>% 
          gtsummary::tbl_regression(exp = TRUE) %>%
          as_gt()
        
        tab_df <- as.data.frame(tab)
        tab_df[is.na(tab_df)] <- ""
        tab_df <- tab_df %>%
          select(label,n_obs,estimate,std.error,ci,p.value)
        colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
        
        tab_df
      }
      
    }
    
  })
  
  output$BiFeatureHRtab <- renderTable({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
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
      
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarAddNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarAddNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf) %>% 
        gtsummary::tbl_regression(exp = TRUE) %>%
        as_gt()
      
      tab_df <- as.data.frame(tab)
      tab_df[is.na(tab_df)] <- ""
      tab_df <- tab_df %>%
        select(label,n_obs,estimate,std.error,ci,p.value)
      colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
      
      tab_df
    }
    
    
  })
  
  output$BiFeatureHRtabInter <- renderTable({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      Feat1Var <- input$SurvFeatVariableBi1Inter
      Feat2Var <- input$SurvFeatVariableBi2Inter
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarIntNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarIntNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarIntContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarIntContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarIntContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarIntContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf) %>% 
        gtsummary::tbl_regression(exp = TRUE) %>%
        as_gt()
      
      tab_df <- as.data.frame(tab)
      tab_df[is.na(tab_df)] <- ""
      tab_df <- tab_df %>%
        select(label,n_obs,estimate,std.error,ci,p.value)
      colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
      
      tab_df
    }
    
    
  })
  
  output$SFeatureHRtab <- renderTable({
    
    if (length(input$SurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SurvivalFeature
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$MultiVarNAcheck == TRUE) {
        
        for (i in Feature){
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
        }
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      
      #if (input$MultiVarNAcheck == TRUE) {
      #  
      #  for (i in Feature){
      #    
      #    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,i]) == FALSE),]
      #    meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,i],ignore.case = T, invert = T),]
      #  }
      #  
      #}
      
      for (i in Feature){
        meta_ssgsea_sdf[,i] <- as.factor(meta_ssgsea_sdf[,i])
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
      }
      
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature,collapse = "+"),sep = "")),
                   data = meta_ssgsea_sdf) %>% 
        gtsummary::tbl_regression(exp = TRUE) %>%
        as_gt()
      
      #tab_df <- as.data.frame(tab)
      #
      #tab_df <- tab_df %>%
      #  select(label,estimate,ci,p.value)
      #colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
      
      tab_df <- as.data.frame(tab)
      tab_df[is.na(tab_df)] <- ""
      tab_df <- tab_df %>%
        select(label,n_obs,estimate,std.error,ci,p.value)
      colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
      
      tab_df
    }
    
    
  })
  
  output$SFeatureHRtabCont <- renderTable({
    
    if (length(input$SurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SurvivalFeature
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$MultiVarNAcheck == TRUE) {
        
        for (i in Feature){
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
        }
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$MultiVarNAcheck == TRUE) {
      #  
      #  for (i in Feature){
      #    
      #    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,i]) == FALSE),]
      #    meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,i],ignore.case = T, invert = T),]
      #  }
      #  
      #}
      
      for (i in Feature) {
        
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
        
      }
      
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature, collapse = "+"),sep = "")),
                   data = meta_ssgsea_sdf) %>% 
        gtsummary::tbl_regression(exp = TRUE) %>%
        as_gt()
      
      #tab_df <- as.data.frame(tab)
      #
      #tab_df <- tab_df %>%
      #  select(label,estimate,ci,p.value)
      #colnames(tab_df) <- c("Characteristic","Hazard Ratio","95% Confidence Interval","P.Value")
      
      tab_df <- as.data.frame(tab)
      tab_df[is.na(tab_df)] <- ""
      tab_df <- tab_df %>%
        select(label,n_obs,estimate,std.error,ci,p.value)
      colnames(tab_df) <- c("Variable","N","Hazard Ratio","Std. Error","95% Confidence Interval","P.Value")
      
      tab_df
    }
    
    
  })
  
  output$multivarSummary <- renderPrint({
    
    if (length(input$SurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SurvivalFeature
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$MultiVarNAcheck == TRUE) {
        
        for (i in Feature){
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
        }
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      for (i in Feature){
        meta_ssgsea_sdf[,i] <- as.factor(meta_ssgsea_sdf[,i])
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
      }
      
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature,collapse = "+"),sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary (Categorical):",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
    }
    
  })
  
  output$multivarSummaryCont <- renderPrint({
    
    if (length(input$SurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SurvivalFeature
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$MultiVarNAcheck == TRUE) {
        
        for (i in Feature){
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
        }
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      for (i in Feature) {
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
      }
      
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature, collapse = "+"),sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary (Continuous):",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
    }
    
  })
  
  output$bivarSummary <- renderPrint({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarAddNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarAddNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
    }
    
  })
  
  output$bivarSummaryCont <- renderPrint({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
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
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
    }
    
  })
  
  output$bivarAnova1 <- renderPrint({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarAddNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarAddNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      meta_ssgsea_sdf <- meta_ssgsea_sdf[complete.cases(meta_ssgsea_sdf),]
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #only rows with both features
      meta_ssgsea_sdf <- meta_ssgsea_sdf[complete.cases(meta_ssgsea_sdf),]
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab1 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                    data = meta_ssgsea_sdf)
      tab2 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature1,sep = "")),
                    data = meta_ssgsea_sdf)
      
      annova_res <- anova(tab1,tab2)
      
      out <- capture.output(annova_res)
      
      line1 <- out[3]
      line2 <- out[4]
      line3 <- out[5]
      line4 <- out[6]
      line5 <- out[7]
      
      text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
      cat(text)
    }
    
  })
  
  output$bivarAnova2 <- renderPrint({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarAddNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarAddNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      meta_ssgsea_sdf <- meta_ssgsea_sdf[complete.cases(meta_ssgsea_sdf),]
      
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #only rows with both features
      meta_ssgsea_sdf <- meta_ssgsea_sdf[complete.cases(meta_ssgsea_sdf),]
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab1 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                    data = meta_ssgsea_sdf)
      tab2 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature2,sep = "")),
                    data = meta_ssgsea_sdf)
      
      annova_res <- anova(tab1,tab2)
      
      out <- capture.output(annova_res)
      
      line1 <- out[3]
      line2 <- out[4]
      line3 <- out[5]
      line4 <- out[6]
      line5 <- out[7]
      
      text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
      cat(text)
    }
    
  })
  
  output$bivarAnovaInter1 <- renderPrint({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      Feat1Var <- input$SurvFeatVariableBi1Inter
      Feat2Var <- input$SurvFeatVariableBi2Inter
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarIntNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarIntNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      meta_ssgsea_sdf <- meta_ssgsea_sdf[complete.cases(meta_ssgsea_sdf),]
      
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarIntContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarIntContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarIntContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarIntContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab1 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")),
                    data = meta_ssgsea_sdf)
      tab2 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                    data = meta_ssgsea_sdf)
      
      annova_res <- anova(tab1,tab2)
      
      out <- capture.output(annova_res)
      
      line1 <- out[3]
      line2 <- out[4]
      line3 <- out[5]
      line4 <- out[6]
      line5 <- out[7]
      
      text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
      cat(text)
    }
    
  })
  
  #output$bivarAnovaInter2 <- renderPrint({
  #  
  #  if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
  #    ## Assign variables
  #    geneset <- gs_react()
  #    geneset_name <- names(geneset)
  #    SampleType <- input$SampleTypeSelection
  #    Feature1 <- input$SurvivalFeatureBi1Inter
  #    Feature2 <- input$SurvivalFeatureBi2Inter
  #    Feat1Var <- input$SurvFeatVariableBi1Inter
  #    Feat2Var <- input$SurvFeatVariableBi2Inter
  #    scoreMethod <- input$ScoreMethod
  #    if (is.null(input$SurvivalType_time) == TRUE & is.null(input$SurvivalType_id) == TRUE) {
  #      surv_time_col <- metacol_survtime[1]
  #      surv_id_col <- metacol_survid[1]
  #    }
  #    if (is.null(input$SurvivalType_time) == FALSE & is.null(input$SurvivalType_id) == FALSE) {
  #      surv_time_col <- input$SurvivalType_time
  #      surv_id_col <- input$SurvivalType_id
  #    }
  #    expr <- exprSub()
  #    meta_ssgsea <- ssGSEAmeta()
  #    
  #    ## Determine type of survival data - OS/EFS/PFS?
  #    SurvDateType <- sub("\\..*","",surv_time_col)
  #    
  #    ## Remove rows with NA in survival column
  #    meta_ssgsea <- meta_ssgsea[!is.na(meta_ssgsea[,surv_time_col]),]
  #    
  #    ## Subset columns needed for plot and rename for surv function
  #    select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
  #    meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
  #    
  #    meta_ssgsea_sdf <- meta_ssgsea_sdf[complete.cases(meta_ssgsea_sdf),]
  #    
  #    colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
  #    meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
  #    Feature1 <- gsub("[[:punct:]]","_",Feature1)
  #    Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
  #    colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
  #    meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
  #    Feature2 <- gsub("[[:punct:]]","_",Feature2)
  #    Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
  #    
  #    #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
  #    #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
  #    #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
  #    #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
  #    
  #    if (input$BiVarIntContCheck1 == FALSE) {
  #      meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
  #      meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
  #    }
  #    else if (input$BiVarIntContCheck1 == TRUE) {
  #      meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
  #    }
  #    if (input$BiVarIntContCheck2 == FALSE) {
  #      meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
  #      meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
  #    }
  #    else if (input$BiVarIntContCheck2 == TRUE) {
  #      meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
  #    }
  #    
  #    ## Survival Function
  #    tab1 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")),
  #                  data = meta_ssgsea_sdf)
  #    tab2 <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature2,"+",Feature1,sep = ""),sep = "")),
  #                  data = meta_ssgsea_sdf)
  #    
  #    annova_res <- anova(tab1,tab2)
  #    
  #    out <- capture.output(annova_res)
  #    
  #    line1 <- out[3]
  #    line2 <- out[4]
  #    line3 <- out[5]
  #    line4 <- out[6]
  #    line5 <- out[7]
  #    
  #    text <- paste("Model Comparison:",line1,line2,line3,line4,line5,sep = "\n")
  #    cat(text)
  #  }
  #  
  #})
  
  output$bivarSummaryInter <- renderPrint({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      Feat1Var <- input$SurvFeatVariableBi1Inter
      Feat2Var <- input$SurvFeatVariableBi2Inter
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarIntNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarIntNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarIntContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarIntContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarIntContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarIntContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary:",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
    }
    
  })
  
  output$UnivarSummary <- renderPrint({
    
    if (input$UniVarContCheck == TRUE) {
      
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SingleSurvivalFeature
      ref_Feature <- input$SurvFeatVariableUni
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$UniVarNAcheck == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$UniVarNAcheck == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature,sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary (Continuous):",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
      
    }
    
    if (input$UniVarContCheck == FALSE) {
      
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SingleSurvivalFeature
      ref_Feature <- input$SurvFeatVariableUni
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
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      meta_ssgsea_sdf[,Feature] <- as.factor(meta_ssgsea_sdf[,Feature])
      meta_ssgsea_sdf[,Feature] <- relevel(meta_ssgsea_sdf[,Feature], ref = ref_Feature)
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature,sep = "")),
                   data = meta_ssgsea_sdf)
      
      out <- capture.output(summary(tab))
      
      con_line <- grep("^Concordance=",out,value = T)
      lik_line <- grep("^Likelihood ratio test=",out,value = T)
      wal_line <- grep("^Wald test",out,value = T)
      sco_line <- grep("^Score ",out,value = T)
      
      text <- paste("Coxh Summary (Categorical):",con_line,lik_line,wal_line,sco_line,sep = "\n")
      cat(text)
      
    }
    
  })
  
  
  ####----Other Plots----####
  
  output$BivarLinearityPlotInter <- renderPlot({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      Feat1Var <- input$SurvFeatVariableBi1Inter
      Feat2Var <- input$SurvFeatVariableBi2Inter
      residType <- input$ResidualTypeInter
      AxisFont <- input$linAxisFont
      MainFont <- input$linMainFont
      TickFont <- input$linTickFont
      linpredict <- input$linPredict3
      #if (linpredict_choice == TRUE) {
      #  linpredict <- "linear.predictions"
      #}
      #else if (linpredict_choice == FALSE) {
      #  linpredict <- "observation.id"
      #}
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarIntNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarIntNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarIntContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarIntContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarIntContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarIntContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      p <- ggcoxdiagnostics(tab,
                            type = residType,
                            sline = T,
                            sline.se = T,
                            ggtheme = theme_minimal(),
                            ox.scale = linpredict)
      p <- ggpar(p,
                 font.x = AxisFont,
                 font.y = AxisFont,
                 font.main = MainFont,
                 font.tickslab = TickFont,
                 main = paste("Linearity Plot Featuring: ",Feature1," * ",Feature2, sep = ""),
                 ylab = paste(str_to_title(residType)," Residuals", sep = "")
      )
      p
      
    }
    
  })
  
  output$BivarLinearityPlot <- renderPlot({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
      residType <- input$ResidualTypeBi
      AxisFont <- input$linAxisFont
      MainFont <- input$linMainFont
      TickFont <- input$linTickFont
      linpredict <- input$linPredict2
      #if (linpredict_choice == TRUE) {
      #  linpredict <- "linear.predictions"
      #}
      #else if (linpredict_choice == FALSE) {
      #  linpredict <- "observation.id"
      #}
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
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
        meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
        meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
        
      }
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      p <- ggcoxdiagnostics(tab,
                            type = residType,
                            sline = T,
                            sline.se = T,
                            ggtheme = theme_minimal(),
                            ox.scale = linpredict)
      p <- ggpar(p,
                 font.x = AxisFont,
                 font.y = AxisFont,
                 font.main = MainFont,
                 font.tickslab = TickFont,
                 main = paste("Linearity Plot Featuring: ",Feature1," + ",Feature2, sep = ""),
                 ylab = paste(str_to_title(residType)," Residuals", sep = "")
      )
      p
      
    }
    
  })
  
  output$UnivarLinearityPlot <- renderPlot({
    
    if (length(input$SingleSurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SingleSurvivalFeature
      ref_Feature <- input$SurvFeatVariableUni
      scoreMethod <- input$ScoreMethod
      residType <- input$ResidualTypeUni
      AxisFont <- input$linAxisFont
      MainFont <- input$linMainFont
      TickFont <- input$linTickFont
      linpredict <- input$linPredict1
      #if (linpredict_choice == TRUE) {
      #  linpredict <- "linear.predictions"
      #}
      #else if (linpredict_choice == FALSE) {
      #  linpredict <- "observation.id"
      #}
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$UniVarNAcheck == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$UniVarNAcheck == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature],ignore.case = T, invert = T),]
      #  
      #}
      
      if (input$UniVarContCheck == FALSE) {
        
        meta_ssgsea_sdf[,Feature] <- as.factor(meta_ssgsea_sdf[,Feature])
        meta_ssgsea_sdf[,Feature] <- relevel(meta_ssgsea_sdf[,Feature], ref = ref_Feature)
        
      }
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature,sep = "")),
                   data = meta_ssgsea_sdf) 
      
      p <- ggcoxdiagnostics(tab,
                            type = residType,
                            sline = T,
                            sline.se = T,
                            ggtheme = theme_minimal(),
                            ox.scale = linpredict)
      p <- ggpar(p,
                 font.x = AxisFont,
                 font.y = AxisFont,
                 font.main = MainFont,
                 font.tickslab = TickFont,
                 main = paste("Linearity Plot Featuring: ",Feature, sep = ""),
                 ylab = paste(str_to_title(residType)," Residuals", sep = "")
      )
      p
      
    }
    
  })
  
  output$MultivarForestPlot <- renderPlot({
    
    if (length(input$SurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SurvivalFeature
      scoreMethod <- input$ScoreMethod
      forextFont <- input$ForestFontSize
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$MultiVarNAcheck == TRUE) {
        
        for (i in Feature){
          # Remove NA_unknown
          meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,i]) == FALSE),]
          meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,i],ignore.case = T, invert = T),]
        }
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$MultiVarNAcheck == TRUE) {
      #  
      #  for (i in Feature){
      #    
      #    meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,i]) == FALSE),]
      #    meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,i],ignore.case = T, invert = T),]
      #  }
      #  
      #}
      
      for (i in Feature){
        meta_ssgsea_sdf[,i] <- as.factor(meta_ssgsea_sdf[,i])
        colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[which(colnames(meta_ssgsea_sdf) == i)])
      }
      
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature,collapse = "+"),sep = "")),
                   data = meta_ssgsea_sdf)
      
      ggforest(tab,
               data = meta_ssgsea_sdf,
               main = paste("Hazard Ratio Modeling: ",paste(Feature,collapse = ", "),sep = ""),
               fontsize = forextFont)
    }
    
  })
  
  output$BivarForestPlot <- renderPlot({
    
    if (length(input$SurvivalFeatureBi1 > 0) & length(input$SurvivalFeatureBi2 > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1
      Feature2 <- input$SurvivalFeatureBi2
      Feat1Var <- input$SurvFeatVariableBi1
      Feat2Var <- input$SurvFeatVariableBi2
      forextFont <- input$ForestFontSize
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarAddNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarAddNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarAddNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarAddNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarAddContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarAddContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarAddContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarAddContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"+",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      
      ggforest(tab,
               data = meta_ssgsea_sdf,
               main = paste("Hazard Ratio Modeling: ",paste(Feature1,"+",Feature2,sep = ""),sep = ""),
               fontsize = forextFont)
    }
    
  })
  
  output$BivarForestPlotInter <- renderPlot({
    
    if (length(input$SurvivalFeatureBi1Inter > 0) & length(input$SurvivalFeatureBi2Inter > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      Feat1Var <- input$SurvFeatVariableBi1Inter
      Feat2Var <- input$SurvFeatVariableBi2Inter
      forextFont <- input$ForestFontSize
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$BiVarIntNAcheck1 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature1]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature1],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      if (input$BiVarIntNAcheck2 == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature2]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature2],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature1,Feature2)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$BiVarIntNAcheck1 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature1]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature1],ignore.case = T, invert = T),]
      #  
      #}
      #if (input$BiVarIntNAcheck2 == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature2]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature2],ignore.case = T, invert = T),]
      #  
      #}
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      meta_ssgsea_sdf[,4] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,4])
      Feature1 <- gsub("[[:punct:]]","_",Feature1)
      Feat1Var <- gsub("[[:punct:]]","_",Feat1Var)
      colnames(meta_ssgsea_sdf)[5] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[5])
      meta_ssgsea_sdf[,5] <- gsub("[[:punct:]]","_",meta_ssgsea_sdf[,5])
      Feature2 <- gsub("[[:punct:]]","_",Feature2)
      Feat2Var <- gsub("[[:punct:]]","_",Feat2Var)
      
      #meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
      #meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      #meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
      #meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      
      if (input$BiVarIntContCheck1 == FALSE) {
        meta_ssgsea_sdf[,Feature1] <- factor(meta_ssgsea_sdf[,Feature1])
        meta_ssgsea_sdf[,Feature1] <- relevel(meta_ssgsea_sdf[,Feature1], ref = Feat1Var)
      }
      else if (input$BiVarIntContCheck1 == TRUE) {
        meta_ssgsea_sdf[,Feature1] <- as.numeric(meta_ssgsea_sdf[,Feature1])
      }
      if (input$BiVarIntContCheck2 == FALSE) {
        meta_ssgsea_sdf[,Feature2] <- factor(meta_ssgsea_sdf[,Feature2])
        meta_ssgsea_sdf[,Feature2] <- relevel(meta_ssgsea_sdf[,Feature2], ref = Feat2Var)
      }
      else if (input$BiVarIntContCheck2 == TRUE) {
        meta_ssgsea_sdf[,Feature2] <- as.numeric(meta_ssgsea_sdf[,Feature2])
      }
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",paste(Feature1,"*",Feature2,sep = ""),sep = "")),
                   data = meta_ssgsea_sdf)
      
      
      ggforest(tab,
               data = meta_ssgsea_sdf,
               main = paste("Hazard Ratio Modeling: ",paste(Feature1,"*",Feature2,sep = ""),sep = ""),
               fontsize = forextFont)
    }
    
  })
  
  output$SinglevarForestPlot <- renderPlot({
    
    if (length(input$SingleSurvivalFeature > 0)) {
      ## Assign variables
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$SingleSurvivalFeature
      ref_Feature <- input$SurvFeatVariableUni
      scoreMethod <- input$ScoreMethod
      forextFont <- input$ForestFontSize
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
      
      quantCutoff <- input$QuantPercent/100 #Quantile cutoff given by user
      if (input$UniVarNAcheck == TRUE) {
        
        # Remove NA_unknown
        meta_ssgsea <- meta_ssgsea[which(is.na(meta_ssgsea[,Feature]) == FALSE),]
        meta_ssgsea <- meta_ssgsea[grep("unknown",meta_ssgsea[,Feature],ignore.case = T, invert = T),]
        ## Re-Perform Stat functions
        meta_ssgsea$VAR_Q <- quartile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuartilePScore <- paste("", meta_ssgsea$VAR_Q, sep="")
        meta_ssgsea$HiLoPScore <- highlow(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)])
        meta_ssgsea$QuantCutoffPScore <- quantile_conversion(meta_ssgsea[, which(colnames(meta_ssgsea) == geneset_name)], quantCutoff)
        
      }
      
      ## Subset columns needed for plot and rename for surv function
      select_cols <- c("SampleName",surv_time_col,surv_id_col,Feature)
      meta_ssgsea_sdf <- meta_ssgsea[,select_cols]
      
      #if (input$UniVarNAcheck == TRUE) {
      #  
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[which(is.na(meta_ssgsea_sdf[,Feature]) == FALSE),]
      #  meta_ssgsea_sdf <- meta_ssgsea_sdf[grep("unknown",meta_ssgsea_sdf[,Feature],ignore.case = T, invert = T),]
      #  
      #}
      
      if (input$UniVarContCheck == FALSE) {
        
        meta_ssgsea_sdf[,Feature] <- as.factor(meta_ssgsea_sdf[,Feature])
        meta_ssgsea_sdf[,Feature] <- relevel(meta_ssgsea_sdf[,Feature], ref = ref_Feature)
        
      }
      
      colnames(meta_ssgsea_sdf)[4] <- gsub("[[:punct:]]","_",colnames(meta_ssgsea_sdf)[4])
      Feature <- gsub("[[:punct:]]","_",Feature)
      
      
      ## Survival Function
      tab <- coxph(as.formula(paste("Surv(",surv_time_col,",",surv_id_col,") ~ ",Feature,sep = "")),
                   data = meta_ssgsea_sdf)
      
      
      ggforest(tab,
               data = meta_ssgsea_sdf,
               main = paste("Hazard Ratio Modeling: ",paste(Feature,collapse = ", "),sep = ""),
               fontsize = forextFont)
    }
    
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
  
  output$Sheatmap <- renderPlot({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    heatgenes <- geneset[[GeneSet]]
    meta <- SboxplotReact()
    expr_start <- exprSub()
    samples <- meta$SampleName
    expr <- expr_start[which(rownames(expr_start) %in% heatgenes),colnames(expr_start) %in% samples, drop = F]
    meta <- meta[which(meta$SampleName %in% colnames(expr)),]
    clmethod <- input$ClusterMethod
    rowfont <- input$heatmapFontR
    colfont <- input$heatmapFontC
    color_choice <- input$ColorPaletteHeat
    
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
    #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
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
    boxplotang <- input$boxplotTextAngle
    meta_ssGSEA <- ssGSEAmeta()
    boxTab <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    boxTab[,FeatureSelec] <- as.factor(boxTab[,FeatureSelec])
    
    ## Adjust 'Sample Type' for label 
    if (length(unique(meta[,metacol_sampletype])) > 1) {
      SampleTypeLab <- paste(" (",SampleType,") ",sep = "")
    }
    if (length(unique(meta[,metacol_sampletype])) <= 1) {
      SampleTypeLab <- " "
    }
    
    if (is.na(as.numeric(boxplotang)) == TRUE) {
      ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
        geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
        geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
        labs(x = "Group", y = paste(GeneSet," ",scoreMethod," Score", sep = ""),
             title = paste(GeneSet," Gene Set ",scoreMethod," Score: ","\n",
                           Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
             fill = FeatureSelec) +
        theme_bw() +
        stat_compare_means(method = StatMethod) +
        theme(text = element_text(size = font)) +
        scale_x_discrete(guide = guide_axis(n.dodge = 2))
    }
    
    else if (is.na(as.numeric(boxplotang)) == FALSE) {
      if (as.numeric(boxplotang) == 0) {
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
      }
      else if (as.numeric(boxplotang) == 45) {
        ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
          geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
          geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
          labs(x = "Group", y = paste(GeneSet," ",scoreMethod," Score", sep = ""),
               title = paste(GeneSet," Gene Set ",scoreMethod," Score: ","\n",
                             Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
               fill = FeatureSelec) +
          theme_bw() +
          stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                axis.text.x = element_text(angle = as.numeric(boxplotang), hjust = 1))
      }
      else if (as.numeric(boxplotang) == 90) {
        ggplot(boxTab, aes(factor(boxTab[,FeatureSelec]), boxTab[,GeneSet], fill = boxTab[,FeatureSelec])) +
          geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
          geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = dot) +
          labs(x = "Group", y = paste(GeneSet," ",scoreMethod," Score", sep = ""),
               title = paste(GeneSet," Gene Set ",scoreMethod," Score: ","\n",
                             Feature,SampleTypeLab,"Patients Featuring ",FeatureSelec,sep = ""),
               fill = FeatureSelec) +
          theme_bw() +
          stat_compare_means(method = StatMethod) +
          theme(text = element_text(size = font),
                axis.text.x = element_text(angle = as.numeric(boxplotang)))
      }
      
    }
    
  })
  
  output$FeatureHeatmap <- renderPlot({
    
    geneset <- gs_react()
    GeneSet <- names(geneset)
    FeatureSelec <- input$HeatmapFeature
    meta_ssGSEA <- ssGSEAmeta()
    meta <- meta_ssGSEA[,c("SampleName",FeatureSelec,GeneSet)]
    
    heatgenes <- geneset[[GeneSet]]
    expr_start <- exprSub()
    samples <- meta$SampleName
    expr <- expr_start[which(rownames(expr_start) %in% heatgenes),colnames(expr_start) %in% samples, drop = F]
    meta <- meta[which(meta$SampleName %in% colnames(expr)),]
    clmethod <- input$ClusterMethod
    rowfont <- input$heatmapFontR
    colfont <- input$heatmapFontC
    color_choice <- input$ColorPaletteHeat
    
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
    #hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
    #Heatmap color
    col_sets <- c("OrRd","PuBu","Greens","YlGnBu")
    if (color_choice == "original") {
      HeatMap_Colors <- c("dark blue","blue","white","red", "dark red")
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice %in% col_sets) {
      HeatMap_Colors <- brewer.pal(n = 5, color_choice)
      hmcols <- colorRampPalette(HeatMap_Colors)(length(bk)-1)
    }
    else if (color_choice == "Inferno") {
      hmcols <- inferno(500)
    }
    else if (color_choice == "Viridis") {
      hmcols <- viridis(500)
    }
    else if (color_choice == "Plasma") {
      hmcols <- plasma(500)
    }
    else if (color_choice == "OmniBlueRed") {
      hmcols<- colorRampPalette(c("#1984C5", "#22A7F0", "#63BFF0", "#A7D5ED", "#E2E2E2", "#E1A692", "#DE6E56", "#E14B31", "#C23728"))(length(bk)-1)
    }
    else if (color_choice == "LightBlueBlackRed") {
      hmcols<- colorRampPalette(c("#34C5FD","black","red"))(length(bk)-1)
    }
    else if (color_choice == "GreenBlackRed") {
      hmcols<- colorRampPalette(c("green","black","red"))(length(bk)-1)
    }
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
  
  ssgseaDensity_react <- reactive({
    
    geneset <- gs_react()
    geneset_name <- names(geneset)
    ssgsea_meta <- ssGSEAmeta()
    user_quant <- input$densityPercent/100
    ShowQuartile <- input$QuartileLinesCheck
    scoreMethod <- input$ScoreMethod
    
    cols_selec <- c("SampleName",geneset_name)
    ssgsea_scores <- ssgsea_meta[,cols_selec]
    
    quant_df <- data.frame(quantile(ssgsea_scores[,geneset_name]))
    quant_df2 <- quant_df[c(2,3,4),,drop = F]
    colnames(quant_df2)[1] <- "Quantile"
    
    user_vline <- quantile(ssgsea_scores[,geneset_name],probs = user_quant)
    
    if (input$GeneSetTabs == 2) {
      if (input$RawOrSS == "Raw Gene Expression") {
        scoreMethodLab <- "Raw Gene Expression Density"
      }
      else if (input$RawOrSS == "ssGSEA Rank Normalized") {
        scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
      }
    }
    else if (input$GeneSetTabs != 2) {
      scoreMethodLab <- paste(scoreMethod, " Score Density", sep = "")
    }
    
    p <- ggplot(ssgsea_scores, aes(x=ssgsea_scores[,geneset_name])) + 
      geom_density(color="darkblue", fill="lightblue", alpha = 0.4) +
      xlab(paste(scoreMethod,"Score")) +
      ylab(scoreMethodLab) +
      ggtitle(paste(colnames(ssgsea_scores)[2],scoreMethodLab)) +
      theme(axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 20))
    if (ShowQuartile == TRUE) {
      p <- p + geom_vline(data = quant_df2, aes(xintercept = Quantile), linetype = "dashed", color = "darkblue", size = 1)
    }
    if (user_quant != 0) {
      p <- p + geom_vline(xintercept = user_vline, linetype = "dashed", color = "darkred", size = 1)
    }
    p
    
  })
  
  output$ssgseaDensity <- renderPlot({
    
    p <- ssgseaDensity_react()
    p
    
    
  })
  
  ####----Text Output----####
  
  output$timewarnmessage1 <- renderUI({
    
    if (input$linPredict1 == "time") {
      p("Residual type must be schoenfeld or scaledsch.")
    }
    
  })
  
  output$timewarnmessage2 <- renderUI({
    
    if (input$linPredict2 == "time") {
      p("Residual type must be schoenfeld or scaledsch.")
    }
    
  })
  
  output$timewarnmessage3 <- renderUI({
    
    if (input$linPredict3 == "time") {
      p("Residual type must be schoenfeld or scaledsch.")
    }
    
  })
  
  
  ####----Downloaders----####
  
  
  
  ## quartile
  output$dnldSplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- Splot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldSplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuartileSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SplotBIN_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## binary
  output$dnldSplotBIN_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_HiLoPscoreSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_HiLoPscoreSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SplotBIN_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldSplotBIN_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_HiLoPscoreSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_HiLoPscoreSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SplotBIN_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## Quantile
  output$dnldSquantPlot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldSquantPlot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_QuantileSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## cutoff
  output$dnldSquantPlot2_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot2_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldSquantPlot2_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_AboveBelowCutoffSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- SquantPlot2_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ## univariate
  output$dnldfeatSplot_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- featSplot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldfeatSplot_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SingleSurvivalFeature
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,"_",geneset_name,"_",scoreMethodLab,"_UnivariateSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- featSplot_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  ##bivariate
  output$dnldfeatSplotBi_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.svg",sep = "")
      }
    },
    content = function(file) {
      p <- featSplotBi_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldfeatSplotBi_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      Feature1 <- input$SurvivalFeatureBi1Inter
      Feature2 <- input$SurvivalFeatureBi2Inter
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        SampleType <- input$SampleTypeSelection
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,Feature1,Feature2,"_",geneset_name,"_",scoreMethodLab,"_BivariateSurvival.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- featSplotBi_react()
      ggsave(file,p$plot,width = 10, height = 8)
      
    }
  )
  
  output$dnldssgseaDensity_SVG <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.svg",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.svg",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaDensity_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  
  output$dnldssgseaDensity_PDF <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      scoreMethod <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      if (input$GeneSetTabs == 2) {
        if (input$RawOrSS == "Raw Gene Expression") {
          scoreMethodLab <- "RawGeneExpression"
        }
        else if (input$RawOrSS == "ssGSEA Rank Normalized") {
          scoreMethodLab <- scoreMethod
        }
      }
      else if (input$GeneSetTabs != 2) {
        scoreMethodLab <- scoreMethod
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.pdf",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",scoreMethodLab,"_DensityPlot.pdf",sep = "")
      }
    },
    content = function(file) {
      p <- ssgseaDensity_react()
      ggsave(file,p,width = 10, height = 8)
    }
  )
  
  output$dnldssgseaDensityTable <- downloadHandler(
    filename = function() {
      geneset <- gs_react()
      geneset_name <- names(geneset)
      SampleType <- input$SampleTypeSelection
      Feature <- input$FeatureSelection
      score_method <- input$ScoreMethod
      if (Feature != "All_Features") {
        SubFeature <- paste("_",input$subFeatureSelection,"_",sep = "")
      }
      if (Feature == "All_Features") {
        SubFeature <- "_"
      }
      # If more than one sample type
      if (length(unique(meta[,metacol_sampletype])) > 1) {
        paste(gsub(" ","",ProjectName),"_",SampleType,"_",Feature,SubFeature,"_",geneset_name,"_",score_method,".txt",sep = "")
      }
      # If only one sample type
      else if (length(unique(meta[,metacol_sampletype])) <= 1) {
        paste(gsub(" ","",ProjectName),"_",Feature,SubFeature,"_",geneset_name,"_",score_method,".txt",sep = "")
      }
    },
    content = function(file) {
      geneset <- gs_react()
      GeneSet <- names(geneset)
      ssgsea_meta <- ssGSEAmeta()
      table <- ssgsea_meta[,c("SampleName",GeneSet)]
      write_delim(table,file,delim = '\t')
    }
  )
  
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
  
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)






