# DRPPM - SURVIVE Shiny App

[Insert Intro]

# Installation

## Via Download

## Via Git Clone

1. Clone the [GitHub Repository](https://github.com/shawlab-moffitt/DRPPM-SURVIVE.git) into the destination of your choice.
   * Can be done in R Studio Terminal or a terminal of your choice
```bash
git clone https://github.com/shawlab-moffitt/DRPPM-SURVIVE.git
```
2. Set your working directory in R to the cloned repository
   * This can be done through the "More" settings in the bottom-right box in R Stuido
   * You may also use the `setwd()` function in R Console.

# Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.7.1 | shinythemes_1.2.0 | shinyjqui_0.4.1 | shinycssloaders_1.0.0 | gtsummary_1.6.0 |
| dplyr_1.0.9 | tidyr_1.1.3 | readr_2.1.2 | tibble_3.1.7 | DT_0.23 |
| ggplot2_3.3.6 | survival_3.2-11 | survminer_0.4.9 | pheatmap_1.0.12 | gridExtra_2.3 |
| GSVA_1.40.1 | clusterProfiler_4.0.5 | ggpubr_0.4.0 |  |  |

# Required Files

* **Expression Matrix (.tsv/.txt):**
  * Must be tab delimited with gene names as symbols located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column.
  *  The App expects lowly expressed genes filtered out and normalized data either to FPKM or TMM.
     * Larger files might inflict memory issues for you local computer.
  *  An example file can be found via this [link](http://shawlab.science/shiny/DRPPM_SURVIVE/Pan_ICI_iAtlas_ExpressionMatrix.zip) due to the size being too large to store in the github repository (157MB zipped).

* **Meta Data (.tsv/.txt):**
  * This should be a tab delimited file with each row depicting a sample by the same name as in the expression matrix followed by informative columns containing survival data and other features to analyze the samples by.
  * Required columns:
    * A sample name column
    * A survival time and ID column (can be more than one type of survival time)
    * Feature column(s) that allow for grouping of samples for analysis
  * Optional Column
    * A sample type column that allows for an initial subsetting of samples followed by grouping by feature (ex. Tissue or Disease Type)
    * Description column(s) that give additional information on the samples
  * An example file can be found here: [Pan_ICI_iAtlas_MetaData.txt](https://github.com/shawlab-moffitt/DRPPM-SURVIVE/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData.txt)

* **Meta Data Parameters (.tsv/.txt):**
  * This should be a two-column tab-delimited file with the first column containing the column names of the meta file and the second column containing the column type of that meta column
  * The first column containing the users meta column names can be named however the user chooses, the **second column must have names matching the format provided below**.
    * For example, the survival ID column containing the 0/1 for events in the meta data could be nammed "OS" or "ID" or "OS.ID", but the column type in the second column meta parameter file must say "SurvivalID" for that specific column.
  * Column types and exmplinations:
    * **SampleName:** Contains sample names matching exprssion data (ONLY ONE ALLOWED)
    * **SampleType:** Contains a way to group and subset samples for further analysis (ONLY ONE ALLOWED and OPTIONAL)
    * **SurvivalTime:** Contains the overall survival time for the samples (can be other types of survival)
    * **SurvivalID:** Contains the survival ID for the samples, should be in a 0/1 format, 0 for alive/no event or 1 for dead/event (can be other types of survival)
    * **Feature:** Contains a feature that allows the samples to be grouped for analysis (More than one feature column allowed)
    * **Description:** Contains descriptions for the samples that may be viewed in the app (OPTIONAL)
  * Below is an example of these catagories and a full example file can be found here: [Pan_ICI_iAtlas_MetaData_Params.txt](https://github.com/shawlab-moffitt/DRPPM-SURVIVE/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt)

|  |  |
| --- | --- |
| sample | SampleName |
| Cancer_Tissue | SampleType |
| OS.time | SurvivalTime |
| OS.ID | SurivalID |
| EFS.time | SurvivalTime |
| EFS.ID | SurivalID |
| Responder | Feature |
| Race | Feature |
| Age | Description |
| Center | Description |

* **Gene Set File (.RData/.gmt/.tsv/.txt):**
  * This is the file that contains the gene set names and genes for each gene set.
  * This file is provided but can be replaced for a file of the users choice.
  * An .RData list is the preferred format which is a named list of gene sets and genes. A script to generate this list is provided here: [GeneSetRDataListGen.R](https://github.com/shawlab-moffitt/DRPPM-SURVIVE/blob/main/GeneSet_Data/GeneSetRDataListGen.R)
    * The app also accepts gene sets in .gmt format or two-column tab-delimited .tsv/.txt format with the first column being the gene set name repeating for every gene symbol that would be placed in the second column. If either of these three formats are given athe app with automatically convert them to an RData list.
  * The RData list provided ([GeneSet_List_HS.RData](https://github.com/shawlab-moffitt/Survival_Analysis_Shiny_App/blob/main/SurvivalAnalysis_ExampleApp/GeneSet_Data/GeneSet_List_HS.RData)) contains over 94K gene sets from MSigDB, LINCS L1000 Cell Perturbations, and Cell Marker databases.
 
* **Gene Set Master Table (Optional):**
  * This is a optional three-column tab-delimited table the catagorizes and subcatagorizes the gene sets of the provided [GeneSet_List_HS.RData file](https://github.com/shawlab-moffitt/DRPPM-SURVIVE/blob/main/GeneSet_Data/GeneSet_List_HS.RData)
  * It allows for organization of the large gene set list in the UI of the gene set selection for the Shiny App.
  * If not provided or using a user-provided gene set file, the gene set selection table only contains the gene set names.
  * The gene set master table that is provided can be found here: [GeneSet_CatTable.tsv](https://github.com/shawlab-moffitt/DRPPM-SURVIVE/blob/main/GeneSet_Data/GeneSet_CatTable.tsv)

# App Set-Up


# App Features

## Sidebar Panel

### Sample Selection and Parameters


### Survival Parameters


### Figure Parameters


### Meta Data


## Main Panel

### Survival Plot


### Survival Box Plot


### Survival Heatmap


### Feature Box Plot


### Feature Heatmap


# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.
