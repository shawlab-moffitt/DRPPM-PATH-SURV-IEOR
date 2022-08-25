# Guide to Immune Deconvolution Pre-Processing


# Intro
why immune deconvolution is powerful
why pre-rocess

# R Requirments

* `R` - https://cran.r-project.org/src/base/R-4/
   * R version >= 4.1 is required for this script
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

|  |  |  |  |
| --- | --- | --- | --- |
| dplyr_1.0.9 | readr_2.1.2 | stringr_1.4.0 | immunedeconv_2.1.0 |

## Immune Deconvolution Package

* Please go to this link for further information on package installation: https://github.com/omnideconv/immunedeconv
* Installation of this package could take up to 30 minutes depending on how many dependency packages need to be compiled
* The immune deconvolution methods CIBERSORT and CIBERSORT_abs require additional files which can be obtained from https://cibersortx.stanford.edu/
  * CIBERSORT.R and LM22.txt files are required
  * A license is required to download these files and may be obtained for free, upon approval
  * The script may still run without these files, but the CIBERSORT deconvolution method will not process
  
# Required Files

* **Expression Martix (.txt/.tsv):**
  * Must be tab delimited with gene names as symbols located in the first column with subsequent columns consiting of the sample name as the header and expression data down the column.
  * Example file here: [Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip)

* **Meta Data (.txt/.tsv):**
  * This file is optional, but recommended if the user plans to view the immune deconvolution results in the [DRPPM-PATH-SURVEIOR App](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR)
    * Including this file allows the script to output an updated meta file with the immune deconvolution result columns added
  * This should be a tab delimited file with each row depicting a sample by the same name as in the expression matrix followed by informative columns containing survival data and other features to analyze the samples by.
  * An example file here [Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt)

* **Meta Data Parameters (.txt/.tsv):**
  * This file is optional, but recommended if the user plans to view the immune deconvolution results in the [DRPPM-PATH-SURVEIOR App](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR)
    * Including this file allows the script to output an updated meta parameter file with additional immune deconvolution feature rows appended to the bottom
    * If this is not included, the script will output the feature rows and the user could append them to the parameter file manually
  * This should be a two-column tab-delimited file with the first column containing the column names of the meta file and the second column containing the column type of that meta column
    * Described in further detail [here](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR#required-files---user-provided)
  * And example file here: [Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR/blob/main/Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt)

# Set-Up

* Input desired Project name, file names and paths, and an output path for the processed data to be held
```{r}
####----User Input----####
ProjectName <- "PAN_ICI_Skin_Kidney"
Expression_Matrix_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Expression.zip"
Meta_Data_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_Skin_Kidney_Meta.txt"
Meta_Data_Param_File <- "Pan_ICI_Example_Data/Pan_ICI_iAtlas_MetaData_Params.txt"
Output_Path <- "Pan_ICI_Example_Data/"
```
* Denote `TRUE` or `FALSE` for which methods you want to run through in the script
  * `mcp_counter` and `estimate` are ran in the [DRPPM-PATH-SURVEIOR App](https://github.com/shawlab-moffitt/DRPPM-PATH-SURVEIOR) upon start up if it is not included in the app sart-up files
```{r}
quantiseq <- TRUE
mcp_counter <- TRUE
xcell <- TRUE
epic <- TRUE
abis <- TRUE
estimate <- TRUE
```
* If the user would like to run the CIBERSORT methods please follow the instructions to obtain the proper files 










