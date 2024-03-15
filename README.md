# Bull-et-al

## **NOTE** 
This repository consists of `dualGSEA()` R function created for the preprint BioRxiv submission:
"**Evaluation of Gene Set Enrichment Analysis (GSEA) tools highlights the value of single sample approaches over pairwise for robust biological discovery.**"


## **IMPORTANT**

Before you run the function `dualGSEA()`, please make sure you have all the dependencies
packages installed and the packages has been loaded. See below:


These are the list of required packages.
```r
### Required R packages ----
bio_pkgs <- c("limma", "msigdbr", "DOSE", "devtools","fgsea", "enrichplot","GSVA","circlize",
              "ComplexHeatmap","ggridges","pROC","cutpointr","ggalluvial","waterfalls",
              "randomForest","devtools","tibble","msigdbr","data.table","tidyverse","dplyr",
              "ggplot2","grid","gridExtra","ROCR","reshape2","stats","RColorBrewer","ggpubr",
              "ggbeeswarm","tidyr","caret","escape", "doRNG")
github_pkgs <- c("kassambara/easyGgplot2")
```

If not already installed, please do so as follows:
```r
### If not already installed ----
suppressMessages({
  if (!require(bio_pkgs, quietly = TRUE))
    install.packages(bio_pkgs)
  if (!require(github_pkgs, quietly = TRUE))
    devtools::install_github(github_pkgs, dependencies = TRUE, force=TRUE)  
})
```

Lastly, load the R packages before you have your input data ready and run the function!
```r
### Load R packages ----
suppressMessages(
  lapply(bio_pkgs, library, character.only = TRUE)
)
suppressMessages(
  library(easyGgplot2)
)
```

