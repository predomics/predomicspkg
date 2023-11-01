# Installation

Follow these steps to install the Predomics package:

1.  Install R from [CRAN](https://cran.r-project.org/).
2.  Install the required dependencies: 

```
devtools::install_github("predomics/predomicspkg")

## install dependencies
install.packages(c("doSNOW", "foreach", "snow", "doRNG", "gtools", "glmnet", "pROC", "viridis", "kernlab", "randomForest","effsize"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioQC")
# install.packages("testthat")
# install.packages("roxygen2")
```

Once eveything is installed un the check ...