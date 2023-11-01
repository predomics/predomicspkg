---
output: pdf_document
---

The **predomics** package contains three methods for suppervised learning based on ternary coefficients.
============================================================

##UPDATES
* 17/05/2016: package creation and compilation
* 18/05/2016: git project creation and merging of Lucas and Edi's version
* 1/06/2016: New population (denseVect) for terGA and different rewritten operators.
* 3/06/2016: Digesting and comparative plot capability added.

```{r}
## install dependencies
install.packages(c("doSNOW", "foreach", "snow", "doRNG", "gtools", "glmnet", "pROC", "viridis", "kernlab", "randomForest","effsize"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioQC")
# install.packages("testthat")
# install.packages("roxygen2")
```

