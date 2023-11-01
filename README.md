# Predomics - Interpretable machine learning for omics data

The **predomics** package offers access to an original Machine Learning framework implementing several heuristics that allow discovering sparse and interpretable models in large datasets. These models are efficient and adapted for classification and regression tasks in metagenomics and other datasets with commensurable variables. We introduce the custom BTR (BIN, TER, RATIO) languages that describe different types of associations between variables. Moreover, in the same framework we implemented several state-of-the-art methods (SOTA) including RF, ENET and SVM. The `predomics` package started in 2015 and has evolved quickly since. A major improvement came in 2023. The package comes also with **predomicsapp**, a R Shiny application for easy training and exploration of results.

## Badges

![R Package](https://github.com/predomics/predomicspkg/workflows/R/badge.svg) ![License](https://img.shields.io/github/license/predomics/predomicspkg) ![GitHub issues](https://img.shields.io/github/issues/predomics/predomicspkg) [![GitHub forks](https://img.shields.io/github/forks/predomics/predomicspkg?style=social)](https://github.com/predomics/predomicspkg/network/members) [![GitHub stars](https://img.shields.io/github/stars/predomics/predomicspkg?style=social)](https://github.com/predomics/predomicspkg/stargazers)

## Table of Contents

1.  [Installation](docs/Installation.md)
2.  [Usage](docs/Usage.md)
3.  [Features](docs/Features.md)
4.  [Screenshots](docs/Screenshots.md)
5.  [Technologies](docs/Technologies.md)
6.  [Contributing](CONTRIBUTING.md)
7.  [License](LICENSE)
8.  [Authors](docs/Authors.md)
9.  [Contact](docs/Contact.md)
10. [FAQs](docs/FAQs.md)

## Installation {#installation}

Steps to install the project.

## Usage {#usage}

Examples and use cases.

## Features

-   Feature 1
-   Feature 2
-   ...

## Screenshots/Demo

![Screenshot](path/to/screenshot.png)

## Technologies Used

-   Technology 1
-   Technology 2
-   ...

## Contributing

How to contribute to this project.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors and Acknowledgment

Authors, contributors, and acknowledgments.

## Contact

Contact information.

## FAQs

Answers to common questions.

# Introduction

## Overview

We introduce here the *predomics* package, which is designed to search for simple and interpretable predictive models from omics data and more specifically metagenomics. These models, called BTR (for Bin/Ter/Ratio) are based on a novel family of languages designed to represent the microbial interactions in microbial ecosystems. Moreover, in this package we have proposed four different optimization heuristics to find some of the best predictive models. A model in *predomics* is a set of indexes from the dataset (i.e. variables) along with the respective coefficients belonging to the ternary set `{-1, 0, 1}` and an intercept of the form `(A + B + C - K - L - M < intercept)`. The number of variables in a model, also known as *model size*, *sparsity* or *parsimony*, can vary in a range provided as a parameters to the experiment.

In *predomics* we have impemented the following types of object:

-   `model`, which can be tested with `isModel()`, is a list which contains information on the features, the languages, the fitting scores etc.
-   `population of models`, which can be tested with `isPopulation()` is a list of `model` objects.
-   `model collection`, which can be tested with `isModelCollection()` is a list of `population` objects. They are grouped by model size.
-   `classifier`, which can be tested with `isClf()` is a set of parameters which defines a learner ready to be run.
-   `experiment`, which can be tested with `isExperiment()` is a top level object which contains a `classifier` object along with the learned models organized as a `model collection` object.

All these objects can be viewed with the `printy()` function. Other existing functions allow conversion from one object type to the other as for instance `modelCollectionToPopulation()`. An experiment can be explored using the `digest()` routine along with many other functions implemented in more than 18K lines of code.

## Heuristics

In this package we have proposed four different heuristics to search for the best predictive models.

-   **Terga1 and terga2**. These algorithms are based on genetic algorithms. Here we introduce the notion of a population of models, which is a set of individuals/models that can be mutated, crossed, evolved and selected for many generations (epochs). The main difference between *terga1* and *terga2* is that the former will evolve individuals of a fixed model-size, while the latter can evolve models of different model-size in larger populations. There are also other differences in the core of these algorithms but that we will not discuss here.
-   **Terbeam**. This algorithm consists of using a beam search approach. In computer science, beam search is a heuristic search algorithm that explores a graph by expanding a subset of the most promising node. Beam search is an optimization of best-first search that reduces its memory requirements. Here we use a window of a given size `(a, b)` for two consecutive model-size `(k, k+1)`. The best models of model-size `k` are used to generate combinations of model-size `k+1` and the best ones are kept for the next round. The size of the windows is fixed in parameters at the beginning of the experiment. The results can be also considered as a population of models.
-   **Terda**. This algorithm is based on a standard logistic regression approach. After learning a linear model with real-valued weights with the GlmNet package, these weights are rounded to either `{-1, 0, 1}`.
-   **Metal**. The four algorithms presented above complement each other in terms of model space epxloration. For instance, *terbeam* is specialized in small models, *terda* in bigger ones and *terga* will make sure to explore random combinations that are not only composed of features which in smaller-sparisty settings yield good results. We devised another algorithm named *metal* (fusion), which is a meta-heuristic based on the *terga2* engine. The initial population is seeded by sets of models found by either *terbeam*, *terda* or *terga2* and is next combined and evolved together. Metal has also the capacity to evolve models with multiple languages in the same time. This is particularly of interest when we wish to discover automatically the set of rules that predicts best the studied condition.

## Predomics languages

A *predomics* model is coded in R as a S3 object, which contains a certain number of attributes among which the `learner` (algorithm) that generated it but also the `language` that is used. The languages we have proposed in the current version are the following.

-   **Bin/bininter**. Let us take the following exemple `(A + B + C < intercept)`. In a *bin/bininter* language we have only two coefficients from the binary set `{0, 1}`. Features that do not appear in the model have the coefficient `0` while the others have the coefficient `1`. The difference between *bin* and *bininter* is that the intercept for the former is set to zero. In this tutorial we will consider *bin* and *bininter* as the same language and will search for *bininter*, which encompasses *bin*.
-   **Ter/teriter**. If we take the example above `(A + B + C - K - L - M < intercept)` we can see that in this model of size `k=6`, we have coefficients from the ternary set `{-1, 0, 1}`. Features that do not appear in the model have the coefficient `0` while the others have `1` or `-1`. Models that have only positive or only negative coefficients are not considered as they would be *bin* models. The difference between *ter* and *terinter* is that the intercept for the former is set to zero. In this tutorial we will consider *ter* and *terinter* as the same language and will search for *terinter*, which encompasses *ter*.
-   **Ratio**. In the ratio language the intercept plays the role of a multiplication factor and the model is the form `(A+B+C)/(K+L+M) < intercept`.
-   **State of the art (SOTA)**. We have implemnted in *predomics* three state-of-the-art algorithms used today in the field of metagenomics, random forest *(rf)*, support vector machines *(svmlin, svmrad)*, respectively using a linear and a radial kernel and finally the Elastic-Net Regularized Generalized Linear Models *(glmnet)*. Most of the parameteres of these algorithms are left by default and some are optimized using internal cross-validation techniques.

##UPDATES \* 17/05/2016: package creation and compilation \* 18/05/2016: git project creation and merging of Lucas and Edi's version \* 1/06/2016: New population (denseVect) for terGA and different rewritten operators. \* 3/06/2016: Digesting and comparative plot capability added.

```{r}
## install dependencies
install.packages(c("doSNOW", "foreach", "snow", "doRNG", "gtools", "glmnet", "pROC", "viridis", "kernlab", "randomForest","effsize"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BioQC")
# install.packages("testthat")
# install.packages("roxygen2")
```
