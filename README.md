# hRUV

`hRUV` is a package for normalisation of multiple batches of metabolomics data in a hierarchical strategy with use of samples replicates in a large-scale studies. The tool utilises 2 types of replicates: intra-batch and inter-batch replicates to estimate the unwatned variation within and between batches with RUV-III. We have designed the replicate embedding arrangements within and between batches from http://shiny.maths.usyd.edu.au/hRUV. Our novel tool is a novel hierarchical approach to removing unwanted variation by harnessing information from sample replicates embedded in the seequence of experimental runs/batches and applying signal drift correction with robust linear or non-linear smoothers.

## hRUV overview

<img src="https://github.com/SydneyBioX/hRUV/blob/main/inst/overview_v1.png" align="center" width=100% />

## Software requirements

This package has been tested for *macOS Big Sur (11.1)* and *Linux Debian 10 (buster)* with R version 4.0.3.

### R package dependecies

```
SummarizedExperiment
dplyr
tidyr
rlang
ggplot2
plotly
S4Vectors
tibble
DMwR2
MASS
impute
```

## Installation

Install the R package from GitHub using the `devtools` package:

```r
if (!("devtools" %in% rownames(installed.packages())))
    install.packages("devtools")

library(devtools)
devtools::install_github("SydneyBioX/hRUV", build_vignettes=TRUE)
```

## Vignette

The vignette of the package is available under the `vignettes` folder for detailed demontration to hRUV package.

You can also open the compiled vignette from R console with `browseVignette("hRUV")`.

## Contact us

If you have any enquiries, especially about performing hRUV to integrate your metabolomics data, please contact taiyun.kim@sydney.edu.au. We are also happy to receive any suggestions and comments.
