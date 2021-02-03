# hRUV

`hRUV` is a package for normalisation of multiple batches of metabolomics data in a hierarchical strategy with use of samples replicates in a large-scale studies. The tool utilises 2 types of replicates: intra-batch and inter-batch replicates to estimate the unwatned variation within and between batches with RUV-III. We have designed the replicate embedding arrangements within and between batches from http://shiny.maths.usyd.edu.au/hRUV. Our novel tool is a novel hierarchical approach to removing unwanted variation by harnessing information from sample replicates embedded in the seequence of experimental runs/batches and applying signal drift correction with robust linear or non-linear smoothers.

## hRUV overview

<img src="https://github.com/SydneyBioX/hRUV/blob/main/inst/overview_v1.png" align="center" />

## Installation

Install the R package from GitHub using the `devtools` package:

```r
if (!("devtools" %in% rownames(installed.packages())))
    install.packages("devtools")

library(devtools)
devtools::install_github("SydneyBioX/hRUV")
```



## Contact us

If you have any enquiries, especially about performing hRUV to integrate your metabolomics data, please contact taiyun.kim@sydney.edu.au. We are also happy to receive any suggestions and comments.
