% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotPCA}
\alias{plotPCA}
\title{Principal Component Analysis (PCA) plot}
\usage{
plotPCA(dat, assay, colour, ...)
}
\arguments{
\item{dat}{SummarizedExperiment object data}

\item{assay}{An assay to use in a SummarizedExperiment object}

\item{colour}{A group vector to colour the points in the PCA plot}

\item{...}{An additional parameter to be used for PCA (\code{prcomp} from
stats package)}
}
\value{
PCA plot
}
\description{
Principal Component Analysis (PCA) plot
}
