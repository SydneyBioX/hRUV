#' Construct SummarizedExperiment object from a matrix
#'
#' @param dat A p*n data frame or a matrix object where p is the number of
#' metabolites and n is the number of samples
#' @param assay A name of assay to be labelled
#' @param name Name of the data
#' @param colData A DataFrame object. This should include type of replicates,
#' sample name and etc.
#' @param rowData A DataFrame object. This should include details of the metabolites.
#'
#' @return A SummarizedExperiment object
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export constructSe
constructSe = function(dat, assay = "raw", name = NULL,
    colData = NULL, rowData = NULL) {

    assays = list(as.matrix(dat))
    names(assays) = assay

    if (is.null(name)) {
        name = "data"
    }
    SummarizedExperiment::SummarizedExperiment(
        assays = assays,
        colData = colData,
        rowData = rowData,
        metadata = list(name = name)
    )
}
