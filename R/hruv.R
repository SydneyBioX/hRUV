#' @title hRUV - intra- and hierarchical inter-batch normalisation.
#'
#' @description Hirarchical approach to remove unwanted variation across
#' multiple batches of data.
#'
#' For \code{intra} method, \code{rlm} or \code{loess} methods fits robust
#' linear model or loess respectively to samples selected from \code{pCtName}.
#' \code{rlmShort} and \code{loessShort} both performs robust smoother \code{rlm}
#' and \code{loess} respectively then performs RUV-III utilising the replicates
#' identified from \code{intra_rep}. Parameters \code{negCtl} and \code{intra_k}
#' will be used for RUV-III.
#'
#'
#' @param dat_list A list of SummarizedExperiment data where each object
#' represents a batch.
#' @param assay An assay name to measure the missingness of the signals.
#' @param intra An intra-batch normalisation method.
#' See Description for details.
#' @param inter An inter-batch normalisation method. Either "balanced" or
#' "concatenate" method for hierarchical approach.
#' @param intra_k An intra-batch RUV-III parameter k value.
#' @param inter_k An inter-batch RUV-III parameter k value.
#' @param pCtlName A name of the variable in the colData of each
#' SummarizedExperiment data in the \code{dat_list}. This variable in colData
#' should be a numeric ID of which samples to use for robust smoother fitting.
#' @param negCtl A vector of row IDs for selection of stable matbolites to be
#' used as negative control in RUV.
#' @param intra_rep A name of the variable in the colData of each
#' SummarizedExperiment data in the \code{dat_list}. This variable in colData
#' should be a logical vector of which samples are a intra-batch replicate
#' @param inter_rep A name of the variable in the colData of each
#' SummarizedExperiment data in the \code{dat_list}. This variable in colData
#' should be a logical vector of which samples are a inter-batch replicate
#' @param hOrder A vector of batch names from names of list \code{dat_list}.
#' @param newAssay A name of the new assay for cleaned (preprocessed) data.
#'
#' @return A normalised data as SummarizedExperiment object.
#'
#'
#' @export
hruv = function(dat_list, assay,
    intra = c("rlm", "loess", "rlmShort", "loessShort"),
    inter = c("balanced", "concatenate"), intra_k = 5, inter_k = intra_k, pCtlName,
    negCtl, intra_rep, inter_rep, hOrder = NULL, newAssay = NULL) {

    # dat_list = lapply(dat_list, function(dat) {
    #     intraNorm(dat, assay = assay, method = intra, pCtl = which(dat[[pCtlName]]), rep = intra_rep, k = intra_k, negCtl = negCtl)
    # })
    dat_list = intraNorm(dat_list, assay = assay, method = intra,
        pCtlName = pCtlName, rep = intra_rep, k = intra_k, negCtl = negCtl)

    balanced = NULL
    if (inter == "balanced") {
        balanced = TRUE
    } else if (inter == "concatenate") {
        balanced = FALSE
    }
    if (is.null(balanced)) {
        stop("Incorrect inter normalisation method is selected")
    }

    dat = hierarchy(dat_list, assay = intra, k = inter_k, newAssay =
            newAssay, balanced = balanced,
        negCtl = negCtl, hOrder = hOrder, rep = inter_rep)

    dat


}

