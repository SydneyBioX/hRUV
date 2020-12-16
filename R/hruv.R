#' @export
hruv = function(dat_list, assay,
    intra = c("rlm", "loess", "rlmShort", "loessShort"),
    inter = c("balanced", "concatenate"), intra_k = 5, inter_k = intra_k, pCtlName,
    negCtl, inter_rep, intra_rep, hOrder = NULL, newAssay = NULL) {

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

