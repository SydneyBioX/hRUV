#'
#' @title Clean each data
#'
#' @description Filter and order metabolites across all batches
#'
#' @param dat_list A list of SummarizedExperiment data where each object
#' represents a batch.
#' @param threshold A threshold value between 0 and 1 of missingness to filter.
#' @param method A method to select metabolites across al batches.
#' @param assay An assay name to measure the missingness of the signals.
#'
#' @export clean
clean = function(dat_list, threshold = 0.5, method = c("intersect", "union"),
    assay = "raw", newAssay = assay) {
    # Order the rows of all data
    dat_list = lapply(dat_list, function(dat) {
        dat[order(rownames(dat)),]
    })

    if ("intersect" %in% method) {
        # Find the metabolites with missingness <= threshold
        keep_names = .select_mets(dat_list, threshold, assay)
        dat_list = .filter_mets(dat_list, keep_names)
        # Get a list of all metabolites across all batches
        met_names = .get_rownames(dat_list)
        met_names = Reduce(intersect, met_names)
        dat_list = .filter_mets(dat_list, met_names)
        dat_list = .imputeList(dat_list, assay = assay, newAssay = newAssay)
    } else if (method == "union") {
        # Get a list of all metabolites across all batches
        met_names = .get_rownames(dat_list)
        met_names = Reduce(union, met_names)
        dat_list = .add_NArow(dat_list, met_names)

        # Find which metabolites in which dat has missingness more than the threshold
        dat_list = .map_na(dat_list, assay, threshold, newAssay)

    }

    dat_list
}



.get_rownames = function(dat_list) {
    lapply(dat_list, function(dat) {
        rownames(dat)
    })
}

#' @importFrom SummarizedExperiment assay
.select_mets = function(dat_list, threshold, assay) {
    keep_mets = lapply(dat_list, function(dat) {
        mat = SummarizedExperiment::assay(dat, assay)
        rownames(dat)[-which(apply(mat, 1, function(i) {
            sum(is.na(i))/length(i) > threshold
        }))]
    })
    Reduce(intersect, keep_mets)
}

.filter_mets = function(dat_list, mets) {
    lapply(dat_list, function(dat) {
        dat = dat[which(rownames(dat) %in% mets),]
        dat[order(rownames(dat)),]
    })
}

#' @importFrom SummarizedExperiment SummarizedExperiment assay assayNames colData
#' @importFrom S4Vectors DataFrame
.add_NArow = function(dat_list, met_names) {
    lapply(dat_list, function(dat) {
        # For metabolites not in the data, add the rows with NAs
        if (any(!(met_names %in% rownames(dat)))) {
            n = length(which(!(met_names %in% rownames(dat))))
            mat = matrix(NA, ncol = ncol(dat), nrow = n)
            colnames(mat) = colnames(dat)
            rownames(mat) = met_names[which(!(met_names %in% rownames(dat)))]
            mat_list = rep(list(mat), length(assayNames(dat)))
            names(mat_list) = assayNames(dat)
            tmp = SummarizedExperiment::SummarizedExperiment(
                assay = mat_list,
                colData = colData(dat),
                rowData = DataFrame(
                    metabolite = met_names[which(!(met_names %in% rownames(dat)))]
                )
            )
            dat = rbind(dat, tmp)
        }

        dat[order(rownames(dat)),]
    })
}


#' @importFrom SummarizedExperiment assay
#' @importFrom DMwR2 knnImputation
.imputeList = function(dat_list, assay, newAssay) {
    lapply(dat_list, function(dat) {
        .impute(dat, assay, newAssay)
    })
}

.impute = function(dat, assay, newAssay) {
    SummarizedExperiment::assay(dat, newAssay, withDimnames = FALSE) =
        SummarizedExperiment::assay(dat, assay) %>%
        t() %>%
        as.data.frame() %>%
        DMwR2::knnImputation() %>%
        t()
    dat
}


.map_na = function(dat_list, assay, threshold, newAssay) {
    lapply(dat_list, function(dat) {
        mat = SummarizedExperiment::assay(dat, assay)
        # Separate rows with missingness > threshold
        na_rows = which(apply(mat, 1, function(i) {
            sum(is.na(i))/length(i) > threshold
        }))
        na_mat = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))
        if (length(na_rows)) {
            cur_dat = dat[-na_rows,]
            na_dat = dat[na_rows,]
            # Impute the rows with missingness <= threshold
            cur_dat = .impute(cur_dat, assay, newAssay)
            SummarizedExperiment::assay(na_dat, newAssay) = SummarizedExperiment::assay(na_dat, assay)
            # rbind the split data
            dat = rbind(cur_dat, na_dat)
            dat = dat[order(rownames(dat)),]

            # Map the position of NAs
            cur_map = lapply(na_rows, function(row) {
                which(is.na(mat[as.numeric(row),]))
            })
            mat = SummarizedExperiment::assay(dat, newAssay)
            for (met in names(cur_map)) {
                for (s in cur_map[[met]]) {
                    mat[which(rownames(mat) == met), s] = mean(mat[,s], na.rm = TRUE)
                    na_mat[which(rownames(mat) == met), s] = TRUE
                }
            }
            SummarizedExperiment::assay(dat, newAssay, withDimnames = FALSE) = mat
            SummarizedExperiment::assay(dat, "na_id", withDimnames = FALSE) = na_mat
            dat
        } else {
            SummarizedExperiment::assay(dat, "na_id", withDimnames = FALSE) = na_mat
            dat
        }

    })
}
