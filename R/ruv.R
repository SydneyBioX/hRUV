#' @importFrom SummarizedExperiment assay
getRUV = function(dat, assay, k = 5, rep, negCtl = NULL, newAssay = NULL, replicate_mat = NULL) {
    # Negative control metabolites
    if (is.null(negCtl)) {
        negCtl = seq(nrow(dat))
    }
    if (is.null(replicate_mat)) {
        if (rep == "pool") {
            replication = dat$sample_name
            replication[which(dat[[rep]])] = "Pool"
            replicate_mat = replicate.matrix(replication)
        } else if (rep == "batch_replicate") {
            replication = dat$sample_name
            replication = gsub("\\*{3}", "",replication)

            replication[dat$pool] = paste("Pool", 1:sum(dat$pool))
            replicate_mat = replicate.matrix(replication)

        } else if (rep == "short_replicate") {
            replication = dat$sample_name
            replicates = dat$sample_name[dat$short_replicate]

            replicates = gsub("\\*", "", replicates[which(grepl("^[0-9]+\\*$", replicates))])

            replicates = paste("(^", paste(replicates, collapse = "\\*{0,1}$)|(^"), "\\*{0,1}$)", sep = "")
            replication[which(grepl(replicates, replication))] = gsub("\\*", "", replication[which(grepl(replicates, replication))])

            # replicates = dat$sample_name[dat$short_replicate]
            # replicates = paste("(^", paste(replicates, collapse = "\\*{1}$)|(^"), "\\**$)", sep = "")
            # replication[which(grepl(replicates, replication))] = gsub("\\*", "", replication[which(grepl(replicates, replication))])


            replicate_mat = replicate.matrix(replication)

        }
    }

    dataRUV = t(SummarizedExperiment::assay(dat, assay))
    # if (adjust == "mean") {
    #     sample_grand_mean = mean(dataRUV)
    #     dataRUV = dataRUV - sample_grand_mean
    # } else if (adjust == "rowMean") {
    #     sample_grand_mean = colMeans(dataRUV)
    #     dataRUV = t(t(dataRUV) - sample_grand_mean)
    # }
    means = colMeans(dataRUV)
    data_stand <- apply(dataRUV, 1, function(x) x - means)


    ruvCorrected = RUVIII(Y = data_stand, M = replicate_mat, ctl = negCtl, return.info = TRUE, k = k)

    ruvCorrected$newY <- apply(ruvCorrected$newY, 1, function(x) x + means)


    # if (adjust == "mean") {
    #     ruvCorrected$newY = ruvCorrected$newY + sample_grand_mean
    # } else if (adjust == "rowMean") {
    #     ruvCorrected$newY = t(t(ruvCorrected$newY) + sample_grand_mean)
    # }


    if (is.null(newAssay)) {
        newAssay = assay
    }
    SummarizedExperiment::assay(dat, newAssay) = t(ruvCorrected$newY)
    dat
}


#' @title Hierarchical RUV
#'
#' @param dat_list A list of SummarizedExperiment data where each object
#' represents a batch.
#' @param assay An assay name to measure the missingness of the signals.
#' @param k Parameter \code{k} for RUV-III.
#' @param rep A name of the variable in the colData of each SummarizedExperiment
#'  data in the \code{dat_list}. This variable in colData
#' should be a logical vector of whic samples are a intra-batch replicate
#' @param hOrder A vector of batch names from names of list \code{dat_list}.
#' @param negCtl A vector of row IDs for selection of stable matbolites to be
#' used as negative control in RUV.
#' @param newAssay A name of the new assay for cleaned (preprocessed) data.
#' @param balanced A type of hierarchical approach, "concatenate" or "balanced".
#' Default is "concatenate" approach.
#'
#' @return A normalised data as SummarizedExperiment object.
#'
#' @importFrom SummarizedExperiment assay
#'
#' @export
hierarchy = function(dat_list, assay, k = 5, rep, newAssay = NULL, balanced = FALSE, negCtl = NULL, hOrder = NULL) {
    # Check the order of batches
    if (is.null(hOrder)) {
        if (is.null(names(dat_list))) {
            hOrder = seq(length(dat_list))
        } else {
            hOrder = names(dat_list)
        }
    }

    if (length(hOrder) != length(dat_list)) {
        stop("Error: length of hOrder must be identical to length of dat_list")
    }

    if (is.numeric(hOrder)) {
        hOrder = paste("B", hOrder, sep = "")
        names(dat_list) = hOrder
    }


    if (is.null(newAssay)) {
        if (balanced) {
            newTag = "balanced"
        } else {
            newTag = "concatenate"
        }
        newAssay = paste(assay, newTag, sep = "_")
    }
    # Creat a new assay
    dat_list = lapply(dat_list, function(dat) {

        SummarizedExperiment::assay(dat, newAssay) = SummarizedExperiment::assay(dat, assay)
        dat
    })

    if (balanced) {
        .balanced_tree(dat_list, assay = newAssay, k=k, rep = rep, negCtl = negCtl, hOrder = hOrder)
    } else {
        .concatenating_tree(dat_list, assay = newAssay, k=k, rep = rep, negCtl = negCtl, hOrder = hOrder)
    }
}





#' @importFrom SummarizedExperiment cbind
.balanced_tree = function(dat_list, assay, k, rep, negCtl, hOrder) {
    comp = hOrder
    names(comp) = hOrder

    odd = FALSE
    dat_comp = dat_list
    while(length(comp) > 1) {
        if (length(comp) %% 2 == 1) {
            odd = TRUE
        }
        comp1 = names(comp)[c(TRUE, FALSE)]
        comp2 = names(comp)[c(FALSE, TRUE)]

        comp = tryCatch({
            mapply(c, comp1, comp2, SIMPLIFY = FALSE)
        }, warning = function(w) {
            suppressWarnings({
                mapply(c, comp1, comp2, SIMPLIFY = FALSE)
            })
        })

        if (odd) {
            comp[[tail(comp1, 1)]] = comp[[tail(comp1, 1)]][1]
            odd = FALSE
        }

        dat_comp = mapply(function(comparisons) {
            dat_comp[comparisons]
        }, comp, USE.NAMES = TRUE, SIMPLIFY = FALSE)

        dat_comp = lapply(dat_comp, function(dat_list_sub) {
            if (length(dat_list_sub) == 1) {
                return(dat_list_sub[[1]])
            }
            # dat_concat = combine_se(
            #     dat_list_sub,
            #     assay_name = assay,
            #     data_name = "combined"
            # )
            batch = names(dat_list_sub)
            prev_sample = dat_list_sub[[batch[1]]]$sample_names
            next_sample = dat_list_sub[[batch[2]]]$sample_names
            # next_sample = gsub("\\*{3}", "", next_sample[which(gsub("\\*{3}", "", next_sample) %in% prev_sample)])
            next_sample[which(gsub("\\*{3}", "", next_sample) %in% prev_sample)] = gsub("\\*{3}", "", next_sample[which(gsub("\\*{3}", "", next_sample) %in% prev_sample)])

            replication = c(prev_sample, next_sample)
            replication[which(grepl("pool", replication, ignore.case = TRUE))] = paste("Pool ", seq(length(which(grepl("pool", replication, ignore.case = TRUE)))))
            replicate_mat = replicate.matrix(replication)

            dat_concat = do.call(SummarizedExperiment::cbind, dat_list_sub)
            # batch = names(dat_list_sub)
            if ((is.null(dat_list_sub[[batch[1]]]$batch_info)) |
                    (is.null(dat_list_sub[[batch[1]]]$batch_info))) {
                dat_concat[["batch_info"]] = c(
                    rep(batch[1], ncol(dat_list_sub[[batch[1]]])),
                    rep(batch[2], ncol(dat_list_sub[[batch[2]]]))
                )
            }
            dat = dat_concat
            dat = getRUV(dat, assay = assay, k = k, rep = rep, negCtl = negCtl, replicate_mat = replicate_mat)
            dat
        })
    }
    dat = dat_comp[[1]]
}

.concatenating_tree = function(dat_list, assay, k, rep, negCtl, hOrder) {
    dat = dat_list[[hOrder[1]]]
    dat$batch_info = rep(hOrder[1], ncol(dat))
    # alpha_list = list()
    for (batch in hOrder[-1]) {
        next_batch = dat_list[[batch]]
        next_batch$batch_info = rep(batch, ncol(next_batch))

        # dat_concat = combine_se(
        #     list(
        #         prev = dat,
        #         next_batch = next_batch
        #     ),
        #     assay_name = new_assay,
        #     data_name = "combined"
        # )

        prev_sample = dat$sample_name
        next_sample = next_batch$sample_name
        # next_sample = gsub("\\*{3}", "", next_sample[which(gsub("\\*{3}", "", next_sample) %in% prev_sample)])
        next_sample[which(gsub("\\*{3}", "", next_sample) %in% prev_sample)] = gsub("\\*{3}", "", next_sample[which(gsub("\\*{3}", "", next_sample) %in% prev_sample)])

        replication = c(prev_sample, next_sample)
        replication[which(grepl("pool", replication, ignore.case = TRUE))] = paste("Pool ", seq(length(which(grepl("pool", replication, ignore.case = TRUE)))))
        replicate_mat = replicate.matrix(replication)

        dat_concat = SummarizedExperiment::cbind(dat, next_batch)

        dat = dat_concat

        # if (ctrl == "batch") {
        #     replication = dat_concat$sample_name
        #     next_batch$sample_name[next_batch$batch_replicate] = gsub("(\\*{3})|(\\^)", "", next_batch$sample_name[next_batch$batch_replicate])
        #     replication[tail(1:length(replication), n =ncol(next_batch))] = next_batch$sample_name
        #     replication[dat_concat$pool] = paste("Pool ", 1:sum(dat_concat$pool), sep = "")
        # } else if (ctrl == "pool") {
        #     dat$sample_name[dat$pool] = gsub("(Pool\\ *[0-9]+)|(HPP\\ *[0-9]+)", "Pool", dat$sample_name[dat$pool])
        #     replication = dat$sample_name
        # }
        # replication_matrix = replicate.matrix(replication)

        # dat = get_ruv(dat, assay_name = new_assay, k = k, ctrl = ctrl,
        #     new_assay = new_assay, replicate_mat = replication_matrix,
        #     neg_mets = neg_mets, alpha = alpha)
        dat = getRUV(dat, assay = assay, k = k, rep = rep, negCtl = negCtl, replicate_mat = replicate_mat)
    }
    dat


}
