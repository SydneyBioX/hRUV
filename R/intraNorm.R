#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom dplyr select %>% all_of arrange
#'
smoothAdj = function(dat, assay, method = c("rlm", "loess"), pCtl = NULL) {
    run_order = seq(ncol(dat))

    # If ctl is not provided, all samples will be used for fitting the smoother
    if (is.null(pCtl)) {
        pCtl = run_order
    }
    dat_assay = SummarizedExperiment::assay(dat, assay)

    ctl_dat = dat_assay %>%
        as.data.frame() %>%
        dplyr::select(all_of(pCtl))
    ctl_meta = data.frame(
        run = run_order[pCtl],
        name = colnames(ctl_dat)
    )

    if (identical(run_order, pCtl)) {
        sample_dat = dat_assay %>%
            as.data.frame()
        sample_meta = data.frame(
            run = run_order,
            name = colnames(sample_dat)
        )
    } else {
        sample_dat = dat_assay %>%
            as.data.frame() %>%
            dplyr::select(-all_of(pCtl))
        sample_meta = data.frame(
            run = run_order[-pCtl],
            name = colnames(sample_dat)
        )
    }

    dat_norm = .fit_model(sample_dat, ctl_dat, sample_meta, ctl_meta, method = method)

    # Form to an original data frame form
    sample_dat = dat_norm$norm_mat
    ctl_dat = dat_norm$norm_mat_control

    sample_run = rbind(ctl_meta, sample_meta) %>%
        arrange(run)

    dat_result = cbind(sample_dat, ctl_dat) %>%
        as.data.frame()
    dat_result = dat_result %>%
        dplyr::select(sample_run$name)

    SummarizedExperiment::assay(dat, method) = dat_result

    dat
}

#' @importFrom MASS rlm
.fit_model = function(sample, pCtl, sample_run, ctl_run, method) {
    for (i in seq_len(nrow(sample))) {
        norm_dat = data.frame(
            x = as.numeric(as.character(sample_run$run)),
            y = as.numeric(as.character(sample[i,]))
        )
        norm_control_dat = data.frame(
            x = as.numeric(as.character(ctl_run$run)),
            y = as.numeric(as.character(pCtl[i,]))
        )

        if (method == "rlm") {
            # Fit rlm
            fit = MASS::rlm(y~x, data = norm_control_dat, maxit=100)
        } else if (method == "loess") {
            # Fit loess
            fit = loess(y~x, degree = 2, family = "symmetric", surface = "direct", data = norm_control_dat)
        }
        median_control = median(norm_control_dat$y, na.rm = T)

        norm_control_run = data.frame(
            x = norm_control_dat$x
        )
        norm_run = data.frame(
            x = norm_dat$x
        )
        y_run_hat = predict(fit, newdata = norm_run)
        y_control_run_hat = predict(fit, newdata = norm_control_run)

        adj = median_control - y_run_hat
        adj_control = median_control - y_control_run_hat

        pCtl[i,] = norm_control_dat$y + adj_control
        sample[i,] = norm_dat$y + adj
    }

    list(
        norm_mat = sample,
        norm_mat_control = pCtl
    )
}


#' @export
intraNorm = function(dat_list, assay,
    method = c("rlm", "loess", "rlmShort", "loessShort"),
    pCtlName = NULL, rep = NULL, k, negCtl = NULL) {

    dat_list = lapply(dat_list, function(dat) {
        if (grepl("rlm", method)) {
            dat = smoothAdj(dat, assay, method = "rlm", pCtl = which(dat[[pCtlName]]))
        } else if (grepl("loess", method)) {
            dat = smoothAdj(dat, assay, method = "loess", pCtl = which(dat[[pCtlName]]))
        }

        if (grepl("Short", method)) {
            method1 = gsub("Short", "", method)
            dat = getRUV(dat, assay = method1, k = k, rep, negCtl = negCtl, newAssay = method)
        }
        dat
    })
}

