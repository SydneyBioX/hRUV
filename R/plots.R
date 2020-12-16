#' @export
plotPCA = function(dat, assay, colour, ...) {
    mat = SummarizedExperiment::assay(dat, assay)
    mat = prcomp(t(mat), ...)
    pc_var = mat$sdev^2
    pc_var = signif((pc_var[1:2] / sum(pc_var)) * 100, 4)
    plt_dat = data.frame(
        PC1 = mat$x[,1],
        PC2 = mat$x[,2],
        samples = rownames(mat$x)
    )
    plt_dat$colour = colour
    p = ggplot2::ggplot(plt_dat, aes(x = PC1, y = PC2, color = colour)) +
        geom_point() +
        ggtitle(paste("PCA, ", assay, sep = "")) +
        labs(x = paste("PC1 (", pc_var[1], "% Var)", sep = ""), y = paste("PC2 (", pc_var[2], "% Var)", sep = ""))
    p
}


#' @export
plotRLE = function(dat, assay, group = NULL) {
    # Assumes the assay to use is logged and no missing values
    if (is.null(group)) {
        group = rep("samples", ncol(dat))
    }
    group = factor(group)
    cur_dat = SummarizedExperiment::assay(dat, assay)
    rla_dat = NULL
    for (g in levels(group)) {
        sub_dat = cur_dat[,group == g, drop = FALSE]
        sub_med = apply(sub_dat, 1, median) # Row median
        sub_rla = sweep(sub_dat, 1, sub_med, "-")
        sub_rla = as.matrix(sub_rla)
        rla_dat = cbind(rla_dat, sub_rla)
    }
    rla_dat = rla_dat[, order(group)]
    sample_named_order = colnames(rla_dat)

    rla_dat = rla_dat %>%
        as.data.frame()
    rla_dat$Metabolite = rownames(rla_dat)
    rla_dat_long = rla_dat %>%
        tidyr::pivot_longer(-Metabolite, names_to = "Sample", values_to = "value")
    rla_dat_long$Colour = rep(group[order(group)], nrow(rla_dat))
    rla_dat_long$Sample = factor(rla_dat_long$Sample, levels = sample_named_order)

    p = ggplot2::ggplot(rla_dat_long, aes(x = Sample, y = value, fill = Colour)) +
        geom_boxplot() +
        geom_hline(yintercept = 0, color = "red")
    p
}


#' @export
plotRun = function(dat, assay, colour) {
    mat = wide_to_long(dat, assay)
    metabolites = unique(mat$Metabolite)
    sapply(metabolites, function(j) {
        tmp = mat %>%
            dplyr::filter(grepl(paste("^", j, "$", sep = ""), Metabolite))
        p = ggplot2::ggplot(tmp, ggplot2::aes(x = run, y = value, color = colour)) +
            ggplot2::geom_point() +
            ggplot2::ggtitle(paste(j, ", ", comp))
        p
    }, USE.NAMES = TRUE, simplify = FALSE)
}

metabolite_plot = function(dat, comp, highlight, pkg, colour_by, colour_vec) {
    mat = wide_to_long(dat, comp, highlight)
    metabolites = unique(mat$Metabolite)
    mapply(function(j) {
        tmp = mat %>%
            dplyr::filter(grepl(paste("^", j, "$", sep = ""), Metabolite))
        if (colour_by == "batch") {
            p = ggplot2::ggplot(tmp, ggplot2::aes(x = run, y = value, color = dat$batch_info)) +
                ggplot2::geom_point() +
                ggplot2::ggtitle(paste(j, ", ", comp))
        } else if ((colour_by == "custom") && (!is.null(colour_vec))) {
            p = ggplot2::ggplot(tmp, ggplot2::aes(x = run, y = value, color = colour_vec)) +
                ggplot2::geom_point() +
                ggplot2::ggtitle(paste(j, ", ", comp))
        } else {
            p = ggplot2::ggplot(tmp, ggplot2::aes(x = run, y = value)) +
                ggplot2::geom_point() +
                ggplot2::ggtitle(paste(j, ", ", comp))
        }
        if (highlight != "none") {
            p = p +
                ggplot2::geom_point(data = tmp[which(!grepl("Other", tmp$replicate)),], size = 3, shape = 1, color = "black")
        }
        p = p +
            ggplot2::ggtitle(paste(j, ", ", comp))
        # ggpubr::ggarrange(p1, p2, nrow = 2)
    }, metabolites, USE.NAMES = TRUE, SIMPLIFY = FALSE)
}


wide_to_long = function(i, assay, colour) {
    mat = SummarizedExperiment::assay(i, assay) %>%
        as.data.frame()
    mat$Metabolite = rownames(mat)
    mat = mat %>%
        tidyr::pivot_longer(-Metabolite, names_to = "Sample", values_to = "value")
    mat$run = rep(seq_len(ncol(i)), nrow(i))
    mat$colour = rep(colour, nrow(i))
    mat
}
