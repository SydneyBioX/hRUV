# This is an implementation derived from ruv package
RUVIII <-
    function (Y, M, ctl, k = NULL, eta = NULL, include.intercept = TRUE,
        average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE)
    {
        if (is.data.frame(Y))
            Y = data.matrix(Y)
        m = nrow(Y)
        n = ncol(Y)
        M = replicate.matrix(M)
        ctl = tological(ctl, n)
        if (inputcheck) {
            if (sum(is.na(Y)) > 0)
                warning("Y contains missing values.  This is not supported.")
            if (sum(Y == Inf, na.rm = TRUE) + sum(Y == -Inf, na.rm = TRUE) >
                    0)
                warning("Y contains infinities.  This is not supported.")
        }
        Y = RUV1(Y, eta, ctl, include.intercept = include.intercept)
        if (ncol(M) >= m)
            newY = Y
        else if (is.null(k)) {
            ycyctinv = solve(Y[, ctl] %*% t(Y[, ctl]))
            newY = (M %*% solve(t(M) %*% ycyctinv %*% M) %*% (t(M) %*%
                    ycyctinv)) %*% Y
            fullalpha = NULL
        }
        else if (k == 0) {
            newY = Y
            fullalpha = NULL
        }
        else {
            if (is.null(fullalpha)) {
                Y0 = residop(Y, M)
                fullalpha = t(svd(Y0 %*% t(Y0))$u[, 1:min(m - ncol(M),
                    sum(ctl)), drop = FALSE]) %*% Y
            }
            alpha = fullalpha[1:min(k, nrow(fullalpha)), , drop = FALSE]
            ac = alpha[, ctl, drop = FALSE]
            W = Y[, ctl] %*% t(ac) %*% solve(ac %*% t(ac))
            newY = Y - W %*% alpha
        }
        if (average)
            newY = ((1/apply(M, 2, sum)) * t(M)) %*% newY
        if (!return.info)
            return(newY)
        else return(list(newY = newY, M = M, fullalpha = fullalpha))
    }


RUV1 <-
    function (Y, eta, ctl, include.intercept = TRUE)
    {
        if (is.null(eta))
            return(Y)
        m = nrow(Y)
        n = ncol(Y)
        ctl = tological(ctl, n)
        if (is.numeric(eta))
            if (length(eta) == 1)
                eta = matrix(1, n, 1)
        if (is.matrix(eta))
            if (nrow(eta) != n)
                eta = t(eta)
        eta = design.matrix(eta, name = "eta", include.intercept = include.intercept)
        eta = t(eta)
        Yc = Y[, ctl, drop = FALSE]
        etac = eta[, ctl, drop = FALSE]
        if (sum(is.na(Yc)) == 0)
            return(Y - Yc %*% t(etac) %*% solve(etac %*% t(etac)) %*%
                    eta)
        else {
            for (i in 1:m) {
                keep = !is.na(Yc[i, ])
                Yci = Yc[i, keep, drop = FALSE]
                etaci = etac[, keep, drop = FALSE]
                Y[i, ] = Y[i, ] - Yci %*% t(etaci) %*% solve(etaci %*%
                        t(etaci)) %*% eta
            }
            return(Y)
        }
    }


design.matrix <-
    function (a, name = "X", remove.collinear = TRUE, include.intercept = TRUE)
    {
        if (is.vector(a))
            a = matrix(a)
        if (is.matrix(a)) {
            if (is.null(colnames(a)))
                colnames(a) = paste(name, 1:ncol(a), sep = "")
            for (i in 1:ncol(a)) {
                if (is.na(colnames(a)[i]))
                    colnames(a)[i] = paste(name, i, sep = "")
                if (colnames(a)[i] == "")
                    colnames(a)[i] = paste(name, i, sep = "")
            }
        }
        if (is.numeric(a))
            A = a
        else {
            if (is.factor(a)) {
                a = data.frame(a)
                names(a)[1] = paste(name, 1, sep = "")
            }
            a = data.frame(a)
            varnames = colnames(a)
            for (i in 1:length(varnames)) if (varnames[i] == paste("X",
                i, sep = ""))
                varnames[i] = paste(name, i, sep = "")
            A = design.column(a[, 1], name = varnames[1])
            if (ncol(a) > 1)
                for (i in 2:ncol(a)) A = cbind(A, design.column(a[,
                    i], name = varnames[i]))
        }
        if (remove.collinear)
            if (ncol(A) > 1) {
                if (ncol(A) > nrow(A))
                    A = A[, 1:nrow(A)]
                A0 = scale(A, center = FALSE, scale = TRUE)
                d = svd(A)$d
                if (d[1]/d[length(d)] > 10^9) {
                    warning("Collinearity detected.  Removing some columns.")
                    toremove = NULL
                    for (i in 2:ncol(A0)) {
                        A1 = A0[, 1:(i - 1)]
                        if (!is.null(toremove))
                            A1 = A1[, -toremove]
                        if (mean(residop(A0[, i, drop = FALSE], A1)^2) <
                                10^(-9))
                            toremove = c(toremove, i)
                    }
                    A = A[, -toremove, drop = FALSE]
                }
            }
        if (include.intercept) {
            if (sum(residop(matrix(1, nrow(A), 1), A)^2) > 10^(-8)) {
                A = cbind(1, A)
                colnames(A)[1] = paste(name, "0", sep = "")
            }
        }
        return(A)
    }


tological <-
    function (ctl, n)
    {
        ctl2 = rep(FALSE, n)
        ctl2[ctl] = TRUE
        return(ctl2)
    }


residop <-
    function (A, B)
    {
        return(A - B %*% solve(t(B) %*% B) %*% t(B) %*% A)
    }


replicate.matrix <-
    function (a, burst = NULL, return.factor = FALSE, name = "M",
        sep = "_", burstsep = "_")
    {
        if (is.matrix(a) & is.numeric(a)) {
            m = nrow(a)
            b = abs(a) > 10^(-8)
            if ((sum(abs(apply(a, 1, sum) - rep(1, m))) < 10^(-8)) &
                    (sum(abs(apply(b, 1, sum) - rep(1, m))) < 10^(-8))) {
                M = a
                if (is.null(colnames(M)))
                    colnames(M) = paste(name, 1:ncol(M), sep = "")
                for (i in 1:ncol(M)) {
                    if (is.na(colnames(M)[i]))
                        colnames(M)[i] = paste(name, i, sep = "")
                    if (colnames(M)[i] == "")
                        colnames(M)[i] = paste(name, i, sep = "")
                }
                colnames(M) = make.names(colnames(M), unique = TRUE)
                if (is.null(burst) & (!return.factor))
                    return(M)
                a = rep("", m)
                for (i in 1:ncol(M)) a[M[, i] > 0.5] = colnames(M)[i]
                a = as.factor(a)
            }
        }
        a = data.frame(a)
        A = matrix("", nrow(a), ncol(a))
        for (i in 1:ncol(a)) A[, i] = as.character(a[, i])
        a = make.names(apply(A, 1, paste, collapse = sep))
        if (!is.null(burst)) {
            for (b in burst) {
                bn = sum(a == b)
                if (bn == 0)
                    warning(paste("Unable to burst factor level: ",
                        b, "  -- not found.", sep = ""))
                a[a == b] = paste(b, 1:bn, sep = burstsep)
            }
        }
        a = as.factor(a)
        if (return.factor)
            return(a)
        M = model.matrix(~a - 1)
        colnames(M) = levels(a)
        return(M)
    }
