pairwise_simpleLM <- function (dat) {
    ## matrix and its dimension (n: numbeta.ser of data; p: numbeta.ser of variables)
    dat <- as.matrix(dat)
    n <- nrow(dat)
    p <- ncol(dat)
    ## variable summary: mean, (unscaled) covariance and (unscaled) variance
    m <- colMeans(dat)
    V <- crossprod(dat) - tcrossprod(m * sqrt(n))
    d <- diag(V)
    ## R-squared (explained variance) and its complement
    R2 <- (V ^ 2) * tcrossprod(1 / d)
    R2_complement <- 1 - R2
    R2_complement[seq.int(from = 1, by = p + 1, length = p)] <- 0
    ## slope and intercept
    beta <- V * rep(1 / d, each = p)
    alpha <- m - beta * rep(m, each = p)
    ## residual sum of squares and standard error
    RSS <- R2_complement * d
    sig <- sqrt(RSS * (1 / (n - 2)))
    ## statistics for slope
    beta.se <- sig * rep(1 / sqrt(d), each = p)
    beta.tv <- beta / beta.se
    beta.pv <- 2 * pt(abs(beta.tv), n - 2, lower.tail = FALSE)
    ## F-statistic and p-value
    F.fv <- (n - 2) * R2 / R2_complement
    F.pv <- pf(F.fv, 1, n - 2, lower.tail = FALSE)
    ## export
    data.frame(LHS = rep(colnames(dat), times = p),
    RHS = rep(colnames(dat), each = p),
    alpha = c(alpha),
    beta = c(beta),
    beta.se = c(beta.se),
    beta.tv = c(beta.tv),
    beta.pv = c(beta.pv),
    sig = c(sig),
    R2 = c(R2),
    F.fv = c(F.fv),
    F.pv = c(F.pv),
    stringsAsFactors = FALSE)
}

general_paired_simpleLM <- function (dat_LHS, dat_RHS) {
    ## matrix and its dimension (n: numbeta.ser of data; p: numbeta.ser of variables)
    dat_LHS <- as.matrix(dat_LHS)
    dat_RHS <- as.matrix(dat_RHS)
    if (nrow(dat_LHS) != nrow(dat_RHS)) stop("'dat_LHS' and 'dat_RHS' don't have same number of rows!")
    n <- nrow(dat_LHS)
    pl <- ncol(dat_LHS)
    pr <- ncol(dat_RHS)
    ## variable summary: mean, (unscaled) covariance and (unscaled) variance
    ml <- colMeans(dat_LHS)
    mr <- colMeans(dat_RHS)
    vl <- colSums(dat_LHS ^ 2) - ml * ml * n
    vr <- colSums(dat_RHS ^ 2) - mr * mr * n
    ##V <- crossprod(dat - rep(m, each = n))  ## cov(u, v) = E[(u - E[u])(v - E[v])]
    V <- crossprod(dat_LHS, dat_RHS) - tcrossprod(ml * sqrt(n), mr * sqrt(n))  ## cov(u, v) = E[uv] - E{u]E[v]
    ## R-squared (explained variance) and its complement
    R2 <- (V ^ 2) * tcrossprod(1 / vl, 1 / vr)
    R2_complement <- 1 - R2
    ## slope and intercept
    beta <- V * rep(1 / vr, each = pl)
    alpha <- ml - beta * rep(mr, each = pl)
    ## residual sum of squares and standard error
    RSS <- R2_complement * vl
    sig <- sqrt(RSS * (1 / (n - 2)))
    ## statistics for slope
    beta.se <- sig * rep(1 / sqrt(vr), each = pl)
    beta.tv <- beta / beta.se
    beta.pv <- 2 * pt(abs(beta.tv), n - 2, lower.tail = FALSE)
    ## F-statistic and p-value
    F.fv <- (n - 2) * R2 / R2_complement
    F.pv <- pf(F.fv, 1, n - 2, lower.tail = FALSE)
    ## export
    data.frame(LHS = rep(colnames(dat_LHS), times = pr),
    RHS = rep(colnames(dat_RHS), each = pl),
    alpha = c(alpha),
    beta = c(beta),
    beta.se = c(beta.se),
    beta.tv = c(beta.tv),
    beta.pv = c(beta.pv),
    sig = c(sig),
    R2 = c(R2),
    F.fv = c(F.fv),
    F.pv = c(F.pv),
    stringsAsFactors = FALSE)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    
    test <- cor.test(x,y)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " "))
    
    text(0.5, 0.5, txt, cex = cex * r)
    text(.8, .8, Signif, cex=cex, col=2)
}

reg <- function(x, y, col) abline(lm(y~x), col=col)

panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"),
cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) reg(x[ok], y[ok], col.smooth)
}
