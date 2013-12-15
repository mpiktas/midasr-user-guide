graph <- function(x,pv0hac=NULL,titlestr=TRUE) {
    cfur <- coef(x$unrestricted)
    cfur <- cfur[grep("mls",names(cfur))]
    cfre <- weight_coef(x)
    k <- length(cfur)
    sdval <- sqrt(diag(sandwich(x$unrestricted)))
    sdval <- sdval[grep("mls",names(sdval))]
    if(is.null(pv0hac))pv0hac <- hAhr.test(x)$p.value
    plot(c(0:(k - 1)), c(cfur), col = "black", ylab = "Beta coefficients", xlab = "h")

    if(titlestr) title(main = sprintf("p-val.(hAh_HAC) < %.2f", max(pv0hac, 0.01)), cex.main = 1, font.main = 4, col.main = "black")
    
    points(c(0:(k - 1)), cfre[1:k], type = "l", col = "blue")
    points(c(0:(k - 1)), cfur[1:k] + 2 * sdval[1:k], type = "l", col = "red", lty = 2)
    points(c(0:(k - 1)), cfur[1:k] - 2 * sdval[1:k], type = "l", col = "red", lty = 2)
}

harstep <- 
function(p,d,m) {
    if(d!=20) stop("HAR(3)-RV process requires 20 lags")
    out <- rep(0,20)
    out[1] <- p[1]+p[2]/5+p[3]/20
    out[2:5] <- p[2]/5+p[3]/20
    out[6:20] <- p[3]/20
    out
}
