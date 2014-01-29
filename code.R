harstep <- 
function(p,d,m) {
    if(d!=20) stop("HAR(3)-RV process requires 20 lags")
    out <- rep(0,20)
    out[1] <- p[1]+p[2]/5+p[3]/20
    out[2:5] <- p[2]/5+p[3]/20
    out[6:20] <- p[3]/20
    out
}
