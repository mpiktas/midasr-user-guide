# Preliminaries -----------------------------------------------------------
# rm(list=ls) #Uncomment to clean up the workspace
library(midasr)  # loads the midasr package


# Example of simulated MIDAS regression  ---------------------------------------------------

## Sets a seed for RNG ###
set.seed(1001)  # just for comparability of results
## Number of low-frequency observations
n <- 250
## Linear trend and higher-frequency explanatory variables (e.g. quarterly and monthly)
trend <- c(1:n)
x <- rnorm(4 * n)
z <- rnorm(12 * n)
## Exponential Almon polynomial constraint-consistent coefficients
fn_x <- nealmon(p = c(1, -0.5), d = 8)
fn_z <- nealmon(p = c(2, 0.5, -0.1), d = 17)
## Simulated low-frequency series (e.g. yearly)
y <- 2 + 0.1 * trend + mls(x, 0:7, 4) %*% fn_x + mls(z, 0:16, 12) %*% fn_z + rnorm(n)
## Figure 1 (coefficients)
fn_xNA <- c(fn_x, rep(NA, length(fn_z) - length(fn_x)))
plot(fn_z, col = "red")
points(fn_xNA)


# Examples of MIDAS regression specification in midasr --------------------

## OLS using lm
eq_u <- lm(y ~ trend + mls(x, k = 0:7, m = 4) + mls(z, k = 0:16, m = 12))
eq_u <- midas_u(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12))
summary(eq_u)

## NLS using midas_r
eq_r <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, 
    -0.5), z = c(2, 0.5, -0.1)))
summary(eq_r)
deriv_tests(eq_r, tol = 1e-06)
coef(eq_r)
coef(eq_r, midas = TRUE)
amweights(p = c(1, -0.5), d = 8, m = 4, weight = nealmon, type = "C")
nealmon(p = c(1, -0.5), d = 4)
## NLS using midas_r with aggregates
eq_r <- midas_r(y ~ trend + mls(x, 0:7, 4, mmweights, nealmon, "C") + mls(z, 0:16, 12, nealmon), 
    start = list(x = c(1, -0.5), z = c(2, 0.5, -0.1)))
summary(eq_r)
## Table 3-related NLS variations using midas_r (trend dropped in the table and can be
## omitted here)
fn <- gompertzp
eq_r1 <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, fn), start = list(x = c(1, 
    -0.5), z = c(1, 0.5, 0.1)))
summary(eq_r1)
eq_r2 <- midas_r(y ~ trend + mls(x, 0:7, 4) + mls(z, 0:16, 12, nealmon), start = list(z = c(1, 
    -0.5)))
summary(eq_r2)
eq_r3 <- midas_r(y ~ trend + mls(y, 1:2, 1) + mls(x, 0:7, 4, nealmon), start = list(x = c(1, 
    -0.5)))
summary(eq_r3)
eq_r4 <- midas_r(y ~ trend + mls(y, 1:2, 1, "*") + mls(x, 0:7, 4, nealmon), start = list(x = c(1, 
    -0.5)))
summary(eq_r4)
eq_r5 <- midas_r(y ~ trend + mls(y, 1:4, 1, nealmon) + mls(x, 0:7, 4, nealmon), start = list(y = c(1, 
    -0.5), x = c(1, -0.5)))
summary(eq_r5)
eq_r6 <- midas_r(y ~ trend + mls(x, 0:7, 4, amweights, nealmon, "A"), start = list(x = c(1, 
    1, 1, -0.5)))
summary(eq_r6)
eq_r7 <- midas_r(y ~ trend + mls(x, 0:7, 4, amweights, nealmon, "B"), start = list(x = c(1, 
    1, -0.5)))
summary(eq_r7)
eq_r8 <- midas_r(y ~ trend + mls(x, 0:7, 4, amweights, nealmon, "C"), start = list(x = c(1, 
    -0.5)))
summary(eq_r8)
fn <- function(p, d) {
    p[1] * c(1:d)^p[2]
}
eq_r9 <- midas_r(y ~ trend + mls(x, 0:101, 4, fn), start = list(x = rep(0, 2)))
summary(eq_r9)


# Testing the adequacy of MIDAS regression --------------------------------

## DGP-consistent specification
eq_r <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:16, 12, nealmon), start = list(x = c(1, 
    -0.5), z = c(2, 0.5, -0.1)))
summary(eq_r)
hAh_test(eq_r)
hAhr_test(eq_r)
## Mis-specification of constraint on z coefficients
eq_rb <- midas_r(y ~ trend + mls(x, 0:7, 4, nealmon) + mls(z, 0:12, 12, nealmon), start = list(x = c(1, 
    -0.5), z = c(2, -0.1)))
hAh_test(eq_rb)
hAhr_test(eq_rb)
summary(eq_rb)


# Model selection ---------------------------------------------------------

## Potential sets of models
set_x <- expand_weights_lags(weights = c("nealmon", "almonp"), from = 0, to = c(5, 10), m = 1, 
    start = list(nealmon = rep(1, 2), almonp = rep(1, 3)))
set_z <- expand_weights_lags(c("nealmon", "nealmon"), 0, c(10, 20), 1, start = list(nealmon = rep(1, 
    2), nealmon = rep(1, 3)))
expand_weights_lags(weights = c("nealmon", "nbeta"), from = 1, to = c(2, 3), m = 1, start = list(nealmon = rep(0, 
    2), nbeta = rep(0.5, 3)))
## Estimation and selection of models
eqs.ic <- midas_r_ic_table(y ~ trend + mls(x, 0, m = 4) + fmls(z, 0, m = 12), table = list(z = set_z, 
    x = set_x), start = c(`(Intercept)` = 0, trend = 0))
mod <- modsel(eqs.ic, IC = "AIC", type = "restricted")


# Forecast combinations ---------------------------------------------------

## With model selection and horizon 1:3
cbfc <- select_and_forecast(y ~ trend + mls(x, 0, m = 4) + mls(z, 0, m = 12), from = list(x = c(4, 
    8, 12), z = c(12, 24, 36)), to = list(x = rbind(c(14, 19), c(18, 23), c(22, 27)), z = rbind(c(22, 
    27), c(34, 39), c(46, 51))), insample = 1:200, outsample = 201:250, weights = list(x = c("nealmon", 
    "almonp"), z = c("nealmon", "almonp")), wstart = list(nealmon = rep(1, 3), almonp = rep(1, 
    3)), IC = "AIC", seltype = "restricted", ftype = "fixed", measures = c("MSE", "MAPE", "MASE"), 
    fweights = c("EW", "BICW", "MSFE", "DMSFE"))
cbfc$accuracy$individual
cbfc$accuracy$average

## With a given specification one period ahead
nealmon2 <- nealmon
nealmon3 <- nealmon
cbfc1 <- select_and_forecast(y ~ trend + mls(x, 0, 4) + mls(z, 0, 12), from = list(x = c(4), 
    z = c(12)), to = list(x = rbind(c(14, 14)), z = rbind(c(22, 22))), insample = 1:200, outsample = 201:250, 
    weights = list(x = c("nealmon3"), z = c("nealmon2")), wstart = list(nealmon3 = c(10, 1, 
        -0.1), nealmon2 = c(2, -0.1)), IC = "AIC", seltype = "restricted", ftype = "fixed", 
    measures = c("MSE", "MAPE", "MASE"), fweights = c("EW", "BICW", "MSFE", "DMSFE"))
cbfc1$accuracy$individual
cbfc1$accuracy$average


# Manual model selection --------------------------------------------------

## First split data into in-sample and out-of-sample data

datasplit <- split_data(list(y = y, x = x, z = z, trend = trend), insample = 1:200, outsample = 201:250)

## Fit two models

mod1 <- midas_r(y ~ trend + mls(x, 4:14, 4, nealmon3) + mls(z, 12:22, 12, nealmon2), data = datasplit$indata, 
    start = list(x = c(10, 1, -0.1), z = c(2, -0.1)))

mod2 <- midas_r(y ~ trend + mls(x, 4:20, 4, nealmon3) + mls(z, 12:25, 12, nealmon2), data = datasplit$indata, 
    start = list(x = c(10, 1, -0.1), z = c(2, -0.1)))

## Calculate average forecasts

avgf <- average_forecast(list(mod1, mod2), data = list(y = y, x = x, z = z, trend = trend), 
    insample = 1:200, outsample = 201:250, type = "fixed", measures = c("MSE", "MAPE", "MASE"), 
    fweights = c("EW", "BICW", "MSFE", "DMSFE"))

avgf$accuracy
avgf$forecast
avgf$avgforecast

## Produce rolling forecasts, where for each forecast models are reestimated using rolling
## window
avgrollf <- average_forecast(list(mod1, mod2), data = list(y = y, x = x, z = z, trend = trend), 
    insample = 1:200, outsample = 201:210, type = "rolling", measures = c("MSE", "MAPE", "MASE"), 
    fweights = c("EW", "BICW", "MSFE", "DMSFE"))

avgrollf$accuracy

## Produce recursive forecasts where for each forecast models are reestimated by recursively
## increasing estimation sample
avgrecf <- average_forecast(list(mod1, mod2), data = list(y = y, x = x, z = z, trend = trend), 
    insample = 1:200, outsample = 201:210, type = "recursive", measures = c("MSE", "MAPE", "MASE"), 
    fweights = c("EW", "BICW", "MSFE", "DMSFE"))

avgrecf$accuracy


# Inspect objects produced by midasr --------------------------------------

objects(eq_r)
objects(cbfc)
## Accessed e.g. by
eq_r$opt
cbfc$bestlist
############ Info ### A specific function
`?`(select_and_forecast)
## On midasr package
`?`(`?`(midasr))



# MIDAS Matlab toolbox example ----------------------------------------------------------------------

##  Get the data
library(quantmod)
gdp <- getSymbols("GDP", src = "FRED", auto.assign = FALSE)
payems <- getSymbols("PAYEMS", src = "FRED", auto.assign = FALSE)

## Convert to ts with the exact sample used in MIDAS Matlab toolbox
y <- window(ts(gdp, start = c(1947, 1), frequency = 4), end = c(2011, 2))
x <- window(ts(payems, start = c(1939, 1), frequency = 12), end = c(2011, 7))

## Calculate the log differences
yg <- log(y/lag(y, -1)) * 100
xg <- log(x/lag(x, -1)) * 100

## Align data
nx <- ts(c(NA, xg, NA, NA), start = start(x), frequency = 12)
ny <- ts(c(rep(NA, 33), yg, NA), start = start(x), frequency = 4)


## Replicate the graph of the data
plot.ts(nx, xlab = "Time", ylab = "Percentages", col = 4, ylim = c(-5, 6))
lines(ny, col = 2)


## Use the same sample as in MIDAS Matlab toolbox
xx <- window(nx, start = c(1985, 1), end = c(2009, 3))
yy <- window(ny, start = c(1985, 1), end = c(2009, 1))

## Estimate the models
beta0 <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbeta), start = list(xx = c(1.7, 1, 5)))

coef(beta0)

## Note the nbetaMT, which is different form nbeta
betan <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3, nbetaMT), start = list(xx = c(2, 1, 5, 
    0)))
coef(betan)

um <- midas_r(yy ~ mls(yy, 1, 1) + mls(xx, 3:11, 3), start = NULL)
coef(um)


## Split the data into in-sample and out-of-sample
fulldata <- list(xx = window(nx, start = c(1985, 1), end = c(2011, 6)), yy = window(ny, start = c(1985, 
    1), end = c(2011, 2)))
insample <- 1:length(yy)
outsample <- (1:length(fulldata$yy))[-insample]

## Calculate the individual forecasts of each of the model and their weighted average
avgf <- average_forecast(list(beta0, betan, um), data = fulldata, insample = insample, outsample = outsample)
sqrt(avgf$accuracy$individual$MSE.out.of.sample) 
