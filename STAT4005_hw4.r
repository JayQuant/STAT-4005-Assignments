## You can paste the code to R to obtain outputs (refer to global environment).

# Import libraries

library(fpp3)


### Q1 

##a)
# convert to tsibble
Y = c(1.33, -0.56, -1.31, -0.37, 0.05, 0.46, 2.00,
	-0.19, -0.25, 1.07, -0.17, 1.14, 0.63, -0.75, 0.15, 0.71, 0.45, -0.14, 0.57, 1.43)
n = seq.int(length(Y))

dt <- tibble(
	index =n,
	value = Y
	)
	
dt <- as_tsibble(dt, index = index)

# fit MA(2), AR(1), ARMA(1,1), and ARIMA(1,1,0) models to Y_t (no constant)

fit <- dt %>% model(
	ma2 = ARIMA(value ~ 0 + pdq(0,0,2)),
	ar1 = ARIMA(value ~ 0 + pdq(1,0,0)),
	arma11 = ARIMA(value ~ 0 + pdq(1,0,1)),
	arima110 = ARIMA(value ~ 0 + pdq(1,1,0))
	)


# fit same models using stats::arima

fitma <- arima(dt$value, order = c(0,0,2), include.mean = F)


# forecast values using built in forecast()

fit %>% forecast(h=5) %>% filter(.model == 'ma2')

# forecast values by manually iterating
# Since p = 0, we let Z_1 = Y_1, and Z_2 = Y_2 - theta1*Z_1

ma_coef = fit %>% select(ma2) %>% coef()
ma_coef = ma_coef$estimate
ma_sigmasq <- fit %>% select(ma2) %>% glance
ma_sigmasq <- ma_sigmasq$sigma2

n <- length(Y)
Z <- rep(0, n)
Z[1] = Y[1]
Z[2] = Y[2] -ma_coef[1]*Y[1]

for (i in 3:n) {
	Z[i] = Y[i] - ma_coef[1]*Z[i-1] -ma_coef[2]*Z[i-2]
}

# Generate k-step ahead forecast for MA(2) model
ma_pred1 <- ma_coef[1]*Z[n] + ma_coef[2]*Z[n-1]
ma_pred2 <- ma_coef[2]*Z[n]

# Generate 95% prediction interval
ma_p1 <- ma_sigmasq
ma_p2 <- ma_sigmasq*(1+ma_coef[1]^2)
ma_pk <- ma_sigmasq*(1+ma_coef[1]^2 + ma_coef[2]^2)

ma_pi1 <- c(ma_pred1 - 1.96*sqrt(ma_sigmasq), ma_pred1 + 1.96*sqrt(ma_sigmasq))
ma_pi2 <- c(ma_pred2 - 1.96*sqrt(ma_sigmasq*(1+ma_coef[1]^2)),ma_pred2 + 1.96*sqrt(ma_sigmasq*(1+ma_coef[1]^2)))
ma_pik <- c(-1.96*sqrt(ma_sigmasq*(1+ma_coef[1]^2+ma_coef[2]^2)), 1.96*sqrt(ma_sigmasq*(1+ma_coef[1]^2+ma_coef[2]^2)))

##b) find partial autocorrelations using first principle


ma_pacf1 = (ma_coef[1] + ma_coef[1] * ma_coef[2]) / (ma_coef[1]^2 + ma_coef[2]^2 +1)

ma_pacf2 = (acf(dt)$acf[3] - acf(dt)$acf[2]^2)/(1-acf(dt)$acf[2]^2)

# Cramer's Rule to find ma_pacf3
# First, find matrix B, whose third column is replaced with solution vector

A <- toeplitz(acf(dt)$acf[1:3])

B <- A
B[,3] <- acf(dt)$acf[2:4]

ma_pacf3 = det(B)/det(A)


## c) AR(1)

fit %>% select(ar1) %>% report()

ar1_coef <- fit %>% select(ar1) %>% coef()
ar1_coef <- ar1_coef$estimate

ar1_sigmasq <- fit %>% select(ar1) %>% glance()
ar1_sigmasq <- ar1_sigmasq$sigma2

# k-step ahead forecast with 95% PI

ar1_pred1 <- ar1_coef * dt[20,2]
ar1_pred2 <- ar1_coef^2 *dt[20,2]

ar1_p1 <- ar1_sigmasq
ar1_p2 <- ar1_sigmasq * ((1-ar1_coef^4)/(1-ar1_coef^2))

# 95% P.I of 1-step and 2-step ahead forecast

ar1_pi1 <- c(ar1_pred1 - 1.96*sqrt(ar1_sigmasq),ar1_pred1 + 1.96*sqrt(ar1_sigmasq))
ar1_pi2 <- c(
	ar1_pred2 - 1.96*sqrt(ar1_sigmasq * ((1-ar1_coef^4)/(1-ar1_coef^2))),
	ar1_pred2 + 1.96*sqrt(ar1_sigmasq * ((1-ar1_coef^4)/(1-ar1_coef^2)))
	)
	
## e) ARMA(1,1)

fit %>% select(arma11) %>% report()

arma_coef <- fit %>% select(arma11) %>% coef()
arma_coef <- arma_coef$estimate

arma_sigmasq <- fit %>% select(arma11) %>% glance()
arma_sigmasq <- arma_sigmasq$sigma2

# generate MA terms

Z2 <- rep(0,n)

Z2[2] <- Y[2] - arma_coef[1]*Y[1]

for (i in 3:n) {
	Z2[i] = Y[i] - arma_coef[1]*Y[i-1] - arma_coef[2]*Z2[i-1]
}

# 1-step and 2-step prediction 

arma_pred1 <- arma_coef[1]*Y[20] + arma_coef[2] * Z[20]
arma_pred2 <- arma_coef[1]^2 * Y[20] + arma_coef[1] * arma_coef[2] * Z2[20]

arma_p1 <- arma_sigmasq
arma_p2 <- arma_sigmasq * (1 + (arma_coef[1] + arma_coef[2])^2)

# 95% PI of 1-step and 2-step prediction

arma_pi1 <- c(arma_pred1 - 1.96 * sqrt(arma_p1), arma_pred1 + 1.96 * sqrt(arma_p1))
arma_pi2 <- c(arma_pred2 - 1.96 * sqrt(arma_p2) , arma_pred2 +1.96 * sqrt(arma_p2))

## f) ARIMA(1,1,0)

fit %>% select(arima110) %>% report()

arima110_coef <- fit %>% select(arima110) %>% coef()
arima110_coef <- arima110_coef$estimate

arima110_sigmasq <- fit %>% select(arima110) %>% glance()
arima110_sigmasq <- arima110_sigmasq$sigma2

# 1-step and 2-step ahead forecast
arima110_pred1 <- (arima110_coef +1)*Y[n] - arima110_coef * Y[n-1]
arima110_pred2 <- (arima110_coef + 1) * arima110_pred1 - arima110_coef * Y[20]

arima110_p1 <- arima110_sigmasq
arima110_p2 <- arima110_sigmasq * (2+2*arima110_coef + arima110_coef^2)

# 95% P.I. for 1-step and 2-step ahead forecast

arima110_pi2 <- c(arima110_pred2 - 1.96 *sqrt(arima110_p2),  arima110_pred2 + 1.96*sqrt(arima110_p2))
arima110_pi1 <- c(arima110_pred1 - 1.96* sqrt(arima110_p1),  arima110_pred1 + 1.96*sqrt(arima110_p1))

