###Q1) Draw Time series plot, ACF and PACF for the data.

y = c(1.33, -0.56, -1.31, -0.37, 0.05, 0.46, 2.00, -0.19,
	-0.25, 1.07, -0.17, 1.14, 0.63, -0.75, 0.15, 0.71, 0.45,
	-0.14, 0.57, 1.43)

#create time series
Y = tsibble( 
	Y = y,
	Index = seq_along(y),
	index = Index)
	
#Plot time series, ACF, and PACF

Y %>% gg_tsdisplay(Y, plot_type = "partial")

##Answer: Please check attached PNG image for plot

###Q2) Find moment estimates of sigma^2 and theta for MA(1) model

#Find sample autocovariance and sample autocorrelation
C1 <- acf(Y, type = "covariance")$acf[2]
r1 <- acf(Y)$acf[2]

theta_hat <- c(1-sqrt(1-4*r1^2), 1+sqrt(1-4*r1^2))/(2*r1)

#Since 27.5594216 >1, we reject since we need a stationary solution.
#Thus, the theta_hat is 0.03628

sigma_squared <- C1/0.03628507 

##Answer: sigma^2 = 0.6414742, theta_hat = 0.03628


###Q3) Find least squares estimates of phi1,phi2, sigma^2 of an AR(2) model.  Find 95% CI for each of phi1 and phi2.

#Create data frame for fitting AR(2) model using lm() function

df <- data.frame(
	y_regressand <- y[3:20],
	y_lag1 <- y[2:19],
	y_lag2 <- y[1:18]
)

#Regress y_regressand on y_lag1 and y_lag2, omit constant

fit.ls <- lm(y_regressand ~ y_lag1 + y_lag2 -1, data = df)
summary(fit.ls)

#LSE of sigma^2
sigma_squared.ls <- (sum(resid(fit.ls)^2))/16

##Finding Confidence Interval for phi1 and phi2
#Create Gamma Matrix
p <- acf(Y, 1, type = 'covariance')$acf
GAMMA <- matrix(c(p[1],p[2],p[2],p[1]), nrow =2)

#estimate standard error of phi1 and phi2
GAM.inv <- solve(GAMMA)
var.phi <- (sigma_squared.ls/20) * GAM.inv[1,1]

#Construct 95% CI for phi1 and phi2
cbind(fit.ls$coefficients - 2*sqrt(var.phi), fit.ls$coefficients + 2*sqrt(var.phi))

## Answer: phi1 = 0.2354, phi2 = -0.2230, sigma^2 = 0.7236343
## 95% CI of phi1: [-0.2396240, 0.7103563], 95% CI of phi2: [-0.6980286, 0.2519517]

###Q4) Find Yule-Walker Estimates of phi1 and phi2 for fitting AR(2) model

r1 <- acf(Y)$acf[2]
r2 <- acf(Y)$acf[3]

#Yule-Walker matrix
YW <- matrix(c(1,r1,r1,1),2)
rho <- c(r1,r2)

# Solve for phi
phi.yw <- solve(YW,rho)

##Answer: phi1 = 0.04718125, phi2 = -0.30200553


###Q5) Find CSS estimates of phi, theta, and sigma^2 for fitting ARIMA(1,1) model 

#Create function get.z

get.z <- function(arg) {
	phi = arg[1]
	theta = arg[2]
	z = rep(NA,length((y)))
	z[1] = y[1]
	for (j in seq.int(2,length(y))) {
		z[j] = y[j] - phi*y[j-1] - theta*z[j-1]
	}
	return(sum(z^2))
}

#Optimize get.z

results <- optim(c(0.1,0.1),get.z)

#find sigma2

sigma_squared.css <- results$value/(length(y)-1)

##Answer: phi = -0.5108904, theta = 0.8439972, sigma^2 = 0.6810138


###Q6) Find the MLE of phi, theta, sigma^2 of ARMA(1,1) model, w/o constant
fit.MLE <- arima(y, order=c(1,0,1), include.mean = FALSE)

#call loglikelihood
fit.MLE$loglik

##Answer: phi = -0.3883, theta = 1, max value of loglikelihood = -23.08444


###Q7) Among AR(p), p=1,2,3,4,5, which model is best in terms of FPE?
FPE<-rep(0,5)

#fit five AR models from p=1 to p=5, omit constant
for (p in 1:5) {
	fit.fpe<-arima(y,order=c(p,0,0),include.mean=F)
	sigma2<-fit.fpe$sigma2
	FPE[p]<-sigma2*(20+p)/(20-p)
}

bestmodel <- (1:5)[min(FPE)==FPE] 

##Answer: AR(1) is the best model based on FPE. FPE(AR(1)) = 0.8041636


###Q8) Among MA(q), q=1,2,3,4,5, which model is best in terms of AICC?

# fit various MA() models, omit constants

fit.aicc <- Y %>% model(
	arima001 = ARIMA(Y ~ 0 + pdq(0,0,1)),
	arima002 = ARIMA(Y ~ 0 + pdq(0,0,2)),
	arima003 = ARIMA(Y ~ 0 + pdq(0,0,3)),
	arima004 = ARIMA(Y ~ 0 + pdq(0,0,4)),
	arima005 = ARIMA(Y ~ 0 + pdq(0,0,5)),
	)

#Show AIC
glance(fit.aicc) %>% select(AICc) 

##Answer: The model of order 2 has the smallest AICc. Thus, MA(2) is the best model based on AICc 


###Q9) Fit MA(1) model Find Ljung_box test statistic Q(10). State H0 and H1, and then draw conclusion.

#Fit MA(1), exclude constant

fit.lj <- Y %>% model(
	arima001 <- ARIMA(Y ~ 0 + pdq(0,0,1))
	)

# extract residuals
fit.lj %>% augment() %>% select(.innov) -> ma_res
ma_res <- ma_res$.innov

#Conduct ljung_box test statistic for h=10

ljung_box(ma_res, lag = 10)
qchisq(0.95,9)

##Answer: H0: The series of residuals exhibit no autocorrelation for lag of 10.
## H1: Some autocorrelation coefficient p(k) for k=1,2,...,h are non-zero.
## Since test statistic is 11.7302485 and qchisq value for 95% percentile is 16.91898, we do not
## reject the null hypothesis. Thus, we can carefully speculate that the model is a good fit.


###Q10) Which model would I fit?


# Check ACF and PACF
Y %>% gg_tsdisplay(Y, plot_type = 'partial')

# Hard to know which model to use from visal inspection.

# We use brute force to iterate all ARIMA(p,0,q) models, calculating the IC values for each model
# p,q = 0,1,2,3,4,5. Total of 6*6 = 36 iterations

  
# Build function to calculate IC values
IC_values=function(x,order.vec){
	p=order.vec[1]
	q=order.vec[2]
	n=length(x)
	fit.q10=arima(x,order=c(order.vec[1],0,order.vec[2]),include.mean=F)
	sigsq=fit.q10$sigma2
	FPE=sigsq*(n+p)/(n-p)
	AIC=fit.q10$aic
	BIC=(n-p-q)*log(n*sigsq/(n-p-1))+n*(1+log(sqrt(2*pi)))+(p+q)*log((sum(x^2)-n*sigsq)/(p+q))
	return(c(AIC,BIC,FPE))
}

# Create empty vectors to populate
AR<-rep(0,35)
MA<-rep(0,35)
AIC<-rep(0,35)
BIC<-rep(0,35)
FPE<-rep(0,35)

# Create dataframe that lists orders of ARMA models along with IC scores
finaldf<-data.frame(AR, MA, AIC, BIC, FPE)

# Populate dataframe by iterating over 36 models
i<-1

for (p in 0:5) {
	for (q in 0:5) {
		finaldf[i,]<-c(p,q,IC_values(y,c(p,q)))
		i<-i+1
	}
}

# Remove NA values

finaldf<- drop_na(finaldf)

# Best model based on AIC

finaldf[finaldf$AIC == min(finaldf$AIC),]

# Best model based on BIC

finaldf[finaldf$BIC == min(finaldf$BIC),]

# Best model based on FPE

finaldf[finaldf$FPE == min(finaldf$FPE),]

##Answer: FPE: ARMA(2,4), AIC: ARMA(0,2), BIC: ARMA(1,4)
##Add-in:: Judging from the acf and pacf diagrams, we see lack of autocorrelation to actually suggest
## appropriate models beside an iid model. 