setwd("C:/Users/idkup/Desktop/Course studies/2022 TERM 2/STAT4005") #change wd


library(tidyverse) #load tidyverse
library(fpp3) #laod library to use tsibble

mkdata <- read_csv("monthly milk.csv") #read data
sum(is.na(mkdata)) #check for missing values and error

mkdata <- mkdata[1:168,1:2] #cut off last two irrelevant rows that return NA values

mkdata <- mkdata %>% 
	mutate(Month = yearmonth(Month))%>%
	as_tsibble(index = Month)	 #change tibble to tsibble
	
colnames(mkdata)<- c("Month", "Pounds Per Cow") #change colnames

autoplot(mkdata)

gg_season(mkdata)

gg_subseries(mkdata) +
	labs(title = "Monthly Milk Production", y="pounds per cow"
)

##STL Decomposition
	
stldcmp <- stl(mkdata, s.window = "periodic")
stldcmp #Conduct STL decomposition to obtain trend, season, and residual

plot(stldcmp)

##Decomposition using moving average filter

# Let d = 12. Set filter length to be 12. Apply 2X12ma
# Estimate T_hat

wgts = c(.5,rep(1,11),.5)/12 #construct vector of filter

mkdata_f <- stats::filter(mkdata, sides = 2, filter = wgts)
plot(mkdata_f)

#Estimating Seasonal Effect using T_hat

d = as.ts(mkdata) - as.ts(mkdata_f) #subtract trend from original time series
d.bar = mean(d, na.rm = TRUE)		#take mean of d

S.mat <- matrix(d-d.bar, ncol = 12, byrow = TRUE)  #Construct matrix to account for seasonability

S <- apply(S.mat, MARGIN = 2, FUN = mean, na.rm = TRUE) #Construct vector that contains seasonality effect for d=12
S <- rep(S, length(as.ts(mkdata))/12)  #Construct vector of seasonality length equal to observations

ts.plot(S, ylab = "Seasonality") 

#Re estimate trend using de-seasoned data, and 9-ma filter`

Q <- mkdata$"Pounds Per Cow" - S
mkdata2 <- mkdata #create one more tsibble, mkdata2
mkdata2[,2] = Q

wgts2 = rep(1/9,9) #Create filter vector for 9-ma

mkdata2.f <- stats::filter(mkdata2,  sides = 2, filter =wgts2)
plot(mkdata2.f)

#Investigating the Structure of the Noise

N = as.ts(mkdata) - mkdata2.f - S  #Subtract Trend and Seasonal Components from Original Data
qqnorm(N)
qqline(N)

# We may safely conclude that the residual is distributed approximately by a normal distribution





