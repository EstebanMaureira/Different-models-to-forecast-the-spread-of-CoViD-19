####################################
# The CYO project
#" Different models to forecast the spread of covid "
# for edx caption
#
#created by: Esteban Maureira 07/01/2021
####################################

####################################
# Install libraries load and  #
####################################
# first we install and open the corresponding libraries

if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

if (!require("caret")) {
  install.packages("caret")
  library(caret)
}

if (!require("lubridate")) {
  install.packages("lubridate")
  library(lubridate)
}

if (!require("covid19.analytics")) {
  install.packages("covid19.analytics")
  library(covid19.analytics)
} # Install and analyze updated time series worldwide data of reported cases for the Novel CoronaVirus Disease (CoViD-19)

if (!require("tseries")) {
  install.packages("tseries")
  library(tseries)
}# install tseries tools

if (!require("fpp2")) {
  install.packages("fpp2")
  library(fpp2)} # install fpp2 for forecasting and modeling

if (!require("forecast")) {
  install.packages("forecast")
  library(forecast)
}# install forecast

if (!require("curl")) {
  install.packages("curl")
  library(curl)
}# install curl


####################################
# Load and prepare data #
####################################


# Data
tsc <- covid19.data(case = 'ts-confirmed') #download of updated covid19 data in confirmed cases

# Let's filter through Chile
tsc_chile <- tsc %>% filter(Country.Region == 'Chile')

# we convert data into analyzable data

tsc_chile <- data.frame(t(tsc_chile)) # we convert to dataframe and transpose
tsc_chile <- cbind(rownames(tsc_chile), data.frame(tsc_chile, row.names = NULL)) #swap the columns and add index
colnames(tsc_chile) <- c('Date', 'Confirmed') # we write titles to the columns
tsc_chile$Date <- ymd(tsc_chile$Date) # we adapt the date in ymd format
tsc_chile <- na.omit(tsc_chile) # we delete na the data
tsc_chile$Confirmed <- as.numeric(tsc_chile$Confirmed) # we convert the confirmed data into numeric
tsc_chile <- na.omit(tsc_chile)

# we convert to df
Date <- tsc_chile$Date
Confirmed_cases <- tsc_chile$Confirmed
df <- data.frame(date=Date,Confirmed_cases)
df <- df %>% filter(Confirmed_cases >0 )

df  %>% as_tibble() %>% head() # Head function to data

# declare this like a time series data with ts function from fpp2 package
Y <- ts(df[,2], frequency = 1) 

plot(df,type = "l", col = "grey", xlab = "Days", ylab = "Confirmed cases",
     main = "Time plot : Chile confirmed cases COVID-19")

####################################
# preliminary Analysis
####################################

autoplot(Y) + 
  ggtitle("Time plot : Chile confirmed cases COVID-19") + 
  ylab("Confirmed cases") # We can observe a clear increasing trend of confirmed cases

# As we saw, the data has a great trend. We will investigate the transformation.
# We will take the first difference to remove the trend.

DY <- diff(Y) # A non-stationary process with a deterministic trend becomes stationary after removing the trend, or detrending


autoplot(DY) + 
  ggtitle("Time plot :Change Chile confirmed cases COVID-19") + 
  ylab("Confirmed cases") # We still observe a trend, then we will apply a difference again to remove it definitively

DYY <- diff(DY) # A non-stationary process with a deterministic trend becomes stationary after removing the trend, or detrending

autoplot(DYY) + 
  ggtitle("Time plot :Change double difference Chile confirmed cases COVID-19") + 
  ylab("Confirmed cases") # Now we can see like stationary time-serie
# in the series begins with a burst of cases, but then normalizes, much like a trend-stationary.

#We can see the correlation with ACF graph.
ggAcf(Y, lag=48) # When data have a trend, the autocorrelations for small lags tend to be large and 
#positive because observations nearby in time are also nearby in size. So the ACF of trended time series tend 
#to have positive values that slowly decrease as the lags increase(2)



####################################
# Train and Test sets
####################################
# Let's create a train and test set to assess the accuracy of the models we implement

train <- subset(Y, end=length(Y)-31)
test <- subset(Y, start=length(Y)-30)
h <- length(test)


####################################
# Models
####################################
# we will analyze each model separately

################
# naive model
################
# For fit with Naive method wee need stationary time-serie. 
fit_naive <- naive(train,h)
print(summary(fit_naive))
checkresiduals(fit_naive) # We can see from the residual that there is no correlation throughout 

autoplot(subset(Y, start = 280)) +
  autolayer(fit_naive, series="naive model") +
  ggtitle("Forecasts from naive method") + xlab("Days") +
  ylab("Confirmed cases") +
  guides(colour=guide_legend(title="Forecast"))

################
# ETS model
################

fit_ets <- ets(train)
print(summary(fit_ets))
checkresiduals(fit_ets)


autoplot(subset(Y, start = 280)) +
  autolayer(forecast(fit_ets, h), series="ETS model") +
  ggtitle("Forecasts from ETS method") + xlab("Days") +
  ylab("Confirmed cases") +
  guides(colour=guide_legend(title="Forecast"))


################
# ARIMA model
################

fit_arima <- auto.arima(train, d = NA ,D = NA, stepwise = FALSE, approximation = FALSE, trace = TRUE , seasonal=FALSE, stationary = FALSE) 
print(summary(fit_arima)) 
checkresiduals(fit_arima)

autoplot(subset(Y, start = 280)) +
  autolayer(forecast(fit_arima, h), series="ARIMA model") +
  ggtitle("Forecasts from ARIMA method") + xlab("Days") +
  ylab("Confirmed cases") +
  guides(colour=guide_legend(title="Forecast"))

################
# Holt’s linear trend method
################
fit_holt <- holt(train, h) 
print(summary(fit_holt))
checkresiduals(fit_holt)

autoplot(subset(Y, start = 280)) +
  autolayer(fit_holt, series="Holt's method") +
  ggtitle("Forecasts from Holt's method") + xlab("Days") +
  ylab("Confirmed cases") +
  guides(colour=guide_legend(title="Forecast"))

################
# Drift method
################
fit_rwf <- rwf(train, h, drift=TRUE)
print(summary(fit_rwf)) 
checkresiduals(fit_rwf)

autoplot(subset(Y, start = 280)) +
  autolayer(fit_rwf, series="Drift method") +
  ggtitle("Forecasts from Drift method") + xlab("Days") +
  ylab("Confirmed cases") +
  guides(colour=guide_legend(title="Forecast"))

################
# TBATS  method
################

fit_tbats <- tbats(train)
print(summary(fit_tbats)) 
checkresiduals(fit_tbats)

autoplot(subset(Y, start = 280)) +
  autolayer(forecast(fit_tbats, h), series="TBATS method") +
  ggtitle("Forecasts from TBATS method") + xlab("Days") +
  ylab("Confirmed cases") +
  guides(colour=guide_legend(title="Forecast"))


################
# Combination
################
# We know that unity is strength

Combination <- (forecast(fit_ets, h)[["mean"]] + forecast(fit_arima, h)[["mean"]] +
                  fit_holt[["mean"]] + fit_rwf[["mean"]])/4


c(ETS = accuracy(forecast(fit_ets, h), Y)["Test set","RMSE"],
  ARIMA = accuracy(forecast(fit_arima, h), Y)["Test set","RMSE"],
  `HOLT` = accuracy(fit_holt, Y)["Test set","RMSE"],
  Drift = accuracy(fit_rwf, Y)["Test set","RMSE"],
  Combination =
    accuracy(Combination, Y)["Test set","RMSE"])

autoplot(subset(Y, start = 280)) +
  autolayer(forecast(fit_ets, h), series="ETS", PI=FALSE) +
  autolayer(forecast(fit_arima, h), series="ARIMA", PI=FALSE) +
  autolayer(fit_holt, series="HOLT", PI=FALSE) +
  autolayer(fit_rwf, series="Drift", PI=FALSE) +
  autolayer(Combination, series="Combination") +
  xlab("Days") + ylab("Cases") +
  ggtitle("Forecast cases COVID19 confirmed in Chile")

################
# Resume
################

autoplot(subset(Y, start = 280)) +
  autolayer(fit_naive, series="Naive", PI=FALSE) +
  autolayer(forecast(fit_ets, h), series="ETS", PI=FALSE) +
  autolayer(forecast(fit_arima, h), series="ARIMA", PI=FALSE) +
  autolayer(fit_holt, series="HOLT", PI=FALSE) +
  autolayer(fit_rwf, series="Drift", PI=FALSE) +
  autolayer(forecast(fit_tbats, h), series="TBATS", PI=FALSE) +
  autolayer(Combination, series="Combination") +
  xlab("Days") + ylab("Cases") +
  ggtitle("Forecast cases COVID19 confirmed in Chile")


rmse_results <- bind_rows(method = c("Naive","ETS","ARIMA","HOLT","Drift","TBATS","Combination"),
                          accuracy = c(accuracy(fit_naive, Y)["Test set","RMSE"],
                                       accuracy(forecast(fit_ets, h), Y)["Test set","RMSE"],
                                       accuracy(forecast(fit_arima, h), Y)["Test set","RMSE"],
                                       accuracy(fit_holt, Y)["Test set","RMSE"],
                                       accuracy(fit_rwf, Y)["Test set","RMSE"],
                                       accuracy(forecast(fit_tbats, h), Y)["Test set","RMSE"],
                                       accuracy(Combination, Y)["Test set","RMSE"]))
rmse_results



c(RMSE= rmse_results, caption = 'Summary RMSE')