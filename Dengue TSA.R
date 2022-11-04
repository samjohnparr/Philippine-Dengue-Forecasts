library(readxl)
library(tseries)
library(forecast)
library(lmtest)
library(mice)
library(aTSA)
library(fGarch)
library(Metrics)

Dengue <- read_excel("D:/Data for Research/DOH/Diseases Statistics [for TSA].xlsx", sheet = "Will Be Used") #Check if there's Week 52

Dengue_TS <-ts(Dengue$Treated, frequency = 52, start = c(2017,1), end = c(2022, 40))
Dengue_df <- data.frame(Cases = as.matrix(Dengue_TS), Year = time(Dengue_TS))
Dengue_Imputed <- mice(Dengue_df, m = 5, maxit = 50, method = 'pmm')
Complete <- complete(Dengue_Imputed)
#write.csv(Complete,"D:/Data for Research/DOH/Imputed_Dengue_Statistics.csv")
From_CSV <- read.csv("D:/Data for Research/DOH/Imputed_Dengue_Statistics.csv")
Dengue_Complete <- ts(From_CSV$Cases, frequency = 52, start = c(2017,1), end = c(2022, 40))
ts.plot(Dengue_Complete, xlab = "Year", ylab = "Number of Cases", main = "Dengue Cases in the Philippines", col = "blue", lwd = 2)

#Data_Splitting; Model Building Set ###### REDO AND USE BOX-COX TRANSFORMATION

Dengue_Train_Set <-ts(Dengue_Complete, frequency = 52, start = c(2017,1), end = c(2021, 32))
tseries::adf.test(Dengue_Train_Set)
Acf(Dengue_Train_Set)
Pacf(Dengue_Train_Set)


#Differencing

Dengue_First_Diff <- diff(Dengue_Train_Set,1)
ts.plot(Dengue_First_Diff)
tseries::adf.test(Dengue_First_Diff)
Acf(Dengue_First_Diff) #q = 1,2,3 <> Q = 1
Pacf(Dengue_First_Diff) #p = 1,2 <> P = 1

#Data Validation Set

Excel_Dengue_Validation_Set <- read_excel("D:/Data for Research/DOH/Diseases Statistics [for TSA].xlsx", sheet = "Will Be Used", range = "B242:D301", col_names = FALSE)
Dengue_Validation_Set <- ts(Excel_Dengue_Validation_Set$...3, frequency = 52, start = c(2021,33), end = c(2022, 40))

#ARIMA Model

model_arima = auto.arima(Dengue_Train_Set, d = 1, ic = "aic", test = "adf", trace = TRUE, seasonal = FALSE) #4584.405
#auto.arima(Dengue_Train_Set, ic = "aic", test = "adf", trace = TRUE, seasonal = FALSE, lambda = 0, allowdrift = FALSE, allowmean = FALSE, d = 1)

checkresiduals(model_arima)
coeftest(model_arima)

#Residuals vs Time
plot(model_arima$residuals,type="p",ylab="Residuals",main="Plot of Residuals vs. Time")
abline(0,0) #Independent Against Time

#Residuals vs Fitted
plot(as.numeric(model_arima$fitted),model_arima$residuals,xlab="Fitted", ylab="Residuals", main="Plot of Residuals vs. Fitted Value")
abline(0,0) #Nonconstant Variance

#Normality
qqnorm(model_arima$residuals, main="Normal Q-Q Plot of Residuals")
qqline(model_arima$residuals)
shapiro.test(model_arima$residuals)

#ADF
tseries::adf.test(model_arima$residuals)

#Acf
Acf(model_arima$residuals)
Pacf(model_arima$residuals)

Box.test(model_arima$residuals,type = "Lj")

#ARCH test
arch.test(arima(Dengue_Train_Set, order = c(1,1,2))) #The Lagrange-Multiplier Test Revealed the presence of Heteroskedsticity

#Building the GARCH model on ARIMA

arima_r_squared = (model_arima$residuals)^2

tsdisplay(arima_r_squared)

Acf(arima_r_squared)

Result=data.frame(Model="m",AIC=0)
q=0
for (i in 1:3){
  for (j in 1:3){
    q=q+1
    fit=garchFit(substitute(~garch(p,q),list(p=i,q=j)),data=model_arima$residuals,trace=F)
    
    Result=rbind(Result,data.frame(Model=paste("m-",i,"-",j),AIC=fit@fit$ics[1]))
    
  }
}

Result=Result[2:nrow(Result),]
Result[which.min(Result$AIC),]

garch_fit <- garchFit(~garch(1,1), model_arima$residuals, trace=F)
summary(garch_fit)


garch_bound <- garch_fit@sigma.t  # conditional sd
arima_residualsplus <- model_arima$residuals + 1.96 * garch_bound
arima_residualsmin <- model_arima$residuals - 1.96 * garch_bound
plot(model_arima$residuals, main="GARCH(1,1)")
lines(arima_residualsplus, lty= 2, col = 3)
lines(arima_residualsmin, lty= 2, col = 3)

# Forecast for 60 periods
garch_predict <- predict(garch_fit, n.ahead = 60)

#ARIMA-GARCH

#Data Validation
newmodel <- Arima(c(Dengue_Complete),model = model_arima)
onestep_from_arima <- fitted(newmodel)[241:300]
onestep <- onestep_from_arima+garch_predict$meanForecast

#Actual vs Forecasted

argarch_forecasted <- ts(onestep,frequency = 52, start = c(2021,33), end = c(2022,40))
ts.plot(Dengue_Validation_Set, argarch_forecasted, xlab="Year",ylab="Cases",col=c("blue","red"),lwd=2)+legend("topleft",bty="o",lty=c(1,1),col=c("blue","red"),legend=c("Actual","Forecasted"),cex=0.7,inset=0.025,lwd=2)

#Forecast Errors
FE = Dengue_Validation_Set-onestep

eval=cbind(Dengue_Validation_Set,onestep,FE)
colnames(eval)=c("Actual", "Forecast", "Forecast Error")
eval #THE FORECAST ERRORS ARE HIGH

#ACF and PACF of forecast errors
Acf(FE,main="")
Pacf(FE,main="")

#ARIMA-GARCH Forecast

arima_forecast_12 <- forecast::forecast(Dengue_Complete, model = newmodel, h = 12)

garch_fit_forecast_12 <- garchFit(~garch(1,1), newmodel$residuals, trace=F)
garch_predict_12 <- predict(garch_fit_forecast_12, n.ahead = 12)

arima_garch_forecast_12 <- arima_forecast_12$mean + garch_predict_12$meanForecast

#RMSE

rmse(Dengue_Validation_Set, onestep)
mape(Dengue_Validation_Set, onestep)

#Holt-Winters

holt_fit_add<-HoltWinters(Dengue_Train_Set, seasonal = "additive")
checkresiduals(holt_fit_add)

holt_fit_mult<-HoltWinters(Dengue_Train_Set, seasonal = "multiplicative")
checkresiduals(holt_fit_mult)

forecast_holt_fit_add <- forecast::forecast(holt_fit_add, h = 60) #Better Model

Acf(forecast_holt_fit_add$residuals, lag.max=20, na.action=na.pass)
Box.test(forecast_holt_fit_add$residuals, type="Ljung-Box")
hist(forecast_holt_fit_add$residuals, xlab = "Residuals")
shapiro.test(forecast_holt_fit_add$residuals)

rmse(Dengue_Validation_Set,forecast_holt_fit_add$mean) #4486.413
mape(Dengue_Validation_Set,forecast_holt_fit_add$mean) #2.506484

forecast_holt_fit_mult <- forecast::forecast(holt_fit_mult, h = 60)

Acf(forecast_holt_fit_mult$residuals, lag.max=20, na.action=na.pass)
Box.test(forecast_holt_fit_mult$residuals, type="Ljung-Box")
hist(forecast_holt_fit_mult$residuals, xlab = "Residuals")
shapiro.test(forecast_holt_fit_mult$residuals)

rmse(Dengue_Validation_Set,forecast_holt_fit_add$mean) 

#Holt-Winters Forecasts 12-week period

holt_add_complete <- HoltWinters(Dengue_Complete, seasonal = "additive", alpha = 0.4040651, beta = 0, gamma = 0.6998076)
holt_add_complete_forecasts <- forecast::forecast(holt_add_complete, h = 12)

#Actual vs ARIMA vs SARIMA vs ARIMA-GARCH vs SARIMA-GARCH

ts.plot(Dengue_Complete, arima_garch_forecast_12 , holt_add_complete_forecasts$mean, col = c("black","blue","red"),lwd = c(1,2,2))
legend("topleft",bty="n",lty=c(1,1),col=c("black","blue","red"),legend=c("Actual","ARIMA-GARCH", "HW"),cex=0.7,inset=0.025,lwd=2)


