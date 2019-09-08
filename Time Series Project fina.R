```{r}
library ( astsa ) 
library( tseries ) 
library (MASS) 
library( forecast ) 
library (TSA) 
library ( ggplot2 )
import0 <- read.csv("C:\\Users\\zoreee\\Desktop\\monthly-australian-imports-from-.csv")
import0 <-  monthly.australian.imports.from.
import<-ts(import0[,2], start=c(1965,7), frequency=12)
# plot
par(mfrow=c(1,1))
plot(import,xlab="Time",ylab="import(Thousands of dollars)", main="Monthly Australian Imports From Japan")
# seasonal Plot
seasonplot(import,12,col=rainbow(3),year.labels = TRUE, main = "Seasonal Plot")
#decomposition plot
decom<-decompose(import1)
autoplot(decom, main="Decomposition plot")
#stablize variance 
#using box−cox
bxTransform <- boxcox(import ~ as.numeric(1:length(import)))
lamda <- bxTransform$x[which.max(bxTransform$y) ]
lamda # lamda=0.2222222
Transimport <- import^lamda #model transformation
#make sample set: drop last 10 data for later comparison of forecasting
tryimport <- ts (Transimport[1:(length(import) - 10)])
plot(tryimport, xlab="Time" , ylab="", main=expression(Import^0.2222))
acf(tryimport ,main="ACF of bos-cox transformed data")
#De−seasonalize
diff12 <- diff(tryimport ,lag=12)
plot(diff12 ,xlab="Time", ylab="", main=expression(nabla[12]~"Transformed Data"))
abline(lm(diff12~as.numeric(1:length(diff12))))
var(diff12)

#De−trend
diff12diff1<-diff(diff12,lag=1)
plot(diff12diff1 ,xlab="Time", ylab="", main=expression(nabla~nabla[12]~"Transformed Data"))
abline(lm(diff12diff1~as.numeric(1:length(diff12diff1))))
var(diff12diff1)
diff12diff2 <- diff(diff12diff1,lag=1)
var(diff12diff2)
# D-F test
adf.test(diff12diff1)
```
```{r}
#Model identifications 
#identify P, Q
op <- par ( mfrow=c ( 1 , 2 ) )
acf(diff12diff1 ,lag.max=60,main="")
pacf ( diff12diff1 , lag.max=60, main="")
title (main = "ACF and PACF of Deseasonalized Tansformed Data", outer = TRUE,
       line = -1) 
par(op)
#check value of acf , pacf , at lag=12, 24, 36 ,48...
#ACF cuts off after lag=24, then Q=2 
#PACF cuts off after lag=36, then P=3
```
```{r}
#identify p,q
op <- par ( mfrow=c ( 1 , 2 ) )
acf(diff12diff1 ,lag.max=11,main="")
pacf ( diff12diff1 , lag.max=11,main="")
title (main="ACF and PACF Plots for Lag Less Than 12" ,outer=TRUE, line=-1) 
par(op)
#at lag 1 ,2 ,3 ,... ,11
#acf cuts off after lag=1 or tails off , q=1 or q=0
#pacf cut off after lag=2 or tails off , p=2 or p=0
#test all combination of p in 0 to 2 and q in 0 to 1
```
```{r}
#Model selection by AICc
AICc<-numeric()
for (p in 0:2) {
  for (q in 0:1){
    AICc<-c(AICc,sarima(tryimport,p,1,q,2,1,3,12,details=FALSE)$AICc)
  }
}
AICc<-matrix(AICc,nrow=3,byrow=TRUE)
rownames(AICc)<-c("p=0", "p=1", "p=2")
colnames(AICc)<-c("q=0", "q=1")
AICc
AICc <- data.frame(AICc)
aicc <- setNames(AICc,c("q=0", "q=1"))
aicc
#smallest
#p=2,q=1 smallest 
#p=2,q=0 2nd smallest 
#AICc very close

#Model selection by BIC
BIC<-numeric()
for (p in 0:2) {
  for (q in 0:1){
    BIC<-c(BIC,sarima(tryimport,p,1,q,2,1,3,12,details=FALSE)$BIC)
  }
}
BIC<-matrix(BIC,nrow=3,byrow=TRUE)
rownames(BIC)<-c("p=0", "p=1", "p=2")
colnames(BIC)<-c("q=0", "q=1")
BIC
BIC <- data.frame(BIC)
bic <- setNames(BIC,c("q=0", "q=1"))
bic
#smallest
#p=2,q=1 smallest 
#p=2,q=0 2nd smallest 
#AICc very close
#Based on AICc and BIC, select two models 
#Model 1, SARIMA(2 ,1 ,1 ,2 ,1 ,3)12
#Model 2, SARIMA(2 ,1 ,0 ,2 ,1 ,3)12
```
```{r}
#Fit and Estimation based on MLE method 
# MODEL 1: p=2 and q =1 < SARIMA(2 ,1 ,1 ,2 ,1 ,3)12 >
fit1 <- arima(tryimport, order=c(2,1,1), seasonal=list(order=c(2,1,3), period=12) ,method="ML")
fit1

# Model 2 : p=2 and q=0 < SARIMA(2 ,1 ,0 ,2 ,1 ,3)12 >
fit2 <- arima(tryimport, order=c(2,1,0), seasonal=list(order=c(2,1,3), period=12) ,method="ML")
fit2

##Normality
resid1<-residuals(fit1) #resid for M1 
resid2<-residuals(fit2) #resid for M2
op<-par ( mfrow=c ( 2 , 2 ) )
hist(resid1 ,main="Histogram of Residuals Under Model 1") #not good enough
qqnorm( resid1 , main="Normal Q−Q Plot for Model 1")
qqline( resid1 )
hist(resid2 ,main="Histogram of Residuals Under Model 2") #good
qqnorm( resid2 , main="Normal Q−Q Plot for Model 2") 
qqline( resid2 )
qqline( resid1 )
par(op)
#Shapiro Test for Model 1 and 2
Shap<-matrix(c(shapiro.test(resid1)$statistic, shapiro.test(resid1)$p.value,shapiro.test(resid1)$statistic ,shapiro.test(resid2)$p.value),
             nrow=2,byrow=T)
#greater than 0.05 , then good
rownames(Shap)<-c("Model1" ,"Model2") 
colnames(Shap)<-c("W Statisttic","P−value")
Shap<-data.frame ( Shap )
Shap

```
```{r}
#Independence/Correlation diagnostics
b1<-Box.test(resid1, lag = 12, type = "Box-Pierce", fitdf = 2)$p.value
# Cor
b2<-Box.test(resid1 , lag = 12, type = "Ljung-Box", fitdf = 2)$p.value #Cor
b1 #>0.05 good
b2 #>0.05 good
b3<-Box.test(resid2, lag = 12, type = "Box-Pierce", fitdf = 2)$p.value
# Cor
b4<-Box.test(resid2 , lag = 12, type = "Ljung-Box", fitdf = 2)$p.value #Cor 
b3 #>0.05 
b4 #>0.05
boxT<-matrix(c(b1,b2,b3,b4) ,nrow=2,byrow=FALSE) 
rownames(boxT)<-c("Box−Pierce" ,"Ljung−Box") 
colnames(boxT)<-c("Model1 P−value" , "Model2 P−value")
boxT

#Test for constant variance of residuals
par(mfrow=c(2 ,2) ) # acf
acf ( resid1 , main = "ACF Plot of Residuals for Model 1" , lag.max=40) # pacf
pacf ( resid1 , main="" , lag.max=40)
title(main="PACF Plots of Residuals for Model 1",outer=FALSE,line=1) # acf
acf ( resid2 , main = "ACF Plot of Residuals for Model 2" , lag.max=40) # pacf
pacf ( resid2 , main="" , lag.max=40)
title (main="PACF Plot of Residuals for Model 2" ,outer=FALSE, line=1)
par(op)

```
```{r}
##Model 2 since AIC and BIC smaller
#forcasting based on Final model
pred.tr <- predict(fit2 ,n.ahead = 10)
U.tr= pred.tr$pred + 2*pred.tr$se # upper bound for the C. I . for transformed data
L.tr= pred.tr$pred - 2*pred.tr$se # lower bound for the C. I . for transformed data
ts.plot(tryimport , xlim=c(1,length(tryimport)+10),main="Forcasting Based on Transform Data",ylab="")
lines (U.tr , col="blue" , lty="dashed")
lines (L.tr , col="blue" , lty="dashed") 
points((length(tryimport)+1):(length(tryimport)+10), pred.tr$pred, col="red")
pred.orig <- pred.tr$pred^(1/lamda) # back−transform to get predictions of original time series
U= U.tr^(1/lamda) # bounds of the confidence intervals 
L= L.tr^(1/lamda) # Plot forecasts with original data
import2<-ts(import0[ , 2 ] )
ts.plot(import2 , xlim=c(1,length(import2)) ,main="Forcasting Based on Original
        Data",ylab="import")
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed") 
points((length(tryimport)+1):(length(tryimport)+10), pred.orig ,col="red")
#zoom effect
ts.plot(import2 , xlim=c(length(import2)-20,length(import2)) ,main="Comparison between Observed Values and Forcasted Values",ylab="import")
points((length(tryimport)+1):(length(tryimport)+10),import2[331:340], col="dark green")
points((length(tryimport)+1):(length(tryimport)+10),pred.orig , col="red") 
lines((length(tryimport)+1):(length(tryimport)+10),U, lty=2, col="blue") 
lines((length(tryimport)+1):(length(tryimport)+10),L, lty=2, col="blue")
#close to observed value . within confidence interval , good forcasting

```

```
