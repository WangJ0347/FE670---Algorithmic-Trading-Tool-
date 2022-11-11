rm(list = ls())
library(quantmod)
library(TTR)
library(PerformanceAnalytics)
library(tseries)
library(e1071) 
library(forecast)
library(rugarch)
library(PeerPerformance)
library(xts)
library(readr)
library(HighFreq)
library(nortest)
library(boot)

#Getting stock prices of BIO
getSymbols("BIO",from="2017-01-01",to="2021-12-31")

#Getting return
Ad_justed <- BIO$BIO.Adjusted
returns <- periodReturn(BIO, period="daily",type='log')
returns <- returns[-1]

#Augmented Dickey-Fuller Test
L <- length(returns)
adf.test(returns, k = trunc((L-1)^(1/3)))

#TS model: AR(1)
model <- ar(returns, aic=FALSE, order.max = 1)
ar.beta <- model$ar
res.model <- residuals(model)[-1]

#Function to simulate SMA model
sma <- function(price.adjusted,sma.short,sma.long){
  #Simple Moving Average trategy
  sma20 <- SMA(price.adjusted, n = sma.short)
  sma100 <- SMA(price.adjusted, n = sma.long)
  
  #Trading signals (-1: Sell, 1: buy, 0: nothing)
  #SMA 20 and SMA 100 crossover signal
  sma_ts <- Lag(
    ifelse(Lag(sma20)<Lag(sma100)&sma20>=sma100,1,
           ifelse(Lag(sma20)>Lag(sma100)&sma20<=sma100,-1,0)))
  sma_ts[is.na(sma_ts)]<- 0
  
  #Security position (1: hold stock, 0: don't own stock)
  sma_strat <- ifelse(sma_ts>1,0,1)
  for (i in 1:length(Ad_justed)) {
    sma_strat[i]<-ifelse(sma_ts[i]==1,1,ifelse(sma_ts[i]==-1,0,sma_strat[i-1]))
  }
  sma_strat[is.na(sma_strat)]<-1
  
  #daily profits and losses
  hb.ret <- periodReturn(price.adjusted,period = "daily",type = "log")
  sma.ret <- hb.ret*sma_strat
  
  out_put <- cbind(hb.ret,sma.ret,sma_ts)
  colnames(out_put) <- c("h&b_returns","sma_returns","position")
  out_put
 # Return.cumulative(sma.ret,geometric = FALSE)
}

sma20 <- SMA(BIO$BIO.Adjusted, n = 20)
sma100 <- SMA(BIO$BIO.Adjusted, n = 100)
lineChart(BIO,theme = 'black')
addSMA(n = 20, col = 'blue')
addSMA(n = 100, col = 'red')
legend('left', col = c('green','blue','red'),
       legend = c('BIO','SMA20','SMA100'),lty = 1,
       bty = 'n', text.col = 'white',cex = 0.8)



pn_l <- sma(Ad_justed,20,100)[,1:2]
charts.PerformanceSummary(pn_l, main = 'Performance of SMA Strategy')

#draw pnl distribution
po_sitions <- sma(Ad_justed,20,100)[,3]
daily_return <- pn_l[,2]
plot(density(daily_return), main = "sample P/L distribution")

#bootstrapping#####################
bt_statistics <- function(beta,data,indices,cl_ose,model_mean){
  tmp_res <- data[indices]
  price <- array(as.numeric(cl_ose[1]))
  price <- append(price, as.numeric((cl_ose[2])))
  re_turn <- as.numeric(log(cl_ose[2])) - as.numeric(log(cl_ose[1]))
  for(i in tmp_res){
    re_turn <- beta*(re_turn - model_mean) + i + model_mean
    price <- append(price, price[length(price)] * exp(re_turn))
  }
  price <- xts(price,order.by = index(cl_ose))
  #sma trading strategy
  tmp_sma <- sma(price,20,100)
  #get result of trading
  po_sitions <- tmp_sma[,3]
  daily_return <- tmp_sma[,2]
  hb.ret <- tmp_sma[,1]
  pn_l <- cumsum(daily_return)
  
  #criterias##########
  #trade_num
  trade_num <- sum(sign(abs(diff.xts(po_sitions)[-1])))
  #total return
  total_return <- Return.cumulative(daily_return,geometric = FALSE)
  #B&H return
  hb.cumret <- Return.cumulative(hb.ret,geometric = FALSE)
  #winning rate
  returns <- pn_l[diff.xts(po_sitions)!=0]
  returns_lag <- lag.xts(returns)
  winning <- sum(((returns-returns_lag)>0)[-1])
  winning_rate <- winning/trade_num
  #max drawdown
  max_drawdown <- maxDrawdown(daily_return)
  #sharpe ratio
  sharpe_ratio <- sqrt(252)*sum(daily_return)/sd(daily_return)/NROW(daily_return)
  
  return(cbind(tradeNum = trade_num,total_return=total_return,
               bh_return = hb.cumret,winning_rate=winning_rate,
               max_drawdown = max_drawdown, sharpe_ratio=sharpe_ratio))
}

result <- boot(data=res.model,statistic = bt_statistics,R=150,beta = ar.beta,cl_ose=Ad_justed,
               model_mean=model$x.mean)

result$t0
data <- result$t
data <- as.data.frame(data)
colnames(data) <- colnames(result$t0)
apply(data, 2, mean)
data <- rbind(result$t0, data, apply(data, 2, mean))
write.csv(data,file = 'bootstrap_result.csv')

###Draw Pnl distribution##################
plot(density(result$t[,2]), main = "P/L distribution of bootstrap resampling")

###Check if average total return differs from zero
t.test(result$t[,2])
###################################################
#########  TS Momentum strategy  ##################
###################################################

r_daily <- periodReturn(BIO, period="daily",type='log')
dayRlistAll <- NULL
dayRlistCond <- NULL
dayCntAll <-0
dayCntCond <- 0
sumLagD <- 252
sumLagD1 <- sumLagD + 1
Ldm <- L-21
for(i in sumLagD1:Ldm) {
  im <- i-sumLagD
  im1 <- i-1
  im21 <- i - 21
  ip21 <- i+21
  dayRlistAll <- c(dayRlistAll, sum(r_daily[i:ip21]))
  dayCntAll <- dayCntAll + 1
  if(sum(r_daily[im:im1]) > 0) { 
    dayRlistCond <- c(dayRlistCond, sum(r_daily[i:ip21]))
    dayCntCond <- dayCntCond +1  
  } else { dayRlistCond <- c(dayRlistCond, 0) } 
}

td <- t.test(dayRlistAll, dayRlistCond, paired=T)
td

sharpe.bh <- sqrt(12)*sharpe(dayRlistAll)
sharpe.tsm <- sqrt(12)*sharpe(dayRlistCond)

print("/")
print("daily:")
print(" R All/Cond pvalue")
print("Ann Sharpe")
print(c(mean(dayRlistAll)*12, mean(dayRlistCond)*12))
print(td$p.value)
print(c(sharpe.bh,sharpe.tsm))
print(c(dayCntAll, dayCntCond))

###################################################
###############  Statistics  ######################
###################################################
mean_d <- mean(daily_return)
sigma_d <- sd(daily_return)
skew_d <- skewness(daily_return)
kurt_d <- kurtosis(daily_return)
plot(density(daily_return),
     main = "sample daily returns distribution")

###################################################
#############  ARMA+GARCH(1,1)  ###################
###################################################
#Augmented Dickey-Fuller Test
L1 <- length(daily_return)
adf.test(daily_return, k = trunc((L-1)^(1/3)))

#Estimate ARMA coefficients
auto.arima(daily_return)

#estimates for all ARMA+GARCH coefficients
GARCHspec <- ugarchspec( variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(1, 1), include.mean = TRUE))
BIOfit <- ugarchfit(GARCHspec, returns)

