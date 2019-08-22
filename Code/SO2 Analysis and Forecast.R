library(ggplot2)
library(Rssa)
library(lattice)
library(plyr)
library(dplyr)
library(forecast)

# Import Data
SO2 <- read.csv("F:/study/myuoft/STA457/Zekun_Liu_1003594149_project/code/SO2.csv")
SO2_log <- log(as.vector(SO2[['SO2.AQI']]))

time <- as.vector(SO2[['Date.Local']])
time <- as.Date(time)

N <- length(SO2_log) # total length of the dataset
testing <- 3*365 # remove the last 3 year as validation testing
future <- 2*365 # add a new year for prediction
len <- testing + future # total length for forecast model
r <- 1
L <- 4*365 #window length

# separate training data index
training_date <- window(time, end = N-testing)

# Embedded & Decomposition

## Determination of Window Length
training_de <- ssa(window(SO2_log, end = N-testing), L= 730, neig=730, kind="1d-ssa", svd.method = "svd")
training_de_1 <- ssa(window(SO2_log, end = N-testing), L= 365,kind="1d-ssa", svd.method = "svd")
plot(wcor(training_de_1, groups = 1:30))
training_de_2 <- ssa(window(SO2_log, end = N-testing), L= 1095, kind="1d-ssa", svd.method = "svd")
plot(wcor(training_de_2, groups = 1:30))
training_de_4 <- ssa(window(SO2_log, end = N-testing), L= 5*365, kind="1d-ssa", svd.method = "svd")
plot(wcor(training_de_4, groups = 1:30))

## From W-correlation, window length can be determined as 1460
## Visualize the decomposed eigenvectors
training_de_3 <- ssa(window(SO2_log, end = N-testing), L= 1460, kind="1d-ssa", svd.method = "svd")
plot(training_de_3, type = "vectors", idx = 1:20)
plot(training_de_3, type="pair", idx=1:20)
plot(wcor(training_de_3, groups = 1:30))

## Further decompose the weakly separate eigenvectors
training_fos <- fossa(training_de_3, nested.groups = list(c(1:15)), gamma = 1000,
                      normalize = FALSE)
plot(training_fos, type = "vectors", idx = 1:13)
plot(training_fos, type="pair", idx=1:20)
plot(wcor(training_fos, groups = 3:13))

## Frequency estimation of the extracted seasonalities
est_train <- parestimate(training_fos, groups = list(c(1:16)),
                        method = "esprit")
plot(est_train)

# Grouping and Reconstruction
training_re <- reconstruct(training_de_3, groups = list(c(1:15)))
## Extraction of total tendency & random residuals
total_re <- training_re$F1
total.resid <- residuals(training_re)

## Visualize the extracted result
theme1 <- simpleTheme(col = c("grey", "blue","red"),
                      lwd = c(2, 1, 1),
                      lty = c("solid", "solid", "solid"))

xyplot(window(SO2_log, end = N-testing) + total_re + total.resid ~ time(window(SO2_log, end = N-testing)),
       xlab = "", ylab = "", type = "l", lwd = c(2, 1, 1),
       col = c("grey", "blue","red","black", "yellow"),
       auto.key = list(text = c("Full SO2 AQI series",
                                "Total Trend",
                                "Residuals"),
                       type = c("l", "l", "l"),
                       lines = TRUE, points = FALSE,
                       space = "top",
                       title = "Decompositon of Total Trend and Residuals"),
       par.settings = theme1)

# Further extracting slow-varying trend and seasonality
## Extraction of slow-varying trend
trend_extract <- reconstruct(training_de_3, groups = list(1))
trend <- trend_extract$F1
trend.resid <- residuals(trend_extract)

theme2 <- simpleTheme(col = c("grey", "blue","red"),
                      lwd = c(1, 2, 1),
                      lty = c("solid", "solid", "solid"))

xyplot(window(SO2_log, end = N-testing) + trend + trend.resid ~ time(window(SO2_log, end = N-testing)),
       xlab = "", ylab = "", type = "l", lwd = c(2, 1, 1),
       col = c("grey", "blue","red","black", "yellow"),
       auto.key = list(text = c("Full SO2 AQI series",
                                "Slow-varying Trend",
                                "Residuals"),
                       type = c("l", "l", "l"),
                       lines = TRUE, points = FALSE,
                       space = "top"),
       par.settings = theme2)

## Extraction of Seasonality
## Window Length = (total-length - 4*period)/2
resid_de <- ssa(trend.resid, L = 1701, kind="1d-ssa", svd.method = "svd")
plot(resid_de, xlim = 1:50)
plot(resid_de, type = "vectors", idx = 1:30)
plot(wcor(resid_de, groups = 1:30))

### Further decompose the weakly separate eigenvectors
trend_fos <- fossa(resid_de, nested.groups = list(c(1:16)), gamma = 1000,
                   normalize = FALSE)
plot(trend_fos, type = "vectors", idx = 1:30)

### Frequency Estimation for seasonalities
est_seas <- parestimate(trend_fos, groups = list(c(1:16)),
                   method = "esprit")

## Grouping and Reconstruction
trend_re <- reconstruct(trend_fos, groups = list(c(1:16)))
resid_seas <- trend_re$F1
resid_resid <- residuals(trend_re)

theme3 <- simpleTheme(col = c("grey", "blue","red","black"),
                      lwd = c(1, 2, 1, 1),
                      lty = c("solid", "solid", "solid","solid"))

xyplot(window(SO2_log, end = N-testing) + trend + resid_resid + resid_seas  ~ time(window(SO2_log, end = N-testing)),
       xlab = training_date, ylab = "", type = "l", lwd = c(1, 2, 1, 1),
       col = c("grey", "blue","red","black", "yellow"),
       auto.key = list(text = c("Full SO2 AQI series",
                                "Slow-varying Trend",
                                "Seasonality",
                                "Residuals"),
                       type = c("l", "l", "l","l"),
                       lines = TRUE, points = FALSE,
                       space = "top",
                       title="Full Decomposition of Trend, Seasonality and Trend"),
       par.settings = theme3)

xyplot(window(SO2_ori, end = N-testing) + exp(trend) + exp(resid_seas)  ~ time(window(SO2_log, end = N-testing)),
       xlab = "", ylab = "", type = "l", lwd = c(1, 2, 1),
       col = c("grey", "blue","red","black", "yellow"),
       auto.key = list(text = c("Full SO2 AQI series",
                                "Slow-varying Trend",
                                "Seasonality"),
                       type = c("l", "l", "l"),
                       lines = TRUE, points = FALSE,
                       space = "top"),
       par.settings = theme3)

# Construct LRR
lrr.coef <- lrr(training_fos, groups = list(c(1:16)))
lrr.roots <- roots(lrr.coef)
plot(lrr.roots)
len = 3*365

# Forecasting
fore_fos <- rforecast(training_fos, groups = list(trend = c(1:15)), len = len, only.new = TRUE)
fore_fos_trend <- c(rep(NA, length(window(SO2_log, end = N-testing))), fore_fos)

fore_fos_vec <- vforecast(training_fos, groups = list(trend = c(1:15)), len = len, only.new = TRUE)
fore_fos_vec_trend <- c(rep(NA, length(window(SO2_log, end = N-testing))), fore_fos_vec) 

# Visualize the prediction
time_index <- seq(as.Date('2000-01-01'),as.Date('2016-04-23'),by = 1)
df_1 = data.frame(time_index, SO2_log, fore_fos_trend, fore_fos_vec_trend, total_re_trend$F1)

ggplot(aes(x = time), data = df_1) + 
  geom_line(aes(y = SO2_log),col = "grey") +
  geom_line(aes(y = fore_fos_trend), col="blue")+
  geom_line(aes(y = fore_fos_vec_trend), col="red")+
  geom_line(aes(y = total_re_trend$F1), col="black")+
  xlab("") +
  scale_y_continuous("Log(SO2 AQI)", limits = c(-1,5))+
  ggtitle("Forecast Result with Original Extracted Trend")



# Tuning the best parameter for window length and # of leading components
# By finding the best combination of least validation error
# With a sliding window 
forecast.mse <- function(x, F.check, forecast.len = 1, ...){
  stopifnot(length(F.check) == forecast.len)
  f <- forecast(x, h = forecast.len, ...)
  mean((f$mean - F.check)^2)
}

forecast.sliding.mse <- function(F, L, ncomp, forecast.len = 1,
                                 K.sliding = N %% 4,
                                 .progress = "none",
                                 .parallel = FALSE,
                                 ...) {
  N <- length(F)
  sliding.len <- N - K.sliding - forecast.len + 1
  L.max = max(L); L.min = min(L); ncomp.max = max(ncomp)
  stopifnot(sliding.len > L.max)
  stopifnot(ncomp.max + 1 < min(L.min, N - L.max + 1))
  g <- expand.grid(L = L, i = 1:K.sliding)
  aaply(g, 1, 
        splat(function(L, i){
          F.train <- F[seq(from = i, len = sliding.len)]
          F.check <- F[seq(from = sliding.len + i,
                           len = forecast.len)]
          s <- ssa(F.train, L = L)
          sapply(ncomp, function(ncomp){
            res <- forecast.mse(s, F.check,
                                forecast.len = forecast.len,
                                groups = list(1:ncomp),
                                ...)
            names(res) <- as.character(ncomp)
            res
          })
        }),
        .progress = .progress, .parallel = .parallel)
}


optim.par <- function(m0){
  m <- apply(m0, c(1, 3), mean)
  mpos <- which(m == min(m), arr.ind = TRUE)
  L.opt <- Ls[mpos[1]]
  ncomp.opt <- ncomp[mpos[2]]
  list(L.opt = L.opt, ncomp.opt = ncomp.opt, m = m)
}

# minimize RMSE errors of two forecasts for 365 steps ahead
total <- length(SO2_log)
K.sliding <- 2

# set validation set for 3 years
forecast.base.len <- 1095
base.len <- length(SO2_log)
sliding.len <- base.len - K.sliding-forecast.base.len + 1
# number of leading composition
ncomp <- 1:100
L.min <- 365
Ls <- seq(L.min, 6*L.min, by = 365)
m0 <- forecast.sliding.mse(SO2_log, K.sliding = K.sliding,
                           L = Ls, ncomp = ncomp, 
                           method = "recurrent",
                           forecast.len = forecast.base.len,
                           .progress = "none")

m1 <- forecast.sliding.mse(SO2_log, K.sliding = 3,
                           L = Ls, ncomp = 1:50, 
                           method = "recurrent",
                           forecast.len = forecast.base.len,
                           .progress = "none")
# By setting k =2 AND 3, same result is returned
p1 <- optim.par(m1)
print(c(p1$L.opt, p1$ncomp.opt, sqrt(min(p1$m))))

p <- optim.par(m0)

print(c(p$L.opt, p$ncomp.opt, sqrt(min(p$m))))
# 1460, 17, .542
matplot(Ls, sqrt(p$m), ylab = "", xlab = "Window lengths",
        type = "l", col = topo.colors(100))

# Predict the 5 years value in the future
forecast.len = 5*365
ssa.obj <- ssa(window(SO2_log, end = total - 1095), L = p$L.opt)
plot(ssa.obj)
plot(ssa.obj, type = "vector", idx = 1:20)
plot(wcor(ssa.obj, groups = 1:30))

ssa.re <- reconstruct(ssa.obj, groups = list(1:17))

ssa.for <- rforecast(ssa.obj, groups = list(1:17),
                     len = forecast.len)

ssa.for.total <- rforecast(ssa.obj, groups = list(1:17),
                     len = forecast.len, only.new = FALSE)

# Forecast with 95% Confidence Interval
ssa.for.con <- forecast(ssa.obj, groups = list(1:17),len=forecast.len,
                        method = "recurrent", interval = "prediction",
                        level = .95, R=100)

plot(ssa.for.con, include = 36, shadecols = "grey", type = "l",
     main = "Confidence intervals",title = "Five-Year Forecast with 95% Confidence Interval")

original_data <- c(SO2_log, rep(NA, 2*365))
original_trend <- c(total_re_trend$F1, rep(NA, 2*365))
training_fore <- c(ssa.re$F1, rep(NA, length(ssa.for)))
forecast_data <- c(rep(NA, length(ssa.re$F1)), ssa.for)
upper <- c(rep(NA, length(ssa.re$F1)), ssa.for.con$upper)
lower <- c(rep(NA, length(ssa.re$F1)), ssa.for.con$lower)

time_index_2 <- seq(as.Date('2000-01-01'),as.Date('2018-04-23'),by = 1)

df2 <- data.frame(original_data, original_trend,training_fore,forecast_data,upper, lower)

# Plot the predictions and original data
ggplot(aes(x = time_index_2), data = df2) + 
  geom_line(aes(y = original_data),col = "grey") +
  geom_line(aes(y = original_trend), col="red")+
  geom_line(aes(y = training_fore), col="black")+
  geom_line(aes(y = forecast_data), col="blue")+
  xlab("") +
  ggtitle("Forecast of Optimal Model and Original Data") +
  scale_y_continuous("Log(SO2 AQI)", limits = c(-1,5)) 

# Transform back to original data
SO2_ori <- as.vector(SO2[['SO2.AQI']])
trend_ori <- exp(total_re_trend$F1)
fitting_ori = exp(ssa.re$F1)
forecast_ori = exp(ssa.for.total)

SO2_ori_data <- c(SO2_ori, rep(NA, 2*365))
ori_trend <- c(trend_ori, rep(NA, 2*365))
fitting <- c(fitting_ori, rep(NA, length(ssa.for)))


df3 = data.frame(SO2_ori_data,ori_trend, fitting, forecast_ori)

ggplot(aes(x = time_index_2), data = df3) + 
  geom_line(aes(y = SO2_ori_data),col = "grey") +
  geom_line(aes(y = ori_trend), col="red")+
  geom_line(aes(y = forecast_ori), col="blue")+
  geom_line(aes(y = fitting), col="black")+
  ylab("SO2 AQI in New York") +
  ggtitle("Fitting and Forecast with Original Data Without Transformation")  



