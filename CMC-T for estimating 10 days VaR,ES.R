rm(list = ls())
require(dplyr)
require(readxl)
require(rugarch)
library(parallel)
################################################################################################
#  1990.01.01~2023.12.31 10-days Value at Risk
###############################################################################################

sp_dat <- read_xlsx("C:/Users/seoyeji/Desktop/spx_100years.xlsx") %>% data.frame(stringsAsFactors = FALSE)
price <- sp_dat$Spx_index[13337:26701]

# 로그 수익률 계산
log_price <- log(price)
log_return <- diff(log_price, differences = 1)

for (i in 2:length(log_return)) {
  if (log_return[i] == 0) {
    log_return[i] <- log_return[i - 1]
  }
}

# simulation setting
length_period <- 10 
num_periods <- floor((length(log_return)-750)/length_period)
num_sim_multi_r <- 10^4
confidence_level <- 0.975
delta <- 0.25 
L <- 10^4

# 처음 750일을 제외한 실제 10일간의 수익률 계산
actual_logr10 = rep(0,num_periods) 
for (i in 1:num_periods) {
  start_ind = length_period*(i-1) + 751
  end_ind = start_ind + length_period 
  actual_logr10[i] = log_price[end_ind] - log_price[start_ind]
}

# GjrGARCH 모형 정의
spec.gjrGARCH <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), 
                            mean.model = list(armaOrder = c(0, 0), include.mean=FALSE), 
                            distribution.model = "std")

# 다기간 가상 수익률을 생성하는 함수 정의
simulated_return_norm <- function(par, length_period){
  sam_z <- rt(length_period, par[7])
  sam_z<- sam_z/sqrt(par[7]/(par[7]-2))
  sim_r <- rep(0, length_period)
  
  old_r <- par[1]
  old_sigma <- par[2]
  
  for(j in 1:length_period){
    sigma_sq <- par[3] + ((old_r < 0) * par[6] + par[4]) * old_r^2 + par[5] * old_sigma^2
    new_sigma <- sqrt(sigma_sq)
    new_r <- new_sigma * sam_z[j]
    sim_r[j] <- new_r
    
    old_r <- new_r
    old_sigma <- new_sigma
  }
  sum(sim_r)
}

# VaR과 Standard error 
compute_VaR_ES_cmc <- function(returns, confidence_level) {
  n <- length(returns)
  m <- floor(n/10)
  sub_VaR <- rep(0,10)
  sub_ES <- rep(0,10)
  for (k in 1:10) {
    sub_loss <- -returns[(1 + (k-1)*m):(k*m)]
    sub_VaR[k] <- quantile(sub_loss, confidence_level)
  }
  VaR <- mean(sub_VaR)
  se_VaR <- sd(sub_VaR)/sqrt(10)
  for (k in 1:10) {
    sub_loss <- -returns[(1 + (k-1)*m):(k*m)]
    sub_ES[k] <- sum(sub_loss*(sub_loss > VaR)) / sum(sub_loss > VaR)
  }
  ES <- mean(sub_ES)
  se_ES <- sd(sub_ES)/sqrt(10)
  c(VaR, se_VaR, ES, se_ES)
}

# setting for parallel computing 
numCores <- detectCores() - 1
cl <- makeCluster(numCores, type = "PSOCK")
clusterExport(cl, c("simulated_return_norm", "length_period"))

# 다기간 VaR 및 ES 계산
result <- matrix(rep(0,num_periods*4),num_row <- num_periods)
time_start <- Sys.time()
for(i in 1:num_periods){
  start_ind <- length_period * (i-1) + 1
  end_ind <- start_ind + 749 
  sub_vec <- log_return[start_ind:end_ind]
  
  # GjrGARCH모형 적합 
  gjrGARCH <- ugarchfit(data = sub_vec, spec = spec.gjrGARCH, solver = "hybrid")
  coef_value <- gjrGARCH@fit$coef
  
  # 모수 결정
  par <- c(sub_vec[length(sub_vec)], gjrGARCH@fit$sigma[length(gjrGARCH@fit$sigma)], coef_value)
  clusterExport(cl, c("par"))
  
  # 다기간 수익률 중요도 샘플 추출
  sim_multi_return <-  parSapply(cl, 1:num_sim_multi_r, function(x){simulated_return_norm(par, length_period)})
  
  # VaR과 ES계산 
  result[i,] <- compute_VaR_ES_cmc(sim_multi_return,confidence_level)
  cat("VaR_fhs[",i,"] =",result[i,1],"se_VaR[",i,"] =",result[i,2],"\n")
  cat(" ES_fhs[",i,"] =",result[i,3]," se_ES[",i,"] =",result[i,4],"\n")
}
time_end <- Sys.time()

# 총 걸린 시간
total_time <- time_end - time_start
cat("time:", total_time, "\n")
print(difftime(time_end,time_start,units="secs"))

stopCluster(cl)

# VaR과 ES의 standard error, relative error
VaR_cmc <- result[,1]
se_VaR <- result[,2]
ES_cmc <- result[,3]
se_ES <- result[,4]

cat("mean of s.e. for VaR estimation:",mean(se_VaR),"\n")
cat("mean of r.e. for VaR estimation:",mean(se_VaR/VaR_cmc),"\n")
cat("mean of s.e. for ES estimation:",mean(se_ES),"\n")
cat("mean of r.e. for ES estimation:",mean(se_ES/ES_cmc),"\n")

# 초과손실 계산
exceptions_cmc <- (-actual_logr10 > VaR_cmc)
cat("#### number of violations ####","\n")
cat("expected:", num_periods*(1-confidence_level), "\n")
cat("cmc with t errors:", sum(exceptions_cmc),"\n")

plot(-VaR_cmc, type="l")
lines(actual_logr10,type="p")

E_cmc_0.975 <- cbind(actual_logr10, ES_cmc, se_ES)
write.csv(E_cmc_0.975, "C:/Users/seoyeji/Desktop/E_cmc_0.975.csv")