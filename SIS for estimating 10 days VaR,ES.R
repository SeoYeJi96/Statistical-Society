rm(list = ls())
require(dplyr)
require(readxl)
require(rugarch)
library(parallel)

################################################################################################
#  1990.01.01~2023.12.31 10-days Value at Risk
###############################################################################################

sp_dat = read_xlsx("./spx_100years.xlsx") %>% data.frame(stringsAsFactors = F) 
price <- sp_dat$Spx_index[13337:26701]

# 로그 수익률 계산
log_price = log(price)
log_return = diff(log_price,  differences = 1)

# 0인 로그 수익률을 이전 로그 수익률로 대체
for (i in 2:length(log_return)) {
  if (log_return[i] == 0) {
    log_return[i] <- log_return[i - 1]
  }
}

# simulation setting
length_period <- 10 
num_periods <- floor((length(log_return)-750)/length_period)
num_sim_multi_r <- 10^4
confidence_level <- 0.99
delta <- 0.25
#L <- 5*10^3

# 처음 750일을 제외한 실제 10일간의 수익률 계산
actual_logr10 = rep(0,num_periods) 
for (i in 1:num_periods) {
  start_ind = length_period*(i-1) + 751
  end_ind = start_ind + length_period 
  actual_logr10[i] = log_price[end_ind] - log_price[start_ind]
}

spec.gjrGARCH = ugarchspec(variance.model = list(model="gjrGARCH", garchOrder = c(1, 1)), 
                           mean.model = list(armaOrder = c(0, 0),include.mean=FALSE), distribution.model = "norm")

# 다기간 가상 수익률을 생성하는 함수 정의:  SIS method
simulated_return_sis <- function(par, length_period, lambda_hat, delta) {
  index <- sample(1:length(z_vec), size = length_period, replace = TRUE, prob = twisted_pmf)
  norm_z <- rnorm(length_period, mean = z_vec[index] + lambda_hat * delta^2, sd = delta)
  norm_r <- rep(0, length_period)
  
  # 초기 수익률과 표준편차 지정 
  old_r <- par[1]
  old_sigma <- par[2]
  
  for (j in 1:length_period) {
    sigma_sq <- par[3] + ((old_r < 0) * par[6] + par[4]) * old_r^2 + par[5] * old_sigma^2
    new_sigma <- sqrt(sigma_sq)
    new_r <- new_sigma * norm_z[j]
    norm_r[j] <- new_r  
    # 수익률과 표준편차 업데이트      
    old_r <- new_r
    old_sigma <- new_sigma
  }
  c(sum(norm_r), exp(-lambda_hat * sum(norm_z)))
}

# 목적함수를 최소화하는 lambda 값 찾기 
find_lambda <- function(par, var_q, z_vec, lambda) {
  diff <- 1
  step <- 1
  
  while (diff > 0.001) {
    clusterExport(cl, c("lambda","lambda_hat"))
    
    sim_result <- parSapply(cl, 1:1000, function(x) {
      simulated_return_sis(par, length_period, lambda_hat, delta)
    })
    
    r <- sim_result[1,]
    R <- sim_result[1, r < -var_q]
    weight <- sim_result[2, r < -var_q]
    weight <- weight / sum(weight)
    sum_Z <- -log(weight) / lambda
    
    if (length(R) > 0) {
      c_lambda <- mean(exp(lambda_hat * z_vec + delta^2 * lambda_hat^2 / 2))
      derivative_c_lambda <- mean((z_vec + delta^2 * lambda_hat) * exp(lambda_hat * z_vec + delta^2 * lambda_hat^2 / 2))
      derivative_log_g <- sum_Z - length_period * derivative_c_lambda / c_lambda
      gradient <- mean(R * weight * derivative_log_g)
      lambda_next <- lambda_hat - 100 * gradient / (step + 10)
      
      cat("step:", step, " lambda:", lambda_next, "\n")
      diff <- abs(lambda_next - lambda_hat)
      step <- step + 1
      
      # lambda_hat 값을 업데이트하여 동일한 분포에서 샘플링
      lambda_hat <- lambda
    }
    lambda <- lambda_next
  }
  
  return(lambda)
}


# VaR과 Standard error 
compute_VaR_ES <- function(returns, weight, confidence_level) {
  n <- length(returns)
  m <- floor(n/10)
  sub_VaR <- rep(0,10)
  sub_ES <- rep(0,10)
  for (k in 1:10) {
    sub_return <- returns[(1 + (k-1)*m):(k*m)]
    sub_weight <- weight[(1 + (k-1)*m):(k*m)]/sum(weight[(1 + (k-1)*m):(k*m)])
    ordered_return <- sort(sub_return)
    ordered_weight <- sub_weight[order(sub_return)]
    
    j <- 1
    sum_weight <- ordered_weight[j]
    while(sum_weight < 1 - confidence_level){j <- j+1;sum_weight <- sum_weight + ordered_weight[j]}
    sub_VaR[k] <- -(ordered_return[j-1] + ordered_return[j])/2
  }		
  VaR <- mean(sub_VaR)		
  se_VaR <- sd(sub_VaR)/sqrt(10)
  
  for (k in 1:10) {
    sub_return <- returns[(1 + (k-1)*m):(k*m)]
    sub_weight <- weight[(1 + (k-1)*m):(k*m)]/sum(weight[(1 + (k-1)*m):(k*m)])
    ind_excess <- (sub_return < -VaR) 
    sub_ES[k] <- - sum(sub_return*sub_weight*ind_excess)/sum(sub_weight*ind_excess)
  }		
  ES <- mean(sub_ES)
  se_ES <- sd(sub_ES)/sqrt(10)
  
  c(VaR, se_VaR, ES, se_ES)
}


# setting for parallel computing 
numCores =  detectCores() - 1
cl = makeCluster(numCores,type = "PSOCK")
clusterExport(cl,c("simulated_return_sis","length_period","delta"))
clusterExport(cl,c("length_period","delta"))

# 다기간 VaR 계산
result <- matrix(rep(0,num_periods*4),num_row <- num_periods)
time_start <- Sys.time()
lambda <- -0.8
var_q <- 0.1
for(i in 1:num_periods){
  start_ind = length_period*(i-1) + 1
  end_ind = start_ind + 749 
  sub_vec = log_return[start_ind:end_ind]
  
  # GjrGARCH모형 적합
  gjrGARCH <- ugarchfit(data = sub_vec, spec = spec.gjrGARCH, solver = "hybrid")
  coef_value = gjrGARCH@fit$coef
  
  # 모수 결정
  par <- c(sub_vec[length(sub_vec)], gjrGARCH@fit$sigma[length(gjrGARCH@fit$sigma)],coef_value)
  clusterExport(cl,c("par"))
  
  ## 다기간 VaR 추정: SIS method
  # 표준화된 잔차 계산와 중요도 잔차 분포 정의
  z_vec = sub_vec/gjrGARCH@fit$sigma
  z_vec <- (z_vec > -4)*z_vec - 4*(z_vec < -4)	
  twisted_pmf <- exp(lambda*z_vec)/sum(exp(lambda*z_vec))
  clusterExport(cl,c("z_vec","lambda","twisted_pmf"))
  
  # 교차 엔트로피 최소방법을 이용한 최적 lambda의 추정
  lambda_hat <- lambda
  lambda <- find_lambda(par,var_q,z_vec,lambda_hat)

  # 다기간 수익률 중요도 샘플 추출 
  sim_result <- parSapply(cl,1:num_sim_multi_r, function(x){simulated_return_sis(par,length_period,lambda_hat,delta)})
  sim_multi_return <- sim_result[1,]
  sim_weight <- sim_result[2,]
  
  # VaR과 ES계산 
  result[i,] <- compute_VaR_ES(sim_multi_return,sim_weight,confidence_level)
  var_q <- result[i,1]
  cat("lambda =",lambda,"\n")
  cat("VaR_fhs[",i,"] =",result[i,1],"se_VaR[",i,"] =",result[i,2],"\n")
  cat(" ES_fhs[",i,"] =",result[i,3]," se_ES[",i,"] =",result[i,4],"\n")
}
time_end <- Sys.time()

# 총 걸린 시간
total_time <- time_end - time_start
cat("time:", total_time, "\n")
stopCluster(cl)

VaR_sis <- result[,1]
se_VaR_sis <- result[,2]
ES_sis <- result[,3]
se_ES_sis <- result[,4]

cat("mean of s.e. for VaR estimation:",mean(se_VaR_sis),"\n")
cat("mean of r.e. for VaR estimation:",mean(se_VaR_sis/VaR_sis),"\n")
cat("mean of s.e. for ES estimation:",mean(se_ES_sis),"\n")
cat("mean of r.e. for ES estimation:",mean(se_ES_sis/ES_sis),"\n")

# 초과손실 계산
exceptions_fhs <- (-actual_logr10 > VaR_sis)
cat("#### number of violations ####","\n")
cat("expected:", num_periods*(1-confidence_level), "\n")
cat("fhs method:", sum(exceptions_fhs,na.rm=TRUE),"\n")

plot(-VaR_sis, type="l")
lines(actual_logr10,type="p")

