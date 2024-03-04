# Revisit the Stanford heart data using MPL method

library(survival)
data("heart")
# We need to do some data processing to use the MPL package
n <- 103

# Create baseline data set
stnfrd_hrt_baseline <- NULL

for(i in 1:n){
  ind_long <- which(heart$id == i)
  new_row <- c(heart$id[ind_long[1]], heart$stop[ind_long[length(ind_long)]], 
               ifelse(heart$event[ind_long[length(ind_long)]] == 1, heart$stop[ind_long[length(ind_long)]], Inf), 
               heart$event[ind_long[length(ind_long)]], 
               as.numeric(heart$age)[ind_long[1]], as.numeric(heart$year)[ind_long[1]],
               as.numeric(heart$surgery)[ind_long[1]], length(ind_long))
  stnfrd_hrt_baseline <- rbind(stnfrd_hrt_baseline, new_row)
}

row.names(stnfrd_hrt_baseline) <- NULL
stnfrd_hrt_baseline <- as.data.frame(stnfrd_hrt_baseline)
colnames(stnfrd_hrt_baseline) <- c("i", "TL", "TR", "delta", "age", "year", "surgery", "rep")
table(stnfrd_hrt_baseline$rep)

# Create long data set
stnfrd_hrt_long <- data.frame(i_long = heart$id, 
                              TL_long = rep(stnfrd_hrt_baseline$TL, stnfrd_hrt_baseline$rep), 
                              TR_long = rep(stnfrd_hrt_baseline$TR, stnfrd_hrt_baseline$rep), 
                              delta_long = rep(stnfrd_hrt_baseline$delta, stnfrd_hrt_baseline$rep), 
                              age_long = rep(stnfrd_hrt_baseline$age, stnfrd_hrt_baseline$rep), 
                              year_long = rep(stnfrd_hrt_baseline$year, stnfrd_hrt_baseline$rep), 
                              surgery_long = rep(stnfrd_hrt_baseline$surgery, stnfrd_hrt_baseline$rep), 
                              transplant = as.numeric(heart$transplant == 1),
                              start = heart$start, end = heart$stop)

stnfrd_hrt_long$last_record <- 0
for(i in 1:n){
  ind_long <- which(stnfrd_hrt_long$i_long == i)
  stnfrd_hrt_long$last_record[ind_long[length(ind_long)]] <- 1
}

stnfrd_hrt_long$left_position <- stnfrd_hrt_long$tau <- rep(0, nrow(stnfrd_hrt_long))
stnfrd_hrt_left <- stnfrd_hrt_long
stnfrd_hrt_left[which(stnfrd_hrt_left$tau !=1),-1] <- 0

# Fit the MPL model
library(splines2)
library(dplyr)
ctrl <- tvc_mpl_control(lambda = 0, iter = c(10,1000), n_knots = 4, par_initial = c(0,0,0.1), 
                        range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-6)
test <- tvc_fit((stnfrd_hrt_long), (stnfrd_hrt_baseline), (stnfrd_hrt_left), ctrl)

# Plots of baseline survival quantities estimates
plot(test$func_est$h0_est ~ test$func_est$v, type = "l")
plot(test$func_est$H0_est ~ test$func_est$v, type = "l")
plot(test$func_est$S0_est ~ test$func_est$v, type = "l")

# Predictive survival plots
# year = age = 0, no surgery vs. year = age = 0, surgery at t = 20
St_surg0 <- LL <- UL <- NULL

for(t in 1:1000){
  Ht <- (test$func_est$H0_est[t])
  St <- exp(-Ht)
  Psi <- Psi_f(events = test$func_est$v[t], knts = test$kn$int_knots, Boundary.knots = test$kn$bound_knots)
  
  dSt_beta <- -St*Ht*as.vector(rep(0,3))
  dSt_gamma <- -St*Ht*as.vector(0)
  dSt_theta <- -St * Psi
  dSt_eta <- matrix(c(dSt_beta,dSt_gamma, dSt_theta))
  cov <- test$covar_H
  
  se <- sqrt(t(dSt_eta) %*% cov %*% (dSt_eta))
  St_surg0 <- c(St_surg0, St)
  LL <- c(LL, (St - 1.96*se))
  UL <- c(UL, (St + 1.96*se))
  
}

plot(test$func_est$v, St_surg0, type = "l", ylim = c(0,1), main = "Survival probability by transplant: none",
     ylab = "Estimated survival probability", xlab = "Time (days)")
polygon(x = c(test$func_est$v, rev(test$func_est$v)), y = c(LL, rev(UL)), 
        border = NA, col = adjustcolor("grey60", alpha.f=0.5) )

lines(test$func_est$v, St_surg0, lwd = 2)

St_surg_t10 <- LL1 <- UL1 <- NULL

for(t in 1:1000){
  theta <- test$parameters$theta
  if(test$func_est$v[t] < 10){
    Ht <- (test$func_est$H0_est[t])
    St <- exp(-Ht)
    Psi <- Psi_f(events = test$func_est$v[t], knts = test$kn$int_knots, Boundary.knots = test$kn$bound_knots)
    
    dSt_beta <- -St*Ht*as.vector(rep(0,3))
    dSt_gamma <- -St*Ht*as.vector(0)
    dSt_theta <- -St * Psi
    dSt_eta <- matrix(c(dSt_beta,dSt_gamma, dSt_theta))
    cov <- test$covar_H
    se <- sqrt(t(dSt_eta) %*% cov %*% (dSt_eta))
    
    St_surg_t10 <- c(St_surg_t10, St)
    LL1 <- c(LL1, (St - 1.96*se))
    UL1 <- c(UL1, (St + 1.96*se))
  }else{
    tMinus1 <- min(which(test$func_est$v>=10))
    
    Psi100 <- Psi_f(events = test$func_est$v[tMinus1], knts = test$kn$int_knots, Boundary.knots = test$kn$bound_knots)
    Ht100 <- Psi100 %*% theta
    Ht100_g <- Psi100 %*% theta * exp(test$parameters$gamma[1])
    
    Ht <- Ht100 + (test$func_est$H0_est[t]*exp(test$parameters$gamma[1])) - Ht100_g
    St <- exp(-Ht)
    Psi <- Psi_f(events = test$func_est$v[t], knts = test$kn$int_knots, Boundary.knots = test$kn$bound_knots)
    
    dSt_beta <- -St*Ht*as.vector(rep(0,3))
    dSt_gamma <- -St*(test$func_est$H0_est[t]*exp(test$parameters$gamma[1])) + St*Ht100_g
    dSt_gamma <- as.vector(c(dSt_gamma))
    dSt_theta <- as.vector(-St * exp(test$parameters$gamma[1]) )* Psi -
      as.vector(St)*Psi100 + 
      as.vector(St*exp(test$parameters$gamma[1]))*Psi100
    dSt_eta <- matrix(c(dSt_beta,dSt_gamma,dSt_theta))
    cov <- test$covar_H
    se <- sqrt(t(dSt_eta) %*% cov %*% (dSt_eta))
    
    St_surg_t10 <- c(St_surg_t10, St)
    LL1 <- c(LL1, (St - 1.96*se))
    UL1 <- c(UL1, (St + 1.96*se))
  }
}

LL1[which(LL1<0)] = 0

plot(test$func_est$v, St_surg_t10, type = "l", ylim = c(0,1), main = "Survival probability by transplant: @ t=10",
     ylab = "Estimated survival probability", xlab = "Time (days)")
polygon(x = c(test$func_est$v, rev(test$func_est$v)), y = c(LL1, rev(UL1)), 
        border = NA, col = adjustcolor("grey60", alpha.f=0.5) )
lines(test$func_est$v, St_surg_t10, lwd =2)