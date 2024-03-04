## A simulation study with TVC

save_ch6 <- matrix(0, nrow = 200, ncol = 12)

for(s in 1:200){
  data <- tvc_simulation(2000, 0.5, -0.2, 0.3, 0.1, 0.7, 3.5)
  dat_baseline <- data$dat_baseline
  dat_long <- data$dat_long
  
  #fit using mpl 
  dat_long$left_position <- dat_long$tau = rep(0, nrow(dat_long))
  dat_left <- dat_long
  dat_left[which(dat_left$tau !=1),-1] <- 0
  ctrl <- tvc_mpl_control(lambda = 0, iter = c(1,2000), n_knots = 10, par_initial = c(0,0,0.1), 
                          range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
  test <- tvc_fit((dat_long), (dat_baseline), (dat_left), ctrl)
  
  save_ch6[s,1] <- test$parameters$beta[1]
  save_ch6[s,2] <- test$parameters$beta[2]
  save_ch6[s,3] <- test$parameters$gamma[1]
  
  save_ch6[s,4] <- test$se_Q[1]
  save_ch6[s,5] <- test$se_Q[2]
  save_ch6[s,6] <- test$se_H[3]
  
  #fit using pl
  dat_long$event_long <- dat_long$delta_long * dat_long$last_record
  
  surv.obj <- Surv(time = dat_long$start, time2 = dat_long$end, 
                   type = "counting", event = dat_long$event_long)
  pl.fit <- coxph(surv.obj ~ dat_long$x1_long + dat_long$x2_long + 
                    dat_long$Z, data = dat_long, id = dat_long$i_long )
  
  
  save_ch6[s,7] <- pl.fit$coefficients[1]
  save_ch6[s,8] <- pl.fit$coefficients[2]
  save_ch6[s,9] <- pl.fit$coefficients[3]
  
  save_ch6[s,10] <- sqrt(diag(pl.fit$var))[1]
  save_ch6[s,11] <- sqrt(diag(pl.fit$var))[2]
  save_ch6[s,12] <- sqrt(diag(pl.fit$var))[3]
  
}
