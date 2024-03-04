# Generate extended Cox model data with right censoring

# generate data
tvc_simulation <- function(n, beta1, beta2, gamma, 
                           t_0i_par1, t_0i_par2, c_i_par){
  
  # draw time-fixed covariate values and regression term
  x1 <- rnorm(n)
  x2 <- rbinom(n, 1, 0.5)
  mu_term <- exp(beta1 * x1 + beta2 * x2)
  
  # draw t_0i values
  t_0i <- runif(n, t_0i_par1, t_0i_par2)
  
  # fixed lambda and nu
  lambda <- 1
  nu <- 1.5
  
  # draw -log(u)
  neg_log_u = -log(runif(n))
  
  # draw event times
  t <- NULL
  for(i in 1:n){
    if(neg_log_u[i] < lambda * mu_term[i] * t_0i[i]^nu){
      t_i <- (neg_log_u[i]/( mu_term[i]))^(1/nu)
      t <- c(t, t_i)
    }else{
      t_i <- ((neg_log_u[i] - mu_term[i] * t_0i[i]^nu 
               + (mu_term[i] * exp(gamma)) * t_0i[i]^nu)/
                (mu_term[i] * exp(gamma)))^(1/nu)
      t <- c(t, t_i)
    }
  }
  
  c_i <- runif(n, 0, c_i_par)
  delta <- as.numeric(t < c_i)
  TL <- TR <- rep(0,n)
  TL[which(delta == 1)] <- t[which(delta == 1)]
  TR[which(delta == 1)] <- t[which(delta == 1)]
  TL[which(delta == 0)] <- c_i[which(delta == 0)]
  TR[which(delta == 0)] <- Inf
  
  # create "long" dataset
  i_long <- TL_long <- TR_long <- delta_long <- x1_long <- x2_long <- Z <- start <- end <- last_record <- NULL
  
  for(i in 1:n){
    if(TL[i] < t_0i[i]){
      # no change in z observed
      i_long <- c(i_long, i)
      TL_long <- c(TL_long, TL[i])
      TR_long <- c(TR_long, TR[i])
      delta_long <- c(delta_long, delta[i])
      x1_long <- c(x1_long, x1[i])
      x2_long <- c(x2_long, x2[i])
      Z <- c(Z, 0)
      start <- c(start, 0)
      end <- c(end, TL[i])
      last_record <- c(last_record, 1)
    }else{
      # change in z observed
      i_long <- c(i_long, rep(i,2))
      TL_long <- c(TL_long, rep(TL[i],2))
      TR_long <- c(TR_long, rep(TR[i],2))
      delta_long <- c(delta_long, rep(delta[i],2))
      x1_long <- c(x1_long, rep(x1[i],2))
      x2_long <- c(x2_long, rep(x2[i],2))
      Z <- c(Z, 0, 1)
      start <- c(start, 0, t_0i[i])
      end <- c(end, t_0i[i], TL[i])
      last_record <- c(last_record, 0, 1)
      
    }
  }
  
  dat_long <- data.frame(i_long, TL_long, TR_long, delta_long, 
                         x1_long, x2_long, Z, start, end, last_record)
  
  # create baseline dataset
  dat_baseline <- data.frame(i = c(1:n), TL = TL, TR = TR, 
                             delta = delta, x1 = x1, x2 = x2, rep = c(table(i_long)))
  
  out = list(dat_baseline = dat_baseline, dat_long = dat_long)
  return(out)
  
}

# fit using MPL method 
# access code for tvc_fit() function at Github repo tvc_mpl
source("fitting_diffN.R")
source("functions_diffN.R")

data <- tvc_simulation(1000, 1, 1, 1, 0.1, 0.5, 1.05)
dat_baseline <- data$dat_baseline
dat_long <- data$dat_long
dat_long$left_position <- dat_long$tau <- rep(0, nrow(dat_long))
dat_left <- dat_long
dat_left[which(dat_left$tau !=1),-1] <- 0

table(dat_baseline$delta)
table(dat_baseline$rep)

ctrl <- tvc_mpl_control(lambda = 0, iter = c(10,1000), n_knots = 7, par_initial = c(0,0,0.1), 
                        range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
test <- tvc_fit((dat_long), (dat_baseline), (dat_left), ctrl)

