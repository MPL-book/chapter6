## Stanford heart transplant

library(survival)
data("heart")

surv.obj <- Surv(time = heart$start, time2 = heart$stop, 
                 type = "counting", event = heart$event)
pl.fit <- coxph(surv.obj ~ heart$age + heart$year + heart$surgery + 
                  heart$transplant, data = heart, id = heart$id )
summary(pl.fit)
