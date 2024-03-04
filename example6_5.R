# Breast cancer example

# Read in long format and baseline datasets available on Github repo
bc_longformat <- read.table("bc_longformat.txt", header = TRUE, row.names = NULL)
bc_baseline <- read.table("bc_baseline.txt", header = TRUE, row.names = NULL)

# Need to create a dataset that is cut-off at the left end of any interval censoring intervals

# Note that creating bc_left requires the bc_longformat dataset 
# to have the columns tau and left_position

bc_longformat$left_position <- bc_longformat$tau <- rep(0, nrow(bc_longformat))

for(i in 1:length(unique(bc_longformat$i_long))){
  
  #identify rows in bc_longformat associated with individual i
  ind_long <- which(bc_longformat$i_long == i)
  
  if(bc_longformat$delta_long[ind_long[1]] == 3){ #if individual i has delta == 3
    if(length(ind_long) == 2){ #if z(t) changed for individual i during follow-up
      if(as.numeric(bc_longformat$start)[ind_long[2]] < as.numeric(bc_longformat$TL_long)[ind_long[1]]){
        #if the change time of z(t) for individual i was before t_i^L
        bc_longformat$left_position[ind_long[2]] <- 1
        bc_longformat$tau[ind_long] <- 1
      }else{
        #if the change time of z(t) for individual i was after t_i^L
        bc_longformat$left_position[ind_long[1]] <- 1
        bc_longformat$tau[ind_long[1]] <- 1
      }
      
    }else if(length(ind_long == 1)){ #if z(t) did not change for individual i during follow-up
      bc_longformat$left_position[ind_long] <- 1
      bc_longformat$tau[ind_long] <- 1
    }
  }
}

bc_left <- bc_longformat
bc_left[which(bc_left$tau !=1),-1] <- 0

left.i <- which(bc_baseline$delta==3)
for(li in left.i){
  #replace "end" time with left interval time
  #make everything after left interval time 0
  ind <- which(bc_longformat$i_long == li)
  left.pos.ind <- ind[which(bc_left$left_position[ind] == 1)]
  bc_left$end[left.pos.ind] <- bc_left$TL[left.pos.ind] 
}


# Now we can fit the model

ctrl <- tvc_mpl_control(lambda = 0, iter = c(10,3000), n_knots = 6, par_initial = c(0,0,1), range = c(0.1,0.9), line_search = c(1, 1, 1), reg_conv = 1e-5)
breastcancer_mpl <- tvc_fit(bc_longformat, bc_baseline, bc_left, ctrl)
