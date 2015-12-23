# run these two packages first
library(doSNOW)
library(foreach)

################################################
# hh matrix  for household level exposures ####
hhmatrix_hhexp <- function(boot.iter,
                     Data=q1q4iccctrl,
                     Varlist = c("pplrm_di","defecate_pitl","defsharing_only",
                                     "dist_more10", "q4water2","q11_tube",
                                     "q11_supply","boilnewst","soap")){
  data <- Data[c(Varlist,"dataid","cluster")]   # subsetted dataset
  # all combinations of risk factors
  mcomb<-as.data.frame(combn(1:length(Varlist),2))
  # 
  if (boot.iter!=0){
    # bootstrap on cluster level
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.cluster,funcfunc)
  } else {
    bootlist <- list(data$dataid)
  }
  
  # calculate cooccurrence of rf 
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  prob.hh.boot<- foreach(b = 1:boot.iter,.combine='rbind') %:% 
    foreach (k = 1:ncol(mcomb),.combine='c') %dopar% { # for each combination of variables
      q1q4iccctrl.matrix.b<-data[match(bootlist[[b]],data$dataid),]
      i<-mcomb[1,k] #ith risk factor
      j<-mcomb[2,k] #jth risk factor
      #  subset a dataframe that include: rf1 rf2 dataid cluster
      temp0 <- q1q4iccctrl.matrix.b[ c(i,j) ]    
      # delete na
      temp <- na.omit(temp0)
      # hh concordant: prob = pairs being concordant over total number of pairs
      if (nrow(temp)!=0) {
        logiccount <- temp[1]==temp[2]
        junko <- length(logiccount[logiccount==TRUE]) / nrow(temp)  
        junko 
      }  else {
        NA
      }
    }
  stopCluster(cl)
  prob.hh.boot
}

# cluster matrix for household level exposures ####
clustermatrix_hhexp <- function(boot.iter,
                     Data=q1q4iccctrl,
                     Varlist = c("pplrm_di","defecate_pitl","defsharing_only",
                                 "dist_more10", "q4water2","q11_tube",
                                 "q11_supply","boilnewst","soap")){
  data <- Data[c(Varlist,"dataid","cluster")]   # subsetted dataset
  # all combinations of risk factors
  mcomb<-as.data.frame(combn(1:length(Varlist),2))
  if (boot.iter!=0){
    # bootstrap on cluster level
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.cluster,funcfunc)
  } else {
    bootlist <- data$dataid
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
  }
  # get co-occurrence of exposures
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  prob.cluster.boot <- foreach (b = 1:boot.iter,.combine='rbind') %:%
    foreach (k = 1:ncol(mcomb),.combine='c') %dopar% { # for each combination of variables
      q1q4iccctrl.matrix.b<-data[match(bootlist[[b]], data$dataid),]
      q1q4iccctrl.matrix.b$newcluster <- names(bootlist[[b]])
      i<-mcomb[1,k] #ith risk factor
      j<-mcomb[2,k] #jth risk factor
      #  subset a dataframe that include: rf1 rf2 dataid cluster
      temp0 <- q1q4iccctrl.matrix.b[ c(i,j,
                                       which(colnames(q1q4iccctrl.matrix.b)=="dataid"),
                                       which(colnames(q1q4iccctrl.matrix.b)=="newcluster")) ]    
      temp <- na.omit(temp0)
      fenzi <-0
      fenmu <-0
      for (p in 1:nrow(temp)) { # for each household
        seedhh <- temp[p,]
        # subset to its cluster
        temp_cluster <- temp[which(temp$newcluster==seedhh[1,4]),]
        # number of concordant pairs - self...only count TRUE
        if (nrow(temp_cluster)!=0) {
          logiccount <- seedhh[1,1]==temp_cluster[2]
          fenzi <- fenzi + length(logiccount[logiccount==TRUE])
          # number of pairs without any NA (both with numbers) - self...count TRUE and FALSE..exclude NA
          fenmu <- fenmu + nrow(temp_cluster) 
        } else {
          fenzi <- fenzi+0
          fenmu <- fenmu+0
        }
      }
      fenzi/fenmu
    }
  stopCluster(cl)
  prob.cluster.boot
}

# overall matrix for household level exposures ####
overallmatrix_hhexp <- function(boot.iter,
                     Data=q1q4iccctrl,
                     Varlist = c("pplrm_di","defecate_pitl","defsharing_only",
                                 "dist_more10", "q4water2","q11_tube",
                                 "q11_supply","boilnewst","soap")){
  data <- Data[c(Varlist,"dataid","cluster")]   # subsetted dataset
  # all combinations of risk factors
  mcomb<-as.data.frame(combn(1:length(Varlist),2))
  if (boot.iter !=0){
    # bootstrap on cluster level
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.cluster,funcfunc)
  } else {
    bootlist <- list(data$dataid)
  }
  # get co-occurrence of rf
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  prob.overall.boot <- foreach (b = 1:boot.iter,.combine='rbind') %:% # bootstrap households 
    foreach (k = 1:ncol(mcomb),.combine='c') %dopar% { # for each combination of variables..
      q1q4iccctrl.matrix.b<-data[match(bootlist[[b]], data$dataid),]
      i<-mcomb[1,k] #ith risk factor
      j<-mcomb[2,k] #jth risk factor
      #  subset a dataframe that include: rf1 rf2 
      temp0 <- q1q4iccctrl.matrix.b[ c(i,j)]    
      temp <- na.omit(temp0)
      fenzi <- 0 
      fenmu <- 0
      for (p in 1:nrow(temp)) { # for each hh
        if (nrow(temp)!=0) {
          seedhh <- temp[p,]
          # number of concordant pairs - self...only count TRUE
          logiccount <- seedhh[1,1]==temp[2]
          fenzi <- fenzi + length(logiccount[logiccount==TRUE])
          fenmu <- fenmu + nrow(temp) # number of pairs without any NA (both with numbers)
        } else {
          fenzi <- fenzi + 0
          fenmu <- fenmu + 0
        }
      }
      fenzi/fenmu # after all HHs' been looped through
    }
  stopCluster(cl)
  prob.overall.boot
}


################################################
# ind matrix for individual level exposures ####
indmatrix_indexp <- function(boot.iter,
                             Data=q2icc_nohh,
                             Varlist=c("handfeed_never","prep2hrs_never",
                                       "q40_bc_never","q41_4_bc_never","q41_5_bc_never")){
  ##
  data <- Data[c(Varlist,"datamemberid","dataid","cluster") ]
  # create mcomb
  mcomb <- as.data.frame(combn(1:length(Varlist),2)) 
  ## 
  if (boot.iter!=0) {
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on matched-set level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist.hh <- lapply(bootlist.cluster,funcfunc)
    # bootstrap on household level ...names...new cluster
    funcfunc.ind  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(names(x)[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.hh,funcfunc.ind)  
    # bootstrap on individual level .....names... new dataid
    funcfunc.ind.n  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(x[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist.n <- lapply(bootlist.hh,funcfunc.ind.n)  
  } else {
    bootlist.n <- list(data$datamemberid)
  }
  
  # get co-occurrence of exposures
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  foreach (b = 1:boot.iter,.combine='rbind')  %:% 
    foreach (k = 1:ncol(mcomb),.combine='c') %do% { # for each combination of variables, 231 in total
      q2icc_nohh.matrix.b <- data[match(bootlist.n[[b]], data$datamemberid),]
      
      i<-mcomb[1,k] #ith risk factor
      j<-mcomb[2,k] #jth risk factor
      #  subset a dataframe that include: rf1 rf2 datamemberid
      temp0 <- q2icc_nohh.matrix.b[ c(i,j,which(colnames(q2icc_nohh.matrix.b)=="datamemberid")) ]    
      # delete na
      temp <- na.omit(temp0)
      # hh concordant: prob = pairs being concordant over total number of pairs
      if (nrow(temp)!=0) {
        logiccount <- temp[1]==temp[2]
        junko <- length(logiccount[logiccount==TRUE]) / nrow(temp)
        junko
      } else {
        NA
      }
    }
}
  

# hh matrix for individual level exposures ####
hhmatrix_indexp <- function(boot.iter,
                            Data=q2icc_nohh,
                            Varlist=c("handfeed_never","prep2hrs_never",
                                      "q40_bc_never","q41_4_bc_never","q41_5_bc_never")){
  ##
  data <- Data[c(Varlist,"datamemberid","dataid","cluster") ]
  # create mcomb
  mcomb <- as.data.frame(combn(1:length(Varlist),2)) 
  if (boot.iter!=0) {
    ## 
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on matched-set level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist.hh <- lapply(bootlist.cluster,funcfunc)
    # bootstrap on household level ...names...new cluster
    funcfunc.ind  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(names(x)[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.hh,funcfunc.ind)  
    # bootstrap on individual level .....names... new dataid
    funcfunc.ind.n  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(x[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist.n <- lapply(bootlist.hh,funcfunc.ind.n)  
  } else {
    bootlist.n <- data$datamemberid
    names(bootlist.n) <- substr(bootlist.n,1,5)
    bootlist.n <- list(bootlist.n)
    bootlist <- data$datamemberid
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
  }
  
  # get co-occurrence of exposures
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  foreach (b = 1:boot.iter,.combine='rbind') %:% 
  foreach (k= 1:ncol(mcomb),.combine='c') %do% { # for each combination of variables
    q2icc_nohh.matrix.b <- data[match(bootlist.n[[b]],data$datamemberid),]
    q2icc_nohh.matrix.b$newcluster <- names(bootlist[[b]])
    q2icc_nohh.matrix.b$newdataid <- names(bootlist.n[[b]])
    
    i<-mcomb[1,k] #ith risk factor
    j<-mcomb[2,k] #jth risk factor
    temp0 <- q2icc_nohh.matrix.b[ c(i,j,which(colnames(q2icc_nohh.matrix.b)=="datamemberid"),
                                    which(colnames(q2icc_nohh.matrix.b)=="newdataid")) ]    
    temp <- na.omit(temp0)
    fenzi <-0
    fenmu <-0
    for (p in 1:nrow(temp)) { # for each individual
      seedhh <- temp[p,]
      # subset to its household.....
      temp_hh <- temp[which(temp$newdataid==seedhh[1,4]),]
      # number of concordant pairs - self...only count TRUE
      if (nrow(temp_hh)!=0) {
        logiccount <- seedhh[1,1]==temp_hh[2]
        fenzi <- fenzi + length(logiccount[logiccount==TRUE])
        # number of pairs without any NA (both with numbers) - self...count TRUE and FALSE..exclude NA
        fenmu <- fenmu + nrow(temp_hh) 
      } else {
        fenzi <- fenzi+0
        fenmu <- fenmu+0
      }
    }
    fenzi/fenmu
  }
}

# cluster matrix for individual level exposures ####
clustermatrix_indexp <- function(boot.iter,
                                 Data=q2icc_nohh,
                                 Varlist=c("handfeed_never","prep2hrs_never",
                                              "q40_bc_never","q41_4_bc_never","q41_5_bc_never")){
  ##
  data <- Data[c(Varlist,"datamemberid","dataid","cluster") ]
  # create mcomb
  mcomb <- as.data.frame(combn(1:length(Varlist),2)) 
  ## 
  if (boot.iter!=0) {
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on matched-set level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist.hh <- lapply(bootlist.cluster,funcfunc)
    # bootstrap on household level ...names...new cluster
    funcfunc.ind  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(names(x)[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.hh,funcfunc.ind)
    # bootstrap on individual level .....names... new dataid
    funcfunc.ind.n  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(x[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist.n <- lapply(bootlist.hh,funcfunc.ind.n)  
  } else {
    bootlist.n <- data$datamemberid
    names(bootlist.n) <- substr(bootlist.n,1,5)
    bootlist.n <- list(bootlist.n)
    bootlist <- data$datamemberid
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
  }  
  # get co-occurrence of exposures
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  foreach (b = 1:boot.iter,.combine='rbind') %:%
  foreach (k = 1:ncol(mcomb),.combine='c') %do% { # for each combination of variables...
    q2icc_nohh.matrix.b <- data[match(bootlist.n[[b]], data$datamemberid),]
    
    q2icc_nohh.matrix.b$newcluster <- names(bootlist[[b]])
    q2icc_nohh.matrix.b$newdataid <- names(bootlist.n[[b]])
    
    i<-mcomb[1,k] #ith risk factor
    j<-mcomb[2,k] #jth risk factor
    temp0 <- q2icc_nohh.matrix.b[ c(i,j,which(colnames(q2icc_nohh.matrix.b)=="datamemberid"),
                                    which(colnames(q2icc_nohh.matrix.b)=="newcluster")) ]    
    temp <- na.omit(temp0)
    fenzi <-0
    fenmu <-0
    for (p in 1:nrow(temp)) { # for each individual
      seedhh <- temp[p,]
      # subset to its cluster
      temp_cluster <- temp[which(temp$newcluster==seedhh[1,4]),]
      # number of concordant pairs - self...only count TRUE
      if (nrow(temp_cluster)!=0) {
        logiccount <- seedhh[1,1]==temp_cluster[2]
        fenzi <- fenzi + length(logiccount[logiccount==TRUE])
        # number of pairs without any NA (both with numbers) - self...count TRUE and FALSE..exclude NA
        fenmu <- fenmu + nrow(temp_cluster) 
      } else {
        fenzi <- fenzi+0
        fenmu <- fenmu+0
      }
    }
    fenzi/fenmu
  }
}

# overall matrix for individual level exposures  ####
overallmatrix_indexp  <- function(boot.iter,
                                  Data=q2icc_nohh,
                                  Varlist=c("handfeed_never","prep2hrs_never",
                                               "q40_bc_never","q41_4_bc_never","q41_5_bc_never")){
  ##
  data <- Data[c(Varlist,"datamemberid","dataid","cluster") ]
  # create mcomb
  mcomb <- as.data.frame(combn(1:length(Varlist),2)) 
  ## 
  if (boot.iter!=0) {
    set.seed(3355)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on matched-set level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(8474)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist.hh <- lapply(bootlist.cluster,funcfunc)
    # bootstrap on household level ...names...new cluster
    funcfunc.ind  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(names(x)[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.hh,funcfunc.ind)  
    # bootstrap on individual level .....names... new dataid
    funcfunc.ind.n  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(2848)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(x[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist.n <- lapply(bootlist.hh,funcfunc.ind.n)  
  } else {
    bootlist.n <- list(data$datamemberid)
  }  
  # get co-occurrence of exposures
  boot.iter <- ifelse(boot.iter==0,1,boot.iter)
  foreach (b = 1:boot.iter,.combine='rbind') %:% 
  foreach (k = 1:ncol(mcomb),.combine='c') %do% {  # for each combination of variables
    q2icc_nohh.matrix.b <- data[match(bootlist.n[[b]], data$datamemberid),]
    
    i<-mcomb[1,k] #ith risk factor
    j<-mcomb[2,k] #jth risk factor
    #  subset a dataframe that include: rf1 rf2 datamemberid dataid
    temp0 <- q2icc_nohh.matrix.b[ c(i,j,which(colnames(q2icc_nohh.matrix.b)=="datamemberid")) ]    
    temp <- na.omit(temp0)
    fenzi <- 0 
    fenmu <- 0
    for (p in 1:nrow(temp)) { # for each hh
      # number of concordant pairs - self...only count TRUE
      if (nrow(temp)!=0) {
        seedhh <- temp[p,]
        logiccount <- seedhh[1,1]==temp[2]
        fenzi <- fenzi + length(logiccount[logiccount==TRUE])
        fenmu <- fenmu + nrow(temp) # number of pairs without any NA (both with numbers)
      } else {
        fenzi <- fenzi+0
        fenmu <- fenmu+0
      }
    }
    fenzi/fenmu # after all HHs' been looped through
  }
}