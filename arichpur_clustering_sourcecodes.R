# run these packages 
library(reshape)
library(foreach)
library(plyr)
library(ggplot2)
library(splines)

#### FUNCTION: concordance within matched-sets over spatial extent of matched-sets

### Varlist = exposure(s) of interest
### data = q1q4iccctrl for household level risk factors; = q2icc_nohh for individual level risk factors
### boot.iter = number of bootstrap iterations

# household level risk factors
within_hh_function_mltp <- function(Varlist,data=q1q4iccctrl,
                                    boot.iter=0){
  if (boot.iter !=0){
    # bootstrap on cluster level
    set.seed(1234)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(data$dataid[which(data$cluster==x[j])]))
        aa
      } 
    }
    set.seed(6)
    bootlist <- lapply(bootlist.cluster,funcfunc)
    b1<- rep(NA,boot.iter)
  } else {
    bootlist <- (data$dataid)
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
    b1 <-NA
  }
  
  ######linear coefficient
  for (i in 1:length(bootlist)) {
    dat <- data[match(bootlist[[i]],data$dataid),]
    dat$newcluster <- names(bootlist[[i]])
    #### calculate spatial extent of each cluster  = "medianclustersize"
    q1q4spatialextent <- data.frame() # empty dataframe 
    clusternames <- unique(dat$newcluster)  # all cluster IDs
    for(j in 1:length(clusternames)) {
      q1.onecluster <- dat[dat$newcluster==clusternames[j],]
      q1q4spatialextent[j,1]<-clusternames[j]  
      median.dist.onecluster <- median(dist(cbind(q1.onecluster$X,q1.onecluster$Y))) 
      q1q4spatialextent[j,2]<-median.dist.onecluster   
    }
    names(q1q4spatialextent) <- c("cluster","medianclustersize")
    ###
    for (Var in Varlist) {
      for (clus in unique(dat$newcluster)) {
        # In each matched-set, count household pairs having the same exposure  = "cncdc_count"
        junksubset <- dat[which(dat$newcluster==clus),c("dataid",Var)]
        # If a matched-set contains zero pair of hhs that both answered the question 
        if (  length(junksubset[,2][is.na(junksubset[,2])==FALSE]) <2) {
          cncdc_count = NA
        } else {
          # if at least one pair that both answered the question
          # discard the hh in the matched-set that answered NA for this variable
          junksubset2 <- junksubset[!is.na(junksubset[,2]),]
          junkcombn <- combn(junksubset2[,1],2)
          # and calculate the total pairs of households that have the same exposure
          cncdc_count = ifelse (0 %in% names(table(dist(junksubset2[,2]))),
                                table(dist(junksubset2[,2]))[1],0)
        }
        # get proportion
        q1q4spatialextent [which(q1q4spatialextent$cluster == clus),Var] <- cncdc_count/ncol(junkcombn)
      }
    }
    q1q4spatialextent <- cbind(q1q4spatialextent,
                               total=apply(as.matrix(q1q4spatialextent[,3:ncol(q1q4spatialextent)]),1,mean,na.rm=T))
    # bootstrap linear regression
    Var <- ifelse(length(Varlist)>1,"total",Varlist)
    b1[i] = lm(substitute(kk ~ medianclustersize, list(kk=as.name(Var))),
               data=q1q4spatialextent)$coeff[2]
  }
  quantile(b1,probs=c(0.025,0.975))*100
}
# individual level risk factors
within_ind_function_mltp <- function(Varlist,data=q2icc_nohh,
                                     boot.iter=0){
  if (boot.iter != 0){
    # bootstrap on cluster level
    set.seed(1234)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(6543)
        aa <- sample(unique(data$dataid[which(data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(data$dataid[which(data$cluster==x[j])])))
        aa
      } 
    }
    bootlist.hh <- lapply(bootlist.cluster,funcfunc)
    # bootstrap on individual level 
    funcfunc.ind  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(6543)
        aa <- sample(unique(data$datamemberid[which(data$dataid==x[j])]),replace=T)
        names(aa) <- rep(names(x)[j],
                         length(unique(data$datamemberid[which(data$dataid==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.hh,funcfunc.ind)
    b1<- rep(NA,boot.iter)
  } else {
    bootlist <- (data$datamemberid)
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
    b1 <-NA
  }
  
  ###### find linear coefficient
  for (i in 1:length(bootlist)) {
    dat <- data[match(bootlist[[i]],data$datamemberid),]
    dat$newcluster <- names(bootlist[[i]])
    #### calculate spatial extent of each cluster  = "medianclustersize"
    q1q4spatialextent <- data.frame() # empty dataframe 
    clusternames <- unique(dat$newcluster)  # all cluster IDs
    for(j in 1:length(clusternames)) {
      q1.onecluster <- dat[dat$newcluster==clusternames[j],]
      q1q4spatialextent[j,1]<-clusternames[j]  
      median.dist.onecluster <- median(dist(cbind(q1.onecluster$X,q1.onecluster$Y))) 
      q1q4spatialextent[j,2]<-median.dist.onecluster   
    }
    names(q1q4spatialextent)<- c("cluster","medianclustersize")
    ###
    for (Var in Varlist) {
      for (clus in unique(dat$newcluster)) {
        # In each matched-set, count household pairs having the same exposure  = "cncdc_count"
        junksubset <- dat[which(dat$newcluster==clus),c("dataid",Var)]
        # If a matched-set contains zero pair of hhs that both answered the question 
        if (  length(junksubset[,2][is.na(junksubset[,2])==FALSE]) <2) {
          cncdc_count = NA
        } else {
          # if at least one pair that both answered the question
          # discard the hh in the matched-set that answered NA for this variable
          junksubset2 <- junksubset[!is.na(junksubset[,2]),]
          junkcombn <- combn(junksubset2[,1],2)
          # and calculate the total pairs of households that have the same exposure
          cncdc_count = ifelse (0 %in% names(table(dist(junksubset2[,2]))),
                                table(dist(junksubset2[,2]))[1],0)
        }
        # get proportion
        q1q4spatialextent [which(q1q4spatialextent$cluster == clus),Var] <- cncdc_count/ncol(junkcombn)
      }
    }
    q1q4spatialextent <- cbind(q1q4spatialextent,
                               total=apply(as.matrix(q1q4spatialextent[,3:ncol(q1q4spatialextent)]),1,mean,na.rm=T))
    # bootstrap linear regression
    Var <- ifelse(length(Varlist)>1,"total",Varlist)
    # bootstrap linear regression
    b1[i] = lm(substitute(kk ~ medianclustersize, list(kk=as.name(Var))),
               data=q1q4spatialextent)$coeff[2]
  }
  quantile(b1,probs=c(0.025,0.975))*100
}

#### FUNCTION: visualize between matched-set concordance vs. distance between matched-set with LOESS

### Varlist = exposure of interest
### data = q1q4iccctrl for household level risk factors; = q2icc_nohh for individual level risk factors
### boot.iter = number of bootstrap iterations

between_function <- function(Var,Varlabel,Data,plot=F){
  # format data
  Data <- Data[,c(Var,"cluster","X","Y")]
  Data$cluster <- as.numeric(Data$cluster)
  #### distance between centroid in data.frame
  clustertrend_centroid_m <- ddply(Data,.(cluster),colwise(mean, na.rm=TRUE)  ) 
  ccombnt <- dist(clustertrend_centroid_m,diag=F)
  ccombnt.df <- melt(as.matrix(ccombnt), varnames = c("row", "col"))
  ccombnt.df <- ccombnt.df[ccombnt.df$row > ccombnt.df$col,]
  ccombnt.df$pdiff <- NA
  # find between matched-set concordance
  for (i in 1:nrow(ccombnt.df)) { 
    #subset matched-set 1
    #subset matched-set 2
    junkc = Data[which(Data$cluster==clustertrend_centroid_m[ccombnt.df[i,"row"],"cluster"]),Var]
    junkd = Data[which(Data$cluster==clustertrend_centroid_m[ccombnt.df[i,"col"],"cluster"]),Var]
    # delete NAs
    junkc<-junkc[!is.na(junkc)]
    junkd<-junkd[!is.na(junkd)]
    # calculate concordance
    countz=0
    for (j in 1:length(junkc)) {
      z<-junkd==junkc[j]
      countz<-length(z[z==TRUE])+countz
    }
    zz <- countz/ (length(junkc)*length(junkd))
    zz <- ifelse(zz==Inf, NA,zz)
    ccombnt.df$pdiff[i] <- zz
  }
  ccombnt.df <- rename(ccombnt.df, c(pdiff = Var))
  
  ######################
  # between matched-set concordance
  if (plot==T){
    ccombnt.df$species <- "between matched-set concordance"
    
    plot1 <- ggplot() +
      geom_point(data=ccombnt.df[which(ccombnt.df$value<780),], 
                 aes_string(x="value", y=Var) ,shape=19,size=1) +   
      geom_smooth(data=ccombnt.df[which(ccombnt.df$value<780),], 
                  aes_string(x="value", y=Var, color="species") , method="loess", size=1) +
      ylab(paste("Concordance of\n",Varlabel,sep="")) +
      xlab("Distance between matched-sets (m)") +
      scale_color_manual(values=c("within matched-set concordance"="blue", 
                                  "between matched-set concordance"="red")) +
      ylim(c(0,1.1))+
      scale_y_continuous(labels=c(0,0.2,0.4,0.6,0.8,1.0,1.2),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),
                         limits=c(0,1.12))+
      theme(legend.title=element_blank(),
            legend.position="bottom",#c(.65,-.92),
            legend.direction="vertical",
            legend.background = element_rect(fill=NA),
            legend.text = element_text(size = 20),
            axis.title.y = element_text(size = rel(1.8),vjust=0.8),
            axis.title.x = element_text(size = rel(1.8),vjust=-0.8),
            axis.text.y=element_text(size=18),
            axis.text.x=element_text(size=18))
    print(plot1)
  }
  ccombnt.df
}

####  FUNCTION: for exposures that showed clustering beyond the spatial extent of the 
####  matched-sets (concordance remained static beyond a certain distance threshold), 
####  find the AIC of the pattern characterzed by a linear spline with the knot at distance k (meter).

### k = threshold
### Data = output of the between_function

opfun <- function(k,Data){
  names(Data)[4] <- "Fakename"  # rename it
  spline.lm <- lm(Fakename ~ bs(value,degree=1, knots=c(k)),
                  data=Data)
  AIC(spline.lm)
}


####  obtain the linear association between concordance between matched-sets 
####  and distance between matched set up to a threshold

### threshold = optim
### Var =  exposures
### boot.iter = number of bootstrap iterations
### Data = q1q4iccctrl for household level risk factors; = q2icc_nohh for individual level risk factors

lmknot <- function(optim,Var,boot.iter,Data){
  #### new bootstrap
  if (boot.iter !=0){
    # bootstrap on cluster level
    set.seed(1234)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(Data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(6543)
        aa <- sample(unique(Data$dataid[which(Data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(Data$dataid[which(Data$cluster==x[j])]))
        aa
      } 
    }
    bootlist <- lapply(bootlist.cluster,funcfunc)
    b1<- rep(NA,boot.iter)
  } else {
    bootlist <- (Data$dataid)
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
    b1 <-NA
  }
  ###### find linear coefficient
  for (i in 1:length(bootlist)) {
    dat <- Data[match(bootlist[[i]],Data$dataid),]
    dat$cluster <- names(bootlist[[i]])
    ccombnt.bs <- dat[match(bootlist[[i]],dat$dataid),]
    tmp <- between_function(Var=Var,Data=ccombnt.bs,
                            plot=F)
    tmp2 <- tmp[which(tmp$value<optim),]
    b1[i] = lm(substitute(kk ~ value, list(kk=as.name(Var))),data=tmp2)$coeff[2]
  }
  aa<-quantile(b1,probs=c(0.025,0.975))*100
  return(aa)
}

lmknot.ind <- function(optim,Var,boot.iter,Data){
  #### new bootstrap
  if (boot.iter != 0){
    # bootstrap on cluster level
    set.seed(1234)
    bootlist.cluster <- lapply(1:boot.iter, function(i)  sample(unique(Data$cluster),replace=TRUE))
    # bootstrap on household level
    funcfunc  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(6543)
        aa <- sample(unique(Data$dataid[which(Data$cluster==x[j])]),replace=T)
        names(aa) <- rep(j,length(unique(Data$dataid[which(Data$cluster==x[j])])))
        aa
      } 
    }
    bootlist.hh <- lapply(bootlist.cluster,funcfunc)
    # bootstrap on individual level 
    funcfunc.ind  <- function(x){
      foreach(j=1:length(x), .combine='c') %do% {
        set.seed(6543)
        aa <- sample(unique(Data$datamemberid[which(Data$dataid==x[j])]),replace=T)
        names(aa) <- rep(names(x)[j],
                         length(unique(Data$datamemberid[which(Data$dataid==x[j])])))
        aa
      } 
    }
    bootlist <- lapply(bootlist.hh,funcfunc.ind)
    b1<- rep(NA,boot.iter)
  } else {
    bootlist <- (data$datamemberid)
    names(bootlist) <- substr(bootlist,2,5)
    bootlist <- list(bootlist)
    b1 <-NA
  }
  
  ###### find linear coefficient
  for (i in 1:length(bootlist)) {
    dat <- Data[match(bootlist[[i]],Data$datamemberid),]
    dat$cluster <- names(bootlist[[i]])
    ccombnt.bs <- dat[match(bootlist[[i]],dat$dataid),]
    tmp <- between_function(Var=Var,Data=ccombnt.bs,
                            plot=F)
    tmp2 <- tmp[which(tmp$value<optim),]
    b1[i] = lm(substitute(kk ~ value, list(kk=as.name(Var))),data=tmp2)$coeff[2]
  }
  aa<-quantile(b1,probs=c(0.025,0.975))*100
  return(aa)
}


