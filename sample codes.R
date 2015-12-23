# linear association between concordance within matched-sets and 
# spatial extent of matched-sets (results in Fig 2C)
within_hh_function_mltp(Varlist="q11_supply",data=q1q4iccctrl,boot.iter=0)
# bootstrapped ci
within_hh_function_mltp(Varlist="q11_supply",data=q1q4iccctrl,boot.iter=1000)


# visualize the association between concordance between matched-set and ..
# distance between matched-sets (Fig 3)
between.linear.boot <- between_function(Varlabel="use of supplied water",
                                      Var="q11_supply",
                                      Data=q1q4iccctrl,plot=T)
# where does the association level off?
optim <- optimize(opfun,Data=between.linear.boot,lower=min(between.linear.boot$value),
                  upper=max(between.linear.boot$value),
                  maximum=F)$minimum

# find the linear association between concordance and distance 
# before concordance levels off...(Results in Figure 2E)
lmknot(optim=optim,Var="q11_supply",boot.iter=0,Data=q1q4iccctrl)
# bootstrapped ci
lmknot(optim=optim,Var="q11_supply",boot.iter=1000,Data=q1q4iccctrl)


# co-occurrence between household level exposures within the same households  (Fig 4)
probratio.h.c <- hhmatrix_hhexp(boot.iter=0)/clustermatrix_hhexp(boot.iter=0)
# bootstrapped estimates
probratio.h.c.boot <- hhmatrix_hhexp(boot.iter=10)/clustermatrix_hhexp(boot.iter=10)
probratio.h.c.boot.ci<- apply(probratio.h.c.boot,2,quantile,probs=c(0.025,0.975),na.rm=TRUE) 
# co-occurrence between household level exposures within the same matched-sets
probratio.c.o <- clustermatrix_hhexp(boot.iter=0)/overallmatrix_hhexp(boot.iter=0)
# bootstrapped estimates
probratio.c.o.boot <- clustermatrix_hhexp(boot.iter=10)/overallmatrix_hhexp(boot.iter=10)
probratio.c.o.boot.ci<- apply(probratio.c.o.boot,2,quantile,probs=c(0.025,0.975),na.rm=TRUE) 

Varlist = c("pplrm_di","defecate_pitl","defsharing_only",
            "dist_more10", "q4water2","q11_tube",
            "q11_supply","boilnewst","soap")
mcomb <- as.data.frame(combn(1:length(Varlist),2)) 
mcomb2 <- cbind(Var1=Varlist[as.numeric(mcomb[1,])],   # exposure 1
                Var2=Varlist[as.numeric(mcomb[2,])],   # exposure 2
                hc=probratio.h.c,                      # co-occurrence within households
                hc_95ci_low=probratio.h.c.boot.ci[1,], # and its lower bound of 95%ci
                hc_95ci_high=probratio.h.c.boot.ci[2,],# and its higher bound of 95%ci
                co=probratio.c.o,                      # co-occurrence within matched-sets
                co_95ci_low=probratio.c.o.boot.ci[1,], # and its lower bound of 95%ci
                co_95ci_high=probratio.c.o.boot.ci[2,] # and its higher bound of 95%ci
                )

