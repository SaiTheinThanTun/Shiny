## cut off value assessment tool for clearance half-life




nn <- 200 #default sample size
senmu <- 1.1 #log(3)
sensd <- 0.37 #log(1.45)

prop_resist <- .4
resmu <- 1.9 #log(6.5)
ressd <- .2 #log(1.22)
cutoff <- 5

#1. Choose a "sensitive" 1/2 life mean and SD 
#(Default: 1.1, 0.37 which corresponds to geometric mean of 3.0 hr)
hist(rlnorm(nn*(1-prop_resist),senmu,sensd))

sen_pop <- rlnorm(nn*(1-prop_resist),senmu,sensd) #sensitive population


#2. Choose a "resistant" 1/2 life mean and SD
#(Default: 1.9, 0.2 which corresponds to a geometric mean of 6.5 hr)
hist(rlnorm(nn*prop_resist,resmu,ressd)) #creating a histogram of random normal distribution with mean 1.9 and sd .2

res_pop <- rlnorm(nn*prop_resist,resmu,ressd) #resistant population


#3. Choose a proportion resistant (Default: 10%)

#4. Choose a sample size (Default: 200)

#5. Show a bar chart of a joint distribution
total_pop <- c(sen_pop,res_pop)
hist(total_pop, freq=FALSE,col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
lines(density(total_pop),lwd=5, col="red")
abline(v=cutoff, lwd=3, col="blue")


#6. Choose cut-off point (Default: 5 hr)

#7. List true/false sensitive/resistant, ROC curve, or some other way to visualize diagnostic
true_res <- plnorm(cutoff,resmu,ressd, lower.tail=FALSE) #% identified as resistant from truly resistant pop
fal_res <- plnorm(cutoff,senmu,sensd, lower.tail=FALSE) #% identified as resistant from sensitive pop

fal_sen <- plnorm(cutoff, resmu,ressd) #% wrongly identified as sensitive from truly resistant pop
true_sen <- plnorm(cutoff, senmu, sensd) #% identified as sensitive from truely sensitive pop

#test
mixdat <- as.data.frame(total_pop)  # import dataset
N<-1

# create output matrices
output.mu <- matrix(NA,nrow=M,ncol=N)
output.sigma <- matrix(NA,nrow=M,ncol=N)
output.lambda <- matrix(NA,nrow=M,ncol=N)
output.loglik <- matrix(NA,nrow=M,ncol=N)
output.mu.se <- matrix(NA,nrow=M,ncol=N)
output.sigma.se <- matrix(NA,nrow=M,ncol=N)
output.lambda.se <- matrix(NA,nrow=M,ncol=N)
AIC<-matrix(0,nrow=M,ncol=N)
AICdelta<-matrix(0,nrow=M,ncol=N)

#nb<-na.omit(mixdat[,N+1])

# fit single component model
for (i in 1:N){
  # 1 COMPONENT LOG NORMAL
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  xll<-fitdistr(lmixdat,"normal")
  output.loglik[1,i]<- xll$loglik
  output.mu[1,i]<-xll$estimate[1]
  output.lambda[1,i]<-1
  output.sigma[1,i]<-xll$estimate[2]
  output.mu.se[1,i]<-xll$sd[1]
  output.sigma.se[1,i]<-xll$sd[2]
  output.lambda.se[1,i]<-0
  AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
  AICdelta[1,i]<-0
}

# fit multiple component models sequntially
for (i in 1:N){
  nmixdat<-na.omit(mixdat[,i])
  lmixdat<- log(nmixdat)
  # >=2 COMPONENTS LOG NORMAL
  j<-1
  # stop if j-component model is more parsimonious than (j-1)-compnent model
  while((j<=M-1) && AICdelta[j,i]<=pval){
    j<-j+1
    res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
    resboot <- boot.se(res, B = nboot)
    resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
    output.loglik[j,i]<-res$loglik
    AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
    AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
    if(AICdelta[j,i]<=pval){
      output.mu[1:j,i]<-res$mu
      output.sigma[1:j,i]<-res$sigma
      output.lambda[1:j,i]<-res$lambda
      output.mu.se[1:j,i]<-resboot$mu.se
      output.sigma.se[1:j,i]<-resboot$sigma.se
      output.lambda.se[1:j,i]<-resboot$lambda.se			
    }
  }
}

##########################################################
# PLOTTING

# GRAPHS FOR SUPPORTING INFORMATION FILE 3
Sys.sleep(0.02)
for (ds in 1:N){
  Sys.sleep(0.02)
  nmixdat<-na.omit(mixdat[,ds])
  plam<-na.omit(output.lambda[,ds])
  pmu<-na.omit(output.mu[,ds])
  psig<-na.omit(output.sigma[,ds])
  hist(nmixdat,freq=FALSE,main = paste("Data from Dr. Kyaw, mean HL = ",round(exp(pmu[1]),2)),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
  x <- seq(0.1, max(nmixdat), length=1000)
  
  hx<-plam[1]*dlnorm(x,meanlog=(pmu[1]),sdlog=psig[1])
  if(length(plam)>1){
    for(k in 2:length(plam)){
      hx<-hx+plam[k]*dlnorm(x,meanlog=(pmu[k]),sdlog=psig[k])
    }
  }
  lines(x,hx,col="red", lwd=5)
  Sys.sleep(0.02)
}


#exp(pmu[1])
