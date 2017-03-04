# Y|theta ~ binomial(n,theta)
# theta ~ beta(alpha=1,beta=2/3)
# n=100, y=65

# actual posterior: 
# theta|y ~ beta(alpha+y,beta+n-y) = beta(66,35.667)
# E(theta|y) = 0.6492
# Var(theta|y) = 0.00222

set.seed(23)

n <- 100
y <- 65

theta <- c(0.5) # starting value

burnin <- 1000
iter <- 100000
accept <- 0

for(i in 1:(burnin+iter)) {
  
  # proposal ratio
  ################
  
  # beta proposal (good proposal; fast mixing)
  # 	thetaprime <- rbeta(1,3,2)
  # 	prop_ratio <- dbeta(theta[i],3,2) / dbeta(thetaprime,3,2)	
  
  # beta proposal (bad proposal; slow mixing)
  # 	thetaprime <- rbeta(1,1,5)
  # 	prop_ratio <- dbeta(theta[i],1,5) / dbeta(thetaprime,1,5)	
  
  # Random Walk (bad proposal; slow mixing)
  #  	thetaprime <- runif(1,theta[i]-0.001,theta[i]+0.001)
  # 	if(thetaprime>0 && thetaprime<1) prop_ratio <- 1 
  # 	else prop_ratio <- 0
  
  # Random Walk (good proposal: fast mixing)
  thetaprime <- runif(1,theta[i]-0.1,theta[i]+0.1)
  if(thetaprime>0 && thetaprime<1) prop_ratio <- 1 
  else prop_ratio <- 0
  
  
  # prior ratio
  #############
  prior_ratio <- dbeta(thetaprime,1,2/3) / dbeta(theta[i],1,2/3)
  
  
  # likelihood ratio
  ##################
  like_ratio <- dbinom(y,n,thetaprime) / dbinom(y,n,theta[i])
  
  
  # acceptance probability
  ########################
  acc_prob <- like_ratio * prior_ratio * prop_ratio
  
  
  # accept/reject
  ###############
  if(runif(1) < acc_prob) { # accept
    theta[i+1] <- thetaprime
    accept = accept + 1 
  }
  else { # reject
    theta[i+1] <- theta[i]
  }
}

#Report posterior mean and variance
cat("posterior mean:     ",mean(theta[-burnin]),"\n")
cat("posterior variance: ",var(theta[-burnin]),"\n")

par(mfrow=c(1,2))
#Plot MCMC Output
plot(theta,type='l',main='MCMC Output for Random Walk Proposal with Step 0.1',xlab='Iteration t', ylab='X')

#Plot Prior and Posterior Distributions
colors <- c('black','red','blue')
labels <- c('Estimated Posterior','Actual Posterior','Prior')
plot(density(theta[-burnin]),lwd=2,xlim=c(0,1), main='Prior and Posterior Distributions',xlab = 'x value',ylab='Density')
grd <- seq(0,1,0.001)
lines(grd,dbeta(grd,66,35.6667),col="red",lwd=2)
lines(grd,dbeta(grd,1,2/3),col="blue",lwd=2)
legend('topleft',lwd=2,labels,col=colors)

#acceptance rate
acc_rate <- accept/(burnin+iter)

#95% credible set
quantile(theta,0.025)
quantile(theta,0.975)

#Plot densities of Beta(3,2), Beta(1,5), Beta(66,35.67)
x<-seq(0,1,length=100)
colors <- c('red','blue','black')
labels <-c('Beta(3,2)','Beta(1,5)','Beta(66,35.67')
plot(x,dbeta(x,66,35.67),type='l',xlab='x value',ylab='Density',main='Comparison of Beta distributions')
lines(x,dbeta(x,3,2),col='red')
lines(x,dbeta(x,1,5),col='blue')
legend('topright', lwd=2,labels,col=colors)

#Geweke convergence diagnostic

library(boa)

theta_matrix <- as.matrix(theta)
rownames(theta_matrix) <- c(1:length(theta))
colnames(theta_matrix) <- c("theta")
dimnames(theta_matrix)
boa.geweke(theta_matrix, 0.1, 0.5)
