# Likelihood: Y_i|mu,sigma ~ N(mu,sigma) i=1,...,n (indep)
#    likelihood will be product of n normal dists
# Priors: mu ~ N(0,100), sigma ~ chisq(2) (indep)
# Actual posterior is unknown (we are using non-conjugate priors)

# data
######
n = 10
y <- c(3.659798, 6.203945, 4.200632, 6.837608, 1.433226, 
	6.449043, 3.906239, 1.211519, 2.543834, 2.785943)

# starting values
mu <- c(0)
sigma <- c(1)
sd <- 1

burnin <- 2000
iter <- 10000
move_mu <- 0 
move_sigma <- 0 
accept_mu <- 0 
accept_sigma <- 0 

set.seed(23)

for(i in 1:(burnin + iter)) {

	# proposal ratio
	################

	# choose move type at random
	u <- runif(1)
	if (u <= 0.5) {
	  movetype = "mu"
	  move_mu = move_mu + 1 
	}
	else {
	  movetype = "sigma"
	  move_sigma = move_sigma + 1 
	}

	# RW proposal (1-at-a-time)	
	  # RW proposal for mu
	if (movetype == "mu") {
		muprime <- runif(1, mu[i] - 5, mu[i] + 5)
		sigmaprime <- sigma[i]
		prop_ratio <- 1
	}	
	  #Normal RW proposal for mu
#   if (movetype == 'mu'){
#     muprime <- rnorm(1, mu[i], sd)
#     sigmaprime <- sigma[i]
#     prop_ratio <- dnorm(mu[i], muprime, sd) / dnorm(muprime, mu[i], sd)
#   }
	 	
	if(movetype == "sigma") {
		muprime <- mu[i]
		sigmaprime <- runif(1, sigma[i] - 1, sigma[i] + 1)
		if (sigmaprime > 0) prop_ratio <- 1
		else prop_ratio <- 0
	}
	
	
	# prior ratio
	#############
	if (sigmaprime > 0) {
		prior_ratio <- (dnorm(muprime, 0, 100) * dchisq(sigmaprime, 2)) /
					   (dnorm(mu[i], 0, 100) * dchisq(sigma[i], 2))
	}
  else prior_ratio <- 0
  
	# likelihood ratio
	##################
	if (sigmaprime > 0) {
		like_ratio <- 1
		for (j in 1:n) {
			like_ratio <- like_ratio * dnorm(y[j], muprime, sigmaprime) / 
			  dnorm(y[j], mu[i], sigma[i])
		}
	}
  else like_ratio <- 0 
		
	# acceptance probability
	########################
	acc_prob <- like_ratio * prior_ratio * prop_ratio

	
	# accept/reject
	###############
	if (runif(1) < acc_prob) { # accept
		mu[i+1] <- muprime
		sigma[i+1] <- sigmaprime
		if (movetype == 'mu') accept_mu = accept_mu + 1
		if (movetype == 'sigma') accept_sigma = accept_sigma + 1
	}
	else { # reject
		mu[i+1] <- mu[i]
		sigma[i+1] <- sigma[i]
	}
}

# Report estimated posterior means
cat("posterior mean of mu:    ", mean(mu[-burnin]), "\n")
cat("posterior mean of sigma: ", mean(sigma[-burnin]), "\n")

# Plot MCMC output and estimated posterior distribution
par(mfrow = c(2, 2))
grd <- seq(0, 15, 0.001)
plot(mu, type = 'l', main = 'MCMC Output for Mu \n RW +/-5 Proposal', 
     xlab = 'Iteration t', ylab = 'X')
plot(density(mu[-burnin]), lwd=2, xlab = 'mu value', ylab ='Density', 
     main = 'Estimated Posterior Distribution for Mu \n RW +/-5 Proposal')
grd <- seq(0, 5, 0.001)
plot(sigma, type='l', xlab = 'Iteration t', ylab = 'X', 
     main = 'MCMC Output for Sigma \n RW +/-1 Proposal')
plot(density(sigma[-burnin]), lwd=2, xlab = 'sigma value', ylab = 'Density', 
     main = 'Estimated Posterior Distribution for Sigma \n RW +/-1 Proposal')

#95% credible set
quantile(mu, 0.025)
quantile(mu, 0.975)
quantile(sigma, 0.025)
quantile(sigma, 0.975)

#acceptance rate
accept_mu / move_mu
accept_sigma / move_sigma

#Geweke's convergence diagnostic
theta_matrix <- cbind(mu, sigma)
rownames(theta_matrix) <- c(1:length(mu))
boa.geweke(theta_matrix, 0.1, 0.5)
