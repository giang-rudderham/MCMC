# sampling window W = [0,1]x[0,1]
# n = 50
# Straussian pair potential function phi(s) = h if ||xi-xj|| < b (0 otherwise)
# code takes about 15-30 seconds to run on Matt's computer

set.seed(23)
n <- 50
x.curr <- runif(n,0,1)
y.curr <- runif(n,0,1)

b <- 0.1	# also try 0.05 
h <- 2		# also try 0 (binomial process), 0.5, 1, 5 (very strong repulsion b/w points)

x.prop <- x.curr
y.prop <- y.curr


# plot initial spatial point pattern
par(pty="s")
par(mfrow=c(1,2))
plot(x.curr, y.curr, main="Initial SPP",xlab='x',ylab='y')


# compute total potential energy (exponent in exponential function in likelihood)
# this double for loop is slow (we have 50 choose 2 terms)
tpe <- function(x,y,b,h) {
	tpe.tmp <- 0
	for(i in 1:(n-1)) {
		for(j in (i+1):n) {
			if(sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) < b) {
				tpe.tmp <- tpe.tmp + h
			}
		}
	}
	return(tpe.tmp)
}


# main mcmc loop
for(i in 1:2000) {

	# choose point to move
	i <- sample(1:n, 1)
	
	# propose to move i^th point uniformly in W 
	x.prop[i] <- runif(1,0,1)
	y.prop[i] <- runif(1,0,1)
	proposal.ratio <- 1
	
	# likelihood ratio
	tpe.prop <- tpe(x.prop,y.prop,b,h)
	tpe.curr <- tpe(x.curr,y.curr,b,h)
	likelihood.ratio <- exp(-tpe.prop + tpe.curr) 	# notice I'm not taking exp(...)/exp(...) 
													# this is more stable numerically
	# accept/reject
	acc.prob <- min(1, likelihood.ratio * proposal.ratio)
	if(runif(1,0,1) < acc.prob) { # accept
		x.curr <- x.prop
		y.curr <- y.prop
	}
	else { # reject (i.e. reset x.prop and y.prop)
		x.prop <- x.curr
		y.prop <- y.curr
	}
}


# plot resulting spatial point pattern
plot(x.curr, y.curr, main="Strauss SPP \n b = 0.1, h = 5",xlab='x',ylab='y')



