
> cat("model
+ {
+ 	# Parameters
+ 	# phi: survival
+ 	# psiIO: probability to emigrate
+ 	# psiOI: probability to immigrate
+ 	# p: recapture probabi ..." ... [TRUNCATED] 
model
{
	# Parameters
	# phi: survival
	# psiIO: probability to emigrate
	# psiOI: probability to immigrate
	# p: recapture probability
	
	# states (S)
	# 1 alive and present
	# 2 alive and absent
	# 3 dead
	#
	# observations (O)
	# 1 seen
	# 2 not seen
	
	# priors and constraints
	
	for (t in 1:(n.occasions-1)){
		phi[t] ~ dunif(0.5,1.0) #<- mean.phi
		psiIO[t] <- mean.psiIO
		psiOI[t] <- mean.psiOI
		p[t] <- mean.p
	}
	
	#mean.phi ~ dunif(0, 1)
	mean.psiIO ~ dunif(0, 1)
	mean.psiOI ~ dunif(0, 1)
	mean.p ~ dunif(0, 1)
	
	# define state-transition and observation matrices
	# define probabilities of state S(t+1) given S(t)
	for (t in 1 : (n.occasions-1)){
		ps[1, t, 1] <- phi[t] * (1 - psiIO[t])
		ps[1, t, 2] <- phi[t] * psiIO[t]
		ps[1, t, 3] <- 1 - phi[t]
		ps[2, t, 1] <- phi[t] * psiOI[t]
		ps[2, t, 2] <- phi[t] * (1 - psiOI[t])
		ps[2, t, 3] <- 1 - phi[t]
		ps[3, t, 1] <- 0
		ps[3, t, 2] <- 0
		ps[3, t, 3] <- 1
			
		# define probabilities of O(t) given S(t)
		po[1, t, 1] <- p[t]
		po[1, t, 2] <- 1 - p[t]
		po[2, t, 1] <- 0
		po[2, t, 2] <- 1
		po[3, t, 1] <- 0
		po[3, t, 2] <- 1

	}
	
	# likelihood
	for (i in 1:nind){
		# define latent state at first capture
		z[i, f[i]] <- y[i, f[i]]
		for (t in (f[i]+1):n.occasions){
			# state process: draw S(t) given S(t-1)
			z[i,t] ~ dcat(ps[z[i, t-1], t-1,])
			# observation process: draw O(t) given S(t)
			y[i,t] ~ dcat(po[z[i,t], t-1,])
		}
	}
}

> sink()
