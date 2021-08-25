 model{
  for (i in 1:nIndiv)	 {			# i indexes individual animal
    for (j in 2:n[i])	{				# j indexes length measurement
      # observation - can have negative growth
      L[i, j] ~ dnorm(L_Exp[i, j], tau[i, j])	                       
      L_Exp[i, j] <-  Linf[i] *(1.0 - exp(-k[i]*(A[i]+t[i, j -1])))
      tau[i, j] <- 1 / (pow((L_Exp[i,j] * CV), 2))
      
      # posterior prediction
      #L.pred[i, j] ~ dnorm(L_Exp[i, j], tau)
      #p.value[i, j] <- step(L.pred[i, j] - L[i, j])
      
    }
    
    L[i, 1] ~ dnorm(L_Exp[i, 1], tau[i,1])
    L_Exp[i, 1] <-   Linf[i] *(1.0 - exp(-k[i]*A[i]))	
    tau[i,1] <- 1 / (pow((L_Exp[i,1] * CV), 2))
    #var1[i] <- 1/tau[i,1]
    
    # posterior prediction
    #L.pred[i, 1] ~ dnorm(L_Exp[i, 1], tau)
    #p.value[i, 1] <- step(L.pred[i, 1]- L[i, 1])
    
    Linf[i] ~ dnorm(LinfMu,  LinfTau)T(90,200)		
    k[i] ~ dbeta(kAlpha, kBeta)
    A[i] ~ dgamma(Shape, rate)
    
  }	
  # priors
  LinfSD <- sqrt(1/LinfTau)
  CV ~ dbeta(1, 1)
  
  # hyper-priors
  LinfMu ~ dnorm(100, 0.001)I(90,)  # indicator can be I(0,) but that's too broad?
  LinfTau ~ dgamma(0.001, 0.0001)
  Shape ~ dnorm(0, 0.001)I(0,)
  rate ~ dnorm(0, 0.001)I(0,)
  kAlpha ~ dnorm(0, 0.001)I(0,)
  kBeta ~ dnorm(0, 0.001)I(0,)
}
