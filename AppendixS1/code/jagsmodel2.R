#MODEL 1:
#psi: landuse * shrub
#p: time + timeSq + date

model {
  
  # Priors
  alpha1 ~ dnorm(0, 0.01)
  alpha2 ~ dnorm(0, 0.01)
  alpha3 ~ dnorm(0, 0.01)
  alpha4 ~ dnorm(0, 0.01)
  
  beta ~ dnorm(0, 0.01)              ## intercept
  beta1 ~ dnorm(0, 0.01)             ## shrub cover
  beta.l[1] <- 0                     ## set first land-use type effect to zero
  beta.l[2] ~ dnorm(0, 0.01)         ## remaining 2 land uses
  beta.l[3] ~ dnorm(0, 0.01)
  
 
  #interaction effect landuse 
  beta.l.int[1] <- 0                ## set first land-use by shrub effect to zero
  for (i in 2:3){      
    beta.l.int[i] ~ dnorm(0, 0.01)    
  }
  
  #prior random effects
  sd.grid ~ dunif(0, 10)
  tau.grid <- pow(sd.grid, -2)
  
  for(r in 1:ngrid){
    grid.effect[r] ~ dnorm(0, tau.grid)
  }
  
  # Likelihood
  # Ecological model for the partially observed true state, psi
  for (i in 1:R) {
         
      #psi at transect level
      z[i] ~ dbern(psi[i])                
      logit(psi[i]) <- beta + beta.l[land[i]] + beta1*shrub[i] + 
        beta.l.int[land[i]]*shrub[i] + grid.effect[grid[i]] #effects parameterization
      
      # Observation model for p
      for (j in 1:J) {
        y[i,j] ~ dbern(mu.p[i,j])  # Detection-nondetection at i and j
        mu.p[i,j] <- z[i] * p[i,j]
        logit(p[i,j]) <- alpha1 + alpha2*date[i,j] + alpha3*time[i,j] + alpha4*time[i,j]
       
        #WAIC
        pd[i,j] <- dbin(y[i,j],z[i] * p[i,j], 1)
        log_pd[i,j] <- log(dbin(y[i,j],z[i] * p[i,j], 1))
        
      } #j      
    } #i


}