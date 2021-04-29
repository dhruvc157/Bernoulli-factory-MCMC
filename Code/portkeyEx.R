#########################
## Bernoulli factory MCMC
#########################

## This code illustrates working of Portkey two-coin algo and
## compares it with the usual two-coin to get MCMC samples 
## from target \pi

set.seed(713)
#  Opening the Gringotts vault to the hidden Sorcerer's Stone

mcmc_sampler = function(N = 1e3, beta = 0.95, K, lam_mu, lam_var)
{
  "
  The function will generate MCMC samples from a accept-reject MCMC algo
  based on portkey-barker's acc.prob \alpha_{(\beta)}
  if \beta = 1, then it's Barker's acc. prob., and MCMC samples are based on 2 coin algo.
  
  inputs: N, beta, K (shape param of Weibull), lam_mu, lam_var (mean and variance for normal dist)
  output: list ( MCMC samples, 
                vec. of #loops to decide on the proposed value, 
                acc. prob.)
  
  "
  chain = numeric(N) # MCMC samples to be obtained 
  loops = numeric(N) # no. of loops to simulate 1 event of probability 
                     # \alpha_{(\beta)} default initialization to 0
  acc = 0 # No. of accepted proposals
  
  chain[1] = 0.75
  for ( i in 2:N)
  {
    curr = chain[i-1]
    prop = rnorm(1, mean = curr, sd = sqrt(0.325))
    if (prop < 0)   
    {
      # reject proposal and hop
      # loops[i] = 0 (by default)
      chain[i] = curr
      next
    }  
    out = portkey2coin(beta = beta, curr = curr, prop = prop, K, lam_mu, lam_var)    
    if (out[[1]] == 1) # accept proposal
    {
      chain[i] = prop
      acc = acc + 1
    }  
    else  chain[i] = curr # reject proposal and hop
    
    loops[i] = out[[2]]
    
  }
  return(list("samples" = chain, "loops" = loops, "acceptance prob" = acc/N))
}

c_val = function(val, K)
{
  return (K/ (exp(1)*val))
}
p_cond = function(val, lam, K)
{
  #lam = abs(rnorm(1, mean = lam_mu, sd = sqrt(lam_var)))
  return (dweibull(val, shape = K, scale = lam))
}

portkey2coin = function(beta, curr, prop, K, lam_mu, lam_var)
{
  "
  The function will simulate an event with probabilty \alpha_{(\beta)}
  inputs: 
  output: list ( X ~ Bern (\alpha_{(\beta)}), #loops)
  
  "
  
  looper = 0
  flag = 0
  while (TRUE)
  {
    looper = looper + 1
    S = 1
    
    if(beta != 1)
      S = rbinom(1, size = 1, prob = beta)
    
    if (S == 0)  
      return(list("X" = 0, "loops" = looper))
    
    C1 = rbinom(1, size = 1, 
                prob = c_val(prop, K) / (c_val(curr, K) + c_val(prop, K)))
    
    lam = abs( rnorm(1, mean = lam_mu, sd = sqrt(lam_var)) )
    U = runif(1)
    if (C1 == 1)
    {
      # C2 is a simulated event of probability p_prop
      if (U < ( p_cond(prop, lam, K)/ c_val(prop, K)) )
        return(list("X" = 1, "loops" = looper))
      
      else 
        next
        
      
    } else {
      
      # C2 is a simulated event of probability p_curr
      if (U < ( p_cond(curr, lam, K)/ c_val(curr, K)) )  
        return(list ("X" = 0, "loops" = looper))
      
      else 
        next
        
    }
    
  }
  
}

lam_mu = 0
lam_var = 1
K = 10

chains = list()
times = list()
betas = c(1, 0.99, 0.90, 0.75)
sd.store = matrix(nrow = 1e3, ncol = 4)
ess.store = matrix(nrow = 1e3, ncol = 4)
ess.psec.store = matrix(nrow = 1e3, ncol = 4)
mean.loops = matrix(nrow = 1e3, ncol = 4)
max.loops = matrix(nrow = 1e3, ncol = 4)

#install.packages('mcmcse')
library(mcmcse)
#install.packages("doParallel")
#library(doParallel)
#detectCores()
#registerDoParallel(detectCores()-1)

for(j in 1:1e3)
{
  for(i in 1:4) 
  {
    times[[i]] = system.time({chains[[i]] = mcmc_sampler(N = 1e5, beta = betas[i], K, lam_mu, lam_var)})
    sd.store[j,i] = sd(chains[[i]]$samples)  
    ess.store[j,i] = ess(chains[[i]]$samples)
    ess.psec.store[j,i] = as.numeric(ess.store[j,i]/ times[[i]][3])
    mean.loops[j,i] = mean(chains[[i]]$loops)
    max.loops[j,i] = max(chains[[i]]$loops)
  }
  print(paste("Executed iteration: ", j, "     ", Sys.time()))
}


par(mfrow = c(1,2))
xax = seq(to = 1e5, length.out = 1e3) 
plot(xax, tail(chains[[1]]$loops, n = 1e3), 
     pch = 8, col = "red", xlab = "Chain index", ylab = "#loops to converge")
points(xax, tail(chains[[2]]$loops, n = 1e3), 
       pch = 20, col = "royalblue")
legend( x="topright", pch=c(8,20), cex = 0.65,
        legend=c("2 coin algorithm","Portkey Bernoulli factory"), 
        col=c("red","royalblue"))

plot(xax, tail(chains[[1]]$samples, n = 1e3), 
        col = "red", xlab = "Chain index", type = "l",
        ylab = expression(paste(X, "~" ,pi, "(", theta,") (approx)")))
lines(xax, tail(chains[[2]]$samples, n= 1e3), col = "royalblue")
legend( x="topright", lty = c(1,1), cex = 0.65,
        legend=c("2 coin algorithm","Portkey Bernoulli factory"), 
        col=c("red","royalblue") )

par(mfrow = c(2,2))
for ( i in 1:4)
{
  plot(xax, tail(chains[[i]]$samples, n = 1e3), col = "royalblue",
       xlab = paste("Acceptance prob: ", chains[[i]]$`acceptance prob`), type = "l", main = bquote(beta == .(paste(betas[i]))),
       ylab = "MCMC Samples")
}

par(mfrow = c(1,2))

cols = c("red", "royalblue", "darkgreen", "orange")
plot(density(chains[[1]]$samples), col = cols[1], lty = 1, lwd = 2,
     ylim = c(0, 1) , xlab = "x", main = "Density plot")
for ( i in 2:4)
{
  lines(density(chains[[i]]$samples), col = cols[i], lty = i, lwd = 2)
}
legend( x="topright", lty = c(1,2,3,4), cex = 0.75,
        legend=c("1","0.99", "0.90", "0.75"), title = expression(beta),
        col=cols )


plot(acf(chains[[1]]$samples, plot = FALSE)$lag,
     acf(chains[[1]]$samples, plot = FALSE)$acf, 
     type = "l", lty = 1, lwd = 2, col = cols[1],
     xlab = "Lag", ylab = "ACF", main = "Auto-correlation plot")

for (i in 2:4)
{
  lines(acf(chains[[i]]$samples, plot = FALSE)$lag,
        acf(chains[[i]]$samples, plot = FALSE)$acf, 
        type = "l", lwd = 2, lty = i, col = cols[i])
  
}

legend( x="topright", lty = c(1,2,3,4), cex = 0.75,
        legend=c("1","0.99", "0.90", "0.75"), title = expression(beta),
        col=cols )


apply(ess.store, 2, mean)
apply(ess.psec.store, 2, mean)
apply(mean.loops, 2, mean)
apply(max.loops, 2, mean)

apply(ess.store, 2, sd)
apply(ess.psec.store, 2, sd)
apply(mean.loops, 2, sd)
apply(max.loops, 2, sd)

