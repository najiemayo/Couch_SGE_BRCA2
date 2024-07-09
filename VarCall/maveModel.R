model
{
    for (m in 1:M) {
        f.bv[m] ~ dt(mu[m],tau[m],tdf) ## model for replicate values
        tau[m] <- prec.reps*pow(gamma[batchM[m]],-2)
        mu[m] <- beta[batchM[m]] + eta[variantM[m]]*gamma[batchM[m]]
    }
   
    for (b in 1:B) {
        ## batch-specific random effects
        beta[b] ~ dnorm(mu.beta, prec.beta)
        gamma[b] ~ dlnorm(mu.gamma, prec.gamma)
    }

    for (v in 1:V) {
        ## variant-specific function effects:
        eta[v] ~ dnorm(alpha[del[v]], prec.eta[1])  ##LDA (see prec updates below)
    }

    for (v in 1:V){
        ##  Prior on Pathogenicity Status, Pr(D)
        pi.del[v] ~ dbeta(a.pp,b.pp)
        P[v,1] <- (1.0 - pi.del[v])
        P[v,2] <- pi.del[v]
        del[v] ~ dcat(P[v,1:2])
    }

    ## measurement variability of replicates
    prec.reps <- pow(sigma2.reps, -1.0)
    sigma2.reps <- pow(sigma.reps, 2.0)
    sigma.reps ~ dnorm(0.0, cauch.prec.reps)
    cauch.prec.reps ~ dgamma(0.5, 0.5)
    ## variability of batch coefficients in mean regression
    prec.beta <- pow(sigma2.beta, -1.0)
    sigma2.beta <- pow(sigma.beta, 2.0)
    sigma.beta ~ dnorm(0.0, cauch.prec.b)
    cauch.prec.b ~ dgamma(0.5, 0.5)
    ## variability of batch scaling in mean regression
    prec.gamma <- pow(sigma2.gamma, -1.0)
    sigma2.gamma <- pow(sigma.gamma, 2.0)
    sigma.gamma ~ dnorm(0.0, cauch.prec.g)
    cauch.prec.g ~ dgamma(0.5, 0.5)

    ## variability of variant coefficients in mean regression
    for (k in 1:2){
        prec.eta[k] <- pow(sigma2.eta[k], -1.0)
        sigma2.eta[k] <- pow(sigma.eta[k], 2.0)
    }
 
    ## function component variances
    sigma.eta[1] ~ dnorm(0.0, prec.se1)
    sigma.eta[2] <- sigma.eta[1] ##equal
    prec.se1 ~ dgamma(0.5, 0.5)
    prec.se2 ~ dgamma(0.5, 0.5)

    ## Prior on batch location and scale random effects means
    mu.beta ~ dunif(-10,10)
    mu.gamma ~ dnorm(0.0,1.0)
    
    prob <- mean(del[])-1
        
}

