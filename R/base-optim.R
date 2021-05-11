################################################################################
################################################################################
### R script w/ optimization routines, objective functions, and derivatives
################################################################################
################################################################################

################################################################################
### Optimizer for glmNB and glmZINB
################################################################################

optimNB <- function(par,y,x,offset,weights,control,dist){
    model <- NULL
    if (dist == 'nb') {
        tryCatch({
            lowerLimit <- c(rep(-Inf,ncol(x)),1/control$maxDisp)
            assign('model',{
                stats::optim(par = invertDispersion(par, model = dist),
                             fn = glmNB, gr = derivNB, method = 'L-BFGS-B',
                             lower = lowerLimit, y = y, x = x, offset = offset,
                             weights = weights)   
            })},
            error = function(e) {
                assign('model',{
                    list('par' = invertDispersion(par, model = dist),
                         'convergence' = 99)
                },inherits = TRUE)
            })
        model$par <- invertDispersion(model$par, model = dist)
    } else{
        tryCatch({
            lowerLimit <- c(rep(-Inf,ncol(x)),1/control$maxDisp,
                            rep(-Inf,ncol(x)))
            assign('model',{
                stats::optim(par = invertDispersion(par, model = dist),
                             fn = glmZINB,gr = derivZINB,method = 'L-BFGS-B',
                             lower = lowerLimit,y = y,x = x,offset = offset,
                             weights = weights,minZero = control$minZero)
            })},
            error = function(e) {
                assign('model',{
                    list('par' = invertDispersion(par, model = dist),
                         'convergence' = 99)   
                },inherits = TRUE)
            })
        if (ncol(x) == 1) {
            model$par <- c(model$par[c(3,1)],1/model$par[2])
        } else{
            model$par <- c(model$par[c(4,5)],model$par[c(1,2)],1/model$par[3])
        }
    }
    return(model)
}

################################################################################
### Optimizer for glmMixNB and glmMixZINB
################################################################################

optimDifferential <- function(par,rcWeights,control,dist){
    model = NULL
    if (dist == 'nb') {
        tryCatch({
            lowerLimit <- rep(log(control[['minZero']]),4)
            upperLimit <- c(Inf,Inf,rep(log(control[['maxDisp']]),2))
            assign('model',stats::optim(par = par[c(1,3,2,4)],
                                        fn = glmMixNB,
                                        gr = derivMixNB,
                                        method = 'L-BFGS-B',
                                        lower = lowerLimit,
                                        upper = upperLimit,
                                        dt = rcWeights,
                                        minZero = control[['minZero']],
                                        maxid = length(unique(rcWeights$id))))
        },error = function(e){
            assign('model',list('par' = par, 'convergence' = 99),
                   inherits = TRUE)
        })
        model$par <- model$par[c(1,3,2,4)]
    } else{
        newpar <- c(par[c(1, 2, 3)],
                    par[2] + par[4],
                    par[3] + par[5])
        tryCatch({
            lowerLimit <- rep(log(control[['minZero']]),5)
            upperLimit <- c(Inf,rep(c(Inf,log(control[['maxDisp']])),2))
            assign('model',stats::optim(par = newpar,
                                        fn = glmMixZINB,
                                        gr = derivMixZINB,
                                        method = 'L-BFGS-B',
                                        lower = lowerLimit,
                                        upper = upperLimit,
                                        dt = rcWeights,
                                        minZero = control[['minZero']],
                                        maxid = length(unique(rcWeights$id))))
        },error = function(e){
            assign('model',list('par' = par, 'convergence' = 99),
                   inherits = TRUE)
        })
        model$par <- c(model$par[c(1, 2, 3)],
                       model$par[4] - model$par[2],
                       model$par[5] - model$par[3])
    }
    return(model)
}

################################################################################
### Log-likelihood function of a GLM NB distribution
################################################################################

glmNB = function(par,y,x,offset,weights){
    mu.vec <- exp(x %*% par[seq_len(ncol(x))] + offset)
    l <- stats::dnbinom(y,mu=mu.vec,size=1/par[length(par)],log=TRUE)
    l[is.infinite(l)] <- log(1e-300)
    return(-sum(weights*l))
}

################################################################################
### Derivative a GLM NB distribution and its derivative
################################################################################

derivNB = function(par,y,x,offset,weights){
    mu.vec <- exp(x %*% par[seq_len(ncol(x))] + offset)
    dldbeta <- colSums(as.numeric(weights)*
                           as.numeric((y-mu.vec)/(1+par[length(par)]*mu.vec))*x)
    dldphi <- sum(as.numeric(weights) *
                      (log(1+par[length(par)]*mu.vec) + 
                           par[length(par)]*(y-mu.vec)/(1+par[length(par)]*mu.vec) - 
                           digamma(y+1/par[length(par)]) + 
                           digamma(1/par[length(par)]))/(par[length(par)]^2))
    return(-c(dldbeta,dldphi))
}

################################################################################
### Log-likelihood function of a GLM ZINB distribution
################################################################################

glmZINB = function(par,y,x,offset,weights,minZero){
    mu.vec <- exp(x%*%par[seq_len(ncol(x))]+offset)
    zeroinfl.vec <- 1/(1+exp(-(x%*%par[(ncol(x)+2):(2*ncol(x)+1)]+offset)))
    idx.Y0 <- (y==0)
    d0 <- log(stats::dnbinom(0,mu=mu.vec,size=1/par[(ncol(x)+1)],
                             log=FALSE)+minZero)
    d1 <- log(stats::dnbinom(y,mu=mu.vec,size=1/par[(ncol(x)+1)],
                             log=FALSE)+minZero)
    l0 <- log(zeroinfl.vec + exp(log(1-zeroinfl.vec) + d0))
    l1 <- log(1-zeroinfl.vec) + d1
    return(-sum(weights*(idx.Y0*l0+(1-idx.Y0)*l1)))
}

################################################################################
### Derivative of a GLM ZINB distribution
################################################################################

derivZINB = function(par,y,x,offset,weights,minZero){
    mu.vec <- exp(x%*%par[seq_len(ncol(x))]+offset)
    zeroinfl.vec <- 1/(1+exp(-(x%*%par[(ncol(x)+2):(2*ncol(x)+1)]+offset)))
    idx.Y0 <- (y==0)
    l0 <- log(stats::dnbinom(0,mu=mu.vec,size=1/par[(ncol(x)+1)],
                             log=FALSE)+minZero)
    l1 <- log(stats::dnbinom(y,mu=mu.vec,size=1/par[(ncol(x)+1)],
                             log=FALSE)+minZero)
    P0 <- zeroinfl.vec+exp(log(1-zeroinfl.vec)+l0)
    P1 <- exp(log(1-zeroinfl.vec)+l1)
    phi <- par[(ncol(x)+1)]
    aux <- (1+phi*mu.vec)
    dldb1 <- as.numeric((idx.Y0)*(weights/P0)*(1-zeroinfl.vec)*
                            (-mu.vec*((1+phi*mu.vec)^(-(1/phi+1)))))*x
    dldb2 <- as.numeric((1-idx.Y0)*(weights/P1)*(1-zeroinfl.vec)*
                            exp(l1)*(y-mu.vec)/(1+phi*mu.vec))*x
    dldphi1 <- (idx.Y0)*(weights/P0)*(1-zeroinfl.vec)*(aux^(-1/phi))*
        (log(aux)/((phi)^2)-mu.vec/(phi*aux))
    dldphi2 <- (1-idx.Y0)*(weights/P1)*(1-zeroinfl.vec)*exp(l1)*
        (log(aux)/(phi^2)-mu.vec*(y+1/phi)/aux-digamma(y+1/phi)/(phi^2)+
             digamma(1/phi)/(phi^2)+y/phi)
    dldlambda1 <- as.numeric((idx.Y0)*(weights/P0)*(1-exp(l0))*
                                 zeroinfl.vec*(1-zeroinfl.vec))*x
    dldlambda2 <- as.numeric((1-idx.Y0)*(weights/P1)*(-exp(l1))*
                                 zeroinfl.vec*(1-zeroinfl.vec))*x
    return(-c(colSums(dldb1+dldb2),
              sum(dldphi1+dldphi2),
              colSums(dldlambda1+dldlambda2)))
}

################################################################################
### Log-likelihood function of 3-state HMM with mixture of NBs
################################################################################

glmMixNB = function(par,dt,maxid,minZero){
    Intercept = weights = ChIP = Dsg.Mix = offsets = id = NULL
    llBackground <- sum(dt[(id == 1),][,{
        -weights*log(stats::dnbinom(x = ChIP,
                                    mu = exp(par[1]*Intercept + offsets),
                                    size = exp(par[3]*Intercept),
                                    log = FALSE) + minZero)
    }])
    llDifferential <- do.call(sum,lapply(2:(maxid - 1),FUN = function(x){
        sum(dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,{
            -weights*log(stats::dnbinom(x = ChIP,
                                        mu = exp(par[1]*Intercept + par[2]*Dsg.Mix + offsets),
                                        size = exp(par[3]*Intercept + par[4]*Dsg.Mix),
                                        log = FALSE) + minZero)
        }])
    }))
    llEnrichment <- sum(dt[(id == maxid),][,{
        -weights*log(stats::dnbinom(x = ChIP,
                                    mu = exp((par[1]+par[2])*Intercept+offsets),
                                    size = exp((par[3] + par[4])*Intercept),
                                    log = FALSE) + minZero)
    }])
    return(llBackground + llDifferential + llEnrichment)
}

################################################################################
### Derivative of 3-state HMM with mixture of NBs
################################################################################

derivMixNB = function(par,dt,maxid,minZero){
    Intercept = offsets = weights = ChIP = mu = phi = Dsg.Mix = id = NULL
    derivBackground <- lapply(1,FUN = function(x){
        subdt <- dt[(id == x),][,c('mu','phi') := {
            list(exp(par[1] * Intercept + offsets),
                 exp(par[3] * Intercept))
        }]
        colSums(cbind(subdt[,{
            -weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept,0)
        }],
        subdt[,{
            -weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept,0)
        }]))
    })
    derivDifferential <- lapply(2:(maxid - 1),FUN = function(x){
        subdt <- dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,c('mu','phi') := {
            list(exp(par[1] * Intercept + par[2] * Dsg.Mix + offsets),
                 exp(par[3] * Intercept + par[4] * Dsg.Mix))
        }]
        colSums(cbind(subdt[,{
            -weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept,Dsg.Mix)
        }],
        subdt[,{
            -weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept,Dsg.Mix)
        }]))
    })
    derivEnrichment <- lapply(maxid,FUN = function(x){
        subdt <- dt[(id == x),][,c('mu','phi') := {
            list(exp((par[1] + par[2]) * Intercept + offsets),
                 exp((par[3] + par[4]) * Intercept))
        }]
        colSums(cbind(subdt[,{
            -weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept,Intercept)
        }],
        subdt[,{
            -weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept,Intercept)
        }]))
    })
    return(colSums(do.call(rbind,c(derivBackground,derivDifferential,derivEnrichment))))
}

################################################################################
### Log-likelihood function of 3-state HMM with mixture of ZINB and NB
################################################################################

glmMixZINB = function(par,dt,maxid,minZero){
    Dsg.Mix = weights = ChIP = Intercept = offsets = id = NULL
    llBackground <- sum(dt[(id == 1),][,{
        -weights*dzinb(x = ChIP,
                       mu = exp(par[2]*Intercept + offsets),
                       size = exp(par[3]*Intercept),
                       zip = exp(par[1])/(1 + exp(par[1])),
                       minZero = minZero)
    }])
    llDifferential <- do.call(sum,lapply(2:(maxid - 1),FUN = function(x){
        sum(dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,{
            -weights*((Dsg.Mix == 1)*log(stats::dnbinom(x = ChIP,
                                             mu = exp(par[4]*Intercept+offsets),
                                             size = exp(par[5]*Intercept),
                                             log = FALSE) + minZero) + 
                          (Dsg.Mix == 0)*dzinb(x = ChIP,
                                             mu = exp(par[2]*Intercept+offsets),
                                             size = exp(par[3]*Intercept),
                                             zip = exp(par[1])/(1+exp(par[1])),
                                             minZero = minZero))
        }])
    }))
    llEnrichment <- sum(dt[(id == maxid),][,{
        -weights*log(stats::dnbinom(x = ChIP,mu = exp(par[4]*Intercept+offsets),
                                    size = exp(par[5]*Intercept),
                                    log = FALSE) + minZero)   
    }])
    return(llBackground + llDifferential + llEnrichment)
}

################################################################################
### Derivative of 3-state HMM with mixture of ZINB and NB
################################################################################

derivMixZINB = function(par,dt,maxid,minZero){
    Intercept = offsets = weights = ChIP = mu = phi = zip = Dsg.Mix = muzinb = phizinb = munb = phinb = id = NULL
    derivBackground <- lapply(1,FUN = function(x){
        subdt <- dt[(id == x),][,c('zip','mu','phi') := {
            list(exp(par[1] * Intercept) / (1 + exp(par[1] * Intercept)),
                 exp(par[2] * Intercept + offsets),
                 exp(par[3] * Intercept))
        }]
        colSums(cbind(subdt[,{
            -weights*((ChIP==0)*(((1-stats::dnbinom(x=0,mu=mu,size=phi,log=FALSE))/dzinb(0,mu=mu,size=phi,zip=zip,log=FALSE,minZero = minZero))*zip*(1-zip)) + 
                          (!ChIP==0)*(-zip))*cbind(Intercept)
        }],
        subdt[,{
            -weights*((ChIP==0)*((-(1-zip)*((phi/(mu+phi))^(phi+1))/dzinb(0,mu=mu,size=phi,zip=zip,log=FALSE,minZero = minZero))*mu) + 
                          (!ChIP==0)*((ChIP/mu-(ChIP+phi)/(mu+phi))*mu))*cbind(Intercept)  
        }],
        subdt[,{
            -weights*((ChIP==0)*((1-zip)*((phi/(mu+phi))^phi)*(log(phi/(mu+phi))+mu/(mu+phi))*(1/dzinb(0,mu=mu,size=phi,zip=zip,log=FALSE,minZero = minZero))*phi) + 
                          (!ChIP==0)*((digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi))*cbind(Intercept)
        }],
        0,0))
    })
    
    derivDifferential <- lapply(2:(maxid - 1),FUN = function(x){
        subdt <- dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))]
        subdt <- subdt[,c('zip','muzinb','phizinb','munb','phinb') := {
            list(exp(par[1] * Intercept) / (1 + exp(par[1] * Intercept)),
                 exp(par[2] * Intercept + offsets),
                 exp(par[3] * Intercept),
                 exp(par[4] * Intercept + offsets),
                 exp(par[5] * Intercept))
        }]
        colSums(cbind(subdt[,{
            -weights*(Dsg.Mix==0)*((ChIP==0)*(((1-stats::dnbinom(x=0,mu=muzinb,size=phizinb,log=FALSE))/dzinb(0,mu=muzinb,size=phizinb,zip=zip,log=FALSE,minZero = minZero))*zip*(1-zip)) + 
                                       (!ChIP==0)*(-zip))*cbind(Intercept)
        }],
        subdt[,{
            -weights*(Dsg.Mix==0)*((ChIP==0)*((-(1-zip)*((phizinb/(muzinb+phizinb))^(phizinb+1))/dzinb(0,mu=muzinb,size=phizinb,zip=zip,log=FALSE,minZero = minZero))*muzinb) + 
                                       (!ChIP==0)*((ChIP/muzinb-(ChIP+phizinb)/(muzinb+phizinb))*muzinb))*cbind(Intercept)
        }],
        subdt[,{
            -weights*(Dsg.Mix==0)*((ChIP==0)*((1-zip)*((phizinb/(muzinb+phizinb))^phizinb)*(log(phizinb/(muzinb+phizinb))+muzinb/(muzinb+phizinb))*(1/dzinb(0,mu=muzinb,size=phizinb,zip=zip,log=FALSE,minZero = minZero))*phizinb) + 
                                       (!ChIP==0)*((digamma(ChIP+phizinb)-digamma(phizinb)+1+log(phizinb)-(ChIP+phizinb)/(muzinb+phizinb)-log(muzinb+phizinb))*phizinb))*cbind(Intercept)
        }],
        subdt[,{
            -weights*(Dsg.Mix==1)*(ChIP/munb-((ChIP+phinb)/(munb+phinb)))*munb*cbind(Intercept)
        }],
        subdt[,{
            -weights*(Dsg.Mix==1)*(digamma(ChIP+phinb)-digamma(phinb)+1+log(phinb)-(ChIP+phinb)/(munb+phinb)-log(munb+phinb))*phinb*cbind(Intercept)
        }]))
    })
    derivEnrichment <- lapply(maxid,FUN = function(x){
        subdt <- dt[(id == x),][,c('mu','phi') := {
            list(exp(par[4]*Intercept+offsets),
                 exp(par[5]*Intercept))
        }]
        colSums(cbind(0,0,0,
                      subdt[,{
                          -weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept)
                      }],
                      subdt[,{
                          -weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept)
                      }]))
    })
    
    return(colSums(do.call(rbind,c(derivBackground,derivDifferential,
                                   derivEnrichment))))
}
