################################################################################
################################################################################
### R script w/ the log-likelihood functions used in epigraHMM
################################################################################
################################################################################

################################################################################
### Log-likelihood function of a mixture of NB with constrained parameters
################################################################################

loglikDifferential = function(dt,delta,psi,N,M,dist,minZero){
    Dsg.Mix = ChIP = offsets = .SD = NULL
    if (dist == 'nb') {
        ll <- lapply(seq_len(length(delta)),FUN = function(b){
            log(delta[b]) + rowSums(matrix(dt[,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',b)][,{
                stats::dnbinom(x = ChIP,
                               mu = exp(psi[1] + psi[3] * Dsg.Mix + offsets),
                               size = exp(psi[2] + psi[4] * Dsg.Mix),
                               log = TRUE)
            }], nrow = M, ncol = N, byrow = FALSE))
        })
        dt[,Dsg.Mix := NULL]
    } else{
        ll <- lapply(seq_len(length(delta)),FUN = function(b){
            log(delta[b]) + rowSums(matrix(dt[,{
                (.SD == 1)*stats::dnbinom(x = ChIP,
                                        mu = exp(psi[2] + psi[4] + offsets),
                                        size = exp(psi[3] + psi[5]),
                                        log = TRUE) + 
                    (.SD == 0)*dzinb(x = ChIP,mu = exp(psi[2] + offsets),
                                   size = exp(psi[3]),
                                   zip = exp(psi[1])/(1 + exp(psi[1])),
                                   minZero = minZero) 
            },.SDcols = paste0('Dsg.Mix',b)],
            nrow = M, ncol = N, byrow = FALSE))
        })
    }
    return(ll)
}

################################################################################
### Log-likelihood function with mixtures sharing unique set of parameters
################################################################################

loglikDifferentialHMM = function(dt,theta,N,M,minZero,dist){
    ChIP = offsets = NULL
    psi <- unlist(theta$psi)
    delta <- theta$delta
    LL <- matrix(0, nrow = M, ncol = 3)
    if (dist == 'nb') {
        #LL from Background
        LL[,1] <- rowSums(matrix(dt[,{
            stats::dnbinom(x = ChIP,
                           mu = exp(psi[1] + offsets),
                           size = exp(psi[2]),
                           log = TRUE)
        }], nrow = M, ncol = N, byrow = FALSE))
        #LL from Mixture
        LL[,2] <- log(base::Reduce(`+`,{
            llDiferential <- 
                loglikDifferential(dt = dt,delta = delta,psi = psi,N = N,
                                   M = M,dist = dist,minZero = minZero)
            lapply(llDiferential,FUN = exp)
        }))
        #LL from Enrichment
        LL[,3] <- rowSums(matrix(dt[,{
            stats::dnbinom(x = ChIP,
                           mu = exp(psi[1] + psi[3] + offsets),
                           size = exp(psi[2] + psi[4]),
                           log = TRUE)
        }], nrow = M, ncol = N, byrow = FALSE))
    } else{
        #LL from Background
        LL[,1] = rowSums(matrix(dt[,{
            dzinb(x = ChIP,
                  mu = exp(psi[2] + offsets),
                  size = exp(psi[3]),
                  zip = exp(psi[1]) / (1 + exp(psi[1])),
                  minZero = minZero)
        }], nrow = M, ncol = N, byrow = FALSE))
        #LL from Mixture
        LL[,2] <- log(base::Reduce(`+`,{
            llDifferential <-
                loglikDifferential(dt = dt,delta = delta,psi = psi,N = N,
                                   M = M,dist = dist,minZero = minZero)
            lapply(llDifferential,FUN = exp) 
        }))
        #LL from Enrichment
        LL[,3] = rowSums(matrix(dt[,{
            stats::dnbinom(x = ChIP,
                           mu = exp(psi[2] + psi[4] + offsets),
                           size = exp(psi[3] + psi[5]),
                           log = TRUE)   
        }], nrow = M, ncol = N, byrow = FALSE))
    }
    LL[is.infinite(LL)] <- log(minZero)
    return(LL)
}

################################################################################
### Log-likelihood function of HMM for consensus peak calling
################################################################################

loglikConsensusHMM = function(dt,theta,N,M,minZero,dist){
    ChIP = offsets = controls = yvec0 = zip = ll00 = ll01 = .SD = NULL
    LL <- matrix(0, nrow = M, ncol = 2)
    useControl <- ('controls' %in% names(dt))
    theta1 <- theta$psi[[1]]
    theta2 <- theta$psi[[2]]
    if (dist == 'nb') {
        #LL from Background
        LL[,1] <- rowSums(matrix(dt[,{
            mu <- exp(theta1[1]+ifelse(useControl,theta1[2]*controls,0)+offsets)
            size <- ifelse(useControl,theta1[3],theta1[2])
            log(stats::dnbinom(x = ChIP,mu = mu,size = size,log=FALSE)+minZero)
        }], nrow = M, ncol = N, byrow = FALSE))
        #LL from Enrichment
        LL[,2] <- rowSums(matrix(dt[,{
            mu <- exp(theta2[1]+ifelse(useControl,theta2[2]*controls,0)+offsets)
            size <- ifelse(useControl,theta2[3],theta2[2])
            log(stats::dnbinom(x = ChIP,mu = mu,size = size,log=FALSE)+minZero)
        }], nrow = M, ncol = N, byrow = FALSE))
    } else{
        #LL from Background
        LL[,1] <- rowSums(matrix(dt[,{
            xbeta <- ifelse(useControl,theta1[1] + theta1[2]*controls,theta1[1])+offsets
            mu00 <- exp(ifelse(useControl,theta1[3] + theta1[4]*controls,theta1[2]) + offsets)
            size00 <- ifelse(useControl,theta1[5],theta1[3])
            mu01 <- exp(ifelse(useControl,theta1[3] + theta1[4]*controls,theta1[2]) + offsets)
            size01 <- ifelse(useControl,theta1[5],theta1[3])
            ll00 <- log(stats::dnbinom(x = 0,mu = mu00,size = size00,log=FALSE)+minZero)
            ll01 <- log(stats::dnbinom(x = ChIP,mu = mu01,size = size01,log=FALSE)+minZero)
            cbind('zip' = 1/(1+exp(-xbeta)),'yvec0' = 1*(ChIP == 0),'ll00' = ll00,'ll01' = ll01,.SD)
        }][,yvec0*log(zip+exp(log(1-zip)+ll00)) + (1-yvec0)*(log(1-zip)+ll01)],nrow=M,ncol=N,byrow=FALSE))
        #LL from Enrichment
        LL[,2] <- rowSums(matrix(dt[,{
            mu <- exp(ifelse(useControl,theta2[1]+theta2[2]*controls,theta2[1])+offsets)
            size <- ifelse(useControl,theta2[3],theta2[2])
            log(stats::dnbinom(x = ChIP,mu = mu,size = size,log=FALSE)+minZero)
        }], nrow = M, ncol = N, byrow = FALSE))
    }
    return(LL)
}

################################################################################
### Function to calculate the log-transformed probability dist. of a ZINB model
### It agrees with the function VGAM::dzinegbin, except that a small quantity
### is added to stats::dnbinom to avoid having -Inf when passing this function
### to optim-'L-BFGS-B'.
################################################################################

dzinb = function(x, size, mu, zip, log = TRUE, minZero) {
    if (log == TRUE) {
        log(zip * (x == 0) + (1 - zip) * (stats::dnbinom(x,mu = mu,size = size,log = FALSE)+minZero))
    } else{
        zip*(x==0)+(1-zip)*(stats::dnbinom(x,mu = mu,size = size,log = FALSE)+minZero)
    }
}