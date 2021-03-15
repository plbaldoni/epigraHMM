################################################################################
### Estimate transition probability from a sequence of integers
################################################################################

estimateTransitionProb = function(chain,numStates){
    if(max(chain)>(numStates-1)){
        chain=round((numStates-1)*(chain-max(chain))/(max(chain)-min(chain))+(numStates-1))
    }
    MC = matrix(chain,nrow = 1,ncol=length(chain))
    MC = table(c(MC[,-ncol(MC)]),c(MC[,-1]))
    MC = as.matrix(MC/rowSums(MC))
    MC = matrix(MC,ncol=ncol(MC),nrow=nrow(MC),byrow=FALSE)
    return(checkProbabilities(MC))
}

################################################################################
### Check if rows of a matrix sum up to 1
################################################################################

checkProbabilities = function(P){
    P = pmax(pmin(P,1),0)
    if(sum(P==0)>0){P[P==0] = .Machine$double.xmin}
    return(P/rowSums(P))
}

################################################################################
### Log-likelihood function of a GLM NB distribution and its derivative
################################################################################

glmNB = function(par,y,x,offset,weights){
    l = stats::dnbinom(y,mu=exp(x%*%par[seq_len(ncol(x))]+offset),size=1/par[length(par)],log=TRUE);l[is.infinite(l)] = log(1e-300)
    return(-sum(weights*l))
}

derivNB = function(par,y,x,offset,weights){
    mu.vec = exp(x%*%par[seq_len(ncol(x))]+offset)
    return(-c(colSums(as.numeric(weights)*as.numeric((y-mu.vec)/(1+par[length(par)]*mu.vec))*x),
              sum(as.numeric(weights)*(log(1+par[length(par)]*mu.vec)+par[length(par)]*(y-mu.vec)/(1+par[length(par)]*mu.vec)-digamma(y+1/par[length(par)])+digamma(1/par[length(par)]))/(par[length(par)]^2))))
}

################################################################################
### Log-likelihood function of a GLM ZINB distribution and its derivative
################################################################################

glmZINB = function(par,y,x,offset,weights,minZero){
    mu.vec <- exp(x%*%par[seq_len(ncol(x))]+offset)
    zeroinfl.vec <- 1/(1+exp(-(x%*%par[(ncol(x)+2):(2*ncol(x)+1)]+offset)))
    idx.Y0 <- (y==0)
    
    return(-sum(weights*(idx.Y0*log(zeroinfl.vec+exp(log(1-zeroinfl.vec)+log(stats::dnbinom(0,mu=mu.vec,size=1/par[(ncol(x)+1)],log=FALSE)+minZero)))+(1-idx.Y0)*(log(1-zeroinfl.vec)+log(stats::dnbinom(y,mu=mu.vec,size=1/par[(ncol(x)+1)],log=FALSE)+minZero)))))
}

derivZINB = function(par,y,x,offset,weights,minZero){
    # This function is correct (Checked with grad() from NumDeriv on 11/03/18)
    mu.vec <- exp(x%*%par[seq_len(ncol(x))]+offset)
    zeroinfl.vec <- 1/(1+exp(-(x%*%par[(ncol(x)+2):(2*ncol(x)+1)]+offset)))
    idx.Y0 <- (y==0)
    
    l0 = log(stats::dnbinom(0,mu=mu.vec,size=1/par[(ncol(x)+1)],log=FALSE)+minZero)
    l1 = log(stats::dnbinom(y,mu=mu.vec,size=1/par[(ncol(x)+1)],log=FALSE)+minZero)
    
    P0 = zeroinfl.vec+exp(log(1-zeroinfl.vec)+l0)
    P1 = exp(log(1-zeroinfl.vec)+l1)
    
    phi = par[(ncol(x)+1)]
    aux = (1+phi*mu.vec)
    
    return(c(-colSums(as.numeric((idx.Y0)*(weights/P0)*(1-zeroinfl.vec)*(-mu.vec*((1+phi*mu.vec)^(-(1/phi+1)))))*x+as.numeric((1-idx.Y0)*(weights/P1)*(1-zeroinfl.vec)*exp(l1)*(y-mu.vec)/(1+phi*mu.vec))*x),
             -sum((idx.Y0)*(weights/P0)*(1-zeroinfl.vec)*(aux^(-1/phi))*(log(aux)/((phi)^2)-mu.vec/(phi*aux))+(1-idx.Y0)*(weights/P1)*(1-zeroinfl.vec)*exp(l1)*(log(aux)/(phi^2)-mu.vec*(y+1/phi)/aux-digamma(y+1/phi)/(phi^2)+digamma(1/phi)/(phi^2)+y/phi)),
             -colSums(as.numeric((idx.Y0)*(weights/P0)*(1-exp(l0))*zeroinfl.vec*(1-zeroinfl.vec))*x+as.numeric((1-idx.Y0)*(weights/P1)*(-exp(l1))*zeroinfl.vec*(1-zeroinfl.vec))*x)))
}

################################################################################
### Invert parameters for glm.nb and glm.zinb
################################################################################

invertDispersion = function(par,model){
    if(model=='nb'){return(c(par[seq_len((length(par)-1))],1/par[length(par)]))}
    if(model=='zinb'){
        if(length(par)==3){
            return(c(par[2],1/par[3],par[1]))
        } else{
            return(c(par[c(3,4)],1/par[5],par[c(1,2)]))
        }
    }
}

################################################################################
### Print message during EM algorithm
################################################################################

verbose = function(level,control,controlHist = NULL,theta = NULL){
    currentTime <- Sys.time()
    if(!control[['quiet']]){
        if(level == 1){
            message(paste0(c(rep('#',80))));message(currentTime);message("Starting the EM algorithm")
            return(currentTime)
        }
        if(level == 2){
            message(paste0(c(rep('#',80))))
            message('\rIteration: ',controlHist[[length(controlHist)]][['iteration']])
            message('\rBIC: ',formatC(controlHist[[length(controlHist)]][['bic']], format = "e", digits = 2))
            message("\rError: ",paste(paste(names(controlHist[[length(controlHist)]][['error']]),'='),formatC(controlHist[[length(controlHist)]][['error']], format = "e", digits = 2),collapse = ', '))
            message("\rConvergence: ",paste(paste(names(controlHist[[length(controlHist)]][['count']]),'='),controlHist[[length(controlHist)]][['count']],collapse = ', '))
            message("\r",paste('Initial probabilities: '),paste(formatC(theta[['pi']], format = "e", digits = 2),collapse = ' '))
            message("\r",paste('Transition probabilities: '),paste(formatC(theta[['gamma']], format = "e", digits = 2),collapse = ' '))
            if('delta' %in% names(theta)){
                message("\r",paste('Mixture probabilities: '),paste(formatC(unlist(theta[['delta']]), format = "e", digits = 2),collapse = ' '))
            }
            message("\r",paste('Model coefficients: '),paste(formatC(unlist(theta[['psi']]), format = "e", digits = 2),collapse = ' '))
            message("\rTime elapsed (mins): ",formatC(controlHist[[length(controlHist)]][['time']], format = "f", digits = 2))
        }
        if(level == 3){
            message(paste0(c(rep('#',80))))
            message("EM algorithm converged!")
            message(currentTime)
            message(paste0(c(rep('#',80))))
        }   
    }
}

################################################################################
### Compute error for EM algorithm convergence check
################################################################################

getError <- function(controlHist,parHist,control){
    
    iteration <- controlHist[[length(controlHist)]][['iteration']]
    
    if(iteration>1){
        gap <- ifelse(iteration>control[['minIterEM']],control[['gapIterEM']],1)
        
        parlist.old <- unlist(parHist[[iteration - gap]])
        parlist.new <- unlist(parHist[[iteration]])
        
        q.old <- controlHist[[iteration - gap]][['q']]
        q.new <- controlHist[[iteration]][['q']]
        
        return(c('MRCPE' = max(abs((parlist.new-parlist.old)/parlist.old)),
                 'MACPE' = max(abs(parlist.new-parlist.old)),
                 'ARCEL' = max(abs((q.new-q.old)/q.old))))
    } else{
        return(c('MRCPE' = 0,'MACPE' = 0,'ARCEL' = 0))
    }
}

################################################################################
### Optimizer for glmNB
################################################################################

optimNB <- function(par,y,x,offset,weights,control,dist){
    
    if(dist == 'nb'){
        tryCatch({assign('model',stats::optim(par=invertDispersion(par,model=dist),
                                              fn=glmNB,
                                              gr=derivNB,
                                              method='L-BFGS-B',
                                              lower=c(rep(-Inf,ncol(x)),1/control$maxDisp),
                                              y = y,
                                              x = x,
                                              offset = offset,
                                              weights = weights))},
                 error=function(e){assign('model',list('par' = invertDispersion(par,model=dist),'convergence' = 99),inherits = TRUE)})
        
        model$par <- invertDispersion(model$par,model=dist)   
    } else{
        tryCatch({assign('model',stats::optim(par = invertDispersion(par,model=dist),
                                              fn = glmZINB,
                                              gr = derivZINB,
                                              method = 'L-BFGS-B',
                                              lower = c(rep(-Inf,ncol(x)),1/control$maxDisp,rep(-Inf,ncol(x))),
                                              y = y,
                                              x = x,
                                              offset = offset,
                                              weights = weights,
                                              minZero = control$minZero))},
                 error=function(e){assign('model',list('par' = invertDispersion(par,model=dist), 'convergence' = 99),inherits = TRUE)})
        if(ncol(x)==1){
            model$par <- c(model$par[c(3,1)],1/model$par[2])
        } else{
            model$par <- c(model$par[c(4,5)],model$par[c(1,2)],1/model$par[3])
        }
    }
    
    return(model)
}

################################################################################
### Function to enumerate combinatorial patterns from initializing peaks
################################################################################

enumeratePatterns = function(object,group){
    
    Window = ref = NULL
    
    # Peaks by group
    chain <- data.table::setDT(lapply(unique(group),function(x){1*(Matrix::rowSums(SummarizedExperiment::assay(object,'peaks')[,which(group==x),drop = FALSE])>0)}))
    data.table::setnames(chain,paste0('ChIP',seq_len(ncol(chain))))
    chain[,Window := seq_len(.N)]
    
    # Creating reference table
    ref[paste0('ChIP',seq_len(length(unique(group))))] <- list(NULL)
    for(i in names(ref)){ref[[i]] <- c(0,1)}
    ref <- as.data.table(expand.grid(ref))
    ref <- ref[order(rowSums(ref)),]
    ref$Group <- seq_len(nrow(ref))
    
    chain <- merge(chain,ref,by=paste0('ChIP',seq_len(length(unique(group)))), all.x = TRUE)
    return(list('group' = chain[order(Window)][,c('Group'), with = FALSE],
                'ref' = ref))
}

################################################################################
### Function to dassociate mixture components with combinatorial patterns
################################################################################

determineMixtures <- function(pattern,group,nGroup,ref){
    if(is.list(pattern)){
        B <- length(pattern)
        z.seq <- lapply(seq_len(B),FUN = function(x){
            aux <- rep(FALSE,nGroup);aux[pattern[[x]]]<-TRUE
            which(apply(as.matrix(ref[,seq_len(nGroup),with = FALSE])==1,1,FUN = function(x){all(x==aux)}))
        })
    } else{
        if(is.null(pattern)){
            B <- 2^(length(unique(group)))-2
            z.seq <- lapply(1+seq_len(B),FUN = function(x){x})
        } else {
            stop('Error: pattern should be either NULL or a list. Check controlEM()')
        }
    }
    return(list('B' = B,'zSeq' = z.seq))
}

################################################################################
### Generate differential (mixture) intercepts
################################################################################

generateMixtureIntercepts <- function(pattern,group,nGroup,ref,M,B){
    if(is.list(pattern)){
        model <- lapply(unlist(lapply(pattern,FUN = function(x){
            aux <- rep(FALSE,nGroup);aux[x]<-TRUE
            which(apply(as.matrix(ref[,seq_len(nGroup),with = FALSE])==1,1,FUN = function(x){all(x==aux)}))
        })),FUN = function(x){cbind(rep(t(ref[x,seq_len(nGroup),with = FALSE]),times=M*table(group)))})
    }
    if(is.null(pattern)){
        model <- lapply(1+seq_len(B),FUN = function(x) cbind(rep(t(ref[x,seq_len(nGroup),with = FALSE]),times=M*table(group))))
    }
    return(model)
}

################################################################################
### Estimate mixture probabilities for initialization
################################################################################

estimateMixtureProb = function(zDifferential,zSeq,B){
    return(unlist(lapply(seq_len(B),FUN = function(i){
        subz.diff <- zDifferential[zDifferential%in%unlist(zSeq)]
        sum(subz.diff%in%zSeq[[i]]/length(subz.diff))
    })))
}

################################################################################
### Estimate model coefficients for differential model initialization
################################################################################

estimateCoefficients <- function(z,dt,dist,type,control){
    Window = ChIP = offsets = mu = sigma2 = NULL
    
    parList <- lapply(range(z),function(x){
        subpar <- dt[Window %in% which(z == x),list(mu = mean((ChIP+1)/exp(offsets)),sigma2 = stats::var((ChIP+1)/exp(offsets)))]
        subpar[,c(zip = (sigma2-mu)/(sigma2+mu^2-mu),mu = mu,disp = min((mu^2)/max(0,sigma2-mu),control[['maxDisp']]))]
    })
    
    par <- list()
    if(type == 'differential'){
        if(dist == 'nb'){
            par[[1]] <- log(parList[[1]][-1])
            par[[2]] <- log(parList[[2]][-1]) - log(parList[[1]][-1])
        } else{
            par[[1]] <- log(parList[[1]])
            par[[2]] <- log(parList[[2]][-1]) - log(parList[[1]][-1])
        }
    } else{
        if(dist == 'nb'){
            par[[1]] <- c(log(parList[[1]]['mu']),ifelse('controls' %in% names(dt),0,NA),parList[[1]]['disp'])
            par[[2]] <- c(log(parList[[2]]['mu']),ifelse('controls' %in% names(dt),0,NA),parList[[2]]['disp'])
        } else{
            par[[1]] <- c(log(parList[[1]]['zip']),ifelse('controls' %in% names(dt),0,NA),
                          log(parList[[1]]['mu']),ifelse('controls' %in% names(dt),0,NA),parList[[1]]['disp'])
            par[[2]] <- c(log(parList[[2]]['mu']),ifelse('controls' %in% names(dt),0,NA),parList[[2]]['disp'])
        }
    }
    
    return(lapply(par,function(x){unname(x[!is.na(x)])}))
}

################################################################################
### Function for log-likelihood function of mixture of NB
### components sharing unique set of parameters
################################################################################

loglikDifferential = function(dt,delta,psi,N,M,dist,minZero){
    Dsg.Mix = ChIP = offsets = NULL
    
    if(dist == 'nb'){
        ll <- lapply(seq_len(length(delta)),FUN = function(b){
            log(delta[b])+rowSums(matrix(dt[,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',b)][,stats::dnbinom(x = ChIP,mu = exp(psi[1]+psi[3]*Dsg.Mix+offsets),size = exp(psi[2]+psi[4]*Dsg.Mix),log=TRUE)],nrow=M,ncol=N,byrow=FALSE))
        })
        dt[,Dsg.Mix := NULL]
    } else{
        ll <- lapply(seq_len(length(delta)),FUN = function(b){
            log(delta[b])+rowSums(matrix(dt[,((.SD==1)*stats::dnbinom(ChIP,mu = exp(psi[2]+psi[4]+offsets),size = exp(psi[3]+psi[5]),log = TRUE) + (.SD==0)*dzinb(x = ChIP,mu = exp(psi[2]+offsets),size = exp(psi[3]),zip = exp(psi[1])/(1+exp(psi[1])),minZero = minZero)),.SDcols = paste0('Dsg.Mix',b)],nrow=M,ncol=N,byrow=FALSE))
        })
    }
    
    return(ll)
}

################################################################################
### Function to calculate the log-transformed probability dist. of a ZINB model
### It agrees with the function VGAM::dzinegbin, except that I add a small quantity to
### stats::dnbinom to avoid having -Inf when passing this function to optim-'L-BFGS-B'.
################################################################################

dzinb = function(x,size,mu,zip,log=TRUE,minZero){
    if(log==TRUE){
        return(log(zip*(x==0)+(1-zip)*(stats::dnbinom(x,mu = mu,size = size,log = FALSE)+minZero)))
    } else{
        return(zip*(x==0)+(1-zip)*(stats::dnbinom(x,mu = mu,size = size,log = FALSE)+minZero))
    }
}

################################################################################
### Function for log-likelihood function of HMM with NB distributions with mixture
### components sharing unique set of parameters
################################################################################

loglikDifferentialHMM = function(dt,theta,N,M,minZero,dist){
    
    ChIP = offsets = Dsg.Mix = NULL
    
    psi <- unlist(theta$psi)
    delta <- theta$delta
    B <- length(delta)
    LL <- matrix(0,nrow=M,ncol=3)
    
    if(dist == 'nb'){
        #LL from Background
        LL[,1] <- rowSums(matrix(dt[,stats::dnbinom(x = ChIP,mu = exp(psi[1]+offsets),size = exp(psi[2]),log=TRUE)],nrow=M,ncol=N,byrow=FALSE))
        #LL from Mixture
        LL[,2] <- log(base::Reduce(`+`,lapply(loglikDifferential(dt = dt,delta = delta,psi = psi,N = N,M = M,dist = dist,minZero = minZero),exp)))
        #LL from Enrichment
        LL[,3] <- rowSums(matrix(dt[,stats::dnbinom(x = ChIP,mu = exp((psi[1]+psi[3])+offsets),size = exp((psi[2]+psi[4])),log=TRUE)],nrow=M,ncol=N,byrow=FALSE))
    } else{
        #LL from Background
        LL[,1] = rowSums(matrix(dt[,dzinb(x = ChIP,mu = exp(psi[2]+offsets),size = exp(psi[3]),zip = exp(psi[1])/(1+exp(psi[1])),minZero = minZero)],nrow=M,ncol=N,byrow=FALSE))
        #LL from Mixture
        LL[,2] <- log(base::Reduce(`+`,lapply(loglikDifferential(dt = dt,delta = delta,psi = psi,N = N,M = M,dist = dist,minZero = minZero),exp)))
        #LL from Enrichment
        LL[,3] = rowSums(matrix(dt[,stats::dnbinom(x = ChIP,mu = exp(psi[2]+psi[4]+offsets),size = exp(psi[3] + psi[5]),log = TRUE)],nrow=M,ncol=N,byrow=FALSE))
    }
    
    LL[is.infinite(LL)] <- log(minZero)
    return(LL)
}

################################################################################
### Function for log-likelihood function of HMM with NB distributions for
### consensus peak calling
################################################################################

loglikConsensusHMM = function(dt,theta,N,M,minZero,dist){
    
    ChIP = offsets = controls = yvec0 = zip = ll00 = ll01 = NULL
    
    psi <- unlist(theta$psi)
    LL <- matrix(0,nrow=M,ncol=2)
    useControl <- ('controls' %in% names(dt))
    
    if(dist == 'nb'){
        #LL from Background
        LL[,1] <- rowSums(matrix(dt[,log(stats::dnbinom(x = ChIP,mu = exp(theta$psi[[1]][1]+ifelse(useControl,theta$psi[[1]][2]*controls,0)+offsets),size = ifelse(useControl,theta$psi[[1]][3],theta$psi[[1]][2]),log=FALSE)+minZero)],nrow=M,ncol=N,byrow=FALSE))
        #LL from Enrichment
        LL[,2] <- rowSums(matrix(dt[,log(stats::dnbinom(x = ChIP,mu = exp(theta$psi[[2]][1]+ifelse(useControl,theta$psi[[2]][2]*controls,0)+offsets),size = ifelse(useControl,theta$psi[[2]][3],theta$psi[[2]][2]),log=FALSE)+minZero)],nrow=M,ncol=N,byrow=FALSE))
    } else{
        #LL from Background
        LL[,1] = rowSums(matrix(dt[,cbind('zip' = 1/(1+exp(-(ifelse(useControl,theta$psi[[1]][1] + theta$psi[[1]][2]*controls,theta$psi[[1]][1])+offsets))),
                                          'yvec0' = 1*(ChIP == 0),
                                          'll00' = log(stats::dnbinom(0,mu = exp(ifelse(useControl,theta$psi[[1]][3] + theta$psi[[1]][4]*controls,theta$psi[[1]][2]) + offsets),size = ifelse(useControl,theta$psi[[1]][5],theta$psi[[1]][3]),log=FALSE)+minZero),
                                          'll01' = log(stats::dnbinom(ChIP,mu = exp(ifelse(useControl,theta$psi[[1]][3] + theta$psi[[1]][4]*controls,theta$psi[[1]][2]) + offsets),size = ifelse(useControl,theta$psi[[1]][5],theta$psi[[1]][3]),log=FALSE)+minZero),.SD)][,yvec0*log(zip+exp(log(1-zip)+ll00)) + (1-yvec0)*(log(1-zip)+ll01)],nrow=M,ncol=N,byrow=FALSE))
        
        #LL from Enrichment
        LL[,2] = rowSums(matrix(dt[,log(stats::dnbinom(x = ChIP,mu = exp(ifelse(useControl,theta$psi[[2]][1]+theta$psi[[2]][2]*controls,theta$psi[[2]][1])+offsets),size = ifelse(useControl,theta$psi[[2]][3],theta$psi[[2]][2]),log=FALSE)+minZero)],nrow=M,ncol=N,byrow=FALSE))
    }
    
    return(LL)
}

################################################################################
### Function that implements the inner E-step of the EM algorithm for differential
### calling
################################################################################

innerExpStep = function(dt,theta,N,M,model,hdf5,minZero,dist){
    
    ll <- loglikDifferential(dt = dt,delta = theta$delta,psi = unlist(theta$psi),N = N,M = M,dist = dist,minZero = minZero)
    ll <- rapply(ll,function(x){ifelse(exp(x)==0,log(minZero),x)}, how = "replace")
    sumExpLL <- Reduce(`+`,lapply(ll,exp))
    
    rhdf5::h5write(obj = vapply(seq_len(length(theta$delta)),FUN = function(b){exp(ll[[b]])/sumExpLL},FUN.VALUE = vector('double',length = M)),
                   file = hdf5,
                   name = 'mixtureProb')
}

################################################################################
### Log-likelihood function of 3-state HMM with constrained parameters with
### mixture of NB distributions and its derivative
################################################################################

glmMixNB = function(par,dt,maxid,minZero){
    # This function is the -loglikelihood function implements a mixture model in the differential component with NB distributions
    # NB distributions from the mixture model share the same set of parameters for both mean and dispersion models.
    # Each NB distribution represents one out of all possible combinations of differential patterns
    # exp(Int) is the mean (dispersion) of the group with background counts in that particular NB distribution
    # exp(Int + Slope) is the mean (dispersion) of the group with enrichment counts in that particular NB distribution
    # Additionally, the mean (dispersion) for the first and third HMM components are, respectively, exp(Int) and exp(Int + Slope).
    # Hence, only 4 parameters drive the entire distributions of this HMM
    Intercept = weights = ChIP = Dsg.Mix = offsets = id = NULL
    
    return(sum(dt[(id==1),][,-weights*log(stats::dnbinom(ChIP,mu = exp(par[1]*Intercept+offsets),size = exp(par[3]*Intercept),log = FALSE)+minZero)])+
               do.call(sum,lapply(2:(maxid-1),FUN = function(x){sum(dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,-weights*log(stats::dnbinom(ChIP,mu = exp(par[1]*Intercept+par[2]*Dsg.Mix+offsets),size = exp(par[3]*Intercept+par[4]*Dsg.Mix),log = FALSE)+minZero)])}))+
               sum(dt[(id==maxid),][,-weights*log(stats::dnbinom(ChIP,mu = exp((par[1]+par[2])*Intercept+offsets),size = exp((par[3]+par[4])*Intercept),log = FALSE)+minZero)]))
}

derivMixNB = function(par,dt,maxid,minZero){
    # This function is the -loglikelihood of the mixHMMConstr that implements a mixture model in the differential component with NB distributions
    # NB distributions from the mixture model share the same set of parameters for both mean and dispersion models.
    # Each NB distribution represents one out of all possible combinations of differential patterns
    # exp(Int) is the mean (dispersion) of the group with background counts in that particular NB distribution
    # exp(Int + Slope) is the mean (dispersion) of the group with enrichment counts in that particular NB distribution
    # Additionally, the mean (dispersion) for the first and third HMM components are, respectively, exp(Int) and exp(Int + Slope).
    # Hence, only 4 parameters drive the entire distributions of this HMM
    
    # This derivative was checked with numDeriv::grad and it is correct for glmMixNB (April 10th 2019)
    Intercept = offsets = weights = ChIP = mu = phi = Dsg.Mix = id = NULL
    
    return(colSums(do.call(rbind,c(lapply(1,FUN = function(x){
        subdt = dt[(id == x),][,c('mu','phi') := list(exp(par[1]*Intercept+offsets),exp(par[3]*Intercept))]
        colSums(cbind(subdt[,-weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept,0)],
                      subdt[,-weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept,0)]))
    }),
    lapply(2:(maxid-1),FUN = function(x){
        subdt = dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,c('mu','phi') := list(exp(par[1]*Intercept+par[2]*Dsg.Mix+offsets),exp(par[3]*Intercept+par[4]*Dsg.Mix))]
        colSums(cbind(subdt[,-weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept,Dsg.Mix)],
                      subdt[,-weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept,Dsg.Mix)]))
    }),
    lapply(maxid,FUN = function(x){
        subdt = dt[(id == x),][,c('mu','phi') := list(exp((par[1]+par[2])*Intercept+offsets),exp((par[3]+par[4])*Intercept))]
        colSums(cbind(subdt[,-weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept,Intercept)],
                      subdt[,-weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept,Intercept)]))
    })))))
}

################################################################################
### Log-likelihood function of 3-state HMM with constrained parameters with
### mixture of NB distributions and its derivative
################################################################################

glmMixZINB = function(par,dt,maxid,minZero){
    # This function is the -loglikelihood of the mixZIHMMConstr that implements a mixture model in the differential component with ZINB and NB distributions
    # ZINB and NB distributions from the mixture model share the same set of parameters for both mean and dispersion models.
    # Each NB distribution in the differential component represents the upward shift in the background counts from ZINB model
    # exp(Int) is the mean (dispersion) of the group with background counts in that particular ZINB distribution, which has zero-inflation probability exp(Int)/(1+exp(Int))
    # exp(Int + Slope) is the mean (dispersion) of the group with enrichment counts in that particular NB distribution
    # Additionally, the mean (dispersion) for the first and third HMM components are, respectively, exp(Int) and exp(Int + Slope).
    # Hence, only 5 parameters drive the entire distributions of this HMM
    Dsg.Mix = weights = ChIP = Intercept = offsets = id = NULL
    
    return(sum(dt[(id==1),][,-weights*dzinb(ChIP,mu = exp(par[2]*Intercept+offsets),size = exp(par[3]*Intercept),zip = exp(par[1])/(1+exp(par[1])),minZero = minZero)])+
               do.call(sum,lapply(2:(maxid-1),FUN = function(x){
                   sum(dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,-weights*((Dsg.Mix==1)*log(stats::dnbinom(ChIP,mu = exp(par[4]*Intercept+offsets),size = exp(par[5]*Intercept),log = FALSE)+minZero) + (Dsg.Mix==0)*dzinb(ChIP,mu = exp(par[2]*Intercept+offsets),size = exp(par[3]*Intercept),zip = exp(par[1])/(1+exp(par[1])),minZero = minZero))])
               }))+
               sum(dt[(id==maxid),][,-weights*log(stats::dnbinom(ChIP,mu = exp(par[4]*Intercept+offsets),size = exp(par[5]*Intercept),log = FALSE)+minZero)]))
    
}

derivMixZINB = function(par,dt,maxid,minZero){
    # Derivatives of glm.zinb2_Constr
    # This derivative was checked with numDeriv::grad and it is correct for glm.zinb2_Constr (April 20th 2019)
    Intercept = offsets = weights = ChIP = mu = phi = zip = Dsg.Mix = muzinb = phizinb = munb = phinb = id = NULL
    
    return(colSums(do.call(rbind,c(lapply(1,FUN = function(x){
        subdt = dt[(id == x),][,c('zip','mu','phi') := list(exp(par[1]*Intercept)/(1+exp(par[1]*Intercept)),exp(par[2]*Intercept+offsets),exp(par[3]*Intercept))]
        colSums(cbind(subdt[,-weights*((ChIP==0)*(((1-stats::dnbinom(x=0,mu=mu,size=phi,log=FALSE))/dzinb(0,mu=mu,size=phi,zip=zip,log=FALSE,minZero = minZero))*zip*(1-zip))+(!ChIP==0)*(-zip))*cbind(Intercept)], #zip
                      subdt[,-weights*((ChIP==0)*((-(1-zip)*((phi/(mu+phi))^(phi+1))/dzinb(0,mu=mu,size=phi,zip=zip,log=FALSE,minZero = minZero))*mu)+(!ChIP==0)*((ChIP/mu-(ChIP+phi)/(mu+phi))*mu))*cbind(Intercept)], #int mean background
                      subdt[,-weights*((ChIP==0)*((1-zip)*((phi/(mu+phi))^phi)*(log(phi/(mu+phi))+mu/(mu+phi))*(1/dzinb(0,mu=mu,size=phi,zip=zip,log=FALSE,minZero = minZero))*phi)+(!ChIP==0)*((digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi))*cbind(Intercept)],#int disp background
                      0, #int mean enrichment
                      0)) #int disp enrichment
    }),
    lapply(2:(maxid-1),FUN = function(x){
        subdt = dt[(id == x),][,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',(x-1))][,c('zip','muzinb','phizinb','munb','phinb') := list(exp(par[1]*Intercept)/(1+exp(par[1]*Intercept)),exp(par[2]*Intercept+offsets),exp(par[3]*Intercept),exp(par[4]*Intercept+offsets),exp(par[5]*Intercept))]
        colSums(cbind(subdt[,-weights*(Dsg.Mix==0)*((ChIP==0)*(((1-stats::dnbinom(x=0,mu=muzinb,size=phizinb,log=FALSE))/dzinb(0,mu=muzinb,size=phizinb,zip=zip,log=FALSE,minZero = minZero))*zip*(1-zip))+(!ChIP==0)*(-zip))*cbind(Intercept)],
                      subdt[,-weights*(Dsg.Mix==0)*((ChIP==0)*((-(1-zip)*((phizinb/(muzinb+phizinb))^(phizinb+1))/dzinb(0,mu=muzinb,size=phizinb,zip=zip,log=FALSE,minZero = minZero))*muzinb)+(!ChIP==0)*((ChIP/muzinb-(ChIP+phizinb)/(muzinb+phizinb))*muzinb))*cbind(Intercept)],
                      subdt[,-weights*(Dsg.Mix==0)*((ChIP==0)*((1-zip)*((phizinb/(muzinb+phizinb))^phizinb)*(log(phizinb/(muzinb+phizinb))+muzinb/(muzinb+phizinb))*(1/dzinb(0,mu=muzinb,size=phizinb,zip=zip,log=FALSE,minZero = minZero))*phizinb)+(!ChIP==0)*((digamma(ChIP+phizinb)-digamma(phizinb)+1+log(phizinb)-(ChIP+phizinb)/(muzinb+phizinb)-log(muzinb+phizinb))*phizinb))*cbind(Intercept)],
                      subdt[,-weights*(Dsg.Mix==1)*(ChIP/munb-((ChIP+phinb)/(munb+phinb)))*munb*cbind(Intercept)],
                      subdt[,-weights*(Dsg.Mix==1)*(digamma(ChIP+phinb)-digamma(phinb)+1+log(phinb)-(ChIP+phinb)/(munb+phinb)-log(munb+phinb))*phinb*cbind(Intercept)]))
    }),
    lapply(maxid,FUN = function(x){
        subdt = dt[(id == x),][,c('mu','phi') := list(exp(par[4]*Intercept+offsets),exp(par[5]*Intercept))]
        colSums(cbind(0,
                      0,
                      0,
                      subdt[,-weights*(ChIP/mu-((ChIP+phi)/(mu+phi)))*mu*cbind(Intercept)],
                      subdt[,-weights*(digamma(ChIP+phi)-digamma(phi)+1+log(phi)-(ChIP+phi)/(mu+phi)-log(mu+phi))*phi*cbind(Intercept)]))
    })))))
}

################################################################################
### Optimizer for glmNB
################################################################################

optimDifferential <- function(par,rcWeights,control,dist){
    
    if(dist == 'nb'){
        tryCatch({assign('model',stats::optim(par = par[c(1,3,2,4)],
                                              fn = glmMixNB,
                                              gr = derivMixNB,
                                              method = 'L-BFGS-B',
                                              lower = rep(log(control[['minZero']]),4),
                                              upper = c(Inf,Inf,rep(log(control[['maxDisp']]),2)),
                                              dt = rcWeights,
                                              minZero = control[['minZero']],
                                              maxid = length(unique(rcWeights$id))))},
                 error=function(e){assign('model',list('par' = par, 'convergence' = 99),inherits = TRUE)})
        
        model$par <- model$par[c(1,3,2,4)]
    } else{
        newpar <- c(par[c(1,2,3)],par[2]+par[4],par[3]+par[5])
        tryCatch({assign('model',stats::optim(par = newpar,
                                              fn = glmMixZINB,
                                              gr = derivMixZINB,
                                              method = 'L-BFGS-B',
                                              lower = rep(log(control[['minZero']]),5),
                                              upper = c(Inf,rep(c(Inf,log(control[['maxDisp']])),2)),
                                              dt = rcWeights,
                                              minZero = control[['minZero']],
                                              maxid = length(unique(rcWeights$id))))},
                 error=function(e){assign('model',list('par' = par, 'convergence' = 99),inherits = TRUE)})
        
        model$par <- c(model$par[c(1,2,3)],model$par[4]-model$par[2],model$par[5]-model$par[3])
    }
    
    return(model)
}

################################################################################
### Create control list to track EM algorithm
################################################################################

controlList <- function(time = 0,
                        iteration = 0,
                        convergence = 0,
                        bic = 0,
                        q = 0,
                        MRCPE = c('error' = 0,'count' = 0),
                        MACPE = c('error' = 0,'count' = 0),
                        ARCEL = c('error' = 0,'count' = 0)){
    return(list('time' = time,
                'iteration' = iteration,
                'convergence' = convergence,
                'bic' = bic,
                'q' = q,
                'error' = c('MRCPE' = MRCPE[['error']], 'MACPE' = MACPE[['error']], 'ARCEL' = ARCEL[['error']]),
                'count' = c('MRCPE' = MRCPE[['count']], 'MACPE' = MACPE[['count']], 'ARCEL' = ARCEL[['count']])))
}

################################################################################
### Function to rename file output if it already exists
################################################################################

checkPath = function(path){
    if(!file.exists(path)){return(path)}
    i <- 1
    splitPath <- strsplit(path,"\\.")[[1]]
    repeat {
        f = paste0(paste(splitPath[seq_len(length(splitPath))-1],collapse = ''),i,".",splitPath[length(splitPath)])
        if(!file.exists(f)){return(f)}
        i <- i + 1
    }
}

################################################################################
### Function to check counts for convergence
################################################################################

checkConvergence <- function(controlHist,control){
    if(control[['criterion']] == 'all'){
        return(max(controlHist[[length(controlHist)]][['count']]))
    } else{
        return(controlHist[[length(controlHist)]][['count']][[control[['criterion']]]])
    }
}

################################################################################
### FDR control method to call peaks based on posterior probabilities
################################################################################

fdrControl <- function(prob,fdr = 0.05){
    if(!all(prob>=0 & prob<=1)){stop('Posterior probabilities must be between 0 and 1')}
    if(!(fdr>0 & fdr <1)){stop('fdr must be between 0 and 1')}
    
    notpp = FDR = Window = NULL
    
    return(data.table::data.table(notpp = 1-prob,Window = seq_len(length(prob)),key = 'Window')[order(notpp),][,FDR := ((cumsum(notpp)/seq_len(.N))<fdr)][order(Window),]$FDR)
}

################################################################################
### Initial peak caller
################################################################################

initializerHMM = function(object,control){
    
    Group = ChIP = mu = sigma2 = Intercept = offsets = NULL
    
    # Checking input
    if(!ncol(object)==1){
        stop('The initializer supports only one sample at a time')
    }
    
    # Creating subdirectories
    hdf5File <- checkPath(file.path(path.expand(control[['tempDir']]),paste0(control[['fileName']],'.h5')))
    invisible(lapply(hdf5File,function(x){if(!dir.exists(dirname(x))){dir.create(dirname(x),recursive = TRUE)}}))
    
    # General parameters
    M <- nrow(object)
    K <- 2
    controlHist <- list(controlList())
    parHist <- list()
    theta.old <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    theta.new <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    
    # Transforming data into data.table and calculating scores
    dt <- data.table::data.table(ChIP = as.numeric(assay(object)),offsets = as.numeric(assay(object,'offsets')))
    dt[,Group := .GRP,by = c('ChIP','offsets')]
    
    # Creating Unique data.table
    dtUnique <- unique(dt,by='Group')[,c('ChIP','offsets','Group'),with = FALSE]
    setkey(dtUnique,Group)
    
    # Naive split
    score <- dt[,scale(log1p(ChIP/exp(offsets)))]
    z <- as.numeric(cut(score,breaks=c(-Inf,stats::quantile(score,0.75),Inf)))
    
    # Parameter initializations
    theta.old[['pi']] <- c(0.999,0.001)
    theta.old[['gamma']] <- estimateTransitionProb(chain = z,numStates = K)
    theta.old[['psi']] <- lapply(c(1,2),function(x){
        muvar <- dt[which(z == x),list(mu = mean((ChIP+1)/exp(offsets)),sigma2 = stats::var((ChIP+1)/exp(offsets)))]
        muvar[,c(log(mu),min((mu^2)/max(0,sigma2-mu),control[['maxDisp']]))]
    })
    rm(z,score)
    
    # EM algorithm begins
    ## Verbose
    init.time <- verbose(level = 1,control = control)
    
    while(checkConvergence(controlHist,control)<control[['maxCountEM']] & controlHist[[length(controlHist)]][['iteration']]<control[['maxIterEM']]){
        
        # Update iteration
        controlHist[[length(controlHist)]][['iteration']] = controlHist[[length(controlHist)]][['iteration']] + 1
        
        # E-step
        expStep(pi = theta.old[['pi']],
                gamma = theta.old[['gamma']],
                logf = do.call(cbind,lapply(c(1,2),function(x){stats::dnbinom(x = assay(object,'counts'),mu = exp(theta.old[['psi']][[x]][1]+as.numeric(assay(object,'offsets'))),size = theta.old[['psi']][[x]][2],log = TRUE)})),
                hdf5 = hdf5File)
        
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(hdf5 = hdf5File)
        
        ## Model parameters
        ### Rejection-controlled posterior probabilities
        probCut <-  (control[['probCut']]>0)*max(c(0.9^controlHist[[length(controlHist)]][['iteration']],control[['probCut']]))
        rcWeights <- lapply(consensusRejectionControlled(hdf5 = hdf5File,f = dt$Group,p = probCut)[,1],
                            function(x){colnames(x) <- c('weights','Group');return(as.matrix(cbind(x,dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1])))})
        
        ### Calculating MLEs
        optimMLE <- lapply(seq_len(length(rcWeights)),function(x){
            optimNB(par = theta.old$psi[[x]],
                    y = rcWeights[[x]][,'ChIP'],
                    x = rcWeights[[x]][,'Intercept',drop = FALSE],
                    offset = rcWeights[[x]][,'offsets'], 
                    weights = rcWeights[[x]][,'weights'],
                    control = control,
                    dist = 'nb')
        })
        theta.new[['psi']] <- lapply(optimMLE,function(x){x[['par']]})
        
        # Updating parHist
        parHist[[controlHist[[length(controlHist)]][['iteration']]]] <- theta.new
        
        # Checking convergence
        ## Q-function (evaluated on the old parameters to avoid calculating the log-likelihood again)
        controlHist[[length(controlHist)]][['time']] <- difftime(time1 = Sys.time(),time2 = init.time,units = 'mins') 
        controlHist[[length(controlHist)]][['bic']] <- computeBIC(hdf5 = hdf5File,
                                                                  numPar = length(unlist(theta.old)) - 3,
                                                                  numSamples = ncol(object))
        controlHist[[length(controlHist)]][['q']] <- computeQFunction(hdf5 = hdf5File,pi = theta.old[['pi']],gamma = theta.old[['gamma']])
        controlHist[[length(controlHist)]][['convergence']] <-  1*any(!unname(unlist(lapply(optimMLE,function(x){x[['convergence']]}))) == 0)
        controlHist[[length(controlHist)]][['error']] <- getError(controlHist = controlHist,parHist = parHist,control = control)
        controlHist[[length(controlHist)]][['count']] <- (controlHist[[length(controlHist)]][['error']]<=control[['epsilonEM']])*(controlHist[[length(controlHist)]][['iteration']]>control[['minIterEM']])*(controlHist[[length(controlHist)]][['count']]+1) + 0
        
        # Updating controlHist
        controlHist[[length(controlHist)+1]] <- controlHist[[length(controlHist)]]
        
        # Updating parameter history
        theta.old <- theta.new
        
        #Verbose
        verbose(level = 2,control = control,controlHist = controlHist,theta = theta.new)
    }
    
    # Verbose
    verbose(level = 3,control = control)
    
    # Getting viterbi sequence
    viterbi <- computeViterbiSequence(hdf5 = hdf5File,pi = c(theta.new$pi),gamma = theta.new$gamma)
    
    # Removing temporary file
    system2('rm',hdf5File)
    
    return(viterbi)
}

################################################################################
### Function for differential peak calling
################################################################################

differentialHMM = function(object,control,dist){
    
    Group = Intercept = NULL
    
    # Creating subdirectories
    hdf5File <- checkPath(file.path(path.expand(control[['tempDir']]),paste0(control[['fileName']],'.h5')))
    invisible(lapply(hdf5File,function(x){if(!dir.exists(dirname(x))){dir.create(dirname(x),recursive = TRUE)}}))
    
    # General parameters
    M <- nrow(object)
    N <- ncol(object)
    K <- 3
    group <- as.numeric(factor(colData(object)$condition,levels = unique(colData(object)$condition)))
    nGroup <- length(unique(group))
    parHist <- list()
    controlHist <- list(controlList())
    theta.old <- list('pi' = NULL,'gamma' = NULL, 'delta' = NULL,'psi' = NULL)
    theta.new <- list('pi' = NULL,'gamma' = NULL, 'delta' = NULL,'psi' = NULL)
    
    # Initializing patterns & mixtures
    zPatterns <- enumeratePatterns(object = object,group = group)
    zMixtures <- determineMixtures(pattern = control[['pattern']],group = group,nGroup = nGroup,ref = zPatterns$ref)
    
    # Transforming data into data.table
    dt <- data.table::data.table(Window = rep(seq_len(M),N),ChIP = as.numeric(assay(object)),offsets = as.numeric(assay(object,'offsets')))
    dt[,paste0('Dsg.Mix',seq_len(zMixtures$B)) := generateMixtureIntercepts(pattern = control[['pattern']],group = group,nGroup = nGroup,ref = zPatterns$ref,M = M,B = zMixtures$B)]
    dt[,Group := .GRP,by=c('ChIP','offsets',paste0('Dsg.Mix',seq_len(zMixtures$B)))]
    
    # Creating unique data.table
    dtUnique <- unique(dt,by='Group')[,c('ChIP','offsets',paste0('Dsg.Mix',seq_len(zMixtures$B)),'Group'),with = FALSE]
    setkey(dtUnique,Group)
    
    # Parameter initializations
    theta.old[['pi']] <- c(0.999,0.0005,0.0005)
    theta.old[['gamma']] <- estimateTransitionProb(chain = 0*(zPatterns$group$Group==min(zPatterns$group$Group))+1*(!zPatterns$group$Group%in%range(zPatterns$group$Group))+2*(zPatterns$group$Group==max(zPatterns$group$Group)),
                                                   numStates = K)
    theta.old[['delta']] <- estimateMixtureProb(zDifferential = zPatterns$group$Group[!(zPatterns$group$Group%in%c(1,max(zPatterns$group$Group)))],
                                                zSeq = zMixtures$zSeq,
                                                B = zMixtures$B)
    theta.old[['psi']] <- estimateCoefficients(z = zPatterns$group$Group,dt = dt,dist = dist,type = 'differential',control = control)
    
    # EM algorithm begins
    ## Verbose
    init.time <- verbose(level = 1,control = control)
    
    while(checkConvergence(controlHist,control)<control[['maxCountEM']] & controlHist[[length(controlHist)]][['iteration']]<control[['maxIterEM']]){
        
        # Update iteration
        controlHist[[length(controlHist)]][['iteration']] = controlHist[[length(controlHist)]][['iteration']] + 1
        
        # E-step
        expStep(pi = theta.old[['pi']],
                gamma = theta.old[['gamma']],
                logf = loglikDifferentialHMM(dt = dt,theta = theta.old,N = N,M = M,minZero = control[['minZero']],dist = dist),
                hdf5 = hdf5File)
        
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(hdf5 = hdf5File)
        
        ## Model parameters
        ### Inner EM algorithm
        probCut <-  (control[['probCut']]>0)*max(c(0.9^controlHist[[length(controlHist)]][['iteration']],control[['probCut']]))
        
        count.innerem <- 0
        error.innerem <- 1
        theta.tmp.old <- theta.old
        theta.tmp.new <- theta.old
        
        while(count.innerem<control[['maxIterInnerEM']] & error.innerem>control[['epsilonInnerEM']]){
            #### Inner EM: E-step
            suppressMessages(innerExpStep(dt = dt,theta = theta.tmp.old,N = N,M = M,model = 'nb',hdf5 = hdf5File,minZero = control[['minZero']],dist = dist))
            
            #### Inner EM: M-step
            ##### Mixing probability
            theta.tmp.new[['delta']] <- c(innerMaxStepProb(hdf5 = hdf5File))
            
            ##### Model parameters
            rcWeights <- rbindlist(lapply(differentialRejectionControlled(hdf5 = hdf5File,f = dt$Group,p = probCut,N = N)[,1],
                                          function(x){colnames(x) <- c('weights','Group');return(cbind(x,dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1]))}),idcol = 'id')
            
            ###### Calculating MLEs
            optimMLE <- optimDifferential(par = unlist(theta.tmp.old$psi),rcWeights = rcWeights,control = control,dist = dist)
            theta.tmp.new[['psi']] <- unname(split(optimMLE$par, (seq_along(optimMLE$par)-1) %/% ifelse(dist == 'nb',2,3)))
            
            ##### Checking convergence
            error.innerem <- max(abs((c(theta.tmp.new[['delta']],unlist(theta.tmp.new[['psi']]))-c(theta.tmp.old[['delta']],unlist(theta.tmp.old[['psi']])))/c(theta.tmp.old[['delta']],unlist(theta.tmp.old[['psi']]))))
            
            ##### Updating paramaters
            theta.tmp.old <- theta.tmp.new
            count.innerem <- count.innerem + 1
        }
        ### Updating parameters
        theta.new[c('delta','psi')] <- theta.tmp.new[c('delta','psi')]
        
        # Updating parHist
        parHist[[controlHist[[length(controlHist)]][['iteration']]]] <- theta.new
        
        # Updating controlHist
        ## Q-function (evaluated on the old parameters to avoid calculating the log-likelihood again)
        controlHist[[length(controlHist)]][['time']] <- difftime(time1 = Sys.time(),time2 = init.time,units = 'mins') 
        controlHist[[length(controlHist)]][['bic']] <- computeBIC(hdf5 = hdf5File,
                                                                  numPar = length(unlist(theta.old)) - 4 - (length(control$pattern)-1),
                                                                  numSamples = ncol(object))
        controlHist[[length(controlHist)]][['q']] <- computeQFunction(hdf5 = hdf5File,pi = theta.old[['pi']],gamma = theta.old[['gamma']])
        controlHist[[length(controlHist)]][['convergence']] <- optimMLE$convergence
        controlHist[[length(controlHist)]][['error']] <- getError(controlHist = controlHist,parHist = parHist,control = control)
        controlHist[[length(controlHist)]][['count']] <- (controlHist[[length(controlHist)]][['error']]<=control[['epsilonEM']])*(controlHist[[length(controlHist)]][['iteration']]>control[['minIterEM']])*(controlHist[[length(controlHist)]][['count']]+1) + 0
        
        # Updating controlHist
        controlHist[[length(controlHist)+1]] <- controlHist[[length(controlHist)]]
        
        # Updating parameter history
        theta.old <- theta.new
        
        #Verbose
        verbose(level = 2,control = control,controlHist = controlHist,theta = theta.new)
    }
    
    # Verbose
    verbose(level = 3,control = control)
    
    # Saving mixture patterns
    rhdf5::h5write(obj = gsub('1','E',gsub('0','B',do.call(paste0,zPatterns$ref[unlist(zMixtures$zSeq),seq_len(nGroup),with = FALSE]))),
                   file = hdf5File,
                   name = 'mixturePatterns')
    
    # Saving viterbi sequence
    invisible(computeViterbiSequence(hdf5 = hdf5File,pi = theta.new$pi,gamma = theta.new$gamma))
    
    # Saving output path to file
    S4Vectors::metadata(object) <- list('output' = hdf5File,
                                        'control' = control,
                                        'history' = list('control' = data.table::rbindlist(lapply(controlHist,function(x){data.table::data.table(t(unlist(x)))})),
                                                         'parameter' = data.table::rbindlist(lapply(parHist,function(x){data.table::data.table(t(unlist(x)))}))))
    
    return(object)
}

################################################################################
### Function for consensus peak calling
################################################################################

consensusHMM = function(object,control,dist)
{
    
    controls = Group = Intercept = NULL
    
    # Creating subdirectories
    hdf5File <- checkPath(file.path(path.expand(control[['tempDir']]),paste0(control[['fileName']],'.h5')))
    invisible(lapply(hdf5File,function(x){if(!dir.exists(dirname(x))){dir.create(dirname(x),recursive = TRUE)}}))
    
    # General parameters
    M <- nrow(object)
    N <- ncol(object)
    K <- 2
    parHist <- list()
    controlHist <- list(controlList())
    theta.old <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    theta.new <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    
    # Initializing patterns & mixtures
    zPatterns <- enumeratePatterns(object = object,group = SummarizedExperiment::colData(object)$condition)
    
    # Transforming data into data.table
    dt <- data.table::data.table(Window = rep(seq_len(M),N),ChIP = as.numeric(assay(object)),offsets = as.numeric(assay(object,'offsets')))
    if('controls' %in% SummarizedExperiment::assayNames(object)){
        # Adding control & keying
        dt[,controls := as.numeric(assay(object,'controls'))]
        dt[,Group := .GRP,by = c('ChIP','controls','offsets')]
        
        # Creating unique data.table
        dtUnique <- unique(dt,by='Group')[,c('ChIP','controls','offsets','Group'),with = FALSE]
        setkey(dtUnique,Group)
    } else{
        #Keying
        dt[,Group := .GRP,by = c('ChIP','offsets')]
        
        # Creating unique data.table
        dtUnique <- unique(dt,by='Group')[,c('ChIP','offsets','Group'),with = FALSE]
        setkey(dtUnique,Group)
    }
    
    # Parameter initializations
    theta.old[['pi']] <- c(0.999,0.001)
    theta.old[['gamma']] <- estimateTransitionProb(chain = zPatterns$group$Group,numStates = K)
    theta.old[['psi']] <- estimateCoefficients(z = zPatterns$group$Group,dt = dt,dist = dist,type = 'consensus',control = control)
    
    # EM algorithm begins
    ## Verbose
    init.time <- verbose(level = 1,control = control)
    
    while(checkConvergence(controlHist,control)<control[['maxCountEM']] & controlHist[[length(controlHist)]][['iteration']]<control[['maxIterEM']]){
        
        # Update iteration
        controlHist[[length(controlHist)]][['iteration']] = controlHist[[length(controlHist)]][['iteration']] + 1
        
        # E-step
        expStep(pi = theta.old[['pi']],
                gamma = theta.old[['gamma']],
                logf = loglikConsensusHMM(dt = dt,theta = theta.old,N = N,M = M,minZero = control[['minZero']],dist = dist),
                hdf5 = hdf5File)
        
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(hdf5 = hdf5File)
        
        ## Model parameters
        ### Rejection-controlled posterior probabilities
        probCut <-  (control[['probCut']]>0)*max(c(0.9^controlHist[[length(controlHist)]][['iteration']],control[['probCut']]))
        rcWeights <- lapply(consensusRejectionControlled(hdf5 = hdf5File,f = dt$Group,p = probCut)[,1],
                            function(x){colnames(x) <- c('weights','Group');return(as.matrix(cbind(x,dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1])))})
        
        ### Calculating MLEs
        optimMLE <- lapply(seq_len(length(rcWeights)),function(i){
            optimNB(par = theta.old$psi[[i]],
                    y = rcWeights[[i]][,'ChIP'],
                    x = rcWeights[[i]][,if ('controls' %in% names(dt)) c('Intercept','controls') else 'Intercept',drop = FALSE],
                    offset = rcWeights[[i]][,'offsets'],
                    weights = rcWeights[[i]][,'weights'],
                    control = control,
                    dist = ifelse(dist=='zinb' & i==1,'zinb','nb'))
        })
        theta.new[['psi']] <- lapply(optimMLE,function(x){x[['par']]})
        
        # Updating parHist
        parHist[[controlHist[[length(controlHist)]][['iteration']]]] <- theta.new
        
        # Checking convergence
        ## Q-function (evaluated on the old parameters to avoid calculating the log-likelihood again)
        controlHist[[length(controlHist)]][['time']] <- difftime(time1 = Sys.time(),time2 = init.time,units = 'mins') 
        controlHist[[length(controlHist)]][['bic']] <- computeBIC(hdf5 = hdf5File,
                                                                  numPar = length(unlist(theta.old)) - 3,
                                                                  numSamples = ncol(object))
        controlHist[[length(controlHist)]][['q']] <- computeQFunction(hdf5 = hdf5File,pi = theta.old[['pi']],gamma = theta.old[['gamma']])
        controlHist[[length(controlHist)]][['convergence']] <-  1*any(!unname(unlist(lapply(optimMLE,function(x){x[['convergence']]}))) == 0)
        controlHist[[length(controlHist)]][['error']] <- getError(controlHist = controlHist,parHist = parHist,control = control)
        controlHist[[length(controlHist)]][['count']] <- (controlHist[[length(controlHist)]][['error']]<=control[['epsilonEM']])*(controlHist[[length(controlHist)]][['iteration']]>control[['minIterEM']])*(controlHist[[length(controlHist)]][['count']]+1) + 0
        
        # Updating controlHist
        controlHist[[length(controlHist)+1]] <- controlHist[[length(controlHist)]]
        
        # Updating parameter history
        theta.old <- theta.new
        
        #Verbose
        verbose(level = 2,control = control,controlHist = controlHist,theta = theta.new)
    }
    
    
    # Verbose
    verbose(level = 3,control = control)
    
    # Saving viterbi sequence
    invisible(computeViterbiSequence(hdf5 = hdf5File,pi = theta.new$pi,gamma = theta.new$gamma))
    
    # Saving output path to file
    S4Vectors::metadata(object) <- list('output' = hdf5File,
                                        'control' = control,
                                        'history' = list('control' = data.table::rbindlist(lapply(controlHist,function(x){data.table::data.table(t(unlist(x)))})),
                                                         'parameter' = data.table::rbindlist(lapply(parHist,function(x){data.table::data.table(t(unlist(x)))}))))
    
    return(object)
}