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
### Invert parameters for glm.nb and glm.zinb
################################################################################

invertDispersion = function(par,model){
    if(model=='nb'){return(c(par[seq_len((length(par)-1))],1/par[length(par)]))}
    if(model=='zinb'){
        n = length(par)
        return(c(par[seq_len(((n-1)/2))],1/par[(n-1)/2+1],par[seq(from = ((n-1)/2+2),to = n)]))}
}

################################################################################
### Random generator for directory path name
################################################################################

generatePath = function(){
    return(paste(sample(c(letters,LETTERS,0:9),10),collapse = ''))
}

################################################################################
### Print message during EM algorithm
################################################################################

verbose = function(level,control,controlHist = NULL,theta = NULL){
    if(!control[['quiet']]){
        if(level == 1){
            message(paste0(c(rep('#',80))));message(Sys.time());message("Starting the EM algorithm")
        }
        if(level == 2){
            message(paste0(c(rep('#',80))))
            message('\rIteration: ',controlHist[[length(controlHist)]][['iteration']])
            message("\rError: ",paste(paste(names(controlHist[[length(controlHist)]][['error']]),'='),formatC(controlHist[[length(controlHist)]][['error']], format = "e", digits = 2),collapse = ', '))
            message("\rConvergence: ",paste(paste(names(controlHist[[length(controlHist)]][['count']]),'='),controlHist[[length(controlHist)]][['count']],collapse = ', '))
            message("\r",paste('Initial probabilities: '),paste(formatC(theta[['pi']], format = "e", digits = 2),collapse = ' '))
            message("\r",paste('Transition probabilities: '),paste(formatC(theta[['gamma']], format = "e", digits = 2),collapse = ' '))
            if('delta' %in% names(theta)){
                message("\r",paste('Mixture probabilities: '),paste(formatC(unlist(theta[['delta']]), format = "e", digits = 2),collapse = ' '))
            }
            message("\r",paste('Model coefficients: '),paste(formatC(unlist(theta[['psi']]), format = "e", digits = 2),collapse = ' '))
        }
        if(level == 3){
            message(paste0(c(rep('#',80))))
            message("EM algorithm converged!")
            message(Sys.time())
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

optimNB <- function(par,y,x,offset,weights,control){
    
    tryCatch({assign('model',stats::optim(par=invertDispersion(par,model='nb'),
                                          fn=glmNB,
                                          gr=derivNB,
                                          method='L-BFGS-B',
                                          lower=c(rep(-Inf,ncol(x)),1/control$maxDisp),
                                          y = y,
                                          x = x,
                                          offset = offset,
                                          weights = weights))},
             error=function(e){assign('model',list('par' = invertDispersion(par,model='nb'),'convergence' = 99),inherits = TRUE)})
    
    model$par <- invertDispersion(model$par,model='nb')
    
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

estimateCoefficients <- function(z,dt){
    Window = ChIP = offsets = mu = sigma2 = NULL
    
    par <- lapply(range(z),function(x){
        dt[Window %in% which(z == x),list(mu = mean(ChIP/exp(offsets)),sigma2 = stats::var(ChIP/exp(offsets)))][,c(mu,(mu^2)/(sigma2-mu))]
    })
    par[[1]] <- log(par[[1]])
    par[[2]] <- log(par[[2]]) - par[[1]]
    return(par)
}

################################################################################
### Function for log-likelihood function of mixture of NB
### components sharing unique set of parameters
################################################################################

loglikMixNB = function(dt,delta,psi,N,M){
    Dsg.Mix = ChIP = offsets = NULL
    
    ll <- lapply(seq_len(length(delta)),FUN = function(b){
        log(delta[b])+rowSums(matrix(dt[,Dsg.Mix := .SD,.SDcols = paste0('Dsg.Mix',b)][,stats::dnbinom(x = ChIP,mu = exp(psi[1]+psi[3]*Dsg.Mix+offsets),size = exp(psi[2]+psi[4]*Dsg.Mix),log=TRUE)],nrow=M,ncol=N,byrow=FALSE))
    })
    
    dt[,Dsg.Mix := NULL]
    return(ll)
}

################################################################################
### Function for log-likelihood function of HMM with NB distributions with mixture
### components sharing unique set of parameters
################################################################################

loglikMixNBHMM = function(dt,theta,N,M,minZero){
    
    ChIP = offsets = Dsg.Mix = NULL
    
    psi <- unlist(theta$psi)
    delta <- theta$delta
    B <- length(delta)
    LL <- matrix(0,nrow=M,ncol=3)
    
    #LL from Background
    LL[,1] <- rowSums(matrix(dt[,stats::dnbinom(x = ChIP,mu = exp(psi[1]+offsets),size = exp(psi[2]),log=TRUE)],nrow=M,ncol=N,byrow=FALSE))
    #LL from Mixture
    LL[,2] <- log(base::Reduce(`+`,lapply(loglikMixNB(dt = dt,delta = delta,psi = psi,N = N,M = M),exp)))
    #LL from Enrichment
    LL[,3] <- rowSums(matrix(dt[,stats::dnbinom(x = ChIP,mu = exp((psi[1]+psi[3])+offsets),size = exp((psi[2]+psi[4])),log=TRUE)],nrow=M,ncol=N,byrow=FALSE))
    
    LL[is.infinite(LL)] <- log(minZero)
    return(LL)
}

################################################################################
### Function that implements the inner E-step of the EM algorithm for differential
### calling
################################################################################

innerExpStep = function(dt,theta,N,M,model,hdf5,minZero){
    
    if(model == 'nb'){
        ll <- loglikMixNB(dt = dt,delta = theta$delta,psi = unlist(theta$psi),N = N,M = M)
    }
    if(model == 'zinb'){
        message('Needs implementation')
    }
    ll <- rapply(ll,function(x){ifelse(exp(x)==0,log(minZero),x)}, how = "replace")
    sumExpLL <- Reduce(`+`,lapply(ll,exp))
    
    saveMixtureProb(eta = vapply(seq_len(length(theta$delta)),FUN = function(b){exp(ll[[b]])/sumExpLL},FUN.VALUE = vector('double',length = M)),
                    hdf5 = hdf5)
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
### Optimizer for glmNB
################################################################################

optimMixNB <- function(par,rcWeights,control){
    
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
    
    return(model)
}

################################################################################
### Create control list to track EM algorithm
################################################################################

controlList <- function(iteration = 0,
                        convergence = 0,
                        q = 0,
                        MRCPE = c('error' = 0,'count' = 0),
                        MACPE = c('error' = 0,'count' = 0),
                        ARCEL = c('error' = 0,'count' = 0)){
    return(list('iteration' = iteration,
                'convergence' = convergence,
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
    splitPath <- strsplit(path,"\\.h5")[[1]]
    repeat {
        f = paste0(splitPath[1],i,".h5")
        if(!file.exists(f)){return(f)}
        i <- i + 1
    }
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
    hdf5File <- checkPath(file.path(path.expand(control$tempDir),generatePath(),'output.h5'))
    invisible(lapply(hdf5File,function(x){if(!dir.exists(dirname(x))){dir.create(dirname(x),recursive = TRUE)}}))
    
    # General parameters
    M <- nrow(object)
    K <- 2
    controlHist <- list(controlList())
    parHist <- list()
    theta.old <- sapply(c('pi','gamma','psi'),function(x) NULL)
    theta.new <- sapply(c('pi','gamma','psi'),function(x) NULL)
    
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
    theta.old[['psi']] <- lapply(1:2,function(x){
        muvar <- dt[which(z == x),list(mu = mean(ChIP/exp(offsets)),sigma2 = stats::var(ChIP/exp(offsets)))]
        muvar[,c(log(mu),min((mu^2)/max(0,sigma2-mu),control[['maxDisp']]))]
    })
    rm(z,score)
    
    # EM algorithm begins
    ## Verbose
    verbose(level = 1,control = control)
    
    while(max(controlHist[[length(controlHist)]][['count']])<control[['maxCountEM']] & controlHist[[length(controlHist)]][['iteration']]<control[['maxIterEM']]){
        
        # Update iteration
        controlHist[[length(controlHist)]][['iteration']] = controlHist[[length(controlHist)]][['iteration']] + 1
        
        # E-step
        expStep(pi = theta.old[['pi']],
                gamma = theta.old[['gamma']],
                logf = do.call(cbind,lapply(1:2,function(x){stats::dnbinom(x = assay(object,'counts'),mu = exp(theta.old[['psi']][[x]][1]+assay(object,'offsets')),size = theta.old[['psi']][[x]][2],log = TRUE)})),
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
            optimNB(par = theta.old$psi[[x]],y = rcWeights[[x]][,'ChIP'],
                    x = rcWeights[[x]][,'Intercept',drop = FALSE],offset = rcWeights[[x]][,'offsets'],
                    weights = rcWeights[[x]][,'weights'],control = control)
        })
        theta.new[['psi']] <- lapply(optimMLE,function(x){x[['par']]})
        
        # Updating parHist
        parHist[[controlHist[[length(controlHist)]][['iteration']]]] <- theta.new
        
        # Checking convergence
        ## Q-function (evaluated on the old parameters to avoid calculating the log-likelihood again)
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
### mixNBHMM function
################################################################################

mixNBHMM = function(object,control){
    
    Group = Intercept = NULL
    
    # Creating subdirectories
    hdf5File <- checkPath(file.path(path.expand(control$tempDir),generatePath(),'output.h5'))
    invisible(lapply(hdf5File,function(x){if(!dir.exists(dirname(x))){dir.create(dirname(x),recursive = TRUE)}}))
    
    # General parameters
    M <- nrow(object)
    N <- ncol(object)
    K <- 3
    group <- as.numeric(factor(colData(object)$condition,levels = unique(colData(object)$condition)))
    nGroup <- length(unique(group))
    parHist <- list()
    controlHist <- list(controlList())
    theta.old <- sapply(c('pi','gamma','delta','psi'),function(x) NULL)
    theta.new <- sapply(c('pi','gamma','delta','psi'),function(x) NULL)
    
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
    theta.old[['psi']] <- estimateCoefficients(z = zPatterns$group$Group,dt = dt)
    
    # EM algorithm begins
    ## Verbose
    verbose(level = 1,control = control)
    
    while(max(controlHist[[length(controlHist)]][['count']])<control[['maxCountEM']] & controlHist[[length(controlHist)]][['iteration']]<control[['maxIterEM']]){
        
        # Update iteration
        controlHist[[length(controlHist)]][['iteration']] = controlHist[[length(controlHist)]][['iteration']] + 1
        
        # E-step
        expStep(pi = theta.old[['pi']],
                gamma = theta.old[['gamma']],
                logf = loglikMixNBHMM(dt = dt,theta = theta.old,N = N,M = M,minZero = control[['minZero']]),
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
            innerExpStep(dt = dt,theta = theta.tmp.old,N = N,M = M,model = 'nb',hdf5 = hdf5File,minZero = control[['minZero']])
            
            #### Inner EM: M-step
            ##### Mixing probability
            theta.tmp.new[['delta']] <- c(innerMaxStepProb(hdf5 = hdf5File))
            
            ##### Model parameters
            rcWeights <- rbindlist(lapply(differentialRejectionControlled(hdf5 = hdf5File,f = dt$Group,p = probCut,N = N)[,1],
                                          function(x){colnames(x) <- c('weights','Group');return(cbind(x,dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1]))}),idcol = 'id')
            
            ###### Calculating MLEs
            optimMLE <- optimMixNB(par = unlist(theta.tmp.old$psi),rcWeights = rcWeights,control = control)
            theta.tmp.new[['psi']] <- unname(split(optimMLE$par, (seq_along(optimMLE$par)-1) %/% 2))
            
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
    
    # Saving viterbi sequence
    invisible(computeViterbiSequence(hdf5 = hdf5File,pi = theta.new$pi,gamma = theta.new$gamma))
    
    # Saving output path to file
    S4Vectors::metadata(object) <- list('output' = hdf5File,
                                        'control' = control,
                                        'history' = list('control' = data.table::rbindlist(lapply(controlHist,function(x){data.table::data.table(t(unlist(x)))})),
                                                         'parameter' = data.table::rbindlist(lapply(parHist,function(x){data.table::data.table(t(unlist(x)))}))))
    
    return(object)
}