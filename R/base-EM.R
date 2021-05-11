################################################################################
################################################################################
### R script with the EM algorithm loops for initialization, consensus, and
### differential peak calling algorithms
################################################################################
################################################################################

################################################################################
### EM algorithm for epigraHMM initialization
################################################################################

emInitialization <- function(theta.old,theta.new,object,dt,dtUnique,hdf5File,parHist,controlHist,control,init.time){
    Group = Intercept = NULL
    itNumber <- length(controlHist)
    numPar <- length(unlist(theta.old)) - 3
    logicalWhile <-
        all(checkConvergence(controlHist, control) < control[['maxCountEM']],
            controlHist[[itNumber]][['iteration']] < control[['maxIterEM']])
    while (logicalWhile) {
        # Update iteration
        controlHist[[itNumber]][['iteration']] <- 
            controlHist[[itNumber]][['iteration']] + 1
        # E-step
        logf <- do.call(cbind,lapply(c(1,2),function(x){
            mu <- exp(theta.old[['psi']][[x]][1] + 
                          as.numeric(assay(object,'offsets')))
            stats::dnbinom(x = assay(object,'counts'),
                           mu = mu,
                           size = theta.old[['psi']][[x]][2],log = TRUE)
        }))
        expStep(pi = theta.old[['pi']],gamma = theta.old[['gamma']],logf = logf,
                hdf5 = hdf5File)
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(hdf5 = hdf5File)
        ## Model parameters
        ### Rejection-controlled posterior probabilities
        rejectionCut <- max(c(0.9^controlHist[[itNumber]][['iteration']],
                              control[['probCut']]))
        probCut <-  (control[['probCut']] > 0) * rejectionCut
        rejectionControl <- 
            consensusRejectionControlled(hdf5 = hdf5File,
                                         f = dt$Group,p = probCut)
        rcWeights <- lapply(rejectionControl[,1],function(x){
            colnames(x) <- c('weights','Group')
            out <- dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1]
            return(as.matrix(cbind(x,out)))
        })
        ### Calculating MLEs
        optimMLE <- lapply(seq_len(length(rcWeights)),function(x){
            optimNB(par = theta.old$psi[[x]],y = rcWeights[[x]][,'ChIP'],
                    x = rcWeights[[x]][,'Intercept',drop = FALSE],
                    offset = rcWeights[[x]][,'offsets'], 
                    weights = rcWeights[[x]][,'weights'],
                    control = control,dist = 'nb')
        })
        theta.new[['psi']] <- lapply(optimMLE,function(x){x[['par']]})
        
        # Updating parHist & controlHist
        parHist[[controlHist[[itNumber]][['iteration']]]] <- theta.new
        ctNum <- (controlHist[[itNumber]][['error']] <= control[['epsilonEM']])*
            (controlHist[[itNumber]][['iteration']] > control[['minIterEM']])*
            (controlHist[[itNumber]][['count']] + 1) + 0
        controlHist[[itNumber]][['count']] <- ctNum
        convNum <- 1*any(!unname(unlist(lapply(optimMLE,function(x){
            x[['convergence']]
        }))) == 0)
        controlHist[[itNumber]][['convergence']] <-  convNum
        controlHist[[itNumber]][['error']] <- 
            getError(controlHist = controlHist,parHist = parHist,
                     control = control)
        controlHist[[itNumber]][['time']] <- difftime(time1 = Sys.time(),
                                                      time2 = init.time,
                                                      units = 'mins') 
        controlHist[[itNumber]][['bic']] <- 
            computeBIC(hdf5 = hdf5File,numPar = numPar,
                       numSamples = ncol(object))
        controlHist[[itNumber]][['q']] <- 
            computeQFunction(hdf5 = hdf5File,pi = theta.old[['pi']],
                             gamma = theta.old[['gamma']])
        # Updating controlHist and parameters
        controlHist[[itNumber+1]] <- controlHist[[itNumber]]
        itNumber <- length(controlHist)
        logicalWhile <-
            all(checkConvergence(controlHist, control) < control[['maxCountEM']],
                controlHist[[itNumber]][['iteration']] < control[['maxIterEM']])
        theta.old <- theta.new
        verbose(level = 2,control = control,controlHist = controlHist,
                theta = theta.new)
    }
    return(theta.new)
}

################################################################################
### EM algorithm for epigraHMM consensus peak caller
################################################################################

emConsensus <- function(theta.old,theta.new,dt,dtUnique,M,N,dist,
                        init.time,hdf5File,controlHist,parHist,control){
    Group = Intercept = NULL
    itNumber <- length(controlHist)
    numPar <- length(unlist(theta.old)) - 3
    logicalWhile <-
        all(checkConvergence(controlHist, control) < control[['maxCountEM']],
            controlHist[[itNumber]][['iteration']] < control[['maxIterEM']])
    while (logicalWhile) {
        # Update iteration
        controlHist[[itNumber]][['iteration']] <- 
            controlHist[[itNumber]][['iteration']] + 1
        # E-step
        logf <- loglikConsensusHMM(dt = dt,theta = theta.old,N = N,M = M,
                                   minZero = control[['minZero']],dist = dist)
        expStep(pi = theta.old[['pi']],gamma = theta.old[['gamma']],
                hdf5 = hdf5File,logf = logf)
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(hdf5 = hdf5File)
        ## Model parameters
        rejectionCut <- max(c(0.9 ^ controlHist[[itNumber]][['iteration']],
                              control[['probCut']]))
        probCut <-  (control[['probCut']] > 0) * rejectionCut
        rejectionControl <- consensusRejectionControlled(hdf5 = hdf5File,
                                                         f = dt$Group,
                                                         p = probCut)
        rcWeights <- lapply(rejectionControl[,1],function(x){
            colnames(x) <- c('weights','Group')
            out <- dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1]
            return(as.matrix(cbind(x,out)))
        })
        ## Calculating MLEs
        optimMLE <- lapply(seq_len(length(rcWeights)),function(i){
            x <- rcWeights[[i]][,{
                if ('controls' %in% names(dt)) 
                    c('Intercept','controls') 
                else 'Intercept'
            },drop = FALSE]
            optimNB(par = theta.old$psi[[i]],y = rcWeights[[i]][,'ChIP'],
                    x = x,offset = rcWeights[[i]][,'offsets'],
                    weights = rcWeights[[i]][,'weights'],control = control,
                    dist = ifelse(dist == 'zinb' & i == 1, 'zinb', 'nb'))
        })
        theta.new[['psi']] <- lapply(optimMLE,function(x){x[['par']]})
        # Updating parHist & controlHist
        parHist[[controlHist[[itNumber]][['iteration']]]] <- theta.new
        ctNum <- (controlHist[[itNumber]][['error']] <= control[['epsilonEM']])*
            (controlHist[[itNumber]][['iteration']] > control[['minIterEM']])*
            (controlHist[[itNumber]][['count']] + 1) + 0
        controlHist[[itNumber]][['count']] <- ctNum
        
        controlHist[[itNumber]][['time']] <- difftime(time1 = Sys.time(),
                                                      time2 = init.time,
                                                      units = 'mins') 
        controlHist[[itNumber]][['bic']] <- computeBIC(hdf5 = hdf5File,
                                                       numPar = numPar,
                                                       numSamples = N)
        controlHist[[itNumber]][['q']] <- 
            computeQFunction(hdf5 = hdf5File,pi = theta.old[['pi']],
                             gamma = theta.old[['gamma']])
        convNum <- 1*any(!unname(unlist(lapply(optimMLE,function(x){
            x[['convergence']]
        }))) == 0)
        controlHist[[itNumber]][['convergence']] <- convNum
        controlHist[[itNumber]][['error']] <- getError(controlHist = controlHist,
                                                       parHist = parHist,
                                                       control = control)
        # Updating controlHist and parameters
        controlHist[[itNumber + 1]] <- controlHist[[itNumber]]
        itNumber <- length(controlHist)
        theta.old <- theta.new
        logicalWhile <-
            all(checkConvergence(controlHist, control) < control[['maxCountEM']],
                controlHist[[itNumber]][['iteration']] < control[['maxIterEM']])
        verbose(level = 2,control = control,controlHist = controlHist,
                theta = theta.new)
    }
    return(list('theta.new' = theta.new,'controlHist' = controlHist,
                'parHist' = parHist))
}

################################################################################
### Inner EM algorithm for epigraHMM differential peak caller
################################################################################

innerEMDifferential = function(theta.old,dt,dtUnique,probCut,N,M,hdf5File,
                               dist,control){
    Group = Intercept = NULL
    count.innerem <- 0; error.innerem <- 1
    theta.tmp.old <- theta.old; theta.tmp.new <- theta.old
    logicalWhile <- all(count.innerem < control[['maxIterInnerEM']],
                        error.innerem > control[['epsilonInnerEM']])
    while (logicalWhile) {
        # E-step
        suppressMessages(innerExpStep(dt = dt,theta = theta.tmp.old,N = N,M = M,
                                      model = 'nb',hdf5 = hdf5File,
                                      minZero = control[['minZero']],
                                      dist = dist))
        # M-step
        ## Mixing probabilities
        theta.tmp.new[['delta']] <- c(innerMaxStepProb(hdf5 = hdf5File))
        ## Model parameters
        rejectionControl <- differentialRejectionControlled(hdf5 = hdf5File,
                                                            f = dt$Group,
                                                            p = probCut,N = N)
        rcWeights <- 
            data.table::rbindlist(lapply(rejectionControl[,1],FUN = function(x){
                colnames(x) <- c('weights','Group')
                out <- dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1]
                return(cbind(x,out))
            }),idcol = 'id')
        optimMLE <- optimDifferential(par = unlist(theta.tmp.old$psi),
                                      rcWeights = rcWeights,
                                      control = control,dist = dist)
        par <- split(optimMLE$par,
                     (seq_along(optimMLE$par) - 1) %/% ifelse(dist == 'nb',2,3))
        theta.tmp.new[['psi']] <- unname(par)
        # Checking convergence
        parOld <- c(theta.tmp.new[['delta']],unlist(theta.tmp.new[['psi']]))
        parNew <- c(theta.tmp.old[['delta']],unlist(theta.tmp.old[['psi']]))
        error.innerem <- max(abs((parOld - parNew) / parNew))
        # Updating parameters
        theta.tmp.old <- theta.tmp.new
        count.innerem <- count.innerem + 1
        logicalWhile <- all(count.innerem < control[['maxIterInnerEM']],
                            error.innerem > control[['epsilonInnerEM']])
    }
    return(list('theta.new' = theta.tmp.new[c('delta','psi')],
                'optimMLE' = optimMLE))
}

################################################################################
### Outer EM algorithm for epigraHMM differential peak caller
################################################################################

outerEMDifferential = function(theta.old,theta.new,dt,dtUnique,N,M,dist,
                               init.time,hdf5File,controlHist,parHist,
                               control) {
    itNumber <- length(controlHist)
    numPar <- length(unlist(theta.old)) - length(control$pattern) - 3
    logicalWhile <- 
        all(checkConvergence(controlHist,control) < control[['maxCountEM']],
            controlHist[[itNumber]][['iteration']] < control[['maxIterEM']])
    while (logicalWhile) {
        # Update iteration
        controlHist[[itNumber]][['iteration']] <- 
            controlHist[[itNumber]][['iteration']] + 1
        # E-step
        logf <- loglikDifferentialHMM(dt = dt,theta = theta.old,N = N,M = M,
                                      minZero = control[['minZero']],
                                      dist = dist)
        expStep(pi = theta.old[['pi']],gamma = theta.old[['gamma']],
                hdf5 = hdf5File,logf = logf)
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(hdf5 = hdf5File)
        ## Model parameters
        rejectionCut <- max(c(0.9^controlHist[[itNumber]][['iteration']],
                              control[['probCut']]))
        probCut <-  (control[['probCut']] > 0) * rejectionCut
        innerEMDifferential <- 
            innerEMDifferential(theta.old,dt,dtUnique,
                                probCut,N,M,hdf5File,dist,control)
        theta.new[c('delta','psi')] <- innerEMDifferential[['theta.new']]
        # Updating parHist & controlHist
        parHist[[controlHist[[itNumber]][['iteration']]]] <- theta.new
        ctNum <-
            (controlHist[[itNumber]][['error']] <= control[['epsilonEM']]) *
            (controlHist[[itNumber]][['iteration']] > control[['minIterEM']]) *
            (controlHist[[itNumber]][['count']] + 1) + 0
        controlHist[[itNumber]][['count']] <- ctNum
        controlHist[[itNumber]][['time']] <- difftime(time1 = Sys.time(),
                                                      time2 = init.time,
                                                      units = 'mins') 
        controlHist[[itNumber]][['bic']] <- 
            computeBIC(hdf5 = hdf5File,numPar = numPar,numSamples = N)
        controlHist[[itNumber]][['q']] <- 
            computeQFunction(hdf5 = hdf5File,pi = theta.old[['pi']],
                             gamma = theta.old[['gamma']])
        controlHist[[itNumber]][['convergence']] <- 
            innerEMDifferential[['optimMLE']]$convergence
        controlHist[[itNumber]][['error']] <- 
            getError(controlHist = controlHist,parHist = parHist,
                     control = control)
        
        # Updating controlHist and parameters
        controlHist[[itNumber + 1]] <- controlHist[[itNumber]]
        itNumber <- length(controlHist)
        theta.old <- theta.new
        logicalWhile <- 
            all(checkConvergence(controlHist,control) < control[['maxCountEM']],
                controlHist[[itNumber]][['iteration']] < control[['maxIterEM']])
        verbose(level = 2,control = control,controlHist = controlHist,
                theta = theta.new)
    }
    return(list('theta.new' = theta.new,'parHist' = parHist,
                'controlHist' = controlHist))
}

################################################################################
### Inner E-step of the EM algorithm for differential calling
################################################################################

innerExpStep = function(dt,theta,N,M,model,hdf5,minZero,dist){
    ll <- loglikDifferential(dt = dt,delta = theta$delta,
                             psi = unlist(theta$psi),N = N,M = M,
                             dist = dist,minZero = minZero)
    ll <- rapply(ll,function(x){ifelse(exp(x) == 0,log(minZero),x)},
                 how = "replace")
    sumExpLL <- Reduce(`+`,lapply(ll,exp))
    rhdf5::h5write(obj = vapply(seq_len(length(theta$delta)),
                                FUN = function(b){exp(ll[[b]])/sumExpLL},
                                FUN.VALUE = vector('double',length = M)),
                   file = hdf5,
                   name = 'mixtureProb')
}