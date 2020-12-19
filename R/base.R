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
### Optimizer for glm.nb
################################################################################

optimNB <- function(par,y,x,offset,weights,control){
    
    x <- as.matrix(x)
    
    tryCatch({assign('model',optim(par=invertDispersion(par,model='nb'),
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
### Initial peak caller
################################################################################

initializerHMM = function(object,control){
    
    # Checking input
    if(!ncol(object)==1){
        stop('The initializer supports only one sample at a time')
    }
    
    # Creating subdirectories
    tempDir <- path.expand(control$tempDir)
    paths <- list(logFP = file.path(tempDir,'logFP'),
                  logBP = file.path(tempDir,'logBP'),
                  logP1 = file.path(tempDir,'logP1'),
                  logP2 = file.path(tempDir,'logP2'))
    lapply(unlist(paths),function(x){if(!dir.exists(x)){dir.create(x,recursive = TRUE)}})
    
    # General parameters
    M <- nrow(object)
    K <- 2
    errorEM <- c('error' = 1)
    iterEM <- 0
    countEM <- 0
    parList <- list()
    theta.old <- sapply(c('pi','gamma','psi'),function(x) NULL)
    theta.new <- sapply(c('pi','gamma','psi'),function(x) NULL)

    # Transforming data into data.table and calculating scores
    dt <- data.table::data.table(ChIP = as.numeric(assay(object)))[,Group := .GRP,by = 'ChIP']
    
    # Creating Unique data.table
    dtUnique <- unique(dt,by='Group')[,c('ChIP','Group')]
    setkey(dtUnique,Group)
    
    # Naive split
    score <- scale(log1p(dt$ChIP))
    z <- as.numeric(cut(score,breaks=c(-Inf,quantile(score,0.75),Inf)))
    
    # Parameter initializations
    theta.old[['pi']] <- c(0.999,1-0.999)
    theta.old[['gamma']] <- estimateTransitionProb(chain = z,numStates = 2)
    theta.old[['psi']] <- lapply(1:2,function(x){
        dt[which(z == x),.(mu = mean(log1p(ChIP)),sigma2 = var(log1p(ChIP)))][,c(mu,(mu^2)/(sigma2+mu))]
    })
    
    # EM algorithm begins
    message(paste0(c(rep('#',80))));message(Sys.time());message("Starting the EM algorithm")
    
    while(countEM<control[['maxCountEM']] & iterEM<control[['maxIterEM']]){
        iterEM = iterEM + 1
        
        # E-step
        expStep(counts = assay(object,'counts'),
                pi = theta.old[['pi']],
                gamma = theta.old[['gamma']],
                logf = do.call(cbind,lapply(1:2,function(x){dnbinom(x = assay(object,'counts'),mu = exp(theta.old[['psi']][[x]][1]),size = theta.old[['psi']][[x]][2],log = TRUE)})),
                nameForwardProb = file.path(paths[['logFP']],paste0('logFP.bin')),
                nameBackwardProb = file.path(paths[['logBP']],paste0('logBP.bin')),
                nameMarginalProb = file.path(paths[['logP1']],paste0('logP1.bin')),
                nameJointProb = file.path(paths[['logP2']],paste0('logP2.bin')))
        
        # M-step
        ## Initial and transition probabilities
        theta.new[c('pi','gamma')] <- maxStepProb(nameMarginalProb = file.path(paths[['logP1']],paste0('logP1.bin')),
                                                  nameJointProb = file.path(paths[['logP2']],paste0('logP2.bin')))
        
        ## Model parameters
        ### Rejection-controlled posterior probabilities
        probCut <-  (control[['probCut']]>0)*ifelse((0.9^iterEM)>=control[['probCut']],(0.9^iterEM),control[['probCut']])
        rcWeights <- lapply(rejectionControlled(nameMarginalProb = file.path(paths[['logP1']],paste0('logP1.bin')),f = dt$Group,p = probCut)[,1],
                            function(x){colnames(x) <- c('weights','Group');return(as.matrix(cbind(x,dtUnique[match(x[,2],Group),][,-c('Group')][,Intercept := 1])))})
        
        ### Calculating MLEs
        optimMLE <- lapply(seq_len(length(rcWeights)),function(x){
            optimNB(par = theta.old$psi[[x]],y = rcWeights[[x]][,'ChIP'],
                    x = rcWeights[[x]][,'Intercept'],offset = 0,
                    weights = rcWeights[[x]][,'weights'],control = control)
        })
        theta.new[['psi']] <- lapply(optimMLE,function(x){x[['par']]})
        
        # Updating parameter history
        errorEM['error'] <- sum(unlist(lapply(names(theta.old),function(x){sqrt(sum((unlist(theta.old[[x]])-unlist(theta.new[[x]]))^2))})))
        
        parList[[iterEM]] <- c(it=iterEM,errorEM,m=1*any(!unname(unlist(lapply(optimMLE,function(x){x[['convergence']]}))) == 0))
        
        theta.old <- theta.new
        
        # Computing EM error
        countEM <- as.numeric(errorEM<=control[['epsilonEM']][1])*(iterEM>control[['minIterEM']])*(countEM+1) + 0
        
        #Outputing history
        if(!control[['quiet']]){
            message(paste0(c(rep('#',80))))
            message('\rIteration: ',iterEM,'. Error: ',paste(formatC(errorEM, format = "e", digits = 2)),sep='')
            message("\r",paste('Initial prob. estimates: '),paste(formatC(theta.new[['pi']], format = "e", digits = 2),collapse = ' '))
            message("\r",paste('Transition prob. estimates: '),paste(formatC(theta.new[['gamma']], format = "e", digits = 2),collapse = ' '))
            message("\r",paste('Model paramater estimates: '),paste(formatC(unlist(theta.new[['psi']]), format = "e", digits = 2),collapse = ' '))
            message(paste0(c(rep('#',80))))
        }
    }
    
    getViterbiSequence(nameForwardProb = file.path(paths[['logFP']],paste0('logFP.bin')),pi = c(theta.new$pi),gamma = theta.new$gamma)
    
    # Organizing output
    z = hmm2_Viterbi(LOGF=loglik,P=pi.k1,GAMMA=gamma.k1)
    logF <- setnames(as.data.table(logF),c('Background','Enrichment'))
    logB <- setnames(as.data.table(logB),c('Background','Enrichment'))
    loglik <- setnames(as.data.table(loglik),c('Background','Enrichment'))
    mu <- as.data.table(mu)
    return(list('Pi'=pi.k1,'Gamma'=gamma.k1,'Psi'=psi.k1,'Prob'=dt[,list(PostProb1,PostProb2)],
                'LogF'=logF,'LogB'=logB,'Loglik'=loglik,'Parhist'=as.data.table(do.call(rbind,parlist)),'Mean'=mu,'Viterbi'=z))
}