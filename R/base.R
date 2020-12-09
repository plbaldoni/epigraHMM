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
### Initial peak caller
################################################################################

initializerHMM = function(object,control){
    
    # Checking input
    if(!ncol(object)==1){
        stop('The initializer supports only one sample at a time')
    }
    
    # Creating subdirectories
    tempDir <- path.expand(control$tempDir)
    paths <- list(logF = file.path(tempDir,'logF'),
                  logB = file.path(tempDir,'logB'),
                  logP = file.path(tempDir,'logP'))
    lapply(unlist(paths),function(x){if(!dir.exists(x)){system2('mkdir',paste('-p',x))}})
    
    # General parameters
    M <- nrow(object)
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
    
    while(count.em<maxcount.em & it.em<maxit.em){
        it.em = it.em + 1
        
        # Updating parameters
        # pi.k = theta.k[paste0('pi',seq_len(K))]
        # gamma.k = matrix(theta.k[paste0('gamma',as.character(transform(expand.grid(seq_len(K),seq_len(K)),idx=paste0(Var1,Var2))$idx))],nrow=K,ncol=K,byrow=FALSE);k=(K+1);for(i in seq_len(K)){for(j in seq_len(K)){assign(paste0('gamma',j,i,'.k'),theta.k[k]);k=k+1}}
        # psi1.k = theta.k[paste0('HMM1.',c(namesControl,'Disp'))]
        # psi2.k = theta.k[paste0('HMM2.',c(namesControl,'Disp'))]
        
        # E-step
        # mu <- HMM.mean.consensus(X.mat=as.matrix(dt[,grepl('Dsg',names(dt)),with=FALSE]),offset.vec=dt[,offset],psi=as.matrix(rbind(psi1.k,psi2.k)[,seq_len(ncolControl)]),N=N,M=M)
        # loglik <- HMM.LL(Y.vec=dt[,ChIP],mu=mu,disp=rbind(psi1.k,psi2.k)[,(ncolControl+1)],N=N,M=M,K=K,model='nb')
        # 
        # Forward-Backward probabilities
        logForwardBackward(counts = assay(object,'counts'),
                        pi = theta.old[['pi']],
                        gamma = theta.old[['gamma']],
                        logf = do.call(cbind,lapply(1:2,function(x){dnbinom(x = assay(object,'counts'),mu = exp(theta.old[['psi']][[x]][1]),size = theta.old[['psi']][[x]][2],log = TRUE)})),
                        nameF = file.path(paths[['logF']],paste0('logF.bin')),
                        nameB = file.path(paths[['logB']],paste0('logB.bin')))
        
        # Posterior probabilities
        dt[,paste0('PostProb',seq_len(2)):=as.data.table(check.prob(hmm2_P1(logF=logF,logB=logB)))]
        dt[,paste0('JoinProb',c('11','12','21','22')):=as.data.table(check.prob(hmm2_P2(logF=logF,logB=logB,logf1=loglik[,1],logf2=loglik[,2],gamma=gamma.k)))]
        
        # M-step
        ## Initial and transition probabilities
        PostProb = HMM.prob.consensus(dt = dt)
        pi.k1 = PostProb$pi
        gamma.k1 = PostProb$gamma
        
        ## Model parameters
        ### Aggregating data
        rejection = (pcut>0)*ifelse((0.9^it.em)>=pcut,(0.9^it.em),pcut)
        dt1 <- agg(dt[,Rejection1 := PostProb1][PostProb1<rejection,Rejection1 := rbinom(.N,1,prob=PostProb1/rejection)*rejection],data.unique = dt.unique,rows = '(Rejection1>0)',agg = 'Rejection1')
        dt2 <- agg(dt[,Rejection2 := PostProb2][PostProb2<rejection,Rejection2 := rbinom(.N,1,prob=PostProb2/rejection)*rejection],data.unique = dt.unique,rows = '(Rejection2>0)',agg = 'Rejection2')
        
        ### Calculating MLEs
        tryCatch({assign('model1',optim(par=inv.par(psi1.k,model='nb'),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                        Y.vec=dt1[,ChIP],X.mat=as.matrix(dt1[,grepl('Dsg',names(dt1)),with=FALSE]),offset.vec=dt1[,offset],weights.vec=dt1[,weights]))},
                 error=function(e){assign('model1',list('par' = inv.par(psi1.k,model='nb'),'convergence' = 99),inherits = TRUE)})
        tryCatch({assign('model2',optim(par=inv.par(psi2.k,model='nb'),fn=glm.nb,gr=deriv.nb,method='L-BFGS-B',lower=c(rep(-Inf,ncolControl),1/max.phi),
                                        Y.vec=dt2[,ChIP],X.mat=as.matrix(dt2[,grepl('Dsg',names(dt2)),with=FALSE]),offset.vec=dt2[,offset],weights.vec=dt2[,weights]))},
                 error=function(e){assign('model2',list('par' = inv.par(psi2.k,model='nb'),'convergence' = 99),inherits = TRUE)})
        
        rm(dt1);rm(dt2)
        
        ### Saving parameters
        psi1.k1 = inv.par(model1$par,model='nb')
        psi2.k1 = inv.par(model2$par,model='nb')
        psi.k1 = c(psi1.k1,psi2.k1)
        
        # Updating parameter history
        theta.k1 = c(pi.k1,gamma.k1,psi.k1)
        names(theta.k1) = names(theta.k)
        theta.k = theta.k1
        parlist[[it.em]] = c(it=it.em,error=error.em,theta.k1,m1=model1$convergence,m2=model2$convergence)
        
        # Computing EM error
        gap = ifelse(it.em>minit.em,gap.em,1)
        if(it.em>1){
            parlist.old = parlist[[(it.em-gap)]][names(psi.k1)]
            parlist.new = parlist[[it.em]][names(psi.k1)]
        } else{
            parlist.old = rep(1,length(names(psi.k1)))
            parlist.new = rep(1,length(names(psi.k1)))
        }
        MRCPE = max(abs((parlist.new-parlist.old)/parlist.old)) #Max. Abs. Rel. Change. of par. estimates
        error.em = ifelse(it.em>=2,MRCPE,1)
        count.em = as.numeric(any(error.em<=epsilon.em))*(it.em>minit.em)*(count.em+1) + 0
    }
    
    # Organizing output
    z = hmm2_Viterbi(LOGF=loglik,P=pi.k1,GAMMA=gamma.k1)
    logF <- setnames(as.data.table(logF),c('Background','Enrichment'))
    logB <- setnames(as.data.table(logB),c('Background','Enrichment'))
    loglik <- setnames(as.data.table(loglik),c('Background','Enrichment'))
    mu <- as.data.table(mu)
    return(list('Pi'=pi.k1,'Gamma'=gamma.k1,'Psi'=psi.k1,'Prob'=dt[,list(PostProb1,PostProb2)],
                'LogF'=logF,'LogB'=logB,'Loglik'=loglik,'Parhist'=as.data.table(do.call(rbind,parlist)),'Mean'=mu,'Viterbi'=z))
}