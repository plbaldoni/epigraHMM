################################################################################
################################################################################
### R script w/ the core consensus and differential peak callers
################################################################################
################################################################################

################################################################################
### Initial peak caller
################################################################################

initializerHMM = function(object,control){
    Group = ChIP = mu = sigma2 = offsets = .GRP = NULL
    # Checking input
    if (!ncol(object) == 1)
        stop('The initializer supports only one sample at a time')
    # Creating subdirectories
    hdf5File <- checkPath(file.path(path.expand(control[['tempDir']]),
                                    paste0(control[['fileName']],'.h5')))
    invisible(lapply(hdf5File,function(x){
        if (!dir.exists(dirname(x))) dir.create(dirname(x),recursive = TRUE)
    }))
    # General parameters
    K <- 2
    controlHist <- list(controlList()); parHist <- list()
    theta.old <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    theta.new <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    # Transforming data into data.table and calculating scores
    dt <- data.table::data.table(ChIP = as.numeric(assay(object)),
                                 offsets = as.numeric(assay(object,'offsets')))
    dt[,Group := .GRP,by = c('ChIP','offsets')]
    # Creating Unique data.table
    dtUnique <-
        unique(dt, by = 'Group')[, c('ChIP', 'offsets', 'Group'), with = FALSE]
    data.table::setkey(dtUnique,Group)
    # Naive split
    score <- dt[,scale(log1p(ChIP/exp(offsets)))]
    quantBreaks <- c(-Inf,stats::quantile(score,0.75),Inf)
    z <- as.numeric(cut(score, breaks = quantBreaks))
    # Parameter initializations
    theta.old[['pi']] <- c(0.999,0.001)
    theta.old[['gamma']] <- estimateTransitionProb(chain = z,numStates = K)
    theta.old[['psi']] <- lapply(c(1,2),function(x){
        muvar <- dt[which(z == x),{
            list(mu = mean((ChIP + 1) / exp(offsets)),
                 sigma2 = stats::var((ChIP + 1) / exp(offsets)))
        }]
        muvar[, c(log(mu),
                  min((mu ^ 2) / max(0, sigma2 - mu), control[['maxDisp']]))]
    })
    rm(z,score)
    # EM algorithm begins
    init.time <- verbose(level = 1,control = control)
    theta.new <- emInitialization(theta.old,theta.new,object,dt,dtUnique,
                                  hdf5File,parHist,controlHist,control,
                                  init.time)
    verbose(level = 3,control = control)
    # Getting viterbi sequence
    viterbi <- computeViterbiSequence(hdf5 = hdf5File,pi = c(theta.new$pi),
                                      gamma = theta.new$gamma)
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
    hdf5File <- checkPath(file.path(path.expand(control[['tempDir']]),
                                    paste0(control[['fileName']],'.h5')))
    invisible(lapply(hdf5File,function(x){
        if (!dir.exists(dirname(x))) dir.create(dirname(x),recursive = TRUE)
    }))
    # General parameters
    M <- nrow(object);N <- ncol(object);K <- 3
    group <- as.numeric(factor(colData(object)$condition,
                               levels = unique(colData(object)$condition)))
    nGroup <- length(unique(group)); parHist <- list(); controlHist <- list(controlList())
    theta.old <- list('pi' = NULL,'gamma' = NULL, 'delta' = NULL,'psi' = NULL)
    theta.new <- list('pi' = NULL,'gamma' = NULL, 'delta' = NULL,'psi' = NULL)
    # Initializing patterns & mixtures
    zPatterns <- enumeratePatterns(object = object,group = group)
    zMixtures <- determineMixtures(pattern = control[['pattern']],group = group,
                                   nGroup = nGroup,ref = zPatterns$ref)
    # Transforming data into data.table
    dt <- data.table::data.table(Window = rep(seq_len(M),N),
                                 ChIP = as.numeric(assay(object)),
                                 offsets = as.numeric(assay(object,'offsets')))
    dt[,paste0('Dsg.Mix',seq_len(zMixtures$B)) := {
        generateMixtureIntercepts(pattern = control[['pattern']],group = group,
                                  nGroup = nGroup,ref = zPatterns$ref,
                                  M = M,B = zMixtures$B)
    }]
    dt[, Group := .GRP, by = c('ChIP', 'offsets',
                               paste0('Dsg.Mix', seq_len(zMixtures$B)))]
    # Creating unique data.table
    dtUnique <- unique(dt, by = 'Group')
    dtUnique <- dtUnique[, c('ChIP', 'offsets', 
                             paste0('Dsg.Mix', seq_len(zMixtures$B)), 'Group'),
                         with = FALSE]
    data.table::setkey(dtUnique,Group)
    # Parameter initialization
    theta.old[['pi']] <- c(0.999,0.0005,0.0005)
    initialChain <- 0 * (zPatterns$group$Group == min(zPatterns$group$Group)) +
        1 * (!zPatterns$group$Group %in% range(zPatterns$group$Group)) +
        2 * (zPatterns$group$Group == max(zPatterns$group$Group))
    theta.old[['gamma']] <- estimateTransitionProb(chain = initialChain,
                               numStates = K)
    diffPatterns <- zPatterns$group$Group[!(zPatterns$group$Group %in% 
                                                c(1,max(zPatterns$group$Group)))]
    theta.old[['delta']] <- 
        estimateMixtureProb(zDifferential = diffPatterns,
                            zSeq = zMixtures$zSeq,B = zMixtures$B)
    theta.old[['psi']] <- 
        estimateCoefficients(z = zPatterns$group$Group,dt = dt,dist = dist,
                             type = 'differential',control = control)
    # EM algorithm begins
    init.time <- verbose(level = 1,control = control)
    outerEMDifferential <- 
        outerEMDifferential(theta.old,theta.new,dt,dtUnique,N,M,dist,init.time,
                            hdf5File,controlHist,parHist,control)
    verbose(level = 3,control = control)
    # Saving output and returning object
    mixturePatterns <- 
        do.call(paste0,zPatterns$ref[unlist(zMixtures$zSeq),seq_len(nGroup),
                                     with = FALSE])
    rhdf5::h5write(obj = gsub('1','E',gsub('0','B',mixturePatterns)),
                   file = hdf5File,name = 'mixturePatterns')
    invisible(computeViterbiSequence(hdf5 = hdf5File,
                                     pi = outerEMDifferential[['theta.new']]$pi,
                                     gamma = outerEMDifferential[['theta.new']]$gamma))
    controlHistory <- 
        data.table::rbindlist(lapply(outerEMDifferential[['controlHist']],function(x){
            data.table::data.table(t(unlist(x)))
        }))
    parHistory <- 
        data.table::rbindlist(lapply(outerEMDifferential[['parHist']],function(x){
            data.table::data.table(t(unlist(x)))
        }))
    S4Vectors::metadata(object) <- list('output' = hdf5File,'control' = control,
                                        'history' = list('control' = controlHistory,
                                                         'parameter' = parHistory))
    return(object)
}

################################################################################
### Function for consensus peak calling
################################################################################

consensusHMM = function(object,control,dist){
    controls = Group = Intercept = NULL
    # Creating subdirectories
    hdf5File <- checkPath(file.path(path.expand(control[['tempDir']]),
                                    paste0(control[['fileName']],'.h5')))
    invisible(lapply(hdf5File,function(x){
        if (!dir.exists(dirname(x))) dir.create(dirname(x),recursive = TRUE)
    }))
    # General parameters
    M <- nrow(object);N <- ncol(object);K <- 2
    parHist <- list(); controlHist <- list(controlList())
    theta.old <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    theta.new <- list('pi' = NULL,'gamma' = NULL,'psi' = NULL)
    # Initializing patterns & mixtures
    zPatterns <- 
        enumeratePatterns(object = object,
                          group = SummarizedExperiment::colData(object)$condition)
    # Transforming data into data.table
    dt <- data.table::data.table(Window = rep(seq_len(M),N),
                                 ChIP = as.numeric(assay(object)),
                                 offsets = as.numeric(assay(object,'offsets')))
    if ('controls' %in% SummarizedExperiment::assayNames(object)) {
        # Adding control & keying
        dt[,controls := as.numeric(log1p(assay(object,'controls')))]
        dt[,Group := .GRP,by = c('ChIP','controls','offsets')]
        # Creating unique data.table
        dtUnique <- 
            unique(dt, by = 'Group')[, c('ChIP', 'controls', 'offsets', 'Group'),with = FALSE]
    } else{
        # Keying
        dt[,Group := .GRP,by = c('ChIP','offsets')]
        # Creating unique data.table
        dtUnique <- 
            unique(dt, by = 'Group')[, c('ChIP', 'offsets', 'Group'),with = FALSE]
    }
    data.table::setkey(dtUnique,Group)
    # Parameter initialization
    theta.old[['pi']] <- c(0.999,0.001)
    theta.old[['gamma']] <- 
        estimateTransitionProb(chain = zPatterns$group$Group,numStates = K)
    theta.old[['psi']] <- 
        estimateCoefficients(z = zPatterns$group$Group,dt = dt,dist = dist,
                             type = 'consensus',control = control)
    # EM algorithm begins
    init.time <- verbose(level = 1,control = control)
    emConsensus <- emConsensus(theta.old,theta.new,dt,dtUnique,M,N,dist,
                               init.time,hdf5File,controlHist,parHist,control)
    verbose(level = 3,control = control)
    # Saving viterbi sequence & output file
    invisible(computeViterbiSequence(hdf5 = hdf5File,
                                     pi = emConsensus[['theta.new']]$pi,
                                     gamma = emConsensus[['theta.new']]$gamma))
    controlHistory <- 
        data.table::rbindlist(lapply(emConsensus[['controlHist']],function(x){
            data.table::data.table(t(unlist(x)))
        }))
    parHistory <- 
        data.table::rbindlist(lapply(emConsensus[['parHist']],function(x){
            data.table::data.table(t(unlist(x)))
        }))
    S4Vectors::metadata(object) <- 
        list('output' = hdf5File,'control' = control,
             'history' = list('control' = controlHistory,
                              'parameter' = parHistory))
    return(object)
}