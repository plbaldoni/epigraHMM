################################################################################
################################################################################
### R functions utilized in various epigraHMM scripts (brief desc. in header)
################################################################################
################################################################################

################################################################################
### FDR control method to call peaks based on posterior probabilities
################################################################################

fdrControl <- function(prob,fdr = 0.05){
    notpp = FDR = Window = .N = NULL
    if (!all(prob >= 0 & prob <= 1))
        stop('Posterior probabilities must be between 0 and 1')
    if (!(fdr > 0 & fdr < 1))
        stop('fdr must be between 0 and 1')
    
    probDT <- data.table::data.table(notpp = 1 - prob,
                                     Window = seq_len(length(prob)),
                                     key = 'Window')
    probDT <- probDT[order(notpp),]
    probDT[,FDR := ((cumsum(notpp) / seq_len(.N)) < fdr)]
    probDT <- probDT[order(Window),]
    return(probDT$FDR)
}

################################################################################
### Compute error for EM algorithm convergence check
################################################################################

getError <- function(controlHist,parHist,control){
    iteration <- controlHist[[length(controlHist)]][['iteration']]
    if (iteration > 1) {
        gap <- ifelse(iteration > control[['minIterEM']],control[['gapIterEM']], 1)
        parlist.old <- unlist(parHist[[iteration - gap]])
        parlist.new <- unlist(parHist[[iteration]])
        q.old <- controlHist[[iteration - gap]][['q']]
        q.new <- controlHist[[iteration]][['q']]
        return(c('MRCPE' = max(abs((parlist.new - parlist.old) / parlist.old)),
                 'MACPE' = max(abs(parlist.new - parlist.old)), 
                 'ARCEL' = max(abs((q.new - q.old) / q.old))))
    } else{
        return(c('MRCPE' = 0,'MACPE' = 0,'ARCEL' = 0))
    }
}

################################################################################
### Print message during EM algorithm
################################################################################

verbose = function(level,control,controlHist = NULL, theta = NULL){
    currentTime <- Sys.time()
    if (!control[['quiet']]) {
        if (level == 1) {
            message(rep('#',80))
            message(currentTime)
            message("Starting the EM algorithm")
            return(currentTime)
        }
        if (level == 2) {
            cLen <- length(controlHist)
            message(rep('#',80))
            message('\rIteration: ',controlHist[[cLen]][['iteration']])
            message('\rBIC: ',formatC(controlHist[[cLen]][['bic']], format = "e", digits = 2))
            message("\rError: ",paste(paste(names(controlHist[[cLen]][['error']]),'='),formatC(controlHist[[cLen]][['error']], format = "e", digits = 2),collapse = ', '))
            message("\rConvergence: ",paste(paste(names(controlHist[[cLen]][['count']]),'='),controlHist[[cLen]][['count']],collapse = ', '))
            message("\r",paste('Initial probabilities: '),paste(formatC(theta[['pi']], format = "e", digits = 2),collapse = ' '))
            message("\r",paste('Transition probabilities: '),paste(formatC(theta[['gamma']], format = "e", digits = 2),collapse = ' '))
            if ('delta' %in% names(theta)) {
                message("\r",paste('Mixture probabilities: '),paste(formatC(unlist(theta[['delta']]), format = "e", digits = 2),collapse = ' '))
            }
            message("\r",paste('Model coefficients: '),paste(formatC(unlist(theta[['psi']]), format = "e", digits = 2),collapse = ' '))
            message("\rTime elapsed (mins): ",formatC(controlHist[[cLen]][['time']], format = "f", digits = 2))
        }
        if (level == 3) {
            message(rep('#',80))
            message("EM algorithm converged!")
            message(currentTime)
            message(rep('#',80))
        }   
    }
}

################################################################################
### Melt data.table for plotting counts
################################################################################

meltCounts <- function(DT,ranges,object,subobject,peaks,annotation,subsetIdx) {
    Sample = Window = .N = NULL
    if (methods::is(ranges)[1] == "GRanges") {
        vars <- c('seqnames', 'start', 'end', 'width', 'strand')
        DT <- cbind(DT, as.data.table(rowRanges(object))[subsetIdx, vars, 
                                                         with = FALSE])
        DT[, Window := seq_len(.N)]
        DTmelt <- data.table::melt(DT,id.vars = c('Window', 'start'),
                                   measure.vars = seq_len(ncol(DT) - 6),
                                   value.name = 'Counts',
                                   variable.name = 'Sample')
        if (!is.null(peaks)) {
            DTmelt[Sample == levels(Sample)[1], {
                peaks := overlapsAny(subobject, peaks) 
            }]
        }
        if (!is.null(annotation)) {
            DTmelt[Sample == levels(Sample)[1], {
                annotation := overlapsAny(subobject, annotation)
            }]
        }
    } else{
        DT[, start := seq_len(nrow(object))[subsetIdx]]
        DT[, Window := seq_len(.N)]
        DTmelt <- data.table::melt(DT, id.vars = c('Window', 'start'),
                                   measure.vars = seq_len(ncol(DT) - 2),
                                   value.name = 'Counts',
                                   variable.name = 'Sample')
        if (!is.null(peaks)) {
            DTmelt[Sample == levels(Sample)[1], {
                peaks := peaks[subsetIdx] 
            }]
        }
        if (!is.null(annotation)) {
            DTmelt[Sample == levels(Sample)[1], {
                annotation := annotation[subsetIdx] 
            }]
        }
    }
    return(list('DT' = DT, 'DTmelt' = DTmelt))
}

################################################################################
### Transform data.table for plotting counts
################################################################################

aggregateCounts <- function(DT,ifDifferential,object) {
    .SD = NULL
    if (ifDifferential) {
        cond <- SummarizedExperiment::colData(object)$condition
        nameCol <- paste0(unique(cond))
        for (i in nameCol) {
            DT[, paste0(i) := rowSums(.SD), .SDcols = which(cond == i)]
        }
        DT <- DT[, nameCol, with = FALSE]
    } else{
        nameCol <- paste0('Replicate ',
                          SummarizedExperiment::colData(object)$replicate)
        data.table::setnames(DT, nameCol)
    }
    return(DT)
}

################################################################################
### Sort epigraHMM object
################################################################################

sortObject = function(epigraHMMDataSet){
    cData <- SummarizedExperiment::colData(epigraHMMDataSet)
    oData <- base::order(cData[,c('condition','replicate')],decreasing = FALSE)
    if (!all(oData == seq_len(nrow(cData)))) {
        epigraHMMDataSet <- epigraHMMDataSet[,oData]
        message("Rows of colData have been sorted with respect to conditions (and replicates). The resulting colData is:")
        message(paste0(utils::capture.output(SummarizedExperiment::colData(epigraHMMDataSet)), collapse = "\n"))
    }   
    return(epigraHMMDataSet)
}

################################################################################
### Invert parameters for glm.nb and glm.zinb
################################################################################

invertDispersion = function(par,model){
    if (model == 'nb') {
        return(c(par[seq_len((length(par) - 1))], 1 / par[length(par)]))
    }
    if (model == 'zinb') {
        if (length(par) == 3) {
            return(c(par[2], 1 / par[3], par[1]))
        } else{
            return(c(par[c(3, 4)], 1 / par[5], par[c(1, 2)]))
        }
    }
}

################################################################################
### Function to enumerate combinatorial patterns from initializing peaks
################################################################################

enumeratePatterns = function(object,group){
    Window = ref = .N = NULL
    uniGroup <- unique(group)
    # Peaks by group
    chain <- data.table::setDT(lapply(uniGroup,function(x){
        peakVec <- SummarizedExperiment::assay(object, 'peaks')[, which(group == x), drop = FALSE]
        1*(Matrix::rowSums(peakVec) > 0)
    }))
    data.table::setnames(chain,paste0('ChIP',seq_len(ncol(chain))))
    chain[,Window := seq_len(.N)]
    # Creating reference table
    ref[paste0('ChIP',seq_len(length(uniGroup)))] <- list(NULL)
    for (i in names(ref)) {
        ref[[i]] <- c(0, 1)
    }
    ref <- as.data.table(expand.grid(ref))
    ref <- ref[order(rowSums(ref)),]
    ref$Group <- seq_len(nrow(ref))
    chain <- merge(chain, ref,
                   by = paste0('ChIP', seq_len(length(uniGroup))),
                   all.x = TRUE)
    return(list('group' = chain[order(Window)][, c('Group'), with = FALSE], 
                'ref' = ref))
}

################################################################################
### Function to associate mixture components with combinatorial patterns
################################################################################

determineMixtures <- function(pattern,group,nGroup,ref){
    if (is.list(pattern)) {
        B <- length(pattern)
        z.seq <- lapply(seq_len(B),FUN = function(x){
            aux <- rep(FALSE,nGroup)
            aux[pattern[[x]]] <- TRUE
            refMat <- (as.matrix(ref[,seq_len(nGroup),with = FALSE]) == 1)
            which(apply(refMat,1,FUN = function(x){
                all(x == aux)
            }))
        })
    } else{
        if (is.null(pattern)) {
            B <- 2 ^ (length(unique(group))) - 2
            z.seq <- lapply(1 + seq_len(B),FUN = function(x){x})
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
    if (is.list(pattern)) {
        patList <- lapply(pattern,FUN = function(x){
            aux <- rep(FALSE, nGroup)
            aux[x] <- TRUE
            refMat <- as.matrix(ref[,seq_len(nGroup),with = FALSE]) == 1
            which(apply(refMat,1,FUN = function(x){
                all(x == aux)
            }))
        })
        model <- lapply(unlist(patList),FUN = function(x){
            cbind(rep(t(ref[x,seq_len(nGroup),with = FALSE]),
                      times = M*table(group)))
        })
    }
    if (is.null(pattern)) {
        model <- lapply(1 + seq_len(B),FUN = function(x){
            cbind(rep(t(ref[x,seq_len(nGroup),with = FALSE]),
                      times = M*table(group)))
        })
    }
    return(model)
}

################################################################################
### Estimate mixture probabilities for initialization
################################################################################

estimateMixtureProb = function(zDifferential,zSeq,B){
    mixList <- lapply(seq_len(B),FUN = function(i){
        subz.diff <- zDifferential[zDifferential %in% unlist(zSeq)]
        sum(subz.diff %in% zSeq[[i]] / length(subz.diff))
    })
    return(unlist(mixList))
}

################################################################################
### Estimate model coefficients for differential model initialization
################################################################################

estimateCoefficients <- function(z,dt,dist,type,control){
    Window = ChIP = offsets = mu = sigma2 = NULL
    # Naive estimation
    parList <- lapply(range(z),function(x){
        subpar <- dt[Window %in% which(z == x),{
            list(mu = mean((ChIP + 1) / exp(offsets)),
                 sigma2 = stats::var((ChIP + 1) / exp(offsets)))
        }]
        subpar[, c(zip = max(0.01, (sigma2 - mu) / (sigma2 + mu ^ 2 - mu)),
                   mu = mu,
                   disp = min((mu ^ 2) / max(0, sigma2 - mu), control[['maxDisp']]))]
    })
    # Adjusting parameters for the chosen distribution
    par <- list()
    if (type == 'differential') {
        if (dist == 'nb') {
            par[[1]] <- log(parList[[1]][-1])
            par[[2]] <- log(parList[[2]][-1]) - log(parList[[1]][-1])
        } else{
            par[[1]] <- log(parList[[1]])
            par[[2]] <- log(parList[[2]][-1]) - log(parList[[1]][-1])
        }
    } else{
        ifControls <- ifelse('controls' %in% names(dt),0,NA)
        if (dist == 'nb') {
            par[[1]] <- c(log(parList[[1]]['mu']),ifControls,
                          parList[[1]]['disp'])
            par[[2]] <- c(log(parList[[2]]['mu']),ifControls,
                          parList[[2]]['disp'])
        } else{
            par[[1]] <- c(log(parList[[1]]['zip']),ifControls,
                          log(parList[[1]]['mu']),ifControls,
                          parList[[1]]['disp'])
            par[[2]] <- c(log(parList[[2]]['mu']),ifControls,
                          parList[[2]]['disp'])
        }
    }
    return(lapply(par,function(x){unname(x[!is.na(x)])}))
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
### Get combinatorial patterns
################################################################################

getPatterns <- function(x,hdf5){
    patterns <- rhdf5::h5read(hdf5,'mixturePatterns')
    conditions <- unique(x$condition)
    out <- unlist(lapply(patterns,function(x){
        paste(conditions[as.numeric(gregexpr('E',x)[[1]])],collapse = '-')
    }))
    return(out)
}

################################################################################
### Print message during combinatorial pattern prunning
################################################################################

verbosePrunning <- function(level,control, info = NULL,
                            patternTable = NULL, idxReduced = NULL, BIC = NULL){
    if (isFALSE(control[['quietPruning']])) {
        if (level == 1) {
            message(rep('#',80))
            message('Full model (BIC = ',
                    formatC(info$BIC,digits = 2,format = 'f'),')\n')
            message(paste0(utils::capture.output(info$Components),
                           collapse = "\n"))
            message(rep('#',80))
            message('Prunning enrichment patterns\n')
        }
        if (level == 2) {
            p <- formatC(patternTable$PosteriorProportion[idxReduced],
                         digits = 2,format = 'e')
            message(paste0('- ',patternTable$Enrichment[idxReduced],
                           ' (p = ',p,')',' ... '),appendLF = FALSE)
        }
        if (level == 3) {
            bicChange <- formatC(100*(BIC - info$BIC)/info$BIC,
                                 digits = 2,format = 'f')
            message('% BIC rel. change = ',bicChange,'.')
        }
        if (level == 4) {
            message(rep('#',80))
            message('Reduced model (BIC = ',
                    formatC(info$BIC,digits = 2,format = 'f'),')\n')
            message(paste0(utils::capture.output(info$Components),
                           collapse = "\n"))
            message(rep('#',80))
        }   
    }
}

################################################################################
### Prune differential combinartorial patterns of enrichment
################################################################################

prunePatterns <- function(object,control,dist){
    Enrichment = NULL
    
    threshold <- control[['pruningThreshold']]
    control['pruningThreshold'] <- list(NULL)
    
    # Fit full model
    fitFull <- epigraHMM(object = object,control = control,type = 'differential',dist = dist)
    info_fitFull <- info(fitFull)
    patternTable <- info_fitFull$Components
    conditionTable <- 
        as.data.table(patternTable[!grepl("-",patternTable$Enrichment),-3])
    # Verbose
    verbosePrunning(level = 1,control = control,info = info_fitFull)
    
    # Pruning states
    while (any(patternTable$PosteriorProportion < threshold)) {
        # Removing meta
        unlink(metadata(fitFull)$output)
        # Find rarest pattern
        idx_Reduced <- which.min(patternTable$PosteriorProportion)
        subPatternTable <- patternTable[-idx_Reduced,]
        # Verbose
        verbosePrunning(level = 2,control = control,patternTable = patternTable,
                        idxReduced = idx_Reduced)
        # Setting combinatorial patterns
        control[['pattern']] <- 
            lapply(seq_len(nrow(subPatternTable)),function(x){
                comp <- as.character(subPatternTable[x,]$Enrichment)
                comp <- strsplit(comp,'-')[[1]]
                comp <- lapply(comp,function(z){
                    conditionTable[Enrichment == z,]$Component
                })
                unlist(comp)
            })
        # Fit reduced model
        fitReduced <- 
            epigraHMM(object = object,control = control,type = 'differential',dist = dist)
        info_fitReduced <- info(fitReduced)
        patternTable <- info_fitReduced$Components
        verbosePrunning(level = 3,control = control,BIC = info_fitReduced$BIC,
                        info = info_fitFull)
    }
    
    # Verbose
    verbosePrunning(level = 4,control = control,info = info_fitReduced)
    
    S4Vectors::metadata(fitReduced)$control[['pruningThreshold']] <- threshold
    
    return(fitReduced)
}
