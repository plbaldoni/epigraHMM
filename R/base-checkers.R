################################################################################
################################################################################
### Functions to check input, output, and intermediate epigraHMM objects
################################################################################
################################################################################

################################################################################
### Function to check counts for convergence
################################################################################

checkConvergence <- function(controlHist,control){
    if (control[['criterion']] == 'all') {
        return(max(controlHist[[length(controlHist)]][['count']]))
    } else{
        return(controlHist[[length(controlHist)]][['count']][[control[['criterion']]]])
    }
}

################################################################################
### Function to rename file output if it already exists
################################################################################

checkPath <- function(path){
    if(!(substr(path,nchar(path)-2,nchar(path)) == '.h5')){
        stop('The file name must have a .h5 extension')
    }
    
    if (!file.exists(path)) {
        return(path)
    }
    
    i <- 1
    dirName <- dirname(path)
    baseName <- basename(path)
    repeat {
        newBaseName <- paste0(substr(baseName,1,nchar(baseName)-3),i,'.h5')
        newPath <- file.path(dirName,newBaseName)
        if (!file.exists(newPath)) {
            return(newPath)
        }
        i <- i + 1
    }
}

################################################################################
### Check if rows of a matrix sum up to 1
################################################################################

checkProbabilities = function(P){
    P <- pmax(pmin(P,1),0)
    if (sum(P == 0) > 0) {
        P[P == 0] <- .Machine$double.xmin
    }
    return(P/rowSums(P))
}

################################################################################
### Check consistency of arguments in plotCounts.R
################################################################################

checkPlot = function(x,ranges){
    if (!is.null(x)) {
        grRanges <- methods::is(ranges)[1] == "GRanges" & 
            methods::is(x)[1] == "GRanges"
        cRanges <- methods::is(ranges)[1] %in% c("integer", "numeric") & 
                        methods::is(x)[1] == "logical"
        if (!(grRanges | cRanges)) {
            stop('If ranges is a GRanges object, then peaks/annotation must also be a GRanges object. If ranges is a numeric vector object, then peaks/annotation must be a vector of logicals')
        }
    }
}

################################################################################
### Check input matrix
################################################################################

checkInputMatrix <- function(countData,colData,rowRanges){
    # Checking if colData has the correct format
    if (!(is.data.frame(colData) &
          all(c('condition', 'replicate') %in% names(colData)))) {
        stop("The argument colData must be a data.frame with the columns 'condition' and 'replicate'")
    }
    # Checking whether replicates are unique
    if (any(table(colData$condition) > 1)) {
        uniqueN <- nrow(unique(colData[, c('condition', 'replicate')]))
        N <- nrow(colData[, c('condition', 'replicate')])
        if (uniqueN < N) {
            stop('The columns "condition" and "replicate" must uniquely represent your data')
        }
    }
    # Checking whether countData is a matrix or a list of matrices
    if (!is.matrix(countData)) {
        if (!(all(unlist(lapply(countData,is.matrix))) & 
             !is.null(names(countData)) & 
             (nrow(unique(do.call(rbind,lapply(countData,dim)))) == 1) & 
             ('counts' %in% names(countData)))) {
            stop("countData is not a proper argument, check the help manual.")
        } 
    } else{
        countData <- list('counts' = countData)
    }
    # Checking rowRanges
    if (!(methods::is(rowRanges)[1] == "GRanges" | is.null(rowRanges))) {
        stop("rowRanges must be a GRanges object") 
    }
    # Checking dimensions
    if (is.null(rowRanges)) {
        if (is.list(countData)) {
            if (!(nrow(colData) == unique(unlist(lapply(countData, ncol))))) {
                stop('Distinct dimensions of countData and colData are not allowed')
            }
        } else{
            if (!(nrow(colData) == ncol(countData))) {
                stop('Distinct dimensions of countData and colData are not allowed')
            }
        }
    } else{
        if (is.list(countData)) {
            if (!(nrow(colData) == unique(unlist(lapply(countData, ncol))) &
                  unique(unlist(lapply(countData, nrow))) == length(rowRanges))) {
                stop('Distinct dimensions of countData, colData, and rowRanges are not allowed')
            }
        } else{
            if (!(nrow(colData) == ncol(countData) &
                  nrow(countData) == length(rowRanges))) {
                stop('Distinct dimensions of countData, colData, and rowRanges are not allowed')
            }
        }
    }
    return(countData)
}

################################################################################
### Check input bam files
################################################################################

checkInputBam <- function(bamFiles,colData,genome,windowSize,
                          gapTrack,blackList) {
    # Checking if colData has the correct format
    if (!(is.data.frame(colData) &
          all(c('condition', 'replicate') %in% names(colData)))) {
        stop("The argument colData must be a data.frame with the columns 'condition' and 'replicate'")
    }
    # Checking whether replicates are unique
    if (any(table(colData$condition) > 1)) {
        uniqueN <- nrow(unique(colData[, c('condition', 'replicate')]))
        N <- nrow(colData[, c('condition', 'replicate')])
        if (uniqueN < N) {
            stop('The columns "condition" and "replicate" must uniquely represent your data')
        }
    }
    # Checking whether bamFiles is correct
    if (is.list(bamFiles) & ('counts' %in% names(bamFiles))) {
        bamList <- lapply(bamFiles, function(x) {
            return(!(is.character(x) &
                         length(x) == nrow(colData) &
                         all(file.exists(x)) &
                         all(file.exists(paste0(x, '.bai')))))
        })
        if (any(unlist(bamList))) {
            stop("bamFiles is not a proper argument, check the help manual.")
        }
    } else{
        if (!(is.character(bamFiles) &
              length(bamFiles) == nrow(colData) &
              all(file.exists(bamFiles)) &
              all(file.exists(paste0(bamFiles, '.bai'))))) {
            stop("bamFiles is not a proper argument, check the help manual.")
        } else{
            bamFiles <- list('counts' = bamFiles)
        }
    }
    # Checking whether windowSize is an integer
    if (!(is.numeric(windowSize) & windowSize %% 1 == 0)) {
        stop('The argument windowSize must be an integer number')
    }
    return(bamFiles)
}