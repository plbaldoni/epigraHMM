#' Control parameters for the EM algorithm from epigraHMM
#'
#' This function passes controlling parameters for the EM algorithm implemented in the epigraHMM package.
#'
#' @param epsilonEM a named vector of positive values specifying up to four possible convergence criterion tolerances for the EM algorithm (see 'criterion' below). Default is c('MRCPE' = 1e-3, 'MACPE' = 1e-3,'ARCEL' = 1e-3).
#' @param maxIterEM a positive integer giving the maximum number of EM iterations. Default is 500.
#' @param minIterEM a positive integer giving the minimum number of EM iterations to start evaluating the convergence. Default is 3.
#' @param gapIterEM a positive integer giving the number of EM iterations apart to compute the convergence criterion. Default is 3.
#' @param maxCountEM a positive integer giving the number of consecutive EM iterations satisfying the convergence criterion in order to stop the algorithm. Default is 3.
#' @param maxDisp a positive value for the upper limit constraint of the dispersion parameters. Default is 1000.
#' @param criterion a character specifying the convergence criterion. Either "MRCPE" (maximum absolute relative change in parameter estimates), "MACPE" (maximum absolute change of parameter estimates),
#' "ARCEL" (absolute relative change of the Q-function), or "all" (simultaneously check for MRCPE, MACPE, and ARCEL). Default is "all".
#' @param minZero a positive value for the minimum positive value allowed in computations to avoid having zeros. Default is .Machine$double.xmin.
#' @param probCut a number between 0 and 1 for the cutoff of the rejection controlled EM algorithm. Default 0.05.
#' @param quiet a logical indicating whether to print messages. Default is TRUE.
#' @param maxIterInnerEM a positive integer giving the maximum number of inner EM iterations. Default is 5.
#' @param epsilonInnerEM a positive value with the convergence tolerance value for the inner EM algorithm. The criterion for the inner EM is "MRCPE". Default is 1e-3.
#' @param trimOffset either NULL or a positive integer indicating the number of decimal places to be used in the offset. Default is 3.
#' @param pattern either NULL (the default) or a list with length equal to the number of differential patterns to be modeled by the differential HMM state. See Details section below.
#' @param tempDir a string where results will be saved. Default is `tempdir()`.
#' @param fileName a string with the name of the result files. Default is `epigraHMM`.
#' @param pruningThreshold a numeric value between 0 and 1 to consider when pruning rare combinatorial patterns. Default is NULL (see Details).
#' @param quietPruning a logical indicating whether to print messages during the pruning step. Default is TRUE.
#'
#' @details
#' If \code{pattern} is NULL, every possible combinatorial pattern will be considered. If \code{pattern} is a list, elements of it should specify the differential patterns to be modeled by each mixture component.
#' For instance, if pattern = list(2,c(1,3)) the mixture model will have two components that will represent the enrichment of condition 2 alone and the enrichment of conditions 1 and 3 together.
#' 
#' If \code{pruningThreshold} is a value between 0 and 1, say 0.05, epigraHMM 
#' will sequentially remove differential combinatorial patterns of enrichment 
#' from any mixture model component with associated posterior mixture proportion
#' less than 0.05.
#' 
#' @return A list with components equal to the arguments
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @examples
#' # No more than 100 EM iterations
#' control <- controlEM(maxIterEM = 100)
#'
#' @export
controlEM = function(epsilonEM=c('MRCPE' = 1e-3, 'MACPE' = 1e-3,'ARCEL' = 1e-3),maxIterEM=500,
                     minIterEM=3,gapIterEM=3,maxCountEM=3,maxDisp=1000,criterion='all',
                     minZero=.Machine$double.xmin,probCut=0.05,quiet=TRUE,maxIterInnerEM = 5,
                     epsilonInnerEM = 1e-3,trimOffset = 3,pattern = NULL,tempDir = tempdir(),fileName = 'epigraHMM',
                     pruningThreshold = NULL, quietPruning = TRUE) {
    
    # Check names for epsilonEM
    internalEpsilonEM <- c('MRCPE' = 1e-3, 'MACPE' = 1e-3,'ARCEL' = 1e-3)
    
    if(!is.null(names(epsilonEM)))
        for(i in names(epsilonEM)[names(epsilonEM) %in% names(internalEpsilonEM)])
            internalEpsilonEM[i] <- epsilonEM[i]
    
    epsilonEM <- internalEpsilonEM
    
    # Checking input
    if (!(is.numeric(epsilonEM) & all(epsilonEM>0))){stop("epsilonEM must be a positive numerical value (or vector)")}
    
    if ((!maxIterEM%%1==0 || maxIterEM <= 0) | (!minIterEM%%1==0 || minIterEM <= 0)){stop("maxIterEM and minIterEM must be a positive integer")}
    
    if (!gapIterEM%%1==0 || gapIterEM <= 0 || gapIterEM>minIterEM){stop("value of 'gapIterEM' must be a positive integer <= minIterEM")}
    
    if (!maxCountEM%%1==0 || maxCountEM <= 0){stop("value of 'maxCountEM' must be a positive integer")}
    
    if (!is.numeric(maxDisp) || maxDisp <= 0){stop("value of 'maxDisp' must be > 0")}
    
    if (any((length(criterion)==1 & criterion %in% c('MRCPE','MACPE','ARCEL','ACC','all'))==FALSE)){stop("value of 'criterion' must be 'MRCPE', 'MACPE', 'ARCEL', 'ACC', or 'all'")}
    
    if (!(is.numeric(minZero) & minZero>0 & length(minZero)==1)){stop("'minZero' must be a positive numeric value")}
    
    if (!(length(probCut)==1 & is.numeric(probCut) & probCut>0 & probCut<1)){stop("'probCut' must be a value between 0 and 1")}
    
    if (!(length(quiet)==1 & is.logical(quiet))){stop("'quiet' must be a logical value")}
    
    if (!maxIterInnerEM%%1==0 || maxIterInnerEM <= 0){stop("value of 'maxIterInnerEM' must be a positive integer")}
    
    if (!(length(epsilonInnerEM)==1 & is.numeric(epsilonInnerEM) & epsilonInnerEM>0)){stop("'epsilonInnerEM' must be a positive value")}
    
    if (!is.null(trimOffset)){if (!trimOffset%%1==0 || trimOffset <= 0){stop("value of 'trimOffset' must be a positive integer")}}
    
    if (!is.null(pattern)){if (!is.list(pattern)){stop("'pattern' must be NULL or a list'")}}
    
    if (!dir.exists(tempDir)){stop("Temporary directory tempDir does not exist. Create it prior to execution.")} else{tempDir <- normalizePath(tempDir)}
    
    if (!is.null(pruningThreshold)){if (pruningThreshold <= 0 || pruningThreshold >= 1){stop("value of 'pruningThreshold' must be between (0,1)")}}
    
    if (!(length(quietPruning)==1 & is.logical(quietPruning))){stop("'quietPruning' must be a logical value")}
    
    # Return
    list(epsilonEM=epsilonEM,minZero=minZero,maxIterEM=maxIterEM,minIterEM=minIterEM,
         gapIterEM=gapIterEM,maxCountEM=maxCountEM,maxDisp=maxDisp,criterion=criterion,
         probCut=probCut,quiet=quiet,maxIterInnerEM=maxIterInnerEM,epsilonInnerEM=epsilonInnerEM,
         trimOffset=trimOffset,pattern = pattern,tempDir = tempDir,fileName = fileName,
         pruningThreshold=pruningThreshold,quietPruning=quietPruning)
}
