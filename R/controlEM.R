#' Control parameters for the EM algorithm from epigraHMM
#'
#' This function passes controlling parameters for the EM algorithm implemented in the epigraHMM package.
#'
#' @param epsilonEM a named vector of positive values specifying up to four possible convergence criterion tolerances for the EM algorithm (see 'criterion' below). Default is c('MRCPE' = 1e-4, 'MACPE' = 1e-4,'ARCEL' = 1e-4).
#' @param maxIterEM a positive integer giving the maximum number of EM iterations. Default is 500.
#' @param minIterEM a positive integer giving the minimum number of EM iterations to start evaluating the convergence. Default is 3.
#' @param gapIterEM a positive integer giving the number of EM iterations apart to compute the convergence criterion. Default is 3.
#' @param maxCountEM a positive integer giving the number of consecutive EM iterations satisfying the convergence criterion in order to stop the algorithm. Default is 3.
#' @param maxDisp a positive value for the upper limit constraint of the dispersion parameters. Default is 1000.
#' @param criterion a character specifying the convergence criterion. Either "MRCPE" (maximum absolute relative change in parameter estimates), "MACPE" (maximum absolute change of parameter estimates),
#' "ARCEL" (absolute relative change of the Q-function), "ACC" (agreement of Viterbi peak calls), or "all" (simultaneously check for MRCPE, MACPE, ARCEL, and ACC).
#' For ACC, it computes the window-based percentage of Viterbi predictions unchanged 'gapIterEM' iterations apart. Default is "all".
#' @param minZero a positive value for the minimum positive value allowed in computations to avoid having zeros. Default is .Machine$double.xmin.
#' @param probCut a number between 0 and 1 for the cutoff of the rejection controlled EM algorithm. Default 0.05.
#' @param quiet a logical indicating whether to print messages. Default is TRUE.
#' @param maxIterInnerEM a positive integer giving the maximum number of inner EM iterations. Default is 5.
#' @param epsilonInnerEM a positive value with the convergence tolerance value for the inner EM algorithm. The criterion for the inner EM is "MRCPE". Default is 1e-4.
#' @param trimOffset either NULL or a positive integer indicating the number of decimal places to be used in the offset. Default is 3.
#' @param random either 'intercept' or 'slope'. It specifies the type of random effects model for consensus peak calling. Default is 'intercept'.
#' @param maxSigma2 a positive value for the maximum positive value of the variance component allowed in computations. Default is 10.
#' @param minSigma2 a positive value for the minimum positive value of the variance component allowed in computations. Default is 1e-8.
#' @param pattern either NULL (the default) or a list with length equal to the number of differential patterns to be modeled by the differential HMM state. See Details section below.
#' @param tempDir a string where intermediate results will be saved. Default is `tempdir()`.
#'
#' @details
#' If \code{pattern} is NULL, every possible combinatorial pattern will be considered. If \code{pattern} is a list, elements of it should specify the differential patterns to be modeled by each mixture component.
#' For instance, if pattern = list(2,c(1,3)) the mixture model will have two components that will represent the enrichment of condition 2 alone and the enrichment of conditions 1 and 3 together.
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
controlEM = function(epsilonEM=c('MRCPE' = 1e-4, 'MACPE' = 1e-4,'ARCEL' = 1e-4),
                     maxIterEM=500,
                     minIterEM=3,
                     gapIterEM=3,
                     maxCountEM=3,
                     maxDisp=1000,
                     criterion='all',
                     minZero=.Machine$double.xmin,
                     probCut=0.05,
                     quiet=TRUE,
                     maxIterInnerEM = 5,
                     epsilonInnerEM = 1e-4,
                     trimOffset = 3,
                     random='intercept',
                     maxSigma2 = 10,
                     minSigma2 = 1e-8,
                     pattern = NULL,
                     tempDir = tempdir()) {
    
    # Checks for epsilonEM
    
    if (!(is.numeric(epsilonEM) & all(epsilonEM>0))){stop("epsilonEM must be a positive numerical value (or vector)")}
    
    # Check names for epsilonEM
    internalEpsilonEM <- c('MRCPE' = 1e-4, 'MACPE' = 1e-4,'ARCEL' = 1e-4)
    
    if(!is.null(names(epsilonEM))){
        for(i in names(epsilonEM)[names(epsilonEM) %in% names(internalEpsilonEM)]){
            internalEpsilonEM[i] <- epsilonEM[i]
        }
    }
    epsilonEM <- internalEpsilonEM
    
    # Checks for maxit.em
    
    if (!maxIterEM%%1==0 || maxIterEM <= 0){stop("value of 'maxIterEM' must be a positive integer")}
    
    # Checks for minIterEM
    
    if (!minIterEM%%1==0 || minIterEM <= 0){stop("value of 'minIterEM' must be a positive integer")}
    
    # Checks for gapIterEM
    
    if (!gapIterEM%%1==0 || gapIterEM <= 0 || gapIterEM>minIterEM){stop("value of 'gapIterEM' must be a positive integer <= minIterEM")}
    
    # Checks for maxCountEM
    
    if (!maxCountEM%%1==0 || maxCountEM <= 0){stop("value of 'maxCountEM' must be a positive integer")}
    
    # Checks for maxDisp
    
    if (!is.numeric(maxDisp) || maxDisp <= 0){stop("value of 'maxDisp' must be > 0")}
    
    # Checks for criterion
    
    if (any((length(criterion)==1 & criterion %in% c('MRCPE','MACPE','ARCEL','ACC','all'))==FALSE)){stop("value of 'criterion' must be 'MRCPE', 'MACPE', 'ARCEL', 'ACC', or 'all'")}
    
    # Checks for minZero
    
    if (!(is.numeric(minZero) & minZero>0 & length(minZero)==1)){stop("'minZero' must be a positive numeric value")}
    
    # Checks for probCut
    
    if (!(length(probCut)==1 & is.numeric(probCut) & probCut>0 & probCut<1)){stop("'probCut' must be a value between 0 and 1")}
    
    # Checks for quiet
    
    if (!(length(quiet)==1 & is.logical(quiet))){stop("'quiet' must be a logical value")}
    
    # Checks for maxIterInnerEM
    
    if (!maxIterInnerEM%%1==0 || maxIterInnerEM <= 0){stop("value of 'maxIterInnerEM' must be a positive integer")}
    
    # Checks for epsilonInnerEM
    
    if (!(length(epsilonInnerEM)==1 & is.numeric(epsilonInnerEM) & epsilonInnerEM>0)){stop("'epsilonInnerEM' must be a positive value")}
    
    # Checks for minSigma2
    
    if (!(length(minSigma2)==1 & is.numeric(minSigma2) & minSigma2>0 & minSigma2 < maxSigma2)){stop("'minSigma2' must be a positive value")}
    
    # Checks for maxSigma2
    
    if (!(length(maxSigma2)==1 & is.numeric(maxSigma2) & maxSigma2>0 & maxSigma2 > minSigma2)){stop("'maxSigma2' must be a positive value")}
    
    # Checks for random
    
    if (!(random %in% c('intercept','slope'))){stop("random should be 'intercept' or 'slope'.")}
    
    # Checks for trimOffset
    
    if (!is.null(trimOffset)){if (!trimOffset%%1==0 || trimOffset <= 0){stop("value of 'trimOffset' must be a positive integer")}}
    
    # Checks for pattern
    
    if (!is.null(pattern)){if (!is.list(pattern)){stop("'pattern' must be NULL or a list'")}}
    
    # Checks for maxIterEM
    
    if (!dir.exists(tempDir)){stop("Temporary directory tempDir does not exist. Create it prior to execution.")} else{tempDir <- normalizePath(tempDir)}
    
    # Return
    
    list(epsilonEM=epsilonEM,
         minZero=minZero,
         maxIterEM=maxIterEM,
         minIterEM=minIterEM,
         gapIterEM=gapIterEM,
         maxCountEM=maxCountEM,
         maxDisp=maxDisp,
         criterion=criterion,
         probCut=probCut,
         quiet=quiet,
         maxIterInnerEM=maxIterInnerEM,
         epsilonInnerEM=epsilonInnerEM,
         trimOffset=trimOffset,
         random = random,
         maxSigma2 = maxSigma2,
         minSigma2 = minSigma2,
         pattern = pattern,
         tempDir = tempDir)
}
