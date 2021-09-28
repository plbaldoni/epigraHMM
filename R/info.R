#' Get information about peak calling results
#'
#' This function returns the BIC and expected log-likelihood function of the model,
#' with respect to the last conditional distribution of unknown enrichment peaks given the data.
#' The latter is also known as 'Q-function' in the EM context.
#'
#' @param object an epigraHMMDataSet
#'
#' @return
#'
#' A list with BIC, and expected log-likelihood function of the model. If
#' the input object contains results from a differential analysis, `info` will
#' also output the enrichment patterns associated with each mixture component
#' used in the mixture model.
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#' 
#' @examples 
#' 
#' # Creating dummy object
#' countData <- rbind(matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1),
#'                    matrix(rnbinom(2e3,mu = 7.5,size = 5),ncol = 1),
#'                    matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1))
#' 
#' colData <- data.frame(condition = 'A', replicate = 1)
#' object <- epigraHMMDataSetFromMatrix(countData,colData)
#' 
#' # Initializing
#' object <- initializer(object,controlEM())
#' 
#' # Running epigraHMM
#' object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'nb')
#' 
#' # Get info
#' info(object)
#'
#' @export
info <- function(object){
    
    out <- list('BIC' = utils::tail(metadata(object)$history$control$bic,1),
                'logLikelihood' = utils::tail(metadata(object)$history$control$q,1))
    
    if(length(unique(object$condition))>1){
        out[['Components']] = metadata(object)$components
    }
    return(out)
}
