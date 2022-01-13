#' Estimate transition probability from a sequence of integers
#'
#' This function estimates the transition probabilities for a k-state Markov
#' chain based on a sequence of integers that represent states of the chain
#'
#' @param chain a vector of integers
#' @param numStates an integer, the number of states in the Markov chain
#' 
#' @return A k-by-k matrix of transition probabilities, such that k is the 
#' number of states of the chain
#'
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom SummarizedExperiment assayNames assay
#'
#' @examples 
#' 
#' trueMat <- matrix(c(0.9,0.1,0.1,0.9),2,2)
#' simChain <- simulateMarkovChain(trueMat,1e3)
#' estMat <- estimateTransitionProb(simChain,2)
#' 
#' # estMat should be close to trueMat
#' estMat
#'
#' @export
estimateTransitionProb = function(chain,numStates){
  nm1 <- (numStates - 1)
  if (max(chain) > nm1) {
    chain <- round(nm1*(chain-max(chain))/(max(chain)-min(chain))+nm1)
  }
  MC <- matrix(chain, nrow = 1, ncol = length(chain))
  MC <- table(c(MC[, -ncol(MC)]), c(MC[, -1]))
  MC <- as.matrix(MC / rowSums(MC))
  MC <- matrix(MC,ncol = ncol(MC),nrow = nrow(MC),byrow = FALSE)
  return(checkProbabilities(MC))
}