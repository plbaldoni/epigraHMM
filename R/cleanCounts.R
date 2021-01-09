#' Remove effects from covariates of interest
#'
#' This function removes the effect from covariates of interest (such as GC content) from experimental counts
#'
#' @param object an epigraHMMDataSet
#' @param effectNames a character vector with the names of assays for which the effect will be removed from the experimental counts.
#' Names in `effectNames` must be assays stored in the epigraHMMDataSet `object`.
#' @param byNames a character vector with the name of an assay containing stratification variables which will be used to define stratum-specific effects.
#' Examples of byNames assays include the 'peaks' assay from `initializer()`. In this case, models will be fit separately for peaks and non-peaks regions.
#' This can be usefor for effects such as GC content, which are known to have a differential effect between peaks and non-peak regions.
#' Default is NULL, i.e., effects will be removed without stratification.
#'
#' @details
#'
#' To be added
#'
#' @return An epigraHMMDataSet with an 'offset' assay filled in.
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#'
#' @importFrom stats predict
#' @importFrom MASS glm.nb
#'
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @export
cleanCounts <- function(object,effectNames,byNames = NULL){

    # Checking input
    if (!((methods::is(object)[1]=='RangedSummarizedExperiment') &
        all(effectNames %in% assayNames(object)) & (byNames %in% assayNames(object)))){
        stop('Check argments')
    }

    # Looping through samples
    offsets <- do.call(cbind,lapply(seq_len(ncol(assay(object,'counts'))),function(x){

        offset <- rep(0,nrow(object))

        ## Putting counts and covariates in a data.frame
        mat <- data.frame(y = assay(object,'counts')[,x],
                          X = do.call(cbind,lapply(seq_len(length(effectNames)),function(w){assay(object,effectNames[w])[,x]})),
                          by = if (is.null(byNames)) rep(TRUE,nrow(object)) else assay(object,byNames)[,x],
                          offsets = if ('offsets' %in% assayNames(object)) assay(object,'offsets')[,x] else rep(0,nrow(object)))

        ## Writing model formula
        form <- stats::as.formula(paste('y ~',paste0('splines::ns(',colnames(mat)[grep('^X.*',colnames(mat))],', df = 2)',collapse = ' + ')))

        ## Looping by
        for(index in unique(mat$by)){
            fit <- MASS::glm.nb(form,data = mat[mat$by == index,])
            pred <- stats::predict(fit,type = "terms")
            offset[mat$by == index] <- rowSums(pred[,grep('s\\(',colnames(pred)),drop = FALSE])
        }

        ## Returning offsets
        return(offset)
    }))

    # Adding offsets
    object <- addOffsets(object,offsets)

    return(object)
}
