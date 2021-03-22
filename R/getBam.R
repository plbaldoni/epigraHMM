#' Download and rename bam files from chipseqDBData
#'
#' This function downloads and renames bam files from chipseqDBData to be used 
#' in the vignette of epigraHMM
#'
#' @param dir a (temporary) path where bam files will be stored
#' @param download a logical indicating whether or not files should be downloaded
#'
#' @return A GRanges object with differential peak calls in BED 6+3 format
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/epigraHMM}
#' 
#' @examples 
#' 
#' getBam(download = FALSE)
#'
#' @export
getBam <- function(dir = tempdir(),download = FALSE){
    
    if(download){
        # Downloading the data
        h3k9ac.paths <- chipseqDBData::H3K9acData()
        h3k4me3.paths <- chipseqDBData::H3K4me3Data()
        h3k27me3.paths <- chipseqDBData::H3K27me3Data()
        cbpd.paths <- chipseqDBData::CBPData()
        
        # Output table
        paths <- rbind(h3k9ac.paths,h3k4me3.paths,h3k27me3.paths,cbpd.paths)
        paths[['bam']] = NA
        paths[['index']] = NA
        
        # Looping datasets
        for(x in seq_len(nrow(paths))){
            bam.file <- paths$Path[[x]]
            
            bam.name <- file.path(dir,paste0(gsub(' ','_',paths$Name[x]),'.bam'))
            index.name <- file.path(dir,paste0(gsub(' ','_',paths$Name[x]),'.bam.bai'))
            
            file.copy(from = path.expand(bam.file$path),bam.name)
            file.copy(from = path.expand(bam.file$index),index.name)
            new.bam.file <- Rsamtools::BamFile(file = bam.name,index = index.name)
            
            paths$bam[[x]] <- new.bam.file$path
            paths$index[[x]] <- new.bam.file$index
        }
        
        return(paths)   
    }
}
