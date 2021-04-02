
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epigraHMM

<!-- badges: start -->

<!-- badges: end -->

A fundamental task in the analysis of data resulting from epigenomic
sequencing assays is the detection of genomic regions with significant
or differential sequencing read enrichment. epigraHMM provides set of
tools to flexibly analyze data from a wide range of high-throughput
epigenomic assays (ChIP-seq, ATAC-seq, DNase-seq, etc.) in an end-to-end
pipeline. It includes functionalities for data pre-processing,
normalization, consensus and differential peak detection, as well as
data visualization. In differential analyses, epigraHMM can detect
differential peaks across either multiple conditions of a single
epigenomic mark (differential peak calling) or across multiple
epigenomic marks from a single condition (genomic segmentation). The
data pre-processing steps are heavily integrated with other Bioconductor
packages and allow the user to transform sequencing/alignment files into
count matrices that are suitable for the final analysis of the data. The
current implementation is optimized for genome-wide analyses of
epigenomic data and is efficient for the analysis under multi-sample
multiple-condition settings, as well as consensus peak calling in
multi-sample single-condition settings. epigraHMM uses two modified
versions of hidden Markov models (HMM) that are robust to the diversity
of peak profiles found in epigenomic data and are particularly useful
for epigenomic marks that exhibit short and broad peaks. Analyses can be
adjusted for possible technical artifacts present in the data and for
input control experiments, if available. Results from the peak calling
algorithms can be assessed using some of the available tools for
visualization that allow the inspection of detected peaks and read
counts.

## Installation

You can install the development version of epigraHMM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("plbaldoni/epigraHMM")
```

A
[vignette](https://github.com/plbaldoni/epigraHMM/blob/main/vignettes/epigraHMM.Rmd)
of epigraHMM with an overview of the package and its functionalities is
available. Users can build the vignette during the installation process
with the option `build_vignettes = TRUE` in the command above.
