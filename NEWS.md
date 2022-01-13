# epigraHMM 1.0.0

* First release of epigraHMM on Bioconductor.

* It is now possible to add normalizing offsets via `addOffsets`.

* epigraHMM now uses hdf5 files to store all intermediate data during computation
of the EM algorithm. Intermediate data include window-based HMM and mixture model 
posterior probabilities, and forward-backward probabilities. This change leads to
a better memory utilization of the package.

# epigraHMM 1.0.1

* Minor fix in the package DESCRIPTION file and version numbers

# epigraHMM 1.0.2

* epigraHMM now exports a function called `segmentGenome` that segments a given
genome (e.g. 'mm10') into non-overlapping genomic windows while considering 
gap tracks and blacklisted regions.

# epigraHMM 1.0.3

* Minor updates in the NEWS file as well as the README page.

# epigraHMM 1.0.4

* Adding function `callPatterns` to exp[ort] combinatorial patterns (or posterior 
probabilities) associated with a given set of genomic regions. 

* Adding function `info` to print summary statistics from epigraHMM output. This
function will print the model's BIC, log-likelihood, and combinatorial patterns
associated with mixture model components.

* Adding new example dataset `helas3` with ENCODE ChIP-seq data from broad 
epigenomic marks H3K27me3, H3K36me3, and EZH2.

* Adding option to prune combinatorial patterns associated with rare states. See
vignette for details.

* In differential peak calling, epigraHMM now exports combinatorial pattern 
table. See vignette for details.

* Improvement of the vignette to clarify epigraHMM's use of blacklisted regions
and gap tracks.

# epigraHMM 1.0.5

* Minor bug fix in callPatterns and info function (explict import of
S4Vectors::mcols and utils::tail).

* Exporting expStep function, which implements the E-step of EM algorithm
(forward-backward & Viterbi algorithm) for a K-state HMM.

# epigraHMM 1.0.6

* Minor bug fix in controlEM documentation

# epigraHMM 1.0.7

* Exporting maxStepProb, which compute the MLE of initial and transition probabilities of a K-state HMM, as well as simulateMarkovChain, which simulates a Markov chain of length 'n' given a matrix of transition probabilities

# epigraHMM 1.0.8

* Minor bug fix in maxStepProb documentations

# epigraHMM 1.2.0

* Second release of epigraHMM on Bioconductor.

# epigraHMM 1.2.1

* Bug fix of output paths to handle paths with '.'

# epigraHMM 1.3.2

* Exporting function estimateTransitionProb to estimate transition probabilities
from a sequence of states of a Markov chain
