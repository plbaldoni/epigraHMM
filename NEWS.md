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
