# StringDecomposer.Benchmarking scripts

The repository contains benchmarking scripts for testing String Decomposition Problem solutions described in the manuscript "The String Decomposition Problem and its Applications to Centromere Assembly".


## Solvers

### alpha-CENTAURI HMMer (AC)

AC uses `chop_to_monomer.py` script from the original alpha-CENTAURI pipeline to generate monomer positions in given long reads, after that it applies `convert_identites.py` script to choose best monomer for each place. 

### StringDecomposer (SD)

SD uses dynamic programming to generate best decomposition of each read into the given monomers. 


## Scripts description

All scripts, except `AC_HMMer.snake` and `convert_identites.py`, use Python3. Required libraries are listed in `requirement.txt`
	
	ls scripts
	
	AC_HMMer.snake                     snakemake file, runs AC pipeline (uses Python2)
	convert_identities.py              identifies best monomers for the decomposition (part of AC pipeline)

	compare_tools_monomerfree.py       generates statistics for Appendix "Monomer-free benchmarking"

	generate_stats.py                  generates monoread-to-monocentromere alignment and calculates overall statistics (Table 1)
	process_mismatches.py              generates mismatches statistics based on generate_stats.py output
	process_abnormal.py                generates statistics for abnormal HORs based on generate_stats.py output
	draw_identities.py                 MinAlignmentIdentity statistics
	hor_analysis/					   script for HOR analysis with HiFi reads from T2T project
	parameters_tuning/				   scripts for building logistic regression

## Datasets

All data used in this work can be found on [T2T github page](https://github.com/nanopore-wgs-consortium/chm13) or on [Figshare](https://figshare.com/s/076674f298ce7d67701b).


## Citation

The String Decomposition Problem and its Applications to Centromere Assembly. *Tatiana Dvorkina, Andrey V. Bzikadze, Pavel A. Pevzner* bioRxiv 2019.12.26.888685; doi: [https://doi.org/10.1101/2019.12.26.888685](https://doi.org/10.1101/2019.12.26.888685)

## Contact

In case of any issues please use [issue tracker](https://github.com/ablab/stringdecomposer/issues) or email directly to [t.dvorkina@spbu.ru](mailto:t.dvorkina@spbu.ru)


