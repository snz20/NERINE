# NERINE
This code is associated with the following manuscript. If you use any part of the source code, please cite us:
Sumaiya Nazeen _et. al._ **NERINE reveals rare variant associations in gene networks across multiple phenotypes and implicates an _SNCA-PRL-LRRK2_ subnetwork in Parkinson’s disease** (submitted to _Nature Genetics_)

An abstract on this work was presented at ESHG 2024.

## Requirments
Python >= 3.12.4

R >= 4.3.2

Python packages: numpy, pandas, scipy.stats, bz2, gzip, seaborn, networkx

R packages: MASS, Matrix, matrixcalc, mvtnorm, sinib, prodlim, foreach, zeallot, plyr, ggplot2, parallel, foreach, doParallel, stringr, dplyr, emdbook, eva, DescTools

## Directory structure
./: contains the main scripts required to run NERINE

util/: python script containing utility functions

resources/: database files to facilitate network generation from 
protein interactions, co-expression, and co-essentiality

example_data/: toy data for doing a test run of NERINE

## Usage
Modes:
./NERINE_main.py [-h] {create_freqtable, generate_network, generate_lookup, run_NERINE}

1. Create freqtables for NERINE's input from a given gzipped VCF file (Mode = create_freqtable)

   ./NERINE_main.py create_freqtable [-h] [--glist GENELIST_ARG] [--mincutoff_arg MINCUTOFF_ARG] [--maxcutoff_arg MAXCUTOFF_ARG] in_vcf in_fam ftable_dir

2. Generate network adjacency matrix from a list of genes using database edges in absence of user-defined topology (Mode = generate_network)

   ./NERINE_main.py generate_network [-h] [--glist GENELIST_ARG] network_dir network_type resource_dir resource_prefix

4. Generate lookup table for NERINE's likelihood calculation (Mode = generate_lookup)

   ./NERINE_main.py generate_lookup [-h] [-k TESTTYPE_ARG] network_file lookup_dir

6. Run network-level rare variant association test with NERINE (Mode = run_NERINE)

   ./NERINE_main.py run_NERINE [-h] [--glist GENELIST_ARG] [-k TESTTYPE_ARG] [-n NUMCORE_ARG] [--mincutoff_arg MINCUTOFF_ARG] [--maxcutoff_arg MAXCUTOFF_ARG] ftable_dir network_file lt_file in_fam out_dir

Steps to be followed in a typical workflow is given in the workflow.txt file.

## Contact
Sumaiya Nazeen, sumaiya_nazeen@hms.harvard.edu
