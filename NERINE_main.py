#!/usr/bin/env python

'''
Runs NERINE to assess rare variant burden in gene networks for dichotomous traits. Based of the paper "NERINE reveals rare variant associations in gene networks across multiple phenotypes and implicates an SNCA-PRL-LRRK2 subnetwork in Parkinson’s disease" by Sumaiya Nazeen et al.

This is written by Sumaiya Nazeen <sumaiya_nazeen@hms.harvard.edu>.

This software requires working installations of Python (version >= 3.12.4) and R (version >= 4.3.2) to be available. Check the README for detailed information on package dependencies.
'''


# import required packages
from __future__ import print_function
__version__ = "0.1.0"

import os
from os import path, mkdir
from os.path import isdir
import glob
import argparse
import sys
import subprocess
import random
import time
import threading
import pandas as pd
import numpy as np
import csv
import shutil
import operator
from shutil import copyfile
from multiprocessing.dummy import Pool
from datetime import datetime
from collections import Counter
import util.utility_functions as uf

script_loc = os.path.realpath(__file__)
sys.path.append(os.path.join(os.path.dirname(script_loc),'util'))

# Setting up environment variables
my_env = os.environ.copy()

# Utility classes and functions
class ArgClass:
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    

def safe_makedirs(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
	else:
		print("Directory already exists!!")
		return 1
	return 0
    

# run count table generation
def create_ftable(in_vcf, in_fam, ftable_dir, args):
    '''Creates case control mutation counts tables from input gzvcf file for six variant categories: 

    in_vcf (string): Path to annotated gzipped vcf file with genotypes
    in_fam (string): Path to fam file with case-control status
    ftable_dir (string): output directory for case control mutation counts tables
    Unpacking args:
        genelist (string): Path to the list of genes in the network with one gene symbol per line
        mincutoff (float): minimum MAF cutoff for selecting qualifying variants. Default value 0.
        maxcutoff (float): maximum MAF cutoff for selecting qualifying variants. Default value 0.05.
    Returns case and control freqtables per category of variants for input to NERINE.
    '''
    print("inside_create_ftable",in_vcf, in_fam, ftable_dir, args.genelist_arg, args.mincutoff_arg, args.maxcutoff_arg)
    safe_makedirs(ftable_dir)
    #convert to mutations file
    tmp_bname = os.path.basename(in_vcf).strip("vcf.gz")
    tmp_mut = os.path.join(ftable_dir,tmp_bname)
    val = uf.convert_gzvcf_to_mutations(in_vcf, tmp_mut)
    if val != 0:
            print("conversion from gz.vcf to mutations.tsv file failed\n")
            exit(1)
    mutfile = tmp_mut + "_mutations.tsv"
    # subset to coding_mutations if the mutations file is too large
    # coding_mutfile = tmp_mut + "_coding_mutations.tsv"
    # subset_coding_mutations(mutfile, coding_mutfile)
    # create_freq_table(coding_mutfile, in_fam, genefile, min_cutoff, max_cutoff, out_prefix)
    genefile = args.genelist_arg
    min_cutoff = args.mincutoff_arg
    max_cutoff = args.maxcutoff_arg
    fam_bname = os.path.basename(in_fam).split(".")[0]
    out_prefix = '_'.join([tmp_mut, fam_bname, str(min_cutoff), str(max_cutoff)])
    val = uf.create_freq_table(mutfile, in_fam, genefile, min_cutoff, max_cutoff, out_prefix)
    if val != 0:
        print("generation of frequency tables from mutations.tsv failed\n")

# run network generation
def generate_network(network_dir, network_type, resource_dir, resource_prefix, args):
    '''Generates network topology for a gene set when bespoke network topology is not available: 

    network_dir (string): Path to output directory for networks
    network_type (int): 1 (dafault) = Physical & genetic interactions 
    from PPI database, 2 = co-expression in GTEx tissue, and 
    3 = co-essentiality in DepMap
    resource_dir (string): Path to directory containing resource files 
    needed for network generation
    resource_prefix (string): Prefix for resource file
    Unpacking args:
        genelist (string): Path to file containing the list of genes

    Returns the adjacency matrix of the gene network.
    '''
    print("inside_generate_network",network_dir, network_type, resource_dir, resource_prefix, args.genelist_arg)
    safe_makedirs(network_dir)
    if network_type == 1:
        genefile = args.genelist_arg
        bname = os.path.splitext(os.path.basename(genefile))[0]
        outfile = os.path.join(network_dir,bname+"_phy.tsv")
        uf.generate_network_phy(genefile, resource_dir, resource_prefix, outfile)
    elif network_type == 2:
        genefile = args.genelist_arg
        gctfile = os.path.join(resource_dir,resource_prefix+".gct.gz")
        bname = os.path.splitext(os.path.basename(genefile))[0]
        outfile = os.path.join(network_dir,bname+"_coexpression.tsv")
        uf.generate_network_coexpression(gctfile, genefile, outfile)
    elif network_type == 3:
        genefile = args.genelist_arg
        infile = os.path.join(resource_dir,resource_prefix+".tsv")
        bname = os.path.splitext(os.path.basename(genefile))[0]
        outfile = os.path.join(network_dir,bname+"_coessentiality.tsv")
        uf.generate_network_coessentiality(infile, genefile, outfile)
    else:
        print("Invalid network type\n")
        exit(1)
    return 0

# run lookup table generation
def generate_lookup(network_file, lookup_dir, args):
    '''Generates lookup table for a gene network: 

    network_file (string): Path to network file
    lookup_dir (string): Output directory for lookup table
    
    Unpacking args:
        testtype (int): 1 = genes can have only trait-increasing effect, 2 (default) = genes can have effects in both directions

    Returns R object containing the lookup table.
    '''
    print("inside_generate_lookup",network_file, lookup_dir, args.testtype_arg)
    safe_makedirs(lookup_dir)
    bname = os.path.splitext(os.path.basename(network_file))[0]
    ttype = args.testtype_arg
    alpha_levels = 9
    if ttype == 1:
        alpha_levels = 4
    elif ttype == 2:
        alpha_levels = 9
    else:
        print("Invalid test type\n")
        exit(1)
    out_prefix = os.path.join(lookup_dir,bname+"_l"+str(alpha_levels))
    command = network_file + ' ' + str(alpha_levels) + ' ' + str(10000) + ' ' + out_prefix
    rscript_path = os.path.dirname(script_loc)+'/genLookup.R'
    os.system('Rscript '+rscript_path+ ' ' + command)
    return 0

# run NERINE test
def run_nerine(ftable_dir, network_file, lt_file, in_fam, out_dir, args):
    '''Run NERINE to assess rare variant burden in a network: 

    ftable_dir (string): Path to directory containing case-control 
    mutation counts tables
    network_file (string): Path to network file
    lt_file (string): Path to lookup table file
    in_fam (string): Path to fam file with case-control status
    out_dir (string): Path to output directory
    
    Unpacking args:
        genelist (string): Path to file containing the list of genes
        testtype (int): 1 = genes can have only trait-increasing effect, 2 (default) = genes can have effects in both directions
        num_cores (int): number of parallel processors to use
        mincutoff (float): minimum MAF cutoff for selecting qualifying variants. Default value 0.
        maxcutoff (float): maximum MAF cutoff for selecting qualifying variants. Default value 0.05.
        

    Returns estimated network effect, log-likelihood ratio, and significance p-value as well as individual gene effects in .RDS and .txt files.
    '''
    print("running NERINE", ftable_dir, network_file, lt_file, in_fam, out_dir, args.testtype_arg, args.genelist_arg, 
args.numcore_arg, args.mincutoff_arg, args.maxcutoff_arg)
    safe_makedirs(out_dir)
    categories = ['LoF','damaging','damaging_missense','missense','neutral','synonymous']
    categories.sort()
    case_ftable_files = glob.glob(os.path.join(ftable_dir,"*case_freqtable*.tsv"))
    case_ftable_files.sort()
    control_ftable_files = glob.glob(os.path.join(ftable_dir,"*control_freqtable*.tsv"))
    control_ftable_files.sort()
    ttype = args.testtype_arg
    if ttype == 1:
        alpha_levels = 4
    elif ttype == 2:
        alpha_levels = 9
    else:
        print("Invalid test type\n")
        exit(1)
    genelist_file = args.genelist_arg
    num_cores = args.numcore_arg
    min_cutoff = args.mincutoff_arg
    max_cutoff = args.maxcutoff_arg
    for i in range(len(categories)):
        print("------ Analyzing "+ categories[i] + " ------")
        bname = os.path.splitext(os.path.basename(case_ftable_files[i]))[0].split("case")[0]
        out_prefix = os.path.join(out_dir, bname+categories[i])
        command = case_ftable_files[i] + ' ' + control_ftable_files[i] + ' ' + in_fam + ' ' + network_file + ' ' + lt_file + ' ' + str(alpha_levels) + ' ' + out_prefix + ' ' + str(num_cores)
        rscript_path = os.path.dirname(script_loc)+'/run_NERINE.R'
        os.system('Rscript '+rscript_path+ ' ' + command)
        print("------ Finished "+ categories[i] + " ------")
    
    return 0

# NERINE interface
def main(argv):
    parser = argparse.ArgumentParser(description='Assess rare variant burden in gene networks')
	# Shared arguments
    numcore_arg = ArgClass('-n', dest='numcore_arg', default=1, help='Number of parallel processors to be used', type=int)
    testtype_arg = ArgClass('-k', dest='testtype_arg', default=2, help='Test type: 1 = pos-only and 2 (default) = pos-neg', type=int)
    genelist_arg = ArgClass('--glist', dest='genelist_arg', help='path to input genelist')
    mincutoff_arg = ArgClass('--mincutoff_arg', help='minimum MAF cutoff for qualifying rare variants', type=float, default=0)
    maxcutoff_arg = ArgClass('--maxcutoff_arg', help='maximum MAF cutoff for qualifying rare variants', type=float, default=0.05)

    # Subparsers
    subparsers = parser.add_subparsers(help='sub-commands', dest='mode')

    # Count table args
    parser_count = subparsers.add_parser('create_freqtable', help='Generate case and control mutation count tables from input vcf for six variant categories: damaging, damaging_missense, LoF, missense, neutral, and synonymous', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_count.add_argument('in_vcf', help='Path to input gzvcf file')
    parser_count.add_argument('in_fam', help='Path to input fam file')
    parser_count.add_argument('ftable_dir', help='Output directory for case and control mutation count tables')
    parser_count.add_argument(*genelist_arg.args, **genelist_arg.kwargs)
    parser_count.add_argument(*mincutoff_arg.args, **mincutoff_arg.kwargs)
    parser_count.add_argument(*maxcutoff_arg.args, **maxcutoff_arg.kwargs)

	# Network generation args
    parser_network = subparsers.add_parser('generate_network', help='Prepare network file for test', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_network.add_argument('network_dir', help='Output directory for network file')
    parser_network.add_argument('network_type', type=int, default=1, help='1 (default): physical/genetic, 2: co-expression, 3: co-essentiality')
    parser_network.add_argument('resource_dir', help='Path to directory containing database files')
    parser_network.add_argument('resource_prefix', help='Prefix for resource files')
    parser_network.add_argument(*genelist_arg.args, **genelist_arg.kwargs)
	
    # Lookup table generation args
    parser_lookup = subparsers.add_parser('generate_lookup', help='Prepare lookp table for test', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_lookup.add_argument('network_file', help='Path for network file')
    parser_lookup.add_argument('lookup_dir', help='output directory for lookup table')
    parser_lookup.add_argument(*testtype_arg.args, **testtype_arg.kwargs)
    
    #  Test network for rare variant burden
    parser_nerine = subparsers.add_parser('run_NERINE', help='run NERINE to assess rare variant network effect', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_nerine.add_argument('ftable_dir', help='Path to directory containing case-control mutation counts tables')
    parser_nerine.add_argument('network_file', help='Path to network adjacency matrix file')
    parser_nerine.add_argument('lt_file', help='Path to lookup table file')
    parser_nerine.add_argument('in_fam', help='Path to tab-separated .fam file')
    parser_nerine.add_argument('out_dir', help='Path to output directory')
    parser_nerine.add_argument(*genelist_arg.args, **genelist_arg.kwargs)
    parser_nerine.add_argument(*testtype_arg.args, **testtype_arg.kwargs)
    parser_nerine.add_argument(*numcore_arg.args, **numcore_arg.kwargs)
    parser_nerine.add_argument(*mincutoff_arg.args, **mincutoff_arg.kwargs)
    parser_nerine.add_argument(*maxcutoff_arg.args, **maxcutoff_arg.kwargs)

    args=parser.parse_args(argv)
    print(args)
    sys.stdout.flush()

    mode = args.mode

    if mode == 'run_NERINE':
        st_time = datetime.now()
        print('starting NERINE')
        print("{:%Y-%m-%d %H:%M:%S}".format(st_time))
        run_nerine(args.ftable_dir, args.network_file, args.lt_file, args.in_fam, args.out_dir, args)
        print("Total full annot wall clock runtime (sec): {}".format((datetime.now() - st_time).total_seconds()))
        
    elif mode == 'create_freqtable':
        create_ftable(args.in_vcf, args.in_fam, args.ftable_dir, args)

    elif mode == 'generate_network':
        generate_network(args.network_dir, args.network_type, args.resource_dir, args.resource_prefix, args)

    elif mode == 'generate_lookup':
        generate_lookup(args.network_file, args.lookup_dir, args)


if __name__ == "__main__":
    main(sys.argv[1:])
