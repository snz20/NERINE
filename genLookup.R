#!/usr/bin/env Rscript

source("essential_functions.R")
library(dplyr)

#' Author: Sumaiya Nazeen (sumaiya_nazeen@hms.harvard.edu)
#' This program generates the pairwise lookup table for 
#' a given gene network
#'
#' @param Network_matrix_file a .tsv file containing the network adjacency matrix.
#' @param num_alpha_levels 4 if posonly test, 9 if posneg test
#' @param Niter number of random samples to draw from the multivariate normal distribution
#' @return outfile_prefix.RDS a .RDS file containing the lookup table for input to NERINE

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
        stop("Insufficient arguments.
        Required parameters: 
        Network_matrix_file
        num_alpha_levels
        Niter
        outfile_prefix"
        , call.=FALSE)
} else {
	nm_file <- args[1]
	num_levels <- as.numeric(args[2])
	N <- as.numeric(args[3])
	outfile_prefix <- args[4]
	rsense <- 1
	
	Sigma <- read.table(nm_file, sep='\t', header=F, row.names=1)
	colnames(Sigma) <- rownames(Sigma)
	Sigma <- as.matrix(Sigma)
	Ngene <- nrow(Sigma)
	diag(Sigma) <- rep(2, Ngene) 
	Sigma[is.na(Sigma)] <- 0
	
	theta_try <- c(0.0000, 0.0001, 0.0005, 0.0010, 0.0020, 0.0040, 0.0080, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.3500, 0.4000, 0.4500, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000)
	lt <- make_alph_lookup_bin3(Sigma, theta_try, N, rsense, num_levels=num_levels)
	outrds <- sprintf('%s%s', outfile_prefix, ".RDS")
	saveRDS(lt, file=outrds)
}

