#!/usr/bin/env Rscript

#' Author: Sumaiya Nazeen (sumaiya_nazeen@hms.harvard.edu)
#' This program runs NERINE for a given network and observed mutation counts
#' in a specific variant category for a case-control cohort
#'
#' @param case_count_file file containing case frequency table
#' @param control_count_file file containing control frequency table
#' @param fam_file tab-separated .fam file with case-control status of individuals
#' @param Network_matrix_file file containing the network adjacency matrix
#' @param lookup_table_file file containing the lookup table for network
#' @param num_alpha_levels 4 if pos-only test, 9 if pos-neg test
#' @param num_cores number of cpu cores to use
#' @return outfile_prefix.RDS a .RDS file containing the input and results snapshot of NERINE
#' @return outfile_prefix.txt a .txt file containing the results of NERINE

source("essential_functions.R")
library(emdbook)
library(eva)
library(DescTools)
library(dplyr)
library(doParallel)
library(pracma)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 8){
        stop("Insufficient arguments.
        Required parameters: 
        case_count_file
	control_count_file
	fam_file
        Network_matrix_file
	lookup_table_file
        num_alpha_levels
        outfile_prefix
        num_cores"
        , call.=FALSE)
} else {
        case_file <- args[1]
        control_file <- args[2]
	fam_file <-  args[3]
        nm_file <- args[4]
	lt_file <- args[5]
        rsense <- 1
	num_levels <- as.numeric(args[6])
	N <- 10000
	outfile_prefix <- args[7]
	num_cores <- as.numeric(args[8])
		
	Xcount <- read.table(case_file, sep='\t', header=F, row.names=1)
	colnames(Xcount) <- seq(1, length(colnames(Xcount)), 1)
	Ycount <- read.table(control_file, sep='\t', header=F, row.names=1)
	colnames(Ycount) <- seq(1, length(colnames(Ycount)), 1)
		
	Sigma <- read.table(nm_file, sep='\t', header=F, row.names=1)
	colnames(Sigma) <- rownames(Sigma)
	Sigma <- as.matrix(Sigma)
	Ngene <- nrow(Sigma)
	diag(Sigma) <- rep(2, Ngene) 
		
	fm <- read.table(fam_file, sep='\t', header=F)
	Ncase <- table(fm$V6)[[2]]
	Ncontrol <- table(fm$V6)[[1]]
		
	theta_try <- c(0.0000, 0.0001, 0.0005, 0.0010, 0.0020, 0.0040, 0.0080, 0.0100, 0.0200, 0.0300, 0.0400, 0.0500, 0.0600, 0.0700, 0.0800, 0.0900, 0.1000, 0.1500, 0.2000, 0.2500, 0.3000, 0.3500, 0.4000, 0.4500, 0.5000, 0.6000, 0.7000, 0.8000, 0.9000, 1.0000)
	lt <- readRDS(lt_file)
	print(head(lt))
	print(case_file)
	print(control_file)
	print(fam_file)
	print(nm_file)
	print("case mutations")
	print(Xcount)
	print("control mutations")
	print(Ycount)
	cl <- makeCluster(num_cores, type="FORK")
	registerDoParallel(cl)
	start <- proc.time()
	print(start)
	r <- foreach(th=theta_try, .combine=rbind) %dopar% {
		l <- calc_likelihood_pois_bin(Xcount, Ycount, Ncase, Ncontrol, Sigma, th, N, num_levels, lt, rsense)
		cbind(l, th)
	}
	stopCluster(cl)
	total_time <- proc.time() - start
	print(total_time)
	print(r)
	L <- as.vector(r[,1])
	ind <- which.max(L)
	theta_hat <-  theta_try[ind]
	LLR <- 2*(log(max(L))-log(L[1]))
	llrp <- pchibarsq(LLR, df=1, mix=0.5, lower.tail=F)
	c(alph, rel, maxl, l0, lralph, zl) %<-% find_mle_alpha_pois_bin_beta(Sigma, theta_hat, N, Xcount, Ycount, Ncase, Ncontrol, 
num_levels, lt)
	res_alph <- data.frame(cbind(rownames(Sigma), alph, rel))
	names(res_alph) <- c("gene", "alpha", "relative_effect")
		
	outrds <- sprintf('%s%s', outfile_prefix, ".RDS")
	outtxt <- sprintf('%s%s', outfile_prefix, ".txt")
	outtxt_genes <- sprintf('%s%s', outfile_prefix, "_gene_effects.txt")
	saveRDS(list(Xcount, Ycount, Sigma, theta_hat, LLR, llrp, res_alph, lralph), file=outrds)
		
	fprintf('Network: %s\n\n', nm_file, file=outtxt)
	fprintf('Estimated theta %0.4f\n', theta_hat, file=outtxt, append=T)
	fprintf('Log-likelihood ratio %0.4f\n', LLR, file=outtxt, append=T)
	fprintf('Asymptotic p-value %0.4e\n\n', llrp, file=outtxt, append=T)
	if(num_levels %% 2 == 1){
		fprintf('Individual gene efects (no-effect = %0.4f; pos-effect > %0.4f; neg-effect < %0.4f)\n', zl, zl, 
zl, file=outtxt, 
append=T)
	} else{
		fprintf('Individual gene efects (no-effect = %0.4f; pos-effect > %0.4f)\n', zl, zl, file=outtxt, append=T)
	}
	write.table(res_alph, outtxt, sep='\t', quote=F, row.names=F, append=T)
	write.table(res_alph, outtxt_genes, sep='\t', quote=F, row.names=F, append=F)
	fprintf('Log-likelihood ratio for gene alphas %0.2f\n', lralph, file=outtxt, append=T)
}

