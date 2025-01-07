library(MASS)
library(Matrix)
library(matrixcalc)
library(mvtnorm)
library(sinib)
library(prodlim)
library(foreach)
library(zeallot)
library(plyr)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
library(stringr)
library(dplyr)

find_interval <- function(val, num_levels){
	breakpoints <- seq(0,1,1/num_levels)
	nb <- length(breakpoints)
	for(b in 2:(nb-1)){
		if(val <= breakpoints[b]){
			return(b-1)
		}
	}
	return(nb-1)
}

map_custom3 <- function(val, th, num_levels){
	if(th == 0){
		if(num_levels%%2 == 1){
			return(round(rep(0.5, length(val)),2))
		} else{
			return(round(rep(0.1, length(val)),2))
		}
	}
		
	al <- 0.5/th
	#x <- rnorm(10000, 0, th)
	y <- rbeta(10000, al, al)
	ty <- cut_interval(y, num_levels)
	xp <- quantile(val, probs=cumsum(as.numeric(table(ty)/10000)))
	lower <- findInterval(val,xp)
	upper <- lower+1
	upper[upper>num_levels] <- num_levels
	gval <- upper*0.1
	gval[gval==1] <- num_levels*0.1
	return(round(gval,1))
}


map_custom4 <- function(val, th, num_levels){
	if(th == 0){
		if(num_levels%%2 == 1){
			return(round(rep((floor(num_levels/2)+1), length(val)),2))
		} else{
			return(round(rep(1, length(val)),2))
		}
	}
		
	al <- 0.5/th
	bl <- 0.5/th
	#x <- rnorm(10000, 0, th)
	y <- rbeta(10000, al, bl)
	v <- get_levels(num_levels)
	if(num_levels%%2 == 1){
		#v <- c(seq(0,0.5,1/num_levels),seq(0.5,1,1/num_levels))
		tty <- cut_interval(y, num_levels)
		xp <- quantile(val, probs=cumsum(as.numeric(table(tty)/10000)))
		#lower <- findInterval(val,xp)
		#upper <- lower +1
		#upper[upper>num_levels] <- num_levels
		#ggval <- v[upper]
		lower <- findInterval(val,xp)
		upper <- lower+1
		upper[upper>num_levels] <- num_levels
		lvly <- round(seq(min(y),max(y),(max(y)-min(y))/(num_levels-1)),1)
		gval <- lvly[upper]
		l2 <- findInterval(gval, v)
		l2[l2==0] <- 1
		l2[l2>num_levels] <- num_levels
		
	} else{
		#v <- seq(0.5,1,1/(2*num_levels))
		y[y<0.5] <- 1-y[y<0.5]
		tty <- cut_interval(y, num_levels)
		xp <- quantile(abs(val), probs=cumsum(as.numeric(table(tty)/10000)))
		lower <- findInterval(abs(val),xp)
		lower[lower==0] <- 1
		lower[lower>num_levels] <- num_levels
		lvly <- round(seq(min(y),max(y),(max(y)-min(y))/(num_levels-1)),1)
		gval <- lvly[lower]
		l2 <- findInterval(gval, v)
		l2[l2==0] <- 1
		l2[l2>num_levels] <- num_levels
		#ggval <- 0.1*l2
		#upper <- lower +1
		#upper[upper>num_levels] <- num_levels
		#ggval <- v[upper]
	}
	return(l2)
	#return(lower)
}


gen_alphas_bin4 <- function(Sigma, theta, N, num_levels){
	val_levels <- get_levels(num_levels)
	Ngene <- nrow(Sigma)
	Sigma <- as.matrix(Sigma)
	if(theta > 0){
		realS <- as.matrix(nearPD(Sigma*theta)$mat)
	} else {
		realS <- as.matrix(nearPD(Sigma)$mat)*theta
	}
	
	#print(realS, Ngene, N)
	alphs <- mvrnorm(N, rep(0, Ngene), realS)
	if(num_levels%%2 == 0){
		alphs <- abs(alphs)
	}
	alphs2 <- data.frame(matrix(map_custom4(alphs, theta, num_levels), ncol=Ngene))
	
	#alphs_ntile <- round(alphs2/0.1)
	alphs_ntile <- alphs2
	print(table(alphs2))
	#alphs_ntile <- alphs2
	if(num_levels %% 2 ==1){
		xx <- apply(alphs_ntile, 1, function(x) sum(x!=floor(num_levels/2)+1)/length(x))
	}else{
		xx <- apply(alphs_ntile, 1, function(x) sum(x>1)/length(x))
	}
	alphx <- data.frame(cbind(apply(alphs_ntile, 1, function(x) paste(x[1:Ngene], collapse='.')),xx))
	names(alphx) <- c("key","pcnt")
	#alphx <- data.frame(apply(alphs_ntile, 1, function(x) paste(x, collapse='.')))
	#names(alphx) <- "key"
	dalphs <- ddply(alphx, colnames(alphx), nrow)
	return(list(dalphs,val_levels))
}

map_custom2 <- function(val, th, num_levels){
	if(th == 0){
		if(num_levels%%2 == 1){
			return(round(rep(0.5, length(val)),2))
		} else{
			return(round(rep(0.1, length(val)),2))
		}
	}
		
	al <- 0.5/th
	#x <- rnorm(10000, 0, th)
	y <- rbeta(10000, al, al)
	if(num_levels%%2 == 0){
		val <- abs(val)
		y[y<0.5] <- 1-y[y<0.5]
	}
	ty <- cut_interval(y, num_levels)
	xp <- quantile(val, probs=cumsum(as.numeric(table(ty)/10000)))
	lower <- findInterval(val,xp)
	upper <- lower+1
	upper[upper>num_levels] <- num_levels
	gval <- upper*0.1
	gval[gval==1] <- num_levels*0.1
	return(round(gval,1))
}

map_custom_beta <- function(val, th, num_levels, Ncase, Ncontrol){
	mu = Ncase / (Ncase+Ncontrol)
	if(th == 0){
		if(num_levels%%2 == 1){
			return(round(rep(0.5, length(val)),2))
		} else{
			return(round(rep(0.1, length(val)),2))
		}
	}
	lm <- Ncase/Ncontrol	
	al <- 0.5/th
	#x <- rnorm(10000, 0, th)
	y <- rbeta(10000, lm*al, al)
	if(num_levels%%2 == 0){
		val <- abs(val)
		y[y<mu] <- 2*mu-y[y<mu]
	}
	ty <- cut_interval(y, num_levels)
	xp <- quantile(val, probs=cumsum(as.numeric(table(ty)/10000)))
	lower <- findInterval(val,xp)
	upper <- lower+1
	upper[upper>num_levels] <- num_levels
	gval <- upper*0.1
	gval[gval==1] <- num_levels*0.1
	return(round(gval,1))
}

genDataPois <- function(Ncase, Ncontrol, Sigma, theta, low, cutoff, N, num_levels){
	case <- sample(c(rep(0,Ncontrol),rep(1,Ncase)))
	Ngene <- nrow(Sigma)
	realS <- theta*matrix(nearPD(Sigma)[[1]],nrow=Ngene)
	c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
	#print(dalphs)
	res <- rep(c(1:nrow(dalphs)), round(N*dalphs[,'V1']/N))
	#print(res)
	#reala <- unlist(dalphs[which.max(dalphs[,ncol(dalphs)]), -ncol(dalphs)])
	rkey <- as.character(unlist(dalphs[sample(res,1),-ncol(dalphs)]))
	#print(rkey)
	reala <- read_levels(rkey, alevels)
	#c(pcase, pcon) %<-% genProb(Ncase, Ncontrol, reala, realp)
	X <- list()
	Y <- list()
	Xcount <- list()
	Ycount <- list()
	for (i in 1:Ngene){
		y_sim <- rpois(n = Ncase+Ncontrol, lambda = exp(-2 + log(reala[i]/(1-reala[i])) * (case == 1))) 
		X[[i]] <- y_sim[case==1]
		Y[[i]] <- y_sim[case==0]
		h <- data.frame(table(y_sim, case))
		#print(h)
		Xcount[[i]] <- h[which(h[,"case"]==1 & h[,"y_sim"]!=0),'Freq']
		Ycount[[i]] <- h[which(h[,"case"]==0 & h[,"y_sim"]!=0),'Freq']
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	Xcount <- data.frame(plyr::ldply(Xcount, rbind))
	Ycount <- data.frame(plyr::ldply(Ycount, rbind))
	Xcount[is.na(Xcount)] <- 0
	Ycount[is.na(Ycount)] <- 0
	colnames(X) <- c(1:Ncase)
	colnames(Y) <- c(1:Ncontrol)
	colnames(Xcount) <- c(1:ncol(Xcount))
	colnames(Ycount) <- c(1:ncol(Ycount))
	
	return(list(X, Y, Xcount, Ycount, reala))
}


genDataPoisNull <- function(Ncase, Ncontrol, Sigma, theta, low, cutoff, N, num_levels){
	case <- sample(c(rep(0,Ncontrol),rep(1,Ncase)))
	Ngene <- nrow(Sigma)
	realS <- theta*matrix(nearPD(Sigma)[[1]],nrow=Ngene)
	#c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
	#print(dalphs)
	#res <- rep(c(1:nrow(dalphs)), round(N*dalphs[,'V1']/N))
	#print(res)
	#reala <- unlist(dalphs[which.max(dalphs[,ncol(dalphs)]), -ncol(dalphs)])
	#rkey <- as.character(unlist(dalphs[sample(res,1),-ncol(dalphs)]))
	#print(rkey)
	#reala <- read_levels(rkey, alevels)
	#c(pcase, pcon) %<-% genProb(Ncase, Ncontrol, reala, realp)
	reala <- rep(low, Ngene)
	X <- list()
	Y <- list()
	Xcount <- list()
	Ycount <- list()
	for (i in 1:Ngene){
		y_sim <- rpois(n = Ncase+Ncontrol, lambda = exp(-2 + log(reala[i]/(1-reala[i])) * (case == 1))) 
		X[[i]] <- y_sim[case==1]
		Y[[i]] <- y_sim[case==0]
		h <- data.frame(table(y_sim, case))
		#print(h)
		Xcount[[i]] <- h[which(h[,"case"]==1 & h[,"y_sim"]!=0),'Freq']
		Ycount[[i]] <- h[which(h[,"case"]==0 & h[,"y_sim"]!=0),'Freq']
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	Xcount <- data.frame(plyr::ldply(Xcount, rbind))
	Ycount <- data.frame(plyr::ldply(Ycount, rbind))
	Xcount[is.na(Xcount)] <- 0
	Ycount[is.na(Ycount)] <- 0
	colnames(X) <- c(1:Ncase)
	colnames(Y) <- c(1:Ncontrol)
	colnames(Xcount) <- c(1:ncol(Xcount))
	colnames(Ycount) <- c(1:ncol(Ycount))
	
	return(list(X, Y, Xcount, Ycount, reala))
}


genDataPoisPositive <- function(Ncase, Ncontrol, Sigma, theta, low, cutoff, N, num_levels){
	case <- sample(c(rep(0,Ncontrol),rep(1,Ncase)))
	Ngene <- nrow(Sigma)
	realS <- theta*matrix(nearPD(Sigma)[[1]],nrow=Ngene)
	c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
	#print(dalphs)
	res <- c(1:nrow(dalphs))
	print(res)
	rk <- sample(res,1)
	print(rk)
	while(rk == 1){
		rk <- sample(res,1)
	}
	#reala <- unlist(dalphs[which.max(dalphs[,ncol(dalphs)]), -ncol(dalphs)])
	rkey <- as.character(unlist(dalphs[rk,-ncol(dalphs)]))
	print(rkey)
	reala <- read_levels(rkey, alevels)
	#c(pcase, pcon) %<-% genProb(Ncase, Ncontrol, reala, realp)
	#lvl <- get_levels(low, theta, num_levels)
	#pos_lvl <- lvl[2:length(lvl)]
	#reala <- rep(low, Ngene)
	X <- list()
	Y <- list()
	Xcount <- list()
	Ycount <- list()
	for (i in 1:Ngene){
		y_sim <- rpois(n = Ncase+Ncontrol, lambda = exp(-2 + log(reala[i]/(1-reala[i])) * (case == 1))) 
		X[[i]] <- y_sim[case==1]
		Y[[i]] <- y_sim[case==0]
		h <- data.frame(table(y_sim, case))
		#print(h)
		Xcount[[i]] <- h[which(h[,"case"]==1 & h[,"y_sim"]!=0),'Freq']
		Ycount[[i]] <- h[which(h[,"case"]==0 & h[,"y_sim"]!=0),'Freq']
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	Xcount <- data.frame(plyr::ldply(Xcount, rbind))
	Ycount <- data.frame(plyr::ldply(Ycount, rbind))
	Xcount[is.na(Xcount)] <- 0
	Ycount[is.na(Ycount)] <- 0
	colnames(X) <- c(1:Ncase)
	colnames(Y) <- c(1:Ncontrol)
	colnames(Xcount) <- c(1:ncol(Xcount))
	colnames(Ycount) <- c(1:ncol(Ycount))
	
	return(list(X, Y, Xcount, Ycount, reala))
}

genDataPoisPercent2 <- function(Ncase, Ncontrol, Sigma, theta, low, cutoff, N, num_levels, pcnt){
	case <- sample(c(rep(0,Ncontrol),rep(1,Ncase)))
	Ngene <- nrow(Sigma)
	realS <- theta*matrix(nearPD(Sigma)[[1]],nrow=Ngene)
	c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
	print(dalphs)
	numPos <- round(pcnt*Ngene,0)
	numZero <- Ngene - numPos
	v <- c(sample(seq(1,num_levels-1,1),numPos, replace=T),rep(0,numZero))
	v <- sample(v)
	k <- paste(v,collapse='.')
	#rkey <- as.character(unlist(dalphs[rk,-ncol(dalphs)]))
	while(!(k %in% dalphs$key)){
		print("here")
		print(k)
		v <- c(sample(seq(1,num_levels-1,1),numPos, replace=T),rep(0,numZero))
		v <- sample(v)
		k <- paste(v,collapse='.')
	}
	#reala <- unlist(dalphs[which.max(dalphs[,ncol(dalphs)]), -ncol(dalphs)])
	print(k)
	reala <- read_levels(k, alevels)
	#c(pcase, pcon) %<-% genProb(Ncase, Ncontrol, reala, realp)
	#lvl <- get_levels(low, theta, num_levels)
	#pos_lvl <- lvl[2:length(lvl)]
	#reala <- rep(low, Ngene)
	X <- list()
	Y <- list()
	Xcount <- list()
	Ycount <- list()
	for (i in 1:Ngene){
		y_sim <- rpois(n = Ncase+Ncontrol, lambda = exp(-2 + log(reala[i]/(1-reala[i])) * (case == 1))) 
		X[[i]] <- y_sim[case==1]
		Y[[i]] <- y_sim[case==0]
		h <- data.frame(table(y_sim, case))
		#print(h)
		Xcount[[i]] <- h[which(h[,"case"]==1 & h[,"y_sim"]!=0),'Freq']
		Ycount[[i]] <- h[which(h[,"case"]==0 & h[,"y_sim"]!=0),'Freq']
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	Xcount <- data.frame(plyr::ldply(Xcount, rbind))
	Ycount <- data.frame(plyr::ldply(Ycount, rbind))
	Xcount[is.na(Xcount)] <- 0
	Ycount[is.na(Ycount)] <- 0
	colnames(X) <- c(1:Ncase)
	colnames(Y) <- c(1:Ncontrol)
	colnames(Xcount) <- c(1:ncol(Xcount))
	colnames(Ycount) <- c(1:ncol(Ycount))
	
	return(list(X, Y, Xcount, Ycount, reala))
}

genDataPoisPercent <- function(Ncase, Ncontrol, Sigma, theta, low, cutoff, N, num_levels, pcnt){
	case <- sample(c(rep(0,Ncontrol),rep(1,Ncase)))
	Ngene <- nrow(Sigma)
	realS <- theta*matrix(nearPD(Sigma)[[1]],nrow=Ngene)
	c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
	numPos <- round(pcnt*Ngene,0)
	numZero <- Ngene - numPos
	print(c(numPos, numZero))
	zk <- c()
	repeat {
		vc <- dalphs$key
		for(i in 1:length(vc)){
			v <- vc[i]
			vn <- as.numeric(str_split(v, '\\.')[[1]])
			if(sum(vn==0)==numZero){
				zk <- c(zk, v)
			}
		}
		if(length(zk)>0){
			break
		} else{
			c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
			print("regenerating dalphs")
		}
	}
	k <- sample(zk,1)
	print(k)
	reala <- read_levels(k, alevels)
	print(reala)
	X <- list()
	Y <- list()
	Xcount <- list()
	Ycount <- list()
	for (i in 1:Ngene){
		y_sim <- rpois(n = Ncase+Ncontrol, lambda = exp(-2 + log(reala[i]/(1-reala[i])) * (case == 1))) 
		X[[i]] <- y_sim[case==1]
		Y[[i]] <- y_sim[case==0]
		h <- data.frame(table(y_sim, case))
		#print(h)
		Xcount[[i]] <- h[which(h[,"case"]==1 & h[,"y_sim"]!=0),'Freq']
		Ycount[[i]] <- h[which(h[,"case"]==0 & h[,"y_sim"]!=0),'Freq']
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	Xcount <- data.frame(plyr::ldply(Xcount, rbind))
	Ycount <- data.frame(plyr::ldply(Ycount, rbind))
	Xcount[is.na(Xcount)] <- 0
	Ycount[is.na(Ycount)] <- 0
	colnames(X) <- c(1:Ncase)
	colnames(Y) <- c(1:Ncontrol)
	colnames(Xcount) <- c(1:ncol(Xcount))
	colnames(Ycount) <- c(1:ncol(Ycount))
	
	return(list(X, Y, Xcount, Ycount, reala))
}


genDataPoisPercent3 <- function(Ncase, Ncontrol, Sigma, theta, low, cutoff, N, num_levels, pcnt){
	case <- sample(c(rep(0,Ncontrol),rep(1,Ncase)))
	Ngene <- nrow(Sigma)
	realS <- theta*matrix(nearPD(Sigma)[[1]],nrow=Ngene)
	#c(dalphs,alevels) %<-% gen_alphas_bin(realS, theta, N, low, cutoff, num_levels)
	numPos <- round(pcnt*Ngene,0)
	numZero <- Ngene - numPos
	print(c(numPos, numZero))
	k <- paste0(sample(c(rep(0,numZero),sample(1:(num_levels-1),numPos, replace=T))),collapse=".")
	print(k)
	high <- theta/(1+theta)
	step <- (high-low)/(num_levels-1)
	alevels <- c(low, low+step, low+2*step, high)
	reala <- read_levels(k, alevels)
	print(reala)
	X <- list()
	Y <- list()
	Xcount <- list()
	Ycount <- list()
	for (i in 1:Ngene){
		y_sim <- rpois(n = Ncase+Ncontrol, lambda = exp(-2 + log(reala[i]/(1-reala[i])) * (case == 1))) 
		X[[i]] <- y_sim[case==1]
		Y[[i]] <- y_sim[case==0]
		h <- data.frame(table(y_sim, case))
		#print(h)
		Xcount[[i]] <- h[which(h[,"case"]==1 & h[,"y_sim"]!=0),'Freq']
		Ycount[[i]] <- h[which(h[,"case"]==0 & h[,"y_sim"]!=0),'Freq']
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	Xcount <- data.frame(plyr::ldply(Xcount, rbind))
	Ycount <- data.frame(plyr::ldply(Ycount, rbind))
	Xcount[is.na(Xcount)] <- 0
	Ycount[is.na(Ycount)] <- 0
	colnames(X) <- c(1:Ncase)
	colnames(Y) <- c(1:Ncontrol)
	colnames(Xcount) <- c(1:ncol(Xcount))
	colnames(Ycount) <- c(1:ncol(Ycount))
	
	return(list(X, Y, Xcount, Ycount, reala))
}

gen_alphas_bin <- function(Sigma, theta, N, num_levels){
	val_levels <- get_levels(num_levels)
	Ngene <- nrow(Sigma)
	Sigma <- as.matrix(Sigma)
	if(theta > 0){
		realS <- as.matrix(nearPD(Sigma*theta)$mat)
	} else {
		realS <- as.matrix(nearPD(Sigma)$mat)*theta
	}
	
	#print(realS, Ngene, N)
	alphs <- mvrnorm(N, rep(0, Ngene), realS)
	if(num_levels%%2 == 0){
		alphs <- abs(alphs)
	}
	alphs2 <- apply(alphs, 2, function(x) map_custom2(x, theta, num_levels))
	alphs_ntile <- round(alphs2/0.1)
	if(num_levels %% 2 ==1){
		xx <- apply(alphs_ntile, 1, function(x) sum(x!=floor(num_levels/2)+1)/length(x))
	}else{
		xx <- apply(alphs_ntile, 1, function(x) sum(x>1)/length(x))
	}
	alphx <- data.frame(cbind(apply(alphs_ntile, 1, function(x) paste(x[1:Ngene], collapse='.')),xx))
	names(alphx) <- c("key","pcnt")
	#alphx <- data.frame(apply(alphs_ntile, 1, function(x) paste(x, collapse='.')))
	#names(alphx) <- "key"
	dalphs <- ddply(alphx, colnames(alphx), nrow)
	return(list(dalphs,val_levels))
}


gen_alphas_bin2 <- function(Sigma, theta, N, num_levels){
	val_levels <- get_levels2(num_levels, theta)
	Ngene <- nrow(Sigma)
	Sigma <- as.matrix(Sigma)
	if(theta > 0){
		realS <- as.matrix(nearPD(Sigma*theta)$mat)
	} else {
		realS <- as.matrix(nearPD(Sigma)$mat)*theta
	}
	
	#print(realS, Ngene, N)
	alphs <- mvrnorm(N, rep(0, Ngene), realS)
	if(num_levels%%2 == 0){
		alphs <- abs(alphs)
	}
	alphs2 <- apply(alphs, 2, function(x) map_custom2(x, theta, num_levels))
	alphs_ntile <- round(alphs2/0.1)
	if(num_levels %% 2 ==1){
		xx <- apply(alphs_ntile, 1, function(x) sum(x!=floor(num_levels/2)+1)/length(x))
	}else{
		xx <- apply(alphs_ntile, 1, function(x) sum(x>1)/length(x))
	}
	alphx <- data.frame(cbind(apply(alphs_ntile, 1, function(x) paste(x[1:Ngene], collapse='.')),xx))
	names(alphx) <- c("key","pcnt")
	#alphx <- data.frame(apply(alphs_ntile, 1, function(x) paste(x, collapse='.')))
	#names(alphx) <- "key"
	dalphs <- ddply(alphx, colnames(alphx), nrow)
	return(list(dalphs,val_levels))
}

gen_alphas_bin_beta <- function(Sigma, theta, N, num_levels, Ncase, Ncontrol){
	val_levels <- get_levels_imb(num_levels, theta, Ncase, Ncontrol)
	Ngene <- nrow(Sigma)
	Sigma <- as.matrix(Sigma)
	if(theta > 0){
		realS <- as.matrix(nearPD(Sigma*theta)$mat)
	} else {
		realS <- as.matrix(nearPD(Sigma)$mat)*theta
	}
	
	#print(realS, Ngene, N)
	alphs <- mvrnorm(N, rep(0, Ngene), realS)
	if(num_levels%%2 == 0){
		alphs <- abs(alphs)
	}
	alphs2 <- apply(alphs, 2, function(x) map_custom2(x, theta, num_levels))
	alphs_ntile <- round(alphs2/0.1)
	if(num_levels %% 2 ==1){
		xx <- apply(alphs_ntile, 1, function(x) sum(x!=floor(num_levels/2)+1)/length(x))
	}else{
		xx <- apply(alphs_ntile, 1, function(x) sum(x>1)/length(x))
	}
	alphx <- data.frame(cbind(apply(alphs_ntile, 1, function(x) paste(x[1:Ngene], collapse='.')),xx))
	names(alphx) <- c("key","pcnt")
	#alphx <- data.frame(apply(alphs_ntile, 1, function(x) paste(x, collapse='.')))
	#names(alphx) <- "key"
	dalphs <- ddply(alphx, colnames(alphx), nrow)
	return(list(dalphs,val_levels))
}

make_alph_lookup_bin_beta <- function(Sigma, theta_try, N, rsense, num_levels, Ncase, Ncontrol){
	Ngene <- nrow(Sigma)
	dalphs <- data.frame()
	print("populating look up table")
	print(num_levels)
	for (theta in theta_try){
		c(x,val_levels) %<-% gen_alphas_bin_beta(Sigma, theta, N, num_levels, Ncase, Ncontrol)
		x[ncol(x)] <- x[ncol(x)]/N
		x <- x[x[ncol(x)]>=0.0001,]
		if(nrow(x)>0){
			x['theta'] <- theta
			dalphs <- rbind(dalphs,x)
		} else {
			print(theta)
		}
	}
	print("done!")
	return(dalphs)
}

make_alph_lookup_bin <- function(Sigma, theta_low, theta_high, N, rsense, num_levels=2){
	Ngene <- nrow(Sigma)
	dalphs <- data.frame()
	print("populating look up table")
	print(num_levels)
	for (theta in seq(theta_low,theta_high,10**(-rsense))){
		c(x,val_levels) %<-% gen_alphas_bin(Sigma, theta, N, num_levels)
		x[ncol(x)] <- x[ncol(x)]/N
		x <- x[x[ncol(x)]>=0.0001,]
		if(nrow(x)>0){
			x['theta'] <- theta
			dalphs <- rbind(dalphs,x)
		} else {
			print(theta)
		}
	}
	print("done!")
	return(dalphs)
}

make_alph_lookup_bin2 <- function(Sigma, theta_try, N, rsense, num_levels=2){
	Ngene <- nrow(Sigma)
	dalphs <- data.frame()
	print("populating look up table")
	print(num_levels)
	for (theta in theta_try){
		c(x,val_levels) %<-% gen_alphas_bin4(Sigma, theta, N, num_levels)
		x[ncol(x)] <- x[ncol(x)]/N
		x <- x[x[ncol(x)]>=0.0001,]
		if(nrow(x)>0){
			x['theta'] <- theta
			dalphs <- rbind(dalphs,x)
		} else {
			print(theta)
		}
	}
	print("done!")
	return(dalphs)
}

make_alph_lookup_bin3 <- function(Sigma, theta_try, N, rsense, num_levels=2){
	Ngene <- nrow(Sigma)
	dalphs <- data.frame()
	print("populating look up table")
	print(num_levels)
	for (theta in theta_try){
		c(x,val_levels) %<-% gen_alphas_bin2(Sigma, theta, N, num_levels)
		x[ncol(x)] <- x[ncol(x)]/N
		x <- x[x[ncol(x)]>=0.0001,]
		if(nrow(x)>0){
			x['theta'] <- theta
			dalphs <- rbind(dalphs,x)
		} else {
			print(theta)
		}
	}
	print("done!")
	return(dalphs)
}


read_levels <- function(key, alevels){
	vkey <- as.integer(unlist(strsplit(key,"\\.")))
	reala <- alevels[vkey]
	return(reala)
}

genDataNull <- function(Ncase, Ncontrol, Sigma, realp, theta, nvargene, N, num_levels){
	Ngene <- nrow(Sigma)
	zero_level <- 0
	if(num_levels%%2 == 1){
	zero_level <- floor(num_levels/2)+1
	} else{
	zero_level <- 1
	}
	k <- paste0(rep(as.character(zero_level),Ngene), collapse=".") 	
	print(k)
	alevels <- get_levels(num_levels)
	reala <- read_levels(k, alevels)
	print(reala)
	c(pcase, pcon) %<-% genProb2(Ncase, Ncontrol, reala, realp)
	X <- list()
	Y <- list()
	for (i in 1:Ngene){
		X[[i]] <- rbinom(Ncase, nvargene[i], pcase[i])
		Y[[i]] <- rbinom(Ncontrol, nvargene[i], pcon[i])
	}
	print("x generated")
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	max_count <- max(max(X), max(Y))
	while (rowSums(X==1) == 0 && rowSums(Y==1) == 0){
		X <- list()
		Y <- list()
		for (i in 1:Ngene){
			X[[i]] <- rbinom(Ncase, nvargene[i], pcase[i])
			Y[[i]] <- rbinom(Ncontrol, nvargene[i], pcon[i])
		}
		print("x generated")
		X <- data.frame(do.call(rbind, X))
		Y <- data.frame(do.call(rbind, Y))
		max_count <- max(max(X), max(Y))
	}
	gcase_count <- list()
	gcon_count <- list()
	for (j in 1:max_count){
		gcase_count[[j]] <- rowSums(X==j)
		gcon_count[[j]] <- rowSums(Y==j)
	}
	gcase_count <- data.frame(do.call(cbind, gcase_count))
	colnames(gcase_count) <- 1:max_count
	gcon_count <- data.frame(do.call(cbind, gcon_count))
	colnames(gcon_count) <- 1:max_count
	
	return(list(X, Y, gcase_count, gcon_count, reala))
}

genDataPcnt <- function(Ncase, Ncontrol, Sigma, realp, theta, nvargene, N, num_levels, pcnt){
	Ngene <- nrow(Sigma)
	c(dalphs,alevels) %<-% gen_alphas_bin(Sigma, theta, N, num_levels)
	while(!(pcnt %in% as.numeric(names(table(round(unlist(as.numeric(dalphs[,"pcnt"])),1)))))){
		pcnt = pcnt + 0.1
	}
	k <- sample(dalphs[round(as.numeric(dalphs[,"pcnt"]),1)==pcnt,]$key,1)
	print(k)
	alevels <- get_levels(num_levels)
	reala <- read_levels(k, alevels)
	print(reala)
	c(pcase, pcon) %<-% genProb2(Ncase, Ncontrol, reala, realp)
	X <- list()
	Y <- list()
	for (i in 1:Ngene){
		X[[i]] <- rbinom(Ncase, nvargene[i], pcase[i])
		Y[[i]] <- rbinom(Ncontrol, nvargene[i], pcon[i])
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	max_count <- max(max(X), max(Y))
	gcase_count <- list()
	gcon_count <- list()
	for (j in 1:max_count){
		gcase_count[[j]] <- rowSums(X==j)
		gcon_count[[j]] <- rowSums(Y==j)
	}
	gcase_count <- data.frame(do.call(cbind, gcase_count))
	colnames(gcase_count) <- 1:max_count
	gcon_count <- data.frame(do.call(cbind, gcon_count))
	colnames(gcon_count) <- 1:max_count
	
	return(list(X, Y, gcase_count, gcon_count, reala))
}

genData <- function(Ncase, Ncontrol, Sigma, realp, theta, nvargene, N=10000, num_levels=2){
	Ngene <- nrow(Sigma)
	c(dalphs,alevels) %<-% gen_alphas_bin(Sigma, theta, N, num_levels)
	print(dalphs)
	res <- rep(c(1:nrow(dalphs)), dalphs[,'V1'])
	rkey <- as.character(unlist(dalphs[sample(res,1),-ncol(dalphs)]))
	#print(rkey)
	reala <- read_levels(rkey, alevels)
	c(pcase, pcon) %<-% genProb2(Ncase, Ncontrol, reala, realp)
	X <- list()
	Y <- list()
	for (i in 1:Ngene){
		X[[i]] <- rbinom(Ncase, nvargene[i], pcase[i])
		Y[[i]] <- rbinom(Ncontrol, nvargene[i], pcon[i])
	}
	X <- data.frame(do.call(rbind, X))
	Y <- data.frame(do.call(rbind, Y))
	max_count <- max(max(X), max(Y))
	gcase_count <- list()
	gcon_count <- list()
	for (j in 1:max_count){
		gcase_count[[j]] <- rowSums(X==j)
		gcon_count[[j]] <- rowSums(Y==j)
	}
	gcase_count <- data.frame(do.call(cbind, gcase_count))
	colnames(gcase_count) <- 1:max_count
	gcon_count <- data.frame(do.call(cbind, gcon_count))
	colnames(gcon_count) <- 1:max_count
	
	return(list(X, Y, gcase_count, gcon_count, reala))
}

get_levels <- function(num_levels){
	#stopifnot(num_levels>=2 & num_levels<=5)
	val_levels <- c()
	if(num_levels%%2 == 1){
		p <- 1/num_levels
		ul <- floor(num_levels/2)
		ll <- -ul
		for(i in seq(ll, ul, 1)){
			val_levels <- c(val_levels, 0.5+i*p)
		}
	} else {
		p <- 1/(2*num_levels)
		for(i in seq(1,num_levels,1))
			val_levels <- c(val_levels, 0.5+(i-1)*p)
	}
	return(val_levels)
}

get_levels_corr_old <- function(num_levels, Ncase, Ncontrol){
	#stopifnot(num_levels>=2 & num_levels<=5)
	val_levels <- c()
	null_rate <- Ncase/(Ncase+Ncontrol)
	if(num_levels%%2 == 1){
		p <- 1/num_levels
		ul <- floor(num_levels/2)
		ll <- -ul
		for(i in seq(ll, ul, 1)){
			val_levels <- c(val_levels, max(0.0000001,min(null_rate+i*p,0.999999)))
		}
	} else {
		p <- 1/(2*num_levels)
		for(i in seq(1,num_levels,1))
			val_levels <- c(val_levels, max(0.0000001,min(null_rate+(i-1)*p,0.999999)))
	}
	return(val_levels)
}

get_levels2 <- function(num_levels, theta){
	#stopifnot(num_levels>=2 & num_levels<=5)
	val_levels <- c()
	if(num_levels%%2 == 1){
		p <- theta
		ul <- floor(num_levels/2)
		ll <- -ul
		for(i in seq(ll, ul, 1)){
			val_levels <- c(val_levels, min(0.5+i*p,1))
		}
	} else {
		p <- theta
		for(i in seq(1,num_levels,1))
			val_levels <- c(val_levels, min(0.5+(i-1)*p,1))
	}
	return(val_levels)
}

get_levels_imb <- function(num_levels, theta, Ncase, Ncontrol){
	#stopifnot(num_levels>=2 & num_levels<=5)
	val_levels <- c()
	null_rate <- Ncase/(Ncase+Ncontrol)
	if(theta == 0){
		return(rep(null_rate, num_levels))
	}
	x <- rbeta(10000, Ncase/Ncontrol*0.5/theta,0.5/theta)
	if(num_levels%%2 == 1){
		ll <- min(x)+1e-9
		ul <- max(x)-1e-9
		step_left <- (null_rate-ll)/floor(num_levels/2)
		step_right <- (ul-null_rate)/floor(num_levels/2)
		ulim <- floor(num_levels/2)
		llim <- -ulim
		for(i in seq(llim, 0, 1)){
			val_levels <- c(val_levels, null_rate+i*step_left)
		}
		for(i in seq(1,ulim,1)){
			val_levels <- c(val_levels, null_rate+i*step_right)
		}
	} else {
		x[x<null_rate] - 2*null_rate - x[x<null_rate]
		ul <- max(x)-1e-9
		step <- (ul-null_rate)/(num_levels-1)
		for(i in seq(1,num_levels,1))
			val_levels <- c(val_levels, null_rate+(i-1)*step)
	}
	return(val_levels)
}

get_levels_corr <- function(num_levels, Ncase, Ncontrol){
	#stopifnot(num_levels>=2 & num_levels<=5)
	val_levels <- c()
	null_rate <- Ncase/(Ncase+Ncontrol)
	ul <- min(2*null_rate,1)
	ll <- 2*null_rate - ul
	if(num_levels%%2 == 1){
		step <- (ul-ll)/(num_levels)
		ulim <- floor(num_levels/2)
		llim <- -ulim
		for(i in seq(llim, ulim, 1)){
			val_levels <- c(val_levels, null_rate+i*step)
		}
	} else {
		step <- (ul-ll)/(2*num_levels)
		for(i in seq(1,num_levels,1))
			val_levels <- c(val_levels, null_rate+(i-1)*step)
	}
	return(val_levels)
}

get_levels2 <- function(num_levels, theta){
	#stopifnot(num_levels>=2 & num_levels<=5)
	val_levels <- c()
	if(num_levels%%2 == 1){
		p <- theta
		ul <- floor(num_levels/2)
		ll <- -ul
		for(i in seq(ll, ul, 1)){
			val_levels <- c(val_levels, min(0.5+i*p,1))
		}
	} else {
		p <- theta
		for(i in seq(1,num_levels,1))
			val_levels <- c(val_levels, min(0.5+(i-1)*p,1))
	}
	return(val_levels)
}

calc_likelihood_prodterm_pois2 <- function(Xcount, Ycount, Ncase, Ncontrol, alpha){
	max_count <- ncol(Xcount)
	Ngene <- nrow(Xcount)
	prod_term = 1.0
	r <- Ncontrol/Ncase
	for (g in 1:Ngene){
		cond_l <- 1
		case_count <- sum(as.numeric(names(Xcount))*Xcount[g,])*r
		control_count <- sum(as.numeric(names(Ycount))*Ycount[g,])
		cond_l <- dbinom(round(case_count),round(case_count+control_count),alpha[[g]])
		prod_term = prod_term * cond_l
	}
	return(prod_term)
}

calc_likelihood_pois_bin2 <- function(Xcount, Ycount, Ncase, Ncontrol, Sigma, theta, N, num_levels, lt, rsense){
	L <- 0
	Ngene <- nrow(Sigma)
	#theta <- round(theta,rsense)
	alphas <- lt[abs(lt['theta']-theta)<10**-9,-ncol(lt)]
	#print(alphas)
	alevels <- get_levels2(num_levels, theta)
	for (i in 1:nrow(alphas)){
		akey <- as.character(unlist(alphas[i, "key"]))
		#print(akey)
		#print(alevels)
		alph <- read_levels(akey, alevels)
		#print(alph)
		prod <- calc_likelihood_prodterm_pois2(Xcount, Ycount, Ncase, Ncontrol, alph)
		L <- L + prod * alphas[i,ncol(alphas)]
	}
	return(L)
}

calc_likelihood_prodterm_pois_no_scaling <- function(Xcount, Ycount, Ncase, Ncontrol, alpha){
	max_count <- ncol(Xcount)
	Ngene <- nrow(Xcount)
	prod_term = 1.0
	r <- 1
	#r <- Ncontrol/Ncase
	for (g in 1:Ngene){
		cond_l <- 1
		case_count <- sum(as.numeric(names(Xcount))*Xcount[g,])*r
		control_count <- sum(as.numeric(names(Ycount))*Ycount[g,])
		if(case_count != 0 | control_count != 0){
			cond_l <- dbinom(round(case_count),round(case_count+control_count),alpha[[g]])
		}
		#print(c(rownames(Xcount)[g], case_count, control_count, cond_l))
		prod_term = prod_term * cond_l
	}
	return(prod_term)
}

calc_likelihood_prodterm_pois <- function(Xcount, Ycount, Ncase, Ncontrol, alpha){
	max_count <- ncol(Xcount)
	Ngene <- nrow(Xcount)
	prod_term = 1.0
	#r <- Ncontrol/Ncase
	for (g in 1:Ngene){
		cond_l <- 1
		#case_count <- sum(as.numeric(names(Xcount))*Xcount[g,])*r
		case_count <- sum(as.numeric(names(Xcount))*Xcount[g,])
		control_count <- sum(as.numeric(names(Ycount))*Ycount[g,])
		if(case_count != 0 | control_count != 0){
			cond_l <- dbinom(round(case_count),round(case_count+control_count),alpha[[g]])
		}
		#print(c(rownames(Xcount)[g], case_count, control_count, cond_l))
		prod_term = prod_term * cond_l
	}
	return(prod_term)
}

calc_likelihood_pois_bin <- function(Xcount, Ycount, Ncase, Ncontrol, Sigma, theta, N, num_levels, lt, rsense){
	L <- 0
	Ngene <- nrow(Sigma)
	theta <- round(theta,rsense)
	alphas <- lt[abs(lt['theta']-theta)<10**-9,-ncol(lt)]
	#print(alphas)
	#alevels <- get_levels(num_levels)
	alevels <- get_levels_corr(num_levels, Ncase, Ncontrol)
	for (i in 1:nrow(alphas)){
		akey <- as.character(unlist(alphas[i, "key"]))
		#print(akey)
		#print(alevels)
		alph <- read_levels(akey, alevels)
		#print(alph)
		prod <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol, alph)
		L <- L + prod * alphas[i,ncol(alphas)]
	}
	return(L)
}

calc_likelihood_pois_bin4 <- function(Xcount, Ycount, Ncase, Ncontrol, Sigma, theta, N, num_levels, lt, rsense){
	L <- 0
	Ngene <- nrow(Sigma)
	alphas <- lt[abs(lt['theta']-theta)<10**-9,-ncol(lt)]
	#print(alphas)
	alevels <- get_levels(num_levels)
	for (i in 1:nrow(alphas)){
		akey <- as.character(unlist(alphas[i, "key"]))
		#print(akey)
		#print(alevels)
		alph <- read_levels(akey, alevels)
		#print(alph)
		prod <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol, alph)
		L <- L + prod * alphas[i,ncol(alphas)]
	}
	return(L)
}

calc_likelihood_pois_bin_par <- function(Xcount, Ycount, Ncase, Ncontrol, Sigma, theta, N, num_levels, lt, rsense){
	L <- 0
	theta <- round(theta,rsense)
	alphas <- lt[abs(lt['theta']-theta)<10**-9,-ncol(lt)]
	alevels <- get_levels(num_levels)
	#print(alphas)
	Lprod <- foreach (i=1:nrow(alphas), .combine=c)%dopar%{
		akey <- as.character(unlist(alphas[i,"key"]))
		#print(akey)
		#print(alevels)
		alph <- read_levels(akey, alevels)
		#print(alph)
		prod <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol, alph) * alphas[i,ncol(alphas)]
		prod
	}
	L <- sum(Lprod)
	return(L)
}

find_mle_alpha_pois_bin <- function(Sigma, theta, N, Xcount, Ycount, Ncase, Ncontrol, num_levels){
	Ngene <- nrow(Sigma)
	c(dalphs, alevels) %<-% gen_alphas_bin(Sigma, theta, N, num_levels)
	#new
	alevels <- get_levels_corr(num_levels, Ncase, Ncontrol)
	#print(alevels)
	dalphs['V1'] <- dalphs['V1']/N
	zero_level <- 0
	if(num_levels%%2 == 1){
		zero_level <- round(num_levels/2)+1
	} else{
		zero_level <- 1
	}
	print(alevels)
	zerokey <- paste(as.character(rep(zero_level,nrow(Sigma))), collapse='.')
	alph <- read_levels(zerokey, alevels)
	print(alph)
	if(zerokey %in% dalphs[,"key"]){
		ind <- which(dalphs[,"key"]==zerokey)
		l0 <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol,alph )*dalphs[ind,ncol(dalphs)]
	} else{
		l0 <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol,alph )*1/N
	}
	l <- c()
	for (i in 1:nrow(dalphs)){
		akey <- as.character(unlist(dalphs[i, "key"]))
		alph <- read_levels(akey, alevels)
		lik <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol,alph )*dalphs[i,ncol(dalphs)]
		l <- c(l, lik)
	}
	maxl <- max(l)
	max_ind <- which.max(l)
	lr <- 2*(log(maxl) - log(l0))
	akey <- as.character(unlist(dalphs[max_ind,"key"]))
	#alevels <- get_levels(num_levels)
	#print(akey)
	#print(alevels)
	alph <- read_levels(akey, alevels)
	zl <- alevels[zero_level]
	return(list(alph, maxl, l0, lr, zl))
}

find_mle_alpha_pois_bin_beta <- function(Sigma, theta, N, Xcount, Ycount, Ncase, Ncontrol, num_levels, lt){
	Ngene <- nrow(Sigma)
	dalphs <- lt[abs(lt['theta']-theta)<10**-9,-ncol(lt)]
	alevels <- get_levels_corr(num_levels, Ncase, Ncontrol)
	zero_level <- 0
	if(num_levels%%2 == 1){
		zero_level <- round(num_levels/2)+1
	} else{
		zero_level <- 1
	}
	print(alevels)
	zerokey <- paste(as.character(rep(zero_level,nrow(Sigma))), collapse='.')
	alph <- read_levels(zerokey, alevels)
	print(alph)
	if(zerokey %in% dalphs[,"key"]){
		ind <- which(dalphs[,"key"]==zerokey)
		l0 <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol,alph )*dalphs[ind,ncol(dalphs)]
	} else{
		l0 <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol,alph )*1/N
	}
	l <- c()
	for (i in 1:nrow(dalphs)){
		akey <- as.character(unlist(dalphs[i, "key"]))
		alph <- read_levels(akey, alevels)
		lik <- calc_likelihood_prodterm_pois(Xcount, Ycount, Ncase, Ncontrol,alph )*dalphs[i,ncol(dalphs)]
		l <- c(l, lik)
	}
	maxl <- max(l)
	max_ind <- which.max(l)
	lr <- 2*(log(maxl) - log(l0))
	akey <- as.character(unlist(dalphs[max_ind,"key"]))
	#alevels <- get_levels(num_levels)
	#print(akey)
	#print(alevels)
	rel <- as.integer(unlist(strsplit(akey,"\\.")))-zero_level
	alph <- read_levels(akey, alevels)
	zl <- alevels[zero_level]
	return(list(alph, rel, maxl, l0, lr, zl))
}

# modified from the .gpdFittedPValue function from the waddR package
calc_gpd_pvalue <- function(val, null_distr) {
    
    # list of possible exceedance thresholds (decreasing)
	distr.ordered <- sort(null_distr, decreasing=TRUE)
	st = .025*length(distr.ordered)
    poss.exc.num <- seq(from=st, to=10, by=-10)
    bsn <- length(distr.ordered)
    r <- 1
    repeat {
        
        # set threshold for exceedance according to paper
        N.exc <- poss.exc.num[r]
        
        # compute set of N.exc exceedances
        exceedances <- distr.ordered[seq_len(N.exc)]
        
        # check whether the N.exc largest permutation values follow 
        # a GPD using an Anderson-Darling test
        gpd.ad.check <- gpdAd(exceedances)
        ad.pval <- gpd.ad.check$p.value
        
        r <- r + 1
        if (ad.pval > 0.05) {break}
    }
    
    # calculate exceedance threshold for so-obtained N.exc
    t.exc <- (distr.ordered[N.exc] 
              + distr.ordered[N.exc+1]) / 2
    
    # fit GPD distribution to the exceedances using maximum 
    # likelihood estimation
    gpd.fit <- gpdFit(  data=distr.ordered,
                        threshold=t.exc,
                        method="mle")
    
    # extract fitted parameters
    fit.scale <- as.numeric(gpd.fit$par.ests[1])
    fit.shape <- as.numeric(gpd.fit$par.ests[2])
    
    # compute GPD p-value (see paper)
    pvalue.gpd <- (N.exc / bsn) * (1 - pgpd(q=val-t.exc,
                                            loc=0,
                                            scale=fit.scale,
                                            shape=fit.shape))
   
    pvalue.gpd <- as.numeric(pvalue.gpd)
    pvalue.wass <- c("pvalue.gpd"=pvalue.gpd,
                     "ad.pval"=ad.pval,
                     "N.exc"=N.exc)
    return(pvalue.wass)
}


proDist<-function(currentP,sd){
  p<-rnorm(n=1,mean=currentP,sd=sd)
  if(p>0 & p<1) return(p)
  proDist(currentP,sd)
}
priorBeta<-function(p,a,b){
  prob<-dbeta(p,shape1=a,shape2 = b)
  # nllik<- -sum(dbinom(s,size=n,prob = parVec,log=TRUE))
  # cat("nllik= ",nllik,sep=" ",fill=T);cat(" ",sep=" ",fill=T)
  return(prob)
}

calc_realps <- function(max_count, nvargene, alpha, p, Ncase, Ncontrol){
	c(pcase, pcon) %<-% genProb(Ncase, Ncontrol, alpha, p)
	Ngene <- length(alpha)
	pcase_real <- list()
	pcon_real <- list()
	for (i in 1:Ngene){
		pcase_real[[i]] <- dbinom(1:max_count, nvargene[i], pcase[i])
		pcon_real[[i]] <- dbinom(1:max_count, nvargene[i], pcon[i])
	}
	pcase_real <- data.frame(do.call(rbind, pcase_real))
	colnames(pcase_real) <- 1:max_count
	pcon_real <- data.frame(do.call(rbind, pcon_real))
	colnames(pcon_real) <- 1:max_count
	return(list(pcase_real, pcon_real))
}

genProb2 <- function(Ncase, Ncontrol, reala, realp){
	n <- length(reala)
	pcase <- c()
	pcon <- c()
	r <- Ncase/Ncontrol
	for (i in 1:n){
		rat <- (1-reala[[i]])/reala[[i]]
		pcase[i] <- (Ncase+Ncontrol*r)*realp[i]/(Ncase+Ncontrol*rat*r)
		pcon[i] <- pcase[i]*rat
	}
	return(list(pcase, pcon))
}

calcPermPValue <- function(X, Y, Ncase, Ncontrol, LLR, theta_try, Sigma, N, num_levels, lt, rsense, num_cores){
	l <- cbind(X,Y)
	registerDoParallel(num_cores)
	stats <- list()
	
	z <- foreach (i=1:N, .combine=c)%dopar%{
		indX <- sample(1:(Ncase+Ncontrol), Ncase, replace=F)
		Xprime <- l[,indX]
		Yprime <- l[,-indX]
		max_count <- max(max(Xprime), max(Yprime))
		Xp_count <- list()
		Yp_count <- list()
		for (j in 1:max_count){
			Xp_count[[j]] <- rowSums(Xprime==j)
			Yp_count[[j]] <- rowSums(Yprime==j)
		}
		Xp_count <- data.frame(do.call(cbind, Xp_count))
		colnames(Xp_count) <- 1:max_count
		Yp_count <- data.frame(do.call(cbind, Yp_count))
		colnames(Yp_count) <- 1:max_count
		print(Xp_count)
		print(Yp_count)
		Lp <- c()
		for(th in theta_try){
			lp <- calc_likelihood_pois_bin(Xp_count, Yp_count, Ncase, Ncontrol, Sigma, th, N, num_levels, lt, rsense)
			Lp <- c(Lp,lp)
		}
		LLR_p <- 2*(log(max(Lp))-log(Lp[1]))
		LLR_p
	}
	perm_p <- (sum(z>=LLR)+1)/(length(z)+1)
	return(perm_p)
}

calcPermPValueEarlyStop <- function(X, Y, Ncase, Ncontrol, LLR, theta_try, Sigma, N, num_levels, lt, rsense, num_cores){
	l <- cbind(X,Y)
	registerDoParallel(num_cores)
	Nbatch <- N/100
	stats <- c()
	for(b in seq(1,Nbatch)){
		print(paste("batch",b))
		z <- foreach (i=1:100, .combine=c)%dopar%{
			indX <- sample(1:(Ncase+Ncontrol), Ncase, replace=F)
			Xprime <- l[,indX]
			Yprime <- l[,-indX]
			max_count <- max(max(Xprime), max(Yprime))
			Xp_count <- list()
			Yp_count <- list()
			for (j in 1:max_count){
				Xp_count[[j]] <- rowSums(Xprime==j)
				Yp_count[[j]] <- rowSums(Yprime==j)
			}
			Xp_count <- data.frame(do.call(cbind, Xp_count))
			colnames(Xp_count) <- 1:max_count
			Yp_count <- data.frame(do.call(cbind, Yp_count))
			colnames(Yp_count) <- 1:max_count
			print(Xp_count)
			print(Yp_count)
			Lp <- c()
			for(th in theta_try){
				lp <- calc_likelihood_pois_bin(Xp_count, Yp_count, Ncase, Ncontrol, Sigma, th, N, num_levels, lt, rsense)
				Lp <- c(Lp,lp)
			}
			LLR_p <- 2*(log(max(Lp))-log(Lp[1]))
			LLR_p
		}
		interm_p <- (sum(z>=LLR)+1)
		if(interm_p > 11){
			print(paste("early stopping at batch",b))
			return(1)
		} else {
			stats <- c(stats, z)
		}
	}
	perm_p <- (sum(stats>=LLR)+1)/(length(stats)+1)
	return(perm_p)
}

#source: https://rdrr.io/github/xihaoli/STAAR/src/R/CCT.R
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }

  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}


