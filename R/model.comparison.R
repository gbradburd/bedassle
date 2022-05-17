#' Run a BEDASSLE cross-validation analysis
#' 
#' \code{xValidation} runs a BEDASSLE cross-validation analysis
#' 
#' This function initiates a k-fold cross-validation analysis 
#' to determine the statistical support for the specified model.
#' 
#' @param partsFile A filename (in quotes, 
#'		with the full file path) to the data partitions object to be used 
#'		in the k-fold cross-validation procedure. This object can be created 
#'		using the \code{\link{makePartitions}} function in this package.
#' @param nReplicates An \code{integer} giving the number of cross-validation
#'		replicates to be run. This should be the same as the length of the list 
#'		specified in the \code{partsFile}.
#' @param nPartitions An \code{integer} giving the number of data folds 
#'		within each run. This should be the same as the length of each the 
#'		list specified for each replicate in the \code{partsFile}.
#' @param genDist A \code{matrix} of pairwise pi measured between 
#'		all pairs of samples.
#' @param geoDist A \code{matrix} of pairwise geographic distances 
#'		measured between all pairs of samples. A value of 
#'		\code{NULL} runs a model without geographic distance 
#'		as a predictor of genetic differentiation.
#' @param envDist A \code{matrix} of pairwise environmental distances 
#'		measured between all pairs of samples. If there 
#'		are multiple environmental distance measures, this 
#'		argument should be a \code{list} of distance matrices. 
#'		A value of \code{NULL} runs a model without geographic 
#'		distance as a predictor of genetic differentiation.
#' @param nLoci The total number of independent loci used to calculate 
#'		pairwise pi (\code{genDist}) in the dataset.
#' @param prefix A character vector giving the prefix to be attached 
#'		to all output files.
#' @param nIter An \code{integer} giving the number of iterations each MCMC 
#'		chain is run. Default is 2e3.  If the number of iterations 
#'		is greater than 500, the MCMC is thinned so that the number 
#'		of retained iterations is 500 (before burn-in).
#' @param parallel A \code{logical} value indicating whether or not to run the 
#'		different cross-validation replicates in parallel. Default is \code{FALSE}.
#'		For more details on how to set up runs in parallel, see the model 
#'		comparison vignette.
#' @param nNodes Number of nodes to run parallel analyses on. Default is 
#'		\code{NULL}. Ignored if \code{parallel} is \code{FALSE}. For more details 
#'		in how to set up runs in parallel, see the model comparison vignette. 
#' @param saveFiles A \code{logical} value indicating whether to automatically 
#'		save the output files from each cross-validation replicate. 
#'		Default is \code{FALSE}.
#' @param ... Further options to be passed to rstan::sampling (e.g., adapt_delta).
#'	
#' @return This function returns a matrix with \code{nReplicates} columns and 
#'		\code{nPartitions} columns giving the likelihood of each data partition 
#'		(averaged over the posterior distribution of the MCMC) in each replicate 
#'		analysis. The mean of these values gives an estimate of the predictive 
#'		accuracy of the specified model given the data provided. The mean and 
#'		standard error of the data partition likelihoods across replicates can be 
#' 		used for comparing models (e.g., with a t-test).
#'
#' In addition, this function saves a text file ("..._xval_results.txt"), 
#' containing the returned likelihoods for each replicate and data partition.
#'
#'@export
xValidation <- function(partsFile,nReplicates,nPartitions,genDist,geoDist=NULL,envDist=NULL,nLoci,prefix,nIter=2e3,parallel=FALSE,nNodes=1,saveFiles=FALSE,...){
	check.xval.call(args <- as.list(environment()))
	parts <- load.partitions.file(partsFile)
	check.partitions.arg(parts,nReplicates,nPartitions,genDist)
	message(announce.xval.procedure(nReplicates,nPartitions,genDist,geoDist,envDist))
	prespecified <- parallel.prespecify.check(args <- as.list(environment()))
	`%d%` <- parallelizing(args <- as.list(environment()))
	parts <- unlist(parts,recursive=FALSE)
	N <- nReplicates*nPartitions
    xvals <- foreach::foreach(n=1:N) %d% {
    				r <- (n-1) %/% nPartitions + 1;
    				p <- n - (r-1) * nPartitions;
				    message(sprintf("\nk-fold cross-validation: analyzing replicate %s/%s, partition %s/%s\n",
				    				r,nReplicates,p,nPartitions));
       			tryCatch(
	       			xvalBedassle(test = parts[[n]],
	       							 genDist = genDist,
	     							 geoDist = geoDist,
	       							 envDist = envDist,
	       							 nLoci = nLoci,
	       							 prefix = paste0(prefix,"rep",n),
	       							 nIter = nIter,
	       							 saveFiles,
	       							 ...),
	       			error=function(e){
	       				sprintf("replicate %s/%s, partition %s/%s experienced an error and could not complete\n",
	       					r,nReplicates,p,nPartitions)
	       				return(NA)
	       			})
       		 }
	xvals <- round(unlist(xvals),4)
	xvals <- matrix(xvals,nrow=nPartitions,ncol=nReplicates)
	colnames(xvals) <- paste0("rep_",1:nReplicates)
	write.xvals(xvals,nReplicates,nPartitions,prefix)
	tmp <- end.parallelization(prespecified)
    return(xvals)
}

utils::globalVariables("n")

#' Make the partitions used to run a BEDASSLE cross-validation analysis
#' 
#' \code{makePartitions} makes the partitions used to run a BEDASSLE cross-validation analysis
#' 
#' This function initiates a k-fold cross-validation analysis 
#' to determine the statistical support for the specified model.
#' 
#' @param prefix A \code{character} object giving the prefix to be attached 
#'		to the output data partitions object to be used in the k-fold 
#'		cross-validation procedure.
#' @param nReplicates An \code{integer} giving the number of cross-validation
#'		replicates to be run. The default value is 10.
#' @param nPartitions An \code{integer} giving the number of data folds 
#'		within each run. The default value is 5.
#' @param N An \code{integer} giving the number of samples in the dataset.
#'		This should have the same value as \code{nrow(geoDist)}.
#'
#' @details The number of samples \code{N} must be divisible by the specified 
#'		number of partitions (\code{nPartitions}).
#' @return This function generates and saves a nested list R object (".Robj" file) 
#'		that can be used to run a cross-validation analysis. This R object is a list 
#'		of length \code{nReplicates}, each element of which is a list of length 
#'		\code{nPartitions}. Each of those \code{nPartitions} element 
#'		is a vector of indices (between 1 and \code{N}) of the random subset 
#'		of individuals assigned to that cross-validation partition. 
#'		Each of the \code{N} individuals appears in exactly one partition per replicate.
#'
#'@export
makePartitions <- function(prefix,nReplicates=10,nPartitions=5,N){
	if(N%%nPartitions != 0){
		stop("\n\nThe number of samples (N) must be divisible by the number of partitions\n\n")
	}
	partitions <- lapply(1:nReplicates,
				function(n){
					reorderedSamples <- sample(1:N,N,replace=FALSE)
					nInPart <- N/nPartitions
					repParts <- lapply(1:nPartitions,
									function(i){
										reorderedSamples[((i-1)*nInPart+1):(i*nInPart)]
									})
					stats::setNames(repParts,paste0("partition",1:nPartitions))
	})
	partitions <- stats::setNames(partitions,paste0("rep",1:nReplicates))
	save(partitions,file=paste0(prefix,"_partitions.Robj"))
	message("\npartitions file saved\n")
	return(invisible("done"))
}

#' Compare the cross-validation output of different models
#' 
#' \code{compare.model.xvals} compares the outputs of different BEDASSLE cross-validation analyses
#' 
#' This function compares outputs of different BEDASSLE cross-validation analyses 
#' to determine the relative statistical support for each model and identify the best model.
#' 
#' @param xval.files A \code{character} vector giving the filenames (in quotes, 
#'		with the full file path) to the output files of a \code{\link{xValidation}} 
#'		analysis.
#' @param mod.order An \code{integer} vector giving the ordering of the model hypotheses 
#'		to be compared. Lower numbers indicate less complex or more parsimonious models.
#'		No ties.
#' @param mod.names An \code{character} vector giving the names associated with the models 
#'		run in the different cross-validation analyses.
#' @param mod.cols A vector containing the colors to be used in plotting the output of 
#'		the different cross-validation analyses. Default is black.
#'
#' @details This function compares the cross-validation predictive accuracy of different models 
#'		applied to the same data partitions (generated using \code{\link{makePartitions}}). 
#'		Going in the order of model complexity specified by \code{mod.order}, starting from the 
#'		least complex model, this function compares the performance of a model across all 
#'		cross-validation data partitions against that of the next most complex model. The model  
#'		that is determined to be "best" is the least complicated model whose performance is 
#'		statistically indistinguishable from adjacent model (ie, the model that is one rung 
#'		higher in complexity, as specified by \code{mod.order}.
#'
#' @return This function generates a figure comparing the cross-validation predictive accuracy 
#'		of different models used to analyze the same data partitions. The "best" model among the 
#'		set being compared is indicated with a golden arrow.  The significance of the pairwise 
#'		comparisons is indicated using letter groupings, as in a Tukey post-hoc test. The 
#'		function also returns the name of the "best" model. 
#'
#'@export
compare.model.xvals <- function(xval.files,mod.order,mod.names,mod.cols=NULL){
	check.compare.mod.xvals.call(args <- as.list(environment()))
	n.models <- length(xval.files)
	if(is.null(mod.cols)){
		mod.cols <- rep(1,n.models)
	}
	xval.results <- lapply(xval.files,function(n){read.xval.results(n)})
	nReplicates <- ncol(xval.results[[1]])
	xval.results <- lapply(xval.results,function(x){colMeans(x,na.rm=TRUE)})
	yLim <- range(unlist(xval.results),na.rm=TRUE)+c(0,diff(range(unlist(xval.results),na.rm=TRUE))/5)
	plot(0,xlim=c(0.5,length(xval.files)+0.5),
		   ylim=yLim,
		   xlab="models",
		   xaxt="n",
		   ylab="standardized predictive accuracy",
		   main="cross-validation model comparison",
		   type="n")
	graphics::axis(side=1,at=1:n.models,labels=mod.names)
	pairwise.sigs <- ttest.all.mod.xvals(xval.results,n.models)[mod.order,mod.order]
	group.labels <- get.sig.groups(pairwise.sigs,n.models)[mod.order]
	lab.ycoord <- graphics::par("usr")[4] - diff(range(unlist(xval.results),na.rm=TRUE))/10
	invisible(lapply(1:n.models,
		function(i){
			j <- mod.order[i]
			plot.mod.xval.summary(j,summarize.mod.xval(xval.results[[j]]),mod.cols[j])
			graphics::points(x=jitter(rep(j,nReplicates)),
					y=xval.results[[j]],
					pch=20,col=grDevices::adjustcolor(mod.cols[j],0.5))
			graphics::text(x=j,y=lab.ycoord,labels= group.labels[[j]],col=mod.cols[j])
		}))
	bestMod <- whichMod(xval.files,mod.order,mod.names)
	bestMod <- which(mod.names==bestMod)
	if(!is.null(bestMod)){
		graphics::arrows(x0=bestMod,
			  			 x1=bestMod,
						 y0=min(xval.results[[bestMod]])-2*diff(yLim)/10,
						 y1=min(xval.results[[bestMod]])-diff(yLim)/10,
						 col="black",lwd=2,length=0.03)
		graphics::arrows(x0=bestMod,
						 x1=bestMod,
						 y0=min(xval.results[[bestMod]])-2*diff(yLim)/10,
						 y1=min(xval.results[[bestMod]])-diff(yLim)/10,
						 col="goldenrod1",lwd=1.5,length=0.03)
	}
	return(mod.names[bestMod])
}

check.compare.mod.xvals.call <- function(args){
	check.xval.files.arg(args[["xval.files"]])
	check.mod.order.arg(args)
	check.mod.names.arg(args)
	check.mod.cols.arg(args)
	return(invisible("compare model xval args checked"))
}

check.xval.files.arg <- function(xval.files){
	if(length(xval.files) < 2){
		stop("\nyou must specify more than 1 xValidation output file to compare\n")
	}
	if(any(!is.character(xval.files))){
		stop("\nyou must specify a character vector for the \"xval.files\" argument\n")
	}
	if(any(!file.exists(xval.files))){
		stop("\nfunction is unable to find an \"xval.file\" specified\n")
	}
	return(invisible("xval.files checked"))
}

check.mod.order.arg <- function(args){
	nModels <- length(args[["xval.files"]])
	modOrder <- args[["mod.order"]]
	if(length(modOrder) != nModels){
		stop("\nthe length of \"mod.order\" does not match the number of files in \"xval.files\"\n")
	}
	if(!is.numeric(modOrder)){
		stop("\nthe \"mod.order\" argument must be of type \"numeric\"\n")	
	}
	if(length(modOrder) != length(unique(modOrder))){
		stop("\nthe \"mod.order\" argument cannot contain duplicate values\n")	
	}
	if(sum(abs(sort(modOrder)-1:(max(modOrder)))) != 0){
		stop("\nthe \"mod.order\" argument must contain all values between 1 and max(mod.order)\n")	
	}
	return(invisible("mod.order checked"))
}

check.mod.names.arg <- function(args){
	nModels <- length(args[["xval.files"]])
	modNames <- args[["mod.names"]]
	if(length(modNames) != nModels){
		stop("\nthe length of \"mod.names\" does not match the number of files in \"xval.files\"\n")
	}
	if(!is.character(modNames)){
		stop("\nthe \"mod.names\" must be of type \"character\"\n")	
	}
	if(length(modNames) != length(unique(modNames))){
		stop("\nthe \"mod.names\" argument cannot contain duplicate values\n")	
	}
	return(invisible("mod.names checked"))
}

check.mod.cols.arg <- function(args){
	if(!is.null(args[["mod.cols"]])){
		if(length(args[["mod.cols"]]) != length(args[["xval.files"]])){
			stop("\nyou must specify one plotting color for each model\n")
		}
	}
	return(invisible("mod.cols checked"))
}

read.xval.results <- function(xval.file){
	xval.result <- data.matrix(utils::read.table(xval.file,stringsAsFactors=FALSE,header=TRUE))
	if(any(!is.finite(xval.result))){
		xval.result[which(!is.finite(xval.result))] <- NA
	}
	return(xval.result)
}

standardize.xval.results <- function(xval.results){
	n.models <- length(xval.results)
	nReplicates <- ncol(xval.results[[1]])
	nPartitions <- nrow(xval.results[[1]])
	std.xval.array <- array(NA,dim=c(nPartitions,nReplicates,n.models))
	for(i in 1:n.models){
		std.xval.array[,,i] <- xval.results[[i]]
	}
	for(i in 1:nPartitions){
		for(j in 1:nReplicates){
			std.xval.array[i,j,] <- std.xval.array[i,j,] - max(std.xval.array[i,j,]) 
		}
	}
	std.xval.results <- lapply(1:n.models,function(i){std.xval.array[,,i]})
	return(std.xval.results)
}

summarize.mod.xval <- function(mod.xval.results){
	xval.mean <- mean(mod.xval.results,na.rm=TRUE)
	xval.std.err <- stats::sd(mod.xval.results)/sqrt(length(mod.xval.results))
	xval.CI <- xval.mean + c(-1.96,1.96) * xval.std.err
	mod.xval.summary <- list("mean" = xval.mean,
							 "CI" = xval.CI)
	return(mod.xval.summary)
}

plot.mod.xval.summary <- function(i,mod.xval.summary,mod.col=1){
	graphics::segments(x0=i,x1=i,y0=mod.xval.summary$CI[1],y1=mod.xval.summary$CI[2],lwd=1.8,col=mod.col)
	graphics::points(i,mod.xval.summary$mean,pch=19,cex=1.8,col=mod.col)
	return(invisible("plotted"))
}

passDist <- function(a,B=NULL){
	X <- NULL
	if(!is.null(B)){
		if(inherits(B,"matrix")){
			X <- B[a,a]
		} else if(inherits(B,"array")){
			X <- array(NA,dim=c(dim(B)[1],length(a),length(a)))
			for(i in 1:dim(B)[1]){
				X[i,,] <- B[i,a,a]
			}
		}
	}
	return(X)
}

xvalBedassle <- function(test,genDist,geoDist,envDist,nLoci,prefix,nIter,saveFiles,...){
	trainBlock <- makeXvalBlock(inSamples=c(1:nrow(genDist))[-test],genDist,geoDist,envDist,nLoci,silent=TRUE)
	testBlock <- makeXvalBlock(inSamples=test,genDist,geoDist,envDist,nLoci,silent=TRUE)
	stan.model <- pick.stan.model(trainBlock$geoDist,trainBlock$envDist)
	trainFit <- rstan::sampling(object = stanmodels[[stan.model]],
							 	 refresh = min(nIter/10,500),
							 	 data = trainBlock,
							 	 iter = nIter,
							 	 chains = 1,
							 	 init=getInits(trainBlock,stan.model),
							 	 thin = ifelse(nIter/500 > 1,nIter/500,1),
							 	 save_warmup = FALSE,
							 	 ...)
	if(saveFiles){
		trainBlock <- unstandardize.distances(trainBlock)
		save(trainBlock,file=paste(prefix,"data.block.Robj",sep="_"))
		save(trainFit,file=paste(prefix,"model.fit.Robj",sep="_"))
		bedassleResults <- get.bedassle.results(trainBlock,trainFit,nChains=1)
		write.bedassle.results(trainBlock,bedassleResults,prefix,nChains=1)
		make.all.bedassle.plots(paste0(prefix,"_posterior_chain",1,".txt"),paste(prefix,"data.block.Robj",sep="_"),prefix)
	}
	testFit <- fitToTest(trainFit,stan.model,testBlock)
	return(testFit)
}

makeXvalBlock <- function(inSamples,genDist,geoDist,envDist,nLoci,silent=TRUE){
	trainCov <- pwp2allelicCov(genDist[inSamples,inSamples])
	if(any(eigen(trainCov)$values < 0)){
		stop("\nthe allelic covariance calculated from the specified genetic distance matrix is not positive definite\n")
	}
	envDist <- make.envDist.list(envDist)
	sd.geoDist.list <- standardize.distances(geoDist)
	sd.envDist.list <- standardize.envDist(envDist)
	envDist.array <- make.envDist.array(sd.envDist.list)
	envDist.array <- passDist(a=inSamples,B=envDist.array)
	trainBlock <- list("N" = nrow(trainCov),
					   "L" = nLoci,
					   "obsCov" = trainCov,
					   "geoDist" = sd.geoDist.list$std.X[inSamples,inSamples],
					   "sd.geoDist" = sd.geoDist.list$stdev.X,
					   "envDist" = envDist.array,
					   "sd.envDist" = lapply(sd.envDist.list,"[[", "stdev.X"),
					   "nE" = ifelse(!is.null(envDist),length(envDist),0))
	trainBlock <- validate.data.block(trainBlock,silent=silent)
	return(trainBlock)
}

checkCov <- function(inits,dataBlock,modName){
	if(modName=="NULL"){
		parCov <- matrix(inits$alpha0,nrow=dataBlock$N,ncol=dataBlock$N) + diag(inits$nugget,dataBlock$N)
	}
	if(modName=="GEO"){
		parCov <- inits$alpha0 * exp(-(inits$alphaD*dataBlock$geoDist)^inits$alpha2) + diag(inits$nugget,dataBlock$N)
	}
	if(modName=="ENV"){
		eDist <- Reduce("+",
					lapply(1:dataBlock$nE,
						function(e){
							inits$alphaE[e]*dataBlock$envDist[e,,]
					}))
		parCov <- inits$alpha0 * exp(-eDist^inits$alpha2) + diag(inits$nugget,dataBlock$N)
	}
	if(modName=="GEOENV"){
		eDist <- Reduce("+",
					lapply(1:dataBlock$nE,
						function(e){
							inits$alphaE[e]*dataBlock$envDist[e,,]
					}))
		parCov <- inits$alpha0 * exp(-(eDist+inits$alphaD*dataBlock$geoDist)^inits$alpha2) + diag(inits$nugget,dataBlock$N)
	}
	posdef <- all(eigen(parCov)$values > 0)
	return(posdef)
}

getInits <- function(dataBlock,modName){
	posdef <- FALSE
	while(!posdef){
		if(modName=="NULL"){
			inits <- list("alpha0" = abs(stats::rnorm(1)),
						  "nugget"=abs(stats::rnorm(1,mean=0)))
		}
		if(modName=="GEO"){
			inits <- list("alpha0" = abs(stats::rnorm(1)),
						  "alphaD" = abs(stats::rnorm(1)),
						  "alpha2" = stats::runif(1,0,2),
						  "nugget"=abs(stats::rnorm(1,mean=0)))	
		}
		if(modName=="ENV"){
			inits <- list("alpha0" = abs(stats::rnorm(1)),
						  "alphaE" = as.array(abs(stats::rnorm(dataBlock$nE))),
						  "alpha2" = stats::runif(1,0,2),
						  "nugget"=abs(stats::rnorm(1,mean=0)))		
		}
		if(modName=="GEOENV"){
			inits <- list("alpha0" = abs(stats::rnorm(1)),
						  "alphaD" = abs(stats::rnorm(1)),
						  "alphaE" = as.array(abs(stats::rnorm(dataBlock$nE))),
						  "alpha2" = stats::runif(1,0,2),
						  "nugget"=abs(stats::rnorm(1,mean=0)))		
		}
		posdef <- checkCov(inits,dataBlock,modName)
	}
	inits <- list(inits)
	return(inits)
}

fitToTest <- function(trainFit,modName,testBlock){
	pars <- list("a0" = rstan::extract(trainFit,"alpha0",permute=FALSE),
				 "nugget" = rstan::extract(trainFit,"nugget",permute=FALSE))
	if(modName=="GEO"){
		pars[["aD"]] <- rstan::extract(trainFit,"alphaD",permute=FALSE)
		pars[["a2"]] <- rstan::extract(trainFit,"alpha2",permute=FALSE)
	}
	if(modName =="ENV"){
		pars[["aE"]] <- rstan::extract(trainFit,"alphaE",permute=FALSE)
		pars[["a2"]] <- rstan::extract(trainFit,"alpha2",permute=FALSE)
	}
	if(modName =="GEOENV"){
		pars[["aD"]] <- rstan::extract(trainFit,"alphaD",permute=FALSE)
		pars[["aE"]] <- rstan::extract(trainFit,"alphaE",permute=FALSE)
		pars[["a2"]] <- rstan::extract(trainFit,"alpha2",permute=FALSE)
	}
	ppc <- makePostParCov(nIter=length(pars$a0),
						  N=testBlock$N,
						  modName=modName,
						  pars=pars,
						  nE=testBlock$nE,
						  D=testBlock$geoDist,
						  E=testBlock$envDist)
	lnL <- mean(unlist(
			lapply(ppc,
				function(x){
					wishLnL(df=testBlock$L,parCov=x,obsCov=testBlock$obsCov)})))
	return(lnL)
}

makePostParCov <- function(nIter,N,modName,pars,nE,D=NULL,E=NULL){
	if(modName=="NULL"){
		pcPost <- lapply(1:nIter,
					function(i){
						pars$a0[i] + diag(pars$nugget[i],N)
					})
	} else if (modName=="GEO"){
		pcPost <- lapply(1:nIter,
					function(i){
						pars$a0[i] * exp(-(pars$aD[i]*D)^pars$a2[i]) + diag(pars$nugget[i],N)
					})
	} else if (modName=="ENV"){
		pcPost <- lapply(1:nIter,
					function(i){
						eDist <- Reduce("+",
									lapply(1:nE,
										function(e){
											pars$aE[i,,e]*E[e,,]
									}))
						pars$a0[i] * exp(-(eDist)^pars$a2[i]) + diag(pars$nugget[i],N)
					})	
	} else if (modName=="GEOENV"){
		pcPost <- lapply(1:nIter,
					function(i){
						eDist <- Reduce("+",
									lapply(1:nE,
										function(e){
											pars$aE[i,,e]*E[e,,]
									}))
						pars$a0[i] * exp(-(pars$aD[i]*D + eDist)^pars$a2[i]) + diag(pars$nugget[i],N)
					})
	}
	return(pcPost)
}

wishLnL <- function(df,parCov,obsCov){
	invParCov <- chol2inv(chol(parCov))
	logDet <- determinant(parCov)$modulus[[1]]
	lnL <- -0.5 * (sum(invParCov * obsCov) + df * logDet)
	return(lnL)
}

load.partitions.file <- function(partsFile){
	tmpenv <- environment()
	tmp <- load(partsFile,envir=tmpenv)
	rep.partitions <- lapply(tmp,get,envir=tmpenv)
	names(rep.partitions) <- tmp
	return(unlist(rep.partitions,recursive=FALSE))
}

announce.xval.procedure <- function(nReplicates,nPartitions,genDist,geoDist,envDist){
	N <- nrow(genDist)
	xval.proc <- sprintf("running %s %s of k-fold cross-validation on:\n\n\t",nReplicates, ifelse(nReplicates>1,"replicates","replicate"))
	if(is.null(geoDist) & is.null(envDist)){
		xval.proc <- paste0(xval.proc,"the \"NULL\" model (neither geographic or environmental distance)")
	}
	if(!is.null(geoDist) & is.null(envDist)){
		xval.proc <- paste0(xval.proc,"the \"GEO\" model (geographic, but not environmental distance)")
	}
	if(is.null(geoDist) & !is.null(envDist)){
		xval.proc <- paste0(xval.proc,"the \"ENV\" model (environmental, but not geographic distance)")
	}
	if(!is.null(geoDist) & !is.null(envDist)){
		xval.proc <- paste0(xval.proc,"the \"GEOENV\" model (geographic and environmental distance)")
	}
	if(!is.null(envDist)){
		nE <- ifelse(is.list(envDist),length(envDist),1)
		vars <- ifelse(nE == 1,"matrix","matrices")
		envAddendum <- sprintf("\n\twith %s environmental distance %s\n",nE,vars)
		xval.proc <- paste0(xval.proc,envAddendum)
	}
	xval.proc <- paste0(xval.proc,
					sprintf("\n\nwith %s data partitions, each comprised of %s individuals",
							nPartitions,N/nPartitions))
	return(xval.proc)
}

write.xvals <- function(xvals,nReplicates,nPartitions,prefix){
	utils::write.table(xvals,
						file=paste0(prefix,"_xval_results.txt"),
						row.names=FALSE,
						quote=FALSE)
	return(invisible("xval results written"))
}

parallel.prespecify.check <- function(args){
	prespecified <- FALSE
	if(args[["parallel"]] & foreach::getDoParRegistered()){
		prespecified <- TRUE
	}
	return(prespecified)
}

end.parallelization <- function(prespecified){
	if(!prespecified){
		doParallel::stopImplicitCluster()
		message("\nParallel workers terminated\n\n")
	}
	return(invisible("if not prespecified, parallelization ended"))
}

parallelizing <- function(args){
	if(args[["parallel"]]){
		if(!foreach::getDoParRegistered()){
			if(is.null(args[["nNodes"]])){
				nNodes <- parallel::detectCores()-1
			} else {
				nNodes <- args[["nNodes"]]
			}
			cl <- parallel::makeCluster(nNodes)
			doParallel::registerDoParallel(cl)
			message("\nRegistered doParallel with ",nNodes," workers\n")
		} else {
			message("\nUsing ",foreach::getDoParName()," with ", foreach::getDoParWorkers(), " workers")
		}
		d <- foreach::`%dopar%`
	} else {
		message("\nRunning sequentially with a single worker\n")
		d <- foreach::`%do%`
	}
	return(d)
}

get.par.cov <- function(model.fit,chain.no,N){
	par.cov <- array(NA,dim=c(model.fit@sim$n_save[chain.no],N,N))
	for(i in 1:N){
		for(j in 1:N){
			my.par <- sprintf("parCov[%s,%s]",i,j)
			par.cov[,i,j] <- rstan::extract(model.fit,pars=my.par,inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
		}
	}
	return(par.cov)
}

post.process.par.cov <- function(parCov){
	pp.cov.list <- lapply(1:dim(parCov)[1],
							function(i){
								list("inv" = chol2inv(chol(parCov[i,,])),
									 "log.det" = determinant(parCov[i,,])$modulus[[1]])
							})
	return(pp.cov.list)
}

check.xval.call <- function(args){
	check.partitions.file(args)
	check.xval.misc.args(args)
	check.xval.dist.args(args)
	check.parallel.args(args)
	check.for.output.files(args)
	return(invisible("args checked"))
}

check.xval.misc.args <- function(args){
	if(!is.whole.number(args[["nReplicates"]])){
		stop("\nyou must specify an integer value for the \"nReplicates\" argument\n")	
	}
	if(!is.whole.number(args[["nPartitions"]])){
		stop("\nyou must specify an integer value for the \"nPartitions\" argument\n")	
	}
	if(!is.whole.number(args[["nLoci"]])){
		stop("\nyou must specify an integer value for the \"nLoci\" argument\n")	
	}
	if(!is.character(args[["prefix"]])){
		stop("\nyou must specify a character value for the \"prefix\" argument\n")		
	}
	if(!is.whole.number(args[["nIter"]])){
		stop("\nyou must specify an integer value for the \"nIter\" argument\n")	
	}
	if(!is.logical(args[["saveFiles"]])){
		stop("\nyou must specify a logical value (TRUE or FALSE) for the \"saveFiles\" argument\n")
	}
	return(invisible("misc args checked"))
}

check.partitions.file <- function(args){
	if(!is.character(args[["partsFile"]])){
		stop("\nyou must specify a character vector for the \"partsFile\" argument\n")
	}
	if(!file.exists(args[["partsFile"]])){
		stop("\nfunction is unable to find the \"partsFile\" specified\n")
	}
	return(invisible("partition files checked"))
}

check.xval.dist.args <- function(args){
	check.geoDist.arg(args)
	check.envDist.arg(args)
	return(invisible("dist args checked"))
}

check.parallel.args <- function(args){
	if(args[["parallel"]] & args[["nNodes"]]==1){
		stop("\nyou have specified the \"parallel\" option with \"nNodes\" set to 1.\n\n")
	}
	if(!args[["parallel"]] & args[["nNodes"]] > 1){
		stop("\nyou have are running with \"parallel\" set to FALSE but with \"nNodes\" greater than 1.\n\n")
	}
	if(!args[["parallel"]] & foreach::getDoParWorkers() > 1){
		stop("\nyou are running with more than 1 worker but you have set the \"parallel\" option to FALSE\n\n")
	}
	return(invisible("parallel args checked"))
}

check.for.output.files <- function(args){
	if(file.exists(paste0(args[["prefix"]],"_xval_results.txt"))){
		stop("\noutput files will be overwritten if you proceed with this analysis\n\n")
	}
	return(invisible("files checked for"))
}

check.partitions.arg <- function(parts,nReplicates,nPartitions,genDist){
	N <- nrow(genDist)
	n <- N/nPartitions
	if(length(parts) != nReplicates){
		stop("\nthe number of replicates in the \"partitions\" argument does not match the number of replicates specified\n")
	}
	if(!inherits(parts,"list")){
		stop("\nthe partitions object must be of class \"list\"\\n")
	}
	lapply(1:length(nPartitions),
		function(i){
			check.rep(parts[[i]],sprintf("rep%s",i),nPartitions,n,genDist)
		})
	return(invisible("rep partitions arg checked"))
}

check.rep <- function(rep,rep.name,nPartitions,n,genDist){
	if(length(rep) != nPartitions){
		stop(sprintf("\nthe number of partitions in %s of the \"partsFile\" argument does not match the number of partitions specified\n",rep.name))
	}
	if(!inherits(rep,"list")){
		stop("\nthe data partitions object within each cross-validation replicate must be of class \"list\"\n")
	}
	lapply(1:length(rep),
		function(k){
			check.partition(rep[[k]],rep.name,sprintf("partition_%s",k),n,genDist)
		})
	return(invisible("rep checked"))
}

check.partition <- function(partition,rep.name,partition.name,n,genDist){
	if(length(partition)!=n){
		stop(
			sprintf("\nyou have specified a data partition (%s in %s) with the incorrect dimensions\n",
					partition.name,rep.name))
	}
	allCov <- pwp2allelicCov(genDist[partition,partition])
	if(any(eigen(allCov)$values <= 0)){
		stop(sprintf("\nthe allelic covariance calculated from the specified genetic distance matrix by %s in %s is not positive definite\n",partition.name,rep.name))
	}
	return(invisible("partition checked"))
}

get.nrow <- function(geoDist,envDist){
	n.row <- NULL
	if(!is.null(geoDist)){
		n.row <- nrow(geoDist)
	}
	if(!is.null(envDist)){
		if(inherits(envDist,"matrix")){
			n.row <- nrow(envDist)
		}
		if(inherits(envDist,"list")){
			n.row <- nrow(envDist[[1]])
		}
	}
	return(n.row)
}

whichMod <- function(xval.files,mod.order,mod.names){
	n.models <- length(xval.files)
	xval.results <- lapply(xval.files,function(n){read.xval.results(n)})
	xval.results <- xval.results
	nReplicates <- ncol(xval.results[[1]])
	xval.results <- lapply(xval.results,
						function(x){
							if(any(!is.finite(x))){
								z <- x
								z[which(!is.finite(x))] <- NA
								colMeans(z,na.rm=TRUE)
							} else {
								colMeans(x)
							}
					})
	modComps <- lapply(1:n.models,function(i){
					sapply(c(1:n.models)[-i],function(j){
						if(mod.order[i] < mod.order[j]){
							c(mod.names[i],mod.names[j])[head2headMod(xval.results[[i]],xval.results[[j]])]
						} else {
							c(mod.names[j],mod.names[i])[head2headMod(xval.results[[j]],xval.results[[i]])]
						}
					})
				})
	z <- 0
	halt <- FALSE
	while(!halt){
		z <- z + 1
		halt <- !any(modComps[[z]] != mod.names[z],na.rm=TRUE) | z == n.models
	}
	return(mod.names[z])
}

head2headMod <- function(x1,x2){
	p <- stats::t.test(x=x1,y=x2,paired=TRUE)$p.value
	meanDiff <- mean(x1) - mean(x2)
	if (p > 0.05){
		betterMod <- 1
	} else if(p < 0.05 & meanDiff > 0){
		betterMod <- 1
	} else if (p < 0.05 & meanDiff < 0){
		betterMod <- 2
	}
	return(betterMod) 
}

ttest.all.mod.xvals <- function(xval.results,n.models){
	pairwise.sigs <- matrix(NA,nrow=n.models,ncol=n.models)
	for(i in 1:n.models){
		for(j in 1:n.models){
			if(i != j){
				pairwise.sigs[i,j] <- pairwise.mod.xval.ttest(xval.results[[i]],xval.results[[j]])
			}
		}
	}
	return(pairwise.sigs)
}

get.sig.groups <- function(pairwise.sigs,n.models){
	groupLetters <- c("A",rep("",n.models-1))
	groupLetters[which(pairwise.sigs[1,] > 0.05)] <- "A"
	for(i in 2:n.models){
		if(groupLetters[i]==""){
			groupLetters[i] <- paste0(groupLetters[i],LETTERS[i],collapse="")
			if(any(pairwise.sigs[i,] > 0.05,na.rm=TRUE)){
				inGroup <- which(pairwise.sigs[i,] > 0.05)
				for(j in inGroup){
					groupLetters[j] <- paste0(groupLetters[j],LETTERS[i],collapse="")
				}
			}
		}
	}
	return(groupLetters)
}

pairwise.mod.xval.ttest <- function(xvals1,xvals2){
	pval <- stats::t.test(x=xvals2,y=xvals1,paired=TRUE,alternative="two.sided")$p.value
	return(pval)
}