#' Run a BEDASSLE cross-validation analysis
#' 
#' \code{x.validation} runs a BEDASSLE cross-validation analysis
#' 
#' This function initiates a k-fold cross-validation analysis 
#' to determine the statistical support for the specified model.
#' 
#' @param partitions.file A filename (in quotes, 
#'		with the full file path) to the data partitions object to be used 
#'		in the k-fold cross-validation procedure.
#' @param n.replicates An \code{integer} giving the number of cross-validation
#'		replicates to be run. This should be the same as the length of the list 
#'		specified in the \code{partitions.file}.
#' @param n.partitions An \code{integer} giving the number of data folds 
#'		within each run. This should be the same as the length of each the 
#'		list specified for each replicate in the \code{partitions.file}.
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
#' @param nLoci The total number of loci across all data partitions		
#'		(the sum across data partitions of the number of loci used in 
#'		calculating pairwise pi specified in the \code{partitions.file} argument.
#' @param prefix A character \code{vector} giving the prefix to be attached 
#'		to all output files.
#' @param n.iter An \code{integer} giving the number of iterations each MCMC 
#'		chain is run. Default is 2e3.  If the number of iterations 
#'		is greater than 500, the MCMC is thinned so that the number 
#'		of retained iterations is 500 (before burn-in).
#' @param n.chains An integer indicating the number of MCMC chains to be run 
#'		in the analysis. Default is 2.
#' @param parallel A \code{logical} value indicating whether or not to run the 
#'		different cross-validation replicates in parallel. Default is \code{FALSE}.
#'		For more details on how to set up runs in parallel, see the model 
#'		comparison vignette.
#' @param n.nodes Number of nodes to run parallel analyses on. Default is 
#'		\code{NULL}. Ignored if \code{parallel} is \code{FALSE}. For more details 
#'		in how to set up runs in parallel, see the model comparison vignette. 
#' @param save.files A \code{logical} value indicating whether to automatically 
#'		save the output files from each cross-validation replicate. 
#'		Default is \code{FALSE}.
#' @param ... Further options to be passed to rstan::sampling (e.g., adapt_delta).
#'	
#' @return This function returns a matrix with \code{n.replicates} columns and 
#'		\code{n.partitions} columns giving the likelihood of each data partition 
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
x.validation <- function(partitions.file,n.replicates,n.partitions,geoDist=NULL,envDist=NULL,nLoci,prefix,n.iter=2e3,n.chains=2,parallel=FALSE,n.nodes=1,save.files=FALSE,...){
	check.xval.call(args <- as.list(environment()))
	rep.partitions <- load.partitions.file(partitions.file)
	check.rep.partitions.arg(rep.partitions,n.replicates,n.partitions,geoDist,envDist)
	message(announce.xval.procedure(n.replicates,n.partitions,nLoci,geoDist,envDist))
	prespecified <- parallel.prespecify.check(args <- as.list(environment()))
	`%d%` <- parallelizing(args <- as.list(environment()))
	n <- 1
    x.val <- foreach::foreach(n=1:n.replicates) %d% {
				    message(sprintf("\nk-fold cross-validation: analyzing replicate %s/%s\n",n,n.replicates));
        				xval.bedassle.rep(rep = rep.partitions[[n]],
        								  n.partitions = n.partitions,
        								  geoDist = geoDist,
        								  envDist = envDist,
        								  nLoci = nLoci,
        								  prefix = paste0(prefix,"rep",n),
        								  n.chains = n.chains,
        								  n.iter = n.iter,
        								  save.files,
        								  ...)
    			 }
	x.val <- lapply(x.val,round,4)
	x.val <- Reduce("cbind",x.val)
	colnames(x.val) <- paste0("rep_",1:n.replicates)
	write.xvals(x.val,n.replicates,n.partitions,prefix)
	tmp <- end.parallelization(prespecified)
    return(unlist(x.val))
}

#' Compare cross-validated BEDASSLE models
#' 
#' \code{compare.model.xvals} compares cross-validated BEDASSLE models
#' 
#' This function visually and statistically compares BEDASSLE models that 
#' have been evaluated using cross-validation.
#' 
#' @param xval.files A \code{vector} of filenames (each in quotes, 
#'		with the full file path), each element of which points to the 
#'		cross-validation results file output by a call to 
#'		\code{\link{x.validation}} for one of the models the user wants 
#'		to compare.
#' @param n.predictors A \code{vector} of \code{integer} values giving 
#'		the number of predictors included in each of the BEDASSLE models 
#'		tested with cross-validation. This argument should be the same as 
#'		the length of the \code{vector} of filenames specified in 
#'		\code{xval.files}.
#' @param mod.cols An \code{vector} of colors to be used in plotting 
#'		the results of the different BEDASSLE models. This argument should be 
#'		the same length as the \code{vector} of filenames specified in 
#'		\code{xval.files}. If \code{NULL}, all output will be plotted in 
#'		black.
#'
#' @return This function creates a plot showing the standardized predictive 
#'		accuracies of the different models evaluated using cross-validation. 
#'		A higher predictive accuracy means the model is able to describe the 
#'		data better. The plot shows mean standardized predictive accuracy for 
#'		each replicate, the mean across replicates, and the 95% confidence 
#'		letter over a pair of models indicates that the predictive accuracies 
#'		interval. A shared of those models are not significantly different 
#'		using a paired, two-tailed t-test with a significance level of 0.05.
#'
#'		The predictive accuracies are standardized by, for each partition in 
#'		each replicate, subtracting the highest predictive accuracy from the  
#'		predictive accuracies of the other models. Therefore, a predictive 
#'		accuracy of 0 is the best score.
#'
#'		The function returns a matrix giving the significance of the difference 
#'		between the predictive accuracies of the different models evaluted using 
#'		k-fold cross-validation in the function \code{\link{x.validation}}
#'		analysis.
#'
#' Note that \code{n.predictors} is the number of predictor variables in the 
#' relevant model. The null model, with neither geographic nor environmental 
#' distance, has 0 predictors, a model with just geographic distance has 1 
#' predictor, as does a model with just a single environmental distance 
#' variable, and a model with both geographic distance and an environmental 
#' distance predictor has 2 predictors.
#' 
#'@export

compare.model.xvals <- function(xval.files,n.predictors,mod.cols=NULL){
	check.compare.mod.xvals.call(args <- as.list(environment()))
	n.models <- length(xval.files)
	if(is.null(mod.cols)){
		mod.cols <- rep(1,n.models)
	}
	xval.results <- lapply(xval.files,function(n){read.xval.results(n)})
	#xval.results <- standardize.xval.results(xval.results)
	n.replicates <- ncol(xval.results[[1]])
	xval.results <- lapply(xval.results,function(x){colMeans(x)})
	plot(0,xlim=c(0.5,length(xval.files)+0.5),
		   ylim=range(unlist(xval.results))+c(0,diff(range(unlist(xval.results)))/5),
		   xlab="models (# predictors)",
		   xaxt="n",
		   ylab="standardized predictive accuracy",
		   main="cross-validation model comparison",
		   type="n")
	graphics::axis(side=1,at=1:n.models,labels=sapply(1:n.models,function(i){paste0("mod",i," (",n.predictors[i],")")}))
	pairwise.sigs <- ttest.all.mod.xvals(xval.results,n.predictors,n.models)
	groups <- get.sig.groups(pairwise.sigs,n.models)
	group.labels <- groups2labels(groups,n.models)
	lab.ycoord <- graphics::par("usr")[4] - diff(range(unlist(xval.results)))/10
	invisible(lapply(1:n.models,
		function(i){
			plot.mod.xval.summary(i,summarize.mod.xval(xval.results[[i]]),mod.cols[i])
			graphics::points(x=jitter(rep(i,n.replicates)),
					y=xval.results[[i]],
					pch=20,col=grDevices::adjustcolor(mod.cols[i],0.5))
			graphics::text(x=i,y=lab.ycoord,labels= group.labels[[i]],col=mod.cols[i])
		}))
	return(pairwise.sigs)
}

ttest.all.mod.xvals <- function(xval.results,n.predictors,n.models){
	pairwise.sigs <- matrix(NA,nrow=n.models,ncol=n.models)
	for(i in 1:n.models){
		for(j in 1:n.models){
			if(i != j){
				pairwise.sigs[i,j] <- pairwise.mod.xval.ttest(xval.results[[i]],xval.results[[j]],n.predictors[i],n.predictors[j])
			}
		}
	}
	return(pairwise.sigs)
}

groups2labels <- function(groups,n.models){
	mod.letts <- LETTERS[1:length(groups)]
	mod.labs <- character(n.models)
	for(i in 1:length(groups)){
		for(j in 1:length(groups[[i]])){
			mod.labs[groups[[i]][j]] <- paste0(mod.labs[groups[[i]][j]],mod.letts[i],collapse="")
		}
	}
	return(mod.labs)
}

get.sig.groups <- function(pairwise.sigs,n.models){
	groups <- list("group1"=1)
	groups[[1]] <- c(groups[[1]],which(pairwise.sigs[1,] > 0.05))
	n.groups <- 1
	for(i in 2:n.models){
		if(any(pairwise.sigs[i,] > 0.05,na.rm=TRUE)){
			groups <- c(groups, i)
			n.groups <- n.groups + 1
			names(groups)[length(groups)] <- sprintf("group%s",n.groups)
			groups[[n.groups]] <- c(groups[[n.groups]],which(pairwise.sigs[i,] > 0.05))
		} else {
			groups <- c(groups, i)
			n.groups <- n.groups + 1
			names(groups)[length(groups)] <- sprintf("group%s",n.groups)
		}
	}
	groups <- collapse.ident.groups(groups)
	return(groups)
}

collapse.ident.groups <- function(groups){
	n.groups <- length(groups)
	ident <- matrix(0,n.groups,n.groups)
	for(i in 1:(n.groups-1)){
		for(j in (i+1):n.groups){
			if(setequal(sort(groups[[i]]),sort(groups[[j]]))){
				ident[i,j] <- 1
			}
		}
	}
	if(any(ident==1)){
		groups[unique(which(ident==1,arr.ind=TRUE)[,2])] <- NULL
	}
	return(groups)
}

pairwise.mod.xval.ttest <- function(xvals1,xvals2,np1,np2){
	pval <- stats::t.test(x=xvals2,y=xvals1,paired=TRUE,alternative="two.sided")$p.value	
	return(pval)
}

check.compare.mod.xvals.call <- function(args){
	check.xval.files.arg(args[["xval.files"]])
	check.n.predictors.arg(args)
	check.xval.cols.arg(args)
	return(invisible("compare model xval args checked"))
}

check.xval.files.arg <- function(xval.files){
	if(length(xval.files) < 2){
		stop("\nyou must specify more than 1 x.validation output file to compare\n")
	}
	if(any(!is.character(xval.files))){
		stop("\nyou must specify a character vector for the \"xval.files\" argument\n")
	}
	if(any(!file.exists(xval.files))){
		stop("\nfunction is unable to find an \"xval.file\" specified\n")
	}
	return(invisible("xval.files checked"))
}

check.n.predictors.arg <- function(args){
	if(length(args[["n.predictors"]]) != length(args[["xval.files"]])){
		stop("\nyou must specify a number of predictors for each model\n")
	}
	if(any(!is.numeric(args[["n.predictors"]]))){
		stop("\nthe values of the \"n.predictors\" argument must all be numeric\n")
	}
	if(any(args[["n.predictors"]] < 0)){
		stop("\nthe values of the \"n.predictors\" argument must all be greater than or equal to 0\n")
	}
	if(any(!sapply(args[["n.predictors"]],is.whole.number))){
		stop("\nthe values of the \"n.predictors\" argument must all be integers\n")
	}
	return(invisible("n.predictor arg checked"))
}

check.xval.cols.arg <- function(args){
	if(!is.null(args[["mod.cols"]])){
		if(length(args[["mod.cols"]]) != length(args[["xval.files"]])){
			stop("\nyou must specify one plotting color for each model\n")
		}
	}
	return(invisible("xval.mod.cols checked"))
}

read.xval.results <- function(xval.file){
	xval.result <- data.matrix(utils::read.table(xval.file,stringsAsFactors=FALSE,header=TRUE))
	return(xval.result)
}

standardize.xval.results <- function(xval.results){
	n.models <- length(xval.results)
	n.replicates <- ncol(xval.results[[1]])
	n.partitions <- nrow(xval.results[[1]])
	std.xval.array <- array(NA,dim=c(n.partitions,n.replicates,n.models))
	for(i in 1:n.models){
		std.xval.array[,,i] <- xval.results[[i]]
	}
	for(i in 1:n.partitions){
		for(j in 1:n.replicates){
			std.xval.array[i,j,] <- std.xval.array[i,j,] - max(std.xval.array[i,j,]) 
		}
	}
	std.xval.results <- lapply(1:n.models,function(i){std.xval.array[,,i]})
	return(std.xval.results)
}

summarize.mod.xval <- function(mod.xval.results){
	xval.mean <- mean(mod.xval.results)
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

xval.bedassle.rep <- function(rep,n.partitions,geoDist,envDist,nLoci,prefix,n.chains,n.iter,save.files,...){
	kfold.partitions <- make.kfold.partitions(rep,n.partitions)
    x.val <- lapply(1:n.partitions,function(k){
				    message(sprintf("\n\tk-fold cross-validation: analyzing partition %s/%s\n",k,n.partitions));
        				xval.bedassle.fold(test.partition = rep[[k]],
        								   nLoci.test = nLoci/n.partitions,
        								   genDist = kfold.partitions[[k]],
        								   geoDist = geoDist,
        								   envDist = envDist,
        								   nLoci = (n.partitions-1)*(nLoci/n.partitions),
        								   prefix = paste0(prefix,"part",k),
        								   n.chains = n.chains,
        								   n.iter = n.iter,
        								   save.files,
        								   ...)
    			 })
    names(x.val) <- paste0("partition_", 1:n.partitions)
    return(unlist(x.val))
}

xval.bedassle.fold <- function(test.partition,nLoci.test,genDist,geoDist=NULL,envDist=NULL,nLoci,prefix,n.chains=2,n.iter=2e3,save.files=FALSE,...) {
	post.par.cov <- xval.bedassle(genDist = genDist,
								  geoDist = geoDist,
								  envDist = envDist,
								  nLoci = nLoci,
								  prefix = prefix,
								  n.chains = n.chains,
								  n.iter = 2e3,
								  save.files,
								  ...)
	lnL.test.partition <- mean(
							unlist(
								lapply(post.par.cov,function(ppc){
									fit.to.test(test.partition,ppc,nLoci=nLoci.test)
								})
							))
	return(lnL.test.partition)
}

xval.bedassle <- function(genDist,geoDist=NULL,envDist=NULL,nLoci,prefix,n.chains=2,n.iter=2e3,save.files=FALSE,...) {
	data.block <- make.data.block(genDist,geoDist,envDist,nLoci,silent=TRUE)
	stan.model <- pick.stan.model(geoDist,envDist)
	model.fit <- rstan::sampling(object = stanmodels[[stan.model]],
							 	 refresh = min(n.iter/10,500),
							 	 data = data.block,
							 	 iter = n.iter,
							 	 chains = n.chains,
							 	 thin = ifelse(n.iter/500 > 1,n.iter/500,1),
							 	 save_warmup = FALSE,
							 	 ...)
	if(save.files){
		data.block <- unstandardize.distances(data.block)
		save(data.block,file=paste(prefix,"data.block.Robj",sep="_"))
		save(model.fit,file=paste(prefix,"model.fit.Robj",sep="_"))
		bedassle.results <- get.bedassle.results(data.block,model.fit,n.chains)
		write.bedassle.results(data.block,bedassle.results,prefix,n.chains)
		make.all.bedassle.plots(paste0(prefix,"_posterior_chain",1:n.chains,".txt"),paste(prefix,"data.block.Robj",sep="_"),prefix)
	}
	posterior.par.cov <-lapply(1:n.chains,function(i){get.par.cov(model.fit,i,data.block$N)})
	return(posterior.par.cov)
}

load.partitions.file <- function(partitions.file){
	tmpenv <- environment()
	tmp <- load(partitions.file,envir=tmpenv)
	rep.partitions <- lapply(tmp,get,envir=tmpenv)
	names(rep.partitions) <- tmp
	return(rep.partitions[[1]])
	return(rep.partitions)
}

make.kfold.partitions <- function(data.partitions,n.partitions){
	kfold.partitions <- lapply(1:n.partitions,
							function(k){
								make.kfold.partition(k,data.partitions,n.partitions)
							})
	return(kfold.partitions)
}

make.kfold.partition <- function(k,data.partitions,n.partitions){
	data.partitions[[k]] <- NULL
	kfold.partition <- Reduce("+",data.partitions)/(n.partitions-1)
	return(kfold.partition)
}

announce.xval.procedure <- function(n.replicates,n.partitions,L,geoDist,envDist){
	xval.proc <- sprintf("running %s %s of k-fold cross-validation on:\n\n\t",n.replicates, ifelse(n.replicates>1,"replicates","replicate"))
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
					sprintf("\n\nwith %s data partitions, each comprised of %s loci",
							n.partitions,L/n.partitions))
	return(xval.proc)
}

write.xvals <- function(xvals,n.replicates,n.partitions,prefix){
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
			if(is.null(args[["n.nodes"]])){
				n.nodes <- parallel::detectCores()-1
			} else {
				n.nodes <- args[["n.nodes"]]
			}
			cl <- parallel::makeCluster(n.nodes)
			doParallel::registerDoParallel(cl)
			message("\nRegistered doParallel with ",n.nodes," workers\n")
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


log.likelihood <- function(obsCov,inv.par.cov,log.det,nLoci){
	lnL <- -0.5 * (sum( inv.par.cov * obsCov) + nLoci * log.det)
	return(lnL)
}

calc.lnl.x.MCMC <- function(cov.chunk,pp.par.cov,chunk.size){
	lnl.x.mcmc <- lapply(pp.par.cov,
						function(x){
							log.likelihood(cov.chunk,x$inv,x$log.det,nLoci=chunk.size)
					})
	return(unlist(lnl.x.mcmc))
}

fit.to.test <- function(test.partition,parCov,nLoci){
	pp.par.cov <- post.process.par.cov(parCov)
	test.lnl <- lapply(pp.par.cov,
						function(x){
							log.likelihood(test.partition,x$inv,x$log.det,nLoci)
				})
	return(test.lnl)
}

check.xval.call <- function(args){
	check.partitions.file(args)
	check.xval.misc.args(args)
	check.xval.dist.args(args)
	check.nLoci.args(args)
	check.parallel.args(args)
	check.for.output.files(args)
	return(invisible("args checked"))
}

check.xval.misc.args <- function(args){
	if(!is.whole.number(args[["n.replicates"]])){
		stop("\nyou must specify an integer value for the \"n.replicates\" argument\n")	
	}
	if(!is.whole.number(args[["n.partitions"]])){
		stop("\nyou must specify an integer value for the \"n.partitions\" argument\n")	
	}
	if(!is.whole.number(args[["nLoci"]])){
		stop("\nyou must specify an integer value for the \"nLoci\" argument\n")	
	}
	if(!is.character(args[["prefix"]])){
		stop("\nyou must specify a character value for the \"prefix\" argument\n")		
	}
	if(!is.whole.number(args[["n.iter"]])){
		stop("\nyou must specify an integer value for the \"n.iter\" argument\n")	
	}
	if(!is.whole.number(args[["n.chains"]])){
		stop("\nyou must specify an integer value for the \"n.chains\" argument\n")	
	}
	if(!is.logical(args[["save.files"]])){
		stop("\nyou must specify a logical value (TRUE or FALSE) for the \"save.files\" argument\n")
	}
	return(invisible("misc args checked"))
}

check.partitions.file <- function(args){
	if(!is.character(args[["partitions.file"]])){
		stop("\nyou must specify a character vector for the \"partitions.file\" argument\n")
	}
	if(!file.exists(args[["partitions.file"]])){
		stop("\nfunction is unable to find the \"partitions.file\" specified\n")
	}
	return(invisible("partition files checked"))
}

check.xval.dist.args <- function(args){
	check.geoDist.arg(args)
	check.envDist.arg(args)
	return(invisible("dist args checked"))
}

check.nLoci.args <- function(args){
	if(args[["nLoci"]] %% args[["n.partitions"]] != 0){
		stop("\n\"nLoci\" is not cleanly divisible by \"n.partitions\",\nbut you must have an equal number of loci in each partition\n")
	}
	partition.nLoci <- args[["nLoci"]]/args[["n.partitions"]]
	if(!is.null(args[["geoDist"]])){
		if(nrow(args[["geoDist"]]) >= partition.nLoci){
			stop("\nyour data must have a greater number of loci than there are samples\n")
		}
	}
	if(!is.null(args[["envDist"]])){
		if(class(args[["envDist"]]) == "matrix"){
			N <- nrow(args[["envDist"]])
		}
		if(class(args[["envDist"]]) == "list"){
			N <- nrow(args[["envDist"]][[1]])
		}
		if(N >= partition.nLoci){
			stop("\nyour data must have a greater number of loci than there are samples\n")
		}
	}
	return("nLoci checked")
}

check.parallel.args <- function(args){
	if(args[["parallel"]] & args[["n.nodes"]]==1){
		stop("\nyou have specified the \"parallel\" option with \"n.nodes\" set to 1.\n\n")
	}
	if(!args[["parallel"]] & args[["n.nodes"]] > 1){
		stop("\nyou have are running with \"parallel\" set to FALSE but with \"n.nodes\" greater than 1.\n\n")
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

check.rep.partitions.arg <- function(rep.partitions,n.replicates,n.partitions,geoDist,envDist){
	if(length(rep.partitions) != n.replicates){
		stop("\nthe number of replicates in the \"partitions.file\" argument does not match the number of replicates specified\n")
	}
	if(class(rep.partitions) != "list"){
		stop("\nthe rep.partitions file must be of class \"list\"\\n")
	}
	n.row <- get.nrow(geoDist,envDist)
	lapply(1:length(rep.partitions),
		function(i){
			check.rep(rep.partitions[[i]],sprintf("rep_%s",i),n.partitions,n.row)
		})
	return(invisible("rep partitions arg checked"))
}

check.rep <- function(rep,rep.name,n.partitions,n.row){
	if(length(rep) != n.partitions){
		stop(sprintf("\nthe number of partitions in %s of the \"partitions.file\" argument does not match the number of partitions specified\n",rep.name))
	}
	if(class(rep) != "list"){
		stop("\nthe data partitions object within each cross-validation replicate must be of class \"list\"\n")
	}
	lapply(1:length(rep),
		function(k){
			check.partition(rep[[k]],rep.name,sprintf("partition_%s",k),n.row)
		})
	return(invisible("rep checked"))
}

check.partition <- function(partition,rep.name,partition.name,n.row){
	check.any.dist(partition,paste0(partition.name," in ",rep.name))
	if(!is.null(n.row)){
		if(nrow(partition) != n.row){
			stop(sprintf("\nyou have specified a data partition (%s in %s) that does not match the dimensions of your geographic and/or environmental distances\n",partition.name,rep.name))
		}
	}
	allCov <- pwp2allelicCov(partition)
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
		if(class(envDist)=="matrix"){
			n.row <- nrow(envDist)
		}
		if(class(envDist)=="list"){
			n.row <- nrow(envDist[[1]])
		}
	}
	return(n.row)
}