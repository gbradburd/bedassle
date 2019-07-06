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