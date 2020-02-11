#' Run a BEDASSLE analysis.
#'
#' \code{run.bedassle} runs a BEDASSLE analysis
#'
#' This function runs an analysis that 
#' estimates the relative contributions of geographic 
#' and environmental/ecological distances to patterns of genetic 
#' differentiation between samples.
#'
#' @param genDist A \code{matrix} of pairwise pi measured between 
#'					all pairs of samples.
#' @param geoDist A \code{matrix} of pairwise geographic distances 
#'					measured between all pairs of samples. A value of 
#'					\code{NULL} runs a model without geographic distance 
#'					as a predictor of genetic differentiation.
#' @param envDist A \code{matrix} of pairwise environmental distances 
#'					measured between all pairs of samples. If there 
#'					are multiple environmental distance measures, this 
#'					argument should be a \code{list} of distance matrices. 
#'				 	A value of \code{NULL} runs a model without geographic 
#'					distance as a predictor of genetic differentiation.
#' @param nLoci The number of loci used in the calculation of the pairwise 
#'					pi \code{matrix} specified in the \code{genDist} argument.
#' @param prefix A character \code{vector} giving the prefix to be attached 
#'					 to all output files. An underscore is automatically 
#'					 added between the prefix and the file names.
#' @param n.chains An integer indicating the number of MCMC chains to be run 
#'					in the analysis. Default is 4.
#' @param n.iter An \code{integer} giving the number of iterations each MCMC 
#'				 chain is run. Default is 2e3.  If the number of iterations 
#'				 is greater than 500, the MCMC is thinned so that the number 
#'				 of retained iterations is 500 (before burn-in).
#' @param make.figs A \code{logical} value indicating whether to automatically 
#'					make figures once the analysis is complete. Default is 
#'					\code{TRUE}.
#' @param save.stanfit A \code{logical} value indicating whether to automatically 
#'						save the full \code{stanfit} output once the analysis is
#'						complete. Default is \code{TRUE}.
#' @param ... Further options to be passed to rstan::sampling (e.g., adapt_delta).
#'
#' @return This function writes all output to three tab-delimited text files 
#'			for each chain run (specified with \code{n.chains}). The output files 
#'			associated with each chain have "chain_X" for the Xth chain appended 
#'			to the file name, and begin with the prefix specified with the 
#'			\code{prefix} argument. The three output files associated with each 
#'			chain are described below:
#'			\itemize{
#'				\item \code{posterior} contains the posterior probability and 
#'						parameter estimates over the sampled posterior distribution
#'						of the MCMC. Each entry described below has its own named 
#'						column in the output text file. Columns are separated by tabs.
#'					\itemize{
#'						\item \code{lpd} log posterior density over the retained 
#'								MCMC iterations.
#'						\item \code{alpha0} posterior draws for alpha0 parameter.
#'						\item \code{alphaD} posterior draws for alphaD parameter.
#'						\item \code{alphaE} posterior draws for alphaE parameter(s).
#'								Named "alphaE_E" for the Eth environmental distance 
#'								matrix specified.
#'						\item \code{alpha2} posterior draws for alpha2 parameter.
#'						\item \code{nuggets} posterior draws for nugget parameters. 
#'								Named "nugget_N" for the Nth sample in the dataset.
#'					}
#'				\item \code{MAP} contains point estimates of the parameters listed in 
#'								the \code{posterior} file described above. Values are 
#'								indexed at the MCMC iteration with the greatest 
#'								posterior probability.
#'				\item \code{parCov} contains the covariance matrix parameterized by the 
#'								MAP parameter point estimates.
#'			}
#'
#' @details This function acts as a wrapper around a STAN model block determined 
#'			by the user-specified model (e.g., just geographic distance, or geographic distance 
#'			plus 2 environmental/ecological distance variables).
#'			User-specified data are checked for appropriate format and consistent dimensions,
#'			then formatted into a \code{data.block},
#'			which is then passed to the STAN model block.
#'			The model output is written to three tab-delimited text files described above.
#'			The parameter values are rounded to 4 decimal places
#'			The full \code{stanfit} model output is also saved if \code{save.stanfit=TRUE}.
#'			If \code{make.figs=TRUE}, running \code{run.bedassle} will also generate figures
#'			depicting different aspects of model output; these are detailed in the function 
#'			\code{make.all.bedassle.plots} in this package.
#'
#' @import rstan
#' @export
run.bedassle <- function(genDist,geoDist=NULL,envDist=NULL,nLoci,prefix,n.chains=4,n.iter=2e3,make.figs=TRUE,save.stanfit=TRUE,...){
	call.check <- check.call(args <- as.list(environment()))
	data.block <- make.data.block(genDist,geoDist,envDist,nLoci,silent=FALSE)
	stan.model <- pick.stan.model(geoDist,envDist)
	model.fit <- rstan::sampling(object = stanmodels[[stan.model]],
							 	 refresh = min(n.iter/10,500),
							 	 data = data.block,
							 	 iter = n.iter,
							 	 chains = n.chains,
							 	 thin = ifelse(n.iter/500 > 1,n.iter/500,1),
							 	 save_warmup = FALSE,
							 	 ...)
	data.block <- unstandardize.distances(data.block)
		save(data.block,file=paste(prefix,"data.block.Robj",sep="_"))
	bedassle.results <- get.bedassle.results(data.block,model.fit,n.chains)
	write.bedassle.results(data.block,bedassle.results,prefix,n.chains)
	if(save.stanfit){
		save(model.fit,file=paste(prefix,"model.fit.Robj",sep="_"))
	}
	if(make.figs){
		make.all.bedassle.plots(paste0(prefix,"_posterior_chain",1:n.chains,".txt"),paste(prefix,"data.block.Robj",sep="_"),prefix)
	}
	return("analysis complete")
}

check.call <- function(args){
	check.genDist.arg(args)
	check.misc.args(args)
	check.geoDist.arg(args)
	check.envDist.arg(args)
	check.dist.args(args)
	return(invisible("args checked"))
}

check.any.dist <- function(distArg,argName){
	if(class(distArg) != "matrix"){
		stop(sprintf("\nthe \"%s\" argument must be of class \"matrix\"\n",argName))
	}
	if(length(unique(dim(distArg))) > 1){
		stop(sprintf("\nyou have specified a \"%s\" argument with an unequal number of rows and columns\n",argName))	
	}
	if(any(!is.finite(distArg))){
		stop(sprintf("\nall values of the \"%s\" argument must be finite\n",argName))
	}
	if(any(distArg < 0)){
		stop(sprintf("\nall values of the \"%s\" argument must be greater than 0\n",argName))
	}
	tmpDist <- distArg
	row.names(tmpDist) <- NULL
	colnames(tmpDist) <- NULL
	if(!isSymmetric(tmpDist)){	
		stop(sprintf("\nyou must specify a symmetric matrix for the \"%s\" argument \n",argName))
	}
	return(invisible("this dist arg checked"))
}

check.genDist.arg <- function(args){
	check.any.dist(distArg = args[["genDist"]],argName="genDist")
	if(any(args[["genDist"]] > 1)){
		stop("\n all values of the genetic distance matrix specified by the \"genDist\" argument must be less than 1\n")
	}
	return(invisible("genDist arg checked"))
}

check.geoDist.arg <- function(args){
	if(!is.null(args[["geoDist"]])){
		check.any.dist(distArg = args[["geoDist"]],argName="geoDist")
	}
	return(invisible("geoDist arg checked"))
}

check.envDist.arg <- function(args){
	if(!is.null(args[["envDist"]])){
		if(class(args[["envDist"]]) != "list" & class(args[["envDist"]]) != "matrix"){
			stop("\nthe \"envDist\" argument must either be a matrix or a list of matrices\n")
		}
		if(class(args[["envDist"]]) == "matrix"){
			check.any.dist(distArg = args[["envDist"]],argName="envDist")
		}
		if(class(args[["envDist"]]) == "list"){
			lapply(args[["envDist"]],function(E){
				if(class(E) != "matrix"){
					stop("\neach element of the \"envDist\" list must be of class \"matrix\"\n")
				}
				if(length(unique(dim(E))) > 1){
					stop("\nyou have specified a \"envDist\" list that contains a matrix with an unequal number of rows and columns\n")	
				}
				if(any(!is.finite(E))){
					stop("\nall values of each \"envDist\" distance matrix must be finite\n")
				}
				if(any(E < 0)){
					stop("\nall values of each \"envDist\" distance matrix must be greater than 0\n")
				}
				tmp.envDist <- E
				row.names(tmp.envDist) <- NULL
				colnames(tmp.envDist) <- NULL
				if(!isSymmetric(tmp.envDist)){	
					stop("\neach \"envDist\" distance matrix must be symmetric\n")
				}
			})
			if(length(unique(lapply(args[["envDist"]],dim))) > 1){
				stop("\neach \"envDist\" distance matrix must have the same dimensions\n")
			}
		}
	}
	return(invisible("envDist arg checked"))
}

check.dist.args <- function(args){
	if(!is.null(args[["geoDist"]])){
		if(nrow(args[["geoDist"]]) != nrow(args[["genDist"]])){
			stop("\nthe dimensions of the \"geoDist\" argument must match those of the \"genDist\" argument\n")
		}
	}
	if(!is.null(args[["envDist"]])){
		if(class(args[["envDist"]]) == "matrix"){
			if(nrow(args[["envDist"]]) != nrow(args[["genDist"]])){
				stop("\nthe dimensions of the \"envDist\" argument must match those of the \"genDist\" argument\n")
			}
		}
		if(class(args[["envDist"]]) == "list"){
			if(nrow(args[["envDist"]][[1]]) != nrow(args[["genDist"]])){
				stop("\nthe dimensions of the \"envDist\" distance matrices must match those of the \"genDist\" argument\n")
			}
		}
	}
	return(invisible("all dist args checked"))
}

check.misc.args <- function(args){
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
	if(!is.logical(args[["make.figs"]])){
		stop("\nyou must specify a logical value (TRUE or FALSE) for the \"make.figs\" argument\n")
	}
	if(!is.logical(args[["save.stanfit"]])){
		stop("\nyou must specify a logical value (TRUE or FALSE) for the \"save.stanfit\" argument\n")
	}
	return(invisible("misc args checked"))
}

pick.stan.model <- function(geoDist,envDist){
	if(is.null(geoDist) & is.null(envDist)){
		name <- "NULL"
	}
	if(!is.null(geoDist) & is.null(envDist)){
		name <- "GEO"
	}
	if(is.null(geoDist) & !is.null(envDist)){
		name <- "ENV"
	}
	if(!is.null(geoDist) & !is.null(envDist)){
		name <- "GEOENV"
	}
	return(name)
}

make.envDist.list <- function(envDist){
	if(!is.null(envDist)){
		if(class(envDist)=="matrix"){
			envDist <- list(envDist)
		}
	}
	return(envDist)
}

make.envDist.array <- function(std.envDist){
	if(!is.null(std.envDist[[1]]$std.X)){
		nE <- length(std.envDist)
		N <- nrow(std.envDist[[1]]$std.X)
		E <- array(0,dim=c(length(std.envDist),N,N))
		for(e in 1:nE){
			E[e,,] <- std.envDist[[e]]$std.X
		}
	} else {
		E <- NULL
	}
	return(E)
}

standardize.distances <- function(X){
	if(!is.null(X)){
		stdev.X <- stats::sd(X[upper.tri(X)])
		std.X <- X/stdev.X
	} else {
		std.X <- NULL
		stdev.X <- NULL
	}
	sd.dist.list <- list("std.X" = std.X,
						"stdev.X" = stdev.X)
	return(sd.dist.list)
}

describe.picked.model <- function(data.block){
	if(is.null(data.block$geoDist) & is.null(data.block$envDist)){
		chosenMod <- "the \"NULL\" model (neither geographic or environmental distance)"
	}
	if(!is.null(data.block$geoDist) & is.null(data.block$envDist)){
		chosenMod <- "the \"GEO\" model (geographic, but not environmental distance)"
	}
	if(is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		chosenMod <- "the \"ENV\" model (environmental, but not geographic distance)"
	}
	if(!is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		chosenMod <- "the \"GEOENV\" model (geographic and environmental distance)"
	}
	if(!is.null(data.block$envDist)){
		nE <- dim(data.block$envDist)[1]
		vars <- ifelse(nE == 1,"matrix","matrices")
		envAddendum <- sprintf("\n\twith %s environmental distance %s\n",nE,vars)
		chosenMod <- paste0(chosenMod,envAddendum)
	}
	return(chosenMod)
}

validate.data.block <- function(data.block,silent=FALSE){
	if(!silent){
		message("\nchecking data\n")
		message(sprintf("\treading %s samples",data.block$N))
		message(sprintf("\treading %s loci",data.block$L))
		message("\nchecking specified model\n")
		chosenMod <- describe.picked.model(data.block)
		message(sprintf("\tuser has specified %s \n",chosenMod))
	}
	if(!data.block$L > data.block$N){
		stop("\nyour data must have a greater number of loci than there are samples\n")
	}
	data.block <- make.data.block.S3(data.block)
	return(data.block)
}

make.data.block.S3 <- function(data.block){
	data.block <- data.block
	class(data.block) <- "data.block"
	return(data.block)
}

print.data.block <- function(data.block){
	print(utils::str(data.block,max.level=1))
}

standardize.envDist <- function(envDist){
	if(!is.null(envDist)){
		std.envDist <- lapply(envDist,standardize.distances)
	} else {
		std.envDist <- list(list("std.X" = NULL,
							"stdev.X" = NULL))
	}
	return(std.envDist)
}

make.data.block <- function(genDist,geoDist=NULL,envDist=NULL,nLoci,silent=FALSE){
	allCov <- pwp2allelicCov(genDist)
	if(any(eigen(allCov)$values < 0)){
		stop("\nthe allelic covariance calculated from the specified genetic distance matrix is not positive definite\n")
	}
	envDist <- make.envDist.list(envDist)
	sd.geoDist.list <- standardize.distances(geoDist)
	sd.envDist.list <- standardize.envDist(envDist)
	envDist.array <- make.envDist.array(sd.envDist.list)
	data.block <- list("N" = nrow(genDist),
					   "L" = nLoci,
					   "obsCov" = allCov,
					   "geoDist" = sd.geoDist.list$std.X,
					   "sd.geoDist" = sd.geoDist.list$stdev.X,
					   "envDist" = envDist.array,
					   "sd.envDist" = lapply(sd.envDist.list,"[[", "stdev.X"),
					   "nE" = ifelse(!is.null(envDist),length(envDist),0))
	data.block <- validate.data.block(data.block,silent=silent)
	return(data.block)
}

unstandardize.distances <- function(data.block){
	if(!is.null(data.block$sd.geoDist)){
		data.block$geoDist <- data.block$geoDist*data.block$sd.geoDist
	}
	if(data.block$nE!=0){
		for(e in 1:data.block$nE){
			data.block$envDist[e,,] <- data.block$envDist[e,,] * data.block$sd.envDist[[e]]
		}
	}
	return(data.block)
}

is.whole.number <- function(x){
	val <- FALSE
	if(is.numeric(x)){
		if(x %% 1 == 0){
			val <- TRUE
		}
	}
	return(val)
}