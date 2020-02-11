#' Make output plots
#'
#' \code{make.all.bedassle.plots} makes figures from the output from a 
#' 	BEDASSLE analysis.
#'
#' This function takes the file output from a BEDASSLE analysis and 
#' generates a number of plots for visualizing results and 
#' diagnosing MCMC performance.
#'
#' @param results.files A \code{vector} of the filenames (in quotes, 
#'			with the full file path) of the posterior results files 
#'			output by the different chains of a \code{BEDASSLE} run. 
#'			Can also be a single filename if user wants to visualize 
#'			only a single run.
#' @param data.block.file The filename (in quotes, with the full file 
#'			path) of the data.block R object (.Robj) file output by a 
#'			\code{BEDASSLE} run.
#' @param prefix A character vector to be prepended to all figures.
#' @param chain.cols A \code{vector} of colors to be used in plotting 
#'			results from different MCMC chains. There should be one 
#'			color specified for each chain. If \code{NULL}, the plots 
#'			will use the required set or subset of 12 pre-specified 
#'			colors. If there are more than 12 chains, users must supply 
#'			their own colors.
#' @return This function has only invisible return values.
#'
#'	@details This function produces a variety of plots that can be 
#'	useful for visualizing results or diagnosing MCMC performance. 
#'  The plots made are by no means exhaustive, and users are 
#' 	encouraged to make further plots, or customize these plots as they 
#'	see fit. The plots generated (as .pdf files) are:
#'	\itemize{
#'		\item model.fit.CIs - A plot of the sample allelic covariance 
#'			shown with the 95\% credible interval of the parametric 
#'			covariance for each entry in the matrix. Only generated 
#'			if either the \code{geoDist} or \code{envDist} arguments 
#'			in the \code{run.bedassle} function call are specified. One 
#'			plot is produced for each chain.
#'		\item Trace plots - Plots of parameter values over the MCMC.
#'		\itemize{
#'			\item lpd - A plot of the log posterior probability over the MCMC.
#'			\item nuggets - A plot of estimates of the nugget parameters 
#'				over the MCMC.
#'			\item alpha parameters - Plots of estimates of the 
#'				various parameters (all or some of {alpha0,alphaD,alphaE,alpha2}, 
#'				depending on the model specified) over the MCMC.
#'		}
#'	}
#' 
#'@export
make.all.bedassle.plots <- function(results.files,data.block.file,prefix,chain.cols=NULL){
	check.results.args(args <- as.list(environment()))
	if(is.null(chain.cols)){
		chain.cols <- c("#332288E6","#882255E6","#117733E6",
						"#88CCEEE6","#AA4466E6","#CC6677E6",
						"#44AA99E6","#661100E6","#6699CCE6",
						"#999933E6","#AA4499E6","#DDCC77E6")[1:length(results.files)] 
	}
	if(length(results.files) > length(chain.cols)){
		stop("\nthere are only 12 default chain colors, but you have specified > 12 chains,\nso you must specify your own chain colors\n")
	}
	bedassle.results <- lapply(results.files,
							function(x){
								utils::read.table(x,stringsAsFactors=FALSE,header=TRUE)
						})
	data.block <- load.data.block(data.block.file)
	grDevices::pdf(file=paste0(prefix,"_trace.plots.pdf"))
		plot.lpd(bedassle.results,chain.cols)
		plot.nuggets(bedassle.results,chain.cols)
		if(is.null(data.block$geoDist) & is.null(data.block$envDist)){
			plot.alpha0(bedassle.results,chain.cols)
		}
		if(!is.null(data.block$geoDist) | !is.null(data.block$envDist)){
			plot.alpha.params(data.block,bedassle.results,chain.cols)
		}
	grDevices::dev.off()
	if(!is.null(data.block$geoDist) | !is.null(data.block$envDist)){
		grDevices::pdf(file=paste0(prefix,"_model.fit.CIs.pdf"))
			plot.model.fit.CIs(data.block,bedassle.results)
		grDevices::dev.off()
	}
	return(invisible("made plots!"))
}

check.results.args <- function(args){
	if(!is.vector(args[["results.files"]])){
		stop("\nyou must specify a vector of filenames for the \"results.files\" argument\n")
	}
	lapply(args[["results.files"]],
		function(f){
			if(!file.exists(f)){
				stop("\nfunction is unable to find a specified \"results.file\"\n")
			}
		})
	if(!file.exists(args[["data.block.file"]])){
		stop("\nfunction is unable to find the specified \"data.block\" file\n")
	}
	if(!is.null(args[["chain.cols"]])){
		if(length(args[["results.files"]]) != length(args[["chain.cols"]])){
			stop("\nyou must specify one plotting color for each chain\n")	
		}
	}
	return(invisible("plot args checked"))
}

extract.x <- function(bedassle.results,x){
	x <- Reduce("rbind",lapply(bedassle.results,"[[",x))
	if(!is.matrix(x)){
		x <- matrix(x,1,length(x))
	}
	return(x)
}

extract.aE <- function(data.block,bedassle.results,x){
	aE <- lapply(1:data.block$nE,
				function(e){
					Reduce("rbind",
						lapply(bedassle.results,
							function(x){
								x[[grep(sprintf("alphaE_%s",e),names(x))]]
						})
					)
				}
		  )
	aE <- lapply(aE,
			function(z){
				if(!is.matrix(z)){
					matrix(z,1,length(z))
				} else {
					z
				}
			}
		  )
	return(aE)
}

plot.lpd <- function(bedassle.results,chain.cols){
	lpd <- extract.x(bedassle.results,"lpd")
	graphics::matplot(t(lpd),
			ylab="posterior probability",
			main="Posterior probability",
			type='l',lty=1,lwd=1.2,
			xlab="MCMC iterations",col=chain.cols)
	return(invisible(0))
}

plot.alpha0 <- function(bedassle.results,chain.cols){
	alpha0 <- extract.x(bedassle.results,"alpha0")
	graphics::matplot(t(alpha0),
			ylab="alpha0",
			xlab="MCMC iterations",
			type='l',lty=1,lwd=1.2,
			main="alpha0",col=chain.cols)
	return(invisible(0))
}

plot.nuggets <- function(bedassle.results,chain.cols){
	nuggets <- lapply(bedassle.results,function(x){x[grep("nugget",names(bedassle.results[[1]]))]})
	graphics::matplot(0,type='n',
				main="sample nuggets",
				ylab="nugget value",
				xlab="MCMC iterations",
				xlim=c(0,nrow(nuggets[[1]])),ylim=range(unlist(nuggets)))
	lapply(1:length(nuggets),function(n){
		graphics::matplot(nuggets[[n]],type='l',lty=1,lwd=0.8,add=TRUE,col=chain.cols[n])
	})
	return(invisible("nuggets"))
}

plot.alpha.params <- function(data.block,bedassle.results,chain.cols){
	if(!is.null(data.block$geoDist) & is.null(data.block$envDist)){
		a0 <- extract.x(bedassle.results,"alpha0")
		aD <- extract.x(bedassle.results,"alphaD")
		a2 <- extract.x(bedassle.results,"alpha2")
		plot.alpha.param(a0,"alpha0",chain.cols)
		plot.alpha.param(aD,"alphaD",chain.cols)
		plot.alpha.param(a2,"alpha2",chain.cols)
	} else if(is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		a0 <- extract.x(bedassle.results,"alpha0")
		aE <- extract.aE(data.block,bedassle.results)
		a2 <- extract.x(bedassle.results,"alpha2")
		plot.alpha.param(a0,"alpha0",chain.cols)
		lapply(1:data.block$nE,
			function(e){
				plot.alpha.param(aE[[e]],sprintf("alphaE_%s",e),chain.cols)
		})
		plot.alpha.param(a2,"alpha2",chain.cols)
	} else if(!is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		a0 <- extract.x(bedassle.results,"alpha0")
		aD <- extract.x(bedassle.results,"alphaD")
		aE <- extract.aE(data.block,bedassle.results)
		a2 <- extract.x(bedassle.results,"alpha2")
		plot.alpha.param(a0,"alpha0",chain.cols)
		plot.alpha.param(aD,"alphaD",chain.cols)
		lapply(1:data.block$nE,
			function(e){
				plot.alpha.param(aE[[e]],sprintf("alphaE_%s",e),chain.cols)
		})
		plot.alpha.param(a2,"alpha2",chain.cols)
	}
	return(invisible("alpha params")) 
}

plot.alpha.param <- function(param,param.name,chain.cols){
	graphics::matplot(t(param),
			ylab=param.name,
			xlab="MCMC iterations",
			type='l',lty=1,lwd=1.2,
			main=param.name,col=chain.cols)
	return(invisible(0))
}

plot.model.fit.CIs <- function(data.block,bedassle.results){
	lapply(1:length(bedassle.results),function(i){
		plot.chain.model.fit.CIs(data.block,bedassle.results[[i]],i)
	})
	return(invisible("model fit CIs plotted"))
}

plot.chain.model.fit.CIs <- function(data.block,chain.bedassle.results,chain.no){
	cov.func <- get.cov.function(data.block)
	post.par.cov <- lapply(1:nrow(chain.bedassle.results),
						function(i){
							cov.func(data.block,chain.bedassle.results[i,])
						})
	cov.range <- range(c(data.block$obsCov,
						post.par.cov))
	if(!is.null(data.block$geoDist)){
		plotDist <- data.block$geoDist
	} else {
		plotDist <- data.block$envDist[1,,]
	}
	xlab <- ifelse(!is.null(data.block$geoDist),"geographic distance","environmental distance")
	graphics::plot(plotDist,data.block$obsCov,
    	xlab = xlab, 
        ylab = "covariance",
        main=sprintf("Model fit (chain_%s)",chain.no),
        ylim = cov.range, type = "n")
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- get.par.cov.CI(data.block, post.par.cov)
	lapply(1:nrow(combns),
			function(i){
				graphics::segments(x0 = plotDist[combns[i,1],combns[i,2]],
						 y0 = CIs[[i]][1],
						 x1 = plotDist[combns[i,1],combns[i,2]],
						 y1 = CIs[[i]][2],
						 col = grDevices::adjustcolor(1,0.1),
						 lwd=1.5)
			})
	graphics::points(plotDist,data.block$obsCov,col=2,pch=20,cex=0.8)
	graphics::legend(x="topright",legend=c("observed","95% CI"),pch=c(19,NA),lty=c(NA,1),col=c(2,"gray"))
	return(invisible("plotted"))
}

get.par.cov.CI <- function(data.block,post.par.cov){
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- lapply(1:nrow(combns),
				function(i){
					stats::quantile(unlist(lapply(post.par.cov,function(x){x[combns[i,1],combns[i,2]]})),c(0.025,0.975))
				})
	return(CIs)
}

load.data.block <- function(data.block.file){
	tmpenv <- environment()
	tmp <- load(data.block.file,envir=tmpenv)
	data.block <- lapply(tmp,get,envir=tmpenv)
	names(data.block) <- tmp
	return(data.block[[1]])
}