get.bedassle.results <- function(data.block,model.fit,nChains){
	bedassle.results <- stats::setNames(
							lapply(1:nChains,
								function(i){
									get.bedassle.chain.results(data.block,model.fit,i)
								}),
						  paste0("chain_",1:nChains))
	return(bedassle.results)
}

get.bedassle.chain.results <- function(data.block,model.fit,chain.no){
	posterior <- get.posterior(data.block,model.fit,chain.no)
	names(posterior) <- get.bedassle.result.names(data.block)
	if(!is.null(data.block$sd.geoDist)){
		posterior$alphaD <- posterior$alphaD/data.block$sd.geoDist
	}
	if(data.block$nE != 0){
		for(e in 1:data.block$nE){
			posterior[[sprintf("alphaE_%s",e)]] <- posterior[[sprintf("alphaE_%s",e)]]/data.block$sd.envDist[[e]]
		}
	}
	return(posterior)
}

get.posterior <- function(data.block,model.fit,chain.no){
	if(is.null(data.block$geoDist) & is.null(data.block$envDist)){
		posterior <- data.frame("lpd" = rstan::get_logposterior(model.fit)[[chain.no]],
								"alpha0" = get.alpha0(model.fit,chain.no),
							    "nugget" = get.nugget(model.fit,chain.no))
	}
	if(!is.null(data.block$geoDist) & is.null(data.block$envDist)){
		posterior <- data.frame("lpd" = rstan::get_logposterior(model.fit)[[chain.no]],
								"alpha0" = get.alpha0(model.fit,chain.no),
							    "alphaD" = get.alphaD(model.fit,chain.no),
							    "alpha2" = get.alpha2(model.fit,chain.no),
							    "nugget" = get.nugget(model.fit,chain.no))
	}
	if(is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		posterior <- data.frame("lpd" = rstan::get_logposterior(model.fit)[[chain.no]],
								"alpha0" = get.alpha0(model.fit,chain.no),
							    "alphaE" = get.alphaE(model.fit,chain.no),
							    "alpha2" = get.alpha2(model.fit,chain.no),
							    "nugget" = get.nugget(model.fit,chain.no))
	}
	if(!is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		posterior <- data.frame("lpd" = rstan::get_logposterior(model.fit)[[chain.no]],
								"alpha0" = get.alpha0(model.fit,chain.no),
							    "alphaD" = get.alphaD(model.fit,chain.no),
							    "alphaE" = get.alphaE(model.fit,chain.no),
							    "alpha2" = get.alpha2(model.fit,chain.no),
							    "nugget" = get.nugget(model.fit,chain.no))
	}
	return(posterior)
}

write.bedassle.results <- function(data.block,bedassle.results,prefix,nChains){
	lapply(1:nChains,
		function(i){
			write.bedassle.chain.results(data.block,bedassle.results[[i]],prefix,i)
		})
	return(invisible("results printed"))
}

write.bedassle.chain.results <- function(data.block,chain.bedassle.results,prefix,chain.no){
	utils::write.table(round(chain.bedassle.results,4),file=paste0(prefix,"_posterior_chain",chain.no,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
	MAP <- chain.bedassle.results[which.max(chain.bedassle.results$lpd),]
	utils::write.table(round(MAP,4),file=paste0(prefix,"_MAP_chain",chain.no,".txt"),quote=FALSE,row.names=FALSE,sep="\t")
	parCov <- get.cov.function(data.block)(data.block,MAP)
	utils::write.table(round(parCov,4),file=paste0(prefix,"_parCov_chain",chain.no,".txt"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
	return(invisible("chain results printed"))
}

get.bedassle.result.names <- function(data.block){
	if(is.null(data.block$geoDist) & is.null(data.block$envDist)){
		post.names <- c("lpd","alpha0","nugget")	
	}
	if(!is.null(data.block$geoDist) & is.null(data.block$envDist)){
		post.names <- c("lpd","alpha0","alphaD","alpha2","nugget")	
	}
	if(is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		post.names <- c("lpd","alpha0",paste0("alphaE","_",1:data.block$nE),"alpha2","nugget")
	}
	if(!is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		post.names <- c("lpd","alpha0","alphaD",paste0("alphaE","_",1:data.block$nE),"alpha2","nugget")	
	}
	return(post.names)
}

get.par <- function(model.fit,par,chain.no){
	par <- rstan::extract(model.fit,pars=par,inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(par)
}

get.nugget <- function(model.fit,chain.no){
	nugget <- get.par(model.fit,"nugget",chain.no)
	return(nugget)
}

get.alpha0 <- function(model.fit,chain.no){
	alpha0 <- NULL
	if(any(grepl("alpha0",model.fit@model_pars))){
		alpha0 <- get.par(model.fit,"alpha0",chain.no)
	}
	return(alpha0)
}

get.alphaD <- function(model.fit,chain.no){
	alphaD <- NULL
	if(any(grepl("alphaD",model.fit@model_pars))){
		alphaD <- get.par(model.fit,"alphaD",chain.no)
	}
	return(alphaD)
}

get.alphaE <- function(model.fit,chain.no){
	alphaE <- NULL
	if(any(grepl("alphaE",model.fit@model_pars))){
		alphaE <- get.par(model.fit,"alphaE",chain.no)
	}
	return(alphaE)
}

get.alpha2 <- function(model.fit,chain.no){
	alpha2 <- NULL
	if(any(grepl("alpha2",model.fit@model_pars))){
		alpha2 <- get.par(model.fit,"alpha2",chain.no)
	}
	return(alpha2)
}

get.MAP.par.cov <- function(data.block,MAP){
	cov.func <- get.cov.function(data.block)
	MAP.par.cov <- cov.func(data.block,MAP)
	return(MAP.par.cov)
}

get.cov.function <- function(data.block){
	if(is.null(data.block$geoDist) & is.null(data.block$envDist)){
		cov.func <- function(data.block,MAP){
			return(MAP$alpha0 + diag(MAP[grep("nugget",names(MAP))],data.block$N))
		}
	}
	if(!is.null(data.block$geoDist) & is.null(data.block$envDist)){
		cov.func <- function(data.block,MAP){
			return(MAP$alpha0 * exp(-(MAP$alphaD*data.block$geoDist)^MAP$alpha2) + 
					diag(MAP[grep("nugget",names(MAP))],data.block$N))
		}
	}
	if(is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		cov.func <- function(data.block,MAP){
			return(MAP$alpha0 * 
						exp(-sqrt(Reduce("+",
								lapply(1:data.block$nE,function(e){
									(unlist(MAP[grep("alphaE",names(MAP))])[e]*data.block$envDist[e,,])^2
								}))
							)^MAP$alpha2) + 
					diag(MAP[grep("nugget",names(MAP))],data.block$N))
		}
	}
	if(!is.null(data.block$geoDist) & !is.null(data.block$envDist)){
		cov.func <- function(data.block,MAP){
			return(MAP$alpha0 * 
						exp(-sqrt(
								(MAP$alphaD*data.block$geoDist)^2 + 
								Reduce("+",
									lapply(1:data.block$nE,function(e){
										(unlist(MAP[grep("alphaE",names(MAP))])[e]*data.block$envDist[e,,])^2
									}))
							)^MAP$alpha2) + 
					diag(MAP[grep("nugget",names(MAP))],data.block$N))
		}
	}
	return(cov.func)
}