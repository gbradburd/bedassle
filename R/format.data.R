#' Convert pairwise pi to allelic covariance
#' 
#' \code{pwp2allelicCov} converts a measures of 
#'	pairwise pi to allelic covariance
#'
#' This function takes measurements of pairwise pi 
#' and returns the allelic covariance, which can
#' then be used in a BEDASSLE analysis.
#' 
#' @param pi A measurement of nucleotide diversity between 
#' a pair of samples; it is the proportion of sites at which 
#' the samples differ, out of the total number of loci 
#' in the sample. Can be a single value, a \code{vector}, 
#' or a \code{matrix}.
#' 
#' @details This function takes a measurement or set of measurements 
#'		of pairwise pi between samples and converts it to a an 
#'		allelic covariance, which can then be used to run a 
#'		\code{\link{run.bedassle}} analysis. The \code{pi} argument 
#'		specified will most commonly be a \code{matrix} of pairwise pi 
#'		between a set of samples, for which the \emph{ij}th cell gives 
#'		the pairwise pi between samples \emph{i} and \emph{j}. However, 
#'		it will also convert a single value (i.e., pairwise pi between 
#'		a single pair of samples) or a vector of values.
#' 
#' 		Pairwise pi is calculated as the proportion of sites at which a 
#'		pair of samples differs out of the total number of loci in the 
#'		dataset. If it is calculated using only loci that are polymorphic 
#'		in the global dataset, it is called pi at polymorphic sites. If it 
#'		is calculated using all genotyped base-pairs, it is simply pi. Either 
#'		statistic can be converted to allelic covariance and used for running 
#'		BEDASSLE.
#'	
#'	@return This function returns the allelic covariance that corresponds to 
#'		a particular value of pairwise pi. The allelic covariance can then be 
#'		used to run a BEDASSLE analysis run with \code{\link{run.bedassle}}.
#'		
#' @export

pwp2allelicCov <- function(pi){
	check.pi(pi)
	allelic.cov <- (1-2*pi)/4
	return(allelic.cov)
}

#' Convert a dataset from STRUCTURE to BEDASSLE format
#'
#' \code{structure2bedassle} converts a STRUCTURE dataset 
#' to BEDASSLE format
#' 
#' This function takes a population genetics dataset in 
#' STRUCTURE format and converts it to an allele frequency 
#' data table, then calculates pairwise pi between all samples.
#' The matrix of pairwise pi can be used as the \code{genDist}
#' argument in \code{\link{run.bedassle}}.
#' The STRUCTURE file can have one row per individual 
#' and two columns per locus, or one column and two rows 
#' per individual. It can only contain bi-allelic SNPs.
#' Missing data is acceptable, but must be indicated with 
#' a single value throughout the dataset.
#' 
#' @param infile The name and path of the file in STRUCTURE format 
#' 			to be converted to the format used in a \code{bedassle} 
#'			analysis. 
#' @param onerowperind Indicates whether the file format has 
#'		one row per individual (\code{TRUE}) or two rows per 
#'		individual (\code{FALSE}).
#' @param start.loci The index of the first column in the dataset 
#'			that contains genotype data.
#' @param start.samples The index of the first row in the dataset 
#'			that contains genotype data (e.g., after any headers). 
#'			Default value is 1.
#' @param missing.datum The character or value used to denote 
#' 			missing data in the STRUCTURE dataset (often 0 or -9).
#' @param prefix A character \code{vector} giving the prefix (including 
#'			desired directory path) to be attached to output files.
#' @param save.freqs A logical value indicating whether or not 
#'			to save the allele frequency data matrix generated 
#' 			by this function as an R object.
#'
#' @details This function takes a STRUCTURE format data file and 
#'		converts it to a \code{bedassle} format data file.
#'		This function can only be applied to diploid organisms.
#'		The STRUCTURE data file must be a plain text file. 
#'		If there are extraneous lines of text or column headers 
#'		before the data start, those extra lines should be deleted 
#'		by hand or taken into account via the \code{start.samples} 
#'		argument.
#'		
#' 		The STRUCTURE dataset can either be in the ONEROWPERIND=1 
#' 		file format, with one row per individual and two columns 
#' 		per locus, or the ONEROWPERIND=0 format, with two rows and 
#'		one column per individual. The first column of the STRUCTURE 
#' 		dataset should be individual names. There may be any number 
#' 		of other columns that contain non-genotype information before 
#'		the first column that contains genotype data, but there can 
#' 		be no extraneous columns at the end of the dataset, after the 
#' 		genotype data.
#'		
#'		The genotype data must be bi-allelic single nucleotide 
#'		polymorphisms (SNPs). Applying this function to datasets
#'		with more than two alleles per locus may result in cryptic 
#'		failure.
#'	
#'	@return This function returns a matrix of pairwise pi 
#'		that can be used as the \code{genDist} argument in a BEDASSLE 
#'		analysis (\code{\link{run.bedassle}}).  It also saves 
#'		this matrix as a text file ("yourprefix_pwp.txt") so that it can 
#'		be used in future analyses. If the \code{save.freqs} is \code{TRUE},
#'		the allele frequency data matrix generated from the STRUCTURE 
#'		data file is saved as an R data (.RData) object.
#'		
#' @export

structure2bedassle <- function(infile,onerowperind,start.loci,start.samples=1,missing.datum,prefix,save.freqs=TRUE){
	#add checks
	outfile.freqs <- paste0(prefix,"_freqs.RData")
	outfile.pwp <- paste0(prefix,"_pwp.txt")
	if(file.exists(outfile.freqs)){
		stop("\nallele frequency data outfile already exists\n\n")
	}
	if(file.exists(outfile.pwp)){
		stop("\npairwise pi outfile already exists\n\n")
	}
	structure.data <- utils::read.table(infile,header=FALSE,skip=start.samples-1,stringsAsFactors=FALSE)
	sample.names <- get.sample.names(structure.data,onerowperind)
	genos <- structure.data[,start.loci:ncol(structure.data)]
	rm(structure.data)
	if(onerowperind & ncol(genos) %% 2 != 0){
		stop("\nyou have mis-specified the genotype matrix\nplease check documentation\n\n")
	}
	if(!onerowperind & nrow(genos) %% 2 != 0){
		stop("\nyou have mis-specified the genotype matrix\nplease check documentation\n\n")	
	}
	freqs <- get.freqs(genos,onerowperind,missing.datum)
	row.names(freqs) <- sample.names
	if(save.freqs){
		save(freqs,file=outfile.freqs)
	}
	pwp <- freqs2pairwisePi(freqs)
	utils::write.table(pwp,file=outfile.pwp)
	return(pwp)
}

#' Calculate pairwise pi from allele frequency data
#' 
#' \code{freqs2pairwisePi} calculates pairwise pi 
#'	between all samples included in an allele frequency 
#'  data matrix
#'
#' This function takes an allele frequency data matrix
#' for a sample of diploid individuals and returns 
#' pairwise pi between all samples, which can then be 
#' used in a BEDASSLE analysis.
#' 
#' @param freqs A \code{matrix} of allele frequencies with 
#'	one column per locus and one row per sample.
#' 	Missing data should be indicated with \code{NA}.
#' @param nLoci An integer giving the total number of loci in the 
#'  dataset. If left as \code{NULL}, it defaults to the number of 
#'	columns in the \code{freqs} argument.
#' @param quiet An Boolean (\code{TRUE}/\code{FALSE}) indicating 
#'  whether to suppress printing a progress update bar. Default is 
#'	FALSE.
#' 
#' @details This function calculates pairwise pi (the proportion of 
#'		sites at which each pair of samples differs, out of the total
#'		number of loci in the dataset) between a set of diploid 
#'		individuals from a matrix of allele frequency data.
#'		The matrix of pairwise pi that is returned can then be 
#'		used to run a \code{bedassle} analysis with 
#'		\code{\link{run.bedassle}}.
#'		
#' 		Pairwise pi is calculated as the proportion of sites at which a 
#'		pair of individuals differs out of the total number of loci in the 
#'		dataset. If it is calculated using only loci that are polymorphic 
#'		in the global dataset, it is called pi at polymorphic sites. If it 
#'		is calculated using all genotyped base-pairs, it is simply pi. Either 
#'		statistic can be converted to allelic covariance and used for running 
#'		BEDASSLE. If the \code{freqs} matrix specified consists only of 
#'		polymorphic loci, but the user wishes to calculate pi (rather than pi 
#'		at polymorphic sites), she must specify the total number of loci in the 
#'		dataset (polymorphic and invariant) using the \code{nLoci} command.
#'		
#'		Missing data is handled in a pairwise fashion in the calculation of pi 
#'		for each pair of individuals. That is, for each pair of individuals, 
#'		the function goes through each locus at which they were both genotyped 
#'		and calculates the number of sites at which they differ, then 
#'		divides that total by (\code{nLoci} - \emph{Mij}), where \emph{Mij} is 
#'		the number of loci at which either individual in the comparison is 
#'		missing data.
#'		
#'	
#'	@return This function returns the pairwise pi matrix that can then be 
#'		used to run a BEDASSLE analysis run with \code{\link{run.bedassle}}.
#'		
#' @export

freqs2pairwisePi <- function(freqs,nLoci=NULL,quiet=FALSE){
	check.freqs(freqs)
	if(is.null(nLoci)){
		nLoci <- ncol(freqs)
	}
	n <- nrow(freqs)
	pwp <- matrix(NA,n,n)
	if(!quiet){
		prog <- utils::txtProgressBar(min=0,max=n+(n*(n-1))/2,char="*",style=3)
	}
	for(i in 1:n){
		for(j in i:n){
			if(!quiet){
				utils::setTxtProgressBar(prog,(i-1)*n+j-(i*(i-1)/2))
			}
			pwp[i,j] <- calcPWP(freqs[i,],freqs[j,],nLoci)
			if(i != j){
				pwp[j,i] <- pwp[i,j]
			}
		}
	}
	return(pwp)
}

#' Create replicate data partitions for k-fold cross-validation
#' 
#' \code{freqs2xval} creates a list of data partitions to be 
#'	used in an n-replicate k-fold cross-validation analysis.
#'	divides the allele frequency dataset into partitions, and, within each, calculates 
#'	pairwise pi between samples for use in k-fold cross-validation.
#'
#' This function takes an allele frequency data matrix
#' for a sample of diploid individuals, and, for each cross-validation 
#' replicate, divides it randomly into \emph{k} partitions and 
#' calculates pairwise pi between all samples within each partition.
#' These replicate partitions can then be used in a cross-validation 
#' analysis.
#' 
#' 
#' @param freqs A \code{matrix} of allele frequencies with 
#'	one column per locus and one row per sample.
#' 	Missing data should be indicated with \code{NA}.
#' @param n.replicates An integer giving the number (n) of cross-validation 
#'  replicates to be performed. Default is 10.
#' @param n.partitions An integer giving the number (k) of partitions (folds)
#'  into which the dataset should be divided for k-fold cross-validation. 
#'	Default is 10.
#' @param prefix A character \code{vector} giving the prefix to be attached 
#'	to the output file of replicated data partitions. An underscore is 
#'	automatically added between the prefix and the file names.
#' 
#' @details For each of \code{n.replicates}, this function divides an 
#'		allele frequency dataset randomly into \emph{k} equal partitions 
#'		and, within each, calculates pairwise pi (the proportion of 
#'		sites at which each pair of samples differs, out of the total 
#'		number of loci in the dataset) between a set of diploid 
#'		individuals. For details on how pairwise pi is calculated, see 
#'		\code{\link{freqs2pairwisePi}}. Pairwise pi within each data partition 
#'		is saved as a text file, and these can then be used to run a cross-validation 
#'		analysis with \code{\link{x.validation}}. 
#'		
#' 		In each replicate, the entire dataset is divided evenly into the \emph{k}, 
#'		the number of partitions specified in \code{n.partitions}. If the total number of 
#'		loci \emph{L} cannot be evenly divided into \emph{k} partitions, the remainder of 
#'		\emph{L}/\emph{k} loci are dropped from the dataset. Within each partition, 
#'		there must be more loci than there are samples.
#'			
#'	@return This function generates and saves a list of length \code{n.replicates}, 
#'		each element of which is a list of length \code{n.partitions}. text files, each containing the 
#'		pairwise pi matrix calculated from its partition of the allele frequency data 
#'		matrix. These files can then be used for running a k-fold cross-validation 
#'		analysis using \code{\link{x.validation}}.
#'		
#' @export

freqs2xval <- function(freqs,n.replicates=10,n.partitions=10,prefix){
	check.freqs(freqs)
	nLoci <- ncol(freqs)
	loci.per.partition <- nLoci %/% n.partitions
	if(loci.per.partition < nrow(freqs)){
		stop("\nyou must have more loci per partition than samples\n")
	}
	n.to.drop <- nLoci %% n.partitions
	if(n.to.drop > 0){
		freqs <- freqs[,-sample(1:ncol(freqs),n.to.drop)]
	}
	announce.partitions(freqs,loci.per.partition,n.replicates,n.partitions,n.to.drop)
	rep.partitions <- vector("list",length=n.replicates)
	prog <- utils::txtProgressBar(min=0,max=n.replicates,char="*",style=3)
	for(n in 1:n.replicates){
		utils::setTxtProgressBar(prog,n)
		rep.partitions[[n]] <- freqs2partitions(freqs,n.partitions)
	}
	names(rep.partitions) <- paste0("rep_",1:n.replicates)
	save(rep.partitions,file=paste0(prefix,"_data_partitions.Robj"))
	return(invisible("data partitions written"))
}

freqs2partitions <- function(freqs,n.partitions){
	nLoci <- ncol(freqs)
	freqs <- freqs[,sample(1:nLoci,nLoci)]
	partitions <- vector("list",length=n.partitions)
	loci.per.partition <- nLoci %/% n.partitions
	partitions.pos <- get.partitions(nLoci,n.partitions,loci.per.partition)
	for(k in 1:n.partitions){
		partitions[[k]] <- freqs2pairwisePi(freqs[,partitions.pos[k,1]:partitions.pos[k,2]],loci.per.partition,quiet=TRUE)
	}
	names(partitions) <- paste0("fold_",1:n.partitions)
	return(partitions)
}

get.partitions <- function(nLoci,n.partitions,loci.per.chunk){
	starts <- seq(1,nLoci,by=loci.per.chunk)
	stops <- seq(starts[2]-1,nLoci,by=loci.per.chunk)
	return(cbind(starts,stops))
}

announce.partitions <- function(freqs,loci.per.chunk,n.replicates,n.partitions,n.to.drop){
	message("\nreading data:\n")
	message(sprintf("\t%s samples\n\tgenotyped at %s loci\n",nrow(freqs),ncol(freqs)))
	message("\ncreating cross-validation dataset:\n")
	if(n.to.drop > 0){
		if(n.to.drop == 1){
			message("\tdropping 1 locus to get an even number of loci per partition\n")
		} else {
			message(sprintf("\tdropping %s loci to get an even number of loci per partition\n",n.to.drop,loci.per.chunk))
		}
	}
	message(sprintf("\tmaking %s replicates, each with:\n",n.replicates))
	message(sprintf("\t\t%s partitions and %s loci per partition\n",n.partitions,loci.per.chunk))
	return(invisible("partitions announced"))
}

calcPWP <- function(ind1,ind2,L){
	nMD <- ifelse((any(is.na(ind1)) | any(is.na(ind2))),
				length(unique(which(is.na(ind1)),which(is.na(ind2)))),
				0)
	if(is.null(L)){
		L <- length(ind1) - nMD
	}
	diff.homs = sum(abs(ind1-ind2)==1,na.rm=TRUE)
	hets = sum(ind1==0.5 | ind2==0.5,na.rm=TRUE)
	return((diff.homs + hets/2)/L)
}

get.sample.names <- function(structure.data,onerowperind){
	sample.names <- structure.data[,1]
	if(!onerowperind){
		sample.names <- sample.names[seq(1,length(sample.names),by=2)]
	}
	return(sample.names)
}

get.counted.allele <- function(genos,missing.datum){
	alleles <- unique(genos)
	if(identical(alleles, missing.datum)){
		stop("\nyour dataset contains loci with all data missing. please remove and re-try.\n\n")
	}
	alleles <- alleles[!alleles==missing.datum]
	counted <- sample(alleles,1)
	return(counted)
}

get.freqs <- function(genos,onerowperind,missing.datum){
	n.loci <- ifelse(onerowperind,ncol(genos)/2,ncol(genos))
	if(onerowperind){
		freqs <- get.freqs.onerowperind(genos,n.loci,missing.datum)
	} else {
		freqs <- get.freqs.tworowperind(genos,n.loci,missing.datum)
	}
	colnames(freqs) <- NULL
	return(freqs)
}

get.freqs.onerowperind <- function(genos,n.loci,missing.datum){
	if(any(genos > 1)){
		counted.alleles <- apply(genos,2,get.counted.allele,missing.datum)
	} else {
		counted.alleles <- rep(1,n.loci)
	}
	freqs <- Reduce("cbind",
				lapply(1:n.loci,
							function(l){
								(genos[,seq(1,2*n.loci,by=2)[l]] == counted.alleles[l]) + 
								(genos[,seq(2,2*n.loci,by=2)[l]] == counted.alleles[l])
							}))
	freqs <- freqs/2
	missing.data <- Reduce("cbind",
						lapply(1:n.loci,
							function(l){
								(genos[,seq(1,2*n.loci,by=2)[l]] == missing.datum) + 
								(genos[,seq(2,2*n.loci,by=2)[l]] == missing.datum)
							}))
	freqs[missing.data==2] <- NA
	return(freqs)
}

get.freqs.tworowperind <- function(genos,n.loci,missing.datum){
	if(any(genos > 1)){
		counted.alleles <- apply(genos,2,get.counted.allele,missing.datum)
	} else {
		counted.alleles <- rep(1,n.loci)
	}
	freqs <- Reduce("cbind",
				lapply(1:n.loci,
							function(l){
								(genos[seq(1,nrow(genos),by=2),l] == counted.alleles[l]) + 
								(genos[seq(2,nrow(genos),by=2),l] == counted.alleles[l])
							}))
	freqs <- freqs/2
	missing.data <- Reduce("cbind",
						lapply(1:n.loci,
							function(l){
								(genos[seq(1,nrow(genos),by=2),l] == missing.datum) + 
								(genos[seq(2,nrow(genos),by=2),l] == missing.datum)
							}))
	freqs[missing.data==2] <- NA
	return(freqs)
}

check.pi <- function(pi){
	if(!is.numeric(pi)){
		stop("\nyou must specify a numeric value of \"pi\"\n")
	}
	if(any(pi < 0) | any(pi > 1)){
		stop("\nvalues of \"pi\" must fall between 0 and 1\n")		
	}
	return(invisible("checked pi"))
}

check.freqs <- function(freqs){
	if(class(freqs) != "matrix"){
		stop("\nthe \"freqs\" argument must be of class \"matrix\"\n")
	}
	if(any(freqs > 1,na.rm=TRUE)){
		stop("\nall values of the the \"freqs\" argument must be less than 1\n")	
	}
	if(any(freqs < 0,na.rm=TRUE)){	
		stop("\nall values of the the \"freqs\" argument must be greater than 0\n")
	}
	return(invisible("freqs arg checked"))
}