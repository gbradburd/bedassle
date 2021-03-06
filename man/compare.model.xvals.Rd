% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.comparison.R
\name{compare.model.xvals}
\alias{compare.model.xvals}
\title{Compare cross-validated BEDASSLE models}
\usage{
compare.model.xvals(xval.files, n.predictors, mod.cols = NULL)
}
\arguments{
\item{xval.files}{A \code{vector} of filenames (each in quotes, 
with the full file path), each element of which points to the 
cross-validation results file output by a call to 
\code{\link{x.validation}} for one of the models the user wants 
to compare.}

\item{n.predictors}{A \code{vector} of \code{integer} values giving 
the number of predictors included in each of the BEDASSLE models 
tested with cross-validation. This argument should be the same as 
the length of the \code{vector} of filenames specified in 
\code{xval.files}.}

\item{mod.cols}{An \code{vector} of colors to be used in plotting 
the results of the different BEDASSLE models. This argument should be 
the same length as the \code{vector} of filenames specified in 
\code{xval.files}. If \code{NULL}, all output will be plotted in 
black.}
}
\value{
This function creates a plot showing the standardized predictive 
	accuracies of the different models evaluated using cross-validation. 
	A higher predictive accuracy means the model is able to describe the 
	data better. The plot shows mean standardized predictive accuracy for 
	each replicate, the mean across replicates, and the 95% confidence 
	letter over a pair of models indicates that the predictive accuracies 
	interval. A shared of those models are not significantly different 
	using a paired, two-tailed t-test with a significance level of 0.05.

	The predictive accuracies are standardized by, for each partition in 
	each replicate, subtracting the highest predictive accuracy from the  
	predictive accuracies of the other models. Therefore, a predictive 
	accuracy of 0 is the best score.

	The function returns a matrix giving the significance of the difference 
	between the predictive accuracies of the different models evaluted using 
	k-fold cross-validation in the function \code{\link{x.validation}}
	analysis.

Note that \code{n.predictors} is the number of predictor variables in the 
relevant model. The null model, with neither geographic nor environmental 
distance, has 0 predictors, a model with just geographic distance has 1 
predictor, as does a model with just a single environmental distance 
variable, and a model with both geographic distance and an environmental 
distance predictor has 2 predictors.
}
\description{
\code{compare.model.xvals} compares cross-validated BEDASSLE models
}
\details{
This function visually and statistically compares BEDASSLE models that 
have been evaluated using cross-validation.
}
