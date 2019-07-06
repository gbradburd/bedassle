functions {
	matrix spCov(int N, vector nugget,  real gamma) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(gamma,N,N); 	// make N x N matrix of gammas
		Nug_mat = diag_matrix(nugget);		// diagonal N x N matrix of nugget
		parCov += Nug_mat;					// parametric matrix is sum of gammas matrix + diagonal nugget matrix
		return parCov;	
	}
}
data {
	int<lower=1> N;		  			// number of populations or individuals
	int<lower=1> L;	    			// number of loci
	matrix[N,N] obsCov; 			// observed sample allelic covariance
}
transformed data {
	matrix[N,N] LobsCov;				// sample covariance multiplied by the n.loci
	LobsCov  = L * obsCov;
}
parameters {
	real<lower=0> gamma;				// global covariance between all samples
  	vector<lower=0>[N] nugget; 			// sample-specific variance (allele sampling noise + sample-specific drift/inbreeding)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric covariance matrix
	parCov = spCov(N, nugget, gamma);
}
model {
	nugget ~ normal(0,1);				// prior on nugget
	gamma ~ normal(0.15,10);				// prior on global covariance in allele frequencies
	LobsCov ~ wishart(L,parCov);		// wishart likelihood function
}
