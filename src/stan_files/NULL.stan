functions {
	matrix spCov(int N, vector nugget,  real alpha0) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(alpha0,N,N); 	// make N x N matrix of alpha0
		Nug_mat = diag_matrix(nugget);		// diagonal N x N matrix of nugget
		parCov += Nug_mat;					// parametric matrix is sum of alpha0 matrix + diagonal nugget matrix
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
	real<lower=0> alpha0;									// sill of the parametric covariance 
  	vector<lower=0>[N] nugget; 			// sample-specific variance (allele sampling noise + sample-specific drift/inbreeding)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric covariance matrix
	parCov = spCov(N, nugget, alpha0);
}
model {
	nugget ~ normal(0,1);				// prior on nugget
	alpha0 ~ normal(0,1);				// prior on alpha0
	LobsCov ~ wishart(L,parCov);		// wishart likelihood function
}
