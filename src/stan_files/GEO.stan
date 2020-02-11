functions {
	matrix spCov(int N, real a0, real a2,  vector nugget, matrix D, real aD) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		Nug_mat = diag_matrix(nugget); 					// set up N x N diagonal matrix of nuggets
		for(i in 1:(N-1)){
			for(j in (i+1):N){
				parCov[i,j] = D[i,j] == 0 ? a0 : a0 * exp(-(aD * D[i,j])^a2); 	// calculate geographic effects on parametric covariances, upper triangle and diagonal;
				parCov[j,i] = parCov[i,j]; 		// reflect below the diagonal
			}
		}
		for(i in 1:N) parCov[i,i] = a0;
		parCov += Nug_mat;					// parametric covariance is is prior matrix plus nugget matrix
		return parCov;	
	}
}
data {
	int<lower=1> N;		  			// number of populations or individuals
	int<lower=1> L;	    			// number of loci
	matrix[N,N] obsCov; 			// observed projected covariance
	matrix<lower=0>[N,N] geoDist;	// matrix of pairwise geographic distance 
}
transformed data {
	matrix[N,N] LobsCov;				// the sample covariance multiplied by n.loci
	LobsCov  = L * obsCov;
}
parameters {
	real<lower=0> alpha0;				// sill of the parametric covariance 
	real<lower=0> alphaD;				// effect of geographic distance in the parametric covariance 
	real<lower=0, upper=2>  alpha2;		// exponential slope parameter in the parametric covariance 
  	vector<lower=0>[N] nugget; 			// sample-specific variance (allele sampling noise + sample-specific drift/inbreeding)
}
transformed parameters {
	matrix[N,N] parCov;				// this specifies the parametric covariance matrix
	parCov = spCov(N, alpha0,  alpha2, nugget, geoDist, alphaD);
}
model {
	alpha0 ~ normal(0,1);			// prior on alpha0
	alphaD ~ normal(0,1);			// prior on alphaD
	alpha2 ~ uniform(0,2);			// prior on alpha2
	nugget ~ normal(0,1);			// prior on nugget
	LobsCov ~ wishart(L,parCov);	// likelihood function
}
