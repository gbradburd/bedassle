functions {
	matrix spCov(int N, real nugget,  real a0) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(a0,N,N); 	// make N x N matrix of alpha0
		for(i in 1:N) parCov[i,i] += nugget;
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
  	real<lower=0> nugget; 			// sample-specific variance (allele sampling noise + sample-specific drift/inbreeding)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric covariance matrix
	parCov = spCov(N, nugget, alpha0);
}
model {
	alpha0 ~ normal(0,1);
	nugget ~ std_normal();
	LobsCov ~ wishart(L,parCov);
}
