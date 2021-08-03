functions {
	matrix getEdist(int N, int nE, vector aE, matrix[] envDist){
		matrix[N,N] Edist;
		Edist = rep_matrix(0,N,N);
		for(e in 1:nE){
			Edist += aE[e] * envDist[e];
		}
		return Edist;
	}
	matrix spCov(int N, real a0, real a2,  real nugget, vector aE, matrix[] envDist,  int nE) {
		matrix[N,N] parCov;
		matrix[N,N] Edist;
		Edist = getEdist(N,nE,aE,envDist);
		for(i in 1:(N-1)){
			for(j in (i+1):N){
				parCov[i,j] = Edist[i,j] == 0 ? a0 : a0 * exp(-(Edist[i,j])^a2); 	// calculate ecological effects on parametric covariances, upper triangle and diagonal
				parCov[j,i] = parCov[i,j]; 	// reflect below the diagonal
			}
		}
		for(i in 1:N) parCov[i,i] = a0 + nugget;
		return parCov;	
	}
}
data {
	int<lower=1> N;		  				// number of populations or individuals
	int<lower=1> L;	    				// number of loci
	int<lower=1> nE;	    				// number of ecological variables
	matrix[N,N] obsCov; 				// observed projected covariance
	matrix<lower=0>[N,N] envDist[nE]; 					// array of pairwise ecological distances
}
transformed data {
	matrix[N,N] LobsCov;				// sample covariance multiplied by n.loci
	LobsCov  = L * obsCov;
}
parameters {
	real<lower=0> alpha0;				// sill of the parametric covariance 
	vector[nE] alphaE;			// vector of effects of ecological distance in the parametric covariance
	real<lower=0,upper=2>  alpha2;		// exponential slope parameter in the parametric covariance
  	real<lower=0> nugget; 			// variance (allele sampling error + sample-specific drift)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric covariance matrix
	parCov = spCov(N, alpha0,  alpha2, nugget, alphaE, envDist, nE);
}
model {
	alpha0 ~ normal(0,1);
	alphaE ~ normal(0,1);
	alpha2 ~ uniform(0,2);
	nugget ~ std_normal();
	LobsCov ~ wishart(L,parCov);
}
