functions {
	matrix[] sqArray(int N, int nE, matrix[] X){
		matrix[N,N] sqX[nE];
		for(e in 1:nE){
			for(i in 1:N){
				for(j in i:N){
					sqX[e,i,j] = X[e,i,j]^2;
					if(i != j ){
						sqX[e,j,i] = X[e,i,j]^2;
					}
				}
			}
		}
		return sqX;
	}
	matrix getEdist(int N, int nE, vector aE, matrix[] sqE){
		matrix[N,N] Edist;
		Edist = rep_matrix(0,N,N);
		for(e in 1:nE){
			Edist += aE[e]^2 * sqE[e];
		}
		return Edist;
	}
	matrix spCov(int N, real a0, real a2,  vector nugget, vector aE, matrix[] sqE,  int nE) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		matrix[N,N] Edist;
		Nug_mat = diag_matrix(nugget); 					// set up N x N diagonal matrix of nuggets
		Edist = getEdist(N,nE,aE,sqE);
		for(i in 1:(N-1)){
			for(j in (i+1):N){
				parCov[i,j] = Edist[i,j] == 0 ? a0 : a0 * exp(-sqrt(Edist[i,j])^a2); 	// calculate ecological effects on parametric covariances, upper triangle and diagonal
				parCov[j,i] = parCov[i,j]; 	// reflect below the diagonal
			}
		}
		for(i in 1:N) parCov[i,i] = a0;
		parCov += Nug_mat;				// parametric covariance is is prior matrix plus nugget matrix
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
	matrix[N,N] sqE[nE];				// array of squared ecological distance matrices
	LobsCov  = L * obsCov;
	sqE = sqArray(N, nE, envDist);
}
parameters {
	real<lower=0> alpha0;				// sill of the parametric covariance 
	vector<lower=0>[nE] alphaE;			// vector of effects of ecological distance in the parametric covariance
	real<lower=0, upper=2>  alpha2;		// exponential slope parameter in the parametric covariance
  	vector<lower=0>[N] nugget; 			// sample-specific variance (allele sampling error + sample-specific drift)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric covariance matrix
	parCov = spCov(N, alpha0,  alpha2, nugget, alphaE, sqE, nE);
}
model {
	alpha0 ~ normal(0,1);				// prior on alpha0
	alphaE ~ normal(0,1);				// prior on alphaE
	alpha2 ~ uniform(0,2);				// prior on alpha2
	nugget ~ normal(0,1);				// prior on nugget
	LobsCov ~ wishart(L,parCov);		// likelihood function
}
