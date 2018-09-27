\dontrun{

## Unless MS-Windows
# library(parallel)
# options(simLife.mc=2L)

# primary particles and secondary phase (ferrit) 
# which is already clustered and densified 
data(AL2MC_20p_k10_F2p_S)	

# simulation parameters 
opt <- list("vickers"=107,"distTol"=1.0,"Tmax"=10^11,
			"inAreafactor"=1.56, "outAreafactor"=1.43,
			"pointsConvHull"=10, "scale"=1e+06,"pl"=1L)

# lifetimes parameters
par <- list("P"=c(0.01,6,0.5,75,-15,1.5),
			"F"=c(0,0,0,105,-12,1),
			"const"=NULL)
	
nsim <- 10
stress <- as.list(seq(from=90,to=140,by=10))

cl <- NULL
## MPI/SOCKS/PSOCKS cluster object (even on Windows)
## must initialize RNG stream (rlecuyer) for reproducible results
# RNGkind("L'Ecuyer-CMRG")
# cl <- makeCluster(8)
# clusterSetRNGStream(cl)
	
# the following code may take some time
W <- woehler(S, CL=NULL, par, opt, stress=rep(stress,each=nsim),cores=1L,cl=cl)
woehlerDiagram(W, yrange=c(70,145))

## do not forget to stop cluster if used 
if(!is.null(cl)) stopCluster(cl)	
}
