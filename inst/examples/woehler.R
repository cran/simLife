\dontrun{

## For use with parallel package
# library("parallel")
# options(mc.cores=detectCores())
# options(parallel.option="mclapply")

# primary particles and secondary phase (ferrit) 
# which is already clustered and densified 
data(AL2MC_20p_k10_F2p_S)	

# simulation parameters 
opt <- list("vickers"=107,"distTol"=0.5,
			"inAreafactor"=1.56, "outAreafactor"=1.43, 
			"pointsConvHull"=10, "scale"=1e+06,"pl"=0)

# lifetimes parameters
par <- list("P"=c(0.01,10^12,10,105,-13,0.05),
			"F"=c(0.01,10^12,10,98,-12,0.05))
	
nsim <- 15
stress <- as.list(seq(from=90,to=140,by=10))

## For use of SOCKS or MPI cluster do the following
# library("snow")
# cl <- makeCluster(6)
## or	
# cl <- makeSOCKcluster(6)
# clusterEvalQ(cl,library("simLife"))

cl <- NULL	
# the following code may take some time
W <- woehler(S, CL=NULL, par, opt, stress=rep(stress,each=nsim))

woehlerDiagram(W, yrange=c(65,145), cl=cl)

## do not forget to stop cluster if used 
# stopCluster(cl)
	
}
