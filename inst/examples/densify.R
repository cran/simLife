\dontrun{

## Simulate and densify particle system
## densify and include a secondary phase

library(unfoldr)

## Not on MS-Windows	
# library("parallel")
# options(mc.cores=2L)
 
# simulation box either this kind of list
# or use object spatstat::box3 
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

# parameter for spheroids simulation 
theta <- list("size"=list(0.25),"shape"=list("s"=0.25), "orientation"=list("kappa"=1))

# for ease of use: constant size particles
S <- simPoissonSystem(theta,lam=15,size="const",type="prolate",
			orientation="rbetaiso",box=box,mu=c(0,1,0),pl=1,label="P")
	
## 2nd. phase (ferrit)	
param <- list("size"=list(0.075), "shape"=list("s"=0.75))
F <- simPoissonSystem(param,lam=2,size="const",type="prolate",
		box=box,mu=c(0,1,0),pl=1,label="P")

# apply RSA, this may take some 
RSA <- rsa(S,F,verbose=TRUE)

## show 3D spheroids (non-overlapping)
# library(rgl)
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# spheroids3d(RSA[1:length(S)], box, col=cols)	
# spheroids3d(RSA[length(S):(length(S)+length(F))], box, col="gray")

## construct clusters
CL <- simPoissonSystem(list("size"=list(0.35)), lam=0.1, size="const",
			type="spheres", box=box, pl=1, label="F")

CLUST <- simCluster(RSA, CL, verbose=TRUE)
cat("cluster length: ",length(CLUST),"\n")

## show cluster regions
# library(rgl)
# open3d()
# X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
# invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# invisible(lapply(CLUST, function(x) spheroids3d(RSA[x$id],box,col=cols)))

# some controls for 'GenSA'
ctrl <- list(threshold.stop=0.01,max.call=10000)
# densify region to clustered particles
RET <- densifyCluster(RSA, CLUST, ctrl, weight=100,cl = NULL)	

S <- RET$S
CL <- RET$cluster

## show densified clusters 
# library(rgl)
# open3d()
# X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
# invisible(lapply(CLUST, function(x)
# rgl::spheres3d(x=X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
# invisible(lapply(CL, function(x) { spheroids3d(x, box, col=cols) }))
	
}	
