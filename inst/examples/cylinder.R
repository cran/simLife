\dontrun{

## Not on MS-Windows	
#library("parallel")
#options(mc.cores=detectCores())
#options(parallel.option="mclapply")

lam <- 35
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

## Spheroids of constant sizes
theta <- list("size"=list(.5),
		"shape"=list("radius"=0.1),
		"orientation"=list("kappa"=0.1))

S <- simCylinderSystem(theta,lam,size="const",
		orientation="rbetaiso",box=box,pl=101,label="P")

## secondary phase: particles as spheres
F <- simSphereSystem(list(0.075),5, rdist="const", box=box, pl=101, label="F")

## apply RSA
S2 <- rsa(S,box,F,pl=101,verbose=TRUE)

## draw some projection
#require("rgl")
#id <- c(1,5,9,32,10)
#cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
#cylinders3d(S2[id], box, col=cols)	
#P <- getCylinderProjection(S2[id], B=c(0,1,0,1,1), draw=TRUE, conv = TRUE, np=20)

## construct clusters
param <- list("r"=0.35)
CLUST <- simCluster(S2, param, 0.1, box,verbose=TRUE)
	
## densify
ctrl <- list(threshold.stop=0.01,max.call=5000,verbose=FALSE)
RET <- densifyCluster(S2, CLUST, box, ctrl, weight=10, cl = NULL)	
G <- RET$cluster

## draw original cluster
#open3d()
#lapply(CLUST,function(x) cylinders3d(S2[x$id],box=box,col=cols))	
#X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
#invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))

## draw densified cluster
#open3d()
#invisible(lapply(G,function(x) { cylinders3d(x,box=box,col=cols) }))
#invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))

}
