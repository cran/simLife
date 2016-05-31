## Simulate a particle system by spheres
## simulate densified clusters 
## and generate a non-overlapping configuration of spheres by RSA

## MS-Windows only
#library("parallel")
#options(mc.cores=detectCores())
#options(parallel.option="mclapply")	
	
lam <-15
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))
theta <- list("r"=0.1)
S <- simSphereSystem(theta,lam, rdist="const", box=box,pl=101,label="P")
# rsa
S2 <- rsa(S,box,pl=101,verbose=TRUE)
# project some spheres
id <- c(1,5,9,32,10)
# get matrix of border points of sphere projections
P <- getSphereProjection(S2[id],draw=FALSE)
	
# densify and construct clusters
# with radius 0.35 based on hardcore particle system
param <- list("r"=0.35)
CLUST <- simCluster(S2, param, 0.1, box, verbose=TRUE)
	
# densify
ctrl <- list(threshold.stop=0.01,max.call=5000,verbose=FALSE)
RET <- densifyCluster(S2, CLUST, box, ctrl, weight=10,cl = NULL)

## Uncomment the following lines for 3d visualization

## get the densified cluster
#G <- RET$cluster

#require("rgl")
#cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
#drawSpheres <- function(S, box, ...) {
#	X <- do.call(rbind,lapply(S,function(x) c(x$center,x$r) ))
#	rgl::spheres3d(X[,1:3],radius=X[,4],...)
#	
#	x <- box$xrange[2];	
#	y <- box$yrange[2];	
#	z <- box$zrange[2]
#	c3d.origin <- translate3d(scale3d(cube3d(col="darkgray", alpha=0.1),x/2,y/2,z/2),x/2,y/2,z/2)
#	shade3d(c3d.origin)
#	
#	axes3d(edges = "bbox",labels=TRUE,tick=FALSE,box=TRUE,nticks=0,
#			expand=1.0,xlen=0,xunit=0,ylen=0,yunit=0,zlen=0,zunit=0)
#}	

## draw original clusters
#X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
#
#open3d()
#invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
#lapply(CLUST,function(x) drawSpheres(S2[x$id],box=box,col=cols))	
#
## draw densified clusters and their projections
#open3d()
#invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
#invisible(lapply(G,function(x) {
#	drawSpheres(x,box=box,col=cols)
#	invisible(getSphereProjection(x,draw=TRUE))	
#}))
#
