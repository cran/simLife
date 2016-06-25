\dontrun{

## Simulate and densify particle system
## densify and include a secondary phase
	
# Not on MS-Windows 	
#library("parallel")
#options(mc.cores=detectCores())
#options(parallel.option="mclapply")
 
# simulation box either this kind of list
# or use object spatstat::box3 
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

# parameter for spheroids simulation 
theta <- list("size"=list(0.25),"shape"=list("s"=0.25), "orientation"=list("kappa"=1))

# for ease of use: constant size particles
S <- simParticle(theta,lam=15,box=box,mu=c(0,1,0))

param <- list("size"=list(0.075), "shape"=list("s"=0.75),"orientation"=list("kappa"=1))
F  <- simFerrit(param,lam=2,box=box,mu=c(0,1,0)) 

# apply RSA, this may take some 
RSA <- rsa(S, box, F=F, pl = 10)

## show 3d spheroids (non-overlapping)
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# spheroids3d(RSA[1:length(S)], box, col=cols)	
# spheroids3d(RSA[length(S):(length(S)+length(F))], box, col="gray")

# construct clusters
param <- list("r"=0.35)
CLUST <- simCluster(RSA, param, 0.1, box, verbose=TRUE)
cat("cluster length: ",length(CLUST),"\n")

# X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
# invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# invisible(lapply(CLUST, function(x) spheroids3d(RSA[x$id],box,col=cols)))

# controls for 'GenSA'
ctrl <- list(threshold.stop=0.01,max.call=10000,verbose=FALSE)
# densify
RET <- densifyCluster(RSA, CLUST, box, ctrl, weight=100,cl = NULL)	

S <- RET$S
CL <- RET$cluster

## show densified clusters 
# X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
# invisible(lapply(CLUST, function(x)
# rgl::spheres3d(x=X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
# invisible(lapply(CL, function(x) { spheroids3d(x, box, col=cols) }))

## save results	
# path <- "/"
# save(S,file=paste0(path,"test_particle_S.rda"))
# save(CL,file=paste0(path,"test_particle_CL.rda"))
	
}	
