\dontrun{
	
#library("parallel")
#options(mc.cores=detectCores())
#options(parallel.option="mclapply")
	
## construct a fiber system
## with volume fraction 15Vol% and second phase 2Vol%

# box and intensities
lam <- 15
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

## Spheroids of constant sizes,
## random planar within xz plane
theta <- list("size"=list(.5),
		"shape"=list("radius"=0.1),
		"orientation"=list("kappa"=10))

## fibers
S <- simCylinderSystem(theta,lam,size="const",
		orientation="rbetaiso",box=box,pl=101,label="P")

## secondary phase: particles as spheres
F <- simSphereSystem(list(0.075),10, rdist="const", box=box, pl=101, label="F")

## apply RSA
S2 <- rsa(S,box,F=F,pl=101,verbose=TRUE)

## construct clusters
param <- list("r"=0.35)
CLUST <- simCluster(S2, param, 0.2, box,verbose=TRUE)

## densify
ctrl <- list(threshold.stop=0.01,max.call=5000)
RET <- densifyCluster(S2, CLUST, box, ctrl, weight=100)

## whole densified system  
SF <- RET$S

## clustered areas
CLF <- RET$cluster

## 'rgl' plot of cylinders
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# cylinders3d(SF, box, col=cols)	

#path <- "/"
#save(SF,file=paste0(path,"test_S.rda"))
#save(CLF,file=paste0(path,"test_CL.rda"))
	
}