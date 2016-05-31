## Simulate a fiber system and apply RSA to
## obtain a non-overlapping configuration
	
lam <- 5
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

## Spheroids of constant sizes
theta <- list("size"=list(0.95),
			  "shape"=list("r"=0.05),
			  "orientation"=list("kappa"=1))

## primary phase: fibers
S <- simCylinderSystem(theta,lam,size="const",shape="radius",
		orientation="rbetaiso",box=box,pl=101,label="P")

## secondary phase: particles as spheres
F <- simSphereSystem(list(0.075),5, rdist="const", box=box, pl=101, label="F")
## apply RSA
S2 <- rsa(S,box,F=F,pl=101,verbose=TRUE)

#require("rgl")
#cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
#open3d()
#cylinders3d(S2, box, col=cols)
