\dontrun{
	
## Simulate a fiber system (with cylinders),
## apply RSA to get a non-overlapping (hardcore) configuration
	
library(unfoldr)

lam <- 5
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

## Spheroids of constant sizes
theta <- list("size"=list(0.95),
			  "shape"=list("r"=0.05),
			  "orientation"=list("kappa"=1))

## primary phase: fibers
S <- simPoissonSystem(theta,lam,size="const",shape="const",
		type="cylinders",box=box,pl=1,label="P")

## secondary phase: particles as spheres
F <- simPoissonSystem(theta=list("size"=list(0.075)),
		type="spheres",lam=5,size="const",box=box, pl=1, label="F")

## apply RSA
S2 <- rsa(S,F,verbose=TRUE)

## Ferrit (2nd. phase shown as gray spheres)
# library(rgl)
# open3d()
# cylinders3d(S2, box, col=c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF"))
}
