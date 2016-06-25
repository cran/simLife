\dontrun{
	
## random planar (in xz plane) orientated fibers
data(AL2MC_15p_k10_F2p_S)
data(AL2MC_15p_k10_F2p_CL)

## simulation box	
box <- attr(SF,"box")

## adjust according to material constants
opt <- list("vickers"=107,
			"distTol"=0.75,
			"inAreafactor"=1.56, "outAreafactor"=1.43, 
			"pointsConvHull"=20, "scale"=1e+06, "pl"=101)

par <- list("P"=c(0.01,10^12,10,110,-13,0.01),
			"F"=c(0.01,10^10,10,125,-15,0.01))	

## fiber system in 3d (with 'rgl')	
# require("rgl")
# id <- 1:1000
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# cylinders3d(SF[id], box, col=cols)	

stress <- 100
# generate individual (particles') failure times
CLT <- simTimes(SF,par,vickers=opt$vickers,stress=stress)

## generated random failure times
T <- unlist(sapply(CLT,`[[`,"T"))
V <- unlist(sapply(CLT,`[[`,"V"))
U <- unlist(sapply(CLT,`[[`,"U"))
## plot estimated densities
showDensity(list("Delamination"=log10(V),"Crack"=log10(U),"Time"=log10(T)))

## accumulation: may take some time
RET <- simDefect(stress,SF,CLT,opt)

## draw all accumulation paths until failure
dev.new()
L <- plotDefectAcc(RET,last.path=FALSE)
## draw last accumulation path until failure
dev.new()
L <- plotDefectAcc(RET,last.path=TRUE)

## largest area defect (inner/outer)
#require("rgl")
#eid <- which.max(sapply(RET,function(x) max(x$A)))
#qid <- unlist(lapply(RET[eid],function(x) x$id))
#LR <- RET[eid]
#open3d()
#cylinders3d(SF[qid],box=box, col=c("#0000FF","#00FF00","#FF0000","#FF00FF"))
#drawDefectProjections(SF,LR)

## last defect leading to failure
## might be the same as the one before
#require("rgl")
#LR <- RET[[length(RET)]]
#qid <- LR$id
#open3d()
#cylinders3d(SF[qid],box=box, col=c("#0000FF","#00FF00","#FF0000","#FF00FF"))
#drawDefectProjections(SF,list(LR))

# check projected area to last defect cluster convex hull area
P <- getCylinderProjection(SF[qid], B=LR$B, draw=FALSE, conv=TRUE)
c("projected"=P$convH[[1]], "defect"=max(LR$A))

}
