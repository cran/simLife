\dontrun{

## Simulate a defect accumulation
## of a particle system
library(simLife)

## primary particles (as prolate spheroids)
## and secondary phase (as spheres)
data(AL2MC_20p_k10_F2p_S)

## which is clustered and densified according to CL 
## additional predefined clustered regions
data(AL2MC_20p_k10_F2p_CL)

## the box is stored together with the geometry system 
box <- attr(S,"box")

## distTol=0.25, so use 25% of accumulation distance
opt <- list("vickers"=107, "distTol"=0.25, "Tmax"=10^12,
		"inAreafactor"=1.56, "outAreafactor"=1.43, 
		"pointsConvHull"=10, "scale"=1e+06,"pl"=101)

## constants 'const' are set to the following internally
## and can partially be overwritten if needed
# const <- list("Em"=68.9,"Ef"=400,"nc"=28.2,"nu"=0.33,
#		    "pf"=0.138,"nE"=NULL,"sigref"=276,"Vref"=5891)

par <- list("P"=c(0.01,5,0.5,50,-11,3),
		"F"=c(0,0,0,105,-12,1),"const"=NULL)

# stress amplitude applied	
stress <- 110
## generate individual (particles') failure times
CLT <- simTimes(S,par,vickers=opt$vickers,stress=stress)

## generated random failure times
T <- unlist(sapply(CLT,`[[`,"T"))
V <- unlist(sapply(CLT,`[[`,"V"))
U <- unlist(sapply(CLT,`[[`,"U"))

dev.new()
showDensity(list("Debonding"=log10(V),"Fracture"=log10(U),"Failure time"=log10(T)),xlim=c(-2,15))

# Area max
(inA <- areaCrit(opt$vickers, stress=stress, factor=opt$inAreafactor, scale=opt$scale))
(outA <- areaCrit(opt$vickers, stress=stress, factor=opt$outAreafactor, scale=opt$scale))

## run accumulation
RET <- simDefect(stress,S,CLT,opt)
length(RET)

#### alternatively, includes simulating times by 'simTimes'
# SIM <- simFracture(stress,S,opt,par,last.defect=FALSE,CL=CL)	
# SIM$cl_info

## some simple analysis of accumulation paths until failure,
## that is, while the critical area is exceeded
LR <- RET[[length(RET)]]
isInCluster <- any(unlist(lapply(CL,function(x,y)
					any(y$id %in% x$id) , y=LR)))
cat("Broken cluster: ", isInCluster,"\t Ferrit: ",
	any("F" %in% LR$label),"\t Acc.size",length(LR$id),"\n")
	
## select only clusters of size larger than 'msize'	
msize <- 1
id <- sapply(RET,function(x) ifelse(length(x$id)>msize,TRUE,FALSE))
cat("Number of defect projections in last cluster: ",length(RET[[length(RET)]]$id),"\n")

## draw all accumulation paths until failure
dev.new()
L <- plotDefectAcc(RET,last.path=FALSE)
## draw last accumulation path until failure
dev.new()
L <- plotDefectAcc(RET,last.path=TRUE)

## 3D visualization of final defect projection area
#library(rgl)
#library(unfoldr)
#qid <- LR$id
#open3d()
#spheroids3d(S[qid],box=box, col=c("#0000FF","#00FF00","#FF0000","#FF00FF"))
#
### drawing only last cluster leading to failure
#drawDefectProjections(S,list(LR))
}
