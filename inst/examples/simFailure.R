## Load a particle system
## and simulate its defect accumulation

## primary particles and secondary phase (ferrit)
data(AL2MC_20p_k10_F2p_S)
## which is clustered and densified according to CL 
## additional predefined clustered regions
data(AL2MC_20p_k10_F2p_CL)

## the box is stored together with the geometry system 
box <- attr(S,"box")

## distTol=0.25, so use 25% of accumulation distance
opt <- list("vickers"=107,"distTol"=0.25,
		"inAreafactor"=1.56, "outAreafactor"=1.43, 
		"pointsConvHull"=10, "scale"=1e+06,"pl"=101)

par <- list("P"=c(0.01,10^12,10,105,-13,0.01),
			"F"=c(0.01,10^12,10,98,-12,0.01))

# stress amplitude applied	
stress <- 110
## generate individual (particles') failure times
CLT <- simTimes(S,par,vickers=opt$vickers,stress=stress)
## run accumulation
RET <- simDefect(stress,S,CLT,opt)

#### alternatively run
#SIM <- simFracture(stress,S,opt,par,last.defect=FALSE,CL=CL)	
#SIM$cl_info
####

## some simple analysis
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

## plot last 
## get particle id numbers of last cluster
# qid <- LR$id
# spheroids3d(S[qid],box=box, col=c("#0000FF","#00FF00","#FF0000","#FF00FF"))
## drawing only last cluster leading to failure
# drawDefectProjections(S,list(LR))
