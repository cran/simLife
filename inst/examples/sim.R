## Simulation of individual defect times	
data(AL2MC_20p_k10_F2p_S)
data(AL2MC_20p_k10_F2p_CL)

## generate individual failure times
opt <- list("vickers"=107,"distTol"=0.001,
			"inAreafactor"=1.56, "outAreafactor"=1.43, 
			"pointsConvHull"=10, "scale"=1e+06,"pl"=0)	

## simulation parameter	
par <- list("P"=c(0.01,10^12,10,108,-14,0.01),
			"F"=c(0.01,10^10,10,100,-13,0.01))	

## simulate times	
CLT <- simTimes(S,par,vickers=opt$vickers,stress=125)

## times
T <- unlist(sapply(CLT,`[[`,"T"))
V <- unlist(sapply(CLT,`[[`,"V"))
U <- unlist(sapply(CLT,`[[`,"U"))

## show estimated densities
showDensity(list("Delamination"=log10(V),"Crack"=log10(U),"Time"=log10(T)))
