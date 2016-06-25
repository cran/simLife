###############################################################################
# Author:  M. Baaske
# Date:	   03.11.2015	
# File:    times.R: 
# 
# Comment: Simulation of individual failure times
# 
###############################################################################

# param[1]= probability for already cracked fibers
# param[2]= scale factor
# param[3]= shape -> basically controls for the scattering of the (positive) log-times, 
#			roughly proportional to 1/shape^2
# param[4]= shift of the log-time
# param[5]= slope
# param[6]= stddev of the random error of the log-time


#' Defect failure times
#' 
#' Simulation of individual defect failure times. For a secondary phase only 
#' the defect type "delamination" is considered.   
#' 
#' @param S				   geometry objects system
#' @param stress		   the applied stress level  
#' @param vickers   	   Vickers hardness
#' @param param	 		   simulation parameter list of parameter vectors for both phases
#' @param parallel.option  optional, in case of \code{mclapply} the function 
#' 							 \code{mclapply} is used from the parallel package
#' 
#' @return  a list with the following elements:
#' 			\itemize{
#' 				\item{id}{ id of particle }
#' 				\item{U}{ crack failure time }
#' 				\item{V}{ delamination failure time }
#' 				\item{T}{ the minimum of both failure times}
#' 				\item{B}{ failure type, either 0  for particle crack or 1 for particle delamination }
#' 				\item{A}{ projection area set to zero for later use }
#' 				\item{label}{ either \code{label="P"} for primary phase or \code{label="F"} 
#' 							  for secondary phase }
#' 			}
#' @author	Felix Ballani, Markus Baaske 
simCrackTime <- function(S,stress,vickers,param,
					parallel.option = getOption("parallel.option","lapply")) UseMethod("simCrackTime",S)

#' @rdname simCrackTime
#' @method simCrackTime oblate
#' @S3method simCrackTime oblate
simCrackTime.oblate <- function(S,stress,vickers,param,
						parallel.option = getOption("parallel.option","lapply")) 
{ simCrackTime.prolate(S,stress,vickers,param) }
# because the lengths c and a are already switched in E$ab at generation at C-level

#' @rdname simCrackTime
#' @method simCrackTime prolate
#' @S3method simCrackTime prolate
simCrackTime.prolate <- function(S,stress,vickers,param,
		parallel.option = getOption("parallel.option","lapply")) {
  simT <- function(E) {
	uv <- numeric(2) # [u,v]
	label <- attr(E,"label")
	if(label == "P") {
		uv[1] <- getCrackTime(E$angles[1],E$ab[1],E$ab[2],stress,vickers,param$P)
		uv[2] <- getDelamTime(E,stress,param$P)
		list("id"=E$id,"U"=uv[1],"V"=uv[2],
			 "T"=min(uv[1],uv[2]),"B"=ifelse(uv[1]<uv[2],0,1),"A"=0,"label"=label)
	} else {
		## always delamination for ferrit phase
		uv[2] <- getDelamTime(E,stress,param$F)
		list("id"=E$id,"U"=Inf,"V"=uv[2],
			 "T"=uv[2],"B"=1,"A"=0,"label"=label)
	}		
 }
 if(parallel.option=="mclapply" && !requireNamespace("parallel", quietly=TRUE))
  stop("package 'parallel' is requiered to run this function in parallel mode.")
 fun <- get(parallel.option)
 fun(S,simT) 
}

#' @rdname simCrackTime
#' @method simCrackTime cylinder
#' @S3method simCrackTime cylinder
simCrackTime.cylinder <- function(S,stress,vickers,param,
		parallel.option = getOption("parallel.option","lapply")) {
 simT <- function(E) {
	uv <- numeric(2) # [u,v]
	label <- attr(E,"label")
	if(label == "P") {
		uv[1] <- getCrackTime(E$angles[1],E$r,0.5*E$length,stress,vickers,param$P)
		uv[2] <- getDelamTime(E,stress,param$P)
		list("id"=E$id,"U"=uv[1],"V"=uv[2],
				"T"=min(uv[1],uv[2]),"B"=ifelse(uv[1]<uv[2],0,1),"A"=0,"label"=label)
	} else {
		## always delamination for ferrit phase
		uv <- getDelamTime(E,stress,param$F)
		list("id"=E$id,"U"=Inf,"V"=uv,
				"T"=uv,"B"=1,"A"=0,"label"=label)
	}	
  }
  
  if(parallel.option=="mclapply" && !requireNamespace("parallel", quietly=TRUE))
	  stop("package 'parallel' is requiered to run this function in parallel mode.")
  fun <- get(parallel.option)
  fun(S,simT) 
}

#' @rdname simCrackTime
#' @method simCrackTime sphere
#' @S3method simCrackTime sphere
simCrackTime.sphere <- function(S,stress,vickers,param,
		parallel.option = getOption("parallel.option","lapply")) {	
 simT <- function(E) {
	## always delamination for spheres
	label <- attr(E,"label")
	uv <- if(label == "P")	getDelamTime(E,stress,param$P)
		  else getDelamTime(E,stress,param$F)
	list("id"=E$id,"U"=Inf,"V"=uv,"T"=uv,"B"=1,"A"=0,"label"=label)		
 }	
 if(parallel.option=="mclapply" && !requireNamespace("parallel", quietly=TRUE))
	 stop("package 'parallel' is requiered to run this function in parallel mode.")
 fun <- get(parallel.option)
 fun(S,simT) 
}


## Disc: disc projection
#' Generate individual fracture time
#' 
#' Generate individual defect time for particle fracture
#' 
#' The particle fracture (crack) is assumed to happen orthogonal to the major axis direction along 
#' the maximum minor axis length. Thus the projection area can be easily computed for the purpose
#' of defect accumulation. The parameter set is made up of six parameters. Here only the second and 
#' third parameters are used to simulate the defect \code{crack} times. The failure times follow a 
#' Weibull distribution with scale parameter \eqn{p2*a^2/(b*\sigma*cos\theta*Hv)} and shape parameter
#' \eqn{p3} where \eqn{\sigma} denotes the stress, \eqn{a} the minor axis length and
#' \eqn{Hv} the Vickers hardness. The angle \eqn{\theta} is measured between the rotational axis 
#' and the axis of main load direction. In this way we account for the orientation of particles (spheroids)
#' when generating fracture times dependent on their tendency to be more or less oriented towards the 
#' main load direction.
#' 
#' @param theta		colatitude angle
#' @param a			axis length (axis orthogonal to rotational axis)
#' @param b			rotational axis length
#' @param stress	stress level
#' @param vickers   Vickers hardness
#' @param param	    simulation parameter set
#' 
#' @return  numeric, the individual fracture time
#' @author	Felix Ballani  
#' @seealso \code{\link{getDelamTime}}
getCrackTime<-function(theta,a,b,stress,vickers,param){
	theta <- ifelse(theta<pi/2,theta,theta-pi/2)
	U <- ifelse(runif(1)<param[1],1+0.0001*runif(1),
			rweibull(1,scale=param[2]*(a*1000)^2/(b*1000)/stress/cos(theta)/vickers,
					shape=param[3]))
	
	# if an already cracked fiber is critical then at least 1 cycle is needed until failure,
	#	the bit randomness is only to avoid ties
	
	return(U)
}

#getCrackTime2<-function(E,stress,vickers,param){
#	theta <- ifelse(E$angle[1]<pi/2,E$angle[1],E$angle[1]-pi/2)
#	U <- ifelse(runif(1)<param[1],1+0.0001*runif(1),
#			rweibull(1,scale=param[2]*(E$ab[1]*1000)^2/(E$ab[2]*1000)/stress/cos(theta)/vickers,
#					   shape=param[3]))
#	
#	# if an already cracked fiber is critical then at least 1 cycle is needed until failure,
#	#	the bit randomness is only to avoid ties
#	
#	return(U)
#}


## Delamination (spheroid projection)
#' Generate delamination times
#' 
#' Generate individual defect times for particle delamination
#' 
#' This kind of failure time (time of delamination from the metal matrix structure) roughly depends on the 
#' projected area of the object, the applied overall stress level and whether the object lies in the interior 
#' of the simulation box or hits one of the box boundaries. The parameter vector is made up of six parameters where 
#' the second and third parameters are used to simulate the defect type \code{crack}, see \code{\link{getCrackTime}}. 
#' The order is as follows: \eqn{p1} probability of already materialized defects, scale factor \eqn{p2}, shape factor
#' \eqn{p3}, shift parameter \eqn{p4} of log times, the slope \eqn{p5} and \eqn{p6} denoting the standard deviation
#' of the random errors of the log times. Only the last three parameters \eqn{p4,p5,p6} are used to determine the defect
#' time for delamination of the considered object.
#' 
#' @param stress	stress level
#' @param E			the object, i.e. spheroid, cylinder, sphere 
#' @param param	    simulation parameter vector
#' @param inF		weightening factor for inner defect projection
#' @param outF		weightening factor for outer defect projection
#' 
#' @return  numeric, the individual delamination time 
#' @author  Felix Ballani 
getDelamTime <- function(E,stress,param,inF=0.5, outF=0.65){
   exp(param[4]+param[5]*log(stress*sqrt(pi)*(attr(E,"area")*1e+06)^0.25*ifelse(attr(E,"interior"),inF,outF))+param[6]*rnorm(1))	
}

#' Plot estimated densities
#' 
#' Plot the estimated densities which result from the randomly 
#' generated individual failure times.
#' 
#' @param dv	  	named list of individual failure times
#' @param main 		title of the plot (optional)
#' @param ...		arguments passed to \code{\link[lattice]{densityplot}}
#' 
#' @example inst/examples/sim.R
showDensity <- function(dv,main="Crack time density estimation", ...) {
	if (!requireNamespace("lattice", quietly=TRUE))
	 stop("Please install package 'lattice' from CRAN repositories before running this function.")
		
	data=data.frame();
	for(i in names(dv)){
		idf=data.frame(x=dv[[i]], "type" = rep(i,length(dv[[i]])))
		data=rbind(data,idf)
	}	
	
	lattice::densityplot(~x, data = data, groups = data$type,
			plot.points = FALSE, ref = FALSE, auto.key = list(space = "top"), main = main,...)	
}
