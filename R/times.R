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
#' @param stress		   stress level for generation of failure times
#' @param vickers   	   Vickers hardness, see details
#' @param param	 		   list of parameter vectors for simulation of failure times for both phases
#' @param fun 			   optional, either \code{lapply} (default) or parllel processing by \code{mclapply}
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
#' @rdname  simCrackTime
#' @export
simCrackTime <- function(S,stress,vickers,param,fun=lapply) UseMethod("simCrackTime",S)

#' @method simCrackTime oblate
#' @export
simCrackTime.oblate <- function(S,stress,vickers,param,fun=lapply)
{ simCrackTime.prolate(S,stress,vickers,param) }
# because the lengths c and a are already switched in E$ab at generation at C-level

#' @method simCrackTime prolate
#' @export
simCrackTime.prolate <- function(S,stress,vickers,param,fun=lapply) {
  simT <- function(E) {
	uv <- numeric(2) # [u,v]
	label <- attr(E,"label")
	
	if(label == "P")
	{
		theta <- try(.getAngle(acos(E$u[3])),silent=TRUE)
		stopifnot(is.numeric(theta))
		
		uv[1] <- getCrackTime(theta,E$acb[1],E$acb[3],stress,vickers,param$P,param$const)
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
 fun(S,simT)
}

#' @method simCrackTime cylinders
#' @export
simCrackTime.cylinders <- function(S,stress,vickers,param,fun=lapply) {
 simT <- function(E) {
	uv <- numeric(2) # [u,v]
	label <- attr(E,"label")
	if(label == "P") {
		theta <- try(.getAngle(acos(E$u[3])),silent=TRUE)
		stopifnot(is.numeric(theta))
		
		uv[1] <- getCrackTime(theta,E$r,0.5*E$h,stress,vickers,param$P,param$const)
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
  fun(S,simT)
}

#' @method simCrackTime spheres
#' @export
simCrackTime.spheres <- function(S,stress,vickers,param,fun=lapply) {
 simT <- function(E) {
	## always delamination for spheres
	label <- attr(E,"label")
	uv <- if(label == "P")	getDelamTime(E,stress,param$P)
		  else getDelamTime(E,stress,param$F)
	list("id"=E$id,"U"=Inf,"V"=uv,"T"=uv,"B"=1,"A"=0,"label"=label)
 }
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
#' @param theta		polar angle
#' @param a			axis length (axis orthogonal to rotational axis)
#' @param b			rotational axis length
#' @param stress	stress level
#' @param vickers   Vickers hardness
#' @param p		    simulation parameter set
#' @param const		constant parameters for crack
#'
#' @return  numeric, the individual fracture time
#' @seealso \code{\link{getDelamTime}}
#' 
#' @author Felix Ballani, M. Baaske
#' @rdname getCrackTime
#' @export
getCrackTime<-function(theta,a,b,stress,vickers,p,const){
	# if an already cracked fiber is critical then at least 1 cycle is needed until failure,
	# the bit randomness is only to avoid ties
	sigf <- stress*cos(theta)*(1-(1/cosh(const$nE*b/a)))
	preconst <- const$sigref^(const$nc-2)*const$Vref^(1/p[3])*(2*const$Ef/const$Em)^(-const$nc)
	#cat("sigf: ",sigf," preconst: ",preconst)
	ifelse(
	 runif(1) < p[1], 1+0.0001*runif(1),
	 rweibull(1,shape=p[3],scale=10^p[2]*preconst*sigf^(-const$nc)*(a)^(-3/p[3]))
    )
}


## Delamination (spheroid projection)
#' Generate interface debonding times
#'
#' Generate individual defect times for particle debonding
#'
#' This kind of failure time (time of debonding from the metal matrix structure) roughly depends on the
#' projected area of the object, the applied overall stress level and whether the object lies in the interior
#' of the simulation box or hits one of the box boundaries. The parameter vector is made up of six parameters where
#' the second and third parameters are used to simulate the defect type \code{crack}, see \code{\link{getCrackTime}}.
#' The order is as follows: \eqn{p1} probability of already materialized defects, scale factor \eqn{p2}, shape factor
#' \eqn{p3}, shift parameter \eqn{p4} of log times, the slope \eqn{p5} and \eqn{p6} denoting the standard deviation
#' of the random errors of the log times. Only the last three parameters \eqn{p4,p5,p6} are used to determine the defect
#' time for debonding of the considered object.
#'
#' @param stress	stress level
#' @param E			object of class "\code{oblate}", "\code{prolate}", "\code{cylinders}", "\code{spheres}"
#' @param param	    simulation parameter vector, see details
#' @param inF		weightening factor for inner defect projection
#' @param outF		weightening factor for outer defect projection
#'
#' @return  numeric, the individual debonding time
#' 
#' @author  Felix Ballani
#' @rdname  getDelamTime
#' @export
getDelamTime <- function(E,stress,param,inF=0.5, outF=0.65){
   exp(param[4]+param[5]*log(stress*sqrt(pi)*(attr(E,"area"))^0.25*ifelse(attr(E,"interior"),inF,outF))+param[6]*rnorm(1))
}

#' Plot estimated densities
#'
#' Plot the estimated densities which result from the randomly
#' generated individual failure times.
#'
#' @param dv	  	named list of individual failure times
#' @param main 		title of the plot (optional)
#' @param ...		optional, additional graphics parameters
#'
#' @example inst/examples/sim.R
#' @author  M. Baaske
#' @rdname  showDensity
#' @export
showDensity <- function(dv,main="Failure time density estimation", ...) {
	if (!requireNamespace("lattice", quietly=TRUE))
	 stop("Please install package 'lattice' from CRAN repositories before running this function.")

	data=data.frame();
	for(i in names(dv)){
		idf=data.frame(x=dv[[i]], "type" = rep(i,length(dv[[i]])))
		data=rbind(data,idf)
	}

	lattice::densityplot(~x, data = data, groups = data$type, kernel="gaussian",
	         plot.points = FALSE, ref = FALSE, xlab=expression(paste("log(T)")),
	          auto.key = list(x = 0, y = 0.9, corner = c(0, 0)), main = main,
	          scales = list(at=c(5,8,10)), ...)


}

## not exported
#.plotDensity <- function(dv,legend.text,main="Failure time density estimation", ...) {
#	.mylegend <- function(plot = TRUE) {
#		legend("topleft",legend.text, cex=1.8, lwd=1.5,lty=c(1,5),..., horiz=FALSE, bty='y',plot = plot)
#	}
#	if (!requireNamespace("lattice", quietly=TRUE))
#		stop("Please install package 'lattice' from CRAN repositories before running this function.")
#
#	data=data.frame();
#	for(i in names(dv)){
#		idf=data.frame(x=dv[[i]], "type" = rep(i,length(dv[[i]])))
#		data=rbind(data,idf)
#	}
#
#	minT <- max(0,min(log10(TL),log10(TQ)))
#	maxT <- min(15,max(log10(TL),log10(TQ)))
#
#	plot(density(dv[[1]]),lty=1,lwd=1.5,cex.lab=1.8,xlab="Individual lifetimes T",
#			ylim=c(0,0.18),cex.axis=1.8, xlim=c(0,maxT),xaxs="i",yaxs="i",xaxt="n",main="")
#	lines(density(dv[[2]]),lty=5,lwd=1.5)
#
#	low <- floor(abs(minT))
#	ticks <- seq(low,ceiling(maxT),by=2)
#	labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
#	axis(1, at=sapply(ticks,function(i) i), labels=labels,cex.axis=1.8)
#
#	.mylegend()
#}
