#' A packed spheroids configuration (particle system)
#'
#' A dataset containing non-overlapping spheroids
#' where each list element consists of the following variables:
#' 
#' \itemize{
#' 	 \item{id}{ the current cluster id}
#' 	 \item{center}{ the center of ball which defines the cluser region}
#' 	 \item{u}{ the orientation vector of spheroid's rotational symmetry axis}
#' 	 \item{ab}{ the axes length}
#' 	 \item{angles}{ colatitude and longitude angle}
#' 	 \item{rotM}{ rotational matrix, center of rotation is the center of the spheroid}
#' 	 \item{radius}{ attribute, random radius of perfect simulation, if enabled}
#' 	 \item{interior}{ attribute, whether the spheroid lies in the interior of the simulation box}
#' 	 \item{label}{ attribute, character, user defined label, here: 'P' (particle) or 'F' (ferrit)}
#'   \item{area}{ attribute, the projected area}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name S
#' @usage data(AL2MC_20p_k10_F2p_S)
#' @format A list of length 4905 containing lists each of length six
#' @author M. Baaske
NULL

#' A cluster configuration of spheroids
#'
#' A dataset containing cluster regions of spheroids
#' where each list element consists of a spheroids list
#' with additional attributes \code{class="prolate"},
#' \code{interior} (logical), indicating whether the cluster
#' fully lies in the interior of the simulation box
#' 
#' @docType data
#' @keywords datasets
#' @name CL
#' @usage data(AL2MC_20p_k10_F2p_CL)
#' @format A list of length 9. Each element contains a list of spheroids  
#' @author M. Baaske
NULL

#' A packed (sphero)cylinder configuration (fiber system)
#'
#' A dataset containing non-overlapping spherocylinders
#' where each list element consists of the following variables:
#' 
#' \itemize{
#' 	 \item{id}{ the current cluster id}
#' 	 \item{center}{ the center of ball which defines the cluser region}
#' 	 \item{origin0}{ the center of the fist cylinder cap}
#' 	 \item{origin1}{ the center of the second cylinder cap}
#'   \item{length}{ length/height of cylinder without radii of caps}
#' 	 \item{u}{ orientation vector of spheroid's rotational symmetry axis}
#' 	 \item{r}{ radius of cylinder}
#' 	 \item{angles}{ colatitude and longitude angle}   
#' 	 \item{rotM}{ rotational matrix, center of rotation is the center of the spheroid}
#' 	 \item{radius}{ attribute, random radius of perfect simulation, if enabled}
#' 	 \item{interior}{ attribute, whether the spheroid lies in the interior of the simulation box}
#' 	 \item{label}{ attribute, character, user defined label, here: 'P' (particle) or 'F' (ferrit)}
#'   \item{area}{ attribute, the projected area}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name SF
#' @usage data(AL2MC_15p_k10_F2p_S)
#' @format A list of length 3550 containing lists each of length nine
#' @author M. Baaske
NULL

#' A cluster configuration of fibers
#'
#' A dataset containing clusters of regions of fibers,
#' including a second phase as spheres where each list element consists
#' of a list of fibers with additional attributes \code{class="cylinder"},
#' \code{interior} (logical), indicating whether the cluster fully lies in
#' the interior of the simulation box.
#' 
#' @docType data
#' @keywords datasets
#' @name CLF
#' @usage data(AL2MC_15p_k10_F2p_CL)
#' @format A list of length 11 containing lists of fibers belonging to the 
#' 		   corresponding cluster
#' @author M. Baaske
NULL

#' Projection of defect types
#' 
#' In case one wishes to analyze the defect projections directly the function returns 
#' the projected objects itself. The two types of defects lead to different projected objects.
#' In case of defect type \code{crack} the circle of maximum radius of
#' the particle is orthogonally projected otherwise the whole particle
#' is projected in the main load direction, i.e. the xy plane.  
#' 
#' Get the projected objects (ellipses) of particles dependent on their defect types
#' 
#' @param S 	list of (non-overlapping) particles (spheroids)
#' @param type	either \code{crack} or \code{delam}
#' @return	    returns a list of 2d ellipse objects. 
getSpheroidProjection <- function(S,type) {
	if(length(S)!=length(type))
		stop(paste0("Crack type (integer) length and number of spheroids do not match"))
	.Call(C_GetSpheroidProjection,S,type)	
}

#' Fiber defect projections
#' 
#' The function draws the defect projections after \code{\link[unfoldr]{cylinders3d}} has
#' been called and returns the points for the convex hull with its points of each
#' projected fiber (spherocylinder). 
#' 
#' @param S    fibers to be projected
#' @param B	   integer vector of length equal to \code{S} of defect types, 
#' 			   \code{B=0} for crack and \code{B=1} for delamination.
#' @param draw logical, \code{draw=TRUE} (default) draw the projected fibers
#' @param conv logical, \code{conv=TRUE} (default) draw the convex hull of projected fibers

#' @param np   number of points used for sampling the convex hull of projections
#' 
#' @return 	   draw defect projections and the convex hull (if enabled) in a 3d plot
#'			   returning its polygonal area and points of the convex hull
#'   
#' @example    inst/examples/cylinder.R
getCylinderProjection <- function(S, B=rep(1,length(S)), draw=TRUE, conv = TRUE, np=20) {
	stopifnot(length(S)==length(B))
	convH <- function(M, fill=TRUE) {
		x <- M[,1]; 
		y <- M[,2]		
		m <- chull(x,y);
		M <- cbind(x[m],y[m])			
		list(splancs::areapl(M),M)
	}	
	
	P <- .Call(C_GetCylinderProjection,S,B,np)
	if(draw) {
		invisible(lapply(P,
				function(x) {
					if(attr(x,"type")) {
						M <- convH(x)[[2]]
						rgl::polygon3d(M[,1],M[,2],rep(0,NROW(M)),fill=TRUE)
					} else rgl::polygon3d(x[,1],x[,2],rep(0,NROW(x)),fill=TRUE)
				}
			)
		)
	}
	if(conv) {		
		A <- convH(do.call(rbind,P),fill=FALSE)
		if(draw) 
		 rgl::polygon3d(A[[2]][,1],A[[2]][,2],rep(0,NROW(A[[2]])),fill=FALSE)
		return (list("points"=P,"convH"=A))
	}
	return (P)
}

#' Sphere projections
#' 
#' The function draws the defect projections in a 3d plot
#' and returns the points of its  the convex hull together
#' with points of each projected particle (here a sphere). 
#' 
#' @param S    particles to be projected
#' @param draw logical, \code{draw=TRUE} (default) draw the projected spheres
#' @param conv logical, \code{conv=TRUE} (default) draw the convex hull

#' @param np   number of points used for sampling the convex hull of projections
#' 
#' @return 	   draw projections and the convex hull (if enabled) in a 3d plot
#'			   returning its polygonal area and points of the convex hull
getSphereProjection <- function(S, draw=TRUE, conv = TRUE, np=20) {	
	convH <- function(M, fill=TRUE) {
		x <- M[,1]; 
		y <- M[,2]		
		m <- chull(x,y);
		M <- cbind(x[m],y[m])			
		list(splancs::areapl(M),M)
	}	
	
	P <- .Call(C_GetSphereProjection,S,np)
	if(draw) {
	  invisible(lapply(P,function(x) { rgl::polygon3d(x[,1],x[,2],rep(0,NROW(x)),fill=TRUE) }))
	}
	if(conv) {		
		A <- convH(do.call(rbind,P),fill=FALSE)
		if(draw) 
			rgl::polygon3d(A[[2]][,1],A[[2]][,2],rep(0,NROW(A[[2]])),fill=FALSE)
		return (list("points"=P,"convH"=A))
	}
	return (P)
}


#' Crack time simulation
#'
#' Simulate crack times of particles
#' 
#' The function randomly generates phase dependent failure times of defect types \code{crack} and \code{delam}.
#' The accumulation structure of defects is initialized containing the failure times of objects in ascending order.
#' The failure times of the defect type \code{crack} follow a Weibull distribution, see \code{\link{simCrackTime}}.
#' The failure times of the defect type \code{delam} roughly depend on the projected area of the object, the 
#' applied overall stress amplitude and whether the object lies in the interior of the simulation box or hits
#' one of the box boundaries. The argument \code{param} consists of a distribution parameter list which contains 
#' numeric vectors for both reinforcement objects (labled by \code{P}) and an optional second phase (labled by \code{F}).
#' If no second phase is considered the corresponding parameter set is simply ignored. Each parameter vector is made up of
#' six parameters in the following order: \eqn{p1} probability of already materialized defects, scale factor \eqn{p2},
#' shape factor \eqn{p3}, shift parameter \eqn{p4} of log times, the slope \eqn{p5} and \eqn{p6} as the standard deviation
#' of the random error of log times. 
#'  
#' @param S 			   (non-overlapping) geometry system
#' @param param 		   parameters for generating failure times
#' @param vickers		   vickers hardness
#' @param stress    	   stress level
#' @param verbose 		   logical: ignored
#' @param parallel.option  optional: In case of \code{mclapply} the function 
#' 									 \code{mclapply} is used from the parallel package
#' 
#' @return list of increasing failure times
#' @seealso \code{\link{getCrackTime}}, \code{\link{getDelamTime}}
#' @example inst/examples/sim.R
simTimes <- function(S, param, vickers, stress, verbose = FALSE,
						parallel.option = getOption("parallel.option","lapply")) {	
	## sim times			    
	CLT <- simCrackTime(S,stress,vickers,param)
	
	## sort ascending by times
	CLT <- CLT[order(sapply(CLT, `[[`, "T"))]
	return ( CLT )		
}

.VV <- function(S,box) {
	sum(unlist(lapply(S,
	 function(x) 4/3*pi*x$ab[1]^2*x$ab[2])))/(box$xrange[2] * box$yrange[2] * box$zrange[2])
}


#' Simulate particle system (primary phase) 
#' 
#' Poisson spheroid system
#' 
#' The function generates a constant size Poisson spheroid system with intensity parameter \code{lam}
#' and random planar (with respect to the xz plane) orientation distribution. The spheroids are labeled
#' by \eqn{P} to denote the primary particle phase which usually plays the role of some reinforcement 
#' material type in real life specimen. In order to simulate the fatigue lifetime model, see \code{\link{simDefect}},
#' a non-overlapping configuration of particles is required which could be generated for instance by application of
#' the the well-known random sequential adsorption (RSA) method, see \code{\link{rsa}}. Alternatively the Force-biased 
#' algorithm could be used as well, see reference below.
#' 
#' @param theta   simulation parameter list
#' @param lam 	  intensity parameter of the underlying Poisson point process
#' @param box	  the simulation box
#' @param mu	  reference direction of particles, here \code{mu=c(0,1,0)} (default)
#' @param verbose logical, not ues yet 
#' 
#' @return 	list of spheroids  
#' @references  \itemize{
#' 				  \item{}{A. Bezrukov and D. Stoyan. Simulation and statistical analysis of random packings of ellipsoids. Particle & Particle Systems Characterization, 23(5):388 - 398, 2006.}
#' 				  \item{}{J.W. Evans. Random and cooperative sequential adsorption. Rev. Mod. Phys., 65: 1281-1304. 1993.}
#' 				}
#' @seealso \code{\link{simSpheroidSystem}}
simParticle <- function(theta,lam,box,mu=c(0,1,0),verbose=FALSE) {	
	S <- unfoldr::simSpheroidSystem(theta,lam,size="const",orientation="rbetaiso",
			        box=box,mu=mu,pl=101,label="P")	
	if(verbose) {
      cat("V_V: ",.VV(S,box),"\n")
	}
	return (S)
}

#' Simulate a second phase
#' 
#' Poisson spheroid system
#' 
#' The function generates a system of constant size Poisson spheroidal objects with intensity \code{lam} and random planar 
#' orientation distribution. The spheroidal objects are labeled by \eqn{F} to denote the secondary 
#' particle phase, denoting for instance some disturbing ferrit inclusions in real life specimen. 
#' 
#' @param param    simulation parameter list
#' @param lam 	   the intensity parameter of the underlying Poisson point process
#' @param box	   the simulation box
#' @param mu	   reference direction of particles, here \code{mu=c(0,1,0)} (default)
#' @param verbose  logical, ignored 
#' 
#' @return 	list of spheroids  
#' @seealso \code{\link[unfoldr]{simSpheroidSystem}}
simFerrit <- function(param, lam, box, mu=c(0,1,0), verbose = FALSE) {
	F <- unfoldr::simSpheroidSystem(param,lam=lam,size="const",orientation="rbetaiso",
			        box=box,mu=mu,pl=101,label="F")	
	if(verbose) {
    	cat("V_V: ",.VV(F,box),"\n")
	}
	return (F)
}

#' Critical area
#' 
#' Calculate critical area 
#' 
#' The critical area (see reference below) is roughly defined by the projected defect area
#' above which the specific material defect becomes particularly hazardous for the failure
#' of the whole specimen.
#' 
#' @param vickers 	specific Vickers hardness
#' @param stress 	stress level
#' @param factor    specific material dependent factor, default set to \code{1.56} for inner defects
#' 
#' @return 	critical area in \eqn{\mu m^2}
#' @references  Y. Murakami (2002). Metal Fatigue: Effects of Small Defects and Nonmetallic Inclusions. Elsevier, Amsterdam.
areaCrit <- function(vickers, stress, factor=1.56) {
	((vickers + 120)/stress*factor)^12	
}

#' Extract defects
#' 
#' Extract all defects with some area value larger than \code{area}
#' 
#' The function extracts all defects (convex hulls of projected defect areas)
#' which have a larger area value than the predefined \code{area} value for some
#' further analysis.
#' 
#' @param clust     list of defects, see \code{\link{simDefect}}
#' @param area  	the area value lower bound
#' @return      	list of clusters
extractClusters <- function(clust,area) {
	id <- lapply(clust, function(x) if(length(x$id)>1 && x$A>area) TRUE else FALSE)
	return ( clust[unlist(id)])
}

#' Defect accumulation 
#' 
#' Simulate fatigue lifetime of the geometry system
#' 
#' The function simulates the fatigue lifetime model for the given geometry system. The overall fatigue lifetime 
#' is determined by the first accumulated convex area of projected defects which exceeds the 
#' predefined critical area. The process of defect accumulation is based on the general results for VHCF fatigue behaviour of
#' material structures published by Murakami. Here we treat two nearby defects as an enlarged 
#' common one if the weighted (single-linkage) Euclidean distance inbetween (in 3D) is less than some portion (possibly including
#' some material specific constant) of the minimum of the corresponding square roots of projected defect areas. The weight is chosen 
#' as the reciprocal value of the cosine of the angle between the main load direction and the rotational axis of the chosen object
#' having minimum distance to the compared defect .  
#' Further, the argument \code{opt} defines the list of controls of fixed parameters with the following elements:
#' \code{vickers} Vickers hardness, \code{distTol} as a (proportion) factor of the minium the corresponding square roots of 
#' projected defect areas (might be less than 0.95 but in general depends on the inverstigated material structure),
#' \code{inAreafactor} default set to value 1.56, \code{outAreafactor} default set to value 1.43, \code{pointsConvHull} number of
#' sampling points of each particle projection border for approximating the defect projection convex area, \code{scale} scaling factor
#' (should be set to \code{1e+06} for \eqn{\mu m^2}), \code{pl} print level of information output with some verbose messages for any value
#' larger than 100. In this case only the accumulated defect projections after adding the last individual failed particle are returned.
#' Here the last element contains the defect with the highest projected area which leads to the overall failure. 
#' 
#' @param stress	stress level
#' @param S   		non-overlapping particle system
#' @param CLT  		the start (initialized) cluster of particles
#' @param opt	    control list of fixed parameters, see \code{\link{simTimes}}
#'
#' @return 		    list of clusters, the following values are returned for further analysis:
#' 					\itemize{
#' 						\item{id}{ the particle id numbers in the accumulated defect}
#' 						\item{n}{ internal: accumulated defect id, zero means last element of cluster}
#' 						\item{B}{ corresponding defect types of projected objects}
#' 						\item{interior}{ whether each object is interior or not}
#' 						\item{A}{ projected areas, possibly convex area}
#' 						\item{inner}{ whether the defect lies in the interior or not}
#' 						\item{T}{ random individual times}
#' 						\item{label}{ phase labels }
#' 						\item{areaMax}{ attribute, the maximum area of defect projections}
#' 						\item{aIn}{ attribute, critical area for inner defect projections}
#' 						\item{aOut}{ attribute, critical area for outer defect projections}
#' 						\item{maxSize}{ attribute, maximum number of projected objects found in a defect}
#' 					}
#' @author			Markus Baaske 
simDefect <- function(stress,S,CLT,opt) {	
	# areas
	aIn <- (((opt$vickers + 120)/stress*opt$inAreafactor)^12)/opt$scale  
	aOut <- (((opt$vickers + 120)/stress*opt$outAreafactor)^12)/opt$scale
	
	env <- environment()
	tryCatch(				
		.Call(C_SimDefect,as.character(substitute(S,S)),CLT,
				opt$distTol,aIn,aOut,opt$pl,env)
		, error = function(e) {						
					structure(e,class=c("error","condition"),
					 error=simpleError(.makeMessage("Error in cluster construction.\n")))
				  }
	)	
}


#' Fatigue lifetime simulation
#' 
#' Simulation of fatigue lifetime model
#' 
#' We provide a wrapper function for \code{\link{simDefect}} to simulate the proposed fatigue lifetime model including
#' the generation of individual failure times with additional information of the responsible accumulated defect projection
#' area (simply called defect area) which leads to the failure of the system. This information is returned in the following 
#' list with elements: \code{crack} the responsible defect itself, \code{T} the individual failure times of particles aggregated
#' as a cluster of particles, \code{ferrit} whether any secondary phase (here called \code{ferrit}) is involved for overall failure,
#' \code{interior} whether the defect occurs in the interior or at the boundaries of the simulation box. The optional
#' list argument \code{CL} defines clustered regions, see \code{\link{simCluster}}, of particles sticking together more close than
#' others. This is useful in case one also wishes to model densely lieing agglomerations of objects (i.e. clusters of particles or 
#' fibers) as an obvious phenomenon of some specimen in mind. This list is made up of the following elements: \code{id}, the id of
#' the region, \code{center} the center point of the corresponding region, \code{r} the radius (as the Euclidean distance from the 
#' center point) which within all particles belong to this region, \code{interior} whether any particle within radius \code{r} hits 
#' the box boundaries of the simulation box.
#' 
#' 
#' @param stress		applied stress level
#' @param S 			non-overlapping particle system
#' @param opt			control list for simulation of individual failure times
#' @param param 		parameter list for generating individual failure times
#' @param last.defect 	logical, \code{last.defect=FALSE} (default) return all defect 
#' 						accumulations or only the last one
#' @param CL			optional, NULL (default) predefined clustered regions of particles
#' 
#' @return 		a list made up of the following objects if \code{last.defect=FALSE} (default):
#' 				\itemize{
#' 					\item{fracture}{ return value of \code{\link{simDefect}} }
#' 				    \item{times}{ return value of \code{\link{simTimes}} }
#' 					\item{cl_info}{ additional cluster information of responsible defect cluster}
#' 				}
#' 				otherwise only \code{cl_info} is returned. 
#' 
#' @example inst/examples/simFailure.R	
simFracture <- function(stress, S, opt, param, last.defect = FALSE, CL = NULL) {
	## simulate times
	CLT <- simTimes(S,param,opt$vickers,stress)
	
	RET <- simDefect(stress,S,CLT,opt)	
	if(!is.null(attr(RET,"error")))
	  return (RET)
	
	N <- length(RET)
	crack <- RET[[N]]
	cl_info <- if(is.null(CL))
			      list("T"= crack$T[1], "count"=length(crack$label), "ferrit"=any("F" %in% crack$label),
						"num.ferrit"=sum(crack$label!="P"), "interior"=all(crack$interior))  
				else
				  list("T"= crack$T[1], "count"=length(crack$label), "ferrit"=any("F" %in% crack$label),
					   "num.ferrit"=sum(crack$label!="P"), "interior"=all(crack$interior),
					   "inCluster"=any(unlist(lapply(CL,function(x,y) any(y$id %in% x$id), y=crack))))
	
	if(last.defect)
	 return (list("cl_info"=cl_info))	
    return(list("fracture"=RET,"times"=CLT,"cl_info"=cl_info))
		 
}

#' Draw accumulation path defect areas
#' 
#' Plot accumulation of defect projection areas
#' 
#' Simple plotting function to visualize the process of accumulation of defect areas separately
#' for interior particles and particles hitting the boundaries of the simulation box (specimen)
#' 
#' @param CL 			cluster of aggregated (densified objects), see \code{\link{simFracture}}
#' @param last.path		logical, \code{last.path=FALSE} (default), show the critical area accumulation path
#' 					    leading to the overall failure
#' @param log.axis		character, \code{log.axis='x'} (default), switch to logarithmic scale
#' @param main			optional, main title of the plot
#' 
#' @return 				the accumulation path of defect areas as a list containing the
#' 						starting and end points of each horizontal line (the current convex hull area value)
#' 					    together with the starting point of the next convex hull of defect projections and 
#' 					    whether the defect occured in the interior or at the boundaries of the simulation box	  
plotDefectAcc <- function(CL,last.path=FALSE,log.axis="x",
		main="Accumulation of particle's projection areas") 
{	
	.axTexpr <- function(side, at = axTicks(side, axp=axp, usr=usr, log=log),
			axp = NULL, usr = NULL, log = NULL)	{
		eT <- floor(log10(abs(at)))# at == 0 case is dealt with below
		mT <- at / 10^eT
		ss <- lapply(seq(along = at),
				function(i) if(at[i] == 0) quote(0) else
						substitute(A %*% 10^E, list(A=mT[i], E=eT[i])))
		do.call("expression", ss)
	}	
	
	aIn <- attr(CL,"aIn")
	aOut <- attr(CL,"aOut")
	
	N <- length(CL)	
	maxT <- max(CL[[N]]$T)
	minT <- min(CL[[1]]$T)
	maxA <- CL[[N]]$A[1]	
	
	A <- ifelse(CL[[N]]$inner,aIn,aOut)
	yr <- c(0,maxA)
	plot(c(minT,maxT),yr,type="n",log=log.axis)
	legend.size <- legend("toplef",c("Cluster interior","Cluster at bounds","Critical area limit"),cex=.8, 
			col=c("darkgreen","red",col),pch=c(21,24,NA), lty=c(1,1,5),lwd=c(0.5,0.5,1),
			horiz=FALSE, bty='y',plot = FALSE)
	mrange <- range(yr)
	mrange[2] <- 1.05*(mrange[2]+legend.size$rect$h)
	
	plot(c(minT,maxT),mrange,type="n",log=log.axis, xaxt=ifelse(log.axis=="x","n","s"),
			ylab=expression(paste("Area (",mu,"m)")^2), 
			xlab=ifelse(log.axis=="x","log(N)","N"),main=main )
	
	if(log.axis=="x") {
		aX <- axTicks(1)
		axis(1, at=aX, labels= .axTexpr(1, aX))
	}
	## draw critical area line
	col <- ifelse(CL[[N]]$inner,"darkgreen","red")		
	abline(h=A,col=col,lty=5,lwd=1)
			
	# get the accumulation paths
	if(last.path)
	 CL <- list(CL[[N]])
	
	L <- lapply(1:(length(CL)),
					     function(i) {							
							 L <- .ROW2LIST(cbind(CL[[i]]$T,CL[[i]]$A,CL[[i]]$n),arg.names=c("T","A","n"))
							 ok <- sapply(L,function(x) x$n>0 )
							 L <- L[ok]
							 L <- L[order(sapply(L, `[[`, "T"))]
							 T <- sapply(L,'[[', 'T')
							 A <- sapply(L,'[[', 'A')
							 G <- if(length(L)>1)
							       c(lapply(1:(length(L)-1),
									   function(k) { 
										   Anew <- NULL
										   t1 <- ifelse(A[k+1]<A[k], maxT, {Anew <- A[k+1]; T[k+1] })
										   
										   list("p0" = c("T"=T[k],"A"=A[k]),
									            "p1" = c("T"=t1,"A"=A[k]), 
										  	 	"p2" = if(!is.null(Anew)) c("T"=t1,"A"=Anew) else NULL) }), 
					 		 		 	 list(list("p0"=unlist(L[[length(L)]]),"p1"=c(maxT,L[[length(L)]]$A),"p2"=NULL)) )
						 		  else
									 list(list("p0"=unlist(L),"p1"=c(maxT,L[[1]]$A),"p2"=NULL))
							 attr(G,"inner") <- CL[[i]]$inner
							 G
						 }
	 	  )				  
		  
	  if(last.path) {
		  ok <- sapply(L[[1]], function(x) is.null(x$p2))
		  ok[length(ok)] <- FALSE
		  for(i in length(ok):1) {
			  if(ok[i])
			    break
		  }
		  L <-L[[1]][(i+1):length(L[[1]])]
		  attr(L,"inner") <- CL[[1]]$inner
		  L <- list(L)
	  }		 

	 func <-  function(x, inner) {
		 col <- ifelse(inner,"darkgreen","red")  
		 points(x$p0[1],x$p0[2],bg=col, col=col, pch=ifelse(inner,21,24))
		 
		 if(!is.null(x$p2))	{
			 # plot vertical line to next cluster
			 points(x$p1[1],x$p1[2],col=col, pch=ifelse(inner,21,24))
			 segments(x$p1[1], x$p1[2], x$p2[1], x$p2[2], col="black", lty=3, lwd=1)
		 } else if(x$p1[1]<maxT)
			 points(x$p1[1],x$p1[2],col=col, pch=ifelse(inner,21,24))
		 
		 segments(x$p0[1], x$p0[2], x$p1[1], x$p1[2], col=col,lwd=0.5)			  
		 
	 }		  
	
	 ## draw lines 
	 invisible(lapply(L,function(x) lapply(x,func,inner=as.logical(attr(x,"inner")))))
	 	 
	 ## the legend
	 legend("topleft",c("Cluster interior","Cluster at bounds","Critical area limit"),cex=.8, 
		 col=c("darkgreen","red",col),pch=c(21,24,NA), lty=c(1,1,5),lwd=c(0.5,0.5,1),
		 horiz=FALSE, bty='y')
	 
	 return (L)
}

## Convert matrix rows to list
.ROW2LIST  <- function(x, arg.names=NULL) {
	lapply(seq_len(nrow(x)), function(i) {
				lst <- as.list(x[i,])
				if(!is.null(arg.names))
					names(lst) <- arg.names
				lst
			})
}

## Splitting a list in nearly equal chunks
## adopted from snow package with slight modifications
## to handle special inputs, used for ballancing cluster 
## loads for MPI or SOCKS connections
.splitList <- function(x,n) {
	nx <- length(x)
	if( nx > n ) {
		i <- 1:nx
		if(n == 0)
			list()
		else if(n == 1 || length(x) == 1) 
			list(x)
		else {	 
			q <- quantile(i, seq(0, 1, length.out = n + 1))	
			lapply(structure(split(1:length(x),
								cut(i, round(q), include.lowest = TRUE))
							,names=NULL),function(i) x[i])
		}
	} else {
		x
	}
}

#' Woehler experiment
#' 
#' Simulate a Woehler experiment 
#' 
#' As done in a real life scenario the Woehler diagram measures the applied stress amplitude versus load cycles.
#' For each given stress level the function is intended to be an all-in-one wrapper function for fatigue
#' lifetime model simulations for different stress levels which returns a matrix of failure times where 
#' each row corresponds to one stress level. The Woehler experiments can be simulated in parallel based
#' on using the package \code{snow} possibly together with \code{Rmpi} for an \code{MPI} cluster object.
#' 
#' @param S 			   non-overlapping object system
#' @param CL			   predefined clustered regions of objects, see \code{\link{simFracture}}
#' @param param	    	   parameter list for random generation of individual failure times
#' @param opt			   control parameters, see \code{\link{simTimes}}
#' @param stress		   list of stress levels
#' @param cl			   optional, parallel cluster object (see \code{\link[snow]{makeCluster}})   
#' @param parallel.option  optional, if \code{parallel.option=mclapply} the function \code{\link[parallel]{mclapply}}
#' 								 is used, otherwise the function \code{lapply}.
#' 
#' @return				   matrix of failure times, first colunm corresponds to the times
#' 						   and the second to the stress level
#'  
#' @example inst/examples/woehler.R
#' @author  Markus Baaske 
woehler <- function(S, CL, param, opt, stress, cl = NULL,
						parallel.option = getOption("parallel.option","lapply"))
{			
	simF <- function(s,...) simFracture(s,...,last.defect=TRUE)$cl_info					
	func <- get(parallel.option)
	
	p <- tryCatch({
			if(!is.null(cl)) {
				m <- min(length(stress),length(cl))
				if(any(class(cl) %in% c("MPIcluster"))) {
					cat("Using MPI cluster...\n")
					parLapplyLB(cl[seq_len(m)], stress, simF, opt=opt, param=param, CL=CL, S=S)
				} else {		
					cat("Using some other type of cluster...\n")
					do.call(c,clusterApply(cl[seq_len(m)], .splitList(stress, length(cl)), func, simF,
									opt=opt, param=param, CL=CL, S=S))
				}
			} else {
	 			func(stress, simF, opt=opt, param=param, CL=CL, S=S)		
			}						
		  },error=function(e) {
			  structure(e,class=c("error","condition"),
					  error=simpleError(.makeMessage("Woehler experiment failed.\n")))
		    }
		)
	if(!is.null(attr(p,"error")))
	  return (p)
	as.data.frame(
	  do.call(rbind,
	     func(seq_len(length(p)),
		   function(i) c("stress"=stress[[i]],as.vector(p[[i]])))
 	  )
 	)	
}

#' Plot Woehler diagram
#' 
#' Plot a Woehler diagram for simulated Woehler experiments
#' 
#' The simulated results of the Woehler experiments are summarized in 
#' the load-cycle diagram as usually done in the Woehler diagram.
#' The arguments of \code{syms} and \code{cols} are recycled if \code{nsim} 
#' exceeds one of their lengths.
#' 
#' @param W 				results from function \code{\link{woehler}}
#' @param yrange			range of stress levels, the min/max values of stress levels (default)
#' @param cols 	   			vector of colors, possibly being replicated
#' @param syms				vector of symbol ids, see \code{pch} from \code{\link{par}} 
#' @param main 				main title of the plot
#' @param ...				graphic parameters passed to function \code{\link{points}}
#' 
#' @return 					\code{NULL}
woehlerDiagram <- function(W, yrange = NULL,
		                    cols = c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF"),
							syms = c(0,1,2,5,6), main = "Woehler experiment",...) {	

	u <- sort(as.numeric(unique(W[,1])))
	m <- length(u)
	nsim <- NROW(W)/as.numeric(m)		  
	
	if( !(nsim>0) || NROW(W)%%nsim != 0)
		stop("nsim does not match rows of Woehler matrix.")
	
	matsplit <- function(M, nr, nc) {
		simplify2array(lapply(
				split(M, interaction((row(M)-1)%/%nr+1,(col(M)-1)%/%nc+1)),
					function(x) {dim(x) <- c(nr,nc); x;}
				))
	} 

	w <- cbind( log10(as.numeric(W[,2])), as.numeric(W[,1]))
	yrange <- if(is.null(yrange)) {
				range(c(min(w[,2]),max(w[,2])))
			  } else {
				  yr <- range(c(min(min(w[,2]),yrange[1]),max(max(w[,2]),yrange[2])))
				  if(min(yr)==max(yr))
					  yr <- c(yr[1]-20,yr[2]+20)
				  yr
			  }
	  
	minT <- min(w[,1])
	maxT <- max(w[,1]) 
	
	plot(NULL,xlim=c(minT,maxT), ylim=yrange,xlab=as.expression(paste(substitute(log10),"(N)")),
			log="", xaxt="n",yaxt="n",ylab="Stress amplitude (MPa)", main=main)
		
	syms <- rep(syms, length.out=m)
	cols <- rep(cols, length.out=m)
	
	if(m>1) {
		A <- matsplit(w,nsim,2)
		for(i in 1:m )
		  points(A[,,i], pch=syms[i], col=cols[i])
	} else points(w, pch=syms[1], col=cols[1])
				
	# x values
	low <- floor(abs(minT))
	ticks <- seq(low,ceiling(maxT),by=1)
	labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
	axis(1, at=sapply(ticks,function(i) i), labels=labels)
	invisible(sapply(ticks,	function(i) abline(v=i, col="gray",lty=3, lwd=1)))
	
	# y values
	d <- if(m>1) abs(u[2]-u[1]) else 10
	low <- floor(abs(min(u)-yrange[1])/d)
	high <- ceiling(abs(max(u)-yrange[2])/d)
	u <- c(sapply(low:1,function(i) min(u)-i*d),u)
	u <- c(sapply(1:high,function(i) max(u)+i*d),u)
		
	axis(2, at=sapply(u,function(i) i), labels=u)
	invisible(sapply(u,	function(i) abline(h=i, col="gray",lty=3, lwd=1)))
	NULL
}



#' Draw defect projections of particles
#' 
#' The function draws the defect projections in a 3d plot after \code{spheroids3d}
#' has been called. 
#' 
#' @param E    a list of ellipses (the projected particles)
#' @param conv logical, if \code{conv=TRUE} (default) the convex hull is drawn
#' @param np   number of points used for sampling the convex hull of projections
#' 
#' @return 	   draw projections and the convex hull (if enabled) in a 3d plot
#'			   return area of polygon and points of convex hull  
drawProjection <- function(E, conv = TRUE, np=25) {
	.pointsOnEllipse <- function(E,t) {		
		xt <- E$center[1] + E$ab[1]*cos(t)*cos(E$phi)-E$ab[2]*sin(t)*sin(E$phi)
		yt <- E$center[2] + E$ab[1]*cos(t)*sin(E$phi)+E$ab[2]*sin(t)*cos(E$phi)
		zt <- rep(0,length(t))
		cbind(xt,yt,zt)
	}
	
	.plotEllipse <- function(x) {
		M <- .pointsOnEllipse(x,t)
		rgl::polygon3d(M[,1],M[,2],M[,3],fill=TRUE)
	}	
	s <- 2*pi/np
	t <- seq(from=0,to=2*pi,by=s)	
	invisible(lapply(E,function(x) .plotEllipse(x) ))
	
	# draw convex hulls
	if(conv) {
		M <- .Call(C_GetPointsForConvexHull,E,np)	
		x <- M[,1]; y <- M[,2]		
		m <- chull(x,y);
		M <- cbind(x[m],y[m])		
		rgl::polygon3d(x[m],y[m],rep(0,length(x[m])),fill=FALSE);
		list(splancs::areapl(M),M)
	}
}


#' Draw defect projections
#' 
#' Draw defect projections for spheroids, spherocylinders or spheres.
#' The method first constructs the projected objects sampling some points
#' on the borders and calculates the convex hull of such points.
#' The defect object as returned by \code{\link{simFracture}} 
#' and used to show the accumulated defects.
#'  
#' @param F    		particle system
#' @param D 		defects object
#' @param ...  		further arguments passed to 
#'			   		 \code{\link{getSphereProjection}},					 
#' 				     \code{\link{getSpheroidProjection}},
#' 					 \code{\link{getCylinderProjection}}
#' 
#' @return 			\code{NULL}					 
drawDefectProjections <- function(F,D,...) UseMethod("drawDefectProjections",F)

#' @rdname drawDefectProjections
#' @method drawDefectProjections oblate
#' @S3method drawDefectProjections oblate
drawDefectProjections.oblate <- function(F,D,...) drawDefectProjections.prolate(F,D,...)

#' @rdname drawDefectProjections
#' @method drawDefectProjections prolate
#' @S3method drawDefectProjections prolate
drawDefectProjections.prolate <- function(F,D,...) {
	invisible(lapply(D,
			function(x) 				
			  drawProjection( getSpheroidProjection(F[x$id],as.integer(x$B)),...)			
		)
	)
}

#' @rdname drawDefectProjections
#' @method drawDefectProjections cylinder
#' @S3method drawDefectProjections cylinder
drawDefectProjections.cylinder <- function(F,D,...) {
	invisible(lapply(D, function(x)  getCylinderProjection(F[x$id],as.integer(x$B),...)))
}

#' @rdname drawDefectProjections
#' @method drawDefectProjections sphere
#' @S3method drawDefectProjections sphere
drawDefectProjections.sphere <- function(F,D,...) {
	invisible(lapply(D, function(x)  getSphereProjection(F[x$id],...)))
}