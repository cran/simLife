# translate angle to [0,pi/2]
.getAngle <- function(phi) {
	if(phi<=pi/2) phi
	else {
		if( phi <=pi) pi-phi
		else if(phi<1.5*pi) phi%%(pi)
		else 2*pi-phi
	}
}

# get the correct spheroid (rotation-size) matrix
.getA <- function(s) {
	#A <- matrix(c(1.0/s$ab[1]^2,0,0,0,1.0/s$ab[1]^2,0,0,0,1.0/s$ab[2]^2),nrow=3)
	A <- diag(c(1.0/s$acb[1]^2,1.0/s$acb[2]^2,1.0/s$acb[3]^2))
	return (t(s$rotM) %*% A %*% s$rotM)
}

# slice some vector object
.slice <- function(x,n) {
	N <- length(x);
	lapply(seq(1,N,n),function(i) x[i:min(i+n-1,N)])
}

# calculate spheroid's touch distance coefficient
.spheroidTouchDistCoeff <- function(S1,A1,S2,A2) {
	A <- try(solve(A1) + solve(A2),silent=TRUE)
	if(inherits(A,"try-error")) {
		stop("Inversion of A1+A2 failed.")
	}
	invA <- try(solve(A),silent="TRUE")
	if(inherits(invA,"try-error"))
		stop("Inversion of matrix A failed")
	r <- S1$center-S2$center
	p <- t(r)%*%invA%*%r
	as.vector(sqrt(2.0/max(p,1e-10)))
}

# check whether two spheroids overlap
.checkOverlap <- function(S1,S2,A1=.getA(S1),A2=.getA(S2)) {
	# TRUE in case of overlap
	return(1.0 < .spheroidTouchDistCoeff(S1,A1,S2,A2))
}

# check whether two spherocylinders overlap
.checkOverlapCylinder <- function(q,p)
{
	d <- .C(C_sdm,as.numeric(q$center-p$center),
			as.numeric(q$u),
			as.numeric(p$u),
			as.numeric(q$h*0.5),
			as.numeric(p$h*0.5),
			d=numeric(1))$d

	return ( d < ((q$r+p$r)^2 - 1e-10))	
}

.checkOverlapSphere <- function(q,p) {
	return ( sum( (q$center-p$center)^2 ) < (q$r+p$r)^2 - 1e-10 )
}

# check energy
# no overlapp in case of (t < 1)
.checkEn <- function(S,FUNCTION) {
	ret <- list()
	for(i in 1:(length(S)-1)) {
		for(k in (i+1):length(S)) {
			if(FUNCTION(S[[i]],S[[k]])) # if overlerapping then add
			  ret <- c(ret,k)
		}
	}
	ret
}

# Energy SA algorithm
.EnSpheroid <- function(x,S,clust,weight) {
	Snew <- S
	n <- length(Snew)
	for(i in 1:n) {
		if(attr(Snew[[i]],"label")!="F")
		 Snew[[i]]$center <- Snew[[i]]$center-x[i]*((Snew[[i]]$center-clust$center)/sqrt(sum((Snew[[i]]$center-clust$center)^2)))
	}

 	pot <- 0
	for(i in 1:(n-1)) {
		S1 <-  Snew[[i]]
		A1 <- .getA(S1)
		for(k in (i+1):n) {
			t <- .spheroidTouchDistCoeff(S1,A1,Snew[[k]],.getA(Snew[[k]]))
			if(1<t)
			 pot <- pot + t*sqrt(sum((S1$center-Snew[[k]]$center)^2))
		}
	}
	sum(unlist(lapply(S,function(s) sqrt(sum((s$center-clust$center)^2))))-x) + weight*pot
}


# minimum touch distance in centers directions of cylinders
.contactRadius <- function(p,q)  { # (ci,cj) (c1,c2)
	invR <- try(solve(q$rotM), silent=TRUE)
	if(inherits(invR,"try-error"))
		stop("Inversion of rotation matrix failed... -> exiting")
	d <- p$center-q$center
	rmax <- .C(C_ContactRadius,
			    as.numeric(p$u),
				as.numeric(p$h*0.5),
				as.numeric(q$h*0.5),
				as.numeric(p$r),
				as.numeric(q$r),
				as.numeric(invR),
				as.numeric(d),
				rmax=numeric(1))$rmax
		
		return ( abs(sqrt(sum(d^2)) - rmax) )
}

.EnCylinder <- function(x,S,clust,weight) {
	Snew <- S
	n <- length(Snew)
	for(i in 1:n)
	 if(attr(Snew[[i]],"label")!="F")
      Snew[[i]]$center <- Snew[[i]]$center-x[i]*((Snew[[i]]$center-clust$center)/sqrt(sum((Snew[[i]]$center-clust$center)^2)))
	pot <- 0
	for(i in 1:(n-1)) {
		S1 <-  Snew[[i]]
		for(k in (i+1):n) {
			if(.checkOverlapCylinder(S1,Snew[[k]]))
			 pot <- pot + .contactRadius(S1,Snew[[k]])
		}
	}
	sum(unlist(lapply(S,function(s) sqrt(sum((s$center-clust$center)^2))))-x) + weight*pot
}

.EnSphere <- function(x,S,clust,weight) {
	Snew <- S
	n <- length(Snew)
	for(i in 1:n)
	  if(attr(Snew[[i]],"label")!="F")
		Snew[[i]]$center <- Snew[[i]]$center-x[i]*((Snew[[i]]$center-clust$center)/sqrt(sum((Snew[[i]]$center-clust$center)^2)))
	pot <- 0
	for(i in 1:(n-1)) {
		S1 <-  Snew[[i]]
		for(k in (i+1):n) {
			if(.checkOverlapSphere(S1,Snew[[k]]))
			  pot <- pot + abs( (S1$r+Snew[[k]]$r) - sqrt(sum((S1$center-Snew[[k]]$center)^2)))

		}
	}
	sum(unlist(lapply(S,function(s) sqrt(sum((s$center-clust$center)^2))))-x) + weight*pot
}


#' Generate clustered regions
#'
#' Construct regions of clustered objects
#'
#' The function takes a non-overlapping system of spheres, cylinders, spheroids of type "\code{prolate}" or "\code{oblate}
#' and a random ball configuration \code{CL} of class "\code{spheres}" which defines clustered regions of the objects given
#' in \code{S} according to the minimum required distance \code{eps} between these objects and a number of objects required
#' in that region by \code{minSize} in the list \code{cond}. The function \code{\link{rsa}} is internally called if
#' \code{check_overlap=TRUE} in order to make these cluster regions non-overlapping.
#'
#' @param S				(non-overlapping) geometry objects system
#' @param CL			cluster regions, i.e. objects of class '\code{spheres}', possibly overlapping
#' @param cond			conditioning object for cluster algorithm, see details
#' @param check_overlap	logical, whether spheres are allowed to overlapp, \code{check_overlap=FALSE} (default)
#' @param verbose 		logical, verbose output, \code{verbose=FALSE} (default)
#' @param pl		    integer, \code{pl>0} for printing information
#'
#' @return 				a list of the following element:
#' 					    \itemize{
#'							\item{id}{ numeric vector of ids of including spheroids }
#' 							\item{center}{ the center of enclosing ball }
#' 							\item{r}{ the radius}
#' 							\item{interior}{ equals to \code{0} if any of the enclosed ball objects
#' 								hit the boundaries of the simulation box and equals \code{1} otherwise}
#' 						}
#' 
#' @example inst/examples/densify.R
#' 
#' @author M. Baaske 
#' @rdname simCLuster
#' @export
simCluster <- function(S, CL, cond = list("eps" = 1e-6, "minSize" = 1L),
				check_overlap = FALSE, verbose = FALSE, pl = 0L)
{
	# check arguments for 'cond'
	it <- match(names(cond), c("eps","minSize"))
	if (anyNA(it))
	 stop("Expected 'cond' as named list of arguments: 'eps','minSize'")
		
	#sp <- unfoldr::simPoissonSystem(theta,lam,size="const",type="spheres",
	#		box=box,pl=1,label=as.character("C"))
	
	# apply rsa
	if(check_overlap){
	  CL <- rsa(CL,pl=pl,verbose=verbose)
    }
	
	# call clustering algorithm
	.Call(C_Cluster,CL,S,cond)
}

.update <- function(A) UseMethod(".update",A)
.update.cylinders <- function(A) {
	i <- .checkEn(A,.checkOverlapCylinder)
	if(length(i)>0)
		warning(paste("Overlaps detected in cluster:",paste0("",unlist(i),".",collapse=","),"Consider to increase weight parameter and re-run the densification."))
	
	## update interior
	## this now has to be called by the user!
	
	#ids <- unfoldr::updateIntersections(A)
	
	for(k in 1:length(A)) {
		# update origin0/origin1
		m <- 0.5*A[[k]]$h*A[[k]]$u
		A[[k]]$origin0 <- A[[k]]$center + m
		A[[k]]$origin1 <- A[[k]]$center - m
		#attr(A[[k]],"interior") <- as.logical(ids[k])
	}
	
	#attr(A,"interior") <- all(as.logical(ids))
	return (A)
}

.update.spheroid <- function(A) {
	return ( .update.prolate(A) )
}

.update.oblate <- function(A) {
	return ( .update.prolate(A) )
}

.update.prolate <- function(A) {
	i <- .checkEn(A,.checkOverlap)
	if(length(i)>0)
		warning(paste("Overlaps detected in cluster:",paste0("",unlist(i),".",collapse=","),"Consider to increase weight parameter and re-run the densification."))
	
    ## changes in current version:
	## test on intersection of objects with lateral planes 
	## now only provided by package 'unfoldr' as update interior
	# ids <- unfoldr::updateIntersections(A)
	
	#for(k in 1:length(A))
	# attr(A[[k]],"interior") <- as.logical(ids[k])
	# attr(A,"interior") <- all(as.logical(ids))
	return (A)
}

.update.spheres <- function(A) {
	i <- .checkEn(A,.checkOverlapSphere)
	if(length(i)>0)
	 warning(paste("Overlaps detected in cluster:",paste0("",unlist(i),".",collapse=","),"Consider to increase weight parameter and re-run the densification."))
			
 
    ## see above: update interior
	# ids <- unfoldr::updateIntersections(A)
	# for(k in 1:length(A))
	#  attr(A[[k]],"interior") <- as.logical(ids[k])
	# attr(A,"interior") <- all(as.logical(ids))
	return (A)
}

#' Densify clusters
#'
#' Densification of objects within a predefined ball
#'
#' This densification method requires package\code{\link[GenSA]{GenSA}} to be installed. Additionally for
#' performance reasons it is recommended to make use of parallization by package 'parallel' installed and
#' available or some 'cluster' object derived from it.
#'
#' Parts of the (reinforcement) geometry system are densified in order to account for regions where some
#' objects stick together more close than in other regions of the simulation box. These so-called
#' clusters are generated by a densification algorithm minimizing some energy
#' of object positions and their (overlap) distances within an enclosing ball of objects with
#' a predefined cluster radius. In order to prevent overlaps between objects one possibly has to increase
#' the \code{weight} factor iteratively and run it again if some objects happen to overlap.
#'
#' @param F					non-overlapping particle system
#' @param CL				cluster list, defining the different cluster regions
#' @param ctrl  			control arguments passed to function \code{GenSA}
#' @param weight   			optional, weight factor \code{weight=10} (default)
#' @param info			    optional, logical, if TRUE return result from \code{GenSA} call
#' @param fun 			    optional, either \code{lapply} (default) or parllel processing by \code{mclapply}
#' @param cl 				optional, parallel cluster object
#'
#' @return					a list of all objects (including newly densified regions) named \code{S}
#' 						    and a named list of clustered regions \code{cluster}
#'
#' @example inst/examples/densify.R
#' 
#' @author M. Baaske 
#' @rdname densifyCluster
#' @export
densifyCluster <- function(F, CL, ctrl, weight = 10, info = FALSE, fun = lapply, cl = NULL) {
	if(!requireNamespace("GenSA", quietly=TRUE))
	  stop("package 'GenSA' is required to run this function.")
  	
    # get simulation box
  	box <- attr(F,"box")
  	if(is.null(box))
	  stop("Could not find attribute 'box' in `",as.character(substitute(F)),"`.")
    
  	densify <- function(L,box,ctrl,weight,FUNCTION)
	{
		S <- L$S
		dim <- length(S)
		clust <- L$clust
		maxD <- sapply(S,function(s) sqrt(sum((s$center-clust$center)^2)))
		# apply GenSA
		tryCatch({
			out <- GenSA::GenSA(par=0.5*maxD, lower=rep(0,length(S)), upper=0.9*maxD,
				 fn = FUNCTION, control = ctrl, S = S, clust = clust, weight = weight)
	 			}, error = function(e) {
					structure(e, error=e) 
					}
				)
		
		# do nothing in case of error
		if(inherits(out,"error")) {
			warning("'GenSA' optimization routine failed with error.")		 
		} else {
			for(i in 1:dim)
				# do not move (densify) 2nd. phase
				if(attr(S[[i]],"label") != "F")
					S[[i]]$center <- S[[i]]$center-out$par[i]*((S[[i]]$center-clust$center)/sqrt(sum((S[[i]]$center-clust$center)^2)))
		}
		
		structure(S, "interior"=TRUE,
			"info"=if(info) list("out"=out,"val"=out$value,"maxD"=maxD) else NULL)
	}

	FUN <- NULL
	if(attr(F,"class")=="prolate" || attr(F,"class")=="oblate") {
		FUN <- .EnSpheroid
	}
	else if(attr(F,"class")=="cylinders") {
		FUN <- .EnCylinder
	}
	else if(attr(F,"class")=="spheres") {
		FUN <- .EnSphere
	}
	else { stop("Unknow object type. --> exiting.")}

	N <- length(CL)
	cat("densify clusters of length: ",length(CL), " \n")

	A <- tryCatch(
			{
				L <- lapply(seq_len(N),function(i) list("S"=F[CL[[i]]$id],"clust"=CL[[i]])	)
				if(!is.null(cl)) {
					m <- min(N,length(cl))
					parallel::clusterExport(cl[1:m],list(".getA",".slice",".spheroidTouchDistCoeff",".En",".checkEn"))
					if(any(class(cl) %in% c("MPIcluster"))) {
						parallel::parLapplyLB(cl[seq_len(m)],L,densify,box=box,ctrl=ctrl,weight=weight,FUNCTION=FUN)
					} else {
						do.call(c, parallel::clusterApply(cl[seq_len(m)], .splitList(L, length(cl)), fun, densify,
										box=box,ctrl=ctrl,weight=weight,FUNCTION=FUN))
					}
				} else {
					fun(L,densify,box=box,ctrl=ctrl,weight=weight,FUNCTION=FUN)
				}

			}
			, error = function(e) {
				structure(e,class=c("error","condition"),
						error=simpleError(.makeMessage("Cluster construction failed.\n")))
		 }
	)
	
	if(length(A)>0) {
		# check overlap and update
		cl.name <- class(F)
		for(i in 1:length(A)) {
			class(A[[i]]) <- cl.name		
			A[[i]] <- .update(A[[i]])
		}
		# copy to original system objects
		invisible(lapply(A,
			function(x) {
				# loop over all spheroids
				for(k in 1:length(x)) {
					F[[x[[k]]$id]] <<- x[[k]]
				}
			})
		)
	} else {
		message("There are no particles in the randomly generated clusters. Consider to increase the cluster radius.\n")
	}
	return ( list("S" = F, "cluster" = A ) )
}
