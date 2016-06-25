# get the correct spheroid (rotation-size) matrix
.getA <- function(s) {	
	A <- matrix(c(1.0/s$ab[1]^2,0,0,0,1.0/s$ab[1]^2,0,0,0,1.0/s$ab[2]^2),nrow=3)
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
	return(1 < .spheroidTouchDistCoeff(S1,A1,S2,A2))   
}

# check whether two spherocylinders overlap
.checkOverlapCylinder <- function(q,p) {		
	d <- .C(C_sdm,as.numeric(q$center-p$center),
			as.numeric(q$u),
			as.numeric(p$u),
			as.numeric(q$length*0.5),
			as.numeric(p$length*0.5),
			d=numeric(1))$d
	
	return ( d < ((q$r+p$r)^2 - 1e-3))
}

.checkOverlapSphere <- function(q,p) {		
	return ( sum( (q$center-p$center)^2 ) < (q$r+p$r)^2 - 1e-4 )
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
				as.numeric(p$length*0.5),
				as.numeric(q$length*0.5),
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
#' Simulate regions of clustered objects
#' 
#' The functions generates a random non-overlapping ball configuration
#' by the RSA method with some radius defined by the parameter list \code{theta}.
#' 
#' @param S				(non-overlapping) geometry objects system
#' @param theta 		list of simulation parameter
#' @param lam			the intensity parameter of the underlying Poisson process
#' @param box			the simulation box
#' @param cond			conditioning object for cluster algorithm, see details
#' @param verbose 		logical, \code{verbose=FALSE} (default) additional information
#' @param pl		    integer, \code{pl=0} for no additional information
#' 
#' @return 				a list of the following element:
#' 					    \itemize{
#'							\item{id}{ numeric vector of ids of including spheroids }
#' 							\item{center}{ the center of enclosing ball }
#' 							\item{r}{ the radius}
#' 							\item{interior}{ equals to \code{0} if any of the enclosed ball objects
#' 								hit the boundaries of the simulation box and equals \code{1} otherwise} 
#' 						}
#' @example inst/examples/densify.R
simCluster <- function(S, theta, lam, box, cond = list("eps"=0,"minSize"=1),
							verbose = FALSE, pl = 0)
{		
	# check arguments for 'cond'
	it <- match(names(cond), c("eps","minSize"))
	if (anyNA(it))
		stop("Expected 'cond' as named list of arguments: 'eps','minSize'")
	it <- match(names(theta), c("r"))
	if (anyNA(it))
		stop("Expected 'theta' as list of named arguments.")

	sp <- unfoldr::simSphereSystem(theta,lam,rdist="const",box=box,pl=101,label="C")				
	# apply rsa
	sp <- rsa(sp,box,pl=pl,verbose=verbose)
	# call clustering algorithm
	.Call(C_Cluster,sp,S,cond)
}

.update <- function(A, box) UseMethod(".update",A)
.update.cylinder <- function(A, box) {
	i <- .checkEn(A,.checkOverlapCylinder)
	if(length(i)>0)
		warning(paste("Overlaps detected in cluster ",unlist(i),"
				 after call to 'densify'. Consider to increase the weight parameter as a penalty.",sep=""))
	# update interior	
	ids <- unfoldr::updateIntersections(A,box)	
	for(k in 1:length(A)) {
		# update origin0/origin1
		m <- 0.5*A[[k]]$length*A[[k]]$u
		A[[k]]$origin0 <- A[[k]]$center + m
		A[[k]]$origin1 <- A[[k]]$center - m
		attr(A[[k]],"interior") <- as.logical(ids[k])
	}
	attr(A,"interior") <- all(as.logical(ids))
	return (A)
}

.update.spheroid <- function(A,box) {
	return ( .update.prolate(A,box) ) 
}

.update.oblate <- function(A,box) { 
	return ( .update.prolate(A,box) ) 
}

.update.prolate <- function(A,box) {
	i <- .checkEn(A,.checkOverlap)
	if(length(i)>0)
		warning(paste("Overlaps detected in cluster ",unlist(i),"
  		after call to 'densify'. Consider to increase the weight parameter as a penalty.",sep=""))
	# update interior	
	ids <- unfoldr::updateIntersections(A,box)	
	for(k in 1:length(A))		
	 attr(A[[k]],"interior") <- as.logical(ids[k])	
	attr(A,"interior") <- all(as.logical(ids))
	return (A)
}

.update.sphere <- function(A,box) {
	i <- .checkEn(A,.checkOverlapSphere)
	if(length(i)>0)
		warning(paste("Overlaps detected in cluster ",unlist(i),"
		 after call to 'densify'. Consider to increase the weight parameter as a penalty.",sep=""))
	# update interior	
	ids <- unfoldr::updateIntersections(A,box)	
	for(k in 1:length(A))		
		attr(A[[k]],"interior") <- as.logical(ids[k])	
	attr(A,"interior") <- all(as.logical(ids))
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
#' @param box				simulation box
#' @param ctrl  			control arguments passed to function \code{GenSA}
#' @param weight   			optional: weight factor \code{weight=10} (default)
#' @param info			    optional: logical, if TRUE return result from \code{GenSA} call
#' @param cl 				optional: cluster object, see package \code{snow}
#' @param parallel.option   optional: In case of \code{mclapply} the function \code{mclapply} 
#' 								      is used from the \code{parallel} package
#' 
#' @return					a list of all objects (including newly densified regions) named \code{S}
#' 						    and a named list of clustered regions \code{cluster}
#' 
#' @example inst/examples/densify.R
densifyCluster <- function(F, CL, box, ctrl, weight = 10, info = FALSE, cl = NULL,
					parallel.option=getOption("parallel.option","lapply"))
{	
	if(!requireNamespace("GenSA", quietly=TRUE))
	  stop("package 'GenSA' is required to run this function.")
	densify <- function(L,box,ctrl,weight,FUNCTION) {		
		S <- L$S
		dim <- length(S)
		clust <- L$clust		
		maxD <- sapply(S,function(s) sqrt(sum((s$center-clust$center)^2)))
		# apply GenSA	
		out <- tryCatch({
					GenSA::GenSA(par=0.5*maxD, lower=rep(0,length(S)), upper=0.9*maxD,
				    	 fn = FUNCTION, control = ctrl, S = S, clust = clust, weight = weight)				
					},
				  error=function(e) structure(e, error=e) )
		# do nothing in case of error
		if(is.null(attr(out,"error"))) {
		 for(i in 1:dim)
			 if(attr(S[[i]],"label")!="F")
			   S[[i]]$center <- S[[i]]$center-out$par[i]*((S[[i]]$center-clust$center)/sqrt(sum((S[[i]]$center-clust$center)^2)))
		}
		
		msg <- if(info) list("out"=out,"val"=out$value,"maxD"=maxD) else NULL			
		structure(S, "interior"=TRUE, "info"=msg)
	}
			
	FUN <- NULL
	if(attr(F,"class")=="prolate" || attr(F,"class")=="oblate") {
		FUN <- .EnSpheroid						
	}
	else if(attr(F,"class")=="cylinder") {
		FUN <- .EnCylinder				
	}
	else if(attr(F,"class")=="sphere") {
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
						fun <- get(parallel.option)
						do.call(c, parallel::clusterApply(cl[seq_len(m)], .splitList(L, length(cl)), fun, densify,
										box=box,ctrl=ctrl,weight=weight,FUNCTION=FUN))	
					}
				} else {				
					if(parallel.option=="mclapply" && !requireNamespace("parallel", quietly=TRUE))
					  stop("package 'parallel' is required to run this function in parallel mode.")
				  	fun <- get(parallel.option)
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
			A[[i]] <- .update(A[[i]],box)
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
		message("There are no particles in the randomly generated  clusters.\n
                  Consider to increase the cluster radius.\n")
	}
	return ( list("S" = F, "cluster" = A ) )		
}
