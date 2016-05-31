# check periodic overlap
.checkOverlapPeriodic <- function(S1,S2,bx,by,bz,D,FUN) {
  Q <- S2
  overlap <- FALSE
  for(dx in c(-bx,0,bx)) {
    for(dy in c(-by,0,by)) {
      for(dz in c(-bz,0,bz)) {
        Q$center <- S2$center+c(dx,dy,dz)
        if(sum((S1$center-Q$center)^2)<D^2) {
          overlap <- FUN(S1,Q)
        }
      }
    }
  }
  return(overlap)
}

#' Random sequential adsorption (RSA)
#' 
#' A simple RSA algorithm
#' 
#' The function generates a non-overlapping configuration of spheres, spheroids or spherocylinders
#' with respect to periodic boundary conditions by sequentially adding the given objects at new random 
#' positions in the simulation box while keeping their sizes and orientations. If there 
#' is an overlap detected the position is rejected and some new random position is 
#' generated until all particles have put inside or the maximum number of allowed 
#' iterations is reached. This function is most suited for volume fractions less than 0.15. 
#' 
#' @param S 			overlapping geometry system
#' @param box			simulation box
#' @param F 			secondary phase objects as list
#' @param pl			integer: if \code{pl>0} give some intermediate results 
#' @param verbose       logical: if \code{verbose=TRUE} (default) show additional information
#' 
#' @return			    list of non-overlapping particles
#' @author				Felix Ballani, Markus Baaske
#' @references
#'	 J.W. Evans. Random and cooperative sequential adsorption. Rev. Mod. Phys., 65: 1281-1304, 1993.
rsa <- function(S, box, F = NULL, pl = 0, verbose = TRUE) UseMethod("rsa", S)

#' @rdname rsa
#' @method rsa oblate
#' @S3method rsa oblate
rsa.oblate <- function(S, box, F = NULL, pl = 0, verbose = TRUE) {
	rsa.prolate(S,box,F,pl,verbose)
}

#' @rdname rsa
#' @method rsa prolate
#' @S3method rsa prolate
rsa.prolate <- function(S, box, F = NULL, pl = 0, verbose = TRUE) {
	# combine spheroids and ferrit particles	
	if(!is.null(F)) {
		if(class(F)!=class(S))
			stop(paste0("Trying to pack objects of different types -> exiting."))
		n <- length(S)
		for(i in 1:length(F))
			F[[i]]$id <- i+n
		attrs <- attributes(S)
		S <- c(S,F)
		attributes(S) <- attrs		
	}
	# check of volume fraction
	v <- sum(sapply(S,function(x) x$ab[1]^2*x$ab[2])) 
	p <- 4*pi/3*v/((box$xrange[2]-box$xrange[1])*(box$yrange[2]-box$yrange[1])*(box$zrange[2]-box$zrange[1]))
	if(verbose) 
	  cat(paste("Volume fraction: ",p,"\n"))
	if(p > 0.15)
	  message(paste("Target volume fraction is ",p," which is quite large for RSA algorithm.\n
		Algorithm may terminate successfully (but very slowly) or even fail totally.\n",sep=""))
	
	D <- 2.0*max(unlist(lapply(S, function(s) max(s$ab))))           # overall maximum axis length
	FUN <- .checkOverlap
	
	return (.rsaPeriodic(S,box,FUN,D,pl))
}

#' @rdname rsa
#' @method rsa cylinder
#' @S3method rsa cylinder
rsa.cylinder <- function(S, box, F = NULL, pl = 0, verbose = TRUE) {	
	v <- sum(sapply(S, function(x) pi*x$r^2*x$length + 4/3*pi*x$r^3 ))	
	if(!is.null(F)) {
		if(!(class(F) %in% c("cylinder","sphere")))			
		 stop(paste0("Trying to pack objects of different types -> exiting."))
	 	n <- length(S)	 		
	 	if(class(F) == "sphere") {
		  for(i in 1:length(F)) {
			 F[[i]]$id <- i+n			
			 # dummies to treat sphere as cylinder
			 F[[i]]$u <- c(0,0,1)
			 F[[i]]$length <- 0
			 F[[i]]$rotM <- as.matrix(diag(rep(1,3)))
		 }	
		} else {
			for(i in 1:length(F)) 
			  F[[i]]$id <- i+n
		}
		
		v <- v + sum(sapply(F, function(x) 4/3*pi*x$r^3 ))
		attrs <- attributes(S)
		S <- c(S,F)
		attributes(S) <- attrs
	}	
	p <- v/((box$xrange[2]-box$xrange[1])*(box$yrange[2]-box$yrange[1])*(box$zrange[2]-box$zrange[1]))	
	if(verbose)
	  cat(paste("Volume fraction: ",p,"\n"))
	if(p > 0.15)
		message(paste("Target volume fraction is ",p," which is quite large for RSA method.\n
			Algorithm may terminate successfully (but very slowly) or even fail totally.\n",sep=""))
	
	D <- max(unlist(lapply(S, "[[","length")))           	
	FUN <- .checkOverlapCylinder
	 
	S <- .rsaPeriodic(S,box,FUN,D,pl)
	for(i in seq_along(S)) {		
		x <- S[[i]]
		S[[i]]$origin0 <- x$center + 0.5*x$length*x$u
		S[[i]]$origin1 <- x$center - 0.5*x$length*x$u
		
	}
	return ( S )		
}

#' @rdname rsa
#' @method rsa sphere
#' @S3method rsa sphere
rsa.sphere <- function(S, box, F = NULL, pl = 0, verbose = TRUE) {
	if(!is.null(F)) {
		if(class(F)!=class(S))
			stop(paste0("Trying to pack objects of different types -> exiting."))
		n <- length(S)
		for(i in 1:length(F))
			F[[i]]$id <- i+n
		attrs <- attributes(S)
		S <- c(S,F)
		attributes(S) <- attrs
	}
	v <- sum(sapply(S, function(x) 4/3*pi*x$r^3 ))
	p <- v/((box$xrange[2]-box$xrange[1])*(box$yrange[2]-box$yrange[1])*(box$zrange[2]-box$zrange[1]))	
	if(verbose)
		cat(paste("Volume fraction: ",p,"\n"))
	if(p > 0.15)
		message(paste("Target volume fraction is ",p," which is quite large for RSA method.\n
			Alorithm may terminate successfully (but very slowly) or even fail totally.\n",sep=""))
	
	D <- 2.0*max(unlist(lapply(S, "[[","r")))
	FUN <- .checkOverlapSphere
	
	return (.rsaPeriodic(S,box,FUN,D,pl) )
}


.rsaPeriodic <- function(S, box, FUN, D, pl = 0) {  
  n <- length(S)
  xlim <- box[[1]]
  ylim <- box[[2]]
  zlim <- box[[3]]
  bx <- xlim[2]-xlim[1]
  by <- ylim[2]-ylim[1]
  bz <- zlim[2]-zlim[1]
      
  nx <- ceiling(bx/D)		                                       # lattices dimensions
  ny <- ceiling(by/D)
  nz <- ceiling(bz/D)
  lattice <- array(vector("list",nx*ny*nz),c(nx,ny,nz))
  
  for(i in 1:n) {
    add <- FALSE
    nIter <- 0
    nIterMax <- 10000 # could also depend increasingly on i, e.g. nIterMax <- 2*i+100
    
	while(!add) {
      if(nIter > nIterMax)
		stop(paste0("Maximum number of iterations reached.
              Placed ",i," of ", n," spheroids."))	
      nIter <- nIter+1
      S[[i]]$center[1] <- runif(1,xlim[1],xlim[2])		             # random centre coordinates
      S[[i]]$center[2] <- runif(1,ylim[1],ylim[2])
      S[[i]]$center[3] <- runif(1,zlim[1],zlim[2])
      
	  cx <- floor(S[[i]]$center[1]/D+1)		                         # lattice address of centre coordinates
      cy <- floor(S[[i]]$center[2]/D+1)
      cz <- floor(S[[i]]$center[3]/D+1)
      
	  if(is.null(lattice[cx,cy,cz][[1]])) {
	    add <- TRUE
	  } else {
        id <- lattice[cx,cy,cz][[1]]                                 # ids of candidates which may overlap
        overlap <- FALSE
        for(j in id) {
          if(.checkOverlapPeriodic(S[[i]],S[[j]],bx,by,bz,D,FUN)) {
            overlap <- TRUE
            break
          }
        }
        add <- !overlap
      }
	  
      if(add) { # insert in lattice
        if(cx==1) hx<-c(nx,1,2) else if (cx==nx) hx<-c(nx-1,nx,1) else hx<-c(cx-1,cx,cx+1)
        if(cy==1) hy<-c(ny,1,2) else if (cy==ny) hy<-c(ny-1,ny,1) else hy<-c(cy-1,cy,cy+1)
        if(cz==1) hz<-c(nz,1,2) else if (cz==nz) hz<-c(nz-1,nz,1) else hz<-c(cz-1,cz,cz+1)
        for(jx in hx)
	      for(jy in hy)
		    for(jz in hz)
			 lattice[jx,jy,jz][[1]] <- c(lattice[jx,jy,jz][[1]],i)
        
 		
		if(pl>0) message(paste(i," out of ",n," successfully placed.",sep=""))
      }
    }
  }
  return(S)
}
