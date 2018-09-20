#ifndef SRC_CLUSTER_H_
#define SRC_CLUSTER_H_

#include <R.h>
#include <Rdefines.h>


#ifdef __cplusplus
extern "C" {
#endif

  SEXP Cluster(SEXP R_spheres, SEXP R_spheroids, SEXP R_cond);

#ifdef __cplusplus
}
#endif


#endif
