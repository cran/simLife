/**
 * @file cluster.cpp
 * @date 04/26/2016
 *
 * @brief Collect all objects (spheroids,etc.) which lie
 *        in a predefined ball, i.e. the cluster region
 *
 * @author M. Baaske
 */

#include <vector>
#include "cluster.h"


#define NDIM 3
#define DISTANCE_VEC(P1,P2,R)                      \
do {                                               \
  int _i;                                          \
  for (_i = 0; _i < NDIM; _i++)                    \
    (R)[_i] = (P1)[_i]-(P2)[_i];                   \
} while(0)

#define Vec3Op(A,assign_op,B,op,C)                 \
(   A[0] assign_op B[0] op C[0],                   \
    A[1] assign_op B[1] op C[1]                    \
)

#define VecDot(A,B)  ((A[0]*B[0]) + (A[1]*B[1]) + (A[2]*B[2]))
#define VecNorm(V)    sqrt(VecDot(V,V))

extern SEXP getListElement (SEXP list, const char *str);

namespace STGM {

  class CCluster {
   public:
    double *m_p;
    double m_radius;
    int m_interior, m_ncount;

    CCluster(double *center, double radius) :
      m_p(center),  m_radius(radius), m_interior(1), m_ncount(0)
    {};

    void add(int id, int interior, int ncount) {
      m_ncount += ncount;
      id_array.push_back(id);
      if(m_interior)
        m_interior = interior;
    }

    std::vector<int> & getArrayId() { return id_array; }
    unsigned size() const { return id_array.size(); }

  private:
   std::vector<int> id_array;

  };

typedef std::vector<CCluster> cluster_vector;

} // namespace

/**
 * \brief Construct cluster of particles.
 *        The ferrit phase objects are also included in
 *        the return list, but not densified later, because of
 *        being treated as fixed.
 *
 * @param R_spheres        cluster balls
 * @param R_spheroids      spheroids to be clustered
 * @param R_cond           condition object
 *
 * @return R list of clustered objects
 */
SEXP Cluster(SEXP R_spheres, SEXP R_s, SEXP R_cond) {
    int dim=3, ncount=0,
        nspheres=length(R_spheres),
        N=length(R_s);

    double eps = REAL(AS_NUMERIC(getListElement( R_cond, "eps")))[0];
    int minSize = INTEGER(AS_INTEGER(getListElement( R_cond, "minSize")))[0];

    STGM::cluster_vector cluster;
    cluster.reserve(nspheres);

    double r=0.0;
    for(int k=0; k<nspheres; k++) {
        r = REAL(getListElement( VECTOR_ELT(R_spheres, k),"r"))[0];
        cluster.push_back(STGM::CCluster(REAL(getListElement(VECTOR_ELT(R_spheres, k), "center")),std::max(r-eps,0.0)));
    }

    double q[3], *p=0;
    const char *label = "N";

    for(int i=0; i<N; ++i) {
        p = REAL(AS_NUMERIC(getListElement(VECTOR_ELT(R_s,i),"center")));
        label = translateChar(asChar(getAttrib(VECTOR_ELT(R_s,i), install("label"))));
        ncount = (!std::strcmp(label,"P") ? 1 : 0);
        for(int j=0; j<nspheres; ++j) {
            DISTANCE_VEC(p,cluster[j].m_p,q);
            if(VecNorm(q) < cluster[j].m_radius) {
               cluster[j].add(INTEGER(getListElement(VECTOR_ELT(R_s,i),"id"))[0],
                              INTEGER(getAttrib(VECTOR_ELT(R_s,i),install("interior")))[0],
                              ncount);
            }
        }
    }

    ncount = 0;
    size_t nclust = cluster.size();
    for(size_t k=0; k<nclust; k++) {
        if(cluster[k].m_ncount > minSize)
          ++ncount;
    }

    SEXP R_cluster = R_NilValue;
    PROTECT(R_cluster = allocVector(VECSXP,ncount));

    SEXP names;
    PROTECT(names = allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("id"));
    SET_STRING_ELT(names, 1, mkChar("center"));
    SET_STRING_ELT(names, 2, mkChar("r"));
    SET_STRING_ELT(names, 3, mkChar("interior"));

    int i=0;
    SEXP R_tmp, R_ctr, R_ids = R_NilValue;
    for(size_t k=0; k<nclust; k++)
    {
        if(cluster[k].m_ncount>minSize)
        {
          std::vector<int> &ids = cluster[k].getArrayId();
          PROTECT(R_tmp = allocVector(VECSXP,4));
          PROTECT(R_ids = allocVector(INTSXP,ids.size()));
          PROTECT(R_ctr = allocVector(REALSXP,dim));
          Memcpy(INTEGER(R_ids),&ids[0],ids.size());
          Memcpy(REAL(R_ctr),cluster[k].m_p,dim);

          SET_VECTOR_ELT(R_tmp,0,R_ids);
          SET_VECTOR_ELT(R_tmp,1,R_ctr);
          SET_VECTOR_ELT(R_tmp,2,ScalarReal(cluster[k].m_radius));
          SET_VECTOR_ELT(R_tmp,3,ScalarInteger(cluster[k].m_interior));
          setAttrib(R_tmp, R_NamesSymbol, names);
          SET_VECTOR_ELT(R_cluster,i,R_tmp);

          ++i;
          UNPROTECT(3);
        }
    }

    UNPROTECT(2);
    return R_cluster;
  }
