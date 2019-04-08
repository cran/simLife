/**
*  @file sim2Life.cpp
 * @date 04-20-2016
 *
 * @brief R interface to object projection methods,
 *        convex hull and the simulation of defect accumulation
 *
 * @author M. Baaske
 */

#include <R_ext/Rdynload.h>

#include "simLife.h"
#include "cluster.h"
#include "GeometricPrimitives.h"

static int PL = 0;

/// extern declarations
extern SEXP getVar(SEXP name, SEXP rho);

extern SEXP convert_R_Ellipse2(STGM::CEllipse2 &ellipse);
extern STGM::Ellipses2 convert_C_Ellipses2(SEXP R_ellipses);
extern STGM::Spheroids convert_C_Spheroids(SEXP R_spheroid);
extern STGM::Spheres convert_C_Spheres(SEXP R_spheres);
extern STGM::Cylinders convert_C_Cylinders(SEXP R_cylinders);


/**
 * @brief Calculate a vector of points on the convex hull
 *        Adopted from
 *        https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
 *
 * @param P     point vector of points to construct the convex hull
 * @return      points that are the convex hull
 */
inline double cross3d(const STGM::CPoint2d &x, const STGM::CPoint2d &y, const STGM::CPoint2d &z) {
  return (y[0]-x[0]) * (z[1]-x[1]) - (y[1]-x[1]) * (z[0]-x[0]);
}

STGM::PointVector2d convexHull2d(STGM::PointVector2d P) {
    int n = P.size(), k = 0;
    STGM::PointVector2d H(2*n);
    std::sort(P.begin(), P.end());
    // construct the lower convex hull
    for (int i = 0; i < n; ++i) {
            while (k >= 2 && cross3d(H[k-2], H[k-1], P[i]) <= 0) k--;
            H[k++] = P[i];
    }
    // and the upper part
    for (int i = n-2, t = k+1; i >= 0; i--) {
       while (k >= t && cross3d(H[k-2], H[k-1], P[i]) <= 0) k--;
       H[k++] = P[i];
    }
    H.resize(k-1);
    return H;
}

/**
 * @brief Area of convex hull as the area of
 *        its approximating polygon
 *
 * @param P     point vector of convex hull
 * @return      area of convex hull
 */

double convHArea(const STGM::PointVector2d &P) {
  double area = 0.0;
  int np = P.size(), i = 0, j = np-1;

  for (i=0; i<np; i++) {
     area += (P[j][0]+P[i][0])*(P[j][1]-P[i][1]);
     j = i;
  }
  return (area < 0 ? -area : area) * .5;
}

/**
 * @brief  Construct convex hull and
 *         calculate area of polygon
 *
 * @param R_points      points to construct convex hull for
 * @return              convex hull points and the area
 */
extern "C"
SEXP convexHull(SEXP R_points) {
  int n = LENGTH(R_points);

  STGM::PointVector2d P;
  P.reserve(n);

  for(int i=0;i<n; i++) {
     P.push_back(STGM::CPoint2d(REAL(VECTOR_ELT(R_points,i))));
  }

  STGM::PointVector2d H = convexHull2d(P);
  double area = convHArea(H);

  SEXP R_H, R_p;
  size_t m = H.size();

  PROTECT(R_H = allocVector(VECSXP,m));
  for(size_t i=0; i<H.size(); ++i) {
      PROTECT(R_p = allocVector(REALSXP,2));
      REAL(R_p)[0] = P[i][0];
      REAL(R_p)[1] = P[i][1];
      SET_VECTOR_ELT(R_H,i,R_p);
      UNPROTECT(1);
  }
  setAttrib(R_H, install("area"), ScalarReal(area));

  UNPROTECT(1);
  return R_H;
}

/**
 * @brief Get the points on the border of
 *        each projected spheroid i.e. an ellipse.
 *
 * @param R_ellipses
 * @param R_n
 * @return matrix of points
 */
SEXP GetPointsForConvexHull(SEXP R_ellipses, SEXP R_n) {
  int i=0, k=0, N = length(R_ellipses),
      n = INTEGER(AS_INTEGER(R_n))[0];
  int m = N*n;

  SEXP R_points = R_NilValue;
  PROTECT(R_points = allocMatrix(REALSXP,m,2));
  double *mat = REAL(R_points);

  STGM::Ellipses2 ellipses = convert_C_Ellipses2(R_ellipses);
  STGM::CPoint2d p;
  double t=0.0, s=2.0*M_PI / (double)n;
    for(i=0; i<N; i++) {
          for(k=0,t=0; k<n; k++) {
              p = ellipses[i].PointOnEllipse(t);
              mat[k+i*n] = p[0];
              mat[k+i*n+m] = p[1];
              t += s;
          }
  }

  UNPROTECT(1);
  return R_points;
}

/**
 * @brief Calculate spheroid projection according
 *        to the given crack types
 *
 * @param R_spheroids
 * @param R_crack_type
 * @return Spheroid projections
 */
SEXP GetSpheroidProjection(SEXP R_spheroids, SEXP R_crack_type) {
  STGM::Spheroids spheroids = convert_C_Spheroids(R_spheroids);

  int n = spheroids.size();
  SEXP R_ret = R_NilValue;
  PROTECT(R_ret = allocVector(VECSXP,n));

  STGM::CEllipse2 ellipse;
  for(int i=0;i<n; i++) {
        spheroids[i].setCrackType( INTEGER(AS_INTEGER(R_crack_type))[i] );
        ellipse = spheroids[i].spheroidProjection();
        SET_VECTOR_ELT(R_ret,i,convert_R_Ellipse2(ellipse));
  }
  UNPROTECT(1);
  return R_ret;
}

/**
 * @brief Get the points of the projected cylinders
 *
 * @param R_cylinders           list cylinders
 * @param R_crack_type          vector of crack types {0=crack,1=delam}
 * @param R_np                  sample points for each cylinder
 * @return                      list of matrices, each with attribute area
 */
SEXP GetSphereProjection(SEXP R_s, SEXP R_np) {
  STGM::Spheres spheres = convert_C_Spheres(R_s);
  int n = spheres.size(),
      np = INTEGER(AS_INTEGER(R_np))[0];

  double area = 0;
  SEXP R_ret, R_p;
  PROTECT(R_ret = allocVector(VECSXP,n));
  for(int i=0; i<n; ++i) {
      PROTECT(R_p = allocMatrix(REALSXP,np,2));
      STGM::PointVector2d points;
      points.reserve(np);

      area = spheres[i].projectedPointsWithArea(points,np);
      for(int k=0; k<np; ++k) {
          REAL(R_p)[k] = points[k][0];
          REAL(R_p)[k+np] = points[k][1];
      }
      SET_VECTOR_ELT(R_ret,i,R_p);
      setAttrib(R_p, install("area"), ScalarReal(area));
      UNPROTECT(1);
  }

  UNPROTECT(1);
  return R_ret;
}

/**
 * @brief Get the points of the projected cylinders
 *
 * @param R_cylinders           list cylinders
 * @param R_crack_type          vector of crack types {0=crack,1=delam}
 * @param R_np                  sample points for each cylinder
 * @return                      list of matrices, each with attribute area
 */
SEXP GetCylinderProjection(SEXP R_cylinders, SEXP R_crack_type,SEXP R_np) {
  STGM::Cylinders cylinders = convert_C_Cylinders(R_cylinders);
  int n = cylinders.size(), type = 0, m = 0,
	  np = INTEGER(AS_INTEGER(R_np))[0];

  double area = 0;
  SEXP R_ret, R_p;
  PROTECT(R_ret = allocVector(VECSXP,n));
  for(int i=0; i<n; ++i) {
      type=INTEGER(AS_INTEGER(R_crack_type))[i];
      m = (type>0 ? MIN(20,np) : np);
      PROTECT(R_p = allocMatrix(REALSXP,m,2));
      STGM::PointVector2d points;
      points.reserve(m);

      cylinders[i].setCrackType(type);
      area = cylinders[i].projectedPointsWithArea(points,m);
      for(int k=0; k<m; ++k) {
          REAL(R_p)[k] = points[k][0];
          REAL(R_p)[k+m] = points[k][1];
      }
      SET_VECTOR_ELT(R_ret,i,R_p);
      setAttrib(R_p, install("area"), ScalarReal(area));
      setAttrib(R_p, install("type"), ScalarReal(type));
      UNPROTECT(1);
  }

  UNPROTECT(1);
  return R_ret;
}

#define PRINT_HEAD_INFO {                               \
  Rprintf("\n");                                        \
  Rprintf(" %4s \t %4s \t %8s \t %12s" ,                \
          "Size", "Iter", "Interior", "Area" );         \
  Rprintf("\n");                                        \
}

#define PRINT_INFO(SIZE,K,I,A) {                                   \
  Rprintf(" %4d \t %4d \t %8d \t %12.4e",(SIZE),(K),(I),(A));      \
  Rprintf("\n");                                                   \
}


typedef struct {
  double distTol,
         areaMax,
         areaIn,
         areaOut,
		 Tmax;
} siminfo_t;



template<typename T>
void intern_simDefect( typename STGM::ClusterList<T>::Type &cl,
                          STGM::Converter< STGM::ConverterFunction<T> > &converter, siminfo_t &info ) {

  typename STGM::ClusterList<T >::iterator_t  jt, endit;

  Rboolean stopit = FALSE;
  STGM::CDefect<T> *head, *last;

  double minDist=0.0;
  /* // to compare to minDist (distTol is weight factor < 1) */
  double MPI4 = info.distTol*std::sqrt(M_PI_4);
  int i=1, k=0, nclust=0, N=converter.N;

  head = converter(0);
  head->project();
  cl.push_back(head);

  while(i < N && !stopit)
  {
      head = converter(i);

      // check for maximum time i.e. treat object as a runout
      if(head->m_time > info.Tmax) {
    	  stopit=TRUE;
    	  break;
      }

      cl.push_back(head);
      // project single object, store points
      head->project();
      // check if projected defect of single
      // object is already large enough
	  if(head->m_inner && head->m_area > info.areaIn) {
		  stopit = TRUE;
		  head->m_broken = 1;
		  Rprintf("Single particle exceeds critical area (internal).\n");
		  break;
	  } else if(!head->m_inner && head->m_area > info.areaOut) {
		  stopit = TRUE;
		  head->m_broken = 1;
		  Rprintf("Single particle exceeds critical area (surface).\n");
		  break;
	  }

      //T &scmp = head->m_object;
      jt = cl.begin(); endit = cl.end(); --endit;

      //TODO: probably this loop could be run in parallel!
      for(k = 0; jt != endit; ++k )
      {
         last = *jt;
         minDist = last->descent(head);
          //Rprintf("d: %f, head: %f, last: %f \n",minDist,std::sqrt(head->m_area),std::sqrt(last->m_area));
          // check minimum distance
          if(minDist < MPI4*MIN(std::sqrt(head->m_area),std::sqrt(last->m_area))) {
              // update points
              head->update(last);
              // append nodes
              head->append(last);
              // erase current element
              jt = cl.erase(jt);

              // get convex hull points and polygon area
              STGM::PointVector2d H = convexHull2d(head->m_points);
              head->m_area = convHArea(H);

              if(head->m_area > info.areaMax)
                info.areaMax = head->m_area;
              if( (head->m_inner && head->m_area > info.areaIn) ||
                  (!head->m_inner && head->m_area > info.areaOut) ) {
                    stopit = TRUE;
                    break;
              }

          } else { // minDist
              ++jt;
         }
     } // end for

     // give new cluster id to current head
     if(head->m_size > 1 || stopit) {
         head->m_num = ++nclust;
     }

    if(PL > 100) {
        PRINT_HEAD_INFO
        PRINT_INFO(head->m_size,k,head->m_inner,head->m_area)
    }
    ++i;
  } // end while

}

template<typename T>
SEXP convert_R_result(typename STGM::ClusterList<T>::Type &cl, const siminfo_t &info) {
  typename STGM::ClusterList<T >::iterator_t  jt;
  int nclust=0, inner=0;
  // construct R elements
  const char *nms[] = {"id", "n", "B", "interior", "A", "inner", "T", "label", ""};

  STGM::CDefect<T> *p, *head;
  int m = 0, j = 0, l = 0, maxSize = 0;

  if(PL>100) {
	  for(jt = cl.begin(); jt != cl.end(); ++jt) {
	          if( (*jt)->m_size > 1 || (*jt)->m_broken)
	            ++nclust;
	  }
  } else { nclust=1; }

  SEXP R_cl = R_NilValue;
  PROTECT(R_cl = allocVector(VECSXP,nclust));

  SEXP R_tmp, R_id, R_type, R_interior,
       R_label, R_num, R_area, R_time;

  if(PL > 100)							/* return all accumulated clusters */
  {
    for(jt = cl.begin(); jt != cl.end(); ++jt)
    {
        p = *jt;
        m = p->m_size;
        if(!(m>1) && !p->m_broken) {
            delete p;  					/*  calling first destructor p->~CDefect() automatically */
            continue;
        }
        if(m > maxSize)
          maxSize = m;

        PROTECT(R_id = allocVector(INTSXP,m));
        PROTECT(R_num = allocVector(INTSXP,m));
        PROTECT(R_type = allocVector(INTSXP,m));
        PROTECT(R_interior = allocVector(INTSXP,m));
        PROTECT(R_label = allocVector(STRSXP,m));
        PROTECT(R_area = allocVector(REALSXP,m));
        PROTECT(R_time = allocVector(REALSXP,m));

        head = p;
        inner = head->m_inner;

        l=0;
        while(p != 0) {
            INTEGER(R_id)[l] = p->m_id;
            INTEGER(R_num)[l] = p->m_num;
            INTEGER(R_type)[l] = p->m_type;
            INTEGER(R_interior)[l] = p->m_interior;
            REAL(R_area)[l] = p->m_area;
            REAL(R_time)[l] = p->m_time;
            SET_STRING_ELT(R_label,l,mkChar(p->m_label));
            p = p->next;
            ++l;
        }
        delete head;

        PROTECT(R_tmp = mkNamed(VECSXP, nms));
        SET_VECTOR_ELT(R_tmp,0,R_id);
        SET_VECTOR_ELT(R_tmp,1,R_num);
        SET_VECTOR_ELT(R_tmp,2,R_type);
        SET_VECTOR_ELT(R_tmp,3,R_interior);
        SET_VECTOR_ELT(R_tmp,4,R_area);
        SET_VECTOR_ELT(R_tmp,5,ScalarInteger(inner));
        SET_VECTOR_ELT(R_tmp,6,R_time);
        SET_VECTOR_ELT(R_tmp,7,R_label);

        SET_VECTOR_ELT(R_cl,j,R_tmp);
        ++j; 								/* index for R result vector (SET_VECTOR_ELT) */
        UNPROTECT(8);
    }

  } else {									/* return only last (final) cluster */

      p = cl.back();
      m = p->m_size;
      maxSize = m;

      // further list objects
      PROTECT(R_id = allocVector(INTSXP,m));
      PROTECT(R_num = allocVector(INTSXP,m));
      PROTECT(R_type = allocVector(INTSXP,m));
      PROTECT(R_interior = allocVector(INTSXP,m));
      PROTECT(R_label = allocVector(STRSXP,m));
      PROTECT(R_area = allocVector(REALSXP,m));
      PROTECT(R_time = allocVector(REALSXP,m));

      l = 0; 								/* p is last defect */
      inner = p->m_inner;
      while(p != 0) {
          INTEGER(R_id)[l] = p->m_id;
          INTEGER(R_num)[l] = p->m_num;
          INTEGER(R_type)[l] = p->m_type;
          INTEGER(R_interior)[l] = p->m_interior;
          REAL(R_area)[l] = p->m_area;
          REAL(R_time)[l] = p->m_time;
          SET_STRING_ELT(R_label,l,mkChar(p->m_label));
          p = p->next;
          ++l;
      }

      for(jt = cl.begin(); jt != cl.end(); ++jt)
        delete (*jt);

      PROTECT(R_tmp = mkNamed(VECSXP, nms));
      SET_VECTOR_ELT(R_tmp,0,R_id);
      SET_VECTOR_ELT(R_tmp,1,R_num);
      SET_VECTOR_ELT(R_tmp,2,R_type);
      SET_VECTOR_ELT(R_tmp,3,R_interior);
      SET_VECTOR_ELT(R_tmp,4,R_area);
      SET_VECTOR_ELT(R_tmp,5,ScalarInteger(inner));
      SET_VECTOR_ELT(R_tmp,6,R_time);
      SET_VECTOR_ELT(R_tmp,7,R_label);

      SET_VECTOR_ELT(R_cl,0,R_tmp);
      UNPROTECT(8);
  }

  SEXP R_maxSize, R_areaMax;
  PROTECT(R_maxSize = allocVector(REALSXP,1));
  PROTECT(R_areaMax = allocVector(REALSXP,1));
  REAL(R_maxSize)[0] = maxSize;
  REAL(R_areaMax)[0] = info.areaMax;
  setAttrib(R_cl, install("maxSize"),R_maxSize);
  setAttrib(R_cl, install("areaMax"),R_areaMax);
  UNPROTECT(2);

  //setAttrib(R_cl, install("aIn"), ScalarReal(info.areaIn));
  //setAttrib(R_cl, install("aOut"), ScalarReal(info.areaOut));
  //setAttrib(R_cl, install("areaMax"), ScalarReal(info.areaMax));

  UNPROTECT(1);
  return R_cl;
}


SEXP SimDefect(SEXP R_vname, SEXP R_clust, SEXP R_dist, SEXP R_areaIn, SEXP R_areaOut,
		 SEXP R_Tmax, SEXP R_print_level, SEXP R_env)
{
    if(isNull(R_env) || !isEnvironment(R_env))
      error(_("Should provide environment for function evaluation."));

    SEXP Rs = R_NilValue;
    PROTECT(Rs = getVar(AS_CHARACTER(R_vname),R_env));
    PL = INTEGER(AS_INTEGER(R_print_level))[0];

    if (TYPEOF(Rs) == PROMSXP) {
      Rs = eval(Rs, R_env);
    } else {
    	error(_("Expression does not evaluate to a promise."));
    }

    siminfo_t info = {REAL(R_dist)[0],0,REAL(R_areaIn)[0],
                      REAL(R_areaOut)[0],REAL(R_Tmax)[0]};

    /* get class */
	SEXP Rclass = R_NilValue;
	PROTECT(Rclass = getAttrib( Rs, R_ClassSymbol));
	const char* name = CHAR( STRING_ELT(Rclass,0));
	UNPROTECT(1);

    if( !std::strcmp(name, "prolate" ) ||
    	!std::strcmp(name, "oblate" ))
    {
        STGM::ClusterList<STGM::CSpheroid>::Type cl;
        STGM::Converter<STGM::ConverterFunction<STGM::CSpheroid> > converter(Rs, R_clust);
        intern_simDefect<STGM::CSpheroid>(cl,converter,info);

        UNPROTECT(1);
        return convert_R_result<STGM::CSpheroid>(cl,info);

    } else if(!std::strcmp(name, "cylinder" )) {
        STGM::ClusterList<STGM::CCylinder>::Type cl;
        STGM::Converter<STGM::ConverterFunction<STGM::CCylinder> > converter(Rs, R_clust);
        intern_simDefect<STGM::CCylinder>(cl,converter,info);

        UNPROTECT(1);
        return convert_R_result<STGM::CCylinder>(cl,info);

    } else if(!std::strcmp(name, "sphere" )) {
        STGM::ClusterList<STGM::CSphere>::Type cl;
        STGM::Converter<STGM::ConverterFunction<STGM::CSphere> > converter(Rs, R_clust);
        intern_simDefect<STGM::CSphere>(cl,converter,info);

        UNPROTECT(1);
        return convert_R_result<STGM::CSphere>(cl,info);

    } else {
    	error(_("Unknown class object."));
    }

    return R_NilValue;
}


/**\brief Calculate spheroid projection only the area
 *
 * @param R_spheroids
 * @return numeric vector, the areas
 */
SEXP GetSpheroidOnlyProjectionArea(SEXP R_spheroids) {
  STGM::Spheroids spheroids = convert_C_Spheroids(R_spheroids);

  SEXP R_ret;
  size_t n = spheroids.size();
  PROTECT(R_ret = allocVector(REALSXP,n));

  for(size_t i=0; i<n; i++) {
    STGM::CEllipse2 ellipse=spheroids[i].delamProjection();
    REAL(R_ret)[i]=M_PI*ellipse.a()*ellipse.b();
  }

  UNPROTECT(1);
  return R_ret;
}


SEXP GetSpheroidBothProjection(SEXP R_spheroids) {
  STGM::Spheroids spheroids = convert_C_Spheroids(R_spheroids);

  SEXP R_ret = R_NilValue;
  PROTECT(R_ret = allocVector(VECSXP,2));

  SEXP R_tmp0, R_tmp1;
  size_t n = spheroids.size();
  PROTECT(R_tmp0 = allocVector(VECSXP,n));
  PROTECT(R_tmp1 = allocVector(VECSXP,n));

  STGM::CEllipse2 ellipse;
  for(size_t i=0;i<n; i++) {
    ellipse = spheroids[i].delamProjection();
    SET_VECTOR_ELT(R_tmp0,i,convert_R_Ellipse2(ellipse));

    ellipse = spheroids[i].crackProjection();
    SET_VECTOR_ELT(R_tmp1,i,convert_R_Ellipse2(ellipse));

  }

  SET_VECTOR_ELT(R_ret,0,R_tmp0);
  SET_VECTOR_ELT(R_ret,1,R_tmp1);

  UNPROTECT(3);
  return R_ret;
}

/* R Interface functions  */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

R_NativePrimitiveArgType myC_t[] = { REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP };
R_NativePrimitiveArgType myC_t2[] = { REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP };

static R_CMethodDef CEntries[]  = {
	{"sdm", (DL_FUNC) &sdm, 6, myC_t},
	{"ContactRadius", (DL_FUNC) &ContactRadius, 8, myC_t2},
    {NULL, NULL, 0, NULL}
};

static R_CallMethodDef CallEntries[] = {
        CALLDEF(GetPointsForConvexHull,2),
        CALLDEF(GetSpheroidProjection,2),
        CALLDEF(GetCylinderProjection,3),
        CALLDEF(GetSphereProjection,2),
        CALLDEF(Cluster,3),
        CALLDEF(SimDefect,8),
        CALLDEF(GetSpheroidBothProjection,1),
        CALLDEF(GetSpheroidOnlyProjectionArea,1),
      {NULL, NULL, 0}
};


void R_init_simLife(DllInfo *info) {
  R_registerRoutines(info, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

/*
void R_unload_simLife(DllInfo *info){
  // nothing
}
*/
