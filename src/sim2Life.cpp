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

#include "sim2Life.h"
#include "cluster.h"

static int PL = 0;

/// extern declarations
extern SEXP getVar(SEXP name, SEXP rho);
extern SEXP convert_C2R_ellipses(STGM::Ellipses &ellipses);
extern STGM::Ellipses convert_C_Ellipses(SEXP R_ellipses);
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
    sort(P.begin(), P.end());
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
  SEXP Rp;
  int n = length(R_points);
  STGM::PointVector2d P;
  P.reserve(n);

  for(int i=0;i<n; i++) {
      Rp = VECTOR_ELT(R_points,i);
      P.push_back( STGM::CPoint2d(REAL(Rp)[0],REAL(Rp)[1]) );
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
      n = asInteger(AS_INTEGER(R_n));
  int m = N*n;

  SEXP R_points = R_NilValue;
  PROTECT(R_points = allocMatrix(REALSXP,m,2));
  double *mat = REAL(R_points);

  STGM::Ellipses ellipses = convert_C_Ellipses(R_ellipses);
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
  STGM::Ellipses ellipses;
  for(int i=0;i<n; i++) {
        spheroids[i].setCrackType( INTEGER(AS_INTEGER(R_crack_type))[i] );
        ellipses.push_back(spheroids[i].spheroidProjection());
  }
  return convert_C2R_ellipses(ellipses);
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
      np = asInteger(AS_INTEGER(R_np));

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
  int n = cylinders.size(),
      type = 0, m = 0, np = asInteger(AS_INTEGER(R_np));

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

/**\brief Get the minimum distance of two spheroids
 *        approximated by sphero-cylinder distance
 *
 *        Comment: Used for UpdateIntersections
 *
 * @param R_spheroids
 * @param R_cmp
 * @param R_crack_type
 * @param R_crack_cmp
 * @return approximate minimum distance
 *
SEXP GetMinimumDistance(SEXP R_spheroids, SEXP R_cmp, SEXP R_crack_type, SEXP R_crack_cmp) {
  double dist=0, mindist=HUGE_VAL;

  STGM::CSpheroid scmp = convert_C_Spheroid(R_cmp);
  scmp.setCrackType(asInteger(AS_INTEGER(R_crack_cmp)));

  STGM::Spheroids spheroids = convert_C_Spheroids(R_spheroids);

  int n = spheroids.size();
  for(int i=0;i<n; i++) {
       spheroids[i].setCrackType(INTEGER(AS_INTEGER(R_crack_type))[i]);
       dist = spheroids[i].spheroidDistanceAsCylinder(scmp);

       if(dist<mindist)
         mindist=dist;
  }
  return ScalarReal(mindist);
}
*/

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
         areaOut;
} siminfo_t;



template<typename T>
void intern_simDefect( typename STGM::ClusterList<T>::Type &cl,
                          STGM::Converter< STGM::ConverterFunction<T> > &converter, siminfo_t &info ) {

  typename STGM::ClusterList<T >::iterator_t  jt, endit;

  Rboolean stopit = FALSE;
  STGM::CDefect<T> *head, *last;

  double minDist=0,
         MPI4 = info.distTol*sqrt(M_PI_4); // to compare to minDist (distTol is weight factor < 1)
  int i=1, k=0, nclust=0, N=converter.N;

  head = converter(0);
  head->project();
  cl.push_back(head);

  while(i < N && !stopit)
  {
      head = converter(i);
      cl.push_back(head);
      // project single object, store points
      head->project();
      // check if projected defect of single
      // object is already large enough
      if(head->m_inner && head->m_area > info.areaIn) {
          stopit = TRUE;
          break;
      } else if(!head->m_inner && head->m_area > info.areaOut) {
          stopit = TRUE;
          break;
      }

      //T &scmp = head->m_object;
      jt = cl.begin(); endit = cl.end(); --endit;

      for(k = 0; jt != endit; ++k )
      {
         last = *jt;
         minDist = last->descent(head);
          //Rprintf("d: %f, head: %f, last: %f \n",minDist,sqrt(head->m_area),sqrt(last->m_area));
          // check minimum distance
          if(minDist < MPI4*MIN(sqrt(head->m_area),sqrt(last->m_area))) {
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
  typename STGM::ClusterList<T >::iterator_t  jt,kt;
  int nProtected=0;
  // construct R elements
  SEXP R_names;
  PROTECT(R_names = allocVector(STRSXP, 8));
  ++nProtected;

  SET_STRING_ELT(R_names, 0, mkChar("id"));
  SET_STRING_ELT(R_names, 1, mkChar("n"));
  SET_STRING_ELT(R_names, 2, mkChar("B"));
  SET_STRING_ELT(R_names, 3, mkChar("interior"));
  SET_STRING_ELT(R_names, 4, mkChar("A"));
  SET_STRING_ELT(R_names, 5, mkChar("inner"));
  SET_STRING_ELT(R_names, 6, mkChar("T"));
  SET_STRING_ELT(R_names, 7, mkChar("label"));

  STGM::CDefect<T> *p, *q;
  int m = 0, j = 0, l = 0, maxSize = 0;

  SEXP R_cl = R_NilValue;
  SEXP R_tmp, R_id, R_type, R_interior, R_label, R_num, R_area, R_time;

  if(PL>100) {
    for(kt = cl.begin(); kt != cl.end(); ) {
        if( (*kt)->m_size > 1) {
            ++kt;
            continue;
        }
        // free non cluster object
        Free(*kt);
        kt = cl.erase(kt);
    }

    PROTECT(R_cl = allocVector(VECSXP,cl.size()));   ++nProtected;
    for(j = 0, jt = cl.begin(); jt != kt; ++jt, ++j )
    {
        m = (*jt)->m_size;
        if(m > maxSize)
          maxSize = m;
        PROTECT(R_id = allocVector(INTSXP,m));
        PROTECT(R_num = allocVector(INTSXP,m));
        PROTECT(R_type = allocVector(INTSXP,m));
        PROTECT(R_interior = allocVector(INTSXP,m));
        PROTECT(R_label = allocVector(STRSXP,m));
        PROTECT(R_area = allocVector(REALSXP,m));
        PROTECT(R_time = allocVector(REALSXP,m));

        l = 0; p = *jt;
        while(p != 0) {
            INTEGER(R_id)[l] = p->m_id;
            INTEGER(R_num)[l] = p->m_num;
            INTEGER(R_type)[l] = p->m_type;
            INTEGER(R_interior)[l] = p->m_interior;
            REAL(R_area)[l] = p->m_area;
            REAL(R_time)[l] = p->m_time;
            SET_STRING_ELT(R_label,l,mkChar(p->m_label));
            q = p;
            p = p->next;
            // free accumulated clusters
            Free(q);
            ++l;
        }

        PROTECT(R_tmp = allocVector(VECSXP,8));
        SET_VECTOR_ELT(R_tmp,0,R_id);
        SET_VECTOR_ELT(R_tmp,1,R_num);
        SET_VECTOR_ELT(R_tmp,2,R_type);
        SET_VECTOR_ELT(R_tmp,3,R_interior);
        SET_VECTOR_ELT(R_tmp,4,R_area);
        SET_VECTOR_ELT(R_tmp,5,ScalarInteger((*jt)->m_inner));
        SET_VECTOR_ELT(R_tmp,6,R_time);
        SET_VECTOR_ELT(R_tmp,7,R_label);
        setAttrib(R_tmp, R_NamesSymbol, R_names);

        SET_VECTOR_ELT(R_cl,j,R_tmp);
        UNPROTECT(8);
    }
  } else {
      for(kt = cl.begin(); kt != cl.end(); ++kt ) {
        if( (*kt)->m_size > 1) {
            jt = kt;  // need to know only the last defect
            continue;
        }
        Free(*kt);
      }
      m = (*jt)->m_size;
      maxSize = m;
      PROTECT(R_id = allocVector(INTSXP,m));
      PROTECT(R_num = allocVector(INTSXP,m));
      PROTECT(R_type = allocVector(INTSXP,m));
      PROTECT(R_interior = allocVector(INTSXP,m));
      PROTECT(R_label = allocVector(STRSXP,m));
      PROTECT(R_area = allocVector(REALSXP,m));
      PROTECT(R_time = allocVector(REALSXP,m));

      l = 0; p = *jt; // p is last defect
      while(p != 0) {
          INTEGER(R_id)[l] = p->m_id;
          INTEGER(R_num)[l] = p->m_num;
          INTEGER(R_type)[l] = p->m_type;
          INTEGER(R_interior)[l] = p->m_interior;
          REAL(R_area)[l] = p->m_area;
          REAL(R_time)[l] = p->m_time;
          SET_STRING_ELT(R_label,l,mkChar(p->m_label));
          q = p;
          p = p->next;
          // free accumulated clusters
          Free(q);
          ++l;
      }


      PROTECT(R_tmp = allocVector(VECSXP,8));
      SET_VECTOR_ELT(R_tmp,0,R_id);
      SET_VECTOR_ELT(R_tmp,1,R_num);
      SET_VECTOR_ELT(R_tmp,2,R_type);
      SET_VECTOR_ELT(R_tmp,3,R_interior);
      SET_VECTOR_ELT(R_tmp,4,R_area);
      SET_VECTOR_ELT(R_tmp,5,ScalarInteger((*jt)->m_inner));
      SET_VECTOR_ELT(R_tmp,6,R_time);
      SET_VECTOR_ELT(R_tmp,7,R_label);
      setAttrib(R_tmp, R_NamesSymbol, R_names);

      PROTECT(R_cl = allocVector(VECSXP,1)); ++nProtected;
      SET_VECTOR_ELT(R_cl,0,R_tmp);
      UNPROTECT(8);
  }

  setAttrib(R_cl, install("aIn"), ScalarReal(info.areaIn));
  setAttrib(R_cl, install("aOut"), ScalarReal(info.areaOut));
  setAttrib(R_cl, install("areaMax"), ScalarReal(info.areaMax));
  setAttrib(R_cl, install("maxSize"), ScalarInteger(maxSize));

  UNPROTECT(nProtected);
  return R_cl;
}


SEXP SimDefect(SEXP R_vname, SEXP R_clust, SEXP R_dist, SEXP R_areaIn, SEXP R_areaOut,
               SEXP R_print_level, SEXP R_env)
{
    SEXP Rs = R_NilValue;
    if(isNull(R_env) || !isEnvironment(R_env))
      error("Should provide environment for function evaluation.");

    PROTECT(Rs = getVar(AS_CHARACTER(R_vname),R_env));
    PL = asInteger(AS_INTEGER(R_print_level));

    if (TYPEOF(Rs) == PROMSXP)
      Rs = eval(Rs, R_env);
    else error("Expression does not evaluate to a promise.");

    siminfo_t info = {asReal(AS_NUMERIC(R_dist)),0,asReal(AS_NUMERIC(R_areaIn)),
                      asReal(AS_NUMERIC(R_areaOut))};

    const char * name = GET_OBJECT_CLASS(Rs);
    if( !std::strcmp(name, "prolate" ) || !std::strcmp(name, "oblate" ) || !std::strcmp(name, "spheroid" )) {
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
    }
    UNPROTECT(1);
    return R_NilValue;
}


/**\brief Calculate spheroid projection only the area
 *
 * @param R_spheroids
 * @return numeric vector, the areas
 */
SEXP GetSpheroidOnlyProjectionArea(SEXP R_spheroids) {
  STGM::Spheroids spheroids = convert_C_Spheroids(R_spheroids);

  int n=spheroids.size();

  SEXP R_ret;
  PROTECT(R_ret = allocVector(REALSXP,n));
  //STGM::CEllipse2 ellipse;
  for(int i=0;i<n; i++) {
    STGM::CEllipse2 ellipse=spheroids[i].delamProjection();
    REAL(R_ret)[i]=M_PI*ellipse.a()*ellipse.b();
  }
  UNPROTECT(1);
  return R_ret;
}


SEXP GetSpheroidBothProjection(SEXP R_spheroids) {
  STGM::Spheroids spheroids = convert_C_Spheroids(R_spheroids);

  SEXP R_ret;
  PROTECT(R_ret = allocVector(VECSXP,2));

  int n=spheroids.size();
  STGM::Ellipses ellipses1, ellipses2;
  for(int i=0;i<n; i++) {
    ellipses1.push_back(spheroids[i].delamProjection());
    ellipses2.push_back(spheroids[i].crackProjection());
  }

  SET_VECTOR_ELT(R_ret,0,convert_C2R_ellipses(ellipses1));
  SET_VECTOR_ELT(R_ret,1,convert_C2R_ellipses(ellipses2));

  UNPROTECT(1);
  return R_ret;
}

/* R Interface functions  */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
static R_CMethodDef CEntries[]  = {
      CALLDEF(sdm,6),
      CALLDEF(ContactRadius,8),
      {NULL, NULL, 0}
};
static R_CallMethodDef CallEntries[] = {
        CALLDEF(GetPointsForConvexHull,2),
        CALLDEF(GetSpheroidProjection,2),
        CALLDEF(GetCylinderProjection,3),
        CALLDEF(GetSphereProjection,2),
        CALLDEF(Cluster,3),
        CALLDEF(SimDefect,7),
        CALLDEF(GetSpheroidBothProjection,1),
        CALLDEF(GetSpheroidOnlyProjectionArea,1),
      {NULL, NULL, 0}
};


void R_init_simLife(DllInfo *info) {
  R_registerRoutines(info, CEntries, CallEntries, NULL, NULL);
}

void R_unload_simLife(DllInfo *info){
  /* Release resources. */
}
