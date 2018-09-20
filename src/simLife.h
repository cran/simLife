#ifndef SRC_SIMLIFE_H_
#define SRC_SIMLIFE_H_

#include <R.h>
#include <Rdefines.h>

#include <list>
#include "GeometricPrimitives.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("simLife", String)
#else
#define _(String) (String)
#endif

using namespace std;

#define M_PI_4		0.78539816339744830962	/* pi/4 */

extern SEXP getListElement (SEXP list, const char *str);

extern STGM::CSphere convert_C_Sphere(SEXP R_sphere);
extern STGM::CCylinder convert_C_Cylinder(SEXP R_cylinder);
extern STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid);


#ifdef __cplusplus
extern "C" {
#endif

    /* .Call interface */
    void R_init_simLife(DllInfo *info);

    SEXP GetPointsForConvexHull(SEXP R_ellipses, SEXP R_n);

    SEXP SimDefect(SEXP R_vname, SEXP R_clust, SEXP R_dist, SEXP R_areaIn, SEXP R_areaOut,
    		 SEXP R_Tmax,SEXP R_print_level, SEXP R_env);

    SEXP GetSpheroidProjection(SEXP R_spheroids, SEXP R_crack_type);
    SEXP GetSpheroidOnlyProjectionArea(SEXP R_spheroids);
    SEXP GetSpheroidBothProjection(SEXP R_spheroids);

    SEXP GetSphereProjection(SEXP R_s, SEXP R_np);
    SEXP GetCylinderProjection(SEXP R_cylinders, SEXP R_crack_type,SEXP Rnp);

    ///-> exported but not for public use
    SEXP convexHull(SEXP R_points);

#ifdef __cplusplus
}
#endif


namespace STGM {

/**
 * @brief Class to describe the defect accumulation
 */
template<class OBJECT_T >
class CDefect {
 public:

  CDefect(OBJECT_T &object, int type, double time, int np = 20) :
            m_object(object),
            m_type(type),
            m_size(1),
			m_broken(0),
            m_area(0),
            m_time(time), next(0), last(0)
   {
     m_id = m_object.Id();
     m_label = m_object.m_label;
     m_num = 0;
     m_np = np;
     m_interior = m_object.interior();
     m_inner = m_interior;
     m_object.setCrackType(type);
     last = this;
   }


   /**
    * @brief the destructor is recursively called on next linked
    *        list element downwards of current head node
    */
  ~CDefect() {
     if(next!=0)
       delete next;
   };

   inline double distance(OBJECT_T &s) {
     return m_object.distance(s);
   }

   /**
    * @brief set point vector of projections
    */
   void project() {
      m_points.reserve(m_np);
      m_area = m_object.projectedPointsWithArea(m_points,m_np);
   }

   /**
    * \brief Insert newly projected points into head node
    * @param d current head node pointer
    */
   inline void update(CDefect<OBJECT_T> *d) {
       PointVector2d &P = d->m_points;
       m_points.insert(m_points.end(), P.begin(), P.end());
       if(m_inner)
        m_inner = d->m_inner;
       m_size += d->m_size;
   }

   /**
    *
    * @param head       head node
    * @return           minimum (weighted) distance as single linkage distance
    */
   inline double descent(CDefect<OBJECT_T> *head) {
     CDefect<OBJECT_T> *p = head, *q = this, *plast = p;
     double dist=0, minDist = 0;
     // search for two objects of
     // accumulated defects having minimum distance
     while(p!=0) {
       q = this;
       minDist = q->distance(p->m_object);
       while(q->next!=0) {
         q = q->next;
         dist = q->distance(p->m_object);
         if(dist < minDist)
           minDist = dist;
       }
       plast = p;
       p = p->next;
     }
     // distance vector between centers of objects
     CVector3d d(q->m_object.center()), z(0,0,1);
     d -= plast->m_object.center();

     double alpha = fabs(asin(d.dot(z)/d.Length()));
     //Rprintf("\t minDist: %f  alpha: %f  cos: %f (%f) \n",minDist,alpha,cos(alpha),sin(alpha));

     /* fixed correct angle:
      *   minDist/(cos(alpha)+1e-6)
      *   	now to
      *   minDist/(sin(alpha)+1e-6)
      */
     return minDist/(sin(alpha)+1e-6);
   }

   inline void append(CDefect<OBJECT_T> *current) {
      last->next = current;
      last = current->last;
   }

   /** using R's memory management */
   void* operator new (size_t size);
   void operator delete (void *p);

 public:
   OBJECT_T m_object;

   int m_type, m_id, m_inner, m_interior, m_size, m_num, m_np, m_broken;
   double m_area, m_time;
   const char * m_label;

   CDefect<OBJECT_T> *next, *last;
   PointVector2d m_points;

};

template<class T>
void* CDefect<T>::operator new (size_t size) {
    return Calloc(size,STGM::CDefect<T>);
    // ... then calling constructor (placement new)
}

template<class T>
void CDefect<T>::operator delete (void *p) { Free(p); }

/// conversion templates for spheroid, cylinder, sphere
template<class T> struct ClusterList {
  typedef std::list<CDefect<T> *> Type;
  typedef typename Type::iterator iterator_t;
};


/// 3d objects to be converted
template<class T>
struct ConverterFunction { ConverterFunction(); };

template<>
struct ConverterFunction<CSpheroid> {
  typedef CSpheroid Type;
  Type operator() (SEXP Rs) {  return convert_C_Spheroid(Rs); }
};

template<>
struct ConverterFunction<CCylinder> {
  typedef CCylinder Type;
  Type operator() (SEXP Rs) {  return convert_C_Cylinder(Rs); }
};


template<>
struct ConverterFunction<CSphere> {
  typedef CSphere Type;
  Type operator() (SEXP Rs) {  return convert_C_Sphere(Rs); }
};

/// conversion functor
template<class Function>
struct Converter {
  typedef typename Function::Type Type;

  SEXP Rs, Rcl;
  int id, type, N;
  double time;
  Function fun;

  Converter(SEXP _Rs, SEXP _Rcl) :
    Rs(_Rs), Rcl(_Rcl), id(0), type(0), N(length(_Rs)), time(0)
  {
    fun = Function();
  }

  STGM::CDefect<Type> * operator() (int i) {
     id   = INTEGER( AS_INTEGER(getListElement(VECTOR_ELT(Rcl,i), "id")))[0];
     type = INTEGER( AS_INTEGER(getListElement(VECTOR_ELT(Rcl,i), "B") ))[0];
     time =    REAL( AS_NUMERIC(getListElement(VECTOR_ELT(Rcl,i), "T") ))[0];

     // convert object from R list
     Type sp = fun(VECTOR_ELT(Rs,id-1));
     /* using R's memory management as overloaded new/delete ...*/
     STGM::CDefect<Type> *defect = new STGM::CDefect<Type>(sp,type,time);

     /* ...alternatively ... */
     //STGM::CDefect<Type> *defect = Calloc(1,STGM::CDefect<Type>);
     //try {
     //   new(defect)STGM::CDefect<Type>(sp,type,time);
     //} catch(...) {
     //   error("Init 'CDefect' failed: Allocation error.");
     //}
     //STGM::CDefect<Type> *defect = new STGM::CDefect<Type>(sp,type,time);

     return defect;
  };
};


} // end namespace

#endif
