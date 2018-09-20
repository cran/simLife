/**
 * @file convert.cpp
 * @date 04/26/2016
 *
 * @brief methods for conversion of R objects
 *
 * @author M.Baaske
 */

#include <R.h>
#include <Rdefines.h>

#include "GeometricPrimitives.h"

#define COPY_C2R_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      REAL((R))[_i+DIM*_j] = (M)[_i][_j];         \
} while(0)


#define COPY_R2C_MATRIX(M,R,DIM)                  \
do {                                              \
    int _i, _j;                                   \
    for (_i = 0; _i < DIM; _i++)                  \
      for (_j = 0; _j < DIM; _j++)                \
      (M)[_i][_j] = REAL((R))[_i+DIM*_j];         \
} while(0)

#define SET_REAL_VECTOR(Rv, vec)  			\
do {										\
  int _i;		                            \
  for(_i = 0; _i < LENGTH( (Rv) ); _i++) {	\
    REAL((Rv))[_i] = (vec)[_i];				\
  }											\
} while(0)

SEXP getVar(SEXP name, SEXP rho)
{
    if(!isString(name) || length(name) != 1)
        error("name is not a single string");
    if(!isEnvironment(rho))
        error("rho should be an environment");

    return findVar(install(CHAR(STRING_ELT(name, 0))), rho);
}


///* get the elements of a list */
SEXP getListElement (SEXP list, const char *str)
{
     SEXP elmt = R_NilValue;
     SEXP names = getAttrib(list, R_NamesSymbol);

     for (int i = 0; i < length(list); i++)
         if(std::strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
             elmt = VECTOR_ELT(list, i);
             break;
         }
     return elmt;

}

SEXP convert_R_Ellipse2(STGM::CEllipse2 &ellipse) {
  SEXP R_tmp = R_NilValue;
  SEXP R_center, R_minor, R_major, R_ab, R_A;
  const char *nms[] = {"id", "type", "center", "A", "ab", "minor", "major", "phi", "S", ""};

  PROTECT(R_tmp = mkNamed(VECSXP,nms));
  PROTECT(R_center = allocVector(REALSXP, 2));
  PROTECT(R_ab = allocVector(REALSXP, 2));
  PROTECT(R_A = allocMatrix(REALSXP,2,2));
  PROTECT(R_minor = allocVector(REALSXP, 2));
  PROTECT(R_major = allocVector(REALSXP, 2));

  STGM::CVector2d &center = ellipse.center();
  SET_REAL_VECTOR(R_center,center);

  STGM::CVector2d &minor = ellipse.minorAxis();
  SET_REAL_VECTOR(R_minor,minor);

  STGM::CVector2d &major = ellipse.majorAxis();
  SET_REAL_VECTOR(R_major,major);

  REAL(R_ab)[0] = ellipse.a();    // major semi-axis (for both prolate/oblate)
  REAL(R_ab)[1] = ellipse.b();	  // minor semi-axis (for both prolate/oblate)

  const STGM::CMatrix2d &A = ellipse.MatrixA();
  COPY_C2R_MATRIX(A,R_A,2);

  SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
  SET_VECTOR_ELT(R_tmp,1,ScalarInteger(STGM::ELLIPSE_2D));
  SET_VECTOR_ELT(R_tmp,2,R_center);
  SET_VECTOR_ELT(R_tmp,3,R_A);
  SET_VECTOR_ELT(R_tmp,4,R_ab);
  SET_VECTOR_ELT(R_tmp,5,R_minor);
  SET_VECTOR_ELT(R_tmp,6,R_major);

  /* return angle between [0,2pi] */
  SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.phi()));
  SET_VECTOR_ELT(R_tmp,8,ScalarReal(ellipse.b()/ellipse.a()));

  UNPROTECT(6);
  return(R_tmp);

}


STGM::CSphere convert_C_Sphere(SEXP R_sphere)
{
	STGM::CVector3d ctr(REAL(VECTOR_ELT( R_sphere, 1)));

	return STGM::CSphere(ctr,REAL(VECTOR_ELT(R_sphere, 2))[0],
			INTEGER(VECTOR_ELT(R_sphere, 0))[0],
			translateChar(asChar(getAttrib(R_sphere, install("label")))),
			INTEGER(AS_INTEGER(getAttrib(R_sphere, install("interior"))))[0]);
}

STGM::Spheres convert_C_Spheres(SEXP R_spheres)
{
  STGM::Spheres spheres;
  size_t num = (size_t) LENGTH(R_spheres);
  spheres.reserve(num);

  for(size_t i=0; i < num; i++)
      spheres.push_back(convert_C_Sphere( VECTOR_ELT(R_spheres,i) ) );

  return spheres;
}


STGM::CCylinder convert_C_Cylinder(SEXP R_cyl)
{
  int interior = 1;
  const char *label = "N";

  /* Ferrit - actually a sphere as a cylinder because of application of FBA
   * otherwise 'P' for particle and 'N' default (no label)
   * */
  SEXP R_label = R_NilValue;
  PROTECT(R_label = getAttrib(R_cyl, install("label")));
  if(!isNull(R_label))
    label = translateChar(asChar(R_label));
  else{
	error("Undefined attribute `label`.");
  }

  SEXP R_int = R_NilValue;
  PROTECT(R_int = getAttrib(R_cyl, install("interior")));
  if(!isNull(R_int))
    interior = INTEGER(AS_INTEGER(R_int))[0];
  else {
    error("Undefined attribute `interior`.");
  }

  /** just copy from R */
  STGM::CVector3d ctr(REAL(VECTOR_ELT(R_cyl, 1))),
	   		        u(REAL(VECTOR_ELT(R_cyl, 5)));

  UNPROTECT(2);
  return STGM::CCylinder(ctr,u,
    		  	REAL(VECTOR_ELT(R_cyl,4))[0],
                REAL(VECTOR_ELT(R_cyl,6))[0],
				REAL(VECTOR_ELT(R_cyl,7))[0],
				REAL(VECTOR_ELT(R_cyl,7))[1],
                INTEGER(VECTOR_ELT(R_cyl, 0))[0],
				label, interior);
}


STGM::Cylinders convert_C_Cylinders(SEXP R_cyls)
{
  STGM::Cylinders cylinders;
  size_t num = (size_t) LENGTH(R_cyls);
  cylinders.reserve(num);

  for(size_t i=0; i<num; i++)
   cylinders.push_back( convert_C_Cylinder(VECTOR_ELT(R_cyls,i)));

  return cylinders;
}


STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid)
{
  STGM::CVector3d ctr(REAL(VECTOR_ELT( R_spheroid, 1)));
  STGM::CVector3d   u(REAL(VECTOR_ELT( R_spheroid, 2)));

  double *acb = REAL(VECTOR_ELT(R_spheroid, 3));
  double *angles = REAL(VECTOR_ELT(R_spheroid, 4));

  return STGM::CSpheroid (ctr,acb[0],acb[1],acb[2],u,angles[0],angles[1],
			  INTEGER(VECTOR_ELT(R_spheroid, 0))[0],
			  translateChar(asChar(getAttrib(R_spheroid, install("label")))),
			  INTEGER(AS_INTEGER(getAttrib(R_spheroid, install("interior"))))[0]);
}


STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  STGM::Spheroids spheroids;
  size_t num = (size_t) LENGTH(R_spheroids);
  spheroids.reserve(num);

  for(size_t i=0; i < num; i++)  {
      spheroids.push_back( convert_C_Spheroid( VECTOR_ELT(R_spheroids,i) ) );
  }

  return spheroids;
}

STGM::CEllipse2 convert_C_Ellipse2(SEXP R_E)
{
   STGM::CVector2d ctr(REAL(VECTOR_ELT(R_E,2)));
   STGM::CMatrix2d A(REAL(VECTOR_ELT(R_E,3)));

   return STGM::CEllipse2(A,ctr,INTEGER(VECTOR_ELT(R_E,0))[0]);

}

STGM::Ellipses2 convert_C_Ellipses2(SEXP R_E) {
	STGM::Ellipses2 ellipses2;
	size_t num = (size_t) LENGTH(R_E);

	ellipses2.reserve(num);
	for(size_t i=0; i < num; i++)  {
		  ellipses2.push_back( convert_C_Ellipse2( VECTOR_ELT(R_E,i) ) );
	}

	return ellipses2;
}

