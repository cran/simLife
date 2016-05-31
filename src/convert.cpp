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
      REAL((R))[_j+DIM*_i] = (M)[_i][_j];         \
} while(0)


SEXP getSingleCall(SEXP R_fname, SEXP R_arg, SEXP R_rho) {
    SEXP RCallBack = R_NilValue;
    PROTECT(RCallBack = allocVector(LANGSXP,2));
    SETCAR( RCallBack, findFun(install(CHAR(STRING_ELT(R_fname, 0))),R_rho ));
    SETCAR(CDR(RCallBack),R_arg);

    UNPROTECT(1);
    return RCallBack;
}

SEXP getVar(SEXP name, SEXP rho)
{
    SEXP ans;
    if(!isString(name) || length(name) != 1)
        error("name is not a single string");
    if(!isEnvironment(rho))
        error("rho should be an environment");
    ans = findVar(install(CHAR(STRING_ELT(name, 0))), rho);
    return ans;
}


///* get the elements of a list */
SEXP getListElement (SEXP list, const char *str)
{
     SEXP elmt = R_NilValue;
     SEXP names = getAttrib(list, R_NamesSymbol);

     for (int i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
             elmt = VECTOR_ELT(list, i);
             break;
         }
     return elmt;

}


STGM::CSphere convert_C_Sphere(SEXP R_sphere) {
  SEXP R_ctr;
  int interior=1;
  const char *label = "N";

  PROTECT(R_ctr = AS_NUMERIC( getListElement( R_sphere, "center")));

  int id = asInteger(AS_INTEGER( getListElement( R_sphere, "id")));
  double r = asReal(getListElement(R_sphere, "r"));

  if(!isNull(getAttrib(R_sphere, install("label"))))
    label = translateChar(asChar(getAttrib(R_sphere, install("label"))));

  if(!isNull(getAttrib(R_sphere, install("interior"))))
    interior = asLogical(getAttrib(R_sphere, install("interior")));

  return STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],r,id,label,interior);
  UNPROTECT(1);
}


STGM::CCylinder convert_C_Cylinder(SEXP R_cyl)
{
  int interior = 1;
  double radius = 0;
  const char *label = "N";

  if(!isNull(getAttrib(R_cyl, install("label"))))
    label = translateChar(asChar(getAttrib(R_cyl, install("label"))));
  if(!isNull(getAttrib(R_cyl, install("interior"))))
    interior = asLogical(getAttrib(R_cyl, install("interior")));
  if(!isNull(getAttrib(R_cyl, install("radius"))))
    radius = asReal(getAttrib(R_cyl, install("radius")));

  if(!std::strcmp(label,"F")) {
      SEXP R_ctr;
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      STGM::CVector3d ctr(REAL(R_ctr)),u(0,0,1);

      UNPROTECT(1);
      return STGM::CCylinder(ctr,u,0,asReal(getListElement(R_cyl, "r")),0,0,
                 radius, asInteger(getListElement(R_cyl, "id")), label, interior);

  } else {
      SEXP R_ctr, R_u, R_angles;
      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_cyl, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_cyl, "u")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_cyl, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)),u(REAL(R_u));
      UNPROTECT(3);

      return STGM::CCylinder(ctr,u, asReal(getListElement(R_cyl, "length")),
                asReal(getListElement(R_cyl, "r")), REAL(R_angles)[0],REAL(R_angles)[1],
                  radius, asInteger(getListElement(R_cyl, "id")), label, interior);
  }

}


STGM::Cylinders convert_C_Cylinders(SEXP R_cyls)
{
  SEXP R_tmp, R_ctr, R_u, R_angles;

  int interior = 1;
  double radius = 0;
  const char *label = "N";

  STGM::Cylinders cylinders;
  for(int i=0; i<length(R_cyls); i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_cyls,i));
      if(!isNull(getAttrib(R_tmp, install("label"))))
        label = translateChar(asChar(getAttrib(R_tmp, install("label"))));
      if(!isNull(getAttrib(R_tmp, install("interior"))))
        interior = asLogical(getAttrib(R_tmp, install("interior")));
      if(!isNull(getAttrib(R_tmp, install("radius"))))
        radius = asReal(getAttrib(R_tmp, install("radius")));

      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_tmp, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_tmp, "u")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_tmp, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)),u(REAL(R_u));

      cylinders.push_back(STGM::CCylinder(ctr,u,
                            asReal(getListElement( R_tmp, "length")),
                            asReal(getListElement( R_tmp, "r")),
                            REAL(R_angles)[0],REAL(R_angles)[1],radius,
                            asInteger(getListElement( R_tmp, "id")),
                            label,interior));
      UNPROTECT(4);

  }
  return cylinders;
}

STGM::Spheroids convert_C_Spheroids(SEXP R_spheroids)
{
  SEXP R_tmp, R_ctr, R_u, R_ab, R_angles;

  int interior = 1;
  double radius = 0;
  const char *label = "N";

  STGM::Spheroids spheroids;
  for(int i=0; i<length(R_spheroids); i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheroids,i));

      if(!isNull(getAttrib(R_tmp, install("label"))))
        label = translateChar(asChar(getAttrib(R_tmp, install("label"))));

      if(!isNull(getAttrib(R_tmp, install("interior"))))
        interior = asLogical(getAttrib(R_tmp, install("interior")));
      if(!isNull(getAttrib(R_tmp, install("radius"))))
        radius = asReal(getAttrib(R_tmp, install("radius")));

      PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_tmp, "center")));
      PROTECT( R_u      = AS_NUMERIC( getListElement( R_tmp, "u")));
      PROTECT( R_ab     = AS_NUMERIC( getListElement( R_tmp, "ab")));
      PROTECT( R_angles = AS_NUMERIC( getListElement( R_tmp, "angles")));

      STGM::CVector3d ctr(REAL(R_ctr)), u(REAL(R_u));

      spheroids.push_back(STGM::CSpheroid(ctr,REAL(R_ab)[0],REAL(R_ab)[1],u,1.0,
                              REAL(R_angles)[0],REAL(R_angles)[1],radius,
                              asInteger (getListElement( R_tmp, "id")),label,interior));
      UNPROTECT(5);

  }
  return spheroids;
}

STGM::Ellipses convert_C_Ellipses(SEXP R_ellipses)
{
  int id=0;
  STGM::Ellipses ellipses;
  SEXP R_tmp, R_ctr,R_ab, R_minor, R_major, R_A;

  double rot = 0;
  for(int i=0; i<length(R_ellipses); i++) {
     PROTECT(R_tmp   = VECTOR_ELT(R_ellipses,i));
     PROTECT(R_ctr   = AS_NUMERIC( getListElement( R_tmp, "center")));
     PROTECT(R_A     = AS_NUMERIC( getListElement( R_tmp, "A")));
     PROTECT(R_minor = AS_NUMERIC( getListElement( R_tmp, "minor")));
     PROTECT(R_major = AS_NUMERIC( getListElement( R_tmp, "major")));
     PROTECT(R_ab    = AS_NUMERIC( getListElement( R_tmp, "ab")));

     id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
     rot= asReal(AS_NUMERIC(getListElement( R_tmp, "rot")));

     /**
      * BUG: Constructor with matrix A has a bug to determine the correct angle phi
      */

     STGM::CPoint2d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1]);
     STGM::CPoint2d minorAxis(REAL(R_minor)[0],REAL(R_minor)[1]);
     STGM::CPoint2d majorAxis(REAL(R_major)[0],REAL(R_major)[1]);

     ellipses.push_back(STGM::CEllipse2(ctr,majorAxis,minorAxis,REAL(R_ab)[0],REAL(R_ab)[1],id,rot));
     UNPROTECT(6);
  }
  return ellipses;
}

SEXP convert_C2R_ellipses(STGM::Ellipses &ellipses) {
  int nProtected=0, dim=2, ncomps=8;
  int n = ellipses.size();

  SEXP names, R_resultlist;
  PROTECT(names = allocVector(STRSXP, ncomps));   ++nProtected;
  PROTECT(R_resultlist = allocVector(VECSXP,n));  ++nProtected;

  SET_STRING_ELT(names, 0, mkChar("id"));
  SET_STRING_ELT(names, 1, mkChar("center"));
  SET_STRING_ELT(names, 2, mkChar("ab"));
  SET_STRING_ELT(names, 3, mkChar("minor"));
  SET_STRING_ELT(names, 4, mkChar("major"));
  SET_STRING_ELT(names, 5, mkChar("A"));
  SET_STRING_ELT(names, 6, mkChar("phi"));
  SET_STRING_ELT(names, 7, mkChar("rot"));

  SEXP R_tmp,R_minor,R_major,R_A,R_center,R_ab;

  for(int i = 0; i < n; i++) {
      STGM::CEllipse2 & ellipse = ellipses[i];
      PROTECT(R_tmp = allocVector(VECSXP,ncomps));
      PROTECT(R_center = allocVector(REALSXP, dim));
      PROTECT(R_ab = allocVector(REALSXP, dim));
      PROTECT(R_minor = allocVector(REALSXP, dim));
      PROTECT(R_major = allocVector(REALSXP, dim));
      PROTECT(R_A = allocMatrix(REALSXP, dim,dim));

      REAL(R_center)[0] = ellipse.center()[0];
      REAL(R_center)[1] = ellipse.center()[1];

      REAL(R_minor)[0] = ellipse.minorAxis()[0];
      REAL(R_minor)[1] = ellipse.minorAxis()[1];

      REAL(R_major)[0] = ellipse.majorAxis()[0];
      REAL(R_major)[1] = ellipse.majorAxis()[1];

      REAL(R_ab)[0] = ellipse.a();
      REAL(R_ab)[1] = ellipse.b();

      COPY_C2R_MATRIX(ellipse.MatrixA(),R_A,dim);

      setAttrib(R_tmp, R_NamesSymbol, names);
      SET_VECTOR_ELT(R_tmp,0,ScalarInteger(ellipse.Id()));
      SET_VECTOR_ELT(R_tmp,1,R_center);
      SET_VECTOR_ELT(R_tmp,2,R_ab);
      SET_VECTOR_ELT(R_tmp,3,R_minor);
      SET_VECTOR_ELT(R_tmp,4,R_major);
      SET_VECTOR_ELT(R_tmp,5,R_A);
      SET_VECTOR_ELT(R_tmp,6,ScalarReal(ellipse.phi()));
      SET_VECTOR_ELT(R_tmp,7,ScalarReal(ellipse.rot()));

      SET_VECTOR_ELT(R_resultlist,i,R_tmp);
      UNPROTECT(6);
  }

  UNPROTECT(nProtected);
  return R_resultlist;
}

STGM::CSpheroid convert_C_Spheroid(SEXP R_spheroid)
{
  SEXP R_ctr, R_u, R_ab, R_angles;

  int id = asInteger (AS_INTEGER( getListElement( R_spheroid, "id")));
  PROTECT( R_ctr    = AS_NUMERIC( getListElement( R_spheroid, "center")));
  PROTECT( R_u      = AS_NUMERIC( getListElement( R_spheroid, "u")));
  PROTECT( R_ab     = AS_NUMERIC( getListElement( R_spheroid, "ab")));
  PROTECT( R_angles = AS_NUMERIC( getListElement( R_spheroid, "angles")));

  STGM::CVector3d ctr(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2]);
  STGM::CVector3d u(REAL(R_u)[0],REAL(R_u)[1],REAL(R_u)[2]);

  int interior=1;
  double radius=0;
  const char *label = "N";
  if(!isNull(getAttrib(R_spheroid, install("label"))))
    label = translateChar(asChar(getAttrib(R_spheroid, install("label"))));

  if(!isNull(getAttrib(R_spheroid, install("interior"))))
    interior = asLogical(getAttrib(R_spheroid, install("interior")));
  if(!isNull(getAttrib(R_spheroid, install("radius"))))
    radius = asReal(getAttrib(R_spheroid, install("radius")));

  UNPROTECT(4);
  return STGM::CSpheroid(ctr,REAL(R_ab)[0],REAL(R_ab)[1],u,1.0,REAL(R_angles)[0],REAL(R_angles)[1],radius,id,label,interior);
}


STGM::Spheres convert_C_Spheres(SEXP R_spheres) {
  SEXP R_tmp, R_ctr;
  int id=0,
      N=length(R_spheres);
  STGM::Spheres spheres;
  spheres.reserve(N);

  double r;
  for(int i=0; i<N; i++) {
      PROTECT(R_tmp = VECTOR_ELT(R_spheres,i));
      PROTECT(R_ctr = AS_NUMERIC( getListElement( R_tmp, "center")));
      id = asInteger (AS_INTEGER( getListElement( R_tmp, "id")));
      r = asReal(getListElement( R_tmp, "r"));

      spheres.push_back(STGM::CSphere(REAL(R_ctr)[0],REAL(R_ctr)[1],REAL(R_ctr)[2],r,id));
      UNPROTECT(2);
  }

  return spheres;
}

SEXP convert_R_SphereSystem(STGM::Spheres& spheres) {
  int ncomps=3;

  SEXP R_resultlist = R_NilValue;
  PROTECT(R_resultlist = allocVector(VECSXP, spheres.size()) );

  SEXP R_tmp, R_names, R_center;
  for(size_t k=0;k<spheres.size();k++)
  {
    STGM::CSphere &sphere = spheres[k];

    PROTECT(R_tmp = allocVector(VECSXP,ncomps));
    PROTECT(R_center = allocVector(REALSXP, 3));

    REAL(R_center)[0]=sphere.center()[0];
    REAL(R_center)[1]=sphere.center()[1];
    REAL(R_center)[2]=sphere.center()[2];

    SET_VECTOR_ELT(R_tmp,0,ScalarInteger(sphere.Id()));
    SET_VECTOR_ELT(R_tmp,1,R_center);
    SET_VECTOR_ELT(R_tmp,2,ScalarReal(sphere.r()));

    PROTECT(R_names = allocVector(STRSXP, ncomps));
    SET_STRING_ELT(R_names, 0, mkChar("id"));
    SET_STRING_ELT(R_names, 1, mkChar("center"));
    SET_STRING_ELT(R_names, 2, mkChar("r"));

    setAttrib(R_tmp, R_NamesSymbol, R_names);
    SET_VECTOR_ELT(R_resultlist,k,R_tmp);
    UNPROTECT(3);
  }

  UNPROTECT(1);
  return R_resultlist;
}
