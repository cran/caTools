/*===========================================================================*/
/* Written by Jarek Tuszynski. Copyright 2001-2005 SAIC.                     */
/* Software developed in conjunction with the National Cancer Institute, and */
/* distributed under "caBIO Software License" included in "COPYING" file.    */
/*===========================================================================*/

#include <R.h>
#include <Rinternals.h>

#define Swap(a,b,t) {(t)=(a); (a)=(b); (b)=(t);}
void BitSwap(void *v, int size)
{
  int i;
  char *p = (char*) v, tmp;
  if (size == 1) return;
  for (i=0; i<size/2; i++) Swap(p[i], p[size - i - 1], tmp);
}
#undef Swap

SEXP real2raw(SEXP DVec, SEXP Size, SEXP Swap)
{
  char *rvec;
  double *dvec;
  int i, j, dLen, rLen, size, swap;
  SEXP RVec = R_NilValue;

  size = asInteger(Size);       /* size of float type: either 4 or 8 */
  swap = asInteger(Swap);       /* is swap going to be performed */
  dLen = LENGTH(DVec);          /* length of dVec array */
  rLen = dLen*size;             /* length of rVec array */
  PROTECT(RVec = allocVector(RAWSXP, rLen));
  rvec = (void  *) RAW (RVec);  /* get pointer to R's RVec */
  dvec = (double*) REAL(DVec);  /* get pointer to R's DVec */
  if(size==sizeof(float)) {     /* if 4-byte floats are to be used instead of R's standard 8-byte */
    float f;
    for (i=j=0; i<dLen; i++, j+=size) {
      f = (float) dvec[i];
      memcpy(rvec+j, &f, size); /* copy numbers one by one */
    }
  } else memcpy(rvec, dvec, rLen); /* in case of same size floats as R is using - copy them all at once */
  if(swap && size>1)            /* perform byte swap if needed */
    for(i=0; i<dLen; i++) BitSwap(rvec+size*i, size);
  UNPROTECT(1);
  return RVec;
}

SEXP raw2real(SEXP RVec, SEXP Size, SEXP Swap)
{
  char *rvec;
  double *dvec;
  int i, j, dLen, rLen, size, swap;
  SEXP DVec = R_NilValue;

  size = asInteger(Size);       /* size of float type: either 4 or 8 */
  swap = asInteger(Swap);       /* is swap going to be performed */
  rLen = LENGTH(RVec);          /* length of rVec array */
  dLen = LENGTH(RVec)/size;     /* length of dVec array */
  PROTECT(DVec = allocVector(REALSXP, dLen));
  rvec = (void  *) RAW (RVec);  /* get pointer to R's RVec */
  dvec = (double*) REAL(DVec);  /* get pointer to R's DVec */
  if(swap && size>1)            /* perform byte swap if needed */
    for(i=0; i<dLen; i++) BitSwap(rvec+size*i, size);
  if(size==sizeof(float)) {     /* if 4-byte floats are to be used instead of R's standard 8-byte */
    float f;
    for (i=j=0; i<dLen; i++, j+=size) {
      memcpy(&f, rvec+j, size); /* copy numbers one by one */
      dvec[i] = (double) f;
    }
  } else memcpy(dvec, rvec, dLen*size); /* in case of same size floats as R is using - copy them all at once */
  UNPROTECT(1);
  return DVec;
}

