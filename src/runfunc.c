/*===========================================================================*/
/* runfunc - running window functions                                        */
/* Copyright (C) 2005 Jarek Tuszynski                                        */
/* Distributed under GNU General Public License version 3                    */
/*===========================================================================*/

/*==================================================*/
/* Index:                                           */
/*  |------------------+------+------+----------|   */
/*  | function         | NaN  | Edge | Underflow|   */
/*  |------------------+------+------+----------|   */
/*  | sum_exact        | NA   | NA   | 1024     |   */
/*  | cumsum_exact     | NA   | NA   | 1024     |   */
/*  | runmean_exact    | yes  | yes  | 1024     |   */
/*  | runmean          | yes  | yes  |    2     |   */
/*  | runmean_lite     | no   | no   |    1     |   */
/*  | runmin           | yes  | yes  |   NA     |   */
/*  | runmax           | yes  | yes  |   NA     |   */
/*  | runquantile_lite | no   | no   |   NA     |   */
/*  | runquantile      | yes  | yes  |   NA     |   */
/*  | runmad_lite      | no   | no   |   NA     |   */
/*  | runmad           | yes  | yes  |   NA     |   */
/*  | runsd_lite       | no   | no   |    1     |   */
/*  | runsd            | yes  | yes  |    2     |   */
/*  |------------------+------+------+----------|   */
/*  NaN - means support for NaN and possibly Inf    */
/*  edge - means calculations are done all the way  */
/*         to the edges                             */
/*  underflow - means at maximum how many numbers   */
/*    are used to store results of addition in case */
/*    of underflow                                  */
/*==================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <float.h>


/* #define DEBBUG */
#ifdef DEBBUG
  int R_finite(double x) { return ( (x)==(x) ); }
  #define Calloc(b, t)  (t*) calloc(b,sizeof(t))
  #define Free free
  #define PRINT(x) { if ((x)==(x)) printf("%04.1f ",x); else printf("NaN "); }
#else
  #include <R.h>
  #include <Rinternals.h>
#endif

#define notNaN(x)   ((x)==(x))
#define isNaN(x)  (!((x)==(x)))
#define MIN(y,x) ((x)<(y) && (x)==(x) ? (x) : (y))
#define MAX(y,x) ((x)>(y) && (x)==(x) ? (x) : (y))
#define SQR(x) ((x)*(x))

/*============================================================================*/
/* The following macros were inspired by msum from                            */
/* http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393090             */
/* Quote from it:                                                             */
/* "Full precision summation using multiple doubles for intermediate values   */
/* Rounded x+y stored in hi with the round-off stored in lo.  Together        */
/* hi+lo are exactly equal to x+y.  The loop applies hi/lo summation          */
/* to each partial so that the list of partial sums remains exact.            */
/* Depends on IEEE-754 arithmetic guarantees.  See proof of correctness at:   */
/* www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps"  */
/*============================================================================*/

/* SumErr - macro calculating error of the summing operation */
#define SumErr(a,b,ab) ((((a)>(b)) == ((a)>-(b))) ?  (b) - ((ab)-(a)) : (a) - ((ab)-(b)) )
/* SUM_1 - macro for calculating Sum+=x; Num+=n; Which is NaN aware and have minimal (single number) overflow error correction */
#define SUM_1(x,n, Sum, Err, Num)   if (R_finite(x)){ y=Sum; Err+=x; Sum+=Err; Num+=n; Err=SumErr(y,Err,Sum);  } 
#define mpartial 1024	


void SUM_N(double x, int n, double *partial, int *npartial, int *Num) 
{
  if (R_finite(x)){ 
    int j, i;
    double hi, lo, y;
    for (i=j=0; j<*npartial; j++) {
      y  = partial[j];
      hi = y + x;
      lo = SumErr(x,y,hi); 
      if (lo && i<mpartial) partial[i++] = lo;
      x = hi; 
    }
    partial[i] = x; 
    *npartial   = i+1;
    *Num+=n; 
  }
}

/*========================================================================================*/
/* Each iteration of an insertion sort removes an element from the input data, inserting  */
/*  it at the correct position in the already sorted list, until no elements are left     */
/* in the input. For unsorted data is much less efficient than the more advanced          */
/* algorithms such as Quicksort, Heapsort, or Merge sort, but it is very efficient on     */
/* data sets which are already substantially sorted (almost O(n))                         */
/* Referances:                                                                            */
/*   R. Sedgewick: Algorithms. Addison-Wesley (1988) (page 99)                            */
/*   http://en.wikipedia.org/wiki/Insertion_sort                                          */
/*   http://www.cs.ubc.ca/spider/harrison/Java/sorting-demo.html                          */
/* Input:                                                                                 */
/*   V    - data array we operate on remains unchanged (we assume that no number in idx   */
/*          array is longer than V                                                        */
/*   idx  - index numbers of elements in V to be partially sorted                         */
/*   nIdx - length of idx array                                                           */
/* Output:                                                                                */
/*   idx  - index numbers of sorted elements in V w                                       */
/*========================================================================================*/

void insertion_sort(const double *V, int *idx, const int nIdx) 
{ 
  int i, j, id;
  double v;
  for (i=1; i<nIdx; i++) {
    id = idx[i];
    v  = V[id];
    for(j=i; j>0; j--) {
      if (V[idx[j-1]]<v) break;
      idx[j] = idx[j-1];
    }
    idx[j] = id;
  }
}

/*==================================================================*/
/* Array Sum without round-off errors.                              */
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty double                                            */
/*   nIn  - size of In array                                        */
/* Output :                                                         */
/*   Out  - Array sum                                               */
/*==================================================================*/
void sum_exact(double *In, double *Out, const int *nIn)
{
  int i, j, n=*nIn, npartial=0, Num=0;
  double *in, x, partial[mpartial];
  in=In;
  for(i=0; i<n; i++, in++) {
    x = *in;
    SUM_N(x, 1, partial, &npartial, &Num);
  }
  *Out = partial[0];
  for(j=1; j<npartial; j++) *Out += partial[j];
}

/*==================================================================*/
/* Array cumulative sum without round-off errors.                   */
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results              */
/*   nIn  - size of In and Out arrays                               */
/* Output :                                                         */
/*   Out  - results of cumulative sum operation                     */
/*==================================================================*/
void cumsum_exact(double *In, double *Out, const int *nIn)
{
  int i, j, n=*nIn, npartial=0, Num=0;
  double *in, *out, x, partial[mpartial];
  in=In; out=Out;
  for(i=0; i<n; i++, in++, out++) {
    x = *in;
    SUM_N(x, 1, partial, &npartial, &Num);
    *out = partial[0];
    for(j=1; j<npartial; j++) *out += partial[j];
  }
}

//void runsum_exact(double *In, double *Out, int *Size, const int *nIn, const int *nWin)
//{ /* medium size version with NaN's and round-off correction, but edge calculation*/
//  int i, j, k, n=*nIn, m=*nWin, npartial=0, count=0, *size;
//  double *in, *out, partial[mpartial], x;
//  k = m>>1;                              /* half of moving window */                           
//  for(i=0; i<=k; i++) {
//    Out [i]=Out [n-i-1]=0;
//    Size[i]=Size[n-i-1]=0;
//  }
//  if (m>=n) return;
//
//  /* step 1: sum of the first window *out = sum(x[0:(m-1)]) + err1 */
//  in=In; out=Out+k; size=Size+k;
//  for(i=0; i<m; i++, in++) {
//    x = *in;
//    if (R_finite(x)) {
//      add2partials(x, partial, npartial);
//      count++;
//    }
//  }
//  *size = count;
//  *out  = partial[0];
//  for(j=1; j<npartial; j++) *out += partial[j];
//
//  /* step 2: runsum of the rest of the vector. Inside loop is same as:   */
//  /* *out = *(out-1) + *in - *(in-m); but with round of error correction */
//  out++; size++;
//  for(i=m; i<n; i++, out++, in++, size++) { 
//    x = *in;    /* add the new value */
//    if (R_finite(x)) {
//      add2partials(x, partial, npartial);
//      count++;
//    }
//    x = -*(in-m); /* drop the value that goes out of the window */
//    if (R_finite(x)) {
//	    add2partials(x, partial, npartial);
//	    count--;
//    }
//    *size = count;
//    *out  = partial[0];
//    for(j=1; j<npartial; j++) *out += partial[j];
//  }
//}



/*==================================================================================*/
/* Mean function applied to (running) window. The fastest implementation with no    */
/* edge calculations, no NaN support, and no overflow correction                    */  
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runmean_lite(double *In, double *Out, const int *nIn, const int *nWin)
{
  int i, k2, n=*nIn, m=*nWin;
  double *in, *out, Sum, d;
  k2 = m>>1;           /* right half of window size */
  d  = 1.0/m;
  in=In; out=Out; 
  Sum = 0;             /* we need to calculate initial 'Sum' */
  /* step 1 - find mean of elements 0:(k2-1) */      
  for(i=0; i<k2; i++) Sum += in[i];
  /* step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m; i++, out++) {
    Sum += in[i];
    *out = Sum/(i+1);
  }
  /* step 3: runsum of the rest of the vector. Inside loop is same as:   */
  /* *out = *(out-1) - *in + *(in+m);  */
  for(i=m; i<n; i++, out++, in++) {
    Sum += in[m] - *in;
    *out = Sum*d;
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++, out++, in++) {
    Sum -= *in;
    *out = Sum/(m-1-i);
  }
}

/*==================================================================================*/
/* Mean function applied to (running) window. All additions performed using         */
/* addition algorithm which tracks and corrects addition round-off errors (see      */  
/*  http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps)*/
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runmean(double *In, double *Out, const int *nIn, const int *nWin)
{ /* medium size version with NaN's and edge calculation, but only one level of round-off correction*/
  int i, k2, Num, n=*nIn, m=*nWin;
  double *in, y, *out, Err, Sum;
  double NaN = (0.0/0.0);
  k2  = m>>1;         /* right half of window size */
  in=In; out=Out; 
  Sum = 0;           /* we need to calculate initial 'Sum' */
  Err = 0;
  Num = 0;
  /* step 1 - find mean of elements 0:(k2-1) */      
  for(i=0; i<k2; i++) {
    SUM_1(in[i], 1, Sum, Err, Num)
  }
  /* step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m; i++, out++) {
    SUM_1(in[i], 1, Sum, Err, Num)
    *out = (Num ? (Sum+Err)/Num : NaN);  /* save mean and move window */
  }
  /* step 3: runsum of the rest of the vector. Inside loop is same as:   */
  /* *out = *(out-1) - *in + *(in+m); but with round of error correction */
  for(i=m; i<n; i++, out++, in++) {
    SUM_1(in[m] ,  1, Sum, Err, Num)
    SUM_1(-(*in), -1, Sum, Err, Num)
    *out = (Num ? (Sum+Err)/Num : NaN);  /* save mean and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++, out++, in++) {
    SUM_1(-(*in), -1, Sum, Err, Num)
    *out = (Num ? (Sum+Err)/Num : NaN);  /* save mean and move window */
  }
}


/*==================================================================================*/
/* Mean function applied to (running) window. All additions performed using         */
/* addition algorithm which tracks and corrects addition round-off errors (see      */  
/*  http://www-2.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps)*/
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runmean_exact(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation, full round-off correction*/
  int i, j, k2, n=*nIn, m=*nWin, npartial=0, Num=0;
  double *in, *out, partial[mpartial], Sum;
  double NaN = (0.0/0.0);

  k2 = m>>1;         /* right half of window size */
  in=In; out=Out; 
  /* step 1 - find mean of elements 0:(k2-1) */      
  for(i=0; i<k2; i++) {
    SUM_N(in[i], 1, partial, &npartial, &Num);
  }
  /* step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m; i++, out++) {
    SUM_N(in[i], 1, partial, &npartial, &Num);
    for(Sum=j=0; j<npartial; j++) Sum += partial[j];
    *out = (Num ? Sum/Num : NaN);  /* save mean and move window */
  }
  /* step 3: runsum of inner section. Inside loop is same as:   */
  /* *out = *(out-1) - *in + *(in+m); but with round of error correction */
  for(i=m; i<n; i++, out++, in++) {
    SUM_N(in[m] , 1, partial, &npartial, &Num);
    SUM_N(-(*in),-1, partial, &npartial, &Num);
    for(Sum=j=0; j<npartial; j++) Sum += partial[j];
    *out = (Num ? Sum/Num : NaN);  /* save mean and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++, out++, in++) {
    SUM_N(-(*in),-1, partial, &npartial, &Num);
    for(Sum=j=0; j<npartial; j++) Sum += partial[j];
    *out = (Num ? Sum/Num : NaN);  /* save mean and move window */
  }
}


/*==================================================================*/
/* minimum function applied to moving (running) window              */ 
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                               */
/*   nWin - size of the moving window (odd)                         */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/
void runmin(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation */
  int i, j, k2, n=*nIn, m=*nWin;
  double ptOut, Min, *in, *out, CST = DBL_MAX;
  double NaN = (0.0/0.0);

  k2  = m>>1;               /* right half of window size */  
  in  = In;
  out = Out;
  /* --- step 1 - find min of elements 0:(k2-1) */      
  Min=CST;                  /* we need to calculate  initial 'Min' */
  for(i=0; i<k2; i++) Min = MIN(Min,in[i]);  /* find minimum over a window of length k2 */
  /* --- step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m-1; i++) {
    Min=MIN(Min,in[i]);     /* cumulative min */
    *(out++) = (Min==CST ? NaN : Min); /* save 'Min' and move window */
  }
  /* --- step 3 - the inner section - window of constant size is moving  */      
  ptOut=CST;
  for(i=m-1; i<n; i++) {
    if(ptOut==Min) {        /* if point comining out of the window was window's min than ... */
      Min=CST;              /* we need to recalculate 'Min' */
      for(j=0; j<m; j++) 
        Min=MIN(Min,in[j]); /* find minimum over a window of length m */
    } else                  /* if point comining out of the window was NOT window min than min of ... */
      Min=MIN(Min,in[m-1]); /* ... window's first m-1 points is still 'Min', so we have to add a single point */
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Min==CST ? NaN : Min); /* save 'Min' and move window */
  }
  /* --- step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++) {
    if(ptOut==Min) {        /* if point comining out of the window was window's extreme than ... */
      Min=CST;              /* we need to recalculate 'Min' */
      for(j=0; j<m-i-1; j++) 
        Min=MIN(Min,in[j]); /* find minimum over a window of length m */
    } 
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Min==CST ? NaN : Min);  /* and fill the space with window extreme and move window */
  }
}

/*==================================================================*/
/* Maximum function applied to moving (running) window              */ 
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                               */
/*   nWin - size of the moving window (odd)                         */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/
void runmax(double *In, double *Out, const int *nIn, const int *nWin)
{ /* full-blown version with NaN's and edge calculation */
  int i, j, k2, n=*nIn, m=*nWin;
  double ptOut, Max, *in, *out, CST = -DBL_MAX;
  double NaN = (0.0/0.0);

  k2  = m>>1;               /* right half of window size */
  in  = In;
  out = Out;
  /* step 1 - find max of elements 0:(k2-1) */      
  Max= CST;                /* we need to calculate  initial 'Max' */
  for(i=0; i<k2; i++) Max = MAX(Max,in[i]);  /* find maximum over a window of length k2 */
  /* step 2 - left edge - start expanding the moving window to the right */      
  for(i=k2; i<m-1; i++) {
    Max=MAX(Max,in[i]);     /* cumulative max */
    *(out++) = (Max==CST ? NaN : Max); /* save 'Max' and move window */
  }
  /* step 3 - the inner section - window of constant size is moving  */      
  ptOut=CST;
  for(i=m-1; i<n; i++) {
    if(ptOut==Max) {        /* if point comaxing out of the window was window's max than ... */
      Max=CST;              /* we need to recalculate 'Max' */
      for(j=0; j<m; j++) 
        Max=MAX(Max,in[j]); /* find maximum over a window of length m */
    } else                  /* if point comining out of the window was NOT window max than max of ... */
      Max=MAX(Max,in[m-1]); /* ... window's first m-1 points is still 'Max', so we have to add a single point */
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Max==CST ? NaN : Max); /* save 'Max' and move window */
  }
  /* step 4 - right edge - right side reached the end and left is shrinking  */      
  for(i=0; i<k2; i++) {
    if(ptOut==Max) {        /* if point comining out of the window was window's extreme than ... */
      Max=CST;              /* we need to recalculate 'Max' */
      for(j=0; j<m-i-1; j++) 
        Max=MAX(Max,in[j]); /* find maximum over a window of length m */
    } 
    ptOut = *(in++);        /* store point comming out of the window for future use and move window */
    *(out++) = (Max==CST ? NaN : Max); /* and fill the space with window extreme and move window */
  }
}

/*==========================================================================*/
/* Calculate element in the sorted array of nWin elements that corresponds  */
/*          to a quantile of 'type' and 'prob'                              */ 
/* Input :                                                                  */
/*   prob - Quantile probability from 0 to 1                                */
/*   nWin - how many elements in dataset the quantile will be calculated on */
/*   type - integer between 1 and 9 indicating type of quantile             */
/*          See http://mathworld.wolfram.com/Quantile.html                  */
/* Output :                                                                 */
/*   return  - which element in the sorted array of nWin elements correspond*/
/*          to the prob.  If non-integer than split into intger (i) and     */
/*          real(r) parts, then quantile = v[i]*(1-r) + v[i+1]*r            */
/*==========================================================================*/
double QuantilePosition(double prob, int nWin, int type)
{ /* the following code is based on code from R's quantile.default function */
  double a, b, h, nppm, fuzz;
  int j;
  if (type <= 3) {    // Types 1, 2 and 3 are discontinuous sample qs. 
    if (type == 3) nppm = nWin * prob - .5; // n * probs + m; m = -0.5 
    else           nppm = nWin * prob;      // m = 0 
    j = (int) floor(nppm);
    switch(type) {
      case  1: h = (nppm > j ? 1 : 0);   break;                 // type 1 
      case  2: h = (nppm > j ? 1 : 0.5); break;                 // type 2 
      case  3: h = ((nppm==j) && ((j>>1) == 0) ? 0 : 1); break; // type 3 
      default: h=1;                      break;
    }
  } else {            // Types 4 through 9 are continuous sample qs. 
    switch(type) { 
     case  4: a=0; b=1;    break;  
     case  5: a=b=0.5;     break;
     case  6: a=b=0;       break;
     case  7: a=b=1;       break;
     case  8: a=b=1.0/3.0; break; 
     case  9: a=b=3.0/8.0; break;
     default: a=b=1;       break;
    }
    nppm = a + prob * (nWin + 1 - a - b); // n*probs + m 
    fuzz = 4 * DBL_EPSILON;
    j = (int) floor(nppm + fuzz);
    h = nppm - j;
    h = (fabs(h) < fuzz ? 0 : h);
  } 
  nppm = j+h;
  nppm = (nppm<1    ? 1    : nppm);
  nppm = (nppm>nWin ? nWin : nppm);
  return nppm - 1; // C arrays are zero based   
}


/*==================================================================*/
/* quantile function applied to (running) window                    */ 
/* Input :                                                          */
/*   In   - array to run moving window over will remain umchanged   */
/*   Out  - empty space for array to store the results. Out is      */
/*          assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                               */
/*   nWin - size of the moving window                               */
/*   Prob - Array of probabilities from 0 to 1                      */
/*   nProb - How many elements in Probs array?                      */
/*   type - integer between 1 and 9 indicating type of quantile     */
/*          See http://mathworld.wolfram.com/Quantile.html          */
/* Output :                                                         */
/*   Out  - results of runing moving window over array In and       */
/*          colecting window mean                                   */
/*==================================================================*/
void runquantile_lite(double *In, double *Out, const int *nIn, const int *nWin, const double *Prob, const int *nProb, const int *Type)
{ /* internal region only is calculated. Edges, NaN's are not handled */
  int i, j, k, *idx, d, n=*nIn, m=*nWin, nPrb=*nProb;
  double *Win, *in, *out, r, ip, *prob, pointOut, ext;
  k   = m>>1;                          /* half of window size */
  in  = In;
  out = Out+k;

  if (nPrb==1 && (*Prob==1 || *Prob==0)) { /* trivial case shortcut - if prob is 0 or 1 than wind windows min or max */
    d = (*Prob==0 ? -1 : 1);          /* runmin d=-1; runmax d=1*/
    pointOut=ext=0;
    for(i=m-1; i<n; i++) {
      if(pointOut==ext) {             /* if point comining out of the window was window's extreme than ... */
        ext=in[0];                    /* in the 2 lines below we do things diffferent when extreme is a max vs. min */
        if (d==1) { 
          for(j=1; j<m; j++) if (ext<in[j]) ext=in[j];  /* find maximum over a window of length m */
        } else {   
          for(j=1; j<m; j++) if (ext>in[j]) ext=in[j];  /* find minimum over a window of length m */
        } 
      } else                          /* if point comining out of the window was NOT window extreme than we know ... */
        if (ext*d<in[m-1]*d) ext=in[m-1]; /* ... extreme of window's first m-1 points, so we have to add a single point */
      pointOut = *(in++);             /* store point comming out of the window for future use and move window */
      *(out++) = ext;                 /* and fill the space with window extreme and move window */
    }
  } else {                            /* non-trivial case */
    idx = Calloc(m,int   );           /* index will hold partially sorted index numbers of Save array */
    Win = Calloc(m,double);           /* stores all points of the current running window */
    prob = Calloc(nPrb,double);       /* stores all points of the current running window */
    for(i=0; i<m; i++) {
      Win[i] = *(in++);               /* initialize running window */
      idx[i] = i;                     /* and its index */
    }
    in--;                             /* last point of the window will be placed again */
    for(d=0; d<nPrb; d++)             /* for each probability */
      prob[d] = QuantilePosition(Prob[d], m, *Type); /* store common size for speed */
    for(j=i=m-1; i<n; i++) {
      Win[j] = *(in++);               /* Move Win to the right: replace a[i-m] with a[m] point  */
      insertion_sort(Win,idx,m);      /* sort current Win */
      for(d=0; d<nPrb; d++) {
        r = modf( prob[d], &ip );     /* Divide p into its fractional and integer parts */
        k = (int) ip-1;               /* k-1 instead of k because in C arrays are 0 based and in R they are 1 based */
        if (r) r = Win[idx[k]]*(1-r) + Win[idx[k+1]]*r; /* interpolate */
        else   r = Win[idx[k]];  
        *(out+d*n) = r;
      }
      out++;
      j = (j+1)%m;                    /* index goes from 0 to m-1, and back to 0 again  */
    }
    Free(Win);
    Free(idx);
    Free(prob);
  }
}

void runquantile(double *In, double *Out, const int *nIn, const int *nWin, const double *Prob, const int *nProb, const int *Type)
{ /* full-blown version with NaN's and edge calculation */
  int i, j, k1, k2, *idx, d, n=*nIn, m=*nWin, nPrb=*nProb, mm, k, type=*Type, count=0;
  double *Win, *in, *out, r, ip, Max, p, *prob, BIG=DBL_MAX;
  double NaN = (0.0/0.0);

  k2  = m>>1;                      /* right half of window size */
  k1  = m-k2-1;                    /* left half of window size */
  in  = In;
  out = Out;

  if (nPrb==1 && *Prob==0) {       /* trivial case shortcut - if prob is 0 or 1 than find windows min */
    runmin(In, Out, nIn, nWin);
  } else if (nPrb==1 && *Prob==1) {/* trivial case shortcut - if prob is 0 or 1 than find windows max */
    runmax(In, Out, nIn, nWin);
  } else {                         /* non-trivial case */
    idx  = Calloc(m,int   );       /* index will hold partially sorted index numbers of Save array */
    Win  = Calloc(m,double);       /* stores all points of the current running window */
    prob = Calloc(nPrb,double);    /* stores all points of the current running window */
    for(i=0; i<m; i++) idx[i] = i; /* and its index */
    for(i=0; i<k2; i++) {
      Win[i] = *(in++);            /* initialize running window */
      if (isNaN(Win[i])) Win[i]=BIG; else count++;
    }
    /* --- step 1 : left edge -----------------------------------------------------------------*/
    for(j=k2, i=0; i<=k1; i++){
      mm = i+k2+1;                 /* window width */
      j  = mm-1;
      Win[j] = *(in++);            /* Move Win to the right: replace a[i-m] with a[m] point  */
      if (isNaN(Win[j])) Win[j]=BIG; else count++;
      insertion_sort(Win,idx,mm);  /* sort current Win */
      for(d=0; d<nPrb; d++) {      /* for each probability */
        if (count>0) {             /* not all points in the window are NaN*/
          p = QuantilePosition(Prob[d], count, type);
          r = modf( p, &ip );      /* Divide p into its fractional and integer parts */
          k = (int) ip;            /* k-1 instead of k because in C arrays are 0 based and in R they are 1 based */
          if (r) r = Win[idx[k]]*(1-r) + Win[idx[k+1]]*r; /* interpolate */
          else   r = Win[idx[k]];  
        } else r = NaN;            /* all points in the window are NaN*/
        out[d*n] = r;
      }
      out++;
      j = (j+1)%m;                 /* index goes from 0 to m-1, and back to 0 again  */
      //printf("1-------- %3.0f %3.0f %3.0f %3.0f %3.0f - %i\n", Win[idx[0]],Win[idx[1]],Win[idx[2]],Win[idx[3]],Win[idx[4]], count);
    }
    /* --- step 2: inner section ----------------------------------------------------------------*/
    for(d=0; d<nPrb; d++)          /* for each probability */
      prob[d] = QuantilePosition(Prob[d], m, type); /* store common size for speed */
    for(i=m; i<n; i++) {
      if (Win[j]<BIG) count--;     /* Point leaving the window was not a NAN */
      Win[j] = *(in++);            /* Move Win to the right: replace a[i-m] with a[m] point  */
      if (isNaN(Win[j])) Win[j]=BIG; else count++;
      insertion_sort(Win,idx,m);   /* sort current Win */
      for(d=0; d<nPrb; d++) {      /* for each probability */
        if (count>0) {             /* not all points in the window are NaN*/
          p = (count==m ? prob[d] : QuantilePosition(Prob[d], count, type));
          r = modf( p, &ip );      /* Divide p into its fractional and integer parts */
          k = (int) ip;            /* k-1 instead of k because in C arrays are 0 based and in R they are 1 based */
          if (r) r = Win[idx[k]]*(1-r) + Win[idx[k+1]]*r; /* interpolate */
          else   r = Win[idx[k]];  
        } else r = NaN;            /* all points in the window are NaN*/
        out[d*n] = r;
      }
      out++;
      j = (j+1)%m;                 /* index goes from 0 to m-1, and back to 0 again  */
      //printf("2-------- %3.0f %3.0f %3.0f %3.0f %3.0f - %i\n", Win[idx[0]],Win[idx[1]],Win[idx[2]],Win[idx[3]],Win[idx[4]], count);
    }
    /* --- step 3 : right edge ----------------------------------------------------------*/
    Max = Win[idx[m-1]];           /* store window maximum */
    for(i=0; i<k2; i++) {
      if (Win[j]<BIG) count--;     /* Point leaving the window was not a NAN */
      Win[j] = Max;                /* setting to maximum will push it to the end*/
      mm = m-i-1;                  /* window width */
      insertion_sort(Win,idx,mm);  /* sort current Win from 1-mm */
      for(d=0; d<nPrb; d++) {      /* for each probability */
        if (count>0) {             /* not all points in the window are NaN*/
          p = QuantilePosition(Prob[d], count, type);
          r = modf( p, &ip );      /* Divide p into its fractional and integer parts */
          k = (int) ip;            /* k-1 instead of k because in C arrays are 0 based and in R they are 1 based */
          if (r) r = Win[idx[k]]*(1-r) + Win[idx[k+1]]*r; /* interpolate */
          else   r = Win[idx[k]];
        } else r = NaN;            /* all points in the window are NaN*/
        out[d*n] = r;
      }
      out++;
      j = (j+1)%m;                 /* index goes from 0 to m-1, and back to 0 again  */
      //printf("3-------- %3.0f %3.0f %3.0f %3.0f %3.0f - %i\n", Win[idx[0]],Win[idx[1]],Win[idx[2]],Win[idx[3]],Win[idx[4]], count);
    }
    Free(Win);
    Free(idx);
    Free(prob);
  }
}


/*==================================================================================*/
/* MAD function applied to moving (running) window                                  */ 
/* No edge calculations and no NAN support                                          */ 
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Ctr  - array storing results of runmed or other running average function       */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/

void runmad_lite(double *In, double *Ctr, double *Out, const int *nIn, const int *nWin)
{ 
  int i, k2, k1, j, l, *idx, n=*nIn, m=*nWin;
  double *Win1, *Win2, *in, *out, *ctr, med0, med;

  idx  = Calloc(m,int   );       /* index will hold partially sorted index numbers of Save array */
  Win1 = Calloc(m,double);       /* stores all "In" points of the current running window */
  Win2 = Calloc(m,double);       /* stores all "abs(In-Crt) values of the current running window */
  k2   = m>>1;                   /* right half of window size */
  k1   = m-k2-1;                 /* left  half of window size */
  in   = In;                     /* initialize pointer to input In vector */
  out  = Out+k1;                 /* initialize pointer to output Mad vector */
  ctr  = Ctr+k1;                 /* initialize pointer to input  Ctr vector */
  med0 = 0;                      /* med0 - will save previous center (median) so we know it changed */
  for(i=0; i<m; i++) {
    Win1[i] = Win2[i] = *(in++); /* initialize running windows */
    idx[i] = i;                  /* and its indexes */
  }
  in--;                          /* last point of the window will be placed again */
  for(j=i=m-1; i<n; i++) {       /* the main loop */
    Win1[j] = *(in++);           /* Move Win to the right: replace a[i-m] with a[m] point  */    
    med     = *(ctr++);          /* get median of current Win1 */
    if (med==med0)               /* median did not changed */
      Win2[j]=fabs(Win1[j]-med); /* recalculate one point */
    else for(l=0; l<m; l++)      /* median did changed */
      Win2[l]=fabs(Win1[l]-med); /* recalculate whole window Win2 */
    insertion_sort(Win2,idx,m);  /* sort Win2 - if medians did not change than  data should be sorted*/
    *(out++) = (Win2[idx[k1]]+Win2[idx[k2]])*0.5;    /* find mad of current Win1 and store it */
    med0 = med;                  /* save previous median */
    j = (j+1)%m;                 /* index of Win positions goes from 0 to m-1, and back to 0 again  */
  }
  Free(Win2);
  Free(Win1);
  Free(idx);
}

/*==================================================================================*/
/* MAD function applied to moving (running) window                                  */ 
/* with edge calculations and NAN support                                           */ 
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Ctr  - array storing results of runmed or other running average function       */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runmad(double *In, double *Ctr, double *Out, const int *nIn, const int *nWin)
{ 
  int i, k1, k2, kk1, kk2, j, l, mWin, *idx, n=*nIn, m=*nWin, Num=0;
  double *Win1, *Win2, *in, *out, *ctr, med0, med, BIG=DBL_MAX-1;

  idx  = Calloc(m,int   );        /* index will hold partially sorted index numbers of Save array */
  Win1 = Calloc(m,double);        /* stores all points of the current running window: Values*/
  Win2 = Calloc(m,double);        /* stores all points of the current running window: Values - median*/
  k2   = m>>1;                    /* right half of window size */
  k1   = m-k2-1;                  /* left half of window size */
  in   = In;                      /* initialize pointer to input In vector */
  out  = Out;                     /* initialize pointer to output Mad vector */
  ctr  = Ctr;                     /* initialize pointer to output Mad vector */

  /* --- step 1 : left edge -----------------------------------------------------------------*/
  for(i=0; i<m; i++) {
    Win1[i] = *(in++);            /* initialize running windows */
    idx [i] = i;                  /* and its index */
  }
  med0 = BIG;                     /* med0 - will save previous center (median) so we know it changed */
  for(j=k2, i=0; i<=k1; i++, j++){/* start with a small window and make it bigger */
    mWin = j+1;
    med  = *(ctr++);              /* find median of current Win1 */
    if (med==med0) {              /* median did not changed */
      Win2[j]=fabs(Win1[j]-med);  /* recalculate one point */
      if (R_finite(Win2[j])) Num++; else Win2[j]=BIG;
    } else {                      /* recalculate whole window Win2 */
      for(Num=l=0; l<mWin; l++) { /* median did changed */
        Win2[l]=fabs(Win1[l]-med);
        if (R_finite(Win2[l])) Num++; else Win2[l]=BIG;
      }
    }
    insertion_sort(Win2,idx,Num); /* sort Win2 - if medians did not change than  data should be sorted*/
    kk2  = Num>>1;                /* right half of window size */
    kk1  = Num-kk2-1;             /* left half of window size. if nWin is odd than kk1==kk2 */
    *(out++) = (Win2[idx[kk1]]+Win2[idx[kk2]])*0.5;    /* find mad of current Win1 and store it */
//  med0 = med;                   /* save previous median */
//  printf("1-------- "); for(l=0; l<m; l++) PRINT(Win1[idx[l]]); printf(" - %f\n",med); 
  }
  /* --- step 2: inner section ----------------------------------------------------------------*/
  for(j=0, i=m; i<n; i++) {
    Win1[j] = *(in++);            /* Move Win to the right: replace a[i-m] with a[m] point  */    
    med     = *(ctr++);           /* find median of current Win1 */
    if (med==med0) {              /* median did not changed */
      if (Win2[j]<BIG) Num--;
      Win2[j]=fabs(Win1[j]-med);  /* recalculate one point */
      if (R_finite(Win2[j])) Num++; else Win2[j]=BIG;
    } else {                      /* recalculate whole window Win2 */
      for(Num=l=0; l<m; l++) {    /* median did changed */
        Win2[l]=fabs(Win1[l]-med);
        if (R_finite(Win2[l])) Num++; else Win2[l]=BIG;
      }
    }
    insertion_sort(Win2,idx,Num); /* sort Win2 - if medians did not change than  data should be sorted*/
    kk2  = Num>>1;                /* right half of window size */
    kk1  = Num-kk2-1;             /* left half of window size. if nWin is odd than kk1==kk2 */
    *(out++) = (Win2[idx[kk1]]+Win2[idx[kk2]])*0.5;    /* find mad of current Win1 and store it */
    med0 = med;                   /* save previous median */
    j = (j+1)%m;                  /* index goes from 0 to m-1, and back to 0 again  */
//  printf("2-------- "); for(l=0; l<m; l++) PRINT(Win1[idx[l]]); printf(" - %f\n",med); 
  }
  /* --- step 3 : right edge ----------------------------------------------------------*/
  for(i=0; i<m; i++) {
    Win1[i] = In[n-i-1];          /* initialize running windows */
    idx [i] = i;                  /* and its index */
  }
  med0=BIG;                       /* make sure med0!=med so the ELSE is used below */
  for(j=k1, i=1; i<=k2; i++, j++){ /* start with a small window and make it bigger */
    mWin = j+1;
    med = Ctr[n-i];               /* find median of current Win1 */
    if (med==med0) {              /* median did not changed */
      Win2[j]=fabs(Win1[j]-med);  /* recalculate one point */
      if (R_finite(Win2[j])) Num++; else Win2[j]=BIG;
    } else {                      /* recalculate whole window Win2 */
      for(Num=l=0; l<mWin; l++) { /* median did  changed */
        Win2[l]=fabs(Win1[l]-med);
        if (R_finite(Win2[l])) Num++; else Win2[l]=BIG;
      }
    }
    insertion_sort(Win2,idx,Num); /* sort Win2 - if medians did not change than data should be sorted*/
    kk2  = Num>>1;                /* right half of window size */
    kk1  = Num-kk2-1;             /* left half of window size. if nWin is odd than kk1==kk2 */
    Out[n-i] = (Win2[idx[kk1]]+Win2[idx[kk2]])*0.5;    /* find mad of current Win1 and store it */
//  med0 = med;                   /* save previous median */
//  printf("3-------- "); for(l=0; l<m; l++) PRINT(Win1[idx[l]]); printf(" - %f\n",med); 
  }
  Free(Win2);
  Free(Win1);
  Free(idx);
} 

/*==================================================================================*/
/* Standard Deviation function applied to moving (running) window                   */ 
/* No edge calculations and no NAN support                                          */ 
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Ctr  - array storing results of runmed or other running average function       */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runsd_lite(double *In, double *Ctr, double *Out, const int *nIn, const int *nWin)
{ 
  int i, k2, k1, j, l, n=*nIn, m=*nWin;
  double *Win1, *Win2, *in, *out, *ctr, med0, med, Sum=0;

  Win1 = Calloc(m,double);       /* stores all points of the current running window: Values */
  Win2 = Calloc(m,double);       /* stores all points of the current running window: Values - avr */
  k2   = m>>1;                   /* right half of window size */
  k1   = m-k2-1;                 /* left half of window size */
  in   = In;                     /* initialize pointer to input In vector */
  out  = Out+k1;                 /* initialize pointer to output Mad vector */
  ctr  = Ctr+k1;                 /* initialize pointer to output Mad vector */
  med0 = *ctr+1;                 /* med0 - will save previous center (median) sowe know it changed */
  for(i=0; i<m; i++) {
    Win1[i] = Win2[i] = *(in++); /* initialize running windows */
  }
  in--;                          /* last point of the window will be placed again */
  for(j=i=m-1; i<n; i++) {       /* the main loop */
    Win1[j] = *(in++);           /* Move Win to the right: replace a[i-m] with a[m] point  */    
    med     = *(ctr++);          /* find median of current Win1 */
    if (med==med0) {             /* median did not changed */
      Sum -= Win2[j];
      Win2[j]=SQR(Win1[j]-med);  /* recalculate one point */
      Sum += Win2[j];
    } else {
      for(Sum=l=0; l<m; l++) {   /* median did not changed */
        Win2[l]=SQR(Win1[l]-med);/* recalculate whole window Win2 */
        Sum += Win2[l];
      }
    }
    *(out++) = sqrt(Sum/(m-1));  /* save mean and move window */
    med0 = med;                  /* save previous median */
    j = (j+1)%m;                 /* index goes from 0 to m-1, and back to 0 again  */
  }
  Free(Win2);
  Free(Win1);
}

/*==================================================================================*/
/* Standard Deviation function applied to moving (running) window                   */ 
/* With edge calculations and NAN support                                           */ 
/* Input :                                                                          */
/*   In   - array to run moving window over will remain umchanged                   */
/*   Ctr  - array storing results of runmed or other running average function       */
/*   Out  - empty space for array to store the results                              */
/*   nIn  - size of arrays In and Out                                               */
/*   nWin - size of the moving window                                               */
/* Output :                                                                         */
/*   Out  - results of runing moving window over array In and colecting window mean */
/*==================================================================================*/
void runsd(double *In, double *Ctr, double *Out, const int *nIn, const int *nWin)
{ 
  int i, k1, k2, j, l, mWin, n=*nIn, m=*nWin, Num;
  double *Win1, *Win2, *in, *out, *ctr, med0, med, Sum, Err, y, BIG=DBL_MAX;
  double NaN = (0.0/0.0);

  Sum=Err=Num=0;
  Win1 = Calloc(m,double);        /* stores all points of the current running window: Values */
  Win2 = Calloc(m,double);        /* stores all points of the current running window: Values - avr */
  k2   = m>>1;                    /* right half of window size */
  k1   = m-k2-1;                  /* left half of window size */
  in   = In;                      /* initialize pointer to input In vector */
  out  = Out;                     /* initialize pointer to output Mad vector */
  ctr  = Ctr;                     /* initialize pointer to output Mad vector */

  /* --- step 1 : left edge -----------------------------------------------------------------*/
  for(i=0; i<m; i++) Win1[i] = *(in++); /* initialize running windows */
  med0 = BIG;                     /* med0 - will save previous center (median) so we know it changed */
  for(j=k2, i=0; i<=k1; i++, j++){/* start with a small window and make it bigger */
    mWin = j+1;
    med  = *(ctr++);              /* find median of current Win1 */
    if (med==med0) {              /* median did not changed */
      Win2[j]=SQR(Win1[j]-med);   /* recalculate one point */
      SUM_1(Win2[j], 1, Sum, Err, Num);
    } else {
      Sum=Err=Num=0;
      for(l=0; l<mWin; l++) {     /* median did not changed */
        Win2[l]=SQR(Win1[l]-med); /* recalculate whole window Win2 */
        SUM_1(Win2[l], 1, Sum, Err, Num);
      }
    }
    *(out++) = (Num>1 ? sqrt((Sum+Err)/(Num-1)) : NaN);  /* save std and move window */
    med0 = med;                   /* save previous median */
    /*printf("1-------- "); for(l=0; l<m; l++) PRINT(Win1[idx[l]]); printf(" - %f\n",med); */
  }
  /* --- step 2: inner section ----------------------------------------------------------------*/
  for(j=0, i=m; i<n; i++) {
    Win1[j] = *(in++);            /* Move Win to the right: replace a[i-m] with a[m] point  */    
    med     = *(ctr++);           /* find median of current Win1 */
    if (med==med0) {              /* median did not changed */
      SUM_1(-Win2[j],-1, Sum, Err, Num);
      Win2[j]=SQR(Win1[j]-med);   /* recalculate one point */
      SUM_1( Win2[j], 1, Sum, Err, Num);
    } else {
      Sum=Err=Num=0;
      for(l=0; l<m; l++) {        /* median did not changed */
        Win2[l]=SQR(Win1[l]-med); /* recalculate whole window Win2 */
        SUM_1(Win2[l], 1, Sum, Err, Num);
      }
    }
    *(out++) = (Num>1 ? sqrt((Sum+Err)/(Num-1)) : NaN);  /* save std and move window */
    med0 = med;                   /* save previous median */
    j = (j+1)%m;                  /* index goes from 0 to m-1, and back to 0 again  */
    /*printf("2-------- "); for(l=0; l<m; l++) PRINT(Win1[idx[l]]); printf(" - %f\n",med); */
  }
  /* --- step 3 : right edge ----------------------------------------------------------*/
  for(i=0; i<m; i++) Win1[i] = In[n-i-1]; /* initialize running windows */
  med0=BIG;                       /* make sure med0!=med so the ELSE is used below */
  for(j=k1, i=1; i<=k2; i++, j++){ /* start with a small window and make it bigger */
    mWin = j+1;
    med = Ctr[n-i];               /* find median of current Win1 */
    if (med==med0) {              /* median did not changed */
      Win2[j]=SQR(Win1[j]-med);   /* recalculate one point */
      SUM_1(Win2[j], 1, Sum, Err, Num);
    } else {
      Sum=Err=Num=0;
      for(l=0; l<mWin; l++) {       /* median did not changed */
        Win2[l]=SQR(Win1[l]-med);/* recalculate whole window Win2 */
        SUM_1(Win2[l], 1, Sum, Err, Num);
      }
    }
    Out[n-i] = (Num>1 ? sqrt((Sum+Err)/(Num-1)) : NaN);  /* save std and move window */
    med0 = med;                   /* save previous median */
    /*printf("3-------- "); for(l=0; l<m; l++) PRINT(Win1[idx[l]]); printf(" - %f\n",med); */
  }
  Free(Win2);
  Free(Win1);
}

#undef MIN
#undef MAX
#undef SQR
#undef SUM_1
#undef SumErr
#undef mpartial

#ifdef DEBBUG

int main( void ) {
  double NaN = (0.0/0.0);
  char s[] = "%02i "; //"%02i ";

  //double x[] ={1,1.123456e10,1.987654e20,-1.987654e20,-1.123456e10, -1};
  //double y[6];
  //double p[1024];
  //int n = 6;

  //cumsum_exact(x,y,&n);
  //printf("%e %e %e %e %e %e\n", x[0], x[1], x[2], x[3], x[4], x[5]);
  //printf("%e %e %e %e %e %e\n", y[0], y[1], y[2], y[3], y[4], y[5]);
  
  int i, nn = 25, k = 13, np=3, type=7;
  double xx[25], yy[3*25], y1[25], y2[25], y3[25], y4[25], y5[25], y6[25];
  double p[] = {0,0.5,1};
  //for(i=0; i<nn; i++) xx[i]=i;
  //xx[7] = NaN;
  //runmean1(xx, yy, &nn, &k);
  //for(i=0; i<nn; i++) PRINT(xx[i]); printf("\n");
  //for(i=0; i<nn; i++) PRINT(yy[i]); printf("\n");
  //printf("\n");
  for(i=0; i<nn; i++) xx[i]=i;
  for(i=5; i<12; i++) xx[i]=i;
  //runmin(xx, y1, &nn, &k);
  //runmax(xx, y2, &nn, &k);
  runmean(xx, y3, &nn, &k);
  //runmean_lite(xx, y4, &nn, &k);
  //runmean_exact(xx, y5, &nn, &k);
  runquantile(xx, yy, &nn, &k, p, &np, &type);
  for(i=0; i<nn; i++) PRINT(xx[i]); printf("Original\n");
  //for(i=0; i<nn; i++) PRINT(y1[i]); printf("Min\n");
  //for(i=0; i<nn; i++) PRINT(yy[i]); printf("Q1\n\n");
  //for(i=0; i<nn; i++) PRINT(y2[i]); printf("Max\n");
  for(i=0; i<nn; i++) PRINT(yy[i+25]); printf("Q2\n\n");
  //for(i=0; i<nn; i++) PRINT(yy[i+50]); printf("Q5\n\n");
  //for(i=0; i<nn; i++) PRINT(y3[i]); printf("Mean\n");
  //for(i=0; i<nn; i++) PRINT(y4[i]); printf("Mean_lite\n");
  //for(i=0; i<nn; i++) PRINT(y5[i]); printf("Mean_exact\n");
  runmad_lite (xx, yy+25, y5, &nn, &k); 
  for(i=0; 2*i<k; i++) y5[i]=y5[nn-1-i]=0;
  runmad(xx, yy+25, y6, &nn, &k);
  runsd(xx, y3, y4, &nn, &k);
  for(i=0; i<nn; i++) PRINT(y5[i]); printf("Mad lite\n");
  for(i=0; i<nn; i++) PRINT(y6[i]); printf("MAD\n");
  for(i=0; i<nn; i++) PRINT(y4[i]); printf("sd\n");

  getchar();
  return 0;
}

#endif
