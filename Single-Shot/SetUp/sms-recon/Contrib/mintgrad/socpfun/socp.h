/*
 *   socp.h
 *
 *   Second-Order Cone Programming
 *   Header file for socp.c and socp_mex.c
 *
 *   mlobo@isl.stanford.edu -- 96/97
 */

#include "mex.h"
#include "lapack.h"
#include "blas.h"
/* basic macros */

#define SQR(x)    ((x)*(x))
#define MAX(x,y)  ((x) > (y) ? (x) : (y))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))


/* constants for socp algorithm */

#define MAX_ITER_PLANE 20
#define MAX_LAMBDA2 1e-2  /* plane search stopping crit. */
#define DIV_ALPHA 2
#define MIN_ALPHA 1e-6    /* max. of about 20 line search iterations */


/*
 #ifdef nounderscores
 #define ddot_ ddot
 #define dcopy_ dcopy
 #define daxpy_ daxpy
 #define dscal_ dscal
 #define dgemv_ dgemv
 #define dsyr_ dsyr
 #define dsyrk_ dsyrk
 #define dposvx_ dposvx
 #define dgelss_ dgelss
 #endif
 */
/* #ifdef APPLE */
/*
 #define ddot_ cblas_ddot
 #define dcopy_ cblas_dcopy
 #define daxpy_ cblas_daxpy
 #define dscal_ cblas_dscal
 #define dgemv_ cblas_dgemv
 #define dsyr_ cblas_dsyr
 #define dsyrk_ cblas_dsyrk
 #define dposvx_ cblas_dposvx
 #define dgelss_ cblas_dgelss
 */
/* #endif */

/* BLAS 1 */
/*
   double ddot_( );
   void dcopy_( );
   void daxpy_( );
   void dscal_( );
 */
/* BLAS 2
   void dgemv_( );
   void dsyr_( );
 */
/* BLAS 3
    void dsyrk_( );
 */
/* LAPACK */
/*void dposvx_( );*/
/*void dgelss_( );*/



/* socp.c */

void socp_getwork(
								/* input args.:  problem dimensions and max. num. of iterations */
								ptrdiff_t L,
								ptrdiff_t *N,
								ptrdiff_t n,
								int max_iter,
								int out_mode,
								/* output args.:  dimensions of history output matrix, */
								/*   number of doubles, pointers and ints required for workspace */
								int *mhist,
								int *nhist,
								int *ndbl,
								int *nint
								);

int socp(
								/* problem dimensions */
								ptrdiff_t L,
								ptrdiff_t  *N,
								ptrdiff_t n,
								/* problem data */
								double *f,
								double *A,
								double *b,
								/* in:  initial primal and dual strictly feasible points */
								/* out: final points */
								double *x,
								double *z,
								/* stopping criteria, */
								/*   *iter on entry is max. number of iterations, */
								/*   on exit actual number performed */
								double abs_tol,
								double rel_tol,
								double target,
								mwSize *iter,
								/* algorithm parameter */
								double Nu,
								/* reason for exit, output matrix with extra info */
								/* out_mode specifies what will be stored in *hist */
								int *info,
								int out_mode,
								double *hist,
								/* workspace, use socp_getwork() to determine required sizes */
								double *dblwork,
								ptrdiff_t *intwork
								);
