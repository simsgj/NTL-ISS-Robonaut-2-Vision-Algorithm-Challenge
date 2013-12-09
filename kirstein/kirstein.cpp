#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <queue>
#include <string>
#include <algorithm>
#include<valarray>
#include<vector>
#include<ostream>
#include<cmath>
#include<limits>
#include<iostream>
#include<cassert>
#include<math.h>
#include<string>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<string.h>
#include<assert.h>
#include<vector>
#include<string.h>
#include<algorithm>
#include<vector>
#include<iostream>
#include<sstream>
#include<cassert>
#include<math.h>
#include<time.h>
using namespace std;
#define VL_USEFASTMATH 1

//*****lsqr begin


#define OFFSET(N, incX) ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)))

enum CBLAS_ORDER    {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

/*!
  \param[in]     N
  \param[in]     alpha
  \param[in]     X      
  \param[in]     incX
  \param[in,out] Y
  \param[in]     incY
*/
void
cblas_daxpy( const int N, const double alpha, const double *X,
             const int incX, double *Y, const int incY)
{
  int i;

  if (N     <= 0  ) return;
  if (alpha == 0.0) return;

  if (incX == 1 && incY == 1) {
      const int m = N % 4;

      for (i = 0; i < m; i++)
          Y[i] += alpha * X[i];
      
      for (i = m; i + 3 < N; i += 4) {
          Y[i    ] += alpha * X[i    ];
          Y[i + 1] += alpha * X[i + 1];
          Y[i + 2] += alpha * X[i + 2];
          Y[i + 3] += alpha * X[i + 3];
      }
  } else {
      int ix = OFFSET(N, incX);
      int iy = OFFSET(N, incY);

      for (i = 0; i < N; i++) {
          Y[iy] += alpha * X[ix];
          ix    += incX;
          iy    += incY;
      }
  }
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX
  \param[out]    Y
  \param[in]     incY
*/
void
cblas_dcopy( const int N, const double *X,
             const int incX, double *Y, const int incY)
{
  int i;
  int ix = OFFSET(N, incX);
  int iy = OFFSET(N, incY);

  for (i = 0; i < N; i++) {
      Y[iy]  = X[ix];
      ix    += incX;
      iy    += incY;
  }
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX
  \param[in]     Y
  \param[in]     incY
  
  \return  Dot product of X and Y.

*/
double
cblas_ddot( const int N, const double *X,
            const int incX, const double *Y, const int incY)
{
  double r  = 0.0;
  int    i;
  int    ix = OFFSET(N, incX);
  int    iy = OFFSET(N, incY);

  for (i = 0; i < N; i++) {
      r  += X[ix] * Y[iy];
      ix += incX;
      iy += incY;
  }
  
  return r;
}

/*!
  \param[in]     N
  \param[in]     X      
  \param[in]     incX

  \return Two-norm of X.
*/
double
cblas_dnrm2( const int N, const double *X, const int incX) 
{
  double
      scale = 0.0,
      ssq   = 1.0;
  int
      i,
      ix    = 0;

  if (N <= 0 || incX <= 0) return 0;
  else if (N == 1)         return fabs(X[0]);

  for (i = 0; i < N; i++) {
      const double x = X[ix];

      if (x != 0.0) {
          const double ax = fabs(x);

          if (scale < ax) {
              ssq   = 1.0 + ssq * (scale / ax) * (scale / ax);
              scale = ax;
          } else {
              ssq += (ax / scale) * (ax / scale);
          }
      }

      ix += incX;
  }
  
  return scale * sqrt(ssq);
}

/*!
  \param[in]     N
  \param[in]     alpha
  \param[in,out] X
  \param[in]     incX
*/
void
cblas_dscal(const int N, const double alpha, double *X, const int incX)
{
    int i, ix;

    if (incX <= 0) return;

    ix = OFFSET(N, incX);
    
    for (i = 0; i < N; i++) {
        X[ix] *= alpha;
        ix    += incX;
    }
}


#define ZERO   0.0
#define ONE    1.0

// ---------------------------------------------------------------------
// d2norm  returns  sqrt( a**2 + b**2 )  with precautions
// to avoid overflow.
//
// 21 Mar 1990: First version.
// ---------------------------------------------------------------------
static double
d2norm( const double a, const double b )
{
    double scale;
    const double zero = 0.0;

    scale  = fabs( a ) + fabs( b );
    if (scale == zero)
        return zero;
    else
        return scale * sqrt( (a/scale)*(a/scale) + (b/scale)*(b/scale) );
}

static void
dload( const int n, const double alpha, double x[] )
{    
    int i;
    for (i = 0; i < n; i++) x[i] = alpha;
    return;
}

void aprod(int mode, int m, int n, double x[], double y[], void *UsrWrk) {
//                If mode = 1, compute  y = y + A*x.
//                If mode = 2, compute  x = x + A(transpose)*y.
    double *A = (double *)UsrWrk;
    if (mode == 1) {
        for (int i = 0, t = 0; i < m; i++) {
            for (int j = 0; j < n; j++, t++) {
                y[i] += A[t] * x[j];
            }
        }
    }
    else {
        for (int i = 0, t = 0; i < m; i++) {
            for (int j = 0; j < n; j++, t++) {
                x[j] += A[t] * y[i];
            }
        }
    }
}

// ---------------------------------------------------------------------
// LSQR
// ---------------------------------------------------------------------
void lsqr( 
          int m,
          int n,
          void (*aprod)(int mode, int m, int n, double x[], double y[],
                        void *UsrWrk),
          double damp,
          void   *UsrWrk,
          double u[],     // len = m
          double v[],     // len = n
          double w[],     // len = n
          double x[],     // len = n
          double se[],    // len at least n.  May be NULL.
          double atol,
          double btol,
          double conlim,
          int    itnlim,
          FILE   *nout

         )
{
//     ------------------------------------------------------------------
//
//     LSQR  finds a solution x to the following problems:
//
//     1. Unsymmetric equations --    solve  A*x = b
//
//     2. Linear least squares  --    solve  A*x = b
//                                    in the least-squares sense
//
//     3. Damped least squares  --    solve  (   A    )*x = ( b )
//                                           ( damp*I )     ( 0 )
//                                    in the least-squares sense
//
//     where A is a matrix with m rows and n columns, b is an
//     m-vector, and damp is a scalar.  (All quantities are real.)
//     The matrix A is intended to be large and sparse.  It is accessed
//     by means of subroutine calls of the form
//
//                aprod ( mode, m, n, x, y, UsrWrk )
//
//     which must perform the following functions:
//
//                If mode = 1, compute  y = y + A*x.
//                If mode = 2, compute  x = x + A(transpose)*y.
//
//     The vectors x and y are input parameters in both cases.
//     If  mode = 1,  y should be altered without changing x.
//     If  mode = 2,  x should be altered without changing y.
//     The parameter UsrWrk may be used for workspace as described
//     below.
//
//     The rhs vector b is input via u, and subsequently overwritten.
//
//
//     Note:  LSQR uses an iterative method to approximate the solution.
//     The number of iterations required to reach a certain accuracy
//     depends strongly on the scaling of the problem.  Poor scaling of
//     the rows or columns of A should therefore be avoided where
//     possible.
//
//     For example, in problem 1 the solution is unaltered by
//     row-scaling.  If a row of A is very small or large compared to
//     the other rows of A, the corresponding row of ( A  b ) should be
//     scaled up or down.
//
//     In problems 1 and 2, the solution x is easily recovered
//     following column-scaling.  Unless better information is known,
//     the nonzero columns of A should be scaled so that they all have
//     the same Euclidean norm (e.g., 1.0).
//
//     In problem 3, there is no freedom to re-scale if damp is
//     nonzero.  However, the value of damp should be assigned only
//     after attention has been paid to the scaling of A.
//
//     The parameter damp is intended to help regularize
//     ill-conditioned systems, by preventing the true solution from
//     being very large.  Another aid to regularization is provided by
//     the parameter acond, which may be used to terminate iterations
//     before the computed solution becomes very large.
//
//     Note that x is not an input parameter.
//     If some initial estimate x0 is known and if damp = 0,
//     one could proceed as follows:
//
//       1. Compute a residual vector     r0 = b - A*x0.
//       2. Use LSQR to solve the system  A*dx = r0.
//       3. Add the correction dx to obtain a final solution x = x0 + dx.
//
//     This requires that x0 be available before and after the call
//     to LSQR.  To judge the benefits, suppose LSQR takes k1 iterations
//     to solve A*x = b and k2 iterations to solve A*dx = r0.
//     If x0 is "good", norm(r0) will be smaller than norm(b).
//     If the same stopping tolerances atol and btol are used for each
//     system, k1 and k2 will be similar, but the final solution x0 + dx
//     should be more accurate.  The only way to reduce the total work
//     is to use a larger stopping tolerance for the second system.
//     If some value btol is suitable for A*x = b, the larger value
//     btol*norm(b)/norm(r0)  should be suitable for A*dx = r0.
//
//     Preconditioning is another way to reduce the number of iterations.
//     If it is possible to solve a related system M*x = b efficiently,
//     where M approximates A in some helpful way
//     (e.g. M - A has low rank or its elements are small relative to
//     those of A), LSQR may converge more rapidly on the system
//           A*M(inverse)*z = b,
//     after which x can be recovered by solving M*x = z.
//
//     NOTE: If A is symmetric, LSQR should not be used!
//     Alternatives are the symmetric conjugate-gradient method (cg)
//     and/or SYMMLQ.
//     SYMMLQ is an implementation of symmetric cg that applies to
//     any symmetric A and will converge more rapidly than LSQR.
//     If A is positive definite, there are other implementations of
//     symmetric cg that require slightly less work per iteration
//     than SYMMLQ (but will take the same number of iterations).
//
//
//     Notation
//     --------
//
//     The following quantities are used in discussing the subroutine
//     parameters:
//
//     Abar   =  (   A    ),          bbar  =  ( b )
//               ( damp*I )                    ( 0 )
//
//     r      =  b  -  A*x,           rbar  =  bbar  -  Abar*x
//
//     rnorm  =  sqrt( norm(r)**2  +  damp**2 * norm(x)**2 )
//            =  norm( rbar )
//
//     relpr  =  the relative precision of floating-point arithmetic
//               on the machine being used.  On most machines,
//               relpr is about 1.0e-7 and 1.0d-16 in single and double
//               precision respectively.
//
//     LSQR  minimizes the function rnorm with respect to x.
//
//
//     Parameters
//     ----------
//
//     m       input      m, the number of rows in A.
//
//     n       input      n, the number of columns in A.
//
//     aprod   external   See above.
//
//     damp    input      The damping parameter for problem 3 above.
//                        (damp should be 0.0 for problems 1 and 2.)
//                        If the system A*x = b is incompatible, values
//                        of damp in the range 0 to sqrt(relpr)*norm(A)
//                        will probably have a negligible effect.
//                        Larger values of damp will tend to decrease
//                        the norm of x and reduce the number of 
//                        iterations required by LSQR.
//
//                        The work per iteration and the storage needed
//                        by LSQR are the same for all values of damp.
//
//     rw      workspace  Transit pointer to user's workspace.
//                        Note:  LSQR  does not explicitly use this
//                        parameter, but passes it to subroutine aprod for
//                        possible use as workspace.
//
//     u(m)    input      The rhs vector b.  Beware that u is
//                        over-written by LSQR.
//
//     v(n)    workspace
//
//     w(n)    workspace
//
//     x(n)    output     Returns the computed solution x.
//
//     se(*)   output     If m .gt. n  or  damp .gt. 0,  the system is
//             (maybe)    overdetermined and the standard errors may be
//                        useful.  (See the first LSQR reference.)
//                        Otherwise (m .le. n  and  damp = 0) they do not
//                        mean much.  Some time and storage can be saved
//                        by setting  se = NULL.  In that case, se will
//                        not be touched.
//
//                        If se is not NULL, then the dimension of se must
//                        be n or more.  se(1:n) then returns standard error
//                        estimates for the components of x.
//                        For each i, se(i) is set to the value
//                           rnorm * sqrt( sigma(i,i) / t ),
//                        where sigma(i,i) is an estimate of the i-th
//                        diagonal of the inverse of Abar(transpose)*Abar
//                        and  t = 1      if  m .le. n,
//                             t = m - n  if  m .gt. n  and  damp = 0,
//                             t = m      if  damp .ne. 0.
//
//     atol    input      An estimate of the relative error in the data
//                        defining the matrix A.  For example,
//                        if A is accurate to about 6 digits, set
//                        atol = 1.0e-6 .
//
//     btol    input      An estimate of the relative error in the data
//                        defining the rhs vector b.  For example,
//                        if b is accurate to about 6 digits, set
//                        btol = 1.0e-6 .
//
//     conlim  input      An upper limit on cond(Abar), the apparent
//                        condition number of the matrix Abar.
//                        Iterations will be terminated if a computed
//                        estimate of cond(Abar) exceeds conlim.
//                        This is intended to prevent certain small or
//                        zero singular values of A or Abar from
//                        coming into effect and causing unwanted growth
//                        in the computed solution.
//
//                        conlim and damp may be used separately or
//                        together to regularize ill-conditioned systems.
//
//                        Normally, conlim should be in the range
//                        1000 to 1/relpr.
//                        Suggested value:
//                        conlim = 1/(100*relpr)  for compatible systems,
//                        conlim = 1/(10*sqrt(relpr)) for least squares.
//
//             Note:  If the user is not concerned about the parameters
//             atol, btol and conlim, any or all of them may be set
//             to zero.  The effect will be the same as the values
//             relpr, relpr and 1/relpr respectively.
//
//     itnlim  input      An upper limit on the number of iterations.
//                        Suggested value:
//                        itnlim = n/2   for well-conditioned systems
//                                       with clustered singular values,
//                        itnlim = 4*n   otherwise.
//
//     nout    input      File number for printed output.  If positive,
//                        a summary will be printed on file nout.
//
//     istop   output     An integer giving the reason for termination:
//
//                0       x = 0  is the exact solution.
//                        No iterations were performed.
//
//                1       The equations A*x = b are probably
//                        compatible.  Norm(A*x - b) is sufficiently
//                        small, given the values of atol and btol.
//
//                2       damp is zero.  The system A*x = b is probably
//                        not compatible.  A least-squares solution has
//                        been obtained that is sufficiently accurate,
//                        given the value of atol.
//
//                3       damp is nonzero.  A damped least-squares
//                        solution has been obtained that is sufficiently
//                        accurate, given the value of atol.
//
//                4       An estimate of cond(Abar) has exceeded
//                        conlim.  The system A*x = b appears to be
//                        ill-conditioned.  Otherwise, there could be an
//                        error in subroutine aprod.
//
//                5       The iteration limit itnlim was reached.
//
//     itn     output     The number of iterations performed.
//
//     anorm   output     An estimate of the Frobenius norm of  Abar.
//                        This is the square-root of the sum of squares
//                        of the elements of Abar.
//                        If damp is small and if the columns of A
//                        have all been scaled to have length 1.0,
//                        anorm should increase to roughly sqrt(n).
//                        A radically different value for anorm may
//                        indicate an error in subroutine aprod (there
//                        may be an inconsistency between modes 1 and 2).
//
//     acond   output     An estimate of cond(Abar), the condition
//                        number of Abar.  A very high value of acond
//                        may again indicate an error in aprod.
//
//     rnorm   output     An estimate of the final value of norm(rbar),
//                        the function being minimized (see notation
//                        above).  This will be small if A*x = b has
//                        a solution.
//
//     arnorm  output     An estimate of the final value of
//                        norm( Abar(transpose)*rbar ), the norm of
//                        the residual for the usual normal equations.
//                        This should be small in all cases.  (arnorm
//                        will often be smaller than the true value
//                        computed from the output vector x.)
//
//     xnorm   output     An estimate of the norm of the final
//                        solution vector x.
//
//
//     Subroutines and functions used              
//     ------------------------------
//
//     USER               aprod
//     CBLAS              dcopy, dnrm2, dscal (see Lawson et al. below)
//
//
//     References
//     ----------
//
//     C.C. Paige and M.A. Saunders,  LSQR: An algorithm for sparse
//          linear equations and sparse least squares,
//          ACM Transactions on Mathematical Software 8, 1 (March 1982),
//          pp. 43-71.
//
//     C.C. Paige and M.A. Saunders,  Algorithm 583, LSQR: Sparse
//          linear equations and least-squares problems,
//          ACM Transactions on Mathematical Software 8, 2 (June 1982),
//          pp. 195-209.
//
//     C.L. Lawson, R.J. Hanson, D.R. Kincaid and F.T. Krogh,
//          Basic linear algebra subprograms for Fortran usage,
//          ACM Transactions on Mathematical Software 5, 3 (Sept 1979),
//          pp. 308-323 and 324-325.
//     ------------------------------------------------------------------
//
//
//     LSQR development:
//     22 Feb 1982: LSQR sent to ACM TOMS to become Algorithm 583.
//     15 Sep 1985: Final F66 version.  LSQR sent to "misc" in netlib.
//     13 Oct 1987: Bug (Robert Davies, DSIR).  Have to delete
//                     if ( (one + dabs(t)) .le. one ) GO TO 200
//                  from loop 200.  The test was an attempt to reduce
//                  underflows, but caused w(i) not to be updated.
//     17 Mar 1989: First F77 version.
//     04 May 1989: Bug (David Gay, AT&T).  When the second beta is zero,
//                  rnorm = 0 and
//                  test2 = arnorm / (anorm * rnorm) overflows.
//                  Fixed by testing for rnorm = 0.
//     05 May 1989: Sent to "misc" in netlib.
//     14 Mar 1990: Bug (John Tomlin via IBM OSL testing).
//                  Setting rhbar2 = rhobar**2 + dampsq can give zero
//                  if rhobar underflows and damp = 0.
//                  Fixed by testing for damp = 0 specially.
//     15 Mar 1990: Converted to lower case.
//     21 Mar 1990: d2norm introduced to avoid overflow in numerous
//                  items like  c = sqrt( a**2 + b**2 ).
//     04 Sep 1991: wantse added as an argument to LSQR, to make
//                  standard errors optional.  This saves storage and
//                  time when se(*) is not wanted.
//     13 Feb 1992: istop now returns a value in [1,5], not [1,7].
//                  1, 2 or 3 means that x solves one of the problems
//                  Ax = b,  min norm(Ax - b)  or  damped least squares.
//                  4 means the limit on cond(A) was reached.
//                  5 means the limit on iterations was reached.
//     07 Dec 1994: Keep track of dxmax = max_k norm( phi_k * d_k ).
//                  So far, this is just printed at the end.
//                  A large value (relative to norm(x)) indicates
//                  significant cancellation in forming
//                  x  =  D*f  =  sum( phi_k * d_k ).
//                  A large column of D need NOT be serious if the
//                  corresponding phi_k is small.
//     27 Dec 1994: Include estimate of alfa_opt in iteration log.
//                  alfa_opt is the optimal scale factor for the
//                  residual in the "augmented system", as described by
//                  A. Bjorck (1992),
//                  Pivoting and stability in the augmented system method,
//                  in D. F. Griffiths and G. A. Watson (eds.),
//                  "Numerical Analysis 1991",
//                  Proceedings of the 14th Dundee Conference,
//                  Pitman Research Notes in Mathematics 260,
//                  Longman Scientific and Technical, Harlow, Essex, 1992.
//     14 Apr 2006: "Line-by-line" conversion to ISO C by
//                  Michael P. Friedlander.
//
//
//     Michael A. Saunders                  mike@sol-michael.stanford.edu
//     Dept of Operations Research          na.Msaunders@na-net.ornl.gov
//     Stanford University
//     Stanford, CA 94305-4022              (415) 723-1875
//-----------------------------------------------------------------------

//  Local copies of output variables.  Output vars are assigned at exit.
    int
        istop  = 0,
        itn    = 0;
    double
        anorm  = ZERO,
        acond  = ZERO,
        rnorm  = ZERO,
        arnorm = ZERO,
        xnorm  = ZERO;

//  Local variables

    const bool
        extra  = false,       // true for extra printing below.
        damped = damp > ZERO,
        wantse = se != NULL;
    int
        i, maxdx, nconv, nstop;
    double
        alfopt, alpha, arnorm0, beta, bnorm,
        cs, cs1, cs2, ctol,
        delta, dknorm, dnorm, dxk, dxmax,
        gamma, gambar, phi, phibar, psi,
        res2, rho, rhobar, rhbar1,
        rhs, rtol, sn, sn1, sn2,
        t, tau, temp, test1, test2, test3,
        theta, t1, t2, t3, xnorm1, z, zbar;
    char
        enter[] = "Enter LSQR.  ",
        exit[]  = "Exit  LSQR.  ",
        msg[6][100] =
        {
            {"The exact solution is  x = 0"},
            {"A solution to Ax = b was found, given atol, btol"},
            {"A least-squares solution was found, given atol"},
            {"A damped least-squares solution was found, given atol"},
            {"Cond(Abar) seems to be too large, given conlim"},
            {"The iteration limit was reached"}
        };
//-----------------------------------------------------------------------

//  Format strings.
    char fmt_1000[] = 
        " %s        Least-squares solution of  Ax = b\n"
        " The matrix  A  has %7d rows  and %7d columns\n"
        " damp   = %-22.2e    wantse = %10i\n"
        " atol   = %-22.2e    conlim = %10.2e\n"
        " btol   = %-22.2e    itnlim = %10d\n\n";
    char fmt_1200[] =
        "    Itn       x(1)           Function"
        "     Compatible    LS      Norm A   Cond A\n";
    char fmt_1300[] =
        "    Itn       x(1)           Function"
        "     Compatible    LS      Norm Abar   Cond Abar\n";
    char fmt_1400[] =
        "     phi    dknorm  dxk  alfa_opt\n";
    char fmt_1500_extra[] =
        " %6d %16.9e %16.9e %9.2e %9.2e %8.1e %8.1e %8.1e %7.1e %7.1e %7.1e\n";
    char fmt_1500[] =
        " %6d %16.9e %16.9e %9.2e %9.2e %8.1e %8.1e\n";
    char fmt_1550[] =
        " %6d %16.9e %16.9e %9.2e %9.2e\n";
    char fmt_1600[] = 
        "\n";
    char fmt_2000[] =
        "\n"
        " %s       istop  = %-10d      itn    = %-10d\n"
        " %s       anorm  = %11.5e     acond  = %11.5e\n"
        " %s       vnorm  = %11.5e     xnorm  = %11.5e\n"
        " %s       rnorm  = %11.5e     arnorm = %11.5e\n";
    char fmt_2100[] =
        " %s       max dx = %7.1e occured at itn %-9d\n"
        " %s              = %7.1e*xnorm\n";
    char fmt_3000[] =
        " %s       %s\n";

//  Initialize.

    if (nout != NULL)
        fprintf(nout, fmt_1000,
                enter, m, n, damp, wantse,
                atol, conlim, btol, itnlim);

    itn    =   0;
    istop  =   0;
    nstop  =   0;
    maxdx  =   0;
    ctol   =   ZERO;
    if (conlim > ZERO) ctol = ONE / conlim;
    anorm  =   ZERO;
    acond  =   ZERO;
    dnorm  =   ZERO;
    dxmax  =   ZERO;
    res2   =   ZERO;
    psi    =   ZERO;
    xnorm  =   ZERO;
    xnorm1 =   ZERO;
    cs2    = - ONE;
    sn2    =   ZERO;
    z      =   ZERO;

//  ------------------------------------------------------------------
//  Set up the first vectors u and v for the bidiagonalization.
//  These satisfy  beta*u = b,  alpha*v = A(transpose)*u.
//  ------------------------------------------------------------------
    dload( n, 0.0, v );
    dload( n, 0.0, x );

    if ( wantse )
        dload( n, 0.0, se );
    
    alpha  =   ZERO;
    beta   =   cblas_dnrm2 ( m, u, 1 );

    if (beta > ZERO) {
        cblas_dscal ( m, (ONE / beta), u, 1 );
        aprod ( 2, m, n, v, u, UsrWrk );
        alpha  =   cblas_dnrm2 ( n, v, 1 );
    }

    if (alpha > ZERO) {
        cblas_dscal ( n, (ONE / alpha), v, 1 );
        cblas_dcopy ( n, v, 1, w, 1 );
    }

    arnorm = arnorm0 = alpha * beta;
    if (arnorm == ZERO) goto goto_800;
    
    rhobar =   alpha;
    phibar =   beta;
    bnorm  =   beta;
    rnorm  =   beta;

    if (nout != NULL) {
        if ( damped ) 
            fprintf(nout, fmt_1300);
        else
            fprintf(nout, fmt_1200);

        test1  = ONE;
        test2  = alpha / beta;
        
        if ( extra ) 
            fprintf(nout, fmt_1400);

        fprintf(nout, fmt_1550, itn, x[0], rnorm, test1, test2);
        fprintf(nout, fmt_1600);
    }


//  ==================================================================
//  Main iteration loop.
//  ==================================================================
    while (1) {
        itn    = itn + 1;
        
//      ------------------------------------------------------------------
//      Perform the next step of the bidiagonalization to obtain the
//      next  beta, u, alpha, v.  These satisfy the relations
//                 beta*u  =  A*v  -  alpha*u,
//                alpha*v  =  A(transpose)*u  -  beta*v.
//      ------------------------------------------------------------------
        cblas_dscal ( m, (- alpha), u, 1 );
        aprod ( 1, m, n, v, u, UsrWrk );
        beta   =   cblas_dnrm2 ( m, u, 1 );

//      Accumulate  anorm = || Bk ||
//                        =  sqrt( sum of  alpha**2 + beta**2 + damp**2 ).

        temp   =   d2norm( alpha, beta );
        temp   =   d2norm( temp , damp );
        anorm  =   d2norm( anorm, temp );

        if (beta > ZERO) {
            cblas_dscal ( m, (ONE / beta), u, 1 );
            cblas_dscal ( n, (- beta), v, 1 );
            aprod ( 2, m, n, v, u, UsrWrk );
            alpha  =   cblas_dnrm2 ( n, v, 1 );
            if (alpha > ZERO) {
                cblas_dscal ( n, (ONE / alpha), v, 1 );
            }
        }

//      ------------------------------------------------------------------
//      Use a plane rotation to eliminate the damping parameter.
//      This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
//      ------------------------------------------------------------------
        rhbar1 = rhobar;
        if ( damped ) {
            rhbar1 = d2norm( rhobar, damp );
            cs1    = rhobar / rhbar1;
            sn1    = damp   / rhbar1;
            psi    = sn1 * phibar;
            phibar = cs1 * phibar;
        }

//      ------------------------------------------------------------------
//      Use a plane rotation to eliminate the subdiagonal element (beta)
//      of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
//      ------------------------------------------------------------------
        rho    =   d2norm( rhbar1, beta );
        cs     =   rhbar1 / rho;
        sn     =   beta   / rho;
        theta  =   sn * alpha;
        rhobar = - cs * alpha;
        phi    =   cs * phibar;
        phibar =   sn * phibar;
        tau    =   sn * phi;

//      ------------------------------------------------------------------
//      Update  x, w  and (perhaps) the standard error estimates.
//      ------------------------------------------------------------------
        t1     =   phi   / rho;
        t2     = - theta / rho;
        t3     =   ONE   / rho;
        dknorm =   ZERO;

        if ( wantse ) {
            for (i = 0; i < n; i++) {
                t      =  w[i];
                x[i]   =  t1*t  +  x[i];
                w[i]   =  t2*t  +  v[i];
                t      = (t3*t)*(t3*t);
                se[i]  =  t     +  se[i];
                dknorm =  t     +  dknorm;
            }
        }
        else {
            for (i = 0; i < n; i++) {
                t      =  w[i];
                x[i]   =  t1*t  +  x[i];
                w[i]   =  t2*t  +  v[i];
                dknorm = (t3*t)*(t3*t)  +  dknorm;
            }
        }

//      ------------------------------------------------------------------
//      Monitor the norm of d_k, the update to x.
//      dknorm = norm( d_k )
//      dnorm  = norm( D_k ),        where   D_k = (d_1, d_2, ..., d_k )
//      dxk    = norm( phi_k d_k ),  where new x = x_k + phi_k d_k.
//      ------------------------------------------------------------------
        dknorm = sqrt( dknorm );
        dnorm  = d2norm( dnorm, dknorm );
        dxk    = fabs( phi * dknorm );
        if (dxmax < dxk ) {
            dxmax   =  dxk;
            maxdx   =  itn;
        }

//      ------------------------------------------------------------------
//      Use a plane rotation on the right to eliminate the
//      super-diagonal element (theta) of the upper-bidiagonal matrix.
//      Then use the result to estimate  norm(x).
//      ------------------------------------------------------------------
        delta  =   sn2 * rho;
        gambar = - cs2 * rho;
        rhs    =   phi    - delta * z;
        zbar   =   rhs    / gambar;
        xnorm  =   d2norm( xnorm1, zbar  );
        gamma  =   d2norm( gambar, theta );
        cs2    =   gambar / gamma;
        sn2    =   theta  / gamma;
        z      =   rhs    / gamma;
        xnorm1 =   d2norm( xnorm1, z     );

//      ------------------------------------------------------------------
//      Test for convergence.
//      First, estimate the norm and condition of the matrix  Abar,
//      and the norms of  rbar  and  Abar(transpose)*rbar.
//      ------------------------------------------------------------------
        acond  =   anorm * dnorm;
        res2   =   d2norm( res2 , psi    );
        rnorm  =   d2norm( res2 , phibar );
        arnorm =   alpha * fabs( tau );

//      Now use these norms to estimate certain other quantities,
//      some of which will be small near a solution.

        alfopt =   sqrt( rnorm / (dnorm * xnorm) );
        test1  =   rnorm /  bnorm;
        test2  =   ZERO;
        if (rnorm   > ZERO) test2 = arnorm / (anorm * rnorm);
//      if (arnorm0 > ZERO) test2 = arnorm / arnorm0;  //(Michael Friedlander's modification)
        test3  =   ONE   /  acond;
        t1     =   test1 / (ONE  +  anorm * xnorm / bnorm);
        rtol   =   btol  +  atol *  anorm * xnorm / bnorm;

//      The following tests guard against extremely small values of
//      atol, btol  or  ctol.  (The user may have set any or all of
//      the parameters  atol, btol, conlim  to zero.)
//      The effect is equivalent to the normal tests using
//      atol = relpr,  btol = relpr,  conlim = 1/relpr.

        t3     =   ONE + test3;
        t2     =   ONE + test2;
        t1     =   ONE + t1;
        if (itn >= itnlim) istop = 5;
        if (t3  <= ONE   ) istop = 4;
        if (t2  <= ONE   ) istop = 2;
        if (t1  <= ONE   ) istop = 1;

//      Allow for tolerances set by the user.

        if (test3 <= ctol) istop = 4;
        if (test2 <= atol) istop = 2;
        if (test1 <= rtol) istop = 1;   //(Michael Friedlander had this commented out)

//      ------------------------------------------------------------------
//      See if it is time to print something.
//      ------------------------------------------------------------------
        if (nout  == NULL     ) goto goto_600;
        if (n     <= 40       ) goto goto_400;
        if (itn   <= 10       ) goto goto_400;
        if (itn   >= itnlim-10) goto goto_400;
        if (itn % 10 == 0     ) goto goto_400;
        if (test3 <=  2.0*ctol) goto goto_400;
        if (test2 <= 10.0*atol) goto goto_400;
        if (test1 <= 10.0*rtol) goto goto_400;
        if (istop != 0        ) goto goto_400;
        goto goto_600;

//      Print a line for this iteration.
//      "extra" is for experimental purposes.

    goto_400:
        if ( extra ) {
            fprintf(nout, fmt_1500_extra,
                    itn, x[0], rnorm, test1, test2, anorm,
                    acond, phi, dknorm, dxk, alfopt);
        }
        else {
            fprintf(nout, fmt_1500,
                    itn, x[0], rnorm, test1, test2, anorm, acond);
        }
        if (itn % 10 == 0) fprintf(nout, fmt_1600);

//      ------------------------------------------------------------------
//      Stop if appropriate.
//      The convergence criteria are required to be met on  nconv
//      consecutive iterations, where  nconv  is set below.
//      Suggested value:  nconv = 1, 2  or  3.
//      ------------------------------------------------------------------
    goto_600:
        if (istop == 0) {
            nstop  = 0;
        }
        else {
            nconv  = 1;
            nstop  = nstop + 1;
            if (nstop < nconv  &&  itn < itnlim) istop = 0;
        }

        if (istop != 0) break;
        
    }
//  ==================================================================
//  End of iteration loop.
//  ==================================================================

//  Finish off the standard error estimates.

    if ( wantse ) {
        t    =   ONE;
        if (m > n)     t = m - n;
        if ( damped )  t = m;
        t    =   rnorm / sqrt( t );
      
        for (i = 0; i < n; i++)
            se[i]  = t * sqrt( se[i] );
        
    }

//  Decide if istop = 2 or 3.
//  Print the stopping condition.
 goto_800:
    if (damped  &&  istop == 2) istop = 3;
    if (nout != NULL) {
        fprintf(nout, fmt_2000,
                exit, istop, itn,
                exit, anorm, acond,
                exit, bnorm, xnorm,
                exit, rnorm, arnorm);
        fprintf(nout, fmt_2100,
                exit, dxmax, maxdx,
                exit, dxmax/(xnorm + 1.0e-20));
        fprintf(nout, fmt_3000,
                exit, msg[istop]);
    }

//  Assign output variables from local copies.
//    *istop_out  = istop;
//    *itn_out    = itn;
//    *anorm_out  = anorm;
//    *acond_out  = acond;
//    *rnorm_out  = rnorm;
//    *arnorm_out = test2;
//    *xnorm_out  = xnorm;

    return;
}



//=========lsqr end=======



#define LOCAL

//=================================sift begin====================================

// AUTORIGHTS
// Copyright (c) 2006 The Regents of the University of California
// All Rights Reserved.
// 
// Created by Andrea Vedaldi (UCLA VisionLab)
// 
// Permission to use, copy, modify, and distribute this software and its
// documentation for educational, research and non-profit purposes,
// without fee, and without a written agreement is hereby granted,
// provided that the above copyright notice, this paragraph and the
// following three paragraphs appear in all copies.
// 
// This software program and documentation are copyrighted by The Regents
// of the University of California. The software program and
// documentation are supplied "as is", without any accompanying services
// from The Regents. The Regents does not warrant that the operation of
// the program will be uninterrupted or error-free. The end-user
// understands that the program was developed for research purposes and
// is advised not to rely exclusively on the program for any reason.
// 
// This software embodies a method for which the following patent has
// been issued: "Method and apparatus for identifying scale invariant
// features in an image and use of same for locating an object in an
// image," David G. Lowe, US Patent 6,711,293 (March 23,
// 2004). Provisional application filed March 8, 1999. Asignee: The
// University of British Columbia.
// 
// IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
// FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
// INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
// ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
// CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
// BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
// MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


#define M_PI        3.14159265358979323846

#if defined (VL_USEFASTMATH)
#if defined (VL_MAC)
#define VL_FASTFLOAT float
#else
#define VL_FASTFLOAT double
#endif
#else
#define VL_FASTFLOAT float
#endif

#define VL_XEAS(x) #x
#define VL_EXPAND_AND_STRINGIFY(x) VL_XEAS(x)

/** @brief VisionLab namespace */
namespace VL {

    /** @brief Pixel data type */
    typedef float pixel_t ;

    /** @brief Floating point data type 
    **
    ** Although floats are precise enough for this applicatgion, on Intel
    ** based architecture using doubles for floating point computations
    ** turns out to be much faster.
    **/
    typedef VL_FASTFLOAT float_t ;

    /** @brief 32-bit floating data type */
    typedef float float32_t ;

    /** @brief 64-bit floating data type */
    typedef double float64_t ;

    /** @brief 32-bit integer data type */
    typedef int int32_t ;

    /** @brief 64-bit integer data type */
    typedef long long int int64_t ;

    /** @brief 32-bit unsigned integer data type */
    typedef int uint32_t ;

    /** @brief 8-bit unsigned integer data type */
    typedef char unsigned uint8_t ;

    /** @name Fast math
    ** 
    ** We provide approximate mathematical functions. These are usually
    ** rather faster than the corresponding standard library functions.
    **/
    /*@{*/
    float   fast_resqrt(float x) ;
    double  fast_resqrt(double x) ;
    float_t fast_expn(float_t x) ;
    float_t fast_abs(float_t x) ;
    float_t fast_mod_2pi(float_t x) ;
    float_t fast_atan2(float_t y, float_t x) ;
    float_t fast_sqrt(float_t x) ;
    int32_t fast_floor(float_t x) ;
    /*@}*/

    /** @brief Generic exception */
    struct
        Exception
    {
        /** @brief Build generic exception with message
        ** 
        ** The message can be accessed as the Exception::msg data member.
        **
        ** @param _msg message.
        **/
        Exception(std::string _msg) : msg(_msg) { }

        /** Exception message */
        std::string msg ; 
    } ;

    /** @brief Throw generic exception
    **
    ** The macro executes the stream operations @a x to obtain
    ** an error messages. The message is then wrapped in a
    ** generic exception VL::Exception and thrown.
    **
    ** @param x sequence of stream operations.
    **/
#define VL_THROW(x)                             \
    {                                             \
    std::ostringstream oss ;                    \
    oss << x ;                                  \
    throw VL::Exception(oss.str()) ;            \
    }

    /** @name PGM input/output */
    /*@{*/
    /** @brief PGM buffer descriptor
    **
    ** The structure describes a gray scale image and it is used by the
    ** PGM input/output functions. The fileds are self-explanatory.
    **/
    struct PgmBuffer
    {
        int width ;     ///< Image width
        int height ;    ///< Image hegith
        pixel_t* data ; ///< Image data
    } ;
    std::ostream& insertPgm(std::ostream&, pixel_t const* im, int width, int height) ;
    std::istream& extractPgm(std::istream&, PgmBuffer& buffer) ;
    /*@}*/

    /** @brief SIFT filter
    **
    ** This class is a filter computing the Scale Invariant Feature
    ** Transform (SIFT).
    **/
    class Sift
    {

    public:

        /** @brief SIFT keypoint
        **
        ** A SIFT keypoint is charactedized by a location x,y and a scale
        ** @c sigma. The scale is obtained from the level index @c s and
        ** the octave index @c o through a simple formula (see the PDF
        ** documentation).
        **
        ** In addition to the location, scale indexes and scale, we also
        ** store the integer location and level. The integer location is
        ** unnormalized, i.e. relative to the resolution of the octave
        ** containing the keypoint (octaves are downsampled). 
        **/
        struct Keypoint
        {
            int o ;    ///< Keypoint octave index

            int ix ;   ///< Keypoint integer X coordinate (unnormalized)
            int iy ;   ///< Keypoint integer Y coordinate (unnormalized)
            int is ;   ///< Keypoint integer scale indiex

            float_t x  ;  ///< Keypoint fractional X coordinate
            float_t y  ;  ///< Keypoint fractional Y coordinate
            float_t s ;   ///< Keypoint fractional scale index

            float_t sigma ;  ///< Keypoint scale
        } ; 

        typedef std::vector<Keypoint>     Keypoints ;          ///< Keypoint list datatype
        typedef Keypoints::iterator       KeypointsIter ;      ///< Keypoint list iter datatype
        typedef Keypoints::const_iterator KeypointsConstIter ; ///< Keypoint list const iter datatype

        /** @brief Constructors and destructors */
        /*@{*/
        Sift(const pixel_t* _im_pt, int _width, int _height,
            float_t _sigman,
            float_t _sigma0,
            int _O, int _S,
            int _omin, int _smin, int _smax) ;
        ~Sift() ;
        /*@}*/

        void process(const pixel_t* _im_pt, int _width, int _height) ;

        /** @brief Querying the Gaussian scale space */
        /*@{*/
        VL::pixel_t* getOctave(int o) ;
        VL::pixel_t* getLevel(int o, int s) ;
        int          getWidth() const ;
        int          getHeight() const ;
        int          getOctaveWidth(int o) const ;
        int          getOctaveHeight(int o) const ;
        VL::float_t  getOctaveSamplingPeriod(int o) const ;
        VL::float_t  getScaleFromIndex(VL::float_t o, VL::float_t s) const ;
        Keypoint     getKeypoint(VL::float_t x, VL::float_t y, VL::float_t s) const ;
        /*@}*/

        /** @brief Descriptor parameters */
        /*@{*/
        bool getNormalizeDescriptor() const ;
        void setNormalizeDescriptor(bool) ;
        void setMagnification(VL::float_t) ;
        VL::float_t getMagnification() const ;  
        /*@}*/

        /** @brief Detector and descriptor */
        /*@{*/
        void detectKeypoints(VL::float_t threshold, VL::float_t edgeThreshold) ;
        int computeKeypointOrientations(VL::float_t angles [4], Keypoint keypoint) ; 
        void computeKeypointDescriptor(VL::float_t* descr_pt, Keypoint keypoint, VL::float_t angle) ;
        KeypointsIter keypointsBegin() ;
        KeypointsIter keypointsEnd() ;
        /*@}*/

    private:
        void prepareBuffers() ;
        void freeBuffers() ;
        void smooth(VL::pixel_t       * dst, 
            VL::pixel_t       * temp, 
            VL::pixel_t const * src, int width, int height, 
            VL::float_t s) ;

        void prepareGrad(int o) ;

        // scale space parameters
        VL::float_t sigman ;
        VL::float_t sigma0 ;
        VL::float_t sigmak ;

        int O ;
        int S ; 
        int omin ;
        int smin ; 
        int smax ;

        int width ;
        int height ;

        // descriptor parameters
        VL::float_t magnif ;
        bool        normalizeDescriptor ;

        // buffers
        VL::pixel_t*  temp ;
        int           tempReserved ;
        bool          tempIsGrad  ;
        int           tempOctave ;
        VL::pixel_t** octaves ;

        VL::pixel_t*  filter ;
        int           filterReserved ;

        Keypoints keypoints ;  
    } ;


}

// Include inline functions definitions

/** 
** @file
** @brief SIFT class - inline functions and members 
**/

namespace VL
{

    namespace Detail
    {
        extern int const expnTableSize ;
        extern VL::float_t const expnTableMax ;
        extern VL::float_t expnTable [] ;
    } 

    /** @brief Get width of source image
    ** @result width.
    **/
    inline
        int     
        Sift::getWidth() const
    {
        return width ;
    }

    /** @brief Get height of source image
    ** @result height.
    **/
    inline
        int
        Sift::getHeight() const
    {
        return height ;
    }

    /** @brief Get width of an octave
    ** @param o octave index.
    ** @result width of octave @a o.
    **/
    inline
        int     
        Sift::getOctaveWidth(int o) const
    {
        assert( omin <= o && o < omin + O ) ;
        return (o >= 0) ? (width >> o) : (width << -o) ;
    }

    /** @brief Get height of an octave
    ** @param o octave index.
    ** @result height of octave @a o.
    **/
    inline
        int
        Sift::getOctaveHeight(int o) const
    {
        assert( omin <= o && o < omin + O ) ;
        return (o >= 0) ? (height >> o) : (height << -o) ;
    }

    /** @brief Get octave
    ** @param o octave index.
    ** @return pointer to octave @a o.
    **/
    inline
        VL::pixel_t * 
        Sift::getOctave(int o) 
    {
        assert( omin <= o && o < omin + O ) ;
        return octaves[o-omin] ;
    }

    /** @brief Get level
    ** @param o octave index.
    ** @param s level index.
    ** @result pointer to level @c (o,s).
    **/
    inline
        VL::pixel_t * 
        Sift::getLevel(int o, int s) 
    {
        assert( omin <= o && o <  omin + O ) ;
        assert( smin <= s && s <= smax     ) ;  
        return octaves[o - omin] +
            getOctaveWidth(o)*getOctaveHeight(o) * (s-smin) ;
    }

    /** @brief Get octave sampling period
    ** @param o octave index.
    ** @result Octave sampling period (in pixels).
    **/
    inline
        VL::float_t
        Sift::getOctaveSamplingPeriod(int o) const
    {
        return (o >= 0) ? (1 << o) : 1.0f / (1 << -o) ;
    }

    /** @brief Convert index into scale
    ** @param o octave index.
    ** @param s scale index.
    ** @return scale.
    **/
    inline
        VL::float_t
        Sift::getScaleFromIndex(VL::float_t o, VL::float_t s) const
    {
        return sigma0 * powf( 2.0f, o + s / S ) ;
    }

    /** @brief Get keypoint list begin
    ** @return iterator to the beginning.
    **/
    inline
        Sift::KeypointsIter
        Sift::keypointsBegin()
    {
        return keypoints.begin() ;
    }

    /** @brief Get keypoint list end
    ** @return iterator to the end.
    **/
    inline
        Sift::KeypointsIter
        Sift::keypointsEnd()
    {
        return keypoints.end() ;
    }

    /** @brief Set normalize descriptor flag */
    inline
        void
        Sift::setNormalizeDescriptor(bool flag)
    {
        normalizeDescriptor = flag ;
    }

    /** @brief Get normalize descriptor flag */
    inline
        bool
        Sift::getNormalizeDescriptor() const
    {
        return normalizeDescriptor ;
    }

    /** @brief Set descriptor magnification */
    inline
        void
        Sift::setMagnification(VL::float_t _magnif)
    {
        magnif = _magnif ;
    }

    /** @brief Get descriptor magnification */
    inline
        VL::float_t
        Sift::getMagnification() const
    {
        return magnif ;
    }

    /** @brief Fast @ exp(-x)
    **
    ** The argument must be in the range 0-25.0 (bigger arguments may be
    ** truncated to zero).
    **
    ** @param x argument.
    ** @return @c exp(-x)
    **/
    inline
        VL::float_t
        fast_expn(VL::float_t x)
    {
        assert(VL::float_t(0) <= x && x <= Detail::expnTableMax) ;
#ifdef VL_USEFASTMATH
        x *= Detail::expnTableSize / Detail::expnTableMax ;
        VL::int32_t i = fast_floor(x) ;
        VL::float_t r = x - i ;
        VL::float_t a = VL::Detail::expnTable[i] ;
        VL::float_t b = VL::Detail::expnTable[i+1] ;
        return a + r * (b - a) ;
#else
        return exp(-x) ;
#endif
    }

    /** @brief Fast @c mod(x,2pi)
    **
    ** The function quickly computes the value @c mod(x,2pi).
    ** 
    ** @remark The computation is fast only for arguments @a x which are
    ** small in modulus.
    **
    ** @remark For negative arguments, the semantic of the function is
    ** not equivalent to the standard library @c fmod function.
    **
    ** @param x function argument.
    ** @return @c mod(x,2pi)
    **/
    inline
        VL::float_t 
        fast_mod_2pi(VL::float_t x)
    {
#ifdef VL_USEFASTMATH
        while(x < VL::float_t(0)      ) x += VL::float_t(2*M_PI) ;
        while(x > VL::float_t(2*M_PI) ) x -= VL::float_t(2*M_PI) ;
        return x ;
#else
        return (x>=0) ? std::fmod(x, VL::float_t(2*M_PI)) 
            : 2*M_PI + std::fmod(x, VL::float_t(2*M_PI)) ;
#endif
    }

    /** @brief Fast @c (int) floor(x)
    ** @param x argument.
    ** @return @c float(x)
    **/
    inline
        int32_t 
        fast_floor(VL::float_t x)
    {
#ifdef VL_USEFASTMATH
        return (x>=0)? int32_t(x) : std::floor(x) ;
        //  return int32_t( x - ((x>=0)?0:1) ) ; 
#else
        return int32_t( std::floor(x) ) ;
#endif
    }

    /** @brief Fast @c abs(x)
    ** @param x argument.
    ** @return @c abs(x)
    **/
    inline
        VL::float_t
        fast_abs(VL::float_t x)
    {
#ifdef VL_USEFASTMATH
        return (x >= 0) ? x : -x ;
#else
        return std::fabs(x) ; 
#endif
    }

    /** @brief Fast @c atan2
    ** @param x argument.
    ** @param y argument.
    ** @return Approximation of @c atan2(x).
    **/
    inline
        VL::float_t
        fast_atan2(VL::float_t y, VL::float_t x)
    {
#ifdef VL_USEFASTMATH

        /*
        The function f(r)=atan((1-r)/(1+r)) for r in [-1,1] is easier to
        approximate than atan(z) for z in [0,inf]. To approximate f(r) to
        the third degree we may solve the system

        f(+1) = c0 + c1 + c2 + c3 = atan(0) = 0
        f(-1) = c0 - c1 + c2 - c3 = atan(inf) = pi/2
        f(0)  = c0                = atan(1) = pi/4

        which constrains the polynomial to go through the end points and
        the middle point.

        We still miss a constrain, which might be simply a constarint on
        the derivative in 0. Instead we minimize the Linf error in the
        range [0,1] by searching for an optimal value of the free
        parameter. This turns out to correspond to the solution

        c0=pi/4, c1=-0.9675, c2=0, c3=0.1821

        which has maxerr = 0.0061 rad = 0.35 grad.
        */

        VL::float_t angle, r ;
        VL::float_t const c3 = 0.1821 ;
        VL::float_t const c1 = 0.9675 ;
        VL::float_t abs_y    = fast_abs(y) + VL::float_t(1e-10) ;

        if (x >= 0) {
            r = (x - abs_y) / (x + abs_y) ;
            angle = VL::float_t(M_PI/4.0) ;
        } else {
            r = (x + abs_y) / (abs_y - x) ;
            angle = VL::float_t(3*M_PI/4.0) ;
        } 
        angle += (c3*r*r - c1) * r ; 
        return (y < 0) ? -angle : angle ;
#else
        return std::atan2(y,x) ;
#endif
    }

    /** @brief Fast @c resqrt
    ** @param x argument.
    ** @return Approximation to @c resqrt(x).
    **/
    inline
        float
        fast_resqrt(float x)
    {
#ifdef VL_USEFASTMATH
        // Works if VL::float_t is 32 bit ...
        union {
            float x ;
            VL::int32_t i ;
        } u ;
        float xhalf = float(0.5) * x ;
        u.x = x ;                               // get bits for floating value
        u.i = 0x5f3759df - (u.i>>1);            // gives initial guess y0
        //u.i = 0xdf59375f - (u.i>>1);          // gives initial guess y0
        u.x = u.x*(float(1.5) - xhalf*u.x*u.x); // Newton step (may repeat)
        u.x = u.x*(float(1.5) - xhalf*u.x*u.x); // Newton step (may repeat)
        return u.x ;
#else
        return float(1.0) / std::sqrt(x) ;
#endif
    }

    /** @brief Fast @c resqrt
    ** @param x argument.
    ** @return Approximation to @c resqrt(x).
    **/
    inline
        double
        fast_resqrt(double x)
    {
#ifdef VL_USEFASTMATH
        // Works if double is 64 bit ...
        union {
            double x ;
            VL::int64_t i ;
        } u ;
        double xhalf = double(0.5) * x ;
        u.x = x ;                                // get bits for floating value
        u.i = 0x5fe6ec85e7de30daLL - (u.i>>1);   // gives initial guess y0
        u.x = u.x*(double(1.5) - xhalf*u.x*u.x); // Newton step (may repeat)
        u.x = u.x*(double(1.5) - xhalf*u.x*u.x); // Newton step (may repeat)
        return u.x ;
#else
        return double(1.0) / std::sqrt(x) ;
#endif
    }

    /** @brief Fast @c sqrt
    ** @param x argument.
    ** @return Approximation to @c sqrt(x).
    **/
    inline
        VL::float_t
        fast_sqrt(VL::float_t x)
    {
#ifdef VL_USEFASTMATH
        return (x < 1e-8) ? 0 : x * fast_resqrt(x) ;
#else
        return std::sqrt(x) ;
#endif
    }

}





/** @brief Insert descriptor into stream
**
** The function writes a descriptor in ASCII/binary format
** and in integer/floating point format into the stream.
**
** @param os output stream.
** @param descr_pt descriptor (floating point)
** @param binary write binary descriptor?
** @param fp write floating point data?
**/
std::ostream&
    insertDescriptor(std::ostream& os,
    VL::float_t const * descr_pt,
    bool binary,
    bool fp )
{
#define RAW_CONST_PT(x) reinterpret_cast<char const*>(x)
#define RAW_PT(x)       reinterpret_cast<char*>(x)

    if( fp ) {

        /* convert to 32 bits floats (single precision) */
        VL::float32_t fdescr_pt [128] ;
        for(int i = 0 ; i < 128 ; ++i)
            fdescr_pt[i] = VL::float32_t( descr_pt[i]) ;

        if( binary ) {
            /* 
            Test for endianess. Recall: big_endian = the most significant
            byte at lower memory address.
            */
            short int const word = 0x0001 ;
            bool little_endian = RAW_CONST_PT(&word)[0] ;

            /* 
            We save in big-endian (network) order. So if this machine is
            little endiand do the appropriate conversion.
            */
            if( little_endian ) {
                for(int i = 0 ; i < 128 ; ++i) {
                    VL::float32_t tmp = fdescr_pt[ i ] ;        
                    char* pt  = RAW_PT(fdescr_pt + i) ;
                    char* spt = RAW_PT(&tmp) ;
                    pt[0] = spt[3] ;
                    pt[1] = spt[2] ;
                    pt[2] = spt[1] ;
                    pt[3] = spt[0] ;
                }
            }            
            os.write( RAW_PT(fdescr_pt), 128 * sizeof(VL::float32_t) ) ;

        } else {

            for(int i = 0 ; i < 128 ; ++i) 
                os << ' ' 
                << fdescr_pt[i] ;
        }

    } else {

        VL::uint8_t idescr_pt [128] ;

        for(int i = 0 ; i < 128 ; ++i)
            idescr_pt[i] = VL::uint8_t(float(512) * descr_pt[i]) ;

        if( binary ) {

            os.write( RAW_PT(idescr_pt), 128) ;    

        } else { 

            for(int i = 0 ; i < 128 ; ++i) 
                os << ' ' 
                << VL::uint32_t( idescr_pt[i] ) ;
        }
    }
    return os ;
}

/* keypoint list */
typedef vector<pair<VL::Sift::Keypoint,VL::float_t> > Keypoints ;

/* predicate used to order keypoints by increasing scale */
bool cmpKeypoints (Keypoints::value_type const&a,
                   Keypoints::value_type const&b) {
                       return a.first.sigma < b.first.sigma ;
}

double log2(double a){
    return log(a)/log(2);
}


template<typename T>
void
    normalize(T* filter, int W)
{
    T  acc  = 0 ; 
    T* iter = filter ;
    T* end  = filter + 2*W+1 ;
    while(iter != end) acc += *iter++ ;

    iter = filter ;
    while(iter != end) *iter++ /= acc ;
}


template<typename T>
void
    convolve(T*       dst_pt, 
    const T* src_pt, int M, int N,
    const T* filter_pt, int W)
{
    typedef T const TC ;
    // convolve along columns, save transpose
    // image is M by N 
    // buffer is N by M 
    // filter is (2*W+1) by 1
    for(int j = 0 ; j < N ; ++j) {

        int i = 0 ;

        // top
        for(; i <= std::min(W-1, M-1) ; ++i) {
            TC* start = src_pt ;
            TC* stop  = src_pt    + std::min(i+W, M-1) + 1 ;
            TC* g     = filter_pt + W-i ;
            T   acc = 0.0 ;
            while(stop != start) acc += (*g++) * (*start++) ;
            *dst_pt = acc ;
            dst_pt += N ;
        }

        // middle
        // run this for W <= i <= M-1-W, only if M >= 2*W+1
        for(; i <= M-1-W ; ++i) {
            TC* start = src_pt    + i-W ;
            TC* stop  = src_pt    + i+W + 1 ;
            TC* g     = filter_pt ;
            T   acc = 0.0 ;
            while(stop != start) acc += (*g++) * (*start++) ;
            *dst_pt = acc ;
            dst_pt += N ;
        }

        // bottom
        // run this for M-W <= i <= M-1, only if M >= 2*W+1
        for(; i <= M-1 ; ++i) {
            TC* start = src_pt    + i-W ;
            TC* stop  = src_pt    + std::min(i+W, M-1) + 1 ;
            TC* g     = filter_pt ;
            T   acc   = 0.0 ;
            while(stop != start) acc += (*g++) * (*start++) ;
            *dst_pt = acc ;
            dst_pt += N ;
        }

        // next column
        src_pt += M ;
        dst_pt -= M*N - 1 ;
    }
}

// works with symmetric filters only
template<typename T>
void
    nconvolve(T*       dst_pt, 
    const T* src_pt, int M, int N,
    const T* filter_pt, int W,
    T*       scratch_pt )
{
    typedef T const TC ;

    for(int i = 0 ; i <= W ; ++i) {
        T   acc = 0.0 ;
        TC* iter = filter_pt + std::max(W-i,  0) ;
        TC* stop = filter_pt + std::min(M-1-i,W) + W + 1 ;
        while(iter != stop) acc += *iter++ ;
        scratch_pt [i] = acc ;
    }

    for(int j = 0 ; j < N ; ++j) {

        int i = 0 ;
        // top margin
        for(; i <= std::min(W, M-1) ; ++i) {
            TC* start = src_pt ;
            TC* stop  = src_pt    + std::min(i+W, M-1) + 1 ;
            TC* g     = filter_pt + W-i ;
            T   acc = 0.0 ;
            while(stop != start) acc += (*g++) * (*start++) ;
            *dst_pt = acc / scratch_pt [i] ;
            dst_pt += N ;
        }

        // middle
        for(; i <= M-1-W ; ++i) {
            TC* start = src_pt    + i-W ;
            TC* stop  = src_pt    + i+W + 1 ;
            TC* g     = filter_pt ;
            T   acc = 0.0 ;
            while(stop != start) acc += (*g++) * (*start++) ;
            *dst_pt = acc ;
            dst_pt += N ;
        }

        // bottom
        for(; i <= M-1 ; ++i) {
            TC* start = src_pt    + i-W ;
            TC* stop  = src_pt    + std::min(i+W, M-1) + 1 ;
            TC* g     = filter_pt ;
            T   acc   = 0.0 ;
            while(stop != start) acc += (*g++) * (*start++) ;
            *dst_pt = acc / scratch_pt [M-1-i];
            dst_pt += N ;
        }

        // next column
        src_pt += M ;
        dst_pt -= M*N - 1 ;
    }
}

template<typename T>
void
    econvolve(T*       dst_pt, 
    const T* src_pt,    int M, int N,
    const T* filter_pt, int W)
{
    typedef T const TC ;
    // convolve along columns, save transpose
    // image is M by N 
    // buffer is N by M 
    // filter is (2*W+1) by 1
    for(int j = 0 ; j < N ; ++j) {
        for(int i = 0 ; i < M ; ++i) {
            T   acc = 0.0 ;
            TC* g = filter_pt ;
            TC* start = src_pt + (i-W) ;
            TC* stop  ;
            T   x ;

            // beginning
            stop = src_pt + std::max(0, i-W) ;
            x    = *stop ;
            while( start <= stop ) { acc += (*g++) * x ; start++ ; }

            // middle
            stop =  src_pt + std::min(M-1, i+W) ;
            while( start <  stop ) acc += (*g++) * (*start++) ;

            // end
            x  = *start ;
            stop = src_pt + (i+W) ;
            while( start <= stop ) { acc += (*g++) * x ; start++ ; } 

            // save 
            *dst_pt = acc ; 
            dst_pt += N ;

            assert( g - filter_pt == 2*W+1 ) ;

        }
        // next column
        src_pt += M ;
        dst_pt -= M*N - 1 ;
    }
}


using namespace VL ;

// on startup, pre-compute expn(x) = exp(-x)
namespace VL { 
    namespace Detail {

        int const         expnTableSize = 256 ;
        VL::float_t const expnTableMax  = VL::float_t(25.0) ;
        VL::float_t       expnTable [ expnTableSize + 1 ] ;

        struct buildExpnTable
        {
            buildExpnTable() {
                for(int k = 0 ; k < expnTableSize + 1 ; ++k) {
                    expnTable[k] = exp( - VL::float_t(k) / expnTableSize * expnTableMax ) ;
                }
            }
        } _buildExpnTable ;

    } }


namespace VL {
    namespace Detail {

        /** Comment eater istream manipulator */
        class _cmnt {} cmnt ;

        /** @brief Extract a comment from a stream
        **
        ** The function extracts a block of consecutive comments from an
        ** input stream. A comment is a sequence of whitespaces, followed by
        ** a `#' character, other characters and terminated at the next line
        ** ending. A block of comments is just a sequence of comments.
        **/
        std::istream& 
            operator>>(std::istream& is, _cmnt& manip)
        {
            char c ;
            char b [1024] ; 
            is>>c ;
            if( c != '#' ) 
                return is.putback(c) ;
            is.getline(b,1024) ;
            return is ;
        }

    }

    /** @brief Insert PGM file into stream
    **
    ** The function iserts into the stream @a os the grayscale image @a
    ** im encoded as a PGM file. The immage is assumed to be normalized
    ** in the range 0.0 - 1.0.
    **
    ** @param os output stream.
    ** @param im pointer to image data.
    ** @param width image width.
    ** @param height image height.
    ** @return the stream @a os.
    **/
    std::ostream& 
        insertPgm(std::ostream& os, pixel_t const* im, int width, int height)
    {
        os<< "P5"   << "\n"
            << width  << " "
            << height << "\n"
            << "255"  << "\n" ;
        for(int y = 0 ; y < height ; ++y) {
            for(int x = 0 ; x < width ; ++x) {
                unsigned char v = 
                    (unsigned char)
                    (std::max(std::min(*im++, 1.0f),0.f) * 255.0f) ;
                os << v ;
            }
        }
        return os ;
    }

    void extractVec(vector<int> img, PgmBuffer& buffer){
        pixel_t* im_pt ;

        int height = img[0];
        int width = img[1];
        im_pt = new pixel_t [width*height];

        pixel_t* start = im_pt ;
        pixel_t* end   = start + width*height ; 

        for(size_t t = 2; t < img.size(); t++){
            float i = ((img[t]>>16)&255)*0.299 + ((img[t]>>8)&255)*0.587 + (img[t]&255)*0.114;
            *start++ = pixel_t( i );        
        }

        buffer.width  = width ;
        buffer.height = height ;
        buffer.data   = im_pt ;
    }


    /** @brief Extract PGM file from stream.
    **
    ** The function extracts from the stream @a in a grayscale image
    ** encoded as a PGM file. The function fills the structure @a buffer,
    ** containing the image dimensions and a pointer to the image data.
    **
    ** The image data is an array of floats and is owned by the caller,
    ** which should erase it as in
    ** 
    ** @code
    **   delete [] buffer.data.
    ** @endcode
    **
    ** When the function encouters an error it throws a generic instance
    ** of VL::Exception.
    **
    ** @param in input stream.
    ** @param buffer buffer descriptor to be filled.
    ** @return the stream @a in.
    **/
    std::istream& 
        extractPgm(std::istream& in, PgmBuffer& buffer)
    {
        pixel_t* im_pt ;
        int      width ;
        int      height ;
        int      maxval ;

        char c ;
        in>>c ;
        if( c != 'P') VL_THROW("File is not in PGM format") ;

        bool is_ascii ;
        in>>c ;
        switch( c ) {
        case '2' : is_ascii = true ; break ;
        case '5' : is_ascii = false ; break ;
        default  : VL_THROW("File is not in PGM format") ;
        }

        in >> Detail::cmnt
            >> width
            >> Detail::cmnt 
            >> height
            >> Detail::cmnt
            >> maxval ;

        // after maxval no more comments, just a whitespace or newline
        {char trash ; in.get(trash) ;}

        if(maxval > 255)
            VL_THROW("Only <= 8-bit per channel PGM files are supported") ;

        if(! in.good()) 
            VL_THROW("PGM header parsing error") ;

        im_pt = new pixel_t [ width*height ];

        try {
            if( is_ascii ) {
                pixel_t* start = im_pt ;
                pixel_t* end   = start + width*height ; 
                pixel_t  norm  = pixel_t( maxval ) ;

                while( start != end ) {        
                    int i ;
                    in >> i ;    
                    if( ! in.good() ) VL_THROW
                        ("PGM parsing error file (width="<<width
                        <<" height="<<height
                        <<" maxval="<<maxval
                        <<" at pixel="<<start-im_pt<<")") ;    
                    *start++ = pixel_t( i ) / norm ;        
                }
            } else {
                std::streampos beg = in.tellg() ;
                char* buffer = new char [width*height] ;
                in.read(buffer, width*height) ;
                if( ! in.good() ) VL_THROW
                    ("PGM parsing error file (width="<<width
                    <<" height="<<height
                    <<" maxval="<<maxval
                    <<" at pixel="<<in.tellg()-beg<<")") ;

                pixel_t* start = im_pt ;
                pixel_t* end   = start + width*height ; 
                uint8_t* src = reinterpret_cast<uint8_t*>(buffer) ;      
                while( start != end ) *start++ = *src++ / 255.0f ;
            }       
        } catch(...) {
            delete [] im_pt ; 
            throw ;
        }

        buffer.width  = width ;
        buffer.height = height ;
        buffer.data   = im_pt ;

        return in ;
    }

    // ===================================================================
    //                                          Low level image operations
    // -------------------------------------------------------------------

    namespace Detail {

        /** @brief Copy an image
        ** @param dst    output imgage buffer.
        ** @param src    input image buffer.
        ** @param width  input image width.
        ** @param height input image height.
        **/
        void
            copy(pixel_t* dst, pixel_t const* src, int width, int height)
        {
            memcpy(dst, src, sizeof(pixel_t)*width*height)  ;
        }

        /** @brief Copy an image upsampling two times
        **
        ** The destination buffer must be at least as big as two times the
        ** input buffer. Bilinear interpolation is used.
        **
        ** @param dst     output imgage buffer.
        ** @param src     input image buffer.
        ** @param width   input image width.
        ** @param height  input image height.
        **/
        void 
            copyAndUpsampleRows
            (pixel_t* dst, pixel_t const* src, int width, int height)
        {
            for(int y = 0 ; y < height ; ++y) {
                pixel_t b, a ;
                b = a = *src++ ;
                for(int x = 0 ; x < width-1 ; ++x) {
                    b = *src++ ;
                    *dst = a ;         dst += height ;
                    *dst = 0.5*(a+b) ; dst += height ;
                    a = b ;
                }
                *dst = b ; dst += height ;
                *dst = b ; dst += height ;
                dst += 1 - width * 2 * height ;
            }  
        }

        /** @brief Copy and downasample an image
        **
        ** The image is downsampled @a d times, i.e. reduced to @c 1/2^d of
        ** its original size. The parameters @a width and @a height are the
        ** size of the input image. The destination image is assumed to be @c
        ** floor(width/2^d) pixels wide and @c floor(height/2^d) pixels high.
        **
        ** @param dst output imgage buffer.
        ** @param src input image buffer.
        ** @param width input image width.
        ** @param height input image height.
        ** @param d downsampling factor.
        **/
        void 
            copyAndDownsample(pixel_t* dst, pixel_t const* src, 
            int width, int height, int d)
        {
            for(int y = 0 ; y < height ; y+=d) {
                pixel_t const * srcrowp = src + y * width ;    
                for(int x = 0 ; x < width - (d-1) ; x+=d) {     
                    *dst++ = *srcrowp ;
                    srcrowp += d ;
                }
            }
        }

    }

    /** @brief Smooth an image 
    **
    ** The function convolves the image @a src by a Gaussian kernel of
    ** variance @a s and writes the result to @a dst. The function also
    ** needs a scratch buffer @a dst of the same size of @a src and @a
    ** dst.
    **
    ** @param dst output image buffer.
    ** @param temp scratch image buffer.
    ** @param src input image buffer.
    ** @param width width of the buffers.
    ** @param height height of the buffers.
    ** @param s standard deviation of the Gaussian kernel.
    **/
    void
        Sift::smooth
        (pixel_t* dst, pixel_t* temp, 
        pixel_t const* src, int width, int height, 
        VL::float_t s)
    {
        // make sure a buffer larege enough has been allocated
        // to hold the filter
        int W = int( ceil( VL::float_t(4.0) * s ) ) ;
        if( ! filter ) {
            filterReserved = 0 ;
        }

        if( filterReserved < W ) {
            filterReserved = W ;
            if( filter ) delete [] filter ;
            filter = new pixel_t [ 2* filterReserved + 1 ] ;
        }

        // pre-compute filter
        for(int j = 0 ; j < 2*W+1 ; ++j) 
            filter[j] = VL::pixel_t
            (std::exp
            (VL::float_t
            (-0.5 * (j-W) * (j-W) / (s*s) ))) ;

        // normalize to one
        normalize(filter, W) ;

        // convolve
        econvolve(temp, src, width, height, filter, W) ;
        econvolve(dst, temp, height, width, filter, W) ;
    }

    // ===================================================================
    //                                                     Sift(), ~Sift()
    // -------------------------------------------------------------------

    /** @brief Initialize Gaussian scale space parameters
    **
    ** @param _im_pt  Source image data
    ** @param _width  Soruce image width
    ** @param _height Soruce image height
    ** @param _sigman Nominal smoothing value of the input image.
    ** @param _sigma0 Base smoothing level.
    ** @param _O      Number of octaves.
    ** @param _S      Number of levels per octave.
    ** @param _omin   First octave.
    ** @param _smin   First level in each octave.
    ** @param _smax   Last level in each octave.
    **/
    Sift::Sift(const pixel_t* _im_pt, int _width, int _height,
        VL::float_t _sigman,
        VL::float_t _sigma0,
        int _O, int _S,
        int _omin, int _smin, int _smax)
        : sigman( _sigman ), 
        sigma0( _sigma0 ),
        O( _O ),
        S( _S ),
        omin( _omin ),
        smin( _smin ),
        smax( _smax ),

        magnif( 3.0f ),
        normalizeDescriptor( true ),

        temp( NULL ),
        octaves( NULL ),
        filter( NULL )    
    {
        process(_im_pt, _width, _height) ;
    }

    /** @brief Destroy SIFT filter.
    **/
    Sift::~Sift()
    {
        freeBuffers() ;
    }

    /** Allocate buffers. Buffer sizes depend on the image size and the
    ** value of omin.
    **/
    void
        Sift::
        prepareBuffers()
    {
        // compute buffer size
        int w = (omin >= 0) ? (width  >> omin) : (width  << -omin) ;
        int h = (omin >= 0) ? (height >> omin) : (height << -omin) ;
        int size = w*h* std::max
            ((smax - smin), 2*((smax+1) - (smin-2) +1)) ;

        if( temp && tempReserved == size ) return ;

        freeBuffers() ;

        // allocate
        temp           = new pixel_t [ size ] ; 
        tempReserved   = size ;
        tempIsGrad     = false ;
        tempOctave     = 0 ;

        octaves = new pixel_t* [ O ] ;
        for(int o = 0 ; o < O ; ++o) {
            octaves[o] = new pixel_t [ (smax - smin + 1) * w * h ] ;
            w >>= 1 ;
            h >>= 1 ;
        }
    }

    /** @brief Free buffers.
    **
    ** This function releases any buffer allocated by prepareBuffers().
    **
    ** @sa prepareBuffers().
    **/
    void
        Sift::
        freeBuffers()
    {
        if( filter ) {
            delete [] filter ;
        }
        filter = 0 ;

        if( octaves ) {
            for(int o = 0 ; o < O ; ++o) {
                delete [] octaves[ o ] ;
            }
            delete [] octaves ;
        }
        octaves = 0 ;

        if( temp ) {
            delete [] temp ;   
        }
        temp = 0  ; 
    }

    // ===================================================================
    //                                                         getKeypoint
    // -------------------------------------------------------------------

    /** @brief Get keypoint from position and scale
    **
    ** The function returns a keypoint with a given position and
    ** scale. Note that the keypoint structure contains fields that make
    ** sense only in conjunction with a specific scale space. Therefore
    ** the keypoint structure should be re-calculated whenever the filter
    ** is applied to a new image, even if the parameters @a x, @a y and
    ** @a sigma do not change.
    **
    ** @param x x coordinate of the center.
    ** @peram y y coordinate of the center.
    ** @param sigma scale.
    ** @return Corresponing keypoint.
    **/
    Sift::Keypoint
        Sift::getKeypoint(VL::float_t x, VL::float_t y, VL::float_t sigma) const
    {

        /*
        The formula linking the keypoint scale sigma to the octave and
        scale index is

        (1) sigma(o,s) = sigma0 2^(o+s/S)

        for which

        (2) o + s/S = log2 sigma/sigma0 == phi.

        In addition to the scale index s (which can be fractional due to
        scale interpolation) a keypoint has an integer scale index is too
        (which is the index of the scale level where it was detected in
        the DoG scale space). We have the constraints:

        - o and is are integer

        - is is in the range [smin+1, smax-2  ]

        - o  is in the range [omin,   omin+O-1]

        - is = rand(s) most of the times (but not always, due to the way s
        is obtained by quadratic interpolation of the DoG scale space).

        Depending on the values of smin and smax, often (2) has multiple
        solutions is,o that satisfy all constraints.  In this case we
        choose the one with biggest index o (this saves a bit of
        computation).

        DETERMINING THE OCTAVE INDEX O

        From (2) we have o = phi - s/S and we want to pick the biggest
        possible index o in the feasible range. This corresponds to
        selecting the smallest possible index s. We write s = is + ds
        where in most cases |ds|<.5 (but in general |ds|<1). So we have

        o = phi - s/S,   s = is + ds ,   |ds| < .5 (or |ds| < 1).

        Since is is in the range [smin+1,smax-2], s is in the range
        [smin+.5,smax-1.5] (or [smin,smax-1]), the number o is an integer
        in the range phi+[-smax+1.5,-smin-.5] (or
        phi+[-smax+1,-smin]). Thus the maximum value of o is obtained for
        o = floor(phi-smin-.5) (or o = floor(phi-smin)).

        Finally o is clamped to make sure it is contained in the feasible
        range.

        DETERMINING THE SCALE INDEXES S AND IS

        Given o we can derive is by writing (2) as

        s = is + ds = S(phi - o).

        We then take is = round(s) and clamp its value to be in the
        feasible range.
        */

        int o,ix,iy,is ;
        VL::float_t s,phi ;

        phi = log(sigma/sigma0)/log(2) ;
        o   = fast_floor( phi -  (VL::float_t(smin)+.5)/S ) ;
        o   = std::min(o, omin+O-1) ;
        o   = std::max(o, omin    ) ;
        s   = S * (phi - o) ;

        is  = int(s + 0.5) ;
        is  = std::min(is, smax - 2) ;
        is  = std::max(is, smin + 1) ;

        VL::float_t per = getOctaveSamplingPeriod(o) ;
        ix = int(x / per + 0.5) ;
        iy = int(y / per + 0.5) ;

        Keypoint key ;
        key.o  = o ;

        key.ix = ix ;
        key.iy = iy ;
        key.is = is ;

        key.x = x ;
        key.y = y ;
        key.s = s ;

        key.sigma = sigma ;

        return key ;
    }

    // ===================================================================
    //                                                           process()
    // -------------------------------------------------------------------

    /** @brief Compute Gaussian Scale Space
    **
    ** The method computes the Gaussian scale space of the specified
    ** image. The scale space data is managed internally and can be
    ** accessed by means of getOctave() and getLevel().
    **
    ** @remark Calling this method will delete the list of keypoints
    ** constructed by detectKeypoints().
    **
    ** @param _im_pt pointer to image data.
    ** @param _width image width.
    ** @param _height image height .
    **/
    void
        Sift::
        process(const pixel_t* _im_pt, int _width, int _height)
    {
        using namespace Detail ;

        width  = _width ;
        height = _height ;
        prepareBuffers() ;

        VL::float_t sigmak = powf(2.0f, 1.0 / S) ;
        VL::float_t dsigma0 = sigma0 * sqrt (1.0f - 1.0f / (sigmak*sigmak) ) ;

        // -----------------------------------------------------------------
        //                                                 Make pyramid base
        // -----------------------------------------------------------------
        if( omin < 0 ) {
            copyAndUpsampleRows(temp,       _im_pt, width,  height  ) ;
            copyAndUpsampleRows(octaves[0], temp,   height, 2*width ) ;      

            for(int o = -1 ; o > omin ; --o) {
                copyAndUpsampleRows(temp,       octaves[0], width  << -o,    height << -o) ;
                copyAndUpsampleRows(octaves[0], temp,       height << -o, 2*(width  << -o)) ;             }

        } else if( omin > 0 ) {
            copyAndDownsample(octaves[0], _im_pt, width, height, 1 << omin) ;
        } else {
            copy(octaves[0], _im_pt, width, height) ;
        }

        {
            VL::float_t sa = sigma0 * powf(sigmak, smin) ; 
            VL::float_t sb = sigman / powf(2.0f,   omin) ; // review this
            if( sa > sb ) {
                VL::float_t sd = sqrt ( sa*sa - sb*sb ) ;
                smooth( octaves[0], temp, octaves[0], 
                    getOctaveWidth(omin),
                    getOctaveHeight(omin), 
                    sd ) ;
            }
        }

        // -----------------------------------------------------------------
        //                                                      Make octaves
        // -----------------------------------------------------------------
        for(int o = omin ; o < omin+O ; ++o) {
            // Prepare octave base
            if( o > omin ) {
                int sbest = std::min(smin + S, smax) ;
                copyAndDownsample(getLevel(o,   smin ), 
                    getLevel(o-1, sbest),
                    getOctaveWidth(o-1),
                    getOctaveHeight(o-1), 2 ) ;
                VL::float_t sa = sigma0 * powf(sigmak, smin      ) ;
                VL::float_t sb = sigma0 * powf(sigmak, sbest - S ) ;
                if(sa > sb ) {
                    VL::float_t sd = sqrt ( sa*sa - sb*sb ) ;
                    smooth( getLevel(o,0), temp, getLevel(o,0), 
                        getOctaveWidth(o), getOctaveHeight(o),
                        sd ) ;
                }
            }

            // Make other levels
            for(int s = smin+1 ; s <= smax ; ++s) {
                VL::float_t sd = dsigma0 * powf(sigmak, s) ;
                smooth( getLevel(o,s), temp, getLevel(o,s-1),
                    getOctaveWidth(o), getOctaveHeight(o),
                    sd ) ;
            }
        }
    }

    /** @brief Sift detector
    **
    ** The function runs the SIFT detector on the stored Gaussian scale
    ** space (see process()). The detector consists in three steps
    **
    ** - local maxima detection;
    ** - subpixel interpolation;
    ** - rejection of weak keypoints (@a threhsold);
    ** - rejection of keypoints on edge-like structures (@a edgeThreshold).
    **
    ** As they are found, keypoints are added to an internal list.  This
    ** list can be accessed by means of the member functions
    ** getKeypointsBegin() and getKeypointsEnd(). The list is ordered by
    ** octave, which is usefult to speed-up computeKeypointOrientations()
    ** and computeKeypointDescriptor().
    **/
    void
        Sift::detectKeypoints(VL::float_t threshold, VL::float_t edgeThreshold)
    {
        keypoints.clear() ;

        int nValidatedKeypoints = 0 ;

        // Process one octave per time
        for(int o = omin ; o < omin + O ; ++o) {

            int const xo = 1 ;
            int const yo = getOctaveWidth(o) ;
            int const so = getOctaveWidth(o) * getOctaveHeight(o) ;
            int const ow = getOctaveWidth(o) ;
            int const oh = getOctaveHeight(o) ;

            VL::float_t xperiod = getOctaveSamplingPeriod(o) ;

            // -----------------------------------------------------------------
            //                                           Difference of Gaussians
            // -----------------------------------------------------------------
            pixel_t* dog = temp ;
            tempIsGrad = false ;
            {
                pixel_t* pt = dog ;
                for(int s = smin ; s <= smax-1 ; ++s) {
                    pixel_t* srca = getLevel(o, s  ) ;
                    pixel_t* srcb = getLevel(o, s+1) ;
                    pixel_t* enda = srcb ;
                    while( srca != enda ) {
                        *pt++ = *srcb++ - *srca++ ;
                    }
                }
            }

            // -----------------------------------------------------------------
            //                                           Find points of extremum
            // -----------------------------------------------------------------
            {
                pixel_t* pt  = dog + xo + yo + so ;
                for(int s = smin+1 ; s <= smax-2 ; ++s) {
                    for(int y = 1 ; y < oh - 1 ; ++y) {
                        for(int x = 1 ; x < ow - 1 ; ++x) {          
                            pixel_t v = *pt ;

                            // assert( (pt - x*xo - y*yo - (s-smin)*so) - dog == 0 ) ;

#define CHECK_NEIGHBORS(CMP,SGN)                    \
    ( v CMP ## = SGN 0.8 * threshold &&     \
    v CMP *(pt + xo) &&                   \
    v CMP *(pt - xo) &&                   \
    v CMP *(pt + so) &&                   \
    v CMP *(pt - so) &&                   \
    v CMP *(pt + yo) &&                   \
    v CMP *(pt - yo) &&                   \
    \
    v CMP *(pt + yo + xo) &&              \
    v CMP *(pt + yo - xo) &&              \
    v CMP *(pt - yo + xo) &&              \
    v CMP *(pt - yo - xo) &&              \
    \
    v CMP *(pt + xo      + so) &&         \
    v CMP *(pt - xo      + so) &&         \
    v CMP *(pt + yo      + so) &&         \
    v CMP *(pt - yo      + so) &&         \
    v CMP *(pt + yo + xo + so) &&         \
    v CMP *(pt + yo - xo + so) &&         \
    v CMP *(pt - yo + xo + so) &&         \
    v CMP *(pt - yo - xo + so) &&         \
    \
    v CMP *(pt + xo      - so) &&         \
    v CMP *(pt - xo      - so) &&         \
    v CMP *(pt + yo      - so) &&         \
    v CMP *(pt - yo      - so) &&         \
    v CMP *(pt + yo + xo - so) &&         \
    v CMP *(pt + yo - xo - so) &&         \
    v CMP *(pt - yo + xo - so) &&         \
    v CMP *(pt - yo - xo - so) )

                            if( CHECK_NEIGHBORS(>,+) || CHECK_NEIGHBORS(<,-) ) {

                                Keypoint k ;
                                k.ix = x ;
                                k.iy = y ;
                                k.is = s ;
                                keypoints.push_back(k) ;
                            }
                            pt += 1 ;
                        }
                        pt += 2 ;
                    }
                    pt += 2*yo ;
                }
            }

            // -----------------------------------------------------------------
            //                                               Refine local maxima
            // -----------------------------------------------------------------
            { // refine
                KeypointsIter siter ;
                KeypointsIter diter ;

                for(diter = siter = keypointsBegin() + nValidatedKeypoints ; 
                    siter != keypointsEnd() ; 
                    ++siter) {

                        int x = int( siter->ix ) ;
                        int y = int( siter->iy ) ;
                        int s = int( siter->is ) ;

                        VL::float_t Dx=0,Dy=0,Ds=0,Dxx=0,Dyy=0,Dss=0,Dxy=0,Dxs=0,Dys=0 ;
                        VL::float_t  b [3] ;
                        pixel_t* pt ;
                        int dx = 0 ;
                        int dy = 0 ;

                        // must be exec. at least once
                        for(int iter = 0 ; iter < 5 ; ++iter) {

                            VL::float_t A[3*3] ;          

                            x += dx ;
                            y += dy ;

                            pt = dog 
                                + xo * x
                                + yo * y
                                + so * (s - smin) ;

#define at(dx,dy,ds) (*( pt + (dx)*xo + (dy)*yo + (ds)*so))
#define Aat(i,j)     (A[(i)+(j)*3])    

                            /* Compute the gradient. */
                            Dx = 0.5 * (at(+1,0,0) - at(-1,0,0)) ;
                            Dy = 0.5 * (at(0,+1,0) - at(0,-1,0));
                            Ds = 0.5 * (at(0,0,+1) - at(0,0,-1)) ;

                            /* Compute the Hessian. */
                            Dxx = (at(+1,0,0) + at(-1,0,0) - 2.0 * at(0,0,0)) ;
                            Dyy = (at(0,+1,0) + at(0,-1,0) - 2.0 * at(0,0,0)) ;
                            Dss = (at(0,0,+1) + at(0,0,-1) - 2.0 * at(0,0,0)) ;

                            Dxy = 0.25 * ( at(+1,+1,0) + at(-1,-1,0) - at(-1,+1,0) - at(+1,-1,0) ) ;
                            Dxs = 0.25 * ( at(+1,0,+1) + at(-1,0,-1) - at(-1,0,+1) - at(+1,0,-1) ) ;
                            Dys = 0.25 * ( at(0,+1,+1) + at(0,-1,-1) - at(0,-1,+1) - at(0,+1,-1) ) ;

                            /* Solve linear system. */
                            Aat(0,0) = Dxx ;
                            Aat(1,1) = Dyy ;
                            Aat(2,2) = Dss ;
                            Aat(0,1) = Aat(1,0) = Dxy ;
                            Aat(0,2) = Aat(2,0) = Dxs ;
                            Aat(1,2) = Aat(2,1) = Dys ;

                            b[0] = - Dx ;
                            b[1] = - Dy ;
                            b[2] = - Ds ;

                            // Gauss elimination
                            for(int j = 0 ; j < 3 ; ++j) {

                                // look for leading pivot
                                VL::float_t maxa = 0 ;
                                VL::float_t maxabsa = 0 ;
                                int   maxi = -1 ;
                                int i ;
                                for(i = j ; i < 3 ; ++i) {
                                    VL::float_t a    = Aat(i,j) ;
                                    VL::float_t absa = fabsf( a ) ;
                                    if ( absa > maxabsa ) {
                                        maxa    = a ;
                                        maxabsa = absa ;
                                        maxi    = i ;
                                    }
                                }

                                // singular?
                                if( maxabsa < 1e-10f ) {
                                    b[0] = 0 ;
                                    b[1] = 0 ;
                                    b[2] = 0 ;
                                    break ;
                                }

                                i = maxi ;

                                // swap j-th row with i-th row and
                                // normalize j-th row
                                for(int jj = j ; jj < 3 ; ++jj) {
                                    std::swap( Aat(j,jj) , Aat(i,jj) ) ;
                                    Aat(j,jj) /= maxa ;
                                }
                                std::swap( b[j], b[i] ) ;
                                b[j] /= maxa ;

                                // elimination
                                for(int ii = j+1 ; ii < 3 ; ++ii) {
                                    VL::float_t x = Aat(ii,j) ;
                                    for(int jj = j ; jj < 3 ; ++jj) {
                                        Aat(ii,jj) -= x * Aat(j,jj) ;                
                                    }
                                    b[ii] -= x * b[j] ;
                                }
                            }

                            // backward substitution
                            for(int i = 2 ; i > 0 ; --i) {
                                VL::float_t x = b[i] ;
                                for(int ii = i-1 ; ii >= 0 ; --ii) {
                                    b[ii] -= x * Aat(ii,i) ;
                                }
                            }

                            /* If the translation of the keypoint is big, move the keypoint
                            * and re-iterate the computation. Otherwise we are all set.
                            */
                            dx= ((b[0] >  0.6 && x < ow-2) ?  1 : 0 )
                                + ((b[0] < -0.6 && x > 1   ) ? -1 : 0 ) ;

                            dy= ((b[1] >  0.6 && y < oh-2) ?  1 : 0 )
                                + ((b[1] < -0.6 && y > 1   ) ? -1 : 0 ) ;

                            /*          
                            std::cout<<x<<","<<y<<"="<<at(0,0,0)
                            <<"("
                            <<at(0,0,0)+0.5 * (Dx * b[0] + Dy * b[1] + Ds * b[2])<<")"
                            <<" "<<std::flush ; 
                            */

                            if( dx == 0 && dy == 0 ) break ;
                        }

                        /* std::cout<<std::endl ; */

                        // Accept-reject keypoint
                        {
                            VL::float_t val = at(0,0,0) + 0.5 * (Dx * b[0] + Dy * b[1] + Ds * b[2]) ; 
                            VL::float_t score = (Dxx+Dyy)*(Dxx+Dyy) / (Dxx*Dyy - Dxy*Dxy) ; 
                            VL::float_t xn = x + b[0] ;
                            VL::float_t yn = y + b[1] ;
                            VL::float_t sn = s + b[2] ;

                            if(fast_abs(val) > threshold &&
                                score < (edgeThreshold+1)*(edgeThreshold+1)/edgeThreshold && 
                                score >= 0 &&
                                fast_abs(b[0]) < 1.5 &&
                                fast_abs(b[1]) < 1.5 &&
                                fast_abs(b[2]) < 1.5 &&
                                xn >= 0    &&
                                xn <= ow-1 &&
                                yn >= 0    &&
                                yn <= oh-1 &&
                                sn >= smin &&
                                sn <= smax ) {

                                    diter->o  = o ;

                                    diter->ix = x ;
                                    diter->iy = y ;
                                    diter->is = s ;

                                    diter->x = xn * xperiod ; 
                                    diter->y = yn * xperiod ; 
                                    diter->s = sn ;

                                    diter->sigma = getScaleFromIndex(o,sn) ;

                                    ++diter ;
                            }
                        }
                } // next candidate keypoint

                // prepare for next octave
                keypoints.resize( diter - keypoints.begin() ) ;
                nValidatedKeypoints = keypoints.size() ;
            } // refine block

        } // next octave
    }

    // ===================================================================
    //                                       computeKeypointOrientations()
    // -------------------------------------------------------------------

    /** @brief Compute modulus and phase of the gradient
    **
    ** The function computes the modulus and the angle of the gradient of
    ** the specified octave @a o. The result is stored in a temporary
    ** internal buffer accessed by computeKeypointDescriptor() and
    ** computeKeypointOrientations().
    **
    ** The SIFT detector provides keypoint with scale index s in the
    ** range @c smin+1 and @c smax-2. As such, the buffer contains only
    ** these levels.
    **
    ** If called mutliple time on the same data, the function exits
    ** immediately.
    **
    ** @param o octave of interest.
    **/
    void
        Sift::prepareGrad(int o)
    { 
        int const ow = getOctaveWidth(o) ;
        int const oh = getOctaveHeight(o) ;
        int const xo = 1 ;
        int const yo = ow ;
        int const so = oh*ow ;

        if( ! tempIsGrad || tempOctave != o ) {

            // compute dx/dy
            for(int s = smin+1 ; s <= smax-2 ; ++s) {
                for(int y = 1 ; y < oh-1 ; ++y ) {
                    pixel_t* src  = getLevel(o, s) + xo + yo*y ;        
                    pixel_t* end  = src + ow - 1 ;
                    pixel_t* grad = 2 * (xo + yo*y + (s - smin -1)*so) + temp ;
                    while(src != end) {
                        VL::float_t Gx = 0.5 * ( *(src+xo) - *(src-xo) ) ;
                        VL::float_t Gy = 0.5 * ( *(src+yo) - *(src-yo) ) ;
                        VL::float_t m = fast_sqrt( Gx*Gx + Gy*Gy ) ;
                        VL::float_t t = fast_mod_2pi( fast_atan2(Gy, Gx) + VL::float_t(2*M_PI) );
                        *grad++ = pixel_t( m ) ;
                        *grad++ = pixel_t( t ) ;
                        ++src ;
                    }
                }
            }
        }

        tempIsGrad = true ;
        tempOctave = o ;
    }

    /** @brief Compute the orientation(s) of a keypoint
    **
    ** The function computes the orientation of the specified keypoint.
    ** The function returns up to four different orientations, obtained
    ** as strong peaks of the histogram of gradient orientations (a
    ** keypoint can theoretically generate more than four orientations,
    ** but this is very unlikely).
    **
    ** @remark The function needs to compute the gradient modululs and
    ** orientation of the Gaussian scale space octave to which the
    ** keypoint belongs. The result is cached, but discarded if different
    ** octaves are visited. Thereofre it is much quicker to evaluate the
    ** keypoints in their natural octave order.
    **
    ** The keypoint must lie within the scale space. In particular, the
    ** scale index is supposed to be in the range @c smin+1 and @c smax-1
    ** (this is from the SIFT detector). If this is not the case, the
    ** computation is silently aborted and no orientations are returned.
    **
    ** @param angles buffers to store the resulting angles.
    ** @param keypoint keypoint to process.
    ** @return number of orientations found.
    **/
    int
        Sift::computeKeypointOrientations(VL::float_t angles [4], Keypoint keypoint)
    {
        int const   nbins = 36 ;
        VL::float_t const winFactor = 1.5 ;
        VL::float_t hist [nbins] ;

        // octave
        int o = keypoint.o ;
        VL::float_t xperiod = getOctaveSamplingPeriod(o) ;

        // offsets to move in the Gaussian scale space octave
        const int ow = getOctaveWidth(o) ;
        const int oh = getOctaveHeight(o) ;
        const int xo = 2 ;
        const int yo = xo * ow ;
        const int so = yo * oh ;

        // keypoint fractional geometry
        VL::float_t x     = keypoint.x / xperiod ;
        VL::float_t y     = keypoint.y / xperiod ;
        VL::float_t sigma = keypoint.sigma / xperiod ;

        // shall we use keypoints.ix,iy,is here?
        int xi = ((int) (x+0.5)) ; 
        int yi = ((int) (y+0.5)) ;
        int si = keypoint.is ;

        VL::float_t const sigmaw = winFactor * sigma ;
        int W = (int) floor(3.0 * sigmaw) ;

        // skip the keypoint if it is out of bounds
        if(o  < omin   ||
            o  >=omin+O ||
            xi < 0      || 
            xi > ow-1   || 
            yi < 0      || 
            yi > oh-1   || 
            si < smin+1 || 
            si > smax-2 ) {
                std::cerr<<"!"<<std::endl ;
                return 0 ;
        }

        // make sure that the gradient buffer is filled with octave o
        prepareGrad(o) ;

        // clear the SIFT histogram
        std::fill(hist, hist + nbins, 0) ;

        // fill the SIFT histogram
        pixel_t* pt = temp + xi * xo + yi * yo + (si - smin -1) * so ;

#undef at
#define at(dx,dy) (*(pt + (dx)*xo + (dy)*yo))

        for(int ys = std::max(-W, 1-yi) ; ys <= std::min(+W, oh -2 -yi) ; ++ys) {
            for(int xs = std::max(-W, 1-xi) ; xs <= std::min(+W, ow -2 -xi) ; ++xs) {

                VL::float_t dx = xi + xs - x;
                VL::float_t dy = yi + ys - y;
                VL::float_t r2 = dx*dx + dy*dy ;

                // limit to a circular window
                if(r2 >= W*W+0.5) continue ;

                VL::float_t wgt = VL::fast_expn( r2 / (2*sigmaw*sigmaw) ) ;
                VL::float_t mod = *(pt + xs*xo + ys*yo) ;
                VL::float_t ang = *(pt + xs*xo + ys*yo + 1) ;

                //      int bin = (int) floor( nbins * ang / (2*M_PI) ) ;
                int bin = (int) floor( nbins * ang / (2*M_PI) ) ;
                hist[bin] += mod * wgt ;        
            }
        }

        // smooth the histogram
#if defined VL_LOWE_STRICT
        // Lowe's version apparently has a little issue with orientations
        // around + or - pi, which we reproduce here for compatibility
        for (int iter = 0; iter < 6; iter++) {
            VL::float_t prev  = hist[nbins/2] ;
            for (int i = nbins/2-1; i >= -nbins/2 ; --i) {
                int const j  = (i     + nbins) % nbins ;
                int const jp = (i - 1 + nbins) % nbins ;
                VL::float_t newh = (prev + hist[j] + hist[jp]) / 3.0;
                prev = hist[j] ;
                hist[j] = newh ;
            }
        }
#else
        // this is slightly more correct
        for (int iter = 0; iter < 6; iter++) {
            VL::float_t prev  = hist[nbins-1] ;
            VL::float_t first = hist[0] ;
            int i ;
            for (i = 0; i < nbins - 1; i++) {
                VL::float_t newh = (prev + hist[i] + hist[(i+1) % nbins]) / 3.0;
                prev = hist[i] ;
                hist[i] = newh ;
            }
            hist[i] = (prev + hist[i] + first)/3.0 ;
        }
#endif

        // find the histogram maximum
        VL::float_t maxh = * std::max_element(hist, hist + nbins) ;

        // find peaks within 80% from max
        int nangles = 0 ;
        for(int i = 0 ; i < nbins ; ++i) {
            VL::float_t h0 = hist [i] ;
            VL::float_t hm = hist [(i-1+nbins) % nbins] ;
            VL::float_t hp = hist [(i+1+nbins) % nbins] ;

            // is this a peak?
            if( h0 > 0.8*maxh && h0 > hm && h0 > hp ) {

                // quadratic interpolation
                //      VL::float_t di = -0.5 * (hp - hm) / (hp+hm-2*h0) ; 
                VL::float_t di = -0.5 * (hp - hm) / (hp+hm-2*h0) ; 
                VL::float_t th = 2*M_PI * (i+di+0.5) / nbins ;      
                angles [ nangles++ ] = th ;
                if( nangles == 4 )
                    goto enough_angles ;
            }
        }
enough_angles:
        return nangles ;
    }

    // ===================================================================
    //                                         computeKeypointDescriptor()
    // -------------------------------------------------------------------

    namespace Detail {

        /** Normalizes in norm L_2 a descriptor. */
        void
            normalize_histogram(VL::float_t* L_begin, VL::float_t* L_end)
        {
            VL::float_t* L_iter ;
            VL::float_t norm = 0.0 ;

            for(L_iter = L_begin; L_iter != L_end ; ++L_iter)
                norm += (*L_iter) * (*L_iter) ;

            norm = fast_sqrt(norm) ;

            for(L_iter = L_begin; L_iter != L_end ; ++L_iter)
                *L_iter /= (norm + std::numeric_limits<VL::float_t>::epsilon() ) ;
        }

    }

    /** @brief SIFT descriptor
    **
    ** The function computes the descriptor of the keypoint @a keypoint.
    ** The function fills the buffer @a descr_pt which must be large
    ** enough. The funciton uses @a angle0 as rotation of the keypoint.
    ** By calling the function multiple times, different orientations can
    ** be evaluated.
    **
    ** @remark The function needs to compute the gradient modululs and
    ** orientation of the Gaussian scale space octave to which the
    ** keypoint belongs. The result is cached, but discarded if different
    ** octaves are visited. Thereofre it is much quicker to evaluate the
    ** keypoints in their natural octave order.
    **
    ** The function silently abort the computations of keypoints without
    ** the scale space boundaries. See also siftComputeOrientations().
    **/
    void
        Sift::computeKeypointDescriptor
        (VL::float_t* descr_pt,
        Keypoint keypoint, 
        VL::float_t angle0)
    {

        /* The SIFT descriptor is a  three dimensional histogram of the position
        * and orientation of the gradient.  There are NBP bins for each spatial
        * dimesions and NBO  bins for the orientation dimesion,  for a total of
        * NBP x NBP x NBO bins.
        *
        * The support  of each  spatial bin  has an extension  of SBP  = 3sigma
        * pixels, where sigma is the scale  of the keypoint.  Thus all the bins
        * together have a  support SBP x NBP pixels wide  . Since weighting and
        * interpolation of  pixel is used, another  half bin is  needed at both
        * ends of  the extension. Therefore, we  need a square window  of SBP x
        * (NBP + 1) pixels. Finally, since the patch can be arbitrarly rotated,
        * we need to consider  a window 2W += sqrt(2) x SBP  x (NBP + 1) pixels
        * wide.
        */      

        // octave
        int o = keypoint.o ;
        VL::float_t xperiod = getOctaveSamplingPeriod(o) ;

        // offsets to move in Gaussian scale space octave
        const int ow = getOctaveWidth(o) ;
        const int oh = getOctaveHeight(o) ;
        const int xo = 2 ;
        const int yo = xo * ow ;
        const int so = yo * oh ;

        // keypoint fractional geometry
        VL::float_t x     = keypoint.x / xperiod;
        VL::float_t y     = keypoint.y / xperiod ;
        VL::float_t sigma = keypoint.sigma / xperiod ;

        VL::float_t st0   = sinf( angle0 ) ;
        VL::float_t ct0   = cosf( angle0 ) ;

        // shall we use keypoints.ix,iy,is here?
        int xi = ((int) (x+0.5)) ; 
        int yi = ((int) (y+0.5)) ;
        int si = keypoint.is ;

        // const VL::float_t magnif = 3.0f ;
        const int NBO = 8 ;
        const int NBP = 4 ;
        const VL::float_t SBP = magnif * sigma ;
        const int   W = (int) floor (sqrt(2.0) * SBP * (NBP + 1) / 2.0 + 0.5) ;

        /* Offsets to move in the descriptor. */
        /* Use Lowe's convention. */
        const int binto = 1 ;
        const int binyo = NBO * NBP ;
        const int binxo = NBO ;
        // const int bino  = NBO * NBP * NBP ;

        int bin ;

        // check bounds
        if(o  < omin   ||
            o  >=omin+O ||
            xi < 0      || 
            xi > ow-1   || 
            yi < 0      || 
            yi > oh-1   ||
            si < smin+1 ||
            si > smax-2 )
            return ;

        // make sure gradient buffer is up-to-date
        prepareGrad(o) ;

        std::fill( descr_pt, descr_pt + NBO*NBP*NBP, 0 ) ;

        /* Center the scale space and the descriptor on the current keypoint. 
        * Note that dpt is pointing to the bin of center (SBP/2,SBP/2,0).
        */
        pixel_t const * pt = temp + xi*xo + yi*yo + (si - smin - 1)*so ;
        VL::float_t *  dpt = descr_pt + (NBP/2) * binyo + (NBP/2) * binxo ;

#define atd(dbinx,dbiny,dbint) *(dpt + (dbint)*binto + (dbiny)*binyo + (dbinx)*binxo)

        /*
        * Process pixels in the intersection of the image rectangle
        * (1,1)-(M-1,N-1) and the keypoint bounding box.
        */
        for(int dyi = std::max(-W, 1-yi) ; dyi <= std::min(+W, oh-2-yi) ; ++dyi) {
            for(int dxi = std::max(-W, 1-xi) ; dxi <= std::min(+W, ow-2-xi) ; ++dxi) {

                // retrieve 
                VL::float_t mod   = *( pt + dxi*xo + dyi*yo + 0 ) ;
                VL::float_t angle = *( pt + dxi*xo + dyi*yo + 1 ) ;
                VL::float_t theta = fast_mod_2pi(-angle + angle0) ; // lowe compatible ?

                // fractional displacement
                VL::float_t dx = xi + dxi - x;
                VL::float_t dy = yi + dyi - y;

                // get the displacement normalized w.r.t. the keypoint
                // orientation and extension.
                VL::float_t nx = ( ct0 * dx + st0 * dy) / SBP ;
                VL::float_t ny = (-st0 * dx + ct0 * dy) / SBP ; 
                VL::float_t nt = NBO * theta / (2*M_PI) ;

                // Get the gaussian weight of the sample. The gaussian window
                // has a standard deviation equal to NBP/2. Note that dx and dy
                // are in the normalized frame, so that -NBP/2 <= dx <= NBP/2.
                VL::float_t const wsigma = NBP/2 ;
                VL::float_t win = VL::fast_expn((nx*nx + ny*ny)/(2.0 * wsigma * wsigma)) ;

                // The sample will be distributed in 8 adjacent bins.
                // We start from the ``lower-left'' bin.
                int binx = fast_floor( nx - 0.5 ) ;
                int biny = fast_floor( ny - 0.5 ) ;
                int bint = fast_floor( nt ) ;
                VL::float_t rbinx = nx - (binx+0.5) ;
                VL::float_t rbiny = ny - (biny+0.5) ;
                VL::float_t rbint = nt - bint ;
                int dbinx ;
                int dbiny ;
                int dbint ;

                // Distribute the current sample into the 8 adjacent bins
                for(dbinx = 0 ; dbinx < 2 ; ++dbinx) {
                    for(dbiny = 0 ; dbiny < 2 ; ++dbiny) {
                        for(dbint = 0 ; dbint < 2 ; ++dbint) {

                            if( binx+dbinx >= -(NBP/2) &&
                                binx+dbinx <   (NBP/2) &&
                                biny+dbiny >= -(NBP/2) &&
                                biny+dbiny <   (NBP/2) ) {
                                    VL::float_t weight = win 
                                        * mod 
                                        * fast_abs (1 - dbinx - rbinx)
                                        * fast_abs (1 - dbiny - rbiny)
                                        * fast_abs (1 - dbint - rbint) ;

                                    atd(binx+dbinx, biny+dbiny, (bint+dbint) % NBO) += weight ;
                            }
                        }            
                    }
                }
            }  
        }

        /* Standard SIFT descriptors are normalized, truncated and normalized again */
        if( normalizeDescriptor ) {

            /* Normalize the histogram to L2 unit length. */        
            Detail::normalize_histogram(descr_pt, descr_pt + NBO*NBP*NBP) ;

            /* Truncate at 0.2. */
            for(bin = 0; bin < NBO*NBP*NBP ; ++bin) {
                if (descr_pt[bin] > 0.2) descr_pt[bin] = 0.2;
            }

            /* Normalize again. */
            Detail::normalize_histogram(descr_pt, descr_pt + NBO*NBP*NBP) ;
        }

    }

    // namespace VL
}

// -------------------------------------------------------------------
//                                                                main
// -------------------------------------------------------------------
int
    main0(int argc, char** argv)
{
    int    first          = -1 ;
    int    octaves        = -1 ;
    int    levels         = 3 ;
    float  threshold      = 0.04f / levels / 2.0f ;
    float  edgeThreshold  = 10.0f;
    float  magnif         = 3.0 ;
    int    nodescr        = 0 ;
    int    noorient       = 0 ;
    int    stableorder    = 0 ;
    int    savegss        = 0 ;
    int    verbose        = 0 ;
    int    binary         = 0 ;
    int    unnormalized   = 0 ;
    int    fp             = 0 ;
    string outputFilenamePrefix ;
    string outputFilename ;
    string keypointsFilename ;

    // -----------------------------------------------------------------
    //                                            Loop over input images
    // -----------------------------------------------------------------      
    while( argc > 0 ) {

        string name("iss.pgm") ;

        try {
            VL::PgmBuffer buffer ;

            outputFilename = "key.key";


            // ---------------------------------------------------------------
            //                                                  Load PGM image
            // ---------------------------------------------------------------    
            verbose && cout
                << "siftpp: lodaing PGM image '" << name << "' ..."
                << flush;

            ifstream in(name.c_str(), ios::binary) ; 
            extractPgm(in, buffer) ;


            verbose && cout 
                << " read "
                << buffer.width  <<" x "
                << buffer.height <<" pixels" 
                << endl ;

            // ---------------------------------------------------------------
            //                                            Gaussian scale space
            // ---------------------------------------------------------------    
            verbose && cout 
                << "siftpp: computing Gaussian scale space" 
                << endl ;

            int         O      = octaves ;    
            int const   S      = levels ;
            int const   omin   = first ;
            float const sigman = .5 ;
            float const sigma0 = 1.6 * powf(2.0f, 1.0f / S) ;

            // optionally autoselect the number number of octaves
            // we downsample up to 8x8 patches
            if(O < 1) {
                O = max((int)(floor(log2(std::min(buffer.width,buffer.height))) - omin -3), 1) ;
            }

            verbose && cout
                << "siftpp:   number of octaves     : " << O << endl 
                << "siftpp:   first octave          : " << omin << endl 
                << "siftpp:   levels per octave     : " << S 
                << endl ;

            // initialize scalespace
            VL::Sift sift(buffer.data, buffer.width, buffer.height, 
                sigman, sigma0,
                O, S,
                omin, -1, S+1) ;


            // -------------------------------------------------------------
            //                                             Run SIFT detector
            // -------------------------------------------------------------    


            verbose && cout 
                << "siftpp: running detector  "<< endl
                << "siftpp:   threshold             : " << threshold << endl
                << "siftpp:   edge-threshold        : " << edgeThreshold
                << endl ;

            sift.detectKeypoints(threshold, edgeThreshold) ;

            verbose && cout 
                << "siftpp: detector completed with " 
                << sift.keypointsEnd() - sift.keypointsBegin() 
                << " keypoints" 
                << endl ;


            // -------------------------------------------------------------
            //                  Run SIFT orientation detector and descriptor
            // -------------------------------------------------------------    

            /* set descriptor options */
            sift.setNormalizeDescriptor( ! unnormalized ) ;
            sift.setMagnification( magnif ) ;

            if( verbose ) {
                cout << "siftpp: " ;
                if( ! noorient &   nodescr) cout << "computing keypoint orientations" ;
                if(   noorient & ! nodescr) cout << "computing keypoint descriptors" ;
                if( ! noorient & ! nodescr) cout << "computing orientations and descriptors" ;
                if(   noorient &   nodescr) cout << "finalizing" ; 
                cout << endl ;
            }

            {            
                // open output file
                ofstream out(outputFilename.c_str(), ios::binary) ;

                if( ! out.good() ) 
                    VL_THROW("Could not open output file '"
                    << outputFilename
                    << "'.") ;

                verbose && cout
                    << "siftpp:   write keypoints to    : '" << outputFilename << "'"         << endl
                    << "siftpp:   floating point descr. : "  << (fp           ? "yes" : "no") << endl
                    << "siftpp:   binary descr.         : "  << (binary       ? "yes" : "no") << endl
                    << "siftpp:   unnormalized descr.   : "  << (unnormalized ? "yes" : "no") << endl
                    << "siftpp:   descr. magnif.        : "  << setprecision(3) << magnif
                    << endl ;

                out.flags(ios::fixed) ;


                // -------------------------------------------------------------
                //            Run detector, compute orientations and descriptors
                // -------------------------------------------------------------
                for( VL::Sift::KeypointsConstIter iter = sift.keypointsBegin() ;
                    iter != sift.keypointsEnd() ; ++iter ) {

                        // detect orientations
                        VL::float_t angles [4] ;
                        int nangles ;
                        if( ! noorient ) {
                            nangles = sift.computeKeypointOrientations(angles, *iter) ;
                        } else {
                            nangles = 1;
                            angles[0] = VL::float_t(0) ;
                        }

                        // compute descriptors
                        for(int a = 0 ; a < nangles ; ++a) {

                            out << setprecision(2) << iter->x << ' '
                                << setprecision(2) << iter->y << ' '
                                << setprecision(2) << iter->sigma << ' ' 
                                << setprecision(3) << angles[a] ;

                            /* compute descriptor */
                            VL::float_t descr_pt [128] ;
                            sift.computeKeypointDescriptor(descr_pt, *iter, angles[a]) ;

                            /* save descriptor to to appropriate file */          
                            if( ! nodescr ) {
                                insertDescriptor(out, descr_pt, false, fp) ;
                            }
                            /* next line */
                            out << endl ;
                        } // next angle
                } // next keypoint


                out.close() ;
                verbose && cout 
                    << "siftpp: job completed"<<endl ;
            }

            argc-- ;
            argv++ ;
            outputFilename = string("") ;
        }
        catch(VL::Exception &e) {
            cerr<<endl<<"Error processing '"<<name<<"': "<<e.msg<<endl ;
            return 1 ;
        }    
    } // next image

    return 0 ;
}
//=================================sift end====================================



struct siftPoint_save{
    double x, y;
    double sigma, angle;
    char vec[193];
};
siftPoint_save sps1[] = {
{28.83,13.96,3.46,1.04,"############WI-######'$##OG######1)#Oh>######1)#ZP##########4GM######f##|.[######1)#NsF######%'#QG#<5#######2>H.,####'##p4[(c$###J##CH9r@/###1#####$##1,#95#l-&.,#$##95#`j+;n.######.r#b.[######"},
{146.12,98.58,4.09,6.22,"Q*-Kh=.[&O-(]I'.HM{-&0/-sy)SaF$##,-%>e$#:7Nn./,#0y*?:6OxEtL:7KTh07i6(RT*Bn*JI+'7#|5JiC4t[-n,$(i<2/+37@JOC&l#/KT{#%~Q'ZQ$w|F%##*6#CV;n>%~P#D,#N2?@R)*$%o?O.,#^KTkY#M5$###D$Nk>$RG#;5#.##2-']P#E>#"},
{35.43,104.30,3.82,6.01,".,##########p(?######g##GDdZP####a%#9BHkf3bc'1##D>##########:wY######:##(CdZP#%##3d#i16d[,5,#ew%95##########{[W######3##KDd######<6#APM###*##=U%############~19#########^Fd######@##/Cd###$##A~#"},
{109.72,103.63,4.05,1.29,"s7-3&HV,%-c#aM>|J+JR)8><aP#)S)$t<kg8'##+e-?m&mT7g_=;C5-v#)/O~-O{@DR%.US&?>H~X91A.cJ-tc(###wP#5-Opb#M,#Zu:74FvZ(a7&77BFw*e,O$6%xQ%%o(qx5$##:5#q6''##4?B4p4###]c&uh01x,P,$J,O/H%Du#{?$6f3%##:c$Rv&"},
{109.72,103.63,4.05,5.44,"OD>[P#LI&jd&3X;u1=###S##|c#AUN/U7.,#1ECBh&wU;qZ#C-(-H%4M'Tn.2]+[cInx)qQ(yQI;T0/+6y/07*EY,#Jj0339%##9,#x''pm,qc'g,$CW,K{?vIOG~,<B+s14^K1[f*gmJpv)Ql$'Q#dl&######l/#9z86#$:%'or6-%,wY#o$)'0/_Q(%J$"},
{549.98,111.96,3.57,3.10,"Fc#+7)pv#`uQhG$iH%J~$yuQKe.###?##_zQVV?###$##^*/zI'YL5Fc$G6&K_<*J):6%)E2S#Dh,$JT.)q3sSQ1,#gY#h-$TZ&UG#p.$|wQE;:*7)`n$G#J$2:q]*+^3i_5bvQA5#OG#x@$######w##[ED$##>5#,/#q+M4l$D5#:?#3xQXv*.,####>x'"},
{113.83,114.53,4.28,1.34,"9]2z8*~7,q*:uY#''2L=;ZV<'##lD@a-'rC=i##oC=F##kp:>~ST?>s$+VT)U}D=Q;K17v'4c$+###u5$G~S###'##G[#lZSUv)9/(N?@6&-i[S07(h@(pJ*4*E$##~P#8S*###88#(S-Pv)vP$>D/C'0hG$.[S=v&[c$GR%^iC%##o,%Kw(###@%#'nM~Z'"},
{576.01,126.10,3.56,4.85,",`>95#######Tf~M>#7H$??%b[).n*s00?Z%B3;%%+3c$p>#'FG###0,#6,#Ef~iG#m,%z>##&0=n$BE:l#%RXDJc$3-&F&&VEE.,####a,#+f~'Q#u#&DJ&*J.KA%2U30[KqB5J#$~P#U.AQ%/######8Q7W%[###+##DF+,R*###X5#EKP0,####)##^j?"},
{691.74,141.39,3.92,6.13,"pD$SEE7,#@u$%>/^Q(######+R&2u#.,####a#$D>#######ZF31Z%Z1+W_:>>~xY$&##,B'VN<p-+###.##YH$-H&######i,&E>#3N%G:~;}I@H%aR+/)Nqy,_2='%,sG#WI#OaD######5##pP$^Z#Yn/###6,#$g07~-IH$7-&+D0Pn/,?#%K/NI%Rd+"},
{602.90,148.67,4.25,5.96,"Bj&{1?=5#=5#g6]tc(1,#-Q#xq31Z%###+##Kc#;c%#######z6Cu$mY#4##}ahD>####*%#jmSoG$ub#?$$#v#Hl%######Qh>)##@-'p5$^ah95####s##Q7OAH&[P#B,##6%.u####5,#]./%##8-${S.(`h######I.$r7T.,####<##QZ%ZP####(##"},
{688.89,149.78,3.89,0.81,"X##AsDE>####%$#x&2|k####$##s10)R*######uy+>o2|>%[1(eBY&##O,$BCJ9+G-6&}P#g$)N13+C/8A/(##2[(:R)CE5*K/`82N##u*XuBYmP$8,#{}1kn0###$I#9V7###;5#K5##q3ZG#cl&'##NFYxY$rl%=#$1Z:RG#m5%;5#h5$&##iY#[P#95#"},
{827.73,158.21,3.84,5.87,"^e(8H%YUZL5$(YZlG$~-*###(M3BZ%ZP####1,#4Q#@H'###$%+ZT&QSZ)##EUZD,#%p6`##0nV0,#TG#4H#F>#1##v6*;5#pm*Kn%7}?ai?ISZ(##6.)>.DGRMD>#C5#z5$_>#1u#y-):5#qk1c#&H##kU9-(*2Z%cZ#DPHR/)EH&I-&F?'X6##@)a5%###"},
{108.99,169.21,3.30,1.60,"j|3J?$mU4pb#LtZ;5#######XA,4l#.,#######Ic#95####]s@1i)|99b##/sZ~P####|##'6K7u#.,#%##%##nG#I#%###f(>Bv#;.,hs+InZ###%##k*+G#ORG####L##=5#K#$bG$/,#L5$+##`d'dNTxd-,##[R(Nb:@%.G5#ZP#%-#M,$T,$D>#1,#"},
{837.75,181.54,3.69,4.88,"qb#4##m/K2l#dlR0,#+^.Wg%)]3kG#+.?E13ZP#Ju#LR&2LL_l&)W,#*CK>#'mRpc$m7.Hn#~D@`l$Q{89e(8~/%##Zl$k:,yS/~mH+u#2-&-qR&%+OG#HH#8a=.K1/H&R##'z9W5$PG#O##76'v5&###+oRyB5|k####_P9<^6E$)###4n#:R)U'6###0##"},
{816.37,186.28,3.80,1.25,"[13kvS###tY#B*A`B(jH)Uc#Y>$pL#b(=Sm&:5#=Q#]e+(bEX+Ea~0xc#-}?g+GBk2`m+W,#ZvSkC'($(bS$3w-$##>##<a8mG$&v%WE.G[*e@-Sp(tb<f>$nwSk#$,$&l#$nq=###.##Sc$*##rN;%%)###3c$U:&EL4F>#gxSju#$$'Q>#m4H###-##A#$"},
{816.37,186.28,3.80,5.18,"5T2qQ%s%)Pl%XdMoJ3%##t>#cn)7uN;#$sG#=c#/ABnZ(S>#o#&D>#~($A)Aqp1sg:@M%X@-/JTeA2R8(`J*t5&`n&)>5'<8###E#$Mp%`$+r?*'##*Q707+'JT$##&0)5%)3C8$##d1(*:4$##t?&Gu$###F>#.B&kOF.,#/KT,[%d82G>#B)9,.*@m(2,#"},
{554.68,193.38,3.82,1.83,"uU(PQ'7#$ub#T{50Z%'##'9,v(;;-(@5#fd>RZ%80/Pm&$&Tz*CM>#J#%H##Z%TmY#nZ&P<,-{7)`;MHLpC/u['8mI'?&dP#h:=$##.,#s.#O)TJT2AQ&Sv#/{2J.@K$T%c#(j?>-&o#'1,#W^895#$##kd$Ze(LR+;##h7,-o/8(4,c#ke,Ir@_>$###>H#"},
{538.22,197.70,3.63,1.92,"s~R2c$)Q$Eo&}'38(2?~R%/*]'/Tz7`./.,##Q$su%]P#~P#]_R}w0PG#6$$=D-]_RW)BcP#K~Rr[&G[+8##:I*O,$I>#2Z$]e*AH'.##~~-ov)T&/`5#<<=X?QVl%%##AT*yR*#n,/,#NZ$yI)j>%4##]~R(d(###'##y^R6Q%L5$=##b[R?c$@H'$##j/P"},
{538.22,197.70,3.63,6.02,"###0(#fuQ######V.#LuJxm+.,#B6#u7+lS/###&##q;16S/###0(#g{AnG$B6(Qe#OPCMc%+RO|#$8?%pQ<A/-,?&R2)8xQpb#7$#.K1wn)k>%%##&i('C8</Fg,&2J#0j@Z433<A-9.Vp0aH(:l#6@+j#$oP$iP#kZ%'v'56'<5#3-#.U6&=EtG$bu$%2-"},
{473.20,215.25,3.89,1.66,"'~.>5#=5#b5$iA_$##]e/&@#(K40##LD_j$%WA2*##aK,RB*lK7I>#95#+##mA_4##sJ3i##{99]H#{B_Cd%M_=Ac#wA13w%0q:'l#D>#&##zA_V5#a~/-6$:95]-#bfU(a<'o1N>#Sx-QAMTg9A5#E>#'##sA_.,#<l#rQ%*x1Hc%8I%S>H.,#[7+2,#]N>"},
{473.20,215.25,3.89,4.85,"%##$i9###Sn+Ov$,bD;p73Q%zY#i#%AK_###ZP#$###94C5#7/,VfN~~0@5#~fT'a<-g6F-#K~.~5$NK_D,#OG####Pg72u#{A1A.&]V><c#KL_?d%Ig8/?#OT5(##)K_.##95#$##/C9I>#No)Yo)506,##6N_}6%~83,##cf3L?#vJ_###=5#4Z$3&2=5#"},
{745.68,227.41,3.62,5.99,"D>####c-N$##`/4###HTQ9I%V6U###I##L()v$,###}L#so195####UW:3,#kT7###O1NJ?%v6U###Kl#YJ$IC;###+%#+:01,#jY#O[(0,#5v(###g[>W5$=;U####f)`#$ZQF###%##|u#K,$.,#3,#XG#;5####H:.ub#wq0###ji.95#fl8###M,####"},
{652.40,235.30,4.27,4.55,"r@#fU^^>#hv+E2,T7..##BU7q>%###+##{q8######d##|J2j`02xRh%*YR+AZ^k-+*##C$%Si@$##Bu#07*###'##^,$FH'XL5Ml%k.&A33>T^###D,#(O(_M@H>#%l#i6%###%##+l#z#'8936#$<,#0gDDe/###%##]N&xY$###4##w$)######*##kc'"},
{471.67,246.30,4.20,4.80,"r[-######5##K%ZTG####kT$`o4Ig4.,#Y_/wP#'X?eY#?R)Lx3###%##'Z#L%Z0Z#M$*J$#&2<vl?$?QfG#bU5TY@<%-T,#Py5Nc&###B,#s%Zi,#je1j##ef3`S$X%ZnG#l<Cjd%HT4:Z#yT.TB7######P(ZV-%WT6%##WK-^M+E%Z$##S4@ox)(@+'##"},
{498.34,250.33,3.98,2.02,"0Z%(-#j#&DC2uR..##Ld%eiU95####If&KSP######,##&0/pm,$##a[%(fUAdU###f8'/Y4Mm*###pQ:y&3######iZ#~N@~./###oC0O'2EeU###Y@$wo/Xr@###g&Frl&######w9,ax3EJ0###VM0[d'SOF'##q-%c.+s3D###7F=Q>####$##va?pb#"},
{498.34,250.33,3.98,4.53,"h-$hv+######c6RC,$.,#e$#4iAT~#[J2o$#J?'sR$(/1Al#C'+yS3######/xUc?&_5%.$#fD>*<)}_AS##<:3/1-=6(.##5W7s/5###+##NyU3=7g,&Q##-^.ik3,o2$##ki;4K-ZP#y##@/&K}K###'##ah*XxUOG#=u#Lp%7wL6#$(Z$n43C07###a##"},
{818.67,252.48,3.55,1.53,")f&hY#j2.qm,G=/9#$?Z%95#Ej,lY#ZP####$c#eG#a5%###P>#,f#DTV|k#[6;)H$(x.Bc$jWV.,#$##r5#|v+6,#h>$E5#0$(C:$I5MY['i_?zY#g,$^*6tRVrb#%##1K%E[*LH%B,$0##{k#4##/q-#TV$7+###PS'>nL}.TOG#%-%--%p6)(H%D>####"},
{502.18,259.80,3.55,4.48,"/$#e]6######$f&I..######qRU#l#.,#q$#[p9Rn#m[-L$#<$#LK6######*k:Nn/###.##.SUGm%>u$[$#=K3c2(rK8B##4$#,*E######EoLh&5###B##gTUGZ<I#%f##G~)[F0Jm*###z##=bL######gh'hGP###*##41(@TUD>#6u#TT%iFBOG#ub#"},
{522.35,275.33,4.15,3.05,"###W5#j>%###A|D=5#.,####apb######)##b*G######Q##2,#_G#LQ&gY#)>K~P#PG#$##Mpb######2##o#Q######<##0,#hP#zZ&-d)(tI.,#wP#{>%8pb######_R%EjG######Lx(######Y@#*WBGh=###C%$PL5cqbeY#4,##K>'&-W>$%##.tb"},
{181.19,319.17,4.06,2.09,"###w5#eu&###E,#,M15?'###LR&Sf146'1,#:5#`,#@?'RG####m##F'7###s^3|S>aTa8,#-ZaXj=Dx/Iv$@$)<,#Zf/cu$###T,#%;<###ey8Ed#=[aeE<{Uah$'Bp1J68Bm(0[%M~/4,####L5#kH(###D>#*##dC-^g7oZ(Cu#EH%Sh4]P#~,$/u#RG#"},
{365.66,338.51,3.46,5.31,"?Q#_E9~-*/,#iE.z0895#'##Dn+######C##Vl%F>#.,#.##;n(^N3LU^3Z$AZ^{05UJ1d5#FdO######:##=Q&######4##Lv(Ru#1Y^a93IT^'###C19798#O95####<@#0u#'l#95#}Y#=#$K>#GK+A^6iZ(###u>#uQ<nc'{k#&##jc#(l#N#%$##4,#"},
{575.94,360.26,3.69,1.29,"DA+g|FRG#F>#/GcK?(###Q6#TIS]P####L$#^P#t#$D>#(##|~2'l#nP#LQ%RDc######-7#fEckH)###/##L##z$)pb####N:<7l#QG#e6#LDcZP####U&#<~=M/3###+##L8#V%/95####7`BhY####@&#i7V95#####&#`L'P6)###$##{C,CZ&######"},
{200.80,400.93,4.27,3.08,"I>=M%(Zl&###,?$*|;Fl%###cy-bA11u#5c#/l#xH&LH'kY#]8WNH$iz9C$#/S(h:Wz:W(Q$c<WjPLml%sc$yw..H%~G#9c$ecR###O%*Jo#|Z(]u$wpK:D:t7WwP$1/&0.@]L9S>#&Z$&c#Fl%$##_G#RR&D>#3,#F0+;$(hQ($##A@&l7,v@.=5#Bu$eG#"},
{200.80,400.93,4.27,5.65,".Q%eG#q>%ju$E5#Rw&8.-I>#2f.Cd(L5$4##nd,PG#%##8?%_5$~Q$)7+eG$ov'wH7`PLX5$6=]BN;Gc%^#$Pf3(d'ZP#*##%l#Pl#*L3tv*tH)In$y<U-N?o8]yu$>~&?>9Fe/?,1fH)<##xb#2d%rc'S5$.,#1##Y%'A05M$*(%#Iq9yC3bG$*Y-<8].##"},
{269.51,421.97,3.91,2.77,"3##oq/dl&.,#E9+LK5D>####XLS/v(###&Z#HuC95#$##~v)(##Cc${o+|5&Dp4/u#9u#_%,)JSJv)###_?%K47X|E###)########QD'R#LoH)yb#`w&OMSL167`</W<3C6{o'_IS(c$###L-'oP${'&?mLDm(c?&BP4l*B>5#r%&qIS-c$<I&#v&^v*###"},
{385.56,438.27,4.22,6.21,"f-#t~PI6(OG#';'3j?I#%###rE/dQ&pb####:,#$r(=6(###.-#lT2{P$95#)N2GU9:5#;,##_PMv%3e.9##Ou$6($3ZP?5#0,#>5#E[(Au$xy9r>%]>$q8%KZPJ>#@T1j0#`5%?##(^PZ,$sb#.,#~Z$?uKBv)###2,#{AHDe/###`l$dL%D>####ij1V5$"},
{441.74,440.09,3.62,4.13,"l$#[)@q[-###;,#k@(W-R/H&;5#%$$?e*F7-_P#bP#Y-(8H&%U.Fr@[P#8,#]FEB[)#J)O29Eu$/,#`z%k-R######('%J-R*f0G>#D5#nw-$0RE>#pY#:y'e'6h5%E@&9gK3d(}[*>+2#1R9-(.,#,##q_9u-R95####r;*}<@Tx-Bu$a2.`l%G_5cd*Ze+"},
{441.74,440.09,3.62,5.96,"n7.{b#H>#(n'OUS_5%###qJ$]wBm969L0PL0^e(RZ%D],#v'uT6RG#nP#[l$WTSDQ&:5#4/$ha0h$ONu$6-%2K${TSG@*5l$^m+###,H#%5FFGF/2:h?(>A+<R%nVSB,$L>#<h,vGK#########-##^%-daC0,#rm'TRS]6)Tc%sJ.[-)iY#tf/Cm)###$##"},
{378.89,455.30,3.73,1.02,"2I#,TX;c%###CX,z:>N5#{,'SG4D>#8##9h;=c$H>#0,#tp8d5$dd&jl&iG$AKEfG$YG#9l#IXX######X##hn.T,$.,#bG####/,#a5#`..0bJ.,#>##HC)TSXb#&###p0#)~-x8-SG#U[(######}##vWD/B5}Y#~G#OP2=APz~/D>#Z~#s:6Nf-C-'`=>"},
{157.46,456.84,4.17,6.02,"*?&#########Hy`######$##@Fn##########/0.,####&##{G%#########bDj######(##_Fn######$##c9795####$##vG%#########=y`######?##yFn######+##Mq<.,####&##nP$######%##;bI######.##fFn#########Gq;95####$##"},
{375.09,465.54,3.56,1.23,"^D&0S/:5####XACI#%######s~-pY####8,#.,#(7#,$(2,#Y&/,u#<,#yu%S^VeY####+-#VtHGn(.,#m5$PG#EJ$*A+=/,(n-###*##0y.)~ViZ'.,#jS#,&M&'33,#`f*F.+J1)mp6$b;0Z%###?5#3>?r%0dS'.x1%'&k+2Z]VTH(~>#f~'v`V2@+2l#"},
{467.33,471.58,4.05,1.00,"`Q#v6OUR,###_Y3})C(c$###8u4g,&######j,$L,$.,####]5$EK/}6+0,#-oE3?&}5&Z>#W3WOG####M##Ke-sl%.,#$##>#$vP$56%7B2yV@sb#:c#`i/b.W??'###&z#6.,Po.w>$4K+OG####A%#r4Ktw0s>$Nc#yG7zIMW7,D>#wK&G/.er4(?%zF5"},
{464.31,482.40,3.66,1.23,"b{(z6+|k####tR?mP$######/~*p>$.,#######V%$k,&###a/.'?&k>$OQ%2hUM5$$##HH#a`?`w)95#q5#FZ&{[$Wx+YV6L@-###0##AL/'eUN#%###j/$E?J`R*hY#5M't[,b{(tx2di.Yl&###k5#K+B993~e)WR,@/'4?8#fUk>%LZ#,v$B{R_?):l#"},
{373.36,509.75,3.49,2.80,"f,%I7(%6&/?%+E<:$)###p,#QoTxY$###v,#g>B6#$###(m&Nc&%##he*wZ'73D###/,#WR&xnTb-*###Ol#PD.JbK######xY$###82*7,GgH)'##^H%EhO;L5$D68f0r80B[#~nT{k####f?(###M(#IpTE>#d##HY3MvM###~J$CeT(Q%E##jn,^Q(###"},
{702.32,514.17,4.31,1.57,"Z##P<AZP####{T2OI,95#w,#kI.'##Ef,#U0qb#on$W17;o.]8)>5F&m&S5$')Vv6+$##%6#*XE_.%do3S7(###*M#y$Vk>%?J/bQ'W-&=d%2%V.c$3,#a~#S{@]q,;-T4#####dU#U$V###Jd*V,%M,#x*92bJJI'wu'ET$sc'{m?o(>2,####&w%@3D/,#"},
{179.80,522.05,3.68,6.01,"VZ'######%##R^d######+##bNk#########DB5.,####$##fH)######0##C2k######D##)Ok######%##IM<######$##%R)######&##x'c######>##>Ok######-##~;<.,#######ju&#########Y,M#########8Ok#########lz<#########"},
{463.29,527.47,3.93,2.61,"%##^F5XQ&I#%Sx)'D<WG#-c#5oF3S/C#$4?#`q9]P#tI%8_:###=,#$r0.%-V~0###C?$wm*WwSOG#%##xc#$.Do#'###W#$######VO0GHOQu%$##+K&qvSGM8BJ/d-(3e+J)*H2A######+y/:#$}3/txS<5#G##ezPIp85##*8)SFGmP$(&(jw/{k####"},
{407.85,529.92,3.95,0.67,"$W)WV=ub####T#<s5&Yy&T%/@c%###OC#{6V95####)$#Y7Vhr:%[&#l#&##p:V.,#*##Vl$#q:E>#k##Fz:)l#eY#a##ZkL#K2G>#/H$0y'O7V6#$###q:#?L8~7-.,#~5#mP#yU995#<5#,e-?#$H,#0416aEV,%5,#,c5{5%ue.)/*BPF*##cE=B_6t5&"},
{657.98,535.86,4.07,3.59,"jZ?######v>${t^######*##4r<6#$###3##$%(VH(###$##dh9D>####bG#lr^######X##uOJgG$.?$Bc#j,#nI-pl%###JS0D>#$##9C*mo^######e'$1,L.,#}G#$6$0,#hY##w$Bn-o#'###`C%8s^~o^###-?#yt6sU;######&f%L5$###9##Q8,"},
{743.88,542.59,3.84,5.54,"n9Tbv)HH'S,#MM5-12&cM.%([p-,~,>Z%$~)?I*a>$95#<,#g;Txz>ZP#:m$l9T7i>cm)0)/OS/]#%sC-P)9D[*E>#_u$I#$P]26#$5##L9Ig;Tr[-^>#9L0,`.ni=uA07l#B#$r%-[c&k5$;5#C6'T,#9O=Kw*&6&i,#=i:^A,Qd*###y{4=5#HI*###2d?"},
{358.42,558.53,4.45,3.37,"|$#i)=.,####MC#@V>:5####8~#W~'n2C######0'&%4HE,$0##zR(RJ.OG#.a.~]2~#&(##r0M@?%M|?h5$7l$g>#4vCnA.&##&Q#-8D]83hy8|G#i~-c')XQQ###T$$Wh&1H&fY#kL$Xf2Q5#A1.Z%M4e-ZR)1Z%Ud'7R>:93CZ&'##_C%4Q$xY$4$#Sd*"},
{290.72,565.78,4.14,2.14,"I#%9##CS(kL6rkK+u#4##}q%n/VZP####B%#}$*[5$;5#1,#MQ%cy*x:3$^4},<W{@rG$UR$+2VSH(###F##ux4VG#TG#H5#g>#]P4X%/###S0VRJ/.,#7$#A.VeY####r&#R%.;5#/,#>##UH%g/0S#$IS-'.Vb5%###8)&^-V###:5#}'#%Z$:5#:#$<,#"},
{764.33,570.12,3.63,0.22,"Vn)BW=?u$c5#4SY>-'uG%/%#H`CT%#iRYA$#,d)b##L'83.'o~0;?%p#'#&&LSYAH$v6+J.#UD@p]#hRYq##p':I##T7.fQ#P80rm$;x2(Q#=WYK~(2I+=##*V9K;'D6TE##eU;=,#T,%y,#fe.;q)&f2T>#/TY9-&P%-D7$cp9Oc#15Jz#$h]6###G>#|G#"},
{421.73,575.70,3.77,6.08,"TJ$K8SWH(###Ih+,0.Mv)###{#:|Y$.,####4?$qJ+#l####J&&//-@H'###XW1j`BQ5$+##`:S:7)0Z%4##;v'Yi*X83###e5$Tu$>S/=u##i=`m+95#U]#]6Sr,$.y50L$=c%{A#27StG$OG#%##f$&bbHD7-######X578K4^P#>'3%M%.,#0,#Q:S%Z$"},
{441.97,577.07,4.13,3.19,"t##{:;gY####vM)bn0###/,#L`0AQ#0w+Ol%&##Kz$kC<eY#*##>[(5S,rP$%O4t%,u5&?5#-;Ruc$(L8=##T5${I#j6R8Z%###9v#e[I`#&{:=p[%XQ'`U'E6R;5#-&-rL$?u$&##/SB)w*cP#nc$IF@L;=F6(3,#yG$m8GRR,###(Q#M;(95####r^&F6("},
{653.52,580.18,3.95,3.30,"wdL#########EE3`P#+u####3w*kl#D/13o&J3F]##3o24_$AeB{k#######1AC###.,#1##0wR9#$A#$OV&{84kv%|y96:$#5:RG#.,#/,#xyR######e##yxR3v(###4$#Le(#V8b#&0##n.,######R$)>zR######aG#b]QA,$###+##%A'N$*###%##"},
{463.86,586.60,3.81,5.98,".,#R^$gXH@$'a..P$'#l#D^.mgSQu%&##i''XZ=1]28^)Yr6ZP#*##HBNXv(B%-U,%],#6J+{gS7d)$##.-$$*+pdS7Z$;5#######&n?:$)Tu%D>#yR#.bIc94pn.GR)]J/YH$b}AFl%#########E{+P6)######j/%i*H###tG$<$GUd+D5#Fm(,e,###"},
{315.75,587.86,3.73,1.79,"CF*3vS95####q^HVZ'######j.,######,##,u#>,#Vu%3,#/=?Qc%Gu$9[%jY^######y##c5M$##PG#T5#PG#0##L-)2,#DI,###]>#Iu1)T^######eE'GcNT5$95#F##&##q5$xG%###7:3(c$*##7N+tGM#######L%q&4'Z$$##F##^P#?c$###1,#"},
{507.39,592.75,4.31,6.17,"FR#OhW..+{k#RN)E_5IZ&###ZP36#$######uZ$T,%######k[%Y^4/H&###D}5PL3Fu$U,$DjWq5%.,#^5#.g3{k####%##_,$1$'27(G,$G<A6[(`P#q(+)eWG>####*{$Rp8######[>#_P#J>#=A&MeWtJ2F>#6,#XCL[p9###$##>M%=n-pb####0##"},
{650.21,594.51,3.75,2.67,"{MB######V/#t)D######,'#E;>6#$###B$#I7*`.*###%##XlQ######@%#wlQ?5#.,#5$#].-eB/+c$Gl#9l#r+>95#yo+RlQ######V7#IGM*##mY#}q+HZ&0R$(@(zpQ###h%&H,$XqQrVA######M1$Eh;A,$###HW+~Q&E/*7l$I&K###Z(&/d(,f-"},
{756.42,598.31,3.60,0.10,"?H&?D)IA1yl$8AWCw(6Z%FI#?r?=.%o)Bs6$Zf5U,#>6(6-#%g3IG3Ze0|b#TAWIm'IH'^A#U*E~n%Y@WQm#jM?oY#J?(X,#FN>0v&gY#|,%oEW6)<{k#sG#DX87VP_@.&##|a?aQ'eY####Nr>###%##c#$Vf+v-+.##r5%n<;bx295####MZGQ,$6#$###"},
{609.26,604.75,3.52,6.01,"'Q%iI%^f5%##3eY+-%OG#u##9-S%H$###K##V?'$I(###i##9@,'6#F=J/##7eYJ>#+u#O$#75LS6%###0$#CS0_I&###j'#mw0]##kNF*##3fY0,#ZP#x##y+ITm%###g'#fH)fN)###2/#nI..##mU;5##yeYXQ'/,#_%#.(:yY3###u%#0Z%-N&###6/#"},
{500.30,610.51,3.77,0.98,"/I$:$O|%0TG#>$8*D>Y>$/##Jl4ll'######9-$<$)###$##q5%P8,0w,B5#U/FkA2e5%Q,#G+Z&Q%###H##@/0@H%rb#L>#$##}Y#)v&&&1W2@VZ&{,%[|02&Zx>%###<D#%~-<K-.,#*##OG####a,#7&Z`$+###2##,yE_p9.H$D>#yB#sb#z$&`5%$##"},
{385.95,617.59,3.99,3.35,"{##NXG######`##qlQ######/p'Oy8###(##t['FZ&###B,#.$$J[Tpu$:m'L7()~T.m$~@,$_T0q;;,#R$%{B7######Ld$RG#2v'Qn$V~TM?'_I,V8$x~Tc[TL#%]##q7EK{>.,####1H#3##P^TR>#B.-]0&1[T0##3J/(:(HC;&##p[*Dp*k]6###2##"},
{385.95,617.59,3.99,5.01,"U5#*w'h,&###S-$J^7######h5#$*?pb####,##(]'@NB###*##L0)9U9###K-%OVTcg98,#(oGuRTUu%N##YI*`g&EWC%#####`##7HO###F@+Z@(xST]Z$BTTI?'Y7,j@$Ex1c7$wRT*#####P##AXD###(?%Tl#kRTr>$'TT4u#g?(=w$}83%##~SJ>u#"},
{217.63,637.53,4.22,5.99,"'m(######(##L~Z######O##qjq######+##L&2D>####3,#cH(######%##C1b######@##Ckq######.##Ag7F>#.,#A5#kZ(######(##Of]######5$#Mkq######I##Q(;.,####,##WZ'######&##@?P######*##wjq######+##Py8######9##"},
{575.46,638.45,4.06,2.96,"d&+PV6Yc&(l#%8.AZ$4e+u>%_K[.,#%##1,#,}E###>5#3,#j)^&m(c,$_6&q&5###){3;I)s&^###Fl$`G#BnV######*##T(^$##]u%:,#*)?M5#~=Db5$kwZL>#v#&*-$A~[###$##}>#[&^###Il%0##{(?^>#:*DS>#^kLYG#HZ&p?#EYM######R6#"},
{359.67,651.32,4.92,5.19,".?#&vBk>%$##%`*>q;###&##h6'1,#D>####gP#~G#,u####v7*ev:i@YyY#NFY|]58#$(H#]*C3,#'l#B#$@,#4v%FZ&.,#OR+<l#DFYah9v@Y]P#+H$T>3z)D;,#kc',v#$##hQ$;06.,#}P$###1y%C>Jk-+###>##twAT,%+##Q6(VZ$###*##CT1###"},
{679.67,651.93,3.59,1.76,"*A[######B##WGN(c$###?$#~jH+u####v;+&~.&##95#$p*}@[######K##+A[######:##8A[######b5#X~0A,#X>$*c#x@[######:{'#A[######S-$;A[######@##d~0'l#rb#/,#-H&######>BJA7-######YHCxbO######1')Om*nY#95#0##"},
{223.47,652.78,3.64,5.98,"Q?(######,##)JY######q##Ejn######3##f]495####$##LI,######,##y^b######a##Rjn######5##n(>######8##I[+######*##TCe######)##/jn######,##zz?######J,#Bv)######*##k?U######/##yin######.##v':;5#:5#>,#"},
{238.68,669.38,3.96,5.98,"^)A######(##oOs######)##[+L######M,#pb#%##=5#lY#`|G######0##jOs######0##8wXH>#E>#X5#,l#%l#.,#0,#73D######.##GPs######&##`eW-c$.,#&##T>#$Z$.,####{C<######$##^Os######&##H6SPG#.,#'##kP#M5$###%##"},
{738.48,712.62,3.92,1.78,"&`$XNE######4s)dJ2E>#$##5l#4w(C95~P####G%#BaZ+u#&SBUv*###'e$r`ZxY$###yc#Ke.Kc$Xd'<[(###mG#n%?Z[,h7/######FG1q~ZN5$###yg%}%.w%+vP$sY#.,#?@&ru$X>$bG$######va-$-'=5#@5#w1/;#$Z#$Cc$#6&###]5#IZ%pb#"},
{560.58,738.33,3.55,1.49,"7##Y?&o#'###qP#V?%Qu%###z[)$H%}>%'-%B,$Nc#4<6Qd)<+/Hx0pb####a$)90*3l$&##N%W}H(%##w%%|$,ED4^h-i~+8aW2?%OG#%$#i&0dA)~#&%##IvHB]J95#F##Mc%8aWX,%)##=dU1##3l$B'#n5%1.&vw/'##Z5#bv<?E<95####:X,'9/WQ'"},
{298.61,763.91,4.06,1.97,"QG#^5#`c&_H(qZ&Bt25L85[*otDk'49%&qA.)c$###<++u[,###7##C93D>#&D<d&'E9_.c#f;_C%)(q9L~%xZ')?$XD5)-&;,#n>$j&295#Xv*E##%>_oV9z9_Pu$~I)gG5tZ'<w%HH&.##1l#2u#)H%######+##9~(~[+Cu$v%'lP#w@*QG#gk.###&##"},
{298.61,763.91,4.06,5.10,"###&##QG#rG1`G#a.*Cu$le'6~(hd+###,##4Q%###;u#'l#`Z&.##jQ'I.&eR)jG5`'_Qu$g+_:r9Wv*L##$9295#;,#p>$F{5,-&H$()?$=(:ie%;)_T.))'_0c#._;mx&,'3D>####8##nF,d[+4l$###6%&Q&.#5D{y4z98Y-*L?&7k2VZ&lc'F>#~5#"},
{364.51,763.84,3.96,3.64,"###%##|y8###ztL###o^8A,#CAY1##liCg5#$[)D##(AY=,#)##8w&`e/###-OCwc#[1:>5#qAY'##r{9q[&)@+###tDYMc$-Z#&B+;Q&;5#_*FP$%481C5#gBYpZ&;q8J6$M7,V[%9CYtY#^P#D5#PH&7l$RA/tG#Mw.R>#nJHe|>Hv)C,#M[&nn@9U9&##"},
{556.11,797.74,4.10,3.57,"n-%fq>sb#/l#$gXxm,D>#~,#BeXF5#6~.Z##uQ)d,#teXD5#~,%]P#-c#A@PCwV.,#^5$fz/~gXr?%Y7.O##V-(@p$sdX&#####A##cu%(`?{V@j5#)&06.*YfX;I(Hx2S6#tu&x7&3eX8,####&[#YR,95#if2]R&A6(B5#0N5S-M9R)VG#Mc#5%HY7.;5#"},
{351.03,810.23,4.04,5.80,"HO<xY$v&*7,#00Nar>0,#bQ#sG$}.N/,#ZQ%[P#di5###cn%~&0$##0TIuP$i.N47+ml%w%+]m*Gw,K>#%ZB3H&P5$###&<-e~0S#%FP7G#LcuLJl%8),HlJmV2[%/lG$6v&Y'0eY####.Z#}|;qn/]?%6f/(R%:V2P/N_%/7]&pQ(0~-IQ&#M(vG%3,#95#"},
{371.40,814.46,3.59,5.10,"m%-mp(qrUMS*XkC|P$A$%dy*jv*fY#S7#3JTA,#%Q%X##RmUhJ..,4Fg6;6'#pU<u#)l#uZ#Z#P.,#&##DB(Z>$###c5#p83zl'5Z$e27/U5NoU)l#p6)s,$}mU######~##HZ&###<,#/Z$#l#2,#v?ImH(bmUH>#4w){]'0QQ######>$#;c%###$##'Q#"},
{380.12,826.69,4.02,5.19,"F0)UhS.o14Z$9X>#^'cB5F?'VeSOG####hQ#sR-###/##%@(P?'.n,fQ%Ml:{y6|G$h4>D93_fS^P#`,%DZ#Xp8###$##X5#kNC###j//Aq/8R*nY#cpNDn*tdS1u#NH%c&'=1;######F##xeS###e-(T5$^u$([%>xP^u%$hSY[(+x0W>#^V=######Y>#"},
{529.95,827.41,3.80,4.11,"$##=5#7S-%##>,#z#'h>$ZH'W[%R$*<c#<c%)C-(c$.,####r#'&##yK5c>#7WBH>#;I%<@A>6&k,&*`$gjEX_3>u$iI#d6)9w,,##a(;>5#`(['##tl%9e(<>G###2W%OL7W?(###{E*|I.Zc&1['^e,G>#w:ME6'Vc%F>#}DQ95#E8'R5$J-&SG#VUI95#"},
{529.95,827.41,3.80,6.27,"Al$g&%)R*###xIX,d%{k#Z%#QIX{%$_~1h0#+u#|1#PIX?,#cm+B[#T1<dG#~JXM$(Y,%'A#-QMdb2HJ1{##C,$42$++H1Z$G6(L##,A.cP#k&MdV=;u#)Z#AI'yKHpb#$$%fc'_n$D,$fw%$##-##SZ&/,#8c#dc'gP#_Q$^>$VH'###q%%+u#3##(Q%<[$"},
{283.88,835.53,3.76,1.74,"gY#k?%Km(3v&5]-uw,rb#<,#PP?fl%;c%I,#aV?sl#hf26m#E>#bp)%N<}#&H02]=7ru&>5#OoU.@*/,#M7$F6Rc6%MT0Dx%###M/D%Q%###??NR<5D>#-##*pU}h9/,#B6#@L3qrU=6(O5####bw*######Jf-<7+######G<4NoU######(##~1T######"},
{403.93,836.12,4.01,1.89,"fP#X,$rv(/$(1l#^;0@9595#-q2S)<Ll%+Q#x>%_P#0fZuY#_>#%A*D/2~P#2x0@N'EeZ9,#^gZvJ+I81,S$w?)&-%cfZCQ$###%H#>g8###1%-|,#hgZ4I'NgZTQ&kI+[n&A7)S6'qeZ######d##}e2###kl%|8&CRTB#$76JXH&=Z%<l#i..hZ%e4K%##"},
{403.93,836.12,4.01,5.13,"/C8###BT1U$'AR*X,$QYA&@'_/[3u#0l#S%%k7/######C##q7U###'J*vl&S~,|I'b1[8?&>2[p?'97,S,#D'7######Nu#e1[],$nQ'_Z%{J22I$?2[1/+&/[4,#d]3Jh'C06eG$eG#1J*S1[#c#sP$/,#F-(Ol#5;6tL8-05.,#2Q$/*3fI*}u'I>#dG#"},
{411.78,840.33,3.60,1.89,"###2##0d'bG$'##fw+c5%###a5#8J.-H&###&##D5#CJY###O>#.I%?U50H&Vd(i-:OND&##}wNB)<Kc%6I$8H&;5#wJY4Q$+##h.'mNF###K{?0z$wIYlG#2LY57(&@*c@$66&B$&qJYT#$###H##'3C###b./N,#0LY;.(LKYFc%k>$GI%4m&&-&ikL###"},
{392.32,855.37,3.66,4.98,"y97jG$######a(h######|Q#7q<&##XR+sn*OG#'##x)?Hc%GV>######&##d(h######)##e>P&##,~,{P$_5$tZ$PEA<#$ug:.,#/,#$##>)h######%##Kf^I5#ld)y$(@-(Cv#f_e_c%N&3/,#0,#/,#q(h#########meZcc$$R)P>#Kl%6m#?PH{R*"},
{732.05,856.21,3.49,0.09,"jZ'R8%'$($##FnYv>$D>#-y#T^:)6#Fq3-c9/,#&$#PpYs5&=95).$J82,##MoY6[(95#b%#s*FsF/8o2j?##l#h2&YsG<,#(J.O##CI+sY#|zTVtF0,#N5#u6'BiQgY#v>$Tm*9g)OG#>n#N,$]G#~P#`P#w5$HI+bP#F#$Su%+@(<5#g.'8m)&#####>A$"},
{732.05,856.21,3.49,4.24,"}k#G>#X#$H#$%l#dP#l5$PI*+@%y5&TQ#xc(?C.pb#'##;5#KI,Q>#;x/<##9[U<5#xQ#DnA]Q'fY#HD#nvRJB095#{%#Ie.P&34,#b;;YG#t1Z$##Vu#fR'pYF*l#5O+M/09l$8,#vQ8K-)u#'+##fg3cP#A`TI#%Il$ZG#w55Z=E7~+:5#/##Q0-I5G###"},
{755.32,864.89,3.51,0.82,":%*1H&`{+s,&4T-5p48u#G##-y3*?&w?&0A+ZP####<Z#._0D#$b#$Dv`LQ'h#Ns,%WT3u$#g(>.,#%##2f$.c$iG$&##c?&,6%&M&9p`'##Uq`_$%)::M##1=DBu$###B##@5#rH(###$##T?&OH3i/`b>#[JZC##Wv*L6#DD7Ac%###%##+###@)######"},
{489.93,888.64,4.01,4.88,"3R*95#8,#Eu#wDB95#g>#Qv$B3q######C&#%N@?u$###@##9R)A,$/,#-##pAX3l$?5#/Q#?4qE>####`$#5L8'n*eY#5##*$(.,#H>#>5#rKbZP#<5#-H#k3q0,#OG#1$#to57d%#J/%##GZ&######C5#in~######o5#K3q$##D>#,$#?e.?,#J80$##"},
{558.08,896.62,3.31,5.02,"3I%qW:95#FR:*'EVJ1###L8&oT0m5%###.##VR*vb####:,#=B0&v62.,lc$WaZrm+E>#D?#XSX.,####A##t?)95####M,#GK5jQ%96'^A$M~Z######1:#)AXW>$.,#9##sl#5?'###$##ey3nn.@H%AS&I`Z~?)###..#`PK:5#kY#`G#_,$Y>$+l#:5#"},
{919.33,899.32,4.11,3.08,"](>######'##J<q######H##TcR#######$#Hz:eY#)##6R%?*F######/##O<q######%$#N@X######E$#(x+Uv)M5$1##^|G######&##3=q######4##58Y.,####5##e>$+l#|,&=5#9D?######%##X<q######D##(lM95####8H$i,%95#&##=l$"},
{615.06,907.90,3.84,5.03,"m1/~0^6?&RH&`a~nn0]P#RG#-x0###~G#h,%/H&###&##+m&BM8M92/z/5#IO2^#H%hZ%D.(h)B95#8,#Wm&cG$###'##mc&`j?.HHtR+)]+I0^p#%Tl%w[#46Npb#A5##I&q,&,l#+d%JA.Eo0h]-AI(H~PR3^f$(T5$'m$udPmP$###<,#=m&kv)Lc%}b#"},
{542.42,917.69,4.01,5.64,"9d&.,####8,#ZvT######g9#6oW(8-6#$['#mT(,K395#&##gf1######1##<oWrY#W>$b(#U}G|C.XaH.A#+67c809v(VG#$n)#########PsWX?(ZP#<##{|-Yg2G'7XG#Sr<PZ%:$(^5$*###########_B%ZP#######'@;}G%hY#L,$UL0GQ%j-).w*"},
{586.49,918.63,3.98,5.89,"^i4+~,9#$<5#qP47q:95#5##cg)DZ&E5#}P$4Q#@u$G5#3c$fy-GI)Ye,;e*p/XQd(.@)HA%|^1Vh3e-)Iv%&L#SsD95####W2&/r@5$'`5$^fE8K3ku%Zu#C&)Q$:7C7pZ%_:(c;;?$'~T-7'#`$+######<;$go5/,####Ir&_M?Z@)t#%;8*Gp3o$)6w("},
{840.91,923.33,3.57,4.36,"sZ(kG#.o,`4;S..M#$WT0'S'p1>95####]'#SH(######YI#rJ.o-=@U678)T@TIT)1o1CI#:PL######j&#du&######U@$Dv'?x*GDT^Q'(BT'$&1%,p##4@T######d~#Yl&(##D#$Mx)-H&,$#LBTd@)S?T%##m7.tL##`A.##,u#:h#+u#TH#?Q&^Z$"},
{871.09,939.62,3.94,2.41,"Vl$<7&N*F2##j(4xA1rl',0)Ku$gl%PJ,3p1###+##Tp1{k#3I+0##t.JJQ$z?QkY#}R)L=2)I)Wu#+CQ%`9###L##4r:=6(O?'kY#EKIvP$3:NLg8fp1Zl#VL'f@Qu@Qm5%XG#/A+4p.+@+GQ$d#%Kj=###/$#=@,1Q9pb#_I#3cNK[$Ev)r93~]3n>$U6'"},
{871.09,939.62,3.94,5.35,"XR)H-&x%.Se-N?%d#&cK&kkK&#={k#A%#<@,.#G###oP#{G$wC2$m($##zQ&zdRh,%r^(zX@-(7hc#o1QN&2hfRVG#vl&U,$g2?mP$###8$#yfRbX>8[)`Q#:&..t1LdR[G#7gRbl$dc'0##GC5Eu$###,##BR(dL4Id&4.*wm,iB(-p.5w(2sDM>#95#RI#"},
{887.94,973.62,3.61,5.05,"&x-P''K-ScG#K[PP5#hv+1##C,KyG%###-6#'y-JM<3Z${,%/&-wK)05M%##I-S2##,m(bH#sZP8#$###kJ&u^;uc&]a1gk=eJ-lP#x05L>#Z.Spb#;v&9c#G0SSH(###Xc#C0SEA,8#Bxw)BK*:5#;Z$:#$F^3{k#1c#c$'ZE3:$)###l#%*t1p+FYl&sY#"},
{916.31,989.92,3.61,4.43,"2-(######a##?8Z######[##4#M?,#n%-$/+:?O&##Zl%u^'<m).,####E##e]XD>####/##99Z###?##`(3}OF###]##(M7a5%;5####9##=8Z;5####N##_:ZeY#8##k^+*B0Fl%)%#q=G>5#Qu%###%##^^2-H&###$##wgJHU:###=5#l2&Z7Z:##<Q&"},
{66.27,112.00,4.46,0.05,"$z2TI,{k#9##*d(},&Z~);[')6A95#W7'O5$LP5^'5eY####+nO######H##iXI###@7$2y.r@+###-M%D>Lkg,2[)X5#tD>QVR######$##>SR###2##v6&1EB###x$#kRRL5$###C6#)RR,v>#########QVR######;w$=kH###,##H;36v(HQ#UI)m~0"},
{66.27,112.00,4.46,1.99,"<Q$~d$Y;;DR+DLOwl%Hl%/,#pcFbw-0,#*c#}e+`T.R?(x5%ll&M{'K07X,#7TU<7(mP$H$#lj?tIK###(##<n$S{SA,$###fI-]e&N5$f@$[<ApV<###x$#WQ&R1NOG#%##*##SZ8|Z)###$H%IZ%l.,H8.1l#KW;aZ'%#####M+10Z%######0V#:$)###"},
{591.82,127.93,4.72,5.48,"Yk1.u#Su#=l$NB-c>$aT-lH'>/N;5#=l$xw#to3######r##.67.,#$##.,#t3XlY#em)=9)/WA.,#V5$%E'aH(.,####UH#e18=5#-u#V$#l0X######~(#`QO.,####B$#UQ%95####$##2NC$###l#4M#/.XD>####4D#f@,[P####:?#'['#########"},
{751.04,141.68,4.75,3.95,"v..{R,bu%wz/eR*k?)NZ$RyPS-),u#8,#m}1I#%###'##9L*MK.n/30U0cA,|zWqC<6-&d.&%(:K>#jY#1r&###2##Ec%G?%.J/Qd#RyODx-.xW'Q#.(7j-#T6RcP#/l#vQ$###=,#w,&E>#G-)';#A{?HT'nvWa##O$*^h#E827##iQ'NS%###+##%[(Tu$"},
{724.87,144.91,4.06,4.76,"4Z$_w$QFJ&##}%SSQ%Yl&t##FWB######<-#L5$######>H#vm(Pv;>;@$##5yV+m%j>%+$#^lQ###$##o@$95####Y5$[Z%G-EJvNX>$]n#6xVqG$.,#jw#zvV###]G#-$$D>####ZQ%(Z$HwV~H&9c$BW(9yVC8-@Z%kA(J'V3Z%=#$P5#gP#0u#dl%G,$"},
{774.83,149.87,4.20,1.16,"eY#71#y^<&##qh?6r(r@/g$#_}H0J.1,#Ys/}b#HK09.(/h9k>%/0#AHR&##uHRrm$up:z7#i2=S,JhK6T?#R,#eIRXZ'###+c$D##YKR###<KID@,fW<K5#?'*4i9<z;?5#0##8T0EU9######9##n,CE>#]6%|Y$%5:'(9r@(7w)Ew+EJ/kY#H%'N$*`G#"},
{529.98,151.06,4.93,3.50,"}Z$CU5yG%###6:^b[*###p##i8^%##QG#l##V5$XG#gY#/,#i.,I]4D>#9##79^rZ(###X&#Y8^###1,#)%$E>#8,#&$&Ll%=A0G>#gY#x>#g=^######--#s:^###.##e-'95####e5#0R*6l$###%##KGBm&[###$##:-5k9^###1##Rn#}k####J##`#&"},
{758.31,151.96,5.03,3.89,"0,#L;;R,$H?'fm**f/%I(K23g.,/Q%w#%-fNk?(&?&'##+`-rQ)@#$HZ$my.L_3vo4%y/uy/hxR'g6bu$)z-jZ(&##TG#|_*un1@,#ne$kt@G^8H%%AVTlS,SRTl>$#J-&.#8K4gP#?#$--$o5%+N,Og.2045w-kC#Q_=C2)dQT>##/H&h1$wu'0##G6'4v$"},
{70.05,162.39,4.33,5.98,"yS3######*##H_f######O,#QafhY#[u#G~,x#&eI&yv)1$'Wp9######/##(_f######b##;Y]'v'$l#T5#YR#gC.0aE###]p9######%##G_f######-##[wW{k#,##q>$m,#4I('U4pb#V&4######$##|^f#########suQ###&##1,#kY#$##z,$zY$"},
{534.49,164.45,4.79,3.63,"a#%0f*126]-*w*8:R)=l$$##_5e[>$###%##|%0/,#D>#)##G?#JwMeQ(95#FX@9<=pb#(##v2ep5%###[##3(;%##wb#xb#Z,#[~0Mu$6#$VJY`&3.,#D##t3e[P#####$#y^:###&##+$&OG#bG#.-'VG#@_9$##/,#5Q#&8e######*6#oHK###/##jl%"},
{819.33,164.97,4.74,5.92,"Nh:yP$c9-b954_,$$&8^U,6'f1O/c$s5&###0R'1Q$<c%###.t@JB-$%,'##SK3b0'(~UVc%L]U:,#+e-&m$gB62,#AZ%T>#UT*$}<-u#F>#Wh1gJ*T02:bDN~UPG#eu$;IE90.Vc&pG$7u#G6$Ce,S>#n#N$X.qu&]Z#'=EI;&`T5z,%..+YR$mT4^>$[P#"},
{527.76,173.73,4.10,3.49,"fH&l~/jT0dG$jRKdc';5#&##Yrd.,####&##bB6.,####(##:m&F10Vx1G>#?-E9/.:5#(##psdM5$###Q##rC=%##E>#Q,#P@#_iB######_Ua|q;###1$#npdPG####%'#K97'##2u#f>#7H%O[+ZP#0##'qd?6(###C%#ypd######H&#Uv*###(##cu$"},
{612.71,176.20,5.06,6.11,"oPL0,#)6%U5$yvfPG####(##`14.c$###&##,v%95####$##yB81,#mZ$J_7osf###$##5%%:1~95####.##G@+.,####8##x/.$.)<H%?e*=rfD>####=@$&qf######C##8J0######=##q04Xc%Tm)###:rf.,####+##:f^######.##6e.######'##"},
{746.04,181.69,4.24,6.25,"o,#3^.C}DxY$q;(*lKIQ&8u#K%%]RLR*=w5%A'2:8/X$(9@%e5%`I$[mU<,#v_9#]0Y.-hT+i[&T_;wo$ue/=n$.y5o5#y#'Ru%###onUJ5#7nU###Ow,Ye%+C8###t'#WL6K>#fY#B0#y':3l$###tpUtP#eFK###VI)X((IJ1###V$#;#C######0&#uV?"},
{221.77,184.02,4.32,1.67,"095@5#95#&##aTaR5#&@+=##qh?F6#PnZ[Z#b2B(##wI.^Q#%h;F5#.,#$##gTaD,#$%+T##Ei?v-%[GGC[$E|C+?$IB6:?#4M>5,#.,#&##NTa;##./1Z##daC:8&lTa]>#KQLWZ$:@,7##u':1,#[P#$##nTaU5#2e.1##c08cZ$t@XT?$t06yA,%[)(?#"},
{265.81,189.38,4.69,1.65,";e.qG$D>#%##~0c5u#.,#N##YsEVV7?d*,##eT6;e(bsH0##<96M>#D>#)##T0c@l#'v'i##^{>oM46)<UQ$7sC2@(gy9;Q#&C9-l####(##I0cB##~w/S##EOCt/)p@YvG#jcO26%,A0U,#i]6=5#/,#*##K0cH,#tn1O##1NB,-#FJ~4H#9GM*##37,m,#"},
{518.30,212.38,4.08,2.88,"###%##`c%{G%ac'###r#$F7)^g:###^>#lj-'A-###N##GX4G-$1@+tb#yb#]iA^?&9Q%QQ$OUX:5#D,#]W.0q3D>#.%#{SXdC*XA1#?&5##h}H1w*###T&#kUX95####-1#ic71I+-##dm)Z;>.l#I#%`-#xuR.,####r(#BiB######r(#`v'0Z%$##Y##"},
{576.57,241.92,4.99,3.32,"3,#qA*ZP####E@*j'3nP$/##5g6fG$=5#Qi395####B,#XdR5u#^v$[S*L-(iXHF$&Z-&SR)&C_H>#ZR(Jm6]c&4,#}i8;jAWl#Y8-DV4{R,xg9OG#F?##`;{D_6,#.n+sd&/g5.##=eO`G#T8#J03:kB/,#,a,o#'%##sG$rC_iY#3-&q6%W;?1-%twOrH&"},
{837.41,250.39,4.11,0.81,"@&/###n8ARG#<j]###xu%H>#Oi>###Z5#}k#######A8#0v(_-*]^'7h]@5#Xf]Z5#Ce-k5#w.[###.,#_5#95#xb#D5#_>$`d+%#6h?)Im'1f]6Z#cG$^e#TvUg>$.,#Z##;5#AZ$OG#$##Gl%vb#|u#=$B#`A/,#%##j|)<.-d>$fY#?6#`P#uP$######"},
{412.37,253.43,4.77,2.87,"M$*###4,#W,$)[U5,#9,#u]*7w-aR&Yp(a2>+##C$(lB%-^8Zw/$##4R(dP#.[U%Q#7[*~##.]2i-<O[U5,#e%#7~UY6'I#%U%/$##g&2H5#M~U]P#p[+b5#58.Td'Y`U0%-D$#1U7t{)J(=oZ(.##En.O5#IlEj>%H5#hw'C&&q~23&$Kf3&$#<R+:(#`EF"},
{361.12,347.02,4.70,2.06,"'##[,$%Z$+u#1##zB.gc'###J?$N]0+d(.,#0,#:H#{m,D>####C?#?J0,u#U].R]B5f]C5#dk]r>E[~/#$$B$)wm#Z(<uY####I##hL:###ax3@@$<l]R;95g]h6&^L3969[-)Pn$E?O@5####Y##B[*###95#1##rz0V81ol'Q>#Z$'pU3PG#V,#nq:uP$"},
{593.60,355.47,4.10,1.27,"^-(Uk.vA4A##SbAFU.~[+#J)T#M###&##nr/f6*(##95#QQ#k.&YM]Sd+*##2M]oI+:Q&%0#$K]$#####^&#}>&a5#eY#3##Yo/oB84,#{b#$N]r5&###a##(L]=l$###[##J>#yQ%.,#$##:/2>5#@5#]l#GK].,####&%#(@G.d)###*##h##7d(.,####"},
{123.79,357.99,4.77,5.97,"{k#######$##rRY######<##Haq######/##f&4gY####(##B,$######%##GUc######B##iaq######/##7C9fY####'##oP$######'##KK`######Y##:aq######9##{99eY####7,#M5$#########k$V######'##6aq######&##6B5######4,#"},
{527.76,372.69,4.46,1.98,"tj=}c(.,####Q@$R(7A,$###_JC)%,###$##?108H&######O9V[5#$m(*##H02y*13[T_G#5;VlS.ku&z>#W#EZP####*##77V(##1d(q~#Le/e5#VfMM(6u8VlY#-v%iM/IwKOG####2##5?'###4,#,U*95####+?#h;:Fd*###)##1O4OD?######0.$"},
{318.25,418.42,4.49,3.43,"q$#nGJ.,####21#D3EZP####?7#v/*`GO######+%$~<EeY#6##Oe+&~+###k</B&2f5%6,#_xH76&-i7Oc${k#2##Y+4,J/######Zj0,%-ZB60,#n-'00'v#Q###z,#Rh&5?'###21#AM<'##P?%)'JK96}>&###s?%a8G'/1######]q&mP$###K##o?("},
{475.96,453.94,4.36,6.22,"N?#-CWC,$###aN*2;<.,####A-9{k#######+v%#########)S%xK1#l#.,#E#8g(>###6##oEWj>%###?##%91#########d5$i5%Q[&QZ&p<Eu?*<5#1V$k@WD>####f(#ZT4.,####8##OG####o?$_SV>%.###4,#STFiB8######?D$UZ&eY####1##"},
{624.57,495.71,4.26,4.08,"TgJEQ$6&1.,#FWT###*BHPG#,//###-|25~+3v(######U[&SW=?l#_d(pFDfhV$##WK)c9+~%/(##bl7KM<OG####G7)5A.&@+###&##O-EfeV4Z$c,%ry(n/1c[%TeVeP#.,#)##oB6gG$;96######KZ#c5L5u#D>#P,#7R(vc$WB6######'##b7+ZP#"},
{404.16,509.96,4.14,6.25,"Z##_JPm>%###A_&{^9.,####t;+b-(eY####3##(F0CZ&###6$#{o1Hl%###6j1*]3###,##fLP8d%<.-9##tG$G(%*HP%##eP#[5$xd+<l$l)ARc&E>#/B$4HP3,#fe,?1$nP$8##fLPCc$95####kH%k+G9@,###%##dc;)n-###[l#Z(&.,####Hr*V#%"},
{784.02,511.18,4.11,0.70,"^Q%6#$$##Q>#E6'.,####LI&q#'###=Z$Y]'.,#,##7PGtu%X[$Ud+###$##<|7fH)###]5#+-Q+##d@+d''W>$&-#loX9$&Z,#8K4######2nFqw0###2##:rXoR%>q<w?#Zu$RP2-nXG5#7##z{B######qx&vtG6#$###)jUBJ-|_?###5~)h0,gpX/,#"},
{286.77,519.30,3.95,4.32,"*h&6EBuP$###<pIKQ'######6q3m#&######8u#}9-pb#&##Qe,1u#zI&318|0UN#%9,#*/'HOAr+C###)##$###68uR.###(Q%vY#li*7/UK2@Z~(q~+QX7=m([2U~#&*##fP#J/U-H&###?Q%2/,Mh3].,###/.#;/Ulc'$##x'.r&5jY#>-'BJ01,#Q,$"},
{395.94,525.71,4.90,0.87,":##4&RWZ'###0))|>OuQ%mP$5)4nP$Ow#5?Q7#$###3##r$R$##|T-L~.OG#hk<U@*sY#tP#T(RZP####SH%~R+.R(.,#4I*ZP#4,#qA)S'8d080,#@?#d34F$RB,$###5)(4d(]n)>&-8g2uG%$##[$$}*Ck$)#%'q~,8j>v9M%<?-u#BT(_6&ieD:S-.]."},
{491.42,526.58,4.97,0.01,"#.#zKTs~.LQ'z1%($JdQ(###N3./v(######I6%>u$######5$#MSSkG$###o43bYM###+##$NTg,&###9##z@,Qu%###$##6,#26&w]-^H(/*Bd$+bP#_2*BITZP####T1#}[-Z>$?5#4##Vl%-c$;S'U>HHv)95#+##PgK=o2+u####bU$|>%D>#1,#(##"},
{184.75,538.06,4.13,5.94,"{,'######*##fmW######U##shk######D##R(;95####$##M?(#########K;j######0##Hik######2##jV=.,#######wu'######%##rhk######1##uhk######.##e_?######&##}>&######=##<.Y######]##mhk######/##vK7D>####%##"},
{802.54,538.83,4.43,5.62,"YQ'yt@###/m9{&(_xV(##**07R#sxV9-%F/2###wO>3$&>u$T#%?|Vnc%TA+wzVPyV>5#]d$?|6a6R###)##*##4*@######(c$uH%JB&k06uvV1,#U,#v&)8(6qQ)###?,#Jc#<$)######.,#.##)B+C6'IJ1###rG#>&($@*95#&##y>$Y#$xY$######"},
{264.81,543.33,5.06,3.16,"c2$Wi@D>####'++4l$$##D>#U#$Y,%6,#dG$@5#ZP#;,#k>%^{)V92zn0&##O5]c#&;5#I##f(8z6)xH'*H%OI&P-)2l#Bu$T~/zZ$1ZHkF890]wb#$##3s'Rq=1,#R.'HA*3['+$'nw*hw,r8)1$(Q%*>DR{//{k####]t.Ul%95#a5#=w(E>#PG#IZ#,J-"},
{451.48,551.93,4.36,4.25,"hp&rrB(6&OG#W'Kg,&0,#FR*k80r>%###d&._P#VS,`P#O$'%A,vb#Z&+4y4UgS8#$bG#s.'xE@725pb#*##1,#?7?,d)###>u$]>#&X-)fSejF6K0)J(VN4sd*;gSW>$*##$##S2,cK795#%Z$rw&=PDA@)2,#u_,c6RQ#%rP#`.SxY$###YI+<93Cu$?u#"},
{348.83,567.88,4.65,0.37,"`/#IOH:5#$##j$#x:0}kN###WZ%-q&dHRK>#>t24.,:5#=,#{/&oW?eY#&##c@Aef2-i9Y5$%v'T>#&JBsT6u-*###7$#F(7Se-%K1@u#XH%-IR&Z$F##~8%996###@1#z<D######`'#ZHRvG%###s5#Z8M=L:######-<+q6*.,#3##/$%aP#<5#h6$Km*"},
{348.83,567.88,4.65,2.63,"#&(?94vG%###MDT%Q%######=//8H%L5$$##7,#kQ(fY####*8,ee+2I+S5#+CT/v&###|##hg6VoF.,#&##s$#TAT######ql'f,$t6(JqRfjH/.%)z8>t1<?'+t0Tr=JH'H7&zN>>6%iG$3?&###^>#>BT.,#*##rAT0C2###xP#m@Lp[-I#%&##hJ+e.-"},
{727.68,574.47,4.11,0.62,"2H#C9SHl%###DSHJD:OG#+##_IQ+Q%@c%v##a_/_l%i6*,##Qc$j[+Yw,$##4;SK#%SG#E##a9S*##yv+Bd#BT-q7'OaGlP#]>$ac$t6*F5#v7Spb#^P#1@#H:SYXDeI-s%$*%(Rk8B6SxY#(c#U#%$##4,#q9-~#&###:##5y'%7SC,$%##mm#MQK&Q%###"},
{640.07,580.72,4.22,3.71,"hc$cI(Rf4>5#c%[D##Ro2m-$zMB###UuDJ6%uG%###[w$Sv'*U(8A0iZ'Fu#f&[gG$6m'R1%_)Bw5$G6D1w'95#J5#Y[&*&0vk1d%0###%##d+[(I*D>#Q##ey-4GCx$+7,#)##Xo.}>%eY#Ta+^Q(######)BCr5&######4r)Vy8######&@&PI,######"},
{371.79,582.55,3.97,2.15,",[).,#zc$R(22bEZP####y:#^>L###]>#`R$-H&###EQ#5h6zH&G^,uy4I-(Yw>KA1.,#K,#c_YOG####W-#6J0###4,#6f*NZ#jO68~/###K_Y~v*###0$#L~Y/,####V'#?%.3c#iY#3$$#m&XJ-7Z$zS.R~YT,%###tM%8~Y'w,###{&#kl%Qn..,#(##"},
{673.88,583.67,4.52,2.83,"t5#($(<5#0,#vv)D>#&##,Z#FA]######k8#OA].,####Q%#?v%Hl%nY#&##bD9^#&###7##2B]######Z$#]B]fP#E>#bc#/l#{k#qH#Oc&sWC=c%5##?Z%sA]######|$#mA]$##0,#K*('#####/-#2-(mf6###^##qe,CA]######IU#IB]{k####9C#"},
{569.75,596.21,4.35,3.20,"DUS###*H%###,g5$##=]/cP#_QS###TG#],#n_@###%##Jc#AVSE>#.$'D,#F,E)m%>S.(##0RS+##Z>$U##74HE>####g5#YSS95#Y7*ow'ZKSfl%ZP#)?#cQS2,####c&#)g7###.##_n%(%,###7h-1TS}QS###Ud&^x(G,P######@(#E-)###n>#xI&"},
{580.69,602.40,4.34,2.89,"9.-s##qf_N>#Cf_&##cc'n##1M>9##~A2G,#Sy8)##L#%t,#SS1[##hf_Q,#Af_2,#pl'D$#h)ChG#Q]4_,#Rz=$##qP$@%#Bv)M##:g_X.'(g_I,$3u#:S$@W@*[&g,&y##K3F$##95#R%#A,$2##mQG4HMc*G.,#H[%(y-1g7E>####W'#T7.###$##iT$"},
{569.79,608.02,4.23,3.17,"W(W###&v&$##NJ/J>#mx1aP#k%W.,#95#2##9|D###(##Bu#[)W###5Q%###];>%##1q5qY#q$W###K,$/6#vPP###$##[Q#B)WD>#8-']5#7dI(R&Yn.)##=%W:,#7l$d##F#O95####j,#B[O.,#Ad'v.(pdL'?%95#$?#q$W2,####D&#sx5###*##c.%"},
{798.29,610.13,4.16,4.00,"M&Gd#$i6)rb#dEW###^R=|Y$tl'###mQ9T]2uG%###`&#{jEW_:7,#lc$>h7kCW2,#ST*lJ'iw/T5##/Auy795#8,#|w*t&3{n1###C##`p2~CWHv'Hu$N@$`A-aC*4lLA5####*H#9K3ZP#gu&###B##d~/.YDf,%95#kc$f6'UA'<@,###?,#tH'p>%95#"},
{571.87,621.07,4.45,3.09,"8(XGl%Mu#Wu$&T2###]~)gu%:nT###3,#@5#AV<###I>#YG#D)X###QZ&%##kD>6,#Bq8cP#j%X<5#f5%b>#~ZR###)##`#$#(X###qv*####7R'##G:7[G#}$X$##bc&h?#&lO###$##wH#n&X###MJ/6H#=%P.v%|?*,##*%X=,#c5%#$#MND.,####.6#"},
{617.94,647.91,4.72,0.81,"2{).:34/W###O3WhQ%W]4###_w*n>%c,$D>####'l#[v$g,&d;)W0W/w*ru&L2WvK1MQ'}>#m=HBu$3,#/H%###a>$S,$E6(?d)bZ'cg#_HMC.WN5$9##V3.s6T###$##5$&F>#%###Z#mQ(;u#H6'gv$Ue.7.V95#,##wI($uK######,##@l$/,#/,#H>#"},
{666.80,648.20,4.47,2.02,"()Z#########?(Z######6##.xQ0,#OG#F##C#$L>#Z>$I>#N'Z###$##fS$i'Z######zv#T(Z.,####v5#t$*R5$/,#2,#-o2######`.>u~2######%fA)EB######]*.`-*$##&##l5$Yl&######E,395####$##AH@###(##5u#/B-W>$&##[5$Kl$"},
{541.73,688.05,4.33,4.83,"'H#~h9######M^0?R+###/##3@+###%##hu#.,####(##$I&AK(LtZF-)E5#LtZ81;###E##Ri@######WR%.,#######eR%Vq;:f'w^:FS%zoZM5$###D;#jPO######Y%$fY#tb####lQ$#803,#yo)sDSEdV######t`(-S/'##,u#:@#.,#`$&h,&[G#"},
{556.83,693.21,4.39,5.25,"B/%-yRQu%###38@#J/###-##pS0######k###?&E5####G##fh7|g/.|C4c#ijYRu%###>##N.V######@##>u$i5#V6)@,#,x0F>#2H;L+BIeY###(##>+.2RV###$##G-#0Z%)##+A,/l#W>$###xI#efYL@-###&##=R:mw0######B$#AH'###*Q#4Z$"},
{336.83,748.51,3.83,3.55,";m)o##awW`G#?mU;##dg9}?#N`C###k5#_Z$######B6#:f38~/$##?|W#@&|vW###pi?lH#bsG######Y,#[P####%##=H%AK2YD,BxW(Z#?wWLZ#o4KR##z+L######B##uP$jY#95#^G#:%&~|Wvn1###cSGXe+ac'V##;{6D>####7##o>$B,$###%##"},
{359.48,750.74,4.96,0.47,"8##Ed(2k7bG$<{-]B6}@(&9.hB/_-(8M;am'E;9~>$?6'2Z#&Q%;##)>>Lc%:>M+##Z:06hMJ7-.##g2W}<=7n,[P#;m(Hc#*v'7##rjG%##{1Wac$iC;{Z$HL5v$%7/WYc$47,$##3.,f~'qG$|G#NND$##p/WZu#G&3;##{g8cu#a.W.##0V;'##OI,yG#"},
{359.48,750.74,4.96,3.52,"]l&'m$8::'##+wW:##P]3%m#:x2p##+xWu#$4NC+##W?(K?#f-*$])T6))##gwWtl$_T3me${L<)R$gyWCv$wjI&##[$*S##4I*LZ#*S->5#0{W;>Bwm,@##?M3whP$GL6##vOBO,$Fl%T##;@+Gc#ZD;S,$@L8eI)v12^d(yf,NK0@q,A/1(Z=D>#+##-?%"},
{555.37,768.61,4.86,1.75,"F,$Kc#>x+hl&)###y*sh0pm,me'/o08-$Ze/xR&{Z'~v*n,&b/1.T'9I&/d''d'SF0;r37[*JeTlK0-v$`>@)7)9H&=o$LeTT6&ws/ye0OG#Ue/wh&[p7t.)2iTa/2(l#8y(j<-|cT:$#;2=###=Q#v],Z[,pb#9##G%'~mK^$(o#'H##=4<Tw#yS3S##|/5"},
{555.37,768.61,4.86,5.69,"###jm#8&2###X@)FF7g7/###wAV7R*0##Md&DR+###VH#r@,###cA'c'6he-b]2*d$twE*&01BVS-)]u#vQ%K]2Lg5RR'r#&.,#%n'1[%REVNQ'(B'3BFI#II#LC)9QS)K8*Qh6YXD###P,#K>#h6;ml&Rf0rG$d56'x/>d(;#$:n'`v'K18D96Y#%$##`K)"},
{322.37,769.62,5.50,3.71,"uP$^_=u[#i6P%p6(S,uG#j*0'g6|m+G>#7z%'Z#://#v#Q?'FC;4##RJJh<:&6QpP#eD:{n'@~.P*>x&-`K3Y.*{`@1l#HW<P>L9##B7Quu#D6Q@,#`,L>c#5'6$##&I$-'5ZP#$##O>#VM:o6Q~G#_SK=v%{=I8,#T~M.c#ip7.,####Z>#:5#K>#pb#_G#"},
{349.40,777.20,4.93,0.46,"+&05##9WA+##3qW,v$=EB^6%,7*=$$bnW.?$I-)###)Z$4&)oR,6?#IcQ'##MpWPQ$0q;F##|B6ku#enWB,#IEB%##<Q&LH$ce/EQ#qYOV5#dnWXm$(95W$#P96<'$jmW.##t5OiP#8Q&+##eQ&XC+%x1<5#Ht4k+D7m)>5#)T-Gr3&K4R##O^5xm+95#L##"},
{349.40,777.20,4.93,3.70,"|k#E##l'4aR+o99;##e~*9=;<R+&##tj4XXAOg9`P#Sc%:&&VZ')##L5JS,$o6U(##sf4L/$D96j,#x7UaQ$H?R<,#*S-.c#Tu%#Q#}i@$##Y8U/##N043c#WM>(##w:Ual$PcP$##in-f>#B5#rS)Rm*###+8UZl#Km)3c#8mP^-&n:UOH%I_;P5#S(9*##"},
{361.30,780.30,4.67,3.62,"###$##T$'`'8}>&a###lJL-'qrD@##p%,gN/X-*7##/6A_92###/##~w.,u#.o2%##OxZa5#iwZ0##cg9Eo#k%0H-#twZHc#(##U7+.[)###8~.C[$?ZM0,#|xZ@,#K:7kZ$]&4+##.{Znl$3Q#vC8gY#RG#Qv)/K(Wd+###%yZ%-%j@,Q6$OU8#.&FzZ4H$"},
{551.08,784.15,4.92,0.69,"{Z#3GJEm(O5$+0&W-SoZ(Cm&md'Be,6_:8?%/h6*c$5c${G#{Y$Jc%Y0)M:9i1:B6(-m$`1Pq$*#l#hQ=Is=CA0]$)PH&lc#*v'###Eg//v'J/S.,#K%+ce'2U6oY#>2SL$''Q$/;1ei>E>#X>$###cA.jY#M.S###>m'uu#>+F&##t1Sel$.,#]5#t0Sdm+"},
{551.08,784.15,4.92,1.65,"&##2%'d0,H?(O%(HR*6$$=~.;e'/m&|?*&-';5#+w#%tG&l#PI+#k2|d)|P$YmPdw,cl#)u=7m&]c&[C%:mP/,#kP#I*/NlPrI.)A$Td)tf-4qP#T3A,#~o)sV)=4J_$#y:<2y%c]6]/#SlPD>####v>#X{;s5$g,&h##-q7x8$jI.f##PJ1%B#GlP6Z#j>%"},
{551.08,784.15,4.92,3.46,"|S0-H&7##xi<9eVJ#$C,$&e$2sCEp$RdVF##a./=6#R1<k>#kH)BQ#%I(qq<geVBv&av*Tw$-C7T2(HdVT,#'D=bZ#}I/3##~@,W~&T@-I5#uQDbvNfR*Bl#>S'uiV3C9kY#+N47v'~#&,##hd,N5#+T,Q@(p[-{H(tq)z029$%pn/7z$HU96W:L5$;##4I'"},
{742.74,788.73,4.01,4.48,"a-#m+LgY####0g*q~2###*##R@,######i##_5%######y5#|w(+kESbBXu$ihRFe/J>#{#$0+G######B##ic'Q,$###<##=Q&;?#$6GbeRhcR###I,#-VP-EA*##N>#=$&0Q%f(&9@,###UJ14##S,$ThRL5$;##;_.(eRyG%9##I:1.w+E>#C1#'dRQG#"},
{310.36,799.61,4.26,4.29,"=m(9l$f5$c%+[T-OU9<l$R&/^YBI?($##GV9z$,###1##EtG,l#JH$>I':OE&p.0C+D<A`?(<NXRZ%T,%=##aiB(c#iY#4$$###7Z5~[)g.//V<i<+jR-7T'GJXgG#ZP#O/#%<@]x*xc'Cl#`5%E208u#%5>b@.-6$-6#=OXOvI###Xm#Gg15_1kG##[&qZ("},
{568.32,800.07,4.30,4.98,"98#K7XD5#iu&$.#(8Xc-%|Z)E.#g8XlA,FT5I>#g'.L9X:R)@w%c'VIw.$##{[DtvW$##>Z$h7)@-QR$%KR*8Z#$N;)A+8?&'&0~[$.M=-Q#(8XgY#.,#v[#<y7###>@$hd)(##],%@o%k'9uc(,##co+7j<?*F######+E)WI-###*##<v$pb####3##88,"},
{341.89,818.29,4.22,1.52,"Dx#ax3.H$###uu#;K*/061,#<6'{l#eo3WD5j5$;d)FJ&kB5}-%1/+3F695#Kg0E^PcQ(I,#GfQ|f1Q[)v{0~H'.f)WwJj80)c$&G4.;;C##x'9*)/+6'{g#IeQ7B0W?(Fh%%H$:F5%p6'##.S/x%&Oy7Q@#L<DG7+###D:#&m'/X<=-(8I$~#&3B'(J/q>#"},
{753.60,817.48,4.70,0.93,"_y*TM<bG$H,#mn$R8O_m+1c#5,#br,$7O-u#Ep0#n*Un**I&m.*f*8=g62l#,#D6GJ&H$BV.u>%Zw*G&OV]1l:;],$Px/UR&M[+`>#z20BK3>8OKu$_H%=f+Se*=-&a8OE7+#L5Vu%d_296'95####n.(B6(m)=RG#Y&/Cc$y6OSG#'0(3$CFZ%F,$2nB;^7"},
{753.60,817.48,4.70,3.83,"3b:?u#8@)IZ&Nv=o>EFf3F5#{v+Do(U?Q[G#.(8~P#qP$5,#?'5fZ#&`?l,%b@QlQ'6@)27%&7(6])U@Q;Z$Ri-3w*rQ)%##:M85/)_..'-$pvMT04o#%E:4,H$?F2=QCTg8T;:]Q%$w)EC+A^6<u#Xw'l{3eh>ub#e,#'7?95#-c#OD(T5L%Q%%##I%&Vh2"},
{542.09,831.62,4.36,3.95,"_P#RG#pZ(%##D>#%##Cq;S>#T5$~>$<['d02(7&;c%N[$mH)zc$S?(t#'###HA1O5#=r@a5#H.Y&##s%&%$;#I)OG#%+-nK41e+j-)8l$>$&604lw,Gg6D5#W3Y;l$hl$D6%wsA###gw?7$'PG#dY=>u$i>$rm&A0YM5$5,#64YMI,.,#,##zKHq,%[)=$##"},
{569.16,835.99,4.56,1.91,"Z##S[JA,$###hw#Iq;h6)###zl%:-&4A,:$&C6&&H%)H%&$%l8*Q_N6NC8,#*+Ys^:o>%B##514,7'505@5#PK0i#&Z>$mP#&05TH#i&Y#S(o%Y#l#qG$%h$,z9HQ%4p/yA/sn0I>#fP#s~*z,'%##L_--'Y1jC>6(c6#=)Yb~&4'74I#:/1[6#bR-)##(-&"},
{361.97,838.41,4.64,3.88,"'s.K96OG#5d$S%)nE;;Z%c',t-)'z99,#~G8qb#q?)&-#9#E-$(###K##UUOS82###HI#:JHrRObm*oG#N|/<R(,42$<8b`<###y$&#d%a5B3Q%>J&m)3Z@-X<3'W>{6)7.(Ed%NQ=Rp1}x0,Z$ML/v>%,m&D>#?8'M$(O{;RG#(c#|Y$AREQG#'##(J+)M6"},
{753.45,845.06,4.35,3.32,"o>#q%-F>####eY#w0(sm,rP#'$(L@#qv+us-A,$7,#Et=%1K7-&Z$'BR'Q-'JB6o[%i%+L6$,BX7-&o,&d;&8&0@9'WAXzU*###*S%T.+yQ(c1<l@$4J-AH#rEXJp4Y>$c?#mR'|6;n{CH5#bZ&<11<Q&?u#;y6Tc#|G%4^%8]2vQ%7S.~]&TZ&/J%L82T5#"},
{783.15,856.60,5.12,1.72,"hc#fR*Q$)[P#:,#uf)qC<95#+D6iz7+6'K,#zn./Z$L#%@,#nF8qx2Dv)###47)Y{+z?R(l#b@RY/-57+#I$pp:7##2I*p>#9DR`n/jH(7d'6K*]g.HBR{$*$AR(m&$w'AM-s83Tc#^m+;##eL(R%.g_'?WA9g%}==M[&,tHCH:^o17m$K[*(]&8x/r5&%##"},
{783.15,856.60,5.12,4.86,"}>&%##&]&-o/6m$@R*YZ:Sf1M[&,tHDp%r4=g_'?WA[C(Gw-^m+<##r83Tc#$w'BM-#AR(m&GBR{$*6K*]g.jH(7d'9DRUe/2I*q>#pp:7##47+$I$a@RY/-z?R(l#47)Y{+Dv)###zO8px2L#%@,#{n./Z$6?'K,#{:6iz7qC<95#:,#uf)Q$)[P#hc#fR*"},
{283.62,858.23,4.16,5.13,"W>$######:##?d*######L|+pb####bR'?qM######Dt9_/10o2######W##BwV0,#,l#y))A/1a>#nPG0P8.,#*##0zV'7)z&3######'##?|VvP$rb#h,#Tt:sm&J|DxP#.,#.##/j6CQ&ZQ'######*##~xV######i%#}#OL>#ko/?&&6#$5##pg5;6&"},
{742.38,862.74,4.84,4.14,"#######9,Y>$+u####xx,:Z$3?&:5#~H$6p1g^1OG#E6#nZ(######Te(4-'%U8+##lQI:Q#x.Z&##}p&^'Htn..,#5O+D;=###.##{6)t#&~07=##[IN6l#S4ZF.-d.)(d%dL/uM7YXX?Z%3,#+%'XZ'RG#<c$-$%LL8Z6)k`+8e..7)@$)t_/-2:<I*xb#"},
{309.91,866.06,4.27,5.26,"5P>######&##w:N*m'Md*N5#Y-(_5$Lq6i,%OG#*##^f+86'BU6######7##f0X}k#%H$6&$~p7(Z$-i3@/*{G%[>#6lFNv';T1######)##t3X:.-QG#V##~DS7&1=.*gd&x7/o>#O0Xx6'a5$#########88PL5$###a##(0XP7(X>$Z$#:I*8_'d+L(##"},
{296.39,868.04,4.46,5.12,"Ru%######'##;wT######5T#vbIf,$V^6Wr.###4##%~N?6'#H%#########`{T.,####B##`{T9d&Rf3VQ#D,$J,#~r8%?%7#$######$##5wT######E%#6xTt>%x~-Zh%.c$xY#~*>Te+rb##########W$G.,####<##`{Tn6*S#%D?#(:7ru$-=Daw("},
{322.42,874.09,4.51,5.16,"X$*######.##(xUfY#>5#IA#b)?v>%:U1B8(+c$S5#6s?<m'm6)######%##}{UDZ&###S##}{U&e,u6)(n%OS0lZ#twU$7'R5$#########'xU.,####l$#9xU~9-,c$N/#bR+ti*ZYMB5#gY##########:[J.,####5##8zUSJ+o,%w$'Iu$g~%)W:M-("},
{870.28,894.80,4.68,1.92,"D>#+##qu%Uc%Mw(mL0Fp7TG#=+_7~,nZ'5-%Af1X-$2(_jl%###2##I.,0$'5^6cR%h'_Ml$5*_|?':o2#-$401wJ&o&_$#####Y>#X7,n>%H?(J##Z&Mzz9F'_|b#Em'(O3j@,g,$S(_&#####4c#hG$G>#95#&##a#$J7,YI,kY#I>#<I(?H&&6$E(9$##"},
{880.42,898.84,4.35,1.77,"OG#&##%c#[c%t,%/I(fY#I#$D2<(?&*##d],($(Q##K+A@L5###2##om*tc'xw+LE2p}L_P#mr~k~.Uu%%m#rD@/o$go~rc%###e##?95h5%0(;nn#ao~H[&)p~,m%R$*)A$iA0EM&hn~)#####fG#xn.{k#2-(8##F{4P`;jtKyY#{5&by*BH&tc#Ko~###"},
{663.10,920.07,4.30,4.88,"894YJ*oR%[GGZvSXH&Mc#{33a+L######I$#t,&ZP####5##>{<OG5>@(q03AxSf~+Ic$SA)`cQ######?##lQ(######=##ye/Ph0v5%SfR0wSf6(0,#Z42'>K######C##(['######&##mS+NyS2,#`x0){SK<B/,#8-&sWC######p>#d5%;,#D>#T5#"},
{776.22,924.33,4.90,4.92,"l~[/c$###4##*/*MA0######_)6DZ&######i82######/##9][_?&W>$.##xA+9d=k84&##Lb[&K1W>$;##hND######5##@][j?&1?&zb#w98z/&2cIfI(}~[O#$h#&m~$lXJ######G##P~[%##lP#VQ$]%/1,#Xf(u7,v~[###oG#Jm%ENB######(##"},
{639.19,927.51,4.26,5.03,"|&Z3l$###Y['cU6-A/B#$K^3V&Z;#$:5#AZ#)^5[P#0,#@u#x&Zt5%;5#aR';16Br.e$)fC7T(ZD@)iY#}u$j)@{k####3##O&Ze>$'##T@&`p7<S'd[&'TU3&ZSl$g,$/)3u_?######?##d&ZtP$###6Q#zM=`M,9d(Gw*$eNz7)Ku$wQ%z(<######'##"},
{629.34,937.01,4.51,4.99,"OC9#########s;^4l$;5#e5$CM8]~/8?%S@*d~~jY#aP#^c#+_<######&##G9^0c$=5#vy-.);tJ-hl%Z+B::^d,%/,#5Z#lL<######$##=:^+m&L>#%w&vC83i/gv&HN>[9^lu%7,#Ln''g7######'##z8^6u#2,#vw&6^6/o'5.({r=b>K4?%0l#<I'"},
{903.95,940.04,4.40,5.70,"f/)Re&v8RkH)J?;].-QQ'###KR%&10jv*######YR#pl'###zm*l),M6Ot7.+8RS#$fQ'YI(p]3Qd&z,'-#####Y0#)]3###(N9#6%8U(v7R}h;Q-'xN2_3A(A)^C3J?('#####rL$w;C###N*.y2AXZ%Lu$5R$-$ArC8gG$ql#&;.=R+$#####0q#$(;###"},
{735.88,948.42,4.87,4.97,"T6)###=5#zP#|$W$##/l#tY#BFr######$##0'5ZP####C5#*n+~P#F>#6,#EB`=#$S5$:,#>Fr######,##J_<{k####mP#Tv).,#######vf`D>####/##lFrD>####3##>q8[H(######C6(#########A,M######$##HGr1u#######8&/c].r5&###"},
{822.07,968.59,4.47,5.02,"$$'#########IJX######$##>KXT5$###;##%.)5~-A5#lY#EH'######(##pJX######Z##/OXwm+###O##Ho.wz:%##yQ%M$*######@##t=L$#####q(#d[SGm'm>%DV#yJXPI+/u#b7%fS2######-$#wQUL,#rc'k_$if5ne&qd,Cn#>JXj>$[P#@~#"},
{114.00,58.69,5.80,5.45,"=z*/o2Qf%Qe-r&H[5$%$&0,#@$($6$-W1RG#C,$###1?6F>#]05cn'@jVAK/F9Yb5#;_6*c#%V<U>#4h5T#$(@+###?O,^5$Tx0OH7d*DXZ%z7Y'e$Zw.-7$5i>h%&rn0.##[J1/,#yu#ZG#4f1=['sx+#^1rZMUQ$87)vR(G^5w?%m?).R%]n/######M##"},
{124.69,73.50,5.74,5.52,"T94Bn'N:Lnp2QpWr5#K3;+l#~<DK5#F(5=u#Kw.$##oT%,Z$L/-~JB9#K>Q%qnWSx${g9~[%b$QHf&uo50##?q9%c#U,$J>#/M:xJ.%q.WM8soWkH%X9.a14.=>{[&SA0X[&eq<_P#0,#/d#>J-W:3MS-_-'}/0V~+-I)[]'[z2pb#$?#J.)Y/.4S+|k#[5#"},
{111.32,81.95,5.75,2.93,"|5&4d%Ir@j#%6T1:@*9f.PT-Hd&Q}H<J,ov(6v#]3FL6&|#&V?'o#%DT*HA10]-[s@w}=I-(3QK;]16J+><3Dc%&Z$*d?Iq5Q95_5$zJ*[-(.K3E#$A^LXc%|1RXn/H(-SA(8p+&=@,0Rk7+d-)PG#sd%5c$2A)FZ&g{-mY#@A%qw023,zH*:@#RcLCe(g07"},
{111.32,81.95,5.75,5.59,"z#RrG$b?${32ZL4Rc$He@U:7UCV8,#]8-xb#ZS04,#0&,4u#IBV4'4~5$K7$9y1kH:(BVYu$zAVE@$q97,Z#{U8).%0&1$##,]-iU8-~%sm+%z79_2JD4:D9gAVPQ$'o+)^1JU3Bd%[7-)$%_>$(Z#py)>81G93SB0aA-?S,uA.+J*%w*7n&Rx+6#$b>#*[&"},
{73.51,84.66,5.23,5.89,"G7Ot.-95#+##2Z$a%Hpb#C5#)$%4wP}>#8B2~q%Gp8u##`m+Wr9|d)7l$j>#b:;#B)|-+vG#ImR*Z$S[#`xIwe*%H$_>7OFG^n/QG#]Q%J8'BNC?5#}B1(v%upRM~.3l#ZZ$i['=44RnR%H%Cp8###1$$z/(F@,###pz'co2QW9|@-:Q#wm+/R'S(4)w*M6("},
{274.80,115.51,5.12,2.15,"###=##jm+###?##*)2@[*###Kn%J94pP$###cP#3,#lQ(######i##y(>###(y-<TAgT^5l#AZ^CkAkd+gl#Iw-T5#$R)(#####9##&r<###0e-4v#AZ^gkDgT^o#%w8+rBE4.,sQ&zG%F5####<##j-*######/##F0*+K3bG$%##Lm&_q7E>#lY#'d'TG#"},
{274.80,115.51,5.12,5.31,"4m'TG#D>#aP#Bd&v(8bG$%##F0*{A3###.##j-*######;##zG%F5#5.,hH&z8+8^ElT^cu$AZ^XbD1e-0v#m_<######8##%R)(##3e-U5#vm+cl#AZ^xOAiT^5l#7'.RfA*2>######h##bH(###cP#3,#{Y$###Ln%'y35R*###@##422vv+######=##"},
{542.19,136.71,5.82,3.48,"s>#tn0Nu$E>#df~0S/###)$#)f~95####|@#xu'5,#Jl$h,%N5$J#$|5&YQ&Bg~.,#.,#=J',i~###%##FI$bm+###d5#Fd)/,#RG#7,#^h~kf5H>#>5#zk~dg~###%##iW'DJ.###C##EH&4,#D7+B7(bFFH81&d'+l#G3:+>AN>##Z$q8+v04###%##1Z#"},
{481.90,139.48,5.47,2.13,"###Q##Nf3Rv)S,#N(03&16#$VA)OJ/qb#^,$+c$;5#0,#6Q$###m##OFH###Iy.WS>>p]6Z$Hu]&j<td+Ed%3~.kP#x-*Ou$###(##`3@###aR,NZ#Hu]JHJmo]8l#Z8*~CIdI-Y>#).+6#####j##vm+######B##c/+ie0|Y$4H%Gu#%p0RG#$6%B,$(##"},
{481.90,139.48,5.47,5.22,"dG$'##RG#s,%#Z#@90.c$Sc%uS*&x0###P,#07+######u##_[+5##V@-~>#~A*((H[o];l#Hu]'mJTI,MQ#{EA######*##o6)j,%#A/lP#~R+?R%Hu]@a=4p]9Z$fT.Ze=w|G######r##%##du$|Y$/,#,u#[,$=/)I/0#f1{k#N,#._/U&3O?(###~##"},
{118.97,142.26,5.58,5.03,"bw,t*:;c%*Z#Lj1t7Q5l$1,#]S*37*L29###8?%MQ&wk==Z$:l$Op'<7QgP#O8QET-i7-(%$G}D,)5HQ&>##o#%er0-^1Ic$+u#$##%H9]u%ia>a~1Kh&Yo1f@(i8Q4##'93(c#N;3###GD2ZP####47#gl&'##+u#@'#YtL&##?c%[##'7Q:5#/,####*[C"},
{132.43,148.95,5.62,0.63,"kc#b?'~D.tv*x)+bc'kv#+u#Wt3B,$######4l#;#$.,####95#]v'Di)B~-Bh4jc$0KGE,$R(S%##f#&5##9w*<#$.,####2,#2)SU94Bc%'z9Y/Ml(<v%%=$S[G#'Q%3'#sv*Zu%###0##1,#58+Mb5rrB8n-Y$*Q,#f.Ba&4A,$###.B#0l#Kl%###&##"},
{818.82,163.52,5.13,5.92,"Yg6Z>$+f('OC7E,_l%t6HN[+*oCGu$K#%###tc%:Q$=?'###=JPTJ+g@-qY#S'3>p',~UVu$h^U&c#A@,.6$4p5O>#oc'5,#1M/TOAM5$F>#UL2qB.#z6Qr98~UiY#i,%JP7:00$m'f#&P>#=-$|R-;,#y?QQW2VQ'uu#?=Exq(U:;'Q$Kd'$e$hg86#$<5#"},
{607.31,223.73,5.95,6.15,"$##ql#bvNw5&MlG#d&C.,=5#AEf######%##JNC######%#####Sd#mmQ###x5O/,#]:/r,%cDf###(##OZ#bbM######Y##3Q#(x)J96###C_:%?&sv%CT-@FfbG$(##B~&LGFX[*###6##2Q$z,$|c(0$'dx./[)8u#;?%@FflQ(H,$I5#.i=m33(U72,#"},
{807.63,236.87,5.20,3.29,"$##TX8gG$D>#A?$a8UjY#/,#^R*d`@/,#zI).Q#?o3T##}p9I>#$+/l@.Yn->D/-FC6Q#:)>P9UFZ&&##9h-}x595#E$#~cL3<D$7#L6({Y9#@*###x$#J9U37U6#$2##dI?@16Aw-6,#Hu#M7U-##j>%MQ#aS/t5%-##:.+&-&$s7$##:m&Vl#M8U###$##"},
{807.63,236.87,5.20,4.82,"D'68Q#k84Q##4vM[o4D>#~$#*H$yg6gm+:,#$##ju#h.U###dR(qm)@sB/,#q1.,0Uq#''##nR?|-U+u#5##:m&>?&w3<$##.,#@6$?/U_P#c;?6D/[EAIQ#)0U:d)###F%#Z[*(J-w#&/##D>#$##3t8P5$G8/TG#k*.CI+'?::3Eh6#@H&OQ#2.U4,#W>$"},
{673.94,240.01,5.19,5.30,"ig&ip[_5%J,#j,4uR.###D5#Em'###$##Fc#+u#######v5#{_52(4AL7II(it[;c%###m>#y}J######1H#,u#######7Q#-o/oZ'.t73bATo[6#$$##YO,-kJ######T-#xG%######P##ZP####}x#dp[<R+###%##E&ANc&######,[#bG$######`5#"},
{409.15,243.80,5.63,2.76,"~#&###-x*)Z$vlST>#-u#.J$4'6xh0_x+ZS,i##;B3Ky%+%-Bv):##BlLO5#8mSF,#>w,*7#_81GO5/qSLl$^'#K$RiS(8Q&l-+y##5{@X,#qnSA,$0?%kJ'sS)H-)&t0M/2#$#vc(=_$Z|GNx3%?#ll'V##FmS95#G##WC,I6%Qu%:&#-C8{##@H'Z&#7x2"},
{409.15,243.80,5.63,5.45,"mZ<###Wm%###3k3###pd&###Qy(>J*f$)###(##^B%z{A###m`6###or,Y?'rUP@,#kI)J[%>%)b20O3B3,#}k#B:(P6)c,#zn-&##jy$D=B7NBF/#ue/P58L?(YU#e)C9##We0,Q#(c$_&#(=4Z7+jY#b$)5c#b;*+x1nQ'jm+i7#qQ)?##sFIA5#.,#z%#"},
{539.11,247.94,5.42,4.99,"xg&-TV.,####ow%jWV7m)/##a9'<TVZP#;##nl$.gN8Q%:5#2m'(UV###3,#tuNg|C###C?#SA0bl&&##&e&=d(O5$Kc%]c$Y6)Y5$$##+H$HSV######N,#9h;]P#&##Fz0H,$;v(4,#j3;N,$~#%F>#RG#8k<pv+.,#+##_a6zRV$##Y@(Bl#7TV%##|-)"},
{256.82,253.19,4.82,4.99,"k;@######,##-~Q[Q'###]H%wfZil&%##XH&d+I###m%/###BNB######*##oeZ.,#Q,#EL--JY.,#E%%,{7<AY###tm+U>#}L<######,##SgZ&##Lu$U#$*IObQ$DW>su%8ARzb#0m(5,#x.0######O5#DgZ$##C,$%##}uNC,#Z*BE5#]wV'##)m'L5#"},
{343.80,264.45,5.63,2.14,"-Z>0,#%K.###K)795#v{6iY####?5#er8######mP#Y>$###W26###VKK/,#^:R95#U],bu$3v($##,LL#Q$###:,#M@*###ZI+###AoFme-q6R###(A'PX0/S/$##-;R(d$###'##Ma;###9/.###T~'<B2&v']#$M?&aj4J#$k[)3>@^5#2,#;5#>eF###"},
{343.80,264.45,5.63,4.77,"rQ)######L,#P.Z###%##v?#'D>pH(c#$L/-'6&P7,%.,@6%0e.######xG#j.Z3?$mP$B$#a(;HpG*A0bG#@k;VvMSH(&##w%1###%##-H$'/ZiH$^./t##*g5;<(Q.Z2##I/Z`-'hv+x##'/1######oG#b6Uo>#w/5g##To47-#T.ZD##O.Z(##d6*>$#"},
{622.92,346.62,5.28,1.15,"wTJdA3###N##GA.vR*D>#`5#B6%a~/.,#*##Q[*|k####C,#S@S2$$%Q%>(#Tx/,-:mA3IH$zTOq^7sb#`H$hw//,####Cm%C+DDb3~#&b&#:r<N|-+94A'-+@Sd>$'l#w3/7e.*c#.,#x>#np+sRPD>#'##O@SPn*95#@/#m{CB5####N~#pb#CZ#.,#-##"},
{622.92,346.62,5.28,5.02,"ZP####f#8/c$je1###C@>|&/>vT###.^/#X6*c$`5$`{TuZ((c$-##w17&##t$*.n$[{T*@)`xT&-%|/0H7'2@+###`{T0x1j>%)##i.+@5#hG$}>#q7A:C8gjGO>#km%L,<)@+%##1y)oe/nP$;,#kG$*#####<##_J,$?&T,%###y5$m7*(c$###XG#z[)"},
{334.36,389.68,5.60,4.67,"|))vWCJ5#^B-UP30Z%###'%$E$(OG####;,#Ic%pb####D5#`+6K/1VG#;:3c<W8#$95#&8(PU7p%,DZ&-##j%-;%-D>#N##;x22##5m%c<WfuQ8##Z/*c<WWn/ll$Ix.dx.Y82Nl$Iu$xK)*9W38)~R+r_);#$&7$2'2Mo/ZP####J,#wx/Hl%^?#d@-v@'"},
{222.31,402.99,6.04,3.16,"5`%'wT[P####(4,7v(4,#B,$;c#I?%.H%+u#>c$0c$iY#dP#f:*=p3,;9UG#_=Zp,&[P#^##`q7N]-Z%-E5#]o+t~0e,%gP#/.,_>#k:UV,B98Z1,#jY#ra-HL:1l#5y+T@'&f0`>$9I%DJ,mP#-Z$[J)m9Zu#'###%###R;+u#)##p5$66%7#$$##aZ%Fu$"},
{242.96,472.96,5.69,3.25,"xq#e5P.,#/,#^i'|,'I5#pb##Z#HH&rP#3u#;l#b5%/,#lG$^|,h'79T.L>#ZOY:Z%<#$i5#dK2:K.N&.4Q%Ig0@I+wY#fZ%pI-Y>#Q^Iy+=ZJY###^G#IF-bf5=5#Bx&u03^R,F>#}>#013Y9*LK2Ff.'+<6(0{v,$##AO,fQ$P6)<##i6(=5####T,#5v("},
{329.39,502.93,4.71,2.49,"jw&8B2V,%<5#92S}>&######S|B&l#o?'I>#QG#]G#@W4,m(*&.x[(p#'*Q#C/S6H&###4$#V}D(~Qrb#&##7Z#^.S4I*TG#YZ'3l#)?$T1SVEDK~(zn0-40gQ'>2SjR-'##.S$+.Spb####95####Q?%I.S###G##=/S8&/###gc#p4AW>$x-*nG$YQ&6Z$"},
{273.56,562.85,5.19,2.43,"M>E_$*0Z$*/*p{P*6'###)##t]1[P#######J>#M#$X,%###8'-9<4_T5M#%/5^^?)###,##-RRPG####-##(##|b#q>%###r96hT)OJ/i(220^|k####;E%uOJqb#3,#0-#2,#>5#M,$###uG%4,#e?$n4^C&3I>#[l%*.9?n+l5%'6&3Q#Sl%###6l#1,#"},
{789.32,570.13,5.75,0.54,"GP@&T'^-S:5#=0(2{8|<EW?%Ry*bm+###g,$x7-D>####'R&#0Se>$q.SiY#VA1$##'/S_z.Ye-.,#a5#/V5zL5tb#7,#AZ$>2Siv'_cP'##%PEwc$]-S5c#R/0fY#=c#Qm'J6(8,#-d'9',-UK,o+)x0YG#a#Ebv&{y5J,$|):=$).,#z,#5:5sH$Yy6Ix("},
{432.68,590.06,5.75,0.42,"J%#MZM%?&.,#m-#_)1?-Q.,#|&-YB)KYHt,%}i-FV<:5#mP#}A)V.QWu%%##]1QFJ/|i64m'W6)&##:O0}>O.m&sb#|,#euH]w.cQ'^7&eL2~-Q?u$G##*_'2x.}u'4y#C1;(##OG#n-#6-Qpb#P[*56#`+C<n+DL8###:f'Nn$AL:%##=5#{c$`?)8,#zY$"},
{432.68,590.06,5.75,2.57,"y104.,i5%{P$7%BQu%;5#%##Y-)>,#WB41,#e,#<-'uD;D>#iS+s~+)y5O>#0pR(v&SG#i,#|o5i`7SA0S>#d##&L6}y/A,$]$*dm&bA.0LR/mRe](005z*0Wv(uqR;/0A5#XN/25K`G#>,#9d(###]H%pnR.,#X##HoRU02rb#|l#GZ:dR,?z9P5$3Q#bI("},
{745.54,591.36,5.28,0.01,"h1=rY#{k#3K#5lODl#[P#]^#wlO&v%D>#v?#4:9RJ%xd-n##5AMj6*###RU#~295-=3H&RT$JmOvS+5H&KI#t&4)7%z?*=?#~D,vS2n>%1,#_9O2FA.,#*##^pOwL9)Q%76#Qh4Q&0EH'H,#q[,@5#Z6)E,#qlO###2,#2@${r3k7/|b#Ic$dHEG-)###&##"},
{815.06,609.15,4.83,3.94,"###WH#d8+Ve0%yX*Z#vQ%t`:byX(##Ia1BB(4l$,##667r_>######o7#{WF,?Q###u##w6D<{X:o+Dw+pS$DH%{0(8wQ$Z$######h##}B7l/4###Z##Tj>-0VyI)c#&0@$yP#w'*f7.###.,#######4.):$)######Uf)i7/)##H6%6T'0,#5l#FU,Y,%"},
{726.08,635.42,6.03,0.45,"5m$_n/&$&4l#_c;*19FQ%1,#q~I:o2)Z#Ed&AN)0J/gZ%<l$<5#ub#TH&Ku${4G$##G@(Ww*9'V###Dm$p&)co+X,%CP49Q%###RG#H>#F'VwQ)$##2l#.*ViaG5u##[&_xE+@)*@)L{2gP#OG####7##8'V0~-Zl%W,#9[DX5$8I(pA/d]2&##_$(I~+@R+"},
{726.08,635.42,6.03,1.49,"vu#PrWU&3F>#SS)3<<Nu${B4Sc#y&3T>#OnW[G#X[*nZ&7GK'R#DpW*I*=c%DO-emW$##TH&oUL){?A5#je-RZ&4-%kw*ae..l#`H&Ql$Ww.B3@@d*N5#$T,coW)Q%;#$Ko$]x2au$@w*?J*fG$<5#RG#+l#e*H&##hY#tR#QnWN5$.,#|8#ZA0&H%mP$~B&"},
{782.85,677.05,5.60,0.54,"t<Uuw.Uw.###iW9=9/usGsG#T/2.,####86$Vc%&H%###0##s7=i^:c&'xP$zk~_YJX,%u>#yaGjY#Y,%`G#]G#=$(B,$###_%-JQ%wy+,U-wf~bG$@,#:|'E}J4,#@Q&V#####]5$JH'###HJ15,#*?#I^/]e~######Cf#f]6eG#7#$=#####R>#Xc&###"},
{385.95,732.10,5.28,6.20,"-Y3n95<c%%##[^-mn-6-((##ep(MH'D>#%##Av'######B,#qxH(7)MU5Ky3Bd?Ht3,n-k,$|1Rx%,A,$J##`19######)##H.-46#`eFy/Rw/3f.%hr?RT+O.R+c#:#$|%$C29+u####,##?[*Bq$b-RS.)fQ(v##,]2`n)t*9qb#@5#?c#}C.A,$######"},
{764.46,735.64,5.10,5.24,"yQ#2^RT,%###m*-i:=.,#&##W[([P#}-'8H&######s%%@L:,U1xt88[S(##c<WV@-PG#2##iHP###%##U,$n>$###4##=m)[@-2,#QpI-kDI7W###(Z#NH8k*H######k6#7c$~P####%##8l$###4e#BRRnl';5#<,#@eEu5&Z>$###AQ#<5#&Z$###$##"},
{399.74,738.19,5.61,6.15,"^0Kz].wQ)>,#mt45^1-H&###L28######B,#v,&.,####Cl#h6G[_,0z5O;7RCSt9-~Z'rA$y@S.,####[##$@(eY####&##cc'@&#g@Sg9/~BSA5##Z$j]%gfH_5%###6##1%%X-*###$##J#%4$#gI))F@O;38l$km$f~/u]MZP#$##pl$rn,(c$###,##"},
{763.63,754.07,5.93,0.44,"<Z%(T/:5#(m%IQ#C}I###$##aH#FI,fY####&##3,#5Z%###(##*yEeaFeY#A1*4pY7?&~?'/tYjMA###cc#s-*nP#7H&I>#%##$I%-;Ne09x[+lw.o8'ZrYToYu5&###wG4em*:-';#$5v&&##iu%6[$X-*aP#>6%Sc#4aCKC44-'.,#^T+[6'u?'i#%G16"},
{763.63,754.07,5.93,5.42,"d,%WS%Q%/%##>f'^NVrJ3###u.AX082,#;5#dm).,#E,#Ov)N5#kH$?%R###U:5o8*^NVJ%,}MVV,%:l#kx(^1;.,####J,#NR(nP$.x)t%,JZ&###;4/{JVUND/,#r>#FiS)@+H,$###DQ#_80hY#7n'cp-y,&G>#'/)e7-#H$*-'RG#zH('##:H&.,#/,#"},
{592.66,766.90,5.56,5.97,"QVRG]0^H%-&.M[&429s6'av)eg6(c$###:,#Of2######&##QVRpTR97([%+QVR#8RZP#f,#GRRpb####N##iT2######(##sv*'9):2:5?B?RRUc%;5#<W*FlEbG$###9##Mo'pb#######.<+sH)_P#uf/p*?######SR%5C.mP$###$##|I$xY$######"},
{337.77,769.33,5.72,0.41,"###P>#=v(G>#}e17##ibKJZ$A7Rl#$TFC[s5Bc%h5#~7R$7(Q##r7.$m(###oT4R$$|5R*##C9R?R%;kKC##,x0jc#56ROZ#AJ(z83|#'>o0XL8f[$|5R'H#)7R>A&1M>I$#~1<Xn#?bL)##sb#dx$Q$)+J.E.(bt4w@/Z,%[W/wIO?6(e5%~U4q&-:$)K,#"},
{529.74,769.85,5.04,2.97,"`.Pty2pP$|C#=S.>$:s$,Y$##S.ni*ml'`I#eu&&##PG#0e#hj.d8Vdf0E5#?lA>ULYl&)$#Nr@vP#95#d%#PR)######A##;~*(d(l&Ftm)#9V_5%U,#>'(&s?95####T##.%'.,#######D>#[#$<t8g80pm,###X,#;a/~H(######lQ#Z?(######G,#"},
{577.73,773.53,5.70,6.10,"m[ANe+2d(2l#q5?u//WQ%/~,Id&^&0Hl$Ku$i.-OG####%##Oj:6u#RN-u-PTIBf-EzH'<@)z/P2z595#W,#AL6ZP####'##+Q%UU#u/P+>A$I)jI%J;;X`77-P#c#:5##C'{&/xY$###&##'m#<G6E;@0l#0D'o97C,$MI(<U.pb####w,$:w&L5$######"},
{349.06,777.52,4.95,0.44,"0J.@,#i{B8,#ggV5$%?1:m6%*R(a6$LeVBQ$3-(###[c&*])F.,nc#~,Q)##*gVNQ$W08K##B94ul#vdV4##H3D'##ZZ'nZ$s%/7v#)}IH,#ieVE7&le18$#/94jp$GdV*##icP)c#r5&'##[Q%b).J..F>#wN/1IN}>&QG#{&.0X83w-6##8q6T7-95#5##"},
{349.06,777.52,4.95,3.74,"B,$B5#}p2,e,B&32,#L~(6dF[?)/,#p;/alIop:'l#Hl$oS%2-(%##?5G5Z$<-S$##Nf1b/%9]32##-/SjQ%y5P+##Q7,f5$&?&_>#|<B%##(/S+##v81&Z#$D<$##42So#%1QO###}e-|Y#=5#wA*LI,###(/SWc#kQ(|Y#<-NLd&D1SbZ%ig7&Q##`;1,#"},
{775.62,780.73,4.98,1.66,"######w5#t5&Z>$###{P#Xu$X>$###TR&(w*###j@%F{:(Q%RH#.:96H%hG$>w+%R'Dd)R>#CsEaP#k[(0U,cG$pl#l)U?R)v}3V&Uj>%%##W@?m82W>$'##lRJL+Cku&f,#Rc%R[=][PhP#2'UZ.*Iv)/&$i&U*Q%(l#a-#}u#ka<sB6$##Y%)f$ASA0)l#"},
{534.63,787.75,5.48,0.65,"###7l#[D6>u$WZ%3x0`B*th8`8-[%//%'j&F`$*T5$%0,/S+######p[%p?*Qq<###wC/o8.owW###E|:533:Z%/@'jzW@H&######,v%bc'/:9$##oq7Hu#yxW0,#&_4).%Nm)e>#~|Ws5%###$##;6'95#H$'^,$gV?###v)TwQ'In.###1T*Nn':xWrb#"},
{352.25,796.00,5.26,4.06,"Ye+.i2z'7*S+95#W#;jU5+C7)c$);3pG#:SNfiCyG$D##O;0*m%#]/$d&<_0N@-Q5#9N2'FA]PN###iJ(UoMnT3###+59q/2######Xu#qW7<[*###D,9h6)#TN###7B,P#$zB2###qJH?l$###$##xG$gu%:5#>,#yp/.,#3<5N5#&q1:5#(8H`P#mp(nH("},
{519.58,818.25,5.45,0.04,"###a>#K#%###C,$.@#d@.'##nC=6?#_I-K)%[J2B%#T+Jci'###e?#'$(###N]3x:&Ro4+##7xY(A'TH('C#'05*i$TwY;@$####&#uR.###+_:GS#TrB@5#y{Yb#A0H&L,#dH'BY0j(=fl$###D$#D07###4Z$d,#pJ2J>#;S(3]1vb#aH#N#%Sd'hG$cd#"},
{326.35,821.82,5.66,2.05,"5Z$io.-C7iP#3Y;9~*T~+^L/{c&Sd%MKLoI,5QPA5#%L3;L&fl%.20,GC.d)bUP@U.}R.)-##'27T'(HOU5#4SPfH'Vf4],#[v*Vz)1SPQc&WRPLw'-&0{G#`@+%i*)w,-##VE8j?(k>%0##+Q%>u#u;+~e0my3q5%QT%'z8n#%9[%A8.3%*ff**],Od*###"},
{417.94,826.70,5.50,5.34,"m*C*l#92/od)WLXiY#rb#kG#`%/###(##Rc$iY#.,#5?#%.+G<;Zu$jWUlc&]NX*Z$>T/Wv$8D=1u#4l#)w(xG$;5##c#^>${NX}b#Ol?rS/(JXD>#`x+@H8Yz6a6)4u#j,$`5$^>$2u#K>#=OX###E5#6?&_J/###*##>~)8R&OG#%##[G#}Q'6#$###&##"},
{558.39,830.99,5.15,2.13,"ZZ#*$BE-)###}U+^g8Cc%###Qm&Q?'ul'?5#yZ'qb#>5#6u#.~+m8*k.WXG#:1W??'l5%2Z#B'5b5%W-'cR*_$*mP$&##M-%CZ&?5#:?@P5JY?S+u#)Q#j2WDS,@A1X-#^<D}5#$'6Yu#'H%#?$`m+<c#F1WvY#wT3Wl#K1W%##d(.F}2,.W:,#:e(cV7&?&"},
{762.96,831.84,5.17,3.78,"fP#L&+X[+}k#I[N*-#9]10$%-;R%D66D:Q[$^-)fx)C8RZ,$bQ'-v&&R';K,s&/:H#F*?zZ(-;R]o0'@*PZ#M@*G&*zGN{Y#.,#AJ#W@-P@)GK36:(pL<*Q#[7RG7*K?''e&.R'x(,n+HA6'###--#z:4*6'/U6Qv%p0/%)3a;@4c$3H#kZ;=5#:H$&W)eEC"},
{708.31,841.76,5.25,3.90,"e05D>#(##+Z$H.W######~Z#iZ(###$##bf/######)##m@,)f0###4:$:=BW/W###bQ#s#$O97###.##BJ+######D##)~.$kH%c#>_NOJHF.W&##C.>4S)[WD###.##rZ&###0,#:,#tQ)>n@z+??0WY#$T1;D[$-0WY,$bz9?#$;5#K>#UG#vG$;5#E>#"},
{721.36,848.90,5.37,0.08,"#Q#jY####.,#)l#?w#=$)G>#hd,8R#s$+':)$7+ER$sx0&g,yG#H@'W>$###mf4j{'I..2##GJX^~(Z$(G=.'7+6C%hLX|U26##IL**6'###={<AL$mf6&##=OX_hVk>%i##GI*./APK6yH####2J$3w-###hu$J[$M$*###'o+>p5###s##Tn/47)###G/#"},
{491.94,890.37,5.19,4.89,"a$+######>,#zsI######/v#nro#######%#4r>#?&###/##X@-######,##Gp_######X,#]sohY####'$#<U7X&0OG#)##rm,######*##JCf######Q##+so3,#6#$a##,T3@@&BA0###a$+######>,#G7Y######[H#oro###E>#.$#t7/&###A*[P#"},
{883.95,917.52,5.62,5.05,"~].Zu$kL)cf3&'V.,#=m%|u$m7.6e'1lL/c#C;:T[)J?(.##>`?9(.#~As^2.%VM6$P)8*8%'^7V>#)(VRR'n<Fvb#ee-[R'+^3xAE5T.)w&}&V)f+ZH(VI#.].:y.6}BKc%`5%fu$,R'Jm'~[*Cm&uV17y,jT4Y7,%##9A$K#$*bAZP#######fZDD>####"},
{825.19,932.11,5.62,4.78,"uZT95#$##m.#M~/A,$e$%Vv&rZT95#qG#/~%<ED######L##E~Tv6+###4$#<'/%g72##>d'b[T;c%)##I-$wWF######1$#VYH6H&###-x#YXE=w'990O35.^T%$%{H(L?%yWD###SG#07$D}Je>$E>#6T#DD=F%&l47m(2E_TO?%Bq-]6&Q[To>#fQ(3d#"},
{688.62,938.46,5.40,4.79,"UwZ######~$#qHRU%,mY#:11%FEPm(H>#0y+Ay7######G###VZBv)###X##d$C@xZ$##N7)>{ZQV>###(R&)1:$#####H##HS)$m(###(##[~SG3E###c##QxZBH'###~%#wy:)##L5$@##+6'###$##tG#MwZ###%##=R#NwZ######$&#_7.###F,$6##"},
{651.59,941.96,5.57,5.00,"@96######$##N9^zG$SG#jR&h=E1n(/Z$-8+^$QZP####G,#OU:######%##J9^{G$P>#>x(/PH%J(&[%/4>ieY~P#$##@-$uB8#########t:^`Q&F>#oc$<YEbp-f,%<U/TnU[>$###cc#R%.#########|;^;%-###%6$,N76fT0,#c~+P@Jj[,###kG#"},
{676.44,948.13,5.51,3.29,"######$##>#$######$##R%-######1##@]4######0##.7,]?)###.##0%)LjF###O5#wiTvl'###6/&29[1c$###Jw#;8[-R*###c##K8[A9[###X,$OKN-`<bP#5j:tD<*9-Xc$=T-~EB######t##v7[vc(###8##2:[1H&####.%,6NyP$4,#0d$,5L"},
{676.44,948.13,5.51,4.80,":$)######1##xIYO5$###PJ#EFFm&.QZ$A9-~FIPu$R>#jd%H82######9##nJY~Z'####'&R{<VD<oY#,4=tGMl,&###T%&g7-######%##(|TVjG###Y5#&LP[JY$##Oc$ZKYxu'###~,#iG$######&##K.*~?)###>,#9JYZ-*###W$#cIY######K$#"},
{813.41,955.67,5.14,4.84,"me1######2##x~^######F%#?eZ8,#Hm%lm%PL;+##1e&o~*)V=######8##?]^OG####O%#%JVp[-H5#RH#K_:k-+}G#S?&8g8######<%#h^^j>%###4&#45B61;)##)w&s|BC7-8H$D8-an02,#:5#0W#]96ac%ZP#p($l]^o,%q5%<g&Fh7Te&BD9SS+"},
{923.11,975.85,5.54,3.46,"tb#######8u#fY#######iw'95####*##}a2,l#[u%fG#gr9Am)######5##*#N###)##Y*.iH)###5[#D(O+u#;5#?H#8DSRd*######%##KCS###)##|l#AAS.,#<$$b:,/GCmq4P7)%0Hub##########1:Q#########,DS.,####ER$Mf/XS,Y%';o,"},
{88.96,77.84,6.77,1.75,"Sc&-R#Vm*8d'$H$Zq&Ny7Oe,UT,{f-Uv'#T.)_6/m&1H$(w(v-+ff#JD@<c#K5LdC%Qq<C6<5.+>B'$=8M8PX}C$.&^%+%H%pQ(--#qdPdP#agPY*=8`9Q?%s@*AyK+A.vl%o4GaK.k>%m##$m'0,#MU-ZP#4c#.d'8r)g7/b>$.x)+f*4&1q6(}D;[>$_>#"},
{88.96,77.84,6.77,5.78,"/|?B`=###;##qc#G[P4##VH'5D'41;8Z#Pc&:R&XG#>]'#[)n95G9/BZ%Pl#XZP~S.>~%E7CB6=&['5Z9<;=T-E^G#Rd'G,$Xo42,#NT))&+2]PsK2_v',w&}C6l|0.~PD[);]P4$%Nd*IZ%^#&$##:'$Yh8W%+#S,TJ$g08-o0qe,He*Dy3<y0%$&mc&L@*"},
{77.42,103.62,6.24,1.73,"6R(bS&Sm(pt@J&.gd(}S)]&1^dM?#$+Q$]l$7d)5##b{6)8/(GEC`3eS-NF:S-(X3/`*>dd+)8PWm'j5%vP#b;9`7-WZ#s7.Rf+mAK|Z$ZZ'8/0p=6J,$,I%ARM>OE###t##(~&K9P#########m?$UC+:~/Vu$(V17%+SH'o5$V9P###%#####IJE######"},
{77.42,103.62,6.24,5.79,"x,Ppb####%.#^]19B.x>%4c#:ZP=w*Xl#L`3x&.w,&*3.=i=YuQ###+##1p%He/@5#a0+H],-xQfm*:6%(w%>y0yN2?6KV?'WuQ###*##c_''$(###5h$*lBhg4B.,[x#)V;h@+wg3Y7*<8/puQVG####T:'<Z%Fn(;e%NU4<##Vd)Q9&8g88H&pu$)p.jS1"},
{150.33,105.32,6.65,6.12,",J@T6'1_9HT/%F35K3WR+Fc$Q)1[&/tl'sY#4[#>8Eii>E>#AV78K,/91uV9goVI&.I[*1o)P'1zZ(&$&EYAkc%%R&+[:Uh9{v*H7'CoV#J.#oVG?%He/uu#c*EuG%'##8A-+?$1Z%Q-#`d*+S*YV1TOI9m#WmV<##*R*b0#Dy6Y>$###I##Jc$W>$######"},
{133.00,147.78,5.74,0.57,"96$gQ&4)0fm)5*,a#&R7%95#:u6,u#######y>$[>$######95#l[(n`,{w-Dz4(m$K]G=5#,VS%##.Q%4##)o.0u#######2,#nVS8:7],%lg92eH4C8Ue$qQSO>#4l$<'#?S.w5&###3##1,#b@)D,8C;>[@-#d(O,#)SA&05pb####aK##Q$`#&###(##"},
{102.68,149.38,6.37,3.39,"###)d%:l#eY#P,#Z^6xb#D>#SZ$zu'e:.S5$######gJEM5$###'|+Um*],%l&(H9Suu$(]1e:S2I+Jc#to-4I+.##3$;x&3n83h.#8)@k7)O94Ec$G$$+9S'@RG[*/##z?=On.Ng/cr<Oc$by3&##er<8,#qF;e@,TG#$['@H%-.KRG#zu%~G#W/S0m'.u#"},
{102.68,149.38,6.37,4.90,"qb#######~COS&1e6*)##U%BKH%O%.IK,hj?E>#X>#@JRIm)tb#YZ$xc(k'/>B-YLR1I+_u#8[=/[Qr6*5##{u%UZ%bnQTG#D>#R,#S(9lY#2T21B(2KR#m$xJR(y3wG$~m#7m'+F:'o/SG#D>####2R&Fu#>Z%.,#:s-d-*mw)YB6M.#$V<F5#8_7&##)):"},
{796.49,155.07,6.59,2.76,"SG#,A&Si<m-+a-)#(2PV6>]0v5$Yf34V.@^79f$cf5<m'PG#t#'gQ$gJR0l#/IR>['Or:JC+#I)(Q$LLR2x-n5#=l$b_4@d*Ku$|k#L?;oG$Tj<>95}=486&87&}IRM4?dQ(w#'yS)&g/=x-y,#M[+d-&###*##]c&(W)L@-/##Bx1YU$IHRB,$S?&sP#BKR"},
{796.49,155.07,6.59,5.60,"6-%pGG{G%i-'r9%pPOm##s'77+3c6*A##UQ',n)###Q$#|5&`V3YZ&yP$`$&}cLzb#de&6PB>oI[l#)a7894OnKJ>#?l#hG$[N?|k#:d%~?%zRO.12Sd*T$$#cIB)*iGM1@&RRON#$B,$>-#wm(.u#8()<'7/D2iK5A?$o]5O{>&I&'R(%}7&)>=u#3l$1~$"},
{734.08,172.52,6.74,3.16,"H$#VSUll'###qx,h<1giD###zVU<Z$3l$###O[*###zb####<.#02>x5$-w,ry4q'/]L91e,=TU(c#o,&BQ$`05;c%&##.,#sI-]Z$R90O;;U07G#$<o$gTU*TURu%B?#y#I{])LJ1<5#(7*%I*Aw&5w,T>5;.-1,#k5#C'I%l#%##>?#HRU?u#mY#37'h(="},
{795.46,195.87,6.27,1.34,"n@(R`>.R%b<D6V3x<8i7-4H%t4IUH%bQ&7|;I,$B,#uT(tOJ2u#PG=,q2OL9I;:11*FGI1Z$}0U$c#$Z$&-$B@,M$#_XBec&)l#%JFbR-(##gD@=^.*$()o$y-U2,#H>#J_$)/1.##61,^^4*?&g?$bI-tZ#(XF=5#95#['#j-U5l$###|%#f%/#T)Nu$GH%"},
{845.19,200.96,7.14,0.01,"zx/-@(o81tH%dD84o1^P#KQ#J,#|~/SA[######&##-B[###$w&(%+U#$^<9y%)>=E*U3Sd)$|:^06}GN+##qb#&##iB[###*##Wv(II*5h5_7.kQ$-B[2-&7A[WG#_kJp?#xY$###hB[3,#l,#BtCzv+3u#Rv)4U.5'7.x#@^9&##Rz;9T#D>####DB['##"},
{845.19,200.96,7.14,3.08,"#f]$##.,#&##||H.]#yS3&##US1tT%#R):K+Wv)Su$<6#q|A2f]'##6#$%##&HP=7#;wXZG#Vf]3-%En.R[$<%,7926##Jw+5f]$##rb#8,#DkK4##u):Yg7SL6dd*$&)MNAz>$f}@r6&e6)we]######1##@lP###)Q#|81_P#L$$*3;gw/[J0vQ%*&,pI*"},
{824.16,227.86,6.07,3.61,"###eR#(J/eY#D,#'m<xR.S#%=ZA0(7###Z$&HFC[P#&##fK0vb#`/#iGP%##*B4Q8%yq?9=8H;<}Y$###fLPt3By.,###&w$x6'dQ#IHP###>KP|G#*q:Uv$**:,[&D>#T/,`>$q),k-+<l##;/qb#aPE###,k8PG#Ws;Fl$CN/G(4nG$0v&qb#bz&i%.W[)"},
{824.16,227.86,6.07,4.72,"@w(B'3[A2(##XN0.$R95#,##3o*#8-LV<%##K#%Z,#''RP>#Z,%J~(@%RL>#+&RlA0aH'6I#cFCJw-d?(Z##hY#~$'GnIzb#A,$$##956]>$8F9f%0_z%E$(@0(&$R&H#OQ'@R)y%RZm*%?%95#mu&F%$qb#*##f-*9z#TB73##to5OA#b#R^?(+o1P,#+-N"},
{788.53,241.12,6.68,3.44,"###W?&Qu#6#$M>#Lo1K>#zb#}G#t/5(##Cl$A##sm,K#$D>#1I&9qRCu$~P#z5ITr>$##_q2uH(je11.#wmRn##A]4dK%~98S8.R^0[H(boR)oR9Q&###cZ;vh;b04@$$s^8$##3)4P=6L@-mYJ6u#OG#=W8_R+'`1.,#511ul%uqR.,#*##1##ppR4l$###"},
{346.71,247.39,6.73,2.31,"~H$zG$-8C###vk5rb#u&+/,#-T+###m|0:#$###-##tJ*.,#_G####.l6###o<7###7R=z#&wKS###gU(JJ)_5%###Lu7.Z$8H$###+{-###|$+###<|,5W;^|FF,#?n&ru9a5%/?$?I>Pu#,L*###[7'###Se*###<0*(-&b/+E,#QI*n?&&Z#OH%JL1$##"},
{346.71,247.39,6.73,4.70,"%-S######a$#BPMpQ(lP#=B+4m(%~+:J0C[&O]5Y,#B&3P##8-SEu#ZP#[$#m+HI0J#.,:Q#p#@<IN^Q(I5#`;Aq,$Uv*d##4-SBZ#<R+I-#}{BaV'%-Sa##x-S~?&ac'.%#E<E1##P6)p##c)CU,#G82{##LrBA$#--SJ6#%-S'##%Q%@%#EC;4##&@+[##"},
{366.33,268.18,6.13,1.74,"PH&>p)K>OqP#hmG#000e.###)d'eu#pIZ######mP#`$+###?[*L,#|IZ###[JZH5#~^:A,#Ke/0##4JZ9,####A,#%~.###Q$*'##FJZQ>#3JZ'##-M<8.$9J01##]JZ>c####T#$s@/###0Z%###kHG,U4`:=###7y1v?9r5&'###KZzl$###Q5#f%0###"},
{386.53,271.61,6.55,1.59,"^P#*e(hG$bt@zm${06mZ&;cHV5#*@%<oT^%/'##s>$J7-###*$'-E+C?S36%EqT|<;ie1H5#wI+Bf%*mT$##$##F?%ed,###{d-h.#*mT0##;mT3?#X1=r-#17,l-#-mTB,####`5#P6)###r@/I$#QmTtc$+mTG##Q:;Rf$;c%T6#4mTG5####W#$_5%###"},
{386.53,271.61,6.55,4.96,"mn0######,c#/&W=c#xB8&##gaFec#_%W$##{%W/,#c5%###p7/######5##l%W4,#.94'##SjF)##*'W/c#e%W###jZ'$H#c$+######Eu#2%W###Q-'Xc$/L9N>#8w?R}AXbMA%'*m&_r/h,&###%##-Q$u&Weu&'##i5$a02&G@#$$Af1qb#e48`P#Rm'"},
{603.73,326.12,6.80,1.51,"Zu#&9V]P#$##1S@+7V95#'##N90D9,'@+'##bg1+w*6#$8##y-*<2<7J%s-)a9V0J.m>$~~%obI_L)B)Az>#$nRLQ%&?&y?#7H&5l$.-#lQ'uxNC_<Q>#wH%|39CCUKQ'&##H:Vo#&.,#(##@5#5l$%##(l#m>#}H*1,#)l#016Jn.=5#<5#c:V_5%###%##"},
{538.04,359.98,6.28,1.79,"kU3CH&SG#PG#sDTv6+######cD(MEE######-2-Nc&######N&2C5#GQ&3l#tAT[5$ZP#Z,#z4<EV7Qu%'##DDT$-'###(##'$(###-l#1G<X?TG5#16&E|,*N@@.'%:3WB._@T:#$0,#6R$#########M]NeY####;,#bZB0Z%###4Q#IcF[L;######~])"},
{554.67,379.22,6.59,0.98,"D%.D#$]P#y7#C~W.,####i(#;TTmc'###9$#E5#:I)######O..P5$3,#)+.9^W(c$###Q(#5`W+u####E##B[)%Z$/,#3,#~c%I].W7+1s5vJT?80RG#q(&n[W0,####&'#WH(/Q$ZP#9##8#$tP#7x(ha=6&2_P#>5#q*-Co2X,%###K%#VG#Z6(95#$##"},
{225.45,473.42,6.15,3.30,"'z#Vv*?##`#&a>#bu%W5#u#&o>$*c$1l#%$&.,####Iu7gu&Uma%S.^P#,##vM5F;6u@+f,%34;qw/$c#O$&P?(###w6%E%,KjagY#a#$6>-XZPY,$x](}k@vU:;#$VZ#HP=vc(A5#XG#el$H,A{v,3##5Z0:o*K~0B##NT-1c#E>#3H#)f1[P#1,#*c#*H%"},
{225.45,473.42,6.15,5.21,"-B495####7)#`?'o6*%##k>#WZ%lH)=5#;5#7d)G5#3c$d#%Pv(Bc$yY$&?#vZ%>t<c~0dP#],B<C8F,$rG#?f1Vc#/S/$##{5$'Z$/Q%1,#Fc%:7%<i6[@+KE?Bi[{,%`n(5:(If[*~.$##Lm(L>#|k#2##cG$kD$DbK?#$gu#]k[~iCIu#J?6{f[8#$;8%"},
{570.81,551.47,5.73,3.79,"pAEKQ'###y5$z9X######lG#@/2###3##Ue-######d##w6+A<?###P&'c<XS9X###3R%K[(#}H###,##qI*######x##J..;93iY#R#:hN:T7X###WD//*8NHS###<,#2S(######`u#H-)J6Q)##A00'9.^2A###aX0U>DV'7-l#F?&Od%%##V>#_?($l#"},
{743.55,632.57,5.91,0.43,"@L)A$Pqb#$##~,7%r>.,#4,#'r&wOJ$##1,#0c6h5%mP$###=h5Qc&*$&=5#XfT(c$F,#Mv%aV+nZ(^A)b?'DR=p>%Cl$###$%,###x>$QN:`dT###2##g</&p5ub#n`1($%^>8a5%^x*1,#nP$.,####0yKCH'sY#J6'lG=|G%'I'FJ-yG$ku$Cu#Xx,dc%"},
{716.72,652.27,6.34,1.58,"rB((=E~P#w$)}T(Ax2zb#$XAqG$1l#yn,J:9O[)###'Q$r8(~)/7uO$##yl&S4Z#f2<5#_H$lU:,l#P.(>^17m)###_G#vF/tU<.,#(##xv'o/ZfY####G.#LcKm#&+Q%*K#01;=5#95#Ay#*L9######ff+T.Z######H^#0#O######^C#lZ(+l#.,#V-#"},
{769.23,669.97,6.39,0.56,"*<&Rn/H.&Gu$-o@e~/O~/###$y-.J,Fp6z,${-+D>####-?#]S-)Z$SF-Il$M@;`08do'yG$Qb[ZL:qb#u>#&B3S5$A,$)##Du$}$':p0gl&8g5x>$5y+{U.6][ZP#'##d:%m&54,#v5&3#####8,#iS&3iAfS2%##I6#e3;@~[.,####N/$*@+&c#6#$,##"},
{361.05,729.86,6.40,1.60,"###[7%ZH(###'##-k5]g3OG#)m#gmPqc%ZQ'gR*cn+/,#t,%yb#1B&>T4###sc&*?:Fh5z>%FwPyf0n6&`uB^A.n,&$p#ZvPWu%Zm#F#Cvl'qu&b9(m{4$wP?yPXS08?#6]OI1+ouP=f%E'695####*e&FH'######aH#g|FC5#A,$'@#=r@M##6~.XA$281"},
{361.05,729.86,6.40,2.95,"dP#lQ$Tr<###LXD>J-K#%T/#n^8>D-zV@e_(>u$N/#:1;uS%pP$K[#*7PlR*)t<CkE,_6~n)nSJ]9PgI-k?#R/36-$,6'zx#######p{0(02Xl%{#&N:P2^3F6P+-'@]*a)/^q>######k6####3,#?%'1Z%(##^u%wp)e&4L-)pb#,m#yY>g,&###$##CI'"},
{361.05,729.86,6.40,5.67,"&##Uc$cG$.,#H,#iX7o#'###hL,(N>:c#7#$S@(D>#y,#`o2%##M(,RQ']>$BI)?*-;?L#[(v]RmK19-'mc$j{9XH(###.##^H(U12e5$|+5;?'rB(aJE&]RN5Jxp52p-S0/j-LaS1###n?#Z?(L5?H>#}L*F,$ju9uo1&~+Q@-k?$3@*H0+>2@PG####$y#"},
{564.47,754.34,6.19,5.62,"M5#];5g,&###9U1}97#d%26%tl'.,#3'&KASXo4###6##pc;gR*-E-{4DT?(#yTZe-hc%%A'dp8/v'Av$rU3Vw.E>#4,#C6%-6'^7%*w>p[S-yTap7f^-pS+EQ?HFHi>$<c#$B2M5$0,#l5#=-%}$>.y2wH)l$),J*;/.MV6GGJp#'###'C&<o2######l##"},
{400.12,769.22,6.49,1.64,"###fQ#I..###.##=g0ll'###'-#%o1|k#<5#mP#xP$>u$$#####oH;tc(###o?&~^L0Z%(##&BHCsC###2##8n)tc&{k####;m#ZqTTQ%#6&gh8eH<L5$yZ$?oTjZ'.,#F%#S%++w)~#&'##MbDmy21w'sbBO;8FH':^#%oTs=7$q;7J#_C;zI#BZMdQ&A,$"},
{400.12,769.22,6.49,5.83,"K@,+-9pE=v*9d[,Mg%q|A.V15.UW#%QG#RV%pM?######)$#9#$.4/BI%hlC{/3k5#s_(w)?<0U###jG#zv%d*<95####/##Vl%H>#,-#S`/m>%###1;%'lKgiC###)-#y[?o/2ZP#.,#5-#m(/######&Z#F6'###?##7e,B#$D>#,##XS+K,$QG#|k#Q5#"},
{378.51,777.94,6.28,0.69,"mJ-Rv($^+5_3yV1x%19?$u,%+^)j{C###Cd&[z3Qu%###`,#Rr:Bc$t'Ob*9rx4fY#{3=]1,<.,bG$###Rm?,/095####%&%p&JG$'q$O1,#m>GdG#{$Ojl#qi=W,%###-[%DQ%'~-###A,#CeKL].3M=###Nr6lB,=#O###=)41s@.,####,l#&%O###/##"},
{378.51,777.94,6.28,1.51,"0##]u>OG####QB0ra>.,#I##xW=@H'###;##@u#y#%+u####<n&wM=B,#-~-fuLJ$'3,#}@$[3A/H&###6$#Y>#~%&E-)###P|?D6(g'$>AO}_5v?*UL$:OGa('oiDLv#=Q&2$#K35381.,#=t4`>O^%%^J/}&#k>OWz&.B5t$#-L53c9:$)uR'$,C0A.O>#"},
{325.33,802.62,5.96,1.95,"'##V(<X.$l':lP#$>?G^2<c%[[@JU1f6)Kv%Ie*T.'P4D+-&d##V@RWB*(y6C$)?<,aAR?c%t@RlS)(x1%$#+03n9)Ro42##i>#tg7-8&VM@?$)LA')KGC.-F@Pu-(g8)zI,c?&@9,gH(mG$|n/fm&,S,|l'ud,;5#Lh%r./v>%{,#|2+hWCM>#>^%K058m("},
{325.33,802.62,5.96,5.39,"=bGf/2PG#s%#uL'HlOhG$0##/o$'`@d%,Wu$;%,wC<hP#[c$?//fR*86%Pd&K')L)?7j9LH&ngOBh=^u%}%'Fd(kg9$##CeM]&4###=R'39(0x1###^pO{@);oOfH)Un,Dx&GU,=lO###YL3a|CM>#f/0`u#4]1Nl$%c8up5u96xG$2/)pD1Xd)0v(###4i6"},
{781.63,816.34,6.22,4.52,"T&%+z:<e(UI,MW*G*C?e)dG$:f)C.*I6EQq79q9sP$x-'*/M|.0|b#Lq5JS*@31gd)YvJ(d'l:TlZ'Hl$A}5+05aH$g_6:oQz,'###2B.6&)8[)cP#>j5z8TI8Tfl&F5#D,6@U38'-Z%/~G#;c%###?##u$(3l$###5##L:30'5OG####C7&7R(<H&.,#%##"},
{756.79,838.12,6.49,3.40,"@,#2p5VG#4l$HZ&YB(AJ0n#$8x28$#;q5q.CiZ(&##J.CsM9R#$Zn-#d%&92;i?X.&)'27d%G:RT(1Z1:)i(Id)[y(-8Rc~(]P#<:+'J-{?*/i?zT*I@,~%$4RL,;/Kw.<~#*02{f(i7/m,#D,$hZ$n_2.g3nJ2BQ%Zm%b/G6@*l'+.N3rA/$M8@I&al%Lw'"},
{777.09,841.02,6.47,1.07,"=,#@s0puOau%XT*^RPs5&%##,~(kn-D>#$##G,#].*eY####IQ%lF3;SQVc%HSQdw,290rR%(17`h/?6(c5##S+.L3OG#;-#r_;#d'oQ?u.,aN<5m(_=48J.Dp2d/-1w+1~'8I+|G$tb#Hp%rAQO5$Mn'A14P'1I-'4VQ(cJyPGY#%&W9/Q#am)R-&.,#j5#"},
{737.73,847.56,6.21,0.21,"0,#`,#oI-|k#H[+Z,#p]2W9)}6+Vw#){<T9+Ac%9~$Ro3Nn(d-).y%+%-(##v.Wa@&Y~,kF03.,P8$Z2W/_2|S0ZQ%e/3dm$N~,3z%3'7###W3W(]KKv)LZ#/91/J@i,Pfv#<R*@c#W?(F^+oP#W%$9x2###.p,x'7bG$B##'NA5x-###N%#'n+-H&###@-$"},
{737.73,847.56,6.21,3.91,".,####[`;iG$&Q%###`/-bM8Q=H###d$#jA1D2@###$##}Z%G-);##ZPI`l$/JWzb#-`2MABf1<Eu#i'Jqj@qFK###B?%jQ%Qm)B7&y'8,Z${NWs81o%-zc$__0;O7#LW+Q$Q(7|P$*8/[G##Z$FQ$Rn,l[(:Y@wl%gA1$l#**=1A)u_>qY#+0+-e)4f28l#"},
{486.01,875.14,7.33,4.90,"^~1######8##J{m######M$#~'dPG#;5#q,#Wu%95#0,#=,#Cg8######M##n{m.,####s$#5CaYS0###>##E6$5::###'##)h;######-##6|mSG#D>#M##8cN7122I*'##l,#k+HlG$B5#1T4######f?$N{m###:5#Fn#7WBN5#,o*Lu$xl'.H%6%%?I)"},
{467.05,885.60,6.37,4.97,"_Q(######&##lJ^######%##k{n######0##ET5###1,#bG#yH*######)##(Di######8##M|n######0##619a5%/,#'##J[+######'##k9b######S##/}n-u####?##D93[x1######5v(#########X#N######$##T}naP#ZP####Rw-Ow(*d(###"},
{565.94,896.91,6.00,5.00,"{L<L5$###&}+h9~mP$###C<'ROA######T,#0$%Sc$=c$W>$M:~>$)###L@#s:~K?(###@-#ivPD>####*##sZ%w5%1?&/,#S8~.,####67#E8~######P]#k}J.,####:Q#67+I>#X5$|P##9~xY$###y?#f6HP6)###Po#w:;###$##56$rm,###/l#Ml#"},
{524.32,899.51,6.26,5.21,"-%*#########b*Y#########J*YdG#>-'.,#>#$*H#cS*wu'{v*######%##s'Y######r##*+Y###9,#z5#58,OG#},#c$+y.0######~##|GPOG####w(#R'YL5$###@(#yB0+u####%##Yp7######?##f'YJZ&###[%#J%KsZ(###~##s:7OG#######"},
{889.41,943.32,6.40,2.99,".z;###@##9g+PH'8c$Zw$~dNl,%]G#_:&=rAd7'FR+M7$V..pbNQG#1l#@g'a//&>By02Vx-k_=}o0i^*&r8E,$VG#T1&/:9cbN###Ju#_@%fm+%6%as1_V<UeN503Zp&_{<K$&3,;IV01%,qbNE>####R6#f&.oJ2pZ#uo29v#BU8I(&?g8?,#sI*Qo$cWD"},
{889.41,943.32,6.40,5.27,"/lK/,#~5${c&:;<EZ$wp+<GGS[OvP#w&-^v(5S.A5#{s9A?&][OgG$],%L5#FPCdZ<wo0Sd(c]OQ'+b%-YQ$)B./x+%B/_>$>~O;l#t5&:c#i<>[~+UL*]V8%3@}@+ZR&iT+K6'j}7.,#*##[X9sb#.,#$##C<0pM=IZ%xb#},$'AIAl$I>#&##86<######"},
{188.03,81.53,7.57,2.04,"i>$7##iZ(###.##Gp/~l&###C$$SB476'###@5#X>#;Q&###.p/+H#_~1###Z@)#@<LRU3l#MWUc=B/[)_l#P.,9u#rZ(###+)RE5#:f1###_6)3v#MWUY*AMSUg?(~x,3RA,/(Hg1a5%$##MWU###qY####7['`,#m~)7e.>,#Q7%C3A3S.a,#>r3gn/{Y$"},
{188.03,81.53,7.57,5.35,"Qh;gG$?,#4f*w{>~/2%##sZ#.~&-~.kJ)N>#dP####oEW.,#JH'1,#$~&P27bT-b7EaAW:m(oEW`i>3d(~l#&f0###g-;TG#-v'###lv(<Q%)R)Ac#oEWvi=m@W%c#{@)-I;|A4###7d%zG#+?&###>5#6##G?'###P-$V/1q,&.,#3##c'0`?)###$##0##"},
{151.73,107.51,7.27,6.13,"7c6)[&%V:P7)$Y2*//W>$aG#hU,hy/][+G5#,/%`TMu(9G,$7T1Ow'+&.YC4g]Ud@*Uu%{I&M92^Z&5w(=a8g6'^?&?,5iGK+I)f~*z~U|I.C~U]l$4&2[H#]bJA,$&##XR(Cc#T,%.$#fn/wR(c<;l&5)f')[U>##wu'<_#xB7OG####U##su%eY#######"},
{816.37,131.49,7.46,5.32,"DU)On/G##,w+[b;eY#:c#,c$7Z%###=/Bsb#######Q)R###J06*6${.)*U1O&X3c#g6'AH&gJ2$##M(.8u#######h*X###B]0T:.8$(DR&H%Xh6%Kl%KA$Ex2{H&GQ%7#####,x%KoN###RJ.l>$X,$d0-[x2UZ$OG#rT&zY$j|,ZP#0#####q41Rd+###"},
{103.89,150.18,6.84,3.48,"###e#$P#$W>$6##Of0BZ$OG#-H$06'q|.|Y$######i57yY$###-{&Mw.fY#x7'eKO7H&;%,syT{d-mZ#b8,_m+)##S#72x/8v(k%#tbNy#%z]4F6%0-'8US/$P%n+%##eG7j7.L;4ir<U#$z6)'##O,E(###G;*d(h#%L?%C$&M,C9#$)-%SG#;QBM6&%6&"},
{103.89,150.18,6.84,4.70,"###/,####EMRm6(r#'###EMR=[&c?)hd&=LR$##jG#s?IK^9###Vc#[d+-^3h~)B/M(~.M-&g6=@+H&d(#Q#]Z%/m%oHRTG####?##FD;+l#o~0<A(#KRmZ$sJRxo4]5$5R#A$'8M6yC<#l#.,####TR''Q$X#%.,#qE.j$+'n']]52.#RC;Nu#V^7>5#aL8"},
{563.28,178.93,7.16,1.81,"'##$fB+u####}Q%4WUE>#;5#0H%{y.{>%Mx-#Z$t5%'##jr8uG$5X.qp:####`5X[*###5,#.&XY>$%##Cz*U?(+.+C5#+>@@5#vb#~>FJ..CM=T,%<,#pf03'XT@-###}?${w+0U7###}P#######o0(#%XQu%###F%#'&XBI*L5$I##K'X#[(6#$###ccG"},
{845.22,201.64,7.13,3.08,"{/`$########i}K.K$>%.###dR-mK%tl'~8+b6)1H%wG#|XC.0`)##OG#$##ClOQ@#N}JcP#hf^*$%{$,A7%7e,WK3+##Fx.-0`$##:5#?,#_tL0###q2-q8%C5S$)Xe)iD>Pc$C+?_$(RQ'l/`######9##=[V###<c#|T2QG#}>#1=@981;%-y,$>/.V7+"},
{208.55,220.87,8.42,5.83,"fd*######(##nmW###'##f8'294D>#P&#eN?#l####V(#n':I06######;##;nW###/##1<*n]495#;A#^GE,Q$F,$XB$pS2Z1;######)##PsWL5$&##56$F}/fS2D%#Y%.Gd#r5&$:#tJ3D.*#########cG4CZ&######(t.O]5$##$##?3-x6+[##K#%"},
{472.00,229.78,7.34,2.94,"<HQ######DW+?o1###G##DLQN5$###t7#VtJwu'###$I##%+ZHQ######h)+vvK###4##G26`z:###61)iXE8m)###HnB%%,6JQ~#&###f0*V101T4###pp2PYCuK5_A1a$%JH'%S'|HQA5#b:9Fl%###R$#|d#[C<eY####qw#]HQL5$###<Q%Xh3L@-$##"},
{413.56,271.60,7.61,3.08,"###'##7H&RG#######K6'-%,######Od&XI-######^l#CZ&^Q($##d5$gH&iDA.##eu#RTKm>%Od&m(()dQ/##4m(%8#AcQ57,###&C&eL;YeQp,#J%+Td&l.-il:tnJun0j$#.dQ)n$?z<zc(M>#RKKZ$*FeQ%##)^2jP#?(8|b#w;,0|B2##4v(dT#AcQ"},
{413.56,271.60,7.61,4.85,"4f3######J##|#S-##y1<CZ#v#S.##Q%S<~'9$SY>$###S?#ZT6######dQ$j$S,Q%Z%(H9,1p6;9/nM-(:R-(91-D###TS&K.-######_Z%aUO3`@*##j5$=7(D(SB##Km'b5%.4<###Hi1QG#######36&wu$kZ($##-Z$F{@p$+###K-#~*F%[(###:3,"},
{604.39,326.06,6.83,1.44,"(-$C/WQG#%##;7=#.W###&##KU0?y0xY$$##VT.+~,D>#'##kv)zL9+e%Cd'/1WxR-;,#9/%<HLWz+L@-F##-.SJ6&I#%)$#SH'cG$t,#I-'yTO,z:(##0d$]Y;/0WeY#'##~0W$H%###7##E5#L5$###'l#DQ#PR,;5#^P#k4CM~0.,#&##k1WL5$###$##"},
{582.85,343.33,8.97,1.18,"SR)4T3:~(A[&fqYN04Fu#rh%[,I+E6(c$[$#d(;Cc%###*$#s[)(c$~>#GH%Y_NutK2,#q-$[pYiC8ZP#Y%#Dg8NH%###D$#dP#>5#+l#BQ&sYK;-(.u#b[#LqYJ#%###y$#R0/3m'###'##.,#$##?I#)bJ6B52u#Hu#oY7L^Xo>%E>#<.#UE?H>####<##"},
{830.69,575.76,8.01,0.31,"Y@$XBKOK6;5#.?$;U1,u#<l$5,#BJ,a>$K#%###R[&TSP###VH@1_5eK5Ih2RRPRZ&+Z#Sq2;l$wu#_;3pf4###b,#>*?###Hz<nP#h-)rk=b>KKB.Sd)lC23?%P@Ap~/`>$$##`J-Do1###rUPlm(%6&=B*2S+6}2jQPT90B/#bRKVZ'###O##FM=OG####"},
{830.69,575.76,8.01,4.43,"XqQil&4H$Sh4Y|6jx+(wI$K3TQ$.A&LmGS6)q$?}^7%z3i04(PE,u####`t6VnQZv'=&+rK/R5$*m#TZ:6{??EA^G#/A)X0O?05L5$###Hn#]YBj^8N,$26$N?&|11,K0?Z%6g7.##d>$k?%d8&l2C###,##dn%C^9######;8-)v'###T5#eA3######b-#"},
{420.01,627.85,7.94,1.56,"P|EWd(T5$=N+$D>3v%}k#q^#s6)|Z'6Z%U,#Q$*######c6$9X/,<@`G#9X/*~'#Q8Qu%5[#?wHDL3D>#/##+;;&l####H##/_PnK4lZ'wK(paA5C(Y]5hQ#0mQs5$L5$q%#6q:fP#D>#_##$C90,#kl$9X/;n.=##xx+Lt?pNF^c$[P#4:(gm)>=AZP#G##"},
{420.01,627.85,7.94,5.52,"e,%kK*Qp6###JR&Ii6p|A+u#o?'~H%4zOMZMa6K2e*'w)/@DK-)@,#gH@6c#-g1;u#*?:Gc%PvP$##V.&TdLI?(Id$8wP5bBeG$###w0,,H%~l&###y;+a(<7#O$##[v%tD5g6*3,#JyP>,#vl'###`5#Oc%D>####p8(pv+Cl$PG#pd'=d(wG$QG#8M5sb#"},
{607.09,636.82,7.77,0.83,"9(#)mR/c$###t8&_W?uEE[G#h:5w#&dn(]n-5u####.:#qq>P&#;;?IU$X83/i(UmR0.%*n,uqR@^74u#u,${7.###TZ#b6)$##.,#xL$HJ1EI*VZ'6q$'>GDmRmP$;##C1,)1:###)##GH&###1,#$9$7m)T,%###S~$i'7]+L.,#0##lo)7&1&?&###4##"},
{763.41,644.28,7.78,6.08,"pf3UG#Wu$V&(YVYeY#K>#$x$[XYpb####>###D<Wv'{b#|$'ll'2,#@[$eWY{2@H,#.q9d^*bTYJZ$vG%w##J_<G=6Uu%sG####^$&Cf.6A/Oh=$@%3*De?&3UYYZ%j,&gm#b[*bm%hl&C%+###^G#C23W?({u'X>#~v'dQ@T;?/,#:5#Ea+D?'OG####U?$"},
{583.93,662.43,7.86,1.14,"8$#]$S9m#Z{B2)S1FDpP#Z'3G'S.,####Sl$zP$1,#G>#iQ'###ZP#lS#x#S+}FOG#Ad#?dK3'SOG####9-$rm)av(95#0,####>5#pe#4|EBH'<5#]?#C%S}iC@e,$##GU,~Q&)j@###%#####[e.)Z#Nc&###e].El$NI,%c#^BJ=R+Tu$rM-G&SZP#<,#"},
{242.16,729.01,9.58,5.97,"95##########$eX######&##N4s######%##0x0######&##95##########)Uc######*##~4s######'##SK4ZP####%##95##########Qo_######6##V4sD>####,##qo3`5%###*##OG##########<ZQ######:##@4s+l#4u#:,#6].G.,sl%vb#"},
{403.71,738.92,7.43,5.98,"LfFU16a.-dQ'mv=?o3###$##u{;######<##m$*######=##oC6]e*FfGi5KiLTI&3=Z%k~#oIT######n##(d'QG#.,#0##ZZ'|g$m5HAC63KT%##Wc$|J(L:TOG####;##P[$|,'.,####(c$I,#/e#7}B;;=###g$#~5DE|A.,####:7$Jm%0Z%###%##"},
{771.94,795.94,7.96,4.64,"4%%Er<4F8]@->A)NU7@h3t?(ZF@w?(Ml%,7'5m&m-@#.,.##sp&R;?-o,{G%|wNQ6(3T,p[E8/1#?#,Z@88M-m'6;*B1Po_9c@+=#$522:C4J/P9?&Tl$sq.s&0W.)CV6'R(%c#E>#,i*P-P(c$###c5#X91uT5qb#%##a~(/7(;v(95#$##il#>u$%##7#$"},
{538.71,808.81,7.10,3.58,"$-$|Q(>u$)##;?'Tu#y7+{(6@T41,#jx&w&3(3D###(##&?#op/JuKN#%]G#sKYR?(qe*G`.,K4I,#ZOY1p/b:=$##MQ%}Q$4q7+J-c7'Q$LdLY&@&,w,KJ'$q6nL&CJY[G#^K2W,$v#'I5#p5%4c#3]-(L7hP?:{8)m'yZ&H/)bAEjS1D,$G:.NZ&OG#%##"},
{381.19,847.68,7.65,5.03,"g&W+q+Jm*0##YQByT+Kd*1Z#%^1lS-D0/a[(Mx1tb#>u#N_4za7]SQ;c%###5{^OK0;l$sv$^*EsZ'w6%bo+'>M/,#3,#d7&yz8A@,######2z^oP$xb#tG#FnY'##zK3Ad%FtK$##%l#Bd##p4#########$z^###;5#'##xRX###Fv&;/+_e0###3,#?7%"},
{402.09,852.16,7.19,4.99,"x:2KuHL5$###6jjfm)oY#5%%k08Ul%S$%=/*X{B###&##.%%CbG'm(######AijOG#)Z$1H#f}L*##~y4x?%kXJ$##QG#:R#6bI#########Rij###wb#F5#BlP###x[)Vy.+95###1,#*A%SD?#########6ij######*###RU######Gv%hZ'######wP#"},
{370.04,865.80,7.87,4.93,"B&]######c##O']JU+T6)|##l3?t'+Xw.*$$n@,Y6'DM6M.)6@Lj>%###%##$UMBUVQu%&##x(]2L1^,%tm$fg6,e+nZ%Cw'II(W>$######?mLST5###+##M&]?c%DZ%6S$xL=)##p.,uR&HZ&#########^lQ######*###&]###x>%gc#v]7&##N[)Ue("},
{915.56,897.89,7.83,3.07,"`$+######$##<Dk######:-#q=L|G$de''U0Zy5;d'rA&Kr>@S0######&##=Dk######v##Y2i).+O#$Z7)q$&OXF4x%pf5<]4######%##ADk######R##~e[;Q%AH%VH%@H$$o0&e&v&4cA3######&##<Dk######m##'PL###%##Pv$pb####4###@)"},
{670.12,135.74,8.97,1.91,";03N5#bQ(###@5#|-%:FBlG$@T2tu%wv(zA*;Q&P,#95#vH&*oQf,#M~0###q6'8+00@U&c#FDUdC3kQ(Uc#cd**7%95#$##DDU=,#lS1###OH'tc#5EUtg6j6SbZ%sR(|F8oc'po+D>#(##5EU<5#o>%###<?&:,#.'0zc(D>#r6#'f0m$*###3F.r5&###"},
{670.12,135.74,8.97,5.23,"e09$#####(M#E93sd+.,#M##Vx,nZ(>8+F>#Q5$###qrU95#}>&5,#%$'yA&1A)ft:kQP9Q$qrUZp5qZ'JZ#UJ0###qrUK>#3l$###tQ'JI$)d(3Z#qrU.U1%nUrY#kd'wj/:T4###)24S,#mP$Q>#Q5$d##c5$+%&bU7=Q%?O;qZ&v,&V6%Tw,###Hl$J5#"},
{292.57,197.69,8.88,1.62,"Bw-MZ%###%##S^e;5#pb#.$#[kMA##FD>O?%5?'[,#O'/85I|T7Tu$###(##c^eZG#W>$@$#pvV5$$o99Q##t~2R7#Py8-Q$<D>^5$.,#)##}_eyZ(95#G##{HN101ll'/##'f2T-#xPP-##i:=^G#.,#,##}^e1Z%###J##B%VJH'Qu%M##Zw/n%#S5P-##"},
{184.10,213.54,8.90,6.02,"dQ&###$##^[)Xg7######gZ#v(?###Z##M[E9#$###*&#;-S4I&#########'&W######9##7&~###4##&i/%m(###a%#@$Ou?&#########@'~######.##l+~_5%+##wK+;y*{v,>%#9OD*d%#########|16#########PgGRd+######k2%e*H.##A,$"},
{394.95,212.28,9.85,1.86,"L8-iY#L5$$##Mr3'kBPG#T6&2,#y`0u[-.M16##CR)OJ,qs=5v'_5$eY#/##$>JS?(O>#Db2nl'6[%`+FcAIJK,1*:3,MK-&.7*^l$###$##X9UQG#Ju$X?$%r=)##%9UUZ%/W?wY#m8U~H%^#%KH%######|9U&l#ZP####)#I###|}E8l$:$)'##:+;95C"},
{814.87,244.00,9.71,1.27,"6Q$9y-4(5sx4yU)Z@-z@$Q5OS(095#4/#&dSTG####G'#hcSmdSA[(eR)%y0~T/ZZ#'zR_06PhS.c#@I*,d%0I)h#&C##UH(*dSe5$%v%{y%0%-qc#=oE'z4FFH$c#)S(Y:.K6(`H(###*##Q$K</01,#DR#^Z&zz05S,^Q$1Z%U5#dc&|Q$iY#$l#D>#%##"},
{814.87,244.00,9.71,4.21,".,#$##UG#K#%Hu$T6%b5%nY#S[*j-&lY#lB.QG#w6%0V3X^7###0##O6($[)27'A;,XNDG#$vBK&'/hu&aI%.[&JM&cYMWQ&t-#*J/a7-|Y$t@,7@&UTS8Q#dUSZ//R80H7#7d'?p(zQSY-'@(#TQS[>$###'_#kQSN|9.,#NU(zy90^*hm)*/-RR*R.*'_2"},
{637.50,274.31,9.04,6.14,"5/0LZ#$v'@l#PxURc%fY#U,#,h1BIP%Q%$##eL'ax4######V%.V>#%l#;v']wUw>$Bd'A6%TB0^;8Nl86q9AULwm,)H#yq2(7*{l$P?&:7+uwUqG$lG#k7(zXJ###WJ#E@Rg$+###N##&fJI.-1,#nu#bvM*#N###>,#8I?2sC######uu#)-%D>####A#$"},
{618.40,304.82,9.85,6.23,"J>#.n'F6(?Z%>+EjP#:#$.v&U]W<7(5y-#.*cp.z'6,1,>kDW>$7,#Z5$1?%N{<Bd%Sd(&02_^W>#$<H#C'/lp8###~-#+BV&95/,#XG#L&*c2?aP#|l$I`WE~WD>#)##xSCLr<OG####P,#&.,###&##B(SAg4ZP#)##?]PZ04OG#%##s^2_I*ZP####;c#"},
{588.55,364.23,8.89,1.14,"s0Zcm+8l#7U%>fT(`<###^$#^uL1m(###h$#]l&r5%###1##S6<zS3$##,H#53ZTp8###i$#:~YbZ&###A%#[P#qR*###%##hL:O5$lY#Q%%t0ZOG####%&#+SJ.v'###,##cP#/v%+u####6I+###m5#:u>L0Z95#3,#x%$?JTF>####V##r>%Gu#6#$###"},
{545.28,400.68,9.92,1.31,"<<Zfd,###d##_=Zl6)###+##1H$(I&0Z%######CH$-$(###h9ZfY#6,#6j,_=ZA,$###g##cC8L?%(c$'#####[H$&v'###~81qb#6@#kDYL8Zqb#*##4|(q83+7'eY#?#####q[%+u####A,$%##}Q#A[Om%0eu%%##AU*8l$aB1D>#*##2##L&.{k####"},
{704.47,520.66,8.94,5.12,"1$$F7SeY#$##M{/0C8-Q%H5#+8'?mHq./4,#F}6*8095#%##Kp1>}E######qpYBm)+A*>I%zp8t$'CMN-C1LoY#Z$]#$9C*5)?jY####iZ&<pY.c$2Q$AQ#0;8.I)G;.(01#J/###N,#`^,-%-###%##u6K^nYZP####`D.dK4)Q%.,#yZ#'@+######sc#"},
{284.33,541.23,10.06,2.75,":J*FM7.^.En+2=1E_;86%$$'?i*W,%.u####gP#0,#R?'###MR'u?'}N/)^6VU.C7(1fKOp6#_PVl%lu&I[$B%,###AZ%*##RB,79NHS+cu%.i96;6pN@O38$=E1,#|,%zG93m(###GQ%Hv'Rc%;m&Z`,lsGG%CSZ&ep/]~.VnHOG#D,#^$'0u####@-&|c("},
{284.33,541.23,10.06,3.51,">##?=D<5#|k#'$#VaHQ>#OG#)$#Bv)E,#mP$6,#D>#LT%yY$m$%{wO(J,Cf0GC(huOJw&2e+1'K5T4bG#ZQ$f-){k#6##B['vl&HH&k.(HwO*%*~~'8$=v=CgxO4^7[?%E0(eo)qx5###Z-'@8&;vO/d(lv*nH%u`0;=Er81Af&siBGR%Bs>-k36;?4m#8R("},
{367.55,557.69,9.63,3.94,"mw(^v*###fP#j>$wy1X&4)##B#$T=4z|Et6)rl&y94BH$7(1;7'V6M'Q%.l#N-&gAN+%*TR+DFC@'4JT0S+@SQ'4R)fm)38Ft-*.o.%f.8F5N@,ld')W2pBNWHJjN?Zo-A'1m6'QANtu&=$'yK60I*6l#'X7},'xK.|>N+&)]-'5M4ctG0x+Fm(tn-0R%j=="},
{367.55,557.69,9.63,5.78,"tO9v[-T,$yv'K1+T$)9)4J%,2@)Km'UJGkl&8O28U84d(s#%,:/Gm)Gp,#M8B=78HG(f/IR(*r<R3=r2;B):I7+2c$2&+S~MgY#E>#`/'1ZM[$*8$&'^M_80'^M}]6Vn)*],h_*6^8g6':C75##~?&$6&6#$0,#cm$?k?eY#J7'ov*q27wS/i_2v81c?(n@)"},
{451.45,571.51,9.19,3.97,"###IR)HH%`Q(%##gT3o04681c,%(9.VFBaS-LR+lv)&$%|O7yb#`ALoP$5u#?7({7NFd'2%,(vKqw/+n(,t<aI-I@*#~-|O5gv*>&.*J*cs45.,$d%~1/A9NP6Mc<@^@)b]/^H&I8NCQ&@Q%>u$###7d#&8Nj>%^H#U7NjB.NH&:)7P<?~/+_c&NJ-(m$DdF"},
{451.45,571.51,9.19,5.74,"e=7jo5}b#[n,/K'U7-M:3IA0G[)ll%^8Ij5%Hk6i/4s?({G$';1.w,#S*i'2W<5$wL#A.A6%|'6[=Cz;:Gx.fw,,Q%gJ+*xL%J.$##0S*3*;gv*&%(exL0.)yxL9f2'R&4&+d|2~./.##E*=DJ0###)6%xH'{Y$/,#BG8#l#DJ(aR,uh2>~,awLtc(%##n%&"},
{352.73,665.89,9.20,2.21,"^p-xP$ZP####Bu#O]+De/###D##p%,*~TeY#VA-jZ$AdMv91obDQc#k6*x,$wn)''H&~TuG$<`T?cHwT1T7)~$*###xJ)|QLJ7-(##iR*x'/'[(/m$<`T{jCw[T;c$k~(:'I>R+%##,v%m-';-%;c%9c#NQ&$#####QA'[R,ZP#L>#XZ$s//95#0Z#gu&/##"},
{352.73,665.89,9.20,5.28,"ZZ'9,#.,#/Z#)?$Jf.OG#dP#nS(QI,######}>$h5%YH$>u$+m%[[(+R*$##he(goG7.T<c$[2T$kCiH()d$L%+6L0Q6)(##k^.~6MQv)3,#Ny2TR'[2T.-J=.TwG$j~*CBHU$*f?%^YG>Z#i.TJJ,Yw+&I$$IS.,#D##2o,s$,###~u$VT,.,####XD4T5$"},
{546.78,707.09,9.78,5.21,",##h~.######z#$JC7###$##MK/D$).,#g?#e@.I#$###B@#Q,#_(7$?&###af*`]G.>IDc$Y`Unz9LQ&:R$pf5wn$Hd*S##vi<IX>n.,g.&uf3Yv%?&Dd^U{[UI#$M[&i8Eu@/x,#N{2RS,5HLRQ$TK0*'(<_9###;##'W56U9###'##CC*+@+###26#f.,"},
{370.55,804.26,8.24,5.41,"AU3uR-W0,k-(ZTPZ[,)c#KJ'@~)g|D#c#G;06u#sN4(E9H7*T]1&-&,r)Q(6qRPfl&sQ$[3/0S+UT4i?$<t<T,%R@(KV+-&.5>6Dh8C6&}u%#SPlZ&`5$Z%(}>M###=H#28+W>$###CK%jg7rUPqd*QG#?Z$+SP###{>#*n)GRP###)##Lu#Qc&###/##N-("},
{914.50,935.06,9.84,3.18,"8Q&######;##D{A######l59z,'M>#1[#eM`3-'A#$LJ%-tFYd+######%##_K`$##1,#RS%&|CTR&N/)t+AOS/*R%7h+o=Ha-*#########uL`###%##xP#&L`y$+@[$`*A@m'lx07^)~WCHZ&#########%M`#########MT~Wm)3Q$-6&oQ&C~.$6$OJ/"},
{818.52,213.71,11.08,1.37,"###%R')''}%]/##E6%^M.&&]U-%-c$Q7#F&]$Q$###B%#1&]nP#S-GD7*eu&b80#'/p93fy68_.R#%OU0*V;U*9F>#Nu#0z:Tm(mr307*W#%0']1R')R't.';T38,#(w>OB183A#l#$$%F.(9x13%(jR,Kc#d']rQ(@5#&n#<w+@o+R/,V-'Hl%G5#5H%yu$"},
{818.52,213.71,11.08,3.52,"Qf-^Q#el&/,#-R(+E3ul'Y>#;U6um+^>#h~+M>#3@)rn+Ov)@/Et.&ie1###x::if(o84if&DXDN>?|b#wJ&mY#[&GU80O#%]_Rj#%*I*tb#S==7@(+S,Fo-@d'@t3*x+GJ.qb##)5Ce%uf5]_R###%##aP#*_R$##k?$R^2~-)2,#n()AC8&6&iP#J&*cz;"},
{451.31,222.04,10.48,1.61,"HH'.,#S$#c$Tv6+###B$#G'T.,####9z&^$T######Mi)=$T_/3XG#OG#[u#+%T7H&B,#L(0*6&Kx(kPFYjD####$#`vMM$*gd*uY#pb#.##bl?w^:###~%*IQ$)~C9~.s1:'##E$$jU3z)BE~/1c$###]##<;?B,$x5%,v:$c#_n(/;;jtG#I%W^2:04a7-"},
{451.31,222.04,10.48,6.03,"PG#/,#`q(&q;G,$###?d%&w+mP#D>#c#$pv+P5$~##F%)&R)L6$+0./8Orl'Xg9:##P(7Pb<j,&:##@D*9dQ###D%#SeQBv)Op/n7*(|DM?$+UM2A-3y6cZ$PJ)zr5@=D+$'###+q#}cQeY#,.&?W9wR.3Z#m%$:(5Jm*###`-#IgP0e.###$##Ji*AcQ###"},
{639.64,277.62,10.05,6.18,"]07>?#{$,]>#,pUn6)|k#4##Pq-MnUb5%###r_+T7.######uA2mG#fG$>?&4oU>H%.[&c6'nU3i(8xV.VHLYySbc'1##n21m$*<l#DH%z]3hnUaP#Z>#Bz.q,Q###6$#q%Q(d(###(##t<71&2###2##1Z@oiD###&##pZ:RNA######{,#(m$OG####'##"},
{243.79,402.33,11.83,2.51,"nw,Cv'.Y8'ODwD6N@Rop0]H'3kC]@+ao1.HB.T3###c>#Nh2###%##VH;9~/.:4Hc$9DR,e+^CRM,$Dm(I&)v%/9,#1$&u-*###$##|3Apb#4Z%4,#9DRM@R`sE0,#o-&GTH5J,=?#Ow..c#######WA,W>$###(##/~(GS0D>#'##Zv'17+Y>$G%#C)An>$"},
{243.79,402.33,11.83,5.37,"mf6xY#e5%L1#O$'2@*rP$@##~7*zH*###gG#/8-A,$###6##4$((Z#@e-D-#1I&]/FTjF7,#uoHf-QzY$f##bN@(c$###L##,m#7C5[&4%##0e+l(.;~S`#$z_SPv'-w+ym#G>G+u####M##^%,#mBZ08s>#jx0gp3:>A3`;GJIO@+%U/Qz8qc:hx4+u#?##"},
{329.42,417.92,11.40,5.44,"]]0O$&PYHc,$S(QM#%]P#{G#Kp/iZ(###'##7%(JZ&:5#$###&.4v&k245#JY-ASL8tg4&C5+UL3V=,u#<$#2o/}G%&c#$%&~=39@+e8.g~/>x.t7+T&QT-(S(Qjl&Am(q[#tz2A,$###pc#t0.P(:XJ/n?(g&/&p2d2;.p2;Q<,m'GA+E)<2G5$q:dG$G5#"},
{418.56,429.82,11.12,5.40,"W`8>?%%D:DQ$=pK3l$$##MQ#nA+|,'###L5#ee,)c$###.##RJ-d5$+z-&[O^%CJW?=T1_f03KHk(>95#$-#XL8eG$###c6#8Y;hl&x05Ow-Kf0x@*yJQ[m'(MQhl&Pu$Fn#h{=/,#E>#7-#B)58L7fR,L?%*A+(e)a=@h05YKQrG$yZ&fN3taF###I>#e?#"},
{251.89,446.06,10.41,2.30,"$7J[/3<c#86&~v$#i9-nB>c%|bL&v#5F@qd&mjH###:##NJ*T8+k_;(G>oH(ovI83?(&.V'2TK64Q#D&.a^Oc'9$##X>#<8'jl&KQ#_^OPu$1~OX#$I%,P~'j3F###36#:L4nP$N##~$(hR+[P#Q##I$F`y7s`DE5#-I'tH>+C74,#zb#9H#/,#?&#('6pb#"},
{264.44,474.79,10.61,2.57,"o.,|m*c9*Yr>1^2SaA&0+|?)'s@Ll$~6(9#;hd,###E,#.p,B0/.'5*;(81:{o1J%*LhQd/2kfQRG#j6)?~&Xe/$##vu$S6(EU.&eQaw+)H%X|>[_9h`?C+>YWC###vl%:AD:.,0##1I*eH&/@)BQ%-Z8+I*;hQvP$Bd(Pm&T{<###L##He,D>#V##ax3|c'"},
{264.44,474.79,10.61,3.50,"ZU5}H(iA)tn,$RI(n-###(I$H-$=6(4##Dv&eP#pb#0q%Ll%nd&|'.3g4|K4R{+LGN{H#j955wBvA46##6v$N.,(c$L##b?%,[(*$&kQ&|KNV$(0K(U=4TINOJNO06rH%U'*']):K5###jH%ZS'BHN%Z$Av%eZ$R+5d(;_R+Tf&r+ID$&903-j/>NC/Q#%m&"},
{348.33,489.64,10.57,2.17,"Aq:/u#S_4EI*W6&`//a]NBH&J3;6g4IR(58*9e-Mv%wo.g'6b)2UuM7S,V>#3@'x_0G]NMQ&1[NlZ%Sn,e%)^f3{Z'[HHQd)Q):Z80T26,[%5;<EI'KM/RYED{:M94lW6yy6m&+XbF6f0k,%rh=~$$Z3@@m'M[N9[*3%&p$([l$|m(,yKi5%4i?x#&4A,Ce("},
{348.33,489.64,10.57,4.17,"SR+TZ%wH(#30d5%Dv%hW<1{3&.+]//%_7U0*J.-7$'#w)]QA0`:C):2H%2S)K$'jfKGS..7)>vMZ~-0.*aV+.K3:C0wT6kS%tH(BC3|sC6J*y~1cS,o$(dBM>|?b:67g3b2<z?(*8Fg5%/S+x/3Ym(mn,1V0Y7.ov'Fx0v2-#H%`n))wM,/,Qx3aZ&Bc$nF5"},
{348.33,489.64,10.57,5.57,"uO<OH&=B3@Q$li.V..xI*cA1Q6'fv(dSMnu&xwF_H(MI*%Q#0/+--&k//kaDri/1}E+K0}w.Fz54YFsV;cS-of-au%JB.36Mha9#Z$3w)r05XR*Tm'VTM~d*vTMz6*S@*B/+^z-}Y$Xx-83BIh2W&2-R(i-)rd)]17_`>ud*Qz0:-'To.?E?zN5Y964v'2?%"},
{696.78,492.59,11.21,5.11,"dx)Je.{k#3,#l&'-'6)##sP$O@%FI,$##N#$H6(######A##^Y7[mSSG#G5#d35BeMfL5&-%7rS5V;iu%x#$kg7######<##1oSg7.5,#Lc#5nS9H%zk56D5XmS)l#R7&855tC<######K##:mS######xE5NnS(Q%t5$Q%&1A0D>#H5#ln&L[+######j5#"},
{435.52,506.90,10.73,5.51,"R124-&1:6?Q$8[?k,&]c%g.-f?&b?'a~N}Q),-FIH&zn-IH$p[*Pl$gB-oZN6P=;(9o]2ap4YW8:<@,{9(A+<00<H%a:36lE@s1?-(R%+Oo2Q%,&%)[~N/d(D^N$6&7%*(K,'r4OG#k##~[LRK.$C6&w+%['z%-K'4Wq86A-lj:T5$l@*_O<*[M.,####:%$"},
{269.88,515.60,10.44,2.42,"nO?U06_n*Kl$el%pK/jAJ=Z%V&Q~#$J.,Wd%1V:1,#yY$:##8/,ujD2E077+IOC2{9Z&0)39_h>>,#R6'3LMAS0/##q,&,I%gl%sS-LpLt-+r%Q1c#M(7(&(&,L###e##,p/D>#x##B933$(IU3pq9Dp5(K0-{?A,#De,'TG(:9%##U,$em%D>#%&#(>LF>#"},
{722.78,558.67,10.86,5.29,"&##6m:6#$$##Pa:/wM######=gU'R(iK*5.*t(;%?%+]':_2SG#=M&-6'I#$D6M,$&PG#an)9iUb5%-H#p-&T93Kl%`c#fw*DI,C##`-(Oe'?h=######[.@V$TD>####pF1>y5+u####gw#Jm*&##nH$33873E######?2&A$)_,%[G#G^'pR'1m(WG#,-$"},
{324.56,567.88,11.48,1.82,"AN/(^1Bw*kI'Kq:c1,'D7`J)n|F<5#L>#9r-K#%p>#2B1KS-cYMz5$qB1{D-&(2)',[8P3'0N:P.R)iG$;v#m?(%##l7+N5$659BV8,y5;[$G93-{-Q6P5d$o6P/Z$Tl%Zx#n@.~5#[l&/##-p1=7)nJ.$q-_..Wn'<8(2xLZ{BYG#95#Oh'k>%]Z#I#%+##"},
{555.87,566.32,10.89,3.80,"u;*,6'###QG#UU0$##|/,66'sb#Fl#<~B9YI^P#*H$A9.WLRk|/0Z%{$#]S0EMR###<##8v%Me/###V7#XbJ###<##;%E)N@cR*0,#H($LIR}HR###C##@wE:(<###_##n8-######4f%+]3o8)O7,n[&An-1T2#l#{>$&=<BH'###Du#pR(######0##Fn,"},
{408.58,581.88,11.29,4.89,"+8-,n+b>$2##J>#e^(A?J@l$0M0Yp/v]3|B,ES*>q+&V7Ad'f#&hC3G]4Mu#)I)Z</hXEQH&%@OI6$c&-P<0(&0}Q$Mk5;16###[>#%E<]w*[u%dd$nF7o)<JBOnf/Q7)i],UK0aD**5J$?$###=##1~+(R(hG$;0$b>OH6&];7Ox+1o/GV,L8.BA'k_4F]/"},
{754.75,750.76,10.07,5.25,"sH&Il%###$##RI*e/1###2##Id%[~-py7]P#/,#FI#0SMs$,g%&9/.s#'###Y8)P1OmtH&?%bMS%tD<.)mI(}%0.,#5L'[lP&)7p%.1@(.9)Le-P$%weD7KSEIS.Z$Md%}9Lb82uG%'##1[&m#>x$+2c#r@'yy.<[)Sv$%C5Bu#&7+G5#Ao-SZ$ac'###*##"},
{567.99,816.85,10.52,3.53,"###_G#n#%Re/Y#$jR+gY#|m*%$&A['x.-ff2>{=_5$*m$q~.###^5#_p.n?*XA+Fc?e*?C6'fyU]&2>A)JO8@C:%H#;,9:W9###%##i')696O[*96%-d:ZYIsyU]29'/+2`4__3wc>0^6sY#######|##LV?/,#.,#P<+5'7E-%VA/Nh*F5N&g0kT2s.&}06"},
{626.33,903.69,9.78,4.98,"Jyb######`##ayb######D$#Sx3###sY#HQ$xY$###:?$xw+ezb######+##vzb######=##K{A###$##Q-$<R+###%##h[%Ayb######f##Qyb######o##D5N######:##N97######E##mzbD>####8##7|beY####A##{{?Sl%###$##6x+im+###'##"},
{166.16,109.67,13.70,6.24,"cK'6~.B?$6?Qhe$M.,$I#dS~im%mK3*/$PS~Rv&{k#x%#RS~(_/KI)+%*aJ/B8GtR+uZ&V6'^V7J91Ve(t|?YM9i#&I>#)J*tI,w0*>013L4.T~>Q$J?'NS%$U6###7##P02/-&Mu$K,$c$)*z5([%b6'k(6ZD?_P####JC(U)A###+##7u?Y>$###l##zS~"},
{76.67,138.42,13.10,6.03,"])A#########$Zi;$)###kc#f=A@J/tP$A@&z'*`#%a5$T?(|WE######$##qViUG#&##]h*[1;wd'|['0lB{C9al$m-(3-&;r@#########9YiDc$%l#@c#vWA9;0ZS,R&0O97c$(`u%1L)kx4#########vWi###0,#>5#$ZM%##b#$vJ04.,###%##E&*"},
{598.02,210.07,12.73,6.11,"/,#7<%pHQ8Q&'ZAsf/W.*#S,pt^|k####Y>#]S.'d&T[(lu&G>#J6#acIY82VrB6##o38@g5%p^###G,$6Z#rjH###5##t,&+Q#710B~/3u#=M>SS%Bh<pG#0q^n6*N5$V##*}8_B7.,#$##+##08,iH)95#F'7qc$:l$Fl$Hp^v-*fl$NR(I{3H958$%;L3"},
{779.22,225.53,13.67,1.48,"pY#SX4ax.q|H3/._R(Mh3K?P^t:OG#Hw''nS<].kG$P5#1mSsI,cu8F]3'##onSWw)ZS(~o,/W>c,$?=0#s<P94]u%L>#F6%gK69J(RRK;n,fqSe.-tG$Hf*U01sK1v@+Q$&uP$Ru$0u#$Z#66'D#$dB0V?OAA.Xc%)##t[GY$(%K1###4v#V5#0d(######"},
{634.30,228.41,13.28,6.03,"oH)a5#Yh0;K3jx~###i#%>7&(&[M>#ku$uu%s6)M,$IR&&'4im+.8&z:;+?%Lx~|P$F$)-.#sy~4[*###~,#'n*L>#mI-8l$'&1Rw'Ll%o>$Hy~[.-du$fv%f8A'#LoQ%OS,aF=L5$1$'o>#Yd+<5#rP#Tp2|w~I>#Wl#FX32L7J,$W-$l'S{/4###&##Pe&"},
{662.90,258.13,11.90,1.28,"fY#~u%@##X07vb#F8/%##O&/Sc#ee0###iG#>,#3H&######E5#T:5g,&###qQ%5:MNT5]G#&0HFGI.u#h5#R.*F-'######Z$#nISH?(###/%)[LSs}ITR'jLSkFJQu$:x%dT*M{@###%##L$#NHSk,#9I+L$$iHSqP#^P#}20SHS;5#W>#aH;s/5###'##"},
{662.90,258.13,11.90,5.38,".Q%#[(vY#YR$2C8z>%95#=2%0h8D>####8M#+6'###6l#aw#p?*;5#T-$?v9dM6k@,Z_3M+0{NWS-)b#%P7%ET4E##-u#%-#xY$###36#{CNqu&aP#@O-5LWNuNF>#^?$_]Fzd-s##t%0WQ#.7,###$##~[%P6)###H##_e.bG$###+##[A-4l$&##|$*OQ%"},
{283.84,426.04,12.37,5.03,"8-$b[*mP$(H#.c#|r-Sz=*##OG=5|<ml'8$#?):ZP####s##O.'JwJnZ(Tl#Nw*lj2pPJAv&3vOA?$&.+]J$FFH######r$#f,L,A'jR,B8'7L6sc$]59}U5?yOi6'zo/L.)(IE.,####X##v`<N?&6J*yv*r:1{9.%`?fP#/ROOR%ZU:;J#fsD###$##ex$"},
{371.40,440.81,12.50,4.94,"*RNau$pZ%Kd&}%-^~&)cD?$')fN,m%HI+5A$?_7.,#%##`6$SV2%o+,T1I@&my12h+/W;bR(<[P;Z#qQ(8(%^97###)##LT)]ZP&v#wd*<](Y1<)$$.G5n:6R^PLm'#J*&J)EcD######_5#jW>dZ%o.*XI(ry36]'nrB?l#c[P6I&1f2''$L3B######5K$"},
{462.84,456.71,13.50,5.35,"n(8L,$oJ+Vq9Sb;~A06;9@e,sgRDv)/,#8$#mo4j>$.,#p##B?@=J.<o01R(IN;2w(ZeR7S+'hR`P#f>$_J&xT6^l#Ew,w>#*>D5@)Qr?AR&tQL}P$E&)*gRI,N$##]>#2%@xn1$##G7(W?$l^6+Q$.J+=|;ydR###%##]z+Ap7######wH%5~.###(##YG#"},
{303.52,497.04,12.00,4.88,"9$%ho2###Pu#9H$243qQ))##uJ(e@JZc&<d$cv'L;-#U5Zl$BR%fRG37,?$%pS,fZ=}e-7@(_uMrm%,n(QM+=e-kH$}P:y~-kNE<v#Co1af(MK50[$^46pz56yMNS+=S+Mv&Q&-)2*<XE$Q#kz7:['&&-Nm%Co.AV*RD?Yu#(N:=w'^~.m')*o/;.&[U.a%+"},
{389.88,513.21,12.54,5.01,"SkJNu#En+2e'y82$J&H>EOe,]35*y/9&0<n&M:8{.'~^5kv'LC2#x+m^6G7'gK0UM-/+By[)_HN|,$XS+-M-Gy52v$Xk7`^6HkH?v#l&38x(&M;($$sF5]);{JN[@)99/M&,`916B(oGNj,$[177m&}%,jI+/_3OU-Ph<P#$Q>I1%&^96TK'L&1LH$$V-K3="},
{479.27,528.85,12.85,5.28,"Rh78u#jS*>U5gz5zd*{<Auw-@;7g7+/h5w03^#K`G#N%*Lj8GE6?y3X~/a?'R{:wn*2mOG~+HoO=5#iZ$`X:KB6###Q##1oG])=7$'5N@Fd%{bJW5$lR&poOqYM###J5#Zm?Tq=######0-#pe/1,#B'/0s;blObG$$##bM+R~.Yl&###Q%&um*ZP####?##"},
{495.06,596.88,13.45,5.50,"qm(,c#dB1IC7lM47{<Fg5|-)_2:9I)|C7-GA3~O###3,#x_->u$$##a:3PT/'N9#?%r]O]R)/^O[P#Mu#%Q<'`@{k#$###{+OG#/##on+;Z%0I*/,#3`/rRL;HOM[(]>#_u:{]2EcJ###I##uu%{v)mP$.##.(8{5%fY#*.&56'a2+z,'yl$+'3D#9j>%_##"},
{647.90,627.30,12.08,2.81,"0r(o:=######oD':$)###%##PI(fG$&##0l#<5#9$'^H&}k#va1S{?QQ')##^NVK?(###B##~z89_1*6'%##xc%4-I(c$###<-RtG#aR,>?#CIV=5#M>#Vy$KA1`@&_?Mh$'9h:Z.+gx3s,#}_A(##yR$<W0-IV2,#Ju$ID%gc'A5#uSSbm'5J0####K/Df("},
{762.38,718.32,12.84,0.84,";j-qZS'l#h5$t%CXf4,u#2,#Uv'(7(+u#######SR&zY$###,{9%Q%e##QyNw]Sf-)7,#x(0Ty1<+EQu%&##/H#i-(5%-###ZP#&l#;,#^]SAe..U+R3Atj>kc><_SOM8B90L$F]~,V-)Nc####VZ%v#'_,%###;d&4z3Qu%`7-_6&KW5z{=F91RE:4R)3m%"},
{762.38,718.32,12.84,5.28,"Ov'}>&###%##%$'pZ$HNC%#####K%#pBSzG%iZ(M##Bz+~R?8@&KO<=m)%##7G9/>F9q6[H'}Z(|,#dI@+@Sa-*M5#oZ$VDSY~,NB-6|7s2<mAS@S,%x)G/E~^86#$d##o<AY?(.,#2##=E=A=;`#&v5#r]1-v%,@+S,#[:4gR(Yl&###7##Bw+qb####$##"},
{321.67,739.83,12.95,3.42,"###a,#_L*-g6pb#bZ#Y/)2dFZP#$##u52Ha;######[8?95#/S)uA.GJ*]d*b7+t6*:w%qj>>c%###p;**s<######E<QI#%r6%C/1#<-re1svL^/3Ul#Sy/6f2###DH#9J,######Ud70Z%###0##?q*%'6c].Zp2zB2Cg2`254v(K5#{c%gY####?e%zY$"},
{321.67,739.83,12.95,4.56,"/m(###06#d3~Jl%######l4~OG#######<4~#########Q2~>GG###@$&3V.P4FA,$###hh3B7*;c%###n90e>$Il%###a%+h7.fG$36%1+3s/0qp2?R*Sn)J2.d;=2H&'?#iw(w|=?m(^,$4u#Y#%%B/DA,%c#>T,-~,>2:Sw,A?%[Q&+k0c?(7m%mx+7<="},
{288.02,785.33,12.09,5.83,"}I/######D##rTc######&&#XUc{k####&$#FZ$t$'I#%###s:>######[##JUcA,$###x%#cKXoo5###R##mH%AB/v$*4u#B1<######xA$9WcBd*###F''pj7LvT&##=?&0s6ZI,.Z$o>$xlJ######jo%pi9}k#4,#zG<PI)*d(fP#OV8m`<&Q%###Bv&"},
{383.86,799.53,13.37,5.46,"2$Fx%16##Of*Tc?_~0pY#q,$].,s/+8W;vu$47S6Q%hc&T[$Q&-B@+qc#Er=UAV<Z%f>#%U,Iw-h?%Ep-eM6KjG&##}P$-V'4v9?7TH>#>#$pCVbG$(##jc$Q;?###e##`g4Mn.F>#'##Y.(gR>'95######REV######4##_XB{k####?,#)-&nc'95#(##"},
{514.13,837.90,12.46,3.69,"######|##]SW-u####2##,XWHZ&######,XW95####4##WVWG>#O5$F#$OTWIR)sG$x6):APJQMF>#1,#aU/;$)###*##'x-a-'@K/^d*1Z$1vK;[)h['zs?YrA`G#[%(JV73%,D>#&##n#%h?(dZ%'{.d821UWg//#.'2&--9++{5,w*Y#%be'_5%###%##"},
{514.13,837.90,12.46,5.18,"C;Z.,####1##_04SH(###/##^Z%7R*###(##=5#<e'{k####_=Z>c%######yD5$ZLF>#4,#{L4jM@S>#1~'tG$+0+jp2AR)z;ZfY####V##SwRc$)yG$HR))b@P[PvQ(Ow'Qn*s9Zx93O['J8Z######P$#YRRJ>#P#%H#$]5$Ow*F'1g@,.K2cI,p#%dV-"},
{167.17,111.57,14.18,6.27,"Ao'g6*[##rJ]mI%Tm)wQ#+9_u6'J7+:w#v8_xZ'95#<%#|8_*_00I*Pc$xm)*.DC7+/Q$K6'zU7-n*jl$$N8UL7*6&SG#e.+ex1b.%o~.P/--9_zP#Iu$%8%{y:###&##_^1m#&]5$A#$S/1`FAe#$AZ%ZT.**C&#####`0'B}I###$##db:T,%###H##f9_"},
{76.20,139.00,12.89,6.04,"]*E#########fGh3l$###%-#lkC'm(%##eQ$iK*L#%/##IH'W5O######%##=Dh:5####uB%u:='v$zH&D35[g6^,$tl%Pl$9tI#########xFhF,$.,#g>#-kE#](:m'#A,5f3TZ$1c$|x'rg:#########eEh######)##T@V$##]G#MA-MI,###%##L/)"},
{383.71,215.53,15.73,1.79,"s%/@5#h>$N%*]}@%7)/$'29.@u#dm&UlC=K4###I##fU2sQ)Xm)bc%sb#{b#DoX#6&PH&zU+gQ(>A'`.Qe:7###U.%'nT%Q%%n*t$(ZP####^pXZ>$m-)l#$V:;lG#msXag9###vG#VqXlw0.v&BR&{k####:qXQ5$3,#:5#FkEwY#^2:w[-###I##@oX8Q&"},
{395.74,471.62,15.70,5.02,"oJ1r6'o7*v%,:O;SS*[/0(v%M6O-w%wf4d$$sPL######~n$~f0*o(+g4`?&TQN5@&ix.k9-77Oxl$~L/r]0v8O<5#.,#WZ#>U7ol$KB,)^1F7OPv%WK/w.+n}FGS&YOEnd&h6OVG#|k#'9#NB0M%)1e,aQ&Q=HNv$I~-(C,9z;4c#N~'^9OusI(##0c$.W+"},
{333.74,526.79,14.63,2.13,"M92H@*K3:E$'tK-z2<KN>u,%M?NYH%W]2a-%:C4J>#Id%gR'ew*)C2h@N&Q$U?NBw*kV9CK+)=:(.)$6GvI,e0N5Z%3Q%F>#8M:}n-u(5+^0'a>wo2ya9B&/xB48M4`FGRK.x>N=,#_l&:T(6t@7n,zB1oH&_e,@/--ANn@)L1<T#%L$(?u9^Q(R$#bH(S_2"},
{333.74,526.79,14.63,5.17,"t>%Cz0X-*w$#q,%'AE&&1,Z$QvMK]*Yw,r.*Dg3C6%c<9'/.u5&^S&FuMI5#R=H0'/W^5zU1)5>w7+)<;NB2r(32(1T_;CR(-H%]G#KyM/c$],Gv7,<X<1I'&<=]B*KcLPI(*vM}G$_A+Z(4*m%vw&(M:@5#oL9Yv$p|EF@'7D;OZ$Mh/Y2;@GA8m&5o/Km'"},
{414.65,543.25,15.57,5.03,"mc'3Z##%K=u#O`<>e(x050I&-E>/8']g7K.'E;=[H%ox*PT/:?&_9'RjDzb#>=E}-&8_5{^0NZLT-%Z22[C5O4CRw&@{<<['zY$MH#0N7QJ1-6Gmw*O02J%+UaBN9)~r>,&)hsEvu$]J)q-E###I##Es@>6(GD>3$%&)<ky,j(>R#$4R%e]Li{C###0##Pt6"},
{902.03,914.73,15.94,3.17,"D>#######7,#]m+###'##e%Fpb####{%#]B~######/(#`EFI#%######'##KA~###'##Y;/3o2/,#j0#Dp[######%(#2A~cG$#########dB~######$Q#)lMsb#|5#O90OG#%##%[$uND95##########LB~######2##$SY######t7)###$##?v%dd+"},
{48.68,39.67,17.36,1.39,"############l%~######('#.&~(&)-u#z%#^l$`L54u#-l#fG$######$##p'~uQ(D>#^&#JkGL&F8Q&9$#|S'zSXaG#)Q$^,$#########`J@xcR###-##m7'&,~pb#%##&$#p'~######$###########=##WZ'######'##_e/#########AI+######"},
{48.68,39.67,17.36,6.02,"#########$###########*H$#########Q5$############xu'######<##|#R###4##3KBY?(.,#a'#dp~######H'#7IW9K4######(##Zs~ZP#'##ce'm7FI#%1%#WGJz[([P#F(#%o~Fo3#########eo~######vR#gYIX>#-Q$j8(rg4S5$_G#1$&"},
{97.03,194.59,17.66,6.01,"4p7######%##X_h1,#$##]%$3M;rG#}5%5S'LB30u#O>#M$&(_<#########Z`h=5#/,#+Z#tOIDc#(Z$oB-cA3%##0,#%1)6^8#########N`h######6,#chf######]?%E?Q###&##df)BA1#########4`h#########gxY#########A$BL5$###H5#"},
{857.98,260.57,17.40,3.02,"4f3######.##2*q######~$#:[V###,##u$&######L5#dc'&V=######2##N*qOG####D$#hx]W%.aP#Zl#H#$}5&U5$ub#rVA######4##@+qW>$###g##=.P-B4###q##{m*<o0~>$+-$p(?######8##9*q###,##SJ$6(:_d*/e'RI&/6%#FC6$&p5%"},
{643.85,486.30,19.51,4.09,"]0'?T4U-'@_:8`;W>$###r)4N$*###8##TU1###&##|c#%{<}6<>EDE##jg/jWV-H&###_S&#]2###/##A11###9##*~-5S+&o/D>#fI#jWV8TV###Iv#qI?A{?###.##:q1###6##.{=C~.bG$###s##PSVz,'###k,#T]K0Z%###3##$q1###'##l90*@+"},
{632.51,578.72,16.87,3.19,"kK$SdTxY$###&s+6'3e5%###&I'Nq/HR*:5#;x/Sl%9?%B-(=<6Z94@j>jY#MWUHc%~Q%w#%(p5Ou#xM,ee,aL;2,#Xm'[J(?S/###/-9>HMXRU###c,#Y#7'q;###U@#LE;mP$###:v$yV7xY$1##E44<N?~./&##x5$J$;Fl%###y>#SM5######~Z#md,"},
{537.62,612.91,18.56,3.35,"1,#D?$X`2BH'D.,t>$rM22|?zx4sP$:l#08+:Q%R7*3Q$H>#[I&cj1LK2}k#ekBp@+i.*3f1QZ&;:0EPExH)d,#d#GDH&P5$ewNru$xV,)'3'$PaP#f8'&G8jH)yc%H9KJ]0nQ%GJ.8*5XA/b$P###s##eM0Zl&###}##vdGOG#%##Mc:^~.Nu#=Z$j-?y/5"},
{537.62,612.91,18.56,6.10,"6B,Sn,#A*s6'N:PZ>$yb#h[#Q39T1:###)##3[#/5ENc&###gI(zu&ge#U`>N:P5?%IQ$6q-M`58xF1T43##.=13P<e*H###Q5$jP#F[%bC9rV9]w,U$)hB,ZS,E&(z6P<]/qi:#@&LuH1h0qY#y6$]7.95#uK-ti:jn0qP#CRJ`-*>-&.L1=?&G,$>J'Cz3"},
{253.58,807.18,20.25,5.96,"############xlS######&##2zbE>####$##OB45@&j6*0,#############ZwZ######;##0}b;$)###5##FL.r/0mH(RG#D>##########@$Q######f$#P}b7[*###%B#iF8(J/J>#eu#.,##########Yy.W>$###-##m#1DOH###L##vU%eKa###,##"},
{459.01,850.63,22.86,4.95,">3>######+##Z$h###)##E7&q1=wY#Qc$Bq,{Z(z#&j,%zy+>94#########h{g.,####a$#4~W5l$4##+g*GZ%C7,*##L~($95######7##+|gW>$###=$#zySVJ0@u$M,#YZ%/r8{c(MH%5U9######Y##>zg$##:5#C/#2`=`~/_Z&Me(G.,|-*'?&vU'"},
{543.07,879.00,19.12,5.09,"9N@#########(ZgL5$###*##H=<fR-.,#=,#Ec$4o,%?&4,#mz=######;##lXg+u####/$#7KQgx3;l$fG#0/-Hq7Rw.8m#~`C######<##aVg######1$#q>OlZ'll%8J*FB5Em)}b#SD)tK6#########=Xg######$##<T[.,####[c#h1<.,####iH#"},
{736.37,920.14,17.96,4.83,"Q{@######@##7re######D$#gQNeY####k##Y^8######8$#kh<######(##Ase######u,#FxWbG$###P8(o)?|5&###pQ#ez>######9##Aqe######)$#:qe95#-##tS(Xn./$(B,#e%)xo6######.$#bpe######4&#qpe###%##o?$fd+###'##aR'"},
{152.31,48.05,20.12,1.63,"Uv*######'##ige95####X##CeVo6*Z>$)##p19c5%~P#%##4[*######'##ogePG####[##vmQ5f0|k#+##-/&i(>######5?'######&##6he:u####E##SsCk10W>$&##pc#A3@<5#0u#Y>$#########^ra(8.###,##oJ.4rOOG####E,#oke###.,#"},
{754.60,177.28,22.78,1.45,"[kM######UM(=K5du#1d$dWZ.u#|Y#;7$TSZI#%###A##zSZ}RZ1,#}k#i0#(kI9x&}x3R@%)W>[6&0?%f~-Q-)RG#D5#PS-'SZyb#D,$;'#/p69?%9C4:S+.^59R(1Q${J-P6%+%,$##H>#,SZG6&tb#R%#m6)VaAU5$R,$<Z#SRFA,$kG$Fo#*h9.,####"},
{786.46,602.07,21.05,5.92,"$&(3#;281###}{U|-+.,#&##Y{8'm(######I5#4'+[J2###%p4mv$A[K57)FxUtb#cP#+@%O=Cg6)eY#3##8,#|M(?5O###CZ&$##'**&U4&wU.,#U##8h-hz6II(z,'5##)##R_$HvU###n#&cZ%4d%+v&V~/eY#2##:&(jQ'M7%o.07##$##~:#HvU###"},
{618.58,102.42,27.15,1.63,"-S/######(##3(g0,#kY#`##BaFx>$]w*8d'Xw,1@)WG#nc%h&5######*##L(g_5$.,#S##j3Ev4=5u#pG$i>#IpL|k#*l#WT6######*##](gVZ$VQ&}#$H95J-:-`5KT3,l#yN=|u#*N?nJ3######*##,(g###xP#qQ$FT5$##Nx%GW@I#%###*d#_nW"},
{834.46,125.59,25.61,1.90,"#########'###########1~)#########rQC.,#######v2RAd*######R##UcQ###/,#YR;Z$*cG#/H%4bZ~$*<5#/,#ZAAtS1######*##N`Z>5#[P#QI$gGG<v$+91|I*i01|>%P,$+v$$$'#########K`Z#########P9Vb>${Q%x,&VZ&Q.',6&G-&"},
{834.46,125.59,25.61,2.84,"~#&######/##g@Z######b&#tAZRu%###r-#JH%bv)b5%###R~0######A##UBZA,$4,#+h#,=B2S.x@)#o'3S*0&1PZ%<[(d[*######%##RrQ=L:$##;##lN)j@Z]c$;H&8g${OK`5$sm*(###########=$#c6*######7'#sEG######x(#O.Z######"},
{281.92,146.66,28.96,1.50,"[vV######y$#rvV*c$;5#f-#Z]5}k#iu%>@&(c$T##f:8P?']vV######`%#4wV#v'###5$#yvVyG$5.+4m#^#&{-#k[UmP#qvVD>####P%#z[RQQ'.,#V##fSUZy/,d)1##kG$[_&TiC###vvVRl%###>$#do,'19###+##hR%R7PD>####$##@6AI#%###"},
{565.74,156.99,24.29,3.49,"###*-#HRPZP#j;>[#$*~*+Z#'/S,u####s5#`?'Zm*###9-%7,#J>#,P3?d*N&.Gu$sk5&v&'1SAu$pP#q?%e.,Rm)0,#E,#7H%Am'./)>c$Bu$3B/#j,%^5U)@8YEvl#);1?-(30S###I,#D>#i0S_P#F5####<2SY,$eY#3,#%1SJ>#~P####L/S######"},
{565.74,156.99,24.29,5.82,"fY#hY#M.#<J0$##wu%]h$SiC%##@#$v1$cWE######`&#:/2/,#;[$&].U,%kZ%7.>ko1zY$R#EpT3[p$hC;EH'###X(#?HR'##SN*BB5D>#OJ0;i'$|CXc#LwUs5%zP$OJ$/J/|>%^R$PS.###0:%.g6b5%{/5Z@&{5&PA(.wU###95#78$(w+Pl$T[)NH&"},
{475.23,229.26,24.75,3.32,"K%(.?9o#'###9~K5@+/##gl%4v(###&$#yL:######m##=FGrC:U~&lq:dl#seV,u####-x#1'5Pc&*##GF:=Z$TR,+##Z,EcB3lG$sxHdl$IhVL#%B#$Yv#WV7Vp8/,#q,$Ll#sXG[>$yY#|#'Pw,bl7%91y$Tby8il#eE/).+~uMG>#KZ####2?Ls>%###"},
{531.93,557.18,24.39,5.82,"oS-#Z$V5$7Z#J]U8#$:5#V-#H~U###Hc#E~*?7*###t5#D,Jg@,p>%Kv$eR+,~UM#%B5#47$x^U681###C##yx(APJ###&##>Q&UG#bd#+r=KRS`%,Hc$@p+k;.?]UQc&N>#fe%Y`UbR-###>H%A~(Ul%>l$i01r@+{u&2@(m$($~'AD8wp7pJ/{z)t]U'R'"},
{552.07,882.09,21.87,5.07,"W_<#########$4d######(##LmF,?&###E,#&Z#{A/.,#:,#.2?######4##Y2d######d##>pQj$+###^>#b@+1~+OG#j[#Y2A######-##V1d######V##Of^L5$###Mn$@q<pb####t]#*B3#########^3d#########$Ec######I##Lr=95####R##"},
{660.03,909.78,26.98,4.85,"SH(######5##HBb######a&#v8Z3Z%###Sn$}Z&Z$)###Qx&US1######3##iBb######J%#iBb.,####5p%L'7D>####:^$Yw.######$##TDb######@##cCb######<$#o(98Q&###kl$gH)######,##YBb######'$#UBb######RI#:f2eY####d@&"},
{888.67,930.47,27.12,4.89,"uT7######/##Lp`######m##so`######rR#JK4A,$###58((f2######d##Eq`{,'###s&#iq`fdT###7-#NJ+Zq`###>##U-(######(##U_MR/3###A##P~(Iq`###&##&##sq`##################4,##################################"}
};
siftPoint_save sps2[] = {
{114.38,85.71,3.89,1.57,"rb#Gc#g#&0.)m,&]I$k7/LH$'95XQ#S..th)Qu%I$#||@:M4###lo%}u'Sc%=y5{s0I#%[>#I7Svn'{99l1%Ev(Mq$@6S}Z$###tU#hS2E['y/4[A&^l%{GCg9SM49YJ1tR([Z%CW*/lM^P#&@+>$#9m)-{(X>$O##.(66G>sY#IR&@_;~u%E>#P$$t@,E>#"},
{555.90,111.84,3.54,2.89,"^%%LC86v#Y7.fc$tp4vA%RC;=.+~,%{Z#TuKFH'###(##ZAGUA(Qy5Z['Y%,(W=}U8Vl$P1-CQEQ6&$].+^0i7OUG#U#%Z(1zY$###j^$:6Olz74A--h%]p54PE?r0d/1@'-A7OF#$W>$:@#######p.#M/3$##H,$bg#>rAtc(a5$`Q#1:OKx3###$##-p%"},
{92.59,114.25,3.55,1.62,"[f-JR:C&[vb#VO=BI%?_9Vl%WxX######|##>;@###F>#M%#9,#A()7HIpb#Vq9vm#e_6Kc%:'[@5#OG#>$#~+KC,$XG#S7####F;#Ip8###')?se#~?)@6#1'[v5&.,#L$#[N;KEC###0#####u@#Lm(6-(QL:$Z#'l#78(NcKb7-OG#56#Gg4~WA.,#8##"},
{757.00,137.40,3.93,3.95,".eN7*B$##lz${.-)<U(y4~C)8?'S8%uK6cA%.,#6,#jP#p['-B,{kHn_5t$*Y8U7JQx?))$#h'9gc$}k#^7#:5#c#$~P#9,#s7/B%$x9U.8,H7Ub>#'(9^##J)=:6&D>#5##hY#dP#%##UG#-7+?r'_(>L6#a6UE##SH(k&#y&5######,-#ZP####)##X#$"},
{596.07,140.36,4.03,5.32,"?R#b~1o5$.,#~#>KH'Ju#(l#p3;Hl#+F@h03^d*-c#lR+0?6Y%#tS2######Pb:d}F$##1##P3Y/m([,$eK.#8.9J*s>%vD3-l#Yg2###1##P-S0p0###s'#W0Y?$)####'#q//_B4###<##lZ(<I(QG#We$ZYNlY####x(#ZdU95####V&#]/3D>####dQ#"},
{778.65,146.44,3.77,1.11,"D>#B&#XC<&##pm,M%#;~/|0$tI-0K1+Q%AV&7c$n}@(H%U%*O$*]$#JQQ.##XQQ0Q#cS/6M&0w,O009TQES'(##=ZJ'/.###d6)'##UTQ~P#4VQLw,;03l>#vS*,W3NRQ[G#&##?Q%>}B###^P#2##@vHo5%mT,>c%O13Qt=Z$&zI+5y/>K1###Ql$#^3###"},
{760.03,149.37,4.08,2.20,"Eh.;m(.n'JU8eG#~[($i/V;A.,#gG#g%(b&4W>$2##hQ&W[*Tz%n`Dov#me1-x+n8SNn+ln/i6S:I)c8'~-BSu%.##9=0wy7WG#(c$,%#xlRY[+a?)p.#t6S.;Sol'+J%eC4`12`]0J;Sjl&>#$###=,#SZ&9Q$/,#|c#{c(~A(0R*A/'k>%CI$wOA3f/OG#"},
{197.80,185.05,3.87,1.46,"@H'######<##EdV######I'#4dUF@#qZTl%#KV>}-$WXI[-#O&3D>####A##^eVS6'OG#?&#m|E#=-VXI?6#OC9w(3A]4v[#8J,6#$###$##wUM?GK.,#0##OI'jiVac'%##5?%|fVeY#)##p,#bG$######*-#o[-######F,#*K3######$##Sv)######"},
{481.76,194.70,4.59,1.76,"######Sc#Qu%11:<5#P#$`5$*JX%##LPH:v$B,$###$KX#c#######=(#Zq>&C9###,&#eT5'KXCm&WuK|@'X$*BH$YKX.f+######[7#jI.DH'###p'#C+JbE=e~,uJ,]=EtL<E[&Wp6=OX######9Z#OG#/,####tG#.H&f7..,#E##dw-cKX#v'$##mT("},
{481.76,194.70,4.59,3.70,"ed$Z^:###qq1*40Dq<yu$~ZBA*-RoRUl$PZ&B#$/V03l$###K#$h)9###ELM4I(uqRin+07D8nRVx+U[)`u%;d)VG#95#*##}G$61RE>#WH%cK6+#C###T,#|lR,6&###@$#e6*######;##-u#;c$]P#3l#'V=.,####a##_lR######D$#z,'######@c#"},
{263.95,204.86,3.46,4.89,"*w,4,#w*H&##Fp7;?#T4I6,#G.X9,#tJ2J,#CZ&###)/X&##W%/&##`-U'##wy8~,#v.X_G#i/X1##y/4d>#}Q)'##'/X$##P~0$##^dV.##595Cl#n*:q$(6/Xh>$&&.0f'Cm)(##U.X.##:e.%##-,N%##`&3=v$f]2Gm(/|>wu%h[+Ff)xl'(##F.XP5#"},
{363.47,204.81,4.36,1.67,",n-######&##J^e######m##w;C###~l%TZ$-H&+##lHT<5#881.,#<5#4,#M^e######f##J>O###wZ'i,$ll'(##i^e2,#'f2PG#0,#3,#R^e######M###HQ###X6(;l#H?(%##i^e6,#6J0###$##?5#G^e######W##DjG###+v'nG#H?((##MJ~+##"},
{457.53,208.19,3.49,1.70,"u6)+u#K$#005{v(O8,n[+9g1SH(~5#`R+IXXuR.###$##IXX~P####T##)@+L`C###I,#$h0#SXZP####WG3}7/YI+###P7Jxb####0,#0,#(,I.,####8,#IXXJ_=###,##&A*%WX###&########1,#ub#^v'D>####4,#YT'kM@######8##1vJ######"},
{290.34,213.86,3.73,4.95,"Sc&$##n-S###3;==Z##r:>l$;s=v$(^A.#&)}S05,#|-Sh5#*Q%&##&.S###H-S;,#C*Aq7,(]21l#xR)60IB6(###q.SW@&jG$&##r-S###p-S###VsEW#$F*E###8c$/~)######|.S4u#nG$%##rtK###s-S###h3E###_}I###N5$%##D>####@[R###"},
{619.33,240.85,3.79,6.07,"<,#Su$GML=L:<2?zG%Q@(Ml#7&^######o%#w%1###QG#6R#$###Z1I,^}>&3WB*##;X2&.),'^###,##4d#W/2eG#m>%1Q#Lc%{;3(~,$##s~26##PD4-Q$-'^E>#1Z$-Z#/T0Sn&_5%###8H&3I)mP$*##lw/[I&Fe/$##v&^7u#bG$'##3[(Be%fA3###"},
{480.88,251.13,4.36,2.84,"|G%^y%Mo2VG#||EJ$'/,#6R$_I-###,##Q$'95####$##@u#rx4;o&8Q&F##z$U<5#B##9,5z-*D>#/'#m%U######3/#2^8a)BbP#(##:Z#s'UuG%)##>/'MoKSH(G.#m<@T6)###@?9f(=Rf2###O##{Q(l)UT,%$##-I%^&UcG$9#$*J%6.-I##$%UvP#"},
{515.62,256.60,3.72,0.58,":_9?Z%###+##|VV]m+$##p-$P22o&5m8&2e&=Z%M>#:l54l#P(<;5#5##e#%sSV=c%###6^$6q8M{6B5#FN)a5%Ut2a|=(Q#VA1###~J&GQ&jWV'S.'##j5#uf,OzS###I,####EW1+c$###yG%.##*O4P#%w,K+u#[,#%m%e/0lQ'.,#+#####(c#[P#.,#"},
{475.14,259.69,3.64,4.73,"*###l#D>####me1######i,#JRWfY####-`&)S.Hw-###i=8@##L5$######rXIOG####?##rRW)R&Nc&R$#6x.3:S>%.2##;##-d(######yaB|e1###)##>VW1RI:$)8##vZ%XqNAB6###&##A0-lQ(###Q6$gwSpG$95#c.%HTW$#####%##8VW.,####"},
{592.60,271.68,3.74,4.78,"3i&h$W######|W+<)W=L9-##qJ/Tx*gsABT*9r:###Z##y?BN//ll'###C,#d&WCv)4,#yw$Gp1eG$~7'9~)K*W###*c#rm&]J1###.,#@5#n'W###0,#r>#i15D>#yv&$w*||<yY$kG#Xf/B?'###0,#`P#TmC######&##M>8>u$###`P#H?%4e-###nG$"},
{561.22,298.14,4.02,2.87,"#;$/<D######34/lI.######]-&Ic%.,####>5#RG#$##;5#Ky.6?'$##:?#.`Upb####~$#4jA95####*##4c$.,####(##?d*###*##R|,SbJ8Q%1w(Sv;9X4Dv(:Z%I5#H$'.,####&##MeNA80W,$@k1XA,Lc7g]UjPA(69U/**&1.,#Kp0[P####%##"},
{619.63,349.90,3.89,1.00,"5U9oA%}Z)h'#D$%Uc9k'7NZ&]pLz97N>#56&Kq8L5$###|Y#jZ%4H9CA1e##[W<]B.n9136LE$P/u#c>#[hRC_<9#$###gv$VI*.:*1}E{@)ceR)?&@5#5|,bq=4,#/,#Xq4+c$I>####3v$5/1m5$l8*lg+gMA###%##,;&@H'######*@$ZP#######([$"},
{619.63,349.90,3.89,2.38,"iG#'t@W>$###Te&Z}KOG####y-'c5%sb####zH&.,#######^l%GSC/a?Dc$XgRR#Na>$=@$+O>To4J>#E,#s6&bG$%##/,#ou$yG$GH:1f0O|?WtD^^259-db51dR8#$6,#x;*s$,######sV(;c%C^-pb#U#$%I&98D_e0m10Aw,BS&MC6+b2<@,kY#0I'"},
{279.49,401.27,3.45,5.12,"0-SH>#/,#?R#P0a'##Vl%Fx%&q9^5#bV<3@(cB53,#_c&$##9dT###$##~@$f0aiP#TH'tI#a{Apl#,%R|u#Y(<G5#UQ'S,#$xS.,####F##]3aS8/C,$[##5j?0}9|w0A##E:;fG#Zl&L##QED######_w#Y8W204###s##Z%V4K2###m##R97Z>#_5%;##"},
{393.37,435.32,3.64,5.68,"-/C},'######{j0lHFA,$###c(Q-n*95####7w)3-'95####25GZP#$##}v#=PA%y16l#uu$5UTf?($##C,#i8/mX@######&*D###<,#3C'RT5VG#V7$1V1YRT*p*$H%~7%>?&eKIac'###|p:###;##az3A[*$##El#SB,,[)}~#HsFS-%%##{g$|2C95#"},
{484.21,450.22,3.99,5.87,"}r(P95pb####5#3FX:uG%###>*U1?&.,####*d'G,$H>#=5#E)7Uc&###@d#e.K>S,9,#:R'hVY2u#$##|P#5U7$#####O5#v(?<5#wb#`U'.D>/,#<-#`AH>SY###3##Fg)}(<E>####&##o_@0,#Cc#1t;Y08###=##(M+m+M###'##Ef$0K395####.##"},
{725.32,530.97,3.31,5.02,"Ff0[P#######HiZ######<##PyNt-(7w--##.##RQ$ugZ###pg8,6'###3##,fZxY$###4%#8&U)_:g,&X##<##W@)ivQ###DuK=6(###i##bgZ~g8###u##6G7keZ.,#%##Q0-$7)?x1###KtIeY####G6#m'6(90.,#<7'HuInn.###d~&dq9;c$yY$B##"},
{738.15,531.03,3.75,4.89,">2=VZ'###h##WmQfH)###''#R3=V(9xd-_##+##oZ$WnQ0Z%:mQ>c%###m%#goQ4mQ###g6#'hOIYK-H&<##'R'26$)2<n>%tED.c$###:W-vS0,&.###IpLUmQ0?&95#bh*~6(v>$5H&cu$|{A###)##K2.;f/ZP#0##i#Fo,&D>#)##<K-]#&######|>#"},
{667.61,534.87,3.38,3.74,"OI<_5%###6##kt1Zf5###$##wf$-iA###1,#QA#cL<######H*4zY$0,#O>#msXpb####t>#2M9@H'G##?f/1$#yS3r,##.,(A/:5#2,#H0*SnX######8K$z-W###:##rf/95#.,#B$#_{A3l$###l##:pK,nX###F##3s1T,P######?K&pb####3##Z&1"},
{606.89,542.16,3.90,0.75,"/f1$##X>$'##(r>###^16uY#R.,###W<=<5#FR)RG#UG#QG#yJ30,#95#2##m[W###{8+s?&LI+###X'IHZ%$J*[P#/$%D>#im+######l5#6]W###`u#p6&fND$##Pu5WtCkl%###853o#'L#%$#####K##r_WdG$###g>#`)ST4D]+<qx/@,#q&,V_WL5$"},
{606.89,542.16,3.90,3.26,"wl'/##lH)'##,_[&6&4c$ep&f=C6`-&a[1q0###'$#3wQN5$^v*O##dw/6,#N_[-%)}>&z##O<;h,4#ZP(#####&/%S08.,#bn0/-#M~0z,#q~[un&z,';%#}H)z2&/1;%#####4@'-H&###P]5r#$-H&d$#Rg9l%$WI-V$#;5#z%$v6+&##$##vY#~P#jG$"},
{751.42,543.71,3.60,5.97,"eJ#;$?W&4###GSCIV9y6+L##CO>>c%###V###-%Tu%###$##W{)j'69m)&##ReRSd))m(p%#b4HKo2D>##J#gP#q%-95#2,#9>8}v,###n##ihRa/4###{##:l7_cR###:##^f-^I,###6$#AK5######d0#<C7L5$###6:$qaEI#%###pJ#]h=4##|k#'V$"},
{487.31,571.25,3.53,4.95,"4w-U>#U5$9@$hRY%##)l#L},C7-T,#v|@7c<:5#%$#&nUQG#kL<###F5##f'S/^iG#cc'OS#n:=e[#z/^KH%###4##*0M*Q%*tBmP$###,Q#31^0w)I#%b##:'4*q,QC;<Q$r5&1,#zg*Z$(PV=OG####nd%fxV)~-###g##7h:AS-1,#M@#ac'###O$&zu$"},
{772.05,571.03,3.91,0.59,"Z@)(7*OG####eaDD>#e80@5#w?*###QbE3,#Sd*###SG#;#$cB/<$(######R&W95#d.,bG#q@.###p'W[G#:e*###$$'###|'5.,#######'&W###8?%J@'n^;###l/E@`;V-%pb#WD/eY#^S*pb#$##]G#0*WxY$&##a6%q(W#_5dCMb_76##B$%[(WZP#"},
{430.73,574.42,3.51,5.65,"jY5M$*###$##Vj.kmO6#$###Y@I>A._P#/,#'R)###hQ%:5#O|@OG####=S$??CwR.-##r7*JKT:?'$##3c#]81x&0c#&###ROI###:,#$V(tK795#Ov#bY?jITGL.1l#s@'WH&$NTL5$###qy9###Q##_i8)8/UG#@?&NJ)~u%d1&.jEGc$2,#0j,mg9-u#"},
{337.54,580.22,3.93,1.59,"mc&B_2F[&7q:7@'k/2$Z$W.)M1)CJ095#2##|'295#######Y'2Gk@n#&+Z$<o,IzRx5&T>#0zR_'4OG#>##^vR######Z##d>#lSEzw0eY#~s?,%R.,#p##AvR_u%###H&#8>N######f.#K,$dR&GU73?%RvR.?&?#$3(%kuR######:(#zo6>5#95#i$#"},
{615.99,589.38,3.63,6.04,"aE*LKXE-)###`F>=OX95#:##DiAuR'OG#1$#`$+######h,#]H&>)6TB7###EJX?8.###f##9wVj>%###U##tx1T,%###~##([)>6$L;@'v&(JXL>#0,#rc$7bFLm*###V##q:;B-'.,#p&#|Z)-##_83hv'OIX###%##+e%vz>1l####o'#hv+$)%bG$W'#"},
{517.36,591.58,3.94,5.79,"ii(*94OG####AQ6wq9j>%###<_X|Y$###/?%ol'######gn-*;7iu&&###S$mBQR6)$##H8*p^Xpb####~u#Ty7######{>$dL<###7$#g),)NB###G##A'I}[X###+##z0*GD>######'##0&2###c%#g}B,g7###(##6i-m+M###$##rA$uJ295####9##"},
{507.77,612.80,3.86,1.60,"cH%c]P######;@$x5N######$M(.o2######u.++u#######Yq+6uGmP$###w]*PyT###'##`{Tqo4###<##*PF######)##D6#/lFfc'###CvL4>J###}##1wT*Q%###t%#VHS######d##+Z$D-'%A*j,%ivTn>%K>#]S$;vT######,&#zV@######~##"},
{626.13,621.63,3.74,3.87,"YG#J[$5mT2,#S^.XR)YI,f-&lQ$6~-nA+1?&3&#h%0Rv$95#9Q&'##GAUJ5##@URG#b923a2U$*Ql%;2&HD;2##v>$w9(3l$6Z%F>#~/S$##NAU###oq7'$%uFI###9%?(d'###)##?3+X-*:,#yu&/{<###XAU)##R)=&##4AU###tX@)##D>####<:1Q6("},
{811.80,623.61,3.39,3.66,"2##9T1^l&###I#MA[%:.-$##3eV.,#HeScG#OG####|fVuG$$##FH&>w$8~/+HN%-$r7(}m+=fV1,#PvQKu#{G%$##qgVk#%######p$#CM?Hw.###M.#<%OpdV%##`p5Z;+X,%$##NNUsQ(7u#D>#,##K$)@%.###+##H'-f83/c#0R)_A%95#N,#[=?(l#"},
{326.17,810.64,3.68,2.10,"uG%O##y`CT5#5iA.##*7*~q&IZ&.H#7OAt7)A5#;.&h%0###Yv*###HcFa#$~mQH,$rn,~('P[*EH$yoQDI'vb#w>#x::1$(},&5l#(YE=5#u%D>95c%,Bc$N_3O<=GHKEz4.H&'##oF3hmQ'Z#9u#h%,pb#DH$~#&=7%&o1MoQdI-eG#)g-8q4PF4`7Mw93"},
{326.17,810.64,3.68,5.06,"syQ5C6S$([:*g$'RF6:xQ7J-4@$7J._.+_5%I%*xY$qY#M>#L=4y_?95#)##=wQ8_6#^-&)6Fw+GZ$QoGkJ1#=E1,#_l%@c#0;>Q5$wb#sv#4wQGJ(TR+a$$r~/|('SvQ+l#^vQ9c#PQ'<,#gw/.,#'##Y%&((6t^-.m(Pl#/80^z%zS3O##ZNBYG#D>#Y##"},
{362.86,812.85,3.69,5.65,"j8+L-%v~X###V9K>Q&<c%&##af.'Q%/##oY#MZ&###5-$GZ&6A.[~&f^X`/20]XE5#n$(}[)bV>iY#3,#$R&G7,]G#tu&B5#_=?J,$J|*n]XQkI8l$gF/#O@Ig06e-Y5$,l#H@&ll&pb####d40R&2Oe)1c$ed#9F6+1TJ#%3]+;['-[).,#~R&%Q%######"},
{595.57,827.65,4.00,3.67,"sb#EZ$Yu$?m(5l#M/+L6%d..6l#Y*>8m&4?&J/PZo/B6&8u#>c$>~'&[)]#%Hm)|O43OFxb#'.PaeMad'lk8O5>zK2I#9_V;A##.J,xm+.,#=7,Cm#B1P)w,B1PeM>Z/)eYD}B1>6E8kA<2:$##sb#dH#o#'###;5#xR$]m+Du#3Z%86#JOF|,%[v)rY#[17"},
{535.83,839.09,3.99,5.66,"b,$pb#######E#$6Z%$##;,#Il%~P####X))C,$95####)CG<c#5u#.,####G^8;-&eY#[##DU~vG%###wE)'90;M4###0@9?u#nY#,u####Vr@/c#@c%G##bY~N:<###d##<r(&X~###4##L5$###_P#E5#Bd'3l$.Q$;#$pn&7EC9,#cG$5$#%'Z9,#6#$"},
{286.34,846.09,3.60,1.65,"l>$qF5*&0.,#w'~tJ/s5&K##XC;]%#j%~7#####+$#&lO###[x2Js7OG#,##>'~8w)yH*x$#NU9V2&k%~N#####+$#]3G###Y8._H(###%##}MRKROxY$2##a?&N9E^rC$#####x##OR,###5##n>%######a##_/3######%##mw,#########*########"},
{286.34,846.09,3.60,5.21,"#########-##.,####(##P/-######+6#}8395####/##MH%RR,######`##]?T###}Z$PS@LQ'###}fF&uFFl%###)A*4x%HB4######+##%EY###8d'EI%jr>.,#NFY(I(*n-###`:2}{*{c'#########KEY###[P####oJR###Zl>|G%Yu%.u#~^-58*"},
{750.54,869.98,4.11,0.06,")c$:,#X>$*##,mT%##[P#{&#TnTG(..7,O%#P>#Q)(*mT###HR+Nu#ZP#1##+nTCl$eY#f&#A2?/g'QFJP$#pb#b8#*mT%##1['&H$?Q&$##[zOax3sb#9##qe(TrT{v,&##;c%Ci(&aF8##=5#$##'v%D>#=f+E-(Ku#&l#9?$`22DH'G,$>u$2w$(Q%U[%"},
{741.49,879.68,4.09,5.75,"xb#;#$I,$;5#V,$_6)###+##[Z'*6&###[;%<c%.l####Ce=###L5#kH)###L06gv&G-)O,#iU[|c(.,#iD$PB1eF7.,#~s*###)##-7+####:84l#-//$$$EY[gOI/,#+6#b0'gX[95#(##/,####|u&.,#8-&86'zR+I,$}H&d}D9Q$)m(=##G~RXG#@~-"},
{703.15,938.46,3.96,4.73,"U08######|##IeWr5&###P{&ZV;z,'###Rz(Ge.6l#OG#@##OE@######3##DjW=WC###v,#?jWrK8###m5#y&3^P#.,#J,#Y.-######5##.gWiZ(###i$#UeWmP$###N'#y%1######A[#pm,######S##XdW######]&#XdW######x0#@H'######U?$"},
{879.55,938.83,3.94,1.31,"1H&Oz$Jm*B##QR,<z#k-+5A$9m)Bm#Wv&9DRD$)b#$Cw'dXE.?&9h$&V=%##v?R_~%2-(4'#j(=1|,`~0=''X##$b;Wm)}Y$pc'a6#4nOO5$~CRg[*).)dQ&yA,=O4KM?K>#}u$NJ-cw.^5$+u####iB(SA/gR*J#%a~%DARZl$6%+9B+l@.}k#mH'o#%Zl%"},
{690.83,941.64,3.97,5.01,"Gw-######$##c:W.,####q##h:WOG####?$#0B2%Z#W>$'##hn0######/##98W4l$###s0#fdOy.0###IU#3y/Hu$pb#,##C.,#########c<Wy]7###<##H/C<7W###2##'8L0Z%###,##3u##########[B-SH(######|9Wpm,######l7W######8##"},
{892.48,942.49,4.08,2.07,"###L&#;}J/,#X&*wL,By7z.&LE3q%.N6(~%'(l#;c#eK1/,#U,%7$#.RTjP#<>MX,#nC;(4/B@,Ed$Fj@GE5aP#Lp+vm,%##HQ&5,#uRT%##`VTzG$@h<PZ#oy1x[&jST4l#X5$2e(?^1.Q%E#$uP#jRT###p3.6@+8*?ZP#$D'DFCN&0=u#{b#4R&*~*p$*"},
{847.37,952.10,3.82,4.77,"mcS#######%#jdSP..###i##&?%0p3F5#1d';-'e-&Ml#<?'ysH<H&###n'#J>G%<A###_##7eS>w,8e'>SBJv(eQ$aAE@W>m6*17(;Z%Lw$u+Ix7-95#m[$+iSWlHxc'Gm$QR&3UK{sE}P$L5$;,#Kd)H@'#m&|5&$l#A%&Re'_C9w>$Z[))?#oE<'o)`Z'"},
{929.49,964.23,3.61,3.72,"A,$###-##r47{k####?##Cs<0Z%###'##I@@L5$'##z5$ot?^T4###z##J;4nA3###=%#Y(6deQ9n+f5$AI?FQ%e`,Z*?y3C`AF###A##8v'so0###J%#7NA9s7$g3t#$Jq;1Z#HF19*9.S/Y#9#########ev>###:##(@*h$(###G##m`C$#####_R##<C"},
{38.56,20.86,3.94,1.25,"^Q(######K##a?U######R(#AcQ+%#cA3R(#pb#*'#=WC,##ey9######V##H@U&Q$OG#g'#N`B%{%sEGi$####A&#w0:###;^6######(##5EU+J-D>#0##/B+869<]4$##%##_@(^Q(###SZ%.,#######Uh(@%.######W.#,6P95####5##^XI######"},
{103.69,38.83,5.01,4.81,"uR-0##.I)dG#v(?######T##2]_$##FjG,##3l$###5]_###oJ3w,#<$)=$#zh?[>#(c$i,#I]_%##>lP+##%Q%###>]_###Pn/*Z#{k#j$#e:=g5#}>&_##5]_$##O>OG##0Z%###F]_###k,&Ru$###T##gK7C5#.,#R##:]_###``D+##j>%###D]_###"},
{93.27,64.02,4.16,2.13,"ac'(##BA0z#$(c$7$#dmU_u#jB5l,$y82W3,XH&R?&?q1{.+$K3&##^^9M5#{S3&##7nUl,#dmUP5#if2^q'u#'q-#[oU=v&jm+-##8<9Ku$:~.+##hoU)##1rU5.)HT2/##%[%j](AnUD>#/,#,##vLND>#H>#-##O,E95#HJ%h$+zT*ac';d$mf4fR),v'"},
{93.27,64.02,4.16,5.14,"991mu%/d$[J/7]*j>%Tg(}v+/?H95#C#$3##]'U######7##VdT2,#)$&;h%xJ2g##k'UeR'P%U7##gS1A##?4AtG$lZ(A##)fRu7*{,'c-#?'5Xi)G$UO,#Q$U|5#T7.>$#m08b>#Kf4;##Yx-D9-+/+H@)2h;3X-Kw,}H#S$UX#$pb#5%#V@-tu#^Q(4##"},
{128.29,78.21,4.29,5.50,"i@*w_*3YK1,#0TTgl#+R*2##$:5~Q$Y-*$##:9/4u#6#$###[96Ym$vb5;|>/RTP5#3&'/o(gU8'7%<.-D,#-C3Mu$[P#'##6{Rk#%gyKE8/ga@?l#._PkQ&,:0d#%k7,&m&5x/*c$WG#$v$).>|L6/]2%##`7'oL.ojG&Z#+L(VH(<5#-Q$7.$7?'###@5#"},
{856.17,104.22,4.42,2.25,";c%######x,#:$)###%##N)'95####g>#Y=0.,####'##x-%J:;######|##?IV###B?%Os,h6*###a$Db;5D>####X?'c5$&E=######)##^NVI>#:6'>H#pD8#c#xLVyG$eY####Wx+qY#0R'#########^NV###.,####UhO.,#NL*Mu$=u####Vj/]P#"},
{275.23,116.64,4.94,1.84,"###]>#gl&###P5#iy/su&###D%)X&1VQ&/Z$D>#kP#ec%]>$###s,#J96###G14dA?VK_a,$#Q_DG@A02b.'/[)mG#QC4lc&###8##72<###o~1:[#RQ_r?IZL_z6'K&-l8@G$(xZ%=7,RG####F,#vH)######D5#`~'R81fl&vb#h5$-]-vb#}5$(c$$##"},
{695.97,136.84,4.69,5.94,"F2%?AXSZ&###wF0fH)######K.(;c%######;c#iH(3l$###s21+@*yi)EAX1FXeY#;##_M14OA3H&###1##yP$R?'ZP#6##[?%3/.cW)*AXEN@bR)eK-J8I$q,rV@{k#8##lJ(2e..,#7##{5#=f0nl&D>#7##Jv%FV7ml';R#A&2Mm%j7/Y8&OT5/c#4Q$"},
{730.88,153.19,4.57,0.26,"hc&Wl%[P#%##SR)em+###zb#?5#OB1aP#^Q'pb#3;&k_<#d(~5#M`=######?w)5eSQG#7,#3RN;jAeG#lJC=-'|`6TZ7YTWO>#iA03l#C8.vI,j@+Xm'rJ1pVQ&TWJ$).8(=R&,XWew+m.-###F>#Z>#XA0###B5#w'58v()##q&/5,K.,#j##~y4m,&###"},
{730.88,153.19,4.57,4.96,"F>#[>#2&0###,o/4l#@&/P>#2<>######]5#cu%[P#%##QQ&C5#eU6eQ(###^7+fv'fA12~,O%TQG#7,#Bw)dm+###zb#kd)C90DPJ.,#*##)]W)d(YS(yMRL|AZ>#mJC_?M[T1VG#QH'?5#x'6Kl%###n##8aW88,IJ.AR&AE5?H7]]Wwl&TV&^z:ul'ZP#"},
{837.95,153.93,4.82,5.79,"1&+p,$4fN.u#dzQ4,#$v'###3L0#Z$eY####Fc#>Z$%Q%###{R+~g+4vQ1Q%xvQJ,#/.,HH$:OEfY#$##Hc#?5#r#$8-(~P#:)<p#%pK*dwQBvQ%##31%BwQka@;5#pP#Ul$`G#5@#zT7}k#[9Ks5&;?#g~.8D-nL4om?ho49U.xZ%Fe,tb#&?$tT$W_?###"},
{826.65,163.06,5.05,1.36,"1,#K6&V#%qb#H#$J90PG#I>#f?)|o('6&F~)W>$6w'>V44y2###fQ;7#$OG#pd*)qR$##w,&[nRi|:g5%hp'lw/2{'.nR_]-W>$4N*^,%ws9QR+qx/lw)VnRLoR=l@S-(pI*kd*yF3kZ(7,#[/3^S*]u%@r*~P#aD-Dy3un.3,#YJ(s@.~P#Q..y>$L#%8##"},
{210.87,187.71,4.97,1.65,"g7/######5##T/_-##?-('$#grAFm$EYHMH$erAL#$}6+/H#9K5######I,#V/_>,#7e.W$#S>Mev$|/_[##@-O3R(7?'C##>y695####+##82_Se*c6*I##_C8.?:RYMJ,#DS,ILQC,$4##>v&xG%######yN40<@D>####Bc#80ObG$###$##tL4######"},
{228.12,190.00,4.75,1.67,"7f3######&##?0b5##)R*6##/r@{,#_J]g,#IND&##'d(%6#Vp9.,####4##50b3##HI+f##f)AvQ$CvNS-%`sF2l#V@-IQ#/1;######F,#20b2##,J/'$#^*D^[%b0bN,#p>Lp,%AH'7##-g7######,##11bUZ$?R+D##rK6%T)v>PW5#Ex/}`:Y>$5##"},
{812.67,222.77,4.51,2.09,"I#%N,$:5#|C2pP$6#$eI&5YFUJ0M#%*S(/<>W93D>#P-#6wQGL9######g?$a1=###g8((A+A-S%##@2+cE:g-S3,#As-@pN@`8#########F.S###<##}>%SbH###'$#W;6=1SF$'Yg5v7&`z2F>#:5#$##GdM#########e}E######%##.0SSG#@Z%ZG#"},
{812.67,222.77,4.51,4.23,"A.RWz1{Y$sz&cR*:KHA95$g+1##*<+&@+######B8$6[*###C&Fq(RZP#F,#W5C|1R95#c,#jI*sL+*6'%#####n:%%Q%###+25^@*Ou$dZ#F}FJ~+95#(J#ll&0*.pb#'#####d9$Q%/###JK,{k#.%)'##@y(:'5$##*##4##Y0.pb#######a?#s[-zP$"},
{665.46,254.66,4.56,0.37,"###GA'kZ(###5c#C>G<5#+?%+.#yEGcG#,Q%k@#pm,######A$)VK%7+J6##ocNu]/}[)M?=oI'Eg2PgN5<>r5<l?)5@++##J[)Hl#m6KE,$PgNk3A@A*KS+rR(bm@ecNT/-qsGG##z6+KB%37(:5#C%&,$(vG#@?'T@#$XF5,#g?'_G#(mI]#&/,####6p("},
{665.46,254.66,4.56,3.67,"###79(Ru%$##VG#alG,##b?'[R#U|F6Z#'-'Nn'EZ&Nd$P5$Xd+fB%Y<E?##jcNZ/-O.(Wd@/&+ES+PgNO3@(eN#l#Lm(Nu#QR+1##WH@(d(PgNk2=zd'a'3e@):?<pcNQ00TaH2##7v(NT%######Fe#0e..c#{G%J$#NGN3,#&?%}b#kYHfu&###$##qw'"},
{533.58,257.06,4.71,3.16,"j@'cA*V96]>$g`7N5$###`l$(Z$F>#(##Sw-###s>#K?'r#'#g58~)g,&)##ewX<5#$##~E1M[+yH'F$$xxX###In$mmKhB8YT4###[P#$##<zXQ,$T5$7?$ZD=ia5Zi=j[+$##v3,kbL(c$R6(###=5##l#t|X>n.J>#$##dL/eBN5E@%#####xf$^rBW>$"},
{264.75,259.44,4.42,4.96,"0]1.,####%##[*@ju&$##Vm%QJXc#&)##T$'.nU###eQ(###CC:######+##aIX###Dc#4;-yIX$##09)o^3IJX###v#'O>#['7#########fKX5,#^c&9c#YmOL-$/5I9Z$WKX[G#<Q&3,#`R+.,#######gKX&##vG%###pGKzP#^uO$##,KX<5#RQ'%##"},
{578.53,281.47,4.96,3.61,"%##8pR###;5#VD,NmR###%##E6;vw0}b#q5$mH#c_>}>%ub####T@*F>#vb#UlNDc%$##*$%GoR3~+{&.{r.dR&-<8joRx[*95####$##qo*v/4%##qY#e+?:qRs7O/@)Dn)H6%[wB(M;;Z%2-(###4,#eT%L5$=##{Z'G]+>c%Bc$/S'Oj9@5#xP$nB--@+"},
{417.54,291.53,4.02,4.88,".,####>5#-u#W+L######K5#47X###oI+*d%{k####$^Q?n-:5####8,#Y>$87X###&##c5#W7X[P#(H#{|2aZ'U$'Qz&l8XH,#+u#&#####u5KD>#$##+##%=XEXF$##aQ$uR(%=XM,#.w+d>#W>$######RR'pb#######PT&WED######8##9RK95####"},
{182.85,319.89,4.30,1.70,".,#5##H?']P#R5#O(1;d)###GT.uo3c>#e?'Z>$G>#G6%cQ'###(?#xp9bG$Ip24%:JB^7Q$QF^ns<h~.Fo']w.^Q$n|A7R'&##_5$y::95#>w-(d#)H^^SQ-C^<R'|[)&pCM[)@@&1n-Ul$###2,#Sl$Vu%######-I$AV<X>$$##B##7(695#$##A5#;6'"},
{365.41,348.29,4.41,1.87,"###%Q#SZ&eY#&Q#(h3I#%###Bf-4A/H5#mu$D,$U,$Ql$Jl$&##H6$-95###X(54KAgo^$Z#gt^3|;n~.Aw&2@+AI%WV98H%T5#$H%<U6###To3.m#eu^N,Cpp^6w'nS-y-8pd*tU*){:/u#9H$3l$%c#.,#`>#3l$'n$Jy5gZ'Ku$>u#FK.C#$I?%T~*g5%"},
{535.44,373.33,4.56,1.81,"[5#hv*###$##e>#e{:OG####Kx-9A/YG#3-%#$'8c$V,$>c$(H#bJ-2I+$##}:5;:G>A[|Y#iP_f|=,8-gn'l$+GH$:{:N?&W5${Y#A93V#$Gp7D.$RQ_M#C5M_WS+{R*X?6`d))o)^I-4,#95#%##Pc$Ou$>u$###z[%C16l[,C,$-##~y,g#&rP$###E,#"},
{229.03,415.11,5.25,2.87,"@+-./1OG####jMO{k####$##Z~.###$##ZG#U,%######YG#;u@dZ'CJ'#?%it[######Y,#tFFdG$###+##{#&ZP####4,#JS0%##:b.Wk@zn[###;,#2k-k::tc'w6)`Q#q7)q-+9#$*##95#%##;0(Sp[[?)###7,#Am8|Y$I,$/8*7m&A.,{G%ku$}m("},
{294.87,444.54,4.55,1.48,"u5%Dy,]7'6x1A/*M$)*###T,;O9@u$###z5#z7-*?$zP$v5%*m&:ZDy,&8,#>T+qwR###.##&{RTn.###>##rbJ-$$)$'bn'#6%:?@iR-F5#t`@Dr?###.p#uvRz5&###7q#94G7|,[J2eI#|Q)NZ$pR,4~'ruRE>#/,#.;&o'9|H&qd+Yz'_H'~=2(^6B5#"},
{382.04,459.72,4.86,1.39,"&6%8}:vR'fd,WS'lS1mP#Z]1xM21Z%###Dv$]w-?d$fG$N,#'7'h4=Pe.###2:-K@R###(##9DRMn.###?##sT34@&am)Do+OQ%Xs4id,7,#fjCn08.,#*&#;@Rhc'###4(#j:;]E-;'7y9%}c(u,$--'2R$c}L/,#2,#<z#Yq;_01aT4Wy$B$',hO8S/B,#"},
{635.44,510.99,4.49,3.83,"|c#u83Ad)0,#'nN1.(mL:###tgW5l#O@T###wu'&##KPH[>#%$#NXGsP#tc(*vJ-E9~c&V-(lgW1u#jW;FZ#@H'###ynFS6'.,#Y>$a##m*Gap9b>$V5#ETM3eWu>$sw-jU&Y6)dl#C8O$Z$eG$.,#$##,8-`98######''*FuJ^m'T,%O?#.l#GS%;K5###"},
{316.33,513.94,4.64,1.34,"g,$Y|7|['5?'/w&?S.?Z#Dy4ph5Tu%###xw(L~/JH$[P#:?#8R$*wE?d*###`p+|-R###T>#H1Ru~1###L,#&U5t$'}>%{[([%,h<6kv+JH#*j@Ff3%##.H9c-R2H&&##1F01z9b2+lZ(u.$zA4Uu#k,&7h$v,R###xY#0k2.w+if-:^6gT0#$&Y8Fd$+8,#"},
{401.74,529.81,4.39,1.71,"NZ%$.({,%hE>k~)q$+qb#/%%m)2gu&###.##mn.S5#dR);o-ff.T{8:Z%B5#tS*@~I=c%V>#/hP>~-###J##iePZd%A80#x'Pu#G,?|.0###vsC;GG95#9o$NdP#H%>,#73-ndP{10#E@RB'Cu$vZ$<C7Du#@dPEu$D#$xo%x}E|Q(+n)5(/0w(ft5$;:nY#"},
{798.45,537.93,4.53,3.68,"b##+2>0H&###HQJ=o,JB6###lhZ-c#0RR###zY$&##.~Ot@,B##z08f[%8Q&^QK%`7tH&$R(XhZa>$yi?Rl#K?(###38Mk-(###$##{[#DB5Q97A5#q-$w'SXeZ{Y#bw,/3*Gv)l>#HG?V?(T,%###.##3@'?L:######{B)iFH)[&mP$pm#>#$w$%cI-###"},
{489.04,545.38,4.79,1.24,".?#Z/SRu%###5]#3-S######D0&;C:######q[(8m(######bH#nG?w@/###^(**-S###+##>2S0S/###1##aS.</**-&4?&+$&I9-#A/e,$UM;j7/###(($%.SfG$###hL#(J.F&$Wo2#q+8-(rY#:Q%M:-1e.)##9l$n{*&z53a>J-)Kp$A$'nAF.7,>?#"},
{356.39,555.49,4.76,3.65,"p$#3aDhH)###i5#$I'CC4}>&###YZ#1G8mZ(H>#3Z$4x(DW@Xq&.bK0,#C,$4pK>R*_d&RR*=$)T##ZD+{GMa~0###AH#+TPZ@*g,&m7#)RPORPOG#D##{j2og:j>%o$#QkGQH'Jc%16#bGM*Q$s8.8|0L-Pi/3s-*W6&9N15[(wNE.,#S5#B##:i?<5#pb#"},
{601.68,558.68,4.38,3.31,".,####re,pG$xY$###[d&?#?mv+###2.#Tr=8Q&###-##b,${R.1,#Ad)LZ#ZJZ9u#Ye'&TBX?(ou#wOZebE95#&##@x-GQ&Oy6eP#Fd**##_OZ90.-R)f?$FA+E+..KZo>$###|#$c/3{k#me1E##J[+;##*KZW~'5?'M##Xu$OD&5U9######qQ&(c$###"},
{766.73,588.96,4.21,0.21,"*L8######@##jGP###1H%>0%d6*U##QmQ}H%*I*2##d$+)##e|D:5####z>#oJYQ,$]7(C}.Zv)A6$ZOYt]/]v'-##='5###+|;ec'.,#@5#ZOYaq9_c&[c#me'Ww>UJY2l#M6%#m#:f3###D:5+u#######^r0kS1######dQ#g=:g,&###c-'IZ${k####"},
{766.73,588.96,4.21,3.29,"eY####Cd'3Q$%Q%###A-${X9###$##L{/f&3###$##az995#v~2###R-&eZ#a~Y4l#[e'jI=^c&gc#raY62:###4,#}*@[l&e07###<I)0##raY0L0o-*D-$*w'E}.y~YG#$$##JZ#rWC[P#V7..##<m)E##Y,M'n%Rd+n##[5$f'%/tJ######H##^08###"},
{650.67,595.65,4.35,3.84,"C##8A/VA1###fnW}>$D4F&##KfY###*N9XG#.,####X{1=I++Q#f2@u>%Su%CQL}S*y&4/$'dgY=5#Oy/-?$]#&###(<1)w+#V-eG$A,##T1$^6hP#YZ$.AKEeYD5#`$(+W)P$*I,#+t>Al$[E)0e.###pY#7#3'/1###[I'BdB,7*mP$/d#q>$Nn(;/2$##"},
{685.41,600.09,4.27,3.06,"/##L5$$##rb#q`CA,$###V5#M)j######D$#BO<###$##C,#D,#s5&###fY#sRXeY####r>#d(j###+##kw#J03rb#a$%pQ%l5#gH)######jB~L5$###-##b*j###&##1v#8h<;5#JA%z80<,#6#$5##o?*43=OG#(##e#&[.j######8,#30]###A##-y-"},
{685.78,610.69,4.26,3.12,"/$#.XG###$##&t4=bL###&##>Jell'###&##*<4#########;,#W>$$##.u#?.WeY####J5#fCe###$##|,#1a7OG#7,#[>#B,#9Q&###[P#^S[pb####F5#6Ce###+##I$$8B3-u#sm$Fv&M,#TH(######'5FbG$######8Fe###$##fG#xh>###Bd#O//"},
{750.97,624.55,4.88,0.18,"90+/@+}G%###_=ZNc&###$##s[<q<Z0v($##@~)5v7z70###QJ0$##k-(Q>#|8Z95#vY#Yc#m95%]/V~,p5$y&0z-(Q5$tl$f$+###gZ%pJ1i8Z###qY#6]'@2?$##OC1K6%-K3$##u?'57'[P#.,#Y##M8Z0&2###(##1;N%m(###2&'>K-LS/###@A(rG$"},
{712.52,663.40,4.75,2.00,"wO6Cc%BQ&0,#ed(RQ%-x/TG#Z<:~>$ZP#.##0'6######x##B=Ywb##6&<5#{/1#~%.#KC#$%:Y`Z%lu&Wc#dXI######W##B=Y'##b5%###+r9zG#t`<~n+okM5,#Xl$d2+8NB######R##B=Y######$##|xRG>#+l#k?$f%0###%##[A$o@,######O5#"},
{787.10,669.27,4.67,0.68,"g8G*Q%&R&Ju$6$<vg0^ZQ.,#Pl>|c%An.>5#Ln&OG#$##f5%}|5X>$X%%4K3GX2PII&M2&S-o^QYx.(6&4Q#RL4OG####N5#,:4TG#f8%-V9M//SQ%q};bj=WZQiY#+m%,O/Iw.######m%(#H%D>#v6%2V8qb#/,#wY#TRC~#&###%##@s.eY#######&n'"},
{373.89,672.89,4.84,3.32,":5#3l#a$)qP$2,#,V/&e-###VQ$ED8|k#@5#<[&^#&$##;5#.,#3##,h6bc&R.-/=0u7ReP#X8R^O>Kp2vj2{y,X'/3E=Qw+N>#,$'_J.>5#uv+Wu#-;Rh[,Y:RHk>`F2S+E$3=1k2^D=%E7*###&/F>#######M%+C~&r5&1,#@Q%:K$W{B(c$XG#;Q#`kE"},
{373.89,672.89,4.84,6.03,"{G#Mk8O&3x$(x1&mq=QQ$Lf/n(-,u#S>#'d(fY####A##v@/Ig5CS*}E>Mb5f37q44yBReRE9DRyQ%Ye/?$%jw/###A##W%.$@(%n*(0,8?&R7+L=25+>w6KW?RvG#}#'tF3=B5kG$ZG#.H$3Z$rb#[#$r7.Y>$###Y##yaC|%095#&##?f)G?&hH(###(##"},
{361.62,718.53,4.24,1.36,":3=uG%###=##1f1Y-)x?&/S+A5#A[&#CR).,95#1##)KGo7/t|CgZ&OG#V5#400$DRLR+j>${4F8.NF<@cW:K-))##r#:S@R-$'<Z$bu%[u$6v'V8)W`7jd*,CRFI+TZ$q/+kh2`K2oGCf7*###&##Fl$-u#95#&##Q7(Au$WA)$?&O-''Z$aAR@I+1u#>V&"},
{361.62,718.53,4.24,5.94,"(##7Q%D>####h5$gc&Su$Ju$A?&ZP#Rc#dh<eY####c##fXD%##w7+gG$[P#|w(v;6+%*<m(syQ:d)P#$e01C$)z5$c@(b'7k,&.6&L#$G22n-*OZ$gL,|xQw~QJN?}q6A3<jI'WxQ4d)jP#pm*P?'YP6~uGAS+AyQP&/Hp0A##P%BXvQRu%--#9yQr5&###"},
{338.11,791.23,4.79,4.96,"8v%Yx+pi;Kv(;f-Ll%bZ%^G#u5&.##,q38J+*n-'##$6%h9J+2=bA'd+A~@)M>MWZ#fK3tu#6kEE#$0R&SK*oEW2~.cG#oEWIo1P|+O..a5#WwVBd$*6'j##[113%)380(##]c7J@Wcu%=,#5-(^?$}K1X~'1x13,####-J$$6&UG#)e*O5#l##JI,Wn(ZP#"},
{336.86,804.90,4.13,2.38,"###(6#Y'8.,#fT.xn+QGMm>#%b?E&0@e+SS'.,#rl#BOC1c$hP#6,#'q2m,&DM?&##1oIeH%fHRE>#e'-t(,-Q%R>#EMRAI*+]/D>#<R$Je-D~.G>#L?;MQ&EMRig9xx,pu$0V3a)={G@]L526&dZ%.ACN07zP#b[*y:6.,#gQ#'.,xf*Z%/OY?*&16,#DS,"},
{352.55,818.54,4.30,1.43,"eG$Hv&PG#T>#J>#iR(#%,.,#,m(KH$Ew+eK0OG#0Q#B;3*'4~P#=K({$,;5#lw-HeCyY$)c#o[P^B0'w+)C)G$)WB%E[Pb7*###oT#@sEaR(^R,y&(8R)sUPI~P|27x-*Qg07Z$WN-sEG%##3g6*R#]S1P),oQ(#0$1E?v~F8,#~3/-B3mc'Z#%Ff-+R*5Q#"},
{318.62,820.87,4.64,4.13,"###@%#4%-No&h,&%.#7K0{x2VZ'6##/w$/TQ######Ay$;QQ:5#G/%]@,g#&OEBQ5#vP$pn(iTQHV66J*XxGVl$VM0zl>=*DE,$bP#z#$hM;Ke-VQ'[G#:oFQ90z)6|3B7F5C//G0.BJ0f>#=Z$X~*CT1TK1z7)vQQ}k#Ze'Vc%ymGv[-68&MQ'a5#~H($]("},
{318.62,820.87,4.64,4.99,"I#%###VV8-K*9FEvP$Cu#wl:t$*ER&}+?DVQ#6$7?&uc'QQ$jH(E>#wp2se*>2Sxp9I,$.&%1b;:[D=dI5]-e~++d'N11;H&Q#%G>#^W=Y#%|N@MQ'(c#u[%Z.SR&*$p1N/(=^7~/%A#Hw5%vZ(96#*`?VZ%9K4F>#P>#ff(#q7h;.&%+pZ$2:9}8%~?)g##"},
{296.83,837.96,4.50,1.59,"mc&:@%1@BQx/SZ&^#$_K.;#$AYJIQ$O$*:##*v'#0$AcQ%##qA.wBLu7N|n,n[*o59^_;2,#pdQ4e)AB6F##%?&|.#AcQ$##H,#LhQOG####:3@/E795#5##(eQFL.Cp8r##,u#=K$AcQ$#####}H(######`.)A6(######H&*;4A>u$######/A'^Q(###"},
{750.28,843.69,4.67,0.20,"###E,#Y>$###/d)VG#9Q$O-&281###$e%s{VeY####h`.6%U5,#I>#yY$###X>O&##)c$a##NxVp?&J%+V`/fc&y/%}yVvR,[>$&##}k#&##uYO%##-u#0$#gxV-o).7,G$#h5$0{&[vV###hY#sY#95#$##^uNWG#95#R##>6Or[(z70f##D>#z$$[vV###"},
{742.62,861.31,4.16,0.25,",l#:5#8#$###7L9$##Y,%:##8`~;-'D>#l##*r7%b-$.W7,#b5%5,#fY#)##|EG:,#eG$Y##-]~r>$z,'X%#yw/*o%X~~4##kY#}G$95#$##74C4Q$eY#&##db~vx2uG%:##]-'<>52uO$##%##fP#yb#95#=v$>c$ql&.,#T<)#L7VG#:5#w,#4c:($(.,#"},
{603.04,920.69,4.40,4.92,"ed,######&##.&~######/$#I&~######l$#@L:######@##</2######.##@&~######H'#K&~######d(#L97$##tb#~##`~0######$##8+~.,####T##&,~D>####U##w^:###[>$*##o-+######&##p&~######=$#/'~######H$#Pq=$##}k#/##"},
{634.14,935.88,4.46,4.99,"D>##########OlK######-##Za^ZP####s##/'[.,####1##95##########V[U######3##H^^######+$#5]^######>##95##########>mS######O##H_^######W$#P^^######,##.,##########]WB######,##@_^######u$#q]^/,####x##"},
{91.63,75.79,5.28,2.38,"yZ%K-&8=GyY#`*0.N<pV=SH$AH$%J*LTODZ%27$fG$/o.ZP#wG%7,#%8P,H$A6Pd>$xD:Yr+-7+tY#N:PaA-c>$&##jtG{u%|l'###^ZD(l#s8PRZ&j3;c#$9'.Cx-;7P5c$uG%L>#}z7g3=ub####642*c$0m$EZ&DX3g,&lQ$]06'n'Z~0~;=fQ'w#%ey/"},
{91.63,75.79,5.28,5.32,"KS)uYATK4n,%_~,S.,7.(ET2Ci3CZ&uv%ZH(<TI95#.,####Q{<Ch9.,#K##M$PTu$/f-$')gE@WZ#q%Ptc&<%PmY#.6'###v#P7u#h>$V,#.'P5y-8@+gu#^+G6W+XFIgG#A$P^,$Ru%2##516OG#jv%fP#*%Pu#%%p)C^/*;>}Q$&t2XB.W)A^u$fY#E?#"},
{154.39,104.69,5.31,5.76,"t,8A'2e8E*(9uU)@p1DcE;Q%]L(A^9=5#0c#5?#:s=6#$###N|,c8Uzw/`P#LHJ]+?<]3bd%_T,{R.nZ#~y2G$%B;2At>}[*W95Z?(zo4d##W7UI#%%##;v#Z|D###V##c'4ZI*PG#iV)Ig2px5,6#GI,&%#a6U#######%#1OF######i5#qv&Hl%Z>#C#$"},
{680.98,119.58,5.53,5.14,"|P#Aq4>u$3,#A](=o2###$##IZE.,#######`{T###(PA###,@'6_/2@RWG#zzT#~,],%W5#u=G######,##`{T###~R*###:6'3l#gCOi.,LvT$##_?%vW-eL;######i##F{T###E>####Zu$g?&Lw*[6(Rw,PG#M5#zA)$-&{k####L,#4:3[>$###>5#"},
{539.39,121.71,5.47,2.88,"(c#+u#r$#P(:Sc%ZP#<##;r9.,####T##b@U######)##q@I%e)[82Be$ap5GM=QG#p>#5EUEH'###1.%~CUpb####R5#FDUC8/i6)SS#UaA}AU$v&-c#GD,J_<<5#$0-JS,pb####Ul#{R*95#;5#h:#@@U,]3;5#t,#uBJF-)###4m$)[%.,####u>#JQ&"},
{539.39,121.71,5.47,4.01,"8.-######H-#I[V######]0#NHP.,#$##|?$)R(######Cu#oq=~>$###U$#X~Vj>%$##l_%H'4N$*A##P)5Q#%95#&##f@*kf4+?%ub#b'(?`VRh=zG%2[$'-8J[Vl5%KZ%*h2W>$$##E['wb#(&%?-(/1*/^.`~*0$(y8(%uEoc'2H&*4.3bH$#####F_#"},
{599.43,127.27,5.60,5.21,"T(4&-'&Z$;5##i9+l#?q1eu%ABXRG#vm(`J'MEXQG####~##or.DB6###$##@DXNc&}H'p:3[K6lu$}u%{]GSK02$%[P#$I$$<>GB5###3&#MBXz,'###ry#>h7Bf0:5#b-&p#&RZ%.,#oG#If4.u####N)&YlQ######p(#I~/OG####X[$D,$95####|Z$"},
{483.91,138.99,5.39,1.78,"###+##T-&I?(8,#/.(/U$Kf4E?&bQ'qC#_-U.,####P(#1cP###<##Pn,o#'3o-ib4G.Uh,%y0Uay2eS)['2;-(5,#H<+1<C###&##4&+'$(>I+FQ#|(Q.cHK.UGu#h[(,w>^?(_G#C&,6u#######L##-R*######{,#tL:+u####?,#Ey195####E6%_P#"},
{804.42,151.14,4.80,2.33,"%##Kx#~,Q%##'L-u(5Jf4y@%#.'qS-rf0N%*S>#fP#7o/.,#_5%k##]-QM5#x,Q<##{07G{)hm+K,#MfP#(1%##LZ#ssF{k#c5%###}/Q&##}0Qe>$?z5t>#,_29[%M.Qd@-###?##GE7[2B######Gc>###.L.95#/F0,u#UF2PS0[$%He.;6%::439+Tf4"},
{530.10,238.68,4.99,2.77,"vQ)~G#V,%>##A@XD>####~D%`x4P?'-n#^nG###)R$T;2Fw.IB6%##L>#.Q#1BXic&(c$*'#yC;*xAQo2#I%QG#>~?/y3{k#uI,###=,#-l#zyJ=3CfY#}G$Z7&fEV4n-PG#&-'~3,;H&oY#ku&###%##7?$v?'CH'<##)T/gP#,w)Md$P07}>&3##g$$=A0"},
{591.35,237.60,5.03,1.96,"#l#bb<###+##eR,c9/(##/])V,%uY#4m$(B0####c#^.'bG$,f2^@(yb#RQ#pJYiu%R$%-1JA$(@[%=zKT&U;5#@Z$~D4/,#@L7O>#H>#EQ%ZOYt#Mj,%2%&>J(ZOYycJ66&b,%Og*PU6%##p>%######?d'w.)^f3D>#(c#3##z3?7#$###tb#zv(OG#'c#"},
{416.01,249.55,5.20,2.78,"}Z)###<I(9,#2@WdG#L5$m[#}99LW,=g3v6'i##w'5)x'KQ'd%03##hM?;##H@W?,#zv*|-#Ao2n'+@EW*-&U%#<N@Gf(|Z)1&2J##f09Q##rAW.,#Q#$bS%tT/3l$=|*yo4H##3l$Xh#1uOB^9O,#=6(g##@@W###1##}V+m6(6#$U%#7N?[##Fl%#&#1T4"},
{352.04,270.87,4.92,2.32,"`W2###Ol#F>#k-'###u;-F,$######TK(OG#$##;5#qY####tE3###H9*###:gS###}w$X.+fu&###Ck0Fd)######_J)###8~,###/,3~v(NoW###Z?#n_+.y6###bc7^%)######bc7###uc'###+z%rC9wA4'##%$$M@=xu'sH$xq3Zm$(##zG$t$<###"},
{352.04,270.87,4.92,4.71,"Km*######8##mdX######I%#S+KM[*1##qf,uQ'<T3.u#OR'X83######O##AeXHH%D>#|$#zD?CVN^Q(A##b35'gXeY#%##|B9######>##)fXv@(SH(V##tI,:O,2YL*##{eX_%*5?'*$#{&6###$##j5#mEF[c#%x1d##4w-*7#lHU;##kdX,##-H&m$#"},
{374.56,275.13,5.31,1.78,"9w*29*I7X0,#RvIcv)~[,###y>%aG#~vU###&##F>#)$(###mn0&##M8X$##v7X%##?n-^5#([)'##(8XG5#'##nY#L[+####S.###K8X*##k7X$##38/jm$.@+%##19X36%1,#<5#-A/###2-(###5:XUA-XYN###tn,`40z,'###v9XN[&###(##m83###"},
{374.56,275.13,5.31,4.79,"+R*######8,#H8Y3?%OG#@##R(9H53ip:&##_SRvW:H?(+##~[,###=5#O>#-8YNc$J?(N##qn/J]%J7Y,##Q7Yc5#K~0o##L@-###&##zY#W7Ya>#Ud+~##X%/m?#M7Y,##H7Y0##m[-r##A7-###%##L5#XYN###97+4v#c6*,##sJWTd(G7Y*##zZ(K/&"},
{551.60,275.63,5.70,0.57,"bP#$##+3*Q&3K?%,u#o%&c5$:S&X/1g,$/U04l#Sx,Md(76'vv*T8()<UdZ%x.0K,#y'LD%)aeSA['D)8*e?fP#pd&;9U'.+JQ%Z$<^'6Q>#u#'/?#wk6k;AN15Ad'(9-r}HD5#ac$2gIvf6z>%LZ%/,#F5#{k#$##I,#M6(######a##(~.95####'.#B96"},
{551.60,275.63,5.70,2.67,"_>$.,#x,#L10C+GF>#M#$gb6{*DF>#G>#RZ#[,%######1,#KW8T[+Fv%6`/I24-eCR'4UP;n@IK%)bl&{Y#Dm(dG$###$##VI%>A,/M:nP$<i9`tCM5$48#J-QFQ%k,#|5:66'wb#zZ#7S._5%=$#}6+F5#~,Q%##95#IC#+ZPsn-'J%:eB###S^,=b63&2"},
{609.30,340.63,5.52,1.49,"HhO1z;.,#;##`174p*'1:hP#;W4ky4}Q)}?$%00fY####I##KIUkG$a5%]&#7<BT&&ElOgc#zJU|l%Im)vn$'g6######'$#lHU=##2d)c'#UJ15d#{IUR$&OIUoY#Q6'5A&=(<######e##lHU0##gm(MU$`#&w>#G:2ZC55f3###5,#1E+n?*######|c#"},
{609.30,340.63,5.52,4.75,"###wH#Sd+###tY#nr,.06/,#iV6k]14I+p>#NH&+]#e@Z)#####M##>iA###:d'5w&gAZoY#MBZ/v%`08ic#$v'G%#f@Z-#####V##^07###[?(G.$4(YXZ%tmUWu#|2@]n%^#&9$#TBZOZ&###?##{w-PG#Yc&Uv#ma9Wo13{?uY#k~,0(.L5$%##Q~@GU9"},
{630.41,343.18,5.52,0.93,"AVPz70###'$#L*Ab6)###xI(JQ#>$)$##[B5;,#95#'##<B5#eT[?%D>#L(#7i<//Cd?)kd$twB__<]P#al$m%,|k#$##V5$Dc%^W+en0_($hD>kr0<.,;i2j5M+Z$lY#*?8p@.tb####EZ#Bu$HA%cx1_J*{DB6l#k5%^D'<R+######w8'D>#######Av%"},
{630.41,343.18,5.52,5.01,"###3##VQI###HU24v&ymG46&9%N|k#C~Sh,$D>####:_S5H&###1##Li;###.n+h[$r_SAR)g~SNc$2N7S7'#-'###qyL4cM.,#,##c8/9#$#Z$z>#cH;o<Dy(>T>#Jd%j~EYl&E##IC.~x1qb#'##X5$t>%.,#-##Vc#V7-###B$#1$(?e,###U$#PR,=5#"},
{314.73,413.88,4.90,3.88,"0##L'2}A3.,#)##_5$R_0tc(95#%##(y%uiAOc&###4##ICN;'&)vRpb#$##_Y>rZ(bH&Og5aQ(###Ol#UxR}~2######A~BmS+c6*'##%8O+xRZP#%##(k4haI######na=xP$A,$###T(6yP$d#%7J'DzRbR--##&e&%:PGw.###:5#T7$4u#[P#.,#+##"},
{403.38,430.63,5.12,3.67,":$#6=EH?(###P,#u?';*?D>####7##v|5P?(###B5#ed%%YG-q%rbN:5####&ACm6*cI(5.+L?(*##Ig){$R-&1###K,#p(R&S*>Q&??#ClLo%RD>#-##P3.V;A###j##D%N(R).,#0##UuISG#_P#:B%#%R@S0###M##s$?bR-######f[#gG$######6c$"},
{470.69,474.04,4.79,1.41,"|c#BAQD>####P&$XkL######tU)Nn/######jI(x#%]>$QG#0w#@<;wu'###_h(p>Q###)##wCQQ%/###8##PaC5c#uu%w.)1$&BJ+.7+7l#I;;bR-###%9#t?QL5$###l^##FFSA%:^7N_)O-)wb#Ru$VI&p(?######2z#-mN#/.@u$/^#~J-+8DC966v$"},
{422.13,499.96,5.24,3.79,"~##_3Ax$,###:##7R&Ph4r5&###f#$s35L81-u#xb#<Z#.7HW^&+6R.,#E,$3~Elm*f$'lS/J?(1##98'U7REw-ZP#;,#HKJL@)9Q&0$#_7RY7ROG#H##5#9x1>Xc&G##[aCrb#qH&du%C(:^5$Og0%D.k7R+~-)h6c@(TV4Kv(t6R###;##4##s6RA,$###"},
{635.48,534.89,5.53,2.40,"3l$###4Z#e%MU}K###1,#Qp&|K4V6)###A##+##/&,6#$.,#-6')##Pg5vR-L-S~>$*Q%)d#T7-Rw).u#Yc#sl'E7&fY#]E*3Q%3,#6s;j./f;9/$(,6&F|:^6'bB3###1BEe=DYEC###y+20l#oc&xm$`%X.u#l9.~,$7(X|5%`wS/,#Fs7WI(}%X###Id$"},
{420.77,596.25,5.19,1.48,"0##Wo)xM+#p6{Z%^m*hl#p?Fzp2%Q%###9y%KS.eY####a##n##>m<SL:###;9(`dR/,#I>#ihRIA1###4##*z4eH(L#%3##~G#q3/mcR>5#CcJcNBoP$N@#PdR}G%ZP#/'#%94;v%$v'.H#'-'J8%'B4F$%r|H=5#R5$i1%-tJ&##d,%K($<$)S>#zu&u6%"},
{624.97,637.56,5.12,0.75,":b1]$*y7/###uC3gv&7/V%##3fSD5#{?*<-'p>%###N##Ie/e<,:80k5$Ac%{V1.%>TuOhP#x0VSH$t?*pG#oe0%##|P$O#%8@)C,$i'#aYLtI-H$&P&)i6Lt-V8,#Gu#~*.dn0$##K#$7Z#.,#$##IJ#5_=.,####0[#N.V*R*###4,#q(1bc'###U>#Al#"},
{363.53,734.57,5.58,3.82,"2*<fl#H.-*##G+Zm~X^l%su#:T(G+Z>'ZF#$ZP#`$#ky8bG#_A0$[${A2-##O)ZvS1}>%J@$Wd(M|5_)>dl$&-'}G#FH'C##@L7Kf*Zu%}P#i(Z0R(Y5$TI$B7)3L*tK1w@,>@+c5$Y5$|c&CR+~G#jc$vE92_<E>#8,#:_(rc'C?$Rw'_%-SG#.v$qZ'/-'"},
{778.73,736.31,5.10,5.23,"+R#GhVuG%###+=-PC;I>####Iv&###Bl#*?&######3##tH)oy/$>7|-V'##,XWId*D>#9##&5I######^,$D>####0l#R#%zI.9,#.{P,*<nRW###cG#gO/R4K###&##Wm#D>####Uu#`>$.,####?B%e,K6.-###*##+H6n?*###$##d6#95####-Q$K>#"},
{762.87,769.05,5.08,1.99,"f##WdDT,%###Kr*|':######ey.n,%ju&$##uY#wZ$wdTOG#X8+yQ;HlP)##HiTy%0fY#9##U+Ct5%uP$Z>#A5#7?$-F<26'2%-G,#1hT2L2AdT###&c#93+FXH###&##4I#j>%###jH#..+D>####a0'+%R?d*###:,#~7@M$*###0,#L-$=6(###&##{Y#"},
{351.14,775.18,5.32,3.98,"5c#-I*I#$ff1-H%Gm)0,#~K*rb#dG#t$(yT2I>#X$&a%)}%0D>#K5#%c#].-A^777)bG$T6$;~X|Y$h@$CETnu&E5#0$8t~X###SH$FZ&ZP#b|@S-$(@+,##UaXUD>bl$1S'#B*QUJ^_XS-(###rP#56'###@v%cG#38.###rZ83K42u#-##;7(#j:dv)2l#"},
{589.77,776.44,5.35,5.27,"###bH#[>]naIwM:#.'NO<:w*:(R]u%_,%Mu#kl%-Z$cS,M5$$##QiPh;]Pc&S[)_'*^bFE@*B;]Yc&/,#k#$@B3?u#nQ';5#XQ%EJGER*[P#@Q&7d$DJ,DT-fcP#l#%##'y%E//:?&pb#-##O@+$Q$S5$lG#w6+###WG#'w$Vo4######%7#rl&,u####)##"},
{573.21,812.17,5.52,3.54,"5,#zl&hY#/,#K_<f,%qb#U##.8[###5?%[g*8Q&###%/Og[*p?)^6'6Q%e5$(ZLp5%@6&,Q$Q8[$##e-'t'(Bd*&##,<[@w(Ul%]P#j,$7$(BsC-Z$.c$~H&|=[:p5#-&[6$l~+1$=W9[s>$###;5#8Z#nP$h?&-$'>l$_P#;N,n3Eqb#$##zl#xfQq#'0,#"},
{418.39,839.92,5.16,5.02,"tVi######8%#Rx2lc&ZP#$&#tNAzG%###Q##xR-.,####5##AVi%##zG%:-#R:8tx+'XE)?#|ad-.*v5&L,#.T0+Q%.,#-##*Wi###Ic%G5#0V=?,#ORI=T/_nYkY#gQ&I(,9w*6@+###)##<Wi###&##N>#cq>###m?$yJ1l,&VG#Ou#@S+kY#OR'L5$$##"},
{777.32,853.07,5.65,3.55,"###7l#qQ)###+)?6##9J0u##Z7X###r$(0M-I#%###l[Jh$*I>#VA/j>%###w}JWl$q-+/##PAY###3A,8f&F-)$##|CYmH&;Z%=$(uY#*l#P=FHc%'l#B,#NFY+U5@I)O.%*S*V486CYs#%:5#?5#eQ$K#%m['*I(QZ&###'s-uvQuP$OG#iv$)DXcQ(M,$"},
{895.06,899.07,5.12,1.68,"###)##]5$1.*wP#9J*yY$#v%#y/g6)cG$+-$wG%/'#;(<Hu####d##x[-k5%1~)RG4hPOR>#b2Vkp2C7-R,#g%/+M#_-V=,####U##(z9.,#0.,4d#&1Vdg2H.Vo.%-^6ah+:o/~V%f-V(##>5#dG#nQ',u#PG#,?#Dx(E05i,&L/#0L8On*],$@R#F0Vgu&"},
{539.05,902.67,6.13,5.23,"'@)#########,F~#########3E~$##$##VG#d5%/,#A5#x5%w,&#########=D~######D##ZE~######CH#[%/###B##]d)m>%######$##0B~######+&#ZC~######N'#Ny6###YG#R,$.u##########V'U######7##JG~######8##?);###95####"},
{796.24,927.37,5.61,1.90,"eY####*m$od(>c$#v%g/+&d&~m'~u%PD;]H%.,####E3d.,#1Z%7##;7,.Z#37(9h+Qo]dP#_pURn+XV<sG#<Q&###}3d$##_>$YG#>%-%##SQ'q>#>VW1'2-ZOM>#0L4R318Q&###]3d>,#M6#-A0>#$.,#b5#X,%d6%z$+V5$.,#z$+;J*######T2d2,#"},
{796.24,927.37,5.61,5.10,"nh_4,#/,####|c'<J*nG$95#q?%@R*U,#Jl%2u####.6#581Kj_@,#Zl&###gT1M<1JlOK>#C_V6B1IH'h>#l[,&##jG$@5#'k_$##w5&###cM<mG#6^SxI+LT^nY#<7(|p+(n,Cl#Fl%2##*j_#########>j?tu$F['Ol%_S)O$'5H%@6&5d$-%)ZP####"},
{93.03,129.07,6.85,1.55,"8(8cC(4#K|w*''5K@${L;J[$Xp7y5&V,#a{3*7)Cu$?'%JU8q80A<,bWD7Q%{QRk6$4[*M-#]EBjQ'1m&Mg-$y17R&BB1fS,h6*:U#|]5U-(9TRm,%{k#A?$]'N~B3Dl$PJ*%v&zP?ex3B7)pb#t.(/U.N&2C.(/kA<5#3-&EB)~SR######8,#5VR######"},
{505.68,142.52,6.76,3.91,"=)?5J+###2##E(`lG$###@##pq<D<66#$_?$)B)~B4.,#*##m`D######F##y(`LQ%Mu$sP##E=:?8&;4w;<=ZFve.0l#=2/i1>######b%#r&`###@5#$h#m)B###-$$=(6'e-###'##G&+otL.,####r.#r/Z8Q&%##FE(s$*6#$%###x+gc'######?m%"},
{124.50,149.92,5.98,3.50,"###+j395####[T2QZI###46#SA/fS04c#qV6%m%P$)Ad%5i<BQ&UD(eo5###tvQFK-8Q&x.#/bGHN1R'5OB+q*0f_8Av&Hx.'v'v,#^vO###dzQ~..0%+8##=:-<wQw@-W5$Bv&TwQk,%u,%:5####Bt4eY#`n%/v(2p)Mu$Zd#&kG.c#302###7QDS6%^1="},
{869.34,184.07,6.69,3.10,"z1>#########*2g'##,u#(##-&1)$#Gh=6,#kB6d>#,R*y>#]2A#########q1g###Bu#hu#y0:,##G`6]U4.@*.u#Ke&l[QvL<######%##)2g###,##b#$-7UX>$im#Q{7Df-GX?@C-3tB{/5######&##-2g######-##ZJ[~P#P-%4m'_G#kS-4Y<+u#"},
{848.54,193.93,6.91,4.61,"K9,Tq69#$###5l#A9)nB8N>#5?'E,#l/1=<15w*>x0(B.|].fd$7K43@#X-*&02Zz99c$87*WC]Dw**'2gF-EZ%87(?nJp+A'##?T4J7#Jx3l[*m84e##IOC*$E]V>cG$_o(6S+AEB|5&7~'5##NA]'##95#NH#EA]%##oG$&I$EA]###jG#cc#DA]###&##"},
{477.83,219.94,6.13,2.72,"jP#RW0z,%xG%a@Y]m&XG#o$${IZ###^##PE3KQ'###L<-TU8t,%dP#=n#@D>fJZ###x-#*q6IKZ4,#wP$*n$)g7H-#RKZCv%lT&z,'Y##'@**`9^#&e[#|95^NZcn.a>$v-%#f,L+4&lO,##,$%eY####af,O%+M5$/,#s~''U,o7/###_5#(?#:D7_5%###"},
{353.79,253.59,6.65,2.63,"k.#$A.y~&Nc&Ip%du&1%#j6*WH&.,#]'#1C9######Z7#SH(K##ZP#O:#b]6&+3;c%9(#}T6WoR.,#f8#I<<uG%$##_(%*:95,####9(#k84H$)###9(#jdQVPL&##i-&'-:&?&k5#hX5e$))T#/v(o$#H[+>6%eY#)'#M96)9)pb#;R#P7*6[#Y>$<&&>$)"},
{353.79,253.59,6.65,4.65,"}+N######i$#;QQ-$(3,#0h+ve0b7+Cv)}I'H'8KH#~./V##KQQJl$.,#d$#XcK4VQI#%YQ#5Q<]RQT,%X>#HbK.m%o#'g##yYN<$%=6(t##A<B7<,p(?f##=RQ&['+u#9%#:QQeG#du&m##w]7f,#r.0o##LaHD$#.;?N$#+QQ###.,#O%#$NBg5#qQ)]##"},
{402.91,261.04,7.01,5.04,"@j=NZ$o%0###`PG.-$~,N###4TQ###[P####uD;O>#PA1###pFF%##Wv)1,#Q|B&##}7Puu%tSQh5$FZ%rZ$F_4IN8jZ(###Y#N###eP#c,$~M>6d&w]*%IKF`B~X1Mc$?W3*Q$onD6#$P#$8q4uG%###5,#zp0%RGB5#E$(j>%G57=5#ac$Nc&|6&###BE-"},
{394.31,278.13,6.09,1.59,"6,#YS*v5&7/NI7$+E>qH'9@MG,#^e'VxS|.0###5##FJ/###>Q%+j,(vSJ?&C/J-Y@%~.-###7(k9(}uS###3,#'Q#{v,###(@+5%##vS6##0vSX,#J'8z6#I-)R-#&vS>,#B5#Y>#tc(###{v,;$#;vS+$$StL3##Cp6#y%bG$c##;vSrY#D5#N>#Fl%###"},
{394.31,278.13,6.09,4.95,"in/###'##<5#1JRfl$en0###YM=I-$~mU###ynU%##4l$###%e-###+##aP#;nU6,#Ce.'##~q<8,#lnUtY#wnU0,#qu&vG#,R*###+##_P#UmU###JQ%4d%_x4@5#Yk7{ZJLmU(J(JZ$`s2^c&D>#$##F>#boU6?'$##BH$PC4;HD9H#._8B,$g7@P>#nd)"},
{421.37,277.98,6.47,3.04,"eY#######kl#95####(##>8-###%##h5#-S.###WG#W>#T,%O[+###$##IQ#@5N*##n#$>&E<Q&eZ$ri+omS&##Gc%F~#B<ENJ0###F~%_?)@oS]l#:7+U6%P&1@,6{8RAR*_-#4[Rk-%mw00$(###;^O:Q&[oS&##`v(cP#i(94Q$;G40{?F##9I+1y#rlS"},
{602.34,323.72,6.77,1.67,"WQ$ZC~######sABp{A;c%$##>00O],A'7O5#d2;KZ&PG#Z,#:f1mz9i#&%##vB~lc'PQ'~,#cOFJ?$|XGDl#m@SUG#>#$l>#J'8$##`~.o,%2A~%##1?&/w#%D>,##dr;=f,w2C###rY#v8&S$*%l#MB&c7P1A~###=,#<P/Y..95#5-$:z4V#%pb#$##PS("},
{710.75,535.34,5.75,4.86,"F@)vK6######PN;jv+###<##ZJWeY#hu%Ag%DH'-&&/HO?m%1L66#$######l/]{k####A##|/]W>#-f13K%bH(C?#`0]vl%bJ)<R+$##/,#Y0]du&###a##Z2]*m&`$+b##^?&M$&(2]i6*pc&j>%R##tA1BdRZP#####B&z%IO#%OG#yv$@7)Ou$`v(_Z&"},
{440.73,566.87,5.76,3.83,"9##GW=U@-###@,#3[&(k96#$95#*Z#71)Q;:dG$^P#8,#2)SJx%8$SH>#|k#0SGiR-8[%b~.x6+###d-$G'SpJ12I'sb#7l7rd)@$)Z6#0%Sd%S;c%9##Rr0Q::cB5<5#c]2=#$(V2`Q(,6&E>#o8).~'GOE?~-MZJp5$fS+Wc%/&S###'##R##:%SeY####"},
{671.57,592.89,6.94,3.03,"s[(r@/###&##]r^E-)###c##u-Hdu%###7##~#$`S&)m(g6&z~0Ru%###(##jo^###&##/w#5FDlQ(zu#;:0_m*6['G$)Fj)nI,7#$Q,#+]25r^###*##Xv&4p^@#$n>#?R@Bu${#$Bp6:W1L5$###Y##:h:8p^###3##39-kr^0v'###t>#+H$-w'Ev)B5#"},
{566.97,750.98,5.71,4.97,"Bc#i@,pY#D>#i[(n./$##S#$[>$+H#ye/;?&<5#K1'sZS%##am(nM2R(YW,%>'YVJ/2l#E8'h./8,#NIF|J-95#I##D*Yb..Mw-xp&9(Yw?'6'YY#%Zu$R.&O(8)$%5bD3l####$##wr2Pg8GR*IL)q%/4d&g%YE,$1,#4%$ox4#Z#*e+T5####(##Zm&ml'"},
{594.85,757.20,6.59,5.08,"''2xT+G=ZFv(l:Zm?(g$(b@'d6(3l#}9LEd)_5%###=(,]RKhS0Y=.K:ZQ$'H:Zt5%Tc%r6%fB5T#$xz5o,&###$##c$%&04Uv)6f(#S+1z.G8Z0c$<5#B_&Y/1:I'@-(+##*##lP#Dc$xY${v,###&##;9%wx4D>####oI#)6%ru&###%##6##:l$######"},
{341.78,760.66,6.53,0.01,"F>#y,&###)##qb#x>%(##L].Tu%Il#E?%x16U5$-&&MJ0#v&0o.)n+ZP#H##r8Y~c%S[&gCL[-)Mo%B=YV/TW>$o~#r)Bo#%lg4Z$)Bm&;R)B=YedLKQ&)w%&S(.CIP8YM6'3l$7$#:B5<w(?Q%r#&-Z#'I)iw'=U7E#$<Q%n5#r'3OQ&o,&aG#H>#N>#<[)"},
{341.78,760.66,6.53,3.73,",H$}G%^>#px0Lc%e,%cZ$WD8|b#C.)<f)q7.@I(jQ(c5#L?&jx3:[%;#$EQ$hoZ`l%^o'a1LpZ']6%G4YlQMw5&3,#r7*7A,a:7su#jd,,##LtZ:uHo?(|-%Me'dCKWpZc#%D>#b##]e/S#$i@-u?&VH(+##-HHb[*Wu%%6#{5%39+d./F5#4l$/##eY#T5#"},
{796.80,816.33,6.19,5.00,"wc$$v@jWV*SV~&2Y]'1UVhc%{RV?l$k>$`[$:?'$##+C1Gc$b?#jWVPp80,#:S-UE+|uN;m&PUVX5$L,$D[$Dn+_P#V[(^>$eP#i.,WH'8Z$-Q%v#$[m'W=;GM=;5#%##wW+rv*A,$3,#f>#C,$:5#xb#I%&Bv)######B'&+$(95####cv#iP#B,$###'##"},
{796.80,816.33,6.19,5.70,"n.&|&/#N9GZ&{R<oH)uP#|>%>.%Fv)2m#$7+J>#[P#0##)7)=H%(I$25~_=H}P5l~02r467*,KBko5^P#V5#E@):Z%$##F>#<5#=m7l1~^c&l-*o?$VYFP/,e/~rb#TG#>]%5/0vb#.,#+#####(_,:$'zY$T@-Iu#8H%bZ%,bI######.H#od*E>#######"},
{764.75,838.63,6.74,0.31,"###hQ$IJ.+u#uG#AM:4S)E'6$v$~V<%$$3g52f*u?)wP#/Z$]?)B,#M-'$6%g,PN>#{n'wOZ??'j>#~XXOKZr#&3,#f/.Q?'g83_P#pb#+##zLZYd&Z%.|v%x~.rg&JKZ{#&eY#?##2(:oZ'cI,%l####%##sLZiI*[?),##Id(k0)wIZ'##X-*6##)~.86%"},
{764.75,838.63,6.74,3.46,"w@/{#%yH*7##zIZ&##o$){B)/v(*##wLZY@*###%##+w+%l#jU:pZ'pb#=##MKZ{#&7o.bg&Dn-)%&-MZUd&pb#+##1f2_P#p8.P?'[l%>5#~XXQKZ??'i>#oe'wOZ^#PN>#[6':H%R6)B,#wP#/Z$3f*S$)&$$*^5|u$xz;j7)R06%Q#42;US.6#$###uZ$"},
{768.07,885.48,6.30,3.93,"il$P#%L5$Tc$]P#@$%2-(1,#Ju$O[*A5#+R(&Z#z,'8##`-*Pw,.u#x5&Y6$2~.(m#5(;C,#cRW6l#Uc$08AO?(.,#y`*8xVD>#gP#*?&ou$u%/3##lr>5,#,XW%##Q6';e&X&S###]@>e]0(##oR+D>####ZG#le(-%,###,XWu>$^Z'$##,XW<o/`C.dR)"},
{893.08,914.11,6.09,5.17,"^@&PY<R~.6#$m@Nh?)|k#`##px0_'49Q&X##0x*jH)###1,#lL6?],fh)fC9+CY&Z$#x*6[&h2=3x*)AMz?'P.P]u%I#$<$%X`BKo(]j;g.*|@Y($$@p1jZ$R`C%##QF0/'2'x1yP#eu#gT-Q~/U'*&U.sS,g@Y|l$7c$5z%KT5M%&`6(El####[b5=l$0u#"},
{909.73,946.42,5.88,5.56,"Vp3nQ&3d&8g4=3;nP#Pq-IC9.)T>5#g?%ku&V-((##n6%5J.8QJ,Z$lY#wZ#N16v'(#5FR?'zyUOH$4I*N6%/:7&I#A@,%##@*9D6'Qu$KH%;;4b-'^B,dbH7wUh?'Cd%7a:.]1WD&1T4*##K^$k)@$Z$E>#c<,UL:/l#'H%f&(04;g,%Eu$:l#}D(yS3###"},
{675.25,947.74,6.83,4.97,".,##########YGN######0##~Ua######g##NkJf>#;c%,##.,##########G.X######6##*Va######M$#B5K^#$0Z%B,#.,##########V6Q######F##dZaxY$###U$#MlDiY#D>#I,#############aJ-#########=[a(c$###&##:Va95####R##"},
{943.11,972.99,6.43,3.35,">u$######8c#%Q%######ND7#########JoT&?&######S~H8e.######8##}-V###$##2Y2m-+###D##.2V'?&5,#Bl#31VZd*#########P1V######OH$m.V1u#pG#c==l?(o',aJ,5;8U5$#########r1V######$##b1V0,#$H#_v(CH&0?$mq):e."},
{188.07,82.94,7.64,1.80,",Q%'##(V_###(##1.'FR+####-%$7)eY#)##|k#/,#$##T>#o%/6##dV_###b?'1N-([SM>#SoT_o/dl&/Z#WQ'_l$Fu$L>#|_>%##rW_###hu&&H#]3ZqA1gbHXH'zm'5y/*-%O^+{k####T#M###?W_######2####A6?'K5#@H%G_7jQ(.##-9(Le+5H&"},
{188.07,82.94,7.64,5.07,"~o.z6*I>#d[%?h4P-)?5#`5$ddDI?(###'##p2V###M0V####?&###8H$.L*}d'$L2|3BhZ'23VgA1?Q&sG#:3V###XM:$##|G%<5#l#&*m${#'-c#ISN-x-z,PUG#n?'y(-L2V###7-'1##$##;,#}k#0,#fY#6,#pl%0m'pZ(###+##*.'ATT######(##"},
{72.85,88.74,7.79,5.24,"T7.######Z&#03C(Q$j5%[))Qx-E7*96'9:'](PI#%Ml$+##d1=_P#.,#pL#_/PFkGzG%}J$RD9u(8Gm(@~)noLqP$jl%bc&Pc&]G#'%&tO8V`@~e-8J*?Z$iwSKS*]6(a-%zPJ{I';7+j>$}k#qG#Z1,<B3SK6$##C<25H%kvS56%rI&R/*N::dl$F@&hR("},
{57.71,110.32,7.15,5.82,"_@.######5##BBZyY$$##gI#]$NdC2VT2$n%X8-<S,a-)Kv&#95######=##?CZ]l&.,#Z##vO;9.GZH(%##|R=vS.4Q%}b#w0:######c##l@Z0,#/,#aT%DA0x5$GJ)d(9;W1`>$Qv#Xn/5`B######7$#f@Z###$##&},I?(&##%w%g'ZtG$&##o]$-::"},
{826.33,129.45,7.23,5.22,"77'4l$6e*q>%6IFD>#t,&###KVS###O:4###Y6'###nVS###_f2JR&_w.:##:RS}G#&I*K##zUS/,#q5%*##a(3|b#nVS###O[*WR'e-%p/.mFK>l#XG#&'%Hs<nu&###>##+E7:_0HRS###)S+BH&ed%mJ.[$*Wv&lP#0n%x>$i+8###*##Pl#CfFo#'###"},
{672.60,134.61,8.28,1.71,"*n-S>#7La###&##/e&Bx2###(d(V6&|e1_c$du&*##G,$rQ'#'6J##>La$##a$(Hb/m?TgP#|TXdK0~-)Z-%HZ&?,#HI+L#$bU;/##ALa###SQ'MQ#`)[eL7Z5NRQ${6(x`4Ml%^.(dc'&##@:;N,#eKa######M6%>A-5d)95#s8%,?&Rv(###]+9eY####"},
{672.60,134.61,8.28,5.04,"M$*######y;(#@*I[).,#LR#1j>M-)####H#A_V###*&T$##?6((##:Z%P~&pm(@35TbIb#$#`V/^4,-'8Z#A`V###M_:%##J@+`P#v>%N5#y?)mc$%yQxn,p$UoY#67(y`-<_V###cH(<##|G%N,$Au$*##)I)bZ$ZI,Xu$Cp6/,#J>#4I%.vK###95#-c#"},
{848.08,195.00,7.35,4.62,"C~$+cG}k#####Q#'z+%C9YG#0Z%X,#Ro1&W3M7+Oe.D7(,'-Xd#<&24%$VZ'KK2+;:El$dv([^^6$'fw*aj0j,%N6'kX7O+B$##ZL:6e#;n.kI+*(;a##J)<Z~L3r@D,$nm%QA+j+L9Q%N@))##.]^&##.,#)6#{~^$##lY#H[$x~^###+##4?#w~^######"},
{216.90,227.64,8.54,5.73,"n>%OG#$##D5#`}L###$##[T%1h:###''#^i9VG####k1#(n-T%/95####D##vdV###2##=E)`95PG#:]#gE;8?$yG%DM$Td+>f2######+##'jV6#$(##--$8k0)R*e8#Om)F$$X>$]q#}6+PQ%pb#######F+0*6'######?#5A7-)#####tq(fH)Y$#>u$"},
{609.98,230.12,7.63,2.41,"M>#.H%On'zf66l$b[$*T1_h:8#$TJ#.|@9$(&##R##wf)4w-###4J#fu%5d),I*<D%sB8|u&t{B>e%Fo+^%CT5$L5#L=0KD=###qC#.S/###jp4J_$`/42##x2Un=BAn*r(0([%OZ8.nL%n+###Bh#8~/%##bl%4R#zH*Uc#HL,M&195#~##zl&;p-xG%]5#"},
{609.98,230.12,7.63,6.22,"###+##D@++u#tB6.##Av(K,$QIU######C##LC:fY####<##fx4%##Qd%VdGo/4+##D40yIU:IU$##B-%WZ$z':oH&OG#<##XLUN^7@Q$wU/d]/4oFANUp^6?JU|#%101$6$_o1Bq,h,&P5#O%(4q7`.,$w+BQ#QAT1Z%###tbCTw,Z,%>5#*K1H6%9n)Uc%"},
{421.18,277.54,6.65,3.09,"X>$######VQ#OG####'##=90###$##j,#fo4###@5#f>#2-('J/###$##/m#MPM1,#Fl#^/FSu%LH$~`,fdS*##>6'}7$+FHfx3###be$p%0efSr,%8d(|Q%;&0#3.=0S27+#$#R5J8w&nJ3mc'###pKJEd*kfS.,#;w*`P#'OB]G#7E-5*C.##y5&-'$kcS"},
{622.55,306.21,9.06,0.08,"X~)(AP.A.$m'vh4vPD/$(%##yILXQ&nI(NI+'7)?%)GI'+r:/,#bM0YQGLC;G4BQp57##'}EaCXZP#6##ST,O[+###m,#Mg6###~l$=+/@@X>$({m+_H#6BXF#Kdu&.##B]D;v(###lH#l?'.,#ju$iC0tm,vP#0S-fQ$,cLKB1ll'-##f{9Z$*###P6#CQ$"},
{622.55,306.21,9.06,1.49,"Q&0Pe,2J*B~(v:;-%-^P#e6#S?&WZ'###yZ$@u#=c%###=.&Tc%z&Vcv(:R'1',i%VM5$E##5'HmFHmP$D##UU1My3I#%P,#$##v/.H&V=c%ZPKWD=4y5K-#)&VoH'sd,]8#sGNAH$(7+OH#Ew,D@'vj=}f3YM@0,#fy*Q~D%-S###Um'PX.On/###&$&EU,"},
{704.82,518.24,7.00,1.66,"d>#jV5OG####sx.[S/R[(Su$yY$ER#.1Z'Z$(##+T(*|BOG#3].|>6op::,#[1Z8S-~@,[6%@~/]-#i/Z)?%###A##3RObG$<J0{5#i#CN;7B/Z_I&R[+~L)v-*lL%</ZG>####4##j`?($(OG####g%$$*Apl'e>#'8-90/:5#j,#jb<J?(######fI%R4G"},
{704.82,518.24,7.00,4.95,"/w$AXB######ej8rl'.u#^5#9d'9W6]I-R>#n&&xsE{k####{);t7/###(##c]Z-u#*e+Xx%nm+%B(>^Z7R&FeMy8/2]2A$$YtC%Q%######5_Z^,%_R,^5#6x.z#&x`ZQn-4'7$##o8,oQ9U2<uG%'##]v&n_ZqP$ZP#=$#R]-{,&+A+P[)95####pG#Ty/"},
{594.35,756.80,6.76,5.06,"kn-ge);Y[B$)RU[mZ'e-(5R%}H)iY#zKIt6)Qu%###g1-YRK}e1|O.JU[Lc$zT[<u#uc'#$$Yy69u#Jh3i,&######lv$((9Um*6^*ev).^+_S[vb#:5#41%dx2c6'G?',##&##sY#Oc$I#%Z[,###$##@0%,y5D>####r@#ac%IH'###&##2##q>%######"},
{364.99,761.42,8.34,3.40,"G>#9#$&##7I*`>$]Q&9d$OR+]v(v>%?$$$lB5l#>R)5]*aS0K,#W%-D>#$##<D=2K+s5&i,#WUW(T.>V1Z9Hyc&PT*,XWKE>Y,#6^5cH&gY#_5Mb]/t5&[6#,XW,fN.[)<I#PQ$b%@9aD;c$;5#^#%xA'EK4<06>Z%]#$E}6r./Gw'X-'#@%Y>$FQ$_P#[#$"},
{319.54,768.40,7.99,6.01,"%JU[P#######`H'4[)###%##e>$DH'&##f7+G>#sG$iG#=K3)JUuY#9#$###H:7u~.:-'=l#UJU{6)s?&0fD9-'(g'R'P:7OoIU%##+l#oP#]C9#Q$5R&17)ANUEwP@Z%md%g6&YxEf<DT?'&JU3H&eG#:[%od(O.-J@'Y.-F7%IT3yv*W,%@$&Ev'0l#&-'"},
{781.87,795.52,7.66,5.01,"a5#Q~,,%)_m+Io(fU8nP$'##:?%se%``C2,#H].zX7CWA^#$Wu$<M-TrTQ)A%oTE~*zp9-d#n&4B,#v[G3d&ov+a>#hH?JvFYu$^1OumT8Q%|nTt,#9&1{[%k_;.,#4S('v&.,####4m#v;=7#$5d&Sc%kA+Xx3*##aP#H<,6~-eY#%##;H#mG#OG#$##)Z$"},
{559.15,800.45,7.62,0.37,"@5#0?%}l'sb#kP#r06>6&)U5Bl#eK6VH#bJ1Ld$O$*/##^P#KI,rY#zP$FQ$^[VK>#v7(p4[)-'Q5#p4[hdQf,$Mu$X/00u#gg8|k####*##(3[gv&,J.J[%`8.I''70[IZ%`P#yY#yDA+Z$#e+ZP####(##+pKR]1j>%###Fm%w329$T/,####=##Fe.}u&"},
{559.15,800.45,7.62,3.48,"U7-C6'###7##B-T:5#K6%7a0Qu%###iTL?y2###(##g.-eY#+)?oG$lY#3c#y&[JZ%^S-{y&x./N[%8)[1.'###+##kC;qb#$]/1u#$?$,Z$d+[XIR(-'Z5##A(d+[3.W@5#V5$$6$ld,,l#;,#^P#Ym$6[*SH#2o18c#Ep7yc&Cp5cG#W05ZQ'hY#3,#Gc$"},
{292.92,829.75,7.33,5.89,"6#$######%##5vT######a##j%~$##S>#kw$^J1W5$2[%u6(_5%######cl#|%~######X:.^)~=&2CH$s:)dw'Bp6^f*IA/p06######.m#g[HA,$###_[Ew~FwK8=,#-N;7*-Z;@Q?'yl'0I&OG#######*F*{v,###R5$qa,$m(C##T.--L,O5$,n$Md*"},
{548.35,827.65,7.43,0.02,"i#&mP$###-##=dUD>#%##eM(f;ArI$DTLe(S]P#S6##KRr#')S.95####N##eeU(['95#6&#A>I7F.FdUAH#%##ZT#.ZPPG#U[(uG%###%###xBP)>.,#1##x6&4JA*::rb####&$#3A/O-)>,#1Z%.,####}#$1R*###,##B5#fd)PG#=Q%###'##G5#sl'"},
{548.35,827.65,7.43,3.83,".,#+##Fu$[c&.,#$##E#$(S,######E##iR,###$##/Z#c5%jm+Ql#A,$Y,#k*G%##bQ#LSC%Q%###bk2tWB###$##Z%*|k#:o.]?&eY#(##K2WOG#OZ#[e(KjD|b#W3W5o-###'##7?BfY#mQ'Au#c,%`>$W3W/v('##jP#P@@lyNH@M=c$vY#vx(P4F$##"},
{754.90,867.23,7.78,0.03,"vP$7-'.,#&##E[VpY#0c#9s/nC<:n$S{R2MWSl%=Q#B7Ik?)v[,,u####1##fJW>R'mP$O/#=}F;j+SIWG6$:5#U%#j}JU-)M-'$l#0Z$PG#vJESN=jG$6##{$({-<(V=AZ#Qu%Z##8m(AS*?##cG$@u####W,#+[(8Z#Kl%$l#2I'wb#~I(}P$]-)~P#)H#"},
{754.90,867.23,7.78,3.76,"{P$Hc#_5%]5$)Z$_H'XG#d-(oG#uG%A##8v(4,#N,$~G#rb#Om*[,#./1k>#8$SU>#q$&x?<iu&.,#^#6+*@###$##h$'Cu$lm*#7'r[-'##hzW}b#+R&kR&xOH|P#~|WoI)eY#+##%s;-l#:Z$'n)H#$eY#~|W~J1@5#R>#'6=~|WzQK],$f-'g0/qy7*##"},
{522.70,151.86,8.89,4.06,"7^0R_:bu%vc#6]NwV;###+##{LV}d*OG#%##</*V'5[P#D5#(D1&w,######4KVN#%###4##[LVI$'7c$i?'<p4D[(.Z$.p,sXED>####*##EIV######w$#qIV.,####8]#lJ2######'e$'FE###&##q$&WIV###$##n.$nDTNc&###F$#U`<%l####s:#"},
{301.44,204.26,9.43,1.63,"a@.1,#95#*##Ayc###M5$K$#,*E###VFE>m&ZP#E##+{4t|E>96wb####+##SycVG#>u$@$#AYL]u#7M>Rc#5?'J$#P5O%$&>:;;5#1,#6##LzcN6(D>#l##@5Ig&0,R*6##&.,_##%^b1##ay9$##;5#9##Wyc0c$95#k##9ZOn#&>6(~,##.,D%#kRY.##"},
{635.18,235.55,9.00,6.13,"p7/2##Rn)Q~/,%X###<,#`-$~ZPpb####q##q6&mP$###%##c'6&&+xR>V'X<%X&H$z0/uS'P|B1M7###U##[T%4n-###$##oB-c(X3z7;Z$i'X%0.9J.o>#Rg3os9_l%$6$Bi-#o.8l$[#$mR)=GJvG$Bc%-'Xu-+rP#T?%?f2MQ%Gn$Ty5pw.Y6&Du#yg1"},
{829.13,243.54,9.89,3.12,"o%VN5$###%##x#&:n-###$##;Q#V7.###=5#L,#~#&}Y#F>#<'V9@).,####nR*9wHOG#,##~mM7C96##m25,[&cc'#/$u,N;&VLv$J?(###(;<8'(J[+A-#.*V#&1,##eA(`^)PST;o(cT5#%V/##P6'tv(LS.>c$vn**OBwf)HjF`I(DQ&'I#4&V`P#2,#"},
{336.86,264.61,8.39,4.81,"4o2/,####$##,8~######~##48~.##}>&l##XS17##m@.#?$}0:>5#95#-##V8~_#&###+%#j-Q$W?(c$,%#9V:<.*1I+&##DL:###4,#F5#e:~'(4+u#5##[@Cq:~W>$'##BhZz-+wu'###8J0###(##^G#,+F|H%RR,0##m8~g[&{d->##_8~G5#}>&-##"},
{552.00,262.54,9.16,3.09,"t/)^2BSG#|H',0I/7,###D-#z8*4l$###%##P>#g#&.,####fd+,u#Ec#gQJ7?IdD>8Z$MJ&K/Rs[,bu$-C+*v&lI(HC66w*f'24o-vr85HG7tE$]-C?K66'P1R~(2c[+_A'jp523,|5OpZ&hf(S5HM7-N>#c7+,/*|n11,#f>:X(7*c$uY#ce-Rg+KI+yo2"},
{641.57,282.39,9.04,6.20,"LB/-.G6Y@)J,(BUx@)c%-,Z#Cx-]E7g5%-Z#Y1*j%.pP$mP#fR'bAU>c%95#$DUmS0t>$@H%p'7/@'*/&2V:Py4OR'@c$~M3TI+g%.rw$E@Uk@U(c$3##h+7wEG###@-#JN<J#%###~,#)h5'c#cu%}.%2@Ux'8D>#:##6yHUd+###I-$g%(E>####(-$eY#"},
{295.79,468.18,8.79,5.17,"1##K)2=6(###z@':[MwR.&##?6'o$&s$U}b#G%U'##UI,iu#qu&9r(D$U$##,&Uc91`.-8R$ie0xb#aeBjx0nV=Zl%C@'DC1o?)P?%7%US,$j%U'##o%.Rv$y{?4##EN3P6(}c&Lo(^cMIH'/m(I>#5`4:S-,N?g-)@m&/^.';5_04h]6.##sZ(g-$V$U1,#"},
{357.43,495.61,10.37,2.13,"%{=S5$`B1,[(tm(^S+qZK%6&jnKAI)2K2gG$#~-6c#rs?i&1J;/a(;:@*M5#zJ+O[&^JQFl$2JQsY#sJ1)@&x&3Pl$]tB]$**#BYJ-YM<C5#'#Hol$Ss;l3B;{;5K/^g4Fr8L_.OK33h9F>#o3G4,#JmM0R&KIQ<c%6J,j-'A['z$'KJQ)c$oIQQH&Ce-/6&"},
{357.43,495.61,10.37,5.24,"`@,tl&H7QhZ&g7QnP$G@)eR'iA/Gd&57Q(Q%)QJ0.':FH@5#y(<<5#8V.2T11U3sV75`<5K/'F;?aAAYGRQ$?{=J5#|>FXn+zFB3I*GB4u>$4T2w-&U7QaG#'8QRu$-T+XR&Um)[>#B)/:_:k4CI//&~-4Z#Y81;#$#oM.%*amN/?&Rm'tw*fx/fH'k)AhG$"},
{442.69,512.52,10.32,5.17,"Y~.&Z$CdME6&eISZP#<f/Y~*8W=.6%uz9'm'&J-AJ(AcO6c$^3BUG#wW6*R((y0w9/]j@o~.~eJP'4G)>&[$<n-nG#9JSrQ'z`:BB2:'6rG#vg9'x&cHS-##bJSl>$g91;d$Z6(###;+5m3C8<=t&1;e-}u#~~/Xl$s5Dj8/EJS$##i#%5z+)//###xG#eo/"},
{376.14,564.13,9.46,5.21,"BH''##X/TsZ'r.T+u#8v&Vw(np5il%'tEo>%6&1%w&q-TY>####'##}b9I@,9g1'h2f(9VJ-/@E@^6uq<>Q#P~/4,#m/TqH'###<##_PHD>#+C8[d$d-T7,#./TvG$1T.F$$'m'H5#{k9x~0###M##Ui@###@[*+##70Tbm)i/Tb5%AH%b@([n+>v&5*APG#"},
{683.39,572.62,8.92,3.02,"O5#EH%jU9###5,Hbc'.,#4##Ly[+u####0m#7%*E?&xP$tR'&##:6#D)?0d(4SXF>#lY#PJ.Sx[######L9(eI-N>#Xm)+a3###8##Bx(U@TwK8###U##IIMgy[D>####M~%a2:K-'x#'7H$Gl#Fn.eA#.x[1Q%TQ'kL$/y[8|B###A##EX6]V7-u#'##?c#"},
{458.54,574.82,8.35,3.93,"%##@{6A,$###)##h|8G07###B5#Lx'%4?j..D>#Qm'Sc$B.D###*zOeY#`>$4S(XwOZ>$j@,,7MX..F6%DM4.m(BI'a>$Ab6###U~)8R'%c>*$'BR'QH&ZyOKYHB=A1Z$@(3-v&=xOeY#xY####$##kv$vHK###.S#`|AfJ,1,#m=8Rf0g[*@5#D2;$##Fa;"},
{458.54,574.82,8.35,5.26,"?'5###Ei9V7*90QV,%3l#yc$AI(7m&d-QZP#JK4)d$m,Q%##bG$6##;{-8;<9q3ZYCZe/7%)oz045HOECnG#k..yP#8.Q@Q%###?##].Qr5&bK5oS(.-QO5#T.Q:6'X.(F%&tu&###mD,AE?###$##4IO###WR+$##N/QZ#%M/Q###0Q$(I%|x3###K##B^1"},
{732.68,632.43,8.91,1.55,"7e&udI?-D5h:CH$yW8x|DD>#on+r~1###q,#{%01c$###R-#CI#B)U=6(###M:+w&U^Q((##f(UM$*###y,#Y06######36#>u#XK3eG$k>$LtGV[+###dv#F%Upb####JJ#,19.,####E[#D>####D#$E%UAB6###%##MyN7>N######K2',d)######Ve$"},
{746.80,653.26,9.01,1.38,"e>#[Y<{fTQ6)2x%9#Lw5&-##Ne)hz:F>#a##';9Gl%###V###f$7iT'y5.,#2xGMS/W>$0d#vT6X>$###I^&9C9.,####c$#/W7M(<####R#PgTdu&###90#_>M######4R#k%0######?@%5U9###6##JQ<}cT######f;(<81######qU'CZ&%##;5#Kw%"},
{394.17,736.53,9.19,4.68,"3xDqPNcQ''~'%<*nVS(oQh7,:m$L~,yN?C6':u#vZ%Q[))6&a,#'8.Nz71,#l((WRSr::J#$nVS3x1TZ%,@'qm+vY#]U-4]2###{[#fJ2G#$sQ(]c#,C6ZO7ZRSQG#VG#;N*|B2B,$7##v?)95#ml#X>$Ie&wu'%##E>#-D(|I-######}$$R[%Su%######"},
{563.52,772.86,8.36,1.83,"d>$?5#jc'@#$D>#W.$PU:F5#P5$5.&<Z%8C2*l#'c#56$cXH@x+@Q&5%-$##0f2r,#PeZ#?#'_]+v',6&+sRy@+E^57;'E_]w*>RQ&dc'P,#eS0q,$1|@'R(0c]_=KF#$m['b](m^]n?&Tm)qv+;5#:5#5@%]P#1,#6H$Yc&Cl#*m(NZ%#Z$/##Qy7zY$###"},
{76.69,50.96,11.10,4.83,"ZP#######:z.Ru%######$TD)R*###?V.cgT######I)4>6(>x.{k#J>#XC+c;=0Z%###F0,4.W###./WYR&eY#####/WQG#N/1[P#$['Pd$=OFj>%'##f?$w-W###1.WF##A,$###@.W###>w,VG#[l%Lv$EC:{k#A,#o?%y-W###FcP=##L5$###D.W###"},
{76.89,97.70,9.99,5.16,"(c$######g$##:9OG####l1$_^64H&###|z%kyTOG#S&,.~&Fl%###)##Dt1z5Hyu'%Z#h;(D[GQw.$##vZ#szTN5$Cl$Z>#6#$&##md%`pKV~/-u#)3._o08wT$-'-Z#g7'KxT=Q&+##7?%>Z%<6&'e+tv%_H(]u#7+@~,$.a>Dc%?H%%w&H(7Au$9,#KQ%"},
{80.72,149.27,9.35,6.10,"g6L######$##onN.e,=5#{m(5'(O@,@,#i1;KV-I.-N>#>m'7K^######-##`9a###&##mW.JR+$##d7$h:^JH$E6(oB*5B5.~V#########W<a95#$##>Z#v_;dH&P@(jQ'###l6#7|>.$(W2?#########%<a######'##uND(c#3-%o#&###C##8n(&B4"},
{433.43,235.33,10.50,2.95,"(y47?%Bl$3f*WI-###D6#>:2qT78##@I'Tr65,#z>$S8(WS1=SSo5$]c$wB-n)A=c#XAGS03jQSD##4p/Cx,Fc$nc%.V+@OG@USIm)#l#E$%/F;nj9_M?(##XUS`Z&^H'Yx+8]*k-+Y%#<RSn.->u$###ZR$^v&g[+###m>$v|/2&2###>5#tC$QQS'##l>%"},
{460.87,233.05,11.24,2.80,"093<0*W7,`%(O05H7(,$'Qf*?6(###f##W^0'$(&##1$$w_7uM<2=6Mc%U7*gIV,d%}m'o44HrAh>#;/F3)6Y<FMQ#(p.zo,P_;=5#B5#bf,1MVmJ/*-'b@&:.O*s4]|E<?#;JVVG#&8&`38d6)######2V,P~,#v'###*~'2;4j[+###~5$i6;q~2/##Yc%"},
{528.29,251.14,10.28,0.57,"###9##b'2Uv*95#(##px)Ew-bg5ZP#S?#qA-RU2`A-:~*){6tZ(o[$4u7`7.rIX~c$<'(7JILJ,1u#|G5=XC1S'My0ST.x.-wc'*e(I&CAp73LXN$(/30-K-p~,q7*aKJ}e0###&##X*0M`AA##,t9%=?OG#Fm(Qu#AU1$e*1H%yP$e5$Y6'######%$#]w/"},
{528.29,251.14,10.28,2.97,"[*>'$((##~o*2dDFl%###?Q#d6':H&L,$xG$Ll$`l&.,#2##=f/0%*VD8nT38UXH7+VZ&Y0'V05Zz))XA$K+%E?`v'Xl%#q-E8-&i3R^8=#$IXXlM8u5&$H#q5Hv?:G810.)|TXmu$X?%K`1l6)Fc%###T,#pz9t,&/,#@Q#hM;@H%N,$2U0pVX-d(~>$bv%"},
{558.17,252.60,12.07,2.94,"g,#Nz3s`B*m(7[>O97WG#v.)Wl?cG$###B##$H$l>%###$##o$*8l#G`9DjA]r<EI*7(3|tBuJP@@*}Q(S'(4J-'o)sS0ld'u(/xz7gND8['yD8Jk;R*B^6(fLPC}9D[*=@%@IPj|09J/sn(;H%{Z&U|FW>#_7,'m'D>#UZ#,`:V-(SG#<d%HIPgu$xb#^T-"},
{583.53,297.11,11.25,2.10,"88+`(1wP$[P#fB(p0SZP####Uj3y/3###bG#>r8D>####D,#=#$-D(u29#I)|4>1q/)(/0e+91S+u#$##}S&<q;Z5$95#u[#5z5FC07b85_5VC8[Q$P_/m7HD-S###$##@41N%.K,$jY#w5#:'0gi:Rc&&n#I;<SI(36&g.B_f4OQ%wZ&W24u#&D>####1##"},
{425.00,437.01,11.09,4.99,"{J0v6'7*;V['Z%WG,#;90*0,c'WD>####<##*H$6#$wY#:#$7}IL,#sU.C%)iV6F:/'a<we-K*WHI+D>#n##[%*h6*?5#E5#o*<x@(xL<Uc#W`B30&'%WDQ#B&W+Z$Q5$z@#cw.rb#Uu#:Z$c#Ob5$,q4R@%(x14##'ZA,*<F%W###0,#m1*8%-%##K,$Q5#"},
{254.88,452.29,10.81,2.28,"JGD0.,$-%cc%pv$7n,v8R#l#E7RK>#/;93n&j{B###&##k%*2J(~o/9HEi#&tJKg%.rT2b~.np:%##-J*B^LRo4)##OG#9S%S6(@5#-;RFl$A7R0,#`8/oo'(`A###'##b13ZP#D$#cm+BI)OG#G>#l#@Uo2ay9###'7&hH<Yw/6##oP$_-$###3(#I2Arb#"},
{274.17,483.81,11.36,2.39,"c81_5$3g/~2=ng-S2<1U4k,%@V5X?(S_5xC53iA###,c#gq-5P@(e-/(1,m':_4]v(]_RDc$a]RE>#oZ&=f*qw0(##gY#K8+`:-ci>dy36Q%}O=%e,,+=?E:%QP###$$%8-;jI.i##Bv)ml#z%-J#$$_OD#$l~RE>#)m&CJ(J82$##F>#a90OG#j%#'D>:H%"},
{706.08,492.55,11.13,1.65,"3l$)##[P#ec#qY##o+,u#ZG#~#$)h.9x1$##x>#%.>);>0,#xY$N##Cw-(Q#ZI)K#7-6P1l#~0T2N8UQOs,$zG%5~%`.T:5#5H&(##JA-xb#P-)yl#:rS#*?P-RXd$DQHlO;Q5$Q6#]0T[:;_-)###<,#]P#cc&###?n#Gf3QG#%##jO8gJ1ZP####pS*^+B"},
{706.08,492.55,11.13,4.96,"6&)e(56#$###0q*.C84u#:5#g6#7K4{#&###.##QG##m&/,#.UP|i@<#$Y5#*b?1l@#HLMQ$S(Q.E?)d(Yc#aw+yb#Zu%'##J'QU,%^>$5Z#A%Q0Q$h'Qc^3/$Q|b#%S)w458n-rG#T,%I##)HM###9##f-@*}E:5#x>$['*bG$2,#'c#pJ+[P#-H#j>%*##"},
{444.04,514.32,10.78,5.11,"_'5e,%DOAsu$e,OF>#[z6T~*w4H>Q$m'4vc&m7,QJ(2PJ`,%(RQ8,#]C1%$%Df.u/+RQHE@*YSJD&.2_:I6$LA1>Q#Z7R)m&F*9}94ln0f,$3h:ff':6RS>#E7R8l#fx-&7%#R)###pE1J+ErNAoR*_%.D6$(x0_,$,j7fC6o6R###T#$DN.rd,###R,#D^2"},
{277.35,519.70,11.25,2.22,"PkEb6)Fv'yd)e%({7.Q@O[P#u-Q<u#RaB;I%#$Q######&-#+%(DS+&eN9I*$'O28.$:59w+3z;3,#Sx-wfJNy82##pb#;A$s.(.%+0/QBH&i-Q1,#ny54(+vDB$##)##rp08Q&T%#7.-R?%dGHIc%?.*JK1[J2###06%PIA.7,`##%Q%C[$95#((#}_AK>#"},
{294.75,551.59,11.01,2.49,"dn-Uc%9L.cp6?<,`_;ae-`P#:m<xY$######xu&$##:5#yb#$h1vd,hD1E$)FD5yH(9DR;Z$vCR95#-6%]6$v/40##j>%V5#*C,QaBq&1AQ%[E7X%.h39M)87?R###1?$sb6wn1g##Q%/c5#0.'Z?'t~EP5$bCR#l#gd(nv&_q<######$&*;c%o$#xz?0Z#"},
{602.67,588.63,10.56,4.68,"<c$g5%5,#jG$=-%M..(##C.+eQ%*6'hG#hR*Ku$###>Q#_?(K>#iO;7#$qb#2p/d@R&##+d'FBRte1]c#Ud>Fv'+@(C=1NbFMl%fS(b/.tw+7:6Y@+@5#RyP9DRn{A/Q$$|3&A'>BR^:39-'Cl${S+;w'BH&;<.B{>G,$m6'>/(JARy>%Q$(Vv%2.K-u#$##"},
{814.16,604.26,10.32,6.02,"9oYJ6%#x.(51_v*gA%`~UJ'.xR,C@):J)8J,5@+cG$7,#c['/tYJqYlI-29&:#H2$@?z8+J*l-*n>$.q-/02Pw*yP$Xl%]c$3JO&J/>##3q+*qY[K6<5#&I#SS(PK2=y5aP#m7)O,$5.,###@^8###^?#DB..q8Sc&+##3e'v&/6l$@H&U#%3-%`P#G@+###"},
{777.30,615.69,10.55,5.02,"###;##@KRj>%pb#.##Mq4]Y@Wu%###L,#EMRY07###i,#<MR###:##pG>2-(/V;$##e/(3TKMLRqC9[6%EMRBu#/E8pS)i6P###;##Rp4Ne,g[+###2Q$(IF%;4Ie.Z5$,*3S>#8x-Ho16-'.,####Y6%RA0Sm).,#.c#O%)M8/H$'~>$@-%P5$%$%'@*eu%"},
{777.30,615.69,10.55,5.99,"i;AwZ(Wu%B-#%fY`#%o5%7s)W]4H%%P|BQE/(H%Qm&.8,J?'V>O###P$$J{4#kY4OCf5%|T%Db?2xJ|/4x6&r6*I#$Ow'D.+Oc&###n.#MgYvD;'$(w5#qE6Y8V[I-###77#7A*zI-T-)M>#######][#,06=@,###zd#1K0A'6cG$'##5w&RR)M5$*l#Al$"},
{559.67,709.76,9.77,1.65,"6w+J,$^P#u#$xl#9_8OG####;g)^'88#$TG#<$'BA#LXH-l#.,#&6#n&4j>$u?&*pH67UxG$)<UM%P`H(]?$ZaA|e$o6Ucl#E>#0##Qh6rb#s,&V[$)<U@ZJa7Ucc%=@'NBHKB4^Z$MC4xm(o,&$##=$&:5####$##CB*(o1.,####9c#~^4###%##Id(q5%"},
{559.67,709.76,9.77,5.21,"tY#CQ$pb####`H$H:6D>####s%$v841,####vG#N#%)c#9#$Jl#A.(K.-1,#yR'SKH2~TvG$<`T'+D4?&C6$-p0iG$.,#'##Fk@XI+Or6+-&1&0E?$B)Sq}Eg[T[5$`I&@KGow/eP#nG$K$&']T:.$_~TmP#XN@<5#ox$e/395#%##gZ#4j?###%##MQ%4m("},
{769.05,752.10,10.24,2.16,"_#$+Z#RG#D,$9##oq2z>%+u#=$#+p3|,I###%xS{>$Q5El7%9m$Qv$&J/###$e&E8D&wStG$C{ScbF8B1Dd%P;9hG$B;39W=_~&:-'+8/###rl&W-%C{S$5G|vS_5$e@'^gLKJ1###@##my1@h2nP$;5#&##}d,fY#PZ#PK1(c$###I,#Kq5/v(###5,#Qc#"},
{769.05,752.10,10.24,5.29,"5,#Rc#tc(###I,#Jq5(c$###QZ#h]1eR,[P#<5#'##[(3cG$B##JL2*/1###e@'ipL{vS_5$C{SFPG|u&X-%{./###$x&/$'X;4O{<N):hG$KT19[%C{SfbF%wStG$#e&/&Dk7/###Ev$Pv$%QFV%%owSt5$EcH###>$#8y3d,%+u#8##pq2RG#D,$l,$5c#"},
{561.79,102.14,13.38,4.75,"88%hm+######l'3Z>$###$##<-SI##_aI0##(c$t##?@X###D')k>%######%>9ac'###*##SAXbu#?@X3##W>$%$#?@X###,]-L5$###&R%zaFZP####~J#O@Xj##?@X?$#eY#*$#?@X###'Q%.,####iM4506[P####@/$KOIa##?@Xg#####i##?@X###"},
{138.85,122.32,12.00,5.80,"9dA$x1+##pI*B?;m[-&##<Z$jW-6U9.##bG#+(#&HQN&#m[-@(;95#=J$I0O?P<RR,ZK*+o05,;Aq<X5$]#$d')(aCCm%8@*o)C1v#Q].cI'.R)|-$O6JT@)RHQ.u#Gc%go%Iw.###[Q#m/0*A0I,#'~'Xd(L5$C##%/)uYFL@-###U,#R8F%Q%###2I#G]2"},
{791.85,127.54,12.36,4.80,"yw+MR+H>#5##ff1ac'###u##|YNC5#}jJ?##(c$:##lS^###e&,f5%7n+4,#Rj:&@+###+##-U^cG#:[V+##U,%K##lS^###}Q)N5#%&0]I)ONB95#&##hw'%T^p5#^-V>$#cG$}##lS^###CS-qH)kY#27)6%+((:'##3f+:y26$DzB9T##/,#@p%lS^###"},
{499.98,236.85,12.19,2.93,"/$Jsb#:c$t-&%D9eu$5v'Zm%fm)&d'Y,%wl$&?&###$##{e(OTTh7*#e,r,$k@Np6;Lw-2I%s.W4J(1R'PY5g09$##[T)|P@4;7pH(95#`,#%6J]'1iY#|.(82WWI*E?&Ox'rZJ=^/+a>1I'LQ'######tT,dg9######a('=?IpH)###&x#Mz/903pb#hP#"},
{252.19,412.82,12.62,2.55,"1C3[c&tr+(937U+O4G;L1S#%}N@MZ&l8+qs5pL<(##4l$a8&###K>#e|,5T4^V5t6)hTGpc&%VWpb#EH$af(fn0B##7e.}R'###+Q$eA)W>$r#'+##Dw>5<?1vS###k$$/TDBS0D$#^DA=Z####O>#g5$bG$###(##Y?$%&1eY####yP#l'7{k#B%#o6Vq5%"},
{252.19,412.82,12.62,5.04,"[P#B5#.,#mk1T5#A~*pb#sS&gd#4@+###*##Il#OG#J>#`P#O#$>9/0Z%>m#AR'?/AcA3'##5sVp*Bpb#~,#/&.ZP####v>$p$'=CSx$,5,#|(=,t-?/Z_c$q/Zlc%KZ&dm#n]5###-##UZ%Z=H)f*uA21[$oB7JZ$Q+8,&,'2Zz5%~Q&Ev%@z6fY#0,#P>#"},
{387.24,511.98,12.18,1.80,"^,ILd'|R-sZ##f0[$$}-N,~'D4G_G#uA0yK(9e*-%'c<BO6'L#Kp,%f%-8K&_~,Bx'c6NNH%%9NZw(<g7Il#SB5lI$a5NZ$%Dj3LK/oT6t5$ip3)y&}5N%7'Y5NIc#J92v/&(/.yH%J4?TI(jNDNl#`L:ET&&&0)m$'O9Eg0*+<8%)AJ/#%&#J,z%&NbLBc#"},
{387.24,511.98,12.18,5.19,"E}Hs#&5/*vm)g'4Xe,&i7q%-3:10z3s|DLl$MRPl>#f:8j-(evKPv(38.[?&O'4,w'=RP<u#yQP$-&10-C/,(B1H6&wr5u06($Ju[)R;?Q#$$*AX5$$w@~@*WSP8[(vJ,[@)S/-;I&mQPI,$+FEnY#bD9.I(^jA^[(ai=(Q$UXA)/)lFI0l#mEBJ,$a(1O7+"},
{654.32,574.88,12.39,3.09,">h,EOF/6'###UrRS$*###9##XR(zd(6I(ll%Z9/GQ&@?&:5#C*CM,$7K3/['$fX.,####SS$TK5:5#E',b)5[/2###g&+qFEfd,%##]%'MfXDfX###5##PW3KgXY>$Pl#+B,:'4'/)3q1{cJI-&FU7dI%p#QXK44l$[,#e6E?E?v>%~>$00+0%*=g-%v',l#"},
{404.25,583.10,13.48,1.90,"mg2RL3Hg7C6%usF+m${)>5K,UASZP#.,#[H#cR*.,#,##4v&=vQnP#/h7pB)*q4h7)SASJ&.VDSKm(_Q(}>#g*@$%),u#zP#MG?6e)_(;Uc$XN=>.%_@SS$&U@SvY#A%-PS$?U7#v$x$*Gc$>YJ;6$$_:vm&x(=J6%R#Ehd'/CSAR(2J.RZ$5i=###H5#lG$"},
{404.25,583.10,13.48,5.13,"'##wb#OtC.,#Np5BH%]@GV@*<4>Q@)XXF{>$Z|EeZ%l*D>-%Gd(c,%9q8:c#QS/ov$l.S*l#^.S4$&y):>d%7h8qc%.b;/w)(Q%L5#3<<zR)I?(>##>2SVm)P/S;]/{U4_7*@16XT*f-SpY#+##`6''w)sb####[5#*0SL5$6;;E8+w-S>Q$'q9QZ%%'.)h4"},
{401.02,670.26,12.86,0.00,"(##,Z$gH#JHQrY#7#$+(#2HQn5%:5#mI#QKQD5#0u#%_$kHQ+m$g'/<o0*Q%x[K_m**v#6q7vQ)1##gA*ZZM$##Y,$pm&RK5p#&|,%-'/Ew*8IQ)U6Su#/h+Dn,]`8BA/Q,$&##6%'?$)###U[+dc$=;71mGuR'YLQ0sAGd([M(kHQ;5#>#$P$##S-E>####"},
{635.61,146.95,15.02,5.66,"Ad#OvT######^,4681}##D>#8;'ZP#,)#$m(=#####,)#9@,en*.V;ub#C,#_`ZZl&###d##rk1]w.yu#TG#d(%VZ',)#ll'Bn.C#$D[$gg-fRX:5#1,#kU#2]/.$'Z$&A-&+=2V,%j5#h>$102y5&}5$w[&I_795####e6#h?(%-$X#%`P#5@&;//D,$;5#"},
{389.03,220.40,16.21,1.94,"RQ'###G##=|=5'3al%0d'e57G>#<l#k'.0<@###(##18*[v*j$+###B,#q-)XRT2c$bl$dE.2v(0H#beIk{:###L##L=Ar5&<R*###^>$2,#xTT###IZ$O?$csD###YH:zp9######wgM|e2a>$6,#q>%###bUT0,#/,#$##EUTqb#kJ,@6(+##:5#JUT(c$"},
{320.82,465.38,15.60,2.17,"e'6[Q'dV2K$(_J,^e-ROEX#%LbF?Z%FK11v&RK/kl%Y_:Ru$|5%XZ&hyNCu$y#MIR*-<;)&,k:7i.-COBWw,`^6c-)Cg54^.tb#H>#hyNI,$mr;jR+UeGR6&H_7'e*BlK=9-_FJN>#+H%FN-/,#`P#|C4;5#(S-Qu$ZyNAm&.5L]>$uQ':a05?'s$#`?)+L,"},
{320.82,465.38,15.60,5.14,"tb#D^,o#'n.#o,%Om>Sy7bP#r[P>B*of4Cc$W6O.,#$##vY#;Q%h+5H_={b#x>NDg.oL65n'M~PQ-%Ur8-e*]]Pxb#F>#7,#-r>r%+QB1aR)]`=TA-3_76w)#lGV&*>>IX-(f]PaP#A#$yG$[C6+6%jf29Q$EC8wZ%V*B#v$6;=8Q$D'/gJ-EeFkl%x-*Y#%"},
{411.77,479.28,15.35,2.07,"6##z<29x2###gm$a+>F>M95#i(<SR(ZXE;.'lx1Bc$Ky2i%,(?%n.'r.QUG#m.QLn+nW>*d%QV7R?&}}Bev'2z1a?'osD}G$pb#$##_1Qbl&t{;ke.YmC{-(313>f-G-QfZ&dPJu5%;'2zQ&:,#Q5$s{6A,$g%.mc&%1QY#$P*C>?&'C3`&-#'0X@+>^6W?'"},
{411.77,479.28,15.35,5.10,"=K3H6&<g1wm)9C3yx-3_:wl%WnQmQ$<B3kZ&X$MD>#G5#,l#5g2/R%(+EOc$;-Qlc%:q4wA*cdIc$&pD8WA.XqQt>%D>#&##X|DEQ$cL3vH&=5D=7'Tq7?-%ZOD-.%?nQr[)rmQcP#)Z$=I%]x.J8,Jf2z5%r,Pm-%_o2kJ*&2?###Ze&O$E)g7###:,#E30"},
{476.27,531.01,14.22,5.32,"xz9`,%so-&m'%3;|c'Lg0'R(uK0lQ'e5Ht?)4JQI>#6w)D&*%>D/w*Po-2v&5V6T&,yHQ8Q%(JQlG$&e'xh7@]4###)?#i.CD~,(e*C8Q[#%44FVv(~S-8}:?HQ###'c##G4fn0.,####,x%J?(###/@@V8/|*I###h#$HX0@@,.,####(r)gc&ZP#4##{m*"},
{341.47,533.60,15.38,2.08,"rK4=c$(r7Fv&Sz2%['NFEGH%TRN'Q$0h9=Z$,4B9l#1H%Y>$Q].BA,.QM(Q$URN?6&sM;7d&_k@wQ&<5E)['lTNU5$.u#cP#;XBkv)8M8O8,JX@BA.|)>me-[2;X.+LQN0L/1QNI5#&H%{9'r_9.R(fM9U?%D17ov'6QNjS*g3F?#$m,%|=6:$);&#:Q&e0*"},
{341.47,533.60,15.38,5.17,"7#$j0*E-)f/#9Z$Il:~_>zb#~kJL9,2h8mv&|q:xl$0D6s$)~,%.U&0QNS>#DQNO^/2`;s@*br=NS,k<>*x-0D8.K,3OCO$(%l#eP#CTNJ,$MPDp?'p+B_6&F{;%@&>RNfZ%;PJ0Q$?B.kJ,e,%vb#e5KO#$|M?LZ$'RNj#%pWB=6%9)51m&@i9M-&XK3/Q$"},
{494.46,598.96,15.60,1.73,"s[%=r;(S-F#$'##tg1OV6=T38Z#+~,/6%e8/|b#VZ$@<9d6)4##,1MxG%###kv&<`T8[*0,#0`T0bD;#$4##On-###>7'w5&m.*OJCG;@/$$W_90}52;9#])E~T:?%Cc%5K$co4###|>$W6'RT3}l#Wz99o&ky5ld%Hj?a&+5^Ti#%6H&M$$d(:X-)[G##v%"},
{494.46,598.96,15.60,5.74,"UV/Xc&LU+jR,1T)D024F=yH)o0R@l$D8-#1-;@S95#####%$<c$ub#h)*EOHJg04&.]6<Ye.QASG,$t,$:|,K=Ipb####hA$%##Sc$+p(X-*.?&95#9N*KsC`?S$c#uH#=eA(|:jy2uG%.[#K,#%@)>5####H##2&0{,#&[)tu$?%,NQ&Q@*V-'~A)j(=mY#"},
{661.41,638.20,14.40,2.44,";L$9X@:[*###/@;ll'###$##=90######B5#jl&D>#&##jP#'t7'P6o&5,-#UaXFl%###r$#ZkCac'###*##,[%7m'W>$###N'8x,#aI'it1|[Xub#F?&a}0Y15_81A@+46%&p0)6&SQ&hG#M/3###{6#sI>K~0%##HL.8v8md,>,#4*6pB1He/ZG#)l#g[%"},
{854.74,164.80,17.11,3.19,"[M=#########Ur]Qu%###6,#d*917,$##;l#0@(#U1F$(Pu$>h;#########9p]###+##n-%@p3q$+Ln'>J+v5$p(:i5$#8*RT4#########oq]D>#+##??$=}D:$):/%h=F,##iZ()?#k~R[H'#########@hJd%0###(##Yk.zn]2##:R*K$#zn]+##7J/"},
{749.04,196.13,16.87,1.30,"b#R}%#V<EAq$Cc%=+2+4Eeu%9%-5A,MA+z|=8Q&###6H#g&R](>g%#]$Rz%#re1@$${&Ru#%r%RpZ'%d&%N0-[';w)z%+RL6sZ(VD(rh?c5#06'E6$=X:*04,j<U80r,$z&0$d$5*@pb####`c&+r5u/4(c#l,%^j=j>$f5%MQ$O,<p>%O5$MA*fK6###)##"},
{864.06,278.66,19.33,2.99,"uR.######,##v(k######J$#.nY###/6%DQ#*##?#$}m).,#AB6######2##7)k.,####L$#~'`_5%`>#~u#XG#:5#sJ*sb#S^:######6##a*kbG$###p##=NZTd+3,#ic$D.+cd(D/-Fm'HU:######B##/)k###nP#rB&K<;-&1.w(8~)$~(EbETl%j?'"},
{603.73,293.11,19.82,2.32,"G5#P~'AS0###w8/5K3mP$###qB,RA2######{?'[P#.,#*##rP#_9-l|D5,#4S&3_Rs&3*l#]_R|FIxb#/Z#LM7cc'###'##Q,#*W)V}H/,#RK.%r2hO<R/,d[Rdu%MQ$N)).q:B,$###t5#&Q#&V'[$OzQ)g#Fhe.{G$x8+|~RDv)8,#oD)Iv({k#M##/|<"},
{457.58,454.72,17.52,5.15,"#]//n)nx2o?'eBUTZ$a~,-/-5EU{Y$1,#vb#-7'|?*%##(l#E;<3H$c01xm';5CNI&|4GU7*|DU1Z$^c&fQ$/o0pu%w5&5,#uh;b~(zK66@'yGNL#$%]+bBKw?U###Z#$$R;sm,@,#xm)JH$'&1sG$=Q%Ic7Dp8###$##*k0%m(######2p)mP$(##ml$}G$"},
{458.93,858.23,16.33,5.00,"pL6##########Gd######)##t5KBl$xY$'##rc'`?''v'Tc#`L;#########5Cd######T##{mVo>%.,#w>#{~2R5$E>#B6#BU9#########eCd######,##WnZ######/p*%.,###)H#<Y@~83######.##,Cd######<$#ODdCZ&###2L(*^0T_91'.8-G"},
{514.71,915.49,22.78,5.04,".,##########8wT#########U)g}k#.,#Q5#V&3^>$8#$rl#############88Z######(##s)g95####/%%-r<@R)/e(M/I############9wX######=##X)gWu%pb#`-#Mb>nTZG.+2f+############.=F#########M*gN>#+u#'##Z&3j,%yG%JQ#"},
{829.42,942.62,17.30,4.88,"b/4######+##ege#######$#9ie&?&.,#cc#we*f`<fd*|[*gB8######6##wge######%$#T_b`c&X>$8,#(H$6T0Ec%mI'd%0######V##Zhe######i%#$CVT,%###t##if4Ol%####@$16'######0##3NdmP$###M$#z3B*_:###W$#)I)f/U###2##"},
{722.10,962.30,19.95,4.90,".,##########;8]######5###EmA5#D>#p>#uw0X>#N5$F-$############fUd######*##3Fm######*?#c^:0,#sY#ty.############uoa######9##MEm######u##HbA1y3*H%_I(############sPP######G##qEmZP####P##Ho/B$(|Y$#?$"},
{785.34,173.75,25.21,4.55,"###.##wm&(8/E#$0u#Il#K{8p@.XG#5l$MD3)R*f&#:[VX?#`P#p#&s>$:&.u6(FR'{p7T6(kuLQ6$:/24-#;?'_(#:[V>##X>#}~VPG#C5#|%)UdN,@+t>#Ba@vb<ay9F-#P5$52%:[V$##G##R[V######oR#D[V######`$$O^VQu%$#####=z*L@-###"},
{659.17,192.38,25.71,5.84,"~##BKS*[(eY#=`){<G<K#lw0j7%uG%?(#WuQ*#####?(#IOIe#%^'./ZIWc%JLSXI,;u#(~%{482.,7$#=$(7@#j98J1#[U;aQ's,$#wN.I(~ISkY#.l#Sx'':5-7)FQ%am(Hv#dU9km&pv+qb#$##LU+8:8M[+###%H#=`7j]5=5#sY#JI''Z$1Q#H)=V?("}
};
siftPoint_save sps3[] = {
{1048.28,21.99,3.46,1.99,"%Q%######[,#`$+###&##JW)D>####BZ#-Q@######B5#6U6=g8######s##MmU###%?#xN,Fv)###pb4a:3######6n(/u#|)A######,##qrU###iG#kc#RX<###uR?_,$######=Y7###5n*#########gMR#########gBI###qS*###j#%###n?;###"},
{952.60,91.11,4.44,1.78,"Oj.~#&BI&###&z%l1<pK.:5#H5#b[)*`9Au$&##Nx0SZ%}>&D6#8Q&pz%H?(BW/z<G?A$ed+EeSF%-$44,24-Q%5w*AeA?R)/,####8/#~w/fw/95#`1$j2@pgSr#'6M/u~+LS)]Z'6~ABn-?u#S5$6,#D,$)Q$H>#m/-{G%k^)a7.no0%l#pT#@FF]J(j>%"},
{736.10,95.36,3.99,0.55,"###O/%]m+######ub8tH)X>$`>#b0NJ'6*c$Lv&jU48]1Z~*TG#QA#jbN###RB4Yq'kjC_y0)@)_8/$P5@27EK5SQ%GH%85<7c$S,#ueNdG$2X6Q8.e%FIS.616s(3N'4Bx);n-z<12c$}n*8w*Au#sn&:6')c#=o$nO4vp:Ic%lq(E7*_w-jY#|B)Y>$Fd'"},
{736.10,95.36,3.99,3.39,"rb#om*&l#>uAmP#t-(NQ&XuA=r/D.-/l#)t3K0-K,$&R(v'+[c#=a93m(t2:RL7qU,w%-@h36vAoI+5t<[/.GKN7c$KZ&~##?&.a266m&;@(KZ=eB2-I&Bz6bC8^T-jB67)*faI###&Z$OB#S@-?#$kH#L%F&?%eY#u##FJNh5%fY####dG9ac'###B5#5e%"},
{762.48,105.13,4.09,4.77,"S94F5#A~/N5#/SRuP#<o*SH&>R&SD8vD8Kl$)4C0x,Tl%*w#Ax1>5#_94_G#(SR|>#Qo-O6&a-)V&'YKKE6'GZA|I-c[*NJ%gH)###h?&Wm'hQRTZ#}n,sK)t./OC%m/LdEC$,DOu#$d&.G=:x2###&##r9*z&4{k#Nl#Vj/=c$wG%vo%,RR'-&}k#9-#.RR"},
{747.44,107.24,3.70,4.62,"t/4jG#@c%*R%]xVVI'{Q(T,${80Wr5;n*@{6zg*W~/*B*(V:<80i5#Px3AQ#XwV^>#&x*$@&:`6uK5T,>/@'}/.E~L6'2dv&r7/~5#9M<F,#`wV@Z#@A-3R$5A/*g'd1VZu$Qz2r]2rw-BS&n?*###{R*xu$hvVru#ed(L&%dI-/W'AwIMQ%_y21w(MT2VK)"},
{948.48,110.66,4.09,4.35,"L5#<-&gu&F>#mi6}Z(?5#/['*ZF###PN-4r;n?*###hE3pl%###TS%w@/###:6E7z4F,$X>#%2Y{#$Qs=mH$:n-<7$'2Y;Q$0,#?v#d91;$)}cQsZ&Qm$Pg2;1Y'Y/8C:`[#V'/FF-`<F%##S6$'Q%z>#pl'aw-95#$$#0T0[(7jA+LT->e(qn+_v$(<36l$"},
{826.74,148.25,4.22,5.61,"]/.W>$######:?p######&##6f~######yG#Ac%3,#`P#ZH&No1)c$######v<pD>####M##R`imP$###I##$@&tv)gY#lY#uA3Z>$&##1,#o<p95####3##Do[Su%###Q5#gc$=[*###F5#h%/j#&.,####4=pF>####*##w,Q5u#I>#+c#J>#2[)XG#y5&"},
{365.06,160.40,3.69,4.46,"b?#j,Ad:;Au$<O,q98.,#$##zH'###$##wb#95####(##Z>$'K++P/T~Z.l#4bZV@+l>%Y,#S+G######8##OG#.,#jY#TG#Id*[H#^aZ>;7S~Z3,#DZ$m+0UjFD>####z##@#$PG#[P#$##95#>5#41+mFGE-)###B,#H6={u'###$##(?#6l$2,#D>#&##"},
{188.68,164.83,3.71,4.52,"A7)z9*v_>'w&dXA-$'C,$4p&4@(bh8zY$J5#A6%TaB:5#M5#*J(gR(>vOhZ&s]K3J-<R*.[%pv(<%B/mLev(7e(&wOu5%7*.+c$o5#.ZHM:9NuO1,#?R'L:13o2lG#[b6]wOEu$UG#v,$BTIWQ%Rp0Dz8%$'PbH[H(wG$i~%iJ2###.##7T,`>$###$##N@)"},
{156.80,168.22,3.95,5.65,"r[)&H$0VV/,#mUV###'@*7##I*>######3##m5$;5#>A.###MI*Mg+ORV)##kRV:##f~1O@#wV?/,#T>#Me'3Q$3-$b[,SG#l6*p-$*X9m>J0TVI#$}@'yV:ti98/1R,$Jl$4l#:d%Hx0OG#QwCR[+y]#cyUaw)&7*UX.~YK`K.'S.V~'9Z%)@'>]/Ml$fY#"},
{69.91,184.18,3.54,1.86,"0&#;z;(-%###0?$,04u[)b>$$H%.6%?%'Lf1yb#]I*3w')Q%a$#}T8l-#]m+xYBL{AN5#(7)}nXM5$L{+qP<m>%I>#OTEa7,###95#;&#f*HniCZP#_$#]uJoqX*6&H^4B%&km(z#&$sX#c#ZP####o##+]1QQ'###Xw#jA2<AE%Z$1?%t5%^B*L-(EP9mP$"},
{133.18,216.28,4.27,1.08,"rh2CFD/A*^w-Fp-^NVj&0.u#$^3FZ$+]%*.+>5####xB#?m)v,#)AV26&95#HL3eFA,%)tb#+KVNZ%8H#GA,@H'###](#P{?sY#dR,G-$*I*1L8ov)Pu#_o.DIVjY#(##Cp%fS2###,%#Ur5:[(C#$tZ%V.-Mo3###-##kY=$IV######;)&8~/###;##X2("},
{679.66,222.26,3.99,3.30,"fY#--#2&0h6(Y.Z&H#D[*76#[/Z.,####q##Z2ZX&4###'##lH)#W&780c,%o.Z56$bG$[%#Z1Z3l$###p##Qs.].Z######d&48S&pb#?l#~0Z/,####Z##gdE95####+##:x-*6'###2##e@.:,#;#$tG$}.Z######:##gHN######$##NU9######D##"},
{816.86,256.01,4.36,2.65,"###g`&Qu%###a&/P%<eY####'E]x[,###'##o82######'#####sM#`/4###OjEPT%9-(J>#'B]F>####/##dC;OG####&##9l#'(&4$(ZP#^i@oP#~5$r#%7F]T7.###I,#H',]M>######+l#J5#|-$z%197)D>#S##;m('r*1_<Uu%S>#~K#CK~N5$###"},
{570.24,258.39,3.51,3.38,"_P#yh:95#ck@}L0>h8.,#N13]K1:5#.v'2J+/,#$##I%TY>$a>$7/QA5#[~.E&Tlp50,#sn'+%T0,#)$(?##RG#%##&%T###zH&7g4ZG#}J/3&T?7*UG#{e(w$T26&k>%>##SG#f@%O$T###g5$:7+'##|,Ch(=Y04$##FF4$]1e'T.,#N#####|v>%~.###"},
{636.92,259.38,3.89,1.36,"wu'&##r>#9'2_G#^c&t6#)`@S5#vG%Li(h;B######Z,8[?)2-(__-i?(5m&u$&#LK1?&U5$CY8T^87r;~5$T5$<5#{eU$##_Q'Dr+93E9,#1K25^-1o-gc$6eU)H%rh>I7$s5&&##`dU8,#:@+_Z#I?S~,#up9&?$VA.Av#CdUL>#+r?8S$bG$0##MdUVG#"},
{607.37,263.99,3.57,1.48,",R(jZ$2EB>,#zf4|>$y@+l#$}%Z=5#z82E$$%Q%###w&ZWG#2B52,#iRU%?#a`ByY#7);T?$_&ZK>#]C:$[#ol'###p&ZM>#r&3]u#vHQA##k2@7,#5W>dZ#{%Z&##b^9dd#cQ($##'&Z>5#}m,`5#ZXDC-$205A,#*A,8R%M&Z,##A&29?#@?'(##d%Z(##"},
{108.25,272.14,3.72,0.92,"-S,6#$###6##}d+em(~[+8c#d-V<c#0o0)0$*A0[##_.VSd$MZ$bg1OG####.g1p/Vpb#)##30VzT-kI.4$#V~/;;%b-VJ,#tb#iV*281(##dbKo02OG#3/#,/VVy(vn1{@#eH(B1%696>##)c$sZ#:_=Tn%%q;Mc#*?&IV$v6+a/#9WBC~#7m)k##7m)i##"},
{614.02,272.41,4.17,1.50,"<d(3K)Ew+M#$AU_]l%[K3h#$z6+###ZV_>5#PG####xR*}k#zx3|>$Tp4,-$3U_iP#ry5bR$+x1###{V_W#$######)p0pb#EM=S>#yr@MZ#IU_O>#2i>)[#Tp8###gU_I5#######ry5eY#ho3X>#@03hH$vT_*##7o0Kv#Tx2(##gT_<,####;5#X]4ZP#"},
{614.02,272.41,4.17,4.52,"e[+OG####6,#j/`:,#W6)5$#$q:vm#20`C##~~.P[%S&3^>#dK3.,####'##y/`S5#r~2l##n{BPn#00`^>#U)@ZQ#5h;{P#rL895#######I0`;?$=]4*##%U5/o$60`6l#^:93m#Pf3b#$3^5QG#OG####Q0`W5#1]3###Fn,y?$6'_nl%r..xY#tl&W&*"},
{599.65,274.12,4.12,1.52,"ix3'c#iA/$-$y9`L>#<L6V$$Om*###:;`wY#######e%+pb#-;<R>#L2=4Q#,:`5,#U)?zQ#$]2###^:`B5#######}&2{k#Ry6~>#3p2@@%A:`1##.^5Jv##K2)##+:`>,#$##^P#^B6ZP#8J/Su#3g//U5]9`V5#PK27J)/801,#x9`P,#>5#`P#`T6###"},
{599.65,274.12,4.12,4.62,"gS2###>5#cP#6y`L,#od,5,#*^4rw(>y`J,#Fy.Sy5l./Wu#/05ZP#$##_P#9y`@,#M81-##m07[v#ty`4##MB2h[%L06`>#,p3eY#######hy`F5#`&4###<{>Sd#Yy`6,#u_=/H#5V;aG#m~-OG#######Qz`-c#F7-###N94g$$Wy`N>#r/1qu#Do2xY#"},
{331.69,283.16,3.70,4.11,"H?(t&#$(;p##g]6#7%-o1dA#Hd**d$fGH|B.h1>_##YH(RL&h&5c%#.95GA#-B1{r.T;<Dw%_wMqS-u?''1/kd,e&#6U92J'Vp9U##x6(mL&8H&Oc#,k4(vMk)A6d#_7+].Eol'c'#ftMF,#}I/###'##5~#95#>##KZ$'7+###c'#Tq=zH)###c'#qrD###"},
{331.69,283.16,3.70,5.85,"_[#qg;:$#u]7ql#p$+Q&#[kMR5#06'b'#[kM######J'#[kM.{0z99[,$y_9iK,EoMy.,G-'OSFH{>]c#K~->R*###8%#+2?ZR$<U9`Q$_h:.@&(:5q2._W@0mMZ$*Ic#:b6b7-j>%&##{P#r##Ce/;%#^rCw-#Sz=o$#<D?Pf%9aGj##=A/@~#@iB-##4l$"},
{313.13,286.17,3.81,4.64,"&V<&##hD@$##?h;b,#A7UJ>#D7U5,#`7.6,#DZ&###/7U###P2@'##n4K%##f;@Z5##9UyG$$8U2,#o%.6H$Cm)###~7U$##/*D$##wFJ2,#K::k?&c37D.+y5J'$'km(2^+zZ($##:7U6##?::%##oWD###pT4Un(O3AeH(vk=(w*r$*|~*]m)$##v6UC,#"},
{269.81,289.45,3.51,1.58,"2SU@5#^[,###%I)Ol#@RUP>#NSU{b#l]5aZ#Bq:###p(<5,#dSU%c#kR-###[w.fH#~SU7H%uRU9,#/B2^%&%_:.,#okII5#CSU###cl%;5#>e.w#$]$<^e/_N>+n*3A*eU1cT3]P#EcM6##7SU%##[P####-/0In$)uJE>#vN8dg1tB6TG#-q5sb#D;>###"},
{269.81,289.45,3.51,4.42,"B1<^##&C9{,#U%/+Q#&?F'~*3=HZG#X@*j2.@n.D$#QFJ`##F=Je##>y72$#pJ.b(/zrAWH&HmAkR+W?(aw'}u&f7#R5P:##zNFX$#mf6s##P]3nU%R5P7$#b6PFZ$>n.?B#hc'M%#T5PI5#5'7c##8iA~##YWDl?#VWCt@$]OJK##c/4iJ#95#r$#S5P2,#"},
{284.54,290.59,3.43,1.50,"abL###kx4###Yn/$##J.TP5#_-T$##<B4m,#fD@$##f^:$##d-T;,#g]6###8A.J?$Z-T*##P/T[5$DU9=##_3E###e)B###v-T`u#|&6###o%0RI#N-Tiu#^-T@5#wo4sw#RrA###c3EB5#K$TD5#RA/2,#S..[,#QvE4n++uLG>#:A/5e'Oy7####*C'##"},
{284.54,290.59,3.43,4.72,"uNC*##Oo2.,#I.+;J)^uMSG#:~ITI+&%,Q,#3e,/,#<.T@5#@ODN>#zq=###df3a7$y-T5,#Z-TXu#'A/jQ#}&6###].T=l#5r>###PsC###Zg8?,#80TK,$7.T'##@&/t>$MU9###:.T'##{T6###7sB###`B5F,#w-T###9/TmP#J'7###Lz:###9N@###"},
{321.07,293.79,4.00,4.63,"cG$###kbM###^*F5##utL$##[PK{>#9cM+##/cM###Oy7%##M5$%##fbM###kbM(##WcMC5#5WAa5$)6?:f--cM95#/L8ad%,u#$##qbM###%cMD##gPM0,#9q2rL1ZE?od)pX;HH'dT6k-%PG#$##`4K###IGM'##FsDrG$`J1Au#PC1PkE]~/E>#?&2{p/"},
{275.83,299.07,3.98,1.61,"(e,I>#GmO0,#AnO'l#WM=5,#KcK###}h=###|;?###D>####S/1Gc#LlOK>#LmOzY#{jG(7$9mO###7cJgG#FmO###rb####=19h>#h.N36'8mO1l#~q5|].B+FE>#PmOM,#PmO###c#&###M)?x,#xf0(Q%;y-7z1b'4~[)|h7Vu%%kF0##9dO###2m(###"},
{216.26,306.70,3.89,6.19,"0v(######$Z#wlPIZ&gP#H7$=:/t>Pw5$}@+fm#IlP4,#F.-JJ1######r>#(mPAu$:Q#'z+.I(qf53T$vlPS##,WAVH#^lPM&3######W##ToPv6+%##Z1-GI&OlPB##/mP###-PK2##LlPlC=######s##>lP######$O/2H&OG####dmP###ZP#(##GlP"},
{463.55,306.99,4.71,4.59,"h,&R,#0EA$##lh?0##v09C5#$8Z)##g7/)##siD0##V{B###5?'F##pGN$##RNCs5$z(>iP#F8Z9,#/&2,##D5N9##0+J###wu''##)$GgH(fV@>5#Ki<RQ%o7Z+##PA2k5#wbN5##3=I###d%,o#'];.;.,_I-bG$`(*D-P[7Z%##QS0,_03NB+##CiB(##"},
{404.44,314.10,3.79,3.65,"}I':{;OG##?#mM,-/W###.##2Q<-H&###&##d#&######K>#C5#8P8eY#*##MmEG16###[?#C3WOG####.$#~.-UH%###)###l#^$&D>#,p'J7-SG####TP3E>MW[&###<`'qZ(2a/###5Z#.,#K>####@l5.,#u%(###Db1.,#G|.###)x&###V/(###T:+"},
{404.44,314.10,3.79,4.21,"cn-]N;.,#Q##Z^#$jC######i7$^#&.,#$##4,#N-#%~.$##Rl#;I=9@,###K|,<ZIT,%,##EyItH){k#;##_-(XS#D_>.##8#$@K$N:<%##EwT;H$I#%k%#VwTh,#tc(i'#R~0C/#'ECx##fu&}H#gS2>##6vT1,#D>#G(#|+N^$#$x1g&#]x4x$#.B5f$#"},
{394.43,319.94,3.80,3.28,"nL/ZB7###Q&,HD-eo5###r91[d)######&f-###1##^I*D6't8*JnW###*##PsWj1>###2##?-Q######V,#OG####$Q#9c$k[,Y6)###=dC$nWeY####0k/[uQd5%###;-#.,#<[)$##TG#.,#2d(###crWam+2(7###?%>eQ(PmL###%-$###H.O###lZ'"},
{260.15,338.15,3.49,4.64,"=7'p*A;l$1u#ua`dm)rb#E5#YOG###.,#8,#.,#$##W{A/,#}w.vZ'S?&v.+C^`7,#k>$.R%X]`###95#,##95####B|D$##9m)###[x(5s@{]`###3d&hc$W]`######-##D>####Q)A0,#0n-'##)oMd&1._`+##>18cG#8AZ######)##D>####i98/,#"},
{238.43,342.43,3.43,4.60,"@/.hu%ll%;w*LoZ16&Xu$j$%UnZ###%##^e*{k####Qe,[~0aH'%w+Sm#*JRIoZ^Q'-[#O,6~nZ######m$#V,%6,#8f36,#B.&'eRG#$j#&7tZ%V7D#$NZ$zoZ######=##L5$(##Kf43,#KR(rm*fl%6['.oZGu#u5%$I%SnZ######A##6#$###YS11,#"},
{377.78,418.58,3.78,4.60,"YH#B}<^WAj5%7O+pz<vP$###uu$vb#eu%.,#&##]P#rY#bG$Co,jq-LC`v>$cH`jw,JJ0%?#4j@G,$kY#;,#D,##l#nY#{k#rZ'Qm%%H`6(7=B`0,#P/+W/?}z?3,#%Z$H6#.##:#$8#$###$##WG##K)mf53l$%##IQ$bF?A,$3##7Z%;c#;5#N>#X>$%##"},
{124.36,451.91,3.89,4.53,"4?#Qs;3:8z>%gV)1)>pb####F?%L#%.,####%##5u#E>#.,#4f,s*0fp`EH%Dv`&90}R.X,#^*C95#$##uG$/,#/,#YG##Z$s$*&@%6u`3U3(p`6,#YA,-I61EC:5####$[#<5#8#$8,#8#$*##~5$<C.A:92-(###z>#J><a5%vb#gY#g5#$##7l$0,#ZP#"},
{827.95,488.57,3.43,2.84,"ZP####m>#/H&[&Y###L5#|P$f&Y###uST'##pl$n5%o'Y#########=##|d,kdU###C##JI)g&Y###9QI8##3c$nG$U*Y|Y$95####*##}-)~_?###/##,9)&&YPG#y]3B$#PZ&QG#*+YL[+Qu%######(H#iy7######Q?#w>@&?&.c$:##]Z#7d(.25{k#"},
{452.17,491.66,4.19,4.75,"NHC95#,##nR)%5;HI+RG#,I'(H>u>%;5#WQ&Fe+5,#w>%###?B1G>#r?#6O=QwSO5$@6$*69hi@;5#pU/bOA}C:-##K$)H>#yD,aU::c#cH'C{S&o0U?'/7&7j>pA-)wSMZ%0M:YG#'d(###@A'Q]50,#SG#m[?;2@95#%##<R??-N:$)2##03B(l#eY#,##"},
{431.19,499.45,4.12,4.33,"$&--'34l$5,#LtGoY####sB-rpU######(T-{sD&##cP#`P#N&2fP#4v&Zy5+nUD>#-##m:/'nU^>$]5#|`:$nU'###Z$>,#Jd*###*8%X#I8oUbG$x6$n_/+nOXm*l>#zA.;nU?5#oP$N5#EW8}c(@R%g[(S=JZ5$|d&@H8w/5?##,{0eO>G{A/##`>$%m#"},
{758.07,506.85,4.32,4.37,"k,&######G,#pRVD>####V%#uUVq.)xd-h$#S>#))%0RV)##Ev)###3##M8'WeNc6*V##/M/jWVRf..K4#6$y,%ax$'SV###*o1###?$#Q>9u/2{k#VS#bZD'kG7H$Be,|f,`>$}G#~DTFl%U[*4,#=$%Z6($L2/,#X>#Z#%gl&0,#H#$Ol$###$##~d$Fd*"},
{336.56,510.24,3.60,4.74,"rq3&Q%A-&}l'Dy+{C6Q?'A[(:.)I[%Kv'km+Oo.###%##;5#9A.}k#K],7$%snX{5&tc$3yF5m(Q>#%S>V/XGK2###]%+Q5$Tg4]v'M6'U>#wqXLZ%+[(Q%'[y4>R%VpXJ$(;e..##s{A1,#aw)491######KiRr1=qb####~o(.qTJ2@$##l@.E$%9w--##"},
{348.67,509.36,3.61,4.34,"%&/v5%y#'?5#:jFuY#;m(y|.|@/?c#cO4|KRrm,I5#Cq1/6&#M7=p6.,#7##GfU4m'3l$,h#MD=G)'3eUGw&Bv)t$#7XEU,#8;70v(###+##DNUB_<pb#B##TB)YqOty:7##6[*mc#5[*u##a$'I#%$##c@)cp0Yl&###^8%@e+,7+###Me#bc';5####P6#"},
{190.55,515.26,3.66,4.27,"[;=###o'&v&2B,h###Fc#yY#Pr@aG#7#$*##/,#-c#yY$$##bz<RH%e>8Lw)r(h4##QS*=v$}lSgG#f6*:##kG$O,$xY$'##]h7v`5_S.QZ%3)h<d%E?'+6#YPM_#$AH',##qP#wP$N5$.,#$&1,?$868#27E(hWG#JH%<K$Jy7Iu$D>#=##cP#7u#:5####"},
{315.70,515.43,3.95,4.32,"y/5iG$m5#<m&EL:###I##2F8(,F###;-#$'3*3=$##.,#$##j09$##t7%)z4h-Ra,%l$%qE2zT6^H%@Y4j]2ttJ=5#[5$dP#On-F,#Y)=b@(l1RM944%*j_3vf1lX;8p0bm*0+H[G#9#$9,#k:<)$&~%.)f%]kLZu#m(/|1RI7-_Q$t/H(M;$2?###6u#i>#"},
{699.12,517.64,3.33,4.81,"P%#J3F######x0#4q<######Yf&[?)###1,#ol#z,'$##G>#TU0l2C%##[7'KN,Cp8C##C939CQ;c%*##0~'@[*Z>$+%&nd)g?=T,%%##O.'qY>95#S##r7M6?Q###.##WQ;~./###o_/se)J*2bG$###6Q#wCQ######il$5CQ######8T'o%0###O(1zH%"},
{699.12,517.64,3.33,5.99,"###P7(######@$&c.,95#;l#8e)W,%WG#~z6BZ%###q>#^vR###D*0######)&';LT###1##?/DN(=###WZ$A=5U7.###vI,YG#7r'Bv)###W18s{2eY#*%#YKT|uO###x$#MK-5JT###$##5Z#Ch&s/5###p:.Oj?L5$7##3-%NLT###,##&-$aIT######"},
{743.14,517.10,3.26,4.50,"OG#######r>$xkN######06#90a######y%#k_@:##[~0o@$zY$######O5#:0a######h##M2a+u####[$#9j=`8*h^;j,#/Q%######1##swZ95#####%#h6atn1###p##_O,WL4(y6###UH(###+##,v#%g5{k#2##Kh-f*9ie1N##WJ.}?H^Q&iv+b>#"},
{330.01,521.06,3.92,5.04,"vq2######D$'.v>qQ(Fu#Zx.rH%ze*/z4(J..&&f>$+Q%,c$HA/###M,#og3Z<0{&5]@&q[+#,FC15dH'<d>:I(B5#pz+d@RQ04Z>$k#$*$&Mf.Du$<o&oQ'9DR<l$D$&hU,W*:QH&9DR(01UQ#Sv'sd,.,#4A#+x-}G%###y6>%^1*Q%###FO0FAR'K1###"},
{661.61,527.20,3.65,2.77,"$m(######Hz'i}L######eT)Ro[vP#sDBN##Ym*J,#bn[%##eY####,##5)0ojH###T5#9A+lo[2,#T*F6##2['Au#+o[#########W,#;S,Kq=###QH#Gd'hn[###pV>aH#.6&X,$-p[A5#######1##l?)'/1###A##^J*gn[1,#$93KI#yG%&##Dr[jc'"},
{841.04,533.50,3.98,6.04,"Hx/###;8.PG#BnP###c//3,#7l$:5#FGC0,#9l$OG##Z$*##>.,'##xU9/,#8]X###5&.rP#D%.###r^X]G#Nc%###t$+&##%m(,##w/2N>#1~X###=w*oK*lA3###MjW/_6hP####ap0###.$((##.R)R>#I^Xv>%Zc%Y2+n2<0g)UaXA3:/,#m##GwQ###"},
{344.71,540.41,4.12,4.05,"6Z%f6%6z8BZ${[-gv&;I+=m#[P#nH#u`Cs.+uG%+##>K)[z6|>$({-[2@XG#a.*i-I@c%P5#NOFGe*jQ(?S'^%/R6#ON@Zf,5##642ZYN###>;5ojD.,#.##M9VGK4A,$+$#grA$e)R/2<Q7(##h?%I8V###RB*6[)Y,%/##/URCB6###b$#7:VF8-K#%JV$"},
{852.75,538.86,3.98,2.41,"###.##Qu%###iZ(&##A,$&$#;n.###Ou#Tm;95#&##YwCEcG###'##Gc%.,#l+K###.u#5$#*rV<:33I)Hj,J$'?4/ToVte,95#E5#o,&%##EFH&##X>$'%#=pTq*4'/1d$#R-(J=/O]51##(c$L5##l#+##MkLA,#A,$^$#i/4(^#2NC^##ZI-e-#I..[##"},
{324.49,543.63,3.62,1.06,"#XC'@?6q7GZ#Jz;[?$jbLBc#X>$W5#DQJrp07l$2,#`Q%[#C&Y<F^.x6P-H#/h2-d=K[+6##2-O)/)q?)U.&/w,5v#Cg7d.)K8#5F5R5P###{F92lED>#*##H7PAw*6#$z##AJ/Z@'v$,jc$_5#a%&%_<K-(Cr1%?&,Z$Bf2TH<M]4%##A?%xq9VR+;5#<?#"},
{324.49,543.63,3.62,5.23,"+##_;)|L<}m,?~&Xl%P#%=q4rF3CZ&+##/E6Rj/{QQhc%;I)>u$N,#Iw&,oQBA1###n&#ERMr-R.,#n##eP<>g0^w/c1%Iq9NJ/###g6#e:4#%*###$h#M::}=AI[*2x(tw.dU'['8w6%eQ(vu$$c#v8.ju&Oo&qc&T%*o,&f5#[~*Y6OA,$(&#G<Azd-###"},
{455.51,544.48,4.00,1.00,"D>#D5#FW<;n,fz>###WG#Lh+RQR######p%#BWAzG%###*##2y4[Z&ZR(~N>>C8###a%%Ea>2SR3u#*##|H$U94K24nP$%##FTRdG$/c#o#%}&1]>$.J&]:8G#Nh>$I5#%D/X[+QA'Fi/fkBP>F###$##lG#Ps>m>%$##XQ$ZQR}k#$##Fq$>}J*7#lvB'JB"},
{204.86,549.62,3.91,0.91,"sC1:'13f04/-{93]wD705_?$;&2?m#E39)p/######N24?90zA0O@%.[N)B/S,ME.%SZN:8$.C8]B'AC8@@$Qx2|u$VA.Fz,G^,xW1.XG###OP9ji2cWE*##_L2Jl:+u#&##7~N|H'B,$R,#lQ#3<:RS0UG#g-#Yy/z&4Jd*3'*Ax06##t5&8AM}k####B##"},
{204.86,549.62,3.91,4.80,"|P#twETc&I@+ar9lu$P5$Wf/g%L###Q##Xh6l^1?5#r~%Qh<H?(3,#w6#sz:%*C###a?#O6E//FwG%}T)pf4F3/F^2U8.F?'DH'0,#.z$ZC:|#O0u#d@&ee*&^/<Z%=/HEH&M%O1d(<e)-E2=e*bu$$z39Z%v%'OK0#M:)l#)B'rD=on,_P#8&O*@(Jv'p$&"},
{291.37,551.77,3.64,2.13,"NQ#tV<4l$'##{e-8:695#e7&H&TZP#(##:A*.7Q###4##.H%=.-c$)_[%wQ&r$TX,%+Q#?u6>h<$c#pB(I:S^ZR###E5#}?%Ee.1##u'3Lw,O)T_H(U,$OK)_r7Xy2jS,jd'|g9Y,%mY#Bc#eY#&##ul$E%T)x/rP#pG$${Rbm)JC2D,$Gn'*v%D%T###*##"},
{339.42,552.15,3.87,2.05,"`F8#.,###|$%%0I|Z)###<$#Gi*/}HcQ()##6I&pe/VkDB$'Ep8###m$#PuB|HQ-u#xG#l_)?q0Vq<`8'nd(|;6D.-(%+R>#`$*ZP#3z#xtI3@K2m(%&)=91{y*Pv)(L+m,%`o+6n,67)[G#M9%>/0o00WH(i6#H90ZIQ)c$A:#><@uv+####n$+4=Am)###"},
{226.97,555.36,4.28,1.96,"oY#v$)`6&)J-yV5Nc&'##n.%p[@tc(%##'A%IV&_&4Cu$L,$######)r,M_6l09###?##dG<4@R.,#<Q#`z-4K,Q?(-B(V7*Mv(sb#.0IfD8wS-6#$Qn%;AR]t?M@,gJ+i]3M9)mM??7'Vl$I@'hV4sBR|5&pE.%R)112tv+Ym$+91s,IW>$H$#3aCR7-###"},
{226.97,555.36,4.28,3.59,"J,$q:%]080Q%oI*6],hQ(bx&yG$<T'.+GMd&###q##e{>L/1{I%-R=s$,'##OC-oIRD>#*?#_'6Ra;>%-wS+U,%+%%O]1XYCA%$Gj7tc(%##d48<h<###2##eJRz~-bG$A7%Z@-NW+^m*P_51e+_5$p7,#$%Fz3######FU)nW;}Z(ZP#OQ:Pl%*e(;V3EMR"},
{207.76,560.76,3.92,0.69,"r6*>9){lP-n&K07W7(|%0*M&pb#o>#=nP0&.######Bw'kHK6+ApD*D1<q,#;D9:wA.H&qH#<jDUT-pR,5U+C6(;R#LA.~38)g'j11OI,OG#c3:ms?###(##LnPM6'PG#XH#)J-{R%3o-Qg3+?#Z,$&r(AmPCr6A,$Cn$.%*amPH,$0,#E~#@S/Y//tl%BA*"},
{207.76,560.76,3.92,4.95,"j##BgNhS.6m(4$'^',eG$zp274995#G##T*=PQ<x5&v5#&&/95#(##SM)uFKMf4###S%#}PJNcN###b[$i<<b1,s-+HF34w,3R);$%$3/tB6F[)nY#z1(fYLw4B8%,mf,+e+L&*:T28i2gG$a[(Z8(6'19c$cg--7)Yg4h5%p,#f'1JcN95#z6#m5J7?'###"},
{339.16,561.40,3.31,0.86,"ZP#2H#GFDzz7OG#$##4n%oAR<V>###GQ#tQ<2?ReY####4(#*|C5m&l#&z^+l6)on$P)=<_5IAR#$%R?'%B'G@Rl^2.,#w##'bIk7(T,%?%#qH)'/$@]3Rg/e}I:['jG$:N-X]2Uk5%.)nm(?N7(L6.,#J##<i?AH%OG#hv$%S-d%)OG#J~%J?(3m%4d%tF3"},
{339.16,561.40,3.31,4.40,"*6'{6#GQHDm&.H&v5#;C4>6%)R*Q5#Lp2(9.xY$$##&]Q6~,@[*ZS&zy82A%%v'9U&4kF}5$rbMRI'Y%.xl#8?'t>#&~QYQ%-:9`,%py/u2+1H&^l#@_QN?')a@Ei73o-ec%@@(v|00f2C5#o#'###%i6M(+###WI#h~Q95#?,#*BHg2B###{7*@5:OG####"},
{450.92,570.98,3.63,1.54,"}?%(m$VH(E>####?.#g2<e./###fQ#`'-('6&##mY#4Q#L@,e^4D#$,c$x,$.kI5,#-f'.<Z2Z%###xKG|7Y######9^,9C7sL80u####dG#e:Zbu%0n+r$'7]1_@'F;ZzQ(i#%*Q$mi<l&1?&/+c$@5#mG#OAB({=k5%XG#0.&`xCG=H/c$v[%Vw)l@-_P#"},
{450.92,570.98,3.63,3.76,"6:3I#%###Z##jo.sQ(0?&}A&*Q%lP#Qx+g_V'e*i-*l9.uA2^i<xY$###t##*]Vhu%###|M'Sm)<E)jy3s-<9#$~y$d~Smc'8<>D>#0##r~-v`V}^8###Y-$Pw'y~C-H&Z>#$##:y$`-*D>#5l$###W##<lK;%*CQ&&Z$TC4fG#B&*A@,)Z$'Q#hv(a#&E>#"},
{450.92,570.98,3.63,4.30,"bo196'1$&Sl$8L9yP$pv'zb4i?)/Q$C57]IJ-m'al$H9/=#$|2:@g8.,#5###pV-R(HZ&_h$:94w1&|nV4A($##e6#nOH[P#3)4&@+###%##P2R<r>bG$=##~8'50G`|G:,#$##:[%8'6G>#SH%5?'*##P92Gz1{,'###;p(8%*8[)gu%|I%$##n5$p^7nY#"},
{420.34,578.09,3.67,4.25,"kI.tP$yH&>d&M97###S##fa9B%Q###4$#*^4FYD[G#eY#&##?/2###l7%]U5G$Q=l$;R$rO37q:H?%Oa0kp4N>It#%#c#cG#R[*7,#uq4yA.uBPFx0vv'S`6;90@=;d90F-(+h6de'{Y$&##n^9Am(Fe*E:)FL:%Q#Ys3S(Q`?(B-%m0MPA1n?(zA'{-+1,#"},
{316.75,584.96,4.36,4.12,"l$+DH$^SQ[u#.C8mG$ue(Z9.&RRA#$;R$tv',2=b/*.,#(##sm*:f$+$Q|d'kq5+17[/+oA-&SREU2IH%el#P80yc;xY$(##gM@W]+Qn-bW,%e-9.%x/H>X8EaDtA)}]2}?%0l#4`)s(:|,'UIDdK1,u#u$#zc&ID*]o3pv%@u$S[%o82/$%mP$_5#<q6Vn,"},
{433.60,583.17,3.56,5.16,"0*5######NH$I|4'n+HQ%C^2Nc#ZJ+*h7b@-aZ#}6*-R)oH)W~*.,#/##tK/Yi.}:=z#$nR+cmQ)'2zc'=d=c6)QG#Vr*`nQw~1E>#Oc#8n(`80S#%-I$c{=XqQ2Z$^l$V9-yv?iv){>9On+Lv%fd(Dw,aP#p&&Ag5kG$`>$Vr*yz4-Q%PG#@`*wmQI?'###"},
{329.96,588.81,4.29,4.74,"*B/qP$l1+7.,wx,XB0Y-'jm*/X0v~+97+*?&6?$>I%1ZLP#%nv(.Q$$JJBZ%{YL$[']C1>,=Dg3v>$CF3E4En6)$Z#yx-bZ&kv)1f$ZPJGH&lxN6B-cT4*d%Qp4Ng0EvN3m&'I*6Q#'V6TQ%S$&zG<?_=RG#PN4z#MiZ('##)cF2vNJ?(&R#Z]5G>#5')Lw("},
{441.33,590.06,3.70,4.40,"gP#f,$h_9p,&.,#%##3A+c92%XB0,#1.'wS.yj@7e+a$'w$(iv'*V2W$*YG#Lq99c$u-'z8+l~.gP#ID1gM0R].@@*$&&pKLyx/ho.ZP####[mE>K47u#N>#9RQV/1zZ&xC';QQ4-#%~M^6=ZI'z,&zu&###tx*U:<######4VQDS0###d##Y[Bb>@B<E^##"},
{193.62,592.49,3.78,1.78,"L&-bZ&&R(v$&7A)K6&>]1LQ&pV,c.S|H'xc&,~(4F8Ln-e&1C;'m(>X>$###$+:&@*Sd&7J,uI-#-'1f#q<;b7+=w,Ge(B+<Hq0q%0%##s~,s/SxY$0##Ae&ZM?Y>$ag$('1eH%0`>7M+Hu$t[-L>#JL)>2S%-S$##0S&@=5F2A###L-#hR'eY#lY#3M%0~."},
{446.76,602.22,3.85,4.07,"K-&FI(oB7$?%{y0XS-_n-CH&1,#&w%O-QzH)I#%,##AB++03m5$/1+GuO:5#;f*']Th,&$##~r@F91ER*M'.Z6)>,#73:H:52,#iM,sZT###}NA_r=(c$:##i]TQe/K,$D$$eC8/%+_p3Qb3H>#m>#u~T###9:1DQ&mQ(%##WwGgB8###O##c]T}-*M5$t:$"},
{618.52,610.91,3.90,5.72,"L,$em)###I5#/E=kZ(###b##0|>qb#0d$I%)~P#7##H1/Vu%R>#:;;95#'##7i3/nOpb#+##.aWb%,7n-i>#N-(|L(2|C0,####;J){?)+u#06Q9x*yA02v#'~W_I$*sEf.#?R&Z)*M[W'#####*##$R'Un.@QR&##Xg3$I%p[W-9$taJ3$#pc$4L'(4I$##"},
{855.93,610.51,3.88,0.30,"2%,qb#Dz3`?)Q0~.,#qG#?5#u_>$##0u#A5####&##cmU.,#M$*$##d}?E&.:/~###{J.VZ#L5O###ZG#$m$###(##W/~-u#7x2%##KdP9[$4/~###?r>i5#4vPUG#.u#S>#'##e5##/~.,#<z;%c#.=C2Q#5/~###l19R-$B<?Ac%95#-##-##eH%y$X###"},
{345.66,613.08,3.75,2.08,"vG$&VPj,&FI'G8(s#R###D^,)R=/1;###J8,lHCUH(G>#.Z#Pc&x>#*$'45:l;A.,#?,#>0J$'52,#;%%B'RQPH@$'v5%2&*v?*####/$a+@i&RqP$0d#BS(4U1_w)me+Uw,}6&*m&S6'[Z%:w+;5#y2+1D8uFC3Z%&*1Uc$gC/O?'^?%v>%xZ&&Z$au#EI("},
{345.66,613.08,3.75,3.34,"n>#b@A~#&###Fq-l?R$l#c>#n%+SJ.k6(5m'aQ&bd(xb#*?%Dg3YO19K5*##tARK%./,#`B%vJ0ow,o7*&g-'Q$5@(H6&UI(B}?a%&o;A&##9DR@h=###Hc#Jy'+@R#H%V$'|Z#g,E+Z$Qe-R'6)##|47{c'xkE8Q&-##+R&aV;A7-###gu9Ul$][*T$'vAR"},
{435.78,614.06,3.67,1.04,"jA-98*XPI<m(###c5#hRPB)<ow0###`l$V9N:QP~P####6]$%C115AJ-)D5#j|F`-'~S,PRCb82>5#hj1+>DVRPtZ&T#$/@&[dJ9U7.,#5##QTP}w09c$WH$9/)nJ/511Y%-[e/>]+P#$P%*l4?FZ&bG#Xe.rX5Y2B%##>l#jN24~+/H&b5#Pp6$&'Y-*g##"},
{435.78,614.06,3.67,3.53,"Am'?w('R'<17Jl#Od'9-&x(=TG#'7&3cJem){Y#*N8BC4OG#GH&XD*t1=uZ(87N=f-b6(|v&fZ't~$qJN3d',##m&$,nK+u#j%'5'Lwu'&##QRANHN*c$}G#V-&jk7^L8Xd*###q&#X=?:m)1B,0D9OG#{##[a>Q/3###(z$r-*~6&A9.BlA0Z%IH#DX;xo0"},
{307.79,625.65,3.78,0.63,"L7)}8*XZ'%##~cQrl#Sc&N7$^v*Uc#B7*IR&vl&EQ%uG$mc$M'2PA-{k#9l#RdQEQ%yb#^8&m@-?.(}T.e8-+Q%FZ$_d'o8-6w+qb#xT&lgQp/Htg:xp*9S.5t9O.Q&-'.6#ty4B:2O?&.q*bB1)`3o%HH/MBc#vI'LhQ3o1Cr5:m)?c$E_3Hf-:q4KH'rc%"},
{284.61,628.51,3.76,1.86,"Af$GD?nG$3,#}#$_k=]R,[w+I80Xe*hc%CL6]FIeY#<##`@'s81E~+0?&G##@nT:I*gZ$C:M~[+IQ$_Z9HpTPL;$##7-%x])+'5@Z$3[#)?MTrTH94Fl$y_1;]+p5>UOAA['PI,F5#gc&$I$D>#O5#tQ#SmT7[(:&,rG$ur=ql$SXE95#3l#gv$'K4###-##"},
{214.19,637.25,3.76,3.04,"F(5^%-(Q$Cv$hm*sY#Mq2{q7`@.%##wl#FmB/,#`P#KX2F]4lz<5l$Y5#2O.;v'H>#c[%h3<)z7E>#I,#<D5t#'[H$=E.D95iBJM..Yl#_&(|u#=f1Df%|R-uS1*<;pG#id(s,&zY6H6(cP#6E*ay9$##p/.pn#B/T0m%^A2l5%$2TrQ&/Q$8l$$?8-Z$802"},
{214.19,637.25,3.76,5.09,"h=Kqb#1##>=0el&Sc%#C$g-S@l#`,%)j,<ED4.+Y~.=h*i.)|-M2m(###Z?#cd']z3h,%QG#Se-&.*O_*($'Ax/K>#AeA2n*U;.=r?'##bl#$i-RJ.)c$Nc#Vr=cQ'1d(Dm#mu%=w%pdSkP#+//&v&-l#Z8(%])WE;:v({#$Ky'5eSm>%.##T>#b%Bnq?$##"},
{192.55,641.14,4.24,0.50,"5y,ZT-Kl%2,#L4J?o';@+^y$@y7]#%iZ&K)(y>$<]1(-&,$&PC.:S,Il%7c#/8Tf?'95#HK#@`?)x.#A*N.%6,#Uf(Lg3_5%a.-uY#~/*7;T0^JKg8$h1)x/:W)87T8#$)##|,#[Y@`82###bR*iK,9C3Hr7h5#@/([QIsV?[x'%S.*Z$0O>5$%3g/Q-)cZ&"},
{192.55,641.14,4.24,5.82,"Y,%kA&TW37o15_:V6(g>#`e)H],XU:1,#;{6>m'm].5w)7i>B$'3E6r~D`L9.|CTl#l//I?:Fv(Ou$sG$O)TBu#CZ$=7%9%TOd)wX;<i4a>$e(TTZ$7Q%H.%i|5je-7v(ld(O?%GJ+M%.[>$cG$0.+F{(V?(aPEgl&W-$J7)fC,T{;(c$'##Qn*ti?IZ&7,#"},
{182.14,644.67,3.66,6.00,"###n.#2ZGXI-*~-Bm$-S(203zo2#[)%##%8'W#$UJ1###(t6+.)ed*^z,L$OIH'e@'aSE&PH{4K6,#~?'XW,p?)nY#iG$.w@m6&k=AZ'2I#$[d)0?:=s?T5$4VQ:c#l,&??#=<0:w,K#%^#$###W>#C_*6#$uP$Ju#IM'pe1u4D[P#:6#;S-;D0Lw-D>#%##"},
{430.45,648.31,3.79,1.05,"9v'f$'>d&[@+c82vx+U$(%J)-H&.-#/u742<$##yY#wg+IPJ0o+d~'um,[>#^RP{@(ul'90%0T2]8&7TP,U.OG#J,#-s<8.I0.*nu#&RP=8)0oFjU6]n/vZ%gHGrmCBn.2I&uL:|v$S]25i6&c#,c#KRPQm'kd'H[(Y/1a/-x6I&S*pI-S5${V8sB+Y-*<##"},
{655.05,648.90,3.84,3.12,"3x,KHMD?'&##Y18&C1405[c%p'(8M6'05$-&}#<(c$###.l#T3B;?%HkGaP#E[O8U&If4rP#ZiALd#|~0d$CG#B.,#VI$9;:4;>)w#T}K*##(*BCS&H?(+-#I~O9##ec$@f)zr=###.l8T[*U94JI$}I/%##k;@,Z$$##_d#}[O###qP#k[&W@M###BW4_P#"},
{301.61,652.14,3.55,1.91,"-K-AZ$4a=h@*Br+:'6,|=6#$eU*LTP$93@?%tH(4~'g`9J2=0o%DdI505###>p2{y2nB2B7+kq;[-*-e#~M7/M73:7vn(iy2JQ%tr3(/0cI,pRPJI+/##gw(aD@{k#@A#&N</?$.17h_(r.0Cd*o5$UK'0UP$OG$##;o&*c>_V@###E##Re'(c$###X8#~y7"},
{399.57,651.57,3.43,4.12,"kP#*+2@B1-u#'O]|#%o-&y>%sP];5####$##r$)rG#VZ'###,R*hP#;A)x@*3vQ###_/#m&1tM]95#)##g,$2^4j6&_5%###W~0F,$Y~*/B*Zy8wu%co%O'2VK]R#%P5#1@%:T1+V2######zH)al%?]/AS-eU9}S-wm(M#$sK]t%+0,#8##el%:~=.,####"},
{321.63,664.15,3.40,1.27,"Dv'qu&^d(&w+paIL,$[c$@O4^Q(5##9^HT}DZP#(##@s;.n(%B3K.$lz</?%3yXR-'E$)vx#-T22;'@xX%R%vm,46#gkLtw&*e)}['pw.;u#Y/CVFDGQ&Y#$.e'qoEMjFnY#2z7j$&,R*|G#Vc&E>#o.)EI'{d)il&H[&fC0>6'Qd&HX6`[*sg8WG#TZ%$$$"},
{321.63,664.15,3.40,4.81,"`(.Gd*[6'%Z$#~&Uz3A~-n6)jm$`d'EI)-R*&@)<c$N[)]n,bJ/h#&pp2lQ%C%V--'Y$%60GE-(2,#)I;2dP&YI-u#*~)`Q%=]2)c#1L3|5%2{V>l$(@)Ie'W033v%8{V;v'u~.<$(?<>O5$T-'F%&f'9:5#jKHv08Fl%###PK(ipUVD?###nu&L$&38/W-("},
{332.70,663.36,3.71,4.38,"{Y$yG#;o0n,%p1={Y#U.,Wt1y6*n#%BQ7u#EN@+K#$[8,zZ%WA1Yg)I%.hG#?oVv6'9-(hC%fg8>1&GoVn%),L7@-%O{?*%$rS-dU1{k#$##]iT<3BL5$8##Yf))BG><D$H$Yd(mR+)J,9[)_l$R]1oP#yA0uz2b$+####q)R8//m(lP#5g&F,$;u#jQ%y7-"},
{233.17,665.41,3.95,5.04,"a.-###La3jn/G_8yd,_Z%X-'|25&p4+l#mo'-V.<f1A,$du#&0,&F8gZJ<5#,/U/[%<~-X@(WK4sS+Yn(bK,#w#?wK0^04u#T'1o1UT,%8##//UD.*hY#%S#Th>5,#f&%2T.###l>$wr*%~.ID>{n)q#'{h&V-UOG####|L#WI-###7Z#l7'######iS%*6'"},
{188.96,669.86,3.98,1.86,"(##)Z#y;7=y6v?%.-'BX>V,%gg*'*@IA-:c$Rm'i$(PU1C@+M>#0r(UvRA,$ug72&*K_1)w*_K3V,%Qq$SA/3o,dB5?/)[v(XZ&bN,I_=V-)t[SF@(@c#tv&/(;rb#Fh$xA0K,#5'4?`)>%.'B3}H&B9'mnN]ZS'##BS%bW3cL<###6$#2J)ZP####L1#EV>"},
{188.96,669.86,3.98,4.10,"N..1,#.]'|x/vlQ$##V%'&w&MmQhR*###0##hP#Zh7###%##]p5(T,3S,^Q%|lQE.(`[*rZ#qV>=d=eY#4##J<<lfM.,#H##(7+-[%,HE;8,=mQrI(cQ'`d%rl&HW+=9.nM?9V3~~+?6$+QNU,%Vw#aq<^v$Gf1RG#c>$Eq'O6(0,#eZ#$qQG,$:5#^&)8oQ"},
{302.04,669.73,3.91,4.05,"T,%J,$&W6wc&%m(tY#|@)&3:&-O###^-#ED:m-O:#$95#>,#7@+)I$|i?nG#,-OBZ$Jx,g2.D;?_c$[*1Up1%-O3v'M>#,Q#1d(*S$u^:?T,;7JXT1z%)hs<+)9.=<j8/1l#,16}{68#$'##>2@kZ$($')a/S6)Nw#%1O)'LDJ/pH&H|<ad&z,&mS(i/,dv)"},
{314.96,674.74,4.21,5.05,"kJ/4,#B6&{J,UE7,-'k5$ci72w$`o/dp4Q~/[6%N%,4g5ic&4@'hY#[I$q`D[E/S1;:g-%w++}Apz7XT0ovFOo/2Q%9;+#GHR%'LZ&s,%i5%7J.*H%Hi147*VDS}G%|w)sh1:)5sv(VDSrJ06~*@c%ZG#nY#W~+0m&Dg4'-&F56WK2L(;.,#AV*zASJ]2###"},
{431.33,676.17,3.59,4.77,"*6'(##vm%F38Yl&u,#XPAK%*s-+)##RM)7uJ~o/2Z$1x%q@,9S'N[({//Ge,`811$%Q~N6?%MZN,l#(V,DU0,8/Cl#0Y5F'/{[#BK1mW;###qT0K+=mU;E##c[Ngn($~-+[#A:2}2-9V:%c####Ac$J'M###SB)Bp36X?(##(3/auCPm*^>#I-$'s.Ve0###"},
{424.97,676.93,3.41,1.41,"^S0V5#K(7,X1WH(om#vNArp,wu'-$#_`9z18}PK###'?#`]1N^NlA+v$+GH$Rp0U0(?ZO&H$<ZOl,#o80eU-&93>5#d:,)'0;=2so2^A0[P#Bt4iO=8~/Y##^[Ow[*al&j?#OK,q'6+@(NH%^$$;K0oU;$##&0JYd+Fu$c5#th6]<7oP$=H#xy+M{7ZP#&##"},
{424.97,676.93,3.41,4.87,"*w,bZ'nZ%];.vG%Cl#DP6nB3h@,0,#V('jFIbq7/Z$vu#=8.K-)4,#yR&Lb>Ad*V##K&N@@)33D(##0z(0IKI@*$l#Iz%vU7LA&DA-|/15m'b%*//+/oOmG$RnOUu$$K/W-&SB1Y8(KZ<K[)w##fd*KF8###/7&.Y>q&3&##NHE98-[?)*6#Gm&M{*{1?E5#"},
{389.03,682.70,3.71,1.76,"tI%S7-)d'?5#4U*Bx.8l$Bn)FR)k@'|?'/T2[?OOG#B5#P>#8f/AT3_l$U,$+AT/Q%*$%-R<K.-Al#nTI4~NEFI###6%(a&'de/3H&4L#<@T|CTG%,t,%-4;_T.eq1}>I%e)G[+'##X6'F9(;5#q,$%&$guQ|Q(i90|c%-YF;6$$=E[P#'m&eI'v6+###<Z#"},
{212.01,683.76,4.09,1.25,"nc$8{6}6)n,&x?*0^%_z:$W8mP$g##gP7}HQ|Y$f>$R%+a3<bZ&]~%}V;(?&HIQ%w%^..fg()f17K%]JQGU1{7/Hu#6L77+2|Z'g%'fFEiu$C9J_^2ax286$X_4(,6k`CUc$0h4|I+zG%[R'oH)gY#HcFr.'ZB21H%,V6g8)PK2_H%gkCH?&7)<[5$jG$~m'"},
{212.01,683.76,4.09,4.77,"PU12u#+U,9Q%zp-p%+r8/)I),X0W['4f/M#%l&.K-&994,x+s?(`5$~r2K[(.@QFZ%%g**Y82^4|b#h#:o(:/ROpl'3Z$J6${c'IJ%~17;7*'CQ0o+2f/9R&5M9WB-EAQAH%1h:QH'NI*OJ*x,%th.=V=R5$:=7B?JSH(H>#X<>&BQH[+h##HGM/Z$Wl$q~'"},
{336.51,686.09,3.83,2.16,"`e-(~G>u$a##{o(;|E%##0v%ig$bsH%##96'oC$M/3###1,#%V<fu#4v'ko)[|EeY#=##IE4T[(3l$MS$9cJr7&r?*z,$s6)?6(###{o#>'S.&S7#${B$WB1-B-Jm(F0.5v'$##?5#[e)NQ'|P<I,$U;(-?G;@'[?&+qNN$)gK([J0Be,5u#>v%K-),-$Pv)"},
{429.00,686.03,3.57,1.48,"M$*EQ#Jg1(328n-Q##,k7XU6F8O###;##jn,hRJ#########1A,O[$>6OId&6|Dj##z05m)5dW?###^g'Ih6D7O###*##8-$q(-$W59(<.##N7Ox?'dR,oc#=(3R93?x*uv'W7O6#$(##yH$?43Sf3`l&2##Dt?tT.T,%;6#/g(I5@3l$&##gxIg,&###$##"},
{227.44,707.91,3.73,2.19,"g&.eoO.,#~##Br**lO###|P#'g#&lO###8Q%&1#vL=###:5#QlOU['oP$T~#uXDpc'QG#g30;%(#@*|R%KnOnI#|J3Vw%bm+#946,#87$.5;ZoO|G%rZ#mn(V90|'0|@+'%*N-(KR&&J(Cd)z~.###L1%{cE7D:L6(ON-R%+-x*my5]@-M5#gS/p,%y5&#$$"},
{324.34,710.63,3.55,1.28,"1Q%wu#//MPw+_SZ$##;]+$x*:z^######<##T1;~,$A?&lY#fP#^n%IV<[P#fN=>H$'916$&9}^E>#/,#o,#%#LZu$@Z%_G#-##Uc$j<A[P#uL:*-&H%+uw&px^###&##C($`OJXZ#$937?#7m#B]+,`@###:/VAm&H?(F##jx^######,$#j(>zQ%;/-%Q#"},
{189.67,715.67,3.79,0.82,"?~*$f(/q:nu%:ECqc#Df3'~$v@-M$)k>#:x-U-%Zv'kZ'3m'Wg,1z1Nm*}b#qmS>$%/$(C]&X~-)A)Kw*{'3&d%%-&tA.BS-73?Q5$cZ$}t8T#AFx1yJ$$FB|9/N%KP?&1Q$9QJB-(nY#.%&W]39m&=n'HDRd>$dS&7rS+7O;v(ju%;'1:30H:;######EV%"},
{167.53,720.62,3.35,1.79,"sA%kw/9u#(Z#5Z#6(3Fm)-z2]A1~>#/@%0L5B.W###*##|>#J<9`7.UG#|u#2/WiG$T[%^0IWm*_G#e-;@%M0{@####$$/()p*Cu5&x##WmLY2WUn-FH&:p)|7*P:1c#Fi@*Fl%###5?%{])L5$###%$#2/WA~-@u$*##H6F`5$p#&F>#wA.#########kR("},
{213.75,730.01,3.99,1.30,"eA0e?(f0-wE9b}EL#$Qj7c6)`UX###/,#'##~XG######=##aI'J)<Hm)'$%(h5h@)TlJ|u%(VX)l#Ju$#$$|RX######3##V[*=-';-&0;3R~-b$*~:,vD6uRX/,#*Z#Yj,R>O######6$#4y/WD36f1L?&Tc<^^3UQ'b#$/TX######oR#?y7######mZ#"},
{693.26,807.60,3.98,1.49,"_l#3'4t5&###Dw(@v(X,$xP$kY#2u#2H$=n*######&##I|B]c#bn,^M=###?j2b<E|Y$F5#+xV^-*mY#Et0D81iP#KH$?|V0,#bu$AoS%-'2+Em-*(@(jB/?|VE-)###AI$wu>;yVPu$^I(###&H$o00a$+PH&n[,FZ#=sAJD/p]5#Q$GI*om#:xV###}Y$"},
{240.37,826.17,3.75,0.50,"m@#Jz]######o3)RrB95####;[$`u%4Z%###F>#ub##l#%##y^/^.JT,%.##T#^RR,###R##u,D'Q%###)##)Q$9Z%###&##_T5BH%h,&8C$Lx]######X2#l|F######a$#|l'95####Z5#DH'######U66cy8OG####?i#OI*ZP####x##SZ&pb####K,#"},
{648.87,843.87,3.26,1.40,"2i>mY#O$'se+T#$v?&LRN@m)l]S3m'3@+Iu#a[S######Sv'h&/T$'=v&>[)g#%LD'~dQ<#$h~SPv$zm+^#$>^S95####%##;5#A5#G9%Yp9,u#o5#U}2/<C([SE5#ee$,b>C]SD>####%H#######wn%S%/######g%#'+I5?'###f%#i[Sm./###.##m02"},
{648.87,843.87,3.26,3.59,"###]7#Og9###'##Ta-[J2###Y`39W:95####h25Ll%#########FF.7m)###J?&<7=8$(=H&L:WB:50,#@H$@.K95####$##Q>#p<,Py6lv&Wu%u%&8m&HzT]$SmY#Hu#`Y>/%Q######8##0_0B?'U$'}I'O7,)m$;-'0W0zd*P7&XI+Td*27S~P####=##"},
{706.00,846.87,3.59,2.81,"O-)i5$Y,%dG#}cT$##^P#hB't$,&##}~'qC7###%##DSH*6'NJ1wY#ZP#+##1JZ###$##B##H^96##J'6]G####N##4JZ###:~.:5####-##?KZ######pv#c:<%##0M23D5###T##*LZmP$6e.######A##$JZ######L;$.D>?5#v<6#6<###*##ZNZ3l$"},
{710.27,880.04,3.53,3.75,"Lm)v*0jcM4z,F?'[S$6,KNQ%_e+S~0K5#66%},$p{C(##kY#'01EMRJm*5##E3Cp7)y>%uK-.Q%Vc&h@#`,M5##.@+H-#bHRK|E+I%|y:{$#]IR###%##sv$[XF###V##I-CL5$###h$#cKRqg;N##cc'q&$SM=###<##?&+DLR(c$&##Um$Tx(du&`,#qR)"},
{723.32,891.56,4.12,1.60,"301###p:2$##9pJRQ'c@&Y5$?6$`Q'2G4nl'95#$##:e&m7/W,$iG#[a6eY#*P3*~.?9)Mu$(yRNR+fU-_U1+u#D5#nb9Lm)###H,#}y)u]7ef4;,#M;)IT2<zRY6(u90OT+mR)>%)[zREZ%###gw$3M9Yl&j5%Z.&0A,q#&/</X069m'R#%Am#2,D{]6###"},
{229.83,898.18,3.81,1.65,"pI)ES/.,#2##-f/(Q$?I(900#l#i#%Vn,4%);5#Hu$F>#{N@`A*ke.JI*z5&ZY6b7Sf$*46%.7SG'2)n,uz'TC:F>#A5#J;Swb#0m&peQac&ZjALT4/I)b6$J;SZl&###N6#GyJ=3BF>#NH$OG#^5#JPAn>%dn.QG#8.'4V/S`6?m($I(-h.T8&)7SCc%3u#"},
{717.39,902.39,3.64,1.47,"Ea/9Q&m[(/,#eF72w+l8.~m'fY#C5#.349m(qY#~G#>g5_5%Ox0OG#t^(>%*Z8XF>#1o)]F0Aw-C5#%=X6&-###O5#faF###$@*76$lB/oc&%=X_81*7*Cv$,L-Vb<|8XfP####?6#McO###;#$ZQ%<H%&I'NS%`]44%*)c$3R#hIPbv*95#_l#el%G].###"},
{717.39,902.39,3.64,4.37,"q::###J>#/-#&A/###k[$uO@O[*###uo&L2>a%,Lv(oY#yb#eWE######-%#f7Y.c#n8-d}3k[,..$B=Y(g2N%,Tl$_I,|?&$3C######hl#|;YM^.m7/c>#>8+6t.o7Y4,#(s3HR'bQ(I5#e8-~#&###V>#pF@)R'|l&Yl$iH'9%&w9LDH&b;:VG#Ie*F>#"},
{495.95,907.08,3.89,3.24,"H>#%H$.,####,c$9-$=6(xP#WZ'aI),H%I~G;##A@,'y#J6QCc%tb#nP$3##?`B8H#Q~01?#Px]E5#xZ%m&B>m)A,#2S<nuKWQ'###Du$2##@jD+##W//+c#T#^bn.CZ%au#CB+-j5n+Hyb#s5$.,#/,####5$&F>#UH&;5#J8?NQ'4c$/,#bN2l@,]u%###"},
{208.19,910.33,3.74,0.18,"g-*bP#En*y6'ne0?l#*^*=_2/,#e>#+;+~XBbB6:5#N##G]Uvl'.,#/d$tJ/1D7{p7>w(xn)[WC`R*ix+vR>sQ)###aH#PsW=5#ZP#|Y#Av'%~,fZ't7)h#&0sW|Y$zP$t7$5qUu,&}Z&9;/U>#87*N5$N>#j,$W[*M[*kY#^xDoY#1-'}Y$/7<q&4ZP####"},
{198.74,921.96,4.21,3.11,"]5$i5$>|_+?&zR*6&/AR*&Z#V5$vg3&l#6,#Z$)O$)###K##M>#Cu#3$`Q$*'f/cQ'fGE8I'KZ&Fv'Wm$6L4Um*rb#$##R@%I>#4BA'{_|Y$~K1dv%1oXbR+Cv&3%+}P#yZ(|Z(S5$/,#Nu#6M0@#`ev*qQ'J0.p1/lM8%17YS)Zo0pu&H>#=J/95####q##"},
{198.74,921.96,4.21,5.04,"A##c%/OG####(Q#av)kY#/,#2,##n*~l%95#&##OU/9w,ZP#4c#6%,[n(:#$p17o6*e,$o5$B-(v?%0S+9u#|7.wc$sL2mM<E>#4H%j&)Ll%]@,]J/D.*p%()m'9T-{5%h}Z@(6pT,Y<0W8VU,$HA.HQ&Up/^[,f5$L$))T?#S.wb#f-$oQ`###w>#DaS>nT"},
{198.74,921.96,4.21,6.03,"###<v#7?'.,#.,#>##Z8-oH)H01II+>e(_r3d6'lR*S%GYtD1,#WH%Yc&gY#bl#~-)UH&f$*BE@OI*(U4|d%l]`H-(wc&qQ`'##gn&rc'FZ&h7&#f.&H%B.,Q^`x$&j.-x#%Ua`m,&<5#%^(###.##~u$~n.gY#/##k,%mC90XEHl#?%(wd)pb`{$,&##E?#"},
{459.27,925.27,3.56,2.13,"0$&J?#Yl^G[+wZJdG#+A.+##uv+###)##rl&######:##OQ'f?(]{#[f^+##'h^Qv#D&3F##KU8Rc&###sP#;,#@u$cG#F,$_S0*O(ENDt6#hJ[|#&Zu%}o#mo,cx45,#bP#?m#-d)######<u4f12OG#Q##1-9[P#U5#Am$gQ%;m'^v(c#&oR%S]3^,%.,#"},
{248.71,942.94,4.06,2.91,"L5$###a>$~u#=6(C##v^8nZ$;)A6##;6&V%:j>%###/K%_9~qc'###Jl%N5#a:=###=M<}G#mfa###f,%`8&O/3-##,r9~A0Tv(.,#ub#$##Uy7yY#(8+.Q$$ha1,#I#$kw&QT5P##`FBqA/Jl$H>#D>#&##4^7.m$16'm5#{fadP#N>#9P,F827,#MgH&@F"},
{226.77,960.21,3.48,5.79,"yY$:R#@o3L>#M6TbZ$-H&[&#e6T95####%(#oR)Q8+1T44##I%.5.#p(?:##S7T*c#_5%-'#h@R=Q&95#1%#@H&Rn#8}J.##ww-GH#psF%##g;Tp-+:#$8##vi+VlP.,#'##`H#~wA'::###k#&_>$E)9UG#2F?4l$Mc$sY#?o'-<=95####.c#2?Apb#W#$"},
{508.55,979.85,3.58,4.40,"Yl&######+##%JU######Q?#wKU95####.m#+IU######V##`@.######G##@IU######+`'wtK######Mu6EbL######oI#%x0######(##ANU+u####;H#$yGr5&###:7&*j@######`,#|H*######*##wZJ95####6##PMUpb####=##}HR######6##"},
{372.52,991.89,4.00,4.63,"D2P`/4###o#%/rR^Q(###V#$%8,,u#:5####H>#:5#&c#0u#FDZbG$###T##[CZD>####d##?r?###$##z>$T5${b#hY#]P#>CZ.,####e@#vBZD>####.K#zWE###/##_Z%.,#&##*6$oP$pCZD>####HS%xDZeY####GT(aOH###)##Bc#.,####Zl#/u#"},
{337.56,1008.03,4.00,4.69,"h&3#########<yX#########5zX###.,####A?''##:l$###HT4######Z,#BxX######jo*BzX######xc$8A/+##Jl%###vFK######d6#mwXL5$###gUKtwXA,$###J}4Q81~P#gG$,##nXE######(##~|X281###xc%kc<m[-###(-%.8,ZP####$##"},
{349.20,1012.51,4.03,4.59,"Oc&######)##6h<######&e&W.V######%J(:.V$##D>#-##]./######;##%.V95####(Q6=[T6#$###SWUs>P/,#D>#Tm#k~1######$##,2VSH(###PQ$sxF`/4###+S*EN;D>####U5#;v(#########=0V{k####*##)0VI#%###N##i>L######<,#"},
{378.21,1014.25,3.81,4.52,"ZP##########v@W######$##N{Wn?*###8##?^TCZ&###3##D>##########CvP######-##dxW######{?#:wW######q6#95##########0mS######D##JyW######Sx#GxW######Cx#95##########COH######&$#$yW######-E+*AS######P&)"},
{109.95,1022.14,3.72,4.50,"k@*se1######lV<bo3;5#s#$N2@lY#Ii*,O=#~,###~6:]Z'c6)I-(;H$@Z%y/URc&E5#XQ%:'W*##RXAQd%KI,5##l'WbP#<.-2,#gZ#>94]XFdP#q,$FV4v&W9)1m09@.#4%-+{%z$W4##Kf.&-&wb#wu%{@+kQ&]R(mp52f*YAH4e+^u%37,,M)$m(Qd#"},
{97.31,1028.74,3.92,4.43,"xWE1,#j>%uI$,M:4e-###/T'zK5A?78i?@l#fA'$(.E]1OG#t#L0d)pb#1##fx-eaG###7##&M<ef/L(,WE9#A,`G#qY3?J/dU:OH'j5$r,$&`:>e.+c#h>$9yXD>#-f,}S)nn0(##O|XeH'sx5%##Xc#R%+]S0VG#,v#~g4LxX9@'a[,jA$~n/PB%<wXA##"},
{50.27,1063.10,3.65,1.21,"1l#)K-`e-t>%~P#r5#a:0Sd*/R*A,#i$)&|03l$'%#SdU[-%[G#-A*b@(MS0v[+uM)h049Z%>eU27'e1=z6%Fl%e%#;dUuY#1$#}gUb,%eY#xFHTv;DZ&?##odUW.&+FH6$#mP$m.#2dU$##W##Rh=######F+;b..###(##37DLG@=A10#####PS&.B5###"},
{944.33,79.71,4.68,3.89,"RmT.e(#d(`(${%1Fp'xmTY.%r5&)##e@-WU&95#&##2-'#7'=bE]MRgR,qS$b%TV/E^[,46#BB5|Y#lu&Ke$###0##Uw-nY#^R,Bw'hoTj-'pmTGQ%%w+s##X)@pY#S5$>,####*##CR)<5#]m+5&#JmT9?#*mT,##|u'(&#t/55,#A,$Y##eY#.##1H&jd#"},
{880.60,91.92,4.68,6.08,"t##jeY######c''Y2B######w7*95####)##P6(.,####9##T^,@gYiP#E0/#kYj_@###K8%6#Hfc'###4##7%(o-+###.##.-'^l$kI&1jY>>I$.)rR,R0GIL.fPK}k#?,#4-#.fYD>####(##vB/Aw);z;I5#9v$1{7a7.1Q#1x-;A,Ov)@u#-D5`v*2,#"},
{29.21,119.88,4.54,0.91,"7m)######2$#XYN######b(#{<H;/$Jm*m1#M5$A:#o.07##piD######9$#;wVC5#IH'_q#Af3To#tvV7S$5,#JI#;n.###jN?######)##?|V3I(+6'E##|/,1_+^vV8,#>5#Y5#a$+<5#>d%#########wa/5d).,####Yp##;:VZ'###&6#Xu%B,$TG#"},
{715.14,125.00,4.25,4.12,"nx(.,#f##8[*whO###)##fY#spY######3##O:;######3##L,$###o$#h@.g_W###|##RJ/ZqY###%##)A&Co3###E##po.######5##|'/P1<###J##xIH7oYNc&(##zr(MB38I+gm$Ux.[P#M5#}k#g]-u^;P,$SG#`&'%dI~3E#c#*S&Wr)8tJCu#o>%"},
{242.60,141.31,4.27,4.15,"-I#gSHSo2yY$DY0,95&##y>$XI*######AZ#bG$######c,$<x*<b/Y_?'##S4ZC[)eY#C##'@Q######I,#G,$######i5$x84QS#jg:fS#t.ZQ>#OG#')#e+K######O$#]-(.,####H##Ce/<##Rw(sJB}_A######|_$YZ'###&##~Q#E07###$##88#"},
{967.68,157.52,4.36,4.63,"P^%KdLWy8###bk2HI,F>#$##EZ%/,##-$36'###2##I]1Rc&R{4W.)t/Wt5%W3W+u#5,#9?$HL9>,#s,%Q-'###sH#=;>$l#U6)%##Dm:b>J=.W1,#P##0.=o1>}Q#V&3FQ#(##}y&$C9#########zT#0*D8?'O#$(7%S3;b5%:e$V~0:##$##J%%6lO###"},
{837.92,161.93,4.50,5.57,"V.-eY####$##;3p######b##@UaA,$###E##tc$-R)###$##Nx2J#%###&##>3p.,####]##h^eO5$###D##2H$#I*###&##&g5KZ&###$##F3p95####F##GKakY#N>#>$&'##bl&$Z#@R+_o4~P#K##8?&23p######+-#ZuP<H&$##^Q$###=I(zb#1H&"},
{98.34,170.74,4.70,1.06,"95#o&#73E###Oh>;/#7x2w/#eA3c&(GI,H{'&##X2/uv*eQ'eu&D%#G$U*##W$Ug##XB7~/#(U8Z[#R$UX6####]I'LC8###,-'0##V%U###d(UyG$Ng86##5L2Hv%9%U>5####'##i13OG#H>#$##.0QOG#yt;;#$X_5QR&Tg/(I'yJ.ax,9##*m'dA,L5$"},
{128.83,170.16,4.25,2.17,";5#.^$M`C%##>o-5o(j^;YA%3J*6&.MH&(/(.,#~6#l%/qP$VZ'N%#C-S2##5-Sc##_81f)*{6+@##P.K.p.###.$#sg9k>%#v'$##T/S&##81S0,#a]25Q#}(6>5#01S6Q%###$##'D.]J2:5####U$A###da395#mM2###RE,Qu%1n$du&fH#zu'WK#]98"},
{153.17,181.83,4.11,1.25,"6l#QQ&#c#Wu%ub#Z7&So3iG$o7/zZ$)@*XK'`l&N$&Xq8Ax+[5#j&U.,#.,#6p3C)U:Q&$c#s%U0/,[w.wV'Am)z&$C%U'x)###hN<eG#`4Bup8so1gl#h,Ll)UllF4v(`m&E-&z#9(XFUG#>u$e$$eR*:G<lY#,~'#B/dS1-?$J03VH':Z%Kc$=S+o#&###"},
{100.50,208.94,4.85,3.29,">##^/T######f5#^/T95####&I$I.T###&##A##QtE+u####?%$5%O@5#95#YT+GlF[P#?5#BHIQ06-n+BH$sH%TOEH'5iR)}6(sP$rI%`>$?(9+Q$yb#Jn%2dS1,#X@*%^&+tA^..Ag1O/EvG%###}0(-A.O]5###G##CD+&V=3##>?&U=17_3}>AoNCqf-"},
{686.39,232.64,5.05,3.38,"###;5#eZ#mv+K~0###nm$3@*T7Y###O>#:Z#,9YD>#]P#)#####7,#[U(vh=&&1RH#,91kZ'A8YD5#gG$J,#Y<YqT7######_>$.}-^DXb80^L:|B%k/3|b#M9YWG#.,#J##S:.;x2###$##TG#j?9MQ'yY#XK5tl#7#$F#$KKX#########1(6######'##"},
{847.47,236.77,4.50,0.53,"*v'%p-###{-$[04(?:###9$#H'VBv'###.%#j[*Pm)###)##E6'l,9###+##A%Vxe*###Y(#G>N@5#.u#5V#pl'yY#qb#Q,#;A-Eg-`P#jG$.*V|Y$###+-#KyR%##U5$VZ#Il%$##vb#@H$&6&P>#H6$1@+Y)V:5#4,#kY#z)V######2,#0[)######7c#"},
{588.74,240.61,4.20,1.49,"VL;######,##ko_###J/2$##x.0*##Kp_L,#ZT51,#n7-u,$B;@######+##eo_$##1f2%##Zf4?##bp_?c#5^7*##9.*6R%u:>######E5#To_%##_.-X#$qx5.##1)ZLf0o849##ye,hA.P]5###+##tY#ko_###gG#5-%(h;###QB%Z>L-H&$##9K$v-U"},
{599.36,239.24,4.40,1.48,"xB8######%##iyb-##[A2###BS/{c#cyb?,#rT6yY#Id(;H$$q;######2##Pyb###6K4'##DU90##Jzb_5#={?3,#Lf1TH$;:;######.##@yb$##VJ13,#ng:B##Szbdu#;q;,##v[*xv%1T4###&##~G#1yb%##k-){5$696,##CnG3C6[836##-x*.y2"},
{790.91,241.90,4.36,1.38,"V5$en$v#R;u#J#%AH#Bg7cu$-o2###A5#~l7T,%###%###kY&Z#3E1r]60,#u&3%K-IR+zY#seYb5%###8C&w=K######g*5m>%O,$Mf1_G#a07Sc$Uy2/m&#kYo?*&##cu#qm;d09###1##yG%TG#iQ&G>#&Z$###`m$?w-&p&`$+Y##_v*#)#G|F,##D>#"},
{803.80,254.32,4.86,1.72,"&?$Ef3t-*Q>#l5$|q?H>#qm&Ye*VB7###=a+[r=0Z%###R2&AH'###/B1kl$<3D###1,#G+/+.,######9v8W>$######wL(2m'hQ%Y&3K,$AjXn#&OG#$I#weI}I/###Hd#`-&_~1/,#tu#wP$Hl#ke,fG$J$@R[+}G$z,&_W)kdX###PG#YJ#kdX######"},
{803.80,254.32,4.86,2.46,"gG$HY1D>#$##c~XVT/###~$#-~X.,####K$#R6)$##@5#(c#hv*u.''l#8c#{~XlY#:5#hA#M]XB6(###<$#'R(o:2/,#6,#]d*0,#P-#}H)dUKk~1hY#uG#k<+p~X.,#&##F##UaX######G>#$##0h&xQ)EH#%@)ZJ,,Q%#'%O*DOl$vl&c5$Ih5###56$"},
{716.79,265.62,4.18,4.76,"ll'######k>#Oh>######39#JJ1######po%8%-###2##V])lS2###$##@c#6%W;5#;5#of$DA0Ou#X$)C])}%WI>#qb#~{'6~.###(##Em&$'W2u#RG#}8(by4.@&BI+G-%n)WI,$pb#>-#pb####V##K*WGv)J##e*;K*Ww%,`](/<@^I)-,<&v'###Q,#"},
{581.60,275.54,4.14,4.60,"@L:###0,#5,#7g`0c$###n>#cz8k{@]>#k6)-$##'6.-#an0,;?###5,#;,#(g`Lc$yH*,##VC3<@HaPO'##W(&*/[;c%(##d(>###%##cG#pf`@##GS04##Kx1Ov#gf`5##So/tw*M@-7##=:9######$##cf`<,#l%0*##i@-~v#ef`&##gn.OH$4e-.##"},
{218.88,283.71,4.58,1.16,"??'95####0##~ZS######('#]ZSSI#ziEH/#8q;6S$HkL+$#+o1######>##q[S0?%D>#~&#jsF^M'k<G%I#)sEV.#=:;7K#?S-######%##=^KwB7D>#2##l`.NuBb./8I$:q6,9/8.-u&$t#%95#######Uv#4[*######y'%u]7###&##q'(0PJOG#(##"},
{737.91,284.81,4.65,5.99,"{[$%f21,#/,#DT.tA3V#$Gl$.//m5%|J&B3AB[(3c$Nd$d/3*Q#fH)###&##,fT^Q(####H#qfTvv'e],#/L5H%m%&IhT`:;OG####B##IU8$aE######ZN+ueTz7*p@+aV&N#$>/'5$G`g9(##ZP#,$#fdTg,&###%##DhO-H&###FQ#4_*######zd$ko5"},
{215.54,324.19,4.30,0.22,"u`A######.##YBSe-*,##eT1j5$p'41##/uH###4o-###i<=1(8######(##e0Uql'###Xy%18/*7*###DgJ###I,$###cIJQf2######%##`1U######?6#.~Ox[,###f[&&##c;<###6c$_H'E>#######G2U~#&######_j6F/U95####uQ&-0UfY####"},
{439.86,486.30,4.26,4.43,"WvR######R.&ovReY####k'+UwR0,#}k#A,#gG$4##hu&###yvR###h##-`5bxR~#&/##PL/1xRQG#iY#eu$dH(J5#eQ(###4QM3l$'-#o=3d(>A,#0q.}yMuuR/##uv'ge)yQ)kd&vc((##&{RbQ(###+-#,GC^^/`p7~H$i:=p>#W[+p,#e5%4?$>067,#"},
{462.90,491.26,4.09,4.15,"s4BxG%,##%6%qo1E>#/.$5W1:S/###1R$@0HaL;.,#}#$-f+QdC.,####.##M'VtP$%##J_$3z9c~&YU4o`0}I/O##s.,KS&dXC/,#4,#UH%.*VY6)###Y?#f+=kh5_5%E$#f%0_5#Qu%L$#H?'###P##6'Vv_<.,####vc9Og7ZP####mg#Yl&.##-6'_##"},
{809.50,496.49,4.82,5.65,"MZ$a$+###VG#cN<wu'###)##_kHOG#lH$eR(eY#&##GV1g5%T,#vjHZG#eY#bj6dIV1l#oY#MjZ=7*Z.-n#$P$*zJ$4cLmY####uu&bl#qx2_FI0I)k~&PC4heZ$7$9M>F8$gZ'nU'+QQ%#####&##[Z&'q2N4K:##q_;,~)7eZ*$#+>M;7#.c#_[$J>O###"},
{324.90,506.30,4.72,0.26,"###B#$rY#95#tb#f@(ZG#G#$RG#kf.N,$u6(0,#H/0GQ$%B2###(/%,Q%rb#@7,vM-Yu$WQ&,8T8B04,#t~Dw5%cE4~4>{9T]S/o:)%Q%v,$Lg5u]+412|^5QyJ%9T'/.:w&Im%g;Tz~0A$'60,$;45?'(##-],0CMM-)`G#cn+Qq3)m(@I$d>$b8)+%+Ll#"},
{490.19,505.99,4.74,6.19,",##mb8%h8CZ&B;%%SQ-H&###XX.ac'######HZ$fG$######$##pSC&~.###+b4fq2#?&J5#3sV{k####H##e]2SG####$##PG#j~%-@+<Q%Q:;*Q$[>$.h%$nV######l1#*05}k####@#####1##7c#B}B`$+######u}2Kg9######'D#FZ&95####nG#"},
{207.53,519.79,4.72,4.42,"h}K.,#su#GK/gGKL#$_&'ax0q&V5,#*Q$D#$pc'~G#3l$###q);{G$zL9PI'A#DNK-j.)|%,$'V%?$6Z$9l#|[,M,$A,$###$RPPd)^..oz)VWCrc#9-;dRBk$V`G#+S)Vw%G.-#J,{k#'##P(Vi7..,#`,#KYDb<4H06yc#jcSn>#J6(<7#=c%kR)-J*4H&"},
{404.32,520.11,4.58,2.06,"$##qu#U;,*6':@)w>%($$]%)qA03l$N,#Qv<#kG###g5#'h3###P##&7??d*^_>%##,y+Kh5GTS?v'$7'X#8]7-*$&~~*U91###(##s_'QQS.$('##eS%uRSc:7X~(jl&s?G;H%,7L6#$#Z#/,#rY#<Q#Uv*###~$#n/1d]6###kh&yMA)m(9,#<KHSc&8l$"},
{404.32,520.11,4.58,3.49,"Xh5*OD95#O##'3:{d,@H%bf+Il#nH&qNBM#%}k#x,#]oTE>#4u9}@+xY$9##P[=enT]Z&[m&%$GFr<bA,Vc%dc'###`>65~.tm(.-&O#$Ie&Zq6qf6.##Rh0)nTnl'.##Y^%%~.###u/#]f1xY$###,H#;,6n?*###v##~7An`E###%##9D%Yl&/,#Ou#m,$"},
{220.38,526.13,4.22,4.74,"K:2rb#%w'DI*|^-'x,y-)ql&lnGM?%5$'--'+j:]P#.,####Se+IQ&oL5,H$e$Q{>%:%(4u7&tEG5#<|1tODC%QYG#xP$oY#6K*i<DTl%5,#S(Qke.YH'_-%`W=@:-l$QXc%E$QH5#5H&2,#V.),g5######oH=<L9.,#$##YGBpFD~#&V##b#QUG#95#m##"},
{768.62,526.57,4.31,5.27,"E$&95####&##ku7Km*###'##z->.@U###Ox&349?@,v5$/J'3H$D>####.l#/cID>####r@('VVg,&###/;)Kp.|e2###lS(yx0######0h)`[,######}AE,a?mP$3,#)J@/n%Xe0o5$oQ'},'###+##jWVr#'###-##W6:AQ%E>#=5#|T11,#RG#_c$_~0"},
{768.62,526.57,4.31,6.18,"X>#~|WuA3{k#xSMc5C*6'H##g179_3###A##$##wg.F>#}k#@6#4xWOG####J?8:QPD>#/##vyJbZN###x,#Zm*hQ&[G#pm(G5#al&###$##66L1v(###w##]xW>u$###9%#Sx2a#&0,#K?#D5#o09######kPGi7/###M##TxW.,####G7#%o1.,####.8$"},
{197.08,533.31,4.67,4.28,"<R(uEAY5#:e,0L6;$&Q,#YM4}(S###c##x@-PdM&#####$##E-)$##eI&?O<?$S`P#xI&JM1[$S5l#7L)Bf.N$S*##qb#Q5#Ie,)6$%:7o[&s$Lcv'=B2g8*3%S>],f.*7d&?$SM>#eY#A##n_4Rw-bQ(pc#_ND>H$C9.i58VsG>Z#6W2-1/w=LP5$###A-#"},
{446.99,534.26,4.07,0.80,"+p(.{8Sv)###*Z#AS*>nSS$*@H'###$o-)=;#ZP######a1$H:-XaET,%(##6.R+S-d-'G_52e.###'B%(nNxlS.,#5##%C&-)>9Z%i,%,7&inS.,#;,#2$$$;;###z@%(z5L*Fub#U5#w_/[[,###@y&A17UmS###)##_%%Y}H######Ow%Vz=&##M>#'X-"},
{640.17,548.31,4.16,2.52,"|.X^B&d093%#Dw-RM$).X[##[f5T,#Gl%@~$[P#$Q#Gc%$m&4.XJT#).X/%#1_;to#).X#$#Hg7v.-95#0$#]P#QT*J#%'##O(=*g#+.X8-#6D/M^-m|H'##%S',FD###%##6v'-A,###a##=$)xR#j.X:?$u@&n[(zL;Dw+b~)c6*###8['2&,D>####q,$"},
{848.37,551.30,4.67,2.87,"xG%$##}k#H5#T,%7##XA,UW;j>%)##1_'i{>(n-###[>#Sw$oA3###ub#`##h8ZIl$2x(-/@[-)4d$_=ZMF@S6)+##=S,~Q$p&2###Q,$+##_=Zxp44I*>$$8]*Ll4Z8ZK#$IR*}G#4S/4,#u>%$##b5$.,#&SN}P$,6'###J&/%m$m7Z###JA06,#AR+6,#"},
{848.37,551.30,4.67,6.05,"k@.)##&x/4,#&8Z###gx/($%Il%###$.JyP$N#$95#(Z$$##./13,#PR*)Q#b8ZJ#$E]*3H5o-*8$$_=Z@:5F#$*##[o1###s7,-6$jH)*##_=Z|*@#I)-d$J&)K8Ak8ZHl$jY#i,#Lo3###K5#v[$%~.###h^&BE?-H&&##Ao,zW=T,%5##D,$G5#)Q%%##"},
{422.36,558.01,4.16,1.89,"3~)]m+###V?%1w@rQ)###'##vY7}oQ1N?SG#Q5#**26nQZ6)0:&_V@######|7JH[+2##gc%6L2AEAT-$yJ.}~/ko0|g.)z4[g0r5&.##A26)nQ95#.##8o).{?qb#9e#Eg3B-&TC51_++[)?x1P,$%/'npQ8rA(##Q[%W35LD@###O##jm&Y,%%##qx#6o2"},
{666.89,561.00,4.69,2.62,"1o/OG####FH$4&Y;l$0,#yV%gsDvy(o'Y^s1D>##$#;6N2-%{n1######mA$f&YG$'(c$,/#x`>)k.9-T0##r5&e$#~J2FZ#WI-######`9'p%Y>6$,d).I#+8/;9$juR0##/7,:##g,&0H#$m(###%##P~)WFJR,#;n._H#Gm)P?#|{D,##a6)Yl%OG#0##"},
{428.47,568.84,4.28,0.40,"eY#3##@B/DE?r@/###6,#|`6Lm*###-##Vh-Y>$vb#&%'@QA###A##~R(fuP:7,y,#RJ.k_;GwPod+kY#M}3wQ'6,7Gs>TKJ6-'W&%SA1:H&@K2N@&rx/Mg4o7CAwP|w.Cv&Je)GzP0e,^>#x['TL*:/2$##j],]ZA3$(0l#Xx.>M0*.,F,#de-qn+:$(Y>#"},
{832.16,577.52,4.40,3.40,"Jn(sZ'qx5###IRRy5&###>##'D5###/##$Z#-f$X>$Zo'95#6.)f5%`3DsP$1BZ{k#%##6l#qU;###e,#Bc$Z7$B,$Ql#|k#vG%###Fz2#9X.AZ###Y$#z.TP{@###R##Qc%(l####a##[Z'0AZ###6d#kFZ*]26,#VR;ZCZx@+,u#]|5Zc&R#$.,#/##u>%"},
{339.02,588.83,4.60,4.40,"F6'YQ%w]5Z#%n$)vl$&lJ*f)1K./%+jn$G^3zK0TJ()A.P,$17(8<,TK56?&~[QYK)K~.7'(x&4E~&QmE/:1Ol%CZ#vT2G?&~A&0^Q:Q&###w&ILkG+u#G##WW7&F='A.<e%=Q&$v#e_8kc%z-%g/4Ed#][Q|o._@.###UwB`y6(c$%.(H*.rb#z#$tmMml%"},
{754.93,588.93,3.66,0.09,":5#s,&###Al$Jq=8#$###KR#cUg######/*#IQP5['D>#?%#.,#$##RG#87+C&^######9d$}VgZP####]$#.<?fS);Q&R5#.,#3,#|b#u5&le~/,#0,#f,$*Wg95####f##?h:Bc$/u#j['.,#wY#0,#H,$-lO;5####1?#bUg######}%#+06######?K("},
{204.50,602.73,4.51,0.32,"W7.$##%A*3F:E96DQ&h>$5:(s,&)[(*R%fx,_$*hZ'sc%NC/OG#[##?x,&7Q:%-v6%cV:#T.67QI7+=l$ZG7EZ%CC-]tBK8L&%*/('ZJ0&m'<//m$&Ob;y[++'J1M;8J,Z')Vm&E-@/]/w,H(p*E:.fH)%##Vg4626%.+9~$d~0r>$[?(f44>H%Fm)Z,$yt="},
{790.28,606.29,4.17,2.43,"G[+R##TvSHx&?*Fj##NvSS[#$vSJ>#.,#$%#?#$'H%###%##(K40&#&vSP?#EvShd#%vS;@#XaF+T1###'$#'##Y1:######T%.W($}uS5##i%L0K**L9,##>&+}GN###%##F##^aD######}G%|@#>wS--'jf/ac%9@*rA/Z?$J]1=l$=Z%;##ll9ZP####"},
{818.05,616.69,4.06,5.42,"###{$#Lp8Bc%|XKCA#A;@R.#s>Q6R#+%-g'#5l$E%$vA42##{Y$|&#p>Q%##e?Qr7%TB7i%#T?Qv30[?)K$#0,#xt:ac'###z>%.%#L?Q#l#VCQiw,BA0=R&Px+..EUG#w-)###aD2cQ&xu'.,#&##IZA;d)+d'v>$AR&[a9[5$jx'W#%tx/###-K(zP$bT1"},
{692.63,625.79,4.55,1.91,"{>$7l$/,#.,#L1TZP#######zrWYG#PG#QG#Q$(`['SQ'Q5$NZ$?$)###Ll$GoWeY####Cd%=qWQG#]P#},#h/2r#%Co0###vb#ZP#)##)gU[&3OG#'##N:L8pW2p7###e]%*n)lRR@Q&vY####)##vG5anW|b#qb#%S&H:7,[%?%..##JB.)Q$ae0$##ww)"},
{418.35,643.91,4.27,0.41,"U[*8##:;4}2?JuHzP$J>#F:1N8.*93h,%wu${H'vn.Z$(Jw(###f##n'/X6R^H(Zy%M3A|OE'dO50/4Z%JL&u$*D8)vV:ni3B5#Nz$,ZM+%-od)|q'-i>nd+9xMZJ/Uc%@e(vK1y14,A+aT.vu$iq&x5R###mp,Cq4c5%&l#9?M88.F,$1R#RU7sT2-c#@S'"},
{418.35,643.91,4.27,4.08,"pb#A##;G>sv(6-(lY#En%Q27)IS]P#)6#/K+_,LJ8-95#*##zG%Lc$;KSz/1&~,R7*^^.K]18JSre-Jc$nZ$fg5{qR###'##L;<l8/Tx/<(2c?)Dm#TMS)(5i<C8]*8y1T[%UZ$]G6Y92xY$Md=%f.]#&S,#LH%Fo(2M=@e*Qu$f>$So*.x+:5#I5#KY7=R+"},
{310.66,660.20,4.21,0.45,"bU/7n-[%-^Q&IO=;g5aP#fm%:Q%uL2ce,PK+($&'?%b.)Rs;<6'+##7f,p-QS93#/-Y&1V//e.QL'4rb#oN.Kv'gb693={~Fk5%B%$[_;4A0YS/nZ$@;8*i<*.CO.Q4y3D.*)%(_1Q:Q&LR'Y[&E()+h<%##M^-^39Z[+Ml$B[*ZK*)g6,&)O6(vc':#$9/("},
{184.63,670.30,4.80,1.96,"$##Cc#e{4eA3/e%ve0Lh3EZ&>/+}g6D91u-*^R%0v'pJ*|Q)vY#pz&bIRRu%6N=vw(qy-Cf1jS/8I*7V%>19^,$X25bL/fc'.[(SG6qM?~^5?IRY?%8m%C)6]q>TG#H($vNA###XQ%);&6HRcJRLd)(@%sW9'z9%Q$=o'Lk='::###7##xI'xG%###.$#P:;"},
{184.63,670.30,4.80,4.07,"rH)|Y#'K+jo03/S###??#[@)(uG/7*###&##+##}82######n078m&:S*_n)]-SNm'Y['qI&-s@KTL###1##}z7'>F3c$nl$j@-4[&i):|//&.S{J+=d(~v&'~*'v:$J*&{>{wFIO>#R&NT39H&n.&-z8RR&Oy6A5#-?%Jr./m(7,#Lv$t0SFu#IH&n6=zbL"},
{219.64,683.32,4.93,4.44,"_P#BR'8o/uc'H02{G$#q4We(U:4qm&1B,xv(6h4:n&|m+9@)@$']T)`J.5n,vHP5J(-K.{1,;B44[$oKP1g-2PFiY#dl%j&(.R&,%A;m)###_9LgIPQu%B##_i;M5B?96~/&xK2Ll%G#$^25>5#6L/Z'+GjEa/-qg8O>#(gK%z8ml':Z$1|-HQ%9e+Yg5.~("},
{190.33,689.96,4.39,4.24,"j?$6)9tp2c7,P$'t8)%3:me/k:Qsb#HQ#CJ.0:Qk5%###$###H%PI';7Q}c&>D>O>#Wq5>N6^6QbG#D0(P{:E8Q8J*$##N5#Q5$Mm&eW@nH([FC}Q&uT3HK,MlI|8+hA.Wc%icJ`).Y?(Q5$an/Ff+1S,ZI'`&4Jc#:W8j59iT5Ec#(b>y[(FB4RQ$(m&~z2"},
{702.72,787.38,4.42,4.57,"(##=U,qqXh,&BVUJ{;wI,9I(mpX]P#BH&(Z#@Q&C,#8D<VG#Q>#?^'msX*y66,M6,#:**>RLwpXeP#+Q%Km%IR).c#:D<###3m%'['M.(N6'gu&<5#MU(CmORXF:5#:,#3F9g-)cP#Cd(F5#nl'######]-$ll'###'##;/-6#$###&##x7,OG#(##^u%'l#"},
{708.60,832.18,4.55,3.00,"I#$`>$.,#0,#nQNJ>#.,#?,#X3B###Dn+j,%OG####|~TgY#SG#cP#B,$###i?N/,####-##dwU###(@%+V.mP$###}{UBS-/,#$##&H%###8mQ###.,#0##zzUX7-w#%8/':7)Nq0pzUa?'######c#$7#$H~,B#$yb#>l#|k4^wUuG%(##g?$}zU)n-###"},
{717.67,850.40,4.31,2.87,"R#%.,#%##KbEW,%=5#xc$b4DG[+###c,#qH8E>####o`+`18Dm'Z>$$##n5$kq=@?&{P$s#$YeZ###2c#A2,9m)###1KKh&2H>#8u#pP$aP#`NCnG$###$##(gZ###yb#AH$D/2###JgZ7R($##3,#gu&###9C9:5#95#:##seZ###*##r2(.y6&##IzK+=9"},
{448.79,862.06,4.14,6.04,"###+J#yV>PG#1(:|l$P.+CI%2n,7x*Qr6Fo0du$fh:*e'-kG###<##xT-K?(K;=###kx'^I);AY###dI%/VSg,&###@x&_(^###Jx.:?%D6(mg8i/.J'2{b#>)^1u#ub#TS&eV;Yc$W<@KC5L>#-TMxY$###L`;`5KeY#.##K)^6l$###8##Ks@PH%:'6$##"},
{740.67,877.84,4.62,1.21,"fu&0##;I(S7+X>$Yn#;*EW>#rc'sw)#@)'(-[e(C@+uv(v$*~P#%##Jf0pH(mw0=##u-RyG#y,R)##eL3_),t>%`,%)LL~-(M,$;6%,B/y?*]~/o>#T.R&##|1R]#%wK47Z#7.)Vl$_'Qe>$###3H#K}35J06u#X5#2wA)Z$u.RkY#Ko,6W*v7.xY#<dH0]*"},
{749.10,886.32,4.35,1.38,"F>#N/$rx5%##H?&rg/)-'rI'~y)h/2@-&ev(r>#]U/s.//,#DZ&Q##IIRaG#@HRE,#u//f<.|Q)]c&id?}x/m>$KR*,A)2Q%$-'(##;KRRG#9LR[P#T91QQ$fW;al%J'K[c%fY#$##|N.=e.:5#l>#)P@mP$e#;B,$;{5^,$9nD%d'r^2oS).,#+##V4;KH'"},
{749.10,886.32,4.35,3.91,"$f+:w'L5$%##Q~UE>#mY#r$$v06_)0g[S)I%<E@RL,+T3K##_o2$Q$FH$xc&8~U]?%Pl%l$$[^9AD%m[UZQ$n+Mn##p09FS#<@+>#$wU0[]1~JQ9ZDw,&Uc#rA2754px5l$#2QQf>#T,%AB#Sl$xZ%j<@kQ(W5#f/1FI'8V;Sd+3['$##hN.jT5M#%###/7#"},
{689.74,890.02,4.58,2.39,")u5nB8k>#52:}'&+iSQ27<d)+Z$rN-fD;xS.CR*B?'o>#$2-J#$m>$G>BI=EI]2_.,yA/G^3#?&'Z#VI%keSn?)###8##k(/%l#v##+iSqZ(edS>##OV=3v$Rg8###66#;K3%l#%##7R$/~.(c$B(#6dSQ,#X{BT##t.0(/#q./7#$<5#L##*##-u#9c#{k#"},
{689.74,890.02,4.58,4.04,"$r,jI-*Q%?##k'..@*###E,#z@/RG#0,#X-$eY#,##eG$bG#5p/=H&^:*g19&9Ur#'M5#Jv${833V;###i,#T##>K47#$/,#..,l%$)<Ue+@8z5h&2+&.[J-il$,8U@,#Si?C.#fEFb##qg:5U3e}2s|G&Q#v=E},%BQ&sW:yd,UG#A##)<UP>#pb#?'#r7U"},
{740.42,890.55,4.53,1.02,"D>#+$#5FGP5$I..f$$tA3rC+(['Nu$kz5:)7E%'401&d(E>#|,'cQ#;dUE5#hdUO5#CL8ag%-&1)c#ChUlI'`,$du%VL4Pu$0?&$?#xeUVG#)hU8Z$xf4cw%|x00R%SfUG?%Qu%)##xM3<$'[l&-[&%<:pZ%/cN'Z#/g1G40]Z'f?#HfUxJ+wv+q?#L'6;,#"},
{755.14,897.13,4.81,3.93,"VnSRq/_Q(<-#B'0,p(1S/f5#<nSl5%'Z${d#$r:?H9|'9jc#5V:]7+<5#_,$;r:Om'5Q$Ac$5nS6-$Fn-1n$/f2-K$xmSoZ$HV>95#wP#oq3|x4Q5$o9,XC57oS#c=6d)~c#q]6m)*4q<b$#ooSVu$Nc%`R'm?'Yv&/cF;%,PZ$tp5^$'MC7Ee/e?&95#4:("},
{449.63,910.77,4.51,1.52,"Yu#4zZ###oY#<v#n|Z95#?5#Q-)|9HJm*O$$O~0l#$`5%GI#.?$SzZ(Q%###s6K{8Z###t6%.K3<R(2,#9D,Pn,###$##S6#dd)Fx0T_4]?'AyZ/H&'##p.%|2B*H$Rl%)6$,o+=l$PG#$##bG$###Y^&lkHch?%##O>#|G5T@-`u#*n-X##co2M,$U,%E##"},
{646.68,913.81,4.66,1.35,"I,$Z>#=$)&##1Z#vw)&'4/l#3n'A-(4v&gc%E>####wCbqb#wl$S>#yn1###^w)dU-}Bb@5#msaE.+Mw-zG#{u''##tBb8,#B@(<5#[m*###u6*E5#SF]IS-hBb;,#9T.&M/5v(M##QBb-##pu$###S>#rb#Dc#pb#lv$SI,v#'Q5$c$)ST0ZP#B,#JBb>5#"},
{646.68,913.81,4.66,4.37,"S/_%##OG#k%#B80f&/EH'<l#=7%@n-f,$pb#dG#E>#16%###X/_3##XZ'/&#>K0lV-f/_V5#94_7x-Ld*U5#W[+$##nd)4,#_/_=,#dQ(o$#fI-:?#83_g%,p/_D5#N%*HV,Fw.###s6'Q5#^0_|Y$hY#R##0-&Q6&4n'Jv(005`P#V5#nS)DZ&'##Cl$K5#"},
{467.24,916.89,4.88,5.68,"6l#E>#'D)$.,/v(%##V:'+kCR@Y###V>#s{'nq?L,#v@.CK&rG$0u#<-%8.,ZT6[>#G+@67)<AYRG#Pc%1E'GB5>w$>cGt36$##F>#+H$/R*-A/%##XW4Fc%NFYLU9;v$9H$~(0y{StA3%?$$##m5%E>####0,#F#$2@'D>#Nm#jd,,@$%[)k5#v94O>#k?("},
{686.34,937.85,4.68,4.55,"Z7.#########_K_#########dL_###.,####&A.@5#Y%*V,%n[,#########6L_######?##JO_T,%###6##[9.7B4:#$.,##w,######9##6K_######$'#ZM_yl'd>$d.#_A(dbK/R'{5&E82######e##lvT.,#pG#:'&5US/,#Nu#&$$=D<2u#AI%y+B"},
{214.18,953.15,4.89,6.02,"M>#&?$(?%>J/E.,q,%2-%ZS/iE]/Q%*c#Q6$4JIeY8$bI~l$###$##%Z$FS,&]2ZQ#2U7e$)wB]Il$;l$a@#3N;Sy.U#Q=#####%##*m(oc$8(9,?#buQM>#fF]4n-}>&L##PT-v<<)sE)##I>#&-'TQ'V>#.v%'?%^YF95#2u2pe1|H'hY#t$$tvH0Z%8,#"},
{240.57,957.78,4.23,3.36,"|k#<-#o.0%##eR-.H#{u&]~)wY#yG%QJ$zWCG8(~#&p[#g%03-(R##s96~u#1dT.##Gc%iX0h6*###iw$-hTz>%####7=k09,m(###1S'Xu$RhT###Ru#y$%PdG###0w'#o-QG####HiTvG%V,$###xu$F>#NP>###06$Fu$HiT###iP#Qd&S..###Jk2a<?"},
{250.47,966.55,4.54,3.76,"###Y##p[-$##ku&*n#8C:S#$|Q)+.)[?$e*@XQ#P..k-#LK[###D##:A0###)^6T>#PdRN5#(K[1,#a,%,H4*w+###m&#9P[F>#`>#J$)/,#Sw.$##t;4GR*9P[{k##v#EJ(BZ47m)o##=n(x5%.Q$OG#$##5d(.,#j##M6(i&&6#$l$#R[+=j)n[-R>#Pl%"},
{521.19,970.85,3.76,4.16,"pJ3######Z##J~W######)D&E_;######>Y1i&5######RS#.(;######[###`W######X##=_W######Q?#C07######n##ay9######O$#p[W######/'#n[W######g'#US1.,####s5#4q<######X$#M[W######i(#N[W######5'#fu&######vP#"},
{465.50,998.32,4.92,4.52,"############o3F######:##'qYL5$###a%#nnYeY####H$#############ecQ######C##_rY;c%###8A#jpYbG$###'-#.,##########G?R######Z##@pYpb####,)+}nYD>####6n&############<<A######$##3sY######Nu#coY######yP#"},
{405.21,1001.48,4.59,4.55,"b$+######&##]8Y######9$#?9Y######77#BU9###*##;c#e%0######3##$8Y######6z$98Y######VE):q<###*##LQ#c%/######$##'=Y######b,#B=Y######~Z#;V;###$##I5#yc(######%##}8Y######Y$#U:YD>####s?#S|B######$##"},
{112.67,1012.31,4.27,1.65,".,#8%#PU:WG#zy/M~'c$+Nw%Mz-$V61S,[[(;#$?m&v8-b#&Nc&j##,IR_>#HHRA##d~.[3.l$+YZ%3LRNV4%c#Fy318/<5#*v'###AJR2,#EMR:5#_&1YZ#1W7z5%EMR}Y#4,#h,%i39cH(h>$=5#PYD###Yk4+$(dL4###%E*s=G%C5%##El$g?'Z927A-"},
{98.74,1016.23,4.74,0.97,"rQ)j%#O]5:##~./tA#ax40q(:m)_6$NB0BhQi&/38,q%-t$*sl'5T#OsG/##2dQf$$:/2b^#((:%M%ccQ0f'P$$XK,ye1UG#GH&jQ#tcQ/,#LhQan+480@?$iJ+R;-G-QE-(<l#%6$$906-(}k#%##9E:[?'{S/@u#Rn*Y#9k,%$%&ST.eg5ZP#7d#qU7Zu%"},
{98.74,1016.23,4.74,4.46,"|g8U$*###$n#.]07}2$_8MQ%5q,yf,8B/V,%q;6.,####$##Ny0&p6###)##gHP<w+Og*@V4RU75,#0pGiS0ZeW###HZ%xP#}~0Xu%8Q#=n+7fWf,$Ix0Z.&`h>Vd#-gWo5$]dW)##Z@-q,#nI,mY#y5$&&.K5F_r77R)OR'2K3E}/a)C(-#EcQ;,#vG%D$#"},
{55.39,1033.94,4.46,5.92,"+QQ######Z&#6QQ<m'###l&#i5%wy6;5#5,####.;14Z%.,#/QQ,##95#,(#=+I,C-0Z%n$#5`>Bo/,-$T`5K,$y@'*k2D2>(SQAc$W>$+%#mK3a~(2S/G,#GKI8tI8l#Cn(=o%qQQP;/l%/f|6>c%)##Y94%~%?9/MQ&WH'aA$pq:Mc%.,#H1'fB7{k####"},
{162.91,1043.03,4.41,4.45,"_B5######%##uVYwb#=#$sb#.a<I>#1H$RI+jB~TG#SZ%af*r&5######(##MD~######5,#=C~###1##^y.EC~p[*b#&uo+Bo3######*##IA~######i##>A~###'##eL(un->$(MH#~U5281######2##0A~######&$#0A~######</$n?*###~Z##S+"},
{42.71,1045.83,4.89,6.02,"_5%######,##}6W######W&#97Wpw)D>#E%#=v(100(H$tI*ed,######@##S7WA5#OG#)'#b.WjA'en06$#~fNw1<v#%nK.An.######Yc#_<Ww5&95#4J)Ss2Bz4em+Ll$5B&J>GAQ&/,#yK6ZP####K8&sU2PA2###%;W,S'u1=M?%iYH*H#'f/DK*)I*"},
{989.65,81.86,5.11,1.48,"###$8#XT6###^>$z;/VR,7%&ie(r&3>Q%@I&eJ'd?M3Z%&##T,%-$#BHQ:,#34IX##zT3HX0XZ'(6%Ok4'C2'##$T0^K(+u#z5&$##GJQ$##.LQE>#m94sc$XD8i,%PpL+?%###,##p+5%Q%/,####qb?95#7+2Zl&7O8###;|-sg6mB5Q#$lY#]@(/2:lG$"},
{989.65,81.86,5.11,3.99,"v_:TG#UG#<5#fYM###D5#EK't-+W?%WTRe/+}QRXl$P^6J##}V?###:-#3%,MRRVH%K#%cd#fo4-M%sQR$m#>QRg##Zg9~%#o6*'-&/24{I-C$HMZI$l#*-%U%-gG7tn10J#tEG>##uG%O'#0,#HQ%$r0_@.O>#Jx1NQ#0,Gt5&+$'###@X/5J0######b/#"},
{974.37,87.43,4.96,0.85,".,#L%#]q>###~./8/#.o2*'&ac'L,#Ah59mB8?%Iv'o06-[)k,&K%#VcR'##1dRXH#PA2[K#My7od#5dRy%&5##*v%DdR###IQ&K##]dR(l#NgR<['^&3DI&d/..<-U)BgP#$##@Z#(?L###^#&M,#lXC20.Gn-nZ#m`>[BJ$m$Bh-I)?nl&:6$SZ&@n*}k#"},
{953.55,90.04,4.96,1.89,"e43CZ&xZ$>#$.f$]p7sz-U,%2,#uQ&/j9WZ'-c#PS,D-'.6'I7#X-*w:%,d)Dr-=$SZ%&wR-@%S}I-RE/(26J?']6'Z+3We.$##.,#F'#O]5;n-6#$~C#~$S^'S/d)#f%BC2MS(9o1~N,uA33##i,&?,#J#%fG#X,%4%'s?*J%&7~.PL3:l$}-#@$PJI)D>#"},
{715.88,108.32,5.82,4.26,"4[9pb#-##4Z%w`Y#######Q#%z8######XZ${?*######A##4s<###I##626Q_Y%Q%'##Py)]U8Ed*[c#aC4_#&95#H5#hl$O/3'##PG#l_,v^YfV@)##,B(*66S7YEl#kZ'<z3W>$E,#Rm&3p7K,#w#'//%^B5R8-*w+Jx(Y%O'B2[u$P@'%]Y######<A#"},
{810.13,144.31,5.58,5.51,"So/{k####%##r9H+u####&##'.eT,%###&##u82+u####E5#{I,uG%###%##rAS(c$###'##{)e######V##z:>$##G>#6-$/[(T,%###$##<CcZP####T##{'e######4%#8_9IZ&$##)Z#{5%bG$$##$##y'e{k####]##l'e.,####X$#I7*q#'###'##"},
{55.35,167.91,5.53,3.74,"h~0cR)^>$pk4*c$Lq'Gr<(Q7ac'k,#9o.5o'OG####Qs-wc&Uj<lhTpH)=0',6JuqQ/6'8m$480oc%|k#mI#'##.H%X%(tY#%_:W].`dT}#$DeTnu&M5$S$#z&2#$'###@Z#I,$#H%N>#Fu#I97,%#[1=k.#{cT###95#D'#8m)###)##DI%Qu%###1##)C)"},
{849.73,194.40,5.59,5.71,"yY#16'7c#fY#JJU&Q%###/,#KL[######P>#Rm*^5$WG#sc'cG#O#%#x$KT5m?ROG#v##iy6?M[L#%###AH$tw-C&..,#dP#eL%Y83n5#9n-&<4=6(s##$dHcK[nP$###h(+dS/ZQ'###6J%O~&^Q(B##:US1a5$m(;##fKLI=GD>#####W$wR-eY####P@#"},
{855.65,221.09,5.29,0.76,"f/'F$TeY#%##`z'vDAD>####rL&.R*######0H#*e-######x{<DA0###|$#Y}9Z?@###=##O)TsS0###<##TR'-J.###$##C@,K''+u#%I#~uLnE-D>#P$#[%T0Z$$##I0#aR,&l#1,#Z5#/,#&;%f6*cP#2xHNx){k#F,#O)T###$##l,#1e)###%##h>#"},
{89.76,271.44,5.83,1.60,"Z:7oQ'RG#_u$X=2gb=J6'bu%|N-0;8sAT~P#t/15,#iT5###5x.?m&?%-&##i;7AX0f@TfP#7DT&S+;PF;?$519(##l<B@5#D].QG#Od(1,#v987l#4{Rtv&@@TxS(;3@Ax%n09~v#e_?9##zd$x7/F>####]H&'r8Vn,Y,$Bp45E:*H%tu#V]53$%L5$p##"},
{748.68,298.12,6.23,0.45,"YA+(%)_#%ny/|8-bQ'`N/0'PzlQ3,#V0*tg.FsD9#$1,#IH#J=/7^4.,#q,$3sV^]3H^0]f/?A0E#$2f*-yL~Z'OG#N##KJ+GS._>$###R0($oV######yh%f_?QG#Q5#{7*###/,#<6#_Q(ZP#`,$###gm<#013J-###QE*98-)c$###B##PG#E>#F5#/,#"},
{204.79,334.50,4.61,6.23,"{k#######5##b(>######f##;']eY####I:'J@*#d(2##~kDG,$95####(##)&]######{##<&]95####K<'wJ2du&###RQE?#$+u#######M&]######/##b+]zo6###=##XB,/&]###5,#eG$.,#######9YJ######-##nHEuT7###&##AA):QP######"},
{717.26,335.78,5.14,0.86,"{4M######H(#}?T8.*.,#H(#R]*^..###)##>H%.,#$##/,#G,PU##(c$H(#o&2X|-]T648#RCT@]/5l$Q5#BK3.,####A,#R]5g/$lw0S1#&K25('Qp7#]+m?TgP#@l$V;-[/3Fu$.,#R,#Uu%Hg$Vf476$g?TT?$jH)If%H'8######]f'R#%M5$###Z>#"},
{401.85,339.88,5.60,4.44,"]?)######/##1%Y######2%#4?RD>#{G#d/,Xu%0,#)Z$i7)@A1###.,#6##h%Ydu%###r$#WlHOwQ###W,#~o#lM@0Z%2,#Gf4######o>#U&YKJ+EZ&r##*K+lVP796%##?NSWWA0Z%4##nf6######J##LFImu#7o2'-#;w-vd#JdV1##J%YX>#z,'A$#"},
{279.81,355.92,5.03,4.72,";T3###$##3,#T;<WZ'%##Aw'>&Y9Q&+##mI)z%Y###{G%%##71;###%##Q>#G%Y###au#5;/.RS###((.5(5K&Y###:l$2,#6g6#########I'Y###]c%<c#jQQ@5#O[HTc$8&Y/,#Rc%X>#^R+95#######k'Y$##B,$###'cIl#$cPM%##Z8UvY#a?)%##"},
{588.07,401.03,5.08,1.45,"<5#>,#.u#IZ%[l#3L3M5$cP#kA-h$+###rZ$d[,###$##Bv$###Z##dy5#v'd25+UE4C^L#$gG^jp6uv(hf**R*&##X.)4A+###$##~(1du&$p5Ol#)H^htC=C^'R'pA+2w;^.-4I&iA1pP#######w6$=6(D>####5w#]|DWw-tP$B##Q{9)6&#Q$:5#Y#$"},
{377.34,429.94,5.44,1.48,"1?&jP#ZP#2##]m&u(8{k#%##Mo*OS0$##~G#Z>$.,#+##:c$X>$p5#W~0A#$.i5*qG,:`4l#U?`]{<-8-sw';R*>5#48-ZH&N,$RG#:A-sb#kB77Q#W?`LO@i:`e-'(10Q@=RI*c['F06pY#;l$.,#+##)l#OZ&###6I#Wo2jI,/v&g,$)g01Z#E'2Tc%pb#"},
{126.59,463.35,5.67,1.42,"95#-##G>#u>%<c#&136#$pP#nR'an/'##E#$/u#.,#C5#D#$3,#<Q#-/1~>$c_4s0EKp`cG#9v`W<>N8.(J')m(E5#,x,kl%7,#H#$X/1iY#;U8@-$Dv`%4?Jq`kR(k^1fd;F[)lm%WS0}G$$##.,#H5#<Z%OG#.,#{I$w08nI-@#$i>#3(4Ll%pY#95#~c$"},
{728.02,521.89,5.69,4.54,"%I'Ym*1,#ol&P&TxY$###&d$8kF.,#B5#k2,=$(L,$rI+dI'=c%###$##Iy1BB~######]7(MB~%##eQ'kA&^?)C8'=oZCv%9@,######Yc$yRY######Lc#fC~qc%jd,|##'.*0q&=A~4##/v(###*##|?)Ah<###&##<%&xLPoc&-H&'I&~d&;I&^B~SG#"},
{637.17,523.73,5.21,5.66,"#Z#eI+95####xN;@d*###+##SXBZP#2@&i$'ZP#)##W)7E,$+##f*?ml&###za5oRU}P$D5#MWUKI)m@-$Q#}c(Qf$zPM;5####=Q$Kx+`7.G=H|H(Y8&l^1_RU+A$X^:l/$|#&e`*H`C###.,#2##tQ'to1yiE:##I].JB*wQUX.#WQS#]#0##.`0C{A###"},
{415.67,544.28,5.21,2.11,"*B-[p9###VH#GM*W_?###.##:s.YbH;Q&gZ&HR*|x.c01+8.fe%fOH)Z$)Z#@Q;s:>###ze(hSP+u#+##(K,vR-fZ%>/'H(:^B6&l#KS&#C,oQP1u#8?$_UO:QPJ>#WQ$Rk<q1;.,#N##I.,D/2###iS&LRPiRPd$(z5%~{4G.*msCZu$KQ%(]'/*E###&##"},
{303.14,577.54,5.00,2.02,"0S(9d)~G#SG#8S's2=c-(+?&M?$y00Zf2S7+RJ$2^7.R(}5&#B)[aDG>#$##fRM6m(P-${:9ml'l,$3M'^#L(##[7&kO<BR+I{:S^9&##l)-H%S+u#2##ur1#+I$##^g#6#JZP#3##XM'{#S`$SfY#O5#G+24sC8#$]5#R*7`FFOG#+##VH%}Z%O5$c##$B4"},
{672.70,579.50,4.95,2.91,"tc(%##E>#;[#PA2###5##/P1mP$###P%$[FBCZ&###B,#r[&2p6######3K%W~YdP#G[&g,58R*;Z#raY@X=%Q%%##gn,Lm%OOC######S##raYB7)>-(XQ#?B.NV*~~Y)c#*6'=##9&2`5#?g8######K$#Z^Yk>$(c$R##re,f-%uPP###u?**##w6++##"},
{789.32,584.36,5.15,5.78,"Bd*.,####7},:o/%Q%###gG2Py3fY##c#aP3)H%K>#Il$~A-C{@3?&J>#`)%C<7#|>y#&o-$w*~[-)E$)Rl#~-(QQ$t7/_P#^#&UG#+Q#Gk?82?*6%]q+}o18&~_P#^/1&T$Qv)nP#8_<J5#######OZ%-{9ZU:/##jcIjd*f&~n>#CT4z5#~,%7~$8U9%##"},
{185.10,593.52,5.33,2.03,"AJ)h6*g>$#6$@K&[sB`c&r?(Ce&-)5r6)qo-o0$d}JGI*h#&SB'iGOrb#8,#6/LXv*hP#jT0~@-}?&d/)4`9*##dK)r<;A%.1^4Ow,_5$9#8=IP6#$pl#t~F]NE&##g%$E`<eY#0##{C&vGPgIPXl%6?%gj1nB4Kd'L9.d*;'::###?5#nm$VZ'###7##@n,"},
{284.11,613.63,5.91,1.97,"+.%LR)F[)2-'a&&qWCfl&_,$57(o|9{v+q]1Yh<g6(X,$R^1gY#su#,;-0]2J^5mB1QL2h81FfSf~.=7'M]Gt/3-m&`8*ka>d-*%##[0$t2B9v(Y5#c:&ddS|bFpE7D-'CvH{7&;eS/Q%N>#Fe'Nc&(##,c$S,#2o%[A1qm,RG#>=0lWDj?(~~)&cGO5$*y("},
{284.11,613.63,5.91,3.59,"XL2cD?kc&;l#eW;o&38d'~A*@5#8n%.&Te5%<'(cw*x+C*Q%R^3f[)F|73w,HfF1'Tqe.[R'LdIJbC*s4F$()[(E>#i`+75MY#$K/&QOF|#'=o0+z6/C3O^22%T;v(T5#HL&a@.W,#o.&$&0IQ&;%%@@*#@)<f2^P#NH$:q-4XFx?*'##tK$e5%^w(Yl&)##"},
{229.82,630.68,5.39,2.65,"kP#z#DZP####_$'7dP95#w>#d.+p&4ZI)r-&qY#N5$Bf,S%.d/2R58|]7=##LdPVf3$##j_((x/6#$nH$CAE6#$###|d%'uI506s>#h>EyS)JgP6#$.##Tr.D*2u)CAv%PI)<-&nkHRQ%W$'|6*5c$:L,I^-o^2ye0Fd&<Q<Ym%[kDg;<4R)46%|,Lzc&K?&"},
{404.48,636.24,5.48,2.15,"P$&gd)o?)EH%_I%Qa@3.*wc'<Z$$b>0.+g5%w?%eW><[*$##+7&_6If5%D#$etCFA/Ld%dM<6v(.$%di+CvRj?'KR%8;0Q#PcD<5;<zb#,(+;wR7l$1##dt7DGL###w.#.wRl>$<5#]J#5vRRvR;5#)6%.s4U<Bg@-5##823{M74M>###'Q#4%$AA0,u#n,&"},
{404.48,636.24,5.48,3.97,"4Z$)p-jT21%,V8TTQ%Pu#I~)i7TaH%0H&-#####K[#_6T###Du$5u#vC2l&277TzZ(}P#cy+i,KE8SD>#C##~c$hV9C<@###f6)BH$$P<IS.RtFZh4{n-&$%r1.$;T5v'$##>*3M8TzP$'##ZR(fm&vn/f-(wG$BH$SM2{?)`Z#ZA.EdF[P#G7&86N3d(wv&"},
{193.94,640.84,5.00,0.49,"/%%^y0AR*:#$|y7Aw*;e+uK%994|5&O,$8E*)l#@w)Z02JC7d'.0(5:A0p#%iAT)7+###c]#-3@Q8.SS-oA%)##H_+v'9###P[*kP#510(CT*VOo98GU3)^3aW+%ATj,&S>#R,#*c9E82%##%Z$86$y02G166Z#I$'nuEYXG;S'CJ.md*?`<Wm)#p-oZ(:c#"},
{193.94,640.84,5.00,5.76,"F>#W_))#>cQ(m(:qK2L-%]w)pD8;@*@l$<00u./4@&iI,*'0Ou$B_/?[<Z..{DA-6$&y,8SAwv*:u#1?$l)U0d(1l#v@$}$Ue6)|y,D,?p5%,(U)v%36%'&'nW7PV9`u%e~,/I%nB3_Q%Cd*&l#fC7kB-v>%)j8te0@H$$v%E'*>-Q######rm&Q`<DQ&ZP#"},
{724.93,640.10,5.62,1.49,"y_/-iA######{IgM$*###)##?s@######~,#{7+######,##VM>OG#0,#9,#zEg######)##@AX>u$###%?#Bd&eY####gG#%V:{k####B##[Dg######w$#ZSXT?(1##^c$:8-yY$&##hP#fT6.,#,##=r6PCg>5#ZP#_2$rg8Sw*,6'*?#Yp1Zu%.,#.##"},
{647.25,661.81,5.36,2.51,">q%uXH`l&$##}('&ZP###%##[p&Vp9######Ce&5?'###(##EW=n'-Z%/{##Y>>XGDL>#?m&p(R?J/%##@?$x9795####>$#~::]-'95#MB%)g5E@)_R&f^N?$Rvb#+v$(<3e:9(c$###P###95###$##y2+*)@C###XBf2*%5Jd,%vn/5$%-n(e./###%##"},
{183.70,670.36,4.95,2.08,"_G#vP#pM2OR,nv$_::@K-YZ'Fd(Jo/Mx.G7,Dm#*e+}.+gu&R>#m8%-fQH?(`L6O@'tK.kL9=@+uQ&Sh&lsF8,#uo'fM10B5'v&mb6^>O}n.sdQyQ%nH(T3:/aF$##LK#N$ND>#*##3h$PcQ|eQ@//L#$wg2P_;H,$u?$-@H32@###(##Rn&_Q(###4##iJ1"},
{183.70,670.36,4.95,3.99,"n,&h#$Sf.RK2OwS1,#8H#s8,mwSvH)###0##/##k82######lA2ZQ&}-(sS,dvS/[(vl$B](75GHxS###C##L9/|OH.u#dG#EI+wl%pD9rA/6cK{J-tQ'iR(`&-jfJYI(KL:L~ERk@#$&t-*)-'G@&1T2[7(dK74,#BZ$XP836'9,#B/'KxSnu$~7*VeB`1="},
{387.21,669.65,5.29,1.86,"X#$Ef/uQ';l$y7&&`<WZ'(##OT+c//Iu#ve,@s=JH&$H#.%+ZP#WZ#M*,@f3;_9|I(C^)6:8TxUl#&1n&Xk8A]3uG$6U*zz5#m'###<:#yV@|c(Y5#hq$9wU#5GHq2w6'@@L18(o*E/6&6H$D8&%-'(##Y>$`>#3V&pS02I+D#$fd=D95$m'Py()RP.,#cu#"},
{387.21,669.65,5.29,3.65,"~w,:YC|k#a,#G90|A2]Z&)x(hP#Td%RyUBu$@m%cn,ct;av*r(7:K/1J+zl%,O8zzU*$'<%'0ySl{;7U*WI+u?*.,#WD&tvU2##I/(JW?9Q&6C7g3>Ee)p%*ZwUhH(F5#Y]$uJ3)H#F.'um+F>#'Z#v(9%d(7T4qb#K$$Ph.evUHH'%##)1$x#':~(H?(*##"},
{452.89,680.45,5.71,0.99,"0?#h:Uhf-_v*|,78vSP>#Tn,cn@.H&###&##,$&L>#wP$###'l#k5$XG4_FHZ4DQG#.T(G[PZ{W###`G#Ol#E:8###yb#`P#of/o@*X%,#'3a@.GQ#N`7yq8(wW###Nu#^B'K)A###$##c5#?*:^m)ZP#W,#G]1Yl%nu#ex,jXJ###8##oB&5J0######p##"},
{343.82,697.60,5.66,1.07,"I,#pqVxp0;-(lG5$nVBu#FI)E/If#&)l#Dd'^#&$##UG#YB.rb#Rv#3sVKr=r3ChY#.o'8?J<pV###S>#_u#Dg7ZP####dZ#k~0C_+pq=:e),x1pl#7N:se,nmV$##6Z$:.$X(;vc'xP$'Q#MB1@v%om+t>#}.07-$N..EH#SmVA5#,u#}.#Co0eZ&l5%qu#"},
{166.52,704.81,6.05,3.53,"+e*Q]R6l$;,#Pt<%V;zb#-e%M]*YQ&2Q%^5$d>$@6%>x.0c$8'1eo,zo/P%.Y&G/CVT$)IR%REVqGJW5$Gu#ve0E>#Fm$|Z(V,$>U)%]0R[)pU9]:4bq1_^.GAVW?(vY#ah${/4Tu%0,#@,#{l'<m%$&-r&0kK6|l'd@(q{*f5LuR.)##f:$9m%8.-###(##"},
{489.15,767.07,5.28,5.79,"8B#@/]######tq%JU:###$##zu$3l#P5$.,####4Q#*c$###$t+U/]wP#Z8/O5]+06###R7$^W<nQ%OG#.#####ov$eY####a$)7#$p-#g3]n/]95####DfAF95K?&/,#K?$###Xl#/l#h,&P5$###6##v,NZP#######y(5OG#######h7,######'##[7-"},
{470.43,795.59,5.42,1.18,"L?#79N6?'###B3,g1=.,####_n*.,#######_Q&_P#95####v~,t<1M6T(##g;Tw$+pb#J##iaCl5%95#&##Gc$lg2+u####hd,Q,#2(Poq7U6T&##iP#O*-<M=Ru$/$&}7'*##%p*TkE=Z%95####Vf$dPK&@+###2###v<Nx3qP#S-(/x%8#$EH#}7T$?%"},
{420.77,837.68,4.98,3.63,"He*jp7###Q##H?Mk,&###2$#|p;######=&#}>&######;-#q9+6>I9q3UG#fXYP$*UG#q,#7#O######Xw#A,$######h?&U:6bw-])5)ZBvSYD>#%Q#c%?D@X######k@%.,####$##z$*j08a5$aG#*;NcI-;v%36%<WYsT2c5%%##I[%9c#{k#'##V-)"},
{679.79,838.62,5.34,5.85,"8,#0?$%?%|k#+6'D5#Td&77(ed,,##-n'.JG=#$|b#9L'>8W^5#Ul%_G#2v(d{B8Z$Nl$}l$?8WNl$XH%[$;l$*T%%K<WUOD###'##TQ%981}3HB,#0%*:6&x:Ww7+^#&T##vv'.a+&7W%#####1##rv+.u#dNDV,#}70?##:8Wi#%KQ']##Lu$0I$6ZQ###"},
{767.43,839.32,5.21,1.11,"###4##-8/###1##bo+HR*###i6$x..H~$xY$N>####}($r5&###OQ#iB8###xv'6#4fdW5,#DjWz:76v(M,#>I)@8+u@&###&##rl&;'6###Um*ZZ#fiW$x,@eW1n%a8.}L-7$(K~@@H''##P5#Rm*i5%###SG#=#$vC-T~/Ru%a6#MIN=03,6'B))/lOaJ("},
{767.43,839.32,5.21,4.36,"a.Xqd$H[+J8$lYB[z9Gl%M##tU,Gw-I>#Hu$0c$###+##i[,0&28,#:$(SG3U'1Jz/)/XE[$t3X9S,-I*;Z#b]4######Wc%8R%pb#d-(3g+e?)>,#i*X3h6G.X(##+n'>c5)_<######h,#vz#LK6B5#2,#mn#~?)]$$o@-bm).,#6##ef,Tn,######/##"},
{684.28,859.98,5.45,5.87,"@,#|P$Mu#P~0l./Dl$}#%t6*Oh]'H%rY#)]&F}BTV++LZJB0###%##*?%te0Z&4e,#NL7K$)>h]:?&|G%-$#5`9rM-oe])##$##&##(R)<l$|K6^,##YI'##{j]R6(`#&O##Tx-N_.AcQ'##A,#f-*3Q%95#y#&Cm'~h8.,#mt6w..mc&E>#O$$^vB%Q%###"},
{734.56,863.77,5.77,4.25,"|V:_E-R+G?Q$|7JvI&./0Id%ESG3w-&##ea4;~'LV?U##mmQb*@UI%tSRm>$$SRxG#TD=]##%?;?{@N5$x>#t%$EQRZP#&##tV7vs3U(=-Z#zQR]$%#J/A-#l]1rd+M5$C5#tY#)-'N,$zY$Ec%l$(2v%qh5`o4--%`?)^9'WL0K&/]?)6,#'8-.,#'##q>%"},
{438.44,873.76,5.54,4.16,"_c#%Q=@H':,#%j-w'9######,R&sZ&95####R]3QG####%##z-)'_*X3?O7,6LU|Z(iP#][&Kg4uR)>u$&##gL/GI+.,####A,$*##9{&5IUf;@)Q%l.#$KUa$)812%8&^T6g9)8^4P-%tu&QZ#Q'8Qn#euQ($$3A-Nr.4IU%c#/C(fMU;:;J'4pR)fv&}e*"},
{695.35,909.91,5.38,6.09,"@,#q$+dP#QG#4,#78+4p..,#&]'6n*r%+iu%-m#r<4qH&5w,*?%:.*g6),c#>H&-F2}K83,#4&U>p1q,&Es1ae,H{8ch*I'U95#2uE2@+%c#TZ&8s8wiCf#%bIF7%UvP$/m$'T(v&U]c%qu%/##f%U######-%#U%U[P####T%'?wR###N,#aw(uVA###O0."},
{219.75,971.86,5.50,5.38,"gG$###2,#Pm)###;$#UZ&17,5-(je%^?'1&&k>OlG$H#$fh$Y~$`S0nP$I,$Xm)*D&o':|b#dnVDn$}>&a(##PF(f(3l$a(#Oz'xf5y$&Ld)`'2Q~(9i?Xc$3sVdS1QZ&??#5j00sV>u$@##:m%D>#a6#,h9(6&0,#f7%)]0SR'_e03?#?o2RH%Pi@$##;8-"},
{214.72,985.27,5.16,2.60,"c##PV=%##95#Y##J3;??&L5$j-#47Qgv%$.,>e#Nn/K%#`V@MZ#&//Cp-R95#)6wS-O9/sR,eR',:6c:)M@Vgq/$.,h5#k=Exb#H,#REVANBu@VP>#xF@)J'Dh=VG#]6#~,LZP#S>#RJ',<B+u#2(#9BVx5%V*GSZ#v{A}e#x.0KR)0,#a,####9d&ll&###"},
{214.72,985.27,5.16,4.40,"vo5ln&Zf5C##&@K8&,VH(&H#wo5###yG#&@&###$##Wd&{k#hJ29[#^7TWc$*8TPd'HT2H7&I<AtkJN>#E5#L##-OFmc%###M$*'##Is5$7)JK49B-@y/mT/Um(@9TQ5$~O3#d$J<E>##8:T-H&.##^95D5#sH(p5#D)=oH&|7T~G#'l#dX0Oi>2-(N##g;T"},
{139.56,994.64,5.44,5.03,">M03A(B1Pa#&B1P2H#bM<&##e)4}>&V##yQ):l#^>$?l$9-'5p-tX3|5Ef/0a-Pdm#2:5qn&H$E'7+###7##YA,Y5$)c$###4nCsZ&mP8bh86E<.d'M@B5e)6>B'x0###@##IV96Z$.,#(##'R>hJ.Y7-/?$/?$1n(Jh:,B1Ly63l$%##ZT';A0<5#0,#%$$"},
{212.03,999.57,5.04,1.33,":#$-Y-uo5jH&J^6=@)Dc%Z-%zy0GQ&8Q&'##Q>#96#lS^###D>#a$&c-(?AQ/.(C&+,nTYA-ooM]7+?@,9##oZ'~6#lS^%##jG#oY#dS0-~-qZ(Q5#-W^+R'-T^8,#O]1k/'uc(f##oS^1##V,#4I'+6'###pG$1I$K(7*d'jI-^G#F@+QS(6#$$$#lS^'##"},
{212.03,999.57,5.04,4.39,"b.[%##eY#_%#z%0*/(OI,^5#j04Zm(,H%Kd$6?'###V5#sI(d.[4##LQ'e%#6K2.p&p7ZF,#)1[`m'5v(a5#Jf3oH)xP#6l#b.[&##qH(CA#Cd*F##l]OJw*IeWPo/:R'Dx);@)8oW.,#kc%b.[###G#$-7#A,$7##mz5i#&|5&[@'@f1tI*S..{d*^P#A-7"},
{131.40,1009.47,5.41,1.06,"###A$#d%/PG#:5#(n$j/0N5$[$(.R(q$*?8(;m$w6(x-)hv(:5#*T$dC;95#9%,',7:d)$m&#IOP/+@f2|:*fl&oC(DHOtm(TH(8C#2+J?[$o%0P8'Vf0ILO^kB?h3>:6::6Hu#=p*+HO2,#@4HH@#9V>a-#8HOr~$%HO}U*k,#Oz-xaHL#%X?'<S(ix4i>#"},
{378.07,1016.62,5.76,4.52,".,##########*GK######]##&<_######D((*uH###95#z[$.,##########*ZM######&##|;_######x##y9_######q##############g}L######D##Y9_######%&#|8_###&##9I#.,##########[kL######i##%:_######T1#=+I###%##O@#"},
{80.20,1017.25,5.62,0.57,".,#h5#Oe+hQ'zG%B_&wQ)R6%e/4Pz%mP$U_'i,&8x%9~)snI%##D-#'g7###1w,uT#LQP=-'3SPIp+^6)Op'zv)iN.un0S6$0m$F[$PJ1H>#Qe-iH%2i:jeN&T0GB+=p.ASP~G#|e&oD9qb#?D*qH&:m)$##9C-oE=3m(l?&xY#cKKEv(Wl%u#'Dr1Bm)n6%"},
{80.20,1017.25,5.62,1.90,">%JPG#B12dQ$l|4Iw,PR@A6&/,#SZ%ni0}]4C.*:5#;q,hR,/o&Uu%3-<.,#I)*L6PM:.dl&tlL818z90^C1Fu$E,#g8P67*S6'SG#`(,8-'6d(dc'P1%qOHk7P=8.l%*pL3sH'fT+*8PTZ&OG#$##R?#p?)###&c#$e#5w-#H$X~+Me%B$)*##Wh1w6'_-*"},
{125.70,1036.63,5.06,4.39,"<?P},$EM;/7$t1<hu#|q<.-$VE<~[,.##3Q$6?PDu$[],)]+8BPuI(51:uG#/EBJ?$8(:R~'pV?fY#h##@)7n?P?K,Yf4bI$Y~KoZ&bG$2,#a@Po5$nP$*6$nOF^>$LS*t)72x,elAX:6DI)lcL######G5#??P######LJ&Y?Pw$'Xh:p8(k,&~6$([Md7)"},
{884.25,91.08,6.51,5.85,"3-#?vO.,#$##J1&gMA######>q'n`E###%##:1$1uO*##0u#Pe'QwOU9*ruN*zOmf6+##j10p*A~#&###=##2I((c$###??$<u#w?(3M*kuOt94Vd*4o-li9I'*RuOgY#+##:_(_(>###)6$,##U@*mu&###&##Fw$ay6###}@-$r5Fd*M##@4Akx4###-`*"},
{782.35,105.54,6.53,0.14,"7##y'V######Ji2Z%V.,#$##Ce+).)X-*$##oY#Q#$g,&###@Q$h6M#v%'{;8(VTg9.,#&V/}iAJn)Bv)C##}P#*-&b5%&##Bc#c@U^G#eA0}k4:HR###+n(:{'@%VOG#(##,^#Gz<######Em#7B4###iP#Zn$V&4###7,#E%#D-T######Y(#ch?######"},
{782.35,105.54,6.53,5.43,">C)(I(1d&iY#3;+y6)~L+ec&4DM(c$_>#^[%i}/z,'6##A##F48+u#=##BQ%;^X.,#yw$3u>HD<0Q%`5#Q^HcS+n-)###Nv#u`@###$##Fl#hE[L5$###lv$#)4n./###p$''l#Kd)###J>#Z8-OG####$##gD[######+##wlN95####G##b5%Mm%^l&dP#"},
{147.02,120.49,6.35,4.63,"k84>##($(l##=<D######>##'K_###H*F$##W>$###'K_###pm,K5#D>#Qw#$r=0?%95#E##IK_95#]PN-##I#%###;K_###,u#gl$D>#[Q$pz>jl#+u#'$#|J_###aYN_##_5%###-K_###E,$0$'/,#IQ#%g76,####X$#tJ_###H<E:##(c$$##}J_###"},
{635.14,144.11,5.99,1.47,"###*##.m#k-+0##LS*<&&5J0u#%0d(Ie#=jF######l-#JdV###@##{02^Q(^x,(n=]fWP,$1iWjL7-J*}8-'$(H5#[U2Cw,######yB+$m(L.-Ac#DjW`aCFeWr>$P10GfD,.,sP#,mL*c#######>Q#q#'######R.#g1=DZ&/,#v-&NN=ZP#%##L26`P#"},
{368.09,172.48,6.39,1.45,"###Yc#%Q%###C##bq6Eu$.,#@@&x%0eP#lY#+c#iY#3u#:5#####6#c'8###Z11C2Mh'_X,$},_x*Bfn,Z@&`v*?5#-K/Pl$.##T5$X94###^&3@6$},__4B#(_x?'t0/`A?47+N6%/04s5%3##KZ&)l####PG#+l#=e%3/0t@.=#$~5#EB/w>%5u#_P#,Z$"},
{170.77,207.00,6.96,4.58,"S##{[P?m)###X#$u95&tA[c%XH(lY#m%)g59`w,6q6mn+h],;##Gf1|q'an0Y2:m'9gd$)92(8Rq[(_O>Wk4Av(7~)='O;-C######3(#AM?NI,95#-'#B7R{-O)A)<f0[}?kT2Le(1f14m%######($#681######i5#'+G$##*##I{9Ox1z[,'c#Ao/z-'"},
{74.99,257.62,6.51,4.32,"###bA#vv*uH)###k')2&*}S3$##>y0<.%pDBR##)803-#1p7Xd+DB$rL9X*7:i='L)u^3}V=$[)###xz+)GJ###)##c9*h&5>1;_:%Q[Pvw'3[P-y'6U75U(hx4j,#HOD?l#$##WI#g)@###0T3'Q#(9Iis=ed,[>#)L/hwHI#%N/#s`DW,$L,#6^(jn0~P#"},
{649.50,263.93,6.24,3.62,"0|6Pc&######[-(M,$###07#:Q&W5$###|(&O5$<x,###1i(qoS.,####+##gmS###UG#wC$g6*$8(JA-O{(6u#EV3Y,$#S*wnS###$##*H%QqS###L>#,?#X`;_l%K])s]51,#lw(W?$ikIxB6###M5#,d(2rSbG$###$##(|1w}H&##Y>$'##9.CK>#I#%"},
{320.86,272.82,6.04,1.47,"%(;######A5#_ViH,$###(##JjD1M65J0###ML:}-'^(>###{z?######A5#hViX,$Zc&*##R;=b)8=z7x>$p<E?@)G6(pP#m(>D>#######)Wi2##S$*%##VE@G~'+HOiP#T5Kr#&Y,%E5#GB6######,##hVi)##.d)$##{2B/Q#.wW.##&>K1,#*c$1##"},
{259.23,280.39,6.05,1.49,"Cq<######+##%MdY>#6w-###vL;VZ#8ZPvG#?#N$##F>#f>#VV?.,####4,#YLd:,#bm)eG#Oz;.$%_s=Y.(3GKJ,$Jc%'$$|1>[P#######VLd4##l$++##sW@7A(NmSfP#b?Nm#&)c$&##a]4eY#######-Ndo>$S6)$##d)@Lv$kHT6##M5NRG#xY$>##"},
{496.08,295.79,6.82,3.22,"%m(ZG####/##*.XI,#%?&y9')n-<r'72=X15###]b?m5%F~.KJ1/,#$##:,#G.X^G#N@*Sd#/f2Yw%TfI6-P%##H[(BZ$:1Xxn1###(v%+c#:0XSu%<,#%/%iJ.]~/$m#16F2,#,$(&##SJQhS2(##tw/P,#n.XT,%###h'%b%.mo4###_<60,#IW?###yo1"},
{688.80,323.59,6.95,1.51,",(67m'o#'###R2]d5$nZ(###TM5&f-Lg7T>#$N6,c$/,#7##RC:L5#&Z$pP$m0]%##Jd*(##i4Fg,$#dPgG#lRNRG#E,$N,#E7-###Y%#%0]Z/]&##qZ'0b5:)?H5#{uK+$%CSU###Qu$%?$6#$###u7&ye[<R+###7c#).E^#&###9K*UJ.5T3###yP#^/,"},
{489.99,327.05,6.85,0.62,"+c#PV-O?'[v&W6$%f0MADc-)%#####,d?###{k####WG#:,#fR)AI>ZP#HZ#M_S+i:H{73-#l?)XQ#z_S2,#6#$'##)I):,#AR+vg%OG#g(%W[Sm*4~Z'|8#^c&DE+qFJ%##L5$?##zu')##X>$f%%{k#%U&#l#]i*###M##ZP#Gw$OG#$##M5$'##D>#%##"},
{489.99,327.05,6.85,3.62,"OG#'##6#$)##.,#'##{k#{d&###:,#;5#cj3PG#Sn':5#4q,5-(.##{k#A##kjH*##9H&OF0r#'{T$bFG'[?###`s.YZ'K&)pd+=,#eY#(##gCT:,##w+6m#ZC7g7#gCTa)<.,#)R$fv(LULpG$B5#95####'DT###K>#0,#C7E@R)f.&:(8hY#lH&fG#<45"},
{429.68,337.56,6.52,1.53,"xR-&##W~V%##CwSK>#g./$##Ql%4,#HHQ######$##s#'###&e-###{~V$##x[V$##[/2_5#b?)$##L~VI5####$##n6*###Mm*###W^V,H%k[V$##F923](>[*&##.]V46%###)##z[-###}Z)###_b5erA<R+###SK+M`V+u####f^Vq[)/,#$##g~0###"},
{429.68,337.56,6.52,4.55,"#[)###:5#&##<fXRQ&.,#)##}D8]1UNc&###Mj/pXG@H'###w6+######-##PeXJc$VH(B##0'3H:+qdX(##4gX`-(n?*)##EI,######(##,eXV5#d$+:##ix4H6#rdX+##pdX-##&@+D##<R+######)##<#P%##_d+i5#US17##3eX9,#udX+##>6(2##"},
{763.68,558.44,6.49,0.25,".##ZhMmP$###G*2-sB.,####12U{,'###4##:NBE>#3,#+?#F?#L04######%s9{1:###*##CiZ%Q%###;$#yiAHl%8u#<@$FQ#Er</,#;5#veZn'5###x##)gZN7,OG#n$#?F0})B{P$hP#E>#yP$'l#<H&bnXF>####'?$XgZe#%:#$>c#%]*Av'Jl$vJ/"},
{188.61,576.03,6.92,2.75,"cJ+]i@96%2%'He+%l?b~0ol$5f1(R'?&-q]0j#&vY#/~)e]L~~*Wz6,u#~>#Hg-foPcG$G##;qP8e.*##Kg+Nq5###'##ZL1C[(&J)?K0^6&lPI$g,yz86v%;qPwR-&H$&E+_{<GQ&_5$R6&Y>$Vu$NT,Mp5LR)&##{V.s{>W^POQ'k>#i</Vg7g,%4H&;.#"},
{768.30,580.32,6.32,0.21,">l#ig6######NM]X7-###~##{K](c$$##hR#%I(y5&w5%+6$J.-#q6###fG#BK]cZ'###f%#8N]b`@Z>$06#xX-qA3:6%1u#^e0RG#%###m&mL]QG####T##6GFbd(,6&G&-s>7Dc%@Z#?93:m)###0,#f>$'K]######+6#+7W######9'&>6&D>#T7$z8."},
{380.42,587.76,7.89,5.64,"Z6'Cw*PH%bw,o8.66&#A'5*B+c$0##K&'emQ*$((H$Vv'PREum(eE5pL9l>$QpQUo.gG$MI'hK5*j@`u#SU9cG#z]2J{9{c(]l&)e#|_9knP]ZQCB,4'42P4QC/fpQRg76,#nN5eq<%d'*-%X>$D6'3@&7YI###0S%,mL9~,A[&$@$|sDM]/L{5dd+kG$v@)"},
{396.48,596.01,6.12,3.70,"$L03jDpb#],#-')@<B.d)$##X5#V7(6yR3l$tT*nH([g/rm+ym$dp4c_=/,#DX/#vRa5%&c#%AM0y6[(+Xo17v(###)V$OvRYG#N?%)4@<J-x80<w-lc%+yN_vRK#%3##V)-:~/$6#Hn)PK41v(/,#,n$x>Ec$*9Z%Yv%GQBEWAHd*A,#>9)EQ&:&*-S/*##"},
{267.13,607.62,7.18,5.72,"<6&Z8*mS*q96{]0ZH&7-%vjEzG%###w5#3dP.,#=##kJ,eU6rm(bs06dPk>$AgP[S+EZ&2$$IK4aF>|Y#~-)4##6PEr//A,$<Q&te$C5JO#FRcP&R%mw+ft5_[*@/FAG='d'./*vjE=u#`?'dG$5,#-?$auMLZ%K~*p&0t'2yH)oy*gcC-A+]6)>5#ec%2)8"},
{408.75,621.24,6.44,2.39,"i9.eB6yl&j,$|l%'y0a?(.d&#c#`10Q-(cw.=l#Te*i/0Vn/@p3m(8SZ&?R'|S,F[F(Q%/Z#fnQI?'Gl#?tAQd*/,#n[%dmQbP#,L/+l=nl'~mQU:7t6'LJ$-nQmu&###*D&Hh5H(<###/w*#l#L5#*6<B_=$OE'##dS)MpQy_=.}4*80:&).n'WoQL5$O##"},
{712.83,627.33,7.06,1.59,"qu#AM?###%##wAHip:######QW~6#$###-##o.,95####0##Dl$Hl%2,#=5#VB[+u####+##pV~ZP####5##`o/3l$###)##S>#Kl%Nn#VrB46Q+u#.##-W5uS~<c$pb#~J#M'1Wd)L5$(##$##Lm$}/CaS~DS/$##qw'm$K65Enn/pb#Me%;01%o.L5$_>#"},
{818.45,631.06,5.81,2.30,"Iu#@5@I,$kd*5C&5WA###[G#2q&tc(######Uc$OG####%##)//^|=dP#nI%4(-v=B]A-D6&MWU-x0V5$C,#NB0.,####$##sQ)H>#{Z&t%BSS0=?$MWU3E8jRU(l#,~+.p#YB3pb####8##B7-2$#nx45)-y&5zJ#SRU+e%WvSOm*FH'E8#XZ%wA4###*##"},
{297.58,638.39,6.22,3.25,"U-(5T+c;3@QHa~+<$&$0-_91Jg4A,$HH#IO20~-Fd'8l#]@#N.*Ot=c'4~?'.1.^'2{E@wG$xKM;e.f>$be%am&0z+Kn.du$W.(^$)p$*5u#;o.sA)PnLc5$~pN#mN)e+p~&Fx+HX;}lNsl%tu$4B/F,$>5#rP#ID)w'7(?&lM4*`=H%*&1,7QK@@,TJ)gK%"},
{707.56,643.24,6.81,1.62,"-5@}I/###$##:6c+u####)##dq8+u####8##o6%(c$###%##P|A>u$###,##a4c95####L##R-K]H(###-##-8,eY####C##%196#$L##_,Ez0c'Z$ZP#/;&+YC~o.@H'@##<GCZ,%###B##t#'.##6r/tQN.QEsI.6,#=L-=0/G`=(c$&$$%AQEu$###M,#"},
{183.71,655.17,6.44,2.23,"^#%jY#eK*Yh;z>$b?(W*47I+c#%,f.[e+880*Z#HR(&00xQ)WI$W:9YT,Rf3%m&d1+%]P07+0#J4-%~%*EFAml'G,#Gg(`ZP*~#W~PVw.6u#]a>Y0ML^9WT-{ZPF5#Z,%KZ;n':###3##qU6H>#J$&fM6(EBqPFJ6'1$&*9LKjBwb#/,#Gy([w,c5%###F##"},
{183.71,655.17,6.44,5.68,"W,%k5#p@--d':5#BD)zSMM5$#01hN;S9.5y.4p27I)Oc%,C4###Ul#uM>D>#JZ&g{,$BO3c$5?O07&G9.`v?D7+16&<v#^AO_>#PAOdw.###hQ'|O@*7JRu$]BOU&2Tv&Cm%g0+'?O{b#d,%Uo/NZK|,%9H%>~+dh;[u%lu$tS)+J.`u$Q/.Hg-uw/_#%9l$"},
{279.81,685.93,6.38,3.69,"-v%2IL###E##4/&U}H4-(###sG#]J*<LT>u$VS*{81722t./<v#V385R(@Q&>N-YITS?'V5$dLTt^;E:)mJ0@[*/,#V_%:IToP$&$$of1Z7*@o1Ln,RT-4.M~ITXu%F,#>`-ae0V-$Qw*l[,)[)'-$}R*0bA9o._$*oR)CV6q4B[&4.##@e&R$'8&,+80}b#"},
{313.04,729.07,6.51,1.39,"T|A5,#]~*P&02rY;#$TG#XG#ZN<%@*VG#)m#^f/c$*95##6#9D;eH&bf4>,#joYwY#Tc&O##9k<$L3@T4?,#m|4(y5D>#+##Gi;`/0|G%3h0knYeG#q$+n9%53A9R$<,M@Q#87STG#fY#4?#S[+eY#T##W1L](>###am'kO.-%-$##W37'{4.o2###:,#zp+"},
{295.40,730.69,6.46,1.73,"c7C###F6&P5$(xGnQ'5@*/,#Pp+uz<Sn-###DL.&[)######r*<}I+Q7,=,#EMRWG#S.,kG#~{:2Q$E5DYu$3mL95#WG#|Z$XjDK~-D,$z+5jHR/,#T#$4O/wy:###?h+WtC2_=###X5#ui4]:;W~.2Z#Vd>b*>'@+###3h*Zg2ZP#3##Yn,XQ&###'##`@*"},
{703.13,753.58,5.87,1.15,"*6#d6H-H&###Jz'Uq=0,####A~,D>####(##Tu%0,#I>#=5#|.*6u6zQT6,#0WTfA1oP$I##BcL@5#D>#2##2c$zm*pb#%###I*J,#0WTny3#RT&##;Q$AX.%$S=##7I+'-#;Z#wq2&U8&##.,####81&rz<Ce/###J##bY:OK6U,#m&5rQ#Zc$qJ'1RT*l#"},
{685.92,778.97,6.20,4.43,"eu$OR(Hh/b..&9Vcl&U>#m[&eU7m%.95#.##,##UL9######.w(QD-[:VP$);:VHI*s#&~-%yB8;n'no5N5#T##H5B:$)###ku&i8()_Q}e/57V@5#d$%M22cK6.6#a,P1##Xl#ua2'7V###oP$~P#PZ#0f.6[*###b>#lz5/H&'##Kr>hG#Q5$V5#1;VK6("},
{660.26,804.77,6.05,3.96,"R#%&$$_Z'z>${~,|u&###R5#?4B{k#######?81######-##G,$s5$_+28jF;n&HL25'0Y5$YMU8-(K>#Oc#YL:######n5#c>#@>4ANU_J2_p29z/d<=/J)QIUqb#J5#W)(-04Yl&###~Z#2c#n57GdP@Q%YmQ)~(sR+KM'N:;>-&k>%_J#iQ$p&4D>#'##"},
{463.13,840.75,6.63,1.28,"^A.Hv&Au$$##~Q$qPA|,'###8A*Cy4>5#V#%.,#@H$4u#'[)xo5Uu#b02Kn':U1y9.QgYGl$=hYxQ(MQ&c380[)J?''?#~fYB~/Zv#ky7gd%E80^S%_gYVI(6jY#OFjZ&99)rJ+EgY-Q$)T1PG#'6#qL4h5%n5$sm%H59Lw-5S(pS2MQ$so-O5#vOD5E=Uu$"},
{691.09,875.10,5.89,5.47,"$##$##dP#Lv)OG#%$#u7,-6'kjIh,#g$)-p#Y@SIf,b,%w(&:Q#*H%%-'pb#vm,'e#(r@8##5@SCH$-H&>(#1,J0s,^~1B%#3H#?W=I#%$##eJ0Ry1=y76u#EDS:;>`5%`>#DA)VDSk-+(##+w*'11X>$RH#'p2mv*QI%gk=1@&_|CX20Wn/}$#aZQsG#,d)"},
{241.07,877.30,5.98,0.55,"?V=zY$###wJ#Yn+4e.###3##j6&Rw.95#(##'##vH'{k####m5%gQ&;x-P~+-7%<JPt~1&##7*.)*D###Z5#hc#JK5###0,#:5#W+2)(7[Z&rjCTt2Dh<@-$~~TdQ'+m&*q)Nu#?$)lH'*H%KZ&IO6`%,/o.ZV>Me%E%*7{R||GL>#*h'<`T.,#'##YN+>.-"},
{241.07,877.30,5.98,5.64,"b));A0=l#6l$OQ#,L0qn+U,%vI+al%l5$rw.4l$###%##;'3*q%Z6QwU2&7*d;+G8V`P#;S*o/Q17+###;?$je0######n?'s>$3Q%?c6x8V7w)b?(%&(c;V-PJTG#Ou#fL/w&2A,$###/##=l$ZD&/<VjT4hY#3A(by56d(M80J5#B#$<u#37(W>$######"},
{177.21,881.78,6.30,3.58,"bd)5R)E6&<[&gPJW,%###/$#cL<######{%#o#'-##0Q%y>#?o(9&G.W=A#$}{U<S-,u#VQ#,GL###$##hI#pb####h,%S#$jL4G3<gH(}{U-xU_#%N,$D{R=#Kf6*/,#YQ#E,#v#'Q#%SG#i^6U%,H>#;/CTZ%Ls;f,%EN<^J%){?###F5#-n#hv+######"},
{194.91,886.82,6.99,3.49,"%?%I$%BEU2%,??&K6'9&-aH&h96{Y$###t,#Rd+######R$#X.&v`VZ]V,l##o)HvBa^7;6&h]Vb6(fY#yQ$e%0###/,#P?##w)[Q>xWC;X2Bx.>&/f>$v`V7tH/H%oY#x#A$S,Zl&95#M,#?S/L5#Xm&H`6-o.IT1J>#eT.KZ#&2:`P#qQ(PR$?d*###$##"},
{194.91,886.82,6.99,4.75,"JR#z,'######X~$?d*######U?$]H(E>#E>#%##0d$B,$.,#@v%T7,Cc%E>#6e$rdQ6.*kY#'Q:EWBI,$m>#>?&T%&$7+C5#r,&bG#>B+$fSR~,}K*sxOrn/cfS=7*.m'_d#2I*Dv%SiA]>#Tm'XGA$1QQgS/~+X]'+iS0T30MSB.,6@&T7++e)yd+B[)n,$"},
{92.99,969.73,6.03,4.46,"GpMu,&tP#VG#3V1Df1L5$'##ZC+67,$##T5$;Q%###2##MQ'$nU+##y(7E-$IW<N)4|mU$H#qrUW/0tZ(Y5#QB5###3,#K>#AmU-d#;(:U@#W:<L5#qrU]y1amU@5#oR(5=1(1:###%##&Z#fmU$Z#I>#fA#fA1>5#f[$1S-#A-$l#M,#bo,rw-######.##"},
{212.58,984.85,5.46,4.33,"JXH7J&w6+<$#LFA7f,KZ&fl#l/4###sG#$T*###&##Ud'-Q%;1;9.$vlNhH$8KX(n(g~/.@%)YForC=u#tP#A##X07I7*.,#+%-(##L49kv's96(y1QB/6],1e(IJXE>#sO7q6$[iCq>#/KXYl&.##oz;N5#u./%Q#yU96n&ZJX4u#;5#1A@GS-Fl%R$#CNX"},
{97.30,985.46,6.58,4.57,"}8)(V.Cg4SZ&Pl9*?&I>#;5#lJ([&26#$TG#i7'#?&$##|Y$Ei=ZH%%m;6:7%/RrY#w^3Ll$Ci8[U24QL8l#30RKH'=#$X>#(@Pr_.Q.Ru#%.-REm#uWAu>#*M=3,#;6=tK2^>N###;c#99+R%-)1-$26g.(M-Rlu#^>$ce#,92?5#1H#sQ(I6&uP$@5#zu%"},
{446.07,1004.01,6.84,4.50,"############lFI######6##=Fd######8$#V5K###'##l>#############9HQ######8##nDd######3$#eCd.,####;##############/HQ######=##^Dd######g$#`Cd.,####M##############NOH######<##@Dd######$%#mdV95####T##"},
{421.70,1008.41,6.92,4.44,"############2C:######i##9'c######l'#Z,Q###4,#S@#.,##########-IU######_##g)c######s&#wYM###8,#7R#.,##########mIY######`##t(c######m$#Oya###%##D,#95##########K6T######b##%(c######:%#25M95####T##"},
{88.20,1043.25,6.09,4.60,"fZNw#%wG$L>#b;;;['4Z$=$$GjDO5#iZ&;`0|I*vQ&gm)k.-v5P7l#|@)mK):6P9Z$y7*Dx%=7Pp0.6/0We%LJ'(W/h2;iQ(@8PiH%hp6i-%tPBE['iz:c,$A8Iv&47H%{H'$i=S#$PG7SL768PtG$}R+lP#vB5Dc$$C3]6%Q$O4l$###e~)JYGWK+.W>Fu#"},
{856.56,95.06,8.25,1.47,"M$*0$$w~^$##)##6(0Ge/###uY#r[)+{:Ed)pP$jG$GQ#mU9X-*O-#}~^###H.)0Z9tJZlY#lpYP93_$)yR(;c%###Id&p7-|Z)k,#7]^$##>v(xP#MDN1T1Y]^-l#'A(>W5wl'[u#-e,B5#z,'>,#A]^/,####P>#0I&@$)b5%D[%d5$R7+QG#H`,yY$$##"},
{856.56,95.06,8.25,4.30,"}>&###/,#/D6P,$'@(`l&Jf/S-%Im)95#Yu$X96###$##&m#(A.Jl$J#%O5#f%+M=1N-T;l#$2TxK22d)>Z#J-T%##D>#2^#$n%aT4)c$'##d$*9&)mRMd)<d-TmP#-@(=c7:-T^##1I+F(#{o0j6*[5#CS-)^7fY#I,#BV3Jx3###+##jr38-T^%#?*FF(#"},
{782.57,105.72,6.68,5.45,"i8%i.-JZ%D>#<_)&%+pI'zH(UqKL5$)##.7%:Y.6.->##5##bt<OG#+##V>##]W###96#25;T)?dG$*##B[8Z9,NR+###B6#}4H######5##V|~.,####4v#Ws9jH)###Z$%E#$.R)###)##wf0######$##k{~######'##jeV.,#.,#<##4Z%gm%1m(*##"},
{911.36,123.30,7.29,5.99,"RU6PZ&nR'imKdq3D}H;l$.Z#t{-[#Q###:?#&&),M4*~)O8+%##+%&*;:CH''@)4|4O<DI>#E9QE|AdG$1'$Si<k:Q{L1^s4###A%$`f4j>%)m(G6#X}>i7-:%Ou04/@(}v'4[&k:Qhy/>K4###;##b$)NJ.nP$D##RE>~I,g-*8l#cQ'4w%(c#~$*2d#)C6"},
{74.84,147.43,7.50,4.20,"~C8gL1/?&j-$,9+^n*$6&}b#0m'J,$[5#Wc%QL;+##UJ(,26?x-9p/]8.(l#HpRDv'>?&'##Mi@95####q$#&>M~$#^^:oT#Ax2Y##~:92Z##mR/##-n,=$#flR,##{k#Q'#?z<9(#`lRr8#WH(O7&,d)4R#sEG0##I#%>'#_lR=$#>%.9(#n?*9(#_lRo##"},
{970.48,150.48,8.89,1.79,"~[$qz6T7)_g9.M-eu&%##[GGf^)91;B$#l6Qcx)H[+U'$|5Q:A*CA.yo&:T2j%NZl%cZ'Gm%Q1,$IJ,i<pY#igNKK40R'o#&XI*#.)CQ:Z/0wGIw,&od(&7'qI+K6'-5:w~0mZGX.-0-$R]/Af/mx1O3.`f/Nw,z7*#:+Cn,o#%)e,aG#i6(<l$5v')##F&+"},
{970.48,150.48,8.89,4.70,"0##[n*u-)(I*MZ#=~-Ac#^R,y6'zc&ER)|d*bO9+m&Cd&l)<--%Kw)3BQN%-%nGh%-&]/=R'+n*,.&OFFV#%N&Gje+yn.Bv&6.&d6)wCQ,o1mf5F,#G{.Q}AQQ&oA(/6L;u#}+6Z7,e?'1v&@:#%?QnR&z,')%#y?QIM0ll')Q#%L66a5cl&>e(5Z%gI%*bC"},
{939.38,170.29,7.42,3.59,"1##H@&G=6X*Ed##=)<M:5OG#}0(O~/Q5$H5#oZ$Ru#~7-|b#g5$v:OrD=hH(Fs9A&K0U4oP#.r*}u'&f+~P#2^&ce+r#'###rZ'Mb2M0V{`:3/V,f.rl&:b1{R,e>#l(:p5$WR%Y94U$({k#/,#po0|]/1M9lY#Dj>e5%wy/(8-'E7@d)0v$0w+|,&dU*Np1"},
{855.49,175.03,8.19,5.62,"mC9OG####1,#^{f#########'#GM$*,##/u#*?#1//bS(4I+#:7yY$h?$[?(J{f95#$##Nu#:$R$R'lG$Z,$-##bB+NH'D>#hC2i,&yc#?:4uzfZP#)##'{)+tE;.*tb#'A)###|Z#Qe+SK5e'3OG#&##uVU6-S######RY*}S2^>$95#4'(TG#@7$1C9km*"},
{170.37,208.71,7.52,4.66,"###]AJx~2###Q5$:90Z[N)m&*[)bP#X7'*?:/x/V^3dI((K+%##kV9@`,;n.mh:7)>}[&mT3&JQF$'3V0kt77-'=['aF6kcB######LC#p3H=I+95#E0#NIQYGHe-'Y%+&[Fa.*wd(f%-a?&######G-#^J2######i,#?,M###$###h26L787,0,#(y09d%"},
{38.46,212.44,7.21,4.88,"#########udQQG#######GyTxY$######2{Tpb#######wzT)##*%)gG$'%R7]-'.,###8#EA$QeY#-##?R>je/}k#a,$m3:###`Q%5u#K{;ie/a5%Pu#$);OyTqG$-v&~S)ZT0E11QsBPm(sG$.,#jZ#w^8,l#K>#ly1T6)D$;|m*'%,1,#]8(9N3rV>Il%"},
{41.46,259.29,8.04,6.27,"sS~TG#95#$##ePKOJ,Zv*(##oy75]2aP#:@'fY#|l&Jc#:r=KS~###J>#%Q#.JZ*##@U/~?&HT~gY#V#$1R$E.,P//RG#~0/GS~###%##7H#GS~###DU$i'2fXFG,$`@#~039o/}92(n&A{=KS~######A##FS~###;$#Of.jc'###-%#+EA7w)Hv)H##$_:"},
{75.31,256.14,7.13,4.35,"####B#PJ/t5&###Ih(`f-Zw/$##]g0i@&mq?5Q#^R,S?#ux5-$(hw#:V8}02oy4^~&@q1S{>L?(###P{+P5NO5$2##JC,XK6Fx2$p#`7NeA,X6Npe%W(9]V/iA3@$#'?L1H$G$)y%$LNB3,#;C8ku$=+4sN=2-(ml#dE6wmKL5$,C#E5N@c$a?&y0*<n.2,#"},
{409.31,319.78,7.87,3.41,"#########bA-KQ&OG####D,:Nv(.,#'##S(QOG#/,#O6$+]N(##'^2###$V2^(17PG###v0,S(Q{H*)##0U*`I,,7)sZ#u~0$##Wh5###(u>bv)`J/###S(Qe$Qev*###h*0M-(]L5###BA,###0,####*fF###F>####V.DG>#9m'###W|9###=U2###>B/"},
{409.31,319.78,7.87,4.22,"T^:######f$#%IRYl&$##|C%@S+=-'d#&<C0WZ'l6#(K4rl#QkK-l#OG#E%#1gOFJRj>%8$#Cd=EFD}>&B,#$04QI%Z[,O$#EM?qZ#qQ)(%#64F:;+F'8)%#4JR/$'D>#2'#/D>Lc#ac'B%#6.-J-#306]?#3HRE$#lw0:'#3HR######='#yS3Y$#E82J$#"},
{457.64,334.30,6.84,1.34,"*Q#s)0-jEQ5$6+./)96&,VZ'XH$_d%iyX>u$5,#9,#ge0###.-'.%$*wX&##myX6R&%g78##pI+1A%&wX$##1,#<##SS1###c6*?$#(wX0##(wX?##_98A$#]?)%$#*wX2##.,#5##DI,###n?*Y##ywX`l$2YL,##yx42A%3l$S##4wXH5#95#,##=6(###"},
{457.64,334.30,6.84,4.65,"7A0######%##Z.U&Z#XI-###iL9^H$:%Y###-'YRG#0Z%###JS0###.,#$##t%Y/##g[,'##a1<5##y%Y&##/'Y###OQ'###qI.###&#####H%Y###e?'zu$]B7###hCQjd):&Y###e6':I&O-)###,##/,#o'YbG$C,#pu$iC2T7.r:&:q8`C:1R*~Q#s<6"},
{389.44,510.31,7.67,5.51,"W_0wHMCe*<S+sy.+r7Iw+F,$`80######J,#ub#F>#D#$WG#YS(hAS~l%+d&QCS;jBD>#Nd#A4C5)@###1##96#cA3P>#<5#Xu%0Q$n[%7BS/<Bn~+oS1%u5+9(GBSZ&4)##rC$2uO######&##`?(2S*#2>###v/$}YLu?)Q6$[y*DXEw-+bs2AR+>5#XQ&"},
{273.14,527.68,7.56,5.54,"GK*IPJ/,#7##KY4gf4###$##Bx,OG####'##/,#<5#<5#$l#PA+I[P#@)Y?&n:SVx3###9[#{FE.4H###,##.?#1D>.,#0,#:?'fP#(&'+9S^aG.K.Tx1Vu8*h-#9S<J0&##b45dh?.,#.##.,#+H%+A)b3F###(@#M7S][+QJ*UR&f=I>I)f>;bc'###VZ$"},
{172.32,551.00,7.44,2.59,"78-BS+T$*f#$AW:DkB1d)c##@UWo#'$##s]&|^:###$##BK.ZP#E,#`i43/0_nUMu#,<2S&,UWWpb#%##F(&J;=#H%.,#@H####$##.p&HRW_6)###Cf$@WWO<B###TG#&n=s?)+I'9J./R$.,####K##nNE+u####-##7mNOG####B5#EK.###&##]['e6)"},
{172.32,551.00,7.44,3.45,"px(4<=DZ&fG#mR>.g7###O,#TK(3.,~P#&##_,$J&-L$($l#].&<_/Ih9?u$z.A?QQ],$%d'?|V{h@H>#_5#xe0P>#KR&+.+fc&$c#xJ)Hg5Vp2@H'}5#/{VwwV>u$###tD'hA3###$##B?$,6'###V$$rs<uG%###4##dl7Sn/######C;%,H%95####7##"},
{761.35,559.95,7.58,0.20,"*##_=ZO#%7l$u9)xwP+u####-CNWd+###6##](=/,#C5#c?$E@#t`@######H)7Rr?###*##/;Zj>%###/$##=DyG%#l#M6#UQ#Xx2`#%gY#D8Zaf3###1?#q9ZDR+###e6#uF2s84vb#$c#-u#W>#.?&Iu$QdU.,####o5$`8Z~>$.,#/d${x-X,%N5#/K."},
{267.01,606.59,7.27,5.66,"'Q$-/(uw)/y3ix/M-'o6%n*D`#&###pH#pHQOG#.##A/-9PEz6(_N-OHQf5$hKQnJ->c%mH$b'76O?],$x[,+##vr=}^7%Q%v5&En#G}D}QK/GKIv%tn,d>8,7)R:Q@X>2H%s[&{ZPkG$EQ%cG$qY#ul$qGMfY#c~'My3Ry15Z%9A%|GE8A,(6&OG#Su#_V9"},
{688.04,640.46,7.36,1.73,"MW;mP$######g|X#########<E8j>%######z[*OG####:##:U7M5$@##$`9_xX9#$###x%%|YDR&0.,#(##$k>n>%###4##Fl%0##k55exX9wPaQ(D##^i4jq4k(;###TZ$sAPN#%###0#####3R#+0IQy8mG$95#(-#}wX|Z(pb#y5#UYG]D=###6##X~("},
{461.68,654.69,9.52,6.21,"L%-Ku$%(/6g+RI*Y$&)p,st@hB6|P#B?&rAEg~1Q5$###9x%|]6e~)]~,U.'Zr39{5{;@tJ+ANU#e,/c${T*$[)%##$l#hn){09AH%L@*>U*$f03d$8~PIP9{IUI,$#6%oi*hd,2Z#Z>$lG#T,%&##OZ$DnMD>#(##KH%nMUxu'<5#3,#h).fY#xY#mG$9H&"},
{296.36,710.72,7.84,1.80,",e(Fl$HA)QE=yc<D>#[,$du%ps4;w,^l&###}~'+05'Q%###g5#QK-a3@dl&fa:G8.?$(d>#[^PnY#t-)FQ$wL7nY#?n)lm(7c$$J(V%Gy8204Fw%-@R(5O0YZP95#D,#$r(ie0###)?#tr:$H$uc#E2)?FG2S,D5=ow'+HAts1*3C###8R#O(*0Z%###?u#"},
{182.03,730.08,7.14,1.72,"P[([c%8?%sn)]GBE,$I,$$6${yR.,####&##ue.###:5####Fx'KS+2[)2,#VF9tJ1###I@(czR6l$###P6%x(995####)##Xf1WR'[s01`:uZP.?&4##y(R.vR.,#,##T~CV82######>@&%Q%$##GU#2wR4&24u#b>#PqP9Z%D>#&##Y;2pb#######nI("},
{688.83,800.69,7.82,1.33,"K_7H#$W>$$##}Y#z]0Yl&###?l#^7)2Z%G?'###Vu$}b#Y~0+5L2##7S/9v#Ff-gL1'~W#c#R^Xh[+Dl$T8CHR*YR(2/'8^X?d*$##VJ0Mx,Se.Tu#9^X^6'UaXgjGyu%;T)(C-x`Xy$(@f0###$##(J,4[)fP#cG#5$BN-)k7(77,9m$;S,9@*Hy6(l#zY#"},
{453.29,821.36,7.39,4.35,"`R(wn-F-%@d)eP7?n-3l#o>#~)8?u$###C##f>$v$,###3##Pc%Y7&D:Q)r?D7QH6%pe/pn#ID=,K+oZ(z-#o,#=7Q{k####v,&uV2a7QT6(Z7Qy,%~K2W]+]K2%;(~sGyP#}x*%^KN{A3,#;5#nG$7%(hT4cQ(###KH$e4@Cu$Q,#EeOG?&IS/??$sAN.:3"},
{692.75,826.68,7.40,1.38,"mG#)].j>%$##U,$oZ%oP$Ne-$##V5$B5#&`?######,n'KB5e16Gf.o1<-c#0yW@v(cG#-;Okc&`B.)B)mxW###GB,q'.8cOA047l#TFAtZ&~|Wy[VfP#i[&V7'~|W@d'e-):5#jyNFQ&]Q'&c#gP#%;3F$)LI'a[,Bc#~[)zm+pS1'l#yY#)d(S#%###Md("},
{692.75,826.68,7.40,2.82,"Xn/###/,#b6&3RV###77(#1,bQ(###>DXId(/S-c5%|H)%##+C9$##Nu$a>#oAX'v&A&-t^(87+S.''EXE@(,H$8w,495~P#:w+`Z%r-(wb#H1LMBX`H(/Z#Fw&1FX0%WZG#e?'Jd&17,[#$###Sl#,'1o#'DQ$D(7cw)aNCiH&ex1Ru#qi?Vv)G#$P>#0U-"},
{692.75,826.68,7.40,5.99,"Q>#M^.Im)=u#Bl#lN@gH&)92~e*7<B(?$8q6af0-H&###Dc#=@,Ru#9d'Xm&zmVZG#Aw&h*X`H(|P#7qL0'X&I'1u#e@,-6%/B4hY#b,$C%--*X27(gm*P.'mA-dp(Q&Xwl&Bl$~>#mg:$##P6(%##+~,Kl%/)Xc6(C6(###E%)CU,e-V###.,#7[&)&1###"},
{704.29,859.20,8.20,5.58,"0,#^G#:d%q?*&?QiY#&?$);)]=IhA'AR?}gO$6&O-#wF=k&2z5&6Q#Q~/-Q%%@Q;v%tZ(SJ#7#J[<.*YK06#g5%_J#.::)##L6'2y0b$+T,$NeHL*AmZ'r,%*%'tCQvn1&##`G#10/(m(0,#G(3VS.=,#RI'ym'+(8,_,LT2e##`sEp%(q~28e'3B3pc&:#$"},
{466.29,862.40,7.89,1.28,"Yu#qU4i,&###_d'd]2gY#Kc%.,#^c#wP$E~/###0##7l#>]3D~,;S)C7Q7Z$LUZqv*N6(J`5H[)&p5[>#)UZ'd%EK4k6%iUZ|H),I#6/Vcu$&iV{WDj-)GI&^w'yUZV-)w-*~^*}SZwP$0$&4,#B,#{;03d)^l$/?&,J'V]1n,%B91kO;.~,J{<'S.,H$[f)"},
{466.29,862.40,7.89,6.02,"E,#Ti02d)3,#zH#H&Wau$Ad)'/*1RN3Z$I~,P:4<c%###Ll#/?&:e&502ju$Q3D'@)3@&`2X+.+qc'Ii5v'ZS]1gY#Ul%{6':T/kX:Lw,>c$g(Z3[)5R)Fn'?p6T>#d'Z<R)DQ#1c$<A0.,#P7-Uc$a]*5r<M&ZtP$$-%/C+Vp7}H%H(Z?v%-##9l#m]5###"},
{223.34,891.52,7.94,1.39,"@6$0~)W>$###,v$/*9_5%$##X?&Zm'P>#o?).,#M>#D5#Ye/vJ0oH%3&-cQ&0h.G'017N&Q$_iYxl'Uu$f%CiI,}$)JK'^gY)K4>Q#Tf2wv%'B3Mm#4gYrR'#kY%aE6H%<p';z.*iY7w(V&0xY$E,#an*7R'hc'L##/N1kB69J,Yd+4?$Pe,)d'x'9_#%K,$"},
{226.77,919.54,7.88,1.43,"6$%So.,u#RG#>c$f,$C5#Sw-###&##G5#N'7$##[G#}5$hd,D*7p7,F:71c#y0W'v'/Z#JWTw#&>]-Yg)=/W###nI('{/,5N4V;Gc#TOC+.&W3WQ.W|b#iw&c[&K3WG7)jm*F>#[BNq[+D6(bZ'3##hS)Lx1s6(OI,9Q$Ed(Um)y~12Z$|P$1v'0?&0,#7v'"},
{226.77,919.54,7.88,2.99,"<@,$##ZG#|5$,5N###R7(g_/O?(;5#U0O/n+X6(tc'HQ&1,#|]6###&##H5#;/Ww#&:]-fg)%%+d[&<3WT@)qG$k$*>x1&Q$Rw-Kl$q5$D5#4ETw0W'v'1Z#kw&W3WP.W|b#-R(~$(OI,EZ$G>#D-%$9/,u#3c#`E7Ow+D(8-.&2V;Ic##+CLx1bZ'2##>/)"},
{226.77,919.54,7.88,6.05,"+##1J&8x036'0?#Hr=Ew&)h:h[)@N?^u#4|9zT2+c$<5#t6&ym,16%<-'sZ&t[W,l#:/'8aWkc'Lc#yMRB^W=H%@,#~7-WH&'05Z5$k>$f[+M_W6~)X@,dR&b~,SD)k~Wy#&>5#R5#2T3###Id)@5#kv*r,&JLS*x-]?)/,#Z6&ZV-C+K###A5#Su#:m)###"},
{965.16,174.81,8.53,1.27,"B_:~>$6v%en's>G8$(:,#H7(p*5s5&*`#Hz:p5%###2W$&'6V]5qP#xQ(5B&~f,r)5]V>pc#7#]A/076%)m$faF26$4KBUR*:J-v>%,R(yl$II*36%/5?4o-Gz~#.+qu$jg*B):TF7}x59,#(6%xG%=5#(e*y>$Su%d>#)U6_8-7?')##1/*,6%kR+F>#uZ%"},
{561.80,253.42,10.10,0.33,"%##3,#A,$/,####B##gQ(~.+###8##.c$<#BT>#p#&.,#(c>{,'(##R-%xm+zGP&##bc%MQ:S6)g5#sE;3:RyY#IJ)R6'pg6J$)D,#-;RNh<-;RNe.[~.@-$xy2LxF7+@wZ'###7?:)$&1@+N>#Td%98RY>$50)XB6sy7$##+$$a%MZP#%#####nJ*[l&###"},
{825.99,263.41,8.51,0.45,"-['`U3QG#zu%tW>XD8<5#9l##TW(c$###t$#)f/OG####N##(-&M]1)Q#?,=-,IG$'XJ#gSWqTW###5##^w(AT495####I##GRW###M##V-:<n.###h(#cTW>=I###/$#(hO5I+95####f-'KRW###G##E21@x2SG#PQ#z@-###%##9Q#g&3######'##T[*"},
{825.99,263.41,8.51,2.30,"*##6d(pb####K5#He-eY####ZI'nH)######(['######(##5##h/T95####,%([2T###&##?2TbNC###L###K2G>#3,#)Z#(##C^Q:]/95#'-MN/IA$&Jl#&1T('5###q$##T/MIP:5#D,#JQ$t#&L&)>u$y&/;Q%q1/g5%>G4e-T(Z$L>#c]/D/T###.##"},
{733.25,504.47,8.18,1.68,"~,#n*;T,%###?v$O/.TT3QG####/V+5N@95#1##cbD)$(###;f+tW7PFDzb#*1Udd+r&,VR(:$)[,#x2U~Z&####m#$3A###e?)3,#:,5d:9..U9u#(~)+*4K81%-$52UtY####z##=fS@u$######P-#ro5s$+X,$16%&U3V#%kH$`%Ob5%###%##?F5tq>"},
{733.25,504.47,8.18,4.47,"ih7v@/####-$z&X6#$;5#hn#;B2RH&1-'(n%A-#)I*D>#####`A######Z~%A%XFl#sl'c~#KB3`b2U#PtZ$3|,r(9_5%###G'8###/##mC22'X-R&:[*B.$SA-YM)c&X=-'Q'XAv&PR*be'bG$###v##>5G@]3<5#vG$5L(Ml%@5#r(/^2>DI,###Jc#V#>"},
{785.53,563.14,9.19,0.14,"_[$H)X'H%ic'BH=2C8###>##%{;+l#'##5n$-$L.,#:d#g7+t'0OaE###'##FVX-H&####$#@QNR#%H>#oQ#AnH:h9-H%<H%N,MJS0/,#w>#QUXVZ'###Y6#9I;Lv).c$kG#N6%fx-=aF$##+4I######Iv${RXOG####Bf%(h1Tc&-6$:f,Dv$2o.Qy6nY#"},
{847.39,580.00,9.41,4.80,"~v)qm&)$'M6%al%Sm*N,$H(0`l%6l$<##SS,95#5,#eQ$G]3o[*+S-=@'LQ&(K,A6R9c#:v'}SSu./`5$?5:+$'Po'[]/=SSM$#2;?{6$HJ10S+;^80[#JFFpUS_FG2u#Af+;R&8&EM[*?Q%T&#QQS$##95#('#QQS&##*c$bK#RQS###^P#u.#VQS###$##"},
{165.86,599.12,11.27,2.72,"^7T&g2Am'cT)Q@*NR@ZD8,K1W8T|m&Tn-8e&e925d)SG#VQ$):0Z@P~u$`-%KN<R*957(Mk?OaCpu$bl$Iw@td*fZ'zc$^o/Z`;()8n~,Y?'G:TAv(sY#@/'GE=qZ(0,#2w&P$%I#%-##jZ'j?)$##[/&o7T2@Q%##`#$Ev;ZA.qm*Xd(pd'ge,`5%$##s#$"},
{165.86,599.12,11.27,3.66,"pM4?.-d$&to/p<0Xr@fu%J,$>x.Cm'7c$j%){c'gu&###U,#(w(~x)N~OEz4P]OXT3:D:dB)gGGaI,.,#76$k,$r]2eY#(##pq/A20+W@@?&dO3svL.o13I%.a/HGM0u#9,#tA-Fn*?v&;I)5['7a0b]4g6)'J%dM6GU5a%-a^OfL<`G#jz-|XFpb#(##W-%"},
{362.36,604.61,10.03,0.69,">`@ol$sU7i1.icOn..Ee.N0'.l#3j10dOFm)&f.br9@m'q6)&q9>-'qZ'9K)bU0BeO*16tc$HY=(XB}s<,V:f%/(c#U0)qcO`w,X$'qJ0]d(>7*b-(|r1vEACdO06'_-$f>9m@.q6&fdOJK1w]4f?%r]3W/0m94.,#L##aK0KJ.Cz7BT1yo,^],z|D^:95d("},
{362.36,604.61,10.03,5.59,"Dm%_3B<I*oR+xl&f.,Q@)Po1g$)_?'_T)7OE#m',J,/x'(nOrP$B,#*]-MV8;n*-`18XCKR)~oOIC6bu$qy2vI+Y)?q@*AGK4%*}7,=y-b~.G$)IH##H=jmO'`>X(39}FyGAY{-).OJS0`5$9ZHsp4]l%NA+/o/-[K0e*~07S>#tp1,mO[I,=%*%H$a/0y94"},
{385.64,612.14,8.83,2.57,"_d&w]2rG$%Q$<B2DJ*Dl$EJ,F~.Z,$lu$%BO,%+|Q%`2<LM5*7)a90h11#v&@@O=J-6l##@%3}A}?OE>#~A,^m$nPL@(6Pv(8$(/l#2i*DbGm_?#R%Yi:HBOgy4%^LEtFg6&HdI4C9/,#:?#_^6Xx23##CK1^>#g3Ca+Cef3pZ(%['{HE?(45&0Lv(^u%3~("},
{172.84,628.08,9.02,3.58,"CB3P5$_c$XT,?&F+%-(##Wm#]%%.^3$Z$fY#|H'QJ-'Z$$m%A-'Y5$%p+.JLfD)x`A?HE5.,LhQ-z:(Z$8-%DS/vb#oY#{n*i12Q(-t:;Q~(g.-?T&?eQe]-^eQ;?'_5$d((EK0GI+###2##`F7X]-3J/aZ&im&:36`B7V#$==1H<CeY#&H#y'')o1hY####"},
{465.68,627.89,8.53,0.90,"g&'(-K{OD%Z$3gH[o4;#$$##&M7:#$95#/##X-*###Z5#Y,$(.*8f)+:WTR);:W'?&E-&HT(p8P>U9###/##%K*wu'I>#E5#RU40&-:[(&8+;XC}~.xQ'0'*kr,E9W4J+<Z%`d<}e2%##yb#w81SG#?H$J<@t&-R,$cd'GZ%u[)Xu$Q^('93j8WD>#A##0A)"},
{121.13,643.24,9.69,1.74,"r5%<v't|-@0W?M:7f1%B&W3WYR,nY#;~#Q_OZl&8,#7{%L0W0@+###2~%04AO[(6d&.}<oRO3/SXl%Qu$*o&(.,V5#LJ,A&0so5###d5$hc#al&A5#g)1yp6*2>;5#-6$j44<m)###Dc#DA)6[*######RZ#16'###/##c6'Xm)eY#%##Cv%4.*Gl%###2##"},
{121.13,643.24,9.69,4.34,":A(`7-ZP#)##&##&T-uc(###2##rR*x#'D>#$##Y#$gH&CZ&g/.:%)Fl%O?$|#%_@?6B5@5#]a8Ur={Y$P5#o,%}$&OG####)#3zFD-H&,##(.'&8Jjr@Iv&R$I5w+]l%:A&u,&on%Y-*&##|l<fAXpP#k@+h)+1FXA/*E@+{]G3=FcP#Ku#ml#fX@l6*.,#"},
{352.59,672.83,8.98,3.70,"###=0*6.-###]G#V1%Oh>###56&#^%j6(G6'n-).R)Rd&Wy2;v#-#BcJ2###{m(Lr*2nP|l'VnPD7)U$'kV3tv)Iv'tg3__8J(-8HM1w,F5#v|Bd%.v47OvO>}F2S)p8+bnP?6&/p0NK2{b#(z8ev)VA.[7&XmPS%*Bd)18'r#$UX1qXH}u'be$L;:v5&/,#"},
{352.59,672.83,8.98,5.56,"Dg-(L5_Z'Wc$:%+j?'n@*sz3wP$6H$LM4kI,Cc%X$$o*=}WB<00T?'$6$,7'Dh1T'PSI*b6'D&P7i;~v%m&0)d((B1ix)t#PI>#g5$vy$@#P.d(N~)8(,9%P<vO1q3|x2i+:d$'8&Pu~0uG$######j$#eNE###)##.z,PbL###-7#uQI9A/]5$yl%n93/16"},
{243.22,690.42,8.59,5.55,"#&+'C7<#$>6&>r02e,gY#Ic#kL-G)9LR)Cv(*c$[J%vW<ebI4$(m>$`d(q82qK0l+<;o1.-%#:R%o/'Z#MJ+m[,wx1^[&d4I###+##HD*{<H{?*/6$i`5G8Ro6R3D5lf4rF76o)n8R[d+@5#######y7#J3F###1##+y+66R###OZ#:8R>e-G@)B6&GL6Qm&"},
{118.45,716.27,9.30,2.92,"6iNf09###.##z+C_5%###1##Zg3A,$###-##0d'-u#TG#$##)H^0$(TG#W.$d6F4r;1@)PH%.x@#/0TG#06$c[((m(.,#(##oC^8#$nP#364C05/?%>T)44;YaD}Y$`G#z{.9R)Ku${k#<##@HDc[,###=8%QT(]T6*##*?%?9.-H&###Cl#m>#K7+Fl%###"},
{460.37,752.10,8.74,4.61,"g>#MILA,$###]I${bJ######Fc$f$)2Z%&#####>,#R<B###j.'=zLHSWmY#,XW_tHac&Ic#_i?OG#k,%m,$###$##3jA.,#dQ'w,$|NW,s@jSWB#$x7']KG7sESG#5-'2R####%##.g2#########Y/%cx33d'U,%4Z#AN6}l%H6(xY$G##%##6u#qH)###"},
{211.95,868.94,8.54,4.44,"M-$%%+@6&(m(Li*CK4@u$'##.f(vd*###&##3##Yo1######R,$+8&dzQB3EtwQn.)vo5W-##);J](I..#$#3R#~IIac'###Nl$b6=VxQ)v'LwQBQ$g(6#T+#p1qe&ytJ1?$lo-o*52wQTZ%2,#p$*b[(WI*|u'###3m$#5@zG%$##;z.LJ,|,'&##%55o|A"},
{233.73,953.98,8.34,3.43,"###T##0H&M_3,u#6$#km+a%-.//ZG#_,%_S)%?%vY#p5$_6'U,%.H#3@+'?$d'9;##V-&{SGym*~P#H)(]5LIZ$cR)Bd&PG#F.,Q5$wb#Ku#KyS###Q,#RA)46LO,$5?9sM;0u#x#%&KFS5$-$'###'##/H%[~FcG$%##>5#s/GRzS827L6'sY#tc:o-MVG#"},
{233.73,953.98,8.34,5.61,".,#fG#K$&s#'JXH4u#)Z#yq'1GKrI&$k4+{R]P#4@#?8D(W@zu&2m%:d(t,&.JSCm'?l$80$G@QTl8SB69R#:5#1h#>iA###=v%oB3Bu$/H%{[@th<OG#2H%W@'IyJz,'R>#5,#x[%du&###5S$-%-@5#5.*:$$>6()##*K2L#$sl'5,#VI*1,#rP#66%Y>$"},
{803.58,145.29,10.85,5.38,"qD>w5&###m##epZmP$###S,#'pZeY####2V%.?5Yl&9,#[6#[q>###%##cQ#woZ######H##1qZ0,#E>#v##c8.sK1=R*4##q.0###1##6B*WnZ###*##`v#`nZ/,#UG#$/#Dw,,.(%a6Pu$/K-:5#?Q$DI(~RX###e>#`T&WnZV5$###]%#d#%%8*J?&QG#"},
{964.80,147.74,10.32,1.74,";6$Z*=2o(9C:s)55m(1Z#3ATKL,XA2W$#2ATz7(h,&,M$e?T5m&X7+hj/Ee-KePBH&~m'>.'5+4;YDvB4l,%1^J|A1]I)x?)&$&MQ&W=1xU51)>F#$P_,a7+/8.'6&t0.q.-J==K[*fG#7n(wb#Kp,OU39.*Vl%N0*u]/(d'z6*y5&=5#J.&su&/u#$##.x("},
{256.88,153.94,10.44,1.42,";wY3$#9wY?##,##c.)()@###Ul#+I*###$##vG$95####'##:wYP$#:wYD##|Q&BW0[wY@5#*7I#f/}G%g>#M-)###:5#A,#9wYt##LwY?##IH'/Q#:{Y..*Y[SB#$u?&K/*'%+&l#D>#(##OsG2$#;wY7##/,#%~$0_9Z>$2T-ZH%CH&nG$sm(4[*######"},
{256.88,153.94,10.44,4.56,"######rm()R*7?&oG$To-YH%^:9Z>$/,#'~$0nY9##&aF9$#D>#(##&%+&l#t?&B&*nmSC#$.rY/.*IH'0Q#AnYB##.nY{##:5#A,#M-)###}G%h>#5@I/o/PnY@5#|Q&FW0/nYI##/nYW$####'##vG$95####$##Vl#6R*j_@###,##g.).nYD##0nY9$#"},
{25.00,169.31,10.51,6.14,".,##########r.Z95####3r&.tDp':###S4Z@##(+I8##n.Z0Z%######%##}.Z95####I$#5~T%h7Su%A##.##FgSs#'2,#Uv*######*##z/Zpb####k##nYBJ~,s.,,l#J>#o]*=r7Xu$m[-######/##x/ZCZ&###M##xB+kq=(l#2l#C5#6p.J[($v'"},
{796.79,264.59,10.94,0.84,"Fe(<%,D>#G,#NJ)'q.Hl%mZ&2)Sll&.,#bQ&d1695####)##F7(PI+nP#WZ%OtIrR)F,$|@B~,H<5####2)S1V:######)^+wQ)T5#kg/bB2U&S###7,#Cs0j[Q###3,#>k@D>####6,#<o.vx+};5R/1,d%2)S,6'###D|.uh>###'##{C.######-##-x0"},
{796.79,264.59,10.94,2.08,"H##%wR.,####4p,IvR###(##Qz6%Q%###;##rd,######6##'d'NxRWA.5,#SyR$GIC#$F.#+eOOh>###g6#Pf0wu'###g,#qZ&yR)$z42,#LeBt_?h$*3Z#kR@quR###1##IH:NvR###+##sY#58*Yl&###td*x,&hG$<J)B27`80qm&}/1-9)iyR#S-###"},
{562.49,251.92,8.92,0.52,"###(##Sc&UG####,##96'v'/###%##5u#l6=a5$#-&[>$[L+Wv*###U,$K?%W5O###[>#PP5{u''##3i.+nE1,#>v#+8+in//~-###=V)OC:&{R+m(dZ%pH$a:1)H:_`9'$&###g;)Jn,j>%VG#4##&{RW..zU)BT4zC7###tv#xyRZP#######BC*%Q%###"},
{353.87,275.40,11.06,1.43,"9f3######'##-:f###3l$M##zMB###z$Qcu%mP$###PRB8aD@^9###0,#Q>#B:fTG#>u$R##$cN(?$py8:6&`v*yH#S`:NL9xy:######L5#A;f96'95#@##AlI3]/@H'###I~/'I#.cP###T08######,##T:fW5$.,#,##4?OCZ%j>%(##BS0A##L#P$##"},
{725.88,293.58,11.19,0.89,"l%/Uc$$o+Yq3nw.~>#]w*WS*}7S######dK)@.,###%##27)P&.{=<Yx.UQ%|N3lE:8d)*m$J;S}>&###hx&~7-###(##D~+%##hI(e'L{Y$IEA.8+hm&Sz&l8Spb####P'#e16###%##bP####i##J;S1%-G[+$##hA%5]FGz;TG#?#$19$TK5.,####,-#"},
{225.77,332.07,11.94,0.03,"M$*###ZG#Z7&YjH######@405[).,#2,#sCXVG#OG#jP#G&XZ&4######zu$@AXj>%###?q'_A0Dg7%##reKE5#k^:B5#^5Ktw.;5####0,#1FX2M>###t>#i&(MBX###tu%3%#:AX###E,$2l#/,#/,#.,#F.%Z83######S##QlN.,####J##DtI######"},
{225.77,332.07,11.94,4.84,"/,#/,#/,#tY#F>####<5#}..######FH$%1:###(###e$5J0=/2######Q7%tq?###Z5#fXY<c%###*g'%TY######mk2NkL'#J.,####F##DVY###0Q$_8(@M;$##[@E)o/95#@5#fVXJ-(lM>######+##eTY###/,#K$#t3D&##Ia@=,#pb#~G#R-KcP#"},
{247.69,618.81,10.08,2.85,"^g7nm&gu%M6&#U6X]0Tl$Yh0o,%dT1bX6EA.*c$0,#-/&|g7mF>VS+(Z$h.)B)8U->yg7nl$UnO3j=*v&HA)X@,Q&+?U.O2<1f/Q?'>Z$SZ;4I*k8,MAFKp/^mO)q6&n'Gf))w){l=N$)xl&]U1d4>^l&uG#n6%QpOg-)pY#I7&GmOmG$qY#,X:oB5RG#6?%"},
{695.59,735.66,10.61,1.11,"tL<%##95#fQ#AI*f~.<#$fG#Q@(H%.A5#W>#'-'$##xG$c5$-V=:##In.ww$}['qnCP[STl$z_SiODhZ'al#Yx0ow,*c$(##}I/&##O&/^e&&m'D6$z_S&=B^[Soc$6^0t%Cw(;#O7yK8+d#(c$###R,#SJ*######v~#zsH+u#0##^FB0bE4bE@H$W[S6[$"},
{695.59,735.66,10.61,4.27,"f[S%[$8=C{5$%4@f}F+u#,##o~#A4I.,####O,#C8*>u$###~^:)d#92<kN6Op0PSC`[Sjc$z_S.FB&m'A6$Bx.=J&an0&##5l$(##.T0SJ-sc'_l#z_SgODR[STl$.e'MADJn.mw$e(>:##)Q$a5$'-'$##@5#W>#gR(H%.<#$dG#|-*|n.95#~Q#U1<%##"},
{463.87,777.41,11.04,1.09,"T7-)##?5#8c#mY#~~*l6&W5$`/.W.-6l#L>#8R)6,#xG%###@x2^Z#Lm*D##e['k.B8RSm>$nVSguHC-(H6$}8-Yh:W>$'##xu'V5#X~-Nm(F?'AR$nVS|$QTRS'd%4;7.CLUi=I)704H2R$L5$###3##{/0###=5#Zd#HPL###]##uPG#i;R17Lv%(INJ6&"},
{463.87,777.41,11.04,4.54,":I?#@*%c#/,#(~)t}ED>#'##KH#EtD95####B5#[d&U%/###?]SZ6%GM:4v%R^/~VRr[S6Z$z_Sm5LA6'jQ$H~.:#$RJ0Y>#e..H5#bx-Z+>Wm)YQ$+MP2#Gc[SpG$S.'t%C6e.J>#Dw,v5#?$)###{P#&d&;5#/,#+T()n,rZ%`#&:H%,8+1Z#tu&<c%(##"},
{445.59,892.67,11.20,6.07,"###c:4eY####[>$^Y2ZA2#Z#?@+B?J*m&+?EbL5=.,ul$ih9###iw(K/1A,${f1g3.r}IEH%=zX3@)gK00P:Rw-U5$lKTW]2###hH#,x/{u's@.Q?$J#>7ZK?zX@A-Y8.4X6W&0N^,bxXM#$L>#@-(pY#de.OG#EI++Q#&V<PA+YK4}5$v$+tG#yy2oQ(9Q&"},
{668.31,925.12,10.08,4.42,"9K5######:##g9d######m$#}J[3l$###_##bQ$^#&.,####}_A######H##t9d######f$#>/Y3l$###S-$3d(1c$eu%%Q#w;C######^$#{:d######6A#p%J(c$$##?v$`?&MR*6H&0R'/XG######s&#GB`<5#0,#As$3D9ml$U6'sf+(h5]Z'<5#Ai0"},
{189.33,1020.95,10.96,4.55,"_r9#########_uk###$##Y,#pXI*6%Sd(Hv$13D(l#VG#Cv$y{B#########^rk######yQ#<,NaP#xP$J1/<81L,$iG$Jw%[)B#########rrkD>####O5#:fU#?&:,#7v%qR)VH'ru&O.)Kp8#########Frk###$##gG#6-RD>#oc#k-'C<=VZ';,#O_-"},
{996.87,113.94,11.80,4.05,"nL/3I%maET5$Z-)pv%t&1{?'L'6Yv&0u#`?&g%.{k#1##'{8-P1mA0C.,95#Zm&t^,i::2Q%mmSmm(*Z$b6$m}Cr#'wb#2R%9H74H&D>####ffFYc%08)$]0vW@[,%{6$[91[6SAc$O5$5h#ez.#########,2P###J,#hG$m4<.,#p,#sq1^C<wS#k-+OD%"},
{701.46,205.76,14.60,3.22,"d>$),;E~/}b#hQ(Ex*[q4.D:)oY4H&3Z#$4,rD?/R)f,$pb5H~+tS-Mc%;l#$V5i4@`.*CB1.sY}e.TH&gm%jX:EoY<u#Dl#>,#,H%jn*{k#P^8Uo/d@+$$$tnY3Q%/,#U.#OC8h6*###tl#######%]#nw0OR,###$@$DO8H,P.,#%##5E'G%-zY$4,#B7'"},
{492.06,259.57,12.41,1.67,"//-OG#.##C#$fV/Si@yP$^$(4##}{4qU9+~.###=##c?(qm,=[*.,#%##e5#i4JEH'?##78C-6'ac%Kp(PnU###dZ%I'+YHT*A/###$##O>#BpU###I5#yv&esD8##7^Q]A0/c#YJ+:oU0Z%QZ&#########ipU#########FoU###mOA$##~>$/,#7zScH("},
{525.62,274.97,13.86,3.34,"06'######Y>#T/3.,#$##Xc#'r@*##NR(9j6)l#t,%G[(_ZGk>%###@,#(ZA:{@###=K(p*D)7T4Z%:-%|{46$(HC6V#$V;T1y0iI,6Z$hW6GK1{^-:9TIK2g;T,7Skv*PZ#*[&C;T###3H$$$%PuI<5#(##Kc$H+=SH(-##NQ#W81A,$%##L>#2o0###$##"},
{752.71,309.94,13.89,0.81,"P6(G,$yv&O]/PmT[>$_G#@`.z:=######`?BHl%######/e(hf.%+9vf1=w(opTPH'@#$DM%5NB###(##j(1######-##7x0W8/KU1Vv@D'+?pTWc&$##G_$Z6O###,##'R'PG####=##%[)>u$/##vP6+%KE`B###?##{*/m(=95#%##%I#4c$$l#%##VG#"},
{382.25,331.20,13.88,4.53,">n-|k####$##BJ~######w##N4K'##DK4<Q####4##nIUL5$so5PG#%##]G#,K~xY$###-%#WWAcH(sI,cR%N5$T,#KK~PZ&)f2###kY#UG#rN~mK7.,#7##oLLT)A9Q&/##tI-??#<J~%##qm,###|k#-##N*D#-&7Q%5H$=K~m5%ju&9c#=x2o##;J~###"},
{738.45,476.86,11.64,1.86,"4.,y>$95#tG#BZ%%M+P974,#,##pU.ltK###$##c]K[82###lG$Lv#EA1<,#'7'96;C-OJ,$FQGSU4LHA~v'mP$;d$@_QiG$*##'H$mJ0eY#%$'Xu#%]FR)?eZQyb#R(10vALQ'%##@_QCA0OG#/,#~G#cZ'95####6?#7^7A,$###FQ$@N?95#(##HB)BOF"},
{738.45,476.86,11.64,4.47,"BLRE;@###8$#|JR#J,{k#I$#Ow#8T395####I##46';#$###`HR5u#;5#(U&cIR=56tU;|Z#y]J_;:bl&GZ#0o,Qu%###L5#Z2A?5#/H#x%DXB5Od$q*:K-GO,MM#$yQ&,n@?@,###$c#m-%p#'###>$#Or7)$(P>#v>$|G=95#k,$8Q%]z6###1Z$}Y$wY#"},
{479.78,510.25,13.62,0.51,"e5$=c$'H%[G#aQ%WA0xb#3,#v>#rd,Hl$###ol%rb#/,#6,#^P#3_*_$)^P#gw'0*V][,6,#CzMY=F+u#>##'T1RG#$##`5#I>#&?7K;=i5%^81g6;vI-4`+Z9Vq?(dG$sh%7(7[5$95#F##Nf1Z~+[C0T(7{l'2n)5A*U%@_=9W_4+y5r~'%l=*?&.,#J,#"},
{479.78,510.25,13.62,4.52,"C<<W6&5B37S%y.,na7@mRw5$/;XQ%-'Z$Vm#y96uP$@H&X>#&D4d#%zV4D[LSp516$#=XO~,N9XcP#F6&]9(.z8RQ'O>#s#%-u####2&&X;?{5&###h+1&nTe+JD>#lQ#q2St-+qb#fG#M%)@?'QG#Vu#X5$95#0,#Q[%m6*@#$cG$<,#gR*D,#_#&G>#4,#"},
{366.50,528.99,12.48,4.12,"Id'e//xw.P8.U$(K%)?p6L-'']0>(6od+.:*duPoP$0,#V^'#|@>6'T]1UN9.&/jn&RnQq%-2nQ+v%`6(<j.U06_P#]G#vy+E:8h3=qn0=w&r'3fjB:RIG^2`mQiY#ml$m#:$w*O,$/l#>w).6&[D-qWE#?$HnQ'^4$-'hK$hlQ>5#E>#]o&^P#kQ$vm+8?&"},
{303.76,539.51,12.44,4.69,"4%*#@)w^0'0+y@+G-&J8EPf/MSRE,$kG#gz(V97:5#.,#kQ#.M(Q15md(ZQ%4e(p}4hE7R.+QVR(%*5l#'B'V:4#l#eY#C##nm&DS+V.*Nu$//.<K)Xr9.R&HTR/6&[5$ap%7U2Vu%_P#lZ$o~.9R(fd'#^*wm(ml%KwA`%+WzP1v'G7(}e)dL4Bc%iY#9d$"},
{248.82,545.99,13.46,0.95,"95#q>$Z5$lY#.,#b&';<D0,#Lu$&|,q3EC.,8B4+[&bQ'W~*%##ef(v./###fZ%T.B]aHJ>#P[M)_4L03Y3>1.,sY#6)7f)=###;9#16Q###QI+0;'V6Q;.(97QCw+Dm'Mp(*g2x^6l^4i.-###0-#s8Q&@+J@,Fv$Ry.17QhND-Z$+v$%r22T1Kd(j%'LW?"},
{248.82,545.99,13.46,4.66,"Bm$/p/sg1eQ(Wn+4IA%^3q>$uqRE%-0,#b>#e'07?'0,####7d'Fi:*E9m@(JS+jI(#qRi-(%qRPQ&p,%uR%6aB_5%###.v#{B+_&/{w)bf+G@'+U,XAE401ZoR?$'-R%Fj.j^8G>#ZP#V6#5n%y%,}$(iG$4.(@`/|&24Q$QoRBd',Z$*S$%:3Y>$E>#}P#"},
{486.39,552.67,12.38,0.71,"?Q#]m@{v,###(b1uOG6#$###)U3qb####,##{5&1c#*Q$7u#{H(.+0@h=`c#uCTcw+?6(he#6JQu5&###0##4Z$;H$2Z%:5#qI.@J'3A+G>7?YG#90nx3mD(*@=UNA*@+,##Ip/lG$`P#+c#;I+`v&+/,c|>,[(Ja7cM?*n&*BT.7(Lx2uR%z5F~#&$##G,#"},
{469.43,588.50,11.63,0.16,"nY#zA.Ug,pN>zm+iv(D#$pyLw6S]P####nz%=.-~P####0$#9$%n177&'2tA)y,|P7n82kj<nVS:~+A,$xQ#(x.w5%/,#H5#A@&otBb%*&7(*91tr2Y7,Gj3gRSv>$#l#?D&'A/G>#4c$Od%{G%;Q#z9-CJ,Q].<v'2-'>K*/USd5%QG#Fw#|m,&##K6(+6$"},
{778.10,623.47,13.61,0.70,"ea.d;;{k####ZX.2g5lu&###-I'qL0,{9cP#7%,N,$<xI{Y#+-8Of0g5$$##K*W0%-F>#Uu#/]0Vv'iR*7q5($(&##[;4Oy-'C5J>#R;4og*g&W######v($pC=###(##@z.######nv'Ww+/v(-##hK0=I<AB6######{h%'@+###/,#|6$.,####8l#U#%"},
{211.34,636.99,12.87,0.64,"U-)CJ(gA+*z4Nz;<$(kH&[{,pC6XEAB^3b@)FS/l>$v&-zL94H&Y>#'k4(?M=kHy,%jo*Wa60q8R~,tD,OYBWg74]/f(8Po+&Z$n>#GKMWR,zkCFf/kz0sK5)I'#%Bc&0~l%V01+6J^.-y-'.,#K,#?~(G]1<5#h,$3yLL]4ml&`n)|;=pd)qS2{b#a#%;D)"},
{211.34,636.99,12.87,2.78,"7~-ZZ&k6&jW8$m&q8.BP:e]1u6'b/2kv%3f/'bA?7,+l#^o)m;:-}>v$*g-&K>Ha_9j]+hU51~-@R%m2/I|?aFG'l#$c#PN.?~-Yh6g{6xJ,z$NXGG[-'*w&YU6Q=4rn.u0-@$N)H%TG#5y':7&A%Ne6(Q$(mq2V#Nk>$Gl$1%N&d(SG#XB)WS.W>$%##Lo'"},
{373.81,645.68,13.30,1.12,"203nR&Y*ATB1XJ/CB)*D9._6cg4D7+Y$'On+]?HDI+U#$l.)>;8`#@xB4Y-(H[N}@,)p-j?Gie-3v&0d>oV::ZNrb#^G#%;,^&1jc%3S*_p2cZN`e*iw,i(.QO?Uj:X{7@/-EZNZP#)##-<-Ox/8c$l~+Yq7c8.D.(901|.-x+K4w'L^3Hd%CZN######0-#"},
{373.81,645.68,13.30,3.76,"*###H<Yl&###>H&Oi,fGNKu$.197A&e@,/12jD9p[+|H)en*v5$~[A./11,#%GD)s<_q8az8?1:R(/pT3gE=Wv%5r7DU9.l#qd(o]+#C5cl%jJNpi=zI-X?%Jo+p??RcL`&2So-hT1A.+5f/_Z&T@%54>k@-ke+6V5WJ0p91iFEm.-ae,E5=xZ(KK&4=D)L5"},
{262.06,664.98,13.06,3.92,"e,#J;5]x4$##(Q$w()oYMw>%&K1]['3d(GK+1*8SJ1~Q%SJ+4?$'^MAn.TG#C`=({7V*?s+F'^2f.+h~.@'MIx*`B33*<Jp.n..RB)nS14f-+~M1]12R)#p-99+-P?iZMw7.(q)X*>{f3Y]2j:;K6%K%,Ay0}Q&KU3tg91H%D+Ash<`e.im)=o03A'-,Jrx2"},
{262.06,664.98,13.06,5.36,"on*e>EX6)b5$a)6NJ*l7,eh7C%,q_*9s;<T28B3s$*kf&{3DV~*y[+z/*Q[*_]Ma&057(GA,@S.a0+e4;WZM$T/;,FX05cS.9l$uY#]{,hYMo*A^N?TW<<{93p-3~Mv~/|$+M,$:R%z17Up8###K-#5s=$m(c>$jc$x]M&I*i'2Af2Te&om*-00u/4P7(P@+"},
{219.07,819.78,11.82,4.43,":b4wI-6$(###9n&WuH95####Y&&Oz<######*?$[m*OG####zIT6R%@W@Xv&ag/rfHs6P#6%$NTC)=p#&4d$Y]3Yc#Pn/-##%I*0##*o*^#D0d(&-$.]E?ASJIT7u#Qm%R/E0S/D##N2Aw,#-H&###Q##id(F>####xd#O]3)-&SZ&gP#[o,RG#=R(^x3.##"},
{788.92,858.25,12.22,3.55,"0e.######J##6pc######S:#~RT(c$###|J%UC8fu&T>#ZB-Bq<######)##`pc###$##qw%c*F95#8##jM3|8VOG#%##&J%N'7#########5rc###$##uP#e7Y###A##ep/5J,ZP#'##3]*}[-#########%rc######%##epc######mQ$#?&$##%##cv&"},
{241.44,1009.41,11.98,4.54,"Oq=#########>Wm######Z##'4E.,####;7%9?&g>$Zu%?6%Q`C#########EWm######a5#2%T.,#+##MI&L^5qP$K,$S49#NB#########8Wm######S5#WA~###~5##%(eX7G'83,#+.&{0:######'##>Wm######3##yuRG5#h#&cP#CQ$K$)Y,%>5#"},
{998.78,64.04,14.07,4.53,"7@(Ic$:#$W>#Wz2I#%###96#W+HV5#3_=n##>u$D$#%IV'##%.(Uv'Nw-Fl$w+@C,$$##<Q$>JVy,#?5O;##K#%W$#$IV###j_%3h<#Z$c#%Lb7%~.###tR&_eTox+$>M&$#ZP#%n#$IV###<&#ie1######mq$&aF###%##R^&ALVX-*$#####b]*#.,###"},
{756.46,84.38,14.28,4.64,"5y'Nn/###3,#Y=:+u#######M4G###;uO###OG####|@Z###KJ.OG####r-$lmJ######_Q#5BZ$##TdV6##eY####)AZ###wu'######w&$#p6######Gq's@Z'##+IVyQ#A,$%##u@Z###bG$2,#D>#]f%mn.4Q%eY#uQ#[sE5,#[jH0##ZP#&##q@Z###"},
{117.21,207.17,15.97,1.72,"^5$6v$35Jrb#2H#,4>@H'$##/z)(,K6#$P##T')m<G###}P#R@+Ym$8IO###F-(Gi/S95M>#0IO4n*Du$pn#K]/'S-}P$_@&{jB)##ILOL5$^l&Z>#~F8')=)HJb@-M$&S~-3p+5L8.,#%##xKMcS.)F9j>%###rW6ZV2E_<a,$AQ?c'6gc'QS&%HO######"},
{985.35,204.93,14.19,3.56,"n8[#########xih###%##mP#fm+1,#2H#ll&###wb#eR'>u$po]#########RihiP#TG#lm$&B4hJ&j&1^7)n,#P8.&v&.,#(wU#########9kh.d%g,&C,#[K3ai*B1;Au#_6&P%,E,$c5${:;#########ImhbZ&,u#1,#sJ-E^0]$*Y#$CA0v#'?5#8n&"},
{648.23,246.86,14.92,3.39,"]-%mt>I#%###3G='w,###&R#r6*######K+0<Q&1,#-[%feKS05Wg.a~.ae-w8Ws#'$##8{%y::7o/pb#m*-Hx,THHnd*nA+cJ/&~*c?&)p3c<W'J.{b#oQ%sz0B9WOG#eG#b#%O@O###A##<7*/o-Mu$J>#29Wn#&/,#2##Dr<YH(###7##]$*-u####d##"},
{753.27,309.62,13.96,0.87,"Hd)W>$mu#w~-T[W|k#)##aM(0L9######N41U,%######{['tJ.Hz0?S+2e)*^Ws,&VG#~o#f|G######D:-######'##yJ1r&1rB/Vv?*o*U^Wz5&0,#y9%`~W###%##^?&.u####0##u?*T,%+##zv<vcI^jH###4##2l4>|D######&@#X#%xb#$##UG#"},
{426.28,521.43,14.89,4.46,"<x..%)V~-u@(Ee,EZ$qcBw7*jTWzb#eR)jA&@]1VG#/u#U?$V|@D['7%,<@#H90V>8)j@ru$pUW7~,PQ%@&&BU8]G#k[+~H$|L8xP#T_3zD<(q9;%&nTW[$'[SWBl$/?&Dn#}L:b5%1Z$aG#ZP####S&&+<BZ-*###`G4}*DYmV###FZ#76:o6*D>#K>#l6'"},
{185.07,553.13,15.41,2.85,"m@+%z4S7)*8)`*=9}Dn?'&6%CN>x]-5v'H]%/<?PG####+&%ld'r/OJB,L^79N5-JX7Z#<q8DMX/R*%##LB&BM6>u$###DI$G>#C5#D{%qIX/[)gY#J$#}LX/lM95#&##[BG#o.{k####AJ(D>####.-#oe/5l$qY#ZG#s`<EQ$0w+|k#[J+]-&*Q%.,#Y#$"},
{185.07,553.13,15.41,4.17,"3L(I3D3H$f>$=;+$17$#####pd'Pp3###$##$R&cu%]>#eY#_%(t16hD9K?&B2*[&UNI,-##&MO1GL###5##&K/gH%iu&/,#=$&Y(3zXGU.+`NAV)2gaHRR#Y&UYZ'.,#8'#X19:u#;-'pY#'8*}f03aBdQ'*3C;v%W2<P}1*5M######u(%uI,95####g>#"},
{381.60,569.64,14.53,0.91,"9l$::'*r@/##)?&+q$UmO=L4e3Fc6%J7,_</A*5[C7oA3+Q#ou%k),2lOM>#jmO&p,U^6+i81e-j')`FGR27MU7D8(s~1}v'Q$*fD*SEDOd'4mO=^3gR,G.&?18^W:OC8sB1*w,sY#L9-yr=in/%d%=A,o<?d$+sG$SJ+p]1.1:ZQ&hu%PV.4A/m.+o&03@'"},
{270.10,583.66,15.23,0.98,"###|J$I;@;5#]#&7;'.5LI.+~M>7n)?I*_~*>94p&0.[(jn*~c%)+1P3FK>#6wQ;z/DU5b<;9e.P,$SW4gN=;'71l#B?&TU.2@+:C$:vQ'.'ovQV.+r-)eg)#r;P92*0.nh8[~0ll$pp2<:5#?&M##]i2~vQY;>)l#jv%KP=09248)Je)G3?fn.Ph3-&/(v%"},
{660.18,654.54,14.76,1.67,"&|4|C=*c#lB6oEW%Q%$##;Q$F3;d,$PH'$##:[%)D6&x.J,$~d+YG#iE+_BW=BW###d7$A?DmBW'##<m&b?%x$+9Z$P_.T3AzG$Yd(hc91L62d)###>h*hnREFH###&##O{+'o1###4##po.~>$I>#7o.KZ$}>&###50.[R*;K2###&##r-&%f,######H5#"},
{844.68,667.15,17.63,3.38,"`$+######(##[1i######b$#(9`fY####^##aP#,H$Xv(7#$tn1######-##i1iD>####r$#?&[3%-###G##-##up.-u#.,#3'7######2##(2iT#%.,#Z$#r<B#_9@#$>c#gG#I*.<sDcP#;y7######4##_1i3,#;5#W~&]080l#d$'K*:pb#I,#>gUgc&"},
{683.82,696.71,15.77,1.51,"HiTHJ1$##$c#J&RG>####,##Dv&>%,.,#sY#TQ$0-'DH%Cc%XeTD>#Ue$xt=CfTbG#?$(m>#x&-}<:2kA,R(o`3)^3OH&xH&xS2(##k3.K.OZdT###Mc$~i.f./N>#1N/V%Q0w+oY#cK(;HH-H&(##Gf-`S-g97###%##^K(dS0###'##)/+E,$###/.$;@+"},
{266.33,726.92,14.33,1.32,"<%,j#$T~,jR)=cMwu#hS.1-$&^S######Y##A1/3l$###'##xS/uI'';2(x-/OC@v%zJ,L<,VaE######HF0#:4E>#1,#-.'a-(_:*b;@4$'o]SA7(B?'?[$[[SV##@c%Po*:?$9t9.[)X?'|G%/-$E`4p~0oZS3,#Sl$sg+nZSq##W07F-#TR'Ub5q~S|H("},
{503.71,958.83,16.00,4.41,"~./######3##1;k######S$#5A~d>#3R*y>#}>#fA/,@*-c$]98######3##E;k######>$#&Vf6u#O6'nR)*v%wd*ex*bHQSy8######*##q;k######I##cB]fv*PG#zY#UI'iq<_#%i?(px5######<##6;k######4$#JOHbP#jY#5l#<7,fY#/,#$H#"},
{344.94,988.71,18.62,4.61,"T&3#########i4o######;##YRO#S.,u#nQ$:7(?);8x/{$);f3######,##+4o######o##KSZjG$x>%UH&uA3*c$0,#vl$Ux3######&##74o######O##S]_fG$1u#CH%8n*ec':5#J,$sv+#########W5o######&##zSYGc$C%-mG${l#<'4Ld)tm+"},
{845.68,138.87,19.93,5.48,"|0-lbJ1,#>5#C;OHJ1]##fY#7s,Sd+')#s$,l,#6#$')#q~2bq<A#$cu$,m$70Z0u#]P#F##J>9D&18R'vb#<k11Z%')#FI,g7/%##^[$k0/V.Z3,#D,$e]$SB3'm$m7.~>#P`1Ig8tY#D>#A,$###c[$fW<]m+###Nu$=z)`>${Y#ax/-.+2Z#Ax/:R'lZ("},
{661.91,171.46,16.58,3.32,"j7-1@'CI)uS/SzVeQ((##yA&C<9`wV###>.$D?%3yV8Z%-c#nS-:L6-Z$qG$l{Vcv*$##b>#8,B9r?###7##&@*zn0###p##Vd++l#D[$J:3gvV.,#$##5{%lsG######^%$L..######:n'D>####A$#/JI7m).u####S?90Z%Q7,###*Q9qb#Or<###)<7"},
{40.14,183.25,18.32,5.98,"############Pn/######Za6U[+-H&{$#53~&##Fl%['#u.~A,$######$##Y/~Gu$###m~$r4GWoTac%aN<U@#F/~s6#PbMyH*######C##%1~su&###|6#^o.s)3Z.,D$(*0-;L6w#%ju&HJ1######/.#V:<.,#$##E4.#v'K>#{Q$+r5cH(%##;&%|x2"},
{476.02,267.34,20.62,1.60,"%7+###i>#qG$.,FJv'*%*%I&=u#Yv%fGFrH).,####Zc$C6(]-*###dP#[>#QnVPc%E-%+59#d(dZ$vwAm?Q:5#4##3M4UH(lm+RG#~>$1,#qoV###$I(;-%r{AI##3sVr-*$##<##;[Fg,&cQ'#Q$Vu%>5#0pV###aP#1,#lvOH>#{a4Dn.$##1,#arVtc("},
{154.46,592.32,17.47,2.92,"aw*-S-;I%g7-lg8UQ%<I&7r7@7R/,####{~#pS.BZ#X-*F##MB,t5K86&D#$u8RRC7i5%xR#B7R/,####~&#l.,.R#ty:H##Zd'-*@wc$.7R09R7m)###{k6-IKD>####[[#>R((S#x5R0##95####u##S7Rl6*Q#%###RAEg[*O#%95#(~$s>%PA#x5R1##"},
{56.16,1021.04,18.38,6.02,"{k#######%##xxa######@'#dya#w,###Z&#IH%I93.,#%##Uv*######7##gya+u####j&#;{a@=H###y##D%#gjG95####m[-######6##i{a/v(###j##]D3Nya###%##u$#dya######95##########Oc%6#$###$##N5#`/2######2##ytJ######"},
{386.94,1026.51,19.44,4.59,"############9ZP######0##T+s.,####I##p06e?)F>#Au#############59a######<##v*s######+-#OC;qb#/,#oc%############YJ[#########p,s###PG#.##'h:z-)$R)eP#############{)B#########f-s###95####+E?v5%cv(ge/"},
{264.34,381.48,20.23,4.77,"KZ$E>#I>#^P#+$(######/$$du&###s>#,dF###+##.C'|2C^~..,####&##5ZO######--#ce~###[H$W.?KQ'###r^Ig;?he.pb####%##+JT{k#######'h~###LZ%5c#[/3###wi~9u#JS.~P#<5#%##FM9fG$SG#0,#zj~+u#3,#2,#DM64H&S#EmY#"},
{501.58,638.40,22.50,0.31,"[e-a%A;v&qx.EUU8A/###9e$Y81J[&*R(U>#.u#J:(GkFJ,$[%-KB,pm'A15@UU;v'@5#dI%}(=@[$o(:j5$###oJ$;P@QT56$(]Z%`f+$r7pRU~P#1c#8X-Y{AE>#3C/eI'######4e&eX@:#$G,$v,%)zLY-)rA.jn-2JAip2)q8J-'w-%AQ&e>#Yl%Hn&"},
{464.47,726.02,21.32,1.12,"e_'YA2X,$-u#Lw%(c$kv#|x4.06###B##do.a@.###_w*DQ%2>5_I,f#&oH'|1R0[(oG$1$&Ly//`4d-)3v'R'+9x/5d)5,#T(<uY#kQ$5{2;-Rf>#QI);i)sJ,[C/^m>%#Jb|13K3d7)&x.jkM6c#TH';@$B-R.,#6##9:(B6'oY#W{0Pr?Yv%D$)(8EKv)"},
{695.34,285.29,23.48,3.36,"B,#Ml$r>#mc')Z$hp1G>#-?%+[&9g3###5##hR,G,$###nm%F##TD2>u$###z[)q.@[[,*##/zW=z895#4'$Px28-'YG#Hw>$##7J@P$*###DL6p41#A/c#$~|WH@,]P#kI#{V5xwWH>#MI%###).&RN5Qu%;e+cI+@7&iI+CzWg,&$##;e$[_9Vd+###lZ#"},
{196.63,704.59,26.75,1.62,"H.,#Z$GQ&{v%zV<H>#4c#[v&>]R|P$:5#m5#^Q%vA0r?)&##02:qY#57*U~('h8Cl#cw,ge+x'~^?%8%+wm&b6'n{:Dx-JK3;m)###3d$<1X~Z'$##qZ%M(~nf5fY#EQ$CoDC,$nG$4@%$$OgP#;&~+##58.Bl#J&~###Go/qP#%'~###fI)$##-(~K>#/Q%"},
{196.63,704.59,26.75,2.97,"|Y#9T/SZ&[Z&E?&J%*|jAUR*fbLk,%;6&Iy)Ad*g##=wXPQ$Xv#ExX<5#Q>#v%&lwX8$%S-(7]I8T4<,#|v&uw-$J#(wX,##48)_n/K>#El#s]/_&3Jc#^e,NxXXu%QG#gI%(//8y#&wX(##nR(jR+9c#*f2_`9'U6Y>#jS-2xX8u#I#%iQ#:6':T#&wX&##"},
{241.03,764.05,24.84,1.49,"Pg7tb#bG#,%'x]R/,#5,#P?$4W7`,%G['4[((##/c#_X:gu&z~1{P#+I(aJ)=~RY#$#?%M~$s;;icG1v(d##9@(oaAmo4hu%*%,]>#d/-N-Na[R7c#q[)AQ;A&0E25h5>[YHg{8<p6Bv#D}E'Z#&@+L##X[R$-&hv+>##WgOsb#`d+jm#G?OWG#8v'`B(,D>"}
};
siftPoint_save sps4[] = {
{661.01,51.90,3.92,1.34,"###o1&{v,###=Q%V9K_5%>c#b@%DHQVG##H$O.#Ir@9l$###bG$5B#&HQ%##k+LB8&q]5)(*3@*&v'^;+:9/t[&gR-1[$%H%Du$0##aIQ###_LQTG#9W?cu#jh2IZ&3`,#-&R#%fY#%T$C1:$#####KKQ###qr3=5#T,?aQ%1t84%,fZ&K-$JI'G815Q${H)"},
{484.68,60.61,3.90,3.95,"G4@###W##Nd*S4Z######K>#S%R######%##+o*6#$######v%0%##i>#.:5b0Z.,#$##f]%C0ZZP####j##nJ){k########80S,#Ol%7g'v3ZQu%###.T$pCLg,&###R5#n~,eY####'Z#jG$Mc#W6)[h,_q5Ec$xY$s~%s/Z.,####oo#Nw.)##.,#60'"},
{637.73,71.56,3.56,2.23,"#6&iY#+##OuADC6I@*###0]%au:Mv)(##aP#hw)mP$s##)@+-H#H+Cj#%rf4Xm&>,8vA/[$(ypOG03^5$Iv$?E:(/-D>#'##]Y9LlO1Z#Ax0pwM^6'':+7>I4,Kfx/`v%4P<x%-smOOG#+##M~.+u#:%#=mOA4?d@,:-$&M8yP#^s3|R+:?'{5$+bEz5%Zv("},
{431.80,82.91,4.02,1.29,"###%##96$PS0/##7e(nB+ac';H$F[)c.$R/3.,####A$#E*F###1##:U1QT5p~+x?9f^Zy>%/aZf;:4e+_A,H-)###fl#YV=######[f*R..KR+e#$4bZ.XBg~Z:u#M1+2`Qld,$##/d%o-(######g5#A6(######Kg&+x16#$###RJ$S<AD>####rH#tv*"},
{431.80,82.91,4.02,4.07,"z,%x-)95####CT)9U4T,%$##)</JJ0###&##$.)eY####-##)e'<J-CZ&'##Sp-gJA]nZ*H$LtZua?Rm*G$$0M4}>&###%##*J#dT42$(###A~,<'(trZ:6EzoZLc$G@)Cm8ED4R6)###7##S:&[?)0,####b(2ZP##H#'/-p].ll''##Yw's-$je1###%##"},
{87.26,103.48,3.72,0.71,"ZT6&1#8-T%$#0tJR,#^,%tL%wb#du%~n+Em&kI(c#&b,$AZ%=HP;-#y-TPm&x-T2,#P,$J@&Dc%oG#}&.E?'l6)?5#/-$'91md,'##W=:e0T9D<&H%;^/K#ET-(I-&C/1hf2By73,#]>$U=1j$'bw.tI*AV9_,%7u#.{2;.T~P#E>#~G#I0T-H&%##$##{SE"},
{457.92,133.12,3.38,3.29,"~>$:S#z]62u##B_|u#7v(y##CC_(c$###|##Ds?Fh:*Q$=,#TI,gW(*x15,#0B_HQ$OG#K$#{E_DI,###B##f.&K9]######t99o.)OG#R##-B_/,####n##voXbG$###'##},$A@+######`K7######Z,#zA_######V##,PE95####$##4$$`5%######"},
{485.52,131.71,4.22,1.38,")6$iZ(######rv'iZ(###,##am)gY#nG$oJ.$##%##~H%H^8,6#mS2P>#*c$DM6QL9P5$D,##oUxY$1,#H+6uG%4##WV0%,F1,#QG#gf#Az<x<DY5$j7&P.,qrUd-)CZ%hM+Gd%ke(B4DD/-######W(#?mU=-'###W(#CmUOR=al&c8#yrCCU&-L4f$(i,&"},
{226.33,160.00,3.84,1.47,"8J0QG#(##[G#<0aa>$.,#E##JOE(*=m[-$##ay88S)F=J&##h]6gY#%##)c#30a'H$R-)b##jq:MF=^N@46%0tGr7+'A/:c#,U795####YG#K0a7##y@/2##{;<S8(r~[P#$X-PCH%pH)b>#c&4:5####'##60aG5#h%06##SD?bc#70a=##SZQ4,#f5%C##"},
{539.89,160.31,3.96,2.49,"J>#{@@######l>=HSR###.##0_3j>%###.##U>#Q#%.,#$##LI,ul8[P#E,#JfVYA.###(%#mQNQ$(###o##9,#M;8######D/2gH$5c$=R%_hVjn/+v'|[$nN+4gVa5%'##Qu#'jV######_~,J#%Xl#sc&P9.-N7Ye/'$$Zm#h-LBd*WG#p,%-z20,#cP#"},
{186.08,163.55,3.76,1.49,"6B5.,####'##7f~c,$wR.&##rC<)@$w@Z'H#[PM5,#Nu$R##F^9fY####(##ie~`G#Pn,F##(h:<v$xcG1S(4ZOpY#l[,^Z#/h;[P####+##je~<##080N,#1N=u~(=f~hP#S@Sj,%{,'5,#GK5######K,$uf~$Q#)~.(##OM>bQ#se~S,#8lO1,#yG%S,#"},
{438.43,163.02,4.16,3.58,"GKNnZ(######te.;5####{l$###E5####B_/###,%(###Eq/LLW######'##EJW###>#$K_(Yl&tw$|f2c<3=5#kh-2?%lQ&GKW###_,#*d()MW###J#$JZ#{^9;Q#v=2`I,$##&d$xo&Oh:nr?A,#LJ,z5&DNW`5%######vh5Ga;fl$A,$###nE2~Q%j,&"},
{67.88,170.40,4.52,1.83,"_*+lV8sT7&##d0HQu%###%##W80###.,#+##8#$8,#/c$>5#edIRm'u_Ylc%``YD>#(##Gc#:FE6l$###/##K>#=l$E>#%##Fp8y$#mxQ,A'0~Y$##gP#}L$B;<A,$###F##(c#Nl%###$##L@-E@#~-*si(BWCAu$###}(#%@*>.,###;##<,#og7######"},
{383.96,178.32,4.19,4.66,"%##XG#<#$E>#>u$`G#$c#=$&~l&[P#D##*%J,##k>%|-#{uRG>#%##zb#.u#/sDTG##l#R?$Sf~YH(S##/DLn6'a838C#Lf~1,#=5#P5$###7cOH>#^,%;,#|g~iL3}I/Y>#$('rg~IS/&?&%##lG$[P####g)B?#$0,#6,#Df~7?#XjH&##HJ,DB/,FH###"},
{275.96,206.39,3.75,2.42,"Tw%%Q%>.#Jm*-##pb#QB#[[,###'##o6$o#'G>#H>#WG#SG#b;'X83f5#bG$G-C5?'7A#k:9gu&###q:#9bJ.,####b/#nv+o0-I#%/(#KQKm9VOG#8##;b1`f5/##n($v{;###R>#JV$>07p$%(c$b%#0~QP@*,u#I>#9?8Au$&d#(/.*v#H##nT/(L,g,&"},
{275.96,206.39,3.75,4.49,"W6)eY####3##57X######H.#<3EOG#6$$H^19Q&;5#hu%Go,b~1E>#/,#3c#o7X~Z&###n$#[|@]TW$##,?$IA#rsG6?'iG$)K4######*Q#_8X0K,F-)^##%&){WV|p;$##}zOMkHQu%1##;K5######N##[h>~Q#DL:i##y[-;~#77X2##^7XgG#ll'H$#"},
{257.76,266.43,3.79,1.35,"/,#oP#/c$PG#yY#T:3[>$gP#wd%`B5D>#&##d>$tb#=5#1,#1,#GQ#%]2zP$k&,q1KiB^)c#)H^`?Kb-)jZ#CA0U>#F#$|Y#/,#E5#B^52u#'~-=d$)H^D}BNB^5?%q&*&pC%&1T-%Z>$I,#'##u5$|>%###.,#|Y#N8&y83fu&B[%2$&%25*v#se-A,$6,#"},
{257.76,266.43,3.79,4.42,"OG#K>#$v#d@+s#%Tz316'H7'Q8&sx395#)c#f,%###%##HZ$E,$[5#4S/#R%Rx)OfA:9^DH%x>^%>C0e-4R$qy6jY#$##U>#0l#wP#g~0D5#T$)2[#x>^&.MO9^}Y#Fo+OqJ%]2@l$1,#GQ#<5#&##sG$;#$95#'##SI%Rp6gG$sY#`G#9(3hG$E>#%##*c#"},
{95.56,282.60,3.71,1.42,"$##Rc$q>%ZP#C##rD6AH'.,#<7$(19OG####pY#@#$V#%.,#'##qd%OK6###BT+6rN.9]lP#[>]1-LsQ(PQ#w./i5$?u$7##wP#X,$a:9###Bn-o-%[>]LtCA9]Yu$2f)$KAbS1`G##$'X>#=H%$l#-c#.,#0,#hP#+0&S]4g[,pb#fc#i15r,&Rc%p>%B5#"},
{321.56,311.20,4.25,2.19,"###(##3u#D>#95#cI%7?'(##%7*x]1A,$oK&GZ%GH'8##2%J###7##W$*qb#dS1=x#}h@B##r[W'~%Fn.1r%rm(uyRdB5r`7$##H,$vm+###m08~5#H@ON>#I`Wc>#F^24c#?K.sR)q@?)n-P5#oI,],%F>#`Q&Su%q0.kY#^@>M$*v[$dG$x2(QPNzn#k-+"},
{262.01,326.94,3.83,3.33,"Fe$8o1#j30?&v9',S+1H&&##}:6MJ(tc(4,#Eh86.%2A+~7-zJ-%v$?|VQ5$meL,?%1L5T>#7i3M$&cD:6,#,[:Ee/'m$f82C@,x:#.wV|5#%$SD##m82#o#j$*qb#43-#@(3L$Gw.KA%;81>u$81#696UA%@S0*##,u#q'$uG%###9$#]v($##.,#d$#C^9"},
{235.07,336.77,4.24,4.75,":y6(-&mn-4J%R+3*o0/A,Ad(>v<vb#rW-d./z-)<,#t59X>$B7RF>#626~u#:V4^P#JQA(H$-;Rhd%vK4cZ#5R)09$S7R(##zV3{c(rH'#?#j@){g9|1/Qw)iEDN@@s?)3x&Gd*PM'9K51v#@c#c5%|>$D5#F,#rp6a.(4l$N>#cE:wu$S80Hc%nc$D>#|T+"},
{262.58,335.67,3.84,2.88,"o6&j,%Z(NC@,d:.8?&zd,y-)VZ%Qc%###}H(.u#DH&#c#tu$GC7Wg${JWj,#S^Vb$&J4JI##5L4mg0pb#.###D<e.,O5$Jw#o#'(j*$=H(I$~WD56#x84CB$TXE;S't':j>#_JWjc%-f)(8*D>#ZI%~P#y;3_5%%##An+b0)uJ.M>#w=Cvb#7:KU,%{u#C?'"},
{274.38,339.09,3.68,6.10,"hK5[G#_[*X>#J$)sS$K~08##G82;@%(c$L2&7l$3o-ZP#9H7jZ(*##}J1aw']L;u,#o#MVl#H%WU##1x0<D$_-*^F2A5H-8'###<5#Cv&<@(%o14,#nF6^H'K*WJl%8A*_H&W^)JJ0}b6N(=###jR)<5#G>#QG#6g1N$&$l#{x'r]4]#$?Q%gB#h$WM,#|,'"},
{489.81,338.07,3.61,4.55,"=Q%F>#]>#(l#e8N######+##)0Y###<5#;S$fH(`Z%II+Mv$+u####YH&I>#+/Y###$##B,#60Y?,#=c%4H#t[,jK&|ZT-c#ZP####R$(|G%)IV###2,#^%$q/Yn>$vG%u$#PR)lm$Z/Y?5#.,####3R#L.Y#95###9##FgH@V=$##<5#b|+i5%1,#<;.e]1"},
{506.97,342.62,4.74,6.04,"bH$v6UW.%oB84%'+XC`m$i6N14=v5&###wx'8H%F$&Fl%+##7-#<AZ###$##.t-)AZ95#,##kFZjv+###R##>@)o6'#.+*##ZG#Ru%###Mm(6,JCZ&###T&(sAZ-u####E'#sv+Fu$MZ%<Q#L6(95####w-&je1$#####tE*f`DhY#.,#<)$fu&###?5#J[#"},
{187.79,347.85,3.77,2.70,"/v%@Q%t5D|c'[M4;c$Am)&f/)m'$##$##:y3###)##5H%?$(#v&/T$X~TBA*B]Tw#%B(<c5#XU6pv+###}>#zG%Y#%^P#<5:.6'G=0Uy8vQ%MbKqH&&[),w#r[T]w-W>$t8#XsGX1/[7.~k3B,$yQ%2,#oX;s#'pP#5x/%g,._T6?%y[-{##>0Rh[%T$GMv&"},
{199.50,348.99,3.74,5.89,"0H&7$#2wS*?%{#'Rn%(B4t>$Oc&XA)DH'JW*vP$9v(###ECM###i6#`(:zd,ln0}>#-wSUZ%-wS,%(F6(8;&4e,?zSxG%GN,/,#[#%Vd)O#%.n-`5$`RE<Q%PyS{H&[7)H6&(n(cg.s8-B1;$##b;:PG#/,#PG#qaBk$(rP$8(-,193Q$5c$-)'o99;##U~/"},
{140.98,351.57,4.15,1.67,"oU<o-%TL8|@$eC*Ey1z81uP$IW2eY#$O6L5$J?(:##;5GQ,#w-HH>#'t<e5$Mx-H>#FwB}Q(6/PR,#N2>FI&Ve.R##:-P+##mE/Ad*bv&.,#|]$|_?oz(oM>xw0{k7.o0Gr51H%+x$hiAYG#W.#)R*/.$8Q&lA#2-P?,#?v(+6&O$ERG#*w'RQ'06%Ml#P?$"},
{140.98,351.57,4.15,3.90,"+u#%p%lw0_R%:?'sG$Jc$ytAjG#h%.`&%>D?8R#`m+tl#h,&X-*1[$T,%b_(]?%jk:<6'yP?K{34o1Da0Q_>U?&.u#:TL9Q&{^<].#rg;AA#LS/9D(i^;f-#NvOA-&`M?r8#I-)#$#}uOC,#bU5)?$VD;w#'_I-sQ&e8/}<5Pf4(M(*y6Xf#5?'P/#<sFX##"},
{509.77,356.31,3.65,0.17,"@##kAW.,####C[;^+L######Y;Vn>%######jl&]5$[G#?5#L,$?7,###.c#;CW}>&###|G#'EW$#####G##=&13,#95#L##C,$######9R&)mS######|'%_@WQG####f(#PK4W,%###S##>$(bG$###K5#<@W.,####''#HHSO#$i-*Pg#ev$qe.U$*(##"},
{200.92,363.41,3.80,4.09,"H7,###bR'|n+^JT###<-&4u#5(:95#>5#>u#y:5*Z$u,&QG#JA0*v#l96W#$+ITMH%*n,I##[;?b3B.,#)##2?#C}FN$*###L.-/[$~]4G-$[HRPm'-v'020F$(-RGs%(7ITV[#>mR~n$5*E6A/ER(a5%;h#;@,>-$b7-}d?H249/,`@(QJTA02[P#a#8lcP"},
{223.71,365.36,4.02,4.78,"PR'fG$i%)B-(Yr/&?%z@+ub#vK2[m%}d,1`,(0,>f1#$&UTF0d&Fv%Lh6&Z$dsBaZ%gK0Nu$)SW^G#}J.8g%'L,rx1]vH,m$?H#H(*MI,R#%|}GnK1mc'bQ$nUWyb#Q.*:m#6'/^>$(WRxP$&##t#$?-$Oz6,x1],%)H#ul:$lDJ#%dP#}8$]&&Jy68T+7l$"},
{238.19,371.59,3.60,5.40,"9m$gZ'LI&2n-sc8'?%q[*?u$eH9O/3o-#tm,y?#;`@y9$*V=4Q#FT0;|@}k#'`4jg00e-7H%8aWFR+:5#>Z#Wf-]*EyY#@Z%U,%C@$+#NhZ%:B5j#%Zc&bb7'~W&Q%$##R;)fA00e.)##tG$W>$J##{v(FlNd~1f#$Al$G+;.J-CuFeY#Am#yu#p[W###`P#"},
{253.23,391.63,4.03,2.66,"Jo&G~,pn0>H&IeC4l$*##TsA`$)H,$)##{dP###g,%07*nB3{R)a$&8r8yv*xgP:H&`>#7w%5D9pkA%##yb#Ll#,eP/c$F>#Qu%Ed#zv'F|68bK^v&4~-^r+sv*;JD7(;B5#?s@:*?3l$i#####<,#hP#{fP95#Y?#fQM9:1kl&*^'-#LW,$FYCL#$h,&_H&"},
{461.00,394.99,4.75,3.30,"###2rWKu$E>#s?$3qWL>#kY#|L,*dNhd,C5#WR,Dl#Nv(hR(###RqWZP####%h5*[MD>#'##?qWAR+###f5#gtK;5#Z>#,$$###P#$(?&D>#@2@vP#%Z$7-$8oW.u#D5#{(+E<AJQ%Je)#<8.,#'##cu%xY$n93`P#vb#O-$NsW~Z'H>#`6$K|6'L3:7*~Q&"},
{135.01,403.69,3.95,4.06,";@,L$#&K4;I$P-)NH&<c#VO@e>#&8.$f$}3HE@$rQ)nQ$-H&%~.p,#[l&ez&ld'A^-7R*M,>:<8zZ(/G5;<Cm?(95#PTJvG%I{A6$#iB8.p#IK4BT%0*EuQ#/dP.H$LjB}~#I?(^##7ePY>#4W;W5#Cg5p5%7S/qu$&U3;N3+{?tL(dK7Ox#%Q%R8#2cPl5#"},
{298.15,404.87,3.23,5.30,"&/%zS2<A,2$(By41Z%FR'jf,zC:(Q%XG#+S&Y,$l@)Cd)V5$k.*%T-6<>@c%I^W;c$cu%[]-_&1sZ(78#5ZK2##^p,l'+&K4V?'a:/_]5-e*c]W|#%1u#i:*X[W###K%#CnM###*##c_$R[W`/4U,#HI+iW3M[W###$##MX-Z|G###'##-]%,u####],#gS0"},
{245.17,439.78,4.02,2.76,"8/#aYJ-%-###/O/oT7###%##|-)8c$+##a,%###&I(0H%zY$6@'t00>lKO#$dzQfw.YZ&Vv#5;91Q>###)##z#&uxQ.,#1##T,%}H%IL5_`4.,N#v$Z&1YE-8@+rt5)3C.c#EkG$q3Bu$Vl#.,#-l#F>#TyO.,#(##X<>jq0zY$<##PxQr-(H03F>#g>$gz5"},
{432.08,533.10,3.99,2.76,"0I$pKSg./<#$T5I2D8w[,9Q#ZdS$##/,#[K)FQ&###%##Iy1lU4cG>ex0`g3YhSJe/###a@'+iS6#$###RH#me0######S?$j>%&##=~%MeSiJ1+u####+iS?dS95####}V(.n-######T[$%##.,#(w#mq=6.-###'##B48<i@######^~$j?(.,####.##"},
{450.97,546.52,3.83,1.24,"*1#@@PT,%###jU$itM######@?#}YN######-##16M######Kg,J15=W<w@,ZCPK972,#}R&X8.@tD<%,e,$[5$Va@KZ&S#$H?(###}x,bkCxkHB6(2##H.K#.(*W9_6'906^u#0]1Il$2x1W>$###c6$yAPF>#.,#U6$IAP/,#^P#eS&VYLQG#fY#;##)@P"},
{531.96,548.57,3.50,3.45,"6V=#########pMl######.##I@X$##Mu$@5#D>#6##eA,(w,[NE######1##]Ml######s##ILf###%##D,#@u$###4##Zc&6OH######4##[Ml######]##NBb######I##k>%######iG#KrB######+##ZMl######g##+QQ######[##nP$;5####N5#"},
{153.89,553.79,3.98,5.28,")h([l&######+<>:5#=J)TZ%)Q%7,#BQ7C{@B92:5#Y]#_pVp=1Wv*###+##$rV6#$+##h,$1^795#_%#8;=fP#95#D$#'oV5=>cG$5,#b8D@nV######x))Zg7bG$###s>#$?#0Z%###$l#OG####$##QBLSH(######-r,?c%.,#0,#GZ#ZG#95#SG#TG#"},
{280.96,553.38,4.28,2.87,"&n%uyTJJ1'##%O>ccIX&3nH%^XF<#$]>$dP5gc'rb####*o&iJ.;+:pi?@:5(zTa./*##39,O{T######wQ#k[,######p?%bG$)##zx(ewT280$##BI#XzTyvT###&##6N):&0{k####1H#7##}#'E],D&3CZ&)##Ym%1y1-{=###$##Zv$v@+eY####$##"},
{313.93,555.16,3.71,0.10,"XH#1:Y######GY0MOH######y_4vG%###8,#;l$D>####>u#^/-|yOvm+S,$|NZve1###%?#QKZ95####A##&m%uG%###<5#sv+Mu#t]/W$I*JZ95#%##:}2?JZ######I###[%W>$######%&/?5#N?%4[EOsD3,#x6'gr9p,L###.,#/##a5$2c$,u####"},
{466.22,564.72,4.08,5.88,"#91c@$|kNEZ$fa@#x/;-%;_7ZwX0c$bG#CL)en/######P##l/3iZ#TwXZ?$7mK)%,R$&.R&&zX7#$@5#t5#6A08u####^,#2B5+##/0P3](,<CH5#lp,RC.QwX:#$c,$^A$&~-{x)Wv)S>#G6(&J)h+:ou%>]0OOC7[(dI%HY?-wXSG#[I#Cx$^)@9%*iu&"},
{20.47,569.29,3.96,0.07,"Lp8###.,####YFr#########TGLD>####%##/l#gY#E>#%##ZL;#########oFr######$##*]]######2##N,$OG#$##<5#@'7#########2Gr#########.(d.,#%##1,#~,%###&##0,#Xn/###.,####sFr#########-T~95#######S,$OG#######"},
{298.89,574.40,4.18,1.34,"Q^#2ZN{k####,q%YYN###%##{5#5#L&##PG#1,#'QK$##%l#{X<IM<CA-R~,)^N,y6)##b~'zZ&c|B7d$|e/~G#[L8X>#Qn-tm,%##KI'W7L-mIwT6Vv%5a>6T.6[NGR(Tc%`c%qND'H%_v&######i5#v~N~P#jY#n31+|A&0/~7-<=5c7,>y5Fv(_,%uU*"},
{298.89,574.40,4.18,6.04,"###j5#'oN###%##?@'k7MSR,{:;g~-&&-iO<KlN+u####^g#ub#x*2(3A~P#1g6o-&Js>fQIEB6)##i~'?pNokN###%##J(&h@-EO5Z.,TT/GmNQ[(JZ%6T.&=C8d$q~/,d&%uK&##PG#%6#;m(^,%(_*aK54jD'H%h$'LQ%7z8Y>#On-~G#,ZK$##ub#1,#"},
{477.33,581.55,3.92,4.44,"O7&&b4a94%R)<_,pg3Ic%I,$gP#ZS(~m'`&30,#2[&~@(He/%K-@0+IAB@?'m(Wcu%&o**@(zI.Om''/%taG.,#]H$]6'l%W.B+F'+l$M$Q$J%Wu##^(;aR$oOJN5$]##M*@n>$eY#0##&$Po#&-8&D^4q(3g$WQ##GH'VW)g&1Jo-Xn/R?#/y)mv*rb#<5#"},
{277.96,589.19,3.53,3.27,"5r&6g7+c$###$J$t<F?#$ZZ&Ze']?)$##ue*~c&>5####^Z#>x$@(;<K/D-'#*+n.V,y2#-&:3V*K4M5$T,#Af/yP$95#'##Bl$HA,M4;b-(:B3L(/Hr=:W($/Vtb#M#%bz#sJ24?&###A##Vc%p~*:vD#6%_>INH%XR)g_%m`A######e0#)?&H>####~,#"},
{302.95,592.66,3.66,6.06,"47+;[=aZ'+B/bwRt?'<?&5g/mq;-[$%14sA/2vR'##&l#h$$DL6?%+i&(*X>GwRJm(ul%U[)a<?Am$(05)##*yR*##cG$)##Gu$.u#%v:DU9#vRhY#nB-sd&'D=(##;`8dQ%wvR$##`Z&_,#.,####h*/Sn.r'9:5#3{'h947/0f5$G]0uB/MrB###/H%3z%"},
{451.30,611.12,3.42,3.12,"ZP#-##?Q$OwSdG$=5#~[#OxSxc(/,#]$#pySpb#######6zS@R*1,#0B+5wS~m)v?&&@L<'5cwSOc%Z,%?^-Vf2B-&6#$yg3~<2}`@@x0-6%E-([R'rwSlB2avSG>#Jl$N<.9&/9?%cZ&pd(_`:4e-OG#tA%:].<Z%R,$wq27d'95#$##ms1]P#D>#:,#6/-"},
{137.15,620.51,3.80,2.74,"qS0FG9Ud+j,#p'+I{:Se/q>#On?_Q'NQ$lu&|c&%l#Fw#@rAeV42:0Fh=F##yC2,%=%T3@##PsWZw*A,$@##kC2D$)H>#.,#j,$L-&G,Ibv)'{>z8'pXENw#xoWJQ&vb#k&#{^-'/1###%##k$'se-nV9p7+0;?YQ%LX>Z&&k}E[?)D##}U/[A)OR,(##4J)"},
{340.08,637.45,4.41,4.58,"-+I######,##O|p######/##53DF>#_P#V,$$##?5#,K.7#$ScP######$##i|p######+##KQR/,#M>#&Q$$#####oc$L5$i>P######)##%}p######(##udW/,#?5#K,$0,#0,#iG#L5$.|D#########S|p######%##qcS###&##uG$###$##d>#yY$"},
{83.71,650.04,3.82,5.27,"I{*ec&|Q=^Q(1FXmP$&v&5l$.B+?R+)##}H(EH%}k####%c#y@-.@'{CXxd+*CXC#$CT2fd%~dML[)Zl&+##=C795#3##(H%bR&F,$~O./;='AX9c$$U&7k>|cN-e)xY$xc$}XI.,#-##T5#`*/8d)#d$c,%J6'mG$&~%1V:M95G>#_,#$f,NV<######Z5#"},
{129.84,654.94,3.83,1.37,"].-`l&F>#Ay%K7,.~+X>$gH#0[$]7-ll'$##5,#(R#6pc###+c$G,#nv+Ty)Km&g<3&dS:l#}eJ&:5N..G,#VZ&=R#6pc&##,Z$s>$<w--##,[)^c#1sch$([pciP#yT4,9(TH(`##:pc,##5,#Uv$z,'###uc'R7$ET29H%:p5&Z#ed+4R%S5$a6#6pc###"},
{129.84,654.94,3.83,4.51,"{]c###S5$h6#zv+5R%9p5'Z#/B2:H%%m']@$z,'###*##Xv$#^c,##I?(i##3g4/9(E^cjP#s`ci$(7d)_c#;w-,##,Z$t>${]c&##bc&GR#8w-G,#-oJ':50mS:l#Jm&i<3nv+g0*+c$H,#{]c###6,#.R#VZ'$##<d$]7-cG$uQ#K7,#S+F>#P'&Gw,`l&"},
{42.35,663.71,3.75,0.77,"(c$*?#om(qA.+$(~($yiEj5##{??9#1f2]r):@*fG#Z/1|;7###o##M3E.,#Q~/o/#jIT*v&xJTj6%kA2Ex'0S*37%2aEN#$^##.@(/C795#}l'8u#|-DWIOdK5@l#3^.9LTW5$.S#V16K-)Qz#Px2RG####6K(=T3/l#=R(OH$,D+nf3P-(###I(#ZT5oG$"},
{42.35,663.71,3.75,2.09,"B[HAu$ee*^6'bn(y82Wn$Pq<N#%:-'KA%`T2wP$###@X.M$(aB(Wv*ei*Oc&1r(1ITIS$r./|JS,(:QZ$`p0-$(8,#$40lx1`H&7#$,i'wm,M6'dQ(T1#aIT'JT@$(g,$fc@uv)7L+&]._u%Wc&###u5#1S,###,##Tm#I1<tb#q-%Dn+.A.0,#W2.U7.$##"},
{42.35,663.71,3.75,3.89,")K3oG$###J(#NK3N-(??$lz*:u#G[(^f(b&3G>#.,#9V$,]2y^5*m(c>$8~#=g.9UTy]5Al#avC;%P,v':u#;^695#c##FR(tNEGu#,A+..%jA2Y&(pSTp6%rRT*v&}7/%0#4|D######r##]A0$s8n$*hG#S82Ni)a(>mB#l2C'H#YZ'h1$}v(O/.bG$)?#"},
{100.91,675.51,3.85,4.65,"Mk0fY#Z,$4?%<`].,#8,#,7%o*C,l#c_8xY#-6'hG#7a;My+Ta4ZP#######_,d95####P?#eXGnG$rm(jS*|G%=#$?x)5e+#;=######+##2)d######O,#~e[###m-&gd(du&%##&o+jc%;T3######$##<)d######&##='`H,#-~-kP#-%);I'fn/gG#"},
{640.46,53.41,4.44,3.99,"N#$bG98~-&l#FL3Xh9yS,kc&jS1pZ&w%(_A+3{@######%A#cd&M%*jp/Y-)?L4A,>L)<jv((:XAE<x#'#6#-PL######i$#T$Gfu%T7+(l#}x4Ae%59X6Q#o7Xr,%*A/K&#eFK######V-#tq7t8.)[&Fu$2K3#i*^{BC]%47X-##i,&o(#ie1^##L5${$#"},
{189.59,65.62,4.33,4.60,"W@&GV7tR,7?%&E-GI,###5##}%/######[>#VC;O##5vT-##t~)yB0(AUM>#5EUqH)|G%c5#_OF######2##c1=X##a?U8##GH'J5#4EU+2:*@U###ml$pG6+aF######6$#dL<V##a?U=#####$##L7$]07T,%###5##K)4Ce/######w5#h&5Q##a?U+##"},
{627.54,80.73,5.05,0.53,"P,$#q$[A2vb#gv)%n$0PDjYH&X>av'0o/Sg0`H':,#E25JT2%%(L'(4W?XG#XmI]$'P]/O03q/106%yS+L'OWQ'4?#Km&RJNX6)Q>#=JEW?'S$O=5#Em&d@'F}Fe>#x5&d_5|Y$O:%J)Aqe+}k####=6=#l#3}:I#%C'+IZ%o<4rL9EH'g#$}6*tS(Y;?_])"},
{74.35,93.58,4.60,2.29,"Q&+k,&+l#]c%l#%TR,P>#ld(7j@mu&1l#e6'_%Q*H%:5#F,#oJ$yL=ix#Xe0Tx*2wW]5#T&2?yWW~0$l#Zg)[:9|x-EH'A7&(##J#%l(#nvW4v'$I*k.#/xWzRUyP$~,$*{WV$)BI%z]3~n.(#####S/#;y7######T$#JHS95#.##wZ#pbM.,#>S$l6&pB8"},
{74.35,93.58,4.60,3.80,"G,#3fVhG$/,#u@'?r>Eo-n>%5%,s?)y7'803]B7.,#Yw%(7&9d'9;;>?'8l#UC+UfV1/1iY#|^T1eV|P$l,$]FID>#1##Nd#um(el$sm,E5#ao3~]*LdVw5#SeV'm'vc((0#5YL######v$#$?%.^-u5&g>$je1dJ$*y6c1%EdV&##6#$k1#tn1'#####W&$"},
{259.97,94.97,4.22,1.39,"RG#hG#Gu$+l#b>#g);sb#.,#^v#31:R>#{k#[G#$l#I,$:5####,R#<C:0,#Oo*gWSqx['c#u}[(mMxc'pl#JA0D>#K>#W#$/,#6l#mU70u#I7,=m$u}[IYEix[1Z$A8(M9E?T40,#kY#2Q#dP#+u#uG#<Q%$##]P#Xo$sB76?'###Q##6<8^#&######q>#"},
{668.61,99.42,4.81,1.12,"zo*3&1Dl$B5#:@)Fu$i9$Fm)Tu%###E)#Y~0######d(#T,%f&395#+?#7;1_f4.,#5I$5f,OFH###[;#v/2>u$###[;#h$+p[,;5#ol%+}3ZI+Hm$(,E'&+|/]q>$kc%zc$<n.###P2#pH(De.1,#ud,LA$DK4%I$80]nH$:eW<Z$xc'1h(CZ&~G#ZZ$tl%"},
{138.36,110.28,4.50,0.01,"n`0pb#/7#zx0&$8,x1F,$gd&=i4Ge/=#$%%%>%('f2)##8g1:f3&##j,$9</9i:Xv%%&0;0+|]WyP$fY#_8%=K5###&##7X9P6)1##AA.JI=3C9@Q#LRRe@'(~W=5#:l$zJ#zWF######<H####u%#W`AR]2,z;`$#b2BFf#M[W###D>#i(#r@/###.,#F-#"},
{108.81,117.55,4.55,1.40,"ID4(E8^P#il#)V/yl'3-$$]1Ln(=6((%#f}J`(%`$+'%#?n.I8-im*S,$@S&;.>U7.$##+?$zBV.,#b$#s^5_Q'###X2'dK5Y-(I,$qQ$#x/aHOX>$`,$i8)Q@V###1c#6]$xd-F,#'u=-v%@h/1u#Xl$r$(]D?kY#'v%h(')@V/##1Z%]1#1v(d'%O+K8##"},
{89.96,119.20,4.53,1.77,"+w&o/4QG##Z#WTH]m+###rP#uSR###U?#Q7,3Q%=u#fh/:v(DZ%1c$6c#:x/#xN7#$^5#A$&SzS###nG$c5#;@+XH#|vS1,#nL)O~0iG#,v&o<D`,%DH$_'*]vStb#}k#h9#gm+w#$LwS2##?tAQ6)B$#M+8$)>DH&[G##.>O`?FU8###X.$9u#r`0:B5$##"},
{572.63,130.10,4.50,0.63,"+K#YdW######w0#b#R(##ZP#'8$2-(J##TH(I>####Q##8w-9$#0[P+u####P+0pdW###2##DjWX-*###Q##;I*###T>#`P#SG#JA(ym,###6fWb[*###5%#EfW.,####f'#+.,%##&l#V5#0u#Gd%5v&h$*cwIHu$Tu#&Q$DjW######H##/~+###J>#;5#"},
{435.73,134.77,3.89,1.51,"95####6##`19Qu%###j##{@Tg,&###r%#T9W$##;5#f1#-7W#?&######H?$rtM######cv$I7W[>$B##6w@D?':w*DB#d7W>u$######s>#37W######?##D9W|&0A,$:,#[['v:TsG$E,$D>#######k5$`NE######=##[mTR-%-f2$##Y[*-o)~f41,#"},
{535.26,133.91,4.15,1.21,"5l#}#&YS-}Y$;5#$-#UYCdG$We'$w(,80tb#83*DsFdP#fG$Dx$>HOy>$]>$1[&_,;F}Dl>$7J/E?$?x0tGA|Z$/6'UZ#CUVeH$>;?###4,#Te+K28du%_G#jLMv(>=u#5g.dTVU%/<##-KE:[)######H5#Lu$hY#xv(h>$5`6(.,/$&Ul%jWVo#'###Q,#"},
{535.26,133.91,4.15,5.23,"I##an*_>$D>#sP#bZ'^c%2u#J#%/,#'c#D9)^~1######A2)lP#;-'Ke$0tCOR*~l&|I%{sA=RPRc&Z#$^N-<.*cn.o>$BW5eY#Zc#mu8MoTDc%&q%RfN*+@<t5}..-o/O,$-L)cq=E>#1,####nC#FpTHl%(l#;r'W{Ake-:m&q#%k#%-4@nT1D>#V##v.."},
{569.89,144.77,5.07,0.48,"a$#Tr?6#$###7@<:q<###K##$%GeY#+##}P#ZP####uP#7l$U6(dS*MQ'F5#>TWWc&###h(#}RW$##F>#y/#OG#G>#yb#^P#Hd(DZ$&&&|@.)XW]P#@c#L6%hTW###G>#V,#PG#/,#Y5$QG#$7*.,#>0$17KTSW###%0#W6KJSW######e$%.u#.,#I>#rb#"},
{470.57,158.49,4.61,5.29,"Zm$@'U$##6Q$[P#l)UJ#%4,####l)U%m'I#%###FfI/c#(c$[#%l)U######319^'U$##i[&1v'C00I6$VcL%l#yH'oh*,*Du$*/8-1H$###]xQ0H&###v>#hU-An.<5#'I(>d%M.)Fh55$(R,$###-q(oP$&23eY#3##TR&nA+5.,i-'SI'x5#lg5YZ&###"},
{515.23,169.59,4.55,0.69,"RG#-##P$'gG$s#LF>#7u#xm#+wS95#I##|P6v5&###4%#'JU.,#<##=_4fH'f@R&##+l#Vh4_LU###$##&J&;J/###0##E-'{7*=@(mo3LC.z25]#&U,#ANU2-R###3##pc9*@+###$##:Q#Ef'-i5>u$=,#ANU_m+'##vc%b4E.,#0,#p-%_P#$##QG#L>#"},
{337.83,179.07,4.44,3.22,"%[).,#$##H,#*[UF##$?&LT%-%-@`(Z2?u&/###*JH5?&n[,9]2[P#<5#,##V[UnP#8w*f$$7n-V~%/_NHcO###vv(.m%E]Uq/4###L6'@##G]Uu5&=,#&f%b@*lw.V-$i5H$##z#'*##N~U/g7-##Yd+W##'~UCZ&###?1&}Q(W'8###^`9%##=)>###<g5"},
{333.27,200.24,4.47,0.71,",Q$?z.qQ'K['~?$LS/Pc8Xv)$#####cs3bG$/,#.,#E#$#l#WS+7JB+u#S##[2T`M;TL12?#i-*0?#[2T<5#$##nP#mR,###we1HV%T,%d1$H.T=q/ku&@'#2v'W`*|OJ$#####b,#f6*/,#;c%K?#{k#PM&d#&8_'(c$%$####58%(Q%###'##gP#OG####"},
{333.27,200.24,4.47,3.69,".,#$##>5#nY#cG$######+/({k#X5#tb#FW0D>#l0+OG#QJ'9m):5####(?#n|G+##:6'KO/;Q&e1$pdTHX9###@F/.R*hA(9S.###$##oP#/hTI5#kI-^Q#_]1r7#lhTtM<###<m#Ke+jhQCc$Z>$/,####o/Jpb#UG#&##hP<D[)mn&*(8^P#rZ%Al#}|7"},
{294.46,205.81,4.48,1.57,"6J.iP#<8X$##Q9XB#$NI,%##2v'~G#7jD######;5#m#&{k#>e.0,#M8X###%8X0,#@A0C##]7.%##CcM7#####bP#<[*###@R+###m:Xsl&f7X$##/B0)L*Pw.&##W.SfQ%###(##S7-###6[*###9j0K=GG[+###@f*}<XA,$$##+:X6n)$##4,#r./95#"},
{290.41,309.26,4.44,0.97,"###&##g6&yY$*c$L6#s@-w#&]#&%~$&m(i3/K5#;.*tY#v)U###@##oH).,#:f3=@#eV@X,#e$VPH#@S0)i&Aw*2b1D<EnT*###U,#&~.###sq<BH#Q,P&##O(Vxl#%V;R,#r~0%~#M&VI8,###zG$Bm)###]$(,l#zXD###4CK#[)M%*qb#)C(F2??8(L18"},
{213.78,317.75,4.72,0.90,"###+##'n+95#)Z$?m#^.-tG$E,$2B(DZ&[)1,##9%*,?$&yW###D##'e-###Ko3f7#-r@X,##wW-6#s%0U=/tZ(S))2}GEs:###3##R]4###O)>N?#}vW%##R{W%I$-V=~,#O80O/$MwW96$###<@)dm+###G[(d#%u#LZP#ofF[e/&J,)c$q.'8W?|p6MQ%"},
{247.10,321.37,4.43,2.41,"###,wDjG$QC6'A/x%(CH'h0*@U5$K2ZP#[I#6R(=[*rY#Dy-###>5#s$*<4FYA2'%$3kI%n)m?Tol$jS,vh'Z['[(:x)4CA,###$##LS-J#%y/44##0TQsY#sDTgY#Uw'&R&:D2V%/+i'y?T###U#$DZ%.,#A#$Jc$eh1OG#0Y28m)FS%uc(kz'`?Tv>#wQ)"},
{169.71,327.21,4.67,2.12,"*##Xd+>I#wL=fG$ed*eZ%qh:bJ0eP#[P#?z'Xl%Y>$W,#7lE###*##7n*yH*uA3o[$1*ClZ'rRX/I%S7-P)&x7'$$HS]16`7###)##OA0###E1:V,#]AW%##ZWXhG#d90w5$uo,:w*d#6ud,;5#{u$LH'%##^l%B#$1N5>5#<[:)R*%o'-u#7`&05NLJ#Xv*"},
{308.84,328.53,4.96,0.88,"F,$(S&l>%,@'F,#h/0G5#WmIt?%yH*&%#WsEiH'D>#C##]u%Y_?'-#7m)I;)n[)mE47A0-&FFN6+T19j0){>Bu$<5#={5ou&EvOWQ#<o28-#jg8WU&1vO7J(>vO-~%x:=`J${Y$hx((+H%##E.C2H&I7*[P#u]*Q19yo,wM:@d*5W+=@,;])$##(()7w,rP$"},
{308.84,328.53,4.96,4.53,"G*FHl#z#'Px$gJ.u267f1Ts:~L1?n--|,881RA.###Rr0iY#AvQ~G#]I,4.#yM=*7&|xQS?&qwQ%m$A<<T[&L~/G$#*xQI>#sa;M5$X5$ju#N)-o1=}M.$g4TT5^58A&/|V5.n--q#Lg9>?#m>#;l$hY#F>#M-#T*D#[&;c%)##mP?{G$e81-H&<?#~P#Ix("},
{159.31,345.08,4.55,0.93,"###2K)ZP#k['nP#98/,##EmF|b#yG%Y[#B4Gsb#OG#r,#*93Y'9q##`l&8W+wl%u23Rn.@oG3h2Rx1>M,{:=O5$;5#1N6#$'#6P:-#zx54J#?f2Dq&t5PN-$S6Pz6&ap8O7#<c%z%&&6P(##E9P?Z%%&.;5#{m)i7,k`8gz4u~2Uq&Pn/Yo&eY#YB$}e2DS)"},
{159.31,345.08,4.55,4.57,"FYMS,#U..(o#~o,3{4/h;Pw([P>zQ)>a0$-'gv)###M8GE>#&.R<,##T1g##S^6&Z#S0RE?&y.RlQ$#q7bv%>7,-e#6.R+##IP;bG$-Z$i5#G^($q:|C-p&2D/2Vk4~%-J:0-d)iU#rg;kZ#--$+H%|k#>5#<$#1tE:Q&95#gP#7P<Dc%lc&c5%#-#}Y$mJ("},
{271.25,353.87,4.68,4.13,"0]+?$):##Xm'E7MD>#/##a$)EA/lY#zG%dn,Cg.:7*)x07Z$fc'###mm$~o-WIU-c$`l#-$%g&3u]3J[+$##IR#D;<mo3O5$R&1Fd'JJ+fn(Kh<=d'J.)^,=nc%E$KGh'HIUAx%|C<Rx'7.-_^7sv)+u#Ez#J&3QU.tH)<SA09)$bC'j+$JUk~+8#$CfE.%-"},
{295.86,356.02,4.51,4.79,"Xd)95#4##`.,up.rG$uG#v?)Aw+WI%0R)6,:fP#a%,W-$$1WRS-9H%6c#MR)CE?qH&v~-nm)<.W2m#M/2h)(co+ul?2*AdB0s-$jT/<#$=5#g>E3&.0Z%N##91WI5#Zy3G##H'1u5$y2T+d(4l#q,%-R$mg7m&3A,$###82-aoJbG$jP#E6#'M'B96To%b%/"},
{163.20,381.68,4.36,2.18,"###aP#`Z#Bz<-c$M5#y#&@V8y%/UG#E>#Sg(A#$9#$;##S4C###1##r-)d6*8f2Rm#=C9Qc%w6UKc#jv*7{&I[&r:6992@:43##)Q%:0'G>N<'6C,#sH;M-)$;U@5#fw)El#Z]-lQ'kF1&n,NB#'D=uv#l7/#8&&Q%uk4du&354c6*<;1;5#fV&IGNL[#CH'"},
{523.59,384.45,4.58,2.75,"_5%]##v?(J)-bf4%##o}@[n-}TY[P#(H%J>#~?%5&+{k####fH),##PB6,9%`..?##QSYQ,$QTY_u#(I*d>#0m'A`&&V=###LQ'q5#xS2I##s5&i>#qmKu$)t:={:(#8.yw'gG$nh#?@X$##@l$+e%I#%######]z%7A/F,$[P#CI;8-(7u####A`&kS2###"},
{523.59,384.45,4.58,5.80,"v@/######BW&Uu%}b#[P#HLI3%,G,$.,#5+-J#%###kY#b%',*E%##}Y$5{#(%+p~'P1<z_)xHK/.)~#&Vl#re1=##<c%yl#A7-###JR*9|(#7+n>#*fZOl#]eZT,$KI,($#~&4w8%Yl&<########9w(8R(zQ(jP#DfZ:5#/QF]J.9.-1##6-'Zr,OG#e##"},
{151.71,398.84,4.61,4.58,"COH-##0w,-x#7p/B/.+&0`C/MW8GZ&`3.iw/)7)###AQ:lY#@vR<##iR,{-#+*@g5$FxRL?%9xR,d$D|?V6%ln0ER#+xRY,#vPBX>$G,$W?#J)-]]5'C.aJ.N'7qk5~~/)p+<C:X1$vL=u5#jQ%Fu$^P#-##v-##XBm,&:5#qP#CQ?m,&L6'bv*J?#zY$28("},
{543.31,397.76,4.95,3.63,"tb#lU)tn1###D~#/yP-R*Qu$h^#jbM###4u#}##$x.-04###4c$5~#quP_P#E929f,0sAs4<Tq&h&5[&#NE@6q$#1:cm&dG$lQ('##UwP0,#TvP###4|>lv&mq<95#&(#Q&1vR(T,%7C#Zp7?%-###F>EPG#%s?###8*<=J.mO@D>#YZ#($'C5#0,#-T+;m)"},
{543.31,397.76,4.95,5.65,"eY#f##m(>OL&n`Bp[&Mg92-#PaG'Q#OG#e?#56&Dn$a%/pu%#l#j,#%vPS,$CvP<d&;EDT$#.M<.5695#($##&+&3,{k#GR$N>#2Q#XwP+u#Q33l'55}Don-f$'{4=###ki.L#%)m%###hi,/,#0,##s3qQ)J#$rb#Zy)1wPUG#k,$RG#TKJ.,#[q,###<8'"},
{193.25,414.45,4.49,2.22,"'?%#d$f1;Z5$/h2t..Gn+TS.},%ph7MU5HH'F5##g(fB5ub#ru#=|-Q]5F>#jk=xo21-%N]1@-(q&1{B(~ZQ###y#${A,)[Q6m&R1*U3B~Q%B]Qb?';?&AV(1NBeY#,7#(HJ###:##NK-HZQ6#$4c#|w,uvL-GIub#f>$aN-$91k>%D5#&%E%##D5#&N4(PJ"},
{244.69,425.57,4.60,3.74,"Y?#W,L_5%&##r##Vg/iQR###LT),g2TD<5c$e;:YJ0j5%-d%DL(DTRiw)[6)8'J%'4P(4:.+-S.<5#{{+-RRSf1rQ(-Q#aYH[-)k>#XW1:B4@RR'##LQ$8U+v:>D##)f*iXD1c#Zw%TK2<::/R*n6$Y/2ql#uQR2[)tP$Z9#MJ/*&-8U9'6#c,#H9*ftL###"},
{155.81,432.55,4.50,1.99,".,#$##@d#X7.fY#b>#.d%Sp7][,M>#;5#D(+0Z%:5#Z##F<=1,#5,#YR(e@./o1ym$xL;/l#<vSB$%@m)*{&ZI&&E:VS-9;7O-#&v'eV*-h<JU8E,#p(QF,$xySQ5#RB0N##<K+&n)]m=2c$=1'dZ'xw(vP$wd'Q>#^zS3,#~#8eu&6K.[P#rz%Nq=0K&-6'"},
{174.31,432.53,4.63,3.78,"$v#8`@######:,6X:<5u#%##I{(X7S%/)2g4}Y$J;SGc%}T,>/'(*?[,$Zc%1u:eJ2pP#'?%y$'Z7S1~&]D?p#'+PCfG#@fFg#&&##P7&,K1d6S###yP#&~&n&5$##F,#Dy1]6S:##c#&8`-u5&{Y#A8.k[*_rA&7+Vl%fd#`%+Kf3{k#~5#8RCcv)`,%>5#"},
{146.29,451.08,4.62,1.09,"PG#@0'^T6%R$i,%8n+E5#W6>=5#-u#a[#'6N######Zl#w%.gMAx$#m[-Lh%-e+(`1=f3u59()4TT2KN3|K7hY#F>#(}>X>$@@Rs>#Qn/&.#'r=2'&`?Rg$%0@RNm%C)@}m#U,%86#G?R&##$TH$Z$3@*.,#ng-q/3rq5s8/uw0q_-_R,r~(.,#dI#=A1###"},
{146.29,451.08,4.62,4.58,";y7*##Il%mJ%~~+Yp0yf5~h4_M6tH)x*0[%/E.*.,#*u:SG#d[R%##NR+%?#}D=VZ$A]R1$&Q~R}u#UkEYv%j~1O-#E~RN,#?t:D>#H,$1,#$3/Ef3ZV/592oK7Y=2aU8=h0+y6M:#;aGbQ#B?$W5$###H>#ZR#zjD[P#[P#V>#D??wl'7?&O-)'-#Vc&jn'"},
{131.87,455.39,4.55,1.41,"4=I^,#o[-YC#Y&/i11/1:Cz+VV28%,sM1R6)PG#1,#qD8###|8SR>#PA1,-#v)=C-%W8SC?%%8SHQ$N2=V6%l,&;?#v6S###tv>fG$P%+###^N,a:<B`5<8/px4AZ=dw,[y/}k#g8%WT6###rH#M5$']'*6'K]#%tIG##O~/Gv&BcJ.,#Tm&_>$RQ%:5#G5#"},
{131.87,455.39,4.55,4.31,"D>#%##I>#QH%###86%;m%PFC?##$x/G9%a`D^e&M$*+[$6#$DT5'##B,$uf$1n*>1/8f2+l;*`6:A.b}3/z:^e-.,#W/F:#$<$S.##Zv*W6#OV=e.%l$SH$%r$S~6%-}D'~$WR,*$#i%S`>#Ct@{k#%H%/,#J'/{R-;W6kK2Ep8TW+DA0Oy+uR.3^#zMBM$#"},
{204.11,466.68,4.79,1.40,"ke-'Z$y,$olGXTSwG$=c$BR)un+R8+=n.###hT-,?&qb####7S'&SS:c$P6'dTSAg5WQ'F#$KC6HJ&3jF&##_k>Xc%3l$(##/M9[5G][*m8+vRSw-+c-*^1$h'9j5#*RSGQ#{lR6,#~l&A6#}Z)5,#W-%nVSI97###;-&4F/ll'(##k37RJ,w%1###9l#s/("},
{303.18,491.91,4.56,4.68,"[n'-cF%Q%*##h<.$T3###&##U-(######~#%.,####I>#{5%}p-tRBS^[F#$Qb[~82e>$Uc#O+G######Tl#.,####B5#Dl$q6*A5#K(KLRP%][.,#;v#i1L(`A###$##-%$D>#/,#2,#i>$sb####7A#tq>A,$###c,#cYApb####'##0m%###/,#T>#9#$"},
{291.09,539.07,4.36,4.52,"N##f~P######vo&N6R######2I%R^8######tZ&pZ(###'##o-&3.IH7-###w9R=/1###Q##&);uG%###3@%[@,6[(D>#AQ#(w)J@)[t<lq<v7RK#%A,#$M.r*H,##P$)V0*nP#||24bJTZ%95####mT#>6RVe0###`##-;RZ[,$##nw,m^(H>#$v#Q:RLZ&"},
{269.43,549.50,4.62,3.07,"Hl?$6IF%-8-$RiWD>#0,#o0%I%-OG####JI$.,#95#%##:7)sBPS~/mY#0u;DjW######BT'DvOOG####fQ#:l#aP#RG#XQ&du&###)v#PiWEeW###1##iZ8uuD{k####5##9d$Y,%D>####hu#Y6)m%'UA1)94###=##jJ,8;795####N5#^7*+u####&##"},
{123.02,575.17,4.53,2.73,"ad$nBKz%10,#<mBf4I#6&g5#E}Gtb#]P#`o)/Q%.,####I].B@)c^+cF<e4IihRnZ(N,#~29ihRD>####Jm#m[,.,####@m&######8_%~cRb7-###O$#@gR&dR######Ti*Ud+######Od%xH$?I+aJ'?d*p?*###&6#V^33)>######=I$m?(OG####*##"},
{142.49,593.95,4.19,1.37,"*'#xcOpb####}y#wbO######xQ#%GJ###$##4,#oPI######a:-}GHwV>}H'mgOr(?0,#:R$7T/f<@p?)?l#uc&a2;k>%5##.R*(##{U52@AmcO5Z%-##}JHdR*mB.rH&>o2yc$}w-Kc$Nx28Q&###96$z~GmP$###VQ$xdKA,$0,#QJ(lq=&Q%###H,#[}E"},
{294.90,618.49,4.68,1.96,"o[$c7T###6##8|+H820,#1,#'[&|b#2i2QG#######raY###g8.6E9D>#`H#I`YS%,.,#d,#(i=iG$.H$I#$/,#.,#jTFZP#Xi@:E4iG$e-#~]Y%m'###x&#E<8)$(/,#+##DH&ZP#]n(:5#4IOoe.y,$xK.2]YeG$bu#uw%?<7X>${P#Lu$,@'mP$=5####"},
{53.14,629.55,4.48,4.46,"NLO#m'X,$;5#q9.3K1L5$hY#,L,0@+2##(d(Vu%###]##eo3(~U>##?'3AH#*s=5`3.~Us,$Y`U@//{Z(al#W&3######U5#)[U$6#:7+(h#>'7;,#CCLBa>I[U2,#mv&*d:P/3######.6#xq>&##&##pU$:n-###v,#F&05-(###.##;B-#m'.,####Z##"},
{76.85,658.66,4.09,1.32,"SG#a%$nZ(###2##Oy)cd*.,#+',%n*pP$2,#'-%v~&HjG0,#&##2S#z^:###)d&Sb4$J.@,#/yR@g1bQ(c5#RR+d?#4vR?5#T6&vB$ruR%##v7/{z&.sEvp%vwRjv)&6&b]%2.)g0)luR###3~-S?#suR/##D$M|P#z{BHS%rw*MI+wy2e6&~5#9g,gPN95#"},
{76.85,658.66,4.09,5.29,"gH$Qv)c&(YPMVW*[u%=E3=R+CxDeY#w>$F$)1.'_5%*##b@*De'jA2&I'2]-jA/uZ&61Ubv)W1UC,$aR+zl$%s:eu&###%##,w>.,#{>#5]-y['[P#`|-V{?{-UnP$_n$+a:IXC,6'###2Z#_H:>%.&##<u#&=.lZ(mc#JZ&W-(RG#5d#_U8Vo3.,#8##D@)"},
{98.93,686.32,4.92,4.52,"`o/#########sJEEl$Z@,Vc$:Ub=5#P$(wT';&2OQ#J`<(.)K?'#########g1P######K5#0Xb###G>#vd#nK7S>#gN=67(h,&######$##'Ub######^##wTb###P>#V]$+%-%##,&*#w'l>%#########dfa######T##]Tb.##Il%o6#q%/S-%g%/%Q#"},
{679.61,34.23,5.14,5.41,"?:#1_=-##k>%q('wu'zQ####t#%###3{&.,#######WM&###mS'1@*{o*Ke/TrTW>$+Q$[5#^(8###?n#'##7,#/,##=.###d#&_S.>R)XS'>mT=5#wP$o(%S96<[(%##`##1##1_.rB2###*$'Z6)%##'C)Sd*]j@###G8$fG$TrT###3#####cqP{k####"},
{596.28,52.69,5.01,5.83,",$#_%ReY####iD)^U;(##+u#i7'0Z%|##fd,6##*?&q&#TM@yv&P-F%W1I6Qp(R|/5E##OE65`;>'7###1##C##oy9%##.,#{k####Jz*=$RW8100/p_7gW;C7(;&RI#%,##;v#a$RsP$2,#OG#2$#D&15l####O1$/$Rrb#S7+Ra6&e-Ne(M~-AK0Ze,{$&"},
{578.23,53.90,5.70,1.32,"UcRX,#DZ&-%#%##,K-de/###9,#uR+$0QN?(zH%<H&Gg)93DSbMx-#JJ1B$#-w(peD1dR?5#lfRgaB5e+0]+EH'.,#'@#`RPRo4^##iU:_##Z[+J6$ihRYT2]dR~#%&8({|5-/-$%,(##;?%wu'c>#%S.3#####Nu#Je'K@,.m(~_4`5#ho.f>$=gR###(##"},
{529.73,60.78,5.21,5.67,"kH#{}E###1,#/Y1Qx3###YZ#&DOpb#2##VH$~@)-H&#/#(~.|L15I+###$Q#S/V[P####m($a.VZP####Z(#v?''$($##/##Bq;######o6$:3V######uc#_hS:5#.,#_##y5%`%)BH'%##,]3###/##ee(p~R###?##1n*$0V###_P#T5#jZ(?-$-:2'v'"},
{531.14,96.86,5.81,5.87,")y1lZ'###-##A}BbG$###BH#=OX.,####7H#1Z9vb#Q5$0##M07######SZ#e@X###%##*x';(R###9##7g2pNX###cP#&c#gd,###&##G7).JX%##$##i-$LtGe#&2##O[(xLXZP####[>#^H(xG$7u#jZ&}=E0w+M>#AI%7o*#{:95#H,#pNX0Z%###'##"},
{51.76,113.88,5.13,0.03,"[/^######*##i[Tyn)D>####K>#e'&qL6mP$eY#K/$^(3GS+Y/^;5#.,#,##adS70+k-+%##`p3y6({w)LH&C,$$##v|1L~.B/^$##eY#@##(mSp>#0r?R5#(AXM>#[H'Ud%J$)95#dx$YL9@/^###<5#s>#0XF:5#c298H$xO8K?(yv)(##6A'JWC1##Ol%"},
{456.16,128.62,5.34,3.21,"|b#X$%W7+bA/%T^I6$4n*w7&2U^;v(@Q%>o#2d(Zx._`8nu%8I)IAI$-',l#bT^Mw'(c$Q$#^X^TED.,#K##Ol#oU^@Z%95#tf68n+o,&u5#3T^hY#####%#Y6KkZ(###.##l5#)w,######XC<###Sc&c##sS^######}$#`V<95####.##mc%k,&######"},
{466.79,139.29,5.91,3.27,"=#$wQ''g1t05SH('##,~%akHz/a%##w,$Z:(s>QnG#H&.s^,MK1;e-wb#Fo)PS/?M-X'50J,n1awQ&WZ&B[#{oVItEtI*Dc#qv&tjAk,&&##WXGscAu5&N,#}1aOQ&.,#U##}K/Ni@######.,#Kl#X@-###,bKmY#LZ&s,#S0a######m##`~,bG$######"},
{55.65,145.94,5.40,2.50,"],#bf1{34s5&P4/P&3T6$oQ(yY6W>$###$##X#$6$%eY####2$%cP#g*,8m)mp5)##:-:vw/JUX-c#{P$Y>#x6)H2%M$*###s/-18(.n,'##4w--S#~#N06#>RUvy&:-(X-#t>%'r$B&3###Mv(08&tH)%H$%m(ge#ec'-e$W6)'r$uR.b#####}1#hv+###"},
{382.56,151.96,5.81,0.44,"$##1u#}k#%#####1,#yP$l7)###N5#P>#%l<~l$4S*/,#tE4<$))l#Z>#FQ%^>O.,#Wl$&55@6(:[#nlF,AF'c#c9*'@)kh6yv+###ka4YT5Y(Rx#'Y$(_6$oz6DW3cAJY7,###&t2_$'Tg5'l#w,#%'RRc&`*32x1tI,###Aw%igPX,%######6^,]P#:5#"},
{553.46,162.67,5.07,0.49,"P,$od'_Z%xc'8KWN,$oG#QZ#MKW######&%#}c(###`P#3##('595#L7#)a;t#L###kB#&-NpLW###2##h@(7o1###K>#=5#GIW###M##w,7l-+###Q(#WKW>uO###`##f/Ez[-.,####P,#CIW###$##SB(pm,###C##k%-D>####*##AJ-C,$######>u#"},
{553.46,162.67,5.07,2.33,"E>#)l####F5#G>#N%,.,#&##f5#G7,.,####Y,$OG#######2##USJ###%##>I'.*V###&##LBH=5L###<##C~+A,$###&##S,#.*V[P#&##-lKJUK95#<$#/(VFJ0###q$#BA-BOD###.##8#$D[(6v$0$(b2>I[)wY#Kv#=}6V%V###=##UI%w&V######"},
{323.87,189.07,5.47,4.79,"w067,#^f4###/cJ9,#umS$##[oS###fY####(N>$##AB4###m)@###>%+>u#$XC$##Z1QGR)NoS###t#%g-&%6L$##k%/&##AcGA,$[>#Ec$LM5u/5?2'GjA)+FjL:,v#r48v<9Uf2YH(H5#4n&0Z%###$##BS$m+M*##-u#)?#MmS0,#E#$lI#kQR95#sb#"},
{282.11,193.13,5.49,3.17,"$##B,$&##(g5~5$uG%$##N<@AH&95#%##KoS.,#$##Q#$wmS*##umS0,#^S0I9+<mS###UU4hqSn?*8,#;|7/.,J,$s/+?D=###SJ0%##VnS?v'1n-###DqS8nS:.-$##Q57Sm)Ji?;c#p$+###V,%###|#O$##[H(###<.S-l#z7/###nFG$##AT3)##v:<"},
{282.11,193.13,5.49,4.37,"72@######?$#+[RgH)*##N1*0m&[d)*e-#U0Z[,{,#`98gH$NbK}P$eY#;$#+wHx]R2-(h,#;?:LIQfH)]G#r(>Z-%m[-A$#]1=:?#<R+/$#W<C~)+:PM3$#N~R8m'mP$0&#lEF&H#Z[,e$#5J0i##~'8xl#XYNa$#A3EY.#IZR###.,#m&#ay9.$#1T4&$#"},
{186.04,322.27,5.72,3.97,"b]..^7###5,#lE.-vQ>w)P$)DR)dzQg-(Z'.K?(oW-s/5GS$>v${L<_g({f6_V4?sD_6&hR+#6&}aA7c#6eOfY#ec%P5#2j4^1<ZG#dzQ,y/4vQ1,#3PCN-%$C9######}w+eY####0##@m'%05:h*3vQ~m#JD@^##JiAh[#0]3######I##M5$.,#/,#&##"},
{484.10,327.85,5.76,4.51,"iB*_D@0,####7&DSH(###:##hz,cG$###5##=o(D>#######Gg5/u#{>$VG#ZeT7##J#%/-#;V4S24T%/)?#{t4c@-95####9~/###l$%Fp2CdTsG#(R)N0'6p2Vp*MfT{6)+hT~6(/-&@%'A,$###n$#JeTtm,###3c#g@@L#%###aq)h,LW&4###4H#Lc>"},
{296.37,333.16,5.42,1.23,"gd,fZ#S$*-'%V$(t~,]6)@d>~w'Yd*LC*,sCE>#8Z$G0/Au$G7PS>#<J0v-#,4A)])S6Pm[%D8P,e(SjCw?%d5%,'&;6P###Ec:hG$#g3###g=41C9${5QT2=B4&Q9G(:9V1tP$MT#vaGUl%JQ#M6']&,{k#F0%J<EGZ#Sw,qc$UGDv?)BI([#$Q@(_I,lY#"},
{296.37,333.16,5.42,4.38,"v[,lY#Y#$P@(`-)OR(qc$#lDHZ#Gn,9'%`NE_&,{k#JQ#7$'1tG`u%uP$Q]#hC:-M1>B4<c9%{5EK2F}31C9s]3###(H:gG$<6P###d5%zx%v*Dj6%F8Pz[(U6Pl[%hs@0f)1A0x-#H7PR>#D0/6l$:5#Ec$YL*ZNCPn'M[*R-)5[>W$(/o,Hv)&y$gd,gZ#"},
{231.09,338.12,5.47,4.39,"~o3vP$$Z$N.$uS/=H%oZ'#><(y1~,%h_(L+G^L4D>#a(+$Z$-HPC,$Ml%i('%+DMc$ZQPK%&9SPxl$~lAc6&O/2RH#^UP/l#-p28H&`5%3&(YK0%_:,h8im&FuKl&I3e-W7%{m,<;(6D>:6#4n*;I(-6'%###H#%s>W7-qP$Xl$5PB#c#+9/s5&,m%0,#+''"},
{222.38,341.59,5.60,4.27,"E.,[,$el&^l#1x.;w*L#$Yz5`,$Zy5i~$mrBI'0'd(b?$86&t(=|G%^>$>h%e(:8-%5L8C=4rvPh>$%Z8N{<6_9###awCOu$kq<jG$4Z%C)-s)@@m&tlPvI&.vPNB*saDO%%;8168$KwP}P#HI),R(o#''c#gR&(tFu-*%m%k7,BzPf5%3%'*m(RT'dQ(VI$"},
{557.64,354.26,5.28,0.90,"]##^)A.,####j##?NB######)-$vd+Fc$zb#]P#$##&e?kG$f##dmTjG#A7-c~(RmTL,#e%/#2T:NAF>#P-%C.(:o/K&)3l#'-&b-*&?#5lMDe-2.,R-#NnT+HKs{7KT*c{;0i')RG+w((##K`)w%1$##*l#9z.17,8##Pm)3?#b/3:),ol'7S%_Z&JQ:###"},
{557.64,354.26,5.28,2.51,"wP$]P#$##@v:%c#-n'MB2ne''##ND&c6F&J)###'.$[H'>H8,l#26$=7+Rc$9$%tzT{{APG#G<>uGJ-`5s&*2v(>H#TT5C1)###W##3NA###g7.We(i[UD,#T~U$J-5w,,-#u-*cC.*@+1#####P##*{>.,#&@+P##o[UV>#o,O9Q%(%,(H#c>$D{'t/5%##"},
{557.64,354.26,5.28,5.74,"zL=3,#)6&,g%.o0'Q#NHO'Z${eW@5#LQ'/##U:;.,####:##g6*(##sm*Jz,en.I,#VfWW%-7fW7,#VI,Tw(&`>95####;##<1;)]%nv+ac#1V19T)}XDF5H$FC|k#wu$%`VGR*aQ$R,$SQ$ld+<G26#$q##{ZCmo*~P#21$829#S%fP#Y.'.,#H4/1c$^P#"},
{560.98,377.72,5.59,0.09,"TQ'.,#D/&~-'X0-L5$mT(yA+yq8B,$=5#Ze&sb#|?#s<G2u#}f0%l#Z[)Jl#q.)L'3Aj@Sw(H};(q:'Q%em%)$'n8##7W1,##G5'v%3I*###v]3b]'|6W:##v7WQ6&Rd+I$#Gv)~(#|6W)##3](.&/SB+PG#${<Jq1Ke/*I##7Wi##-o2T%#bG$A(#|6W&##"},
{560.98,377.72,5.59,3.46,"~x_0##27,*##m-+p##yx_Nu#Q:;8,#uS+0;4Qe$Fl%3B%=v(Px_###'~-<##q.0###Uz_OQ&(y_###`%*+.(eH(###|:'<d)dx_###<?&7,#'m(3Q#y39}:=q16K8):&-Ah:1S&Sv(An,95#fx_PG#$##*##Xc&iQ$@&/,c$A6#0x'?35sQ)cJ#3w,/-$Yl&"},
{179.82,379.76,5.18,3.90,"v%,t09###2##i57:iAxG$dG$~&(^xT[Q%o]/;?'2yHz706x%[/'?#L;Z$XH'/Y8dq>F5#dH'Vd(TwT16$8X?Th>&^1#?&%|*XZ&_G#uR&Z06uvT%##C5#J7(s{B|k#)##:9/DxTZu#$o078%v>%mQ%(m'c5$_OHQH'oY#Gv$H&.'r?D5#eP##a-LM?R,$D>#"},
{509.60,379.23,5.24,0.29,"ER)O~-###V##+3_mP$###E$#90[{k####Q##c7,OG#KQ$@c$Pf3.,####|,#&0_###F>#{&#B-NJn,=.+U-#5V'Em)0Q%###.K1#l#E>#&##t2_###=5#@##h7VB#$c:+)m&RT&46'Ox$uG%-m(/,#/,#(e(m0_######Eg'|/_###)$#Ue(-##O5$-($#e-"},
{238.31,389.42,5.45,0.68,">I+`c:=80Ym#%.#GxQj/40,#h,#Kx,YwQI#%X[$8B3}d(mR,oG$%@>{I*%6&3E+mN@ed'YS0f@KxZ(4'+3h7IH'}Y#-D.RvQt[+Zu$*$%zK0&.*+u#Tn#|vQ,vQ+c$p>#A)0^%/@@%76QTd)+z*+057u#0H${8'9T4wG#.v&x7+vtJh,%=-%NR%gjD}U;G>#"},
{556.49,393.89,5.14,0.25,"(.)E>#.7(}P$+@'6'4u::]H%T<9Oy7T,%~>#E$(l,#R@Y###K;)BH&d?)###W.+~d%`@Y###NAY7Z$1I+%##W~0:$#R@Y###7L'5n,BJ%{k#@01[<7`C<8##a@Y'Z#F-)J$#cR-0$#U@Y3##0##.?&7i%kI.Pn.Q~-o%%.>?D*FSG#u#'3;&T,%l##]@YT5#"},
{189.83,414.33,5.36,2.37,"KZ%=Z$@h8gZ&'n)vR*Od)K8-%l#3]&1M<*Q%{g5I$&0w+c.(OQ$BC*v994,#~a<7f.u#&R]16-(Xu##f*HmRfm+N>#Zd([fRnl%IT(j)@~6'mnR2-&8-''()K{@###/H#KbD0v&?l#F>?%_;E>#B5#7.)qmRR#PC5#UQ&Bt2W*F3v&EQ&%pI|5##o0[2)[mR"},
{542.99,414.38,5.34,2.19,"0,#=n,###rqTOZ#5h<###+gQ,m#P6)###_M9H>#$##=5#l-)$##=qTR>#?@+7'(nnTb-&ru&7MP'(;?5#I,#,w*###wb#vb####3S)<%*_eS2B3~@(B-BrL3[nTY>$/Q$0C$P]2D>####~G#.,#X/#_08cm)dL<z@#195D&$itJu$+:5#u.#Rn&$e-###$##"},
{453.23,420.40,5.56,1.75,"LSJ%m(######u}[j>%######d',Jm*######,[%,u#######NV;.,#=##KSR$|[A,$###eL*`,Dam*-u#-##bg1-u#######A,$%##3A%^y[gbHH7-UZ$Y$Al<7~PMmG$I5#i6J1Z%###-#####k,#c7)ug;1,#QH$N00)/09f0rZ(*d&4;.X<D95####}I#"},
{302.83,428.91,6.35,0.72,"mr?+@+j$%w'1SfTTh>$##Pm#gd%k-*PG#L5#ku&kY#.,#P,#{8*?GIC$'dm)&+/YIUr,$zA2ANUPA2###7$%Hw-######BQ#g^88m$vU9L'0V_8~,%yQ#cLUyIUOG#A##JP4ww0######Wu#{`3x~/m(:iI)WJFo#'4##>L(RWCOG#+##Tf$D6'],%###C5#"},
{302.83,428.91,6.35,6.21,"kR+eA)%i/~w-Yc;Sd*76&f['E4:+u#M#$.#6RZ&I,$jY#>w$mn/3f*T[&b24NY=Mm)da0E02>2SA,$S5#L~%/J.pP$%##NQ%[-S:5#VQ${k=_v*$##G<,32SVND###S,#H6;LQ'.,#'##ym(^Z'###+##(@KL#%###3,#*}:B,$###+##Qe(xY$###%##|u%"},
{212.81,442.94,5.81,2.69,"N*8AK21,#VG#T5#CYF4I)eY#A##fz;.A'*6'oc#^%/<m#:&2sD;@JHTx2Od%t#IXvMa#&<Q#B7']wM|C='c#<N-B4Gec&.[)KH'J-%qxM`x-{jEPI(P.+5)/Tc%,p',vM?6'8~LBR)p~.<p*-d'AR)ZD7)Q$38)`92hA0)I'=.*xZ%qe.$wME6(=5#UR(%KK"},
{212.81,442.94,5.81,6.11,"l5$4-=S[(95#D~-3C/$T*Tf/o'.`7)$n)n.-h<13I)C#$j%+w/-=X1CM;U5$'JN?[(hZ&}v%K.*Hy((%M$w(IwFL&+Le-r~'YR+}G$#P5-HN+HNkY#|l%fnLoc'O>#W`9xJN)$(;,#>_0u$A)##KC6gR%Xo4*[$>]48##gW?=e)z,'%##>V5.H%O5$n#$(%("},
{244.45,443.34,5.07,2.43,"07#ut=E-)$###4.{^:.,####+8+hu&(##UG#'##7$'f5%95#gQ%5b2jg9*['}{UnJ/{#'1T&DcHksG3,#V>#O##4D<Td*.,#gY#T#$Sw*NxU>h;xp0nL9~eC|U4+yUMQ'u>#z^/~OIO5${5#gP#we0K%,`ND0,#'@%KxUK.-O96,[&hn-H)4505RR&*[)t6%"},
{194.87,467.29,5.76,1.68,"gh1F-(R>#y:0Uj4#I'y6+;u#BT&/S+g%0###4x'AH'######z$'4TO&H%Hv&<`TbL7cm*X>#t2<Hc$&bCqP#VbFOG#`P#e5#bv){y.^p-?a>/~T8v'Vl$f*.+g7###cz,iL3oB8.,#e>#Co)/@+1v$`J'bUM,d)E>#e7#<`TL5$4##qf'|w.l>$/c$zP$Uu$"},
{430.61,548.38,5.38,2.68,"%##X?%(7Q###G$R}l%sH'4?#nKSD>####D###J(bG$$##1,#r$&{^OKC8$##(JS&S-_?'[/$3=D######<i-MJ/######>I$%D5js9g02A|=PLSZl&###Og('IS######Zo&bG$######<@&C,$###|.%@JSI<E######O,781;######_~#yY$######ic%"},
{277.69,569.86,5.70,2.83,"P>#Uu#BXD-c$vI%x.*yL<Zu$7^K).+.,#6%%^A1vb#$##A$%jH$Y(S}p:hY#m`>dU4U07(30FjEqY#$##BF0'n,s>%###v#${&/R>>x;<#FAt'SDd*$##)h,]$S######YS%~P#.,#*##nd*[P#<#$(g'~$Sj3G###I##l[@7OE######~6#a5$J>#2,#_>$"},
{500.60,576.54,5.38,5.16,"5K(^~.9-'56&Cp/<#$9,##h)'d&Z%.95#`I#c6'_5%###>##D'*cZ&zy0=.-(l~T,$k?(Gc#QK/cK2JK-w6)LW:[l&D5#gZ%+w'pZ'w~,ee-+g^.%-0$&2_*@x1#3B77$:]0IR+UGK1##9n)UB$*f^)##Xl%>w(xf^%##zu$=u#Vh^###$#####wh^######"},
{445.93,591.38,5.95,2.41,"fJ*tK75K%(g5.6FY&0Cn('C1_c&HI$$Z6gM8pb####UM&n[Tu,%:o(pV1<g6AcKoZ&A-%$$BB@,/6$@%(o_T=$'2f0XV'n}CW[)>:0}B8@##G~T<d(:5#+e%YV;^P#0@%mA1*6%2[(e'*H^7a/1B.+c~/<w#f[T0,#E>#mK',_5D>#&##xP#7Q%[P#5##26'"},
{56.21,641.99,4.96,4.60,"CS&Ga6J%,l5%@6:ed+I>#c>$V8'+^5iY#pH)A/*DZ&.##O~/w//gm&tt89H&/:T4Q$;K/<u#Gj;:;5IYEml%D:T2$(6u#,6$wd+vA%[|BBv&b6TRH#fw-2T$WV?&##nW0[bCI|E###nG#|)2>6(Q##6k?kB-^uQ.##*l#-{&vT7###1##q@*2$(###$##rH&"},
{539.30,81.32,5.95,5.51,"(g4;c%###U##l&~{k####J%#Q&~######:)##92L5$&##B$#AB6######Bf$.'~###(##10*`)~###.,#*-#NJ/EJ*M-)0##mq?######if$0kJ###2##gC.1&~###I>#2y#AS0`Q$h_2QR)cx0l-+$##n-#$AV5?'###~?#m%~######49#^Q(7##-$%D7("},
{650.03,91.95,6.41,1.86,"v#%f^4N@&04I=92{Y$6,#9JSn12eY#0##<JS7@)###T##vHS.K0mo2$?#VD:>ISM5$N?$)r-B#;76'4;-'J,)%E|#'c5$46&?h-nA3v5#Rp3CKS{c(o>#V~(6):E[)K{(j82ME3PI)n['l-*8%++6&q~%HV6di84V<qG$^['3&(pIS~#%r,%p]4mT-#6&p$$"},
{153.46,106.71,6.59,0.15,";h#)3B-u#.,#*;%X'895#(##5[$'<5xK5b$)c-((V0K92dRJiu8xR.4,#t%'SMTn-+.,#aI$do2RQ',Q$<t1l5%fY#?l#`Y@5S./,#6~*{{0GJT/,#;$'`m#mHT######$w$gG$.,####4?$*6'--#@|A]y1YHT)##Gm);h#,z;###%##?&(E,$###%##SU5"},
{455.46,125.86,6.16,3.18,";Q%pv'$R'jf17o]%m$8@(S])Ap]*R)L$(`y$4Q%=&-?{8I6'iR*z.OK#%^G#1p]he+6#$H$#os].+HfY#H##g>#;0Z5?&OG#)U8R$(U6)k##?o]fY####b%#F-K_Q(###4##JH$n-+######U1=###tu&P?#$o]######L%#(_9;Z%95#1##Xd'Y$*.,#'##"},
{81.37,153.70,6.40,1.49,"s-T###1Z#Bm&_n.ll$[p.,w*DcB#?%H,$uY#ve-.,#######9-TJ,#iH)W]#[A.~E,M-TN,#[2Tn/.R6)<##:sB.,####/##ex0vM9,d)Kw#E[*`D)d-Tf@&(.TO#$9.,';%XD@95####t###6#XW@[G#Ae-^Z%8K-_Z'O,7K'6c6*QG#$)&d$)tc(###?##"},
{81.37,153.70,6.40,4.90,"###N##i@-^{@,Q%b_$liC;/14/1W'*f#&qz.W5$xl%^?';4C###H##*5Kpb#6K1C9%|SW?u#CTWVZ$od+%]&(.,p7#B[Tw$).,####2W<:5#X6)2,#,XWfS.>SW%l#8&*;V/>e-f#$7TW'########Bd&OG#.,####jX1Ll%?e)ZR+`g/)v&<$$+L6$q6.,#"},
{535.01,163.96,7.04,0.83,"V-$[/00v(VG#qh54m%lH(B~.qrU###'##vH(lK0###&##TG#]o0Q?(j>${G$HnUjY#Q5#IyPUnP###)##qrUzV=######.n%d7.%##e/)5f-&pU###8,#P)+[uM######R26?u$######~?%)/*)f*S&0Q@'.rU=Q&&##c{(yL;######N$#],%######&##"},
{186.05,321.65,5.88,4.06,"S{3l0995#@,#MV*RdQ4%)<[(([)LhQ+%,w@(FI,;g%gS248$ru#~3EbS%a&41V3/dQXZ$^[)},'V{7dP#p$FN5$(c#4##}23>06SG#MAEih9(dQ95#T5?oc$BNC######$8)N5$###'##d?'|C=p])~dQ>I$~q>86#gQO([#Rq=######K##m>%######%##"},
{300.70,331.57,6.19,0.98,"Dl$&R%]-){H%f>$)/,3u#o&J6-%-$(|@$-vNGu$.,#]-%yZ(7#M-?#YS1>w#2GHp0+mz=q31P>C_x,MX=?K1(l#Dx)NbEnP$CIGY#$CL9&##L8GG]0@D8Mo//19}{-7_;uB-[G#C<._rA95#~$%c6((A--d(t:'aM?c5$[n,(?#j6?cd+T#%0##97'T/.I#%"},
{300.70,331.57,6.19,4.49,"WJ/qb#'##8A)1@*Pm&5m$m[K5m%&m&kg$Y3F2U1###KQ#IZ$H6NI#%<5#iJ$l+JB0-5q8bs3XE=9S,uG:mA2&C6+##{b<?#$86N95#;5#NL+RmIW7,;t;iw+2OE^;+6O@to-7e-=o#1bHT5#f-'-c$.,#66&8&&`kINc$b5%M>#{vAUc%=n+7?'F[#<l$Ye&"},
{145.27,349.87,5.47,1.18,"3w-6[#$m(]&$E-(X$'ic'm#;Hd&ZH'`o(/FE95#tG#-M4PQ'@vOF,#*n-@7#usDS@&$vO8J&lHM;.'Xr?L-%b5%Uw#euON6%)0Jsb#qJ0.,#&c8r'9:z3+7*`%/w=5HT2x~+5h<S&$BL:df%%@%pb#2o'TH(oK%}rD?##`..ru$?5AeQ'eu%WwOG6$x#'I##"},
{145.27,349.87,5.47,4.35,"NH'9##N.M<-$[H':?&/?$S>BA##`..jB%dNEX/(3-(n-%eY#np:Mx$tV@J&$QK3$f+`%/>Y5@z36@*zk7p'9rJ0###r/Isb#o#Pc6%W,%7n#53AG-%]6M(%'6$Pw7&hsCf[&+n-D7#Q$PD,##;5GZ&.,#,Q#FB)F|CHd&L?'ul',H:=$(l6'$m(`x#I..?d#"},
{265.50,353.37,6.28,4.07,"6S(1=E+m'(Z$Zp4[u%~R&p`?;T3}b#}>%Y?;Jp76$$|&5H')Vw,~u%CH#TM<udSDu$&##iw)Zg3J(6yY$X$&?z0v4IWv)1w%VT0_5%~,#J^1}p9d?(s-$FyLP$'IkA(D'm]6en*.aCTu#=//;`?I>#Q5$'<*%f/811U/-@&E=/-,w+.pIPp7LNB95#~:3iH&"},
{252.49,368.60,6.31,4.04,"?U7lI+H(;XZ&,##Ci9BPFs5&j[)K;<{u%K&(/A/f>$W5$Mk3?K+R.S'I)^,%Ds;.o0dp.Gi>@'6G,$[Q#>2SQp7k6$q'9j*/}Q'oH(67$*sB[.SbG$-##M^0fbH~f0D>#7T+Uq0SjB3S._7&jf/BH'`,#'T/`L:#H%WQ#3JE-e+eX@5^&l'86S(raEm#$Tm)"},
{218.02,393.55,7.33,2.32,".bHi$)A6'Yp.B91@M=)e*?3By5%r%++ROr&2Q&0U~,hl&Bn'W$)1^,+PD5(7oQITXB~H(DQ#t7+9p+YV=OQ%6%*q.)+e+vA0y#%n6(cxH2w+u(9A7)}0,'A,x]2{h,y^9%x,%dPKZ%q5$e<7'm$$C0|-D&g66~.%$$K`7VFF'[)D5#lA+~dPIcPfP#yl&0+2"},
{218.02,393.55,7.33,6.04,"-R$TvCZI*F>#H&+~P=m~-jG$7{*-ZM?l$vd+e]&IkGI@+N7,Qv&b{0s[M=#$QsAl01k-'-/,}&.Bg6e<7an,G$Bh_7Rv(Rv%%6&*p,s~ILv(JQMuG$Ke(vx.k]0iG$`0-(==G6(NR&)=::38V93PL1uu&i5$I&.mI'i81+['e+8*g/66'Qm&*$&B?=Zd)6?%"},
{413.79,397.18,6.83,3.11,"t6$=8IU]5o5%SSB0_9x5&D$'}>7du&######^Q%+?&$##/,#JqUg95JQ&5z*sbG(##H%$E*;qpU.,#/##vc$Mw,oP$###DQ#tnU<Q&;,#*+/0&0tP$pn#YE:#nU:5#]5#Z]*Te+.6'###i>#rqUT6)1,#<R$Jz/y#'sY#c5$iU8rb#1,#qG#y-&*6'###)##"},
{85.89,403.54,6.05,4.47,"a.*dG$$##a,%S>#VI+M5#Cu$7@)EH'&##Ku#Il%###&##[-(Mf0tG$wl'x5#Av%i<8ED=,l#*TNHo1&Z$;Z#d7.2,#D>#a,$ku7nS.*m'5,#Jd&yuG/}@em(zoR=w-Mc$A9(z%/xP#ju&kG#{$E#g0eU3Yc=X[$hTVG2&VTVc#7Yp9ch$R82az,k&10Q%$##"},
{490.85,402.02,7.76,1.22,"+r%=2?OG####W:%;y44Z$em*yu$C.+=Z%Xs?3.,5,#:#$O_0mB&Hd*%Z$O5$yE,}jJ%##|P$K*Wpm,###`d'A.-###$##ec$M>#qb#{&$e|G'QKRu%#I#/K+C&WZP####s8#*&0######2##fY#(##Y{)>EBY/3.,#o]$f@H>RSOG####xT%H~.@5#95#(##"},
{186.96,414.10,6.20,2.53,"w>%uG$-14II(5[)Vm%+S.N[&-['uR&:(9k@+O@QrP$TG#jw%oc%0x([p84,#i;9m7+S6(891mZ()##yv'.AQ_@Qyb#{d)'23:?%fJ(&M;j?(a@QIQ%xc'XB)LsF:5#%$$d~F^/+jG$z),V?QD,$H>#]Z%1RPrbL0,#q-()=9=@Qql%6I(Ra08H$w'74h(x<D"},
{282.88,434.06,6.61,1.03,"he+=V<Q>#?v')B*rvPkY#^l%h$?ZNE###t5$Yp4######.##.R'[n.uJ-'D7kT3l%-kc&Lb<:wPg,&###?#95#L95####E##ZZ&Ml$axP<x0E6<'^6203p-&=yPqQ)###q%&2r?V,%###N##oQ&Gi:b=Bm,&|T2Wv)oT-QgPuOHpb####1<,|$+Ac%.,#;##"},
{282.88,434.06,6.61,2.06,"4##M36`5%###QB%H6O######5&$`J2######?##EQ%E$'###QZ%fs1?rAT>#(9OJ*=v5&[6#:O1G)>###=##`9%9w,D>####'d'ST(<U8($&K<:%g4W-&~.+/`1Y9O[c%=-&1:O=7,###/?#QV=$e&@H&[T*zn0au%=B$/7OkR,_@(@:,17Os5F-94/u#Se&"},
{229.78,434.45,6.95,0.85,"W_:JZ&nG#/}Ax}J6A.Wr@1B/)%(bw*YIQpc'?H=^H(.Z$&8(fB+{1=,c#l%,+/+kJQmT2su%paA<D;hY@k>IHC7eG$$-$(MQJ/,C%,bf1ty3SA)P~.lU-H19vHQ4Z%FZ#jj87^6j,&0,#vn,Z,#]/.UIQ,$(;X4,g6o#%co+IOD]v(I>#}@%wT3yH)###,##"},
{229.78,434.45,6.95,5.54,"wA.aRQO?(`,#2n+&p2K^3xe+*M<-c#7S+%L+5H&tG#V4C|U8$'5Mc$70,_a5q..2q-4z8eJ)jSQbK2$$&$&+de.${=w&/G5Lr-)>5#&A%$FC(Q%IZ#)t8}QQ&{;n>BP>Hus@Ye*YRQoc'em'`6)###,##`B3Il%0,#xJ.6p5qP$NH$4VQ&z6)v'e,%(A(t,="},
{79.57,454.03,6.05,3.61,"eY#$##.,#lv&_>$i,%RG#qc%bZ&+c$yb#h5$.,#&##XT-D,$M.$%q8?m'S#$u6%k{3]e/^5$y4?Ze.G#$gc$&-'###tl#dH(w?6l2=W?`m>%VfK'7({J.3z5:r?(l#Du#t`1L-)###fl#o#%3IR[V#j;`v7$}9`gP#Hl%Cp&66&b%*Y>$BR'AR%h-*G,$1,#"},
{457.47,485.72,6.96,0.90,"zx5=5####C$#5?&Co/###1##yl#IS/;5#0,#J>#(##Cl$PG#?6(4?#t$,R-#_7(QyJu6TO#$g;T2dL}#'{u#Un,c#&Nl%0,#.,#ZQ#i&4,##O$)Nw$g;T~V868T.j<6~)=G55q*K7TZP#&#####u##im+95####h##t);$v'wb#c_-h-NB6'r<1IkA}J1iv("},
{457.47,485.72,6.96,4.34,"%%Fl06ol%wR'9,@U%.D>#/?#%E:fc'###A##>$)$##/,#~5#G/2=5#2](MWU5(2.,;<TU6L.MWUNq8Kv(Um$6]2I##o#'J,#_>$###C7({g8R6(d#$MWU_ZKfRU:u#$/(MLL][,J##]S1LZ#]5$~P#TG#$##6u#fY#qZ#@%-95#(##*c#Pf/###'##hS0VG#"},
{418.13,527.36,5.99,3.13,"'-E~J-Hv&Y5$GMXOG####[$#2p6######H$#T,%######T?$98AgO=pA1ko*=OXZP####[@#zuO######}##]c&&##95#~##Y..]G#wT)=OX<JX######@**b6I######4##m-&<#$######P_4u07ZH$oC6Tg8######iC-l:6######k#$YR'D>####(##"},
{148.40,546.51,6.95,0.91,"P#$RU16#$###8B&[;>PG#/,#x,$AZ$vG%###[5#Du$######9?%f}3aISmc&bMS#ZDC%,f8(:27w<Apb#,##H$#195######|k#|,#OnG*03k{ABB*4K.BG9|x*bMSYS1U>#js/};A[P#F5####D##B.*jZ(###T?#aj:x6+_c%}w%>cD&}E3k><$(k5$u:3"},
{120.79,591.64,6.55,2.63,"$##1K&obN###sJ1Ix-xT5wu#zSS5Z%%##CH#[v%l,&%v%gY#38(8hORg8'l#>USj^9tZ'#&%&`>######oD-`H&+u#/,#ld(a~-jo)|[@,vOTTSxG%A##nu9l4K######.U*#########K80pm'g-*hT&&W?/h<###C##hl9k8395####M.%1,#-u#%##ul'"},
{444.48,591.19,6.55,2.43,"`n(K%.Ux&Me.ME;.I)5x(qK2(l#wG#}j/.tC######D_#GIVFZ%BJ(sE5z/3}tI9?&*-%rY=yc(Ju#H.(|MV^#%/@)gM&z,Mi$*k'-^)BI,#VJV2$(E>#47%E2<=5#x$%MS.BZ%f,%m.&wL97q7NT1l[,ER#r6T###;5#ip&sU7###%##hu#0Q%OG#(##7Z%"},
{434.42,619.66,7.01,4.51,"M:<######'##ehjQG####T##'.RA[*D>#,##/H$F,$QG####qrD######6###ijiY#zb#>-#XeU&7*4d&w?&yw.Xu#(c$[5#&aF######D$#hjjVZ'%##y%$x$?TK6M5#Ev&oR(>r2(c$>,#YYN######Y$#&ijsG$QG#CT%`90Ig0;$(;7)z6'RU1J?&^u%"},
{130.91,639.13,6.36,2.26,"I##^2;QG####X6#PK4xn%eY#0##:e(lG-z,')###l#l<$jI.V%%Z[GVZ'###>~E,93Bm'6l#8~+*'-`45]I,e7+7H&MO&'A/U&0$'07/1Tn%G)feu&Q5#<x#=J.'l#cD0_m*X#%^%+V3>d#&p'8:#$tP$9|0C;<######_9.|v'xP$_@-nY#$##r5$3mK###"},
{76.07,673.72,6.33,4.73,"ZFE###g.+<v$uzWrl%hc&C5#p].A~)Ay.H,$^X>mQ'+['wP$K@E[>$Z%++c#*xWQG#pQ&.~%0U7@6##O>q['5wW:##a?';8(=c5qb#:5#)##W|W95#L#$Qc#Z2@&##&D.zK0*WB###tP#@U(DB/#########/{W######?,#0dR4,#hl#0~*`I-0,#S>#Mn%"},
{613.23,73.43,6.43,0.65,"$##R%',d)M.&$##VC#dHR~P#gd,iH#-6I914-%-###q['o_6kw'CZEe7.T>#TQ&h_'^IR0u#S<C$m%em']<89m)7##j?(|oIG'3)I)IR&pa1=Q&pP#+$;CS.5JRfY#Gf&a@(jW:K%*R82/%'|I(,d)?##e/)Ac#;c%>i)TQ'&e(eY#B_%3V=1]*Em)Bu#Y15"},
{613.23,73.43,6.43,3.14,"6d'w~0}P$EC*Fl#=J'YGHBv&gTP(-%uI.###q::####(.0,#;&&QmRpb#`u#>d&?m<b98E5#=nRCI%TH(>##@PG{m+@6'$##Ep-/*<.,#z5#If.Ww)Ud)^<7J#M%6%oZ&8,7*f,@eO1~+y5%#%+7R'$##9I>}Q).,####Q.A###>5#U.&CGAjY#RH%I2-?;="},
{127.45,140.53,7.19,0.76,".;3_C<###0Z#R6;Nc&&##;c$F7C######+##h[,######n%#hh?###1##<i)I+C###=A&y90i0T###S>#wG#h.T######.&#0Z%-[&l,%x+4`%/f-'&`5zB,X-T###8l#='%6&K######5##E?&q.BC,$:%(AZK4z5X,%N~#QEE######>/#l7*95####0##"},
{80.78,153.94,6.89,1.52,"7[R###Y>#NI(k&2y#%AA*RA/V$Dx>%ub#bP#18,95#######MZRG,#)$(:B#AT0OW+^ZRyG#]_Rr~,yu'?##|WB95####+##cB39;8'$(*J#w6*r_)iZR:J&2[RrP#^$*vL$oDAD>####e##xu#ojFE,$+Q$-e(:V9y5&G{,$M7Gq=:5#<h$II(|p;###9##"},
{80.78,153.94,6.89,4.92,"###:##(%+88VfG$}0$qNC@}Hb[,G~%8m(e;6Ql%h5#W-(m_:###B##e=H_5%+&/gT$A8V1l#)8VPH$vI-i8'p?*K%#A7V>@)######A{:OG#-$()##F<V<8.K7V/c$5]+hC/PR+&-%b8V'########Y-&{k#######D=/5Z%2['A.,,;//$'+Z#;^4SU56#$"},
{336.74,155.96,7.96,1.72,"K8,+u#>5#WG#_q,||FL,$@@)2##?a8-g59S/###=##1-'jR-:v(95#M>#|5$GbJzu'M,#m%B8?'O?&|8)ReU6##s@+aS)?dU*e-###$##-c#0gU###bG#=I&:|Ae5#LgUy%/E.'[:4heUuG%Xc&######*##DgU######$##h.R###r5H###x6*$##KhU<-'"},
{249.07,162.70,7.14,1.41,"un1######F##/]_###6#$G$#/<D0##.TZUl$g,&-##7_Uu)@[U;######E##J]_K>#bG$C$#jOIHQ$LC8e6'~6)([$xL/5cMlz>######>##l^_|Q).,#W##3,ES04r5&&##Yn.;o$;V>1,#?_=######.##m]_Ll%###8##-7Sx#'r5&6##U7.k##aYN)##"},
{164.71,196.47,7.65,0.19,"WT6###/,#L'%5jF######7F0wl'###*##L&M4u#D>#7,#,ZE()?######VQ#(&Wt5&###k:()e+EN@$##T0I3##1PG$##EGD-B1.,####2,#QNUU3F###e5#be&z'W###e,$}$#V&W###/,#)c#;5#D>#95#dm#?A1.,####V##%=E######7##uh=######"},
{164.71,196.47,7.65,4.91,";5#D>#95#)c#.,####2,##91######WQ#JD?###/,#?y$>07?A1.,####dm#`<F###Z,#QNUt5&###]1(*&W######7F05jFv3E######U##'(W###Z#$`e&^`@$##|]H)e+###*##tSLlc'42=######6##Z&W###/,#y$#n}G$##IYC3##D>#7,#_5E4u#"},
{491.82,308.19,7.92,4.43,"sTLC&32,#-##|1Rkc'###2##AF4k>%###$##~d)######-##,-R2,#@?#$i07/RjV4&90Y7'|1RHe-8c$$d$ry595####*##1R*F>#a6#|1RN@,(Z#4r+R0RG-RUG#86$DeB@J0######N?#U$(%v&O>#$f,L[*###:##~p1;c%###'##q&,mP$######Yl#"},
{443.72,351.81,7.51,3.41,"hKY;T1,u#0##^LY(c$###bd#{@+r5&'##z@+$/-+u####=u#)gUN>#JZ&N,#ZOY(c$###U##2F8VZ'###R#$h]1OG####T5#4K3,##],%,7;@JY######EE('EB######e:-~n/###%##K~)ZR(%J-*H%s.CW|DL5$###wi)s/5###)##Ox%{k####:##f6("},
{229.67,434.34,6.96,1.00,"yYJ/m(?5#-10vU:Z@*o[OB]0KA)K@)d%POI*A]M4Z%@5#h8'WK-}:<tb#SJ+e.+i'PM04'Q$E+Fbo0227t?E{^;.,#7,#>d>gx-Wn.Hd'^p2?.(%K0m/)]B4B$Po>%jG#&i0Z`=}>&###=v$7,#=@(Y%Pp@.hL.D;=]@+Dw)B%P??'D>#pm$4C2{Q).,#)##"},
{229.67,434.34,6.96,5.46,"~/,3IM)n,dl#3S+%:4<K1;w(Gf3,l#Pw*oq,j#&jP#w4C7T0n$+O#$|y-'6AYI+#p*o94Kx)f~Q{w-_u$n%*Bx/V&2x'.]}H=Q&###/~$4PG|Y$sG#aY6fZQer;C>Eg4C$V8&K(U[QbH(2$&3l$###*##G^2D>####%*7NT3bP#Av%@_Qn:9Vl$K$)d.(<QC"},
{462.07,504.76,7.33,0.77,"#H%m11Nc&E##^^)kPIeY####bI(E,$1c$#########3u#UG#oY#C)+z[QYZ&w^Q[*=N92FU+$vFK82+u#:##N,#@$)ZP#$#####h##`-JK6(nx2Cb<TV8{9/y30,~Q]#&2##`_'nS2#########1##'L.SH($##qR%@_Q*R*++:rS-0g.A2:D42MK65,#Q#$"},
{462.07,504.76,7.33,5.71,"5##KK05?'###S,#6'O%w+ZP#=+2)%ROc%Jc$eg1mP$)##l#&7?#'oEk-+###lU0BW9YV4tq<%GCK~-Mw'p(R/)=E>####w~%u.*Kk:0kFGJ+X'RbI-D5#03,lDA}k#'##O{5F>#/Z$:5#EQ$:x2O,#b/0;)2ZXI###&##K;(YS095####X##0,#K>#H>#.,#"},
{302.79,506.99,7.08,0.73,",w,L,$###Z##b,$V03######NI$A81/,####)c#1,#W5$.,##?%%~&Vv*-##sv&tLOcaGJ#$+iS[-Nj5%{u#FK/ZZ&fY#%##1,#v$%zK7###/v'he%+iS;)9keSiV9?v%'G3]M-IdS###*#####^Z#yv*D>####OH#{q7oc'Du#[b42`=RH&)u5r#GtQ(#-%"},
{302.79,506.99,7.08,4.48,"G,=LED/,#G,#mL18o1D>#>,#ah4lZ(###;##n,&###:u#i5$s$TwP$_v'EI>DU0G5>J&TA~)O)T6h9yl&Zm$lf4###&Q$:6%W>$###w$&WKPgl&n5$$9Gl.Sk$TJ,$?w%SzO/S/###4l#dv$B,$###2##*$'U5$OG#K$#']2:5#.,#V,#cU6.,####%Z#/Z$"},
{142.09,528.38,7.44,0.66,"D5#>*4%Q%1,#],$IT/+u#/,#FI%?w-$##Lu$K5#gP#&l#fG$)##bS*%m(###?I&}^OciBX,$EMRM[O[,%yu#hJ,^6(.,#$#####eI$1T3###bQ'#]%6MRbz6;JRp*@qH&gj1D;,lHR###'#####e##3%+mP$.,#j,#9`7}c(d>$*55%4Awc'9?:`OAnc&r-("},
{142.09,528.38,7.44,4.61,"^&(],N######^/)>g6######[q4>6(###+##%[)######G,#ZTTAe-bc%cB)PM2B-DMmJ~d'0WTpU9JH&am%d97######'m#0&2###_H$sqQf6)Hu#M[=STTURT)l#q[%0WTo]6###*##Hv$eY####'##:]/~l&###g##oo3QG#J5#k5$Ng5~5$/6%7c$)l#"},
{416.69,527.10,6.96,3.10,"Fa;g7*r804,#pqW+u####a##;1:######4$#6#$###%##X$'8,=IF82kEhx-PsWmP$###7n#S6R######+$#-Q%?5#QG#2Q#.%,rZ&{|1)rWFnW######ej-DdK######5##L6%xP$.,####sOA-D8PH%C:-S07######A)-wB4######U6%g$(95####,##"},
{451.03,529.74,7.64,6.19,"###z-&JA0###a>#+r/bC8:v(~(6n%,*/,IeBHp6aP#.,#IA$###>@#/[R###FJ(]TIisE{5%]_R005J>#h.'Tg4TG####[5####s##,~R###L7,8d%U_Ry'27~R#Z$kP#;;*Qz7IH'###-#####H##_2:pb#.$',##$%FTc%(/O5l$Q#$Yc$Z/-%T/)c$?5#"},
{148.73,546.82,7.08,0.72,"J5#Ti43l$###Fq(FbG95#$##hR&NZ%fY#.,#'Z#V5$.,####{Y#P3.O@Q}c'wCQ#O=Ie,ny*pt@MD?###8##?$#r[-######.,#X##Y5Bv[,EK40)1,T.w_2=i,7CQ#?&3##v30vA4###$#####.##|w)[?)###(n#KeGSH(=B/S9,&V29E>78Go#'&##R[&"},
{148.73,546.82,7.08,5.34,".##eo+=6(###c##DO=Sd*####E+L@PlH([G#+|4uG%###:,#B8%^JJ$m(###j|7D8-k11$&0li:yn-;{-#CP?-MD>#$##QD+9^/~o,+)9G{7RAPmP$7##gN.$L8###N##+RH;Q%PG#.,#[A)j>%###96$s/J_@.######'F02d)######;-#O,$95####3l#"},
{142.84,576.68,7.75,0.12,"###'e%D^5{k#EA-CU3^L2qP?$[ExR.###t0'Q$'xG%###%#####}v$GnQ95#YA)<{S&3=G?&C{S~WCUG#SH$(`7L5$###$#####]##IeK###]-)1d%nJEX;:}wSuP$Ll#J|0%j>nP$###*#####k##aA0bG$N,$aG#mP?pQ(7k@~P#YQ%?%(y2</,####F,#"},
{405.35,623.88,7.64,4.49,"u]7######9##;`n######b$#(0a:5#2u#m,#/,#0,#`H'<5#$3D######;##F`n######U$#cy`*c$.,#?##=,#vb#Ju$###TXI######C##X`neY####Z-#`[Nll'$##Kl#0d'?Z$OG#(##H`C######y##haneY####~-#*F:z#'$##_u#RQ&-'-ZP#+##"},
{633.53,94.05,9.22,1.61,"Mc%*X7TJ,CI*=80qQ([m#UdGgx3###&&#<TUv5&###T(#6RU$^4.c#w)498,aRU8#$0?#KD+6-LOu$o[$9f-6q5PG#%$#'g63-'VG#2?8j7-+l;7S/o?&>U45{6`1:fu$D]*nWAfY####:-#pv)4J'F:4Kv(x7-Sc%p,$'xS|u&tq4:5#g+>3](^jE###B?#"},
{75.62,132.03,9.15,1.89,"en(OU4mP$-##ih5]-*=l$pI)4M(A7+56'*H%V]((c$######4`8AH%o-)bG#tnRaG#%A/du#m28RJ+GmRS#$}oRb5%iY#M,##x/%##KU&]s?GOFFT/xZ(^w&JB4sz5oU:}?$3tFR[PSG#`?#/##aRI%A'`%/;5#uqRbQ''H%.Q$|pR.,#-##T,$opR###'##"},
{75.62,132.03,9.15,5.05,"###&##;u#3(R.,#+##zG$O(RA$(Ul%F>#p(R2A'he0;,#plDSG#U?#pjEQ%R1z9g6$304yz5mQ(Xw&eFHFf.*C&Lj?&f0$##iY#W5#('Rl>%U$RR#$g)8VS+Te/[u#.&RaG#2[(S>#~)8ml%######>K(3l$??'b,%@_'?%,Hu$e@)7M4#I*W>$.##Ye(<:5"},
{440.25,150.22,9.13,3.35,"k@(xfO7?'VG#6+^+f2###H$$w@.###1,#oy0o,&/,#0c#_D8c]4.f+}n.0<=['^Ll%sY#;s)zL;r?'_f0M:/dT0|81@H&|?'BC47_7/l#aI+`,^6/0_P#FH$6V2j.O#I)9,#>Z$)XAgY#%##gU9?A->Q&*##n(^Il%###+##I^/]7.######nd*`5%###1##"},
{357.58,168.22,9.09,3.33,"95####&##Nl${70###3,#*d$3C:;##<J,A|4xY#G~*[e,:4CW>$###D,#xkC8iA###&9(|WB=7V[>$YR'TX:pu&pn.x?&-;VM_7|d+ZQ&U_4-U3(D-0TRRx0F<V;@T'R)sZ$|6'c:V###4I(7d&FuI9l$$##^#$45<SH('##%7$0:8+u#5,#]>#kq<######"},
{264.28,200.74,8.86,4.57,"6&195####.##3%Y######G$#xsI'##y]5IH####)##3vLmP$+;>.,####xG#*&YFl%###3/#i(;7.,@n*_~&%l#S>#4'Y7v(bB7######dG#P)Y&sBA,$9##hUKPvSl>%F5#Mx0eu#7%Y;5#xI.######0##zy8$I'Cn,{,%K&YUQ&}l'G-$YL;h##1%Y+##"},
{513.32,320.73,8.75,5.78,"$##@dFl>%###^U*GJUc5%P>#wr1qQ)###)##H6%D>#######C@#YKUYR+|k#4&MRz8C%+[h+@KU0Z%###C.#ZA.95####*##Nh-JQOB,$r[$3NU=R+0,#]9$)5BhI-uP$y5#j;8A,$###%##=$)(##E>#VE+&NB$##gP#(N(M81bP#8~*_v%+N@cP#H>#wG#"},
{251.38,406.89,9.39,3.96,"&##1V-|w/E>#z,%U)09|Cg-)$n)k$)U7(CT1UM0TEC?R&u#&n#&*a/T6OnY#~6OUy1_S.=(2@R*MB+t)<L;;980bA-U/*se)lG$zL*AtA.n-|+?>sB<7(td+b;9btCTU7LL26w,4m'EQ%EnC###|P#lG9lC=t?(H[(xy/N19TGL=Q&q>#$H:hC50:91H$#'*"},
{300.47,424.21,8.27,6.15,"/-&Gw'Qk=#~+ka@h?'xR*C6&ph<,m&N?&1?9CNB95####Rg$M?'*,54C5)7(OYB2Z#US)v'5PTS###+##%V)-x1######nm#5%-W$'VS&-zOpv+###m$#nVS}PO###3##wv>*R*######`%(1u#ZP#c##.mKyG$%Q%###ilBZ5$ZP####OU..,####^P#)R'"},
{189.20,446.98,9.03,2.68,"NI$&z66#$###[>#z_6r5&%##`x+J.-OG#Z,#X?&Y,$RR*P,$rS)6}G[v'uu%2f)amJT[+*c#xKR3w,gY#5K'-n-2,#wZ#$KR3v&dm'e'3kZ&i&2rH&[8/EGAJKRQ5$'Q#lQ=,IFL]0FR'4W8<~+oI'#8//.*%IRh>$J,$'k4aV8A,$C]#:KRgf')$L>V,wm,"},
{299.25,553.35,8.78,6.27,"$##:%'n[,###Ec$M_.'*>qw-|$M:e,gQ&zs0i&3iY#95#^#####_Q#9vN###v.(&1ObPKqc&>2Sxg:bP#zR&Zq8D>####.#####0m#m-S###}R-xR&80S*L/e.S],%5l#Er)u2<>u$###-##OG#^,#IaA7#$xm,[##~|9-U58.S###-##m~)S05eY####Sc#"},
{453.33,566.55,10.35,2.56,"###*##q.)c8G+v&`P#AL+T<BUd*###'(#o#M######)(#&6Q###:##6p,*'1'I)~v#2%J,w+j6Q-c#`c%[C0mc'###>%#hYJ###FJ#frA6#$_p6W*/zOJPQ$`6Q0Z$xY$=A#x@.95#+##uQ&Rl$?<6(E?:Q%if4r['*M59><d`D95####Ah%A,$###-##$[&"},
{514.84,579.11,9.37,3.57,"<R+######N##>QR###%##:s#u%/ku&@5#*a+YB0yu'/,#Q%([:<######/##4ze###1##Op(_T6]>$bc#QJO`Z%?Q&{n%WC6<p6#########/|e###%##d,$-nW###Hl#kC5A,$###N0)B/.F.-#########j{e######%##Kze######|c&D>####PH$fl&"},
{440.46,590.53,8.24,2.53,";6&jY#6J(Z07*.+###u'#o3E######Y(#k$V######Y(#S$VLQ&2Q#V|4G.,EEBZG#h#$i)26?'###}%#n&V######Y(#r$VV$)>K&[OIP>#w%V[,$U,%)d#oq<###Z##bn-sb#.,#T##J^6@18d(/E%-|Q%q$V]P#95#}/#NA1######Xv#<#$O5$4,#'l#"},
{140.60,616.31,9.25,2.46,".##4d'3Z$95#,##NB.S~%fH)K,#eI,9N)AM?###&##PD&O#Q$##}R#qS1.,#Av%$68{&6'##Cb?7B3n6'Zl%Nl%9,#2E/'f2S>#ya0>EC;5#6F@Ym>sI.xH#M7U?#$SG#4f$.[(J>#xB2h>$5Q$FU-xoR~I,vQQ,['%w&C#6d]6###$##N(+[P#$##z?'_$*"},
{124.05,635.74,8.87,2.33,"+##%],gZ$OG#I##l..[j(<R+$##C5#nj)0p7######e)#yS3~c$2G4<R+&##4=<gT4Q-'Z,$Uc&cG#'Z0p7/###N,#3G-3'7#M7DE7s#'5I#?$R?#$K>#U&%@@+S>#8`;8c$'##y[&6IJ{k#0*CKu$:c#}N->R+###(##;(.E,$*##dJ.ZQ'%##zG#mGCpP$"},
{475.60,124.83,11.03,3.18,"<,#$8@0Z%###@?&b+2K94s[+b8XO7,>u#:3*Ah<YG#.H%)V.Ec#}O;pl&OG#hL26E@Z5$]H'n:X.n-3,#u%&yD9CN?7H&eu#%##@u#zM7A,$W`Bxm*D7+g-$n8X0Z%###4%#;z2jd,###0########eO/2H&/&2###hv%zi-%cO######8C#sn/.,####hH#"},
{42.45,162.79,10.24,0.14,"dh=#########9OhF>#:5#-##qU6|u'G$%q~0uG#;$);H#qsG/D=#########uMh###.##c$&.r<T..|v#ti?oR(]w/kG#y{8#p5##########Oh###%##f,$8^b###h##VXCpb####/##|OBQI,#########%Oh#########{Mh######@5##H%######hP#"},
{309.83,371.69,9.77,4.76,"`%&>C40d'5Q%:j;CB.Wl%Rd&XAPhm({G$=z0}v+%%(Kw$G$P'##c,%^K.[6KMM;Pu$}J'T3=GKW7#$S>#n;.'(9)@%L<5;kD#######~#NJWF6(/,#wq$cJWnKWG@,J,#gY<rT/V.AF%.+Z$###/,#Lc#5?'###:,#EQ$f@..l#O$'#l#8v'sY#Ad'F,$L,$"},
{25.52,642.33,10.59,6.20,"L5$######$##q&a######&%#K'a=$)###v##6c#J[*95####z,'######(##3'a%l####u$#VKZgr;.,#C##Af#U07######8[*######'##n(aeu%D>#T###N;B|82-('##E1*;(;95#'##{#&#########0CE1o2.,####vS&y'a3l$###=$%A(a###$##"},
{27.76,676.50,11.19,5.90,"############oOK######$'#0%V`I,###Y(#:',pn0###>##@c%######&##{%VM#%###~'#m.T?X?3l$1%#t3/JA0.,#,##rQ'#########.*V:~/###;##;*+,%VA,$)##nN.^$V.,#4##1,##########NR%L5$######FR#o;A######9$#.%V######"},
{693.79,26.80,13.89,3.31,"/v(######2##K8[D>#)##(f#)xJu092c#rZ$yn+>^4Z>$0I$dR-######<##O9[&?&d,#I*)AM:/i<Ff$)@F(n+jAV6,#DE,GQ&######$##|=[po5(##%Q#jf*k:[<c#?$(TG#g:[###:u#############I5#.,#######&##.,###################"},
{662.91,129.13,12.97,3.36,"s$,######/##n1j######Z%#|1jM5$###($#<5#U$%9@,###Zf5######<##*2jZP####}%#ryc`x0.,#f##H,#uy3?^7###<h=######<##-4jFl%###'$#Dk<2i<(c$)##B~(-g19v'2,#S^:######Q##j2j3l$S5#[J%To*m'9*Q#^m'f7,O.,H>#=J%"},
{326.78,160.39,13.09,1.67,"`S0###O5#0$'/#D>m'HR)d$'iP#lH%c+Db-*###%##=R'_5%_d+###A,#^#$Q~V7Q%UH%XG7F6(/d$=nB,[P###|G#I25ac'+S.iY#E5#rY#n]V###(@)8v$R{?gc#d`VC~.###cu#]+=R%/CR)oc&0c$cP#W^V95#aP#4,#D7NjY#er/jK7&##@5#S`VnJ3"},
{182.38,222.00,10.89,4.89,"16'######i>#[-*######N@%CZ&###p,#woQ######9_$JtK_-*######3##NnZ######*@#GnZ###8[$gd?0Z%###ew>hU:]-'.,####(##boZ######%##OpZ###>6&9l#Dw-###3pL2u#<@%yY$.,#%##4GBW>$######RrZ.,#>5#A5#mg8yY$x0+bc&"},
{530.43,364.78,13.28,0.12,"3m$V{1qQ)###-2*/&1.,####y[#a-*pb####E,#;n#;?S###+e(mt4K'8(##VDSmn/95#b5#Bq3bQ'T,%l5$+6%l&$;?S###:&2pZ$W@SBJ*n?S}Y$&c#-D(W8/B%,`-)+['3l#6B$<?S=5#xY$&##2o,/fPk-+###W,#%%?.u#,Q#/_7;U795#X%#H?SF,$"},
{443.47,457.32,13.17,1.50,"DzN^DA(##d$(tLT4u#$##}>#t>$Ze-%Q#Tu%5##aI,r6$>%.0JTbG$e$#bLTKKTq,$-v'_I%xd&jP@O`:Kl%0R#tJT*I'+?&4-(###%A#2SRG*E$##LQ%hW.BQ&{b#ja0uI-`G#Sd(=C(+IT_P#.,#4##hK4Y>$###?5#Fw'.,####%-$Ad):5####I,#?05"},
{160.84,655.92,13.83,4.65,"Ng5#########ham######d?#[L;@6&###^J#<-&Up6###;##N:;#########;am###1,#c5#]FJ<%'MS+N@*W,$L)83~(UU90'6#########cam###%##X>#:8~###qH$?g5pb####[Z#IlL#S.#########Mam######(##Myb###9##mZ&'Q%D>#+##fQ%"},
{541.49,30.46,17.38,1.46,"############u[Y$##ZG##H#QaHT>#dT+0n(5B49#$c>#-]+############+~Y3u#0,#L##yZSP;2;$'mQ&|['T:T3,#D$'ZP##########[~YTH%eu%i>#*;<-S>cM7oG$}7+:9L<?%$?%g,&######%##u[Y&##5H$_w%-o2xG#:t.:<@CQ&Nu$^x$)]Y"},
{611.19,30.61,15.43,1.51,"############JbKC##Zc%lf,*e-Z$#hX3yr>bG$&##Kp#g0]############P/]+##&d%YZ$9{@J##IN/BT,.2>$##M,#s;3L5$#########l/]VG#:u#K5#)=G(@'gJ+Xw+|e.1O=9,#&x+#?&#########M0]ZZ%gY#(##'h9AI9,-'N>#o#$AjT95#%##"},
{58.65,62.29,15.83,1.23,"(c$######,##ovW&##95#Y'#6wWPq3;c%p%#NL,4^5OG#*##B%.######@##?xW1u#|5&i'#i|DQ&*RXB'A&k;<%R(d>$Dn$q$*######&##ZNUK97sb#9##~E*1wWQI)<H%>r-AV>6,#W?$&###########B-#P6)######D&#END######PM$N4K###%##"},
{95.76,90.75,16.26,1.51,"27W#########_7W0n)'Z$2,#,:,>o1]P####[G;+[)95####>7W$##.,#$##o7WOL,<&1A5#%O9TU5Ol%n>#|s2|H*###(##09W{k#/,#$##PTS&C7rp116&*^4Rz;c?&yv&Aa=Oe/###4##}o+*6'######t0'27W_P####(I#17W;5#/,#`8&|6W###%##"},
{577.33,89.37,14.27,5.75,"}B0NvG/l#fY#a&Tec'm-#p97)6&$l#(V#DlP/,#.,#_9#xz?`QPPR(Pl$:d%ieT]P#nY#1Z#N80[I*=S+.d(qP#0_7:T&281Co3###k6#HfT:dT&##xG$P41{d-7H$$z1l.-R>#[7*Iq1@H'######Z$#IeTxY$'##j>#pQ?95#m?#^J(Ie.###Tz#V2=3l$"},
{553.18,135.71,15.86,0.79,"n[+Xx3-u#uT#MH$gA2%I'pm)P-(I#%Uu#q9.{k####X#8wd)Sd$yRSw[-4##8<+RQS'##+l#^hQxd-&##|Z$$R)###4^-AQ$w,$<.+*I>v-*FdMbR-G,#'1(sRSpb####@(#)w,###'H$o>#######<=0i^9vQSRG#SU-bB,bQS######C.#-c$B,$.,#0##"},
{48.42,179.45,15.42,0.14,"4(:#########Jjg:u#J,$;u#nJ0,-%wH&8D<994qb#&##YD2:g7#########{hg######=w&4r?######SmDR82######W|,2J/#########Kjg######.Z#?qd######cx-i.,+u####;m%;?'#########Ajg#########=ieAu$######Lc#A6(######"},
{303.20,333.25,16.06,4.91,"o5H;f1,6$Jv'?yVk91P6'Ux,9h36%([H'v?(u{VmP$###:l#34@-'2%.'m81@zVv[*_c$n?B7(8k6&,X:LHJ<,Bqb#X,%Xl#nv+hY#}%$*xV0yV*o-B?%sRG3y1#E0~<BMd&jf29l$l>%=6####|6'O6%106k>$S:4[#%_@,fI'^S/.,#,##P[(y,&?#$nG$"},
{152.29,346.67,15.34,0.22,"|6%kS1TG#SG#,H#?i?$##hY#qG$1O;iY#)l#Yl%'9,.m%sl&<n+An&%6&Ku$aD>Wp1.d%DO7*/,AV5IwB0q7NB.>ASi&/J-'IU5>v#pd+%##rySS$&|R+Gn&?&.VV)>zSi6'[6'R%L43A,Q#]r<O?$~6)###,zSI[(mH(1$$z&/rT+*zS:d&CH&Fc#5QI(%'"},
{152.29,346.67,15.34,3.33,"z=F>7'rl&`u#g^TzQ&-o,<p*,[(f?$u^Tov(|Q)###yq:T?$~D>7Z#O$'WoQo^T^-'4/-_)*e.,`w&Q^TW$&YR+%##<C6M$$bx.{Q'fB/^RPW.C<z708,OV6$[%ba72i>^y1u,&@l$b%-/~&SH%EH&+H%xo,]P#wb#S,$^j=$##iY#DQ#hD?VG#_P##%&d~0"},
{142.83,455.38,15.25,3.01,"Yc%Ie,dG$8,#[,$tl&2~O|m,)R*/##n)8}wJ5w-$##Ml%;T*Kx,w}AD>#6,#0~O:[)'%&eX=Hc%r6$'6=pZOI,$f>#Y81ud+`NC:J,|,&+%$n[OS7)Ce+oM*6[(3;(>[O{Z&06'$6##M=.$$=`BpP#fZ$d35LZOj$%|p7h1)6Z%Yd#pZOT6&~#&f##q*HqZ%"},
{292.68,468.11,16.91,1.05,"rS$Nk5,J.Nl%r2+>6(,##;H&kQ'OG#>l#qI,$##:5#9g*O..fo)8q,0y2.d(wCQ$I)]5$Ql#[_3WZG(m(vb#&/#PiB`G#i,&a]3;['q[(6w(@?QDl#-.+%0$-x*UE8iBQ3c$c(&t?QP9,t?*7/1.u#bg'vr@JFF;5#o>$U~)ru&oY#)z*}7/8$&`#%GM'TXG"},
{406.81,534.26,14.56,3.07,"1.+:5#%d$1zR^*@/,#G>#%{R<-($##5u#]xR###c,#Og6C>M`d(`u$0yR,&/4zR2H%P.-NZ#n)B###Cu#[6&###?##;wRoP$xG%l,#'xRu[(3vR@5#+J,Q))W3D/,#D,$t5#H>#<##]HOZP#M>#N&0Fh4sS.9w(s?*1l#i(+17*[P#95#s5#ub#Yu#tn,bG$"},
{154.31,55.46,16.78,1.46,"k-+######(##o:i######z##Bf[k>%1,#V5#JPKD>#$##A,#X83######-##{:iPG####/$#/]VuI-###5##P[GDZ&######d09######1##B;iQZ%,u#f##]s?ZPB<Q&+##`c4Bn.###'##L07######)##g<ix>%Yl%]5#D'408-6&-{l%T29{,'###tG#"},
{657.39,134.91,16.08,3.35,"3w-######3##gye######b%#ryeU>#ZP#9$#rb#c]'=6(###F'8######A##}ye.,####y%#~ze2J,0Z%e##0,#<C'|-S###2C:######<##t{eL5$###}##|}?Kr?OG#)##6-%1~,lL1(##L@-######D##Zzepb#$##iJ%zB1tI.%##;7'SR*[&1D5#{5#"},
{524.37,188.35,17.52,1.02,"n#$g%X[P#aP#m6<9rA###4c#u:7###mv#Nl$95####C=.OG#Be+=;=CT*~%-S'X)R*&##vI$_%X###(##4-#|k#F>#:K'$##/H&###C2%=wVM&X###r~#%h1n.Vj>$.,#c,#=5#dS-cP####L5$;##$F0Hf0u@/&##o'(qtA)807Q#BR+`[%$##036_o4###"},
{124.46,652.04,21.83,4.67,"@u$######6##=<D######Ci+n[-###-##+tcQG#95#_##Rpcy6*######$##^qc######[-#ncQFJ/vb#Yf(@.'^B73,#o/38%-#########~qc###$##%Q#&IU,R(B-%Ci<dm(o6*lP#@bB?m)#########sqc######1##npc@Q&%##&/+zY#Ed*wb#jn-"},
{248.37,182.81,26.56,1.39,"W%[######0S#HU:###_Q%IC*bG$###La;N.+###&##a.+_>$f%[######y$#h%[&##an+ue%^?);##o([(-%A#$uu$Hx2mY#=([z,'###d##G*[;.)l$+@,#zv)oT')&[2,#&##4?$5h;eY#(@G1S/.,#/##G@;Bx0ZP####}u#9V0I..###3l#x,%t6)[u%"},
{248.37,182.81,26.56,4.48,"^d),?&3l#Hc$3w-###GZ#[1/eY####?n=X&1.,#3##1-C3e.SL:OG#3,#eZ$h%Z3,#DR)w9'0@+K,#g)ZV7)###($#6(Z@H'<o2yb#5u#TZ$0(ZK?%#[)F##;~+@]%N%Z'#####n%#Q%Z###?S,'Z$###&##$+=j%,3l$###hu$s^)$(;######_A#D%Z###"},
{443.24,303.02,21.46,4.03,"AR$%,Cgu&SG#+(&6B5%c#qb#:7(fG$H5#^6'Ac%7u#D#$4v$Z])H,94mMZZ&@_Qi[,]G#L6%6FD95####I##Uo4.,####6$#7A-yY#.X1DdPhZQ0,#-Q#E9JYM>'0/h,&sZ#BO:$x/95#G##z8H###1##kf2J6'(##M6(i`;%Z#Eq,CZQ;Z$|:&$>D;n.###"},
{480.42,372.32,25.01,1.18,"U;=#$%Cc#eP#(##Py/;O+_5%.##lZ'Z<)Uv*######-)#G[+ES//R$CU7%.'p6%hTJ0QH=l$1JJV^8=U&q/0uQ).,#C;#-{>###)##/q5kG$#?%}G$Qw=}D@lRU#H%2I$odBT%.Qx1Lg'eH'TG#BQ$D?&jG$Jc#3$'N$$(&1zG$.u#zu#,h7:5#xQ(G~&0%-"},
{480.42,372.32,25.01,3.45,"`/_:5#)##PQ#k1=`[(/w%3d&90/V@,nP#m5$Ru#Ll#We.ZP#d0_2H%-u#7##TM9K/?z]5hP#w'NVi<s>%zG#@B1aP#mQ'Vl%_/_?5#8Z$_d$DT4pQ$X$=Cr:wQRR,$Dm%y#5z.0.,####Uv#S/_###'##N-$AB68##o$$B18NQ%}n+Kc$=L1J%+RQ'###_H#"},
{285.70,396.21,26.44,0.43,"`5$QVR>u$S5#@H%^F3Fu$f>$>5#Lw&ap7;#$[l$C)+3TRtZ'H,$Ir1DS.>#$Hw'{TRVu%rY#DF=x2Aue,/Q$2H&5,#MJDnWAtR)+R%XSRjG$u?)Lp.zC7Xg+]SR7%,yb#+B%m~*d,%5M5oQ(Ao*C&/lT06#$}l$y6)qo-a'.%f'oS-@;;$[$~3.%o-ed)G>#"},
{285.70,396.21,26.44,5.01,"<q7bf0Cu$M~)0[%2nQ0c$Tn+Wh2x97e>$:Z%XqQ.,#-##_c$0v&6B0gy/P$&2y/<:6S)2<[)WoQQc&sY#1A(L.BTG#C#$jc%k&/|97fc$E2/v6*lP#=g&YnQ_NCon-j,$XX<gJ*P.-^>$oY#D@)gT3$##|_0Nu#x`8[R++w)x>#XqQ`'7+Q%w2+9mQ#l#F5#"}
};
siftPoint_save sps5[] = {
{469.18,68.25,3.71,3.95,"/%J######-?$|X[######3##`V=######D-#@U7######O##B(:.,####m@&-V[######SS#;#K.,####Dp,zm+F>#95#O$&w':pP$###|&$0Y[z,'###SJ#B<S0Z%F>#*6%`m)rG$c?()Q$26'w,$wl&m2+hr8J6'/u#}$$zU[ZP#$##v6#_&3$##{,%1A("},
{15.33,76.51,4.40,0.10,"############7m)######Y##kvVOG####t},}T6zR)_I'*7b############s`E######[##'2bWv*###}B#[+BkD[&m&a~)############@+I#########:2bW>$###'##Ix1?v('l#mY#############w98#########x1b#########'p4fG$PG#Q>#"},
{621.52,75.89,3.96,2.06,"|k#)##Vm$<cOT7-Y,$jl$w3@cdOxP$1,#|,$_x1g>$mA*fG$m8'(}GLS)`C;I@(T&+AwN$J.reO7I(fQ'eZ%Aq7}w-M>#;5#D;1}rD<w#e2>>dJ:R)IN3fwNh|?,02U+:M|AAn&bdO)c$0,#1H&%##vx#xjGxh5:h8]f)w/0Wu#Vs6j4EKH'=-'y$+%m'jR'"},
{621.52,75.89,3.96,4.00,"f6*O#$T6&1?$###{?%a5>0Z%.l#WZ%}<>bf0A7);Q&MQ#Ao.[P#l,:Fc%@u#B-'#_P'//y5%/~P-_5?B06477A/wR(#5=L26ev%Y?CPu$tb#ne.tS+,J)x}@m4C]b?uy8qC/&.M*+?67,h5#eB.AR)P#%G>#W~P+l#WQ&+.'IJ/K%&=~PN6&]ZPNl$_@-K?#"},
{606.15,79.44,3.75,0.59,"w>%AM*z,'&##)2;O~&`(7@b>zc'Y_0.wNs=Cs['%=@K$)3I(e?&22*me0###em@b%O5@+bZ$-p3pP<^#LCX?{Z',H%Ql$i&O&##3Q#Zr0EI,<S-6T1,o(6.+s%O$v'3,#wI&{p31,#/,#%:2######T`/s?*u>%4,#}|9iu&`%OVG##$'b5#{E?,Z$+u#O##"},
{606.15,79.44,3.75,5.27,"###cE22v(###0Q#y{0#.,%##RV*+x0###>$&<),~#&%##V#%6,#OO;3H&j#&`/1N/)tm*t%.j7O|?*VQ$YRA?/'<15XO<wq;K>#-?&m5#p7O.$(2,#t%'w7O16>k,LV|>*y2bz/S@N;5E*m'i>$6#$R##+|=1,#/,#@h5^y2},&5?%]8OuQ'D4B_6)n?'-n'"},
{417.31,90.35,4.15,1.15,"###7##uK,0H&1##0S*RL2OG#GQ$Rv)6@%pf5.,####D$#p_?###,##/$Crl'df/957Qh~bP#Ti~x^6Je,_y.u5&###c6%9i>.,#F5#SW2`l&-80qG#zk~wz:_f~gP#ED/}[Bcd+<5#(.(uR+###6,#*d%aZ'######`2)=o2@H'.,#Do%p`?WH(:5#V,#K%+"},
{417.31,90.35,4.15,3.96,"~,$6o/yY$/,#1e%Gy0^?)F>#gk6]J0D>#&##T.+=c%###(##C8*an.a5%0##P13Mv8LeZ3Z#@kZT3:%80FQ#?*:wG%###X>##w#x1<Zu%###an-pq+*iZ(a;XfZc5#tI+L,3=6KDu$###e,#-&&ie10,####&C2j>%>Q#hI+As@###+##^%(%(0ZH(###.##"},
{70.44,94.61,4.07,4.23,"*;9to.;-&RI*A=4hJMu[)]P#o~CL6(F,$8l#0Q%+l#??%A6'1(6Q>#Jw%hw-Ns@Gc$iI=2R)DUY###kc%dZ$[83####J*P$'d4;gl%`5$B?&}`@F:))TTic${RYM,#AR*`@#.U8lP#Gm)Zc#LH%Ep.;5#J[)wU;Dz)uG%w'%qRYD#$D>#P'#3K4Il$mP$#$#"},
{159.43,101.42,3.54,2.26,"3R=g$SPQ%y6(TJ&{%S_:57H&UH%{,&v>%6,#EQ%w5&###L>#<s/m$S`[#N'52)S(HN;5#OS&my0A14eY#-##7o)*m(###9##n>%R5$F1#b%Suz>=l#WR)f9Khe*x:1n29ZH&}i,]@-cG#qG$/,#0##*-$$x0###$##0%&WJ0~P#G,#TW/4d)?d$k$*7],66'"},
{665.83,111.22,3.45,3.35,"'/1######0##CLf######@%#WLfoH&###:-#(7+r[&Du$n['M:<######=##iLfQ5$###O%#,0~/W2###O##0ZDtC4pb#vl#^rC######2##GNfz6*###T##>*:^X>###%##dM2JA0###)##wL=######(##yMfZZ'###(##@D5Un.######8e)1<=######"},
{76.75,136.32,3.39,1.97,"rv&WQ&PR+Ew(V*WeY####BI%t3XD>####H#$pT1J~'YB4)m'a,%=-(E?%L9/*dGCS01,#5.&t3Xi,&###>##Ey.=r11)>$##/,#Vu$vB'vy:2bI%$')R%&1/Z.XPG####N0#VI,@c$c@)zY#FZ%w8+@K.Qv)Q.X`H$Cn))J'Q.X+c$###~%#d5%4.(rb#*##"},
{656.11,141.39,3.84,3.42,"Kx3######&##Q<q######_##U$U0,#/,#-H$<5#B5#TG#cl&p':######%##V<q######a##{x_F,$3,#Cc#WG#r,%G>#[5$}p;######(##H<q######3$#6&^F5#pY#Jd&###Qm##-&%[(9g8######'##O<q/,####h##|>P1w&ZP#c>#7S.Y@%}d+w]-"},
{81.69,158.02,3.66,5.52,"yp(H82:5#0##8~'/v(g,%l>$qb####=H=;H&###B,#<:-^H(Fq+UA2kI$S6(]9QVZ'/##N7'>e.UG#dV(y%-###+##f)+$m(vB7###}-#3fPf6QdG$<R#}9QUe,[93Kx)tQ(###x~,ff)Fl%UH@Nc&2-#Of.~N.+7QTK)obI&##tP?eN8{,'###H;41H&###"},
{174.67,174.32,3.99,1.49,"GT5%##%##:u#p0c<##hd,*##*W@-$#q5P[5#i4K)##:Z%mG#B^9G>#;5#K>#Y0c5##wm*wP#x;A]H$v+B7e'2lN)l#y$+cZ$]'9###E>#-c#O0c2##0w,>##F4D*~&O~X$c#q$RJl$ec'5,#_A1###%##O>#82cL#$]-*+##?i?WH$ZdUG##.YJ@5#h5%D##"},
{500.35,178.04,3.65,2.68,"%##[-&cG$###|l%$@*[P#&##i39iZ(###'##CU.lz7######@5#n9*Il%PG#j~,I9Opb#4##cZC}098%(8D/<-$.S@cXBlH(A##j7Az#'###J+BQD,xY$@##}{U:4Ck(7GV-TH$w22,yUZH'`,%g/.9S/4##p6(Vu$ub#Ol#F{'uq=U#%;d'>d$1U1H)=Gl$"},
{470.31,207.67,3.95,1.45,"2h'McG2?&yb#1~:ec'###^,#z$'*A,95#.##$H#OH'j,&(##0(3M7,[P#'##/Y`e>$}>&A##d{9th3GM?7,#WW-G7,zY$5,#V^::5#RG#[G#AT`'##uQ)U##I;>GQ$|T`gu$n>FK,$[Q'2w$1B5.,#@##iv)7T`%##NH'JR#aB6wP#.TZ&Q#TYEUG#}Q)/Q#"},
{265.09,216.29,3.84,4.47,"[?)$##/,#+c#1%Y######_.#7GMD>#:Q#E9.Sc&=5#Ic%&J)bJ2qb#.,#7##h%YMc%###m$#@5GdRP###xG#we#OiAdu&H>#pf6.,####o>#O&YSe*AH'+6#Oo+/`QhB8%##*{OG|BFl%0##(U8/,####Q##NGN8H#p84TH#1~.j[#7%YD,#R%YY>#bc';$#"},
{470.52,219.21,3.88,0.81,"R^9######^1#BIUJ$(.,#S(#$9/|w0###O5#Ku#7-&{k#0,#U1=X##W>$^1#d81)=.Kw.ZU#A^Tt93N5$X>#nR,oP#8#$ju$p?*PK#l84M_$5T30q$]T5j],2IUxY#K,$@<-ln0###$##wl#X>$_B#Up7U$&1sERQ#lu&Nz&:/2###<5#RK&Ec%###f>$m5$"},
{470.52,219.21,3.88,4.63,"###Y,#;V;###5S.j#$R/RmY#T']+Z#C::qP#Gc%bQ#%&]$##*##Ec#vC:.,#f$*US%aBUW5$1']y,%_;?n5$M?(-##I&]$##$##E5#X/.pb#Kl%cG#]c4hw.paIL>#<`8g25du&,##q)]8u#^P#xG#]v*###>Q&7,#=$$Rv)ZP#'##dd&ao/###<##;d706'"},
{625.40,243.17,3.68,3.45,"z/5#########U<n######2##HeV{k####>c#AZ$l>%@5#oG$`g:######)##|;n######C$#]1gpb#$##@m#M,$Z,%<Z$c,%w09#########P=n######1##pzbeY####+##bZ&oP$1,#@5#GA1#########L<n######9##k8~######O,#7d(######D5#"},
{618.10,266.95,3.51,3.45,"y]7#########G`n######B##fdR95####bu#/d%Hu$0,#;Z$#M=######%##K`n######3$#E&^.,#$##r~(M?(.,#%##HR'dC<#########Gan######4##'`emP$###Nu#^H&k>%%##C#$of6######'##E`n######u##re~.,####:v#Y5$D,$Au#Qc%"},
{249.23,276.22,3.67,1.31,"gY#0##>#$7m&fG#):/N?(|P$oI'?&07l$-Z$J,$O>#g#%^,%^P#zc#]&4u5$G9.&:FM9_I5#:?_f+Dz$*+@%+w,O>#Vv&H6&/l#Ml$1C7VG#>f2%R$:?_DE;d9_xc%Q0.'n9Pw-LR'e[+Au#Dc$>?'[u$nY#tb#@5#oy)Ly5od,)J*5$&Ch23?$v~.@u$<5#"},
{84.36,295.62,3.57,1.22,";Q%IZ%~>${b#BQ$1Y:^#%F>#I&)hWAN>#j>$RQ$@[*YG#rG$(##t6&x^:95#fC6]8?1p_h>$Kt_s|>.~+WK*rd+sl%[~)1A,$##&n*Uy5j>%^&4rl#'v_4O@^q_^I*q~*r5;Cn*{o/R6(V,$###7,#Lm%ac'95#%##_&%qf5L7-PG#H6#1L2>-(<#$?5#K$%"},
{501.66,348.66,4.12,6.02,"U##CUX<%&)]3Z.&YLR*n(}K6)N/4o2###1v$Km$4~.VG#hP#v##mSX$##:5#{49?SU###*##.WXHl%###=$#8@+~>$66$XJ(=5#.Q$F>#SG#Z-SgG$###V$#*SX######c&#D&3[>#&?%8R$yY$K>####.6#xz?######r(#URX######O1#jv+I5#1u#j5$"},
{256.49,349.60,3.80,2.76,":v%R#%ln@#w,'.&dc&;J-7#$tG$>5#]>#[Z'2,#F,$P>#Gl$@*Amo#VqX^Z#Br;Sm$/mR(##D.,o7)/,#?##yl'f$)###R;+bQ(As*i?U7m$px46I#'^7$%%~aHfw's?*O.#,4IeJ(Jg5|;'###Z~%EQ&a7+D>#(##Dw+J7(TSX)##ze.H?$'qXJZ%Ra4Nd%"},
{180.19,355.61,3.55,3.05,"@-#]v*1^,W>$)n#_?)fY#PG#uG$07)###9,#T6)ru%###$;&1n*Sl$2XYOl#i3:fu%180+##zWDa.(O$*&$#xRYM-$U82G)%3-(T(#]SY+v#GV><,#{~1k$$T[S?,#SJ-EZ#fXYk>$3_5|c#95#[&#C06Nd(5?'%##Q,$3&'6['###2w%)?%O+/k-+>Z#|Y$"},
{191.94,356.93,3.79,6.09,"~v&N#$b93###'?%*/'du&)##<95|v)95#=M#?o2-K15##ui)kv*N-#;YKHc#>L9{Q#Fh<g5#pwWGe$Xd+';#UR+I|V305'9$A,$X##q[+Z7*io4&##^D2fc&~|W(d$701]c%;n%Px(;RAU3F###rP$.##:l$%l#_#&m6#Pc&)1+Yl&eQ#A$)U]#5<DSl#D%."},
{309.38,389.54,4.62,0.45,"Iu$;UIE82wQ#x,#abB}Z)$##^m#FR+jY#G>#5,#1Q$>#$Y>$%$${1RrL:D>#Hp*>EXac'.##iWU781.,#7##[Q'vY#;#$H>#7,#A48rbL+H%/p5$9/-?&uh+:AXO5$I>#gU#dQ(P>#S,$c5$2f1@[)Gu$Iw&c/4OG#=5#>Q6&L8xG%;5#0^$hG$L,$3u#ZG#"},
{20.13,409.82,4.02,0.02,"^B7#########$4o######*##0JU&##D>#*##?#$~G#`>$|k#9(<######$##83o######o##F^cH>####;$#r>%au#7#$A5#J^9#########.4o######-##@qetY#.u#4##>u#dl$6l$$##}S3######$##A3o######.##>&ZSG#P5$%##pY#6u#Yc&.,#"},
{147.12,410.65,4.21,0.82,"OG#hR#iv+Um%.,#IZ%c>#S6KfP#'-']/$iV@[Z&ZP#kZ#YZ')$(*-#'d(D}4:$%?eDR@,HZ?rh1p1<rt7-'6y#'$##@j:L,$,ZKQ/$1|Ds%'q~0){)0,N/@$M-N?I'<WB)&#5l$;?#C,N+##y+@UH%@):_5%Q-)'w'&D5R=9ke1%z$#J/tS$###YI#U7._?&"},
{147.12,410.65,4.21,4.79,"oMAf6%b/3[[#fD.(E:zT5^Z&rQ>CZ%$W.W~0+?&1H#+k?tl'M@N4,#6;7?5#MU2P>#dAHc$)gBN*$%-2>Sw&Um)Lv#C>N(##c<0T,%t#%%##qL&Ag8sU(k;=Lf31Z:2f1g|4IZ&OT$mC=5##Jc#yY$pY#&##-.#tDARu#oP$zY#[lF^G#y5%Nl%{l$OG#U##"},
{270.89,447.02,4.43,4.16,"+u####O>#16$:#$+Z$al#ltAN,#|H)jK$tMALv&A,$8?#C?'LK6;,#3l$6T#3g2+l=$8/cO4Ib;h96HG6Kp7W-)####c;X#%B?QK,#V~0`-#S)@1T&d?Qr$%$@Qvv%e<B3/$.H&;R#Y@Q_G#|j=P5$U?'F>#RU0.S-hD4:j=?g8uD*b83T^,D>#,y$E&3###"},
{117.71,449.36,4.01,0.71,"|A$PuO/Z#j6)pG$;~K5,#(I'uG%Yw%tb#Vs10,#Q.)>#$K9HV>#j5%>5#.n+i'8{c&OG#lZ$w$V`[$$[)LW)*v'2O-'i@%'*:5#L,#CZ&0,#k*C>6$)R*%##z&ViR$r;>nP#lc'B^#5SQ_:3###9,#fu%.,#An)gP#*o..,#xAE%.,ym'5l$^I$?aBDA'K;<"},
{170.77,459.23,3.66,3.04,"U7#EC;###[P#en$7m)###F5#A6&1v(;5#pG#&I*oY#1u#&#0jS(#B4B5#3c$E~G#-'OG#F,#v%Xf6(Il%C$#Wf]rQ#e1:+X(ql'.,#`###E=,+I$##$c#p'''KWYG#(S(%v$<l]k$+6h286$95#$##Tm%nU63-(###6,#~9(IH%OG#_%(*?%sW'yS3TQ${,'"},
{123.58,533.97,3.72,3.81,"n##|ZG.$(.,#A])|d-kP#-c#xP$###bu#Xl%5l$$l#1,#u>$kh,Qb[r5&=##Qb[Lg9%##o,#b'6:5#7?$IZ%95#=5#GQ%ZZ&/)>[e'K//-Y.d][6#$###bM#hU9IH'+l#KZ#eP##6%QZ&f>$8Q&###GI$Qb[J..###$##c<)_,%pb#4,#;H$WG#<#$iP#W5$"},
{430.04,539.84,4.05,2.56,"4VQ53?T,%,-#|?>4VQBQQ0$$,]/_S&rq>;8,:H%^P#$d#uS1}8'*p3D6%FnL0h,B-M+o*g:4*LQ>@*-Z$5v$m$)fH'/Q$p5%~u$eY#3-#dTQR#%tb#gV(ARQ/^7PG#|u#+JHTc&fY#9##AA-HR&c5%dP#zG$d#$ql''?#e$+V>#%-'/c#]A1OG#<5#x>$eQ'"},
{430.04,539.84,4.05,4.54,"=v'sG$J#$]#%D-&Iv(0Z$[P#d8-^#&;5#lP#+Z$R#%Ql%{b#iS.SJ*Hc%5l#wl$p$>gd+}b#wnFn97~P#u>#Bn,>5#P$)K>#VZ$;%)6V-@]SRI(qf*=yMIf0k^SfH(:#$%f$hS/#$%AR+kG#8R%O_SK4?GAR7m$Ca.|]SFv)D[=k6)}u&_u$'D0c?(<c%Z>#"},
{448.59,552.92,3.67,2.72,"($(C,#Po3-##8eZ)##'d(K%$4'7S-#heZLH$-[&#-${o6=5#N/2zZ#=A1J##*gZ+$&#/.gV%on/q%&GgZYI&h,%zP#HB5q7+F%*a&)lw0%###;MmfZ:Q&G,#X8(p*WIFIC5#|&.JH%QQ')6%'##?7%jo1]6)oR&e&30l#w(3<[&uJ1u6&]e+yg3|G%}b#Uc$"},
{448.59,552.92,3.67,6.02,"|5$AR(|//Pl%jZ$P0.Ed'If2Xl%.9,v-&|p7vL:@u$%##A@$9?']G#iT*3I(1bJ@5#)/'EY[pl'/##GqKrT[]T6$##kI*QI&{7.k.*^Q'S>#MW[Kv%s@-<%&991so%EU[--&bf5>##^'8{>#Xd+~5$&d$Lc%JT[Du#@[*e>#We.?-$/T[$##'f0'##8x1%##"},
{460.24,558.86,3.72,2.94,"hY#P>#[S/(##iFI$##ln-_5#(&Y###Zg1$x*#-&=5#=%MCu$.,#B5#6e+PG#wjF&##1179,#f&Y###o_:E@%g#&mP#i'Yw>$&##dZ%.[).,#9FEDZ#.q;9##y)Y_A/FC5t@%h-(hK+W(YFZ$$##S>#*7*eY#h-'F%'kg:###6G1k&Y.H&7,#h6%-8J)J/mP#"},
{442.59,562.36,3.51,5.89,"pY#G,${?'Z#%&Q%5,#-(.mm(Bd*,c#5S%<DVfY#aG#i{*j4Ij5#=]2W#%:5#y;AAe)_v)IQ#0]W8v$5L3uF4uu&Ow#g_W<J.'##)-'J.%0Z%5YK*$&QQ$S$%B]WAn&uy:.$#@6'rf$P[W)#####(##K.*OG#:FF&##v[+xG#w~Wg?&YC<X##Ml%`Z#O[W###"},
{280.03,564.46,3.91,3.49,"zZ&j'.K*Wj#&pn,R@+=@*0?#fZ$r#'tb#B5#}k#$##kY##?$`7$K*We'W_P#Je&q+E992Kc%(?B#80O>#uG#xQ).,####DH$tg5H?JzQ)'7$y088H%}n$C0QPYK|k#_m#(YAn-*_5%'##L#$[y8XG#~H'`|,so5)n&A.(TJB|k#H[#O39@C9(l#^P#/Z$uG$"},
{280.03,564.46,3.91,5.39,"j>#Zl&###$##'Z#_v)eG$###'d'-I)W,%>c#c#&/,#{G#f]2H5#rl%},'kY#N#$R?@ue15,#MCWisBN5$8v#6[G.n-[[#)s;K5#2w*5%-0,#8Z%dd%R7QQ.,@#E#<>Y@)@n)JX-p@W.Z$jY#~v(y7.&/&8N5o-+al$R`2oEWL>#i$'oEWd?Q9d%nw.zR+U02"},
{485.54,583.29,4.24,4.76,"?g)>{3PU8fG$[nA+M=###)$%Sm&V9U###=&+D>#827.,#Cp.+:3:-%.P2CGKT8UeG$/~'e0/{.0UZ&}5%7P>)Q$w#'K>#PJ+(W3L1/eX>jn-m7Uk-$#M8em$EW?;x*A3CSl$up)IC6g@-TG#+A-_m%PS,#f*(bJ|G#Ic%,r&*S.*Z#F*>Rl#Q92Z5$j'0p.+"},
{451.93,587.99,3.54,5.20,">5#lQ(G>#8u#.,#8I$R@-;5#3l$QH#G&2qT&rS2tu%ub#ji)o?%lv*dl&Xu$c[,H]%^1=]G#AnVTR'Fl%2M#G&1'O,tc(EM$0o%JD;yY$R>#m'4Ts4o.0*##6rVl4F+u#I##LR(TxD%Q%A##O~.}?&v:7,]+.%$bGBZ::95#m{(a}J###(##e.,,]-###_@%"},
{133.37,614.93,3.68,6.01,"F>#D,$D6%$l#%Q%;5#Bn&Hn+[?)+c#j.&E(V0u#U5$r()r}IF,#uw.^P#D>#4jDAI*7Q%xu#LoZCc#C^2$c6OZ&r?#jrZle/+##|>%4c#1I+4QPEc%$c#)m%_pZ1%&J97I##lc&9T$InZ$#####%##Y-(j>%G}I:,#-m(J>#HoZ+?$:K5(##;l$PQ#WnZ###"},
{299.15,616.76,3.66,5.42,"xZ%8l$=Q%%?%}7,T>#mx.wP#<1MDw-Ev%Jp&$&X,C-yZ&2`&eY#7,#jZ&h$*ae0`Z#uS1Bc#J&XlH'b>$>M#J|CIa,w0:u&####A###I)95#G943-$ZU:+##Y*X6YE7#$A##L@(YpI)R*(##'##a,#YI-###Hl#/%&pf6###@.%F|AmP$'##Jm)}y2###[##"},
{18.86,618.25,4.32,6.25,"%~.######)##*|p######f##I'a95####G,#KH%`5%%##&l#Q]5######+##/|p######=$#VUdZP#####?#}l&Mu$H>#ZG#$U8######)##V|p######E##yB_*c$###B,#TH&0Q%gP#cP#eo5######5##6|p######E$#ltJsP$PG#u>#cc%Z&.Xu%:5#"},
{135.32,630.69,3.95,5.97,"0##X7*#c#0Z%^J0Nw-(l#y5%]La'$%Id)$K'oo3W/$EMa`m(###5,#?-%0%-Gh=$c#Pm(P$(DLa9$%K~03-#4&0_x$FKa$##/,#RG#>R)sb#[2?C,#6^7&##POa}I)%~.3##X$(`U'EKa$##3##bH'Uu%###'@(h,$(y5###*t4^q6Tu%###Jl#0t-T7.###"},
{405.84,634.30,3.85,4.55,"-L8######8,#nXp######&##(5Ju#%fY#%##bZ%Yu$D>####r:=######*##2Xp######0##umVf>$iY#/##nH(S5$0,#9u#Kh=######-##@Xp######/##P.XBI'-~.;##k/),y2WZ'?5#}/5######3##VXp######)##.[S|5$$bGT>#KW<0m&0:7Zl$"},
{55.26,643.79,3.75,1.28,"RH&2u#Au#Q5$[,$(o.g,%b5%ku$a$*Km&m#&fY#6?#`:3Hm)WG#cZ$h7-vP$4e(gfF'kH$##@j]sh<#v&vd%>n.+$#Xg]>e(Q>#2$'}o34Z$=f2Q$%Ri]O['*g]cu$M&0~B'R)?2R$Bf]}u#(##AI)L-)###<c%N#$F`6H]0@vOmQ';?%yB,)%)7A(kf]$##"},
{139.64,643.68,3.71,5.55,"2,#>5#h#&P5${k#.-#=A/VG#=aG^,#,H%Zq#=7Wbx'H?(p1#vG#N$'.H&###a~13n#R08B##*8WtH'(c$e(#IEB]<+H`Ca%#z-%O19qP$>5#z96}@''/1(##m;W<?KeY#6##fm'A2Q?d*$##@d'ju&YQ%Qn*JZ$G%)OU9gP#gm$0RO1H&/,#Wv)>i7###C##"},
{19.01,700.49,3.80,5.65,"95#######,##>%.######,0#wFL######>(#uND}k#2,#V8$tQ)######B##k?S######>(#h@SbG$###0(#pINMm)tI'W=7B[).,####&##2DS(c$###?##{?=t/5###1##]6<%tGEQ%Pc$Bc%###%##?5#rkG8#$:5#7,#UM/fK6|Q&{G%aq%:rAqG#6#$"},
{623.91,58.01,4.58,4.00,"7c$NQ:JR*L,$e+FF'1w~-v(5be/%7(UU.W`7oJ3###%##~8$QT1x/-:v&~g29a>w3>fM=k'0n]U4*>-7+oZ#7'7$##yb#dR%@^UEH&Nu$N?%&i<0C)^~UK?%][U3H$FJ0_7#9C:c#$Gu$al#sx0N8-]Z%$.*<p5FO0Xv*50'+[U`>#bG$$M#(@+V6#Uu%0v%"},
{181.48,88.72,4.18,1.30,"'##$v$-S/###S>#$y295####z5#uA2pb####7,#|5%SG#E>####uH#@z<###+[&7_M[2B)##6<WD$P<#$/e%)S-^P#u$&yJ.###]##{=ITG#Oo1C&'6CXHu#1FXKm(w$*(B&/39Te)$15(n(.,#C5#4K0kY#zP$&##,VMfT5&vG-c$?I%.m=ay24v&{k#}5#"},
{245.79,104.85,4.29,1.26,"P>#Kl#5?&ZP#F,#Dq2ru&pb#bd%;K3-c$L>#@#$Pl$-Q%F>#0,#x,#]U9.,#U'0xpH%']&c#C,]65Fu[*.w&8R*Q,$0S*%$&$##4,#._6sb#wn0R6$C,]%[L;']ju$vo,xpH1%,6c#{~/eQ'$##M,$~#$Q5$5u#iY#)n#x1=5H&.,#cu#D<<OG#$##ol$Ql%"},
{60.79,105.75,4.50,3.71,"<6$i?SRG#6u#s^8saCbI,Vo(p..QR*T7)}&KU,%>5#tI&co.W&3|Q)AZ$sA+]k?;FE$'2RS*VDSO-RVZ&-[%Pe..Q$zu%Lm'xtM;u#P-'.%&J95zm'yASP%&~@S[u%?[)Nq$WS0;l#8$&Bm%$'282<iG$R6&U&3cg/KK5g:$vWF,###6&-V$xY$.##k[*$e'"},
{601.72,119.35,4.36,2.12,"3$%:aClz(l{A*j2?I+f#$jK3o4>sb#:5#L#$a#&2,#tb#S$(3,#7o.0z)iV@+x-b'(z#J-94R'W|l%'Q%0Z#k@.T5#~$*+-&9##|(W+Q%PG#6WA|m=;x2id#=%W|>$+u#P&#B~/r5#0A03##^>$U^6A5#@i8#f2u-(4,#uZ9Z2B'###l#]_$%Q%v%#9p7#?$"},
{126.49,123.79,4.65,6.19,"7W8sP$I8*|c%#8ES]3b,%3H#17'KL5u#&@%*v?&of5<u#H/0(U,tS21Z$ud%Zs5tOJo6(]m&hpTqv+hY#-^.-R*.,#+H#W=EmQ'|P$-n+yxG^V>]6%PnTRp,pmTjY#&?&39#8J0/,#,Q#.6%###s%#g{@aM;$(;(%#PbLqp$*mT$##[P#X1#hu&YG#F>#=,#"},
{125.26,137.11,5.01,0.23,"kh4S~/fu#DH%lh0=OCQ,$9c#JmEk,&2,#(]*yG%###A##Je,h'52]2o>$S''l*AoR,>k>2T*(;Y7#$4-'<,#c..I>#[G#pP$]-*E5#b.,NzL^98S##L9YFf(~7Y$##am*$A#S~0L>#|k#|Y#%##}5#9#J;_9tn1I$#m7Yfw$M7Y0,#W,%O9#a83QG#.,#.1)"},
{280.16,153.88,4.97,1.49,"95####qY##l#tkN###$##R5#bDl###L#%:##Hw.$##FbK5,#OG#$##4u#jY#`@Y###$##U5#aDl###n>%M,#Y83###@IUF5#:5#%##A5#qb#$.W###%##lP#jDl###Y,%?##Hf4###z$U5##qb####I>#]P#ojI######f>#]Dl###{Y$?##h&5###sFJ5##"},
{264.24,155.63,4.84,1.49,"OG####lY#iY#:QQ######>,#aDl###,c$C,#0e.###T4I9,#F>#$##?5#gY#&JZ###$##mP#vDl###Bu$6##a/4###XvS2##B,$###I>#%l#)7W######m>#cDl###,c$A##*L9###@vR8##.,####6,#/c$naI######eP#wDl###C,$1##O975,#LFH*##"},
{116.38,159.20,4.88,0.69,"yq>~l%3n%#~(yZB.,#Q^*@Z%z/N/,#1,#/,#n6*1,#gY#>-#KNDD>#AJ*(`)|U:###mpNhc%MAQ95#LZ%R5#ecK~P####+$#*c$TG##q1j+5-%-###wCQZR)-?Q###fm(Uw#6L1.,####0##^*>0q9(S*mB.a7.[d$pCQEq3`(>'##q[+6~$'Z$K>#:5#YG#"},
{198.94,162.37,4.01,1.48,"95#3,#H>#XG#o7[###$##q5#-Vhk5$&d(;##s$*LJ,<S.?,#]P#&l####:,#VVh]P####2##LVhA##XA2(##6K1Zm&FiA3##:5#vb#%##A#$'Vh0,####9c#IVhK##Fo3,##Bg7i##7jEE##95#.,#5,#5u#g@Z.,####X>#,Vh;##ym*8Z#|[-vG#ao0B[%"},
{73.71,165.88,4.84,0.50,"=#$<5#-6#aGN@6&Il%5##8KNe5$Zu%~I*p,@.,#$##-PCyQ()e%vo,1,#Uu%no)nKNOG#'m%nKNtaGW-(Gy+NR+JS&-JNC138w#Y{5uu%3%-%A-y(.76'`q1|KNa|B8#$wA'}[*[,<X/3W@%5,#hY#S?#~tL}k#$##?5#`HN2?%U,%###,*7(l#qP$$##gm&"},
{73.71,165.88,4.84,5.51,"qd+.?$y3Fn,##k;am+rb#/%#^d&;c%V93[#$=5#I,$^$GzY$4l#f5%gHF-u#+|-S'80/'Cn,,pOTH(>m$?x*]l&L>#T3/7d(D>####Hs2s,&s{>+u#8A$2:OTmOa.,:o)ypO:H&j9-n^09Z%#l#/,#[v#-c$>24=6(kR#h.,^g*NnOa~*+f0###`48*~+4[*"},
{356.09,167.35,4.16,0.10,"######fm$,%-fH)###rY#H$%|e2###<5#m-<eY#U5#Xv({&V###.##?)VI.._aI###k..g##X$V*##t#'JU&I[+J6#d$SpR*###,8#'%V.,#FbIG5#PJ18##_(VsB7uG%@##5d'x#G;h9W>$####-#N%V95#_v&%?%4I+###wB([jG######8##D~R95####"},
{426.55,170.04,4.26,3.54,"],D-6'###$##4B3######}-&/,#1,####Y<495#.v&###qi2r/W######,##~.W###@5#;|*WZ'xc#DK0t#<F>#~0-mZ&GI(O/W###%##Al$v1W###N>#ec#1<@i,$U{-=A.95#~m%to'+h96r>###<,#%Z$Q2W0Z%######)j8O?Gx>$P5$###F?;CH%xY$"},
{539.69,171.00,4.51,0.64,"`u%0l#s/&@_<C^S:5#X5#Yc%VJL95####S>#/[(ZP####~l$.S)a`6z_S2HL=X?&##OL)nM>5^S###-##_~'+T2###&##:e*=FA$27y@,bh2A$)%##=9$*]Sqz>###:##N%AaQ(###E##h4@?~S###1##'[%gK4###-##v?'nh:######/w%R07###&##^['"},
{539.69,171.00,4.51,1.36,"J<-F82######pK$#U8###)##q%(du&###O,#1~(;c%###Ov$$g%dn0Q>#;f0oM)^rC%##2n*4VQfH)###w-$]U995####o/'X.,f%)]G7rTQvI,<u#qd'4VQFRQeY####IX0{q5=6(###B##V)28UQ+C4'7(rj8c[*($'F%(bn'VZ'###9?%aJ*&@+###9##"},
{539.69,171.00,4.51,2.39,"Ll#^7-.,####xo/xn/###$##HW5Be.###'##i>#GV:######sY#%F6OG####48'h_S###&##PH:24H###0##Am$QM6######6u#T_Q.,#%##X2?c,>###P.%>]S6R*%##Fe$h-(D^3.,#-##Y@-MA-wG$Pl#o[Swc&yI(^nBRh5(3;S?DC*=+$%hoQZP#$##"},
{326.15,188.84,4.74,3.23,"xc(######Sl#dmWH,#V,%)2)iv+}h(H'5:s<###*QD(Z$*&0%B3D>#%##R>#(nWX>#dd)h$$R82M&&48FdlO###a.+s5$KpW9f3$##EI(i#%aoW)c$=,#Ue$vA/97+uZ#rPF$##<H&(##sIR%055##fy86##5oWz,'###5T$i.-5{>###0M30,#x+I###}S0"},
{489.35,197.61,3.94,0.77,"K(3@S/V6%J$><d(8n+h(-5uFO$*###;,#}^20u####$##;K,CX-V[V$l#O5#P`Vuh=sZ'{R'[S1###+##8]R######$##?vGd^5??'.,#'.#D]VeY####v0#slP%Q%###b~*Z>#v,&###tm*g7/$#####j1#_;A#e(OG#_(#.M2MT3###K##Mv&<Q%.,#=5#"},
{321.81,209.25,4.55,0.70,"D#$$U,%.)Uv';Q#/~-wH<&d(######~k5+u####+##Mc%,u#C%*Pd=I#%O##,hSN2;VM4f##4$(~Q#+iS'#####E##KA0###:w-Rq$%Q%}^#_dS~r.GH'S&#u,&Y))_tK$#####h##y$,.,#3l$).#>u$'_%O5$DV%:Q&a#####Q8$X#%95####T>#$l#.,#"},
{321.81,209.25,4.55,3.75,"gY#.,#$##V>#Y5$.,####^J&E,$N,#QG#FE+{k#,^(.,#Nx$37,.,####o##vaI&##~u%*j,iu&c]#DlMpF595#R2(f6*O''&o0######M##1rS,##_6)wZ#4q4P$#(qSBr;.,#%6#A.*<AC4H%M5$###)##p%CZP#0,#$##,SAfH(K-$g/2^#%2v%'c#6D0"},
{451.10,210.41,4.65,1.57,"](7-[(M5$$##d:YP#$}6+###by0y.,sg9jP#qr:L5$###@,#6z:3,#1u#Fu$x8Y###T[*4##{D>cP#S$NuP#R@O;5#gG$T5#<@,###7%#L8YV7Y###,m&/d:g98###`mEBI'v6U###Xc$^v%######aK(i7Yz,'###fl#d&ML5$###z9(Rx0bR-###n5#op2"},
{283.61,214.91,4.45,1.53,"`[+;,#OwV$##>8V{b#xw0$##[c&N>#cIV95#$##0,#gZ'X,%e[,###TwV###;wV$##yx4A##Om*%##nwV3##%##H>#`R,:5#kH)###oyVGH&.wV&##iB3x]*9I+)##FxV]Z%###,##5J/.,#[?)###li/cEC|Z)###xx+R{V{k####tyVC@)###$##J/1ZP#"},
{283.61,214.91,4.45,4.58,"*7+ZP####&##W&Xic&95#&##%r6R(V}>&###>),G}Gll'###^I-.,####=,#g%Xk#%sQ)H##qT3l1-#%X,##e'Xj$)qQ)+##7e.E>#&##M>#;%X@##PR,1##l]6u##$%X0##{$X'##x6+Y##]R,a5%###&##V[W%##yv+wG#qJ3;##^%X`G#'%X(##ql'G##"},
{301.12,213.15,5.10,1.43,"~Q&/1.PL8g..[_*Iz6`U.I97Yu#(6$z^V=6(###&##Hm)95#=$(ev$_dW(##:hWdR(7|E5##An+Om$hdWE>####(##r$+$?&zH*_##~dW*##]dW9##XsGn##8.-6##ndWB,#'##N>#P-)Y>$SH(2##;fWWQ&R[W0##^V<%0*TH(D##)eWY,$0,#7,#;?':5#"},
{301.12,213.15,5.10,4.63,"*[):5#$##6,#,JW6Q$mv+.##r)?@/)LIW'##}JWY#%{,'$##Z$*8#$'##WG#xIW>,#W%/'##@tJB##UIW&##]IW'###[)%##zH)GZ&###$##ZIWOG#5.*jc$9*E&##~MW}d)^IW%##W-(gv%6$(:5####%##m'UKQ'a>#f,$Ch/[U;6((z05g95`o3}u%81/"},
{451.29,223.02,4.94,1.46,"*TW<#$(c$*##nYGxS,0{?jP#ycCOH'Au$s>#(@)-u#;5#oY#DRW###:Z%&I$tmUjP#HSW]c#PTWJ>#f?)VQ#le0###%##GZ$BRW###k?%U_O?(<%##OyR`A)gRW###Vm'%K'r.0######_$%OG#$##=v$n=Cpb####Pq+XB4SS1###BQ#fr7B6(###$##9Z#"},
{221.98,328.09,4.18,1.83,"###'##'-'###T@-<##r-*t.&JT5qZ#/x1hc7'I#8W7=/0}q<.##O#%{#'###;$NgP#V?(RZ#zJX'6#FAOf.%[%&3T,a|WmG$[#%:#$+H%|b#7kHE>#{5%|6%cMXx99PT(:e%l$&31:4H60o1eY####fP#Xl%id)V,%H>#rP#sg*~6T7c$%##V%#rIX:m%~#&"},
{483.10,333.02,4.65,1.66,"S##5j>xY$###i,#40-WT5###&##s58M/1###+##OdLzY$###Ko(bYCcT4$##[BP%J-`T-Sm(-H&A?#O)T;$(###Bc$T&0D>#`I,L>#PG53x.z$T95#f6%+:,lx4}b#O)T/-&###1,#<c7.,#OG####ax#$y4?D<.c$S##7^0=e+K9.Zq2TG####GH%'@=T,%"},
{483.10,333.02,4.65,4.28,",T2E>#$##3u#QfX%##9,#0.(V&3E6$zA0uB4n$'YH&KH'hY#v84/,#&##Oc#(fX9H$X>$..#aV<4+.d2BmP#dE2739z,'###YI-###x##YcIffX0m%cH(UC,rI+Hr'RfXd#%EfX;[&Um)|?%.,#$##5.$PGLBH'$##bc$T`4,u#'##+t0dL9-S/###8$$7Z>"},
{146.48,335.65,4.31,1.74,".,#$##Y>$aP#E-)F##TR,06##D>kc#$?&g+1LH$*39Lm%m[NQ>##l#n>%:5#JjG-##/%,k,#1.XkH#`sAa~#-f(*=7<1XbG#^u%[P#M>#U#$f2A.,#JQ$RQ%(3X5-(JU/~c$c%)M?(ShNNw-7#$###(##,Z$fZ&$Q$Q5$vb#X|+X1<v,&$##g:#=vT0S&pv+"},
{478.76,345.93,4.66,4.32,"vG#2J/###lY#v.Q~-*###5Z#8TY95#Q,#D:1{?)&##>H%Zn,V5$.,#5,#}b##TY/,#$##(c#7UY8Z#|c'v%+M-'@U*h1;1@*nY####e-#R$*)eX###(##u@&2UYse&;n.M$#|m(W}/o|H$#########'#~RXAJ0###i##LLX$X@:-%6f/'L+*-'Iv#L6F-%+"},
{226.30,347.92,5.26,4.50,"e&3>Q%&?&r[$A'-~L3ZA1J<8=k;:v(5h$yM@(-'*##4q(ad+OvQ?5#-~-S-#[4Cs$'GxQB[&8zQAA&X(548(UQ&X:%-QF{G%zC1Cu$2Q%t>#%i+$z:Mr3@(4%{<.fF;y5Tr-{5&$9&=U9^>#y5%%Q$G#$>5#F-#1M<M?%.c$ju$&+CK>#<d&Wu%e#%OG#}G#"},
{590.92,351.78,3.95,3.40,"=y7######'##U)n######~##5~Z###Pc#Fu#I$(.,#,I#u6*X_?######+##U)n######/$##1e###$##Y,#;v(###F##J6'V{B######4##Y)n######=$#Uf_######k##E$)=5#E>#)v#32@######1##T)n######8$#)ZP######-$#_u%`5$|k#eG#"},
{151.25,356.91,4.74,0.92,".,#Bc#|k#8v$###RZ%R5#M#F@5#%H%Z/$wz>d,%pb#FH#T$)$m(C##9-(XM)wl$ju=_~0kc<621|C;.b4@B5Vu%###B38.Z$s6O-[#I(<y.$SA0qq(t5O(w%q6O~d'[h=kS#xG%86#s5ON>#ovAs5%H8.pb#VI(@J-{N8OkA}I/UL&2&21z*OG#+%#&x1hn&"},
{151.25,356.91,4.74,4.65,"crC[Z#L7-DS$q],G,?W'7qn+R#><I+${(=.-#Z$J5#Es4CQ&~ROG,#Gg6.##$V6_Z$jTOr?(XTO}l#cg7`.'M6')S#7QO&##Fj2cG$%?%)##^_(:p73M,oD=^o38u8SA0)47GZ&ao#,z;8##AZ#qP$'Z$%##4.#K<D:Q%7#$8c#JcFF>#R#$8Z%96$D>#^,#"},
{175.80,412.64,4.05,2.78,"dK$Iy795####t/'@R**H%,Q%3c$gY#4,#J7+/,#[P#x>$ku%(o'e/4>,#OH'{:S#-'###rG##B1Pw+###`?#sZ(&6&1,#9n?l,&95#A-#e}G+jD|P$<u#R1,T[URm',Z$Rq#@[UZK(p`AcY4######<H#FV;ZP####]Q$ay1L^U###vH%rZ$Y`U*['c>?p-("},
{130.37,415.31,4.77,4.21,"pb#5-#IjA=Q$PG#IZ$Pc$g*?e5#X-)j~#2NB;[&W>$#6#|Y$qJ3>S#M@-Px#)x,e`5'S.u=63k<kS0`|/YL9}c(###&P4[l%C?Qe5#o%0{-#m;@.K&p?Q]I&)@Q}[&0h5|S%;Q&|d#DoO5Q$M4=eG$;-'G>#1L/[J0eq2:s=uK7a`,zQ)0L,Z7.|]$H[+/@#"},
{455.08,415.56,5.07,1.67,"I5#?6($##.,#[HGk>%######rUVD>####$##ug1XZ'PG####3,#:5#Z$#Rn/qaF$##&##A/.LUVyl'95#wl#up*,;=eY#$#####*##Ej,]RVKd*95#_[$TTVukHeJ1a,%ML*n@(Pg3p?*A,#G>#(9(X9M;WAK>#95#r'(7RVg#$|R.*8#DRV0-'A@+=Q##O;"},
{239.96,438.46,4.35,5.49,"GA'qL7U,%3,#@^/Bc%iG$zP#_Q'?5#o%)Ef3.,#.,#'$$,PK_/(oQ(9l#Nv(,xNZl&###_l#]x1oS2$##W#%(?#m-+Wl%},'Bu$###@Z#mmQHf3D7)i?)n10c-(hh<eY#jKFVQ$gH)y>#BhZ/,#<5#n#&g@-aU;c>#+T16j)[rBx^//n->2O?]##gZ`['zgZ"},
{284.56,442.22,4.52,4.45,"zz?:,#cG$Ow#a'3eE9jn/$M2}`8Vn.yr,NJ1?c%*##Qa4:l$3RQ)##Pv)f##LsA07&iSQN$']SQ6[%n`?sd&BQ&%A#NRQ0,#Sj7B,$<#$%##041K&3|^.m:8Py6kb6>U7|i4qb#o]#s=K1,#.H#+c$UG#.,#FJ#m)C=u#Q5$],#RmG7#$Ou$###wc#yc'P,$"},
{241.02,455.79,4.52,0.44,"26$?94F>#vY#B.)oS,Fo3###d+[TQ$2o0K-#HPCtz'urYmi<Vv%t4AZP#)##$+<me0hv(XZ%OXW[-*0w&[w-Eb/M'[;S&[A1y?*,v%4?&dJ*2$SOG#)##_x$~7,Q5$?c#3%+8R%t98^P#V5$/v(###1,#Tq/;$)RG#:5#4K&K,$.u#N>#mP#1u#;5#3,#[5$"},
{141.70,462.86,4.49,4.59,"LrBS5#.H&De#9'0Oj9EK4Y/,t=?(%,^V(`R,yY$(##j*3<Z%Q[P+##_7-?##s;=PH$K]PBm'~~P<v$l:9c7'^R,@n#)[PI5#7P:ZP#X5$%##Uj0M&3kU,N2<<'6+56y/39X6Bw-?g#ez>+##oc$tP$OG#'##pw#||G{k#gY#?l#-eJ+u#(c#Zu%3v#6#$1Z#"},
{126.41,467.39,4.52,4.33,"{k#&##:5#kZ$[P#>c$Tl#'?F]5#wl'Rf#7sE-R&L5$/6#m>%Wp9*##6#$Pf#n80GF:`@->b5n`7j960|,Pp7w#'###Ot3l#&F6R.##).,>$#];?R/'s6R0m$<7Rs$'8(6X%%1v((@#u8RS#$T>DOG#AZ%;5#W13g-*LE6kD;z:>YD*,~-^M0RR,9g#+06Z##"},
{21.00,531.16,5.57,0.00,"e%0######&##KFs######7##4~WO5$E>#@,#<,#`#%9#$/,#?o3######%##AFs######O##d1gq>$F>#K,#$c#=Z$iY#~P#X&4######$##JFs######9##wf^m#&###)##&Z#D6'/,#0,#6J0#########BFs######+##]uO}Y$###c#%/l#sP$/,#_,%"},
{455.45,575.64,4.29,2.94,"hY#*v$n~'wOH77+H>#x>#&0Ky(;.,#w##w]/)@(]P#mY#>I(Iw.5,#5)9C[(mn[###~e&G6=y%.$##:x@]'6%y.-u#z,&I>#o@.$##9z62,#Lp[###b8/R$%$L7vY#Xr[JH%K7%}l'Q/2p,%sZ(/Q#bx3$##Fq[*m'@w+;I$[q9vx+qr[X$&M?&^,$CT2gI)"},
{455.45,575.64,4.29,5.77,"(c#bI,R$&kY#E[T7m'X5$FA%NZOa1&)]TXB.E-)H$#_+GNv$9?'5,#TB0{Y#g[TWc$am*gS#hh=yz%oZTS##oT7-$#b]6c,#v#&s>$PU8G>#)UJ$C3u[-4##g~+}+4+QQ)##<_=u>#-d)W,#]G#Hp1_R,;5#Am#&>AR~/.,#kw&HwP?u$,##~_0u?*TG#k>#"},
{120.91,588.89,4.17,5.04,"3m#cG$2,#.,#f,#w7/.,####GQ$)@+###fP#R,$?u$(##R6(6c#dl$fQ(}k#XZ$V?Ek@.?5#KKT91:;5#87%=96###WI%z<>t6&QZ$i7.3[('v&V@'NlC7J.$NTQ.-ZQ&O[$(@FZC9w.,Ov&NQ&R,$~c$rqQL>#Z$%$NTEJTk2)Af1N[IFQ&:40n`DD>#Y,$"},
{304.86,605.67,4.06,3.08,"1Z%(##27*TZ#yH)###V##h;5D@'D>#J$#p$+[l#&Q$iY####R]5###D%(d$%H%Y###Iw$6m;HQ&:5#7A@>'5@H$&-'bR(eu&DT3###wd'0l#E(Y###$7'.[%ky8-l#*+YmJ0=,#w$(7X5j:=>H&(##)~,###-)Ysb#^l%1H$$nOg%)])Y9I(tb#QB,Y>L.Z$"},
{304.86,605.67,4.06,5.67,")15Kl$)m'4d%/~QOc&|G#=s5[[Td('KlM:{*+u#S&#u[Tau$Q/.%l#xQ&.c#B]T&m&eG$b7#hEC.i'oZT~##gH)4&#t)DF##)n+A5#N7,K>#y_TD:55l$H##g[(H,6AM?(##</1Z-$Uv*)##fP#&Z#Pe/###TI%zC9-6'###R$)|h395#.##Z]5.l####7##"},
{141.16,607.48,4.43,2.87,"SR,$##wS1ZG#'/~###jm*C@$-S/-##C0~;H$$d'Eu$k~1&##X/3D##;p7K,#F0~#6%(K1aT$IA057%y0~~6%P>#-H%fx3)I*G7*>7&<y7%##0WP~0~aQ(H,#@A'GXVwtMXG#l6%PH&kc'U?(%##^6$H172[)9[$zT5ac%;#EOv%MK4im$}95<].`l&*c#ul$"},
{141.16,607.48,4.43,6.07,"Ru#G$&iA.IH'5I$7;99m%3'5X?'da>a?$ny6#;:36'%##2v#GH'qu&s-%X-(8,NL>#$/'S4Y_?)@,#[{PO'~n':%##ZI*)R%<&1dd+L>#v>%f(~R?%fw.*@%ZT2j~$H'~+?%<g7B,#9066##]R,$##qc&Vu%t'~k5$Td+(##&A-1R$;&~###KA/gP#yw0###"},
{154.28,613.24,4.60,2.99,"######87'(c$[:<###-/),w%KnY###_y/l:,Wc&###adH:Z%OG####mQ'}P$(YI###2q8]G#&oY###c`=8.%GZ&aP#wpYxG$D>#9,#I6(?u#|`BNu#vq>=,#hsY]q8W19S$$7m&mz/+pY=c$###,##Z$)nP$MH%XR&KrAOG#W2+moYic'w?(<m$7HFV-)tc&"},
{385.65,647.14,3.82,4.48,"zY$######(##n=L######]##;{l######Y$#<V=0u####U,#l>%######*##Lpd######V##8{l######c$#Nq==5#/,#=-$k,&######'##|(k######J##I{l.,####W##j1<_l%.,#(##c5%#########Pf`######S##<{l######*$#(]2Q5$###,##"},
{54.47,656.75,4.86,1.34,"E#$Kd(<u#eG$ic%Dm(C7(^u%ZP#0J$UaAxu&`I*79)*@*m7+rl%%C.MNAJ>#y8Sfx0WK3P6%gC:FQ#>8S~:-XT2j?%nwJ>^L?c%$?$$uD{Z'-7S$-%e'49K*Bt;k/-B7S26$HG5rZIw.Si.)>H#?A0[u$Pl%vv(|d,s5%X6'OZ$J.+Dj;06&[/.<n+j0/+{/"},
{54.47,656.75,4.86,4.47,":U/O;02f-Hw+Ma<??&[c$4w*}>%L-'mm(-n,Dc$Pl%EH#TS01ASv7)uk5ocHG@SA?$h4<h&-m04@K*4@S&-%okD,d'?c%/H$t7J9LM&02V-%@AScC-M1:AH#_T3U6%$BSgx0>3BJ>#ql%h0.z6*Gw*R@*]T)d|@ac&eY#,J$5.(Sl%vl%f-)<u#Z>$Q,$Xm("},
{272.93,658.34,4.14,4.51,"ko5######.##g)o######{##ldX$##1l#IZ#ub#3,#8Q$oP$o1=######(##q)o######d##X:gZP#,##Tu#G5#hY#`c$Pc&I_>######R5#.*o######P##>paZP#&##r,$fG#=u#MZ&0u#Wp9######B,#i)o######`##+@V)l#kY#Zl#/,#Su#ru&]P#"},
{296.18,664.55,5.17,4.54,"]P#######$##2[T######,##HXq.,####6##yw/8#$qY#kG$bG$######5,#n/`######b>#?Xq.,####i5#a]5O#$6Q%Lc$8l$#########;0b######*##@XqrG$###0##Io2Yg0pb#%##.u#######%##juQ######/##:YqJ?%9#$*##-R(%;->H&E>#"},
{620.95,43.08,4.96,4.10,"jq7N]0qZ&>A,5M9FJ.3R&JB-5T4###$##4e#G-)I,#A,$_5#eo0U8.)XA>%)bLZ+o/td,M,#%aE^G#c>$46$0H&x&$(~.3,#)z4gV.PD>f>#,JZJQ#5%-9.#siDgQ#jZ(VZ$.,#52#`PO:5#k?)Ie(Jc$AX59vT_>#95#E{$zd-u1#,z;_Z####*)#wFL###"},
{642.08,46.74,5.15,4.49,"kA-U7(pd*BR%4*4zK1I9/E#$ABK;#$ub#>5#,Z$o5#Gp8I>#{</>7+I6&=H$*W6[x)xmL.l##_RnP#HH'B##8&1/$#xPP`G#%])E6'Z#$C:7Vf1kf)z%/))1iZRX>#xY$kC$j~1X7#IZR9##CI'6#$O,#BuL)Z$:5#yP#R^Raw/Q>#(c$MC'6l$j@#IZR,##"},
{503.37,56.41,5.20,4.98,"^C%](>######W_)8OH######;AE<c%7-'###pH(/,#yyO/,#M`9L5$###R5#&PB9?'u>#hS-bxO95#B5#R$&<C7###TfN4,#MuO.,####I7#|FL###@##8:OkND.,#(##&??oB30,#[L6:##2uO######,&#KuOZP####C&''6%W>$$##7v';I':5#Wl%%##"},
{21.09,59.83,5.40,0.52,"############pb#######)H#{k#######TP4###$#####~pG95##########/EC######%.#)0~Yc%nG$5&@#d'WA'%f/,WP95##########6[O######1##25~tn.gY#LH#Bx)3`2D6(wY#############=e*#########;`POG#######Wa5P5$###'##"},
{513.61,67.47,4.92,5.28,"Tp'uf6######F)2`%/I5#*.*SbYpb#'##md'VX:95#Ve'1##n29Fl%######V^_###*##w36>XE95#%##?hI=h5L5$H>#Rm$0i=-u#95####e__######}u#n6O9Q&###`%)Wu#`l&###,l#_B7$##/,#`,#X]_######(&#zXJOu#cZ':m#jG#rN<&[)###"},
{98.55,79.63,5.18,4.85,"B9+|#&Ie(jJ2iUK:,#d,%~P#1V3###Jp5$##?#$###Z2Y###O%)/0-c@-yY#=1Ynu#ZZ'z##kdHHl%am*(##7['.,#g2Y###XH('m#cl&-C+@6Sj,$mP$}h$BPCuI+K#%B##*C595#71Y%##B%-tb#3,#pg,fZ'@7({k#$x$sg:e#$6#$~$#$g5/,#5%R1##"},
{62.67,84.75,5.20,4.26,"E^/Q{>A#$#@%-C+3S.X5$3?%$R)4,#yl$$%*%U7uP#eh*HGBCT-='4c])_H'vhV8v(9?$^l$[_>$##6H#-J*}}Gz_*4:9k-%t]1w@*|/->Q%]dVg,$mv(7d#>OH?##<?&Gd$[Z'1M#FdVy5#a[+H8,D>#M/'LdV_5$95#p9#J_>`$#,o2k##+u#`(#EdV)##"},
{616.16,89.87,5.50,0.46,";d&6)/iS/F?'fe.bx+O,C&)9$S,W_4r<DWe+rc&+?$WC7c:759)3<:Jf.W5$w*?]D8zU8Wa?jZ'Fd(Dd'{9TF>#(-%Aw%X7TJm)R5$TnB<l$B8TZ>$:$'&.'22<yg)17,`p4Xr=:O4$n*/J-^>$$###w>8Q&hQC;Q&a8(MA1f]+aL16f2<@*5q9hQ$R.-m8)"},
{616.16,89.87,5.50,1.80,"<-'~u$To&1sBS-SRG#;l#uC+>;=?/0.r7]S'iw&>bHa7*IQ&PK+@&.QM4pK38.Sxc&+~*)T*o|D(K3h)3%6%GI),e)Rr.1o1Ef0kw*I1.44;N`8QGBF14;&.#0(Z-S.l#l#&k%._dG|Q(C[${b#%w'=b<)Q%aG#0/)4/,cE<vP#g@.2##5/Sgc''d('##qH<"},
{168.35,93.79,5.20,1.16,"s$,%.#Kg9j##iZ%dABNy8%##VZ@&bF'l#Gv$tl&1c$Gc#iH'6#$:$#&,L5,#wy8ey(~gYaQ&yhYf7+f%-SC)h~._-&(K/07(###V5#<{;###zl'_u#S*Vr:97gYu#&4R%9,82&/xl%W>$=c#%##&)5(v'D>#K?'#q+:^7H#$),I=?%zG%w$$)I*.,####R6$"},
{118.20,95.63,5.25,5.63,"vD'+3B&q41H&z[>m]4*?&###4W*A:7?~-pb#X;)zH*(C#~./ko/S%+Q2<0c#%2WeP#'%+0[#m#KbH$-q:6##q>C?Q&sb#'##J)73c$|#'u%)y0WL6([,%r/%xE<'9017,hu$X|<rH)D>#F,#S@*#m#BA-Q`>QeN{>%Ol%[''E&.g,&T$##<;(8'6@)@d&6@*"},
{593.98,106.95,5.79,2.74,"]H%tV5Gv'm]16h')&1t?(L-'TV,{Y$0v&XG#VG#?##ckF###OZ#c7&eI(-?LTf,sp&$7Q`[*k:QPc$~H(I,#O.-I##67Q[G#{/%_,9Y6(`l%jf5jz&>n.S~#96QO5#bG$9:#~e/{Y#17QC,#(.)=w']N;2K/XYM,?$Rm(3p'.jCB,$###XA#a>$tP$Dp4>,#"},
{40.90,129.06,5.64,0.09,"CKY95#######kW>07+~,%aP#|G$(K'BvG)H%8T3g7'/&-|-&oKYvb#95####W+AFtB7[*<5#G/,6mHW^1hR*>c$Ce(A^,>A0#KY&##dG$%##}uPSv%_GM,##cJYnl&Sm(FJ*K?'{G%)B#){>DJY###H>#(##cGN###,O9>u#rPB8@+g@(.u#,n&wXED5#k5%"},
{40.90,129.06,5.64,4.85,"###0,#%##*xW&##B,$###}xW1u#95####xyWD>#######@yW###ON6F#$IvRW$&+lL(##*vOotDd$+;5#y3=wm+/?&xb#?V7g-*gn(Z>$t+?ac&K[)Sw)5xWG-Hzy1,~+H/-L&'+nIn#&Il$t=EA5#?l$l[&@Q&?8#'i?]6(X7(Uy+r/4i>$@%'ew+m-&OT4"},
{369.80,159.04,5.72,0.44,"###$##W5$=#$.,####2,#/x*###+##2,#V6>>u#^6&iY#s+:xu'###lP#Bm&E=I###<u#xt5dc'$$#TE=FBIM>#<x%P.,XK2}v+###Pr.4p6:VRVu%&m&k-$lq6rq0KmGT$(###Ka/9w)Se/1u#V,#HUR7-(&456]308.###%J&AVRm>%######BM-O5$###"},
{369.80,159.04,5.72,1.11,")c$###&##4Z$kNF###$l#Ep#_x42?$mh=nJ%II+Nd%?f3.I$tm,###%##RZ#cdR.,#Xl#vh'D/0YG#+b2U|;ZP#S%#=dRXJ,WR+###+##*6%z:Qv6+4##uc%e;)y[M1h27@+'##pD)(dRk>%BH'###?###}7@R#8Q&-L#mbJ}%#X}E;7$g,&&##jK,JH'rb#"},
{313.55,197.86,5.05,4.83,"GD9(##^B5###qQL$##%4C###?'S###.,####92:$##>A/###P?L.,#R@*B#$*XB###cmA5%+v&S###4c#-@&m}D/,#_6)(##/'SmP$M5#3Z$M)6by9n0%ME?5<A0OE7Q#li2';0d`A*c$=,#C](wu'######nS$x#S&##[P#q>#J$S$##UG#N-#F$S###hY#"},
{270.85,202.78,5.42,3.28,".,####$##RC6?l$eY####A?E=6'.,#$##_1Q95####j#$B/Q(##FE?###g:7#^,&mO###hh4_1Q,d))##w233w,IQ&,w$y07###{K4/,#G/QU-(hw.###_1Q&.Qmv+###kX29R)v_<5,#4n+###95####B-J/,#Xu%###`?H*l#H~.###Is=1,#u^8###.U3"},
{270.85,202.78,5.42,4.28,"Th>######a$#{ZQVZ'$##BV&UI)-R(M?(+g.7m)[6#FT5hc#NZQ5u#OG#,%#n]Qg~QCZ&+$#&6:u>KYl&<,#|y9*%&<R+4$#a3G7?#|Z)0%#f<Duz*~DAd$#5~Q&d'{k#B&#ZrBml#Jm*`$#1T4{##Zf43R#,5NO$#*g7F&#6ZQ######j&#px5|-#9K5w##"},
{615.84,266.39,5.33,3.45,"NK6#########3Xq######-##;@T95####V#$OZ&[P####Iu#_U;######$###Xq######v##5]].,####<@#9Z%kY#PG#$v$Y^:######$##KXq######L##G:dOG#%##{5#|G$Dc%hY#N>#S]5######%##9Xq######T##qmU95#2,#i>#nY#D>#5,#_>$"},
{299.06,338.12,5.52,4.39,"GL7;Z%S5$3~$$8.oB1.m&l[C8S){@.Y(&Q)ABc%###3''hl&3?P'##x5&t$#/>H2L-g@Pt$&HAP9&,A*:hd&2H&NR#u%O;#$H#Apb#ub#@,#jk6ro5CE3uy6?z:C,89g4zr6}Y$>^$$(;%##{#$X>$;5#%##HJ#aEEO>#)c$`>#;[H}k#g,%4l$[Q$X>$B,#"},
{288.97,342.06,5.65,4.29,")8-_,%Gl$t['zY$d#%6-%D#EM,$|,&?T$}{C~Q'W>$u5#ac%msFQ5$+u#oA#^q:S|8305t+6M4>ly6R|/Hg7k,&&##oX3?H&FvO$##Ac%G-#iGI$%&HxO=x,>vO|~&Y{;i0,Jl%Dx#:QKSG#_'.W>$;5#&##o;+r{B*%&@S.R?'G.B;I*R/-8#$OS$)R*'##"},
{135.85,361.78,5.33,4.22,".6')I#vkC{#%fY#*Q$l?&)>DK,#Gm)cx#XrBNd&3l$;?#pP$FL:#U$QR,Ue#:B1kO9,S.e=7@O:gB5m30M^8bH(###`b67Q%XHP[l#Sm**$#/E@z/'IIPx%**IPq%'K17ix&w5&-~#$JPP,$I)4h5%=l$1,#hL,c07)B+f:8/J/zj2NH'if,`>$_o$fH)+##"},
{246.89,387.01,5.11,3.96,"xu'k?'I07vP####gn%H.T},'7A/q6&<-':y,F#LD>####>]#1v$K0T6-('##,:2w&3Nh5:_;%m(###:Z#;1TesH4##%Q%bi*;J-z'2*?&8J%o.Tg,&###/:)n_?:5#F>#yL4,0T-Q#-e-%]&Ch<SG#nG$oE-G81^5$TR':/ESm*Sd)DH&~$&4=1+S,1m([P#"},
{246.89,387.01,5.11,5.43,"wA$?kIE>####3=/]83+##|k#qD)wRT:?#Bw-eG$PxE`-(Kd)+K-q?)Ie(;.+kTTCZ&###Om$:W9a^:/,#F5#M?$h6)OR*q,&uG%###pA'RRT#^5vq6OB4jU4bD,ZRTW>$&##rlB6@+Ac$6e(*Q#gu&$&.T..*##{-&:TTT,%/(*<w+Zn+.Q%O#7+i?ub#N#$"},
{140.60,386.83,5.77,1.81,"],%tb#O-#BV<`5%)##W.(gZ&I..Pu#p#&i+6Jl#AI)+%$[OERL*C16V'+]f4(h7.,#9l9-v&.lNtZ$*r7YA&7z-A002b:+R(P|;M$(g|1Ug2O*Az$$jmN/7&~pN)L7-:2m>$6o)S<C4|.CK2Uv#.mN([$d~.u,&4,8<n.*d%5[$a81HH&7,#z##|XH|>$_#&"},
{551.00,395.35,5.48,3.61,"A|?f,%e-*###II*4M5$jD9##jHB$T2hG$aw%nT%xm,Km'F$(4xU2,#c[,$##v097,#~wU[>#6wU###F6'?K){m*X>$5?#e4CXwU###Y?'k5$II,###2>?29.<wU###RZ$LQ9/)?.,#0l#Qf)zwU###g5#^H&iK6###Yw#Ai6qj@cw-Y.*jeBI%)SI*/8-i-)"},
{431.49,428.20,5.90,3.40,"O>#gz1][,lY#gB#{lOOG####TV&WI-###I>#f?&eY#,##W5$###wf%OFGL5$.(.iC7hf0M#%jWV}Z)H,$ud%@'1eY#E##4@(.H#1VVIn-pb#s`>));;93i,$lSV2,#Nu$R.&OXB###9&#U&1Mn*g|@###&##<TVk6*###LJ%hRV###4##i:.}7RL5$j##Cm("},
{138.07,439.01,5.87,1.76,";l$0,#*$#<2<cG$/###%&uc'hv+sP#1Q$l48+l#bZ&iR#@PG2T(gU4FC-nL9@g3sb#y#<xQ'V5OV$%W17ly&A:.Gz3<=;n[*|h5=%+/Q:?x.gr<&7$g6O5I'1:O*&0;_6bG#Wy-sp8b47j]1Xd#P6O6@&87+au%`}6e@.I$'7n$Z81IH&@5#X$#{5OLu#'-'"},
{558.60,452.45,5.36,3.36,"0e.######/##n1j######x$#Yge.,####,-#yl'95####Ru#gB8######5##)2j95####/%#.Dd5?'###z,#_m'L5$###Wl$1_=######/##43j_5%###W##>6Anh?95#&##^p$NK6###$##HU:######6##'2jE>####p##|o3pw+8w-*##|A&Y+Fwu'<##"},
{438.26,526.47,5.06,4.08,"A[)u%,r0.L96c'0Ys@ul'G##+4B8l$ZP#_,#IQ%2c$pb#H5#=v$L[>XURE]4>SR]o*Y97D&#>tH%-%A,$h%#eZ&Tu$cG$T5#sZ&^l9m|8Y~/qSRD$%c6'8D/pr=Mx,D>#m,#7S&uJ1###$##eY####uR#Y`Au$,1##?m'M@CD?'yX2LiBgZ$9P8OvJ#?&1##"},
{291.23,547.82,5.33,4.18,"`5$`n+FH<5I+U0,d'6_#&Mv$0x+W>$*##jH&xu&eY#$##;5#<I)8n?BUUnc'CSUwH&r-*ia/l7..Q#r{80a9;?%Rc%jH(E>#hd+Lp,EL)D93KTU&7(9?%dy+co0^X/(4GZ,$OT(n-GDv)###eY####$$#Ty3c..vG$_u%_T+-Z$Mp&USU/c$8|>vf-Fp4aA)"},
{486.07,582.41,4.79,4.85,".U(rF=h833u#:oDh=J###z>#ad'(DU###:%&.,#cE;J>#to..y/2S'<=7?XEvBUdc'UH$0(.Rf3%n*|Y$931Tc#@I+6,#u@(f(2jn)50.yA2|@Ue?$FA-=~'ksDpo+$3@7H$3M.<^5<R)O#$`25i6(h,%fv'{cN$7*xP$.D'yS1yR-g14jQ%We.4H&&e'E/,"},
{141.97,589.79,5.68,1.59,"'##0Q#)7*###(##'M4N?(###^G#U(;%l#g,%###[~,'l#gn.D8-GZ#w^5qZ''~+C*:GICu,&F&U^^90Q$bI?>v'Fq4FZ$](UcHOIQ%%J(U8+Z7-:w$l)U9@*l)UrV>en+SS'?A)O)Ut#&?n+#S)y@,b#&'Z#=u#-o':FBX[+Ce-/Z$l?'{f.t>%8?&iv'Y,$"},
{141.97,589.79,5.68,4.17,"c,$^,%`e#vlRsl'jQ%n(/[h=_7(K'2'v'eu#'/+L5$(##O?%%D%NmR=@%?|Dk['h/GzpRDK5YmRA~&?o2Bi)W@,Ml#{n+IE4y?#[:;rH#f%00I)$k7qM1MT3wnR-v%-v%mp,k8.jV+iT6l,$N5$###B##7?&A,$###9$#.:72R*;5#5l#>1/;l$5m#UkFl>$"},
{423.83,628.11,5.59,4.45,"E&3######2##O2m######t$#XnXm,&###M##wc%xY$###4,#Hq=######K##N2m6,#4l$5%#wwUKT.un0u>#'q2@Q&E>#;,#:ED######[##<3mN5$n>%C$#PcM0[(CC5/R&';<tb#WG#0-$O]5######{S#25m:$)###Jn#Gs4VK6wb#dQ%4T-t@-.,#[>#"},
{146.87,631.49,4.91,3.13,"fG$;c#{f3@Z#N%,F,$nG#;E0y@./,#M?#Rf/&R)###%##.e*n[-###=g*X.(K%Y###(n$#R9q$+###[n?(C54~*&Q%Z6'<5#fd+###dA),6&<(Y###XI()R%H`A$##*+YJ6&0c#Q?(Pr;v>%fG$$##NJ.cP#5)Yfc'6@*dG#e$JGW7n(YIc$[>#Sd&=OB+$'"},
{146.87,631.49,4.91,5.77,"e>$2v'w,$Gu$$wVp>%?5#Ny%:GLq]$YyVtW6SH(S$#,mK~m&d$+>,#k.-uG$FwVxl$#-'@8#8<BaD&[vVU##sJ3:7#=L:[##,.*a['OR,$##V{VTj=Qu%6##5S*Jd;;)A%##ML5:-%VZ',##X#$YI*;6'/Z$K-$*4@r-*:5#=/0Wr6;5#]5#`M>$l####A##"},
{162.78,639.88,5.49,3.57,"/,#c#$RJ0**4ER',M9>[*2-$lT0YR,6c$VA*_~,95#B,#.n(###+##}r5t,&IL:J#$y_/%e(KJW95#9$#39Fl@*ZP#'{$<+GeY#$##9S$KZ&O[+###;^#8D={NW###Y6#J]0:LJ^P#hr(VA0Q5$$##gP#],%eY####j,#yv+3B*.,#{~%KZ&t%A`;=H~):l$"},
{339.95,644.18,5.53,4.54,"H^9######)##8=n-Z$E>#*##2z7K20+m']P#[5=.]/`P#F-$Yh>######+##5<n2,#`G#_-%S_>;Z#:y(Pi?}m%8Z%CZ#8-L,M=#########_<n###%##$Z#3.W|k#_##[T4cP#=#$0##yC:]T6######*##|;n######6##o5PcG$###B##<#$$l#H>#k,%"},
{180.51,673.65,5.30,4.62,"y1=######&##*+jW,$-u#$##IT2&C(I80[P#Ks1In++Z$1H$4h;######$##G*j1,#/l#e,$Q;@a>#nC,x97K?$[>$3$$]?M8y6######&##O*j###$##T>#$8Z###i##[B3=#$###2##d3?;e.######.##Q*j######(##j0a|k#$##a,$H?&{k#$##:l#"},
{129.69,681.06,5.99,4.53,"|S2######-##@Wm###:5#W7#sz>3,#3m'<7$21:&##^P#?.#6C:######'##(Wm###9l$a,#n[WB,#SK5,?##tG:l$[P#v##5C:######/##%Wm%##0c$M##OmSVH%G^3kl%2M6n83K>#[m%;96######*##JWmY,$rb#3##7M;Rz+][*xb#T$@E~-_P#eI%"},
{44.71,48.41,6.68,1.50,"95##########|cT######$##dTb%##],%)##Jm*,##h*A%##95##########TKb######1##aTb*##HH&HI$Km*1?$;q6Bu#95##########KwT######0##Z[bs'8d$)(w$Q[&Y+EJL7d5$############}-'OG#######ya)d|G######7%#%UbD>####"},
{560.52,61.07,6.03,1.28,"}GQm5#X>$u%#/u#Gg,9J0'##R>#R/,D[PB$)uG#]d(yw(+NB=lP&@#4&2u$#Vd'4$<^QO+l#jJQx:8_d)3K-gH)###)@#G?Nh&5K##.:8w##,-':Z#o]Jlp6THQeu%=R%8u77e.By1.##w#%.6'NQ%<e.H,#)##3I&L@(%d(~,%l@BvY#5%)~P#(MQ###(##"},
{611.31,107.79,6.32,0.19,"eZ$t%*+>B'l#8|8YA.18-S]-BH&Z@)HT+Z$R%l#~u#N1+-$Rnl'%##}&RhP#^$RCl$<//U]'9(9#}6@-%~=E,T1ax/aK%I2?&d%PG#Uv=kH)ePDrZ'31'E|Cox/7r/3e,=)::5KK-'WZ%m/)N`,;c%R$$eu&'?#a5%[n#kWEQG#;5#:6$GtFl>%L>#]#&{7)"},
{97.30,131.92,6.39,1.50,"E5#Yh/Q@,)c$Y?%i$I[?)Pl$R0&35N3##s/2~q%)]33##/6'zb#)P8Nv)&$%aC1p'01@+AZ$^oRdG$m##{*@,T,{k#4%#YmR*6%--&J:/P/-J06Pu#oU/7~'OmR1H%5,#8h'Q16u;==v$IbCF[%lI-$C.GK2P%.oG$NJ*M`.o4KMI'8.-pU#r%)&.F#V<;,#"},
{638.50,184.33,7.67,3.41,"A7-######(##$tu######E$#<@W$#####r##Ic%fP#rP$Q>#X83######(###tu######6$#;f^;5#5,#3?$<l#}u$)m'Y>$696######*###tu######8$#j@XW>$'##}>$qP#L,$E7)zY$<]4######(##|su######t##+OGGZ$qb#Hu#Ou#-o/XG#~P#"},
{631.41,209.23,7.17,3.44,"@A1######$##;tt######C##o6T.,#$##cG#8Z%.,#$##vY#$95######'##|st######($#.9`######r##QZ&J>#$l#T>#3'7######$##rst######L##]A];5#1,#<c#m5$Lc$`,%8#$cA3#########mst######?##buPpb#&##Au#G5#yb#*m&*c$"},
{470.59,336.24,6.01,4.29,"{0*.cP######1KKS^:###hc#o25`$+$##=-%i/(p-+######$E?Nc&&##/##>]ZiP#kY#9S&&L5F8)zm*X8,{C1kl&PG#$##ey9###.$#LT0I]ZdI&|,'Y[$'i6%7:gf5P>#XBL;n,95#.##Qu%###t.#=]Z2:9Uu#&~)Lk;(S-$%%wj99{9$?O^P#@5#&L)"},
{221.15,350.21,6.29,4.33,"~M<###/,#'[#uH(LH%gu$kF@t-(fc&wT#<HP=w,{k#f$#n6)>HP###7#$W/$W>HJ&+.5K,g)uJPaw)8427_6a#&Um#&.@}H)>s>7#$pb#gd$Rb:,&19|<W0/oHPYZ;T&1VC+H,$(;'a'8J>#g6(=Z$:u#YH'3S$0<A],$sc'Im&+[KcG$9Q#N5$WQ$rP$Y#$"},
{580.52,378.03,6.19,3.50,"`B7#########3ze####l#J,#6T4+##,3@DH%>QR###a-)0&'up:#########+ze###=,#j-%Ah=###x~'.F:L^8d>$x#%.9BiT6#########3{e###%##1c#8x[.,#]##(U2yT/in.m,%D&-%S.#########){e#########vBb###$##:Z$}5&ZP#-##c-)"},
{454.26,393.88,6.39,3.56,"###Pk5Bd*###iH%9GBW$)###%ATsu&$Z$3##A%P###6##C6'$##1H8lZ(###cj?YC6###2##o@T6#$)##h_)4,I.,#n5#CU/.,#t>#D{<###bYI=5#fY#D##sDTxY$8,#<e&&4/$?&w@&A%+######.5H###yR)###qP$8Z#_oGD>####n%$w$>pb#wY#^l#"},
{508.22,402.07,6.33,0.18,"AU9^5$E>#D##;y]{k####/$#&kD+'2+c$(##vY#K;0HT5%##JsF%##O5$5##%y]X>$8,#~l#uU6](;lR']u$f5$zh,A`B&##&95$##:5#.s5Qx]###)##k:(,<C3l#mD/RH%cG$%.$v49}@.{k#######>}]@S0######QKCG[+###'R$].(######DB$]w/"},
{248.70,412.31,7.05,0.92,"lP#v-'Mo12Z$rY#9W++`?c,%KL:K$$nH)nx%5C:?#$;,#K>5`G#P}7,:9###1x.LJ-dU50M8n+JB5#h>$r-=:~T-f(RV=zi+XH'df)[h=0-$(.ST6(gY#yS%5HCW7)5[*}l%A:M>M:$L3b$)%f0K[$td,]m$'^4:@',6&Re%0.()|5dQ(+##4<2>'6D>#7##"},
{248.70,412.31,7.05,2.00,",&&D6(wv&Hl%fLPyu'?#$7,#owDRHPb;-W80CI)~dA;a>b/.4H%o>$hl?rm+QJP:Q&+?$&A(n21_HP&R&]5$Z,#/LPv..qb#D>#^5#eD3;S/q/-PL4y81zR*?).m&47d(I>#7l#[$'Zs@Il%%?#u5%n+8hn0q,%wu$Mt:g6*`z0@e.`?&Ll%X'22[)=R'bm)"},
{539.69,417.58,7.01,2.36,"#########68D]c&]>$###%fEf2@OG####OF1TH(###*##X&,###Qx+}G%Nl?i.*(5C###wo.zJRnH)###fZ$+95###*##V?&.,#b7',(,(HJgH(E.*`$&UKRhHRB6&/c#O<0]/.~05###vG#_,%n,%G8'/OCSu%vZ$nV3Ig5jS1M#8sc&F6$&A#cLR###%##"},
{451.15,535.52,6.90,1.60,"&?#<-(0v%xY$%##$w'%%+.,#T5#>]1;Z$-c$~P#(J).d&9~-}wDkx0I.,EH$R.+C-<%kB@#$ZnQ?=D.Q$K$=Zl%n:1?@(%pQ{M?M>#vv&F5<lI-|5${IB&J,XqQ[V;'^+vU.t.*XqQ&d&?x-OG####%6$t7.C,$###y:0I-(?m(gY#ex(|^3Uc&xb#F-%=.)"},
{132.07,573.53,6.56,4.31,"~$)%d&@$&6n-AT+aS/qP$]5#gR)D>#)##,6%c5%###$##)c#^H&r<34:PrrC*7PyI(/'5sz'Xw-JI&Kp0'{2*d$0o/q5%d5$_$(w6?7*=gJ1n8PZd&?.+A1.o&/k30;6P3H$r5<LRFUf3s,$L5$###7?#0^5=w-0,#Bu##;1*?&9Q#GSH-e,s08Lu#3q/36B"},
{22.55,579.02,7.21,6.22,"7m)######(##XWp######E$#1hfJ#%###h5#4l#Jc%95#0,#Ce/######.##^Wp######e$#A2iCZ&###&H#VG#g?(ZP#<5#]x4######3##]Wp######c$#d&]L5$###:##|P$n#&1c$;5#AB6######3##]Wp######@$#$OD?u$###1##I5#u?)D>####"},
{472.52,600.24,5.86,2.17,"<#$f>#AO@oV>-?&r>$=U(]^8Do1gu%>R$n=@(v&Y%(vdB}e0j$);u#'/IFu#;w)=FB`=<1[(EYF,B2SS*|n+}G%Dv$sSK.v'2Z%@$#X;;Fw)<Z%@6&Z#;8$NNXFX#%~S(W7JMl%E5#X100tEZl&LQ#Q?(;{*+u##.#KZB;@E5?'=##~f-$<9F>#ub#,7%g06"},
{472.52,600.24,5.86,4.17,"[[,###A##ZY7e}D}>%<5#k1)B>=i?)0,#-H$2@(R~/>5#m,%QT0g~-|@)+)3a~,_QBd@+,o't'OuM;6c$h?$NbA*.,R,#Z6CccEY/10u#-H#ST3/J*5|5v~)g$O$?%k$'%O0{t7mq?,##;K0(%O}k#(##{6&8S.G$(h,%ds1}m+vm)j5%X{+l?&u#'D,#$f1"},
{382.63,636.26,6.02,4.50,"So4######,##3*q######s##H8]D>#'##9-&zb#[P#&##sQ(&V=######4##5*q######?$#{'e95####&-#zP$$l#'Z$Vu$1{@######+##F*q######P##8B[A,$###.##1$%xG%PG#0,#dL<######0##=*qF>#.,#^##W}Ez6*pb#)##Td&X>$######"},
{108.69,685.56,6.45,4.64,"Ki*8R)W>$'##`)cuG$B,$u,#OPIce.Fl#)%'TOCml&aP#2v#OL2######$##n)c######7-#XRW###_>#}d%?+J&##qb#.@#8f2######$##V)c###$##4##%(c%##mH(:Z#~}Hf5%,u#O##?[*#########D)c######$##GCaFu#w6(W5$i^2<@+zb#~l#"},
{649.41,78.16,8.21,2.08,"######<,#{mQD>####FZ#,nQSH(###$##XqQuG%:5####0qQ.,#W.#('26nQ~-('9%3U1Y4HQmQ}P$pY#l:,rx3Q&09?%Fx*~$(:.#:oQ=Z%EJ.T[$pA,p2=smQLl%_u#hT1aB3gf1;A(QZ&Z33###/()e6*}Q'+u#GU$:}G_d)>[(Cx-(y1vl%w(;.Q%Rn("},
{324.80,164.77,7.89,1.69,"mT0eY####dP#]r/MFGb>$d-(>,#,=9jC:5S/###8##BR)Q%/B[*.,#%##rZ%ZPLK?(5##WJBNQ'F?&|e'OAW(##E$'ao)-@Wp./###-##E#$(CW###P>#vm&=XD?##<~Mc&2%v$+K-qnUdu&ku&###1##rb#4CW###$##/,#PJU###u*@;5#+-'$##mpUvu&"},
{347.58,168.65,8.23,3.24,"L5$###%##]R*K~0###8##6n*a?U+##7d'ng+uG%xZ%j~,o,L=[*%##l,$MgSZo4B,#7?>StGm@U},':v&yo)%A.i:;vG#=8Q.;:%s?@Z%@@'Ce+W#7aRTGZ%T'JtlP`?)O,#<R&0BU###(c#pP#O#L###$##0l#3)7{k#$##Av%I@,.,#&##,l#d%/###-##"},
{236.45,172.75,7.19,1.44,"9K5###%##a>#Y^f###9#$3$#Q;A###J$L2-&bG$###@dAlED?z<###&##=Q$j^f;5#/u#5$#$GK:l#^f1Fd'-m(I?#*r5biAq(?###'##/Z#p_fc#&.,#~##|GICe,r5&+##X..?R#X+L+##L(=######<,#,_fpb####7##C-Q(c${Y$<##%~.2##_lP2##"},
{475.12,190.74,7.77,1.05,"jl&0,#4Z#tn+C00`5%]>#;K)Bz8D>####(B&SA1######OS$C-(C8*t7(QR(H{0bi>6c$FR%l2[Bv)###_&&g[*L5$###bA&$##jc$7c3+%-X_<3w)dB(C:,(1[>u$:5#FK#N{8xY$###G##(##+.&OtZ^Q(5?'$##Es)6(X+x1.,#[$&m_.<'7###0,#gd#"},
{154.21,208.13,8.88,0.10,">%.&##~d*w.&oDB######7u5Qu%###(##5;Zpb####(##|9ZA96=5#E,$)v$`8Z<c%###wy'cn.2z9'##7/K2,#.D:'##6YGde/=5####)##_=Z|DB###m>#ve'u9Z###&6%a##-:Z###[P#~5$:#$###&##7R$rA3######9##%PH:5#.,#'##.)=QG####"},
{154.21,208.13,8.88,4.85,"/u####&##pG$=5#.,#(##<A/1,#]>$(m$6g7$##3R);S%=A1Gx2######B[$$EB###~5#_=ZSu%###Qp'o8Z######}>7KrB'tE:5#.,#8##T:Z###=Z$he''D;&##i.HCS.###(##Z;Z;c%NC8QG####$##]:Z###95#O##N`=&##v*C2,####'##noT6#$"},
{630.91,209.57,7.25,3.44,"D07#########2jp######5##CIV######0##Vu%/,####=,#~C<######%###jp######g##[pd######T##AQ&A5#M5$9,#d:=#########tip######C##nTb.,#$##iG#@u#9Q$Ol%[P#3p7#########rip######5##g-U.,#$##hP#7,#L>#L-'|k#"},
{619.80,243.18,8.03,3.41,"*A0######'##xrp######6$#2/]######$$#1-'95####2##>07######)##xrp######J$#}0eD>####d##v,&OG#2,#E5#x0:######)##%sp######'$#tfa######W,#~l&gP#gY#B,#F'8######,##wrp######G$#GbL.,####S##)Z$Yu#,6'2,#"},
{483.52,314.89,7.98,4.33,"aBS7x2$##A##XAS+.)R5$xR%l)0r-*E>#K>#Rm&gu&######.@S}k#},#[_5dASEQ95//<d%VDSng3k5%y#$eK/@H'###-##kd,I>#Z%$UBSbw.'d$lh.fBS>mR=u#*v$)&Ds~1######9H#9H&N,$FZ$3]*1c$#Z#:Q$>4>PG#L5#h#%lB4R5$E>####3Z$"},
{25.72,383.88,8.65,0.02,",A0#########Jtt#########klPsP$/,#1,#(##nH'Yu%:5#QA2#########2tt######'##<^c1c$F>#F5#/,#*?%#I(8l$Xe0#########Ktt######%##jT`,H%D>#$#####(I'nP$###;@,#########Ntt#########'-PYu%.,####$##}#%jc'[P#"},
{25.32,413.79,8.47,6.25,"<R+######'##b|s######u##.B`$Z$###N,#/,#CZ%$-&}k#z70######(##f|s######u##/1e)H%95#5##$##/I'(c$###B&3######(##n|s######Z##08ZWc&.,#*##%##VH%jQ(fY#ie1######)##j|s######U##fOILl%.,#)##%##yc&[>$95#"},
{539.76,417.41,7.17,2.40,"#########_R?L6(,u####F.@Gr@95####Nj.2-(######00,.,#xm'Sl%77Bi~+W2;95#='-V]Sq#'###uu#B'7######G-&###?v$pg+ebJ:$(Qv'y#$U]S.[SO$(<,#*`.Y]/ph<###G,#Al$6u#Fe$Z*E+u#b-$F:/#r=>I+:6:uc&pu$aI#=_S######"},
{454.96,492.61,7.09,1.01,"w[Lsb####:##/['?03,c$&##p##bn/dm'95#(##/,#OI'95#P4Fq5$Rd**V+|7)-LM`HP(H$fLP,ZJbc&$$$dd)37(fv*$##OG#'##JB/[C4vu&*d$fLP,+CrHPIw'N90L.B5A+.8K~%/'Q#######+R%}>&######pf',J/###Z,#gQLr82xB6pv']5K2g'"},
{454.96,492.61,7.09,4.26,"bwP1S&@'6r5$?s=qK5###+##v&&FA1######yc$K?(######xp8lu#CJ,87Cg'1inDYvP6@&GzP,kD<6'9v$-f,~i4du&&##Ld*###_m(V8,5-'w#$GzPe,KDvPg5$3S(,(Oi-*]~%9mJWc$TR'ZP#'##0,#q-&A,$n##L~/*H%$##y,%gT2###0##ua=0u#"},
{301.17,517.47,7.40,1.55,"###B,#~Q&_5%9##Df/{>$;c%r##]R,`Z%wu'%##PR&nd+D>####j,#__=###W%'<^KVRS[5$nVScOE,J,Yv%u@-)V.<%K@H%###.##)vJeY#zc'0v$nVSYaC<SSE6'yw'A8D^1:;~(nVS_m'A,#$m$Ng66#$%##TG#)o&Y%/5m'1Q%9$$Ho0.Q$QR'Vp2B6("},
{301.17,517.47,7.40,4.53,"lFCU-)M6'.8({Q%Ye,An+Ev'P&({[-###5,#`-(OG#Z5#fu$6)UkI'J^7~~(BT+dwAv%U-v&l)UtaCJ-(Av$$uGZP####L##R`?q>$<I)7r0/%+5I%l)U~OD5%U_5$l7'bTHq)B######t##/@)6#$&##Ul#|,$z,'5-#,w,:H%6#$8##ao/;R){k####C,#"},
{139.64,541.89,7.87,1.49,"+c$I>#1H$(@'~,$H]2dP#,6%V##dn/v$(yY$%##=.'(~-###$##$?#7FD}k#m7'VDSR@Si>$VDSr6P88.B[%Wn-|(0bmJBH%%##zc#6XD###gQ'V6%VDS54D&ASKv')T)TfGB19HJ'PCS<I')##Jc$f#%,u#'##A#$J@'Rm*ku%VH&_6(L7+F#$jm&Qp2qZ("},
{139.64,541.89,7.87,4.47,":<=]-)8-'|%'T-(VR)m6(I7)K%)c-*%##eG#y>$M5$,##U5$*:T_I&<U75o'B9,^JD{7T9['g;TSOC@$(~-%Y95###*##j?$$N>Fu#,7)JW2hR+*@%g;T)IOO7T%Q$*J'JqOaZQD>####J-#6%,###%##-6$5d&}Y$v,#fw/tb#GH%^#$&]0%f-xZ':5#4,#"},
{302.36,562.88,7.91,1.68,"(##P5$^u$5?'$##]v&,R)D>#@5#R.,Pu$Cl$###<R${?'Ye/ph0Cp5-7*Y5$*J-HD1Gb@>c$bnSaL8}5%~m=N?'Y^*Q25+qS(oSX6(ZI'#22}q<[['@MQ3I(7rSmFI]900p*i.(7rSX@,Qw+(m'O5$/H#Q&0|5%`$''00fQ(;H%Du$oR(c]2R>#A-'Bl#Ld*"},
{285.56,567.76,8.23,1.40,"j<8.]0oH(3?$t?(N(1_p60c#I%*hK5H,$[)7&##^m&a.'j3F|,N96&RS*8a3f{?4/)_yQ0e&dzQK|B%w)]gM5$%L/FOF8YvQ*v&,-&<v%=v'5[(JI&ha>D~.R32[NBh6&=o/8?$1oIOQ$|u'######vl&.,####&##$S+mP$lY#dP#b$)^#&]P#5,#mQ%:Q&"},
{150.05,675.97,7.51,4.50,"X83######,##VLg######=$#?6S######Z-#W]4k>$3l$f##2C:######3##cLgG>####86#QnT]-(?5#u5$U):hl&~P#N^.<h=######0##&Mg[5$mY#VZ#:GC}A-p6'iv)e)LSA2&##<K*S^:######7##VLg###7,#m@%aT6%##5$$#>EA.(~#&'##{_6"},
{676.73,58.31,8.50,3.58,"rK4#########)DSp-+3,#rP#W/*c,P|Y#^5$}e-O06eG$vY#qx4######%##w[VE>#T,#uS(WT4Q#%TT#nh9j'6:Q&@##%q,_w.#########2^V###,##Zc$&-MDR)>%#MtGM#%1j7P,#;G?h5%#########mMRE-)######QO4*^V###}k####>DP###OG#"},
{620.38,97.48,8.81,2.01,"sH)iZ#bJ*1%S&4H%l#1,#P'K9I)G.-95#1(Sv,$3l$###9%Sxe1J[#g%,@_7i$SJc%Nc$KK)y92cy4>o+u5$YN<Jc%/u#ZG#L?($##>p%=cKIW;NC9?e+&B/4g2=HK6v&{S(W+JH.%eu&8I#`d()##CM,Rn/6H%YH&M`4]:4h-)4$&3Z%`Z:DZ&2M%Zf5O@%"},
{461.62,130.57,9.79,3.16,"Nl#,W9*R*###O-()],,C4bd)ix]n,&gG#Ns(PXG)-&Qv&:e@aH%>f,s$,###V(6:U8Du#[Q&q|]P~0B5#De$}5A[y]^Z&r#$.##k-)xY#pb#CECId)dG$rQ#4z]du&###Z$#py-ey9###%##$#####~~#jJ2)]3###H5#Hr*ObLpb####_]#^%-EZ&###%.%"},
{59.26,144.89,8.79,1.94,"w~+mS*4$(cP#K]Qac&.,#Lo(VR)];5Av'L'00.''%*i5%[P#%%+)Z#.~&nJ0k~Qh6'-v'ZR&zK3.b7+*CU5#alKW[*P5$`,#pb#,##tT&<5HUo2RT/LB21g0+K2Z4>>S/d,#Nz:I]Q~P#m#####1->{@*Pw.###gCP.'44Z%)Z$@_Q6l$$##VG#-_Q######"},
{426.61,157.01,8.86,3.34,"*&+(s=cG$&##(gV2-(###]w%2v(###%##EG<q#'$##Nu#t-ERD?Tv'*@(0n)U/Z8l$E5#:*)kB6'~*?n(_=:y'5Hi<r#&7A)'V9<R*8l#(R(S4ZM]4<5#4Q#fT,i1ZAH&D5#zP#o]X:5####iT4Qm)ZP#'##X1Zp#'###%##C~'PB6######@R(nZ(###&##"},
{502.60,162.81,8.77,0.65,"uP$'Q$E,$+c#PR+I#%%##&])oC6j>%/$#QD:}eOxY$2##qu$%-'###I7'J^/A;;###z5#7C1lPG###y?#,<8#V<###A##<40rY#}[*ZH:uf/=h3Oe-Y6%<24z|E######)?8LK6######|<.(##Zd'<`T-Q%dmJq7-Y@*8z(r[T######$O-)?&######G)5"},
{502.60,162.81,8.77,1.49,"oH(W4;###k5#]L-yJ3###Op-c&P*6'###3%'HuH######@##~B3:u####Z,##}5;y7###~/+D'PJm*###[@&u$P######p$#f&26#$3##Q-&G:75v'@l$mg/u'PE)?9#$HR$6(P`$+###F$#-.,.,#r$#F$P=w-J$%)(-B$P&_3dh7P0-~B56(PI?(%##u?%"},
{369.25,167.77,8.75,0.44,"###0##kl&7u#OG#$##C#$p110u####U>#5fGlZ&Il$]^1Mk=-H&###Em#p=F46OmG$=?$:572[(0J&U281b=0$%Hj2+K3Du#dG$tG#'KH_7P8Z>T;=E,C^Q&NR&01Oa@-M>####o),eu&###2Z#04:JGI*H%8##{R*`7P###'##)],5v(######}G#ac'###"},
{369.25,167.77,8.75,4.60,"RE.LK6###e5$T3.]cQ[P#yu$Yd$._PkbNqb#(##F%%l':###C022-(###V##+,Gs84hl#`m?Ic$F13+7?`M>.##|Z%Rf1###Y`@D>####-##ofQ1S.m>$0.'pQ%2_1MeQ[Q')##8?$ZFI###wg5/d)######feQql'L>#/[&*H$'n'.dQ^6($##:d$*)@###"},
{345.78,175.47,8.98,3.32,"mP$######GZ$X83######jm$bh?$##M-&vX5>#$#Q$u-(EbD6#$###@,#4[FH{A###Ex'7GGN[TJ#%KQ$gV2D-(v82_#$J^TQ'56S-s,%.E5Z8/*{.g.Oi/1<`TRRTj?)Bc#:[&Y^T###Wl$;?%C@PhY####$Z#36@5?'###Wm$z97pb#%##E#$ny7###*##"},
{252.78,210.89,9.12,4.55,"Zw.{k####)##FnZ######x##}>Q7##J/2f5####A##3JXD>#[y8.,####S5#0oZxY$###?%#>jCLm)-%+Zd$J>#u>#ToZW,%WT5E>####B5#orZ606.,#6##SfCDsE`5%,##o[*Z$$BnZ###+/0pb####%##6=Fy,&zP$gu$RoZ@Z%5Z%}#$'o14$#AnZ%##"},
{182.07,352.57,10.04,3.02,"XG#{-)g$%A#LKH&U$)*Q#qmO4nO{k#1,#h21fj:d.*r$)|914~.$H$Xm%RnOxR-#?$P7JkXE;mO1Z$*7*(j.OT-cw*sL5Ui4*=43*?e6(6:5O.,/E05mO[A/ypOl%0)@)L.)g]'B/JyA2{A/G6$;n-bP#|G$;u#7R(Mu$[Z&ol#*6'Y>#&v'tc$:o0vb#R,$"},
{194.06,407.98,9.89,4.39,"Bh+V6OZP#N>#S#:TGFO5$2l#MX9jvME$'Im&Q^3@'/'8+/]/_W?Ke+`7.0$#$W5bt>J[+X##.m=X6O|k#8##@#9])@~5$}>%U95hP#fA.Nx+sJ1SI%~tHdl#36OVu$.-'3S%DK1Ln.C$$nkE/u#B##n`1#<?oQ(@v$oE@f$'(6O$l#]#%A]%H=DE@*zZ&ne*"},
{433.27,436.98,10.37,3.22,"P5$aZ$812W[*@_8&n*CR'Th0zx2L5$###d###I)######9l#D>#`U-Kv'l5%Kx'O@A'.+OH&g~ClR-###L##_e/######N##'c#TrTG[+$##^?JTpTdu&]$#KnT0c$8,#aV(_[,###F##Yo.vd,[g-Y=K<##SnTMv).,#OV%rmTZP#G5#X<,Mm(eu%+d$x[,"},
{222.42,450.03,9.01,1.02,"b]3VT*S82cc$88P{H)xb#B'(BBJ[aCpv*u?#A##WO24GK/,#e|F,v#Hx25n$YV8mN7lw/Pc##?:k6PW>$-##EI'[a<_5%###7[*|##to/LY:x7/FK$.6P,I%A6P~,$:Q&p-#9J..,####+#####N##l6B4D=WZ'A$#8a4vbGt]6'##S5#{T,7I*######P5#"},
{188.30,458.09,9.73,2.98,"###}l%=#$:5#_G#yw,{k####ed(FH':5#Fl#%l####t,$j'6###/F8h06$Z$M[%}IQnP#M7,GKQTR,&##zB(]z<%##wc$(MQWG#yA,j;6?K2/%'3f1<6#IJQmJQnP$=##fl:<1O#p19[&hz2F5#L|>qu$%S-e;4d1;<,#Jp3%JQM5$tY#b`,GA,f$(|8.<O9"},
{188.30,458.09,9.73,4.46,"3~-V5$0,#'m%(MQZo4$##ZQ$pW8uLQ'].37'&<9m%)1J+6L2iG#}-(qZ(.,#cy(4KQ37,(##4$<HJQ{k#G##Zt4DIQG>#V,$###E5#Z//ZP#s7.[d%EIQg>#|IQMn'pw/l6#eg3st>8T39##95####?d'h,%e#&###a=;v(;h&2mP#4g1A10B95b>#L`>/?$"},
{513.55,584.14,8.63,3.54,"?d*######I##I]`######TL#/W?Gc%*Q$z]+jz9JH'6u#Nf)A(<######8##Qyb######<)(QT5lY#6l#'qT.M;QG#gP#8?>oy8#########m{b######2Z#>wUg>$QG#.9.Dw+E,$y,%1e)q@.#########%{b######&##kyb######S$&yY$###3c#Fv&"},
{321.24,623.56,8.95,4.08,"U&4######E%#=HQ6#$###+(#jJ1)v%M5$4%%xH'yY#^P#lc'a4L###$##m^#<IQTH(l>#*_'FJQ]H(J>#Yt:/w(VH&eg+dIQ{GQ###&##Q(%Le/D>#Q-#(4<q&I#.,6##Z'1[U(pKQ^91/J-]`D######l0#'?&$##WG#?v&zl'pb#7,#U6'eP#Ev'(m$ov*"},
{343.61,640.34,8.47,4.59,"rq>#########^sfj,&###$##I>@LB5RG#=Z#RQAd08cP#~N/(r?#########:rf}k#%##/v$*}A36'2##zE=yD0G:;###g.(rg:#########8rf######_>#T8]######q~,t#&95#0,#H$&ZA2#########?rf#########vAZ95####(##|G$fY#H>#A5#"},
{163.55,650.75,9.09,4.09,"#95######~%#)6RZP####3(#ZC;[P####n~#>6'_u$?5#6$$UGO###%##0q#^QP,m(.6$6;(#8RHl%L5#_)3id({#&HK(T$Qx5R###(##l(&P..:5#OR#7-FWeHac'?##{_8[g'C6GYx-fS0@iB###$##rL$j>%###+##q7+Bc%95#'##z7,rY#XQ&??$5//"},
{229.63,659.35,9.57,4.58,"`U;#########6=s###&##Y$&mC<.,#7,#U<:O$(^#&###q@'Kh=######$##t=s######B##fvV######Em$96'######-?$PU:######&##r=s######+##eT`:5#;5#&Z#@$(95#/,#$Z#eJ2#########$>s######&##c7XJ#$ql'G5#(d$dl&8l$2['"},
{91.27,128.88,10.06,1.72,"Tu$sv%8B2gG$0I%D%P(c$$###W*C#P###uG#,h&D#P###|Y#LZ&AK#z#PE>#m]1b|1&S,C#$U%P1.,%##.U)Be*?//c5$Eg00T2A##6(PsP$Pe/T5#C;-sU6I$PU/-e#&cS't]1/,B,$(Ml#G*@4e+'E1Wu%7#$;'/4;.D2?*@)A}?z$+:v&>92@*BFl%q##"},
{226.47,396.02,10.88,4.58,"d(5lU8%-'%##xT+9JR^?)###ma4rjG.,#@c$+n$y5L3u#gw.p06lP#c~,8S+Yg6]?%&j>FZ$bRRZ>$mG${R&IV<bZ&FI'?GC#?&(##K2/xRRi@-_?$cvH8%*<SRWc&_5$.7%W;7@j@2?&G-(OS,d/+'[MdR+HK4iI$DRRcP#,RRuP$ub#[p+y6*>[*)n$pdP"},
{27.63,659.40,11.55,6.16,"OG##########YA^######l$#@B^XH(###f##H5#pd)OG####8Q&######'##gA^.,####:%#8]Zc@-H>#S?%~Z$Yo1dG$TG#<[*######%##(E^tc(###O##}s='eP:5#Z5$}5#8z~D>#;5#z5%#########Sm7Zw/######XB%20XL5$###2##'{[W>$###"},
{652.10,122.33,12.12,3.29,")R*######3##]Tb######F&#dTbS5$Ac%z$#[P#T6$4i<$##9K5######D##'Ub95####V&#8y_O]1Eu$}##jl#Ny6iK0###YC<######@##_Xb~#&###8-#su8)2<7H&7,#N&,d/.&R)p>$2C:######+$#:9^A,$M##+s-nC2q-+Y##+11&~,lZ'UG#6z,"},
{315.17,169.77,13.24,1.67,"O95###1##>?&X5Ea?'I~,kc%N>#ac$tGEJ-)######|Z$j,&J@,###J##?6&>eTDZ%x#%>}59-(h,$K[=LIS###*##6*5%m(B@,$##W5#lG$:fT###|#&d?%&<A5##HiT~@,###.##KeE2-(66'Q,$2c$bP#IfT###2,#B5#beT$##Bz+pJ2%#####HiThv+"},
{524.87,372.06,14.22,6.26,"3c#sL,+K4###a|3|~/`#&/,#c6&*.'2$(:5####J&#i{C###gH'eG8Y+L*##,LP7L4?u$P##8y3zB-KQ'(#####v^$hGP###=@,0o#EHP}92;HP.d&SQ&w(*).+k}6iw-cP####>M')4I###F,$>,#3I(P-GZl&1##7.'r39.,#W6#eRKqP$###}$#2*E###"},
{674.05,31.91,14.55,1.83,"#####################ZP#######3##[?)######*##2A06#$######0##7WB####Z#G&CYZ'###Mj,19XI#%###l$$K;XeY#######'##w:X###9Q%A&(LYF1,#-AC$a>4h:D>#2I%W$<############.;X###%##(##y;X3,#n8'le/N[)gG#ge']&/"},
{674.05,31.91,14.55,3.45,"[T5######)##D|Wo#'###[-$eb702=ZP#-@%1T/:[(eG#~A(]m+######3##zxW{u'###ti,?|>~#H.,#z$@$//~|W1,#+8&pb##########hSBCNC###oP#aS'({W.,#e,%&##CBO###$##############/##A,$######'##pb###################"},
{543.48,128.84,14.58,5.86,",##=Y3xg;###bu8E*?RG#hY#rQ@o#&VZ#{$,t,$-R)B%$S:<&##9'$46N*.,MmJ%Q$^u#`95?gSH5#sc'5?$IH&0S+88.w#'$##{5#E[%tdS=n.###~$#:fS{cS$##U#$NW,bc'E##m%)z7.######.7#QZQOG####w6#~dS&@+W##yb#X;3xY$B(#'.+zG$"},
{456.18,183.62,15.24,3.34,"gP#Tu$A,#/?&<6%<02###<Z$pe*6I+###P,#8v(95####^R%6,#rj;.,####bx-dACuG%+##%qWAJ/.,#pU&rv+:Q%pY#i&E###.R>s5&###^XD=|1KQ'Y5#PsWuS295#,x#<],WoWkY#?A(###;e(rL0mP$416t[*%Q#K?&urWo?*###n,#7(2hB7###/##"},
{525.89,371.68,14.52,6.16,"2##_r+U7.###|i1By1A,$###x@);[&IZ&######*0#l':###lH'2P4ziE1##[pQlB3pb#i##Dq7og.Fl%&#####X:$-tJ###1d)^~#SaD2j;tlQ5d'{#&)j-9I*IQ:d..A5####h(&TM@###`>$4,#&-%,vI6#$F,#Zo-,{7###DS#rmQE>####d%#*L9###"},
{138.09,406.36,15.59,3.26,"GL5@x-W,%zG#|xZzP$;@(#Z;`,%:c#?UR.L6{k####,U4[u$*bF<R'dR+L@&%zZVu$[%-Ue%OI+>?$[zZ4-&W>$###wp6<6%5p5bl$fu$R9RCyZ+c#KI*$:)%x/il#pyZf,$W,%%##%17o#%I7*T[)ZS,gr<.Z@1B2Y.,K$&u?&BN6TjErY#iY#]G#8n,2H%"},
{276.23,442.03,16.48,1.48,"tY#L?&`g0Y~02m&rZ'<-&Nv(95#)##fm$+w,&##AQ%Z?%6#$+n*p1.wST>[)nST5%*<U5X34`v)|6&`qPdh9$H$:%)PC62$'0I)r8*XSThQ'0WT:ZMT'7zu#bK*dUT-jAR,$)c#U#$2o,6m(&c#2%'UTTf08o$'M7,2q0i=HNq4f7.]6%p953[(sb#8c$2l#"},
{135.05,463.52,16.54,3.26,"2,#R5$6H%5$'A?&v5%`A-+^3Pm)fY#I,#n]/.,####f{=7H%nn,rK0Ol%Y>#H8Wd#&`$&R(M{5&~>#WKHa-NA,$###@a@2?%^_;OA*DU6DI&3M[;?&&w*MK'Z/1d6%PM[<I(W>$###cW?hZ%48/e-&vZ'tYDdK[sY#wZ'H^*Q]3*H#:L[IZ$oP$$##t^8Nc$"},
{523.45,43.05,18.97,1.47,".,##########Wn[###*c#Sc#uK8UG#WS(n.+m[,[6)M##,01_5%#########fn[}Y#qb#H##o:=rO/wQ(u5%%S(P_X0,#K,$$m(######$##sn[[l#d6(cl#]J11|'<dE>R*[g4s]-|Q&JT3]m+######*##Tn[###`u#E&'r@/$##R`(_@R}u'###0~#Jq["},
{95.71,115.70,19.66,1.48,"/$T######9##i$Tt>%;5#N5#C&TK#%###+##8&Th.*AWC###/$T$##.,#Q##X$TXA*FI+]5#/7>6&0E,$H,#@W-}HLP6)%##e%T>u$0,#A##AlIB183^1?Q$M_4$^62?%LQ$X~'Io3###$##Kp.wu'######Xe&R$T$l####3@#@$T:5####S[#-$T######"},
{294.54,339.25,18.07,4.72,"O16gB7f[%~&3c'Vm7-($&?@'k?'4n&$15fQ'x&VmY#gY#*##O1:u6&~(V$17k&Vtc&<&.IW1sv)x7'~mBOL6Q(V_G#L?'0,#)J,xT/`NAD$'F)VB)>k6(].($S'K'K2M;Nu$Ci4`?'s5&###^d%_I-2,#nY#7R&aI-U5$3l#lc%Xe.ZP#2##[6$5@*OG####"},
{286.47,395.00,17.27,4.75,"zG#IU76/+X]4tq6=jEx?$i/3t'Ugw.~u$7R&I?&6@&cf1$v&lr=p5%D%'q7L7(92$%l)UNOEi&U=Q%hn-]L-<@*av%Z37LC6UQ?NK3-Q%Xu#l~-(E5u%UsH'k(Ulx3d@+1/)@x*]dCX19KH&|Z$w$,###5,#j6%rR-6u#c>$+%(9I+7u#Dc$}#%/]1D>#,##"},
{39.90,672.59,19.57,5.95,"############i+M######U%#r&_v#'###h'#|u&#01zG%A5#gY##########q'_%Q%###8&#*)_M*F###w$##e'U`AoY#Cc$Ru$#########LsSw99###5##Hr+*'_###'##|f#P&_##################;c#@c%######I##I2?######=$#UcR######"},
{447.23,48.06,20.41,1.65,"yZ(#########a&W`l#>c%###SQ&E*+:x1~P#<U.@U3#l#qY#56'#########C&WhG#8v&uu%#n,zZ#1}5v,Nvd*~G#iu#q'Wb5%#########m%W###*##L6%t+M###t##-9Qxq?###1$$>*WfY##########;&W######.##?%W######v%*0cO###3l#1/("},
{475.26,382.11,24.32,3.51,"nJ[:5#$##~5#HC7w$(G#$w>$,P@<6'###A,#W'5;m(H>#Rv%DM[n5%F>#<,#Cr7u45'f/X,$iM[//-(Z$2Q#0sCR5$0,#B%'^J[;5#B5#][$]_>Xu#br-$b?WJ[L>#_H$#F3nsCOG####E##VJ[######X5#y}LU5$f##:'3On,t~-:##i]/Ai9r6*###:$#"}
};
siftPoint_save sps6[] = {
{90.79,16.08,4.04,4.78,"+7TZP####6$#.T~.,####|##qT~,##tc($##4H&*##E-)###.T~K>#tP$uz%{YIn['S%.L/'OV~X#%2v(=##'I*%##e$+$##cT~&m%0H&=]$5<AUS&=p7[B(fS~VG#Z6)US#Z[,###qR--##I]2/u#.,#&q*U'7xb#PG#~f%^S~###XH(>6#Uv*$##|@/(##"},
{379.90,20.77,4.33,4.45,"k'$=WC######R*+)]3%##Cu#55EjP#=c%sR#U,%R,#-d);##Mw%Q;A###$##=xB-S/###4##1FX<c$W>$c###-&0c#0Z%)##pd,{k#9##ru%$AX######g&#'AX%##pb#{1#t5&%##B,$6##_5%###U%#CW@?@X###/##x*-2dU/##6#$$'#.,#(##pb#$##"},
{405.79,22.46,4.21,5.05,"HI&ae/######msX8Q&######7rX###]>$$##RG####F,$###Tm(:5####`%-zpX######Ke(DsX$##qb#g?#wc'%##fG$(##D>#######lpXxf6######msXZnX###&Q$&|(ju&###O6'2c#&H#O7-####28%&-3l$###pB2KAR%##:l$gu#gH''##0@*$##"},
{514.91,29.62,3.74,4.79,"0pZE>#Hu$p%#rtJ&##R>Il5#SnZ###[[+=,#5?'###XJ1%##{qZQG#A,$;[%WoZE##~GNF5#SnZ'##EA1A,#1v(%##k09$##9[S?l#n-(8y/1oZ7R%4q<8###oZ}Y#L@-,##~H(mG#0B5###xg6`%)PH'j@*=N?Mv&>u$~.'EoZ0Q$6#$sQ#%-&PZ$0Z%%##"},
{189.67,60.46,3.92,1.87,"P>#Yl%dP#D>#,##|g5qu%W>$Cv$W/2fZ&UZ%w-*.,#$###-$/,#h?$Fx2.,#/.&h*XN&XQQ&h*X]&Xsl&$e&*<>.,####L,####<,#2:12-(u#&km%h*XF&XE&Xku%1e&h*XWD>###$##GH####&##tQ&bG$###&##S[$kA3[P#8,#EH$b05$.*$l#L>#dG#"},
{189.67,60.46,3.92,3.44,"#-$w-*.,#$##L,#*<>.,####GH#VD>###$##dG#/7*$l#L>#UZ%Bv$W/2qc&$e&h*X]&Xsl&h*XE&Xku%1e&b05fY#8,#EH$W>$,##{g5qu%QQ&/.&h*XN&XF&Xu#&lm%h*X`83###&##S[$D>#[G#Yl%cP#.,#/,#h?$Fx22-(###<,#2:1bG$###&##$[&"},
{189.67,60.46,3.92,4.85,"$##e5#U@,[P#CH$YB37l$C5#][$f~1###)##e?'ZP####(#####S?#7;>###8n&~|W:xWlu%~|W7xWju%bd%Bz1j>%###=,####y>#3i<.,#C6'(~&~|WPxWAxW16&0.&~|Wc82ZP####F$$I>#k?%F-(:5#7R)b,%RQ#BK4Hl$pb#2##uC8QG#&##$c#Tu$"},
{189.67,60.46,3.92,6.15,"~G#@6&95####9##(D5######&m#$p2.,####d5$^P#A5#^P#G5#'e%aw/###=%&h*XI&X7?&h*XV&Xu#&=7&:16mY####L##1,##?$aT0NQ'z#&g@&h*X@&XZ&Xmc&f@&h*XZ96}H)ZP#+-#;5#0,#06#t?*###$##PR#+C9D>#=5#rH%*h7=u#$6&R?&9c$"},
{38.19,64.48,3.67,3.22,"^:6###AH#Y[(%^R######%##sxY###/,#+##wc(###/u#(##W8,C]3&I&7c$yxYD>#$##y###yY###D>#<##.%,&##xG%%##`xCH_>$##1H$*yYeY####n&#:xY######m$#F801,#PG#7,#~/YOG#*##ar,IxY.,####=&#;xY######d,#'J/$##uP$wP#"},
{320.53,64.80,3.79,3.95,"xQ%SH'.,#####A'/L6######*)4w$,###-##OK5######Q##f[#Zf.NQ'###d8(DjWh{?HZ%DjWJ*C>l$e$%2FE######_$#IJ(k-)Dm(2u#$%+?%'AwB|oRreW&H%su$d%@_<F######;%#mt=###$##N5#Z)B###1##o1,-{?###$##7V'x..######h,#"},
{359.61,74.28,3.64,3.41,"LH&@Z$OG#F>#c1Y######'##63Y#########fo.#########Z>$0,#/,#;f.j/Y######mB-Y2Y######@Z#QXD######*##Nd*;#$>5#~_5Z:<######>/B~.Y######j1)cFJ######1-#nC77R*/,#|b#//Y######J.%HkL######^{&G[+###%##|(0"},
{104.68,78.10,3.48,0.44,"dR-###PZ#P#<D`9]>$f-'3['M:VE>####M5#M.-######'$%v6+###ZK*Z`6FT.bQ&|9VQ#%Z;V?l$hu&B##,jC###hY#,Z#OG#RI$AV8lT1ZS1t.#X8V*R$:7V(##V?(B]#MFFjY#+u#5##=5#[I><$)Ac#hv+dB$>A1{^$D6S$##D>#p8#D1:.,####`,#"},
{37.31,112.44,3.84,3.13,"G6&mP$###C,#>g5ZP####d,#xe]###D>#J$#qm,%##/H&5##BA$i1>######^cF}>&###-##Mf]0,#&l#R6#{6+6,#6m(a>#@LGMuO.,#+##Ji]ZP####C##Vf]mY#fY#g##|H)1c#},')##&i]~c&ZG#g_(Qh]######T##Tf]###95#@##*$(%##W>$%##"},
{381.22,112.40,4.08,1.23,"tp#3)@###%##_i)&aF###8Z#Wr(-S/5u#k>$Q]'4I+Uf)D,$/A.rP$E>#L##^fU_l&Q#$:X-1}/YI+up1pq5s1%{*H^c%F>#,o1###D##q-*(iUh-)U5$,[$zp3Cg,@y3|[+5m&-I(F[*;Z%U,%###q&#EZQNeU95#'$#$);e`>'6%vQ%ke.$##6@%G/0GZ&"},
{399.09,126.86,4.00,0.14,"Iv)/c#Ce')TQ+QQ)##sQ%m^+TB7<H$D>#dm#F,#5J-$##9Q%W%,<_4]7+Nw)BQQ;Q#Pv(S8$UV>]&%j&5bw$H04{8-9^.e0Pe?(Xm)i,%ju$kTQE$(a5%=Z#|ZF5@&I2;qH&k7CqKOJC5cL2b,$[j?E>####WTQYJ/$##9,#?jDZP#]-$zS0P,#<7+FH#ANB"},
{417.71,131.33,3.71,3.82,"2##(J)pb####I>#*Y=mY#S?'oH%K(:O/*G$)Uc%OG#&W)a..L?(yS$d(>3##>7Rn'.9B.~fHyg5ZS+pUO^-MtV<J>#wM6w5%3m(W5#*YABu$-;RTkEp].eT/)C1K9Jb[Q~/+p]6M,#|n0F-####E>#N-$1v(<c#lu&c6#2r@?#$`?'dP#r:7bG$sG#eu&9Z#"},
{338.39,141.23,4.00,3.28,"`1_D>#######u&]######>u#T5$/,####@[)8u#.,####A@*F1_###E##&@)z0_'##|k#Qn'Zl&Zx--m(5):(##}}JD>#Do0%0_###O@#ag431_*##?Z%i#$[~0%f*7c9[Z&SG#4W>f/(Vm*V{=rI',R'#d'B3_3Z%%##gP#%g/=^7Mx&]u%###;80,x%bg9"},
{387.53,147.37,4.25,0.01,")6$*.)Ze.%v%B.+_@,BH&b?&I;d@-&#-'cm#O=C_~'aJ148'i%#5SWX>$###;y00h^[P#&##q<do7.Ac$X?$9V6Km)M9)293iQ#q;BEc%;5#iRUn2@###4##5:d6l$$##tm#U&4###'?#}81[P####Xl$Qg9rkN###(##Y^4k9d######kI$}A46,#2u#]G#"},
{152.08,165.42,3.38,1.22,"Fw.C%#w0:z##(n-T'#IZPq@#MEEm$#^T6g0#6#$b%#V{B)##'?&ZU##ZP?##4vO:y%-i@:J#{/4=p$/5NdI####0%#9(<###.c#;h%&ZP;5#Om>;24}v,mP#.x(3}5<]4'#####v,#<~/###^G#6p$IZP]Z'Q&,w#&NQ&Fy0{?'b$)7#$4c#0##)l#)Z$.,#"},
{152.08,165.42,3.38,4.85,"O#%4,#.l#U#%'l#~I%,K.,w,p&X*R'FQ$/e+2&X$##+m()##u$+.,#4,#x#%Z08$##%o&LvI9%X###+s.$XA}%X###>?'_G#{/2######)##^&X3,#-n+_?&N=FH?%}&X.-'t&XG,$36'`P#i[+eY####&##7&X'##h6*###vD>Tu#w%X*##%'X^P#^6),##"},
{257.40,167.49,4.21,0.63,"H?#e[FK#$ep0@x%O]4r;,=~.lY#_P#3K)PG#-Z#~>$#l####gM4B=Y###3A%B=YHg7+o)Jn#oI-Oc#[/IeP#tG$YG#n?)H>#bm+f['###0[9W8Y+9/1c$O2%LI*}<0Y1;2,#oP#0-&$?&0,#H>#oZ&$##Ts2W5$k~+7#$M?#%##cd&s,&###@,#ll&J>#.,#"},
{201.73,176.08,3.96,4.77,"$Z$PG#cP#M,$D*D6#$###&##bn[.##?{<k>$HIV4##O<?F~(jY#mY#0l#~P#MGMPG#2,#XG#{n[^>$>I&V.+3f2>]2/~&pGBwb#V5$1,#=5#s*G.,####:##<q[@<3r5&XG#;S(it[XI-R#%ZP####cP#tP#mw/######r5#sy481+zEG*##?x/8G0}n[(##"},
{383.89,198.42,3.70,1.14,"sDTyH*.,#fH#jn/_l$kQ&]$'S>#au%06%)c$###,l#V>#%Q%y?Tc##z,'H(#}J0W=1dC<(d#RcB<]0wG$x6(K-)###17$~S0p3Hr%#M:<S1#1g5GD&msHU[$I@TuG$]>$y1%mI.[##7&/+[%:$)4&#'}HrQ#IaG1d#}08pe'R?T?##mP$$&#W>$j$#WI-6##"},
{383.89,198.42,3.70,4.83,"~v*(##d}D###7]/-m&C8N&l#gz]5u#5B4pY#NH'(##Oy]###Yn/:##yV<###E7,7d#H|].6&Jy](c#7'4:v%_7.###:|]qb#Vn*V@,VS.###8l$7R'@%?Gf2Z}IW#%uv'$N52I*9Z$qhPPm*x>#A&2.,####2,#;J-(Q$~>$,H$v#'@5#j,%S>#du%#6%k>%"},
{531.51,218.19,4.44,6.26,"6H$`#&RG#G5#^<Br5&/,#E,#H4u.,####,$#%(;######0##xu&K,$vb#lP#`RV9#$;5#'H#=4u######}##&V=######+##|5$Cl$9#$bP#]o^:c$8#$)H#C4u######p##0z;######)##_Z%mG$(l#`P#ZwX<5#95#*?#<4u######7$#vA4######-##"},
{74.33,224.74,3.83,1.75,"iY#XG#?,#NZ%@5#AB/Ol#Ic%zl#e&3j>$A#$ZG#$##Sl$mP$F>#vc#+03tb#tm&UaX=]X&-&UaX~]XKH&=.&'^3=,#AQ&eP####/##,r8|k#ZZ&1w%UaX8]X(]X76&xR&UaXI/2{Y#ZP#)$$*##G,$ec%D>####5,#[e'qS2[P#.,#v,#H18(c$$##>,#D?&"},
{74.33,224.74,3.83,2.99,"sb#E5#~P#26&nP#nB4iP#$?%JH#P82U#$7#$56&|k#SG#'##oG$/?#no4Fc$Bw%UaXH]X8-'q<W2]XVQ&H%'sf1|k#:5#@Q#n5%VG#z]2m5$el%bm%UaX']XF]Xlu%G.&UaX~A2&##6,#Rw%H,$95#L>#'$&OG####SH#bL:###$##E,#vM;.,#J>#^5#}u&"},
{74.33,224.74,3.83,4.49,"(##'-&rb#_P#@Z#xx0-u#E>#Mn%ue1'##3,#NQ&95#>5#f5#C,$D?#/x1Lu#/n&WsW%TXxl&IXX7TXku%b@&G);###$##K,#u5%eG#Q'4rY#}u&M%&IXX:TXySXpu%em%IXX8z:OG####~H#9?&%l#:,#]P#k#%wP$}5#''5s>$n5%XG#N03A-'tb#.,#O>#"},
{74.33,224.74,3.83,5.96,"i5#EZ%D>#5,#L##eV:###$##c?#1z:D>####H5#U6(PG####&##|7&on0&##zR&UaXE]Xku%UaX']Xz#&uv%-g32Z$V#%dP#~P#1Z#090sb#2$'}d&,OW2]XA]X|u&Z.&UaXNK5kl%tP$p,#ub#)##FZ%rb#Pu#9#$B6#]S16u#%?%]>#cB3jY#du%hG$;,#"},
{202.42,229.07,3.99,0.05,"2u#w#%6#$>5#2Z$.L4~>$>5#e-%(x0-u#&##kG#16'###0,#'##z?%Po1pb#d@&~|WexWsl&~|W@xWz#&4n&So0^,%I>#w5%1H#%U6tf0D>#)I'Be(~|WCxWNxWSQ&2%&~|Wg~1L>#tb#nQ$:^24U4Gd*)H#1q7j>%{5#$]08c$gY#A5#&'1?5#xb#'l#B5#"},
{202.42,229.07,3.99,1.62,"0,#lG#16'###w5%So0^,%I>#nQ$~S1L>#tb#B5#?5#xb#'l#&##e-%(x0-u#5n&~|W@xWz#&~|WNxWSQ&2%&&'18c$gY#A5#>5#2Z$.L4~>$sl&d@&~|WexWCxW)I'Be(~|W$]00q7j>%{5#>5#2u#w#%6#$pb#'##y?%Po1D>#2H#%U6tf0*H#Eg2(L4Gd*"},
{202.42,229.07,3.99,2.55,"/l#L,$=5#sG$`G#>K00,#WG#y##'f1Q#%;5#*##Md$@H&PG#W>#&[$f..Y>$j@&t3Xc/X;?&t3Xk/X0m'J@&QA.Yw$8M7n5%###D,#+p1C,$MH&2w%t3XQ/Xe/Xmu%_w&t3X?T4###0J&1%)qb#0l#(Z$Y5$rb#S>#q?$Y82Z>$(##g#$Q'4<H&###A5#Pl$"},
{202.42,229.07,3.99,4.57,"M$(Fu$4n$^_:t08EK11,#fe'V6%x]4OG#(##+v%sb#>5#I>#6#$m,#y82Sw-]I'PsWCoWTd(PsWToW46&X7&+:6gY#7,#o,%I>#M>#hf1iG$ql&I.&PsW<oW(oW}#&P7&PsWoS1eG$2,#U-%YG#lG$A5#2,#L>#bP#Pm%>J0c#&=5#T,$'K/_5%$##bP#?Z$"},
{316.03,232.85,3.84,1.54,"jY#=c#;l$/u#@##Gz37Z%95#xH#y&5)c#rb#a#%M5$?5#gP####@?#E197#$X%'w`U2]XSc%UaXDxWoZ'nm%p95M5$M5#Gc$###W5#]q8PG#8-'G[%UaX5]X&]XGZ%1S'UaXuo5[5$*l#Q6$&##1l#3m&Y>$$##H5#r6%<06sb#J>#R5#up6PG#uY#Hc%d>$"},
{316.03,232.85,3.84,3.17,"YG#d#%yY$3,#f#%Kg5B,$M5#K6$=06f>$3u#pG$E>#]G#3Q%.u#7d#w&5wY#gm%IXXpAXP?'IXXzSXGZ%w@'917]P#WG#eG#95#B##Jz4Dc%]l%q7'R3V5TX)TX8-'G[%IXXI96###:,#9R%F,$^P#?c#Q#%,u#$##E?#bU8[P####G,#Kh7C,$%##(c#Cv&"},
{316.03,232.85,3.84,4.48,"Z,%YG#iY#0Q$sG#&:5.u#G>#r6$&B40,#yY#7-&zY$5,#rY#}Y##I$?f3+l#3e&msXCoXVc%msXLoXGH&kd%f19|k####Cc#F5#7c#gB3M5$#v&1w%msXbJVGoXbl%a@&2ES(V;pb####,6#vY#Bu#a5$gY#]5$95#x-#q&4O#%E>#F,#hC2w>%.u#hY#Eu#"},
{29.46,243.47,3.39,3.22,"pb#$##i5$1c$c9c###$##Mc#:so######8##me1######0##,u#R,$$Z$o5$r1g###$##,H$(to######.##9^8######$##Rc&PG#&##4$%H8]95#(##9m%Uso######F##o^:######$##LZ&Iu$95#9##HQP|Y$###)##;so######U##`T6######;##"},
{248.63,251.08,3.33,6.22,"{P$T['26'0,#1v%MATB,$)##b~'/2=iY####_H&6Z%YG#/,#%##TI%bmTZP#K9+TrT^mT}Y#TrT>_<]>$/H#nK3;c%3,#eP#bP#;l#onTZP#L[*,Z#TrT?IJlmTH>#_,$U56aT695#$##9@&4H$&~+@d&/u#O,$g>$hA%zi@Wm(fv*A5#y:/0I([-*/,#9Q#"},
{109.85,281.49,4.04,3.54,"eU8f?'{e-J@&e&0gm)}G?4Q$E]XA,$c#$)g%be0###%-$rd&-u<+-&m-(j5$^(;]%&$eN.d%/_XJd(>%+D('@S,9R(c//d?&rK+eQ'uJ%,vI5w*2J(7&0d`XwJU<H&Z>$0/ArI+ol&9#$FH#.,####/d#F6H:5#F>#<5#kZGP,$rb#$##)C0#l####$##:c#"},
{130.50,284.15,3.85,5.73,"Sw#R9H=7,yQ(Q)'#'6&##=Z%O&$5f3###0,#rl%Uc&###/##u@-eS*Zc%~4;B*+xT0Er<:T1:_NqK7(c$r##zM;{k####9##AI+,0+|R-LJ(gV>WC&5tJ[5#pfVY>#Fl%p%#~2=:5#)##E6$9u#dr5~I'nv+1U74K($I(Ox,DfV2,#95#u.$/?P$##0##(9'"},
{529.91,285.99,3.76,0.06,"]u$Uu%###&##kmS+u####/##Dsq######Y##I97######2##W5$;l$0,#Z>#DLa[P####,##>tq######)##@(<######%##T,$Au$$##7,#^]`fY####>##@tq######(##^D@#########Bc%G>#%##9,#Y>N######2##-tq#########l(>#########"},
{233.81,287.29,3.64,3.54,"aR+WQ&'q,'B1Jx/=v(l+7$H$m@TeY#p5$(2){m,###57$Hz5]V6oc&gA.$$&&_;Zu#KBTf?&^BTlc&D7*l_)..)sH'w]-(8-Sq+#f0k$')'.HI*M%*-3@ODTGSQ[u%7c$_JDcR)8-'xP$M?$xP$95#xl#IwDP#%h5%iP#/v>ub#.u#}Y#o03{k#$##bP#m#$"},
{69.89,293.54,3.74,3.53,"av$>-':5#%##mZ#dM:95#$##)^)PT5###<5#jd*.,#^P#|b#1/FU?%VI+-##6~)10HydLRu$0WTdW@RZ%,v$Ki>?5#95#E5#,ST###l-&ah);%-:c#|6=<kE?RToY#_.%L'J~'6F>#E>#$[$sNCL$'N*7zD1IJ/>5#1I%kn/yP$Wl%i#$xg4%Q$KH'fY#bZ$"},
{69.89,293.54,3.74,5.97,"C5#p5$Wc&iY#H,#NL4PZ&SG#`[%FK4Ku#1c$Tu$?6$Tu%8u#$##1S#Ze0~P#+7'3S@-|Cf>$^NV<5E9l$k#$rp5|d''$(k5#$##k?#bi>.,#8$(+f$3iRlGJxIVnm&+@(@R>)]/EoD~.//##yY#*l#=n+######(##L/)>T4W>$--#UD>|81.,#'E(?JV$##"},
{174.34,296.57,4.06,3.46,"0%-^@+$##aX7JZ#VL:M7#OdO86#^-*)2(#.,nl#'v%&K.M#%Uy8ia9X966.$|U,GdO*H$UH&ywHd$++`*lc%F,$&##7X0|6+D^5BH%#c9%B3gI,LQ&z6%OgO/ZN{K*d.-0r,^>$i0.4{7F>#4C('95OQ#1@+7J$w:9j5$CW=F>#;{*ro5Nm'(Z$sm(Dc%h,$"},
{271.05,319.33,4.23,0.69,".^$QSYJ>#.,#03(&U8E>####pH%I,$>#$gG$K,$$Z$.,#XG#](-G,N3k?Q,$fXYA%.}k#o5#=FCC#$PG#3##|G%xY#,u#D,#uR-Q['L'TZD0CSY:5#8u#oW)<`AKl$P#%E?#>?&yl&mP$$##Kw.0['36&dH;KM?k$)w,&~W(0I'E<@4Z%O,#~##%y5.,####"},
{453.84,328.85,3.37,4.13,"]P#D5#ZS%/p6II,/,#2$#*qW_5%###V##PsW######p6%WZ>Md&OU-:93:#$[pW<[):5#5n#IL:;,#`>$L^+###,##@Q$uc'^w,|['PsW5o-noWdu%c#%]h((*CcG#7#$bQ$###*##`P#N$)&Z$7Z$cl<h}GyS3###c5#1e>iZ(######Bn%######8##g?("},
{359.00,332.46,3.57,0.78,"###k##Y7.###1,#Z2+J..###qm(:t?CQ&k#$*l#eI%D,MS,$###C##Ii=j>%;%-OD&DwVJ5#9xV[8+@T/n`+s-*pG#[zVUe*###f$#+jD)c$u?*|##8zV9I*$|Vxc'tS,}~-]K+f/.gxVWG####2R#J_>######G##Fz6?6(R$&j,&|v#sv+6H#tw,]L595#"},
{359.00,332.46,3.57,3.91,"S(7.,#8H#(&-.%$GR+x?&&?&U:6K?(###F##)D>######5R#u/WcP#l]+/]-9o,kJ-03Wyc'I1W/@*/R*%$#(jDxY$###q$#a1Wnw*q-*tG#1K/'j+H/W'T+V.WJ5#$n,WD&Rr=j>%###D##DPKH#$6u#4e%,?&_u#~[(#P@?%.###1,#Y2+C%.######l##"},
{232.26,349.65,3.71,3.52,"D//26&zJ+6@*%]0{?*l24^,$]%WX>$@c#}^(A96$##yv$Zg1AX;au%Lm'=Q%C::Ec$R[J_@*H'W1?&>R(yp)W//{c'9A(*J+XC-TQ'M$$>t8C7*'['Tx0K*W;YDiu%W[&k(WtA/36'hY#wZ%07+Ku$T?#3v>i@+N#G/Z$y<3H##cPG^-'yR-`l#Wy8$##D#$"},
{232.26,349.65,3.71,4.79,"1L+i96N5#Hm&1M/IS/zP$xK,mH&T&/xu&aH'Ju#9m&Qu%###y9%YXDUu%`w->*+L@VIu#?n,A9Q{RS{>%Lw%TI+;$$Bx1tQ%Dl$l,%-6%1lGCT0:e,`[%4cJ?CVQ6'3@(D'+QA/)@$k@V2m%$d$h+JY5$Su$<I(YuF?l$cZ$FEV{I*-c$Qv#9'.JR([q:B,#"},
{473.84,357.73,4.16,6.16,"gv(B^/;J+gn./tCs]1bR(ee@a;<fL1ff28G1>-$nd)%v'.,#Cz;_G#4c$'A)D}<VJUFR+W&*4mCx9Ypl'H$$SJ-06%L5$?u#O..=,#U..]v$&E@e/)P7YjP#W8Y&v&2Z%?J&2[)###.##Ux0E-);##B96~Z#I|EM$#M7YY5#d7Y+##I#%TR#H[+###.##RS,"},
{404.23,411.85,4.06,1.73,"Nz-KC;######;9HR[)bG$###>$%|AHH[+mY#Jq({,OE>#~#%H:8D>####,##mJTc>#Rd)xY#Z%.q@%[-DZ~.4(6CH&^u#'{2Rx3###$##U%(hHT###Z>#nK&tK8###Wf$n[Oo>%###~##>MTW>$###B##qwBun1.,####.G3s[-######}S'DZ&######~7)"},
{106.69,421.80,3.72,2.30,"+c$o9,ZA2`V33B&S6N(?%:8/.D/`Q(eC&/d)~#%:,#GE/8v(Q[+9R$?9NGu#~U/'7)4C,ud*K9NOH#We.$.%G6'%T#%6N1,#&~(9p7Oi-bH(2v%su%c^&k6N6o2j0#b'9%O63Z%%0#Uz=Ql#88%NiB?5#}k#Rg&<z1RZ&u?*%##1(%;B5(Q%###g,#3m(L?&"},
{106.69,421.80,3.72,3.59,"T5$+$'E>#DZ:ZG#R%.?f#1$J,Z#JZ&iU$A7-Zv'PG#36#rP$i6*U}6A7-90*Qe&v7PSc#@-'@F8)x1)M&wu&Wu%###GV'Av(V7-|n&)wJ5.,m@,cJ.ZQ%qG>_6PHv&Bd%+{(=c%>~(3Y6S,$of._./$]%F96-$$y$+KH$h8Pju&{g(JB6jB-kY#xV2}6+&##"},
{121.70,421.09,4.04,5.75,"nL<S00_u$uI&EJ$KPIic#M1:&A$4$(|E=sK64u#6c#]r<&Q%fK74'+AcDcu#1)+WN@M~(Fv'%1OO[(4z6F6$#6&'7$QuD.u#3D1>A1*E,%v'zQ'&8-Sy$R.O%;=Q)+&?&E)0-7*o2,Ed*4##L]#e:=Nl#Ll%j-#n.OzY#]@-D>#Zj1b>$,A,K.-%m$b5%rH%"},
{170.48,423.15,4.18,3.58,"Il%K6&~P#U|6N5#[..aR#)JO}>#}Q)J_'u$,Q$'|5%-A,O#%d:=[V/o[-0p(Oo(eJObc$J.+8c?_./f:&bm).c$###3i*M%.lEC`v$OlB.%,)~,|[)9e(|SIGHO<&(lS.>r*}Y$&U,y`8:#$C23<x28m#-.,|7$cq=lu#;#El5$7).d{B.A)d>$QS*|$,jP#"},
{247.83,425.65,3.95,5.91,")&1I^/PH$nB4#7#ZuMj?%&V83v#tv+ku;lZ(2?&H>#o.(K$)_~01W/`ZKpG$.z*)vN2I';d(`TLGx-39-}c%/Z$kx&u^3dc'60,q[-(P5s@/=v'-.)c&%'wNjg9=h&;~/ZV/B~-BD)G[+V##2'$~M@Jl#PQ'g?#cSLs5$NS095#Ti,w#';&-?6(b##+u#KS%"},
{527.99,434.17,3.86,6.20,"_G#@-'fP#PG#9=A,6'###0,#}5r######B##%V<######+##%##:c$WZ$6l$R[VG>#8c$fG#W3r######o$#32@######@##,##AR)|Y$###df`Vc&pP$R##S3r######)%#8g8######9##$##(l#`>$OG#:K_$l#ub#]5#_3r######=$#d%0######/##"},
{220.07,439.28,3.54,5.80,"K#%;=.#l#VA)rC;8-%2u#7n$R7W4Z#W>$5S#Rx2;-?&[&PS,###>l#.Q%x)772@'##wS1^d&=7W|,#&M;-%#1d(6}1b9W6c$'##Ul%3Z%1Z$9C8K,$)T1bG#c<WJw.(o-A,#?A&o-+D.>|&6R>#p>%OG####Eu#Dl$3w+###fg%%i?,.(5l$TV$AjFBQ#V6)"},
{417.87,441.32,3.43,1.09,"<q'+mR###$##]+0)R*###+##y-(eY####6##T#$A,$###O>#OW,(SSVZ$/e,it[X-*-##h,$?jA.,#VG#UH%nP#$6%Z,%vb#Kw-;#$h@&F-FPo[###lP#0=+7tJK##2?&yR$6,#^v%I~/fY#D>####2$#|q[(n-&##ZP#S675?'}?#SS1nQ#?u#$&)/v(&##"},
{114.41,445.10,3.28,4.62,"95#$##bP#Sl%5$(nG#Vu%EZ$;81@,#3c$wI?ZP#PR&i,&Y<Z5l$8##4-(%c#SW@(.&v~2Q5#M&[VH$7m)|9'Z-*UE(TcR:~)eY#;##Hd*+##>|DW5#w96dH%%([p>#(B2YQ#Em(G7#R)[_5%95#$##gu%E>#(w)^P#Fn*mY#b7=h@.XH&4,#2_$J'8sV,#w,"},
{528.39,446.58,3.75,6.20,"$##g,%r5$7#$)<CSG#@l$J5#Zss######T$#()@######B##,##Jm(pP$###KdU16'C,$K##Yss######1%#xz?######A##%##(l#xP$ZP#p0c'l#S5$]5#oss######W$#Y'9######7##1##&Z$6#$###6%WYG#'?&7##xss######f##PA2######+##"},
{394.15,465.04,3.76,1.72,"[}DF##.3?`5$2m%^T*%dM<5#=iVcv),u#A##JK3VG#uP#0v&rdVX>#mE?}$#y?(TS%dhVaQ'DgV6m%k$*oR'sjDNl%XG#Bu#sI._G#}A/aJ)rP$2H#(CJtdVKcPZG#)%%;iV%V<E>#%##H$%l#&1Z#Ul%r5%S#$U,$06$_./fY#V>#Nv%Z/2EZ&H>#fP#5-%"},
{394.15,465.04,3.76,3.47,"=$((%+]>$9c#rY#=29|5&]P#^$&JL9###'##36%?c%%##;u#~G#g;TO-)PG#.w'g;Ttl%.d(C;T>tHUG#f?%Ro3@5#mP#mc%<5#v/+e'1E;;mZ'b%+?R%g;T+7T/c$oP#RhPFA1iG#_5$JH%t>%Rv'f5#R#FX5#O7R=5#y5@`%&9aE`P#fo0v6&[?(D5#c,%"},
{394.15,465.04,3.76,4.86,"fP#5-%EZ&H>#Nv%Z/2fY#V>#06$_./S#$U,$Ul%r5%l#&2Z#%##H$%%V<E>#*%%;iVKcPZG#(CJtdVrP$2H#}A/`J)h@._G#XG#Cu#}sDNl%k$*oR'BgV7m%chVaQ'y?(US%a<?$%#qdVX>#uP#0v&UT3VG#,u#A##<iVcv);vM<5#'d%_T*})?`5$-GEF##"},
{194.47,472.10,3.39,4.85,"PG#a##IoVK>#JnV.##yS05w#=oV###Zl%|?#jA20,##-D-H%95####3sV8Z%tkF###B_0~C20rV###I>#5S%Y'6.,#@E3sY#######=b2PA2nu&###~f%=pV?GK###$##<-9,U2OG#G?$R#$$##.,#z>#;w-###O>#%-%Iy6:6&>#$VG#hH'E%(D,$###-##"},
{161.68,476.13,3.85,1.78,"E@)f5#O]X###nS/wQ#R]X(##/^XYG#dG$N##RT2[P#.,#H##zg9*##.JPqG#ZA/m##;aXQ6(r]XL,#Jd)#e'DL:######eQ$Pn/2,#r(9LH$5Z%3##}yJe7X$IU(##|$%~`XLx31,#:5#^c#X,%k,%NQ&j5$###},#4w),d)pP$_G#(H$q@-Vu%QG#%##;l#"},
{161.68,476.13,3.85,3.68,"L,#o6){k####U,#_]3$##PG#]Z$B/1/,####nY#Zu%UG#.,#,##)gPpb####vm'[2TH,$]l%M2T8+F%##lZ$>%-;5#~>#p-*:5#_?<a#%:927[)v(3El$j1Py-T2Q%&##y%C%J//,#/##sI*L5$1##nP#KzOW>$k?&'##}&H:u#1M:###>D23c##C8cP#FQ%"},
{161.68,476.13,3.85,4.68,"M>#Al#@u$=5#n>#`n,@u$<5#a-&nl'/,#0v#~u%QG#/##4[%95#bv$9~/2,#Y%%@2O]HT,##@'FLkI`5%b##rp4v,$h,&P,####i7%Wz=###$8.NU&'xZ@H#-{Zgc&b$*sJ#J4H87$BB6f#####)$#h4G95#N[+`$#UxZ,$$dxZ'##fG$3T#~;@_>#(OD86$"},
{522.26,503.67,3.42,1.67,"1##MRS######I9)]QS###$##Jd&Dy/{k#$#####&eO######+w*hRS;5#sP#3US/r?I#%|$#kL9FW*_>PR#####5VPZT6###eA1I,$SG#nP#GVSjw-bG$f##c%*{a2^QS+#####H?#~RS###3$($##YG#`G#2N6*m(###,##6-$<y0vQ)######o>#:{=###"},
{403.43,542.84,3.97,3.07,"H$#M8X%##6l$pn#0y6qP$gG$q5%aP#|Y$dG#(H%.,####zY#3W+/jZ3@+>,#@kZip:###U##$;:D>#'##4Z$Hc$ZP#4,#tb#`~/3H$`2/(WXQfZ.,####Xt/1V=%##_P#g$$Q5$2,#%c#kY#D>####).#xgZCZ&3,#:5#L68+u#M>##Z$Sl#;5#1,#xb#PG#"},
{447.84,578.51,3.78,4.24,"###$##P)(X|EL5$###,/&ajXFl%###@Z#ajX9@,$##[G#yx)9v&^m)I-%wS/_HOxY$<l#;k5jH)###.22$p0Yl&%##N$'S6&bB1:bAD>#T##^gX9$(x5&*x#Tw.2@$m1:U#$QG#tc%gH'G,$_[Vh%*D>#,&#qdXS5#7/16M$mP$3@$@M>)c#$##=R&+u####"},
{523.61,592.71,3.85,0.11,"2~$U(9######GZ<lw0######&QkOG#######OA1#########hm&w?)wP#Z$*lnGX97O#$.c#qQk{k####P##Gq<######%#####:5#3@#5C:`#P0c$&%$P:0)Nk######0K#9q<######3########eH##J/=)?###jv$v7-@Pk###;5#yG#Ny6######$##"},
{282.34,597.71,4.53,6.07,"2##j^Ou5&###B_&xQPhY#:5#hV+E7-######Uu$F>#95#TG#)##]`-P)8H+HN[BD(9aA%KRPrUP>u$+##Dd%oJ095####aP#_#%?5#AK*[7Pa?)$Z#h_(&%PFGM.,#@##_i,PJ0fY#'##%Q#hd'JH$[y.Lx14/0{.&B91go,8&2[P#%##%1+uG%###&##VR)"},
{269.74,601.08,3.98,1.56,"###F</fH)###fR'sSP.,#,##tB+q09SG#(Z#$~%~L:r>$Z,%|P$3`*XQPTG#~RPM27|H)8d#nQPrP$[G#d]'>R)XH&~R%4D9B,$>##rUP0[)MRP7,#%p'Rq6X0OvP$###J@%1A,C.'bu%^8*######YC&]./Gc%###30#vQPA7+###.##1H:ZH'[P####m]("},
{99.18,614.30,3.57,1.58,"@pUp5%hY#lw&te.XK'#AR8g/|mUQ##'L5}]$'%T]I$7C:B##Q|<Yl%<v%T7M{K43](`V>z$',oU^5$y#'}S$uN?4f(aB6;c$.,#9##)k=sh>A6(d,$Yz7,^4_9Uf?)G>#qm'HX1b/3-Z$g5%###K0%iy8.,#_c&ZU1--'Ql%Ly35w)7#$hP#J>DFd({k#4##"},
{127.51,614.52,3.85,1.02,"f|5/o2.##'-$N*,5?'(##95#hz,A,$###%##{$*95####o##(:995#mZ#VV0fj.dI-U:/Tm*.*V<c%L,$S-#:&2e,$95#=^$42>3S)mc%Lh1<f2g/&H$G*1.:%VFc$hc&az#B/2>((tc(*e#Fv(`n(Vl%]$'<n.&.%y5&^s-N2AwJ*(c$Kh#y#&=(,uG%,##"},
{280.58,623.67,3.78,1.34,"oe$&95######MB+Uv*######,K/L5$7,#m>$_5%###+Q#,R(W{6Gv)sY#u%&#:(7NArd&i6)rrVZ-*0##^v%GA1###O##Ha=l;?6u#X>#Zy)-n+Cd&nT&Zi;DnV]P#Z,#@s+f09###&##,UI:M8HZ%;5#Ex*B&1+~({P$=W3hmVQG#$##|z${mVD>#'##zb4"},
{422.64,626.31,3.75,1.79,";##|YE:$(###p1&:RQ:d'###.a.)n*nl'###mG#Ti8=&2###i00Ux0Q%*^5$<<7Pw+tvAu,%4VQv>%27+n5#s-(9?%1PDc>$C,;x(=6l$5##Ig5sB/`2@VK%TQQ%##Zu%Gz$3$('##WK-,@(hp/Q_>I5#@h1($O5R*/,#h:%Ve0###$##<q$eY#$##-%(8n+"},
{10.97,639.68,4.56,3.23,"K6GWc&###,R$z.-}e&3jF3##[_ck##nw0'$#bG$######%##V[+FK.$##e%(3q:.t0en0Q##S^c|5$bQ(&&#ml'######*#####|r,Eu$pb#L-OdE.kG$ic&t_c=#$###~>#wm,############9[$nQ&/i?5B4d#$Ac$tUWP^c######6h**n-######)##"},
{291.45,654.10,3.71,2.28,"i[,######7[$Hl%###E##?'0=6(###V##7A,j>%####H#P8.AV6######L5#w}G###A##{92pl'###d##YKU95####<##.KU$8-######w[&SJU######O&%NlP######(zLOG####E$#NKUNK)E'7###JZ#^LU0Z%###/8#?JU95####K<+mB54-&>1+1MU"},
{245.79,658.71,3.78,5.76,"ae-###3##DZ%Cx1###S##=395w,###&o#EzS######ZC$4#O`7-######wc$0xS######Hm#=xS###U-$X1/e6*###{]$:xSSA/A,$###u5#8mMpb####E.#*yS%Q%###:e#[#B9h8xI*CC1/@%;m)######z5A%Q%###$##z7BJx3######[S$=wSOG####"},
{294.37,665.42,3.69,3.68,"D-(###7,#Dc$$x0###[5#h&))%X######/g$FuP######{D*Aw'z,'###4,#q{6|Z)###n,#{'X95####2$#;'XD>####+&#-&+J?(###.##-E>0Z%###m7#^'Xv':###.%#fLLuZT###,##ZJ+Hl%######QL0_T5###3##Jn'9&X###)##Yz&7$R######"},
{55.28,669.11,4.29,0.34,"jl#t{=eY#,%'~K'sJ3&##1f.p>$5,#|Y#&///,#'H%{P#--'lS)KzOBRLF%+J;S|M@~G#(]'zp67r0]>$mY#5,#37J[P#:5#Ll%W>#c>7|7SftK)~'2w*PTHw6*{#:{v),R)AI'I<<bu#3T3###-##gK/Q$*###}-$~W?'@*l6$LK1GA-l1<CJ%8.-U?#;|D"},
{489.56,671.69,3.99,5.22,"|6(g7*:T-6m(I.&mL/*D8J#%C;.}Q)###'##9m'OG#$##fG#C6&_g(}aHQG#:'*^7?bnVh5$3sVuK4du%>v#2z9$##%##:Q#Je+]l#MoVY,%[v)UQ#j:NmoVHnV7u#?[%<8BP96)m'###~##oA,av'Cw-vb#SR(3?&;v&a95bc%XT4-##`x-5c$f19###V$&"},
{249.44,675.01,3.73,5.37,"i6)95####U,#E.-######&S'dWE######_y'qHU###TQ%0r+37&VZ&.,#,##Tb>OG####LZ#^IU######r1%<JUOG#)##Q`)L-%`d);v'###$LU'Q%95#U##TMUbG$###1$#fTH3C:###J##d,$Y5$)h4[P#sF4(Q%5?&G5#td>8Q&###%##'p'>NC######"},
{81.63,678.03,3.81,5.04,"+o0:I%$+Gm@%?LL(7'B%.P5#$JA+c$95#`>$],$BZ$Jv):5#&m(?n#t-TJ^-5mMg>#q&2<M2`0T.,#G>#f>#F7,}P#r6)oY#EH'ze(p'7R4751;y,#v}G;U1N-T&##0u#xo#;@,###bI'^$'o5%Z^)5K4nB-vL=r@#&2?$^#8-T###95#d1$uG%###{H&.d&"},
{93.02,685.12,3.94,5.44,"_J..Z$;$Ig$&97=+u#yP#p6&J&K###L>#d>$'?&%##IR%h?)V%/Gv#wtDuB*d:91,#IO<>B,&yV$##g#&]5#V/3###n,#(@(zC<?7%nn/sV)wy:v6#,xVo$%dvV,##{Q)#%#ty:###mP#)m$e9+euGOc&,R$L:;N3//B5B/#]vV+##{k#<(#%J/C5#TG#h,#"},
{100.48,693.16,3.72,5.65,"^[,-##e'/MK*5i7/,#2'-h?'(zU###SG#7,#sd,###8##N?&KOHjc#UK2.T&LL96##byUwu$NxU&##X?(6##h2A###8,##Z#M{8wq8>I+=/%IU9I*,[vU76#YvUR5#Yl&-&#5;?>5#95#~##KI$*/TZQ'D>#c6I]i8bG$e%#GwU;5####?'#+.,F>####:Z#"},
{88.33,695.45,4.02,2.24,"###i&#i1>###z,'u'#2q;C',k>%p,#K^*?eOh-($l#R7'913)c$Y%#vbO'##CcOTv#11;7o$%x03W+216&0,<m&].)]279@(M#%q,#[cO$##`5G;J)Op77,#,`<cL*i/4a%%FR+^o#`cO[&)=5#+l#y,H###/-%B%(6s9/,#GeO@J(d@-U>#gD>KI#|bOE7&"},
{88.33,695.45,4.02,5.97,"j@)BS/-1)xu&@z.*A0(?#9@(+3-Yl&sG#fG$~q1######0l#ij@9,#(|-Kv(~WBR>#/z+2g.b{4E>#>3.G-(UUO###Y5$$##/TOnI*yB7O>#7RORm'K18=.%9U6kd%TSO`#$#SOjG#VI,2##2(3fv*-Q$K:4j]+9E>f?(_H&MJ-xSH%f23##;SOVR*;c%^##"},
{74.28,707.65,3.73,5.69,"we%mj@A'Kn@.z$Q`Q$_QH:u#)'5/,#h>:9v%5%+.,#Bd%n%)fS%HTQ/S/'##XRQkV6OJ1h##]sFZ@$}|C-$$Zy8W>#LK.0J'n6)Op26.+/~NHSQQ@,oP#]39C{6[S+~%.1n&}'87w(nA28~$Nu$=c%k5#zQJg5$K@+-7%W)?`.'w&0-.)e#&II&%dK;c%.##"},
{487.06,726.48,4.02,4.71,"-S.h8&G820?#5y3Fd)r[(**:;~%L?(Ts1]GM&|.=c%;q1zY$S&3l6%Rh7S.)yM:JRFOG#jl#*KM.i8QO=}Z&~FFvR$VJM(-$:5#eH$v}7$6&3^7Jp*$g*F&,NHM'`+5:8f@%3=H$C%cT5./%###(@(7B*j>%KQ'(R)5p$kg7b$+MQ%~90H~+|[-XG#g6&x.*"},
{5.76,750.30,3.49,1.61,"######A,#+o/pA3ZG#VH%?J+<{k/,#A5#c1$(?R######u##%##/.(QZ$~o4q2=-M/X6&517B#l9?&%##P~&@A[######1##$##Ml$bl#$K^av*XZ$S[%3jiN{kG>#%##-E0A_=######(##################################################"},
{513.16,767.42,3.81,0.88,"c<0###'$#_R,l?6###)##fG$LX)#########V>##########=H&@5#MZ$~n-NA>Nu#(A-lu%.6_;5#95#m##Y7+######$##&v'Tc#hu&_,$5M=@6#;J0&X'30_$##95#Y)#j,&######%##Yl&`G#eG$O_&%Q%(##}k#Hb->u$######]-#############"},
{359.83,34.36,4.52,3.52,"NGID>####5?#2vS######uK%^Q(###%##H:.mP$######z$&oyT######&##{U[###$##B6$wx5###Ul#+~+{k####:,#XR)2lI######4.%KT[######Ip$7U9###pG#+eJ###SG#*##S~N#.,######KrP7S[ZP####Uk/%w,Y-(%##z;5###vQ'$##W%-"},
{74.52,103.40,4.58,0.86,"|99v7*rb#xJ#K?(Z1(&_;CZ#$EBN,#_?)#.#=D=######6##gZ%h^OWm*P,#P@Jt}=W>$D$#Xq=3,#.,#zx%896######@n#ql%w/+/^L2^6wAP6I)2[&2Z9lsG######kg&<.-######,~#######gC%|>Pll'###{,#ZCP%Q%###3,#Dq.hG$0,#F>#|5$"},
{234.44,108.97,4.07,1.58,"fY####A5#Vm)7D?######[7'pTc###<5#ZK%5J0$##o?NgI).c$###$##|u%/cM###3,#JQ$~Uc###{k#v##qB8E##wTc6##.,####`G#;m'nFK###{Y#eZ$+Uc###+u#q##3L99##rTc)##|k####~>#Rl%MEE###Q>#m#%pTc###C,$r##,o25##{fb8##"},
{358.76,112.34,5.26,3.03,"pP$,H#hB4a%/:v(Q,$qd#'M:SwZeP#6[#9;,W+L[Q%C%+7n%$e(PE85J-0/.K92[0QmZ&gu%/zZ,f/}b#)m#olG9/1OJ+i#$~6$LzZ)$(1,#U-R'HDfY#b##;yZ2H&###N$#~(0b/4/,####E>#0H%GH&?m)%SZrb####@~&_wZ######e$#A//OG####%##"},
{290.62,138.21,4.64,3.00,"g>$+u#1##UL:=5#~P#/?#$U8.,#2,#*@'nl'qx5%##?#$QH#DM&kuR<,#2'6ll?~[,58&([RRc&3##`TI}x5{<DK>#Yd*T,#G;4N/3$6#{wRXxRM['J81p(*}d-SS#]vRN#$.C6R>#9w-1B'/~%&V=A##P_;K{5H$L+6'(J'$-&/F15U9'##gQI4?&mP$_-$"},
{290.62,138.21,4.64,6.15,"xY$<v#PHH@H&w0:&##iu%is1z,'B~'Di6A$L@##VD=rR%Rz=Z7.m&'Yg6Q>#u-SM#$Z7.FS#4&1`(*x/Shv'y5#1/SC24So4ov*P,#<ODJ>#bTIB96Rc&2##08&_-S^u>rm,;,#<y6x1&&-S>#$OH#aK7%##)@'nl'.,#'##6H#tK82,#~P#/##dg9f>$pb#"},
{122.84,158.77,4.47,4.83,"J>#7$&[P#UG#@6&tl'###Z>#AQ%s5&%##J].'c#)c$g,#3B3Cl#w#%.Q%yb#~'3Gh6Zd)4Q$4VcXZ'Jl#w@>SH((##&JBZuIm#&^G#UZ%a,%62<4w&F11u-)/Xc8c$K6'M.&;M=%##`Wc)6%LZ&>5#;u#hH''&096$p9,}OEVVc{b#]$&`x.B3D&##rUc###"},
{369.85,198.62,4.32,1.45,"8s[o[-eY#d##[2?hS(Kp8_>#G3;~7+1H&vl%Wm*$##tv'CR*an[)##6?'<'#&3CK@$`n[t##ro[DZ$x$,M-#^o4G##d83O,#Tn[*##'.,:W%A/2U,#ro[=m%-o[4,#_$)Q%&%B4Y##So40##(^7mY#w#'V1+d5%gP#Ih4~R*WB6#H%2c#^f)`Z'mH&>R+2##"},
{369.85,198.62,4.32,4.66,"O$*1##vZ(R$'Hl#S]*WT5O#%m(5Sd)`l&gP#b#&[)-];ATG#Q]5.##WT5D##u-)@.&;B]2,#%C]hQ%AK5=,#/R*#h$CA]%##-954##U&3>##f$+s##/C]{G$[A]R##.3C%I$yu'D%#gA]3,#eI)zH)*[)$##*Q%E6&a|<*.*~^:[G#Rg6mS)bG$E,#SF](~."},
{355.53,213.71,4.40,1.68,"?7M+d&-fY<,#/jY]P#B.-5##(m'%##beY.,#######MR(+u#nEBVG#VjYqR)ViYcP#td*%/'gK52Z#$fYH>####0##6.*6l$,e->u#d@?sg5teYYI+<?$>s1m[+~10R/2>5####x#$-d([P#=#$4l#r?'Dc%Ou$Q#%%Z#{,&=5#2m%u>%fY#:##jo0;l$95#"},
{192.25,236.44,4.33,2.61,"vb#W@+&##lY#P##w@.#$'[P#4,#9f(~-)#l#78'v)?95#XG#Uv$9WT/}Fc>$*MO`JU0d(hZ$4/-~0')}?4R)wn-]J,D6&C:+uu%|w'ANU{#NCKUZv(tw'ANUV^9###Tx&Q^3eY####P?#e;7D>#&##y6$}84>u$$##BZ#~2;hc'.,#.##@-&95####K##2H&"},
{29.46,262.64,3.90,3.23,"vP$PG####_G#)9^D>####7##bjn#########By6#########rY#qb#+##u>%rhhZP#&##lP#8jn######.##VM@######2##F>#H>#`,$9H&df`###A5#Fc#Jjn######/##_M@######,##wb#4u#gY#gG#G=I###.,#o#$Ajn######,##o(>#########"},
{254.48,286.19,4.21,4.73,"su$0-&QH&LIP)h,3y6ev#h5Nd]U8e.gG#fI'B~.^l$j>%R,#3T'16QzY#m~.@J,Q[K>m%;y3*^U.R'l>$y$&ko2j>#En+H#$A&(*_<Lv$z39V.,JB4###kfGxZRH,$###C8CTf4###^5#Tm%OG####M$#T^Ull'######(68uG%######2z1gu&.,#/,#p>#"},
{254.48,286.19,4.21,5.72,"s.$j[Dy?*UG#`V&t99=u#1,#4m#W7.hZ$K7,Bc$95#'##Y.)dA,5>;Fu$|:(u)*.#L|&4>e(n7A-|DW>$W##OU0ZP####D,#z98Cn&ZH(v'&mU8.N)~rB=H$[]Tn,%Qu%&&#)V:######v##^I)%s3m[(/v'D_:9s0EH&wn))]T<#$###DS#Rp6D>####@##"},
{528.36,312.59,4.50,0.05,"X5$7#$%##2,#awZ.,####3###+r######%##Rx3#########&?&lY#0,#-##ZLe######-##C+r#########A(<#########TZ&Y>$###F5#xnY.u####'##t+r#########f_?#########7?#pl'######|=H`l&######g+r#########}U<#########"},
{465.42,335.01,4.56,0.56,"zG#=v(######0##pf0pb####}5#Fw,PG####D5#C7*PG####pP#R.(M$*4,#P%)h9IFz<M>#$*WFOBB-'E$%@$(Z7).K17Z$W6#.(WI..%##}w/G`*|%W1f,t'W*r68v'?](c]+|u@j[,P>#j;'A%WOG####wb4R/1U#%/6%FJ,}6(ZR*&6$>K'0J+3T..,#"},
{108.58,343.00,4.05,3.63,"*w*5I)z.,9?%fz:NI*uV7xY#*TZ{k#dG#7/%Ev)###;-$t@+?z+0Q%+R(%H%mM?y>$psAES,6VZ=H&y#&6K(x.,Jv'jv(z,%D7(#l#O&$i~JCR*NQ%6@)(YZ0}CL-'yZ%(WZL/.0@+.,#:Z#W-)l,%7%%'E6II*>IDr5%;;3*6#cSPml%V$*OQ#j,M###0,#"},
{108.58,343.00,4.05,4.97,"5o&mp9$##bl$+S',93HQ%M%*o5#;&.rl&X>$0##%e$Vo4###Xy#`AYNl%2d(>|)EAY]5$d?(&gQwZOhu%c$%&m'iv#esE~l$Q#$$v&J-'v}HBL5rS0kH$T`@kCY.v'Zm(=x':S-%d$ACYB?%x,$SC6*6&0?&V7,7D46u#-%%vDY#d'C,$20'-o+tQ'S'1q5$"},
{180.37,344.93,4.32,1.81,"h$%Zc%:5#4,#ziC=Z$sY#Jx*~%/U'+z[&rQOYI%c]27@#gPOqy1`>$g$'XZ%JQOhG#am*lB#ln/-i*'SOmZ$K%I%g1g@+,R&Zy50##/+AqP#XgNmP$8?&W,#N205c$#R>;L9'T3B.'UA&m8Nj$+x-$SA1#?#.116#$%##XQ%kF8}(9nG#ju&7,#V#:yP#;v("},
{180.37,344.93,4.32,4.68,"ac#@y2{k#'dF-e$9I+LM7RGH)##r#%kYEOc&^P#6R$&-'8l#;K+#Z<KWB[l$=COkn/`03[c#+R)?H#/CO'##ry9hZ#/94CH#Q7)jH(Hk;.PAt>O1x$6~.FX1M?(vg$EuNS>#3`<s,$xy88,#vg$iiDb5#W.,uI$]k?W/3L[';5#2e$$T1A-'AZ#>l$lB/G,$"},
{164.65,352.01,4.73,1.27,"~S1}Q#|/5?m#`3Eh,#v$,5{)d?'*d<n(=CS'9E8zB69B&C?&n%/0c#D#$X7'=nQf#%n'4@v#a[+g.(XqQ`V9oXHcm&Me-ms0|Z'B6()##%?$i=>rlQIu#/H#mu$>mQV~&.]1hG$XoQ66'2T,V5#M:7(?&/,#D##CTO.,####]&'nHQ###+##iR);%-###<6$"},
{164.65,352.01,4.73,3.80,")l#Ld)R?#&uB,##c$*N1&6D?###kP#}=6qx5?c#G>#g.%xOK+.'i]KRv(|J,q%Hyy8(1$(L616'>?#=s/j:=###+##G8&$z:CI*bw'*{2+QE}uNsJ'h7.dW.Cv(MV(I`>qY#Q>#>$'lQ&A,$7%)gXDcc#w]4mP#5xJm-+yR(,l#pL-Gw.]Z%O>#nA2,Q#$@*"},
{164.65,352.01,4.73,6.05,"mf2+Q%/,#}Y#^l&)##SZ$f?(yP#`,%:.&1K2w,#nx4)[#+p6)SGZP####9##ZlO*Z#*Q#lkAo,%~(6$$$%/Q7?#&80ED(c`B<.Q###'##LS%O.QuR%hS0p#%v[(xUPhE>1%,{Q=o;@qv%0@(q6*###>Q#2x*20Q)?&cw)hG$m@(xe0t30`XH5f1x/'Fv'#}="},
{195.42,353.82,4.42,3.53,"u7-8e')R*(##QvT#e%PI,3J$h$)A(NL/1(e(<Y;Y2@j$$,d%g_;z,%J-(km'fvT%?$tp04@&ER*Zn)kc9xwT?NAQ@'$Q$^7@5_<[P#BQ#~'/XxTyp9:Z#o$'l/'rvT2[#zK6KH#09L$?&Y/.EB2pP$kP#P6&IZ%Q?'v>$qG$%?#rK1.['/u#T#$&i8OG#jR("},
{195.42,353.82,4.42,5.70,"iH%#s>Rd#Gi@(@#GI,.j0K`C.c$$l#-=8U[+4l$8u#wG#pp9z8(L@J=h4rZ'ILO|~1r'1Hm%p$+27%LSEuu&g#&;Z$IQ%{>%Nv%IJ.+r)z>LHaFs,;Jc$v|6.n*k@Cq7-A5#+%*`/-XQ'iP#z6#A6Ja#$vQ).,#k5:$?%F80xZ&'(1x?*W?'{u%(ZA4l$.##"},
{538.74,352.91,5.38,0.02,"%o.<Q&######G3n[P#######%6R#####################f'7bG$###B,#73n.,####I##/B`######+##############6^8######5m%C3n######L@#j~~######M##############4e.###$##Ei6_2n######/}%Y'9######,$#############"},
{436.54,368.62,4.51,0.16,"f~,{/5###7##x>^CZ&###P##hfQM#%;5#2Q#47+Cl$O#%cZ%kS2###%##/-$=9^######S'#u8^######^&#A@,(##=m(d>#[7-###mY#dl$;>^######^##W;^######~##M@-###bQ'T,#vP$###;5#p5%a:^######Xl#4:^######>##<7,o#%.c$2##"},
{23.32,392.18,4.04,3.25,"B,$###D,#f?(8Lb######(-$Z=l######$##%B3#########yP$O5$1,#C5#>Ca######1##D=l######@##_2A######)##gY#.,#tY#v,%Hgd/,#UG#z5#_;l######8$#@)A######3##4u#.,#R5#3c$r#Q###=,#/Q$9=l######(##;h;#########"},
{395.60,399.52,4.59,2.45,">&#r%U######rq%|f5Y&/C6(<o.0,#Wd&47FSV?######-D'?[#4@*.,####*-90Z%###vY#l)U###$###A(.V9eY####'m#S,$PG####t$)v#O.,####_V+/%U3c$###7C#M//0K1###9##Z83######/Z99@,+##$l#l)U=h<B%*|Y$Pp*du%CS,###1%("},
{100.35,410.71,4.60,1.28,"RiA######*0(Y[P%?#V6)}|3x,%&c6X)Axg2X`5b2<.]&e$*WT3###%##nd%EgSdZ%w/26d#q[+5K'+iSsx0/GI<d&f%-+=/t$*OG#=5#Ql#pHCgtLGl$}5#WH%&eS2g*pe/fl&<fS~c%w9.)##>u#rb#$##>##:T0X>$###F5#_h;[Z%0Z%I>#)y4M5#S[*"},
{100.35,410.71,4.60,3.97,"[>$lc&+n#1dNJQ$IZ&Y8#^XI.H%###gQ#[A1cG$###;##J.+@7(7BJgm(TK0kII*^5nT$9FFt#'<-#rO4C_=ZP#=,#d$(j5%kv)F/(d;7`#A5-O0.&pB7^`-;-(:(%1-O_5$?5#JR'>?'###P8(J<CBH#p;<??$8/Omo5=K*(##U)-QOI/,#+##fL1W>$###"},
{100.35,410.71,4.60,6.12,"xb#RG#1,#4,#95#&##Xl#]R+###$##-7%iz8######z7%g,Nb~-yY$###=##V}Je5$l>#_dG)v%RL6Mc#MzQLH#>J0nj2qFDF^4Yc%>u$7,#>xQRA')A.CH%Ce'dzQE(6H@+=u7DkJnn*UH%EI*bG#Zd*J>#?yQ<c$4/+n,&?w)]e.i3/+>Km]3XA(FI*3G="},
{239.39,409.78,4.25,4.55,"nl#VC5D>#EQI{w%5R*gf0_-P'##Yl#gZJ>R+%##|H%/c$###Py/U5;n1=S#$N:P_7-5B1HZ#@d)2m#(:P?5####Tx$WT4###JS,S6(pE9bX>#6PZ^*8v(Ti,lc'fN.CjE3,####wd#1YJ###H(&Sh>N,#*S,HR$1=8b./V6&/,#{%'6~-tu&###s5#=sAPG#"},
{225.27,416.70,4.44,1.23,"yp:######t&&(jBE##X@-I=1w5%Ct2aGN(f*V*;T:7[n&_Q&lA1###(##HJ(JxREc$^04Vm#;[)-%&&{R#V84sDXQ$^7+Zu8Yd)Y,%1,#:Q$Il@2?Q[,$6?#|u$[vRU&'z83VG#`xReQ'1q3K#$#6&sP#TJ,>##KL79]'K~/Gm$J}H]?$^>$OH%s&40,#|l&"},
{225.27,416.70,4.44,3.81,"g>$@I+3$#xuD<##O?()q$iMAp%.UG#=y)=&2Dc@3-($$$=R*7~'9^O-$'*S)rfL5p6Z0$8^5RQ'F?#Ta14::)@'_Q'N@'~?(8@*`7'B'/T?@lZOAS&7/1(2(yl'__(GFE{b#<,#0d&N?'OG#id'sWAy>#aL6yu$4]OB7-xd%*##nN3DT5###5##Z:;.,####"},
{225.27,416.70,4.44,6.04,"el&VG#Bm&a5$Z>$&##z#%BI*###&##(A&;C6######`S$6<CT04De&_/2R5#=lN'Z#|>#@-EUu$Ih74$$RBR],#}84b*/%i=p*Cem$gu&%v#$AR5S%$]1r#%Hm&9DR*|>,R)/b3Q?R47'XH&+~-(##KQ%Gd%WBRru&'B.<5#i-&Rw,+c7y{COx0?e&pH'oGD"},
{184.32,422.45,4.25,5.45,"*i@4e)`5%]0$,8&oZMD[$yo1}~%wQ(cp/=OGlP#q5#~t?r#'AC8uR&T>Gl[#g_-0|<1r5*I(3gMNd'Q1<.d$Hu$&[#,cMeP#2W,mv+6g,O5$(['z_={:(okFrU;m@D@$)~g,bQ'WU-NA1qG#,S$)T3+H#r5&F##CdM_#$%e-cl&ePC4,#~~,=B57-'###38$"},
{107.82,436.18,4.17,1.87,"B?#j^4@?%guQb@%H?(^2*auQ1m%|b#kN6Zl&Ev'0,#C,$###Tg/I<7aaDk>$8zQ(d(2D4pl$-7+],#dzQuG$am*oY#($&-l#gd+>l#rt5OkEwuQj;*X$)?`0?-(<{(hFGO>#Q#%^#$qu&@5#<v%>+8R-%qd,D>#[s/@Z%Fm(Kl%Ox&PG#8d%cG$=,#_#&F,#"},
{107.82,436.18,4.17,4.76,"(c$6,#[P#oG#PG#j?%Fu$,8+D#$0d(/,#{l?Tu$A,$E5#/7Gil&~G#7Z%9Z$e=IdG#6H&]}3oZ':*1RC;`JEuW-XjAsu&(p0AI*kG#<S/B5#?TRSu#/e-%6#HB3.0&{SRwl&USRed(.x.G]*c5%$##`v(@5#nwDp?*Q?&D5#tz(8`B];0~%/Q~,0^6T$%~P="},
{22.26,437.58,4.85,3.08,"3v('##7c$|P1,%,D>####`0E`tL######}4.rK8######k%#uA0mP$2##z',z7K$m(###ym%,|^######H/#o7.######W##mG$.,#:##hn.DIUZP#$##V~%qx^######~&#;~/######?##bG$###@##'R(PmV###)##t[%=x^######Z$#fH)######(##"},
{37.35,451.56,4.72,3.97,"o.0T##n;BK6#.xY###eY#`'#+}YbG$###$$#xs9.7,######d)C$##8E>%m#ZwYOG####|%#*VTFC9###{##D|)6GM######<T4###Gq+UA0)zYa5%2##b7)U01k95###N,#Wr-Am)###/##eY####G-#tiCg~0]P#0##Y*;U[*oP$###GR$4-'######A?$"},
{13.78,526.79,4.30,3.01,"*H%D>####bO?Wm*OG####V=;u'f######XX&ay9######O$#N&)DI,###Xu$G,Fie1###RZ$f(f95####M%##95######O##r-'gu&$##oY#)h`=m)###9##9(f######q$#r@/######4##I>#YU5$##/,#&*f.(6###9##O(f.,####R$#VZ'######*##"},
{137.23,576.90,3.59,5.62,"jh#J@WOG####Xb1ac'###$##}J-+u####@5#Yu$j>%$##4u#ru65T2~n';o/64Y{k####76$X-P######w>#8Z%###/,#g,$|$,(##Gf$[mJe.Y###&##O|*C.Y###$##T-#)c$###%c#%c#6#$C5#I#$>+1eo5######P`'k84###/,#r$#pb#%##&l#}Y#"},
{257.29,597.91,4.19,3.36,"c&0we/o5$3H%y8$cy4q,%M5$*K$.S..,####,7*rb####)##87'u@QP5$?5#i/)GqP-u#%##p(RAw-###6##8z8PG####*##G5#ySQ?_3/I*2bEO$LlP#B),R$RfG$###9:#'93eY####N,#.,#M>#`z'5$RjH)RG#X5#p(R&x1.,####CD'5Z%qb####k,$"},
{407.47,607.97,4.44,2.87,"K,#lf0,w)G81R`1-w,/$%C$);k9ZP####O>#xv'D>####I5#^,##GCR?'p,&R0,_u;5[(^>$[2T5d(###?##'02.,#####v#Su#:&Sb;7h-*$E=U4AL>#'1),.T+l####~:#Vf4eY####h]&SG#xb#p)*e-TAm)OG#C##[2Tp:=######1`(Ew*8Q&###e,#"},
{125.89,628.97,4.61,1.34,"Z_#&r@###$##aG<^Q(###m##Vy8y6&(c$<%#p5$l[)O-)E#$PW/_T6H](nJ0UaX1d(KQ%mQ#iV?j;*Vp9v##X(4qx/I#%,##*C8&~%jV8sf+A~XC@&8@*+z#fK4tW*1T4*##KuLe#%{k#7##nU;j#$qc'c7%C6Olm%%d(+~#_@,59(dw-gH'y5P[>$###hZ#"},
{429.97,636.76,4.35,1.27,"_x#4YIL5$###3+14y6D>#%##v~.mP$###m5#1l#%l#Bl#RQ&5p/&C1GM=z>#&{R{I,1n-z##&t>##D3l$/##_1*,1:%##nG$<x,WS)]@-6t6(vRz5#~d+DX/Cy3q(.=NB0?$dxR&w+]>$-[#f5%###%##,9HL5$/##b8/6*8o#'?##=i<ZR&*bK###?5#~C$"},
{50.64,650.35,3.99,5.50,"/##+x,z[&YS1,H#W6);7#@z<<5####T&#H'8>,#_5%cS%KQ'r.&sDTWz<K>#sDTjrBa>#H?&aJ0###&y#ym+######,E,6#$^[)'x(zLO;RLlATZu%uP#*@='1:###i~&Pe&######|C'Qu%`-%k..I@%nM<Vf)4[*'##jF47H$KZ&fZ&*?####2,#'S(A,$"},
{93.74,706.26,4.53,3.99,"n#$Y0'sVA###%J(T~*=~.Tl$|Q(`,$Gq3LdO/:3M.-$6#A6P/K.KxG^~1o5#W7PJv(vP#~B*Ce)1[)%o$I{:FQ#uL4e%)`83;K2b~.Hv$N:P@?DP6)&%#`}CI~&lkI@[$<06###%0-IS+2T1######-$#Y7PD5#OG#z&#b5P2##/$(+$#f5P###yP$7,#|6P"},
{455.99,712.71,4.42,3.99,"#H##0/Og.;Q&-9)NM:@I)Eu#T~E7w-%##.Z#Og)8.-$##s>#)m&I8(*G=AS,{=JCI%/g.|1/[N.]-*5(&Wo0;}3Fw.K,#$6$ry4U?$E~G7w)L-O)-#BlDT-&s:;6,#HN/r80tK6###xu#u;.bJ0>k7Tz={d$TYM7I#]x4=:$EA1TG#N#$Bz/_$*Uc$nG$SB,"},
{51.01,714.97,4.30,4.27,"|G$Q{2###$##Ec#~}A76&QG#^e%$tG1?%.e)J('v6+D5#}#&O6#;44_/3###>E1R|D]#%dG#s,OQ?'.e)_O3L8*`Q$2s=}~-S5#i$(u/O,u#)0OSx2%o-8['V0OG6E^$*w@%`w)Gj061:G6&k>$T#$5E9.v'(x.cG$t7)UnNB~)Qg4kY#:,GEZ%7f)S#%[S."},
{65.61,717.55,4.51,0.71,"~Z$YA1Fl$Nz6$e%F:;J-$}B4a.%-C4<@'m,&[>#p~(y_<?Z%4-&5y/hW6S]Qen%#]Q_~Q0d)W$H@,JKQ&K6$tQ)a5$I/+QJ*A%+f>8a<A68+w$+O0%d^QC@*I[QY#%sH&SA'/o1?5#/?&DL*9R)a~-9e'MS/Au$*##PN1&`=uC=###|,$)-<i/3L#%$##rS%"},
{65.61,717.55,4.51,4.08,"_P#PJ&,M<,c$tm$tl=?f3####r.%z8DZ%Fc#jn,X?(^v'2a7:#$C/$VsC]P#3.'*/+%KQv>%KLQ;/1rZ'ue%Q7Qb+;0B2}M*rd'.o(=S/1,#($'>H$]~MKaB*cFFd*T-%0JQ.B,LRIl>$V2?.L18A,'$(#H#}w.Xl%^H$H(3{Q$e80U$%Nh<jY#(_/_-'aw."},
{507.69,735.71,4.60,0.49,"oH''g3dG#%.*#g-9l$b,$4H$G+Z95#S#$95#}-B#########(-&ZC2bR(uT4+D88%-ku#tZ%u*Z###*~%El$G+Z###$##/,#^u%8/+EB+h@,k=JVu%L,#r9)d<B###GI$r.,G+Z###:,#V5$X?(u-)Iu#wJ,;x2;5#'##/K%pc'hu#PZ%~v'&(ZN5#_>$J@$"},
{46.85,743.30,4.85,2.01,"10,OJ18,#P#%d(-bm*zY$%##J;SL,$[P#.,#j{6#########HC*#EBb>#'/+'pP,@+3,#h6&@;S]P#]5#=d(J;S#########1-'fG$`A'sX<B6S7#$(##B(/z}I###;H#$7'd:S###,##oY#7?&iv'3165/+{T7L#%/,#}^&z$,1,#H5#UI'P9S_P#=5#8u#"},
{433.73,51.96,5.59,5.02,".,####*v#?xRE5#Kl%$##&{RLe-_5%###&{RMvR###^P#d)*.,#D,#-$&WtDb6%r===B39L0WyK)L7(l#Pw(}xR###E>#~>####)##MR'F6(EH&'v$@F8M:6$#JwP$Nu#~=23xR######M#####_H%9Q%ZP#Ev%q7((@*.Z$G_.S$*###%-$b_6>u$###(##"},
{433.73,51.96,5.59,5.78,"A`Asb####?(#RJSLv)###?(#e3.Iq=&##/##p'##955##nP$GHSoP#?u$Gz#-(6gdD;K3m.$3~A?*CI,$~u#O;*h&5###$##]98*##)6&oM(dH(S-%OF9?_9qQL?d)/$%R{3+o+w(=.,#(##A,$###$##)8,OG#%##<Z#u..yY#_?'TT166'D##_N=LI,###"},
{514.33,68.55,5.39,5.71,"tR&e<AqR-&##U*-<y7###,##]|-:V>2,#P>#kK&)B4/,#)##pw.2Z$`G9w6*7b>&A.J##kI+vpT+d(}Y$zH$FM8/w+|c(&8$rP#j,&v:$-WB%[(,Q%O%#^nTwmTR5$I>#N}17'65u#L5$9D%)D1O81:c#pQ'xc'###q##j:9?z<###+##E<+*L9######L&#"},
{347.33,91.78,6.08,3.14,"?6&2Y<}u&uu&FN]Qz7@5#N#$lqVif5nH&`G#Z,#}T7][&u7/(8./AKzY$J>#JL]MI*###@##82[pm,###%##][%,o2######nI.W5$&Q$EC3ZJ]######p]&x8~0c$95#D##)x+Ov)95#&##1h;{k#&##RU+>FI######Ya*3o2%Q$CH%K],w7,Kd)kc%q&/"},
{389.83,101.24,5.51,0.35,"$##a/)Zo2I#%.dNdm*@5#G-&]0^######S$#i~1######&$#$##c%(5V95l#G8UeZ'$##oS)52^######([$}y:W5$###7$#%Q%###{$E%U*)^7d>$#Z#=DLr/^2l#P#$zW6RNC)[%Q#%`['`5%PG#tb6E{28S.g$)'m&G:*sRLd$)X5$tG$vN9Um($Z$N-'"},
{505.21,113.57,5.69,0.37,"Lu#rP$EH#kTSIc%Od)95#5@>Lg4M?(###(K&FTS######4##%Q%Q##~-)B`7l$&FxEXOC7J+nVS1sARu$K[%CTS######2##.,#+##cd'7/19Q%P-%ad>9RQ'$OMu$-d$h'LTHK######y##2,#aG#mY#Au$,Q$E>#Q##lI-o#&###)##&/,8e,######:##"},
{505.21,113.57,5.69,2.20,"###aoH######,##JyK######j##vL6######2##nZ'.,####h##Nt<,$(###&I%@_Q3*A;c$cBKFHNY#%tZ$z@+^u$D>#(##zg#[p8{?).,#[.'sm&xeE.=CJ[QCZ%b6%5@?Yw.jP#0u#tG#r+5#J-O$)zY$/l7/v(J##|Q)-].;c%,##2I(k?(###(##_P#"},
{58.71,115.20,5.19,1.89,"5(,Il%Bl$###dW6L5$######On*O5$/,#0,#(l#~P#0,#L>#op'XlHyU;WG#2)S$'6/,#M,#2jC###$##4?$pb#$##3,#K-'nv(t@*BY?']R}%S1H&Q>#D$<b2B######Lg*eY#$##B,#Ih9###$##sP#2)SA,$######2)S{k#######P:R######$##?$L"},
{216.97,154.36,5.39,2.93,"E##Uj@DH%Oc&~v#H'8j>#w5&(d&W>$o,#7C8.,####{##J'8t##WdNqg+*WBfN,.$T<$#'i>O)T8m)r5$O26s?*CQ$/g,&i?###95#_^#-$Tc$)I#%e&#2&Tc%T$m(k5$l-=X[)#S-yf+O%-U##g7/>$#::;$Q#:$)e$#d+L~$&pm,K$#[;@H6#4'7W?$_@."},
{216.97,154.36,5.39,4.54,">:;.,####v##wq>#I*Zc$0g->e-|v)z%0x8*lw0&-#nf6`H$>q:0m&L5$~##(r5HiTd6*y>$#7=WeTfH)Eu#YWC.d%ed,9$#ku&Gv#YT6h>#ph=W<,|cT{##%fTZv(0Z%/&#M5OmG#[?)K$#p#';##?N=#-%j_@j##5&TJ0+zcT###PG#f'$.M>G##T7.1$#"},
{247.19,155.22,5.51,4.95,"6~*rw'mC<$##8$F5N3h7V###^;Vj,&QG####4s@###Ro1###B$)1,#hV:mP#/+F$##F<Vbm)s9V###},%+w&%@Lf>$dR,'##sx1{k#'v$@?&4J-l?)f<,F8VlaHB8,(m$~R?Dy+qh9<c%I5#M@)bG$###N5#w?%?C9@,#6H&eG$Up32,#nc%pG$(7).,#nT0"},
{139.81,159.61,5.01,4.80,"~>#m>$hG$<5#g958%*xZ'b>#Wf^7#$kQ%8X4>u$-##)u?>T2$l#(##g,%*6&,sB`u#q&/$7)dg^H>#lw.2I%rI.(##'h^Y5$Iw-95#L>#5d&o:=8,#^S'17LTf^%##m~.<B-]83-##}f^$##;z995####*##n?RG##@x0Q?&Wg^B##]x2@u#uS2###4g^~P#"},
{373.92,212.80,6.05,1.47,">nD`x1bu%wP$YS,$%);6&zd+)Q%%##Fw&/&1######kQ&A,$[D_YS(d3GK%#OD_)K+h1=vQ%](=G##:$L0d'###$##K~,Y>$mA_D,#TB[R<'XB_;c#rsC]9)p(>tZ#t4L<#####)##c[)}k#_w/[G#-K.^B+YK5rZ'fm'F0,YQ'Ow(2%-0##E>#iG#y>%TG#"},
{182.77,248.96,5.45,0.17,"###Xu$oG$hY####9@'(030,#>7(F?;(qQ~p5IpQk95*[&z`3###TI(DZ&###06%x=:F::K>#B*=NT29)0?4C$8/$l#$H#=v?XG#;Z#S..###.?&XH$^bDcp3/*B>Q%Ov&tj85R&*n,>#$D5#E$%O6&y5&=v&M&32v&l#%.eA7H&r81M<-=qQ9I$j7-U3.=$)"},
{182.77,248.96,5.45,3.43,":X/9U9PH#zS0){'dzQOw.B92>6&JF1dT666%u#'B#$FQ#%['BQ%?#$h6&8/00@'6G:qiAl,%e5Fz&0/?&E?$W..###B5#?c#.##xM2ud+Fu$rK*{4DPj=^&15M;J>#26%>}7DZ&######NI'Eu#rr0PxQo[,dzQ9tAGB/`5;@r;$#####5.&d,%(l#95#5l#"},
{529.74,259.92,5.43,0.06,"6l#V#%<5#/,#'/YfY####'##J4q#########d./#########V,$a#&3,#:5#GUaD>####%##94q######9##FB6######(##|P#w#'######~(d{k####'##Z3q######~##I97######.##k>$Y>$###$##/~WOG####$##,4q######&##[&4######$##"},
{202.19,328.38,4.93,1.06,"xG$aR&i39`#&c@,nY#gR'QQGW>$%##VG#@Q;Cc%kw%ie.C%&&K.jv*:8F%-%5$O:5#*p,}$'zn0-u#e,$Tw)###)-$F_03v(U81Gn$m%O1h)*NA#6#Aa?S/(W/0l>%lG#0>7qb#E?$E+5.L3Q#%c45UA2NR&yiEm##:m)iM*6I+}z0P7-&I?C%(h4E631c@*"},
{68.38,356.57,5.49,3.37,"eA&o08QG#/,#xZ&`v)B5#)Z$5e*GZ&$##^#$O&2###.##Vc%-G[:XF[>$Q6#N_2&vAgf2%-%QpL~06nG$bH$(p5######16$PB[k5%Dm&6F+h@-FQ$)O2WILxWC>#$|5$yv:en0######+d#XB[7?&=?%zs1cl%d>$ku$,f.4-(###$##bw(jZ(######FH$"},
{68.38,356.57,5.49,3.88,"Iu$iY#,##>Q%#.(^-).,#U>#gR'6?'###%c#]*2######3##2J#?`>a5%G>#w[%VmA,p4+Z$89FNy6K,$C?$,p2######=##3jVjxW]u%Z>#j,F8R'IV2WtCaPKjG$p,$e#8G81######i,#T|W`$)']*W?%gzWtc'=$&)U14m(F,${b#U7(*[(Dc%###P5#"},
{460.23,394.30,5.18,4.12,"$##CA%xm,###-m%If,h.-1,#U5$N>#?FXYZ'###+##Kb0{1?>#$+A%heYI>#PvNh%*od+ZH$/J.[P#h(%um*######;i%F=JZP#/##].>]v*5`8Qc&d1'I%,^g'ne1}l$1u#[##Qc&^x)6#$######K%#mJ24##+u#.2#_K7^[#$m(p0#G'8Hw&*6'FZ#,7)"},
{380.65,397.13,5.39,2.67,"Gz$5XFB-'3?&W~E%Q%mP#8_+n<@mP$###P$#56$'m'######>|2G[+###c5#.*V?u$###s##GbCcp5###)##fS(Ep7###:##'T3######+G2v$Vz>%Du$&b05:8v//,u#W7$g.-0?&###Bh&Y$V2,#.,#8`(7I+LZ#_@-M(-3p7sY#c#&5U#;n.######1J#"},
{116.93,404.16,5.29,1.30,"@7PcG$+##>Q$+178l$###Bm&0c#h$*qQ$#_5wP#YV;PA&Mm*hkM95#C5#no%/?P>,#+c$>G6AQ%1N0so0#]J_~*;y2iq'.K4[97######5^%VAPTu#s$+17$5T0]_*VAPv7+-#H3A+PB2[g)x#'######}Q%6uC?f3V5$R,#:@(p?PH0-JI*ru&'mJh?(WJ)"},
{116.93,404.16,5.29,2.27,"76%ZI)`#&[G#.o2ql#~@,5><U6%%(2&&0i/O+](1.,qB%):8C8+}I.W>$8##*@Q,Z$YA.8&%Jx.v/-wCQ<p3)RG5e(>T.EC32*-4f3###~G#KoGRL;W,$N,#y:&fGMN')D06wl%%W*Uq:(q7*M4|l&H>#Ol$TZB*6'###%##}o'Xf3######G##|T-$I*###"},
{116.93,404.16,5.29,4.81,"Jc$B6'QG#AP9YQ$wG%Wy-)ZF###D##X)8;c%###k%&LU:###[o,jh1?<B>f+8UNPe-/_0|v'16'L?#8UN3u#%##]@#9QNPG#f],uo4ns:Op2EQN,|2&I(W{-_l&Z('olMXG#$##1Z#wSND>#i%$#+I(?$xu%M%&*vJ,c$C7(H>#=%&??&&$&###+##AU+W,%"},
{362.01,404.25,5.30,3.41,"LI%dR,A,$###L5#I/-P6)###ux,gH(OG#$##NA/######Y,#@2-bD>pb####wc%o]J/v(&##:CT@:895#X##*}E######L-#3+-KYK>u$###&*5Hq-Vp98##nTZw5%W>$2$#nOA######2##(YZaQ(###~$#toR=m&m@*r-%>VZJ>#)Q$^P#rq4#########"},
{177.74,406.52,5.15,1.11,"X/2,[#>?E-R$nI-1Z#4x.mR):?%Vu%yG#215Nv$s:<Kw&|c(p-+YW,|VAcm#6jD5?#=[*C[>HZ%p;0Yx/X]QQI'Xg45*-:7,%w+3w).c$C[&s~Q3?$0~-:%%sS/YW,G]QXJ-kjBM/-(y/A(+&v&=Q&3,#8Z$<k<|MAVc%<Z#e[&dZQxg.SS.8?%VbA(7*_9-"},
{177.74,406.52,5.15,2.24,"0$&:7)N5$eP#cx4vu#lR,`O9,I%pU2ix3O[L|K,6d)WT'a^9RA.t?):,#,-&f-Qe,%@K.eJ%o~-XS+/gK])=AZG9w%?91T{6SJ*bG$kM)rp9O8Fr(?9.&}>$h:&E-QK/&e'81?$gW-Cr?t956`736$)b:cc%P@)Qm'9{8H,$j%'+T1{k#.##SQ&y/*%w+1~*"},
{177.74,406.52,5.15,4.68,"#H$@-'yP$@&KS6%1H&lf,.AO/,#;l#}p2j?)OG#0H#m?(p6'~T.DX622>{A-=CO6]1?h0hI)a6)&m#=COU5$Y+K#H##8/I]'3p*3o0mE;zy3N?ObE3`Q&9i-;d)|;,~4HXG#v@O2$&Df3D##uo$~U;q>#}d*97%1b@_>$SI)0,#iR%@$'.6'qI%07+p7+95#"},
{464.27,407.84,5.36,1.10,"R&#|1?x?#pm,a(#JkLZ5#uG%Y0$Q;?######2$&p?)OG#5,#H6$gG$lJ&?K5T9($&/C[Ie-*3sV}[-=Q&j>#i;=6#$###D##t($Sz<wl':5#~[(XZ%?IM|S-6nV1,#tu%3h){);{k####7##X)SY'9.,#1##Q.Gqb#Lu#;7(V~0###&##:T%GR*######A?#"},
{24.40,417.04,5.62,3.23,"Mx1ZP#P,#AD/XgSQu%###Xv#Xvm######l##Q%.######(##37+95#(##JT0v9c.,#%##0A)Lsm######Q##xT7######$##G7-.,#(##3~'iZR.,#*##/7(Atm######2##i(=#########8l$.,#$##$?$u)C.,#'##1c#rrm######k##m^;######.##"},
{170.12,462.04,5.20,4.68,";##98+~#&###v~(fR,95#X5#`l%r5#C-&x#&A##,7'$U0;]096&nK)h<G%##[^V~%,[?)0##:4GC[$@6(87%N`C0n%8@+v`ViZ(o##^RTtP#G[V-##J[*7R#D~Vo5$eA0wl#J^Vjx(f-R*A&(c$###TD38%,csH###9H$e'-xC6Bc%Jp/1Z#E+0q;BkB/2Z%"},
{183.54,461.40,5.21,4.78,"nH&mf,l':###R$Fad*uG%$##M2??-$xc(r@&wC=#J'(e-xTIfm+j,#A@UQ>#Z@U7,#Jm)y>#dAU$H$88.Tl#vCU0B+5AU1@&}>&###0b:3@*q?U####v$C0,N<<u#'fg.d,$b|-p?U^p-XZ'$#####/@$,J/Id*###J,#2g0q-).Z$Dd$Z,$f,#CR*&0'b]1"},
{524.17,518.93,4.94,5.86,"Dm%~#&)##vu&j~-OG#N##V`<^v)###`%#X8U0{@###4##-a,^9$0{@###$##kQ:eo5%##BZ$j:UOG#g##.i82I+###*##^:-tx)$95@H#BQ&w9UBv)###(/#=7U######9(#n?*######U##z6+.,#T0%Xn,Z6U###:##=D%XYN######1'#+u#######,##"},
{103.62,536.91,5.31,3.04,"&##@A-Cv$/v(3##S~+k6(L~/###5z)UR+(6&###Es?7l#KQ'K@%sDTHo3<5#@]F^dSWZ&|G$SI+ql#af,r`D###_#%sR$^>PT@*?p,EMQ~-L.BT%v'Y>#>$;<T4###m,#OU2######2x)GI,###.,#Do$N@T%Q%?,#ZG#-JDZP#Y5#Lc$J6&###7,#u?&W>$"},
{92.72,552.19,5.04,4.80,"######V*=WZ&(2?###b*AZ5#XcQ######O.%=$)######+f'%###l#[P@U5$Pp8###D-=9m'=dQ###5##Md%H]3D>####M?$OT%/y62e'A,$V?'###~z%H_<arAWy/*I'GC3q',=eQZP#-##LhQkS1gG$~>#kq:###_-&y$*/,#aD*yeQJH'4D.JhQJ`=5['"},
{42.37,661.75,5.69,2.71,"F5#D#$N'&A074R&,-&^R&Y{>ef,$@*R5$&C*{Z(NQ&sP$Ep&16$fZ%^e,d@.;[%.6;lRP~A1rUP`;<FZ%|I&5{><d%###m?#{%#H06bd+###YQ%k&.rUPm2>GRPOe-3R%Yv>C{>EW0###$$#3S##f22,#Z,%a##{,O9Q$?Q&*n*+RP@5#e6&B(;jB1###-1#"},
{42.37,661.75,5.69,6.23,"hY#&m#ecOTr6=c$VH'_P#jc>1H$u5&A##r+I.,####|##<kK1I+t##Zi?l2/Td%tdA^dOte-mgO*E?y5%{~,k..Fu$x>#I96A7-&##8`:8x%DZ%e?%mgONW?KdObB3|['<v=iw,gcOaP#;,#uG%W##A'6Oc#.,#E%#Pi=3[)2Z$>f)Kf-#J.cH&l&3]H$8n,"},
{486.91,694.49,5.28,4.72,"r1/oP$a8(519J,<yY$7Q$%o,y97m$(EH'4m#Wc$d/-Fx2###&z-ie+xQCjm+M~Pwu#_g5a##K^8X=3y(>:,#q{7((PH%-s,$L16*^*GK1S.)tZPsv%~m*6%%H%,FS&c$Gj./#q8y6';:.ITIW>$Cu#*7&?kG)K2PG#K5#I[P>I'[l&_,#(#NAm$E-)A##P|A"},
{115.77,724.12,5.55,1.66,"=l$iP#j5$rb#[>#ow+*Z$95#{m%rd,.,####-aCV,%E>#'##XG#Zp*{,'###`R&$rQx(?D5#M)S'[O6l$Z,#ldPI#%###-#####`c8Z$)SH's7.8*+xfUQ/,SgUKH&/v&t0'0~P(c$###6##QH'=e(XQ%9A+c#&L?%g22KA/i#G/,#lG#D8,TcDmP$###&##"},
{115.77,724.12,5.55,5.44,"V#%RH$^.)QJ02T*VW?6-'zP#REA'y0zu'H.#JD@/,#.,#8%####c>#O<+7q<Zo+O02=WTK7,nqUMm)C?&:T&V(=.,#(##d$$###$##tM)S..t,&J>#T7?3eT8RT;5#3[$zMS?/1J#%%##=H#######4E)+u#######@A$-T2PG####S##NK3~P#PG#N5#jG$"},
{492.57,731.14,5.71,0.89,"0Q%nS%Lm)(d&mx)'g0a5%4Z$5O(f^;/##fl&xt.iZ($##rb#e%0Q:%'/0B^.tn-J0-;U5l_8'@GWI+wY#Bw*m&A###C5#_.,qGDL-&9_;@v#k.+fn'fY@O~-lf_)H$FZ$7'+tL9G5#~u%~S(5U1>d&7I*@[%`7*&e'76&|B1ig:X>#~P#6{'yd-%##:5#rE'"},
{411.41,731.92,5.53,1.05,"###+##7-%*[)&##l6%$1,0S/f,%$6%NC+0)>/m&OG#I##06M###G##*~,:Z%4d%g#8fvP-H%&)U#D7Y?'PB.{84###)##AY;###&##e7*[,%z>%_Q$poGu4H&lK#Q$|?%0gIOaH######Z1#######vG#=c%.,####5?#/8/Ve0###,##L|,Z|G$##95#O(#"},
{411.41,731.92,5.53,3.53,"###QVR######*##QVR######m,#|SR######[5#2J-######@%#AtAdu&###0R%jzQd{@Qu$/LLzbJV#%%d$N8.{k####+##lN.se/hc'X>#P7+qw(huB]mM|uNg,%2-%^wB/J.{k####K##'4Atb#t>$xl%5`9ZP#*##KA-`c%B,$$##l[)k>$vG%###(##"},
{41.52,17.80,6.02,4.51,"3$#TB7######Q##Y83$#####C##Yl&######%###########w(#2%Y######Zr)1%Y###&##x==0e.###m,#-c$######4##g[<+3D###3##*+Y9@,###x##m%YZP####e/#>u$######B##I&Y###2,#-8#r%Y######?'#2%Y######.%#ZP#######%##"},
{119.75,34.96,6.23,4.83,"N,$uG$V+3?m(pZR.,#tZ#Xh+/[R######y$#qZR######P##%##Qu$tn(Lr@}_>|Y$mu#0+1M[ROG####o%##[R.,####E#####:##*U/Y&4>M0?M6az8ER(Z_RVH($##Rc#Q[R$###########&##1f*Ru%_-)$H$+r0%5E#[ROG#)##AV-sZR###95#+##"},
{137.60,41.46,5.98,4.76,"(##'H%vd%d)C{m+`>$zm%D+;P#L+u####B^#7IV######D$####)##L'-'/1en(wz2+ILkv)^NV%x/{P$m-%sIV######/########,]*Su%)H%#?$`[=4-NhIVI,$/6$$TD,IV###.,#x##>5#7c$BQ%>#$>Z$J#%G##bJ0iw/###&##4&,8IV###95#;##"},
{86.03,46.29,6.63,4.97,"4u#q,C~/31,#FT(;ATiH)###TAD;K5######jBT#########Z@-Vu#%TK9H$1M3.v&rF6fZ'ZBT+u#-##J?$hBT######;##t[)tR,Yp,]5$A?'p,%k|.$g3{?T###|5#[q,I@T######z##3##uZ(OH#Vg995#{P$m6&*jC1^3r6)G?&(9,I5C^#&###F##"},
{415.27,48.57,6.90,6.25,"<?#=iYQu%###&Q4$eY###&##|V*$OG$##-u#V##_aI$##yY$uR,Zz1*)@+##/hYnR-OG#?$#d#D|?N:?'E,#w&(r}L95#%##k7/{,$Ry6|#%ZeYK>#^>$bd#dV={[(/g0rw-)^0LR*=Z$a%+^07zc'&l#L$$'eY###$##1:#YI-###6##v['.,#%##R?&M$("},
{117.27,85.14,6.09,0.42,"78I-i@;5#8A%^%Tfl&46%r13?g8###,##e|0Z$*.,####_H#+2=*##*C1^b46'T%##v6)Sl#m5M######N##AA,######Vc#]-)%`2Mr>Rf'F$T?##WI,;B#k}KZP####GT%VH'R-)###*s.`L89x/D>#L;';)A&##.,#e_&Bi?rZ(###fn$'$&>(:$##tT,"},
{60.84,98.00,6.12,2.97,"<,#Zg8###$##.6#KD?######'o(iv+###-##6{?######C$#(##VXAD#$qP$HI%wzTII*L,$`{T.W@SG#+Q#xvT######7$#:Z%z%.,7%g04A6'+f+ys/r(:<wTa,%Ol#VV-=vT######)$#bu$gxThP#]5$'^(gvT}c$hl&NwT]#&F5#&.%6vT######B$#"},
{187.70,121.65,6.20,1.59,"-S/###(##pZ&2dU$##d6*{5#tn1(##feU0$'SH(&##>zMGfUoJ3$##5,#eZ%MdUW,$&@+06#[EER?$HdUvl%Kp8Iw&?f.JaCT]52,#yb#N5#UeUo8-g,&5##QlLC{2r@/-##J`B2d$0[U2##{./R#%.,#$##Y<C4d$d%0$##,eUv6&bR-*##t;BTu#%-R3,#"},
{253.62,123.66,7.05,1.84,"m04eY#$##s#$I^.[M=P>#Zo/1,#SM25w,[S-)##5?%w#&xJ1+~,95#/##Su#?,LfQ(=Z#{s5DZ&[R&ka<P,G0]&c~K<EB|[+kd,###4##dQ%9KU###8H%7m$]L89##wLU=Q%_D87.(qLUY,$Au$###9,#16%-KU###E>#1l#'PE###57LgY#/R*%##LLL$#F"},
{319.96,131.84,6.36,1.69,"8@+######,d`p6*%c#wZ%;c`'##a>$./'q]`Fl%###tG#u``63DF>#.,#KN1ih;r)2lm+J`8Uq5NK0`.*aB7qb####QI(.J/}z?/,#QG#O[$5C95m$LW:.'1*bIX,$'x,lh0bG$&##|e.pG$3]2######jH$-D7*@+6Q#$M7,Z$P'1w@-<;8###'v#%/0D,$"},
{319.96,131.84,6.36,5.66,"TH(jc#?v'cL/-##$J(x(/o%/,##f>$B)0ZP####b##Pp495#.13nw'/94cQ#mR&MX7:GAn5%F*<^6)1_)F%+{k#Q&#o;]|G%ZI$lm+;-'###X$&0w,OD,98/w{A33(}e-q'0(Z$f`&A9]-c$/##[Z&95####V##kB-%?&sb#M5#XP08x26u#s~#?7:^$*eY#"},
{127.53,145.07,6.56,4.87,",Q#+R)######}c%[?)'##tm(dG#V#%mw*t[,.,#$##:M+_I-{-)Td($6&2,#z&~Gc%o6&`473S/'##/;NN*?$f2'##$V1'x,.I*XG#qv(|I,|'~F>#fw,l].R>M&##u)~o#%7EC'##o<?&$$ZV=/,#^G#,f-w&~F5#k?(<V52GJ$##B(~/l#hL<$##f2<Zm%"},
{37.97,147.90,6.29,3.26,".,####9##j?(4x]###7##ER'@ri######3##-@+###D>####Z5$ZP#)##2['U/^###5##e$(:ri###.,#5##P~0###-u#%##h-(:Z%E>#2,#i.XOG#'##O-$(ri######d##)f2###,u#$##]?$z[S.,####CoBIsF###)##gtiD>####(##XJ1$##E,$&##"},
{201.14,148.85,6.32,3.24,"5?$%]2A,#T6(Ce*j>%F5#1,G[P#5,#cZ%b$T'##3c$UG#@NAFT)q$TB5#&=CO)TuR.[#$J)2RJ08m(Ip+IWB4,#0%,U#$3[SGu$,c$8,#P&T->Iuc(%##?;R57*W3ESc#Dz;###?YH+R&A~/$##'{>###XrA0,#VPM###+`?$##hYLgu$>'7###{]3DA*@H'"},
{203.77,164.05,6.12,2.97,")##n,&j#$Kl%<5#1,#y,$nd,95#$##X5#-n-###'##XG#C6(cK%svU*Z#^l&qBLT7.],#dXB>m)*##H0({vU$##3u#Zv#DGN(8(rT71T#^wU?zU$m(9c#2TEz^:#m'Es4Xy3I5#+m'US&tEG0##du&#%#jvUze.:$)w##f.O;7&3M>>R$[~0w##/`AjQ$AB6"},
{530.72,202.54,5.49,6.26,"/c$XG#OG#;5#t$W######+##]{m######B##fA3######&##e#$m#&.,####p^bD>####'##f{m######D##of6######'##<?&*l#######z2m######.##M{m######Z##<]4######(##f5$}Y$######H'bD>####5##K{m######O##}I/######%##"},
{181.08,298.10,6.43,4.69,"95#m##Qp2u?(TG#vY#op,@'NfP#|k#i@'W'NnG$PG#~?&m%)g(<D?#QjDV.&nq7Hj5yK7Wl;O@Esp6h`4T^23c$FZ$4><C,$P(138*?<@2,#+=2~i;9W=_A,}$NXt<L@*6L*^u%OL,wU8@#$=##+.&:z4.,#aM+C:6t@+YG#}:0%_7###A?$$6&($&$##eu#"},
{103.43,300.44,6.17,3.60,"dm+8A-kP#f~L#(4JM>3,#/w(uR$i%0###+##5?$]d*D>#'##s,%n(PIZ&fl%Vo).]U|Q(WG#}2UwL<>#$BQ#Me,+['uQ&@7*N'2zp+$V5iQ(o|Ag6(4C3I$&I~UQ5$CQ%3B%0S.:5#G,#|w,+6&&##gD'qjH{d+RG#W6%1_U<6Q.,#]P#`b40%,95####P##"},
{228.73,304.37,6.28,3.80,")19Kd()##0i+i[*[<@vG%6c#a^7mC.Zf5/##LZL9'/B[(5~)qu&T7A/6'F[%,7'B'RiQ%Hu$8&ROW?lc'q[$S&R(I*--%Ay/mf0~h*><@ol&jD:1R(Hp'X7+=%RN5$`?$R0)B'2ac'B##w%,.m'5,#;*.d4Kvn,ZP#n6$?mDAp6F,$lP#ZE-Vm)eY#$##5m$"},
{179.41,360.58,6.41,1.78,"?d(eY#nu#M~-1Z%M>#`C':ZMOS.G>#XR#v+H<e+U,%f>#<d%dZM-u#Xl$>n%?_:<z+^4A{{:;bA28,?a50XCF<@@5#|{9?@&2~MI,#>o/7Z#|}<jw*4u:D18GV<zE12(40N7^8-TJ&qg9.##K'2X?&1Z%I##WB+z$CzP$P5$2d#'^M0c$7#$r9)K~-###)##"},
{179.41,360.58,6.41,3.85,"Sc$8uCWv*j-%eY#JJ'#e+7i5%##,v$Us4&z:U>#Vc#Wt8|^<(2;ge,/H&~-#|L8#k8d.,>N0wFAty42V(L3@@u$36#Ea2~{BI$IhJ/X>$I,#97Ly//TT-.a?Kh;;b7O~.2X4oP$7;*_K5Y#%y5L`5%'##XH##D5qg7L#%EZ$4$(6mB_Q(]u#q6*bA+],$Z:0"},
{101.32,362.05,6.35,3.63,"fl$l@@%e-Xn*>##0[?Wd+$##_6#n%/~P#:5#MQ$ym*OG#%##b>#)<U0Z%###en(t9ULH'&###hM1h:6c$4?$B%+p-(Bv&~w+c['aq,Vh69I+rq9K@*[&03/,e7UgG$7Z$cB)w$+95#L,#eA.L,$(##h))->L/v'PH&(m#S;U~2?fo2~G#QF<|R'TA2###-##"},
{227.00,367.05,6.25,3.76,"_6)b58>~/{c#+##GZE2R*###2L5>L2MQ';##ttHtS*N[)We'E5#=:N1Z%###z-'`ARYu$<5#,BRC`@G>#4$#^ARmR,:[&uS.#A)w_/F2<Z,%CM7Nm(57%N['+@R8#$W5#'g'Z);4J.(c$},#hu%2,#8E-xiD[x1lY#Y##0eBxA3wA/3,#0a.N6'x2?bG$lP#"},
{82.38,368.83,6.08,2.78,"0i%l[X###%##I2%X~XD>####gc$T14sd+95#.,#%##O[&vR,@b0Rh7-6'@[$UaXgR-###;e#h_9TH(>5#]G#rG$ZP#5,#D6&H[*#$$Ew'dUK_|EDQ%]@)bw?-;5m06B[(h$'I'46#$$##vu#^S.gY#tc#FG63Z%3,#WS(l&0F6(VG#07&Ws>/d)######3S&"},
{384.66,369.89,6.87,3.21,"'[$47MH6(bP#kv=s99###1##}V5hn+eY#6##-<B<l$###j##)7*i-(<g79##m8T3l$$##F9#.@R0c$0H$'e$|,P###$##<##6%*A6$QL:mY#b;TQG#.,#Z?#rdQ=5#,Z$yS$I^8(##vb##%%_>$P>#eo*bvSy6T.,#%##ZI?-XFc,$6u#rh&E$)RH$i5%=%&"},
{264.02,390.33,5.88,0.97,",.#@wU95####XM(Xo495#$##S6';5#95#X>#)l#_,%D>#&##S9+V+GE](kH'}{U)R*/##r$%rV<`G#fY#I5#$##6$%nP$$##E@*t5&R@$kh.`vU###0##<X-R'7D5#T6%.5=?##I%*Eh2Im)IZ&.,#+##0&C$3D.##+Q%cG5pQ(j(.bL7@/J6I(_//&d;]6("},
{65.31,419.05,6.44,3.74,"######xK'#_<Z(=.,#z$$*BOiFAo#'###Ie$$;QmP$###K,#^u#e5%Mw#]:<m-'2d&md)sy6=FBFu$J>#qI%A|9######3$$;o.m..92'l6Sd,%[G#E*-K#Npp9/,#d>#dC1gI,.,####S,##tE@D+[b?z>9G>#[1'r'R/@+(Q$=H&~5$pZ'[6'~Q'###J5#"},
{65.31,419.05,6.44,5.55,"7##zE-uG%###I##P0,6#$###F,#~@,OG####J>#;Z$G.*W>$.$#%j5Bv)###=v#*>B*Q%>5#N_3c'7E>#uG#uQ(O>#D@).8-+34-8U_5%0##1(3rR)k6(F[)E7U*?%A#$SU'}FI<S/%A#)<U[GJ07,###nK'o3DeG#c5%Xu#F7UIZ$g-)_1$qnMi[M|W:?1-"},
{114.62,421.05,6.42,4.66,"xY$aP#/,#(~(xP$'H%WG#s>@7c$E>#a[&vIL95#wG#98,5$(`kKC5#n>%{n',s>0:-@2>AQ>9&K)A,#4;FD6oP$Q-#W8NgG$V>C>Q&(m&I>#T,:uV=.*<2K.A?KR,<*x-S_,+c$Aq&R6N3,#-Q#3Q%,.)3?&VT%(q;)v&Jz4`.&oE>b#$sT2###,[$*y/M#%"},
{177.19,423.04,6.38,4.58,"zQ(|Y#AZ%Y?&DH$ll&6%+|dLIA'w%0Yv%2eL_G#.%)d7*,I*fbL8u#;#$Je%oN>KM0rB7/IC78HYJ.]|82X>2v(c,#feLvl%Nj9??'Od(rG#A+5_^5R;:8K.:dLqG<28-m^+KL8Xq((4EQ,#l,$|G%/%(1.*+U&q.0@?$TkDRR%I28oQ%3T1SH%?7&}[+pb#"},
{239.94,425.65,6.39,4.29,"el&RZ%8#$m5##m'^,%D6&|=B_n'D,$$w%s5O*S&n/2pH%Z~/xh=c,$d?)b7#Bs?;=7RA1-}5_~Iz]2pz8QuK%Z$|v$z6Od$+-i/B?&<x/@,#e45{n*]:7Qp0d6NzF7[o3C_,}Y$n_)B5O%##$H#E>#2S'FZ&.N.kQ(|>#rI+gv%7N4{6*}Y#&##@9&H5L.,#"},
{48.04,555.50,6.10,4.39,"QG#XK&#5L###az1xA07M:###K*W'/1|k####[<+)%W.,####$##ym$3`>B,$HB4h,#&.S)##3(WS#%#?&3##SC4hV7{k#######=Z#BE?W@,Uo4^$#t$W^?#r%W,?&(c$E.#+S+^x1###*6#;,#-Q#xQRS#%Ro4P##ANB_x$M'8.,####d{,Bu$95#$##ux)"},
{504.97,557.53,6.61,1.65,"###S:Y95####Nd$t8Y###$##vl#D9Y.,####~,#<9Y######Nl%zBM=WB%##u.SMwT+Q%F-#)x-eY@Ih;h5#OR*pt@;@+?##6?'D,#R8YpP#^8Y1,#|c'lZ#ux1D5#a4A<,#;J0,##k15X>#L5$0##oNCe,$qL<###7c$)w#ao3###WS-U5#?6(###HV6N#$"},
{474.63,570.94,6.52,4.44,"L>#_Q$g~*]7.Br;--&V#$`q3y-+:##5~+Ew)###n>#&@<C/2zu&Ay%vh@###nnU7%(DH'_##I07'.#e^:M,#TG#?2+lH(7#$r@.<K#3=I&##>nU?Z#lH)M,#MK5/U#0mS8,#Vn.GE*+%-E##P?(J.#O]-5<=amU,Q#dl&ID-en0*z#ZmUa>#gOIlS&P]5R-#"},
{120.94,583.13,6.57,4.88,"QS-K6(=6&fu#Yx+>=H[c%AZ#r$%hR--h+Sm*{x4S6&`[&[x09H$;w*(=DpY#h9*xAScbJ?u#VDSbK6Oc$`v&&~-Ed$J%JC7+gG#fm)+~B(g7f$*E-&R#B)BS$@S^P#AZ$V#8t@/###YB)&A)######F_,Sz=######m##x@S{k####)##}<6bG$###'##zu$"},
{270.88,589.92,5.97,4.63,"$##_g)$y5###z$#_AS_5%###b'%3y6#######Q#)?%.,#######Gh'@6S###&L*J[D]?SCu#VDS(x0Bu$NQ#,y2.m#+%-&#####>~(tO=OR,uH)RZ#OyRzJQs?SH>#$Q$mt4$K3R,#c::M########Ff&_@.######S6#<IP3l$###=5#EV1O5$RG#66&@,#"},
{75.30,605.43,6.34,5.24,"tP8x$,ol$L?'dm#ND9(uAD>#k>%q,$|:2J~,=Q%###.i)k$*@[$fI-mF:m>%E7&$nP(@*2,#*mPTS-+6&%~&o$+1-'a:%le-###=5#7pKmP$W?(.$'x1-{n0([Im1:TZ%t#%T@'HpPJd(1u#0,#=Z#?k8;c%x>%X,$@*.gy9<Q$]n,&#<Rl%#6$y`8H7HEQ&"},
{22.28,664.92,6.77,2.78,"f$'N$(5d%u95rV6K[*Y5$ax$w&5?l$,l#*:$`'9E>####i$#_m&:E1TwPx~0GzPDg3G-'Df(^cMo'1###g##kD@^,%###8##TQ$k::{N9n&2<>IkbH#m$2|39vPIE9###.'#[S1PG####Z##T,#Ii=4R%.~.ZOFImNSZ$P;.duPG6(###<:#ZP#######;##"},
{356.85,42.94,6.86,3.65,"=X@(c#mY#ll#lUY.,####`##[e.1c$=5#o%(P7-PG####Qw$y|B0u#o#&{u#=TY######A%#.::######ng&Fl%######qo'wC5.u#]]-%H$bVY######G@#mQS######'Y2#########2>90H&###@(+-JMvRY######6m8sVA/,####d2(###1,####5['"},
{406.95,87.78,7.24,0.34,"A_6qc%Ei=tb#`,^=u#>l$`P#eL5ul%s-(990$m&=c$cH%o@)b~/>@%dNA~5$f(^sb#'##>Q#o&^######Ym$}>&###Z,#v?&p[+gm%?r<[@'r'^rb####QT%T&^95####W?#kl%mP$4,#3?%bG$###G)1,=51$QPG####{<)xwU/u#bP#]Q$xZ&hu%y?)-?%"},
{498.37,94.15,7.44,1.52,"l,#@0RCZ&###X9'C-R######<x#5-R######mI#n,R######4I'WO1{XK(##|1RG]0Nl%K,#KL1ukA+:4KR)+?@ox4fP#>@&}c(.[#K/RTo-5-RiP#m>$,D(nB87,#gn(fW;4J/:5#F5#bL.hY#?u#A(.3:5xd-###)##jW1an0######Vv#4l$/,#7,#jG#"},
{338.25,125.36,7.71,3.03,":7*dJ-ME=>g5*h]=6(*##5M)M+J$l####En%wb#OG#.##Z$)}?)VR*tJ,mNAqf]L#%@#$BY.T;>i%,SA/J%&HC3-[)i5$8@)N|=qoUdG$R>#ek]Xd+~P#_,#;{1VB6[@*Pc%vB11e-=u#`6&@uN?7+hG$']&xf]pb####$%#Z:1$m(###$##Z7(&m(###'##"},
{156.53,127.94,7.27,1.50,"=V>######SQ#,%S5/-KQ'|##QV=4&*F&2:R(*6'L'&Wh;p83RWA4l$###?,#ORPBq2Yl&-##-%SYS+[U;E##ac'+A$}#S&##0|<A$).,####G$Pw[$.7,7##9$S8.$h<Gt##4c$x44w#S###lV.}6+######A+?X%*bG$$##B/R=)4()@*##]5$ux'($SlY#"},
{272.15,135.34,7.30,3.29,"95####H,#o6)Ie/%##Su#em'Yd+4C'MZI#/-B-%h9R%]-f<D%]+Gv)t>#?q7Kf2$##f1.^K5A;>P,#Fj5NmL~>$Ac$HS'gTVv7JlA2p#%k7)C29m^.]jB1Z$`BKaD?kc'YA);I&vRV###[{=gn'V`C9,#3?%@m$<YBW>$$##*A$ViA95#(##P6$5cO###%##"},
{38.38,149.09,6.60,3.27,".,####1##}Q(;o^###*##Sm&Uze######.##0~.#########S5$.,####W$'S&_######2I&Hze######1##AK5###.,####`H&,?&###5,#P/[{k####,6#?ze######L##/g7###95####CQ#4JT######z8A@tJ###%##S}eD>####%##;95###95####"},
{136.20,302.49,7.68,5.88,"iG#%#6~#&$##r8$'9ROG#%##EO0.S/###'##D6&Bu#sb#9l$qb#hy)@-Ck#&R]-%p-Xe*jd)C:R7?&^P#@w%2_;iQ'P>#>a=S-#+f/Ji*OR,9$'<c$rK&g7RV6RUG#hG#g%EKPGLy0//-/R'^:)k1::u#sP$C[H*6';##z]2/8.{k#S$#j<B{I)#B0Kw&o-+"},
{261.69,308.82,6.89,0.51,"ic&'s8uY#Fl$mJ#d6R######|@#nA395####2Z#26'######d5#Pc9uK8###^;)`7R:f2jP#-;R#'6gY#i5#%].:#$###%##+d$+J&kaD<y5QJ-/6&''3^mB67RPc&lY#)i(sJ03m(###8##OG#+##ER&37ReY#q%&FQ&r8R_@(<8ML#%Re'h0%:NC###$##"},
{413.43,308.69,7.52,4.01,"j/*1%-A$%GR&GDT@H'###-J#s$MA,$.,#`Z#}e0######X,#_5%###=o$<DT.GA'U47(1k,8sDTC&2h>$`m%~2=######3##ac'###+R#tATUl%vG$Uj.$AT/,KF,$--$XBIX]5######H$#Fl%###-##nd'######X##W]3L5$###'##fo-:Q&######W,#"},
{413.43,308.69,7.52,4.94,"X&'2^8C5#H>#%/'DU9###$##L@#.o2######TZ#;c%######j*0P6Q.,####(<,*9Qa$)hu%%^K3XEwb#qZ$a9195####(##c7Q4-'=5#PC-juM*I'y7)h:Q?cL=l$#Z#j@BAK4######_##(J/###,##o0MAr>######G{+UQ'######XJ*26'######6##"},
{134.70,363.95,7.62,5.91,"wG#`O5u5&$##w&$@YD2[${,'Q'.P6)<$#nXF1x0.,#[##a2;>5#y9*:k:^H(,h.:B1fI(,e+E:S2H&G>#w%'hy8Y?')$#VuH$##.,#&2&A6S_l&&##T1&p7Sr6SI>#WG#)[<[tF7y.u.+0%*~R%[/1~9-d1<@l$YZ$V<8?U7u%-{k#O,#*<>FJ)jo2FR%-.,"},
{259.16,369.88,7.92,0.38,"H,$#I&4?%a.R'U0?D;sb#F1+OS&-S.95#/##+##k#$>Q&95#&##{|0:;?Qm*@:)q0Rq7.?c$|1RJ82^P#9Q#67*>#$~P#mP#J>##A$xvLhC<o$*h#%)&,c0RZ-R~P#gP#at4Hd*S#%/,#%H#`5%%##8%'SuJ4l$%6$KQ$|1RtP$ov'&L/`RO4Q$Qd(b[+%l#"},
{523.65,367.96,7.56,0.25,"06$f6*######e<]W>$######^;]#########.,##########28+h$+###4##_:]OG####j##{;]######C##;?'#########RcIbG$E>#h##u~Z######?0#v9]######O%#j7.######*##UdM6,#C,$+##.c:5?'###;##qET;c%###1##A(0#########"},
{426.77,415.64,6.82,0.64,"*H8###:l#&l#C`T.,####'##k90f#$mv*{P$#p3?,#|5&9d$Zw,###Sw'pi8WBU######KM&d@U&##-l#Y~(?d*###?5#To(ZP####.g%WAUdHFt#'T8*@`3$CUOG####f$$e$+######_?$.,####uf$8D=_$*QG#'M'znRD)@jG$&##Pr-=d(B,$###7##"},
{522.73,420.02,8.48,6.05,"+OG?Z$3l$o7$KGK######g]'fwZ######~3'@S0######1%#`A2;5#0,#2M/vqZD>#/,##H#`r_######,%#O[+######<##L,$pl&<5#@$'mo_x5%RG#;R#eo_######''#g,&######(##_P#q?)0,#95#Po_S#%###J%#Eo_######|%#############"},
{409.86,429.17,8.06,0.83,"M'*L5$,/'/.)$NT)c$.,#}v#g>HH>#+Z$7w&KS0###$##c$$)c$###8i)6KTT.NCH')x&zoI$NT<c%sY#C%&=A/D>####c##@H'###X-#(KT1?&fY#ah%?JT8*Ceu$G?%kSD~n-n#&.,#T##o#'###D##Wn)1Z%$##xl%?80qb#(Q#d~-u.-zu'oY#SG#,7%"},
{409.86,429.17,8.06,2.42,"Fm#MS0######_##58/OG####_,#fw-{,&.,#-7%zu'pY#_P#/n&aYFG>#{P$5w%$NTRu%hP#`JD.|Beu$RH%t.-pb#%Q#g~-Yd#$NTJ#%.,#<KHM~Odc'rn&FJT2?&qb#KV%bS0=c%###UQ%]d(<L*(c$(J&EKTU,%###ZD)*KTVZ'###V-#gw)5?'###F##"},
{122.09,530.57,8.14,1.05,"(l#bP#(c#/Q%5l#K~-1,#zl%FQ#bw/%##U,$B5#@c%$##F>#EZ$&$$27,###2%&0WTpRT.H%ODR7STMQ&uH%cw*ce/###%##Lv&EQ%s%/###PH&`S'0WT><?+TT0~)M8*+~@qe.vxGa./9##i>$]f/|k####`5$^i,MS-7Z%qD8Se(h@,9$&#v&*&'?[Fiv("},
{122.09,530.57,8.14,4.45,"<D1s6*ux09c$M]3;6';'/T0+B02DZ&'##<0(QG####H$&E$'(mN1Z$9/.5a-LT+*9GX]TQ.*<`TSjC%$&?e&dn.###($$jl&xY$$##g7*^r5]Z&&[%:2RU~T-~T-H%xv%<`TgR-0,#6l#fl#######1##5-'/,#`P#z>#W&3.,#7Q#b,%2x/`,$`[(C,$?5#"},
{122.09,530.57,8.14,6.04,"*m&QG####Bv%%$'en.###jl#el#s[-<5#6l#?5#_,$~R(Y>$io'-y1Oc&'##=e&<`T8FD&$&<`T-~T8Q%xv%H/0.,#)H#m5%O'+H95E?'+0.~7*X]+''HY]TX~T]Z&%[%tqQd/3/,#TG#y>#X#%_C/~$*,p1UW,E[PFl$</.ni5_5%$##~.*A6'######1##"},
{271.07,535.88,7.98,0.17,"###>l#$l####1##PB2OG####Y##>04{k####*##86&D>#######fH#|w0###%.&F<V<8V16&]2SH8V&-&W@&O/.{>%/,#/#####5##J'2mP$},&0e&F<V98V2.Q#R(I@&F<V~@-gY#R>#)d$###*##|u%X>$###/{.VS,qI.4A%09VmG$08.J/(E-)&##(l#"},
{271.07,535.88,7.98,1.47,"###*##ul&.,#,##x~-@?'###2-$n./1u#mY#~,%tG$.,#~G####+6#)B3###cR&&rQ{wU:?&}{Ud%T9-'>7&>8/Oy(q':9,####4##do0fY#:?&-%&}{UuwU^wUju%*S'}{U~I-$d#(zUTv&###/,#Y,$eG$.,####%[#9/1+u####2##^y3.,####f@%jd+"},
{271.07,535.88,7.98,3.61,"-u#vP#I(*^QMHT4rw-@,#,I<+A)#A/###Q##`>#Dc%######)##=$$]R+%R)<7&5EUbHN%.)5EU:AU]l%A.&g]1{k####9##UG#oP#UJ-}b#%-&J.&WDRIAU+AUw#&Ym%5EUFJ0######TH#W>$###,##i5$bG$###W##b&1W>$###,##290,u#$##$##/c#"},
{271.07,535.88,7.98,4.98,"###%?$@u$###7##s]1bG$###U?#ge0######)Z#.u#.,####5o&I/-.~.###qR&}{UvwUhu%}{UswU26&hv%=02ZP####/##PzU|#'x%.2I$Tv(2.&}{U1RO+xUv#&Q7&.iRRK3A,$###T##_J2$##WG#o9)N5$F5#,d$wd,9c$;5#<,#B&.QZ%ZP####+##"},
{35.32,538.61,7.28,4.12,"Jm$.};ie1###T)Q1T3I#%'##b#5|Z)######R$#P6)######[-*W-#}z?1##W9Z1l#du&($#VAAhbM###*##()#p7Z######`$+_$#JT5eI$QwV)?&ZP#ee#6S*[(;5l$2c#R.#XQI4'7$##uG%4##2v'ez)16'.,####;{3)l#~P#,I)Kn*<5#d>#^uNG>#"},
{420.36,540.47,7.86,0.82,"xY$=,#jG$yP#(##]T*cc'F>#3Z#[9,yn.f-*PG#*?#tu%UE<###F?#l~11,#~m%gBJEwS@Q%C{S[QK8$'`e({/3aG#8#$H&+###g>#rA095#Qc%/@%C{S;[O(%Sqd'7~(9VP?//w6K0Z%]#####$##*Z#+Q%######7.$V@-###Fw#g4G*S,Z7*'d?1jA1Q$"},
{420.36,540.47,7.86,1.74,"I,$I,#~;=|Y$w[,|v*,A,0x(1:3=I++Z$[##_5$QS&@'1;#$D>#F##5;;]P#U%'H(OoeS_H'+iS^?N<v'47&<w*Rp&meS?#$PG#(##KK/N5$.6&;[%+iS&eS,eS;Q%Z@&+iSYd+%##Q*2H//bG$%##Il$G>####$##s?%tI.95####3##~y495####0##cR+"},
{420.36,540.47,7.86,3.61,"B,$2##H(+&[P,T3nR+[u#ISCJf.ov+###y##^Q%gY#######eY#aR&*&/tR+;@&y(R}6Q&n*&{R{-QEZ%NR%ao0OG####W>#(m$,i=i@-D5#/d'<~'&{RjwR(wRiu%@d%YLNQ..######$R$p?'vc(&##uV5:H&~P#K,#r_5eG$;5#'###'.E,$|#&###K5#"},
{420.36,540.47,7.86,4.48,"-s1UJ1######8?$:g5######}c#4/1###%##4Q%95####.##_JT*6$o$+DL(lS*$NT<wSXc%$NTO7Shu%C[%@92E>####1##j$+0##Q&/yC._Z&6%&$NTU,J;JTw#&E7&YBIlK6QG#.,#d##hA2y##,w,7,#^#%i9*#S*'d(~q:'m']G#Zn(iH)]G#{k#N,#"},
{87.60,554.64,7.01,1.04,"[,$5H$ar<###:d(6R%b1Q4?IEIQrQ&{#%-l7F@+e->,d).##B5#N_4bc&###1-%6<,c^3>R*:+;zn*L~.+@'Rc%{L)KKQM$(###(i.cn/eY#/K3`2*<R+d##wHQ>6&|,%VR#ic'I?%&Y3u$*OG#B$#fHQ*Q%>V>#[#'/1j$#U`?$I)K>#:v#D#$-0,av(%l#"},
{87.60,554.64,7.01,4.95,"&##7,#z(1xG%k27D>#Xy/?5##vGn,&###j##qw0E>#+##6I$c?)5u#@9+fn+eT6###N#95m'yZP+u#R5#II%iT3o#'&##o-$*T,~OHyQ&p,$8J(sw0tq*dT5TD<k:3vB3x05>*.^~Pi6*?5#+~&^>G(L7S,$#_Ph83&H$pv&AS/-@%R:P>n,KW9U14KG;PRH"},
{440.54,582.31,6.91,5.48,"J-#4fQ+u####tT&@p7+V3X>$]H'e,$~582eQ~'77u#xu%cdA*0+PxLRx-wS0LhQ+/14##:~'{1;6#$96#GL7&##2,#x01Cn-=c%'##en(37I%dQd5%}b#fO2Q;1cq=|>%Au#W%#Gw.|@'cG$.,####rG#BcCfY#H>#4[%)V1},#r%/~q1sb#}&#%q;|c%###"},
{478.10,588.23,7.18,4.38,"###(##X8)`5%bP#k#$PL*LB64o0Y5$jc#GnHZl&4##&B0WB-###G$#7/1###jv)sC'fq>$##(~S$@'?c%z5#v$,RR#in/K5####s]%;x16#$|y9iC$X{B1##*[S@Q#Km*`##8v(r1%'V=)#####qQ#0a-~ZS)~.u?#H_,o~S)-SU6#Co1D21[6)/($-ZP_>#"},
{418.90,592.53,7.88,4.53,"Z,#MD7D>####}w#qvQD>####(~#ND={k####L,#w?*######@##va2@&07#$FC*ZyQ_D:%6%dzQ0C9wP$K-$T%+v*9?d*&#####L##X|03'7X-)fP#Vs0lxQ9-QP>#5R(y,:`-)=q+PyQb7,###3##kS)?d*1~+rQ%v[)tB4e'6dP#eQ'_~)###hl#Ah/f[,"},
{386.56,600.66,8.52,3.23,"8K*|(8Qd*^5#^->Bv)###K##'J(4e.###:H#1$&d5%###;-%z')0%C;h<8,#1FX994###&$#2aBwc(###J8$5u#Z@+0,#1d%-S-S])oeAcCX*BX#v&S>#^-:cj?a$+###N##n5#'n,######YZ#4cA_/(,FC^i3,:9'##|z1g81OG####2##l>$[u%###)##"},
{16.97,627.80,7.79,3.37,"T;0>'3C90)7)658bC4$##Hu##`:1v&###+##############0d(_'/i^.=B2WZG6Q8<6'&%(aD[Im(.,#;H#W>$#########7m(8W;wP$8w'?:9:S'Um(ayKBC[|Y$/,#;*)C7,######*##W94o6*O>#}$&Rw)0RO###&/*-G[Bq;###:Z#/D4#########"},
{260.22,692.44,7.63,4.83,"D>#%##6H%QQ'dl&J>#g>$B?&o%0/,#'##g'*eY####.n$dxK-c$'l#Qu$=I*780JZ#W'7}5&h~YQ>#W>$^;%]uKS?&C<6gMP<Q$,v'_5%%##{e1^#$GC9^G#raYvc($Z$)6#%gF%`Yse0hQ$x>$/v'D>####[#%AR)MH'SG#tU-b~1^P#VG#R7#,]Upb####"},
{426.75,714.69,7.61,3.34,"O0#}@X######~K#oCXx5&###cL.SDXK6(k5$cz2fWC###1##`w(ZU3jR,uZ'FCRn8.[R)~J-D`=xS,h/-?`URsB{k####MS%Qc&_#$A%*D3<2AX###P##c12_i@###&##Ny.mQ&C,$###o,$~R&%Q%>##7h7fS/Yc&4##U{3yd,.r=###Wm#/Q#2;=95#mY#"},
{79.93,76.36,8.31,1.27,"<,#$E*b98###s',z93Pl%wY#V1(SH(###?5#GT,L5$######l?'dM(K?R5,#uARC%*Ql%:?$5i696&o,%9H%n3<OG####'##Jl%J##xRA1*Bi,P;6%.w&xb?Q);['0Mm'od).$Ppb####`H#######4(#.?RD>#(##j9$Q?RsZ(R>#)/%z@RB%.###=##wV4"},
{261.21,333.28,8.14,0.95,"X##h.Q{k####n:(v<G######h6%QZ%3l$######O,#e..D>#2v%$t:Uf3:H$_1Qc/3yP$sI$:V5{3EZP#*##%-#zf5b#&.,#gG$xG$C91gH?ff5sn&>K4*a.H.(_1Q`&4`G#BQ:njHRG#X>#du&###Od&k*7.,#E$#|.Qo7+qG$fI%%5HSx,*bF}Y$dP#)r("},
{261.21,333.28,8.14,5.46,"+##d93eM7+u#ZG#+A,`KQRG##1.W96N%*j5$;KQT,%###.6#kA#.IQUQ'###=M+0f2y}7GR+V7-(H%(V(4KQ~PNOG####@a/_&+r818:0O7*eKQWH(A,#gB)?L9.,#D##}=Cu,&.,####Z&+?u$(##{0(eIQp':.,#'##+>5F%,OQ'###O##cc#|H*###&##"},
{476.96,513.05,8.05,4.79,"kc'4c#h).e>Hu>Lc>$P>#B(.v./F5#-u#)##8#$K>#6l#X,%f5$3K&,W@5l$]JT>%)X,%@c#Eg8&##A$']u#_c&0,#U5$=5#HI)Cw$|HT$##wJTfQ%Ad*0H#})D$##f8+6w(Xl$lY#=I)/,#`6)1##,KT]bI#IT)##-7*sUN{B9###Ri)JBIjY#$##2^)V,N"},
{121.21,559.58,10.66,1.11,"uY#f&.0Q%E6&g%%Rr?D>#hG#p#$'[)###<5#z#%r5&###)##r,$[`1J%NoQ(W'N8[ImL4m(1R11a&Ndu&2##jV*;iA###'##'c#@U+VV5%v'ii<Ay*4D5LU0L].zb:3JITx-W'Nj./`G#3x&wG%.9'Bz:vP$@>K[.'c?(DQ#|d,#%&GL-hd)TiB)l#W,$,&&"},
{121.21,559.58,10.66,5.15,"9?&PQ&i9*&m'(tA9#$SJ$,@)94@|97H>#,$%IR%/d)-##CI*lA.}_<fR(%?%NE6z_=?9*io2qE8.%IhRCZ|>oxLNJ0p>$eL-,c#&8,@-?D^74wL`e.(J*u#=|B8%##=y';xLI[+&##3H%v`4######7^%/FET,%###q>#PIE>u$###4,#vd'OG#*##=#$ZG#"},
{416.99,565.93,9.66,1.88,"'_;pb#bS,+7&+C6~#%zS)Au#.[(tC17^0q#&][%5ZM.,#%##h-'eE9LaAol%mc<tuJ^J0q,$).(M|/EKO'~-f`:qp65d&oD3l>$E@'PTJwU:edM/.*J]*ILO(%,%##@E-$mLVH(###fZ#0LO95####vQ$t?*eY####M,#by3?u$###%##:w)*[)###-##zH'"},
{416.99,565.93,9.66,4.63,")R#(sApb####($#:(8L5$###3##y@.OG####Cu#r,&95####(&'~0P|S/HZ%~9LdNBYQ';m$$9*B1PR'7pY#zP98h;hY#R>#*$'Ru#`1,l$N%YG}Y#@&-]b6?/.(U,B1Pys@~.PhH(JQ$N33ZH&KH$j&1T@,/<<w5$b$*%[%CH&#x%i'3c$*Zq;^c%sY#q$'"},
{38.90,603.23,8.47,1.77,"}f*=-'qu%V6)`,$_c$By3[R,7}A<5#zG$(m'u{=YI)$?&&##P`2U#%:$'/,#uM;{P$mH%%'1@B1@l$jg/^J0l^/xw*}-*YZ&N&,###gU%1L7#oT###o,$6H80w+{P$(-9+oTWy1Eu$>@&]oTfY####i&#,nTjZ(###7##TrTeY####q$#LoT######R##PZO"},
{38.90,603.23,8.47,4.26,"*%%dX?87,[G#*G2#=EtZ(@,#)X,9`B###%##s##.^8######;$';R&(i<n804JN^-'{S0[e%l=6>{U+%-C##@N)jwUg,&###Cu#t]$O;@zY$H`4JA+]l&wG#46%a(*1}I?5#A7IA['ho5&##f>$<m#0)@M5#@J+Z#%o#&bw&qY#nu#GtHZQ&Rz<E##7J0@d#"},
{233.52,631.64,8.29,1.15,"95#b>#1Y9|I.PA1.l#aR&Sm;4PL###$##4j+e6*7,#FU+oiV######t?#@&2{99###a##]b<AhV981$##A;%7h4{0M#i:[8+>#$mP$%##/l#Zy795#&##/$$K7Bp_@###6##OH#zgVbG$###*v'bG$###F##{6+vY#eY#P-#W6&R[*{k#5##=l#;R)OG#%##"},
{64.70,659.15,9.19,5.04,"$##p)2k/4###~c%LG8T+GjP#?y.qi=KD/W7,C(.S94$'3V5$mH$5~D1n-95#`94Li-%_7o}9#~,v>%WwEZV1_v&*'4[`;iY#Y^/P^36-&_d'HvM1Q$du$69*`824##qwMq?&2c$}b#Od@[>$*z33?&xv)e@'>r@###yP#uC'/v(-##Q+;m@(###(##*}3]~."},
{64.70,659.15,9.19,5.94,"_G#}M+H-O8e-LR&2W;;[&lD<oG#N&*}e($^7ic#`.*^m'bo4(c#/%&b0Ov7/U.O5e*$^0sN=1w+W>$>L#f,O(c$###G~#k,OG>#u@$}.O:Z%0=<2932N;o7)HY6tn1;-#6@*s.-###K@#Z=H:,#na7th=QG#5&+3M5{%/@-&OZAFR+$##A6%PuFOG#+##cH("},
{426.77,676.30,8.61,1.62,"vQ$KH'gv)-u#`>$RG#Go,V*;{k#S##1%Pao-=?%(~+,r;$[(|2+i{C&##YG#ZTQZ[,hZ$RA+Dv)9##zu:aXGrG$^#$xE4m?RAj=P6)'##et5O@Rj>%2##sy*`U;###|@$L'6<u#BQ${fJM(<}>&###/##*I?XM@###$##{N.?n.###8,#Xv$95####4@$9g4"},
{242.60,125.88,12.07,1.79,"@x2###GH#BN4sEBx5%Pu$SM0mG$1-%;:5j@+###Dl#'J,pb#`I-###vI'.7)UnVN,$/8,VM-W~+~|7_pVDA-#y,$^3[80w>$%[)zb#4I%WH&.oV.,##x-V$%k;>%$%3sVg-PpA0)Q$Q+<XeOj7,W@)?%+u,%joVk>%K>#h>#138Yw./T-pS2yY#|G#5oVrI."},
{242.60,125.88,12.07,4.72,"^JWTe/###y##tX<%V<5x,Yd'?c$06#sJWeu&NH&o$&[U8O6'+h4ri;/15[-'{NW+lH>:8hn'1z8M/%:IW'##7v%wZ%[I,v#&Mo1H>#jw(_T1:KW`A-O~+S|4~%,ls-TIWlY#&f,hZ%5?'%##fn-L5$###P,#2]/E:5B$(iH%E5#Si0jEDR5$do)7^3%Q%###"},
{446.70,376.25,11.18,0.11,"dw(0+2dh?###-v9F^8pb####40)&x1######_f/(c$###%##T.+J:,$~U?5#d_Ued+[u%_5#2nMX>$.,#2##1OC95####t##SH(0##p_Uae-8[U###dZ%Ui+7,MN,$:5#ev%A#I######e,#6#$###J*,Nx0Cp8###;-$_X/;$).,####k'-~m)SG####LH#"},
{385.35,414.58,9.87,2.49,"v##z3B######8&*-~R###&##kg3mH)###~c#:6'######?Z#bd&^NVL5$%##MLVW4F###2$#^JVA,$###Z6##m&eY####Y5#bB48r,7l$:^#nJV`5%###ez#OeO95####0.#xd+95####|?$?%.:,#Lu##e>B7-######h<,g%0######D_#T,%###8,#?n'"},
{23.53,484.93,10.79,3.42,".$(3c#K5#I)57*]@u$###>S#ZmF95####*##############;#$.6#Jw*(8/7n@mh2b?)1c#C,]4H&###H##(c$######$#####A$#S`?8#$N}ILw$yB6q;)4']^P####|1#;w-######0########ZK(u[-`:<###=v$&W08)]######K$#Tp5######$##"},
{514.83,538.07,11.61,6.18,"a>#Bd*0,#?#$.w)dQ(%##uw)LB`OG####4r%lw0######V##Q$&iD?7d%$v&pP8=7WVl#j%,{G`H?(###SJ$H7-######9##/Z$zw-PI$<f1O4E2y4,[#1cF*C`{k#$##s_',d)######.##_l%ve-JH$p[,;~VcJ0.H#|B0LB`D>####A8#mP$######&##"},
{270.36,564.18,10.26,1.78,";5#Fl$i5$7#$QQ&;?&QG#2,#p,$W_56#$###Uz1JT4###T##Sl#dk9WM>SG#X,9hYFNw-<c#,J*fG87SPH$(`TP~y5.?%s7''Q$+w&D^MkC:6dL~v(Eo*'9Ilm+:,#_44qvNs>N%##yu#2BI######^H#rI.W>$###?##TU3t5&###2##~.+^u%fY#%##.w)"},
{270.36,564.18,10.26,4.58,"6##Xq4=_:###f6#~(;pb####7##^7+A,$###xG#HH'/,#$##^6%DH;hQP1H%<gLVXBY-)FR%2B+rUP$:8&c#k,:kD?$l#xP#^>$Z?#t@E^sCEM=sP#7/,7wA4n+jJ)D^MN}B|-P?6'SZ$lk6$##&##07'~v*###0##9?&CI*'l#6,#G6$t7.8$(###)##eI)"},
{87.65,582.96,10.32,3.56,"###;W.Cn.###'l#N6>C,$(l#|&3`g2q>%,R%4IG1g5gY#He%g,$^t6rNE%##@a?eh4S5$}R)|~0#y)V$)FS,Ny/k5>{k#l##IZ%rv&pKM&i=/ZKdI(<~*Hc>>-(DL$Y]2PU7~S1Pm#3u#(+3###g##5|5re1&K.{+8[<?pB3Or=/d(]H&JPB7Z%9,#jZ%A<?"},
{515.23,590.00,10.95,0.03,"T5$X~.PI%C&2XX?fp8CZ#=XA0mpxY$###`@%i6*######%##;l#|^2&-$)S.[cH}r=4c#A17Vlp3l$###KI$a@.######&##BR*2d&IQ$WU7Z}JB?&KQ#v5IKjp.,#$##SS%cn0######'##ql'Tl%^P#`m&zh>8-(>##G@*+jp###$##@Q#9~/######$##"},
{416.92,594.37,9.78,3.61,"_G#A',Ny8$##%p#=}DVZ'###es0+A0###(##*J-Gc%###;#####j6$;$Q.,#oJ(U>BAINqQ'S(Q%f1B#$od%Wz31d)###*##J?'q//YK14@'LZ&0Z#t}5L%Qe?NNB2/-%t-?')3TH(###9##Fm)oR,$##}f$/u#IK,VU7H8.`9.R-@vsBKx.U@HN[+0,#27%"},
{416.92,594.37,9.78,5.25,"2##jR+5?&.,#(##*(1Hv)###6##ZM1L[+###Ln&y7J0S.cP#'##k|0We0###OR%/hPN/26u#BR?L?MIy2nc%sJ/HJ*3u:5HJ###^K#IkGg,&Z?'(T(UcEB?J-eP)-'(Z#Pb5;T/a5%`J+kf3$##bG#y'.K97###$##S6$*dPdH&_c&tA-#14KT$.e-A7+$##"},
{264.15,652.84,10.27,1.77,"h$'ql'$##Q5#ev*###)##}M095####~v#3*@{k####b#$~H'$8MA^6D,$S,#OAWM,$5?$o$<|Z(G6$D`S0mK+$(-##N{;|b#v@W,Z$.$%3E4oEWpU8mG$?R$C8)I/CZdS>c$w>%Z6$Eq<$##x#'###/##&'0tq0|m,###kG#^H$:`9bG$###VG#qu$)Q%.,#"},
{110.90,672.37,11.50,5.65,"cJ0qQ$>=E;w(;S-+l#)V$,M9G,#FH%fM'UEES.,/.(.g*#.KC^0cf-a/3m#$pxROG#xd#-S+WZ'###q]#2JTOG####@$#BKT1L2He-1o%Ip5,KTHZ&-Z#r6&I'8###T?#sB5S#%###ac#Qw,Lu$%l#}C$uHTu-Tsb#-S#F7F)8095#@,#XQ%3u#Q5$w>#^Z'"},
{91.66,332.35,12.42,3.72,"H?'&34+J*]m)&w&f:6>R(Y#%J$A(c$###,##(TH#########j#%FR@5r62d(F8OJ{;tY#zI'a5F~#&###{##r&J95####%##fY#(v#0@@-'5^(6hE<$90rK0&&F/g7###K##vX;ZP####F##^Q'6[&<|8wu&E6'f@(oY8Ii=06O6Z%9,#@W/*o0V#%gY#f,#"},
{72.82,617.16,12.39,1.47,"(^.CV3fA/O6'8482,D28,|R+CzOo#'###M5#LU06[(###$##gf-#~(ho0fZ'U(:vY#Ag-a94qIR-l#'l#<6$z95%E1r5&0##9D:(Z#sz0be,*q9WG#OA)y(5PvNFS.hZ%nm'H6'=(0av):x-Yl&###S1$DJROQ'$##JV&BKR[Q&Du$'K%fJQsb#qP$C##aID"},
{370.79,159.72,14.14,3.06,"E>#WG####m#%c5$Bo0###y,$8$$r7/$Q$.Z$b>$jY#J#$bH'3,#VS)/v(1,#+.&$;NIM?,##[sWbOFP5$t>#3B3###*##GI)*##Hx&*L9###;]0^`)rwXIQ#*zX.v&ic&HA#8XBpb#&##J-&nG#ng%Rh>:5#sL7]j,FjG%-#F{XBH&D>#`$#K23r5&###7##"},
{122.61,295.98,17.30,3.96,"Gv(KE-HB2_R*5+8'E7$d'Pu$qG7r5&######UV&{k#######D$&<S'>v;ax3wM6JC4Q#<0(7nVSdx3f>$<d%%G3n>%.,#4,#z;2lq<h%*l#&tI+|Z&)H8$lJm{AN,$-d$cR@9x.rP$WG#RR&v|?gn/RZ%6S)716=I'pZ&J?&8.,xY##c#W7'9Q&###=5#19-"},
{122.61,295.98,17.30,5.98,"###':#pkN###^#$dT-hw/TG#mQ#&/0###$##o6#Np695#######/1$okN###e?%gnE5(:6u#W.E`)?i>$qZ%6v'{d($m%3R)###y^%usH###Zx09`.N>H=/++mNFv&8['ME.Y-)Ti.Tp4[l%Tm)T@%Dp1X.-&q9%.'o^*svLmsG|?$}%.qp+&Q$-^'{`>Rl%"},
{121.49,363.24,16.37,1.61,"Yl%p6)#8%CU5o6(gR,l-$T,E+J-[[*)R%u{:xP$rH$Xq8#Q$?w&HBPqn.m,%`C,L8QgH'ec%MBP@GHVH&1m%nZ'C[$~W<>c$+-&F$&$|0HaAp@,rv(_0(G@MA7Qg.,P6%u(2zH(Jr.^]2zb#######7$#{5Q######,%#(6Q^>$rG$k-$66Q&##Rv%GL*l=L"},
{121.49,363.24,16.37,6.23,"###n$#c>P###Lu#Ur.oFFr5%.]Nc]2#$&Nw&sw.Z%%=8/IH&###G&#p>P###Hm'V:(a-L5w+~@Pa?'>-&kV-~../.$0+?x-(al%p%%(?PV#%y.-fH%2;0'@P>ZMxu%i-&}SMx$+_d$4}@p[+pn)Vz,;y77,#M)0&g3uY#h?(C.%P{;7Z$V?'xH$*C59?%gl&"},
{244.99,369.03,16.90,1.80,"###Du#20%C&3*##Iv'd6$z#MX>#x#';$&#h7$##<5#t#%O-)o,$)eCd..G>#6&*PgN.8-8?&3fNx96uH&JT/HH'$Q#t'/b&38[(r8+r>=&29E04gR)<y*PgNxcN:H%4@&DgNDR*_v$}=AE$']$)o.+7w$Zh;<x,]4C>6%@X<9B,}[MLQ%QC6zd)@w(bV7fZ&"},
{112.94,405.01,13.15,1.61,"O7(O=FEZ$0A-D8(KM<FH%^~-Tc$&8&wi>b#%T5$Qc#.`>'Z$V9,ofQid)ru%IhQ^`A3l#S[%E]2vd%pb>w?'###]Z#O2;O5$e-*Iu#9((hQK)eQrv*9c#ui2Kx.w,:ef3pY#4##u*4,80#########;%#[cQ_,%6c$_6$3dQ>5#p@'0M.73Ej>#p8-2W/[83"},
{447.03,567.85,14.80,5.46,"###t$%;[)]P#4Z#P5:DR+~P#?h*Wh;[c%'Z$oc$.,#t>#Je/'H$I-:ZK6yb#n);<$B[x0'L3%47b[*=S+1`::k8#H%sP$@[&Ju$,S%W7Ll)?MOD`d*I[&-H@hS,)B3Tg2F0./&V/c$_6'$l4-c#]I*Cn'Jf3(-&Hu$j7%/X=<[)Xu%U$%`9ILQ'###%##hj."},
{121.08,679.98,16.02,5.69,"qP$kR$$HO_#%ou&Z6%Q9*VC9Pc%JR'uT(OGKf`Be,%{Y#845TR)Fh(;*@Dc%QeSqG$D-$S^5]?)###~$#GeST,%###W##{wP,6&r5$$*+pcSbdS95#j%#B7Kio5###@##o]1&[)###L##[01######/'#m':%Q%###B(#XdR(c$$##1/#XI+2B195#.##oZ%"},
{71.25,47.59,18.28,3.29,"###v@$Rn*=@,/o*,N4=$(b>$<F~7d)###%##(e,#########$##]D*Wo21Z%|o.DG0MM:cS/CG~E[)UG#vn'];?######.#####Bc#Q/*OA~KZ&{G#Cg+<C~g[W<5#F5#.aTSA2######j########&##/95######%##>f3#########-$'#########%##"},
{71.25,47.59,18.28,4.70,".##ZU:######9##E+J######.-#|e2######(##.,#######D##wM~8w,###Po(+N~oR+ZG#A}X3.X<5#W5#mQ'.,#######iG$/d%b*4]?'s&3;/-M>1xU4sK~'?&g>#)i)27,######1##>?'/,#2U+Em&.7,)##XV-7&-;J~###9##QB'=6(######,##"},
{496.71,66.18,17.63,4.84,"$##$Q$ff%''a][*{5%g^&:)a%(agY#.##[m6Q97######|#####d##LL4$m(EC5d0)Mi<~.,;*a$-&xb#kI%m*G######,#####6##Cy+Zz=k$+3Q#{(/|[R;'aG>#^G#3i.eU;######)########0##cc'######0##$?&D>#######`P#############"},
{426.12,100.38,20.20,0.38,"3##TD2Qu%###CI#?8TZP####l,$#;TS5$G>#Rd'/q9###1,#06%Mb3r6TG>#M;Tj~Q_82:H$eK4,V5BR&+17/N@Wc&<,#}7+zY$36#k1Pk?)~7TA5#OE5P~(aaG###/##LJ+/o0/,#'##67(&##j>$-s4du&t'8/,#:*3R|;F962##,u#&J&q#%Iu#D>#I>#"},
{247.04,302.34,16.93,3.83,";w*5b:zI-r>$w>:aIL#-&~#%cW0zS1L|43Z%<v&ou%~RA:T2&[&+m$Xj47~.3W5+T0ia3qy6UUO8g706%C7'0y2}Q'Qy'd]4`bI###[$$V7)e5%wb#mi-_QO#`?F,$ku#i9MBn-96&U,$PS*fPO$##2,#,v#ou&@u#{5$C~.KZ&95#*##*V78m(F>#Iu#9B0"},
{247.04,302.34,16.93,5.73,"*JOx[,9u#P/&`m)F&.Zv&[Z%?B,'7+K>#eP#<J)#v'###9,#n09L%&$93oi-nd'ACOQq<r>$NoH}z=iY#5Z#eS-%I)@#$'l#GH'f('~tGGd%ey5q<0y*EuJ):dOLH'4u#$M&08/3f1@5#7l#1?$&n'/l8'm(H06&m%XU*1};5y5%^-r5&uv#dP#mgO+u#$##"},
{417.40,687.33,16.82,2.08,"&##kyLZ,%eY#n5$RUMzQ&U6)u&S#(8Y>#dl%gL:######'##qd$6}AWl%f5%]/.2$'>m$_|CO%S###Q##k6Ku'S.,####=c#<g/TW=hm*VS*A&Sac'###Fo&+bI###%##KL1N~Eo#'###du#/6'ZG##6%,RA~08.,####]z(5_<###&##}n'8@){k#6##&~*"},
{417.40,687.33,16.82,3.58,"1##EV=95####Sl#iJR######56$G^OuG%###g-(M%*A,$.##yu%iJRI`?)Z#w?K%IR###&-#dz3SaG###,##Y8)ny9###(##'m(kP#EMRuH&Zr?RK08-':[$^B'#JRQu%$##>()MT5.,####ZP#1##EMR/c$}Y$%/&h3Bn5%b'/Zg1w:8|-*&[B_#&1,#O#$"},
{376.67,215.31,20.57,1.55,"cq4m-+<5#iQ#d%*Ov)d>$&##$##>n*b%.######gk4:@,###C<*0RV)##pA/@KFKg9R>#ER&hm*D##gW9;u####IU&xrUac'K.(qQ)E##MWV.GF<c%c,$jWV$7+7u#vI&-6$###1l#u))]m+xY$6##6m'Hs:###,?#%&.f~0/,#[I%H%-QG#*###r7Q6'qb#"},
{56.81,222.84,21.74,3.29,"R80>[&#7Qeu$ay_$Z$a>$V-$Hy_######*##|k##########W?(qY#`)4Y{;yx_$##gG#)F4sy_######:##L?(#########gY#.,#'?#$HJDq<RG#6Z#j;8My_95####/6#o7/#########-c$9J(Nw)#D8a.)KE08f/W-)3$`dw-.,#*##+y3#########"},
{202.22,56.49,25.78,4.82,"l?'/}Av@*,-&G;9{ZF]i15p1#CXAR+8,#6S&.,##########yjCiY#7I%Qe)Vy7uY#qP5X(6~@X###L#$fJ'cG$#########Z};al%rP$N,#gz5Ay*'C83,#hBXEl$cG$&##]Z'#########%n*sX>:Q&6##h&(1FXZ,%/,#1FXbV>/,#$##N7)#########"},
{385.30,58.29,25.04,4.21,"I6%0M6vQ)5,#sx%VdUqb#D,#^h79L7.,#5$#oJ32,####)&#oP#L8,N[D###c*0#eUp=5Hv'eiU,A0H5#Xe%w]7######.$#[P#+d#`i9xb#<x1=H$xE/2'H`dU.,#G,#NR=+6'######&##eY#Q5#.&-*I'P]5[G#a5$ri*Nc&######=&&############"},
{180.21,137.62,27.75,4.50,"F2.+;?1##`6'UD-hrCwb#1Z#97$c7PiS.[5$PG7y5P:c$Kc#Z'4DZ&@##yK1q4JTG#*-%(B'tQ(7d&?%Brf.]QO;H&CI&Yr+{B5I#%###Q~&HQIsZ&9l$[?#yC:</(|B6Ad$(PL%##A#$}A$D[(z%,8?'<H$oR(N7Po,&C##T9O/-MD>#R##zOI######t##"},
{411.03,366.15,25.97,3.52,"eY#4##`g-d@--S*Ly2[,%4m&Hg,mv+%##.Q$`d)im*oP$~H%###[##/UOpb#6&(s7A]/U-Z$x2U]g726&oQ$Ny4/v&Mx+{,%###$##IfK.,#iZ'?l#x2Ux4H7.UlY#>7%5]E@]4.,#t6$8f*######7S*#########LM-cB6NQ'qH&e#%t254S(}K4oP$_u#"},
{411.03,366.15,25.97,6.10,"fK.D7+6#$$##U#$j{5o,$hu&sm$~r>F5#8v']L9J#%###T5#BC8%n%g%/|u%).'`&EpbMiP#?|VcZHz5&HQ#{wVD>####t##,u#Am$TV<Pc%$$']@$?|VoM:_wV&H$(7'2$9mvV######.%#p#%A[(A-(&##>,#d$'o&,Te.hi@Tu%5##Ng-y_@######Y##"},
{72.11,719.78,25.37,1.61,"at2,z:###kZ%1O3dM<L5$'##smPKl%.,#G##############N[<4`?rP$7?$_8D%TVG>#rc%jWVI.-###WZ#D>##########FA/iY#,-$OUV{C:Iv'ZG#${PeRVhY####'a+pb#######(########8##of5eY####&##De)@H'######ul#############"},
{483.95,723.76,27.07,0.16,"`l%5H$gSO)7*4'[Pu$Pd(AB*W&[######c##############D>#'##+Z3p2BC&Vgl$QP9&4@#)[xb#cP#o-%ZP#######$##OG#'##pB(?'[:&2C5#QO,p&[n%[$##b5#u0,6#$######%##.,####%##<c$######(##6#$########################"}
};
siftPoint_save sps7[] = {
{81.33,51.99,4.04,6.21,"|H%gvP>@'6c$L(%NuP$c#$@)0M&VuPlu%MH'Ay%W&1J7+###=/+Oo/|I+T,$TvP6v&(m&z{7Kg+b'1X92@D;6D+P;=M5$a>#1/+ve/]#%:u#XvBvPJ-c$e#$_%*lQ?4e-DB/(_.+%+`l#9y5sH%BJ/9$$H=G;~+ce.2c$C&K6]2Uc&M5#:-<Bc%###+B'X*C"},
{478.23,106.04,4.05,0.27,"6v'gG#DZ$Tc&DMX0,#/,#<5#`LX###&##Q>#|?*###/##b,%yc(@?$kG$I'/KwQ2,#S##$wMPLX###(##4w'No3###*###7(JIX0,#Lu$Zg+Z=K###a$#ZMX?kK###F##?CK.n-###/,#bH#JIX###K,$@Q#JIX$##*##KS+D,$sG$VG#Rf/Nl$Z>$N5#D/0"},
{478.23,106.04,4.05,1.81,"Nu$Ld*###0##8I(m/4###+##qQ#r[-###/,#$J-Yu$D,$K5#H5#qgX###&##E%(RgX###)###(KHtK###K##{/0/u##Q$WG#=5#:hX0,#/,#JwOwdQ2,#^##|hX?FI###{$#%o+kdX%##,##IZ&>6'kG#Sc$Zy/I-)TQ$kG$zp+kdX$##W#%qZ#kdX###3u#"},
{414.57,135.40,3.57,1.50,"2$(WG#?5#s#&_dGll'###L,#o12kbK###tR$--'qo0###Jp+eY####c##b_?tNBD>#1##Zq:lJWgQ'TG#1*-lQ'(26nG$v'.######n##aIWwc(###P##_KW2*>K.+dm#PJWGu#|`64~&:OH######&##4#H-H&###<##:KMpb####(&#%JW######O$#TIW"},
{414.57,135.40,3.57,3.00,"bf+Mc%?V8###ZM2H6')L3B#$iL<j5$VY>s@&E[V###$##R@#_R$ON5?+I###Kj.P~V2-'I>#S~V}N@#J,4n#b[V95####-'#f>#XQG>u$$##N;=#E@.,#C##z~V_Q(###*$#0BOj>%###Y##K?'Jl%M>#D5#r2CD>####3$#[[V######K$#CkE.,####*##"},
{176.89,139.94,3.76,4.68,"vL<L6#lkMN##meWMR$s^;;R#)h;i##ctK~##%##Q#$v81F>#?U7]-$zdW6##KeW=6#IsEQ%#4h<*$$ReWo##3,#7R&+A/;5#Oo3DH#KeSR?%$fWKZ#){9Z%&/V;G$$;eW2l#<5#[G#mJ/uc'j]6N,#Gs@(/,Y$Tb>#YT1G36?n.+6#ilP`#$hG$J5#?[(7Z%"},
{433.93,152.23,4.04,2.28,"Iu#CI(Gr3u%/###=J$/mF*6'.,#s6$=)*51;###M5#~y&6?'UW/w#OxK2iG$u-){w(t$O%Z$VV<@Z$^')[p5[>$;##,xF)[)_'Orf1T,%-##t1:]`-*::V##4ZMTR%L.-Sv#X#%GB${#OwH'i]'x83###%##.+8&:5OG#A##3`6bA,U82Ym%OG#lH#}844j1"},
{433.93,152.23,4.04,5.79,"$x/P/$WS.De,`l&I,#U=6Z_<&@+###+x&aQF.,#$##dV*pn0)-O9@%5v'V?$MA1lG#g$KqZ%.-O'##uB/BB*Q-(.,#O<1`-O%1Ogd+G#$,Z#m/*12<9E8m,%=l>4w,tQ&fc%;}3Qe/J5#.r;}x$^o4.,####EB$N,OfY#)##mi9u<G95#P$#'|<l,&wb#nd#"},
{11.78,170.25,3.69,6.27,"*A0######,##~{n%##u5&J$#`:==##fZRA#####YG#~R,pb#jT7######'##d{n&##{l'$$#vVA=##(qbF##D>#&##}I.'l#5U9######*##4|n<##fu&`##*i?Y$#6A[-##XG#>,#B.-lG$3'7######/##`{n(##jc'S$#F&3R##U#NC##-c#N>#M6(UG#"},
{251.67,169.81,3.96,4.71,"t>#X0TEZ&PG#R2/i3F,@+:5#j'8[,$_5O8##1iA###$l#k#$j7,020cWE$##2:YOd*Mv)C##s4L*##NHPF##PFIPG#/u#?c#je1;##{cS-##Y7Y'##9[*E##U>O3##NlP6##7+J###;#$d#$md,I5#S3EG5#i7Y3,#C6(N,#M*F<##~EE(##kNE1,#{k#)##"},
{274.78,169.54,3.60,4.71,"696(##VM?3##.SZ'##>6(=##7iA7##/}I+##)=H$##M5$_>#n98oP#erBZ>#QSZ(##V6)F##K4JC##R?S*##JdU0,#mP$*##pe1E5#D)88o-:SZ,##_?)8f%UaHJ##rlRI##%HQ$##X>$h##yY$U,%+&$x~XWT6`u%j?)+?:By7pP#|sGGc#bWCzY$gY#Q##"},
{236.99,197.69,3.76,2.35,"v));c%wb####d8)+u#V$%fZ'rb#h>#f$'_Z&A5#~5${P$@u$0c5A,$O##:5#]`VOG#`,#p@)57,(##2{%>:84,#/Z$0p*.H&jM8###-q#GP=k]V%##4,#,W'Eg8?Z#H:,Tm&$##}b#v42pb#qH%95#JA#W'S1J/<H%{k#g*-w>%qC-,S.E##8I#@L6=~(.,#"},
{445.74,201.34,3.95,1.27,"b'#NPN######uT#gh?.,####~@&+@+95#I5#4l#U,%N>#;l$}&#[>O+u####t3-PmV###(##3sV|e2###:c#Zm)3?&sb#^c%tY##m'qP$F#$w%VH?(###L.#XnVcP#*c$q9#sm+l8*aT6Ru#bG$$##%c#?c$PmV######i9#PmVT,#Pm*y'#AZ%z@%`V?}b#"},
{13.18,267.54,3.24,0.07,"Ow.######+##S(e######g$#M(e######H$#8H%F>#$##'##^98######;##|'e######w&#[8^###<,#9S#B,$/,#7,#bP#<'6######$##s+e######:##<`b$##cP#=,#$l#2,#nY#VG#W[+#########E*e#########E)e$##$##XG#W,%G>#1,#Oc%"},
{213.69,352.36,3.71,5.87,"UG#OnKN>#-{9uw#g,Oe@#*::H@%<c%62)I..$#####%C)]#&7,#W=73T3$##AG7tg:;5#;##b0Ox5$j<>C##F>#Y-#0.O.,#.##XG#W`4;A0:K4{b#Y>#=P6lsGm9$y^<R]#E>#j]#t)D###y-#(::%d#|e1b5$q,I*##ZD57#$F+4'$(RZ#95#c-#3-(<5#"},
{136.04,417.52,3.41,1.37,"Z,Bs^:$##>0)9/,z_S9HI'=>K##Nj.l[ST,%6T#t&4L>#OG#(o'D^SW>$$##[B3z_SYQ&V7)wP$eA&*j0T'7n[%b5%B&&5?'L.,}&',o2$##~~Sfl%R5$AH$h-(+Q$6r8Wl%fP#4,#YF?OG#a6)3##uy60@+WYDi,&:l#LH'M[%I/27-&D>#-##p6'N<A,u#"},
{136.04,417.52,3.41,5.64,"_5$>p,Y@*,K2w?$f6)Bm#`U;ocE+u#<##3I(@3?###0##wf+UQ'fG$B##T(3i5$_,%$&%c5O-IT.c#ET(zD2IIT7l#D10qwBUm$]q<###0Z#'d#aPFq5%U#%SMT3E35I)=c#$NTZ#ED]24m$ER#i]6,##Z7.E,#m2;%##Aw,C&,[.)###0p+2b7qe/###'@$"},
{484.88,465.18,4.69,0.25,"8u#Ul$]Q&XH(m&Vqb#K>#4c$LCZ)##vG%'##%-'$~'MA10,#c>#9d)w,%|-+^,LH#$-?&='/?CZ$##q,&[w%Od*$Q#w`8>#$QQ%[l%a.'C9W3:2c5%F>#kFZr|F###M,#DLIv5&###S[${@-U6%ce)Q`>.7*kFZ0%,gc'K~'CN@eY#C,#S.*PG#2,#4-&f6)"},
{484.88,465.18,4.69,1.33,"`%#'IVuY#.,#[%#~98e[%'@+xG#5?'Zc#`N@4l$###-##}o-XK#-IV$##'##`E+'IV(##SG#^NVy^<&##rZ%d7.^P#3,#(R%3H%9?&S#${Q%{$QPd*_#$j6$bKV;2=###R7#*z*-jD###sP#tb#EH#'Z$<6%w;Bn5$=u#:L&l`@`O>jd')~(aQ;`3Ctc$>Q%"},
{8.61,568.91,4.27,6.22,"############_&`###(##t`,xf5uG%i8$R)`l,#N$*aK$Q08D>##########h&`E>#/,#9J$e$V'$%F137K/###2l#C>3OR,;c%#########.(`n,&.,#J##}|C9tA5-((#####/j>An,95#Zu%#########E,`gc'.,####6h0f03~?)M>#.,#LS.2c$cu%"},
{8.60,587.96,4.30,5.88,"############C{A######i%#jnZsZ'D,$ZM#'-'q$'9g0}5$qb##########MoZD>####W'#1pZ(E>L5$y%#I6#A,LBQ&###*v'######G##LtZ+u####%&#&AGEo.Y,%}R'P##j&34Z$E6(9@,######y%#;%X######ZM#$19D>#$6#VN1######mZ#^%/"},
{579.50,588.71,3.87,1.26,";5#SG#/N(|FL.,#%##`,6K'8OG####G.$XlO######Q$#1XG]G#sH&Q:9+c$Jd*O##@'W8Q%:%W+##N.*CfP]?)###V&#T%W;5#+H#lK11u#&[)6##y(W)c#<%W###5J+@[$:1;tb#O5#H~,###jG#vJ+PG#>u$S>#w-G<5#H%W###(B/pG#_^:<5#Bc#8$'"},
{579.50,588.71,3.87,4.21,"j5$;Q$$d((l#*[(Xw#XdW###1fWO,#r@/O>#g]2######dQ%(%#SM:x@/PG#5%+^U&XdW$##kfW~l#pm,F###)8PG####u>#j(#:eWg,&###Oy/coOBtKM,#NeW(c#&m(y@#1/1H>#6l#XR&Y$#,o2######GD*$:9###&##?zLl-+###T,#I^(q~1Y5$h>$"},
{587.27,596.91,3.92,2.90,"Yl&######5##XdW######m'#hdW7,#}Y$M&#5l$:,#-H%t5%=A1######F##KeWGl%###,&#ZQMu<FE>#f##:Z$4QOSG##c#=]4######w##BgWET4&##$/)Xe)%fW&##s&0Al#CeW###.C0d09######<K##&1N5$+##](NI>#7Z%A##)gWgY#.,#(##)49"},
{307.55,705.10,4.20,4.72,"a7B5K4###5##I<<wc'###R,#Wc%(8,4R)4['###89%&|;G6(@IQ:5#[G#3*-uHQN0,1T3z-#*JQaM7rI*W(/>$(-:)kLQn&1id,###du$`SE7S/`5#M7L%~%zIQlo/]e*Sx+'-&E(,]j5~x3U[)95#`G#.d&9?&95#)B)zm+aL:}k#E-$4z3hH'KZ&bJ%=S/"},
{13.44,31.21,4.17,6.12,"z,'######*##~/_###%##uT%g{B8Z#l8*|+AeG#N[)P-'Up50e.######4##C0_/u####C$#Z6O4:)#~-xY#/##{8,L?(###|e2######<$#.0_SG#~G#RN*V^7,A%TS0tl%&##'L/-H&###PWD######-:#a4L)##J|;tLI6#$Y$##C~'n+###j~$hv+###"},
{212.53,69.65,4.57,0.45,"1Q#r[*D>####g>#Y26:$)###A$#_ECmP$###]I$Pf4OG####Il#d'.d$+###[%'=OX/7R2H%=OX=KXhc&Z[%i/.x5&wb#CZ$n>$Uo-`S/###T-(</'=OXM8T)KX7v&C@&K_N`@-pH%3c$yH%6m&MZ%xb#>5#?~.j>$m#$5/.Y-)9?%K>#Wx+gY#aZ$%c#=v'"},
{212.53,69.65,4.57,4.39,"95#,##m$'FQ&y#&q.,Hd&cu%MJ,fm+5##^-&~d+###.##fv%kH)*H#Jg6p>$Lw&g;RH'XH6'h*X}QO/-&gw&a(<###)##S7(5~.P#$e&1>l#%v&9@&h*XS&XG&X[l%-n%h*XG2>###&##~I&jG$WG#.-$'[)?#$zb#7$%h95<H&###+##d)4E?';5#I>#m@&"},
{212.53,69.65,4.57,5.50,"Y@'R5$###zG#sA(>]2###T>#qx*DS0###C##67(gc'###5##~l$,[%)w+8$&7x&1FX$7O8$'1FXVAX1H%zR&oi:y6+###-##5c$XG#he*aS.mZ&^I&|LM,BXxAXwl&hm%1FX=jCX,%/,#x,#p>$eG$+6%o#&sQ##v'kZ$RJ1)c#'Z$s#%_8/xS2/,#2,#oQ$"},
{370.03,71.68,4.34,4.06,"K5#PU4%Q%###(?$i3AZP####W[%b^8OG#@5#%h13l$###|G#|$#QN>g,&###0%&3_PE198u#-;R@{>Z5$V?$2#M.,####m%$9&#$lL>c%###]K)$y0[a9^x.X7R3c$@?$rW-/6R######7'#F8%;C9######OZ<=6(+##.l#%9R95#%##Sc#q6R######N$#"},
{64.93,74.25,4.33,1.57,"8v&SJ+627lQ(&Q3W'5%8/.,#ROYI#%###@,#5I*%##sv)I5#tP#fr9f/3OG#D,DF@&%V;9,#.a[TG#:5#Z##v83V##[$*6##J>#[6%UJ/_l%9ECzG#[C6=9)j~[fP#H>#xU#YA2`7$O..C##d,%7u#wc$VL35w-###JH%A^G}sInl##6&q3,Vu%lf$CH'q,%"},
{484.16,148.58,4.69,1.41,"VG#95####'##2u#ZP#/##}w-F>#g>$A6$?w-4,#kY#P5$###'/0D>#.,#_##7BYOH&?~'NFYxZ'3~&NFYDBYvP$fP#LA/2##tA/bG$1,#F5#NFYLBYpQ'].'gw'NFY&BY`Q&Ru%k,#9w-9##c8)5?'%##jY#77'7A/0m&g5%#Z#N03<l$0##lP#XH(###,##"},
{484.16,148.58,4.69,4.88,"+u#%##]P#dn-i6(7l$oP#Qf/`G#$6&L7%@81######16$5J0sI..##b#&fu#;BYE?&L~'NFY(d'z[&NFY<BY/,#*##Ox/OG#=A.+##Uc&A5#NFY)0YdH'3n&~e'NFYDBYlc&###~##Ce.qb#TG#'##R5$[P#gQ$aI-%##Q,$^G#VJ.%l#yb####)##`P#OG#"},
{294.25,171.46,4.67,5.59,"W(#A7-t[#I#%j($tc(+J#Fl%tM'k6)j,$X>$D,#3&'+_6D>#u<-9~/y?#8m's:Ow5&_S$>n*oB'/g4tC/+[)@c#Z<7HH'`G#FJ-0,#m8#'rUcnUv$#K$).,46-([f#(E?Cu#A6(o>$Sl%ER$8.)o%,{Y#oj>:5#cC':-(KJ+:5#m6#8I*pP$[n+rb#A5#.Q#"},
{449.39,170.85,3.69,6.23,"@v&'6%I&'66'$q('U6}b#B#$+4/{iD###1##fFBaw+|b#S#$n,#(zZB6&ZP#ax,izZ0u#0,#7{Zi@.###R##(i>95#R##iH&Dc#(AY.,####*&W&i>###h##{wZ6#$###U%#_;A###(##2Q#j#&F,$hP#;l$]wZ.,####p?#LwZ######<&#q.0###&##zc#"},
{292.94,198.42,4.30,0.55,"@T&n'X-##V@'9b1hS2aC*Sl%;#$###)2-$?%OG#:##Y?(/?$x>C}TQ%##Np(h*XrZ(<7(_w#i%/Jc##TMJ#$###K##ie/_P#t$,O>#;5#A0Gg&X9//$Z$<`+j[)PE4@_<&#####,6$3.,ZP#O>#Iv'dG$i21Q,$,8+#v'A$'E##lK5M#%###C##.@+9c$7#$"},
{427.82,229.40,3.74,4.62,"lY#{H#T_>wb#cu$X'*RrB###/oSYQ'OG#Q,#WK5e$#g=K2l#1u#N5#L?JQ5$ZR*d?&UgWoY#1fWT#%8Q%^Z#lB7De#XdW7##e>$Lw$&)<%##Q-)%&$8gWq#%.eW}b#pQ'Sw$5/1TR#YdWC,#;5#^c$W^2/,#Vc&jl#j*3p&3n4LUG#_Z&_y-~?)K##ZeW;,#"},
{442.20,229.63,4.32,0.65,"AK4######Bc3cFITH&2,#V<(Sv'y,&1u#@?%=5#0c#b#&0,#/,#1K(Du$#&?8R&g/B#Z$WA&c6BSB35c$n5$Ul%o>$TQ%el&.Q%0G/G@,Y$&C:9Y)1(6&7P7IRTY5$aG#Jc:Q6)3,#b$%im)Fl%CI$:.*_z0H_>9,#)l#L}-ol'ac$u#&9:22u#Zv'u,&R>#"},
{442.20,229.63,4.32,4.81,"5u#?,#'5AT5$^f0($&keIO,$M_~^>$Xo3-c#TZ&bP#$^~###zc&0d$#D7VG#]d*@%$7a~U-(4^~;u#`J/B@'sd,/,#/a~~P#.6&6v&t%+C5#kG$:-%(#3^z<%)>dP#Dd&`PAzc(S>#6rPD[*:,#1@(gG$###l>$Fm&0Z$i5%]l$hY#6u#fZ''##5u#L6%%[)"},
{75.87,270.83,4.40,0.33,"5m&<Z%}b#6Q$AQ%[:5$##E5#1I$ZB57u#:5#uG$kP#l6(zb#%l#(8&705y>%K.'IXXU$PCQ%IXX7TX]l%2@%Bg2<Q%Lc%Bu#$##:m#5N>ZP#8-''J&b|WKBXxSX<?&`m%BET4x1F>#7,#1%'###P>#[$'$?&4Q%.,#?-$,:6y-).c$*##x^01u#)Z$gP#4H%"},
{75.87,270.83,4.40,4.75,"$##M-&Pc&/,#C6$Zq83l$###5B&*K3RG#I>#l>$tl%A?&B#$E>#@Q#Ah8[v)+~&)3TL9X~Z&%=X{HOfl%wd&,V9;5#Dc#~Q&W,%W>#${7'H%PQ&Pd%%=XT8XZ8X<Q%n[%%=X|'9###<5#xH$F>#mP#Zc%dG$:l$J>#Yv%#T2[#%7m(vY#lU6k#%SH'=5#Z,$"},
{227.39,271.57,4.51,0.13,")e(`,%gY#.##Ol#`;;mP#0c$&[#I1:8Z$bl&{P$c,$,]1J,$Ev(:R$Ie/Y>#wR'h*X?&XWu$h*Xi&Xnc&m-%Up4pY#Pw):Z%`P#jw'J]3gP#X6(+A'h*XIvMN&X'-&6e&,hLHD=;5#W,$1$%vm)jh5nl'$##}07'-'Zc#ue/dv)###7Z#dT1Sf0x5&>5#)c#"},
{227.39,271.57,4.51,4.12,"%##Dc$jd%5?'A$(@e+S7%88,$z3([)#l#Xc#el%#c#06&%l#'##5R%fx/H?(wI&IXX~0Xv#&IXX]mR2-&vm%z]3###Oc#Im(2,#9u#SL-6[*MH&ev&L_NySX'TX4?&Rn%IXX/95X5$yb#+&(A,#N@*W5$|k#D?$Kw-/[#7/1jQ&v5%U#$wp2c>$).)oP$V5#"},
{362.02,272.16,4.45,0.53,"m5#D7*%Q$a5%)-%/L4Q@*sG$8f,5I+H>#UQ#>.,bG$/,#eH#'H$i[&{~2)##se',XWG8Wpu%,XW]cN1H%|d%2L3$[)1,#t5$&Z$M5#y&1w@-tZ'r$&A{QDTWjSW7Q%)@%,XW7x1[P#yG$W%*######O,#ai=95####A6#~GLt5&$##ac%*N;t#'###/Z$IQ%"},
{362.02,272.16,4.45,2.35,"V5#gv)8v(ub#uc&&'1<-(I>#ad'zT6###9u#j-'y5&###T,$^##zB1=@,###z$&{NW]?Pv5%{NWwJW8Q%e-%hN<3l$###u>$b,${6(Qo05v%bl%bS'{NW)JTCKWzZ'.7&JEUdHI95####f,#0$(M5#.d'Q?%&##-Z#3&*f83Y#%+Q$|Y#*g3Zp1ZP####rG#"},
{617.77,324.24,4.33,3.18,"OG##########Kyc######[##hCg95#.,#;$#Y[+fY#PG#B5#OG##########VDg######3##5Fg######Q##8J/95#?5#l#$OG##########WCg######+$#nCg######l%#Hx2=5#<5#Y5#pb##########KCg######?$#MCg###&##e~#<-(:5#6##lu%"},
{617.37,346.61,3.93,3.23,"D>##########QQS######U##ZMl######]%#Z83###1,#cH%eY##########cMl######A###Nl######<$#8060,#gG$/c#eY##########{Ml######;##mNl.,####`##xo4M5$E>#C5#D>##########>QM######%##/Ql######0##q'8.,#$##xG#"},
{539.39,392.97,4.52,1.11,":##2v${R-###h,#up5yY$###`c%]m'l#&]m)>5#kG#fI):v(###)Y1P6).,#I/)MWUF>#;6&MWUX:;R5#D],_n/dS&jN3+f0OG#R`)o*Dnw*4w+}U0$Q$MWU5UUuJ-6l$d`0(8-F#6/`<N>#h$*^#$Bh9J@'r.-7Z%*l#k'.md(m#&U5$_v&'c#F?%*D.}u&"},
{11.94,456.64,4.99,0.12,"Yn-#########9}_T,%###pP#w0)@o3,##=Z%&U#+WB&##qb#FR+######)##py_95#)##GN.y&4'?&GJ#_z_LQ#FI,U-#b~[9?'#########o{_95#%##qc$`[P4-&/[$HK3c5$m,%J?$,V<OG##########ez_95####'##o8Y@w+^,%c5$:Z#,/+v#'D>#"},
{9.77,556.39,4.85,0.19,"s7-##########EX$m(###p>$W[%l:=H,$Cv(7?$Ky7?#$Z>$^Z'######&##]AXbG$&##ok5}e0,d)@/$sBX/##@e.w;'Q;AZP##########-DX95#0,#.R%,TUfm'-y->K2###:.']}24[*############VDXeY#######`/Ie:3%Q%###$##~L/~Z'kY#"},
{63.03,785.09,4.75,1.59,"&l#=R#$L8###|Q(ei/Zm*T>#;%,8*2y[+SS,<V3Dl#fH&bm*dP#Y?$W;;xP$a&2Y<,oyTCc$RxT/g-xg5UV0R90EH%dfTDn,^T1-m',@)]u#9d)AI'`{TP~,hyT^S.zp1,:-<S*8%*65EDH%Fd*yb#G.)gR'&##+R%l37*n-c8&+K1P_,){=<o,+S-Q#%-L+"},
{63.03,785.09,4.75,4.36,"?/0i'0>?&2m&cd'-&-ua3{I,cv=iQ(*@&T%'$p,AH'###oP#0IS-H$Y]3Sn$]y4kW,-JSRR&&LS;],EI)^&)me0C?$N-(i6';j:ZC4xf4j['n80%+0DISYHEfIS76%Y?(f|0Jn-IZ$Fw*c6&wH)Ru$xd,do,}#'-H$J?&^dK+80?,#bu%{h/CI+;5#A5#}@%"},
{532.58,838.09,4.57,4.50,"#############92######I0%X<>*[).,#8.9el&=8-Df/@b?############aeA;Q&###h##^Z_Uo3D>#F7$Id(F<9#n,,l#############M7V]>$1u#jR#nT_WG#V5$HL$p?*tu#-A*8I*.,#######&##(T_###1,#v.#QT_95#@,#Fg%(-'D>#.?#cn,"},
{14.28,16.48,5.41,6.26,"ZP#######'##^Q(######,n'.,#######dT2######-##De.^m+######6##uBb&##]G#g.:X82^Q%=S+%Gb8##iS1Al#O>Kbw/######$##NFbOl$&l#n#$/GD52+d07nu%,##O=AT,%###07,######M##WCbUG#+-%y4,lf3'8%GsA3A*%##^i/]m+###"},
{425.90,161.89,4.95,5.04,"Dc#U/R######0-&l0R$##M,#L6'F/RN-(3$$###-45,C6G>#x$*p/R.,#(##%.R:<AE>#Cy#n7.iv%O06z0(###^##u,Cwe0fS-RH'Pl%3,#|1Ri5%OG#h,#`{53%(I7-~R%|Z)-##.W5_26T#%/,#{6'Yh9CN895####sb?hHKpb####y(-]{@sb#)A%(B,"},
{153.03,193.56,6.46,4.74,"hY#7l#n,%0?&~7,_w*~w*$-&>$R[u%O6&;j7yY$l5#%|:*U60u#<,#96&0@)k06]w%Ys9*n+O1^v#%Q05r[&IR+<6#.0^J#$D$)YG#P?']l$%)>0Q#843ki=t/^U>#K<:`T+8w-[5#h0^:##M[+jY#_P#ju#t*G]5#F3AYv#~0^Dc#baE'$$E@,e5#J0^$c#"},
{14.82,353.88,5.22,6.23,"P6)######,##-Eo######E%#w0e######w##l5%tb####+##o.0######)##1Eo######L$#9Di######X##b,%.l#eY#5,#4f3######,##@Eo######2$#{Ka;5####^##Kc%<#$F>#N>#)]3######8##0Eo######V%#4|DjZ'###r##ec%x,&###$##"},
{125.56,445.01,5.58,1.27,"cQ&zV6bI+xb#2J('R&`@*.,#jh5.,#.,#'##([(OG####NR&(o.C?'cH'G,$fm)au$yU9dP#n?M2,#yY${5#w8_7##ZP#2u+%-%+%+=A.)v&Zl%lQ'8d)lf+N$T###hY#^p'{`g~6#Gq</X#H>#WG#y6'X&1~P#95####yB-#g3&l#xb#1@'Afg}H)Po1K5#"},
{484.02,462.91,6.06,0.30,"'l#zG#ed)*c$$:TmY#^P#~P#W:T`,$mP$###3u#7])m-+###2,#UZ%57(06'y-P^P#I#$IS*69T.##,m'.@$RQ'T@'p`6/6'H>#R5$X{){6T5_7gG$Cc$HCMK#N###<,#xG7U,%###ze%^A11,#|P#7^Sg@.I@KbZ&le.<:-{uN$##=5#q~'.,#%##9n*NZ&"},
{484.02,462.91,6.06,1.11,"{U#SsG.,####K'#9h<Qd).,#c##Wc$aC2T07E>#E5#^,$qbD-B#MC;c>$(l#JW*a4L)##R>#F'Jpm,A##wR*hQ(###?##9y0@5#lP#V-'{c'Y}E]u%M>#ZA%~zS@96###q.#kh0[?)3,#k#$###.##D[%$7Q:I+xb#Kv#nUNJ_3|B7?&'kJ-,pRXH(lY#P[$"},
{26.51,505.07,5.14,1.79,"###K##*jAZP#;u#ML*rEF95#DuIZ/.W>$2##bf3LH%4v(tP####Z-#y|D+u#<I*xR$`eW&##geWvu$Yd+r5#%04#H$/7+`Z%<5#17)0)6.$(J-)*.&UhW,Q$@fWFo2y.-{c$xe-gZHq#'=5#%##phWal&/,#H>#DjW[961,#HI'peW_,%A5#V>#+RN######"},
{26.51,505.07,5.14,4.83,"######rP#F]YR5$F5#6/)K~Y6x1>5#&c#c`Y#?&###&##l^Y+$(gP#A&/<=Ad@,RR$A]Y0806^Y=Q$3d),/)-E;x5&?5#9y2I$)oQ%By6/Q$id,Z-#W~YoH%U~Y)##m?)rS$8N>L5$###<R#?c%[l#dL:o#%I#%=##>aB^g2ssH.,#VG#,()bD?E>####.6#"},
{590.69,523.63,5.42,3.45,"mtA#########|BH/p2###/H#.?&zf,###^('c~1O>####]^#5z7#########$CZ,u####SK$_08RG#0c#-c348/8H&3,#z6%ef4#########dDZ######`,#]?Jkl%HQ&}>$bP#^1-j,&=5#>$(#########PFZZP#######=?>|2;ZP#%#####F=/Bu$,d'"},
{479.50,691.27,5.05,0.09,"^G#iA(}3A^~0mLVwv*oZ&sD3rKV.,####n##em*E#$.,#`G#######Rr(?IVgZM###%9%x+C$NV.,#(##r-$MA-o?'.,#<,####iP#1u4#95ln/$?%gD+6E<QJVJZ&6##KD)_m)3@)?5#|P####wn*Iz:6#$mn-w-(6L6@R(D?K0Z%$##j-$[w-}G%E5#T#$"},
{50.30,746.47,5.75,4.43,"pz,2-(,##u#&h])6%-|##9(;v$%'.,fA#NQReP#E>#>z#CQRmkFDc$y5&y6%cT*:LMzL9qu%QVRz<D|P$=I%l%.RG#lc%A$'|ZPoY#gc&-S&^-)v?%zb:/TR@RREl$3-%v~EmR-aP#1?&cQ$He/###&##+T%,05###,##yK.Lm*H#$n>%Lf)4u#]#%O?([G#"},
{50.30,746.47,5.75,5.32,"Q#%hG$###kAD+Q$K&0I>#pk3Sm&h$+O5#b$(bu%;5#F,#^#%EH%%W<Ru%y~'nI&v`V_D>Ol$V)S2>JY#%fZ$.o.<#$Ac#*Z$/T)$h9el&B5#37)_7(wdDQW<8]Vm#&G6%HH8Y18M>#Nl$'Z#UN,&{==&,Wl%BX2*~._>#aQ&pz5%Q%%##d?%)p-D,$$l#*##"},
{454.52,801.09,5.67,1.99,"###<5#j1%rU;S6(^Z'i?$_*2ku$o#'_p,0x(/,####7rS[P####qY#.q0mZ(E%(K>AafRYH&/{Rrh=^$'_$&n$+###7rSsG$yb#3,#mn)wl'bu%:Q$`nBhi@`mSmG$'S%g&HI%.1,#P{18?%XG#eP#?l$u#&95#fu$]d%]v)J,$jG$aZ$w/1)l#H,$#R'gc&"},
{454.52,801.09,5.67,4.49,"A7X1,#.,#1/#Dw+RS,L#%wP#Aw'v@.$c#9u#--&N5$L5#F-(K7X'-$Pc&'&#$8*r[=?3Bn>$%=X,*>m5%RH$fe.UQ'3,#3H%Y7X@$(9m([$#nH(jS)%LIHC5G8Xc#%6[%FG3903R5$N5$a5#]`A[*CUG#S5#$T,q5L.H$Al$6uGq>%F5#i6%Y$'4l$4,#W,$"},
{478.25,50.90,6.73,4.86,"*H#Mo2###K>#_?EYl&###W##JTY######0/#@_=###@4F]##]Z%-R)P>#T::CASpb####[02oUY######8m#&kH###h;<}P#.,####c>#kSY>%-W,$0w(TUYcTY?l$/u#Fq*+SY###2/->Z#######I$&eK7G#$1H$S|<sw/L5G$?%sZ&9U*RIV###AZ%x,#"},
{507.41,50.87,6.29,5.29,",`3L5$###u@'@_QD>####.@$[~HOG#8##]P#^r7###6i),?&&?&###%H$@_Qrj>W~.uJ.eE5@_QZd+UG#ku#$yI###U-$L>#######.I%bN;YZ&3Z$`|1_U6;~Qq>%>l#]:,z24ac'###/########}c%|c(K>#m#$@D;}Y$l9'Mo0$[)qY#]]#9K5######"},
{590.08,123.78,6.77,0.82,"g;Agc%[P#OM$/(9_m*S>#sZ#`?%d$+Ga,a>$2,####Js+5?':.,R5#E[*$)'SR&cJG;PG+$&$DNQXEbH%Cm%7J.:5#8W'd-*ZQ'&##-I'r&-{>%+v$On?--OZcNNu$Sv$;:Kt7/,c$^c$l5#BQ&)##WG#KH%wb#[G#ol#^7.}k#K>#%d&(8.lG$`,%#l#bG#"},
{590.08,123.78,6.77,2.79,"w#SH>####;(#?~/@J-iY#D@#g#$:~.7Z$vb#+Z$OG#(##HH%x#SQc#Yd+;(#-.(y9LbaE~u$2)SG5JMu$9v$V/1xb#SG#O#$pS2Vv'bw.X.#Pc%YR&.(Nm3B-%Sq#&E-%K7@7T17Z%1,#Dl#`G#C5Dx#'D5#Ae##$RIu#N$)5n$_g6sb#z?(^e&&R)###$##"},
{451.08,128.45,6.37,0.80,"qP$###HR'$l#>bF###zb#:.#StG######H&#z%1###1H#x.B^P#hl%9[(1,#4mH2,##Q$cc#YpT:5####^$#TrT#>Egw*Sm<xY$Y5$Q[#MH%^C;;5#1,#.;&FnT[A0###)f##=.YnTD>#+##6#$###o&#|o3-&05l$G,#i9,A['Ty1hc'GQ%u~#-M7%Q%###"},
{411.75,167.72,6.02,3.30,"3,#((Y#I*<5#Au@A%X###7,#4)Ywu'###9,#d-)cG$###&-%;5#vG$|05caCU#OeY#r##.6K@&Y###+##Ld%eJ2.d&zY$GJ*O,$|k#4J'1SNb7-&-'T'&DZK7&Y]P#ZZ#`7)[M>fQ%EJ)mZ'A5#+1/hT/=m)hQ$_4Cl-)iY#|(YR..]G#/l#`i8<m)g5#s>%"},
{245.45,183.51,6.62,2.51,"jH#R7H=020l#lL#u(>;Z$7#$)w%a5%il$cH(E>#;5#EH$Ml%U##Cm'pX1an0'j,{70wU#9g6nVS@u$4?#Np.%m(_G#aU'_@,,##rb#S($uR.k$)]P#tC$j6NXRRj#%_?&96:iH(I8)*(2W#$4%#0S/>7$bG$:-$PQ'_^%^u%5p%oI.td%r5$2.#Zm*?T)eY#"},
{245.45,183.51,6.62,4.57,"[T6.,#'c#/d#[]5Xv)<m%h14*[)R@*6M=Go,I..K##';>RI&AL75m&nP$t5#PU2jWV][,'6%f?=8TVMJ1|G$qz>k,$,o1c[#qZ'D-$Py8yG##i>YN*/RV]##fSVw6(XZ'}@#'OG8,#+I*_[#lc'K,#Wh9cZ%y^<e##3UV$^.1RV$##uP$2g&b)C3##6I++-#"},
{280.65,183.82,6.78,4.94,"<n)me'Dz;&##f,C-=5nRU###+VU^l&9#$YG#u)@###}w,O>#3[(lY#?W=}Y#w*D&##MWUov)~UU###4$&De'K[OlY#9I*N>#?f0+u#&d&}Z&)e+'R)X41XSU,EBV~*r-&F9H0z1#q5el&Y,$td(uc(###E,#?d$QU8`G#`l&fG$v90hG$:6&|>%CI)95#no/"},
{421.56,214.56,6.36,1.66,"q%&rJ]{k####zBFrg:3l$###vS,g8.=R+%##oo0V,%H,$5l#Uq8#93{w/,l#=N]m>%3$(5##$PF*H$IbH8,#=-JRG#3Q%:,#zd-*##(%*6L]dJ]###W#%7#5ANC;5#=j:P['=<BO5$z>$h[&D>####+6$GK]X%/[P#H>#Nc?V6(Wu%s#%ZZ%0H%BH'*##fu$"},
{547.24,406.84,6.27,0.96,"8##4].;Q&###1c#+8+E?';#$,##Jv$/%,-u#?5#;'$rPL###cv&k:QUu%fZ&.9Q?h:Wc$=e(B$)#y&^>G$d'.Z$hL$f5Q###t$*)S(2@*'t<37Q(g->m(g9+yu&{m@s3A=5#6Z%n)/K6QG>#n{+*f2PG#}>$)K-q>$G[(r5$0,#@-%E/F:5#95#`l$,49QG#"},
{547.24,406.84,6.27,4.97,"em'Pc&oG#x06XG#_#%+&+8q;2J,2/&&OCDd(YSSr,%X>$(-#Wg*UQS/,#<##0TS,T36T)ze+,J/hc#;d<1?EgXI)##zY#?3-emASy85,#E$#nVS5%-yP#z$$[8+8U4f8)S%-dw-T5$rY#C@'NO0Rc&E~#>.-29H6#$.##.c$E)5Su%###vY#?['3Z%###hP#"},
{104.71,427.26,6.54,2.49,"D5#+j@8J)L/26z+9GMH-$4n-BSBI-T1,#UG#%%'A1Tx6*$##vG#Tq4M3=3]2RS-a8+{91{/3/0T^5$f>$>m%%8.`u$}d-&##w5%&J+R+2IX@7m(dZ&/S%|1T-:9DQ%~I(8,?Uv'he.7-(H5#<5#X>#[2TB;>gY#>5#Q$%DW=###Q>#1r4>[*^$)Cl$$~)Rc$"},
{104.71,427.26,6.54,4.23,"###:c$rC3kn.%##V[)_Q&'n,(c#Cm'yw-U$)*H%VZ&aP#[o)%##AQ88nThY#'d%WqTzP$H,$vZD{98q>$PI(q./###&##D'+`v*jE-GmTn>#I95Wf.#J)=L3zqT,A/Ed(:@'4QEs5&<5#zu$B%-(##A@O#?$;M;A5#EK0g:3)V82?&dB1r=<~mTvb####~H:"},
{13.48,483.87,6.65,6.21,"############6S[###&##f1'trB{k#88#ddIBA'Nc&=$#K=GA,$######$##kS[A,$###B7#|S[&&.3$%?7)=I&a$*&H#yy8uG%######&##mV[)x17,#Q?$gE3sS[_>#g#&fg$^#Q###E>#r5&######.##w,Q4Z%U,#=M/Cd&Jo3.J#w?Sj0.<R+_##-18"},
{540.05,740.86,5.83,1.37,"###)##xd%nS2+##`['dK(i;Bgl#fH(ZU#MHS4,#*##vU#FHS###I@&RS/A,$S6%Tl8ZZOCl$&0H#r:Rc%>d&?I):y'Zd(7[*###O6%B~,1u#`,%B-$v&HI<AH<B,H$(/,s@BR6(~M'KHS[>#4,#T5$Nc%A#$###_G#(d$1%-###v,#Mf2{I-`P#T6#nIS/,#"},
{540.05,740.86,5.83,4.78,">mG{Y$###3##5m'C/0###?##Q-%2n-###4,#6Q%@u#:5#>,#L/UOu#'-'cZ#,e(v~DwiAR,$_^LW<Bm5%3-$@x0~P####5d$,}I^P#sH'+:'QH'J6%~AD[_:K.Uc,%|H%>-:pf3'd(###5H#Jv&P-U)##9-$%?#e-U]Z#G$)hR&X-UH5#;I(:-$$.U###&##"},
{132.60,805.31,6.16,4.84,"u7Y;##t$)au#aZ'0m#UHGx>$pFJ%c#MH&|9)=6(###O>#3;499Yk5#T82,##GS,3q'f8Ybl%a:Y4I'SH'Z7&K064,#oG$JR'38Y$##$7((6$Jd*R5#xI>S6MubM`P#XZ$NoCC~.A,$###[Q#Z8Y###&##_>#A.-###@##|7,6c$<c$C,$;%)9$&e5%95#*##"},
{148.71,48.20,7.33,5.06,".H&###q',{U7h-Qp#'2##D2)9.L-H&0Q$T6#Z%-###J}=eY#eY#)##Q0+Sq=l1/z,H)>B_@+&{R4M=~5$1@%g*B%##;wR8########`x)XZ'pl&EH%x#:xwRkvR:l$m#$jnCFD?$##E}A6$#ZP####@$&/c$gI-ZP#;##m~-Z^8D>#$##$/(Yx1$##n?'_G#"},
{596.61,71.83,7.53,1.15,"D?%B=1B1<V##9R#w;=TH&###S##GA0|8'I#%'##95#a/%*6'p%'.Y3jjI###Gq)L9R(c$$##Z[CB*AC?$XG#5.*vG%}U$t848@)He'%7RE>#L_:Ui1c?)f[$q6RXH&)Q%j&#[.+K<;T.-L#$%##=u#W7Rs?)u@/bu$/?%4}1Z)B(##%v&XB$}>%{?&R`6Z]3"},
{596.61,71.83,7.53,5.54,"a>#D8.|r?###Ox+6lN%m'gZ&O*-|p;5##[w+aQ:+u#/$#mH'F>#4,#>c7?d*7{6K>#w)+CT2>8T.,#:##fO6blEOG####lg*iv%gB4M/%&m(wZ'qb##V#WmSc+KgZ&K##l-B$]2SA.###0d#vD:0o.Nv'$p,`e/G8&e5$~~,,u#/6:aP#U[)OG#eR)###)##"},
{55.60,121.69,7.20,1.40,"z#N:,#;Q&:##w-)uv$Og9%##nAU1Q$xY$%##on.[P#.,####-GLMv#T6)g[#0/-$D&%<AKH&EBUz#%fG$/d$G'6D>####G5#j,&]l#{d'u#NYu%}5#$g.dAU?aF1,#~G#6n?2e.######X6$]###@UY%#k?UT##g,Q9.$.@U;##]1<|Z#<.Ta>#7m)7##H_<"},
{55.60,121.69,7.20,2.60,"q##64H.,####V$#4;?######5m#~[,######Qm'OG####4##`##[V=kn'*6'^T%'[Rs.*3u#pIA82@2,#kG#9V:B,#mP$h$#7$#|70we%^Q(>[(3d(Dj.0o-j[RAu$m>#*W)6sEQ%#9q;U1$z1'8iAZ5$$##t^/gc'^[''[%IZRF##(6&Rz$DT57(#MZRnx#"},
{393.08,166.46,7.86,3.29,".&(wf[fY#Hl$#qIm_@.,#G5#5S,eY#####H#`Q(######uZ#-M9%y31['&,Gjg[?u$+##<8)x(>..*N5$VJ'%e('&1.,#{$%_l&u?(h@&bJVZf[KZ&9H#By.a=HiI*A&)_.-(/,(/0Pl#np11?$(ICrZ'8#$]k[^A11,#A5#iG:/S/V,#LZ&/Q#XT6/##'I("},
{229.60,195.46,7.59,2.81,"Y##N~/Vc#zc(3l#W>$a##G.->l$^5$el#}u';5#RH$Yl%W>$A;&dmT8[$?e.TrT8~/?H#oL69I+)Z#'<*M3D###x>$1&%q~2<w'.d))V#AnTQpTA6(uZ%z~D+C75J*BZCbv((##BT0CB'o.0S5#H7-~&#8;?&8*TJ1@/%bT2b@$I>M>&+(Q%d$#0`AF.$Y83"},
{19.03,292.29,6.76,6.21,"@H'######+##D)m######(&#noa######b$#X>$q>#a5%3,#WI-######'##B)m######z##L)m;5#}k#~##/,#W5#'%+F>#SS1######:##D)m######D%#[gc[P#.,#E##gP#H>#4Z$D>#fS2######J##A)m######;&#:jF.,####U##1u#}Y#D>####"},
{130.34,354.33,7.29,1.31,"Cn..,####>d#Yd)yG%Y##7AKg,#e^8ny$]_?###UA,sR'ZP#ykM1,#[P#o.#;aDM;3zo2)F3+^.j'5{l=W/1kY#g#$dnMG>#7oMM5$#l#4###$=H{=VD57//VA/iP9=M;=U0RG#HU%kkM%##&7%ou&.,####^/$ukM'##OG#@,#)oM5l$QG####%I$;e,D>#"},
{130.34,354.33,7.29,4.40,"7S-95#####@$V,%QG#B,#ZJM'##PG#U&$FGM######k$%LH'=GM$##SG#&z%Vq:u90aJ/O>9TD5ge.CQ<[2>}k#7##eJMX>$7JMG>#_P#Bc#2$>,f0{T.N^5'y2p32G|C4{3fY##/#IGM1,#V@'OG#$##$U1r'%IM?9?#cy7Y##,AIdv)fG$###wH#G%..,#"},
{130.34,354.33,7.29,5.62,"/6'/,#$##;I$a5%qb#?##k8M&##|k#o&#*5M###/,#X-#cH(-5M&##[P#Sf#pV=rz5NA/|t8vh4K]2/*/k_=qb#F##Q?>AH'08H<Q&PG#;##EG8FU7tB0/f/AL6V$@OA/g0.FZ&{B%Vr@B%(/o'UQ'###cG#g^&*jE1,#wY#>$']7M.,#B5#.XE2w'eY#DA%"},
{203.63,355.08,7.29,4.01,"Su%###Gl#7/)oZ(.,#(##c#<le,###'.#wGM?E.1A06v#6v(K1<(##O#$LK&jC7Be+zm)>a1>uDZ$(IX6kjC[l#S7)EJMDZ&m5E1,#]5#%d'f??=A/}7(-V65)=H+5@'6JU/^P#`U%S3F###y?%X>$.h&9_=@C%Jr@W7'E,$Q,#8SE8#$G>#1,#>I$$H%###"},
{277.45,354.21,7.23,1.48,"jZ(.,#0,#+~&uG%TG#a##y.N###^P#?g#7kK###$##x6$KZ&8,N%##ZP#ze#4N@e^.%p2'5:&U/&T/lP9;^6]m+:##b7I|Z$l/NM5$[P#C,#i>:%{<kz3P~-D8/@,9GV;Jy.d2@B]#D=II##u$%l>%###ql%}~#5,NT,#fx3<,#]/N3?$0v(Jn+]v$+c$3,#"},
{277.45,354.21,7.23,4.33,"M5$$$#0cL.[$#c#_P#Bu#[~KTd#1I+dJ#gtLh##V(9zH$T,%U*F7e#^~1)]#vg8w].@8/[k7Y14*/.Pb7?V;[P#`>#DxLcG$/~H+Q$PG#7##<#:+T1H9/Xf1to2Z|1<M<j^0ZP#4K#htL$##8I$L5$###H>#l'$zrD###tb#l##]xL0Z%^P####[m#<-(.,#"},
{277.45,354.21,7.23,5.45,"u>%:|3g,&Mv#|k#I6$)6&?eH%##PG#V]#8GL###A5#ld$5l$^+K26&QG#0T#Nz9m'2TJ/0Y6*h1H&1<F3a^8gY#:##&JE,u#aIJZP#cd#nQ(tb9OC9oK/`e.cK4>$?<x/RK.{k#YU%r<E$##=R%k5%KP8Yl&}n#?kJe$)D>#=,#YeL###;5#sb#&e&D>####"},
{274.32,389.43,8.71,1.55,"qb#)l#f6#1vO.,#<5#Ux${n1{k####3c#;U3######--#2)?h03sB/Up.5X>L]13[%ZCPOw+u>PTc#7^4{=5S$'rA-A)2<tEH^+j@P1A*@A-`J0yj2)tIG?%ZCPAf0CB0Dc#.J(#cAP24*x/Jc#/M;~-${*E&Z$2/,+&'jv+kn'3B3DQ&>5#|5#5ZHLl%###"},
{274.32,389.43,8.71,4.84,"C$)$##V,#5W;J?&iY#&7&{[,Bw%rf6hY#RU/[>#+D9N-$T=F[y..8/Jn&6oLj&0B#$ZCP<'46?PqZ%gS/F3.&y-Px.;h0|tD,s68`=07(1B.MU4^|0~?Pu#$11O)T1('3-?%(p+4.OU//cf.*R#MFG$##mY#V>#ti9;c%%##K[#?S/eY#%##CH#v<E###(##"},
{274.21,425.14,7.13,4.07,"$?&###$Q#&/(QQ'.,#<,#5?>_P####OT$xkL+u#/,#tv#IJ-fq>J>#5l#xx%N17>w*An+rE1MuF<m'Gk6.OBA,$n,#aoLF6(UHFPG#<5#E,#[Z=Q~/<J)eg40r=)b5Hg7}9-6#$.U#WkL0&*%I${k#/,#.,#hp$A{@###SG#~R)_wH{k#]5#`C<d-$Su%~(("},
{463.80,490.61,7.56,0.92,"'8'I#%3-%2$&x2U6#$$##dd#MRP###:##RB0[P####BS#{e1(c$###kR$j/Uj]R)c$N.#qHCx2U/u#.##`R%Wd*?c$N-$Hc%Ny5E>#o$#//Uc6)D>#[1#l/U6+H6H#Gc$Kd<jd+-@&}k#?##7;>.,#0##do+vQ)%##Om%o97OG#,?#/'3:~,q#'4c#I,$s?%"},
{126.94,494.94,7.08,4.79,"(c$###%##Rn'8#$D>#7##LoN'##95#-]#2lN###5##(n&z#'JbL$##D>#Bx#:2<cB05J-J>8;s;kI,]>9SW@ZP#@##;pNI#%ASMeY#F>#5##06<;y5?],191r1:/H?I;>21.eY#.]$XlN$##k$%3l$###2,#YK$0tJ###F[(8Z#HnN(H$E/1###l.'Vd)bG$"},
{199.54,494.47,7.09,4.06,"o#'%#####Ld#|Y$.,#1##FQAm%(D>#7A#YPL3V'cp9#d#cl&hDA###:5#hx#HC64A,/e+wj38s;b$)L59UW>|P#g-(a8NJ#%$[E#?&0,#L,#X-?]~/Ze)WC6w'78F2r:<hL0[P#O1&:5N###^@$w-+F>####6h&vDA###QG#4H#'/F@H'G>#OG#O.$bR-$##"},
{290.84,496.64,7.49,0.46,">%)(YZ;f2#[&ES*w40.;<Sw,V>#y$)]v(4l$B5#;#$Ml$E>#j*8jyQ:.-;,#|#8JVZOm)EQ%O4=^PE}#'bG#wI,E,$wY#YG#U1;j5$,_6%w*G04*T,OA1Ax&(xYbc&pP$0:$Bw-2u#/,#A,#:#$J>#tC-`)B9l$vb#|G#pG9lH)E>####xf%zG%c5$.,#?,#"},
{290.84,496.64,7.49,4.98,"G>#d/%kU:;5#d5$Eg/]8/Kx2cE9pe1s,#a`<25~%Z$&Z$6f&F>#<##p#=$l#:o+cR,{1-~e.x0~fZ'DQ$?]EA4X1~,J@*cz3F>####h/(4-(76'qb#]V%86Pt[Qsc'0Z#Lj:9%)B%+`l&=u#X5$%##7,#/c$hG$>5#R,#]R,|Y$oG$G5#nR+pG$|,&###)##"},
{582.95,608.89,8.57,3.14,"@H'######,##*B`######Q%#.C`6v(###E6#3Q$j?)/,#EZ$AS0######4##mB`mP$###e&#}@VT::###*S#<?&}=C~P#'d%RJ1######*##UF`:@,###rv#Mr6z7X###<U)`c&Jh9###Jh0(/1######:##jB`b5%###VC%X&0+K3###Ki,^?&u7/.,#>o+"},
{74.29,829.48,6.79,4.65,"############+^cgP#ZP#n##.sD._&/%-I,#Q6%P10#x,O5$OG##########?^cA5#^Z&06$h|FYI%#BSfu$se.)~(yr?V#%}>&######%##J^c###_u%V?#d}K^P##O`_#$;c$K[(A2:sS*SH(######/##|]c+##0u#)K#ag:([#--KgH$>@,%Z#v/1Vp("},
{390.23,142.16,8.69,3.01,"fu%qQ&N23U5M3&Y3l$H##x}1c|Ff5%###,B&MZ%A,$###,S'TR(Ks>rn({T6@&Yv5&:Z#`3,_::[@+*S(P8*<;7N?(G5#vw)~lGUuL####Z#:*YFl%'##W,#kr2Rd+zl#%%+lA,zd--##je*b=K[P#%##[I&l'Y.,####V##W2,yH*######bR#Yw/######"},
{240.72,353.84,8.69,6.24,"95####WH#ZQPqb####y@$$94/u#/,#wY#4L/95#$##IQ#]uMC@(NU1E8)WQLoR+L[&rUPzn/7QPZl$H%,y<0a[(Le,Bo)@SP?8(C@MSv)R,$O&/v*1p{C7,#rUPk.-s-*Zu#dx)ytA|R+>.*H.'MjE###&##U$(]J,eY#~Q%5')nv+3,#(5>I6#HaAXR$O)@"},
{471.68,354.67,8.63,5.03,"E_1-S-,?&###pn'.q;######$R$cJ295####UZ%6#$###$##HW+YmQOG####dD*aoQxl&:Q%H(PIECJ>#F?$Y8/pb####7##yoQLm)$##%A(zmQII(?I'YJHTmQQ#%_G#.Z8B/0+u####O##</2###=5#*qQK-ORG#RG#8z*x95D>####A.&:Q$pb####(##"},
{163.48,494.00,8.58,3.30,":5#:5#<##i::4l$###h##3e+ll'$##%##LX40,#]P#CR#2@QzZ(c>$97%?AQ<810,#|b7~C8W?QxY#bI)z),,e+cv'9L/5a>y&0YV9[['KU2Y7,=3-]@Qo?&Y:P1r>9w+~l%Be&~BQdd),//F$$r@Q.,#gY#T##/Q?_5%$##A-##V7+u####C##ElHD>#E>#"},
{584.19,564.91,8.68,3.28,"%i=#########9P[$w,###:##7y.:,J###fv#5$(4/.###-x&)C8######'##_L[I#%###J&##j?;o2###R2$gd(3S/###D8'U]4######$##FN[Ru%###z##8b@}C<###o/$r(;kv+###0C#N..######(##jJ[######i.#>aE~6(###cz%n7,?E@###g-#"},
{135.93,608.84,8.77,1.11,"###[$'LQ&G,$1##P'/~>$N,$Y-$|I-C,$###L,$L>#rb#7,#%##`S(TR,###cm%Z(P/6Nu5%YDS>[O9Q%_-%se-`T.###'##|P$wP#|w-xb#fu%SR%}qRaRRA7R%7&7B/bzPQ@+s6=T)A'Z#xd+###.##kY#C,$%##36#q/3.,#L##B<A7T12,#[$$w^SuI."},
{135.93,608.84,8.77,2.37,")##f>$iP#fY#2##`8/u>%95#D##6&.>0R+u#O$)S@'.w@=W=QG#S?$xJ0X>$v6%|1R(IOf,%>LM/.R,f*rR(zI-95#A'%`-R8#$8##(h6TG#&?%E[%|1R,#J1.REZ%XR%}pN7f295#f>#$@&###+##f8+>u$<5#WG#w['fm+L,$^#&2##u~-d5$/H&###2##"},
{135.93,608.84,8.77,3.38,"###(##4E,_K7Zc&'@(LK'T'TA&.V$*$##A=0:u#:#$%##`m'###HH#?A/+6';[%R0J+SS]%.[qPAdO<Q%S.'po/a,%95#9,#.,#E,#192r>%|>%R[%#`S;%T$-M^l%AR%O)T8J/;l#W>$7?#X,$JZ&hP#jc&fP#+u#tH#WK6######MQ#XU8$##Q>#V-'3d)"},
{135.93,608.84,8.77,5.04,"$v%S>#^>$$##-##(9/X#%OG#bc#h7.4l#*c$fI-###=5#pu$(2(Px.h$+###eR&lLO=dO:Q%&{Rk>KGZ%;R%RL8######,6#&{R~m*4R(.w%QA00R%)'I9wRyvRq5%DR%&{RJy7ZG#bG$K?#0sE###'##&O.Sd+###+6#rn.tP$.,#B,#Tf1x#&J,$s>%B#$"},
{304.20,608.30,8.63,0.87,"H>#b-'pb####<##u0595####[H#2]3######zP#c?)<5#95#$##bd${n1###R[%[2Tm6QSu$eLNX.TDZ%M$%U&,/n-###&##&##:Z#-g3.,#ju%*w%[2TG,J(/T#B,88*^9JFo,*0TW>$7##L>#06%/l#.,#%##d[%/~'6I+###uA%W?OD@+@w);P;E`=@[&"},
{304.20,608.30,8.63,3.79,"A?'~c#mh.S~SHaE-I)L>#mb5|n+&[)###kn(C#$###I>#Ud'###:##+C/,kIum&#KF]~So]1z_STPHUc%]@&jf1###UG#BQ$###&##JT-^Q(6Q%Q$%1qNt[SEvPk,%NR%z_Ssw0###$##3@$=5#/,#eu$vG%.,####(d#xw0.,####M,#pp6D>#(##`P#&R&"},
{304.20,608.30,8.63,5.18,"g.+oG$.,#0,#+$&A92$##wb#7@$6S/###%##L$'OG#+##UG#E0*&90zc(###g.'nVS@>H^l%nVSS6Pu5%Y[%zy6.,####K##S_QzA3t?(N?$3T1`v%3xDaSSdRS+H%^-%NzOxw0.,####Av#haH:5#S5#pX1r~2###A##{p1}>&###(##~K.fY#/,#$##Xc%"},
{304.20,608.30,8.63,5.86,"Mu#k5%95#)##X5#q'6.,#$##S@%uw0######pQ&C,$.,####UZ#T7'*~.5,#2%&[2TilMhu%[2T8dPQu$H[%9y2OG####7##6$#%i;Qe.###8%(lB03_OL[Ob.TPc%H$%q9K&U7######%-#V-$pfL<7+~>$D&E]U:]>#UI*}h:###%##,A)e?(D>####6##"},
{470.05,606.23,8.39,1.39,"$##tY#{%)mH)0##G@(M^/{k#uI*@v(FZ#rd*Y//,u#'##HZ$I>#m>#<z:f,%{-&/I=DeSEZ%+iSAs?ZH'-n&>&1|M0[?)Ml####.##*D:z>%Kc%K?$+iSmdSGZN#?$(y-+iS-R*^z'{eSGZ$.,####F$&%Z$######{m$jo5######zA.iB6######]eDsQ)"},
{470.05,606.23,8.39,4.48,"e'MDv)###*##e92c]4######|m$ax4######J$&tP$.,####tISk,$8[*,3+r].bMS^#N7H$bMS77S@Z%P?$+D:_,%###3##~#&1?$&]1_(0C6'E%'bMS=j?pISDZ%j$&$R<op9q5%H>#s>#)##46%x/1OG#'$$47*N@)Cv(H12eY#/##$w'<e(K-)&##$c#"},
{297.33,732.09,8.33,1.03,"seOS@,<l$zn$Gf(e=H<c%5##R@$B{?O5$###]5$DZ%hY#?5#HJ1X,#gy7)O2>S'VhQtdRL?&ihRF-PBZ%jm%uf0-.,=5#H5#_P#>5#$L2H$(|,&$d$ihR%eRS[RTc%hR'ihRq&3E-(QG#dH#Bu$.,#=,#)Q$OG#P>#'m$&y5fY#L5#GA-[o1yu';5#zb#{R%"},
{297.33,732.09,8.33,2.60,"?5#]5$Oc%hY#H5#%p0-.,=5#dH#fx2E-(QG#|R%nl';5#zb####S@$M)@O5$jm%ihRF-PBZ%ihRS[RTc%hR'[o1fY#L5#GA-5##Gf(Z4H<c%L?&>S'K_QtdR%eR|,&$d$ihR&y5OG#P>#'m${n$reOH7,<l$4X2=A1X,#gy7H$(_P#>5#$L24Z$Bu$.,#=,#"},
{297.33,732.09,8.33,4.18,"zb#pI%yu';5#HA-[o1fY#K5#'m$&y5OG#P>#=,#)Q$Bu$.,#QG#dH#q&3F-(hR'ihRHRRTc%ihR%eR|,&$d$$L2I$(_P#>5#=5#H5#uf087,BZ%jm%ihRF-PtdRL?&>S'aqQfy7zE2^~1X,#]P#?5#]5$DZ%O5$###Q@$B{?<c%5##;](pFH<l$vn$ueO^I,"},
{589.55,39.04,10.15,4.85,"26$7B(>'7>5#T9GW-(}k#]5#9:V`>#?x05,#~P#5##4w,###b,$?I%%oX`,$dMQBQ%tu&Ve'pqX1##w^5tl#{,'-###D7C5#5c$4Z%_j6=y1zf6###w,$d'I-nX7##dg6U_)T7.P##NYH=,#2d#j^;r>#iQ'e%,{k#(##S/*SsG;##s$,b-$9Q&X##px5$##"},
{125.79,396.29,11.22,4.85,"Fu$.l#|b#NQ%4//-u#$##Z,#@Q&H>#pY#4_3%l#E>#a?#X+G$Z#t-&L(6:l$es;)d(95#:H#7~T9#$?5#ku8}d+<c$`l8s]TU,%26#a]Ti>$1sA8#$pP#iJ%<`TZw/=5#-7$E_-)SJLD<xI)E>#oQ#Gs=Bu$zS*_v)n#&h5$D^T[?)###K~#U29(8.(K*d49"},
{595.57,402.38,10.39,3.36,"@h<#########(;gYG#D>#Sl#W@-rl#cl&^:2Jm*)##]P#jD0sL;#########J;g######gw$rU<$##RG#VS=ZP#Xc#-%+%i78T3#########m<g######G##q:g###@5#?n&.,#C,#Vw(S81Lv)#########G<g######$##4;g######3Q#{k####'##Fv'"},
{136.85,669.23,11.54,5.35,"-$%nn/###&##o6#$<@:5####I6#/|4~_<6#$:|5iGCwA0|m(:,#OU,e%0###&T']_R|:=o,%]_R&YF$q2)K/-z8au$bB,=^R;5#Y%%ZV7_l&<v'Dn%s?H+^Rb[R?l$;u#wQ<z84###5,#n-(eP#95#I7$Hw-RG#1,#5Z#U+Ace,-u#/,#}n('.(~P#D>#-##"},
{560.78,760.82,11.59,3.24,"#T3######=##7;j)##3u#kr,37,tI&Cg53+>&8.WS,'-&4g-7(;######$##,<j;,#Hl$(6$>V=*@$gi8P3@P-)<l#/H%oe<T^:######'##);j6,#|b#*$$itK;v%t$'M_7+A/%c#d5%zS*_K7######*##y:j###%##H[#t)D###N##P14eY#(##SI*5@*"},
{460.69,22.37,13.01,4.88,"w'%dZM95####xZ7ze->u$###uZ$Ve'll'###############?s2L<DG,#I,$p4[WQ&o-*oG#Rx1G?$r`@^G#############[&1&Q%C]#L0[V/[###E?%bF6=)A###2}4X6&######:,####:5#`>#q9+zPOMJYH#$gH%DV1=+H%##iw$j6&######jG####"},
{195.39,399.22,13.11,4.81,"~P#,v#VmNX5$})@RG#|P$&I#~pNSv)dP#G@%uZ%EK+pYHju%F>#k?$LRLzG%m:4q#&66'Sl#HoN16'7,#'p&ew,B.(b?D+B.J>#XJ)YmNF6'Z..=5#'['J<6>pNlZ(9##PS(ho+Nr3v{?Bv'-h:)x-@Q&xn'b&3gP#(o';pNV]0Ru$Q9*sX7w?(vH'~C.5*;"},
{585.52,512.29,12.15,3.24,"Io3######%##Fha######%$#7fTOG####C%#|n0~P#0,#hJ%ZB7######.##&ha{k####w%#s9~L::###/'#fM25EC0,#6v#GT4######%##|ja_#&###n,#?PA^FDo5$]v'o#<s/4ZG#]5$mI.######)##+ga######y6#B<D>c#/-&fx,xI,B5#i,%c9)"},
{300.80,666.69,12.17,5.82,"[5#/R(xY$###o,#xV7a[,###=T-ZZ=1+C(0/JJOMd*K>#XA'P##vC4L-).,#A[%UUOM(9c,%RUOnD>Tn*lM3EaF.,#'##{K(.##3g,Y81[P#]~+/253l=0i:5RO1Q%Ml#Ok5S{>.m(###T##2l#NU,^kG9m(H.O%C4jv({1/;sBcc'###MI$(m$TK6###2,#"},
{466.22,666.32,12.12,6.05,"###[/,eY####&##332hXF{k#ef-]|<&O=;E;cS,%8,+6%BA,###/W-TA2###H@&0WTD=F@#$0WT.}E:d&(14bw.%##s#$TV;###}['``A95#W[*ce'0WTC]-`STJQ&sH%)s+OM9[>$U>#%H$$##W>#hT3OG#(c$1##OE86U3+q:z5%}Y#ig*Q@&e]2+u#'##"},
{566.51,673.36,15.92,3.26,"Jf4######$##o)o###%##;v#Zq>$##'c#zz/m>$Md'|c(>R%&(;######'##{)o######WI#YZR:5#J-&yB.zP#]-&($J9H&Rg9#########a+oD>####5##RxU.A/#Q$rY#qY#*~)Zh4oc'-80#########:-o+u####(##kHH[S1###c##IQ&zQ'J?&A8,"},
{555.29,804.96,15.17,3.44,"#####################S>##########kd)ZP#######4C.2%-######+##U&^###I>#r`)be0yY#7S,cJCTg9VG#vb#9l1|A2######$##%+^###nY#=v#SoZ6,#`&-P}9:x2J#$K,$(4*:6'#########/*^######&##%(^###$##-0)tc(xP#H%,KS("},
{59.74,56.00,18.37,2.89,"'##0SQ[P####$Z#eiUeY####/}6L5HE%.5##yu%+7$+jE(##$##ci4)a6`$+7/-e^L$)0$A.:gUS05v?'R])Y7-{G%wI+al$$c#Gu$lg'Ao3*v'(H$*j+)eUoISh36.r87dJdl%Vw+[>$2,#######$##OG####&##Uc$p?*###j-#ye0_5%###:##+u####"},
{59.74,56.00,18.37,4.74,"QG#OG#&##SJ/?l$kc&c?(Ru####lP#4~(AS0######.##-H&*l#@a8,q83sAZe*{NWDp4'v$&LWQ4BOy-ykH_?)###6##6~-###6@&P^JEZ&m..jy0<d:o]/tJWoZ'<$%n{)3&1######:########mY7'6&z,'###nF0a8.DXHD5#8d%nV-`5%6,#^>$8,#"},
{429.11,186.58,16.38,3.08,".,#xY####lG$)##jx/.,#sP$u5#X95)c$SG#gu$Zc%@u$2c#~P#WH$Gw+@v'TQ$|nB[OE^>$D]E)PG;l$W,#0p4PG#.,#5I#xo1Y'0.;6)C0>I)-F-m_V;R(D]VV6'_v(3f$]+JD>####k-#V$&Vn@i%02,#?K0)[:nS21c#n_V#v&ZP#M##m{4VZ'###-##"},
{429.11,186.58,16.38,6.17,"###2##a;4o#'D>#W##0VV^u%P..Su#ZC6a[=^e0:5#*v$(fC###F7#Q|E.,#]-(}o$|SV.v&#VVDd'TI*hX/`M8tK/?'2d_4:5#Ud#Q]2~P#/c$f5#FoD0GH-bFiG$^Z$(fC$A,&7)hY#~Q$K#%)c#<Q$2-&X>$iY#m,#tT5.,#qP$)##~o/###J,$:5#3l#"},
{129.93,353.10,19.52,4.76,"b5%###4##&13#?&###%##X=@######V5#_GMoP$$##L>#$_9uA2P#$Sj=fI+oIV(l#oH%?KFj#&%Z#tJDS,KM5$$##%7(=H&b96,v#PJV~#$?MVB6&ov)tn&Df/Oe&mKVn?'SG#5,#zJ2/Q%om+Wd'7=DX(5sJVdl%K6&#*6$x-=8'*KVyZ&Ql%E5#-A/ZG#"},
{279.83,350.27,18.48,1.61,"Xl%X>$D5#h,${Y$###@c#pI+######Aw%^w/######q,#v@/y6+###a>#8m%YPKV5$rl#{fI<$'e?&K{(FJS=u#[c$/b8?ISn@-{b#<Z%@5#=MSFe,p#%1x'89+f0LOa=tS.X@(9MSgNA>7*<[*gP#<Z%A5#nISY,$-d'An%Z7,>/(bJSpQ&tq<*9/S'7,Z#"},
{279.83,350.27,18.48,3.38,"}]6###%##N?$juJ$H$5o+hdEr>$J.*=UMJ2>~P#>@*ol$P~P4q:######_H#nxP-['Z.*~'*X~,#U.Q7B/lB8.).O?I6%exP@n,{k####JH#8pK_tHV5$`l#g[&GzP$S,i>$3?$owPtb#L>#zY#z,&)c$2,#%R%4%-.,#/##@#$OR*TG#l>$QG#?#$[G#kZ'"},
{271.14,421.88,18.14,3.31,"1,#6,#?Z%:#$h6)56%yH'ww/xIPe..eG#511-~)*GCwc%EL4k@..,#?5#:6$$<AWG#;R%eLP<8,~f1z%E&>I:c#EZI?7)N&1SB6######BH#RJPz,%_$(R9)(n*|?%fLPmi<8Q%GZ%]I'_JPz-+######/v#6KP1T2K#$A.%CJ*wdEz]/F7),@&SIPf5$g@+"},
{271.14,421.88,18.14,4.91,"7y1|d(s>Fxl%`T3uP#WkDsm(j8Phu%UQ&5.'07)%@&M7P(Q$xy0y7Ptn/{Y#c}H(8+af1ASD`W<zR+%I%N:P9.)ZJ+O7E`K/zw/]d*??%T-'N:Pw2AWG#8R%(')W8Pw,%Dm')%%u8P/x0&c#D,$%##8,#bu%MH$K%..,#L>#}5#K96######xc#z?*######"},
{271.14,421.88,18.14,6.00,"GK%YuNyl$l[,L^.ZN@~o)8z7zn-;S('*5?[N###]H#hE<8Q&i0-0`Av##Wo3uxN>C9|>%Bg'8y0q31SuNU@'###WU(JT5###w$):n(;a5XuNn7DmbGJB2Ym'*e%_8HR2A/l####MI$Rd+###?m#,>@VA1'l#0[%K8+qS1(v&:R'6.+95#wP####o>#ic&95#"},
{256.84,481.35,19.54,3.26,"D>#&##,m%'v'zG%&Q$<u#(@(7l$_P#Vu$`'495####,Q#LtD3,#I5#MQ&{Y$US/Q['gm(Se.x7QaJ0^u#yvB1/+-tBD@&*9Qhw/95#pY#T$%T;?~G#Dx(]:QjM7%tC*eCn<DhQ$08Q*$&in.1C9######lZ#p7QL?%uA./(,WQ&e6%C:Q*#HCl$fu%Cv$E8Q"},
{256.84,481.35,19.54,4.84,"{}E95####|G#e/PIA+JbDHI&ue/I?$pwPB6&.xPOu$#-&K$%Xg4Bu$TG#Wu$T[B`wPAx0g#$(=Ek_6djBb@CLYF+v&c6%&zPn-(*Q%&Q$0l#I~.me/Q['Zd(>zP`D?[G#jS(lp+YwPG?%ZA-I6(D>#&##8v%{Y$2,#I5#MQ&Cv$2J/95#pY#ZQ#(L8######"},
{126.43,490.71,18.81,1.49,"N?'R11.RSJ>#ymSsA/A$'_U0_R++@$R&OkV<4,#9Z$5/0~P#=e.T%${mST>#inSyZ&;7)B:(pw-_]'XpS.n('##0[%kC;95#0I*h/(U_=}b#7rS7-O:Q%&-$HI&7rS.aC/l#V#$_R'X~/(Q%0##Kz54l$###`##C`=.,####/$%Oq8###(##MR(h?)###B5#"},
{201.68,488.99,17.88,1.58,"Y6%=C3R^0=<@2,#cJ'deQ.n-q?'M?&yu&VR+X5$]G#7$'-u#|v'6gQ_z<1,#Q+H2_2Ky3_uD,-')Q#6RA<$PTG#Pu$V7,hY#K~/,[$ZdQB5#`eQ&v%/['3o'4e+:S%gfQsv(&##MZ$fg8.,#xd,av#t>NK>#xfQ'f0uu%_R%U@)2t7_eQ2$&(##8m$,L7{k#"},
{576.11,59.13,20.14,4.90,"2H%`-'>b2Y~Tu,O?#$K)(<`T|uQ.,####03,95#######`P#^l$bM0W3=XZ&QM7A/&O[Lsx0U_TpG$$l#@w$?-(######9,#'##'/FBv)###KA/vA+6C6xJ*W[T1,#dP#,D&#[)######>##%##lq<######.Z$SQ'###?,#CZ&######A##.,##########"},
{132.58,425.95,26.47,4.82,"######R##pI,ZP####)##rV9L5$###'##6#C######8##%GH<K4;5#al#?H?ce05,#(5;M7MF$S/,#qw*>[=Yu%%##Bt9Xq:'&S2c$v>%.I%,|>_m'5&S^d(k%SYG#Tq:IJ(m~0bG#G%S^5$?p0em+SG#O>#:O=Nh:9A/K~)#s?P]4Km(4P7M?'Ow)K$Nbl$"},
{466.30,424.62,27.01,3.54,"ia<@c%^P#Oc#AA-^f3$##~Q#-w&u@/gP#U5$Ky6UG#7?%4)39pJ:%-<5#hG#hx)%/L4]/Pl$xzR>^7-Z$Ic#+xRD>#*##M.&~wRE>#8,#:~&j/3Eu#cr,'7LAvRSG#iH$N-<OHJ-$(.,#_$#swR######%Z#xbJ###L##ez8g,&-Z#C#$3`6}u%^A-L5$P,#"}
};
siftPoint_save sps8[] = {
{300.12,70.02,3.49,4.22,"fq9######'##TJZ######A$#AM?######<0%pb####$##bT.bz91u####%##JKZ######l$#&+I######cg#_5%######sz'i/3/,####j#%[LZ######g##cJZ######h&#7m)###SG#W(%|Z)######AK%*JZ######@F,zIZ######-{#8m)###%l#vZ#"},
{113.76,74.14,3.75,6.27,"Vh.L#%SS+_l%=h.b{2'29-c$+fF>i;Ax-3w(5/ROG#L>#y5#hI,###X3.bJ/..R)##Y_-+16MV>%##%B&|1Rw|H%##G>#RL''.)E,#P?Kn5%#.R(##j)?Dd&tOJ###.##k2:#-'###%##Rv%^Q&Am#/L9>,#p,RH##.]3(T#mC=######5R#|k#$##H>#Jl#"},
{375.28,83.54,3.85,5.82,"wP$R8'w6+/##N96s>$95#^'#9FF######iq##q8D>####:;$2I+X5#C%.A##KwW.,#95#m'#t{W8Q&###)'#i|,:$)###l5#&~.9##1m'},%ezW8Q&###e,#[a-oDB###&##i|,l':###-##YZ'###-f(5d(lC9eY#G>#$H#b00`$+###Ic#PkC)R*###A0'"},
{423.51,92.92,4.20,1.46,"1u#Dm%M5$eG#(H#6a8P,$wb#XV(2B4~G#zG$c%*######Jc#_>$T.$w:>3,#&;=Q9)F^3YR')6;iQ(K<4C?&&/F###^P#$##AH&g?#<#P$##R$Poc$MM=Bc#3E=M-&='PyY#w&PaP#Hd)$##C#$sw#2YL$##,k9z,Kv]7###jJ*p&P~J1###jZ=#V;v5&###"},
{423.51,92.92,4.20,4.16,"####%#?,K{k#$x1q$#URPS:3@H';##'41/[PBv),##9'+X*4pb#2$#3ZP###>ZPrc#+ZPR$$Y97ay#yZPjR'F'8k##df3r:&ZP#L#$n[P.,#yFFg$')?F+@*_.-f]&>T4iM,-S/I5#6Z%pg#,##{?)FS+D>#A#$+l#19)cV>O5$L>#fP#2==mP$lP#3c$&-$"},
{409.05,115.48,3.31,0.46,"Hl8du&###B##}8'P=:>u$Ic#4[*C9$SH(P^#fS2+##95#3(#NHB(c$###L##h7PQe-###0(&s[-%Z####s;*0Z%######-((THMZP####W##q9R######F-#C9R######8w#(I*######E]&uy9.,####iK%?kC######+I#-;R######-##J|@######x?#"},
{409.05,115.48,3.31,0.93,"<B0`8.7.-?-#-^8ew$AtK}##rlS'##W>$~%#x%1######H$#yh(tmSmP$###FnS/p/[?)'&#nmS.,####}&#QT5######*%#^E4?d*###(##YqPQu%###@##0>4an0###:##S,7CZ&###D##7&*Z[,###&##tfPbG$###8##JO0n?*###8##,2&<R+###$##"},
{446.09,148.02,3.89,1.60,"#95######a##JIX######F$#X(i######J.#<n.######M@'6M>######M##_(i######L##o(i######]##?07######aH&OV=######$##P)i#########o)i#########Kh=######0,#-'4#########R[Q#########N*i#########<q;qb#######"},
{486.00,148.05,3.85,1.66,"o#'######<##*h<######E##`8_######gu#Nc&######b.--^8######'$#(il######C$#Cil######G$$O]5######oy4<h=######f##;il######X##@kl######4##U3D###0,#uu%x]7######V##UlQ######1$#:jl######{##?XF######$[$"},
{406.63,151.67,3.53,0.57,"00.######2##&}:######2##XA0######In#ZP#######2}-)d'OG####hG#,JL.,####6##p8[######7i$37,######Q,3Uc$t?(@c%%##-^ZAu$###K##|=[######]$#0%<######Vv#H#$:u#um'Ac%Zs<QG#H5##$&fZ6######$##Va+#########"},
{307.29,154.71,3.90,0.32,"hH)####I$`.(MQ'######T2817Rpb#$##]K,z/P0M<###7##x7/###rW-&Q$R5Mtb#rP#r{1sx1o#%8w&ikB@4@iH)R>#j@%VQ'###lIArP$H9R(d(wA*&d$K92(T.YQ&6{,5p71,#K>#~E-<d'###dG7`>$9N>.,#`F2OQ$R?Ak-+C5#p:,RB-,o2'##*B)"},
{312.16,164.13,3.64,1.13,"<N,u&495#Fc#x'.__;PG#Nl$hpOwu'###|,#kt9eYC|[,-e&b`?~>$VG#cf'6(6cZ%y@,Do-JnO'l#qb#=-#_h5d06###}R%Va1TJ1V>#b$*0X6.3>pH)&?$mnO-{<OG#[$#[$(r]K###*##s>#0Z%@%#-lOUS&M$*V$#6lOc8'>+Io##o{C3##I|?~,#3I+"},
{93.44,257.55,3.85,2.65,"<d(ec%OG####;%(EL5.,####>~#I06}Y$###Fc$$?%u>%J,$;Q%U$%@w-'H%je(l|WiSY+Q$fXY]TYnl&'R$>L1'd(###/##|>$Mu#-C6~>$a?(Lw%fXYQeQ9TYpu%W7'>zKxT5zY$.,#76#(-&dP#oG$E,$Hu#T5$s?$6(9+l#5c$~5$iU3c5$KZ%1Z%F5#"},
{93.44,257.55,3.85,6.13,"<?'T>#;5#p>$B,#{C3}G$Pl%>.&8B3dP#'l#)Z$$##C?$Jl%-c$<?##040Q%6A'*;NZTYcl%fXY8-LlQ'C.&c]2l,%<c$l>$###-##V^2Yc&&v&A-$fXYGTY$TYb5$88'fXYa./H>#fl$gm%95#2,#@$%ml&Fc%E>#O-#U(;QG####`Q#Ji=+u####gQ$,m'"},
{198.21,263.67,3.67,2.50,"rc&d#&###(Z#wm'Y19###'##5R$;2;bG$%##4,#+q)6?'&##S5#C%(s[-###}I&'jV1eVJl$'jV?fVY6(Pv$,f/'I$|L3j-)0,#Tl#=^4Q#%'v&~7%'jV>fV|eVq,%:I&cqO4z8:5##I$[d(c#%tH(~l$`,%Jc%.,#eQ#[_=fG$0,#,Z#|04WH%~P#QG#cG#"},
{198.21,263.67,3.67,6.11,"Ru$e,%nY#fP#6##:;5b5$95#k[%q06|b#F#$Nl$>c$~,%<c$0,#XQ#iB3<c%i7'oCM$0X^l%t3X#vK5$'|v%C0495#&##ql%zl#2C8)L1bG$BR({-'t3X^/X6/XDl$8~&t3X5e.(##],$~J*`x/g.-6J/$H#Nq:L5$n,#jU3-H&###)?$:N;###%##&@&?-("},
{165.67,320.48,3.80,3.79,"W0$7|@1@+s>%I))&(667+Bl$oe#:[S[c$5l$@5#R>>@C3^P#jS'MS*rx14m(y~Mf$*-v$VU6dm(pz=F&#3[SQG#+q3l$$P[SWH';-$Nu78d(u[SH>#.$$Tp*Z}KfY#:$#SlK$##jY#uH#}ZS$Z$u7*vI*}A*M3F.,#jG#a_&d]4+6'###q##0##$v'K>#D>#"},
{247.06,320.58,3.66,5.47,"yR(+oUH>#1Q#Mi()aF2,#*Z#OG4},'?5#-c$fl%E>#W>#gQ(-J/)Q$aQ$CO.Gh3)l#||9QU0&rUE>#+6&I,#OT4D>#@5#:l#SV?rG#vQ(H)'H&3pH#4oUn.'OmU$##/.+^'$`/4###/c#&-$#3?kn'hm+>R%p[-4A#<B5C1'|{D$##D,$>_#g,&###d>#T6&"},
{174.67,325.23,3.94,3.69,"dS$vK8/H#bg90_%nV<L$)xZ(4r)TT0Sm):c$R8$qsGOZ$[>$:5#[5#v;+6)@qT-B/(aq6xZ(a~PiZ'[6$rU6zZ%C^8^/#p5Pz>#f{;Ic:@H';%-.R%yc<`I(EZP%##=Q#Er*I&3qb#i$#{WBxP#QjAW-%:z:|G%1['~Q&|}67^8fY#5,#8V&I$(p#'$##9##"},
{488.08,328.52,4.22,2.24,"############xm+######$##i6N######1##G~X######q,#4A-######3##-#L######v##M~X######s##g[X######I$#=gH######&##R'T######x##[]X#######$#._X######;##8k8###V,#CQ&ym<######$##+J?######$##Ct/W>$######"},
{347.59,332.80,3.61,3.76,",JJE$&ou&g%#_:U######9(#>/0p,%^P#BH####z#$GZ&.,#@-7%d(I*77.*d+[D>#$##td$gr<eZ%+u#4##0,#iH$A,$###Km)###yr'g([3&[###;,#OMNvJ3H##CH'3d#[P#I,#3Q%F>#OG#(##6-#jWC###T>#2H%8M795#/##.I*;u#95#C5#$H%SG#"},
{100.22,335.45,3.67,3.19,"D'DL5$######&A)wK7###$##sH%0y6######n?'r>%OG#$##.6_-Z$O,$>#$YV4`cFf..B5#,9PBC91u#Y5#?/0Fc%eY#&##m2_WZ%'Z$tH#5lLj,$oy.Iq5IIUP#%t,$'V.Y7+Qv)A,$)##Z/_H>#QG#5W#x`DB,$;##|S,=H&FH'=,#(w)@l#.$(UG#/,#"},
{112.18,337.11,3.70,2.97,"q4:t>%0x(f#&FH_.,#(##_P#j_7Qu%###+##FZ$P6)###$##n~05v&.`1#8/nD_|b#9H$>H%jlI~190Z%$##QT,M0795#*##Cx1b13Cp7_Z#2B_56&`>$*U#d2Av6'of/e%*)*ASZ&{G$P~&L@-uG#dd+<u2T@Y.,#$##U`$(w,D>#DH#i/.3$(|k#?,#B7("},
{165.95,347.40,3.78,2.62,"mw)A['Q_:.,#qrUa>$$Z$0,#@;9u6):5####],#BK2Y7,###Lv)Pu#OuC[,$)oUHc%e?&j6$h*BX-;9#$&##*$#=qU>l$###].-@q6Dn,jl#znU2'5Z#%iv#-p3vO:w7(Nq3Q,#'C1342^7.1$(&-&Wv(HG9_V>BZ%$K1GS(=.Tlu%B[([A(`P#Jc#.9So,&"},
{120.08,356.31,4.03,2.52,"n'$HdSf#%D>#Uc4*~.b##=~/Ve,###6##6<>T,%###%##f~+:8+0S+O-M>5#>>~au%;Z%P5#0&WeY####9##lm*jY#F>#V>#>I+6v$z#G]Z$38~$##2v%m9$4~XOG####+$#l?(;#$.,#:##2x/ZJ.<-'Zo$+/~2H%vP$hM#L:2cR,D>#S##j$&W>$###*##"},
{194.91,370.93,4.03,1.42,"bA'I7+OR*$.)[1,vu&95#7l#Kl$Bh'Su%$##j,&DO2n,%86#J~'1r.m./#l#=COxu%0v(4##/I)`U$P>O>,#6*AW@%>6(HB#)8.[B/[P#=0'-?OY>$*w*?:$~6(j,$8BO1Z#@F28d),e)}Y#Fu$D5#F>#=-=Zl&$##^#$_q'%##=5#&|-Av(>?#0Z%<(%NtJ"},
{194.91,370.93,4.03,4.60,"xg$FlNx5#uG%Rr+,R)/,#2,#6c#z:'_Q(###QG#Sl:{G%T>#[w)[5$aN/u-+wTOHl#+[(_,$,.*+($OQOoP$qb#)f%MS.k9/`Q(.B#94DUd$jPO>,#i-)6C$sQ)/##f'N3-&un08#$s%'-2.s5%Q6#OQ'V307?'$##wG$o1'D>#}Y#*(+5-'+I)rv([8':%,"},
{242.19,369.99,3.72,1.41,"l.$go3uu&?u$JD,xc(###%##gI+qm(.,#p>#x~2'6%9#$63*,f,T',}c([P#KBVlG$I#%5##KAVi['Gw.k##sBV4H$ox,#z%Le,d-'-u#/:)v,P$##[Z&4;%<r4eu$Q>D9c#6TE*NB)n&O#%^>$###;5#^Z9Su%/,#W>#;L(7,#Du#|:*;Z%H%#Io2lU'll'"},
{126.36,379.39,3.55,3.82,"/x%sP$q-;jZ(`j5[P#_>$Q,$S@+.,####>c#A#$D,$%##?5#WI*D,$7#]SG#Ez~95#R]1]>#&h8{k####&##xP#&?&$##1,#*w,oz#my~{,$'HQ&Z#F1:<f$mS0In-###>##T##n./######{k#?3'FU9zm'4-(F9*r#'[.$%c#X^T6#$+##w$#zRY######"},
{222.92,389.25,3.40,3.73,"b[#Xe/..%UH(WW4sI.9#$K>#}R'8/Y5,#(l#4,#E6H][,$##Q?%+u#(iPn?*PW<-6'Vw*?I)BH&L1:s@#0r?$##bp6A?#GNC`w-G##64Y>5#h.Y'##-6IU?%~1=###5$#,;<###$##56#dkMYm*,3(n.Y8l#8D?'?#YrA)8%(^6gu&%##7Z#*##*-&2u#xY$"},
{153.38,394.31,3.76,4.00,"76%i~,q6*%s92U%P&3{Q#LcOUH#Cv)f+5K'8'w'fY#u{395#l@*KH%zE3e$*dg'GA1)x$<e-r#BJ.,JGBE6%*[)[##BKPnP#6V&%(;]R#J-)L$%M}Hd'%H{>vS2+c8s~1^]-k>%lf#ijH7?$6V&FT5=#$###0S&IHPG>#PG####wW9wP$Sc&###2?#3m(=$("},
{153.38,394.31,3.76,6.01,"Du$58+D>#z~&yY#Iv)H##%pR[Z%B,$R.#$|CKl$/,#Z,#}Y$[U;s>#Tc&gz%ZH&(a6F%,m?@DB.cS0J`*FK5=#$_5#t<5D,$*nR3Q#O07%&#L%-W8'rlR7Q#XmRcl%by4#7$Yl&K$#SnRgP#CyJ]#$E%-.,#y@++v$;2:P/-{;C<e%dn0Y~%cG$AJ#alR%##"},
{232.47,392.55,4.03,0.27,"R5#/c$c5$v>%dG$s@%U7,uH)9%-ce+,u#,s/_,$y>%CQ#S7OB8#E/1<c%}k#8q:;B'*WBI##/JZ6$$^J2W<)BJ+R^.#|ANh6Nm#@?'~@%;1;Tr>?d&:z7S5$wOZO,$Fg6:##Y:5q#%k'Y,-&I5#'w(P91'.,P?'Y5#{90vl&L>2m-+j%&X,%|M%OrBm.(B[)"},
{120.10,404.64,3.49,3.02,"j8,.,#gUHub#nAJ.,#l#%/,#{[)V$*######)?#VH(######rm++$#Mc^tG$;^^tY#?h7tY#lA0}%-Gl%###^8,U-)95#%##~v)uD#?]^5##K8]2d#5]3X6#x@/nP#Q[(%7)QA1H>#fP#B@'+6'r2$v(?;~%[q>6##A,$kK#=6(###1##K7)pP$###&##s-("},
{153.64,411.36,3.84,0.21,"q>Qo,#@.,yA$8J*h'/_'5{Q'z<;u5&*V,S-)1u#(##EyN]P#?CQ$Z$6].o>#y7+q5%P9Q/H%DAQBl#Fx0g?%R6(4?#=AQ&##CW+x$,Tv&95#S:$LsFL^0rR,(f0y?GZ?'>0-.Q%oC*<h;4,#1$#j>%@S#?d*Zo#0iAT##>80WH%*QK###?['rc'q$)$##vG#"},
{153.64,411.36,3.84,2.96,"|6)9,#<#$I@(###~Z%Q$%'6N>##'8/qK%%g7`w$n-+`Q$D>#utM6##]#&W_'YH&a1/|w.;YAqx.{[,b<.*]2`.)D>#[?<SG#9$RF,#Vo4Td#n%/9n%s#RVu#]$Rj#%*q3)I$d-*Z##E(R2l#.7F8u#LR+9##Sn)UZ%rW=^A.8WAO/%q7/Cf()-'$1#B,O0##"},
{225.72,417.14,3.39,2.52,"rQ&/c$Ut[,u#(9Q0u#n07###TV;yb#~P#<5#/,#<c$K7+.,#Z..C-#'p[QQ#9o[2l#56Rb##A*BLC4###(##?##EGE`P####-H&a;$[n[{>#[=H|6'U]5TH#.@)T27^,$d%-@##P14Y.(+$(3l$bD'p,&;e(4p6iP#H,$lw%^m*ub#>$%X=<###2,#WO7h$+"},
{344.69,427.81,3.67,3.12,"7v&I8X###6Z#H5#(:X###wb#[>#~9X###|c&S@&K7X###v?&ll%>SS.,#%##Iz9$OC###^p,E$)DA0###V;Qmu%;c%)##%oG6.*v-'|,',##D9Xqb####Em#Q7X######AT(X-*###>,#Dn(u$,+##Zm(SZ$B7X######].#h?U######D$#9m)5,#<5#,Z#"},
{250.93,436.38,3.73,4.74,"{0%4'7<,#~#$G0&]?'W>$>,#C5#/M$je1J5#pb#c_)h,&~2)IT*2-(L8$k7.-0M)H$]6)Xu$G&2S7#^vQJ~%bPOrc#d[+#`'F6($R#ytAy?),[T~G#56&EI$z`D###Lq0j@(P_TLu$}d(%x&'Q%d-#-(:y]-Lx3(##aP#ug&{l'(##p-%Nc$k2,$w,<c#]>$"},
{416.50,445.79,3.76,1.31,">Z$:6$GU~A,$wW~|Q(`o4$##e6G3'7######+$#s|H######cc'Yw$FS~%##_S~],#/1;o##hnX6#$###E##Pu#=?'95#)##r[+5@C1jF8##QS~V##+%-H$#pHR######=##5Q%WG#*c$P>#?r<U.-{k#J$#HS~(##D>#z&#/h;######%6#BZ$[P#:5#:5#"},
{131.46,454.89,3.95,3.54,"7W0)l#1L+Ey7fXYUG#Qu$G,$&T0b5$%Q$0e(6I*pb#C5#A$$wm*C/$ZTYUu%'WY.Q#fm+c5#d[QW>$&##P?%uQ&Qu%&##X5#eQ(`y*V;@)0+xRY,##H6(x2'2>Kpb####{##V-'Vc&###)##RG#Uc%P>#WUJ3R*;D:$##,=,BJ(8)@###O##OQ#;@+###$##"},
{127.95,468.53,3.87,3.22,"P>2;c%d7#*.,}rSP#%7,#]P#?A(A&1'l#jG$48(r#'###>5#z..F-#gJAPw.:rZZ>#g$(tG$j`A~c%Xu$Km(7B1}Y$%##9$%+w,g)&TsFg5#bnZHQ#hu&c8#,tII>#,c#b/,TQ&PG#M>#0J(iH)m6%v6*3G9ovW'##ZP#Yi&xx395####e##xc&|k#.,#-##"},
{161.96,505.91,3.58,1.39,"z0%9N>fp7C5#ym8qm,###)##S8+b5%######(Z$:#$/,#iY#ap1Fz(^+L###5>`N#%95#/##4IQ6#$###)##:v'/u#######MU6.h),d)|d#):`9#$###m1#dPN######OQ#eu$6#$###`P#LQ'?Q$.,#R9Cip:######Q4+TH(###4,#|e*]P####I5#[H'"},
{332.50,632.65,4.54,4.62,"###%##IW~.,#>)A$##^&1SH#LS~######|&#8x2######C6####^$#<T~.,#}YN4##[L:b#$CU~######b##Py5######6#####[##:V~Yl&DL:'##6/-}h6'T~######rS$+:8######d##.,####rO1v]7g,&(##4~%/U~q(;)Z$M,$|0-BU495####+##"},
{206.22,671.45,4.24,1.69,"###R-#y%[###eZ'5I'sJ2/##P&[yY$0,#:%#QXH5,#TG#3?####.9#%&[###cg9G%(kA3E##$&[fY####_'#,:9Bc%%l#D6####QD#X5P.,#<|Byu#gu&)##(+[eY####9##N11vc(.,#%######$#/|?OG#3n*[>#?Q%`,%3*[######b>#=5C#########"},
{207.36,682.36,4.58,1.79,"2Q$}M8qo5$##-E[&800,#G##/eUB5#zb#8c####o,$)6%ZP#tH)}Q%tXH4,#TA[J,$>5#M/#(/[1Q%a>$*R#$##I-%|&3D>#uI-&J$H1<###lF[/u#.,#F##}[CUR,.,#/##/6%f5%]%+:5#Tu$g5#jI,|k#>NSI>#QG#wY#ssWOG####%##z7*95#qG$.,#"},
{181.92,683.43,3.28,1.48,"Y,LK>#D>#A##;=_######~##]5B-H&###%##Cv&e5%M#%$##jFGA5#ZG#Zd%_;_######+%#[:_######1##8T21,#0u#2,#9K5###q##zP>*9_######d]#2:_######1##+U6/,#F>#(##@H'###rm#,V6d8_######7/%k8_######Q##9A0.,####V5#"},
{434.74,703.03,3.41,3.84,"3l$###A##Oi5u%X###W##<8)NJX######r$#Uo4######W$#qQ)###fn#qTGFFD###V<*<f21LX###lG#ju#sMAWG#pP$MH#V04MA0u?&kw&C@,(S*]JENI'[IX###^6&p_%*L9'##w$)<%%jl%'V8S#%}5%LJ0$x.EQ&H:$}_A###UG#ch#Uv*###gu#H-$"},
{108.31,715.10,4.26,6.16,"3R*[##]yd###cyd(##:Q&1?#iB8'##{#%FI)cG$###t>#?A0E.-?Q#hyd###&zdWG#8?'+##@3E'##sZ&}P$O$*###9,#Yl$-n-<##vyd$##9zd&##ol'1##0=G7,#Y#%M>#T[+O>#/u#wP#BR+9##jyd&##$zd&##eu&A##2;>MQ$%Q%1##vZ&OT2eY####"},
{491.23,753.32,4.37,6.10,"0x&5?'###N>#dj?qQ)###|Y#j{j.,####n##</2######.##c/&Gw.###D5#dJNGf4###j,$^|jpb####F##an0######(##I]21Z%###CR$6dSQu%###Vs/'{j######s~$'$(#########s09######kc#T:h######CL$-/]######XK#.,##########"},
{158.35,782.98,3.95,5.19,"2a.l^;%##e>$-'+8v'/9,GI+E[$:&/P_8kG$vj1]H(K,$rv&-,6v6+h##QaBdzQP$)R@)hg.vH(h#%wd@^'6`F3#H%59(E&2880iY#Nx#ewQ^xQ+6&^?%&h0If,Ie)E4A;@*9R((d#LdM.e-_&+i'8=d%l6)C&+(.+6v%z,&3y.~~*<$){Y#([)Do#$YKg>#"},
{479.15,801.87,3.96,1.12,"}B$(y6######5},7m)######[Q6#########N,$#########}x+WZ'######c/BG[+###0##it[95####J##Kw-######'##&o/F>#.,#$##?o[bG$###W%#yn[######+(#7.-######@##y@.######0##]n[######h'#Tn[######Z'#>u$######/##"},
{338.41,807.30,4.18,4.43,"Wc%@D&06S###%TY:?%eY#'##AHP######+##Zv)######(##Fw,u''hRY6##|SY<l$###n##QnVOG####@##|R)OG#######yI.>?#+nQqw)cSY###'##pd#FuG+u####,##Rn){k#######X?&Q/2=j*#SYfIT/?%.o&*j@4i6eY####(##h%+{k####$##"},
{263.93,825.57,4.09,4.54,"a7-Q?#/0a$##g1aJ5#jZ(%##vV>pY#.,#A5#QG#;l#OG#c>$WJ1D$#=0aUG#21a/##/d)Y>#NuNE>####Z5#U5$oP$###Ol$nl'7##z5a-;<-0a$###v$'H?SbM.,#$##AQ#Q5$B,$(##nG$PG####.+4y_@j7/$##px'NjB]08/,#`P#}P#m>%/,#0,#{Y#"},
{307.79,825.58,4.00,2.54,"R))}T8$##|P$x)+mq?######]1&e*H###'##V[$H?(###&##c$(95#s$#8ZOoY6}I/&$#=T3))S@;@8,#tV+b]-y7+Zd)zT/######d%#|+NLR+95#&'#E$S2)S#A/u5#^r:gS'i@Db%.iQ'######K##}Z)######h##;f3.Q$xY$[##nK7M5#IR)=[$/7,"},
{49.58,49.15,4.79,4.03,"K##?:3Du$D>#F)#,dOxY$###-`#88]q,#I#%Q2#r3HnB#~./c@(4h6{c([G#W~'9/-V182u#XfBYd+>u#uG#PY/QWD5,#D##7n,'Q$SH&kn,%8.Hd)Tp+~v'NmTkY#<v$o9&&o0eY####F7$<5#AH${Z'U$*6.,^I+?#$*x*1iA]P#G>#)g#-H&######ZK("},
{399.10,102.40,4.15,6.20,"zy*()@######Uy)h5Q###)##f)4WjH###VJ)@H#LV?+##I(9r]-OR,###*##xvO6.-###X##q8QnZ(###B0*^$('D6g5$_;=Ng6###~,#N?&:,N######N9#r7Q`?)###-41/?$6L5###YmI0H&###@$#iHM^(>######0F0e5Q######vr,Nc&######_lF"},
{471.74,147.65,4.98,1.65,"Ve0######_##gKc######w##QCg######fm%-H&###$##UT3Y'9######s##dCg######l##qDg######yG#}@/###+##j-)?z<######V##SCg######C$#{Cg######r-#;x2######Un)4^8######>##$B^######9##ZDg######?##606######w,%"},
{390.90,158.40,4.71,0.24,"######3##%z8f5%1u#;5#{7+k@VD>####FH#XaG######4X/###$##K?%xu'Om(&I'[S.>#$t9_[,%TG#^A#H:_######BZ2######@Q$7#$PQ'=5#_/*&@)a;_.,#Z>#47':?_######P,#######zY#R5$OG####:l#5Q%CM9###$##[,$:?_#########"},
{298.27,178.63,4.56,2.89,"&###[$XB^###Nv)1u#)s:.Z#/B^OG#S>#1J#'2?A,#B,$^##3##dK0PxN###wA4?Z$M9JI,#oA^gY#u,&z$#~p8=R(W>$;##G>#sAR_6(/u#@C9)=4hXFE>#DC^T$'NQ'+##0/,osA@H'###zQ%7y0zH)@l$5;>g&(21:kY#KB^$Z#J?(3##=m'Xv'Rg7Jl%"},
{317.72,190.46,4.63,5.95,"&?&Z##g?TwG$|D>)##WR,P^/wH(E>#-##oQL,n,&l#_G#z5?Yl%#8#O?T###88ID,#VZ'fP#VAT###T##JC0X90#.+0[##|<9l#ER#?ATWZ's@TN>#O#$(@'i?T###0##ho&LR+=#$yc#{^9<,#2Z%i;*+WBAU995#}$#rBTx?TT,%###S9&mQ&W[+:5#$##"},
{474.23,285.94,4.62,4.57,"95#{##=IW(##UGO1##h/4X7#*#O######^.##,N######X@#cG$2$#IIWS>#EIW)##[..Ne#*7T######y##fIW######s##$-&M&$zcT&##]KWZ#$8Q&J##_YA######*##)JW######x##<?'PR#^g:'##+dRhP#6#$P##E<A######/##QIW######u##"},
{198.21,318.07,4.31,5.36,"S8*z6Sxb#b#$W;(kXJ$##nG#o_,Wu%p,%CQ&rP$$##jP#dS.2v'kQ'PQ$7%@%M0>u$?_05{5y:SOG#s5%tG#Y82&##I>#16&|MB:##87+SE,w@/4$#]8S+f*w6S=$)3@*Y7##'1T%/:5#T>#raCt.(xI.km%8$(%U'MsFpG#by3v7SW,%`,#(U%>6S###,##"},
{216.88,351.51,4.17,2.52,"<e%`D=b'695#`{Twm,PG#zP$FD9Cu$G?'~-(4,#L,$}%P}k#+R)uu#@eSMQ%WxT:l$$H$O.%?#H9}Aqb#*##J##1wR`x-D>#@%,2~'S%-8R%|wTku&.##)&%Hp3K*9p9/Rw*%$#_nQ$h5eY#pu&xY#iH'-,=g4J?5#y?)Qw&eaGIH$^N>N.(N>#)-$(xTPG#"},
{326.28,365.21,4.55,4.37,"$##}Z'lZ'95#mP##q,jv+###{3By@+OG#M##y8495#'##kl#2##g_9Gv)###DS,T=,#ZP$##}U]_%*r5&-$#jT7###`,$:[%6&#NT]8Q&###v)@Z31;J~5##,T]{P#:$)H$#/C9@5#6?&Hl#zE5tMAI>#]5$'QH_A+:]3@,#MX]0/1bG$9##+(,e7.yY$$##"},
{138.14,369.68,4.33,1.57,"G~&rS0k-*1Z$~D-Sd*.,#&##{vQ}Q'PG#r$#;PK8&,;/-2z&Ln)fK(6e.=5#oxQNu$r>%H##dzQ*d're/8##%KG$WA<&)-~-^m*YR%)Q%NU*xuQPG#=u#@q$}@)#-')<0>v&}g#YuQi%%j7/D>#$##QG#<Q:;c%######]g'OG####;$#9@*%##/,#3%$nf6"},
{402.50,376.83,5.03,0.33,"0u#6f,###vb#_P#u$>###vP#-[$%.H$##c9-1v8q-+J##J8.HZ%%jV%l####/kHok<###eg%SR,v,%'##Lm:<f+###qJ&Kg3c.-Ww*RA(Fu$7zWiG$e,$q##)vL######YH#Il%###=R&Su#ZK2####x(tb#YxM###Je(<5#TgJ######&##1I).,#N>#S>#"},
{209.50,387.57,4.51,1.95,"{w)|J3vn%_J1,-#NRPN/'yz?]l&+6Ds#%d(3.J,tR,?-%@p/v7/{G#4I(On+sZ(Cg#.cIRn/6QPK$$fv)+).R@-v#$b[MbC6O.%SB3j7/*##>~-EJ$dQP%##cUP9C8I$)N##y7*D[@)GLO>#?M+O]3ic'$##Mx-^>$ck9=?'BI$=U9YK(g98h#$D{<8H&oX<"},
{209.50,387.57,4.51,5.46,"W5#,nHNQ'~R+Vp%.`Aku$?-(Np->u$z[$p98<Z%<5#56%^aE1x1g,$W6%t&P{S-OG#uP8|1=dnQ.,#tu%fu$5%-+l#fn.H$),oQOQ&7@)Hw(Md);%'/pQsc&[mQTp7y[+qc#4@(_tJic'hu%2i4pJ1.l#-w*Pl$-W1L'7@?%cH'>L5:Q&l?=El#Wo3qQ#RqQ"},
{108.31,403.70,4.74,3.03,".gG(c$)#61,#{i4N@-######Ud#'@+######J,#p,%W>$###B]XiR#UaX-c#YeOMy4+H%W,$EP?G@,G>#sG#X$)D#$-Z$k>$%[T+)$r[X}&$;XGhP#}#%JW4Sq<;5#N>#sB,Y[*^,$n>%[G#t%0'I#`Q(<a+)e-######+A%s,%oP$95#^l$X#$zP$bG$###"},
{225.50,418.82,3.92,2.41,"'-%*$'9XY;Q&a?G$?%z|EiY#I;=<5#H,$o-*.,#:,#w(<v5&NI+9&$OSYo-$HSY&v#vRYn##I+G#_8OG#/##N##:2={R+###>u$Pj+6dUfu#1aC@n'U&4{6%fv(PF@}H&)w+h$#etH-['D>#6l$n0*o,%Hf0~:;O>#'l#ap*56'I>#7M/`_9[>#yl&+lDeY#"},
{409.76,425.74,4.58,4.77,".,#8##f$S3$'G5#B5#e,JZ[+Fc#Dd*R>#4+?%##|$+###/oQ+u#,##*&S?5#|K8(##oq9LRDXZ'###)##2)S######,##|&SM5$###q'S%##w%S###K<='[$7tED>#'##VA,###>#$0c#W6)######;b;{6+%a=.,#LP7~08KO7p&56##Gf3/##s_<H>#q82"},
{409.76,425.74,4.58,5.94,"NA1WI&L@--###eY+##W>$9$#uh@######&R#.,#1##g$+o5%Ge/IQ%3Q#YT2&fY4l$###`R$eaCCe/###a,#h>$in)g%0tY#3Z%###@'#x-UB#?|e2>$#ksE:L'%eY:##206Bo/HAS(c$NT(6#$###=%#*sD+##pb#k$#)eYD,#2-(+$#GeYj3G{,''##DpG"},
{449.74,428.88,4.05,0.08,"^%./##3/0###2}f&##*c$%##IRS###c5%###.,####=o2###xw07/)w6+0##Vzf`G#W>$O$#1zf###a5%F##95####Y97###Eo3D7*###7](2{f;5####F7#){f###4l$-##PG####F96###vc(######Ro+K|f######|>#hh`###{k####E,$###X~0###"},
{201.39,437.35,4.54,1.46,"bn']?(Cw):n)SV*Rm)pb#gQ$Lu#~z(Ud*A5#<@+KU.,m%lv%LR'~V)C;=|k#[TN2R'VH(3##jw-29$]PNA,#X?IDd&($(5T#rv*Jz'lu&Lh-(QNa>$2v&Z('Bv(Z,$8UNR?&=i+pR-[&*<-(E>#8,#A[$-UN$Z$)l#W>#)i*###=I%E1(-v'{G$p?&&C'Uo3"},
{209.38,437.26,4.56,1.17,"DI(-H%h~*mf2op//c$95#)}<o>$SH&,I(b9/eP#wl%qZ@E,$Dn)ui.+U5tb#zKNY$'bG$E##^I*Ot5eu&.##NH&##B&.(,%%%'4}(3=5#'i0aGN2$$n~1yg&l%.;z%vNF)6#Ja@N$':#$M;*7#$1,#mP#'LNI#%7,#Uq5?(/sm'6H$RKN0u#BH<Jn-'?$=Z$"},
{209.38,437.26,4.56,4.95,"@C()3D3,#EZ#^e$lK/2g7&##'Q%ty#$mNtu$^jH~Z$T#%@1&OK*-n-pT(;h<59L(6&i]1E#$`v*;,#H&H~.,OmNMQ'gc#79,^Z'Q[#lnN<S-5lN-##v[+bx$kI.3##sL'u1;jc%L%&Z'-}^;{5%@'/=1:'@&;(3CS0R5$ho$D$&v.'t%0hu$M$*2C#4C9vn#"},
{148.21,439.32,4.49,4.51,"9p#@;@pP#OG#SC'v$,######|P#rw(eY#######837ZP#aP#J8*16'`z-]u%2:Ro>$/Q%gG#~[*ED'Hy7N[%F,$'20L#%Z8FBv)I$#KPIi6%'6R2##5m(S/$:=IX##%;;M;(gWD2/-3S+?+2@l$g,#?7,hM1AH'$##6c#DV*jGJ$Z$&m$pI(<]J~n/{I'FS-"},
{229.83,450.64,4.69,3.79,">L'TL:NH&Gm)3|5E(:U5$C[%HH%}TQOU8Y?&Yd+bn';FFW[%08&fu$?.G?(<AAG+8-Uw+_d*XH(JV7nw%wQQW>$h$&qu&CwFg#&6')vQQGQ%&RQI#$0/0lo',h;Q$*m#$kM=%##%7%BOFWJ0eY#6&+c$)?*;D/2a/-vc'B{-).(LRQY#$'/0F-#OOFZT49Q&"},
{264.63,451.64,4.40,1.78,";c%*##PG#jn&0R*qR(###P-#=-&^M=###*Q$6c#_^9Hu#4[)LZ%SH$:Q&$R$`v)#+0yS3?##[7SN)8#l#vB*-~-J$(Li8)}BJd#i~,i3AD>#Xo1N]&n6S>,#J;S-z9=?'*6#JC*:7D?C:qY#:(4ST'v6SJ,#He,@Q%oAHB$(jg$cq=tB*5B5xe#aQF;Z%J8/"},
{216.88,453.88,4.84,1.94,"pm*u>%ZD-r~1Sv#vx4`~%)QNN?&>A,y,&hG=:.+f@-vQ&}B,Vx3,H$#w'q~,sQ('^%/+B,C9>QNPd&&m&dz0Jv)KQ%QSNZ//(T&~_=TR,6##qw,oR%lQN%##8UNTsE;@*GQ#x7(W8L>QN1Q$u{/r%.h[,'##lC4$l#+s3lH(oZ#<06m^(VrB0-$KC7bl%9;9"},
{243.98,483.23,4.92,2.42,"3,#;q-T7)(v')y#T,L0%)OG#Y{,R..ZP#:5#Iu$B,#?c%WG#o?(ZQ%B|2[16*y.i[*>2S`>$)2SC?'Qm*J,#u[+'S,P#%###3p*5/Ss..fP#^w/ui4).S#7$s-S'~+ld+vJ&2d&}mKkH'M5$Ec#'YFV6%~p6kR,y94wl'w)0618#Q$oZ$a;2=#$yZ&ub87@+"},
{243.98,483.23,4.92,5.48,"jk<ql'Q5#q@+O@&]i5^f3v5%z#'J|0%~,V&0d?%x)@)Z#XNBZl%PG#8I'2SQ_m*''&7@T,A-&@T}R$m7/8=6K.-3,##p*QATd5%$##L~.,I)(.,`,#oCTkl&7DTe>$#B.*w(r`5OL48w*pH'dG$cP#)H%`G#eY#/,#i3.2n-C8,pb#oo#k}G9n(p-*0,#BL,"},
{184.03,489.47,4.80,2.44,"zv#HW>MJ,D>#CJD_[,qb#I>#vS/`P#hY#T#%=5#T5$]P#UG#<-'bl$4VQ)Z$,TQIc%-e,=?#{{>AtD.,#&##R$#:*D######IH''k9y}H%I$NQQO$''w*o:){$*P7F@:-_I,*1$nQQ4I'###t,&+.*~H'+}:)+?~#%8$'.q,f$'_#%m59>f0~d(`c&iHC2,#"},
{98.65,627.13,4.71,6.14,"17,G##3Uc'##2Uc(##=m)6##94G6#$.,#5##.Q%8,#:#$(##pm,###eUc1##xTc###xv+v##2lNgY####E##]c&[5#0Z%&##}v,$##6Uc(##aUc###G7->##XYG-u####&##]l%`#$Au$###Om*r##uTc###.Uc:,#)R*6##CW@6#$###/##}5%:H$,c$###"},
{427.98,632.45,4.85,0.97,"g5$hh,S$V###brXWS/L5$^##jnNtc(###`,#Y%#-z:.,####v?*77#vmX/##OnX9,#uG%v&#q6USG#S5$$~#(##{>$mu&.,#xm+^B#b)C+##2rX5,#eY#R##:8P$Z$eY#3##&##zP$nl&.,#o#']##yv)1p4#oX######mU+1RNX,%.,#M?#7,#YI+|k####"},
{201.18,639.99,3.95,5.17,"3W$ymX######%f<s/5######qM9+u####M,#dG$:#$H##mQ'O|8X-*###$##}Xc{k####E##fUc;#$###@##8c$%v&###(##Q`B###K#$-6#7Uc######Z%#/9Z|k####>##D6%/u####%##[_>###r6(w#$1Xc######y##IF;eY####*##D@&95#######"},
{214.24,648.55,4.50,5.07,"t$)lu&qu$M>#(q]pb####B$#eo]PG####N$#uu&k#&###.##ZI-###/W8HZ#Ao]######.(#Pp]95####b$#Tm&gY####'##x[-,##{mRbl$`r]######8?#fr[######/##XA*.,#######{5&)##3PE6.+|o]###$##F](QBXX>$###g,#ke&mP$######"},
{102.41,660.29,4.54,6.08,"g,&t##tJ_$##}J_2##M/3d##BIT.,####C##I6'xY#.u####Mm*I&#tJ_%##0K_l5#+%-%$#:$PlZ(###M##/[&a@*,u####9~/E$#tJ_%##-L_)Z$$m(4##{(8j_?###$##Q,#Y}B######'~._6#tJ_(##2K_nY#<c%^,#[7-wR,,Z$wb####H&*z?)OG#"},
{473.60,702.44,4.56,5.98,"R5#1g7>5#;5#/q417,###0##`Ta######;&#8}J9##s$,<$#9u#SH'lZ%6c$c.VcG$6,#s5$sTa7Z$###]&#7;?DE.ll'f##v>#E06jY#:5#%7I5[*$##[G#=[a@5E###T##Td(s(J###&##)-&Fl%###-##rI&uG%###N>#8t0h1<###0##0@(Z'2###'##"},
{485.27,706.19,4.18,5.78,"kd+C,$9u#Z5$idW######7'#YdW|Y#eY#j(#Yl&?R%`$+F##A019Q&0,#}Y#YfW]Z'###['#HwUx'O.,#_%#95#5pM{k#$##|J+A,$###@5#l2S3_<###5##.'*DjW###&#####({5######iu%######S>#^(7j>%###*$#<~,>v(###C#####G>#######"},
{305.67,746.09,4.16,4.12,"2H&^P#3,##:3c5$SH&qY#oA1}J05l$###Fe)yd*SG#dP#w$)CR(H9/N5$mG$<%&[SAWH(>#$?3U2B3.,#a-&Rf2/,#5l#~J-I#$4U*V/30,#bx2lh'8h<.w#LTX=-'F>#]C#(M;B#<h5%^5#%c#*6#8i;y?*s84)6#Uw+w#8A`B'7&Qn-Qh&u>%oBHeSX2,#"},
{181.64,754.74,4.22,5.11,"BZ8%H%ow'NH'';M###$##B5#X9.~#&###.##Xc#WI,G,$3,#-].0$$+?ABm)+4[],$`Z'Q>#aYL######=6#Au$eP#c6&R?%*x09:$?cP%?#J/[+H#H?(30#_HTD>####%.#rP$Z,%l>$lu#l-+qm#r7/&N't[Y&##ZP#s_#?~/wP$###D$#T#%7Z%###'##"},
{150.52,755.04,4.66,4.58,"me)/w,[6(3l#vR#%@Nul&cG$at3{$TLl$6m&F$Eg,&###0Q#9&/n7'[)?86%HPFcw*vI*72:pr?ZH'9v$O)Tb$T######EU%xQ(([$?&T{l%J%T?u#bw,=%'u$T######q0/AU3D>####~Z#%d(@1)Ex0g5#7$Ty,$PQ&~T#1)?######@[#7e'Y>$###4,#"},
{157.84,765.73,4.68,4.56,"K~SbP#Mc$Bd%-J*Fn*;/0M#$6@&Z3>B7)r7.AP<):89Q$D12J*EvP#:/,^X05T0q%'.~S]H&/~SBQ%{Z&[{9{T7###S5#RKI%M9SR'1%,>[$iw/KS%{~SXc$&[SR>#I-'t@$gEC######?Q#6[(7`5.u#J#$a83AN-.$(:'%aZSsY#ZP#(_#:d)###:5#/d$"},
{148.29,793.93,4.76,4.62,"~##]oRrb#gY#H2(WFJ###bP#,b2)]3=5#Jl$DS$;f3_P#+H%U>#f]0$H$Hc@lG>pm,4/#J5J1pRh>$}Q'6w&=@'US,%f0G#$`07PG#Rc#T57D/1=#$c)*X-E)mQ:c#pf/,b2VJ.cm%}ZMWc%lx/oP$oG#f6&Up21-'Id'zu$Z_:@I'@@+uH$ER+UR$}W?7l#"},
{125.96,855.59,4.35,5.37,"0,#[D&8Q&###LR+DN)'$(A-#VZ'b##=l$nz'###M##=^69@(96'7r&YBV1u#]@VY~(mc%,A$KT4###=Z$ro%###(##flJ*l#xY$'##Ob2.@+M7OD>#o($|g9vW/-d)2Z#d>$a##5Z%*{7#########8$#0Z%/c#95#~(#>075A#-%-/{%ac']&#,n--d%###"},
{118.23,60.15,5.09,6.16,"B}2Yl&###:,#rK&7*E######zA$BtK######Fx#'PL######[k5<Q&m,#,e(mW-|eQ4v&Sc%LhQd_>cP#B?$2A.K#%hP#I5#Jv(H>#m^+6L3G:9cH&U~(TgP6dQS5$%Z#xv?#w,6,#fQ'NZ#47,#v#nT5su#AcQ'##x#';2+wu'###1,#Vy+D>#G,#0~-U>#"},
{73.14,84.80,5.59,1.11,"8T-*;:###1c#c3?,v'%##jA*nR'lu&7R$m]3;<-wG%:-%W>$sA,}f0mJ0f?'opR%.+dP#,.&C]3On(IA+n@,C16|d$;&0N>#@]1(Q$nn*tG?hmRVv)9c#5SBG.+iu=PT0Qv)?,I^g,)m(N##6?'###s'#MpR~,%i#&jR$WnQ$##`6%M-G_M@><@}?'>8-*$B"},
{410.50,95.86,5.23,6.15,"8F2s:>###7##nU$&aF$##>5#>R#l2C###oG$Qv##.,'##9l$R`1KrB###(##9gPgB8###Af(kw)og8_l$x*Ayy-C[)Ll$(m&o*FOG#&##qI$-ePCH'###-i)=e*-_6f>$M,CO-)^c$J,$,56SS1###5##p5;/cP######U<.DZ&.,####%gP######6##{*A"},
{23.19,105.64,5.24,3.01,"cC5L5$###$##c(`5d'C,$*##k(9}+.R>O%#####M$#}T8###cjD######4Q$S'`?Z$^[)0R${n0DT&[M_>Q$###+H#?[*###xq>95#`P#5u#R(`|6%0R*?,#}..OV%^~[4,####^u#jZ(###{g7]#&.,####-'`V5#x%1(##D6(A6#[,Q######-##T,%###"},
{62.01,116.22,5.52,1.37,"DI*5-%GZ%~c&q-)^v#U7-ZP#p,KM>#D>#$##gI-######8,#LI*:w(wc&#%*|-(UA(]94DH&:qZy,&iY#`?$]o3######A##gY#e>$i?';h:sP$bG#>^0'ARjB7<5#XG#(I9@d*######mv%<H$d5%K6#-an3u#tG$Rn#Gan8#$###{##{an.,####Y$#Y`n"},
{62.01,116.22,5.52,2.37,"=$#FD?I>#.,#^%#px5######4@$Bv)######T6&D>####-##_###%,zn*+Q%Fo%A,Ltv)(l#~AE#B4SG#_>#yU:######m&#[,$T?'FC6jl%C-(Kl$(g,);8G]].u#pG#N`*G'b>u$###v)#py1Un.l>%oG#*S.M>#r,%)A(H,P=5#cP#(a%qkK#_3WK5iX)"},
{321.65,128.26,5.48,0.47,"9l$cG#6l$+##H07######N.#)V=######6m7pb####;##OM[bG$3,#&Z#[?%<J[######OR#(K[95####VB%Do1D>#(##I3;{k####?7#~J+i`D######vU)WO[[?)###mc#nI=7m)%##'H$I#%###Z(#2&0@U8.,#j##gT)B<<xY$'##_&)(xQW,%:5#dI#"},
{282.29,160.12,5.42,2.79,";5#^+>1OA###(K^i-%.C6v6#_K^@B4###5$#_,#)&/.,####e-([><XA2%##TK^.8'n?*s##7mLv4C;/1lP#1s.+L80u#iP#,U4J>B_5%$##@L^_c$@u$P##dW@<l#gT0k8202;E,$Iu#511>N>Ex2###$##(N^OG#######{+<.,#&##Ol%I.(X>$$##uu%"},
{296.28,165.19,5.42,2.81,"6##%F9(=ED>#du&b?%Z_Y[>#y[Y1,#y?(&T#K;@{,%95#z##A5#CzSOv)###wn0C6<w2@*##|~Yod(Z6)l##j02D|B+u#+##HH#_+=fJ2###oNCL32WS1/##b~Y3-%Y6)c,#pw*nA/,o/S6(Y##`F?>u$###M,KlB3###$##/^Y/,#wb#TG#]@*.,#=Q$-S-"},
{285.84,178.35,5.32,2.83,"$l#{u#5(]YG#x&]{Y$Cv&Gx#B&]Q5#ZP#o$#;5#JQ$6#$$##0H&Wp._KS7,#O&]r5$,x.s@#'wT9'495#p##<##l7,95####K%-D+6+95%##G'];w'F-)?##A`8u?Mem+1,#,V**'6/,#'##yT5`(2jZ(%##9'](Z#x5&:##v~/cu$Zo/0J._z8dG$D5#bw)"},
{199.49,185.08,5.39,2.28,"]>$bp($y5%##Aq&]gRVS0###0?9P6)###%##Tv'xY$###)##/,#yZ$K]Prb#sV2'.(GVQS$(bMSOG#'Q#>~&318,u####0##Cc$###U}7}k#Km)###3F/#[K;}H2?&u[(H/E)0+E[*s,&D,#ud&QG#(/*###6J)###,C-,v&w&)>#$5K.Iv'-S$(H%Ay-###"},
{199.49,185.08,5.39,4.66,"###g,#{Q).,#%v';6&&);wZ%co3f6(G96GR$PA2O,#pf6^?$###(S#wm,.,#Zd&ZCNcJ1}b#jv=>~S[-*a>#yDA&-%a$+D-####(A#L-T###yB7h2){-T(I%}.Tt?(>c%&p#V=JT>#CZ&p$####:##K+>~#&E-)0##`??s/T-;?$##gG#8#65T3#[%xG%LQ#"},
{434.19,283.11,5.32,5.32,"eY####%##4u#>u$###%##eR%KvS######O1#NxS######D(#Zl&###$##w5%s(>######Up#KwS######L1#C{S######2##a5%###A-#dr?2mG###~H#tZ&<eI######G##Bv;95####$########Vz&o{?Aw)###y2*cG$b5:###[5#E>#lo'Qu%######"},
{145.39,314.11,5.12,5.21,"x>#4yQQG#;5#a1%ouQ###=5#y9'wd+)?&[P#%l#3Z#SZ&f-)8-'e13JZ$^F?w=4}d-B8'q2>dzQ+u#lP#Xl#I@,'##mG$uZ'2sDz5$`$)sX2Q&3pH#FfQ[^1NvQW,%56&&g#,7+3?%NQ'U,#Eq8]%(8w+$o-q@+>2/^08X5#2u7fGOgY##6#J0&dx3cG$&##"},
{490.98,316.50,5.26,1.87,"############95##########o#'######-##F&3######?##D[*###Xm(0,#5JW###0,#;##`HT######)&#9IW######4(#3n-###sa4Dl$aKW###=c#h5$bJW######{%#XIW######z'#.u#%###u:7#$nSLC5#nS/QG#{NW######*##2MW######2##"},
{109.98,323.93,5.39,3.37,"qj,#i?zl#Q?()v6)A/###0,#Px)SC:###1,#6I',H%9#$###Y^16H$TRB#~-X}Z26&q>%bl$w=F%m'>5#7]*s@,vZ(###G,#3[)?U-huNw9/3xZN>#(-'xr'EaF0v'OG#cm$,Z#{7/95#+##gY#Ml$bl#d{ZSw.C,$###J+.&?&BQ%^#&[##b5#m5%6Z%###"},
{226.35,324.91,5.77,3.98,"E?$bSQsH)sb#yI%;RQnJ/>I)$7%ycO)%&<19/H&U&+t:1z8.VQ$RN=-x*gQQGh-_rBtw%OA0~SQ'B3~>#S*>z,'7$&2l#JTQcu#7/0Y9(i_>w$+x.*/j0jh6fQQ}k#Y#$8{(6J/<Z%1u#vd)&##Ic%7g$F2AOG#q#%O?$dSQo$+@Z%xZ(T2/e6%zd,}G%dP#"},
{252.15,402.04,5.04,0.23,"47,{Q$S?'*o)?R+cP#Dc#poJ}k#<5#vf%/;>######4J(R5$|ZQs-%Qc&lI#*[QU@&?;=AN2f95V?&@>6:8/.,#%##[$CiY#w^Q9u#AR'N#%G~C|/4^^6XH&';:Si6l;>W%'D>#yI$)[Q<5#b['0u#>.$DA.wU$Oh>c$#;S/]8'9~QeG$.?%8H%y@)IQ&|b#"},
{252.15,402.04,5.04,3.34,"nl&|b#n,%|R(pP$~u$[8'KeQn-#4~.3h$vL=0n#$&.B@'tb#1dQ1,#95#`7$Qr>Uw&(M99`5:C6gQ&iwC5]2DR'Z,%|gQ,l#u?CiY#.,#%##mb6`e.8p4dH&U_<z)22dQ]@&v5&c@#(dQ-@%,8)E,$###$##sf%RV>,u#<5#Cc#5'K?R+pY#T?'ww(*.,Dd$"},
{209.88,423.00,4.99,0.84,"p5%Ne.},%W<:N@$NH':y*M^8&K%>>NgG#*|AQl$H?NlY#3<>9U1PbEz5%U0/&w)fu<rL6d$*%@Nx#&lG$~17tc'S6(g$*u{?[H$$@Nq.%nx4,z7_bE?5#JW2)?NTx,,u#cn#eZ&77MU,%/c#/,#X>#w{-T>N`5%3,#Cd(oAN(-'6]%s;AWI'@05SS*3l$59)"},
{336.44,439.58,5.04,1.53,"h>$sQ(9?#`sF-B~{k#>##Q03xA~$##I?$T^3_Q($##xu#~yV######~e#;A~)FF###-t0&sBBC~###S#$-?$az<######mc%$##4,#*/)381V7.%##2BE}d)YA~###C-%,f&eND######F##$##2,#B-&Cm)T,%$##&U,Z~-XuQ###b>#D0([97OG####a##"},
{275.98,445.29,5.01,5.73,"U9)osE8&+g?(MZ:X:<E,#DJ)d$%~0Qkv*c#%bG$l@$8.-4##y^*Y~0X/#_.Q_1Q6#$zQ#2.(Wx2H##]rAvG$Q6(1?$_J2.##|m)<d&.b3O.Qv,Q6##FW-4X3t~2*##K3<kZ$vY#Mv&6D4ZP#1q,.<;dR+/[$k>%e#$Om(C1,eY####p$$Rx0###,##yh(5w-"},
{242.47,453.77,5.69,0.74,"}[$7$P(Z#d'7>7)`xLkZ(&7(0[%>#H~P#z#%GS'On/)##hm(S@,}Q)lQ'cPKqQ(}f(n#P'x.+$P90.*m(EV**]/3['ym([4;fm)%mL<c%mu$f$)}K4N%PKR)s&PN^9Ff*}c&<T,=&Pzu&6?%jU9f%+xG%}g._Z&C7(l&+N}BvZ&L.,~g(Mq<=Q%5'6z5#|$P"},
{242.47,453.77,5.69,3.88,"z5#|$P=Q%@06~g(Mq<vZ&L.,m&+N}B_Z&C7(xG%}g.jU9f%+zu&6?%<T,=&PFf*}c&s&PN^9N%PKR)f$)}K4<c%mu$fm)%mLym([4;*]/3['*m(EV*+$P90.n#P'x.qQ(}f(lQ'cPKS@,}Q))##hm(GS'On/~P#z#%0[%>#HkZ(&7(>7)`xL(Z#d'7}[$7$P"},
{193.64,457.39,5.73,0.62,"T,$3n-%]%x,Nyn0Oe*w,&Bt6,S,u#KeY#U6%5?$Py4J(:+v&g,$+J.)S*.q;bZ'.g'S,Ns@-i,N^K.5d)i*58m'nR)sS*u,I4@%)-NT,%KR&rv(G92C.NpI+U/N+3BE]+R6&U.)R/N1$'.v&ZB4qm*1,#IW-Ll$o&-n^.@q5~v(uR,op'1EB4Q%jw/(6#n-N"},
{193.64,457.39,5.73,3.82,"$6#@RN$H%1w,5C(v2B.@)-w*e0-{V7d#%_/,:5#C2*G'5G.+PQ&<$'1n(?TN=B+?$'mnKjrArRNQ.+@R(S'41Z%6v$[-%hQN7o*mPFW-(q.*6d)Ks3QQN{/-+QN^%-*v'%''Ln*@EC;l#SR+g(:u?(v#%tx26#$xZ%N7+kHLr,&)|0`f4tw*dI$oQNSu$YJ1"},
{232.68,489.78,5.49,2.32,"Y-#Kr:HC2m>%vX2d&4Ml%###hI*:l#vG%@5#6,#pI':@+$##@H%]J*[9RVZ&K9RS?'[h=K?#4j@3(9Uu%(##],#ZR+Nx,X>$TQ&BJK%;>Qw&gPMgm(+&/Y;*IJ*U8RSK0Fc%fj01M=-Z$:#$c6(pH(v[+6F?o=DEZ%G7'lp---&L['-;RGv(O8R%6&lA1c5#"},
{184.13,490.37,4.97,2.38,"^m#F<?;T/OG#=AD}d-}k#jY#='3F>#<5#G-(0,#:#$<5#7c$3$'kQ%_(R|P$b&RDc%#~,PH#O4C-4F###&##F$#;M>######nu&A=9[NBF@%luPX6'bd)9D)B%*v&Rqp/RH'b1%1$R=?&###jl&-m'<m'qs=Hs<yP$F$&1z,96%g#%ySF1.+CU1sQ(Za?/##"},
{96.33,596.87,5.78,6.18,"T7.###Szc*##vyc###l%/1##oL;Du#pb#O>#6#$^G#&l#U>#[[,)##%zc,##byc$##S]4~,#)sD:5#1,#RH%PG#L>#ju%.l#p-+0##jyc%##=zc(##~]5,##m}I:5#.,#C5#P5$A,#qc'$##VZ'###>zc'##byc###)T22##/jC2,#.,#&##2u#m>#h,&###"},
{411.46,603.56,5.67,4.42,"T,$m>9:7,$##TJUjA0.,#r##}HUJ>####~$#{#'md)D>#3##$I*l,$p@FN5#CIU~P#~5$'0#_bMS>#HZ%ve$%l#{d'O~,.m'D6(###0MU>5#+LU###pm*d5#TJU###vP#0H$pb####;d$w@.ZP####,qPkR->@U###_I'1k8x5P######C%$pb####%##PZ%"},
{430.08,688.96,5.39,4.22,"d]Z######W-#,^Z=5#.,#3$#eR+0n&-6'Q,#up2KH&OG#9,#6&Q###A?#|uB4bZ1,#:c#'R'dB6CZ#ko,3A.1_;eP#S#$7x+gD?m7+$n&iC3b~Zj>$SZ$j9'piD.,#j##PL4h5%.,#C,#^(7PI,{l&b%$&^Z/~Ztb#@##V$;}S2+u#%##Dd&Fv'pb####(Q$"},
{446.04,710.18,5.53,4.58,"6o.+u#%##2H$J'O95####[5#ZZ4M5$.,#?5#C7(D5#1-&6H&I..###$##q@'kYJ:5#(##=D-qG<8u#m,$%15;vI###*##n#&Cv)###$##Ho+6K[M,$=5#K@&3h89I&hl%1J+JwT<5####-d$du&###'##2&*#_<###%##j^HYZ'3,#VG#:u]J[+######0r]"},
{136.30,750.57,5.36,3.46,"=u#WC2v6&BD:*D9(~,UG#dR(U}GD>#(x'gm&sv*###QZ9IZ%:Q%;^.t#%*.CqC2d@,rM2Rz4}LSMZ&*[%&8'aw/###~c9gH&.,####)?#4MSo,&<5#c`+nJSEbK:5#$v#{qRG[+###D|3E@&<$).,#(##y.(cc'###9##=r8.,#$##f5$(V6###<5#*e)+Z$"},
{136.30,750.57,5.36,4.83,"3?%16'###@oK&d%T%/###*(NZd%[?)###>AJO,$###/,#b`=Uw%P=Htb#{l%H&'J]RVl%Gv&]_R4aF0,#tQ$B:7###%##,6%9&.$].KT3xY#q`=@B/:[(z4<,~R*Q%2,#ra1Z2={Y$###J##Tn-T#$D4?a?%M]R+l#3e*_Z$DRM.,####zl#Ge&I?(###%##"},
{490.86,753.57,4.75,6.11,"CB*3l$###'##fdJg,&###'##A;cZP####.##_e0######$##aK%[?)###'##&NYk-+###H##g;c######3##F7-#########N94{k####9-#RpbeY####K(#s9c######j%#2Z%#########Zf5######7[#O9c######C(#hGP######i&#############"},
{418.35,760.21,5.70,3.13,"*Z#{n-6#$4,#P%&)w,_P#>-'s2>Fl%$##|J)?*F######V]#=-(ou$hc&^l$0T-0f*C%.+c#pz~8?&ZP#s##'x~######('#,u#:##Xf-P?(sS2oQ#Z1;k$%Gx~@5#-u#QB#xw~######O&####$##;:+~_?VZ'$##4w&mUSN)B######5/?nU<######nt/"},
{185.71,770.37,5.64,5.44,"|%)o6)Tg*@%-nO/e$+;d$j?)F7=0Z%.,#*##S?(qG$HH${#&Z*<<[)TB*r7-/X89m%W9VqH)%=XfP#xu'G##n]5sP$7,#[>#Q~%I^6Ac%V5$<q:WT*_PNo%$X7X3,#DH'5D#gw/w,&OG#@##Ql#H80m6(Jl%t7Xqv)Ac%3T#wvVrP$.,#_0#hG$H6&nP$-##"},
{338.33,806.17,5.00,4.40,"^5$g_$y$X###1&X8.&@H''##lQS######,##F@,######+##O[*<V${$X9,#$&XDQ$j>%Y##y%X######+##MJ-#########w6+|##Z)Xf/1A%X&##=c#f/+.?I######)##8K/#########jP#<a?fa.up:|g7=8/].&+167D1OG####(##Bx*#########"},
{223.71,833.39,5.22,4.46,"E>#0$#NJ]###<&]j,#AtKP5#BU9=5#)y-Sl$0,#.-%:_8/,#tP$U/#OJ]###>K])I#C<E2##I<DYG#T~-K#$/,#[G#Ws:g82Ru%x##jJ]5,#VJ]G,#@WB16#kNF###O5#A?%_P#95##6#S&2y5&s$#QJ]$##yJ]],#j&52##i_?H>####I5#>#$^P####gP#"},
{173.04,844.93,5.26,4.48,"$##A$#|6W###%4E#[#-tJ###>~.T8#=2@######sy%X07eH'qb#9$#)7W###?7W.6#2[U<##4&2[##t~Vzb#4l$5##U(5PU2+u#2##J7W$##/7W2##IIVU5#DS0>##4<;V5$p>%$##iw'4I)|k#N##}6W###57W=##qXJ'Z#</2%##Jp,%6%0,#KZ$^8.95#"},
{66.42,70.32,5.94,1.11,"K>#P[;3e.iY#Ye)m3A###M>#.16)$'0,#5.&JI&mP$:Z#e$)Ul#VS(v)C~c&WT-in,C@+cZ&0_[Wc&?5#xc$3@+J-%D@)cl%q6'vH'#G82:8+~,tG$Id'DH?:>K=H&Cu#Zj4V?'&:-Fe-%Z$O5#I#%bM#V~[ZP####5)#;_[.,#$##kS%W&XH>#uP#D;5|p:"},
{472.00,146.62,5.84,1.68,"{v,######F##E/^######H##t^g######vQ$-H&######P&0u]7######x##+_g######d##u_g######6Q#yn1###$##Me,^g:######S##x^g######k##D_g######s##C^9######a.)005######7##J^a######+##q`g######,##-2>.,####9l#"},
{266.41,166.37,6.74,2.53,"W82'H@hz=D?#KJY{c'^6)WV#uo2)//###o$#AZ#B-(######:JY?`8D>#9&#ALY|(;U6(1o#~Q7ay8_P#tP#i$&$.)L5$###?JYs#'/,#1/#r6S$H$}[(jg2$QIF,$}Y#wp,UR*Im$If2ZG#nJY######1$#J%M###$##$Q$/L1L5$###j,$FJ-?5#5%'~l%"},
{135.57,172.38,6.47,1.64,"K%-yU-#.,5##>RT$&)TA2D##%s@.d${lS0##I#JSG#eY#'##wn0+8&QlQ5##=@V`Q%vK8v##@F?ve&x?V0##fDV2Z$%Q%&##pn0HQ#2-R6$&O@VHc#hf52[$k'8#[#pAV6y0^@V'##Kc%0g(LH'iP#x$*~B0t%04,#[#%j'.Y,%4,#E7(F:6xn1###'##6f("},
{135.57,172.38,6.47,4.87,"'##je(781####7'X93}G%3,#Cl$&(-Bo22,#Oe-fA-^c&H#$^l%g'*lIV%##CKVA01g96=Z#to5,-$5JV)Z#IIV2Q$z%0zP#vG%'##><V{P$9IV'##n*<B.'L1<0##qIV8?%Z#Q+##1f1uw&###$##L+BE>#+HP'##'a>+$%U'8-##vmTNR&1e.-##q@-z'+"},
{183.87,177.61,5.96,3.42,")$#SdT######*.#~_<PG####'##>+>@#$5l$###j*5*H%97,&z)IfT6v&dy4HiTuD@S,$-f)}$+u020d%_FH&##EC7D#$[?Onl&Pc%'~%tfTL6OL-)$##H^K(@);3?%Q#FM=$##?}C%A'*.,[#$bm+8##}sF:5#e/25,#HkH&##p]4#Q#'h:&##ad*K&'B7-"},
{125.94,187.74,6.59,1.74,"uIV.n'7S/)##AF@pZ%O7W(##2:NzP$uG%###W@)pb####5,#v7WEu#}/5/##eXB,R%g7W(Q$c<W7Z$x%0rQ$iU3A,$###>##t7W`G#'~-HA(G]2=c#8yUBD957W###%I'W#64'6######u6#l?)/,#XG#l]+Au$$##fG#(/-J82###$##J.']81######A##"},
{125.94,187.74,6.59,4.81,"###K##A&1###$##;n&0o2###[>#he-cG$$##6,#7]+-v'/,####PR#Oo3###@d'7Q6>dU$##WoS)W;S81Qc#EI+*y*d-TeG####:##Zq5{k#&]2#R$J`WHZ$w[W1Q$b=Fb.%e/4W##l[W-Z####7,#n.,OG#iZ($##i_W6Z$^[W7##74Di6%kv+A##,~Wnn)"},
{138.51,188.66,6.10,1.38,"(D<RV50Z%D$#'{?u7$+;?37#d-TiR%5J0}##N'8.,####z##C-Tru#Cv)G0#8D<iq%8-TW-#90TE.&L@-l##s1<######d5#TuPN$$uI.zJ%1B4d1$R.T_x,3.T:Z#%v'a1'IV:######D##PQ'S>#*m&$J+>u$B,#HB-V=A/&2###5,#N*-C6(######7I%"},
{138.51,188.66,6.10,4.89,"###@-#p/4###Zc$&+5by9###&r7X94w-*eG#PZ&`J(^WBA5####j##0<>D>#:04Vo(bCZ$c#:BZaH'xL9)-$d./1##nAZS>####&##%U-+u#Vv*###cOY~Q&;AZ&##U<<=@'#A/*##pAZ4d&###)##EI*ZP#6#$%##>JP^P#|cRH>#.E<:$%`&4A5#.tFD$%"},
{332.77,193.11,7.05,6.12,"M.-$##$?I}.'_07###$##Fm94>LgY#0,#H**'m(0,#rH&j#%w(576$/8Y7c#%;YG>#D>#cV,TL9hG$iG#.[BQG#;5#DR&ed*c%.?l#91O,GHx8YcG$U5#i,9vNAE]3N#$%n*+##Md)p-(yY$cc$-&.r.(1h:)a:iZ(-##R_1Z.(Ve0###-##DZ#WR,.,####"},
{314.07,235.91,5.98,1.43,"nv#beTUR,$##;6#a;5Ce/###dl#-'06#$###nc%e5%###$##HiTddT;5#E$$;/Ie[@8o2l>#k&NLL4Il%K5#Xe-eY####)##GdT],%xY#HiTWNDXl#$4<3`2DdT$c#]H&](.dA1@u$###?##|b#eY#+##[,GfY#.,#?H##g32Z%=?$96'l~+3e(R'3{#'(c#"},
{456.87,282.58,6.29,4.79,"{k####bn&;A0Ce/###+c#mB*|p;######o42KvU######F}0_5%###9oQe#&IvU###A[)G6#evU######1I#[wUOG####E##L5$*##nwU)##eHR###R7,?##GxU######'##/yU[P#######0,#--#*cO###K(58l#gR-%##WkA#########mxUD>#######"},
{358.67,341.93,6.04,3.59,"*$%9l$dJ,p^5bV:xY$<$%m^1EM:Ke..,#~h#X6)/c$:5#y2(%l####1]'DjWBg5I>#5UIQN@nfW;5#bZ%KM'PB5u,%[P#E.$###*##w%)Dh2l,&$##*67U?Pj=J$##&@$SiSV[+Kl$X>$2Q#2,#bG#c#&E>#&##Lu#;Q#Wv*###-Q#~l%;w,###E,##?&<5#"},
{358.67,341.93,6.04,5.99,"###'##8l#FZ&bG#,.*lY#+Q$I?$A6(%##eZ&<5####)##,v'6,#JX>B?'Ac%W.'BzLp84UG#WMPUEApb#nG#pR+$##`P#<l$&##'(XZ$'uG%6h9kv;c98xQ#l(XT7,(c${8#691$x1###Jl####/A'd>Ho#'1^.u(6em+B##$[9pmWP5$V,#Mg'z$X###}u$"},
{234.09,353.55,6.87,2.27,".,#(S)gG#yG%,e(taE%##5Q$MR?-mQM>#bZ&h?%XqQp~0).*'Z$e|/7[*'##zA*XnJyg;9,#EqQRr@y5%T.+|v*cU3(j<[81hZ%r]*RD<D,$ZK5)g({^;He&{mQ;x0Vc&Zw#;w(J;<b*ChG$RG#p>$i(,_p9S..s5$:&+74A,z9v.*%L7:%%S)8>L4T]4b5$"},
{375.18,354.38,6.67,4.82,"]u#<wM=6(###{{.@:8gu&Z,%~F1lcSG>#aH&4p0Oc&ZG#Mu#kR*EK2`G#q^8q5LV$*C5#tQ@iL7YI-6,#+iS:0595#O,$]x&3v'&~(2?$4dSEfSd>$a5#Vz4TtG###D5#9C3.,#K>#?u#g#%95#%##/J%IdSY#K###Fs0&M;NeS###pP#;-&eY####n,#fm+"},
{375.18,354.38,6.67,6.05,"*##oU-2R*dG$@N0=O@j>%WG#|w+OG#=5#q>%######K>#9Z$/c$oxGR82~G#IpSay.$.,O7#U,HX83###JZ#cl$z70###5?$;H%Lz0EM<T5$'n?P=HO#%x5#<X/zlSrb#xQ$JJ&rlS%##SB2?m)`?$}lSe5#|T6|l'X%-s$$TmSN[+wc'nh&b@-;c%;c#YJK"},
{182.28,371.23,6.20,2.20,"yM?Km)`H(il$Td'w?'x}IP5$H1/u_:z%/2Q$X6)1~(Pq<ZM+3J+;i<$A/tG#Th2j]4)y.GI+.w+3w,2w$ULTwK8xQ&#l#}u8Qv(:9(~;AvG#mJVfZ&PQ'T@%<mM`$+,##S:2KMVin-+H%JK']v*2e(46&5%)&kFKe.$H%[d#W')D>MBu$&##|{)3uN7,#[,%"},
{180.34,396.92,6.51,3.04,"Z,$CB,}R(@d)Cj2=K3I7%y,&>6HQG#hc${_4}b#xP$Rz&GFFL'0*~,cn'.A-mGMJw'*Y7;[(S/XjG#N.-.q)W96yc$wk;cB2*v#<]1_/'=[*?q;G]+e-)%M/t3Xb7.36'om#7^-,/IqA2H?%U5#oS)W7(Qc&*y'kx2+Q$v.*>g)490D>#r5#q,#GnJ+u#0,#"},
{204.55,403.26,7.63,5.34,"iE=[F8=H&JA)F#$Oc9csEU$'bm%_i<m:9$v'+A/nP#-K.&g29f+XI@uz>eG#1}G_F22/.wg1U;;Zv([q1c04ry7C?&WT3.8.0e(JK264Ak?)[,J]v'>g/{L8ax/eo*mB5i20~@)/%,I8*1AOdg5>?'p.)Hw,K/1Q-(^L7kf39-RKQ%1f1^W,Rm)IS+Pp1>=>"},
{355.96,450.65,6.90,1.55,"$x/V{;/v%o>%bfX_A0J,$?,#[P;heXyP$P#$l9-q%/|k#1##d06Co/?Z$`x2)fXjZ&mZ%300ir=*02ed&CjXaq:S5$5Z$*y'zu'/,##f$?eX+fX###$w$?90y@X###%##NU1Zl&.,#&##G-%######ff*Y7.AdU###Ax($T)%dS.,####n##`#&Q5$.,#9,#"},
{388.52,451.63,6.77,1.44,"D5#Y-I0:9B,$F9Wqe*t':%##K5HK-).,####B5#e-*.,#SG#=&/E9W.?%+l#r8WZ94A,$@$#9APoJ3###8##ol#fc'###=5#3jD~w,J?&<Q$@8W`^6lQ&`7(+e>{7W^Z%10.dh3K.,|Y$o5#22?c5%N?#Y_;&7WO>#eZ$S;6?/15u#-6$DpNP[+&##7Z$.w%"},
{163.43,454.28,5.95,1.41,"mZ(Il#el#7xF},':Z#}D1|8*Ln%D[)K2.#e-7-$IC0GZ%IL61,#gS'VwLx96j#&T|/LS.d5$Q)@5K+Xn)@h3,m&2I&oB37a@=6'|D+yrC=Q#2A-BK&`-PiP#~0Pdy1}Z(kl#+.){`.#V<Z#$u5&g$&'$%5z+K#$0,#)G5dJ1rI%Te-E&%VaHqG$Ho)'?%q>D"},
{163.43,454.28,5.95,5.31,",##88M/,#D>#?~#H6OW$*ZP#i7OVH'|S0B$%h9O}x2bP#KR'{.0Ri2gu$,r7Z:.F&2no$w)At8OX>$qG#zQ$Q22da<+A0~,#=7Op>$X6&M/+mo117'[u:<f/k6O':8HR'gd%).(Lz6)OE9,#QC5Ku$)w%FM;+w(>`-AA0$6%yH&/F@k#%GL6C])eT6Gv%im+"},
{142.18,458.40,6.53,3.57,"_##O[E_v({x5rM)bjHp5#q7/|s1Fw.$##*##%7(8Z%%l#1##Dd'JJ+.&+;nQb14IJ&QvEeq=JpQO>#x5&k,#a&2#l#UG#?##,T-2*90.'*K1Ed*?h*_K3t>ATmQR%.,c$oU&VJ*wn1###'##AT2lZ(bQ#{^5Wl%]>$im#,oQ9-%N]3Mo$}mQl~#mT7?&$Hf4"},
{373.43,461.41,6.34,0.16,"E>#p##XW?6*@gu%xb#v>$)9RS.,cP#/,#y9R;x+:I+###ITO###5##TI(_9R/~)$~(d.Lc7R)9R>-'HZ%nmEhd,######(%H######>I&)y3u,&`G#G-<tFG~6Rvb#}Z$^.BAJ0###$##LH#E>####^G#2Z%###6c#o?&ov+}6(~H'=u#^.+rI,.,#/,#a>#"},
{373.43,461.41,6.34,1.64,"y%Qo^2|H*###GBSw6+######R-$oI.###$##^>#em*###0,#rAS]d*|b#1,#`8NmASD-'u#&@w@a6SkY#.d$tI*Q7*kl&]G#HAS(H%`P#?-&G@SM7)a7(t8QF=Es,&bG#_d>9I+$##<Q$-R&*bG###%R#vz;XAS###,##(n(pw/######1.')c$95####&c#"},
{373.43,461.41,6.34,2.63,"TG#2?&QG#/,#(##-x*bQ(###eu#w@,{$*[u%C['J#%'##~u$-&#2cKxY$###5R%nxH*C9fP#XfH2}Ch5%FQ$_w+.,####=u#/(#dmQZv*###`:)(.IpNB`w)xoQn7.v#&{L(39.v6+###/##y?#F]0G81YH'O_&4GLf5%+c#t{+LlQ###pQ%U(%cWE###*Z$"},
{318.92,623.53,6.18,4.71,"'I%Lf(>WB###M|Vt,&OG#,##:+C.##Z,%<,#sb#N8'8W<(Q%g-*vP#CDXF5#CAXG>#;?&t6#^@X###qY#=m#_5%'##-g*zx4?6($##rBXzb#jAX###^R+#$$RAX######:##II,###&##<?%ZP####7/G@o24kJ###m?%>b<0IQ.,####q5#}I-######1##"},
{283.07,626.64,6.03,4.30,"7,#~^&;K5###sFDuw*3l$0##CGM######^##jH)######3##y-'}41g7/###J:YIn+.,#E##P7Y######}##$[)OG####J##4'3o')b}L16#I9YEu$###^&#T7YmG#+u#B$#:c$8N7_5%5##s:>l##oMAGT#G7Y######y(#8~/U,#c%-An$#Z#T{0=?G#d("},
{205.60,684.43,6.66,1.64,",c#k{3%3CF>#</E{=H%Q%A##E=A#Z$xP$I,####[G#eQ'###jc'8p#TmV'##DnV~Z$X-*(&#tmVzG#lZ(~$#uP#47(6d)###qv*3L#*RShY#TqV~>#Z>$%H#;nLpY#K$)1##z5%sY#Ko0rb#qP$3##_k2JK5]nV###w,#EU1(oV###aP#4##o,&###X6%ku&"},
{305.20,725.57,6.31,2.06,"###OH%{7)95#&##~-''`96#$t6&,['Bo/~P#xR+D>#nG#)6&zc$qrU+h;I>#h&F6oUZ.+v>$x~/F>#f^(:945A/(##AQ$y7,u6(d30qrU`%RVpUD8/^H%d0J.p6*##f['Z&-###/##S@'$d(Q#%###*]%}}KA,$*##,$&w?FD>#s5#X7-s5$###8Q#B-(###"},
{117.88,790.59,6.89,3.74,"+-#-4EG.)B,$l;'pu&h7)zY$m*-###3##+u#?x$#########_h/e_>(I$|?*)oQ2Z%xZ#Oe+&j*###5y#4I*hE+###[5#.,#$n,$##BM.f3<~7W$l#L,#+V)wB4###mz#f7-Pv:###59#^d+V5$qd&L$',(61n(U7*,Q$N~)RY=2Z%L?#bv(L~-###+:#-~-"},
{205.73,806.42,5.99,5.66,"F8+mB7/,#nZ$cv(wV9sS2U5#HwX}Y#3R*#K#9aFe5${k#b##K$'zm,2?$Y'3Kg4v/1N&,(8.DzXQ#%J>#m@%<{:V%-.,#2##W-%@H'#?#RoWN-(c5%wx%(yXgwXfY#F,#qd<:x.5%-###;##1z$c4LP@#co3;h#,wX=Q#be/i9'[U;&##jm'?$%{6+ZP#&##"},
{227.01,846.58,5.97,1.43,"###*##D1:.,#AD?5,#rqdQ5#lpd$##s>%S,#y%1######$##3##ec&6U8.,#wh>8,#Oqd8,#-rd&##{c(<##5y6######%##)g.N-)Mn.$##[^:2,#GrdHc#^pd###Jv(;I#of6######0##/V3U5$2H&.##A94p>$c$S,##;qd###(d(8##$o1######$##"},
{113.82,864.16,6.18,2.84,"M##5_=######$##zm,######$##OG###################{?#HK_EH$&T2tA)lL_r&-g5%/P_]83bP#M,#sI.######'##-o#Sd*Gw(MB5V:9{m)dMSlI(QL_L#%4$&up$-h<######N##bc3W[*vG%%##'XCIe&y1>d,#vJ_*##7l$_.#8~/######9##"},
{472.52,144.71,7.25,1.68,"uG%######/##~U;######1##JkL######Kc#I#%######'w*=A1######-$#|Mm######v##MNm######?$$_~1######i93K97######[##&Nm######f##XNm######m##Kq=######;8+J'1mP$###+##<{d######(##iQmD>####'##Ri>}v+###iP#"},
{484.01,161.34,7.61,1.70,"############.,##########.,#######hY##########|Y$<ED######K##6<n######}5#R)B######V14#########j(<?9b######}##%=n######T##(<n######['/#########eM?mlQ######>##o=n######4##f<nD>####)R%E>#9Z%###(.+"},
{283.80,185.54,7.75,2.84,"`>#;+22nWfY#|xF9a@-I*o,$:W=(l####D##vP#Gu$######P#%*{(rnW^G#9oWml%s95r%$9mQEw,###[##D##;%,######pZ'qm<n~1%##hpWd8.eZ'd>#5=7G|EG>#hP#5n'+w*ec'%##3&-xb@~#&###ppWN?'2,#2u#]V9z5&YG#dI(@]/SG#a,%uP#"},
{129.52,188.45,7.55,1.63,"Lr>B_/yH*4##LODH.'ZsGEu#@0U-?%5?'6,#8K2######*##w-U3.&3C:d##8wReo&b-Uc#$x2UCv%+J/x#$E{795####-##@cP<##;o1^L-yU;l##-0UTi9r-U+##t$)@P3TU9###.,#*%#c@-###@5#8y,k>%###pG#-C2T08###&##t/(205$###l#b##"},
{129.52,188.45,7.55,4.89,".,#D##Xy6###(##3]()g7###wP#Q'0t#'###g>$FL,E~.$##.,#g$#WV>###J@)9G3XRU'##hTU:M6%<@C##HU8un(VRU.#####3##}14{k#M7-]u#MWULH%DRU9u##wLIw&Zp9'##ASU3w&###$##TJ-D>#L5$&##*:Rb>$UFIyP$TX?Wm&(J/zb#6{;k&*"},
{133.32,351.87,7.76,2.43,"QmLq6(Vv*P#$n.ErWE###R-%o*0AfQl>%Nl$r#'$$;OU:A['ZsDaH&xZ'~@)%_/m$O51:_G##hQ#ODxY$*@'z6+tH'.,#_:3y$)X.+KM70I)193.g+O2@-A&7dQVG#&Q%$B#Py6######~>#F>#Pu$Z1-mo5s5&l#%Zc#&gM4U8######Vz*p@-7#$###?##"},
{184.91,353.58,7.17,2.24,"_YKG-%uc(2R#SYHHy4DA0[u#rd&-K0gW@G$)}H&7C4+8/{r-Ug8z>$'I(-C,^p27P>|L=k5#$oRdn.Kd'Xh6pmRzY$8##9I>?['eT/Xh8t5%5q:B2,w]6yH%/oR;B4mP$}H#_I@eOJ###)H#95#p>$9;,5B5N%.GZ$ef*7V7G'4dT04T3zH$%k7QV;Km*G,#"},
{387.17,373.10,7.74,0.18,"y&*)M6)R*###Wp-*n,###3,#xR*^_;###2##'jVqv+###r5##C.#^,uJ3###kM0seV$##$##v_<Q/V###bJ&'jVc-*'###q.uu%fd(MX@}k#__;HT3/Q$z5%?fVx5&+##sd$cFE###sG#IS*Eu$&##I7M:5#)^2N,$B[)9u#R;TY5$2,#ZG#6F9O7,(l#=5#"},
{387.17,373.10,7.74,4.79,"/##hRS<5#9H&u?'&,GUG#tZ'Py+Y~0.,#T]-MV36.-###{&,}P$_c&K#$UU5+x0k5$96%#E?AwV$##'##J|3x[-###1,#(*6S5$1,#D5#${V0H&*##Zx%ZwVowV###^A'QL6J1:###8##X$(:I*pG$UG#zi<###wG#SK.N#L2n-&##[g/?|Vcm+####-#4{V"},
{203.08,402.75,7.08,5.38,"B15^b7Dm)VQ$2,#w}3+YI=l$>v%^|<d'6mu&&7+XG#UJ,LC6{.)|xGOU:&##c)A9=/ov){K0~g6r?(yf+_g5I_;rG$qw.|@,0m&|8.7$MV?(rZNHv&|~,U^46x.8&(`f19p)QS,9Q%7@'4eHB03mZ(h?%Xw-7w,,Z$zA-Sz9tuRmP#_.+0W(D/2d5$n%+_6@"},
{187.81,415.85,8.19,3.04,"%R'`?'tv%I/,ZQ%`T._@(.T,1w*{Z'I~&CC/n5%d&.[n'#T/^m%9y/jp1ju%wy,B;6OK*I@*8_V}>&c[#[`6,-&;5#l_%m5KVS-gI(gGIxQ(GD;KK.x)2l[*H]Vh5$9I(BC*+92HI&=@=M02Am$37(ne)=d)zR+kA-{['xw*f_VRo2ZP#NI$7J'wpOUR,.##"},
{412.73,418.59,7.99,6.01,"xQ)ZG#o,&,Q#Yh>pb#?5#-/$t.0.,#]>#pB/me~###8,#tR%qI-$Z$3l#[-)/f~M5$$##ql#(i?s,&X#%XQ#^e~0,#tb#K0#YZ'G>#TQ#j*EOM:g,&D##`r?bZM,J/95#u^)be~M>#mP$A)#9m)0,#oP#`K.,u####7##=NA^e~.,####Ns.~e~?##*6'('#"},
{354.41,452.19,7.37,1.61,"Up4(PA&H%###I2XnK6RG#~P#?>5MdUSG#f,$Z8+Ne.G>#:,#h_<^w+#Q$c&2t/XrQ(R,$)11^aB4e-V,$cKG`f3^P#Q>#un&jH);5#O@#y.XJ/X###U6#=q3t.X######Zx)BH'OG####]l#######-8%f83j}L###&~$mK-K.X95####d##}5&N#%###&##"},
{258.23,468.88,7.53,0.76,"vP#D^NL:8~>$+7&ejC/7)P$)(B/%Q%2##@g/)e,######s7&]H'_}@'jB]n.{[N4E=Gd*vU'd{8<g18-'u9-G(2q>%###%##4R)Xp21X>*K/#JFBYJL9+j95e|:K~NT#%wh1kGFVc&###;@#xo0I3=c['mo.'_8E/.Fe&SbEx@.LQ')##D^N<7,.,####:~$"},
{258.23,468.88,7.53,5.46,"=s=d['R].WT0@^2M|=A]/?[)vF@=|Bgw.hQ'PgNJL7hG$zP#7&.9~&15FR(9LbJ8'+j95|IF)<=Rm*)_'*eND=DEI)9m(|-&LQ')##PgNX%.PeNIu$j_1h|:.^18-'(C-m)90Z%3##]'0y8/.,####:~$0.,Vc&###8@#U5Fq>%###%##W:2######J~&'e,"},
{161.42,472.33,7.59,5.25,"###9B'381.,#;,#=gNeG$:l$h&*xcNB$)i5%+eN_?(p#'bH#vG%OL$&cN0,#[sFZV.5o.$)4.a>3I*Lp'ch9|dNx.,*m(HZ#T,%K##PgN@m)(dNa?&}g-{q9h//|)3rf/6I(<y,wrAQm)r6)###iY#7e$?d*/u#cl#xdK=06P_=0J&h:;AA*ye/1$'A&*>L5"},
{202.39,488.86,7.85,2.03,",x,hN63;>|?%vIOZ~._Q(=[#.~(`f/&e,sb#uo.&[(w5%Xu$Uo0}d*+8+cE5G)5FN;pHOO8*ILO`/2]$)r-(je0.,#/Z#`^4w^.4_8<1567+Qp3[+;$HN@U-oJO,04{?*Bx%I(+;~/###b,%5]2OZ%Gp-QD:So-)%)WT/}276g,wU5oHON.(6LNq[,'$(`,#"},
{295.32,580.37,7.09,6.10,"#l#qY#}b#A#$}>$sT495#5##:&+b./###0##~$*######dG#*##S-$w84###O%'TrT~ITh,%TrTRlL5Q%P$%g_<######F##,##Gm(303###P?'57&TrTUZKPnT.H%;[%K&EVV;%Q%###Z##,##naF@5####u,#-nTHH#r[,Jx&XuQ,##.A*s?$+bK###(##"},
{180.40,592.35,7.90,3.13,"###'##L-#b2A:,#j.+1/$M~U-@&1.,(##WY9AQ&######x$']G#Lc#kR-&Z$TI&5(MV~Ugc&Y`ULdO:-'z$&n(;######B##%H%B5#Qn,~,$cc&t?%Y`U*JU>~U/H%Jn&Y`U?1:######n##Uc$&d(aP#+c#Mu#*H%A-$Uw.oP$4,#/H$vJ0+I*######@,#"},
{180.40,592.35,7.90,5.23,"###Rv'.,#%##)##Xg1eY####+-$S/0ZP#'##'l#jY#TG#2u#4##;S(vm,###o$&WrT~cM66&x2Ue[Q,H%=%&Ao0VG#Sl%Gl$jD)Y(<x%.[P#8R(C.&$;P:/U7@R@?&4d%|(QrI-}l'*c$`,#x2UL@-(##<5#?h6eY#M##Mn.L>#ZG#b>$k7,gG$&l#G>#[u#"},
{180.40,592.35,7.90,6.24,"###B,#A[*###/H$c80*c$)##C-$>e.Bl#e,%VG#*c#^c$]6)###r##<1:###Vw&Y`U@~U/H%Y`U~%UXZ&j6%Rn,j5$$H%B5####G##Mh:###D6'%%&Y`Uw6OP~U~Z&VI&CCL_I-jG$iP#Mc####~@(q,&###9,#RQ>q-&(%,M/$1~U/##a%+<-#`q=###'##"},
{411.53,601.01,7.23,4.47,"sY#+}/I[+###v$N;U2ZP#:##*?P*Z$###L##_6&0I*###$##%d(Ew%i]NR>#@.U0H%+$&@@#`bMy>$,$&'@%B-&CR),6%FR*Su%###x2U>#$K/U###XT/_v%F.U###[5#A['pb####L,#`f2.,####e43O-U[08###eJ&H1UAy7###6##YGB95####0##{o5"},
{412.29,679.32,7.13,4.09,"_S[######P%#x=LOG####(%#Xe/95####n##=6(&##P#%l,#-T[######xw$1PH=-'C#$NH$Kr;EH&;5#4Z##m'?)5qb#A##NS[######ez)qWE=5#9c#|V5cR,Y?$Yf2:p1(l#k8B]S[/,#WT6###U&#JW[je1###S##8U/g,&0##>`=Wv(6#$Mv#YU[S>#"},
{465.68,691.43,8.97,5.99,"MQ%v>%J5#V?'C~.v5&###CZ$#hg0,#{k#k/#;4JC##:/2($#s,%]:8<5#;5#MkE3)@###9##khgt,&I#%1'#|L=|5%DI,b##2[#1=C######ofOJB6###)##Mog.IP.,#m##*R(<_8OG#(##Km%N-)###%##{^-ac'######^[=1^7###$##3l#'d(######"},
{310.15,705.23,7.15,5.83,"o>#t>%(##Zl%8##{x0cw#}Z)n$&2w,M9#k[Vl5%###%'#f[V}l$|$&HI,###h[&__P_>L>c$v`V)ZJo5%Qd%9_5du&?##-d(TQ%.I(1o0###D6'U~(v`VnXBk~VQc%sH%em<)L3)R*###F##/##F8-EH&.,#j6%i81R6%Rm)Uf*=R+)##5@'+m#NK6###%##"},
{117.57,759.50,7.60,5.87,"###B%#ay9###C:;a/#V+Lw##$lMGH'pb#2.#u#%.6'a>#Wu$=#$lg#R5P###G6P>(')]3Q##m3?lmG5?'?##<2)9)?D>####-R*@C#R5P6##.6Pe&+xd-l##r7+Fr.gK7PG#lr,Ng8T#%J#%Bm(&9$R5P)##]A/S5@r5&1##BB+3bDGc%sG$P%(%d(,m%l[,"},
{212.10,782.57,8.41,5.68,"xM.WR+(]+ES/<l]]H(}-*U5$u*?al&7##V?(]Z%PI,'##|#&;h41N=+A.EQ$1h]Wm'5J/H%#mf]J.*D>#Z##FQ$kI,######wd*qR,;7%BdOJh]IQ&Sl#?=7c>E~K5###@##R@%^Z'######OT%1<D%I#:5JiW3pDB<##t39FT)M/3###<##|7%0Z%######"},
{118.03,791.08,6.96,3.89,"R?#og6^n*r#'zh#B]4[l%^>$ND#Gf4###PG#''#iZ(######xB*sIT,H$x#'u7JWH()Z#bS-n{&r@/5%#=T3n{&B&3G##fG$_~.sP$7B#rFDCfZD>#9##)V+q(9OG#6&#n1;'W.+u#>&#fq<-Q%k5$a$%a)=]R)R$'Lc%<y--}6-H&4##Bn+'16###3%#PU7"},
{409.00,810.36,8.51,1.40,"Dc$;c$<5#5l#?p_.,####P5#k)j######0##.%-######%##J,$95#A6#'`>~eW/,#U-%Ld(A+j###$##[>#`&4###############9:#)h<.sE###<0%Ap1z(j###(##f.%yo6######(########TC%z,'(o.|u&s])D14}*jyG%0,#`x&*80######)##"},
{203.79,831.68,7.61,1.08,"3##;bEZu%###b6$y06du&###]u%>/$c#R>,#c#R%-#YT68U#a?'*ZJ-v&n>$vRMq&2uc(`5#j]5#o#q#R..#b#R?##0R*1(#Kd(7['(A.@u#BlG/d%Cv)Y5#|#R~I$OrB>K#b#R###fY#1(#xA.e)>N5$Zl#vaE56%?R+|,#d#Rf,#8.-]/#`:=######C$#"},
{203.79,831.68,7.61,4.48,"95#7##U7Y###i98c>#<%WE5#:Q&0##es:%S+j-*###}P#tM8fY#j##N7Y###>PLa,#08Ypl$<m)V>#J>?.-&]w*/H%3v'4m&vP$O%#H7Y###F8YI$#^,PvP#R/3)##hO:XB4lG$VG#`?&l~I3Z%/-#K7Y###e7Y[5#liDV5#4p7###F##[.,Z-)###(##78+"},
{41.93,53.56,10.51,3.03,"K>#KJ)ku%N#%-5G}I,Fl#b6':+dD.&1T41##sH(9S$a)C###6,#w&-5A*uG%rmHh8/a#%hZ&%+d<H$*L9d5#h?)CQ#Y2B######)##IN'okNZq=%##vv#?`bv)d*R'=~/f8('v&^Q%M/3#########)%#l2C%?&###dJ$[)Ar;9I$)xP#D6&tG#PH'D>####"},
{362.94,59.67,9.43,5.82,"vd#]o4@5#}k#Lc#g'8nY#mY#cc#%m(2c$L>#$#####G,$yb#w1#6wW######-()4yW###+##cyINq=;5#x>#5d(###G#$LZ%q%%lo5I##n820&Tm*E<###@%a{WSu%###)8#(B.95#J,#Z?'OG####r]#zQNbM<.,#*p#h/0-|W###$##x%%cV2OG#&##BH$"},
{388.22,66.57,9.58,5.72,"S&#1IT######3|+pYND>#&##[n)7#$U5$yP$###'##g#&2u#YU-HZP4##z,%xLT11:###v&#x>K###2##,~&######Ed&V,%505A,$)q%hS/|MT.,#/##RS$EKH###-##bQ$|l%###.d#wG%nP$###`))}U7_:;###J##<xDHL6###/##%k<?7+###f>#>I("},
{315.20,76.06,8.77,4.56,"8##n=@M5$###u4<L81###&##{eY######M%#;y7###$##*J#=,#8d'5L+(c$CpS?c%8l#S5#;gY######x$#{B9$##G,$fv####O##v1N^#&4IT&##J7*hR$XhY######%%#?V8yb#M5$L5####z%#0,M<?&Mp30##}c(hM:=jY######OK-`XD0,#$##JZ#"},
{209.20,187.92,9.31,1.88,"bq<0,#6##eO3$5C)f0###*[##L)R~0######H['k,&###'##6`B`>#fd+MT'1B,GgKIRSzc$nVS8+E#$'xu#.i:&c#k5%A5#rcR'##I]24?#2%,QH$nVSa+DpRStG$V['Ym=eiA###?5#q>#,_:Su#IZ&?##jI*rm&s~/4%,Qh2TZ%=w-gd),D2E>####$##"},
{482.99,247.47,8.76,4.59,"SiC%##an.7f$y9e######?+(:IUQu%###'8C?,#=6()##v=Kt:eD##)@+W##Y:e######4$#';c9x2###8##b##V1=######A|C4##L5$%##O:e95####@##Q2=+K4###$##=##[C<#################################OG#######$##L5$######"},
{218.18,368.29,9.06,2.15,"Yl$_'4[>#CH&xeFb(=###~m#o/-;oG$)>aH%KR'>C,BC9Qc%mc%Qc;*6'+##mgO]W@Ql%)H#vJ-4mAQg7t,$O7)bB/r'8q$'@6(Wd$vf+%f-mZO7B0lu&_%%I017M7px3X.+@w+?%+/R){l<3u#o,%t29h6)^d*jg-4w+nl$[dOy[*+c$%e&1fJK..$##QD*"},
{356.67,384.68,8.82,3.03,"###He+40U(c$;[(8XD~T2D5#+/URy8$l#^7#4K1VZ'###g~$###fo$3.UOG#2x1qH#uiAwG$b.U?5#N#%M6#hD:W>$###q5####3o%4,@)#ODv)O,#++2,=DQ/Uqb#Du#V:-@h.z,'###+S)###7,#]],1_=W>$###JAD`q5kH).,#],$-[;Eu$0,#5,#d{:"},
{407.55,384.90,9.16,0.07,"/f)..+ZP####$R&E96###*##$|]4-(###-6#b9[###g>#HV4a@'z/4xG%)##hV<Ug9###g-$t{]1Z%###FS'9#^$##Fu$fP#|$+$l#17+6##6]XD>####g,#T%SK#%$##ru$/|]D>#bG$%##>c$0,#Zv)0~-<V6ZP####~~,AO?(@+###S$&3y]6#$pb#K##"},
{390.71,426.04,9.71,1.30,"w9-'8ZW>$$##e/%)8Z3l$&##^r&M8Z;c%###(r%Z7Z&##OG#pT&AE@=6(###lR%V8ZA,$$##=~P};=T,%5##U):@S0###(##,?#Y@+FH%j7/e&2/^6>##Jd(d9ZDv)$##h[&2{1ll'###Sc$###.,#5-#v<FYI-[P#9n#0{:fq>OG#+##{z.U$)OG####>&+"},
{390.71,426.04,9.71,6.04,"95#36#1=G$##Cw+ac%mI.-?#A*BH?(###'g&1/[A,$$##zT,fY#)e#5i<}v,u97;##Dd(No1j/[W>$$##PI%5/[I#%%##9x$OG#-##w1./r@e6*$##a[&j0[vD<T,%0##z.QT/[;c%###})'OG####8x*Q$)SH(###mu$m_1-S/###&##]D9b.[%##D>#Gq%"},
{226.98,437.10,8.53,5.47,":<8_q4;o2mY#Cp5/N,(kEm[*f1:q6&wR*(y-A&0kP#c7*/t8Ee*4*@'I(Lf/T19im*~/,8a?I5Mfc%1/)7P6+QNXc$V~+lZ=X<@Om%g~,:^1UT/G*=b~-u.-CTN3h:oQ'd[)8UNUGBLQ'H##c-)=?$$i6~B5}(>,@&qJ/m{6['5j>%&w#ZYEN~-=.)5m(^Q%"},
{181.62,441.00,8.66,2.75,",K28R'?m'OA,s.,pf.6`6N?&vE<X.+M~+=d%u~)Q;;1L/+Z$pZ%~T/Y]/@w+HS+r]/Be*Mi9jJ0%]/wg.|C5aI,^13xp2re)Cd'Ap/7_5#$&]9TJR)?A+#{+{]6*l#6a-@mG&w+?c$oK({4;g$)C/.?x-qZ%y6T>H%&['a;(w2Aa-&z%C=48x80s,#b&O$c#"},
{181.62,441.00,8.66,3.88,"Lc#8t=Lf3k#&H[%g.*?p23E>]*=p,$dy4*f-}*==n.$##'6%#R%kx0+B-/N?M=@>7*zI+]]/UB3VQ&c?(Y`UX%-Bv)|Z#G_UPn*]n*(o.wA0xg7u:8R04U//m.A{kLd>$@y.Sv;GC8Fv&$a;?]2jZ&L-'pK2,Q$p['S29rx3^B'lR*C/-<D<Da-.W?SZ&@#$"},
{181.62,441.00,8.66,5.10,"nu#Z%RE5#GQ&n?Kf19(?&7v#i%Rjy4-Z$/r*V#%r59nO?Z@,]96Q.*>n&i:4o%Rpv)>[%CJ)3ZBp(RJw.I,#y?)p(RJ?&<@(#A.M9,8'2oR*pB,eh9{c%x~.pw+AL.sD;`8/-'5'@'q[(2^/?~/XZ#qsBh[(P]3.%(gw)kT0s&0j-)?7(HA,RS0Uu$a[(-{."},
{151.48,491.47,8.45,2.17,"dh9le+|q>x#$iUS{n1[P#?H$0&&DI,###y>$pG$~P#I>#K#$v@,dR)?~+kM2k|31tCT4G(7&nVS?&2Qu%H##3K3######GH#UR%dp5[>F~R+i]1R`7ejFWi0HRSWG##-'Pq$P97ZG#.,#1-#)$(o>$hF@B3@K#%R#%R5#nVSHp82u####,*,A$)uS'95#>##"},
{163.73,639.91,9.59,4.06,"@u#Hh+>u$$##w2<Gg1:5#T>#,;>.,#R##@/.6u#pb#Y%#rT7pm']=0zB9'##^qV/B/D>#Q##_-TD>####x5#8$&ZP#=##Uu%A[)4{%Y]5%>8nnVGc$cG$hF0qr?=E;yu'&6#1B&.05.,####Yl&;##?u#YR=o#',##e?(sO0eH&|y/-$E#f.hj<n@-hu$T0-"},
{324.79,662.37,7.11,1.32,"5Z##;.VrB0,#Kn>:L6L5$+##zD8@5#xY$&##/,#6##@d*###m?)ZU%R@Y$##GBYb$'E-)+$#1AY|G$A,$[##:5#]d'iZ(###Dw-2&#*AYEu#tAY>##'?&{Q#iQNce(5-(9##6Z#[?Bw?*0,#}>&)##_HGSS*d-V###aG#k)'&B4d5#[&/;v#&Q$hv&/;3wn/"},
{205.49,687.89,8.86,1.59,"EQ%/9-=)<pG$2>4S4F0Z%*##mT/u>%)$('#####/##GR+###h#&O{%sRY$##2UY59,o.0e$#lmTNH$eu&q##fG#cc%xu&###Hm)kB#sUY<H&0UY2H#m6(uv&4.RyP#8@*&c#ll%t>$~?'7l$eY#$##$Q48n-yQT###o%$S`6L,O###`>#j?%ZP#2##n.+>?'"},
{345.17,822.87,9.69,1.36,"######__#o.0EL:/,#9g#|'1$f]###'##;B%*A0######-########&F1{k#;D6P@)f.M+w(Uj]2?&Bu$v$$:o2######)##.,#%##=3=<5#7C9'c#>i]4J'Uf]TG#|l&jC%~J2######:#####2##,93###R7+p5%Zz4L(0g-R{k#C5#LV'($(######0##"},
{132.02,844.08,9.23,1.94,"P%$&'4jQ'|k#bX1KR+8['R5$YQ:###)##&l#&###########xv'}W-dK5b5$K1T}.,W#%Om%[38bA-m?)($&,Z#l,&###*##PQ'qm#~r0oOD4.T5['=c$+}2YR*n&+oR-j&GxP$sP$.,#0`.###jH#p7)EY>%Z$:5#2l#[2T;5#.,####sgM#########DT+"},
{309.74,201.59,10.39,0.03,"######eP7SH(jv+###5#3O)3kRY###)##x2'@^9###I6#f^0eY####fL,:_<O195,#E=YDJ,oTYG>#|#&b1&(^6Fc%Hc#=f*.,####%n(,%+nH)###>oBV^5VTYW>$%[$t)2I:2W~0G>#4,#$##V>#~Z&yP$8l$/l#18*}d+_A1.,#nG#pK/De+6#$###`5#"},
{457.26,300.28,10.36,1.45,"###'##JbK{k#pb#&##IiA8,#A7-###^>$G-#o#'######q?####%$#(7Upb#<m)$$##7UA,#h6U+##{5&t/#qq>+u#$##|h&###R-;,#O###RI)53,y^<###)<U(K3pb#D###>:/7U###8H#%###eG95####c5#eB1######Je#CC:OG####/-#&8Ug,&###"},
{408.11,387.27,9.78,0.03,"b%&T@+D>####;n'X83###(##a2_Fl%###fH#Z7U###z>#1E7I@&k..06''##B*?^~1###5$#%3_L5$###m-#s3_&##Tl%/l#*@*G>#Fe.O>#|~W.,####^,#k1_pb####V##<2_&##g,&%##S5$%##W$(|83xT3.,####U9/FAUL5$###,%&00_###uG%O##"},
{331.86,406.72,11.66,3.10,"2m#dfIZ-*###PV''-QD>#'##6/'WlPOG#.##Bl#YA+^m)qb#RZ%NGC6.Q$##T/Q7=HQ#%N$#5uFed,###@%#R:9-Z$)$%A/,<Q&U$#f.Q)w,l.QC5#fQ&<7'@6AuG%###QZ#(?B?v(@5#2?$95#8##V*/OWC7/0###Rw%-#7[.+PG#%##SN4:v&>m&@R+rc%"},
{336.07,640.41,10.56,1.15,"%##Pu#i816#$aP#vZ#uITh6)[J0{Y#+B0F7?qR-YG#QH&;x&$m#I'//&2###[I(@X.jHT.l#XMTk~-[,%qm%KR*M,#l$+%Z#4T-[T)*2?$R&m;Aah%iHT~Z#]IT4I'(c$l.#0.+3g1~#&.##OG#X##:$R<?&-o2h##<IObD.PtK:-%rv)Py$@d&XL4Ie+$R("},
{336.07,640.41,10.56,4.79,"$@&x-)Qz;=5#+(/%g+MnV$##_rVA?&NQ&mY#_R*DA++]-oS/y81:Z%BU2]T)8L8kP#3sVzu%AoV:5#-.(k6&)_;###^l#>z6SI,###o?'Om%q#'-##{~CEy5NoV)v%U8)Q;7-lB)c$###9,#w%+{k#%#####(H$V.$vq:k>%(D1c<0SD=Au$n@H>-(0,#^P#"},
{226.32,669.93,10.61,1.49,"OG####~,Cd,%kC;,c$1V0jg,;V9yG%MQ$6R%0,#E,#g6'3%,###Z,#E+F.,#rS*V>9kcOa5$yKIwy7L5$A##(d'#Q#_m+%##=-%mo%<$T%##<FF{M'KRVb5#USVc?&nP$+%#`m*>$&>Q&-##K[)Yc#EdGhH%q1>>-#yKI}x/SRV)##,c#w.'nc'|Y#6m&`,%"},
{226.32,669.93,10.61,5.15,")##N$'l$K###WV.F5I(cBI>#jWV?n.)6%=l$g7-%##F-%Yp5GH&)Q%s_2;u#9U8SG#jWVc?'GUV95#^v&).'1`395####F>#rH(.,#uZ&r5$Tc&###N[<wiBc{=o#&[]+qdL6`'_/2lc'H>#Yd).,#$##&##b6({5&V,#Hv)T>#fQ'AGCWv*z/&xo1dPF.,#"},
{411.07,705.61,10.25,4.22,"()':$)######KJHpb####(##R;@###$##[c#CZ&*##46'Lv&H)4_?)/,#o5$<UX95####K##`mT######G##w$+E?%^c&D,#4e.###$##D2.zRX######F'$*,K/,#SG#k$%`l%Gy(-19i>$r5&###)##IXXSS1###g$#IXX}I/###O[#(~N###>##EAQk?)"},
{410.29,808.40,9.94,1.41,"<u#.w'fG$J>#ip]V,%###,##wVf######$##781#########rb#VG#~R$|(>ZK_###S##mn-YWf######B,#Ig8###############W(#n`E<rA###a&#gD9cUf###%##an$D1<######)########wp%wu'X~.F#$h1.5S+dWfvb##l#)@$px4######&##"},
{345.53,131.94,13.37,0.69,")##@r(|e2###-A$qmONc&a>#O[Ev:>###>Q#^IU95####G;&0H&4;%7817,#fV>+7'L5$~C(D,JL5$###ia1ANU######q%$*6'i##:9.K81@JU$##B,#>A*mbE######0.'E#6sb####1Q$######M-#MC;z@-###-$#A(:_|6eY#=##z,&$%&+v'###2,#"},
{396.23,189.89,12.72,6.26,"Bv)###2##ax+Rd+###M##iq4O]5###$##L.;.,####q##zf]mw0###'##{?$&>M###1##9e&$f]###%##r#4^m+###QR%-g]3/0{k#$##{P#57T6#$$##P>#^h].,####:,#7L7a5$RI,_P#`u%.,####yY#?(7.,####'##<l]:$)######(g*+cMbG$###"},
{253.14,401.49,13.20,6.14,"m?)2c$-g*_93f~/Cm(S8'Xr:L,Mx>%kP#,1.iu&09/|k#W#%5A/|l#:RN/@(EQN(l#Wm$8UNL[*.$%7xGy$K:5#yT,Y5E$##48/jw(C3B^$%8UN>h:d,%wd$@B*_a:%RN@H%###m$#`RNW>$(~,%v%q6(/7)rg63w+>-'Q~'=T11-';v%r^..,#P##lr/}NF"},
{155.02,402.33,14.04,3.04,"]H'8m%7J-e6(zI,2['w%,]f/.uH8H%s[)I^-}g8&Z#Sy8a##kd+(w$mcK5-'svQGH%V&):IBx^6P-%:nA/ZE&K-gu$Jy7W5#|n+:)61(:4,#dzQ^wQ#R(YQ$4E3dzQq3C,H$_K19u#f6)Su#Q[&/%*h[)zn/ce%72;%w)W~/wlAeU7C6'Od&2=>OG####H##"},
{314.58,609.86,12.96,4.73,"Dw*+I&mq4R&09&G1d'5l$'Z#.S)('3fY#dP#C$&BH'Dc#qP$g7).K)B&V3c#.*Vo7*O'5o5$pN?,g0^8,bi;+16<Q&A,#xB-P?(0,#.*V=e)7%V0,#sW5{g/L%V.,#R,#7W93c$.,#rG#h92[Q&jc#$V1Yu%ia<<1-xq6_7+}KMuQ)###R>#rR'6?'L,$R,$"},
{168.33,628.22,13.45,4.13,"R-('L+Nu$:-$Gi:%'/.?$=w+1x/I#%~$#VOF4$#A7-o%#%g747(uj.bV@q5$OWWo13;c%G[#MHJ[~0:##fQ'LH#Nc&1.#AS0>?'tg#MGHX[<}dVM$&>e+t+2qE8]}C)L1=B.xq2S6)$##XH$,)>]?&ym(~h-e%0###I-%Cp'>Q&(c#P~'T{:lc'%.)bG$|-&"},
{403.81,626.17,12.94,1.00,"g$#[y7ZP#####;#(TWYl&###gD&dRWeY#B,$9M#BRW###~P#'##c-'`%-.,#,['@{+.+I9u#s;Tap5@6(Fc#5W'6SW3l$$##)##*R$[B6###XH'8I$HFE500gy8hP#r#&VD/h5%fv%Bg8B5#n5$0y&A7-2,#Ci<tn(YZ'??$>J/eG#(?&gl$###E##+M=###"},
{403.81,626.17,12.94,4.46,"M@,$##1u#Y,$Z>$sG#.*;(-&jM@}b#}-(ve'R@-###B,#-n*A'3mP##Z$8,#$$&zg(Mz9_P#IJTMw))6&3c#}o1W>$###tY#1T3%##Pc$`|?~Z'|G#,8CaL:tJ^eP#l@(*24,n,OG#%##3f-6#$###I$#_De.,#<5#Sh$0Ce'$(gY#<J#sDe+u#.,#<$#UCe"},
{303.53,643.53,12.96,1.42,"%##Id(DL6+u#~m&>o):TWS$'wvAxv)ll&iI%3H%'Z#LB3?5#ym&P,Fq83h#%to0G;*sRWyY#1WWR~)mJ20Q#Ge+H<8am*B#$/v&tu%c01~186@+m>#jWW-o+/SW}Y#m&-X`+p6)Y$'^/*<y3*##Pc%67+xY$N,$;S&/8TGl$fk9j92{-*&m$nI)@l#VH'}Y$"},
{303.53,643.53,12.96,4.59,"G?'Q5$ES*^,$)n++v$Gk8[^3Y~SEl$3u#/n%v-*`5%(##-H%8o)hB45R)T6&]o,Li+,JWzY#zNWwe+8@+sG#,^0d18V6'vc&}$+qG$Iw*0|7Y828Z#GNWBS)mIWwY#7'1x:*~f4-Q$CI&/lEhT32,#_#%Al#HQ&,%%ov@sv)=KW.I'f['[8*pC8A,$$##kH'"},
{373.33,649.73,12.60,1.01,"(##d['*9295#'/-v;+:ZP:Z$#_PX05%?&dl#}C&D[P}>&###(##Dn$z*G.,#M.,5%%a?Np*;H)@2l#0$'B+3<l$`7&zjI@5#yl%c;*Z'9)##`~Pg'-?$)_[%I'6Gc#Kd*}Z%###a##FZP###=m)4f#eaG@I&;ZPmv'6H&b]$wR+ZC7Il%B,#;,#lQ%AL8###"},
{373.33,649.73,12.60,1.96,"W$#J;SrI,###60&LQEi$*/,#P,$yY#Or,9m)######[1$~'9.d$ce,s7B>H&V:S[I-'w%?A)hT3Hl%#I$A$($##8,#}T'/v(.,#*##};*p@.O92%B1846%o,};*.2=M$)'##9v#KS'he/]P####rQ%p[(3l$5##j$(N:S###mX19v'tN:Lu$aV4:Q$0R(;m("},
{288.67,677.50,13.05,1.66,"gA)Rf.~nT7?%#S@SS-6v(1##bl$,1.xA1###gP#?$'Ao+0H&X~+dv'-pTVu$TrTyZ%6_7tm&$y1q'.)C0SU5ul&nl&SQ%yJ.Rl%~P#/KF+S+YnTsv*f(.V(/rL2*7+O[';J.%##O,#Xn.1-'%##=Z$u`>fY#U9.zN4&nTvb#{qT|S0Fl%###X#%g,$t6*###"},
{141.19,813.80,12.33,1.96,"$##4L%[U:###R6$1E:GR*.,#ga1:7,_P####{I&.,##########5.#6RQ.,#fQ%Qz.$*@nG$%^Oa80iG$x5$PA+A[*###Rc$###.##$UN^Q(7l$&?#sY7$E>1V<9u#o#$?THCc%Q#%$##{[A###Y-$AD7cA0###)##7A+4VQ######,##4VQ#########mV0"},
{372.04,839.12,11.82,1.39,"2S/$##o-#]::}o^###.##Jd%Ko^######+##D>##########fu&###{2%Hz9do^$##Uc#hz*Bo^######S##D>##########yP$K>#MRCIu$qq^,l#Cx/s$%uo^######A##D>##########E>#$##*a<gu%bJ]/,#1e)$i)i*H######7$#############"},
{316.12,45.87,16.32,4.55,"r,$[G>eY####b[Oje/###:$##YK.,####&/#?u$.,####/##H/./(2n*5:l$[L[IH'T,$`-#N5O###>#$od#.,####F>#/,#Wv*c##qN[W8.TK[###zZ'dE0$.T###^P#wl#######/,####.,#2##F[9MILU7.(##~5$YqL_m+:,#}k#+.&############"},
{55.93,116.64,14.90,2.80,"}P$<c$9#$$##~f1J#%###0##:'b0,#KZ&/:#F'8N##aYM,w#I5#Vf-QG####z}Dsq;###1##v(b4[&AR+F8#9f2U[$|)Dp5#&###p*5J+OG#D(bl&0kG$e5#?)bVS,pb#]##pu%B.(xY$######iP#T3)381a[V)Z$X.%9N4d'bxP$###n-#&##_P#######"},
{342.23,195.09,13.50,0.20,"@H'?,#c.&V.)S^:###%##d`.>)A######Y}7/9`######J,0VZ'###MA#2S,CdU###7##B~%}$T.,####fm$J?`######>d#L5$###[w#vH)R<C###y5#FI(W)8(c$###,##W?`OG#############>Z#r>%=Q&###E,#$.(+@).,####n>#n|2{k####*##"},
{440.92,609.53,16.29,6.14,"wP$TZ%###nG$?i>w,&###vY##Xo:,#Nc&4##`e08##Uv*###d,${i>b>$G>#u#EE=G^P#%Q#zYoT,%9#$a,#Vo4$##|Y$G5#1Z$~R*om%?B1W>LW%.>?$P|4(XoOG####m@#`K7######6##5,#%o-_G#Ul%lM7u^8A5#kl%{]odu&###/##V7,D>#######"},
{206.31,623.73,15.27,4.91,"'f(Q92Sm*Qm&Kr(/L6ZP#gm&5@&k}EO5$)-&p6'-%*tR-G>#h.&Gi6WSVaP#jWV=i<2x/I$&Z*@='1@7(e*<s/26?$^R+C[%SQ'@5#jWV2S,zSVub#*q0NE2c?Jw5&)##rS-b,$M-)D-%1v&5-'wZ(o&*jm+bC.lA-OYD@A.u?;<]0Fc%?l$+9/,f0}b#mR("},
{189.28,660.39,14.52,1.58,"$##C#$X00z-+DR(QS)h&Og-($?9.B1^c&_Z%(v%1c#5.,4l$CH$%3;_/3A5#~I*l*,R&YO,$.)YL],mA2>Q$-^51d$`-*4,#uZ&oR*UB.^aA.v'0?$QjWlL7G%YA5#ng)E?<Eo3}5#,sB@-$###fx2#Q#?&2###KT-wV/X6)o?)>e*'f)|//%6&EJ+|g:Yc%"},
{189.28,660.39,14.52,4.77,"xf5T?'SZ&*o,H]*%00J$)<e)9;-6R*###A],xG#p]5###^e.oECL6$BK5f5#<()_??3SY@5#(sTh(9F-(LZ$z]/(a>M$'LS,$%,@5#pf3,[$/f1Qc$IWYnA,ATYX5$M%+I3,RK5>5#e#$I)9fI-A,$BH%.c#DQ&/?%[,8df2R~Jw$):e*L8+L11qm+/,#XG#"},
{102.76,454.55,17.77,3.11,"=@$ibM>,#9-'ay'^97rB80c#BV;>u#xGN,##?OF###E>#(##5I&%HNoG$ZJ.,LNNT4)w,R/,?jC:l#'HNV,#VGN####l#8##FHNel%=w*bHB+r<F%)f:4J?GAC8CZ$-)>9-%}<G###tb#BQ#jH?>HNqu&1Q$(a1kINDw,e>$G;8;#$i>$Vd%]m+###(##d-$"},
{107.03,503.61,18.13,2.45,"7M'Yh>E>####e^,u,&=S)6f/B$)=7#m[+m;1=6(W/#&19P[%;h'nw/Z6$(&.*zO.e+7u#t$%bx1,q$F'8FQ#&(;AA#,d)}&#9w+F['XE.aeIisF<&%z-)[W-8S//($DT5rg%fS26##eY#h0#~W84q547*UQ$wH)_%%f$*Zw'0e.Z##3l$fU$(c$######w##"},
{107.03,503.61,18.13,6.21,"]P#5H#^M@###%:82$%cK5-H$(z2D4@Qg5pn+N.+NO:BRQAv'D,$G##6QQ###~QQc,#|XGSu#Y[+JL14VQ605&?%iK2ER'_PLY>$0##]QQ$##UQQ3##a<B,H$Pg9DZ$[i+EEC@,#K?'t.$XQQOG####L3C$##OQQ1,#Ah:?,#`^9,f(M{@G>#4L0cB2aR+a#%"},
{131.98,339.42,19.06,1.81,"hl%)$%N]4'##%_0x%,yM@.##BY=i^7>u#i;3nu%<I'Qx(3h9Q[*E>#'y2zY#Fr?0,#zeQMm';1TnZ'*$&si/b7*WR)mO9[_:F8/###gc%?%(ZH(###S9+d1TeW@|Y$^>#48CpH'[I+56%ew,l>%.##/l#[2TD>#X,#36%w1T###2##G$%[YG###1,#VZ#QU9"},
{460.18,358.18,19.58,0.01,"ii995####;,#KK?34I###K?%dc#peZ5,#L$)############p;8)e,D>#/##O{d_5%*##D,5?w-1H&z,#nfW############xGHJ$(^-)$H$e|d95#%##G$%a$Sc5$u>$;n+############[83###nc$zL4Mzd/,####qv#)-P=l#uG%###############"},
{135.39,626.98,22.92,3.87,"lL5%Z$(c#zY9p46fG$I5#h~(jU%pb####xb#U##OG#######%$&b|4E@,}i85JN6i;=l#C^-W_.Qu%{@#S'4w&#o?*6H#95#X@,IC*COCdY=[INa6'<6$HwC#j:VH'jx#tm)9w#_[+~:&mP$qEE=v&#?%-X1Iw-&-$hm'U(5h-)0,8Lc%0$&8/)fN>U-%6Z#"},
{144.75,716.95,20.30,4.78,"*##R6P######^5#$$Q######h>$c>M###;,#ub#_L8###:[(b9-#&Q[,%C5#.'Q8uNqb#~u#*(2|V;0H&oH%VC-E/2###DT-np7Z?'On(1D/<$QlG##R)uf#e:2Ei2b$Q3H$WyLA'6_u%II%#J/###,@&uy*o09D>#FQ%Nx#Wm*vG$P?>e@*87OBl$b?&]J("},
{63.95,57.89,24.53,2.79,"QL4c#&D,#{H(O0/nT4.,#%##w`fj,&9,#yx#62?###LH$5~%qG$1Z$:a&mV@'<<q(4>+3(w*&ef8$(oP#$7%LW<###&##:,#######OT#6.-X7.3,#E4'.=AL_f###$6#U)06m(#####################D>####*##ub#'?&###%##fP#############"},
{420.26,293.02,25.70,5.91,"E>####nu#J?'yz?###K5#Bq)@RV95#u##+hL9#$:5#W/#vGMy>%i#&6,#_Z&(xWA,$###)@$~|W%95B,#(8&9d$M..p5#PZ&+I(;S*R7-(c#*xWz5&0,#p)(#F;9@,5$#U@DI##*6'.##6#$:m)H5#z%*r&.uvW/,#3,#h2'^C<###J##Qy,######%##.,#"},
{306.24,327.62,24.21,4.07,":,#gd''T0M#%)?$F$(U~,rR-8l#xA+>A/*H$=v'tc%95#A6#y?$yY?)00$Z$:_.sK6z?)zG$:l$(Z#uS-8n+ZP#5Q#b5%K6%i#%CQ%uO,G]^X'57#$,/$^GD=c%DZ#l~-vu&###|J+9[)95#95#,%#G)NC)@OG#7##fP0E_^###'##c3(.]^###&R&X6B;c%"},
{370.30,686.12,27.18,4.19,"Su$dR$?n.0,#~S/_,$p>%gH%'H%z,$6A/^?&O~.w>$;$(kw&?A*8V)=)A'##TAQSc%5c$5H$xd,YH$q]38I&jB7fP#lG$se'q$*b6$HNA$s5lbKqP#c%&NFY]o4ou#JL-f4BvR..##LZ%z9.######TG#i&-######3:#BCYI#%###a_#6AY(c$0##.3({@Y"}
};
siftPoint_save sps9[] = {
{110.19,79.06,4.37,1.72,"YG#0$%4d)/,#Qu$e$%Ho3+##^x4:Q#4o0j)-(c$W,#sYCxx/F>#0$#`~/}k#`B67M%nI.O5#?nUu[%|rCK_$Q$){A$imUgQ$###3&#vC=0l#mT6go#qe1+_-kpU:O6&B4zI'CH%n2*>mU%##H?(.%#.z;(R#_@.O##18/2}/}P$=~&.}G#.)<5#A.$cA3$##"},
{109.92,89.47,3.63,1.76,"W5$#m#wn1&##K(=U,#T[+vq(:$)K##x8Qt26Y,%T5#zB5dG#F%-&y$ed,)##EK]lZ$Pn.ly#K]4JJ#3K]sd$kR*Sv$|1?$##1o1K&$,%-Xm%iN]m90%m(n?$yw+]b/PJ]*##Y.*r#$;:;###z,'?##gZ'~E/~7+Sw)AT2'],J5#L:)3C:####6%c,#N@-###"},
{109.92,89.47,3.63,5.07,"_83###WG#W##[D@/,#1c#uf(hR*11,3f/oH'Q~0co%iZ(a##Y{A###5v$>l#>o]&##Eo+*b/iH)}5#Hu]c&.?07j>#{@-H/$(_<###-S)Yl$(q];6$x~1GQ#lf328$hp]e#$DT5(##7J.^7$5n,;,#uQ(tY#whY@K/du&<,#.n*Rr,PsG+##:g7/##Uu%W5#"},
{485.03,96.46,3.97,2.95,"em'{;APc#~H''{:'[(DZ#e:5a`:F>#H#$Qk;A#$0,#)R&xv)8@+1Z%KA#mmS6q6H.)5@%M91MqS9H%<#$*@%FI*=5#hH&.l#######S'#slS7l$$c#V($wmS7L9<5#H,#]~C6?'###Mc#+m%######C(#rlS######C(#vlS######S%#>mS######av#P~0"},
{485.03,96.46,3.97,4.59,"9v'dP#/,#r?&g>$q$*;5#WH%{Q&LQ'###>Z#op:######8R#xW5c)7F>#n>$m-%:`U+?%G,$vSDXU:<5#C,#R[U######:%#s;9h)<yu&,Q#'00zz7%e(8[%7~Ue5%(c#py#,[U######((#am*z6(hr@~c#7RQl7.zu'H@#,[U.,####_&#'[U######|'#"},
{110.91,100.12,4.07,5.03,"eY#'##4]*)?&a1=###>Q$:6#$aEyb#M@(I*/tZ'IU*>h7iH'###$##F{2GH'0SZ###k8+tl$N&~n>#GHH$<,V6)&6#d)~R%(###;##Tr@95#BGL+##zS.J6%')~ep+<rAo5#g%/_T%/&~IZ####k##/^8###Rn*[l#G-)'##?j4Z2:~#&5##iI-sn(|Z)Fl#"},
{505.94,117.59,3.99,4.84,",x'rU^######lyF8^6.,#$##EI&9<43?&E,$kB5<.*36%NH%`NAeR,###%##BU^_>$Tl$L?$x~0kA0yD:l#%*mMS.-Xc&`,#4=I######4##uS^%##z#'z5#0L9,##uwVMH%h}L$##0R)G]%<)A######w##lS^###mG$$f$HJ1'##C917tA$m($##>#$~O5"},
{683.45,137.57,3.39,0.98,"ZP#C&#7(<$##;y7/7#=R+t&$.f1;/095#E`*7-&lbITG#BB+Yl&d##3@R8,#J?RyP#t./i7#;95),<PEEx##3##-AR]?)###`5%###]6E:#$M~Q2H%Fo-{$)G%*)J''@RV#%&##4&/s*A###I#%.##KT0/A*)d(A5#SA)9DR4u#s>#V:Ps1<###9##^N<###"},
{669.99,139.46,3.98,2.19,"Mz,:@(eB3mP$6##J6%r|5?h=###'##V$$cHR######R.&,&1Mg#*y69[#He/oB/{}JFI&993'IR`u%&0&EMR_,%*##jH<*;<0##6#$m$#JPMGJ/z,'L##3IREMRS05z6'sS+].'RY<NLRDZ%######k,$t#'@#$###,Z#>?']H$=w,G%+.u#P,#gC40%,###"},
{691.08,144.64,3.79,5.54,"Uc#/10PfQ56'4Z$w80GS,`3Cz5#%K3z.#gEF99$fH)}Z#o#'a?#?fQWR,^P#Iy)ncQEl${I,P4>K~0+K$=U7.R&[P#l#:xZ(4l#kN>###$##:cKQ(;###r##jdQJ#%p,$iR$v./`G#LhQ26%yu'O?%gG$A,#ecQ]P####u##JgQn,&:5#5##o%)/U+O+J%##"},
{660.88,157.53,3.42,4.73,"###ND;######2['C2>###.##d[+}5&w#%527VZ'###eH#TL7B,#g{?######?eS(80###`##V.X}b#9w'k7?Tc&fc#hzOMr=ol&[d+###8##}/XbG$###E##t3X2X=9d)=-$Ld&+$8<%W7u#'$&v#'###1##&(0iZ(###/##`s,zU:.,#$##8%(s00X@(9m("},
{735.01,165.76,3.42,1.59,"J,$R^(bG$###YiB][&7#$F0&ed,%$#-QFjF:.,#h##s?NCu$Gf0J#:###&##goX:%*o#'7/#D^8lM%%nXUm#z5%3p%|MB&##B95B~)0,#-(1GrXNr95l$M$%RI(lc7v4M&##PJ.U%&tm,2,#eY#%##6e)#cFvl%uw+W/1)v&)##_&);81###*e,H5#ou&4,#"},
{735.01,165.76,3.42,5.23,"te/######;##eo0%Z$`,$T%(b#$xL+)y1y,&Xe0Jw#^l&b-#Uq:+u#n#%oG#J~Z###e%&Xu6{,'=##BrPLL5696-##r@(.M(,q;###z6%Mm'9^Z0,#m$)3R%418X>#4bZJ6'vA4$##oy1`&(bm*3,#^H&Z5$t`Z9m(A,$7,#$U/R{1vtL0,#%f2{G#w@.8Q#"},
{450.80,210.64,3.48,2.97,"|Z)%##:5#K,#D[V######I&#[~V######7%#;U7###Dq75,#1y6.,####k##4~V######8D$%]V{k####y;(~q;:#$e]3o##;D<L5$###9##v`V}>&###w,#Ut1raJ###9v#eg&D6RE6(5##;o295####5%#y{>+u####b$#d0)1T4###+##AA#l*H######"},
{460.47,219.32,4.36,3.21,"A,$###'##?v(stH###$##m6&Y*Z###.##.U-J^5###mM*B^7Nc&######P>#UbK###:5#wP#:)Z######r$#B9R###)f+m5#%?&/,#2,#:u#8tGQG#95#<##++ZyH*###;/#Bc@]=KY5$f$#hY#GZ$ne.?R)E02ZH'cG$/##+UJOR,###F##z/'C6T###%##"},
{737.08,225.63,4.10,2.56,"(##Of'BQ&<m(p/*`R+%l#(C.}o*'$(.##Ij<~A*ZP#C##kS1~#&@##dl#@X::sB###W$#lX<#AT.,#z##sDTl7-^P#KE+H@TYM=###Q##sA++g4.,#f$#f,E5DT|Z)/##&(*8*1f1;(Z<Bv(jp0.,####(##*h+(c$###|G$D;*`$+###J5#|z)4o27l#F>#"},
{574.26,231.83,3.62,5.95,"bG$###Zp,zQ)%^c###p>#DJ*}]c######k$#wU<######O6#CH'###+R##M:YmV###S$#aHIp^c###'##,@$}U<M,$$##g6%U;@2,#&##66%|aIfP#@7%X]0;_c^>$oP#_c${L:HW16#$=5#lT46l$###/##hM9G94~>$&##+`cNH&k>%%##S7,Y;(5WB:,#"},
{709.58,237.89,3.77,5.15,"######],$ZP#5_=###|Y#*d#&nX5##/r?E]#(c$9##MoXKm)###&##Al$|k#D6P8,#Cc%P,#&pXPu#ZV<P[$1u#O5#~xOQ#%.,#%##FH&|k#wuO7,#kQ(J##KpX/l#W&,G')o#&KR$jCQ'-&.,#%##4H&/,##9-Jc$@n.###X:LDg4mH)4,#][$Yc69bK.,#"},
{614.94,247.22,4.13,0.41,"X,#;tA(c$###WI$m5O%##|G$k/#D5OsY##Z$BB#)]3######V$*Hy)+OE?,#j6OeB4C['q32+J&~{:MAMY93Pl:5x0>m)(##+7(L,$/*6&-'fRF`U6b@&Bx/Nm&_k6GFE091#HMu>#wc(LL(v''L5$v5#;Q&:v%^>$a6#ph?0,#f5$eu#z6O*6'######ih-"},
{614.94,247.22,4.13,4.00,"###2K(&?&###^P#hj?Mu#i5%m##wV?J-$+R*5I%z,''$#7B5Jv)o_)fq>-##:uLjp3Nv%wj6x$'*K10~CPL7KlAq#'`6(l,$X[+^,#SJF:7,=COMT2)I&h17#.'7W3CcJ4W;hf5'##_-),P5.,####q0$-*EIc$dG$:$#I?OL>#)Z$QG#v@O#########jtC"},
{220.44,311.70,3.79,1.73,"###s5$&Q%###X>#`_1uG%###N'2K^5U#$=v&>l$Eu$>d&5?&.,#u##8T3~P#J02M$6}A_pP#.G_ZF;E8/J%'he/A$%r<@yc&E>#cZ#u2A0,#IA1_w#fG_6,B,D_/w)NI)4[6(n)xe,Q6)2,####x>#sH(95#.,#&##Jx'Vq;0-'ZP#=##Dz4uG$6l$.,#=l#"},
{381.71,338.74,3.85,1.84,"/,#qP#>l$E>#Ou#=D/Fl%###eq7v80Su#_v'cG$A5#)7&VZ&###{5#X]5###~^4X,2=o]#Z#Wt]O{4sS0hm%Dw,)0(0OBJl$###>##Hh:D>#{.0f?#Hu]ibBbp]*7(9n*d$8'I({J)TN>G>####bP#IZ$0Z%###UG#wn%}082Z%]P#x,%=^1:5#Zu#)~+9Z%"},
{533.04,352.98,3.95,1.18,"GD%mf6###$##,c5;c%###3,#:7(j,&###&##jG#iH)######ga:_5%###Cc#F<V:Q%###k##R:5>-G95#-##8^1ty82,#nH$D{A######ux#y6V|u$Gm'b;'O&2:N2TA-7F5?4G&%*Ul%9)&W_?######[(#k84###6Z#+t/^l&###y#$~>C2.+?-'eY#L[$"},
{533.03,363.92,4.06,1.76,"&v$#d(###'##]>#Vf-######`%,<7+7,#$$%z>%I,$3l#1Z$v{&J*D+u####7C0iZ6w]7D5#SV[hy1/@)sv&0m(*Z#e93hu%EY[~v)zY$2,#]-Kz%%UtDh;6rT[Pm'WH'-2+sQ'/[&6d):u#vT[######:,#^S[###S5#[z0FR+95#%##Qg)rP$95####u,$"},
{379.69,399.67,3.95,1.98,"I>#D,#x.0<5#CZED,#2@+uR%QZ@KC*/,M@S&`5%5g#@[TjP#Fc%^P#/@*_,$C}F{b#a?'`l8*&-SA&B^T`3;uG%GR#C_T#&,@m(hY#E>#+##4LNq^6=#$pQ$#%'./CUh:w5%$##Yu<713~?(FZ%<5#E>#B5#_l#$S+VQ&eY#-##q'+my/Jm*$##2w(BT.N@-"},
{379.69,399.67,3.95,2.89,"J]R|8457&F.*Dv%6p4O9Jn@.^@,&##|{+S[RdI,f?#_%.rR*nz0L2AR?$^,%:]R3(:c%(C(1oZ'xQ%]_RTbI.,#j$#e[RS.+&R)hY#rH#I$)y/Hl^:V,$pc%f[%t@HH[RW5$###N##Ls?vy4###xb#yb#95#e5#hc'G>#.,#D,#II*d,%|k#95#*##i>$?v'"},
{379.69,399.67,3.95,5.98,"xG$w6)95#*##Nc%8#$D,#NI*F>#.,#u>#gc'(l#95####(l#[FC]L6###T##1RRc>$m[%;nFV,$`Q%WoG{y9nH#|l'R-)^P#PRR#~+.,#J.#QVRIbH{c')R%5e'G1/:TRUg86d$b,%U(0@L:@7+IS,48/eH#ij0TlM[R+%###_O$.+5w'0(8<$%K['ETRoA3"},
{495.19,413.69,3.79,5.24,"wJ3m##?`]`Q$lg7FZ#:GJ%##HI+U>#bZ''#####,##&?&###Q96(*&.]]s>#x.WFB%D;@-##'95-##q[,G5#$##sP#5-(###3f1C|'d~]3##M]]C%&z,'-$#@K54##7J0.##3,#]G#z6+###mL<9W$@S0n$#n~]fP#.,#x%#u?*=Q#:$).#####H,#|,'###"},
{396.86,416.45,3.63,5.21,"W/-%9+2e.xY#-3Z??'XG#?m#[&/[35LT0?Q$?-'c/'-~-$##}_@46#1M>o##J/Z;u#LZ$zh&(]2:~$SUKv.)g[,DH#Yh6@,#1`>5##P::Fu#S4Z190.u#u5#}A-s?9fx3C5#%y4,f'@m);##E6(9#$oe)p^5;D9`/,yY$j.(e[*0=4.,#.##o2Bz5%###v##"},
{310.16,465.10,3.55,3.22,"(o)4e*W}H1##E-#}f48Q&###&-#oe0.,####9m#uG%###<5#z&0*V)5y[x;8S&GlL2:x22l#GA-.-'/Q#eI,v7)(c$(##pZ'cu%9q%W{[DK2wx[}l$<g3P52<96###p##}{9@u$###U##`'7I?%.z37%LP03Nc&###zc$OrPA,$###2$#p<A######M$#:81"},
{705.86,463.93,4.05,4.97,"#o$l{]ZP####'?3pDB######FS')7+X>$###&##9?%A%.###u2:a36cR-f,#H#^z,'###M##AbD?m)###&##%##ZR(uG%###,aBfu$6o1(J%Jx]######n(#g2AuP$D>#W##$##06$]#&###^m+%##Yw':UF=WC###:5#OW%FH'|b#ZP#V,####W>#u5&###"},
{353.73,472.20,3.67,5.24,"]I(t/+jA3-6%Q^RL$(^-):d&sc'V7#G#INu$o>%XZ#Z@*F#$Lg8>[#mmWd5#knW{Y##%)iK%Se/R?#:sWq6&]#&=##FK/}#%HU8E##anW+c#:rWq%*KZ&8H#ye,Js+h4K7,#[..(R$l6*^>#M?'46'+|6oJ/rz;Yx'f,%D/*vv*Nr(~#&/##?p7K#$/,#U##"},
{385.66,475.05,3.65,5.45,"gN6K41fiD.##^aC+:'DT5&##M_=.,#$?$>u#uY#_P#77*.,#lI,z_'`q>t#$C7TF~&^Q(?##OaB1U8U5$6,#Jl#BR*:S.$##;R*=a3LQ'Sd#rbL4R(^P#*y);f-L<BA5#DX;_w+lZ(KR'?M37Z%X>$&##/Y2]#&iP#;5#CSBt$'dv)J5#g;Tpc%OG#D2+%IL"},
{396.55,476.06,4.11,0.37,"%##N,$%Z#NI,)l#GQ$s$(j~1[u%YG#{,%6.QaP#OG#w,#:ROP6#IT36Z%p,&Co.+L,te1>l$d%TMm(gY#(=5tR)h{?4##~(TV%#23Cqe)###~*<Jy3/J/Y5#>&Te@(V.,_Q#a/*l'TL{4rv*V6&y,&h>:M>#}s:)6&NR*$##O)T$##Qc$~P#M:*7l$%:'g_?"},
{462.62,489.80,3.52,5.25,"*<<nW+t/5_l#x]6>2&@iB>##Z?TtY#ZP#0R#67,1,#TG#kl$Y,%qt3#J/g5$j[Ugy&2-(v##G/Rj7/###n##g7(&~-~P#>5#+-'Yi6###'I&5i;Of2@5#5R'Q<0?aG###*~'+z.@d*F5#;n':5####%##cN>4Z$&6&Q>#X/P2J,xG%###Y`U16'9,#u@(a%E"},
{696.72,502.55,3.75,1.71,"1h3|047#$6,#l^ZmP$$##9Q#]h=Jm%I-Kb$'###g##9]Z###h`DT>#t5%1d#K]Z@6&pb#M$#Lq8FR:KC;%#####/q'H.Y###mf6,##mY#&X+i~Zrl$>6($B#E$()M(:IPz,'###ZI&TJ*uR.{k####.##&5=^H(uY#c5%DI%H>#vG#M008Q&###/##uA(n?*"},
{691.36,511.48,3.83,5.49,"YA1}k#.c$B##.T/nl'$##7,#kl%Ll#i:8Iu$:5#.I#STZ.,#gd%p{@eY#$##G^+QdO###1##eUTtQ).g*ym(h,&###=|TAH&C6#`g7######k|?P]1###1$#GVZB,$'?#);'[S/S5$_Z6(x+I,#(Q%95####j0--H&###%##eDPr5&$##U5#fN,)93ax-|b#"},
{679.31,515.80,3.97,4.00,"L*/+u####4##bO/DT5###&##6v$~YL3,#G,$###R;9-?&k>%H~,95####)6#X'WW>$###3@#5^6cl&8##(27###9Q%~5#8YJs-+###=##];3A^U###i##*4<cQN%##K>#]M.sb#;5#>5#]U7-H&'##.d#K*WF6(9[#%C)5)W1d)iX6X[*(=8Uv#L3CiY#nY#"},
{728.80,520.70,3.59,5.01,"r~,r5&######MNZ.,####<##@TK{.,5^7fG#&##v,$nMZj>%lD6pm,###&##IKZ{k#####%#onSqn.X-*s##6,#nQ%f-QA,$$cJ(c$###$%#(MZfn0.,#T$#S}1v-V+u#)##h?$xw+UA2/,#8sC######`/#[B4nQ(CH&X/(ir9od,G,$/~'RJ,?Z%OG#3##"},
{743.05,521.90,3.82,3.99,">u$###Wc$CQ?uG%###X##3:U#w,######P{Sp-+.,####RCM8Q&&##d$$H)San,lZ%^.&FvH=:U/u#&##Wb3.;8JJ-###&h/f-*H>#vJ'1(3k#&hP#Lg&}+F/mN%-'(##z;.88-:W8#l#C5#$H%/,#I~+4u#%Q%###K5#'.)'d(C,$$##b?%>5#1H%P>#EH'"},
{783.97,522.29,3.71,3.77,"iK7`>#|U7%##xLe&##|n0'##V/3###DNe,Q$[c&J>#(e)-6%N^8zJ(fe0)l#4MeB,#Jw-V5#Mq=###~NepY#yY$###8.*--&L..Ou#su$[kCaLeE5#]Z'Vh)Pz;{G#@h_cP####$##tS,rZ(DI,###$##@q/H#O5l#eY#FT%o-*z6$PC;*#####,##He-ZP#"},
{354.06,529.33,4.15,0.40,"p5#75LN>#kc'Y,$U35]&1q6*z5&.Q#)S*F9X&c#Yc&pZ#u8X1##o#&SQ'7#$6'2-0+X'7<u#79XDn,0l#1@@iB.Q>L[##V;X7[$k6*9H&iY#TT4AQ%dh5q.+R9X{x.3~+bw'EM*Z:X@M2F/1uG$4l$E##47+Y%.+u#;,#'n'6j;QG#'H$]H%R]'yu'se%.`A"},
{595.54,535.32,3.88,0.68,"###/##k,&###%OG'##N5$E##QmV###oaB+v#3I+###AARW>#D>#C##Fl%$##YmV'##OG#9##ZmV###Rh1A&*SQ'###X)S;v'-u#P5#{k#%##n#O######4##]rV*-'{6'2w&O.)vI)3sVLH&j>$K,$######[A-#########(>2c]5D>####;R#8nC}h@###"},
{417.73,542.97,3.75,5.36,"pr:h`+;?S4##HWA>_$7>N.##/$S$#####e##~>K.,#:5#K##|Z(mE-1_=uG#dvSDB&ed,_##|$L=A1###?##Gj5{v,###%##{u'5b6C,$$.&C*BCx/%l#Jm$aD2/}H.,#pR%cy-C%.#6$-o)OG#E>#-##Sq6Lu$4H%[,$wL5+(50?&95#aX0Z%.k,$SS*@*0"},
{389.39,557.70,4.28,5.18,"Z&4###B,$*-#sm,]Q#CYL9m#<[T*?$w7.;q$X7.T##-]TM[%4aD$##@c%0##-;>5##u[T$6$J~T2J'O~/GS$IJ/`V&0[T7Z#R12z5&iY#$##c:.':98C0N~-sz;v#9$m&'K,Ze/9r'>%.H##o>#*d&wc(###O$#gIS..+Y>$0l#dTOUG#`c%'f1}d'.,#26#"},
{743.30,565.93,3.80,0.60,"W##pp6######B<@k%-###)##gnY95#Gn.V,#Su%3,#`?S2,#G#$v#'<5#.,#InY:#$###3##5nY###.%(4T(eu&<5#MhMn$)7#$###Q>#]P#=.T######?,#esYP5$w5$KI&+y/nQ&=sYiH'D>#&##1,#;5#le+######+##rm;ol'.,####3q%?*=kU;###"},
{595.93,568.18,4.06,4.39,"F5#ac'+##SJ0]>#5?'###)02}b#pb#$##AR(###%##1H%iG$&-#n6V###/,#i9+`EF###0##PU5g,&+##wM3.,####<.#LYJ%m#*A0###NJ04%ISH(###m@,F<VOG#+##I-:Py.Tl%'s*'9V.,####$##)vQ[R*###<##tlQjSC$m(<##7~--<)KmS;J+l5%"},
{737.52,583.42,4.04,6.21,"###wb#sb#/,#'$(###SG#gH#Zf5###'##P`%5?'y##P]2N=,.,#.,#+##Ac%=PM######Am#]x~Nc$n5%BE&<~.n'%/y~(g*0,#95#*##oP$O#M.,#$##gP#7#]bL3W>$>##ad'1#1Y=K&##$##cG$)##D>#/x.|k#%##.,#]G;Z~-eY####tG#rT*Fn.2u#"},
{471.12,599.07,4.21,1.95,"Ao%d.SD>####r`3SS.rl%^>$p.-@7%Nf-YZ%]V<_P#/,#k,#.k:)(87o&+.S^1SXu$R6&_w+504f?$b(5sc%MPJRG#_P#fl#6b7TZM}T%;.Sv-S95#k/%fPDe./###o/,_^5^J1###M>#v/,%1S9'22`>8d#)-O'##W]0x6'f?)95#UG#A-'Du$ZP####n,%"},
{559.50,600.43,3.61,3.67,"###*7'######lP#dp5######B~-_-*###gR#DZ&cu%###ra6###p&&?%.###XB1eP@R6)5,#/fWOd)46'o9#U7.M.?`g7tf(###_##xaB###E_:8m%=292H$DjW*-'ld*56#&:+|@)nSDR_=&#####R'(8Q&e,$B,$ZE*5[*7W'^m+07$=@,2r%(::d5#$S-"},
{750.02,609.03,3.69,5.76,"Z#%B,$L>#)##A/[95####7%#p4[T8T###J%#Gw)p4[_5%3##Pm*0u#[G#B,#Z/[[>$###:$#V}=7TLKl%L>#3c$Qv8dd*7I)M$*PG#|##iV=XmU-l#K>#@e'[7-5@%HH'i5$'$(4##t,%a(5PG####~/&S5O0v($##_I&:T/II,)##VH&g,$896######8m#"},
{825.45,613.96,3.88,3.82,"###u5%./+6#$|HT+##&143,#AnW###W~S~>#-u####;pWQ5$###U.+@$'@H'%cMPc#.C6a?(goW&##p?OAl#f5%$##TrWF-)###2,#U?#ZiCRn/%##oQ$o8TEnW'##Lf0sV+3H&%##[fER-)`>#;c%0##x-+K80OG####i~+gD=+?%Fc%2@$3u#YH$Uq8>5#"},
{767.65,621.87,3.71,6.05,"oVVGn)o#'%##s.-}x'd(<|%,BC:(##EQ%}6<;n.###%##JU%oSV%c#`l$}`9UE?vP#ue+hZ%rSV###5,#UR#u~2I,#Q#%II%VZ'0,#=#5=IONJ1###*p%|&.BRV###X,#QJ%Nf4e>#zl&?6%CZ%47)@TV&c#'/1###A9+;w([RV###0?#l@%&x0.,#X5#/C("},
{748.67,653.31,3.62,2.47,"+D9P(/qb#@6#xW.|d*###&##~7B#########ay2#########6L6s&(r#'&[#(+6vl&17')e&GzP[>$eP#Hc#eP?######$##xuP}#%yb#(B%|B8=5#pU&%k;vuP$##ec#T:-RZND>####<##p&IiU9pb#Z5#@F1cuPr6&ZQ&1#Fad+(H$3-%Q{A######~5#"},
{460.66,725.40,4.15,1.82,"$##J>#-M4wQ)%@(c#&Ad%{}Blc&7Z%De$o=G<,#Vl%?$#~3FN,$2n'XtE&I)c{.J=DH^7Bw(}xRF7-%m$SR@%~,e5%Vp#XxRcP#2K);`?zl'nA0H$'VT2`I@&{R=C8[>$`V/r{-7vR;v$3w*###^G#4v&jc'x-)oc'~G#hv'@?$I^49#${b#y##LFIOG####"},
{823.79,811.71,3.54,1.88,"%l#5,#lc'/,#KV>0,#95#z##<)A######;u3,u#.,#y##.CY(l#R5#6e.7,#OAY>5#{k#*R#rBYbG$###Ur&Z:6Bo2pc#IxLE-)?##07,>?#S3EF>#F>#_<*NFY}]72,#/e&(h).AY<5#N>#g7/+##<c%b?#MFF|Y$$##Vd#9-%>.,SG#%c#9##YU:PG####"},
{596.40,821.62,3.36,6.17,"]d'#########o7.######I?#/<D######%)#]x4.6#mP$~V#&%)###*##PG##8W######j##gwYx#&###%)#fK6Q+1xY$1%#p5%###b##,u#/yY.,#&##?5#T|YJr;###:##C6&ZTE###$##95####36#{k#.&+j>%|P#PG#UC)NQP/,####?,#plA.,#;5#"},
{670.80,823.43,3.74,2.06,"R7$5oK}>&$##(>0*94K#%###~I(?$'W?(o?%hZ&q,&VG#[?%R00e20E-T$##4bZ#d([P#H5#]i;M-)0c#).**%'wQ)2,#PZ%;A0.##W_Zr$'^~Z###-l#$2%?sF###m,$A],=#$/,#Bc$4-'I#%$##,G2_ZHqDB95#<,#^,4'@+F#$N-(5d#%##{G$Q?(###"},
{683.20,827.34,3.89,1.51,"@Z$LV:###&Z#Jg,C]4###$m#X@+rY#il&eu#_$**l#+$'g,$=A)it[`$+3##it[`aG###b##Ls@aS*6Z%*H%18,@I)D,$mQ&=(;;U%>q<b9%Vo[0c$N>#*i#n3@il%n,%G~*8@)4c$S,$&K/M$*0##@e*;rO`/4$##|,$j|*77)~P#e?$e[*)w)fG$:5#vY#"},
{538.40,829.62,3.90,1.66,"I#%$##.,#-]&=$($w)|k#mc#S['IJ/F>#C5#F>#=w#FcQ;5####<##Z$*(8*zv)R}/D*F6u#|'YJz4oc'Ql#aZ'ue#M%Y1l####Y##G`@$-'6o2-]#2&Y,n'Z&Y*7'X$*Y%%86&G)&1%Y######k,#/r?###(c$2d%Z*=.K0?n-PU/Ql%#J(I>#Gi&1%Y###"},
{538.40,829.62,3.90,5.20,"zDB(##0c$ZL$:~).g0@.-Lv#l$Ej-*OG#+##,%)yQ(######HoY%##_v'@I&JT1^d&1rYXQ&lqYEZ%v$*:Z#/e+<@'xY$###XrYE>#Hu#]u$)//L,$(VL2T/=-S=5#08)GE1h?)Bd#w6+;,#ksYD>####$##392U5$n?'}6)pb#&##Cm'cJ..,#c##&'6<5#"},
{807.90,848.16,3.68,0.53,"b]2E>#######`|Xc,%QG#*##*^/1g(>wX0,#gY#t##FwX###S]495####2##PxX'l#,u#G$#Y^7~R%ovU=##%##;##{xX###eB6######&##%yXQG#W>$:##<L7|c$p}J)#####/6#dxX###Rm)######<,#zwFDd*######e~((Q?_Q(######M:&r':###"},
{870.41,848.14,3.69,1.79,"r{YFR+######6v$~D:.,#'##tPMyu%.,#M$#UM@######O@#=zY9#$F>#fP#Zp2oR+w?'cd*0yYZ>$6,#zd%D#P######$$#5yY###%##1.$wU<:u#_?$'pR7IUI,$`>#$Q>[sG######B$#cwYOG#$##?0#^p8*7(2w*6](IEB7$'%Q$Bn'hK3-H&###<##"},
{793.66,867.54,3.55,5.85,"K@,###&##mY#nP$###5##i?'Nc&######W`&D>#6,####Xe>*o/.,####%##Iq;pb#$##cZ#Ug[Fl%###J)$['3LM4.,#0|'0~,######$##`OFOG####/##]k[w_@###?##TV',h[.,#)##fZ&D>#3,####ow*%I*/,#$##8n'DlISG#95#d##?g[J>#A,$"},
{765.65,909.33,3.71,4.94,"[qWA[*G>#%w#*oW9H&G>#aQ#B{?WG#OG#9#####*6#wu'###-pW7x.UG#FU%GpW.7)I>#06$.lLXZ$A,$$#####Lw#P6)###7pWkR,###'<-DnW&$'###AO-oNEKQ#SH(D#####7$#Q%/%##3qW<R+###Mw%[]RyH*###|I&2o0A-$Yl&1##7#$r,#wu'.##"},
{799.05,931.42,3.69,5.04,"hq:#########)3Vhv+###W##,&Ns$,###F$#'&.EZ$L5$'##q]5######(##t/VTu%###?&#C7N=06###y&#k8,Yl%.,#+##pw.######$##:3V:r@###<##r9JM.V###,##(BQuG%###0##$c##########8.&z%1.,####60V$80######L.V######B##"},
{789.88,936.54,3.85,5.00,"N5$#########Z.U######V##@0U]m+###aT(;OAk,&.,#Uv$qb##########r/U######4##x2UL@-###<&#(ARHZ&###x$#D>##########8.U######V##l0UuL<###0'#@48kU;###:$#.,##########1L2bG$######Yl:T$U######41UwT7###%##"},
{1006.24,954.43,3.97,3.27,"eu&###WQ#{kE=c%###(-$oh76l$[P#gc#B3Bc>#|,'L,#TFIfE?###0f$~d)9&2###<('i6(e,MW,%s5#Q/MyG$E$'{]')6O1:O###|##:Z%GA/###{0#tf5CZA$;<in%)p2X[%L9OlT./[)1:O###$##/,#tO?###-$#ER+=5#6c$41&IJ1Xl$pm)sD+1v("},
{1006.24,954.43,3.97,5.43,"QN<`#&>l#^v'O^7Su#W57tD;<-S.##Wo*C9*t83L,#|G%7##j@A#?&:5#5,#&l9%yI`2;?c$4.S&f(3f0{Y#U95gK${v,$##E)3]6)<5#F5#lT1BW7###]##O-S=w&(c$:$#lm+6_$h1>&##@9&^x4###%##Ps:aZ'###n##?V<SD)?d*n##PG#hL#d09###"},
{93.72,34.94,4.95,4.84,"/7*3,#{Q'_>$}NF######4##Wx_###j)C&##xY$###ex_###XS1O6#m@.#?#^OI###$##]##hx_###O=J*##0Z%)##]x_###&95Vc#W>$l$#fh?/,#F>#)-#Ux_###8FHJ##Qu%###yx_###>.-.,####+$#`'9######N##]x_###/_<,##Qu%###<y_###"},
{97.81,50.62,4.58,5.17,"[-(^P#wx,~n,4j>1,#=v'###5:Q#########wxO###owK###Px2ed#N6Qql$K6Q#6#zf69##$eP###M>#L5#k:Q###Ig2$##(x.5{(sp:~G##6Q^Q#|Z)r##-|BhY#;5#3##]:Q###Dl$%##:.-2m#6#BI[&-GKC5#H>#,@#4^5.,####)##ak@######$##"},
{78.28,65.38,4.44,2.04,"OG#T##V?SL>#1V=L##lD@q_'p%.<[&,?H[p+###RQ#^?S$##k,&J,#l?SL5#r?S7l#^1;YK$]f4R%$R@SR$%###r##s?S,Q%oP$%###eO2,#n.I-8.+p0[>#Ge(|E3C$Seu%eY#Q##;r8G|B######5ZC###c5#t5&X9(TH(tc&@R)ru$|w,}v*4:/y7,^B0"},
{78.28,65.38,4.44,5.44,"Am)2d%1w+f/)pl#Qf2(Q$'@*[^(c6*G,#%Q%L56#########N01=}>8d)Z,$2dRn5%f.'eN:e0.G,$+nCZn.YZBF>#X,%$##9?Ks7..,#%##xxS'Q$/S-p,$yB3fv$ewS+l#_xSY>#8w-$##>s?gP#cG$###u^R{H'Vl%'?$&o,wz(jvS|P$VwS$H#mZ(+##"},
{738.68,101.68,4.71,4.84,"o99'##ZP#y$#AkK######k##ZwY&##`x4'##?u$.##IwY###1&2######:%#0OD######C##-yY8v#)]33##4c$vJ#9wY###_?(W>$###U##-*BA,$###3$#ZkJ'g&/1;g$#G>#FD#9wY###R>#f%-###$##e]-z3E###0##8n):3TyH*2#####$E&5U9###"},
{254.57,108.82,4.38,1.84,"###k##<$)###F,#.X5}>&###uf.W_9Au#e>${5&P#$$$%X5$###J%#696###IT/8[7<o]}Y#4u],+<E7+J$%M/1ku$Vx.Kl$###z##}jG###k@.u[#Hu]$>ASp]eQ&`[(Nm7y@,J~+%Q%(#####_5#6w+###.,#9,#X0+jK59d)%H%A,#0U/'Z$]v&ZP#&##"},
{129.42,140.25,4.73,3.04,"RG#5b6;c%$##/`8q+H###[$#xs>:$)9##x23Tl#JZ&*S&W]4/%-{~$;{@B##9-P~Z%zP$s^$ec@riCI/+Ih-|D+P=Ig]'Vx2fI-+##M}CO>#B1PFl%}G$Hc#Cr*C^8XK,FQ&2`.le08%'*8/#l####'=:.,#.M-W>$'d&=#$9z$v]7%##sb#*%#0E@3H&{k#"},
{662.59,148.30,4.19,4.70,"Ec$'N?###$##N.,z$+MH&[01P6)###q?$%3<X-*###(##>[$|M<}09###2##JxWtl&8[&<$9@w,XS&23U/,E'R'i#%oH(Hu#xjC0Z%###2##~|Wc|@Ac%/Q#9])S0H`4H7c$nT#(g3VQ&95#p@&Q..###&##w<,3/1###.##f`7z[*}n(Q-(NJ'&-'&M1Y>$"},
{478.19,161.49,4.23,3.52,"{$)kaAG:6lP#CK3p~/7L/K%*:T~eY#$##Hm#$_<######>##T$#[T~95####sj6ekLTG#?5#QV~0Z%###G##V*F######2##&7#n1>.,####MV~~./###s##{S~.,####/&#Pq;:5####=##AQ$K-(pb#/,#yV~mP$###=##RX~######[##C91###.,#$##"},
{513.90,168.59,4.56,1.23,"95#kN5Fl%e5$>n-+^6%##In'Jy6fQ(R6'3F48H&:w&3K237)E7,^5#4l$)6$<r=+c#]u%-e(#0V-c#ST22n&^>$$v#~IJQ?(;$(1,#x5#.p/V%.XG#-u#'#6u-V###`6%CH8CZ&###84/j~/######C2&C+B3w-###&##lZ9Ly7.,#'%$1MOF>####FZ7fU8"},
{731.66,176.23,4.31,5.28,"g,&###i&(e#%*`>###-##CZ#927cm(5d&)A)rZ'aq)oo0VQ%D>####hF1C,$}}J95#)S&{u$4SX'##sz+2Z8)@+5##[SA0z2######m=>###Z;@###O9*&@(UWX*v&0o.#v$Tf0P%(RVX5$&$##I5#J)=###^Q&u#$GT2G>#c?9@'5nP$&##eT.w{:(@+-##"},
{461.33,177.68,5.09,3.59,"Ex),-&Mc%kG$9gGv{C###%##UN5AfY###%##(6%Y93######OI*|Q&[7G|d*zfYl5%Nc$I@%ieY######X##zv*:5#######g?&d%PUV5EH&IhYV..r,$9R&+fY######4##mm*######Y>#O@%WbK######=hYm[-###_##teY.,####C##b-(][)9l$Z>#"},
{231.33,181.72,4.73,1.66,"co4######&##q9d1##*6':##GEDp##bZSg5#Z{B(##u#'i,#,_;ZP####%##n9d1##?$)U##HsD-d$_4JPc#f|ET#$B6(n5#Y;@95####)##j9d2##B7-|##'aBTR%4eY:##@W;,p17#$(##;h;OG####%##5;dUv&^Q(;##SS/I<1P_>,##q#&5.?E>#&##"},
{724.40,233.75,4.70,4.61,".,#1##~,%=5#vA4@,#Y>$2K#@L:q-#qe1BM%.,#z%#:GI(6$Z?(fc#du&.##r[Tju$'H%))&w/4m/&l[TtU)###G7#1cO###Pm'e@%Ce/###<`TqbD1Z%~,#>0,<`T_83J5#|,$2(/rK8###8,#1$$UM@###7]$x4EDS0###'N5K4C######5A+~x,Fl%###"},
{277.60,250.92,4.51,4.98,"fB6######(##-?J2Z%###z?%yJYV,%%##xQ&s=J###_[+%##o1>######@##pIY###2c#^h,3JY7,#|I(sz5?JY###Lm)F5#U1<######6##SKY(##??'{P#z7TFH#fOGa5$MKY5,#bZ'%##SS0######&##OKY###9l$$##O6R&##C*?9l#*JY###eu%[>#"},
{440.46,249.44,4.00,4.98,"2,#T=1M,$L#%J#%7s,.,#P,#ac'R6%###v1%95#1,####]s.2,#8`2######z^;6-%###o##[$TFu$###|U#BI+k:3###8z%6,#XH'######Mh5O.-###'##O)TauK###C##W?&O)T###)#####<p/######T,#XdG######lR#O%T######C,#|.M######"},
{581.72,259.41,4.85,5.87,"VH()##|k#}b#Jr?j>%.##U.*&g_gu%A5#}$'((8HP<>u$<,#NZ&K>#d5%D5#nr:nx2T,%+##Ah_/7)V[+Z,#{&4UM+z}IJx'G>#]/-+u#&##GuMXU.CZ&<##]f_xb#D.)8x$I:6~6)C%AmR&:5#jq3>$%V,%af_O/.&##SQ#@f_.,#f>#$B%<$)D,$%i((B2"},
{314.94,385.15,4.87,2.04,"C,$95#gG$UG#2-(C##'>KRv$XuQP,#4_:8<(a~.BT($8YeA&?u$'##;Z%hP#OM>'##v)A8$$/:YRn%:XF$y#Wy7Q:$78YX$%JH'3,#2u#A5#P06:5#>o-W-%G9Y:B1#z6;%$E[*Nq64pPY^1:l#j,%.,#/,#q-'n#&],%4,#vw*JuD_$)J,$*m%p$Nel%X6'"},
{285.84,391.24,4.51,0.25,"iG$)l#E>#7,#$R'K7+###hP#cw,=?'~5#lT/<5#OJ.Un%}o5i5##@(cc'$##nU4,i1O^8bG#ehZy%/^l#YW5]&.6<B61$LfZ{I&m?)%m'k>$h07|>$3E8G0.efZWq4lB2{T(M2)1iZfN7,7+2Q%###M5#,n*h82###1##0n's$L+##iw(Um&PS'@$(NO,f>P"},
{351.17,400.78,4.32,5.20,"|7WHi'(K4`$$[06wD'~kM3##^aGQ#$bG$7##SG#TG#O5$95#'.+>Z6Rd+C$&f7WV)'pm,($#cYIcc'.,#e##uP#5l$sb#D>#E$)AIC$##]d(Z-HnR+###O[&%Y8>%.(##SB(V6'A,$EZ#n%,PG#mY#`G#lWAy8/&$'E>#JF6t`@+Z$u5$1#4$Z$=6$0F9TL4"},
{392.90,399.70,4.83,1.88,"I5#iY#95####>#$qH#]]5Pl$g,&k##sg6ToHaP#)%(}HL~2:=5#oY#j>%0,#a/2@H#wo5;d%_oTx.'IU7*P2]H'9D%UnTo[)?5#rG$zY$###Nh;H,$4R);D,7W<U9*nHM|>;BQ&~0%qoTE.)EZ%^#%D>#%##O;7hd+gY#Q,#nS+^&F|82u>$:5#^s60C4Vl%"},
{392.90,399.70,4.83,2.67,"X>$;5##*6MH&C[Tqb#08'rj19o08Q%c9J;D8]W8:5#Ps5~#%gY#E,$<&,xb#F&IN..b%)b6%0$@,_9hyKeS+Fm)L##>xDZ`@<5#iY#,Q$$##dJ/#Z$`@&yu%9_Txx1j7+*o(FZ$t^+H]Tlv*.,#2,#,u#$##rY#.l#wP$:5#W7%ne/}k#/,#6##CS*zw.B$)"},
{321.04,402.30,4.65,2.47,"KQ'&##~g1{$*AcQM5#g$)hh%*B+vp2(g4zf-[:+FH'C8*.@+I..###/*;fl$_-U'##0^0j9&801Oc$U2U|u$qpNC#$=]2###oH)###.T,,c#Z2UAl$z@,Xc#@U/il%;2U'6&g,MiG#WL2xn*%Z$###L?%F>#|)9OG#6v$~#%U{7K~)X/*Ac%Dc%l^-0%'.%*"},
{321.04,402.30,4.65,5.38,"Z-&mw)j$+7].Oj9qH(E&0Lf-nl%m#$H-MR#%Kv&yb#Dc%###ux3jd%OwWi5#fxW:?%.91lm$1x08R#DzWQc$Ff0lG#C[*###GK5###KzWBH%dxWjI%Q]2H[%P93kD&vvW)##li@*6$*6'###E%'t'9k9-<I*fn/Nb4v@,(K-]81?{%d09$##<i<'v'OG#%##"},
{363.55,405.97,4.50,5.99,"_U'1?&_n.###d$%<.+:5####>-'oc'*##D$'PG####s5#t$+fU(c6(97T###~20{PJ###)##>8Tn[,Q5#*;0v[,Tc$vG=DE>yv'f8%I6T;##l<@}7.ZP#.%#g;T*+GO5$hQ#s44I);%j>S%-xu'HI#]Z'UB+f$+$c#sQ'@'&*p'cM;tr:Fu#FeA&&17Q$MQ&"},
{473.24,411.44,4.75,2.00,">6($##[P#O5#p#'Q##~<EJH$$3D>##1^6$s(z-+c0)=eY`L&zu'%##.u#?,#y995##.:8'v#HfY4I$v^;]C$)~-tz#HeY@6$P?',u#=5#9,#KK4ZP#z$(|#%xgYP7+2y4Jc#;[(]B/?iYOx.2-$N-(.,####HI&8?'0,#&##j8+q|@OH'2,#M5#SgY[d(Vu%"},
{442.27,418.21,4.40,0.23,"YG#G>#.,#p>$mQ'?H&.,#>c#td+6#$(##fn)[P#EQ%mQ$WT4fS#5{=,u#UG#<'/[X8Nc&5,#2gW]I,2,#q^-@J,:q;[?#MfW)^#BsDQJ.###Y%J$q4#?&_,#>fW=^0sR,M6#HL)KhW,U2dc&`v'7Q$TM6]P#OgW=#$c>$3,#%hW=,#K.)T5$s@(;$'KF.)q;"},
{478.80,429.39,4.26,2.49,"&@+-##;|Bz>#f$Wg>#,@*7g#w7,<|4)i>:d%=M*m~0x6'xY$$x1###/N9kZ$r$W&##u/.BB%S@,Dl#K*Wn,$z?DB#$(y3###C.-###j/-s>$y(W(l#Zd)fu#FK-2Z$~)WR5$V2?*##s<;C8,?l$###hZ%rb#Lj6E>#Tu#&Z$0D,8I)%g,xY$|Y$kA(^x(8x/"},
{478.80,429.39,4.26,5.36,"dw(=y.2I+{$'MRFSQ'`%*).)Wl%uY#wdKH,$jZ%5u#yl'###a:;Im$5-R],#?xXBu#5n+ov#Te/a[#JyXPu#'o/Ml#s%0$##N'8###RdI0R'SxXG7%2n-z$%ke.%N&'wX'##cr?sZ$^Q(%##w6&XR,E]'c(9=81x_)Sf43q/(A/OD$8rA<##S<EmP#ZP#P##"},
{342.06,434.29,4.90,3.60,"Ym#iTTpQ%*/0k7#0nROx/s-*Id)N//ep1hA*r5&###*##&s7H,#QJ*V%,A,$V/-0[C0@+J>#*STC17'Q%l'$sd,MmI-m'q{5-##ud&$_<###,q9B7$9}GB5#0WTy#%W@+$v#x/*O(4k_29>L}@(n./3S-L@)WZ&$##t17Td)C+1j>%l>#By58r'mf6D##?y5"},
{347.66,457.77,4.29,1.88,"PQ$UQ'.,#'##1v'5$%x?*eZ$UA24##A8.JQ9BZ%<7&:nS>'0^.##A/######&s=cm(eY#5?#|oS$o&I1;xA%*n,=_$FmS$d$w:.kv+8##Fn-9$NkG$E>#$0$~lL'9/T`@c8#=m)ae'doSb/,(I)6#$W/#TmSdi@@$(yH$oS-~R)VoS}-*Au#/6%JXChv(#$&"},
{321.54,461.58,4.44,0.13,"3%#9M>###95#?H#mgT9l$*c$l#&qM*z(6&_:tb#@##GT)s4FR##Ag4uG%###qS--?;re1.##UfTFS.L$$mz4A@*o-+VM&m-R*##O.,@$({k#_/3e?%E3=g@)(fTd}D*e*Y.&?W)~eT&8)I&1###:5#*H#dc'Fl%###DZ#}%/iy7s#%)f)]w+]w%`f3J`)a{B"},
{498.94,462.57,4.49,3.71,"wH#;uEI#$mc'u%$MAUz$*{b#A@+g]/5J,5S&r5&M,$2,#I<75##jR(<S-6#$DK1Lk>Dd*,##$~U9J-P?(.D$8m),P5j9691.###P##nh:###+V:F-$EE<R>#Y`U(d(-e+},#iB)D8,oQ?`q=######*~%g,&Z,$95#%{&<R+rM'bR-$~#<~/+2%}_A[5#j7."},
{425.56,470.15,4.35,1.79,"iY#+z.xY$'###l#]S)iH)N5#xd-o##;7,h`*uG%yu#{_=SH9_6#UXF######7h7Nx0{k#N##UnU3[%l,&u1$*8/CD%`mU7e%/^,B81'Q#bZ',oUW,%###-6#noUCR&R@-]##~d)Z(&LnU.7(GI*vG%O-#~T4<tHO5$D5#~v&g*983C5H&*##|l##QNHS,1m'"},
{441.98,480.37,4.77,2.54,"###rQ#B/,zK8r@/9,#^B-3o,nYO66%Ym(0E)#'',p2j7+a]0###*##3],?w-`DA&##|F?Vv%IIUD5#WV/AV+ff,)?%,{Q9?&mY#Y,#i;=###5B3Uc#|p6gP#ANU2R)jn-],$=K)am)ZMU&6&Wl%rd#H1<###Y>#Z^/LQ'###*|*+e,A,$###V1&y]35.(5l$"},
{503.30,485.53,4.63,2.26,"&##Mc%vQ)###-~.6,#Q04vS'zo5F5#'dIiN.3T23H$*SPG,$5##GT2mP$###}{:d$(v#'h,#*0Njl%{$P5I$AA1O>#rUPsv(oZ#or>:5#s>%2ZD}-+95#7$$+SPA/,kR(Z/%Cu$Gi-$l7X.,Tc%-v&J5#RTPt04cG$}-#/RP9L0FU1Nd&9m(%##5`-ye,yH*"},
{681.63,485.03,4.21,2.86,"g1#;@VcZ'###=:$H{?/)1pv+~>#P5$9%%z(;3,#$##M],*?&Q*,M~0u5#86'REVI[+1,#S5#d:88#$s,$ml%######Ef'z$,d$*###C'#e~Uk@V###,##xt3vU<###%##?C)fY####@##mM8######i$#6AV}>&###%##>@=o>%M5$###6o%$?%zc(###Z-%"},
{431.42,487.32,4.86,2.43,"I#%5,#d].|[,e;BU>#y?(u_'4T*}x0De,c;6`f)wG%EB*qn0hv+,##L)<9H$#dS)##:'.FV(j02[Z%+iSPv&,KJU,$|f3###yH)UH$Xn.C5#+iSXQ&H@+4Z#bU.Vv'JgSNc%u}GPQ#/z5`d'I,#T01A,$###v*08m(ZP####hz(Ep3BA+X,%W,%*N0Nm&=m'"},
{431.42,487.32,4.86,5.28,"LQ%7@'d$+W1/});L$(`~){g6V,%%##ueIgl&###G>#$o)#J.4f2@6$A@U4m#z@Udl$^93N8%Tw.x?#@DU4m%Qc&5##V.,k.+U'6G>#PBUKZ%RBUZ8(Tf1ov&P]1a2(h?U.##&L8r>#du&AH#m%&]ED?/+86'Nm)<7?k-(R7*DK3cM'Jm**##D:6[Z'###8##"},
{474.10,491.71,4.56,0.04,"yw&~H(E>#$##qZ#)|>pb#/,#aZ&k4?)d&*S,95#PR$L'+t^:?9#&(9oA.###Ux+'%DOG#(##.&S^w--Z#sf-+~-`$*#Y2f1:ul#M.,0'S$##63@D00t5&|##o%SWGGx?)D6#4_(C%S|[*y#'y5&$B#H$ST#$7sEV>#0Z%;@#E3D_-'ZC-bm'VR#q07{;*=R+"},
{645.65,500.23,4.51,3.94,"%##Kl$.00&Q%b@U:,#eNA$##B:]###J4E%##pb####,dK7u#-##MOA0m(OG#|lOZx(>L96Z$b:]J>#|q9qG#7-(###zUSW,%hY#|v*o>#[V>k09qu$Gl$/wE>9]r>$K.-5p%d$*q#$PfX.,#n>%hY#3,#~p5,/1######f:0K,KEH%3l$gH#m>$dm$)::###"},
{375.65,502.38,4.88,3.65,"G6#Y/QqH%_[,5##50Q3J.-$(Z$)#K.VJ/Gv%u5&~P#&##hT2>,#<x/&K,`$+?d(JG:Z82:#$w-Q7V6w5&uv#W@,c/Q,m'G'2'##|>%@#8z,'<%-{c#&yNR,$i0Q'm%<S+DH%8K+Rp1.'.O.Q7/)YA2sD2cm(:l$###PgO;6'4K)95#k6$U1;K|.fH)2##TB2"},
{451.93,518.09,4.61,3.65,"*7$=IJoc%k-*%8$<vL`7)]m*PR+''/4'.Qx(du&QG#)##@,=>##W/,[7.95#j//E@H`$+)##BwT`'4U,%gy#)A/m~LP6(Yr4$##`u#j*?6#$'_:b6$^3ChP#`{THH%%@)Lc#4p+Ay.9C14QLy[''S.W_6NZ$ll&###=E<'v&,4/uG%@,#n/2m{)O]5.##,f/"},
{714.88,523.52,4.73,4.87,"/$'+c$######tB4W>$###i##ysFI#%4Z$H*)jH)z##o/[48(Bv$Kw.######?TY%Q%###F##}0[}G$dm*6g&[-(@?$G1[Cv&(J%Zf5###(##z0[Yl&###I$#&1[P$(Fl%1%#$d&pd&d.T`5%P.-nP$&##Dd#IeWA,$$##SK#U@Cue1E>#:Q#{H%C7)O[+TG#"},
{381.49,525.11,4.49,1.90,"w5&Zg&6.-)##_5%]w#:96)Q#l':P##(J.Z<+bl&V-%IGJG42Y##^+@95####N_<w%)tc(e,#;xVVw&Xd+C1$KA05{%tvV}Z$0K(3z9O5#Bm):IMy>%+u#RZ#9yV)x+fJ2~##x-*B^)dcMmv's$)q>%w6#@YJ7;<D,$QQ#7x.y00@+Hn#&mY#QH$i|EA-'Ou$"},
{793.06,527.78,4.68,3.83,"###IQ&,T'(c$jRVV>#n(26u#?(`###A4Avb#V,%###lk@Y~-&##x<Cel%A,$rZPm0,af3x#&e(`eP#(r;0Z#(m(###*/OU#%###5-'fc#HD?>g8HZ$EZ$+gND'`gG#5S.|1(Sv)_5#<6JY>$OG####1##*K1zB9######&z*[5N%?$uG%W[#b>$cc#5K4###"},
{609.81,529.73,4.27,0.80,"Y80zY#ZP####p>L###`r?1,#VI,###.KX1,#J/0eY#hY#5,#o839,#a#&$##(JX###Ug2[#$j7/###=OX[,$Y82/,#lQ&cG#Ce.1##Gl%$##kJX###hQ%SH%D(<###KTEBx.B$'rb#2X1^P#/m(,##6#$$##wLXD>#$##%Z#<~H2o.-%Cwm)9##C7'1MX.,#"},
{609.81,529.73,4.27,3.49,"<7,&##'m&{Y#]1[m,%>Q%rn%kL7c9'51[ev'###)$#obL###o[,LQ#BR+6,#22[:R&ac'M##C_7=|(b.[&##.,#}$#ay9-##&f2y?#wR.`###/[07$]m+_%##n,4{#rlS*##mP$[##,d)S,#=L:Y##Z[,w$##,N_-#LK6#&#{k#=%#S)B,l#xY$(##fY#;l#"},
{528.41,534.33,4.48,3.78,"3##XqQBQ%%U7(##AeB4~.Q#%4K1&S+;$(]Q$,%,/7+###9e'.##1|:Ye-:Q&Ym)&G6Yx3'##nmQ1^1Om)H6#SI,JqQ?92|H'###K##e);I#%|7/8d#<@I)Q$rpQ'7*`J,@?$2])tI*q35a-P######Nn$b[,`P####LM'Zw/>x%Cv)vS##1:e|/$m(7##KT1"},
{456.53,540.12,4.26,1.95,"qY#3J(=6(&##Yc&RQ%E[*{Z%o[-*##'o+GKHK,$}l%~KSZK3j-#lp8######3s<H-)95#36#WKSQe'2y4:y'D6'PC%vISJQ%XB+i%0L,#t83>JS(c$$##Gx$~<@@f).i?TB&g5%n8%3JSU?&X?'6l$6/#*ISU<@*?&}%$O:8|.+/H;7e,vY#4,#Ki3me.Z,%"},
{387.02,541.67,4.72,2.40,"xY$'##'D6.w+Y2B-##ym*BN(2%)[%+-A,(k:P-%Au$(0*nv+DI,%##Z}EBH$@%X)##jK/ci(j'4QH%T*XWI'%38D#$,U4###Ud*9H#;1:N>#h*X4-&p7.5Q#|z21e'='X_#%r;=%6$`'6ul%k>#JB0o#'###LH<nu&a5%%##Mh-8z4C^29l$1Q%u9/,I'66&"},
{387.02,541.67,4.72,5.37,"J?&zu%xl'_J,93=m#&Co+^16'Q%(##<BN96'xY$&##ZR&J80Ao2q5$XjD5-$f~WAZ$wU5YK&;S/dH#j`WAv%LJ1G5#y[,0I%bp6###g`<'Q$'_W,](bB4Lm%4'2jr(][W6,#(tHTl#iZ(H##vS($'5?-%:Z%od*9ABp?(+7(t/3x2(HJ1)##QW=O-)95#(##"},
{370.38,544.53,4.58,5.46,"ZP#$##3c#z-*nG$@Z$7$&t);]Z&iG$_7'VZE/,#&?%n;5$-'Rd+###_P#C%&d&.#J*,o1_0.GyUI?'mM1y'2#?&%##42Qr#&OD?###?l$<H#9)?U#$p.RjZ$ZxUGZ$g16L8&Ww-$.$PzU&-%Dz5###RG#$##/z595#539Lu$2zUA~''T/5d&h.+jD'IvU$##"},
{429.14,544.71,4.42,6.08,"r(*o>%wu&###qx$zz=######v$$hnQ5Q$TQ'###V}3eA,AR+I_&;I*qOG###Oy)q[M.,#'##<nQ|]49c#[f-$@*6l#C}2#r<?S)xR%clQ/##Uq7Vw*,u#g$$DqQiL:qP$fQ#qt5.80Bg/R/1:m)`$#4e-kw(27,0,#.l#iA$_.&P/0>10qu$@k3Go3&m$Ku$"},
{794.84,543.08,4.26,0.58,"yI.Re#M~[M5#'Z#,%B6D?Q>#^Z#<*D###X5$`?#K~0###:5#aU9+-#K~[.,#|v+Dm&_~[)B)tc&8Q&6##BtD'6$3l$2##H1;C{=Q,#E~[###-E;W5#G~[Xl#sp0###,##@J.nu%###-##W^8`>K@,#A~[###lA'lu$>~[###vh&wu'######_S%$m(###3,#"},
{794.84,543.08,4.26,1.55,"###V@#EK4###|G%<g)3l$K##..,E6'###o%#2~.pb####=%#W5$-]$aKQX>$i2>#n&t,%tv#F&0iv+###I&#pc&]96###M##pu#NH'QP/0M>4L0ec'Or#ViBs7&p,QGd#Wu%'##6AS/,####<'#Sh>@A#$m(M%#Hf4h;#l2CQ##5~.AW$qrD$##]&1@'$DI,"},
{702.46,551.50,4.93,4.60,"lL%_hU}k####=//eiU4l$A##Y_3FtF###D##_0):`B######)v&RfUD>#%##zdU2'1###^$#`q;$-'###v##}K46#$###'6#L-(S#$U,%###Vv?,d)###0##Yi)Kx3###&##LeUL5$###m.#=5#&l#+7$Q%/$R#k6*]>#3l#-}6qQ)###/.$:eUD>####a1#"},
{531.31,554.54,4.42,1.74,"E##x.-eY####(?$$w'-H&S5#0e.46#X$*)F-xY$g[$Qr@L?9%x$p^:>5#95#(s<7@*eY#S##S7V2R$am+y1$~R,cz#z6V>7%^8-&-'D~#=05_7T6#$.##^6%U:VM%)Bf39##*['N'(I8V7],k/0(m(Rn#?A/pd+<#$Ov$~w+HD-G3DyQ'RG#|?#XmU:['3m'"},
{493.34,556.30,4.29,5.26,"J;St>CKl%RH%5^2=TG-A0-##Hy6*Q$E,$n?&d(<c>#id,<[$97*rF60d)k[)R7SGU1{k#'$#LkFzd,[P#t>#jy0b6(S../,#D6'~8.iG$cL26FA`7.G>#PW:3;-|z?,##6C5jR%218dZ'vP$ZP#$##3,#.^Nnu$uQ'}k##;SQ.*m>%?g&J8SkY#<Z${}<oH)"},
{504.49,561.36,4.23,0.40,"E5#?f0Ol$]#&G>#ge'sS1Hl%Cu$6c#rI*F`<qP$###'##NBTI?#?{=y#':5#-&/cV-yA4(##I@T_[)-u#Ms-uS-muMD>#MDRn-#?PJxn,.,#z<@^z69[*3##ZAT4@'CR)(H#OS+lBTrN5U..O6'r,%KD3]>${*As5%~H'B5#_CTD>#qP#_P#3C),d)$f$U08"},
{408.00,568.27,4.32,3.46,"'##rI+95####y5#2(:######w7*%m(###9m$yY$###'##$)8###{w&tn1###$g1IHC#.,+##Ig[:w,yY$e;%Iv)wL6N?'Sl?L5#.d$UEE###1sBsQ$`aDO,$]k[yY#_-(gd$,M4KT/iC0FSW.-&p5$%M7_@+BR*?,#Lf,k]2-0D###3$#1g6o[;;n.I##UV;"},
{522.37,576.57,4.93,5.22,"I>#O5$Gc#5D9X#%Z>$%##QRMi?&LQ'mc#AdPMu#h?(%A+*-'Jm*&##G,$5D.Tn+Jy*ux5q382&QH7*T:.]D<T,%%##h]KQZ&UND/##w6+:$#o;A].$a6Uj,#e7U*d$gC9Wv#sv+[R#b8UU#$Rh8QG#<-($##Z05P,$k.Psm)l7UO''$S-.n'g@-eq%Y'9=##"},
{482.49,584.25,4.42,3.66,"'l#al%h#$M[(bP#i~-Gl#ed)A%*0T3###iu#8$'Y6)###L#$.,#m6$W:6=l$V-(H6=K]4K>#1@O>`;6v(k,#KI+LBO0J.96&%##lP#-UMKQ'.%,mH$fAOt,%@BOwH(Y8/%v$>p1AS*X27;BO^n(jH)=COmc%@c%2,#[AO?I*T]*vG%c$$%}F=COuG%0##zC3"},
{921.90,631.24,4.38,4.43,"yc(###YU7I5#`tK######Q##f&^######;$#I1<.,####/@$A,$+##@4ADR&1*E######jK)c(^7m)###7I#5M8MtK###B##:w&8p3=&0P.(Gn+0_<###1s8x|:/)^.,#Ev&=-$y'^######n.$TWD###'##yc#N&^###mY#%##Wr]###.,####Y}E######"},
{446.09,673.59,4.35,1.98,"a,#xd>%Q%###Ni+eV=######M;.bc&ZP####Qc#OJ(#18mP$]J+/Z5_,Q'##8aW~S/[P#4##@,GoP$###,##Iu#`l%s#&rb#Vn/{,#n]W:J)r[W&##6u#v_%gsH######T-#fY#%##TG#E5#W>$###AN/4SPcA3###+##JQ7jH)###/,#N?#$##L5#yG%/,#"},
{622.51,713.99,4.64,1.99,"a,#d?=0Z%###(F.uL:######?b2o7-dG$0,#96$*U-5YHl#&Ue*'53UvU'##}{UHS.N5$8##^+C<?'###(##vY#5m'iu%:#$/%-D##ZyU=8+cvU$##Ru$%W'4iA######Kv#nP$###D,$M>#pb####u_-|5JbR-###0##jZ9uG%###&##'I&.,####8u#fG$"},
{486.92,742.69,4.58,1.75,"2,#rL&=6(###,n,W~'pb#Z##>[*^P####>n$Nm(cG$###6?#wQ',D+[y7zP$|$S#d'###ec#I~.=I).,#`##;n(2H&###$##L?'4l#Pi6vHLng8Bc%[7#zfXO-&~/0#@#$2?-9+_?)G##}>&og+W:5?177].m0U7m)S7#biXg-*F>#fh#NgX^,#-H&fh#$IV"},
{429.05,757.63,4.69,0.49,"###2c#bG$###JT2Ic%.,####QD]5?%:l$aQ$X93?U&7B]CQ%-##)v$T6)###(z5/'10v'uY#kC]>$'F-)c,#J&1*T%OA]*##w%#an/bQ'.,#DC726&q[(4S*BD]p@+)$(p6$57*uq)JA]%##8(#vA4######wh&/&2&##mY#Ls-]]395#(##%6$x:.A~/hc%"},
{614.35,787.00,4.37,0.53,"@5#bu$######GN>;#$;,#Hl#AbK###=v#E0S,u#7,#[H?P*E3,#Xl$######NAZgY####1##dBZ?Q$A6'4I&N?'.A$&BZi#&rY#Rc%######'AZOG####}##VAZgZ$DI,w$#Y#%=8$g@Z%##uP$=#$95#$##'AZ.,####;##9AZf>#r~2X##PG#4$#e@Z$##"},
{416.31,799.89,4.20,6.06,",v'#########/(X#########6(X######7##5BQQK5###A$#Wu%#########u&X######9##4(XZP#0##-K(xx-iL7Ol#:K16l$#########3&X######L##h*X?K3Cl#B(+Q'+d=F,[%<286#$######0##Fh8######-.'DBFq83###lZ&iJ$Z>K95#$##"},
{601.73,801.90,4.54,0.54,"%##hu$######Gq;P5$###5##~L_}G$B,$.-#Bo0:T%'K_A5#4,#vu%95####5vRM5$###F##+K_V>#Sd+4%#ov+#I#wJ_5##rP#hc'D>#$##.RP{k####*##?O_>6%v6+,##o6(so$tJ_###U##/v($##fY#]&)*6'###.,#[[8:&0bG$###Cm#2+4<h=###"},
{442.17,810.10,4.30,4.95,"W>$###VS&['1VWAE>#*##wa2+8-h>$2T,A:Q?/.|k#{5%yK)c#&]P#cO=aw)GpP2Z%I>#gI$'9Q}#&2h1Zf-=0);%+(+>1-'2u#2l#VcOtY#P6Q.c$X,%5@#y7Q{(/;8/QR%B-&*r(SZOXl$`5%n,#{q>6H$y98hY#q>$50(r]1==5Ow+uZ&[o4XT'k>%)R#"},
{851.87,832.46,4.65,1.94,":M9######6##4X@95#$##d5$n[,F>#&##6e*aA2######mG#c04O5$E>#)##*:`95####f5#VA1j5%xP$#d%:n,(-'$##X>#*'6/,#7c$+H#e9`6#$###SU+F7)id,o>#YoYQ#$/x/,?#@aG#h:A5#KH'R,#C9`######oH5KT4Ju$iw#W?`###b5$$;%v8`"},
{806.14,833.64,4.71,0.35,"###Z>#:Z$7#$QJ1<5#V,#Q%*vy:###5.$+NWD>#$##l>6_?S&##ju$######y6T#l####F##YKWS$&^7,9])t#&Go$9KWQH'A#$j5%######.RU95####H$#TJWu$&Ce/&%#R,$@f$7IW###,c$95#$##qY#c@W.,####Q##uvR|[)z70[##.,#d[$7IW###"},
{795.93,850.09,4.30,0.39,"'##KQ%######1M<-u####>##ZV~NZ%{k#2$#_p3;_(LS~+##C#$:Z%###1,#^-U95####a##iS~mP#tc(t%#xI.Fm#YS~:##>#$D>#%##Ll$@QM######U5#UY~S'5uG%.##1d&?v74GM###BH#_5%###3,#ed&r5&###)##vq'+19######[##DP8F-)###"},
{508.46,852.63,4.85,5.03,"62YQc&###/##o/KF~?zPOA,#2{7Pn'W08$##N'0wL:eY#D##E0Yz#'###H##YGEf-:-S/M###mHW&,Nc&a-#9x1-c$###yx#Xp/0w,###$##M{:Ku@95#%##q.YKQ%$##yI#fJ2gu$LI+'A%ZG#eY#######gU:cG$###$##u.Y###$##I,#?U95,#%y0Oc$"},
{964.22,889.88,4.91,1.68,"]P#V>#[,%_p1hc%ym)OG#16%r^5)R)cG$=Q#|,'k&#SkL-H$.,#]##8R*m>$Sn*?Z6EjGW>#kEWur8]m+N,#M95nz#-@W<#####:-#~97.,#9A0Z/$RBW^V5>AWRJ'@A0bM+7&/x'%@AW.u#~P#&?#]?(rb#,u#q,#q&-&)<%?&#f#tH)xn+tb#E6#i(4m?)"},
{964.22,889.88,4.91,5.24,"G&3nP#mP$w##&92t:4p#'T##:C0|n/~P#0##Ld*######2##kdK:-(pQ'?-$[g4?q-WUXXH%uVXj@+)n+z#$PS/mY#D>#+##vUX+u#'B+^c$%*C%##%{ODC4eRX&##2B)o599I+Cc#3Z%]>##TXu5#LaC4,#5,KN>#Cg/46'OG#7##ae*s@,95#s0#%$'G,$"},
{718.53,912.46,4.29,5.09,"sB7#########;'~######]##-'~######O$##94$##9#$-##kU;######%##c&~pb####^&#']Xll'###0(#*04.,#5,#&Z#`B6#########n*~$m(###^,#m1L7(<###xu#~F@[P####2,#*%,#########k'~ZP####;H#FfW%Q%###3~&tGN######0##"},
{886.75,918.59,4.42,1.81,"95####?v&Pd%=c$PQ%p02^u$b@*h#&m#&(d$ZP#.##RK`xb#{k#I##Kn.Wu#RR(73,<J[4,#'|_Y8-`d+T5#7v(4##NK`?5#pb#A##3:9*##..,;R#iM`7]-`K`)Z#7~-{:*Y6)I##3K`,##D>#;##Md)/,#OG#&##_B,t::Yl&$##vH(_h3D>#2##EK`J>#"},
{886.75,918.59,4.42,5.20,"-?Hd5${k####Oc#vW7'$(###gC.r82D>#'##T?'%##V,%)##)2Y[>#SR+0,#E8-<z+s/YtY#Q2YR@*J[*=Z#ox3+##_5%(##/4Y/,#ql&%##8//L>#*3Uew,}.Y;5#w@(6V/]x0AQ$xY$(##_3Y95#.,####]B2m,$%e)Vu%1w'^Q%_d)vZ'Je&G6&M5$###"},
{767.56,933.67,4.17,4.93,"ZP##########8?Q######k##E'Y3l$###r%#O&YL5$###K##pb##########<&Y######H##r(Y^Q(###V/$E&Yg,&###b%$OG##########b%Y######C##I(Ypb####RK$@wPqb####ZJ$95##########A4G######2##g(YeY####2%#7IN,u####S$#"},
{967.99,933.62,4.36,2.18,"/,#c[#D?SH>#9`05K,4^8a.%g&+8e*J~-s$'$##u5#<E<QG#I#%=##YRS,c#pQS<,#F{=lW,D@,u-#JSSoL1fP#go%:B5<5#o#&$##iSS&##NVSFZ%}{@-H#_g0ww(tRSF#$zb#.-$Fh7:v(Sm%d5$wPJ###j')uZ(&Y@a5%t'+J:8Fo0WR(gu&mY#uT([;<"},
{967.99,933.62,4.36,5.30,"Fp(HV:Su%VG#p&27%)`y*sy6^X?m>%jg)&d(9QL###j-&~,$]:7S?(+l#O?$.nSG#$~p0VJ(CEA0H#_qS>Q%rnS&##+6&%##_o4<5#ZG#uo%gnSa(2C@,|-#h2?5s,1mS4##jmS+c#>u$A##nr<E>#$##t5#i.-*.')y+[%+s':z@%E20I]+IcQH>#95#SI#"},
{956.38,940.85,4.61,5.33,"cF3f]3###F5#<W+0xM}:;Yc%L?%vz3Ji8880V4<0d)LH#^,%S'2g@)`P#nG$>t9u08%l#Iu#^[PV#%PN/DM7D5H$##5eBE?'U]4zY#$[&Zm(th>Hu$$H$pK'x[Pu8,Qy5@w#6<@K^'>[PVG#fq5lI*Nl%tY#-+C3u#E>#(n#:o.|f,(f,GR'G_=W0'K'1{$'"},
{972.57,958.67,4.92,5.05,"@,OK>#,l#6w#0.,,Q$U$'%n(8=ELl%L,$K9'RFIwG#t+CT]((.Ob-)tb#9Q#s}Ha@+*H%r%'e{?Bm(mG$_D*4i<mu;Lz<j?#w,O5u#h#%^R%j-OP,$YH&R?%A-Oql%U0.UT+6n-EC-7??(L,55C@u#Y>$)##hvMd,$(Q%6##c/On^7C7*;$%:]2n/._15rw'"},
{37.75,16.60,5.31,1.36,"W>$######*##=mU######)(#=mU<%#U1=u&#95#[%#^rC$##hd,######6##1nUf5$6#$A'#;FGf(%BGNM$####Gw#Oy8###,~,######$##qrU_96D>#4##mx*qrUI..%#####YcCQu%###mP##########|n$je1######v##n[U#########6N=######"},
{481.45,115.01,4.99,3.90,"7VNL5$######w=1+u#/,#%##^6'###pG#vc'>u$###-##=~,~~AA,$###OS$,XWmP$###X7#lL8SH('##N7)a>$ZP#$##q?(vm,G5#hY#UI='gUBe.fY#fL'Ks1%kJ$##7c#yg3;c%###'d%$##DR#Gd(pz1,J*&%(qP$My*eA/MQ'###5b0_&3.,####<j+"},
{100.09,119.14,5.54,1.32,"Zp80:'yC:A}1Vu%l0$&%R`R)Oc&VG#9l$7),O#%%H%06#si;Bg4UwA2e.GH#;,LrC&&@+J$#=L7?u#Z>$U]$+S.q-(p5%1^+wl'];)Hl%Nu#mSSs5$.,#b,#keOim+:5#^.%;.,b<9--'0:.9#$~[)pG$0I'[8/-)=###%?#%h3oRS###8##$Q$0US.,#'##"},
{533.24,119.63,5.75,5.23,"BI@sn.###$##|tAQu$V7'nu%I~S~P#$6$_7&{{U######7##x-J;5#VG#0,#KMW###nG#&A)P'8###A,#u</0a=######k5#(|C(##i5%T7#%JW######Z0#h{A95####<$$Nd)######-Z#6D?###V5$Yh#7IW######g(#>x2.Z$###o##'c#nd+###7,#"},
{619.97,125.15,5.69,5.66,"v8&XZ?_:;###]e=LI,E>####lN'I..######Z2#cA3U##ZP#JJ,I?%FzIH8004^.u#5@$KC/ir:g6*###1##/X(H?(###%##i5$37)Y|+>/1Cq8h-(&0+.`5OV2;m)###<##v%-[P####8##7,#)v%[-)###El#?-&4y/$Z$[D7m>%zP#Af)1.*yc&?I'Nn+"},
{437.68,131.07,4.93,1.74,"###`##+A-lG$U5#|/-4v&N%.(0/p@-7##MB2W>$###1$#N<D###f##a:9{k#v/.av9*oZ8Z$HtZm;9_H'q[&,f1L>#:.'S7,###3##I3;xY$'%,Nm#LtZOTUMoZZ,$=d&:BD&f0*l#Wu$Xu#######$[#rZ(######=R#(QL3l$###,##ts<+u#&##=#$x5$"},
{206.08,177.62,5.22,1.49,"yH)+u####*##d[X###.,#I&#u[X)7$wA48%#,:4_Z%AH'F##KT3ZP####5##d~XIH&eY#D&#9,LhO5?z<($#y~,rTVCZ&&##a&.o#'###%##UaX.bH###6##'~(UaX+u#$##$##0`X######K6#{v,######8J$1o2######)6$|y8######7,#Iv(######"},
{473.57,189.27,5.67,3.60,"Hn+######3,#-49sb#F>#=5#I(NQ|F######oe&yfU+c$###?<A###uP#{7'nR)0-&sh5N.*CiU@[*p5%{G#AYE47,###0###fU>A0i>$}?$h@+;X?1t=Cv&qeU?-([?&%I%^uO######1##j@(.OCD>#&##hm$&eUOG####LgU0%-.,#4##cPL.,####,##"},
{748.07,214.37,4.98,3.12,"###5B%go5.,#p81F_5(w,Xv$5SO'$(###^C&}O@&v'Lv#UUOxY$5$#A=IW>#J:;3##_d+#<,iSOOG####rN/(BI^aHtY#w,$(?&'##rsF)##nnNOG#2$(dc#kj2?o3###p6%O]&t4M###-##8#$###rh8#Z$gf.eY#gH'X.,Xe%:L96c#ul':$#uPOTG####"},
{748.07,214.37,4.98,3.72,"8$&,a2y$,%##N:P5&0.,#p##XIAL6PmY#t7'p5#eoI,Z$Uc&%@*--$,GDX#$P8P)Q%'c#RS#gSF&M<###c$#&-$0&P######Ac%###;j7J>#}P<@H'Pv$vG#&()q>O###5##+##a-HhY#.,#E##7#$}n'.,#0?#L5$/y%W>$h%'Q@-V,#pv%o5%7I(|G$BlD"},
{211.75,222.43,5.36,5.83,":#$######3,#/}I######?##IK^eY#-##mm%yJ)Nc&h##ge.wP$.,#######mJ^######Q##uJ^###<##`W-$I*.,#f7#<XAc5%######*##(fZ######/##5Q^.,#+##3-%mE2>u$W%#z~1Z5$######%##C:6#########*=UOG#######=a(PA2(##D>#"},
{605.69,236.68,5.43,0.57,"o%#M:<######<$#mlQ######&$#YlQ######t$#d;B######IB%b~/T,%###W?&XqQ#p6%##dV8'nQuv)H7'WV'9mQTo1Y#%+2'Ax1r5&)##6/*?B.m$N5H&B.Lx]/Ap19'2Oy22y*:jCGU/P8'{n0.,#+##7z#tU;lu$Uu%/C(66'=Z#V2?,-&$c#R,$v}>"},
{605.69,236.68,5.43,5.01,"Ne.{k#$##0o$e?(%Q%/##}E0XZ%bG$###Et4Nc&######Q1((Q%mZ$sH)24298(KPFRH&8y-*1O?A1.##lv%){?######70#+@(qm$s3F[.*y^1+^2oK,sQK'mO:?&TI&|{-QsG######*&#=?$EQ%|cA|H*-21uU:qd']O9;lO($'U5$+q$kI./,####W$#"},
{401.73,238.92,4.95,2.65,"]m+###:i6wG$9IVt5$n,&56#1e,mj6/o0N>#Z-#sf4>n%(c$mw0E##SIV5##?IV1##~w+$I$j$*[S*n1O?6'jB#fA2T/'du&481&$#26SW##qIV###Hu#Jo%&03###ID%kT4i5#pb##D#vL=Zp9S$%ed,a##-IV###-##5_)p?)###@%#?C6Z5#D>#9'#WI-"},
{401.73,238.92,4.95,5.06,"dM<(##FsC###ECV4,#J#%###$`70[(<c%###kH&?)+bn0###'B3%##&HFP6&5BVRl#Ul%Wd$5].fH<DI,$##Q5#i%A+u#[G#;$)2##lf(U+DPh=:q'^#%5{/|5&LE*s$,s5#OG#@8)OG#Rq(h@+:q,(?%{m+yY$yM'O#%#?%7.-Y%%###'C#T,%U,$###)W'"},
{353.21,261.97,5.11,2.23,"3|4=l$]v'.,#^/.###4`0#l#######X_/.,#######Ll$####24###ob8###@gT###B&'e?'nl'###O[=CH&######Cx+###Fd)###qeCeH'*eT###?n%bM+Vo4###OfF+%'######sQ;.,#CH&###{y)`/1US12##}H%HH9r#'6R%v*;<$$%##lY#)fD###"},
{353.21,261.97,5.11,4.73,"=.-.,####0##J7Y######5%#2aFL7+J5#hf,bQ&?'5=c%a-&:K5######R##p7Y--%pb#p$#gD>|fF1I+J,#4ZAKwOA,$/##M:<######E##@8Y:.&4[*Z##IT3YW(%-S5##{7YS-&-H&V$#}T8######U##,YK/6#,o2]##NK6{-#IZRG##H7Y+##>u$V$#"},
{374.38,266.20,5.50,1.82,"Zc$d10,U8Ow&U1.9;:H82/##/l#Lu#oIX######'##`-*###tv*jH$eJXI5#2KX}G$l/3J5#m6*'##`JXB5####*##g./###t?*###QKX(##vIX###>C7q6$T%/###4LXSZ$###%##ox3###s5&###1LXGI).cP###+g11+.X-*###<LXE[%###%##N^8###"},
{374.38,266.20,5.50,4.78,"^m+######4##|7Yu>$eY#Z##^)?552-^8+##Z8U:`7o#'7##p.0######3##&8Ys,$SH(H##dy7Ey$G7Y.##P7Yr>#Z[,0$#=A1######F##R7YW5#J7-v##vx5V-#w7YZu#L7YO5#yc(5I#]J2######>##oZT###e6')J&/v(]G#xN7o;:w2CYe).l#]y)"},
{390.06,269.25,5.80,4.95,"(f1######&##ofWOQ$/S/###@D<*7$pdW###,fW4,#XZ'###U~0######,##eeWF5#X..(##9h;O5#?fW^l$xeW-Z#Bv(I-$3n-######)##}dW###_u$B6%:96sG#u<3kGFh{BU:+{,%bM0mH(#########BgW7-($##$Z#;L3D6=EZ${R,Ru%&k/vb#T?%"},
{606.73,333.59,5.29,1.67,"<$)/,####Z##j6)@-(###k##XJ/vG%###W$#uc(######iA$Rx+2I)5?'6##{d'.??VK6B,#HXZof4fY#K-#+{>######V6#EhLR{@J[+%##Tz6]%'uSZ'-$5VZ=Q%Sl%Z~#_OH######~##>UZ'-'UQ'[##R;>y#$oeV1H$}UZqY#kc&_H$1}D.,####(##"},
{300.71,392.00,5.41,0.23,"Bl$pP$+##c>$GQ&#l#E5#x6'1,#%@)+?$mc'4##87,iP#+u#KI)en*L?('##n@S@d)[c#u14.R'(8/Ag$&@S,Q#J817{'3PLP%-Qc$^w,m$'d@Sk^52e*:A')`((ASo/+3(:*D%e?S'=/?B5jv*###(##N6%C`;kP#R@'d-'Ue%7S.4i(H?SY0)by8yQ$L@S"},
{300.71,392.00,5.41,2.03,"yu'%##^K5a#$q%Ynu#ky7Gi&*)=,q%e%Yuw%d&Ye6%295[,#T@-###m.,'?$,'Y3@',z8*.$$_:Wy)()YEg.b#Qz~(mv*po'd-'Ol%vY#jP#rW;5dJKv(Kc#gd'r(Y*9/;@)N$&[*@J>#`l$WZ%*H%ZP#*##*##=@(j@*Fl%)?%TR)p6'g/16v%I-(,u#(c#"},
{300.71,392.00,5.41,5.06,"95#D5#)d%X-).v%WA/2$&z[+Je*p#''##bd(|k#-##x5%Lu$ZG#&d%l?&:*>{]1{m)tR'39Vr?)Uc#eD7Y.QmY#Q>#Wv&r,&U@,D1(tQTcw&M9V}o,=:9qz(9K3HA$t7V(@'1w*+?$3d)$##}f68?#.8V2~&>7V^:(Ky6sh%Z/2W2'-7VHc#C81yY#_5%+##"},
{334.37,398.56,4.95,2.69,"###+$#r#J<R+xo6S##%8(`8)TQS=-&5v'6W(>&%S'5(7)~B0###N,#iy.U08P}K###[N-#/,adUaG#:D0FT)pT)?R)uLNz>%>5#95#/?$'?&i83###~q*bZ&eiUZu%78*=-&$y*'$':)R*d(iY#.,#%##4u#$Z$###rZ#F,$Wh-L5$iQ#-?&A{.QI)47'@c%"},
{304.78,407.96,5.40,5.33,"OG####N5#n-((l#m#%Km&tV<d$)Jc%zZ&0RCK,$Q>#Iw*Q@+m[-$##[>$]S%yw+#g+q(>QL2aLUS@*v:9xK/V6'h#$NAUb,%psH###Q6(]H#G{@8-$xJUI6$3JU@@%u{?2e$CJ.zx$,JU*c#>a;g,&mG$'##Rp2Fw-?mBU[)9ZM~`+r&18J*+'3hV&5U9###"},
{427.98,410.91,5.16,5.27,"MLRK)4$f/,R$Jf1[s/'YFvY#Mo3H>#WG#T#$###$##N>#bG$+G=3eE.u#Cd&pIR-k:eY#d##Kr?6#$D,#]Z$$##%##%-%)c$gv*W:1-u#{o*<JRfc&###N_.^j9ac'+##l9.{G%###5$$Lw,~P#]P#?5#l&H5.)dH'eY#EMRu/0cG$.8)EMRTG#~G#[g0rI-"},
{412.55,412.96,5.24,3.14,"E>#q##>ITF>#Rq;P5#2T1+I$g7,l'2g:3<:7(C#b:<t$$L~0###x##(-E=y4B3D###Oq%>M9:JTAc%2y%mLTVd&9$(|_'/IT###iP#sI(~K3%I*I>#[1$Zp7e9MBg8Om#nn-pO75NA@0(dK5+##Uc%~P#.,####Au#Ol#{k#)c#W,%oQ#bQ(7U.O[+cG#u,&"},
{320.56,419.64,5.29,1.85,"####6#N2=_d'(##6?&+d'U9RN%+<c%0%#y6RuD?gZ'|G#HA*KQ'U##Lo2/=0#I'6y*;{>q.CoXA2n*-(+6*C%6&PZ%#46^m*=6ReR$,T3PC%}y9d1%76RX@%]6R%%%]A2z-#dG$C-#5FE.,#yU:@J,jT5B?#6%-LS)i|ANg-:B5r~'du&4%$###r-#m[-###"},
{320.56,419.64,5.29,5.39,"0e.$##6Z%z-$2f.nn(G*E,R&,fR-v&Fp3P6%<-'zP#fuIoY#KOG###7v's>#22??,#1gRru$4eR%~%:V:M$%?J-{g%&dR@5#yW5#[)xb#$##+^-YM>Ea5fd*4q9@56<K0D&-]T142'%g7###9~(av*^Q&8A)-6#$JRI'3N#%1,#rLP`c%r?)i[)C%'pb#t?&"},
{381.12,422.24,5.92,5.19,"4##~bGj.)@813/.Zf.(##:L2u3>8#$mP#|W+z,&>5#|.*f/..,#&##E~&PM?4L3R~($_9'_7';W6~)5I):E)K[)yB'/FAtm(2Z%=##,'4E-(6PJ.%$V7W4H$Q:W`B..'3*B%t@-$z*vFFf,$p,%>-%AA0M,$bw.A@)ePDFe,Y$HX=>#Q$37(N{>7h2(c$T##"},
{461.09,434.62,5.52,5.29,"&l#pZ'W>#?M3rd+jY####Fl;a7+A,$Fl#Ac@8#$%##h-&l[+~#&$##mG$pP8^n*,0,dx4o*4)^SJ[)&y/#D4Zl%6l#lZJt,&>(<%##V[+#A$:2?wH$j[Si?$e[Sx-%'h9T%$L%-#T$R[ScG#nV73l$[l%<,#M02RQ'|=;[B2fPMED'1T25y,tx3z(&6962##"},
{419.37,448.99,5.12,1.53,"{$#EvQI?'-u#sQ#,vQ@5#|5%0I#mvQ###%##mZ'/o+W,%w%'.##835/@*Fw-hI%ruQE$$UI,S@AW&41##q#$GvQbc&.,#N1$hG#DS)I281=GO8,9%,C/$_(;0vQZP#n##rV0%xQ3n+~#&E%#4v'7;';vQ@d&pYJyQ'],%lf&+%,kG$Nl$f&,_f((D<-I)B5#"},
{419.37,448.99,5.12,3.60,"Nw#qeOf?%7A0###mgO}@-3R*i,$Vr7R:9[u$96%I-)E5#ZH&T,#s:8PC6&Q%[#%q]KIB6%##:E:`+DHl%4##8d%smN36'fY#&##[5$4~Cj>%vZ(R-$DfOvY#c,Hlu$x%-EQ$37*M/-(p/ZROXe(s[-P*7QQ%_P####O^Nkc'Ro,.,#/w$iC;;IF>u$3##j_6"},
{477.97,447.54,5.51,5.41,"M$*%##Ml%rg)x&0on)5{?y[%>xQul&k]0_Q%nl&{P#4ZHWG#jD@###Iv(dH#]_>@,#}xQR6&;wQg[$}p8Gd%aS.i^%puQ'##^M2iZ(zb#B5#OC-3g7gC,^05Py6rt5H19oh43(:V(%`:=J##+n(AQ&Q[(PT+c6#@HM~93>c%}b#VUQ7-'Cd'Kx3e6$uG%g@#"},
{384.17,468.27,5.16,5.24,"8/PO]+AB3vQ$DR)@2*|-Pg#%2H&@5#yu$Yw-E>#'?#x4J$m(.E;5?;@?%PR(ClK4+3I#%L##[x2Ac$,R)t5$E>#Y##~-PB,$Z?(C4?Q>#9p-e4E0v&###KE8dA+:[(7I)gT2%##%]*<2:a5%PG#?5#%c#4fGSI)#m&kG$B1PSd)E>#`47y-P###nl%qdH0Z%"},
{339.86,477.90,5.20,5.21,"###$##ul&-u#L,$4$&b7+_A0Zn+Hv(jl%|ZA,I%VA/(_8w[*r5&###:Q%fI'+'0V/*WNB_o/7qWv$(?:77227l$`-#zoWxQ(Vq=/##dI-|l#rYM/%$FnW46$KpWk~)>`=/o%yu&sJ&7pWOH&}B3Q?(`H(1,#e]3|%/j&R@I)(eM4t7L-'Sd'm/2/r.2-(D,#"},
{461.21,480.49,5.37,4.99,"suOO]+K~.(o%PS/2E,UuOml#BH'9H#(S-G'0+u#:##w>%=};pI+*zO],$.$&asDqt6L5$E$#__=YG####T/$Kx38$#2)@aA%a-(?2<-##Mg7C6C|>%###_y1f$Eqb#B5#u-'KO?l>$|N>P$(]6'v?&;6'v<DC);|5%b5$;yOvd(sd&-i:E(7.D'-FEA@*<c%"},
{354.05,487.64,5.49,1.99,"###HS+S:8Tc%P>#]].c?(l|?~'2h?)`,#bZNN+B;n-=u#@-'@H'],#=C7==5;w)D~*c@PIICf~P[d)DB0+*A`u%Ic%U007@+7oUP6&]@,tK(~i@<K%HoUg8*(nUTH$t;BxR#c5%|>#RC;I>#TvOS[)9v(h6#PA0+:/joUdR(%M;]S,H$(J@&lY#mZ$g,&###"},
{354.05,487.64,5.49,5.14,"g,&###wb#mZ$a6(TI&Qq;~S,yxUVI(s~0p'/hQ(d6#6IP0@)RC;I>#n>%|>#G`BnR#5wURH$~xUf8*ir@#B%H.,vK(DxUB-&_90+7+Tl%=Z%MK0z{@UJPO[)s[OgmB=w)Ne*jg7ls4*6'[,#I#$4$'2}@#J.],#YQN/L2R-)L-()4@Q>#_].=(8Uc%###vw+"},
{413.98,492.29,4.99,5.32,":#$cv)~>#a12P6(tb####|#>2-&B,$(d$1$I:5#1c#JD5vl'8Q&###sY#Rb4eR)|[(0T3[s2aHJT?'+g+a:5mP$###rTJ0-'RrB###+$'J/$m^:Pl#CnR)[$>nR)$%Up4L.%Td*D.$ZoRo#%3}@6#$>l$2###g2Zu%DIHq#&g@PSo(O~,,R'F@*VM'`:=$##"},
{431.58,503.96,5.23,1.84,"###w>#Sq8<['&##]u%N-'unOmm)1Z%/%#;$QpD?du%n>#-8)M$*V##]%.jX2)R'2~(]<BoSFKs?^-(Bh,5sC7v'qZ&x|6Gm)N$Q).%5~.IU%6^7Bz$y#Q-%&($Q;R$xB8k6#cG$~$#2uJ$##86J8n)5n-)##?n-8',&*@x.)|T7JD0}>&Am#eY#?K&{v,)##"},
{431.58,503.96,5.23,5.20,"hv+###_P#:^(9[)/R%@r@40*JIR`H&D&.RJ(p#'(##_/Okl%F,M###Zu%_6#B_=1##eIRTH$WIRT.'/C79e$c[+3q%4IRPu#It6|u'J,$S#$up-OECCj;c?(HV:ec:}-(3w)5K2T{)v6+1##e6'K-'tT3Ie))$#AIRz/2mP$/l#G]PZG#SH'c&.)7'###jl%"},
{492.78,508.27,5.77,5.09,"2##|w/ET,c83pc&Ye.=,#8PGV(39#$C5#v5@[U/]c&;@)EJ,h,&/##e?&y:6Bo-gg*@D>POA-VTb~*fm*mF:=6&df&B@REe,A&0aQ$~w/}u#ljGUz$,RTqH$OSTtJ*0'5Fx#C?'gJ&<QM(-&p-%q..R@,.u#7f0~8/k|?Bp15uI<QA9l$e@'m&3#U/s5&CQ#"},
{416.90,536.94,5.09,5.01,"vGM;L/2d(k7%J~.9a/>+JWc#M..iZ#~I-Yx&;c%3##(-'1O3t6(/hP@5#r>$}aGDu8ZP#F$#scPZG####8&#5iA[##@{@iS$tu&Ly5X5#pf3HP@Iu$%##aA*VfI|k#M>#hm&Jc?N@,YV<6l#*H%h,$T-'ZJ1zD@E#$`P#`s2-n*tH%V:7X_3:x#(>I^A1OG#"},
{601.25,548.28,5.74,3.53,"5H%###hp#8?&>`;###Z##k(-Y@,###8m#D'2A5#-##IR%rm,LQ''##Yo&@m&5JYM>#'8'-&A5@+t>#ZOY@>C###8##@g5eY#'@*O5#u#&~G#oNY$e'U@,iH$np3:<(cJYn>$###[##w84Y,$YH((?#(.+;,#ZJYK[$DI,x##mR,aB$Y=K,##6#$5##Fl%W,#"},
{493.03,554.54,4.75,5.15,"~zQzL9HH&Ov&'10dzQnf5;,#c..=v%Q6(pm'G7-W##<e.?'+x[*}=?yY#fv(>wQtr=###%$#y=F?6(###76#.^2JR%-`AL,#'6&3~,<,#CkBmFAZu%$###QG-`-[e0bG#P&1|d$Mq:Eo1D>#ZP#cG#1u#oyQfR)}Z'Dl#8xQJm(wb#Ak3EvQC,#m/02|=pb#"},
{480.27,558.12,5.44,3.32,"3,#>$$@%T###KC8zw+LA1PH#S'4A):ld*n$)9S#a3FEu#o-+(l#r6#}$TX-(Lh9+##+p-+L5[%TB6&x%*2(T6~(v@.R_%H%T.c$4n%*2=Qy-R6)###FS#_6LI(TJS0yl#Nh5x(+)%Tb/'V06nn-1S.E>#_w&|Z$~?)qH$?K2sP#m?&Q90;Q&)I$d|<ml''##"},
{480.27,558.12,5.44,5.13,"TWB)[#AYLTl#I8Sm%*]p5NU&Rl%)]&S8SWQ&.,#'##w7'y^9%q8?R(|}D8n+J;SU9S/v'-I'`(9+ICQe/{,$0J/O5#CQ&&])xP#BXB/8,Km)I'5eW:2,#Ww)Y8ST[+###/m#tx-Od)7[*8##=5#]>$[.'`EA2?&F,$0,#{8SZJ*Gl%7$$C7SvQ$:8/bw+?-("},
{682.70,558.56,5.31,2.60,"VG#H>#5##~@,>{@###4##&~&@7X.,####`4/+7*_Q'###J&K######J$#0}A(nX###v%#+(5j0ZGu$%##`4/mI+.y)O5$jP@.,####l##-/-7B5###*%#0)4?1Z~,%)##.U&&90l~(;c%-########]e(E05`$+###&##J}.dy9&##.u#$r#;Q&=Q#aQ(H,#"},
{445.69,560.18,5.20,5.25,"H>#$z4Nm&^T43&.9$'$##@C2><8Qu%wY#p)*@]'_96f]1-w(######`6%h;AH^4C[&XB5Ig4@BW#$%L7+9{(tl'FI$JBW_v'|k#&##U.-u-*eXFhc#V@W#H$@DW&8+aT2T&&@v'XT)cBWsu%G>#Gc$Bm)/,#>~-tu%flHX$)OKK-s<+-&^$(Q]0D+<'$(,##"},
{461.38,570.55,5.72,5.03,"######p6%L[+4Q%/[$U^5>B5x.Wr?&)R)k`-@o2>n$#YF7C.sb#(c#,$'R#%F'7X~#C.Wt#%^0WY]-|97nM&d~/Yh,S/W?I&N5#LJ,I?(###/w+pC7RGEZ7,x1WmATQ-'f%+w:8$dF%[)RZ#Su$ju$p~/~G#o>#<q8.U15~.1A.pT2*##5U3L|>GZ&$##)@&"},
{538.12,570.22,4.97,2.31,"D>#r>#krC]G#ay9A##v83qz%ew+q/.WS/9)2-&*??'$S&ym,vc(.-%8U88Z#1mS&##~B1&i'['5gZ%kqSfd&*oSV5$&B1(##i-)dw,Pc&*##7rS)6&=@+'H#vL1Q6'1qS'm&{aIsQ$$y0@w'~?$_g7/,#Hu$x=32v(95#ZG#(U&jFD}@+9l$X>$Db:EH%Vc%"},
{538.12,570.22,4.97,5.21,"Ac$>v%,6'Ys;m^5(v&Je&RbG?u$ZG##/DVm*.##FR*29-tQ)Ko2{@%VcR^c#rSW(I&V'4bv$=J0[?#{VWQQ%pb#-##U/-Bx1b]4L,#.oVrc&$TW+y(H80*%'S&1^i(,IV8##II,tG#%-'mm'|$%eC;~['<m(Td*TZ7v#&aI)_g9fz&}>&=$#?A10,#)##n6%"},
{539.87,587.93,5.15,5.20,"HJ1$##7#$_J&+A.vn&I2Aqe(ueUkv'_o..~);$).##bUS`l%'bH'##Zl&M##UD?+$#4eU.H$YeUr%'{B70[$G~.')&odUfG#TL0`#&8#$$##K0,D^8o27_n-K{>ob4Td(rx-T96{D()R*q##G##CH&.d(.,#z##O6N{[+*c$.c$C7E4u#qm&{o6=-$95#:8#"},
{745.19,636.40,5.17,1.41,"E.&=^7K?%7H&`8.qV@###&Q#@f.5~./,#l$(zo2tY#kY#W/,87'xUT-%*cG$);~(/00,#[f'Y82###:e'LlB-V<###l,$+1)NC7j{=.,#|5#(;~t5&$##>K#sjG$##LR';s6Jm*###.c#1,3ID@###0,#K9$%8~######sV#_@.###(##aE./$(###%##R7'"},
{786.82,659.96,5.88,1.86,"D5#t'798(bR-6##o7.x;*WI-G5#mP$C1#MkK&##D>#`8$.<Cv##H,K83,<h=<I&-N7RnL:6'l?I>R*y/%3<@dc'###R'#ycRa?#%fR@Q&/,#xm)YV.XV>I~'}dR>#$9l$''&Dg8###.##&f+Yv'n_8###1##rA-R$)E5#R%(xdRA,$,##2v$_q>###$##C-#"},
{786.82,659.96,5.88,5.29,"###U>#rcM###.##fv&%r6^Q(J#%#H#o[(*p.ZP#0,#K5#W_/###E##}jD###r,&?%$[JPr#'7IQSv%g-)hg+W%*On/B,#Y&R/##W%Rmd+###&v$<5K^6D?[*V}>_7.q-&x(9GK(b#R###690_,#*%R######F/#2$R)l#{k#hU#l#R3,#6$(O)*0e.###%w("},
{804.08,667.81,5.68,4.04,"q_5eY#uP$6##|=[_>$###&##?&-B)2ZP####8##<`3<$)###z93nc'&I*5##m:[l,%+c$0##4z7.i.0_9*.)}l$|K4.m&<i:xL6r7-:$)4##f8[pb#N,$i?#FaG%##J9.p8*Vh=###jy+O}>-U8,Z$~#&_/#t7[###1##yD%g'8###te&0w'40Z###JU/n>$"},
{462.15,720.66,5.34,1.94,"###.##3a;L5$WZ&7c#~V4%;7X%-_>$@##x?Lf>#>H&+-#~5O'H%J5#$g-rQ(i.(#0/[$J:T/C)WF.-J#$:^L>V8?80[y$>(W~['`R&(S.yb#4d(R@%bdQ~y.w(W(_;Y5$Kz'?<2j%WYc$M6'3$&>Q%Dc%$##hG$[G#)d&]6(?Q%xR,#Q$Du#)##)U6{>%.,#"},
{986.45,937.66,5.57,5.50,"2#E@R*7H$$7)C*?1Z#?->,e,Q9PE5#kd)Z>$,-&N5#]v'n,&k6Pl#&`P#).$GW<lr,glLIZ%I8P&R$c..p,$4A.L@$n[-$##7H@R5$H>#Ql#hq6z$''E/;`<16PTZ$ke(1p.`..a^#Q;A'##z_-e#&E>#@5#*j/Cr<76&X5$g`=Gp.2d(kY#HH&^L$l2C###"},
{1014.90,983.23,4.93,4.50,"############0*E######n##o%Z|k####D%#fRW}k#$##vy#.,##########O[S######6##^&Z######YI#)hTLQ'###NJ%/,#$########o$S0,#.,#A##f(Zwu'###O0%$:1Ph>###Ra/######;5####>]-E>#VG#PG#am:x*I###-##C~%[%Z###@,#"},
{33.87,71.49,6.66,5.93,"D>##########`f]######)##ht]r5&###,##y<)kNF>##=c%95##########Zo]######/##eq]h5%gY#.6#&T+lC;R5${,%.,##########;8Z######A##pr](f.6Z%gd#3<5]B1yu&q-'############Cz8#########GqS?S.E>####eq(i&2G,$###"},
{172.57,76.86,6.48,1.71,"D>#:##_f`###&##o8*n?*###E$&K'3###$##<u#>?'95#$##|k#{##Tf`###aZ&[t,T>OO>#:_~dz4]l&g>#CR*z7)U,%'##oP$T##uf`###Qc&j?#[r]Mx0JYGan-ZI)co-06%xb895####a5%$##_g`######R##[;64-(?5#cd&}`A'I)%##nB):o,e6)"},
{172.57,76.86,6.48,5.22,"@h5cK/YR,e5#Jr6:e.###?##?d=3I+######5EU###>E7###{&4zY$zY#wz,%T+%;7|_:07*`DU=/0-?&Y5#5EU###r,%$##]./###d,$my'Am)P>#GKGWK0*@U=5#OR''N,hCU######<##1Z%###:5#k-%Y>$$##v$'nI+o#'###P5#`f)p/1######S##"},
{67.34,78.41,5.85,5.23,"*A0######h&#%^7yb#H>#hW+4'3WH'ub#)h&Ib:>u$###6,#|QSD>####x'#4qTA_<^P#bf%D;;C8/)m%BU/_26D,$g.*}6*|&5~5$c5%cT#w&U5w,bG$>##/&UQd&Xe.W.&dh;&8%?=F^#%j>%+##5J)]n)_$UC5#.d&}c$[%U{[&`8-;o(l;@C&&(01-H$"},
{753.74,96.53,6.01,2.40,"Qu%######T##0p7######}^#zB9###$##q_&-H&0##PI,cu#YB7######/$#J7X###.##bW)v$,###Bx&n}9######;y0t>%CM:######.##%=X###8,#;$$},H###:r(%%*######PE,3l$^?&#########)(K#########5pG###L-#|G%Wd$###pr(6#$"},
{508.49,102.60,6.49,4.90,".w#9{?###$##y6%d^3###$##YJS,R)######RJY###ZWB###,Q@&2>###%##@W:YJ.Mu$F5#,KY2c$;5#2##$KY###c;@###|IYZP#/,#[##nFK0,#cl$8p-+JY###@5#$T':KY###pT6'##)sE######u9#y^<######ag*dx4######=7'KJY###QI,&##"},
{67.18,110.02,6.58,5.84,"dA2######(##[)Xj@-D>#f,#957psA5l$*Z#U>6H$(n036$'W~0#########c&Xe5%_P#3,#RE42e,w6(U5$g0I5-'vZ'iH'@.-#########]%X###3,#hQ$QL:###Gd$f#JA5:.,#y,#<}G8m)#########+%X######JI$0vS###Rv${SOcP#/,#3x%^D?"},
{732.60,123.13,6.63,5.28,"Y.+N5$pn)aZ&WwK.,#vl&###&{R###tZ&###SA*_P#&{R###Zp6(A'C19l#$ZvR@Q#jw/6##JzR.,#.,#&##La5Z-&.wR###?J.}~(_~)ud(yuR[l#-Z$78$]h7m.,.,#7##&0,LeB-06###Tx1v5%#^*1S*;f3AZ#:?%:~%qP$pq&iZ(-##%##B*-;c%###"},
{476.04,143.07,6.73,3.39,"Hp*N#Op>%@,#=&Sm'8T6(ic#;J[######Y$#nv+E>####B,#-9+L^8######%M[iZ(###-%#]J[nG$qb#b$#[?'.;9ac&hP#W.-)Z$*##DGIlL[######{)0BK[L>#rG$4@$@#$DR(?&*}n0F>#Cu#<v%3K[%OB'l#2,#C&Bor<###)##AJ&fY####2##9I*"},
{708.62,143.16,6.05,2.35,"7l#M%&$GL7,#d(*bE=Q0756$-6$qw*{?P1Z$?%#v$+M@,###6#$E##8.Q$c#OtIO#$>{;C{*>v(bG#Z0Q}w+W#$'H${`C)6&7#$###,1Q0,#C~K6Q%q{9Ec#gJ*td)o.Q@-(eY#<5#)(2&-O$#####*w@###Dv$?u$EX0r5&[@&/J/l$&GJ/<h:[u%R?$yD9"},
{708.62,143.16,6.05,5.34,"`A&`5Hen/3u#)e(ZS.qn)mm+Bj2Oc&L@&m>%@_Q#########ez5z;BE>#D##)~QhZ'Se*s.'-E=KZ#%]QGQ%E]Q%##X>$###lWD&l#Vl#mQ$>]Q}/->v(OZ#'XC;N+v2B%H#xZQcP#pb#^##5o2###=~$jl%M5MYl$^%&Rg.yK8S?$RC)uq2>|EwY#OG#eI#"},
{732.71,150.74,6.44,5.66,"DPBOd*]Q$&[(}p1x5$:c@FQ&'CNxb#nu&###nh)(%,PG####7'O$[(tb#L>#kL2T_-}tIKv(g%OQl#*@*{l%OB2l%)Rd+%##6F@fY#s##NR*F%O^u$cK(66IU#O3l##0(or;9&1Ee#N)B+##3@*###IT$~d+Zd@Af2c6&fc&#9)d;;Mw*g5%^H%6g%l2C###"},
{490.57,189.83,6.78,2.07,"hd(V5$?v$g$*6Z%###8Z#ur>eY####&##%W9######$##2M32F=(c$)##Ec$?$QU5$~5#'oM1c$aQ%E&(z%V###2,#9l#@:4l2B######Bi5O&V;e+?u#8z00['?<7N7)*{8%-'X5$)l#lI?A,$###$##>(V8?'=5#F5#&(VSu%XG#mG#='V?%VD>####Hl6"},
{490.57,189.83,6.78,3.67,"#D2######%##L^2###3,#8l##@<W6)c>$(l#654qRU95####&W8eY####%##MSU1c$eQ%?/(z_9Sv'q{5F7)cTUk>%M>#cG#ci=Uu%###tG#$]N`ZOU5$W5#AL/+TUU.+1l#|UUA6(=5#B5#pv*^v'J,$P$%x5%Z|;j>%*##rL1niC######GVUbG$######"},
{529.61,189.55,6.03,3.88,"~G#_?(ZP####s~'ll'######QzQbG$###7##$RIBD9n.),I',##VED^>#Z>$)T*^(>###9Z#zyQ&S,+u#I##N]*HyQ(c$###0,#u#'QC#:GL{?)2H&FH#{E7P:9H^.r?*u$$]P#2<0~./#########0(#YuQ1,#3u#dC$tuQ[G#-f,Nz&)}I###uv&l$%Qz<"},
{707.58,199.55,6.66,1.42,"5<=ZY6e-*ju$oB7H:(fw/I$$q~2###JZ#TN1,u####hf&&SS7R)>|-E@,/u#.RSeH%.,#M$#gQS###&##c:&fM@Jw$&^2(H>q./{w(|#%el%(US]P####E,#fTS~m(###X$#7^7zr+L@-l##K#%###Ge*.7*VW>######vI&l8MJ]2###4?#/m&#_1D>#'##"},
{504.65,202.86,6.61,3.56,"H.*aP#eY####mRI######/##=|=######`>#2)S,u#/,#3,#V$'ul'95####*BQD>####3##N%S/,#~6&NC*%4816&Xq4aI)Rn%tn1-##m,%1&S/H&###3?#V(S<ZN=m'g.%H29X-MGW:U-&AQ%pP$Vm#xK2[p6_x,pb#+S%oe'$&SOG#'##6Z?8$S###&##"},
{588.73,232.68,6.50,5.88,"Ve0###'J(V-')wX###(##i7$IxX######X##3+D######$##@n.###x>#T0/MwXR#%%##=^&FyRQC9###X##_SA5?'######;C4Du$.,#[5#SyXt$(h#&x>#c(1RcDw?)<8(OpHH[*vb#zH&FB3$J+.,#(##fwXH,$hH%a[%c%-Y8/>'+e/,,mHEU3@c$f~)"},
{354.10,245.47,6.51,2.53,",R$yw,mD2=#$iz(i,&)~#v5&U-'###U(#Oe/######3A#T,%P,$###LM%~I,Ca5###U(#x[,^]U###`0#9z4>u$###9M$cd+8l####1(#+m(]H(###U(#DW>|4L+##A%'jl8~#&_,#|}3bu$hK#~#&c-#OG#ll$ZP#:B#?H&`J(,u#h6%V$&d,#B#$`9+1,#"},
{354.10,245.47,6.51,4.69,"^ZS######^$#gZS?-(%##+:*wm+l.-EI,0%&PA2zG#Nn/U##xZS`5$.,#]$#AmO,1M-H&]Z#`6CimOI#%S5#bV?7Q$'$(w##%[SZQ$KQ'k##0+G&<*SiCZ##J[SmQ&bG$=%#J|F3##:$)#$#{z?T,#hI-C?#NEE8$#OWBsv$'[SLl%D>#8.#q96Cx,5?':##"},
{627.98,254.57,6.50,6.12,"xf.53:t5&K,#kE+dL:iY#aQ$N+0|k#95####,-&/,#7m%.u#DV9x-'YJ+M@'js7#O787)89*,XWJR(8#$pc#bM=###)##Wl$?SW/u#i[#792Yv)'A*TA%`UW=SW<u#%##@k05cK######O,#-SW###)##6e'll'###/##XM8q^8######[[&1a<######$##"},
{413.70,269.04,7.02,2.88,"&V=######|P#95####4##`w-######R##V%/######E##_5%73E###&##HI%PFI4,#:H#'SDw5&Pu#cz&[QP###F>#nA#}T8l&5###{y(Ib?qRPLH$z7*{6'3f.q>;C4:qm+OS#21:%/$Hp8(U81,#UTP%~(BRP(##PN=*Q$p1;,m&/j.r{AjI#nI.TU#VtL"},
{413.70,269.04,7.02,4.82,"Fw.######:##KQRG,#3/0<Q#EQRk,#m?P+e'+RR9w+95#LH#/'6######3##1RRVl%V$&{/*2B4Ah0n^,<$GCK4R/F/,#2J'HJ.######$##x&J]3C1,#pP#*/-QVR^5$zu${S2NC1###Ns.Mu$P#%###&##;Q#?y4######uo.YOE###;##J=;4RR###im#"},
{556.88,312.40,6.99,1.65,"<[(~oNd5$SG#,VYC05Vu%W,#CtGlG#2h:Dc#E~V/,#SG#h5#gR-Lc$WI$|d*lSYzY$,Q$U8&)+EZv*%o-UK+:jB6?'1,#i.$.v'ZP#:##=$':DP=L:'##26$j{)hRY-##'6%S|*R,P###T5#W>#_5%###fP#|u#c6*###q5$-@&>07###EZ$M0U.d)###q5#"},
{579.91,314.58,6.19,1.65,"VH$zZM######WD+Y_=0Z%###gS-`%)Zo43##y#N[,%PG#Z##QS,-{Xyb#1,#e{X>95Tu%F##+=D2Z#4EA|P#WxX2,#R5$e5#po5bc%*@$<7)CwXF>#3c$$J$h_?WG#Y93PJ)O+J95#VG#^S$hm+###K##x6(&{XuR.'##7$$LU+q-WQ5#N$'2U(.M>###QZ#"},
{585.02,331.74,5.99,1.58,"Vy'ey8pb#%##F.(Cp-*/1/##Z/T_R+nP$&$#:p7######;-#JtZyq>cG$P##jtFzQ$,4HV,#voZ'c#u#'=$#s`D######T##anZ:#$Eu$%]#2aFM,#(FCu7'ynZ&##BZ%px$2W>>u$###2##cpZwu'M>#W.$OU00^8g?&J],aV6L@-%##(9%AI(LI,###/##"},
{395.40,432.84,6.65,5.11,".,#yG#jD7<~/G6'G6%`00Xz;0&Vqc%PZ&Fi)Fe.s[$e^5^f,+n,R?%<n,0v&c]5;n#H%Vx5%5'VEf,?95zA$J]2b:,j3FNc#0A-^Q$$/.W~.ER+f/+k7Rx6)m'V.$Lc?'j[(8}D5|<8Q&h##E6(>,#E8KQg66,#N/.bkA7802I*>&0G5#;}<;<AC,$%##$_1"},
{489.51,490.41,6.18,4.95,"Vc#tQP#?$bI-vN?O^4###.n)wZKE>####Z-$D94N##[D<c-&Vu$Q?%`m'<QOv3B>H%/Z$rUPN(4Z7(D(7xW<%E+0,K#S*J6(Nf2@_&OiBa$'uRPow'[%/ix%am*-q%lQP{l&0$#JE@Y90#~-9f2#~'MV;W%)7RP2a>`?)<~$vQ(,'-)y4gZ&95####UQ#:y5"},
{489.51,490.41,6.18,6.12,"g0**~._G#7#$J[$^lJ4?$4//cG$+a.8N>}-)/^7-$#E@,-e&:](p+K.~.&##c7Q&2=>v$Gx+M?'Qv&Ij/kC9KH'<5#4I$l0/m^9JJ+tx5s?#F9QPkIrd)$I$u(';6QtC-Vm*y-$?o,9PApZ(q.0@##{5%Gp*K/0Ov(4<.FR&jz&Q'7gM*f$*If(`B/RU:Bc$"},
{529.99,507.73,5.97,3.23,"###c##_3@{k#Tn.oK+wM?76#2FC)y2bH(V.#`p'3D=>c$9d&######Nu<xY$n^9'##uSF>$&9[Qo-(BM/PO?bU'zB7UV'jZQ######Y~'.,#dc'###9{(0f1t~Q?6(eI#fj>%N/I1<Kp#uZQ###'##Ju$######3##Y%*T,%}#%&.'<x.JH'KB%g|AGl%6,#"},
{529.99,507.73,5.97,5.51,"60K4~+gq=*c#/~GxR'4%M)l#ZH'A?$9C4i~.`5%######uw)rvJ3Z$NAOa@*CgMk$@U..~c$R@O=L-VH(I,#R@-######ol${G%|N:50.2w*2g6?P6+m(`%&r>OXc$ZP#/&#r@/######a5####3R'i;-VJ1,%-J6'a>##u7p;B######y8#QH'ZP####3##"},
{708.75,508.63,6.59,1.62,")Q$BJ-###%##YA/?v(-m%0d&mP$Tm#dKTQQ&.##e0&HsE95#E&-J7>pm,)##tCY5f0*e,f?$<w-~e#aAYNZ$###u##k[M^Q(UB6m[$ksB/L.xAYpo*..,Ho$Cm(DE(FAY1,####<##lU5Qu%mP$###TJ&sYGOm*JQ#sS0:T*/,#36#lE6M-)###$##`I%)i@"},
{708.75,508.63,6.59,4.98,"Jv#P3BOG#$##?W5LH',c$s5#Sv(a{0)1:X>#xz/)tB~#&###q91u$,###/##ZfXN5$U@*X('j[,Hv#9iXuw)+mPI?%A93E]'UW8z70###%##ZhXDl$5[)IZ#:B1NH%{iX`e.|Z)###H~(xG:?<@ll')##9m$EhX>Q&D>#0-#tA*;%+<%,+Z$1,#=5#Gl$1v&"},
{505.67,517.79,6.49,5.17,"[P#6c#'13of3$d&c-&t..VHJYnMFm&-@*&k4#I'(&*9ZL8n*$^4pl%>Q&8[%Qy5;p%0RTf$'yTTp/+:q;xZ#Nn-N/(TRSCH&C'-3%,Vc%,@*rd*B;9s6NF%+RST-vE4@)`J*8)<5<97?'dZ#^?(h,%5Y9Cy5C,#NM;Oh4ev*e6)dC6G5#@q+Ag8nG$###+9#"},
{747.92,550.43,6.65,0.51,"_P#Mo)3~-O-(PR$S03_5%###):2SH&`/4###r5%>Q#V;AU,$-##AEV)$(D>#gy,r@V{k#%##_CV18/I^9X5#Lv%<v$HXH&##'##f&-Yp9###lRTPK5i,&O##|@VhG$:E:9C({#%GH%OBVVu$###Av#hR,OG#$=B(##Au$_G#REV^x01B2??$f$&7'.{@VmY#"},
{673.77,572.06,6.37,2.71,"p%/xb#D>#U,#XJ~95####Sf#J3Bz.-####q.###f.)8u#H#Dlh?###I##W'3zJ~E>#%##_O,K/1y$%7l$-3R###&K$pc'S)8Rn/###T$#r`?1N~pP$8##oT)p=A=g*p#'~?%###~R$8d)###0Z%###*H#%_1puQ###1,#K`%_&3t#$J-)/-####^Q%#?&1,#"},
{750.18,579.17,5.99,0.26,"LI'*04j>%3##w9~.7,&l#?]#Uz<eR#[NDg$$od'F-$hv+###D(<I5#&Q%&-#,9~#c#s[*YN'i~0*n$H:~&](2~'F6%We0####z7$##1,#xG$>>~V03P?(p5#1g+Qm8@8~aG#]94&-$S%/I,#xG$0,#C#$gY#4E*4@+###/,#xS&FL3tQ)7,#+x/A5#*c$J,#"},
{750.18,579.17,5.99,3.52,"L#%J5#Lw+|b#J?()##*B%Vz7###$##Wz&O[+6l#:#$Q,$$##&f2B,##g3Zu#gJ~[G#(L*H[<<-(h>#VP~OB3F>#S,$GU6&##<~/###3~&OQ%dM~O%(HJ/W$%[n+wi'dK~+l#Hl%_##'`@6u#H?(###%e%g,$m>Ogu#yA3Kd#Rl%1e#2L~x$,6#$.##-S)=D>"},
{763.97,589.30,6.09,3.64,"###O>#cS+###{Y$Su#2J*z,&J#%?,#>y%1M;/,#;5#n.$vc(###&##[&0###_(=gG#y^7AH$onZ=5#HD)BUMY,%5,#V3Ub&3######gS'uG%_g9###('*jH'`sZ[['/y277'tR)ko&{pZ(Z$###$##nw#@H'cc'###ee#yu&NoZt5$:e-J?$>m(g-$woZU6)"},
{683.90,605.09,6.95,2.84,"tD$+V=ZP####bI;ed,###3,#N&/7,#F,$uY#.,#c?'P5$H>#if0W,%n>%^>#r_~4l$###+$#x=I+9+OG#T#$###*`0=Z$%y2lC;K,${k#[5#T]~7l$###do$hT33&+[P#R[9###dd#n.-_%Gay8/,#4##T]1c^~gY#'##D2()91f7%H-)Lt1###3.#|@.l#%"},
{790.15,677.50,6.30,3.75,"eW:95##['F>#zsZ#########K;6N$*###.##eG#*B4###9,#o(<~5#695e,#'qZrP$###~##74AaOD###$##-##T~O.,####~YG=v%5?'&$#qnZA-&.?&=T#7&.ImG<@*fe/ac#o;9EQ%7d)-kI6u#.,#A)$NPN&##&m%X2'L..8,#PA*wPBTG#YG#G)2~82"},
{455.59,742.95,6.63,3.68,"4e-q>%H@(17+>S[L5$m[*(0'RR,Q5#1V[kd((-'^G#88+nd*(95zc#I@,T>#8T[e#%R//:g%2o1c?%(X[#%'xY$$##Ko/iH&U/1|6(8$(&##EY[pkElc'^5#O8)kX[i-StY#Au#^l%*~,(Q$6$&,e*,$%p?)X_30A,Yl%'o*5@'''2,$&qW?)-&$Z$0l#i6("},
{825.03,744.76,6.74,5.68,".,#$2%cf5YG#%x0_R$-.,Z-$uU;{b#E>#>{Ppb#$##XH#Y/YMm)A6#KiARu#6#FNu#JXE:$%x2Y;5#;5#4y([&3###5##y96~7)v>$-_</,#5{=-l#|/1z_-//Y###1,#*`%b98######S##J[&I,$.;;###KcCU,%###f5#'3Y######y##Xe-OG####3##"},
{803.77,753.43,5.97,2.53,"$##X8$pSW###8u#lH%5:5eG$Zh=###VR'Q,$<SW######A,#0Z$w2(gSW0,#S3A6A*y6)P#$L2A)##/x,(6$JRW###<5#R$#?u$x@%xUWL>#<}GEu#gJ/$Q#sNE*6${Q)2.#BRW###.,#G'#{k#XE1VA2wP#O[+S5#s#',8%6C:?5####jV%^rC######C:#"},
{640.79,758.78,5.84,2.22,"95#1##O49OH'?d(bu$.h1hp0wI-OG#K5#UG<nP$.,#F,#5TRnQ'.,#935;Q%R-&uJ/)AI58,MUR%/0U5$/N.guF=cL4Z#gURTy*O$)-A+'l#M6'g-'JVRUR(dSR'K2c?&&T'ao+QSRrb#$##mg'|7/pb####b#$B%*yw.uP$_,%Ru$+v%[#%###vu%,v%###"},
{640.79,758.78,5.84,5.00,"8l$###Z$#-5M--'{Y#I?%dS0O[)Hl$=#$_u$B,$0##V//|u&SZ#Yy0g@@,dTk,%4N)5gTrS2UeTP&,=?&d7'iw/4,#L)3Rw*9A#ZhT6y5N5$(?%Id@lJSv-*LeTvH(X?%rS*@C8#Z#M93XG#zd&1C8}>&###6m&|*@En-o>$rfT#w+TG#Pc#g/0Yv'Vu%$##"},
{464.06,759.27,6.01,3.77,"hl#9-'%6$zT7hi@QS0_u#Pm'$SY/u#W/-g:-Z>$V>#GxPDv(###i,${c&J..,-RD##dA/~#%>SY###z03tw%6?'$##*XYM?&###Iv$XZ'###WPKw,#ad++##(XYiB3q.-7R$5['%|5GUYL#$+##_@,-c$###28,g@*l#&###,$;k);5Z%sb#&[$Z-G#R)sH("},
{667.13,767.08,6.32,4.87,"'D$MxTne+rI.e>#knB$yTA$)#wT$K+5/-^&)Lq<U>#X:557'-:(rmTI#%$##ZR&0yTs/1Lu$VxThf47-%-J'%{=7-%l~0K,#tZIm$)MQ'd5#6@*aR(i95KK+owTw#'K>#xB&4K/CA.+u#&##vJ31,#H>#.o&+R*1,#iY#O0%Q]4.,####sI#%[&I#%###$##"},
{438.93,768.87,6.64,0.33,"1Z$SQ%OG####xRVNu$1,#Ce'JNBHB%)TV1;7W,%6$#_YL)H$O~'bd(m,&dG#KSVEm&~l&1I#}_>D{&/RV7##vm+,R#w0:vu#b)(,6'QG#T>#WLLjT14l$,?#fv'RF/TFG*@*hR'Q7+9A+}u'}9$pb#######|G6'd(.,#/##S%*/`<PH%vI+Iu#F'4@f*:[)"},
{438.93,768.87,6.64,3.92,"R#%f5$'6%D*4pb#I5#rz'p`A######x9#Q?T######@%#kNFKx,Y2>yP$|R)--Q_-(rI'8E5JZ&~G#w/G;D<NH'bP#L%'CV=5~-b#&]m&}7/*BTOG#w$(f?%6z:%##sDTf$(CH'$##|&-M$'Eu$v5#_m*eY#)CT'H%*6&2c#RYCeh3t/O/6%3,#1d&__74u#"},
{843.57,776.09,6.11,5.25,">-(?J#tEG###H~.X#$DQ&zG$B%.^5#L@*^6>lYO###$##.a*Pn-}H#4%X###5=EcG$###f>#G.Y######.a*tn1######7;&p$)y#$1/Y$##IlFgH)###$##H0Y.,####Y##BM=######f##6,#y,$>bI###r:7HH'%l####:3Y#########Vr8OG#######"},
{472.18,781.29,6.50,5.50,"+6E)1:.l#sp*5o)G(~.,#2a8q-&=+~(d'Cf,/,#|;QR:6Ql%[B+JM92PCPv(|S0l'~$x0V.'V,%CGH5-#i19.,#G.+Hm%4A.^)4v'-Lh9s?(Q8Y.v#tw-`d'8x1###Bv$uZ'######cn*I?(b|=/m',8(ae,_19Vv'Tx*3w*ue-pP$cP#%##%###c#NI+.,#"},
{451.49,785.63,6.80,3.91,"ZP#`.(W,$<f/5l$Cc#`H&[6=D>#.##R)(}}H######3'#csHH[&=8*{7,d^5r9.K=E(H%[C/i|F;n*ie(1X8xG%A5#nu8Y<C~>#g6&A$',i?Bx/d[,um&jp8YxT6l$`R(yd&Ip72,#`{Tr-(.,#Q#$PG#vl'8Z%TQ##J-V,%EyT8c$rZ'Gl#W}Dh/-G&LZQ%"},
{668.91,790.11,6.31,4.58,"&d#kQ'xG#<[*5r&Kx2fP#{Y$8)&x1V<T0dG$56&r#8G4C)6%6?%c6<hH)XG#FdD%N=}k#H##{p*?-K;[*C5#Ui9l{?{%,U%'2$(.A%Nj?*U0=QO~Z%gG$rU(O$)8%'tS09p0SuMNH'qG$q),.,#+##xm&2'3(c$###%##>q.A,$$##1,#qM.k[,.,####B/&"},
{814.35,792.36,6.44,2.30,"######hg1ZZ'SkF###)-%o5%msXZ82/,#}h+Vy/[sXN?&bGE###1##cGHxY$JD>###QJ,o>$msX*y5H>#J,#On'oCXtY#tb#9Q%^e(.94`P#|o6###|@,ku#KnX.,#<5#jI#{$,G>#t>#6['~82fl%su&*H#4'7+##w-*1R#vmX###$##b/#KQ'###Z,$Yu#"},
{652.99,814.61,6.47,2.62,"qZ#+i;r5&###Uo(_[,###'##q>$)?&:5#?5####a5$(m'###V.(an++X8gH('MWx5&/Z#GI'SI,sZ$mu&wG$>5#sm)dG$Y%,9d&F,$T/%Vp7Z'8`P#v,#76@^Z'z#$*Q%{NW/,#?l#5Z${NW*00j#&)6%|e'9s@yG%+##I@=~[*1_/k,%{NW###Se@Zl%=?E"},
{420.46,819.52,6.36,5.86,"95##########`|G######M##umUD>####aw#E,Da*C0Q$B@&ZP#######,##sQS######0.$qrUFw.###B]%PV)&(8Y'03/.O$*######?##_-T###'##c+2va:#.,~$#MYBeN+]96Pm'}5&M$)######$##SrU_5%&##I?%e3/~./c@#pn/hS++u#?&$)I*"},
{832.44,845.52,5.98,3.65,"###*##+R).,#J;@9##0@+s##47W###<R'5`.bG$$##4wJP[*###16&%d(###3?QB5#%n+H5#qwY###5&/um$SH((##MyYy>$I##W7-OG####7=DwH(|Y$+##s|Yf06xm*_6$+@($|5@yYM#$<##.p3######-R%k|D######+N*rwY95####&-#uyYac'.,#"},
{657.43,892.21,6.21,5.26,"#d'#########EO_#########9N_######C##)$(###/##RQ%v>%#########aM_######C##<O_######V##c82SG#P>#9u#<c%######'##>K_######0&#zL_######E&#.g4'l#95#&##o>%#########6h~######G##&P_######9##PD;###G>####"},
{966.55,905.45,5.99,5.19,".6$?p-Q03<Q&Y18]-(Wu%d,#jo.n(7rd,D,#pq7NQ'###(##:p0oz1`h1&8-SpY0?&3B++R&>;=V-&t.Af01-HNG,$V>#Te'i6S6S'R+>SB+jnY?-$K{856%.NB:,#_z(Ay3-H&~@#iQ&%&,|T57E/HB0#@&gnY7e&C-(}@$Tg8G0&dZ'2#####R{%Ll%###"},
{31.07,68.87,7.98,5.89,"############5>M######~##t`[xY$(##m^+'50M/3q&#<lL.,##########E7W######+##x_[h,&###C###4+x;C(##,l#.,##########d?S######A##H_[$~-.,#-$#.;2=A0j#%U6&############xg8#########+BN}d-0,#/,#)a,S.._>$SG#"},
{143.18,70.92,7.66,4.82,"QM>(Q$K%+by&o@.fd'Q3@d&)rZT}P$_P#{w#u?V###ZkLr##ks;*J,%?&OH#H(8be,N@,LR$E[T###$##1I#2@V###o6U0##BH&SG#fP#9W0e]0,I)y?)}v(8KT06'D>#;##,@V###,$SE########]Z&oS+#$'Cc$|z9b7*jHSF,$2u#%/%t?V###l|Gg##"},
{596.88,125.98,8.19,1.87,"~P#-##uoX###?5#L[%b]5.,#u6*hH&A-(zc$T,%9,#%l#6H$#[(,##hqX###U-'@z)_7WpY#GnO3S*{l'5v$*?&2$#xv+O,$Q^6%##cqX###n,&q>#'rXt'6X{@,H#xv)}i4p>%SC#q~25,#abG###SpX######3##s|@:I+###'&#0955@*0,#82$DT5###"},
{596.88,125.98,8.19,5.09,"*WB###$##pL#/'6OR*###u$#@ZEId*###+##UyS###]RL###zB97,#Z,%q0#C7*KE4{V?yG#XzSR04o,&f>#C{S###l04$##4~.(l#Fc%=$#7$(cZ$jRLu[)kvSdP#i6'J^)uzS###Kc%+###Z$=l#mP$5##1?&eQ$fd+9?%ko4<5#WG#Xm%gwQ######.##"},
{494.05,132.94,8.24,3.47,"kn,7K2###I,#%r2ItJ9#$(##[/X<Q&E>#N##Zf4-m'OG#+##3,#s.,m[$%w,R#Hdv*&##/K/V0X95####y6$Q810m&Av&K?'######JV#6.X-@+###|Z$h1X6/X###/,#Ha.%f0###-##(%)######7&%WL;,u####6B'5AWJS.###3,#)2XYv).,####MnL"},
{638.16,144.30,7.73,5.96,"5).?16TL0Rn-6$=an0&##*6%5z$yiE######%1#}bO1##W>$Rv(]5$~h)4jCb)8[6)1w(qw,0*8'$'pl%Xd&gc##=<)]/)6&K?%^l%[-(95#Lm(pc%rL4[5$,AFqdOPA.Cn'@I%/'K]jA&d'%w*###9#$%##=5#ru#0cO###v##9wKvkM.,#+I%/~Jdu%nP$"},
{442.64,234.47,8.22,3.20,"_)<######lP#BnR######mc$+cIG>#L>#*9,>.*Hm%?o-h'7Y(=######Ie'|[W###g,#qF38s<;5#Ch(m;<&+<'?$W8+Ho2F`B######Mc#m^Wj>%qG#E6%CO@'o-k/G`Z&oxO,Z#.B0f#%Te.4u#+m'qQ%R`Wc6*qb#V##B'*=vL4R*&##go(qn,fc'###"},
{338.56,256.38,7.83,4.82,"7x2######8##LwZ######Z$#4HR###MQ'-%#K~0$##M82T,#1(;.,####-##mwZIc%###:%#m<Ch`@Oc&p.#M95#%*T%/<##[z=######-##@yZW9-g,&7##jSMBNW^Q()##7yZBm)uG%4##~f5######3##@[SyQ$$&11##ywZfm$KB6g5#;xZqZ'W>$2##"},
{413.36,268.45,7.29,2.91,"crC######{P#eY####2##6f0######]##]83######S##Yl&+#O###$##:e%n3G'##'H#z@D6Z%2Z#j1'p#O###I>#zA#<h=/^8###a8&(nD0%O1Q#i@($S)<8-=E0Nj8[I,d%$aK5Y7$P97ZV?cP#muB+K+7$O)##1X>Kl$%`?,Q$tV*g=Gp6$+$(HU#2#O"},
{413.36,268.45,7.29,4.80,"281######A##,6RU5#5%,gc#WlQ;?#-|@d@'*7RI.-###66#J:;######<##u6R=Z%4$%`K*,94my/U])K5C9y4c]N$##aR&J]1######$##j&H6OE;5#oP#|A/y:RoG$_Z$,T2*12###/=0Su$:$(###%##BQ#/3@######<_0DvQ###=##Hs6O6R###Hd#"},
{453.32,275.17,7.75,2.99,"`U9#########u8Z######$##<T4###5,#Z5$######A5#<Q&uU:######5,#A9Z######'##f*G###$##vd&%Q%###]##k,Ix@/######+],,8Z###$##S1(u7Z###J##{D1XISDZ#:7'{+?=A1######sn(Z7Z###3##}j-uDA###P=0*mC6%O%H#/],UZ&"},
{628.48,581.67,7.72,5.87,"->>iq;{k#4-#P<4$NT*e-e.(3?&.<*G/2Yl$bQ(u5#Bv(Cv%U~-#N3`EF4##6KTC?L4[*4##jV?7Q%###+@$9Q&b##SM;67&kZ(V##wHTX5#wHTF5#WJ136#%uMJ[)###d%#3H&~33M:8iI%bG$'###>AW5$icIM$*N7'{6's]0HvL###Sd%:5#046Lg7mc&"},
{705.72,615.52,8.66,2.55,"kv%''0^H(kY#GG4N&3OG#3##>z1ZP####]G#tG$95####:l#=##N8)|&Ucf3Fc7_w-'@+5##l)U######)H#G$)aP#}k#u>$3u#j*.#%U1-'|,PnG#CZ&,&#H%U,-'###d9$O-)qB3)c#qQ(6l$jR'uP#$y2E$U/,#-##y)*1i@-?%<5#@JB0,#Xd&zR*.|A"},
{837.77,647.41,8.58,4.41,"&^*(H%{]0ZP#u8Hr,%uP$Rc$k96T>#9$$FC03Q%###B9##]2fn-Kl$x^~-##t^~$l#rP$s6$t>J###EZ#XS+]r;###=/#P95+n-9##m^~'Z#o]~###4?&_Q#1#MPG#}Y#{f,l_;<?'S>#Re(=?&s#$:]~A5#q^~'l#Jl%B,#s[SYu%###)$$ac%$v'###,-$"},
{490.60,724.67,7.80,4.81,">`)<lNa7&:aEmv$7&F-VSNL:Gi?g.(SB5<%%ju&2##>)7xH'+(#pcRDw(mP$bS'%VS$_;>u#CUSF&0/[)7v#f.,tb#|d({#'###:6#.^4sb#B-(T['v07>:,WRS9l$I>#`(&hJ-mP$<5#6,#D>#4$#Om*F?$n?*-##D,$+z%b./######xI#a$)######(##"},
{824.33,744.95,7.05,5.75,"###n_$b98###Oe/,A$wQ)>m#^y8fP#<5#D:K{k####t5#6fXB#$Z-#^FI$l#y=E=u#xg8DH$ThX/,#.,#=B(~T5###,##H18`l#<Z$qjH.,#.r=VG#2S,6')heX###$##c($/V=######P##,H#2u#1#G.,#Ls8W>$gY#K,#ajX95####R##Pp0pb####*##"},
{452.04,784.88,7.12,3.89,"OG#TR&1Q$rJ/O5$&Q#5d&AdB###0,#E_#_$S######n'#-tJ}Q&y$'@e+(p2d^/2h;sY#hB-n<E^$*ZJ%tO=(c$###,Q6MjEvY#LZ$=v&rC<6y3|c(v6&.U7WSVB,$B.(o@'Ox3###jWV~v'/,#h,$]P#x#'Z?(,-#&e,xY$-VVCH&/6&Al#v_9`p/.SLRl$"},
{630.79,785.08,7.54,0.59,"/?%Ud)X5$Ol$/?%V4Atc%N]2wG$B_3$9-X7,bZ$LA/,6%~Z&/e,xc&dP#)Q$:vRgP#/A(1<Z1H&3##Mc^$kE.,#k7&:D946'#(8tQ'###$##e^^=Z#^o3^?%7/0PR#E^^r5%fP#T5#&s<iv+1B1T?(######k`^bv'm[-%##+%*TK'w~^###PQ&P5#Ze0^P#"},
{630.79,785.08,7.54,3.63,"L$)gG$6?&C,#j%~&##7?&ty%yS3J##i(~HJ)###)##VD;FZ&=3;zH*UG#LZ#w&~l,%wd,CT#rT6P7%B&~|l####2##cM?zc'`L756'G>#Nw&|+~DYHj,&B##FS(CDV3RU[G#>5#*c#r7.4m'bQ&(6&fQ$uo2*&,g7-6Z$WN;J$&e/2KQ%:a@yb#Z5$'Q$]H'"},
{675.50,810.85,8.20,2.22,"###P##~L;###g?'#S+w6+C,#aR'm-+.,#1##c5$D>#-H%^P####/$#5,N###M7'z|7nFGyb#sfJz/35Q%Ul#zv(1c$8#$9,#Uu$d[$#RRcP#P.,9H$l>:/V7%,L6l#2R%vH@36'-7)rY#;f2.6&l,$'~Bu%,=N?(Z$cc$5*RBR*'7*$6#-F^/,#gR*y.&gA^"},
{821.95,828.47,7.20,0.47,"###jG$9.#hB8###QQ&6?#:tJ###l5$pl#kiD###UG#G5#p96pZ(:,#j$$Qz:(3C&##~R%yO~n>%-##>UG4K~######Pw(-aDC0516%QG#,##hL~y,$1~-5w&<A.I%$>L~wl&;Q&/##LL9ZH'-o0W5$###$##FM~N@*=6($##s$)JL*<J~$##~d+<##fS2:,#"},
{821.95,828.47,7.20,2.11,"}^7###0,#7,#i;?######N%)^l%26'&##MV=,##+I*D##{K8[|F###sP#lu#z8[zG%)##-pEVQ&8f/o?$9:[###7J+iy+'8[.}H###/-&n>#|=[*NA%##'7%xd&P;[q5$#n+###M;[7d(~#&&g7###Al$L-#-aC.$()##<d#O>#;o0[c$ub####j[+&Q$###"},
{821.95,828.47,7.20,3.49,"c6*nP#`m+8##6eZ&##MQ&Dz'.7,D##igZ=o,###/##;19hY#6B3y6*>u$=##/gZmc&}m+)y#xw/@~&<fZP6$.,#6##*V<Yc%KS(.+I###$##g9G#fZK#%7##Xn%RjZpiD'##,I$te/<m)zY#sP#9/1###4c${Z#ez>###RZ%5@#TaH###NZ&b8#=]4###sb#"},
{433.16,839.93,7.46,5.20,".,#######*##OG#######ZA+7.-######Lt0eT695#/,#g`){H*######N##>cP###Q5#Xl6<=F###i%'U1OV^Wyc%C-(m8&'//######)##w_W###/l#q?$`>I$##}W:d$'.lF|v'~R(dI&AZ%#########T^W######;,#)]WCc$n17#v$$r6;I*q7+i?%"},
{473.51,143.08,7.83,3.43,"%23d^:?#$}P#$gZ]..Nu$x5#nS[{k####j##=m(Gl%###%##xy,wrC$##0,#jV[c6*###I##jS[$v&.u#M##$R&yy7_>$TG#JJ/r>%b,#%%P>T[######Vi,^S[%c#Ql$u$&Pc%YQ&iQ%Tf195####5A#<T[}`E###2##9m:(FF###$##,8'eY####&##Sx0"},
{757.44,188.47,8.65,4.54,"6(/@>Gdc'MQ#;,#'a/`QSVG#pb#e##n=HC:1(n+K$)KQ%3_+EJ#15K&7(t5&Fn)a,D6-(6H%FbHou$3n*'i/,-&[c&I^.%j;I5#lS/r%'5w-j-)`h=B,#5B1Z5EpPO^P#G6$.7'/JXJ#$zv)###FKXjY#.,#D5#(KX$##F>#&Q#[JX###$##*##nJX######"},
{297.75,195.95,9.24,1.62,"Nw.;5####*##Uyd$##95#M$#A5OC##rK6A?%mP$_##vL4Z;?yB8N5$###,##iydK>#ZP#L$#YQRE6$r@/L,#?$)Q.#t3Hxb#*r?U,%###1##qzd<?'###e##UGFrS/VZ'%##Y$*0/#y$X###TD@95####2##'zdL5$###V##<,KLZ%iZ(0##4v(W8#oe]$##"},
{465.56,198.81,10.52,3.02,"z@/QG####Y##^(>######Uo&Q[W###--${O0&OC;,#nS)b01?q;.,#$##;n#dvUG>#b,%#K&d_Wy[,KI)iS%mbE+g0PGEfc%oD3{d-F>#xc#MRLl?)9H&;T$;_WTe/ZP#x%#A'*o;=Qu%%##)y2fu&#d#nf/,~W}k#.,#{($jJ0Q?(QG#W$#o#%bG$###*##"},
{415.10,225.86,8.38,2.89,"Q;A###.##5`)=_=###~Q#iV2~z;+$$X7)i8.Y5#A['xn'YI-^ZR;5#~l$-i()(:^G#HREHT.EQP0Z#No+m8/kl$7$'ML&0}HX^RtR-Z>$'$#wp2m|;q)B/##h~R3Z$tc'te)|%-bG$g%#5QMIx-WZ'####$#)w%OK4.,####>Q9|d-###1,#mq&yiE5##@c%"},
{604.23,275.75,8.38,6.16,"2M8D>#$##7?$8Q@#8/###vY#&j,ho4<5#/l#D.@mP$0,#3,#PA*e$(j5%:5#5LSP#%-##^m'JL7s,&XZ#o#F9JS7#$%##xT)F$(/H$[''|HSpIStb#f5#2v>2bK###)##vF7TM<######4A([P#>5#hA%(IS[S0/,#h,#`LQv7/###/##~x'Om)###2##AZ%"},
{492.95,329.95,11.30,2.17,"`>;6jC###eG#_)0+]Z###&##]QDb[,###8,#V6)$##},&g#$SJ)wiB###K,#u&LYD>####$#4bZ6#$###k##aS-xl'$##:,#ZZ'pb#%##c%&[:<######,)#h~Z<5#95#72#'d'??&e5%>,#>u$###Q##v[?hv+###,##i*+>u$&##X#$fo%0,#[G#uP$0,#"},
{504.10,357.15,9.28,1.88,"k)4#vS###m>#EE/Y+L###'##L91TH($##ZG#tb#V5$@u$&Z#0;/Ke[###;##]k[)]3###Z##M@VOG####e?#_5%*##-v'xc$QS0D>####,y'>f[######8_$$f[cG$###f##`S.I#%=5#Uc#eY#######+]A3'7######rb1{R-ic&,u#'-#UI&5m(OG#%##"},
{657.09,534.13,9.48,3.18,"0oQ0,#OG#pl#PpQ)?%dc'.e#uY#Wu#4r:Fu#e?)$##OL.rW:G{<###*##]X5JqQ>l$,u#w6$A'59H$tR-rG#wlQ1I&w?(h<2Wu%;7+YK%FnQ~+IJl%}>#v4<nB5V[)q,&9p&5`@gf(jZ(p@#s],rmQ>~(=&0qz8:Q&g5#8s:'?&TQ$kQ(l44pb#E$#<]3KH$"},
{803.09,564.54,9.44,6.00,"2)Sg~0###/f${U:2I@;5#4H#S[%Y&S>5#(l#R.'v09###$##c;?[j/p97??:Y{>0?;]#%a-'t/1S6'S?#ox1x:;ZP#&##cH$Ci:h~D2%-3R$#%SM-&]P#]6$V:<2,#R$$D<<#7+s>%b,$EA*y09`u%(##tS$4$SL6(###H'#PC7ez:&-&jm&+Q#a7*qc'9#$"},
{605.98,575.98,9.37,4.40,".,#@-'N@+.,#u-$#f0$o+H~.pA(@q;B5#SC1uY#+A*O~*9]2###e.OTu$`Q(V9.wcR3Z#.'3(fR%m(XQ#m.C`Q'iY#FN+#>GT$*B@(yu%TeR<n,Qm*@##peRihR7&2t,$k)66{)_WA-RA1[(+s96H&p,$|.-Lm&m?(FI$m98&m#yd,N6$u?*rC)mI-{Z&+u#"},
{825.54,596.68,8.76,6.22,">|<&),Cd*Yu$318)##z5#Zi>%2>.,#6##$q2P..{>%Sl%0H#?^T6~&SH(I##8]Ts>%W>#jp1m$+>e)v&0x2=su$192t#'O>#4[ToY####c&#<~T[tB:5#[7#V/(Z]Tf[,WG#C-#)].Tm*###`HKeR,*##{$$To.zRL###Y##:](v~T###&##i,#;`T8#$###"},
{868.04,631.39,9.75,4.32,"*7+N>#gf~@,#3f~<5#4H&~[#&W@###D##g7(Ln.###sI#:p5Dn.F6#,f~R5#ue~2,#5H&&R#R`BA,$&##CT*)01sl&|5%>d&?A0.H#8f~9,#Qf~>5#(Q%q##d<D(c$###-I$Jm({Y$=5#46$AI*;L36OF2Z#1h~IU:S5$h?#4z0~e~###<Q#B6$nvW###Q>#"},
{712.41,639.55,8.56,2.01,"o@(e17lQ'sb#LIMw5&###~$#p>P######G$#QH'###Z>#1u#Ug.1k/C3B,u#<l]z#'###V##mHS######Em#m-)######}b#B/1U](Fl%S0&sf]:5####/N$m5Q###$##/x'b5%?5#<5#A6'yY$.,#'##GTBX96~#&<,#zj+`[)/n-<5#0d#:##g>$e5%ZP#"},
{449.97,658.56,9.09,2.15,"pc%iH'95####aZ$]i8`P#B,$QS#,2<-)1^#&=mImZ'{=>98*>90?R'37,$##an(0fC-AWG#$oEW4@L)e,C6$wS.'v&2A,=?&GK4A5#Be,gH#[-)R6$oEWZ}E/AWc5$b/*Q1M=.-@##W6)}>#:n(~#&8,#f5#%?$iY#I7$iJ295#7,#g%+FL6###T##a%/J>#"},
{449.97,658.56,9.09,5.33,"g%/K>####I##'J*81795#B5##w#iJ2*?$hY#7,#n5#Bn)g,&)[)}>#S@-;##h/*c1N)8Wb5$c<W?kEg6)C-$+S,O?#+g5L>#Iw*/6&{n-^Z&+e,?6$c<Wo-L'8WQ,$T~(ExCTR,###x8.Xd'Ys:Re)+.JBd)]22gu&;J#MM<nY#Y>$76$43995####QH%=$'"},
{871.26,692.69,8.83,1.21,";-(S,#w8]3,#Q9]3,#@?'uG#QsD^P#O>#f-(###%##yG#(S.[H(P,#}8]C5#f8]9,#wZ(}?#O-T######6@$F>#.,#zY#XZ&6[)Y6#G8]$##[:]M5#O$*5##B8X+u####.##WG#*-'[P####2Z%Y,#K&U(l#w:]gY#?v&V,$`5C]h;######:##slN######"},
{627.61,699.13,9.31,5.21,"Om(N,$###$##(v#Qp3.,####2x$lJ2######Dl#V?$m1<###-I'xl$W$*.,#s%)xAEO%VV,$.*VEXA]Q'+$$L04/,#6x*`>$DL5Tl$+B.cI*pI-(v#.*VG-JH%V?u#8e'+]D@o2E>#,##qQ$B'VT.(tEBil%0{90Z$VV)CU7OG#-l#kZ$6a=:5#/l#G?'(?%"},
{403.70,742.98,9.20,5.79,"d%0######&##5]_######'##raI######&##@l$;5#*H%$##U&4######'##N]_######)##&HL3l$###=,#Oc$s>%-H$X6'vA4######&##h]_######^##t+H`5%###1|%#m'U.)$L,rJAan0######6##B]_95####,%#[QMAN>eY#5V$OQ&G2R<$(M.&"},
{611.14,813.00,8.93,0.29,"T,$/~.######EQQz5&D5#['*#'6Z?#yoPGHFA5#Q[%ElIK#%<J,*n-.,#+##AeU;v%Fl%V/#SrAXV%HdU'd#~P#<w#2?JYl&N~+U,%qP#x5&eiUBi:4l$pP#aR(]Q9V{B(##]c&EZ#6v'>l$>c$eC<,##4Z%1-$n5M1,#|k#P5#1nRD>####ku#v3H###%##"},
{611.14,813.00,8.93,3.94,"95#######oZ8bG$:5#;##]+5:5#VQ&,m#h7,?l#p,&S>#DZ%n-+###J#$io-S@X###,w$~[=ml'0,#G7=wh<1,#95#,I&`e.o@-Fm)Xu$SG#3CX|k#;@(*v$p)AiP#1FX:H%######0(4'7(|%+Dy4(Z$Tc$1FX-w,&l#4##D{3X6DKQLA5#`P#%m%iA0GH%"},
{808.33,856.64,8.27,0.20,"7Z%Xd)~,$^P#RZR+l#ic#Zv@,/16-#uCM>gX###7##it6`:=uA2_6)###:##|eX'I&<c%WJ#pM?Di&LeXUQ$(c$S$#cEBs-*'g/pb#'##=#$3hM.4AA,$7##6R'>[:C{A6##%.,a##:$)du#>l#q>%l>#i,&OH#uJ1=5#95#:Z$Qx.OG#:##SH&cc'###;##"},
{808.33,856.64,8.27,3.88,"bG$###VG#cf*5l$D>#2##9T-$##s>%sQ#6A0$##-?&uG#K-)x$,+##t-*j#%&%W%##yQ$S$;|,'###^d<y;>######On)pH'9w,xY#jR+%##`(WRG#E$&Km%TjC=H$K*Wd?'H,$e5$zN<RZ$M,$1]/dG$###ILKhWCSG#0,#6<2K*WIN?H>#d>#p[JGx1$##"},
{819.24,936.64,8.67,4.82,">u$#########Cyb######~##syb.,####-$#Vf3eP#D>#/H#{Y$#########g{b######2##w|b######($#~V;95#G,#zl$eY#######$##Qyb######b$#Kyb######M1#*J/.,#P##rR(bG$######%##.yb######*%#_yb%?&###y%#wR*Pc&###Bu#"},
{232.68,220.21,10.25,5.54,"Oc&###$##U#$cL<###%##kq(Uv*###TL#7k=######p1$[l&:x2######U##@[U###&##mh%ZA2###+U#-E4######%i&KQ'm07######*##Y`U###%##KZ#MRB###Y.#y-*(l####cU#|'9fd).,####%##$KE#########n(P###,##qb#vn*###Q##0v'"},
{439.19,219.93,11.44,3.05,"u/4######ZH#1%V###%##;]%#jD(##-H$l{3Vl%1c#)x(@q;pDB######3S#/%ViY#DQ$N*,e5Ma,$5z0@h2T::{Y#fw'HdP[W@L5$1,#kR$J)V3e-lG$fH#h><a4Byn/`G#D&D41;U>#:o/=,GI#%###D@#(XC<c%###J%#97&'w,###)##WJ#Gq=######"},
{753.35,234.22,10.10,3.17,"+@Fd5%######W>#U&20,####=c#=.-ZP#=5#=5#F>#0H%sb#VDS]R*95####7~)NAK0Z%'##QASAy66,#A0+pQ(4u#XA&4`<VDS'6%mP$###+`;Ze'q~1d5#5DSQf4TG#o$%0]*&BSRe)Rd)fCS###:5#4,#'XA/,#wl%HU4Px*.p5v@(L'5I##.BS87(7#$"},
{467.37,254.00,9.37,3.35,"xS*qb#######O:K#########EGC######gl#@J0###ZG#O<8+J+:5####.##_(X######Y##opY###$##<&&G-L<5#:Z$6'-uf5/,####C##NnY######x$#YpY95#0##o9)h`>?c$T'+4q4Yy6OG####,##',L######6?#/tY;c%$##G5#;j21s>9@*G>#"},
{583.41,292.24,11.27,1.55,"--%UI?KQ'%##Kh,BbH######E~'JJ0###)##,T1T,%###}##<H%?56,o2###VM0;fS|k#*##}'MTYM###J##}dSdu&###f##,?$@uB*c$&##L(:lC3{l&IR&hgS7I+$##~~#GKS2-(###U$#NZ%7#$+##-d(sR*D>#'##xu%`G6Jx3###Hl#K/E;n.###_##"},
{704.43,482.83,10.38,1.46,"lP#Ym*###%##-Q#Xn-.,####/d%4e,UR+G>#O>#(G8dz</,#Om)Au#7[)I?$VR(yu9$+F}G$I1TE|:tL;8?$|l'cM'#.TH>#U7.'##g.*M%&uZ(3v#n'L;kE>#Lg/'@}BWa:kG$|8$M0T'J/?J.###/##Gl#KI)F>#%h,<A0QG#G,#MRGi@-vb#G>#|A)R5M"},
{704.43,482.83,10.38,4.98,")T''y3f5%qb#|z)QU9~#$@u$*$#e97Ln*ZP#(##>Q#|T7###Yl@9=E=#$pG#JE8$wG7s?<6%6KH65If6)mu#<['0-%KT4QG#LhQJ-)hP#o#$ndQ7Q$-pJ3q4bOGD#$cn(Y#9I6(C5#e,%-$&UeQZP#K##Ys7y$L-c$,m%lI(3u#.,#o>#CA.######)##m?'"},
{462.43,501.92,11.23,3.31,"fu$6w)h2>kw,9J,h..Y6%T+<U7+5%,G(&=kGaZMW-)e-%-E<>o(?~.N#$Nh8Op1(M4+^14o0|?O2e,N$&+q2]v)n?'S7L:U6oe+f..A$&FF@q/00I*$N+klMg@O9K2>/*ni=}v'nAOCFCT#%$K1h$)=U3^5$X?O]?'*d%L].M,$h7)O@LDJ0,Q#0*=*s@Nv'"},
{417.92,554.33,12.52,3.36,"N6$`6)Tl#z..Rl$MS-<5#zc%0m'JH'Y>#K.*######A^%nm+ww&cJ0:@*K00=y2P`6:A'}&4hIQkv)'[$(a5m5%>Q$L=8G]2<f.Wo1{d)O4;<3?no3'2'uIQYJQPh9yV4^%Q<J)}~Nli?>Z%1e)hn,|A,7(9zIQ:w+EI'XB0Ke-3y/OKQu..K03Z%,1FB5e("},
{832.56,718.98,10.94,0.87,"{k#oL$Z|G&##MZR+7#}e2B&#9}J######&y#fY#%##Fc$n#&DZ&c[#Y[Rf5$^[R3H%U@,ZS#ikF<]4###|##S5#w8/4l$####Z$###W0K+-'dk;*L9]K-N,$PV'NZR%##D,$9$#.^R.,#######$##:b4W>$-##uG%a)*s/5IZ#CR+[6#y:>Y>#0f.&c#%W;"},
{832.56,718.98,10.94,5.91,"###]y#|w-###<,#gQ9O-(8Z%'-#beVVG#kS1o,#SdV+###x0*c$1M#=L:*##;U8@f%4Z%zX7^Z'r#&L>#2hV)##=c%l##0eVp,&6$#I[U5,#5fV*##/Q%)%%F5L######Dy1fY#=5#g,$3f2Y5$VG#HlF%l#v~K+c$k,%P,#CH:=6(###%##r>#}6*ZP####"},
{519.58,753.81,9.61,5.56,"??%v`VT3>0_34B*fQ9K~T#Q$Y;Rv~1kY#hG#3x(X-*######Jw-L&Q9I%>w+|R.@%&NF4PR)/]V/,#~#$=v$5*;D>#######I;7r,&&d&UZ%U?(U#$3>CtG$R~V$##>?&,Z#ys=.,#######&O-Pv)ZP####<d%<A)M..$##SaAwY#nP$'##Rz2######(##"},
{651.25,895.02,9.93,5.16,"Z,%#########lNe#########%Pe95####.##&h,V08[P#_5$9#$#########`Me######6##ANe######_,#+K48,#}?)_u$L5$######$##LLe######:$#%Me######d$#<g7%##fl%'##|k##########)TZ######(##?Ne######(##[074##U,%$##"},
{855.65,923.26,11.45,4.96,"Y_:#########`Zi######0##0IJ(c$+##mQ%EH2xkNL5#f93n]5######%##DWi.,####{##VIUbG$=,#wd(36#D80IS.'m(]98######$##iVi{k####>##?8R/$(D>#b5#IH&CH#wS21u#dA3######&##^Vi.,####P,#@=G86&:#$I6$1v'SQ%5l$d#$"},
{500.60,95.94,12.23,4.81,"}Y#;7+?u$%##Do0xG%95#'##4OG###^NE1##pb#$##0~Z$##YR$+#MD>####2t8$x1###+##(]ZD>#/lO/##L5$$##6~Z###KJ.)n*b5%;I$JRRmP$###k&#7~Z$##h5Q6$#L5$'##-~Z$##Fl%0,#:5#RV)P96C,$###q8#ZFJ###QaHd##eY#&##.~Z$##"},
{694.38,117.95,12.59,4.86,"Tl$X-(}k#Q5#J'7#l#iY#C$$7.X###ZlP2##pb####Po~###)B)-e+{,&;,#t|=6[*###L,#Ro~D>#mdV.##xY$###io~$##sd+I>#^J,d['d4IOG####/.%:o~bQ$PGNm##pP$a?#2o~###y6)?Q%xu&n$(/S-~o1%##F~)EK1kfKT7.B##0,#0)(q}M###"},
{563.05,136.59,13.64,5.56,"xu#$l?7#$###aI=-]22#####?>3###1i%95#Fc####QE).,#]I*S]*{&4`G#ngXCc%TG#p5#V@=Dw,c,$%##;3)U,%]2&D>#sw0)I%+v'=H#tdX<5#`P#o7#9(7Kv'4e)Sl$g55dG$5Q#<u#iS2*##wb#gR$<$T######'%#{m,i,$ul&gP#xI*`R+0u#P>#"},
{353.10,420.04,12.40,3.27,"}5%84<BSUPQ%oh94054J,<z,9m&OS.zx,MB3R5$4u#:9(jy6Wm#-TUY8/_i<Y(2RSU+J%F_;VTU8(:_H%.V3Lm*P5#~h1R)>u6)#Z$[A(IVUi7-j-(A3+ESUbSUEc%kR$nW8U/.v7*j3Cy>%wG$<c$p,%Ic%&?%9%(-7)*H%qv)1c$Iu#eR(_#$ed($[(~>$"},
{487.48,573.30,12.40,0.08,"U7(m3:f6Q(v'{8*?B1<W;vy2*m$gYC&'+Bw,1mApd*[,$(x-_q&P7QE04*.*U&*E:OE/,^S/H7Q/y1cR&lt@):6Em&x:5r7QH6%lZ'$&+R*?US/bo.KL'W03`8QAm(hG#AA(HA+I^-njGW#%=,#z90gc&0m(*o(a06Rv$N@+=:5J#%+##Jw(uQ&.R(ll&r>%"},
{773.01,650.67,12.60,1.99,"###(##tO6bG$`,#B,$NG5zH*dn#k/3Me#q}LF,#so-qP$jYH4,#n49Z//[P#1?&[M4xD6(-'0U8<u#%@%1oSj>%###,##'oS'##uN4xp0Ye07r6h}ABn,J>#|mS~>$###)[#lI.###6##i$)Z5#-rSxg4K%.&@Lq;/f82-R#+nS95####y-#Pm*G>####Al#"},
{926.85,680.48,13.83,2.78,"n?*######/##d~]######Z.#(4I###8##/y,D>####h~$JK5/&2######6##P]]uG%###(%#;#Knp:###'-#Hu#7_;TH''l#$95######C##H_]}A4(##6S)sn)0]],##NU9H##g]]K5#p&5p99######9##j~],u#+##'G9GZ&Qc&F##'^]$##Vu%K-#r~]"},
{778.43,789.84,12.77,2.65,"_LP+u#xH#}#%0WTD};WZ&U2/5,#'E)jEA;I*86&j>#$e+-u#Lq3C,#%IBsP$0WT7^45u#&c#1J.nR'W~'aQ'8Z%VG#kS)ub#5e.[$$(RTL##)RT&##B,$4A#4C:###0w%2I&$##-c#j`:.,#{,'@e('~.&)(~DA###95#ez$/v(###PZ#{m'######:V0q#'"},
{559.66,845.44,12.17,5.01,"PX:pb#######TMQ[?)###4##N5H<l$95#<##eu%;6%zc(U5#1^3#########s0^######^##t4K@c%I>#9m$^@-t,&95#Ov%H]4#########M0^######+##/0^.,#>5#S{1TR*qb#'6#E2^@n.#########;0^######(##?3^KQ'###oQ$c{*->L~u#Lg5"},
{786.12,914.27,13.26,4.82,"c6*######0##Xo`######<%#cS]'H#bG$C-#dc';6$zY$'Q#;/2######1##ko`######]$#yo`cG#N5$t[&{H*.##Su#SvItn0######&##Pq`######b##Yq`ZP#$##vv$$w<EC9`Q$I>D9.-######4##po`######W%#.,KQG#dP#~d%h[#G(7F.-X>$"},
{149.09,142.79,14.92,6.27,"XV(d09###'##JJ%9bK###0,#J9-]NE###*Z#?/(Gf4$##0Q%i&P>u$.##m-%bK0(n,0n&i$(}$P4l$O>#@$$|I.qb#$##g7*6(P###4,#jP#[+A###;%&Cp2b*G###L5#5s.iK7###*##4(P]RD#########6(P######:-&25F######CR%]'P(c$$##;9."},
{374.97,212.52,15.15,1.91,"BH'###$##/6DfJ1%##x>#o'MQG#(##TS'[[T###'##Jy0XI,f~1###$##z$'J[T$##Cc#,V+=$)/##[?<+h6###+##}ZFr5&'S->5#:5#oY#w]T###uP#XQ$@W?###9}0j/3######<`Tn?*;Q%%c#+u####l^TD>#######tHF*Z$rA2-u#&##R>#n[T.,#"},
{715.66,225.28,14.35,1.36,"###1%)M]*yHU{k#>u#Ay+[IUho5###0##}LU]m+######bKUP7+fc>Eo2%l#zlO,e+B.)s@HwJ3]>##T(SKUr.0]P#dP#u-&ap7*%&dYB$e*xLU|[,=u#k~(nR(l;8~e,qH(RG#KT12u#`P#(c$m,$RS)Hz9lA/'r66,#gJ,Cm#eIU######oG#M<@######"},
{421.21,439.77,15.51,6.22,"o-%fe+|v+]5$qv&.7*K,$5$&h6(:Z%)$%o[(95####$M&0tH2f$tC2wFE%e-zh+>dO_Z$.S-nfOYf3Ll$Y/**~,@@&`;;M%,$7)}I'q9/N/2#jCqv)`@'XgOpdO#y0cS/.%GuJ,`h1=sD)l#S%.6?$UR)jh7(dOgB2(~,l`6^.(o|6|,O+7*&~$eh:qv)&H%"},
{901.47,606.19,16.24,2.98,"X(;#########Zi[%p5###sG$--%Fg[###(J-OG#It>&?&C03Zy7#########vf[<?'###XD0[d+'.+(##zi[###b>$Dv%af[7&1#########Eh[######,-$cdTZG#iY#AL3$##aZ$3C6/S.sZ(#########;h[#########|f[Q5#QZ&|P$###UH$:.+1%,"},
{893.02,655.61,15.26,1.26,"######<##?mU######jH#HmU.,####?##1nU######G##DmU###$##+]+?mU=$)###jqU0}HE>N###'@)nSNNc&###66#xmU=u#<6%-19M5$<I+~##-nU$##?nU2,#fw/jG#vn06Q$$-&F6'4?#c:8*R*###k#&]-%;nUPG#6[H/R*8S,F>#0.((}A+u####"},
{893.02,655.61,15.26,4.28,"xY$%##<v&Cs=]H'*##;KUt[-R@TG>#)%,o-$x6+$##W?##14N6&v[+UR,$Z##w,y6$&IU-##qHU.##l7/r$#BU9$##Z,$E.'T%#GIUuG%###~w+]yT4|E&##1MUm_?wu';##/1/lHU0,#.##u##p`E######@6#_#P######E8$oHU######t$#lHU######"},
{632.69,661.07,14.76,5.07,"7C6yY$###@##8%(Lc%|5&'##e#$-g,IC9,##/0H=8.m>%NI%<~+-265I*XG#}+6Ii=2m(rP#[@+n%%6bEUn(yE?CB0&6&TM&^m)Vd&1|3Nz9v+GfH'9e'.wAZ..Vm#sy1D/-%~*0WT3^7hZ%ybBa-):H$G@+Z#$WH'/$&kf2.,#&A#oQTrY#_5%32&,ST=w("},
{839.90,676.84,15.15,1.07,"B#$_5#sA/J.._5%J$#3xT%I*kYO;##|978N54[*###B.#BwToH$9_2en0$##SR,,U#8vT+##VvTn>#s~2(.#9K4EH%RZ%[$'qS&d18dQ(Nu#EI+z[%6yTNc%.zTQ97aw*0##-~(]xTpb#0,#EI*}Y$?5#Q;-SG#/,#.W(`f1yu#;n.72&s/5a>#&M;D##I{="},
{839.90,676.84,15.15,4.38,"]##oPK)##s./F_#))@V5#H-)CB%^_5%Z$%l#{k#G8$y&14.+.u#H,$Ad%7^VwA*pP$mJI*C9H^Vs,&)[(.7&:@,3,#$%$>4C-v'$Q$YS/(A*~7.G,#9]V`>$_[V$##(x0U-#ph?$##W#$iJ+,v#}lQ'$(###,R)Pr7F[V####^Vi7.p[-(##9q5PWD/,#)##"},
{506.14,791.61,16.44,4.79,"ez>######8N,KC7######<&DIh0p..wZ%)<U`5$0n)'T-FT0,9U######|-#n8UG>#E>#xR#fL,6tD1I*`Q$J{<IV9(v&Lv$I8U######s$#b7U<l$95#gd#xm+UQ%9-(0S$hK7<5#QG#,f#W6U######58#kWD95#<##}G@77,95#8v#3|;UR+3m'1c$<Q#"},
{927.42,810.49,16.66,1.37,"######0##AvS+u####)[#fvSz,'###8,#^yS######&##NwSZP#'##'J)IvS]7.T##7zSx3FivSG5#;I*7<6}#'R?%L6(3E@)6&^e&=^7Yu%7J0V0#4vS*H#svS:l#P6)J-#:v':H$Ed)%c#Nu#bQAnl'-Q$Y-*5o$ye17b2hA31,#Nv$}O8a5%D6#C816c#"},
{927.42,810.49,16.66,4.63,"7&0eu$FZ&1H#p?%mQ;{]72,#c(>}g(*R*z]%,d)>5#K6%'#<Zv)N#$>R)k5$EI,*-#5JT)c#0IT/c#9e.d%#j'7q$+UZ&M6$[m*+]0s#&C-%,R)[10RJTZG#9:Q#uJFA0/##JI'JJT95####TG#TKT.,#J>#'##qMT~?)###n,#[ZK>u$###*##|bK######"},
{993.41,904.54,15.90,2.91,"{k#######$##696###2##]N5xY$###d&#S4G######>&#m[-UR,######$##;'`###J,#Wx)h]6###Fr%NN?######x;#AeYO&3#########[(`###I>#^G#),L4u#^S)fI,yb#7u#l?$8]1C.-#########:+`P6)###2,#2|9<'`sP#Cc%-6#]QPpG$7#$"},
{784.25,913.21,13.62,4.84,">%.######5##.yb######E%#eS]6c#OG#E-#+$(Zu#;#$Ll#G'8######2##=yb######D$#EybB5#PG#R.'rQ))##yP#%?EbK6######%##Yzb######L##nzbOG####u6$NP;)U5)[&-s;@A1######5##?yb######'%#mmTsb#}k#yZ#SR$=)8[I-G>#"},
{790.83,116.57,16.71,2.88,"D>##########MYM######*##Kf]95#1,#*@$vL7(n,>l$ol$OG##########oe]######8##oe]###m,$l]'aZ'B$(L:.H//D>##########jIT######+##<l]iZ(G##QA)$(,w%1[K#caF############TR'D>#######|=.kC=######^)$oe]<##%Q%"},
{762.14,155.69,16.08,3.02,"w{?#########4)Zj>%###,##l=A,R*$##.Z#}H(}80[$)Il$P;>#########q&Z###&##W6$cL8@%,_%*G.(lZ%xsB#d&W-(qK6#########U'Z###'##D?$qXH%##vf&'FAfY#I>#zR$x$Sl6)#########G+Z2-(###G5#am:$IVd##&T16%#YHT:$#wVA"},
{686.85,177.33,19.41,1.31,"ip:9##+-&F8DD>#Bm$.y+B%TE>#&Q$7e%W$TFl%###*##I&T-$T9%#*_<`z$J?'#i)(?Nlu%sjE4%*/d&9:3jQ(Gl$wu$3_6{MB`f#lMA_%#(?&|?$`OA:7+E*@_80^l$<~*hc$V;>U5$RG#<$)mp)-::J##;#$juEg>$Au$_5#9eI2,#gG$+x%A(;######"},
{132.19,203.37,18.16,5.89,"4&2######'##)K^###M>#)%$jy4&Q%`Z%q-&s'4uG%###UH#W]5#########7K^###C5#XZ$.;?###aZ$es8hv+###K##=v?9&2#########'L^######J5#EK^######GA'9AX###x##pvB^v*#########NL^#########(M^######$##XyE+u#&##Kl$"},
{761.74,693.51,18.32,3.35,"]g%I{@a[$9.-Tm>~6);l#^Q'_m*###b6#TV<0~-6Q$:6&*D8%I(.%$V@A%7+4/QcG#VR)q##paE}>&###+d%_5$ku$S/0R?'Nc&Di*e7.~z/m-Q*%*&R(F;.?h,?-Q{5&'##pB%X{A5Z$fY#xG%R#$T>#<8EF>#j>$q39L);#6#.;:'c=/~.<8ET&3,A,2l#"},
{672.53,805.85,17.35,3.51,"0,####,##dB]/,#SG#?$#RA]7#$.,#F##BD]######V-$*C]wP#K,$N6%E'6?.+ge+Vv&.^7ZB]/$'Wm'Wv?eH(j,$2AGlq9Sv%`u%Zd&k088d)(c#@x%VcNR9Y]m(3$'K<98R'M&+)(8{l',@<[l&1##c,%|82###6##,C4c5%j>#TA0^91>m'id'9o/<[("},
{790.87,282.89,24.71,2.80,"3w-######*##fVl######u##^RXC,#_Q(7##Iu#m,$sR-#Q$k84######.##|Vl######-$#xwZRG#cv*E##1Z#Jc$K(<###>07######1##jXl.,####h##6LVsu&Q6(g#$mY#wm%&>KaP#>07######@##|Vl/,#;l#kB%t2<nq<D7(Fe)hQ&mo0v_<2?%"},
{692.50,391.53,22.46,2.10,"Ve0%#####))#hH)|6&:5#21#.,#N]&eg7&l####:/$%DZ-f0nd,l,#&Q%`V#e,%g^($m(5##,]1]%(l@*`T4hG$_#$Op*[AZb6)ox$W@-tu$3f/t_*|Z),##%@Rp5%###/$%'v&Pe)7Z$[m*t,$T$9E7-XG#zV*O9UbG$$##GH6-~.###*##)6$U%(3H&|Y$"},
{692.50,391.53,22.46,6.02,"8Q&###v5#vBS=ZP###J##e-=X@SF#$hc&f#$_>$_G#KCS^#&<c%###'w#$dO-~.###Dz&:U6D@S###G,#Xl$rQ()##(@SF>#eS1g5$_u#K-&}>%pS(xU+Iv(4W@Pc%]u#j-&Ud+###/@SP5#d[+-D.Wv*IZ#|/+e8F0d)4,#W6C>n-.,#@##e?(5e+n#QK5#"},
{712.52,561.46,22.96,3.13,"@%-q?&rd)akA{<3s2AU$)t#$:1(zw/9I)^>$[?&`#%9-&D@(Q5$hG#&?A&C7LU/l]3g;TER*U;T;R*m@*Se)TJ.wP$/7%+`<*~-0##vS*nsDBH'####w>F7TW6T###K%'.0Hbx4k>#.T0ah30FA2,#*m'fm)iY#VG#*`7cR+d~.CQ&M6&Oz*.Z$~u$|I-~l#"},
{712.52,561.46,22.96,5.61,".p/bd*^G#`e'-['I=?ub#Tm%*R'(#E=[(~,$][+s>%z5%Vm$<:.T]/tm+O,$.B*PUL8g88##HiTBcG6#$B$#6r9fl&{k#]-#A7+T-&?_9Z-'KK4+)'`gT0/(ngToH'_#%2:%gi2|&1@H'Tc#(v%aH&j7,XZ${l'pw'k-B5I(iQGqw/h,$cJ&j.*0z.<]4F6$"},
{803.04,746.71,23.14,3.09,"Zv$|=H>,#sd,{P97T39$&N]0$R%#%*I25r&4fG$,Z$5A%{84kx2&@*mI$f[O]</k%/5t3J^9m]O_[+QR'`$(Om*nY#|e%tg6+uAz4Dwl%rv)*&1Vs1*YB0T,([O*Z$UZ&%0%l[*+Q$vQ'_A/1YDIo1###>X70H&fv'b.*^~O^L/.18CT0e[(E)+xS2'6&1,#"},
{969.79,846.63,19.93,2.68,"6#$######$##S$V###$##>o%nJ3###2U#=)6######2D$ed,~./######3##H%V{k#$##x7#cNB~v*LR#`%*<,#K#%d1#>]4:96######P##((V[?)*##yA(9D,Z$V#6#|m++N3s~2W,#(m&](>######N$#Z$V95#?##F>4}u%@H'$)%-ZK1ZAuG%pn&oZ'"},
{511.60,147.35,24.94,5.19,"MZ#GR'>HK###E-%Zv'~&2###U6$&Q%GUO#########}}?###=?#5yK1e.###UR&YKKX-*&##fH@cS1p288,#Nl%###UUO/,#&##t.GQQ&XG#gf4p;7$##/8#%ROM#%ZH&V%#Ce-###UUO%##95#fG#b8&u7-0&2G>#9,#:2)VL68H&yb#$@#vx.gY#Sa<$##"},
{607.42,185.74,24.79,5.68,"tG#Tj1fy8###V$:)D9Hm#1Z%~9+D>#s(#5f35,####s(#=A1&d(xh$/B5/##p]XP6%W>$e##Jv=om+F6#O#%c~#1.,2<(#.,P6)K$#iw*(q7x[X3,#2,#TB(Tg5y[+8R&.&.ex&J0/9A-_Z'95####aI#a}JPA2###C##Sr1h%0J5#l[**8*G>#y[#btJ'Q%"},
{790.23,288.43,27.01,2.79,"jI.######,##&2k######1$#=eZF,#U6)Z5#d5$&Z#-~,iZ&9K5######0##82k######=$#k8^,##uI.F##0-$I6&*95$##8g8######3##_3k######q##[9XXZ&j[,N5#E5#@o(yXJ$##Y'9######F##@2kE>#]>#?g%>V8w97>R&@A,yG$o7*p+H0H%"},
{764.65,505.09,26.65,5.33,"J%)nJ1$##U$%CH$N2?xl'n5$l5%6?&w%/0[%ac'g##4Z%Qm#rI(w7Cy$,_>#9(Va+GD>#5$#l]4AH%eG$2~$ZP#C)#h&5Fc#GR+_@%[YG0A'Iy~'Q$aP##U$<z1yO.5?'`5#$##i)%D_>$##u@(j]0y6)P$'+`2KK/%l#x-$[,$Ta*vA4G5####dD#WI-###"},
{747.13,767.37,25.00,3.13,"dK,|m,-6#OsDL6&Xl%Qq*nlQD>#'##k:%SlQeY#$##:1#rlQt{+]S1Yn%5;>C-<xo4[J+Q?'|m,6,#(''`q8oA3f5$N/*56Ci?'#z-oX<Y&25nQPu#f%-'8%:D<T5$Ad&gD:fH(&7'6S-Be+yQ(|A-(6%rlFQ2:kn.^~,fp1J*.#_;YQ'8,#,@'i,%1R(D6&"},
{381.34,777.48,24.61,5.83,"D>##########Ro^#########I)eF,$###0,#C~-Qv(P5$/c#D>##########w'e######&##_(eoY#iY#H6$E93a-'Kc%}6'D>##########ZBa######;###-e|$+OG#P,#S9+<L2k>%A5#.,##########OL3D>####%##HJ:8g8######&+(Xp9###%##"},
{695.27,875.05,24.77,5.07,"tA2#########-ib.,####S##TROcS1)c#hT+<v&m-*bc$ld*3B5######*##{gbD>####o##^fV*1:?5#86$7<<Q81;5#x5$Yx3######%##thb######=##o/]@5#W#%9J&wm+%W2}e16$%M6(#########Ojb######%##sibcd'uH)G$%Cw,;0&+OF>H$"},
{991.56,922.75,24.75,2.97,"##################$##Q#$######2##jG$######%#####D>##########?4I###@##;<2a$+###](#_JV######](#$IVPl%#########sLVpb#(##S$&V-M6f2^.#D5KE-&Hv)](#<IVr,&#########?NV%Q%###m5#pP=;cOH##yp4o=5506p,#7/0"}
};
siftPoint_save sps10[] = {
{1084.46,90.81,3.72,5.92,"dM+EZ&IK*###/hM.,#0,####cC5######*##dc&/,####'##Ko0El#y0U)##VKV###XZ&I##]5O######D$#Fl%+##*c$zG#qZ(C?$R|AZ|??IV###XZ%+JH|uN95#&##_H%~5$.$$=6(:,#c6LfY#sP#<_Pz%,Kc%.r(eKVbd*tY##1-mw.wY#x@$Yw/###"},
{791.41,96.33,3.60,4.78,"Yd'C94QA*#?%Wi?=E<pl&|>#@eTub#.,#r##%6M#########X[*&R$0lCkG$p|95/.R/./6$HiT>u$###/-#feT######)##?'6Y[#Za:EH&.h81.&,)4sw&<fTD>#$##Gy#{dT######;##]x3_m#](54U6Qf1N6$fR)*fTMAS######h<4S|D######E6#"},
{1049.86,99.92,3.99,3.98,"_Q$,-'X>#~#&4ZFeY#&##gY#d6T######sR#I..?##.uF$](iY####h##;i@qkI###h##i::`7T8?%D>#:~$x/4uz%fZSP?#:v%*v&Ef-Wy5t~-a@-P,#a[OzYIqcAD>#fx.4[(/Y2(n-v##|c#>g0m4D;$)a.%.7TS,#Fe.xb#@gS.,#v/-'Q%o.)###oV'"},
{756.30,106.29,3.83,3.01,";5#N$&M%&2mTr$)=5#iu#)oT0y3######TrTxG%(##eH'KP43D<hJ-WQ%wA+Kq4ml$vO@|e-NqT%l#'?&DI$bn-w5$fv*s5$gD;Ev'f6$TjD;e,IJ(_z;UvCMmT0,#,Q%2<*Um*>5#T#$HQ$/,#.,#Q$#;mT.,#/,#9,#'pTbG$###,##9F8.,####vl#`c&"},
{766.88,105.58,3.88,3.09,"kP#il%T#$*.T/,#m#$r-&W-TA.*E>#`G#c0T7//###/,#Dv;]I'^j?D#$Q6(jC:z~+2I's].PD5^u$i{9,A,w0Ttb#dG$EH#IR*9-'VH#$cD-}@[m(pH%n3CLe-xR'2_:d+;FbK###Q5$(2'######f##(/TSG#:5#1$#d-T:5#0,#1,#P0TL5$###&##Y11"},
{766.88,105.58,3.88,6.19,")##-h0L5$###4,#zoS:5#0,#J$#8mSSG#:5#u##.nS######Q5$[q&SYL$##Z(;Xt;Un-.~'0R%IaCx|@Qd(g?#1mIPm)cu%dG$B$$5pSiY#4E:Q~,0;5@Z$i-'}o-d1;d/-.l#ac&Q%(_|?###W8EEe-###jG#VoSp$):5#8%&#mS/,#]Z$Vu#@mS.c#<$'"},
{806.36,106.57,4.00,0.10,"q##ixV95####*.$09W######ma6KK4.,####GH&U>#IZ&###Dm$sN;8d(b[(ny.T+ElH$n7+4:WU$):,#6?$BT4$##7v&S>#Iu$S[*B6$Ja;Zv)n^9K8$;t?w9Wuf6o>#*p+,(0U6($d(5,#zb#N9WP>#Ru$'Z#L8W@5#lY#t2(,7W/,#=5#:Y1<%-pb####"},
{940.59,106.96,3.88,4.45,"5Q#@z11e.###Up(?QLkU7uY#-v&Vc%U&F(wS)R*###=6#C{Sj5$KL-R3EEu#5ZE`7,;m(:1%Aq:D>#[,#>Z=.,####?##SbGJd%9)6[-*N5#ayS?H%)c$L7#h7S######Nd#Gl%###8##TR*gx+<x/iu#n@,?xSCc%1,#Lw#n*H###.##Yp&6#$###QZ$,['"},
{1002.01,110.76,3.96,1.96,"b##q07>T,eY#VG#g?&<90$d(hP#wT0{-'?Q&G~+6m(<5#B5#2B+'HPm#$L.-X/Yq$+oS&<>?C$)6Z$Aa+/4El[,fY#lv#*&,Y$*{k#f%#*/Y`1YOG#,d#Jq2As7<Q&B[:y?)gY#E>#yr(Vy7vP$###W[$FS/n9,_#&b8*#$'5i%:}IZ@'+u#W##EbBH:,v6+"},
{1002.01,110.76,3.96,4.26,"dc%j@+95####ITXzY$F5#W6$jiC,##/:Ss?&M$*&##m^5Od$.f/be+:-'_P#6TXa$(8#$/8#m`C3j*uRXf?#Op6Rd%I]4|,#U7-eP#s[%w%/MVXy_62H%)@':y,%pF18/2,#bRAjR*fY#&##t[(L5$6##,$'>w,Ml$g-&=A))R'M6%t5:{,&AZ@[P#sl%###"},
{980.87,118.32,3.49,0.44,"3u#:_4###$##)05C/-.,#Mo#sm,gd&tY#]/C###s{'`g4:6N5&/}g2ZP#*##oAYH$'pb#='#s^:qX.=A1FJ####fO.{e1TG#Ho/YH&RQ&$Z$NFY_B2mP$/##P%(Rk09YL0,####$.#m]25?'X>#au$d~.6#$Ui.e7-~l%3Z$.&+k@*R%.?5#I>#0u#zl&H6("},
{980.87,118.32,3.49,4.24,"4c#wb#we(xl'+-&e,%P%&S>Lh?'e#&Gm#fS00x/8l$A5#Pc#sK0Wl%`G#su%N7S###Ae#~oSPc&###>u2|V?(/1###km%W8((;88%,###/##1U[.,#b-'b.%eK7-##EY[Q-&9@,$##+h5Ue%IlKTm)iY#~,#^S[O$%4v(C0#u[-Lq%YS[)?#/f2g#$9S/F6#"},
{1007.43,121.03,4.37,4.21,"7#$d?%8#${,$_H&Q]/[P#J>#'gY5Z%f,$;w&l(>(##%:OL7)###=n#:%-QG##16t00bR*[5$ifY2J*x5&xK#VD?`*+MeY}H#E6%C[&L?&;Z%hx2Q,$R&'tf3|hY-u=X$(~%),;21pF]w-3,#+@$z7/:5####XC2xY$4##wc&v~1Gc$n/+So*pJ.>Q$z5>T#%"},
{990.00,122.28,3.76,4.28,"(m$P5$wY#b5$9~+C,$Om$eT3Zm*/,#xe#|.VBQ&tb#~J$-n+C5#T7)}k#.l#3.KU6(.,#Ou#F/V###T)0Ag1Qu%###(gH.I(0##QT,g,&###O/V(w+.,#%$#).VdQ#N_<,T#L?(ed#u.V1Z#UG#FH$Be*0H&~/V'6&]#$?m&t/V2v9Uv*k##`J+5}/8~/&##"},
{153.37,125.27,3.35,4.66,"UH(v,#2-(M##87X&##.,#I$#5S[###pf53##6#$###xT_###jI.F##uG%m$#(eY######i$#?T_###`08-##?u$###vT_###bR-&##95#^%#:YL######v##DT_###'B4&##k>%###+U_###3v(#######$#gh>######+##{.[###pR-&##~#&###*U_$##"},
{1022.13,127.37,3.89,0.01,"UG#FxEM%AoGOh]-R==000yR*OQ$y->~#&$##$##+`(%~.###w6)nVS#&.fZ&ISSTU1-6'MH$nH(Nf$nOJ$#####G%#%+F###j>%a>#?](%<@LL9Pl$*6%mq7Ku$8m$JN>Ad*###$##hE45[*;,#?u$g##/L9gY#%##(6#=3D###(##,7$CbL######O~#REE"},
{1022.13,127.37,3.89,3.34,"g##[YL######Pc#QlOD>#%##1##A}EgY#95#6-#`_>$##-u#Yg'FC:######/{9fw/_#%ul$q5%k:4ziATu${T-6M:.$('~%HGB{k#/,#P##-AW.,#Sc%)%$fm+aG#i[J3L0q;=OH'hd)toF``D###R>#-K#jT7###8,#q<,e-)OH%29.ZbC.8ALAWg-*,_("},
{967.90,130.50,3.91,5.99,"b>#d?'AI'7-(=c#if/jR+lY#/E=4K3xP$k(%sS/5I+###(_M###7$#.D<(c$m7,21+g(<.c#=qVfo2?#$o($9PCQrV###y<-###g##s;B###o6*eu#$X=hg/WoVy//'6&3q(]6&@fG)w*qb####N##vA3R5$lm+yY#Oc%FJ'S-(bw%kZ(t6$ec$Lc$bv&9l$"},
{280.25,153.93,3.74,5.47,";l34[*<,#.,#qe@.,####%##d%,######+##/[(######'##o7?###b.$(.*wOZ######ZZ#.PF######)##t-*######'##b[,###@M#,MZ8JZ###&##{t18EC######5$#ZR,######1#####R#$r'(u$Rr5&###)##sG:Qu%######Fv$lZ(######?##"},
{754.52,178.15,4.10,1.96,"3m%em';Q&_P#QCS,Q%95#<Q#~3=bT0)m%t8/eY#[;-{l&Zu%5R*6?%L5$Z##p@S#l#95#I%#l29dj8<I+{Y####5X.5m(-8,M..######`/-C@S:5####jp*:n+c#$cd*dO>###3##aZ&#DSD>####$##RBSFl%######VDSeY#######4DS#########8AN"},
{43.69,194.36,4.04,1.55,"%#####6D#&K4tm*uG%4'#=cD0i@{k#}K-K-9.,#%##Y&FM7,Y>$###1S$RA0{o6###b:(C{7e@X$##(ZIw^'Au$3##6CXoY#'m&TG#VZ&O>#v%/I>#(:4|b#b_U$00NU8N>#x5$_])EAX###;Z%H>#A,$<5#XG#8l#c~0###%m#.y/td,###,##wx)}EDD>#"},
{1073.49,194.60,3.82,0.71,"KL9{Y$&##lP#O>O###xn)77&'$(###g=:q,%g,&###{x)J>#VRDkn0###&##e/[%##J].'/&D7-,##[3[pZ%Nc&###jd((c#k0[DZ&###:%#d0[UH$*800-#G/0V8$7/[5,#vc(0##Nv)qY##`A###C5#+k-J/[F6$7.--[#h6)'x#vbO$##F6(,l#mP$$##"},
{1073.49,194.60,3.82,3.78,"eY#'##gc'?#$h<G&##M6(@f#z70|H#s6T2@$jP#d+1ip:$##?v()c#ac'0##:A[;,#kw.:'$PJ1g$#2B[yH$###19#8B[(c$Jw*rP#8Q&###(E[8d%dR-6##rJ.]o%`A['#####0##rSJ<J0RK,[G#-H&###8QA>?%2-(###u[(&~%U#Q###3,#%Z#{x3x5&"},
{819.17,197.53,3.57,0.23,"1J.fY#a>#di-UQ%%Z$SU/L%+^_QiY#_R(%Z$D:T######'##_d&`H(JZ8Q.-`g0)7+z~*Mu$bHQ###<8%{//J~T######@H#V5$eG#<`T[c&+i>yb#~o,N7(^aET#$%?$E?%X~T######Y,#=5#p(+W~T,u#L5$9##Ef.1[(WA.[G#@$(<5#^pS###.,####"},
{981.62,196.82,3.99,2.55,"kX3<&.*6'W##3%&0X/[vV###2YBUQ%Q%/J##W%+bG$###$Q##$'6v$>o2cc#V@+aD&~vV###twVQ?$G[+},#$J/$##wG$[m&Ul%vS$z`EtP#gp9I_#[vV8###wV7Z#du&)%#{.0VQ$Jm)ZQ$0H#D*+c]6###zL=8B#`q>2w#$wVz,&95#L%#B6'qx)%Q%%##"},
{981.62,196.82,3.99,6.04,">u$$##<17K@'<&2k5#KBZgG#{@Z6,#.Q%<e#0K42,#id'cn+6v(A5#h_9Md&6T4%##tEZqP#jAZ###v6(y>#6;;U#$d$++##h6({u%pw/'##jZ(@,#{EZ]5$VAZ###b8+Lm%_h=;,#0e(1K/6u#a#%p5%:#$95#S5#Z*=Qu%>BZ###6K+hR+2CZ###`e#~cG"},
{57.10,203.42,3.66,2.10,"^P#xH'dH&J-)|Y$|b#'o,`J03l$*##Kw';o1Jv).,#)##u#%kx1bW@sb#)m&lnYzG%2@#tvDaQ(###rW)S,MOG#:,#fT)86'j@.wG%1##*pYmpY###],#(q./+=D>#5l3qQ($##5,#</A.,#pb####Y##xOIuC4###($#3f0h=/jJ17e'}k#B##jJ,woK###"},
{57.10,203.42,3.66,4.38,"_Q(###,##-A+^IY$##/,#F:$IjG}##d#JRo)l6*&H#iy7A5#)179H&$##c>#vJYXH%6#$?&#{)BR|*StLj##(p.`/*Fw.%##ey3'$'###Pu#UOYKa<,u#;##G%'#rOE[*$##cN)>p3~P####?I)D>####oQ$sR)(7)Bc$P-&mG#$J(NU*Rm*4h*qP$r$%ZP#"},
{160.39,203.29,3.72,1.62,"(l#~>#`95<5#}R,`#%'m&Me([P#$[$4+>L[)###98,*n+xY$P-)#d#Z)B-##9nXlP#dn,g;+6[*(##XsXJ&-###&##[(51S/sw+cE,]m+###IsX<R'En.s>#8'3C~$ZnX~G#1c#x#%Ep696'@%+<(-###@Z#msXF//xY$~5#e%(O>3@*F###Q?#(~):m)###"},
{160.39,203.29,3.72,4.67,"9Q&&##$7$.7(*WB###Sv%bP2;c%e>#FOZ^92###S-%q810B-)p4]6)(c#m,$?JZX>#(f0ff$HS0*6#)NZ:@'-H&$##$p0>k2fi;DI,###D5#JNZUK.yH*8##De,@|*+JZ4c#Kq=/##H$)Y7$8e-OG####zx/D*<4w*;#$]-%+m'p.($e+:Z$fy73,#~P#X,#"},
{794.77,205.47,3.62,3.17,"}5&######@,#pb####%##<&-######{c#1^7######3~)mP$Fo/,$'pb#E##C>LD>#`G#+w>wG%$##y6=+&U5,#wY#2,Bj>%IW<Wm&Fe.K>#2)UCc%1c$GI%f&0-B,i'US6'Gu$}v$8'UN5#pn.^G#L~)s,&$zLWC:1l#_P#-w%l)UF6(###5@'Uc=.~.0##"},
{65.25,212.74,4.19,4.14,"###tP#6$&ml'h5%d,%9Q#Yv)6=I###,##|r)VXI'$#h94(+.###*0)#v'95#Fq8?8/.,#-Z#1TWcn,###h(#2jC5m9eA3<&#eu$0n&>u$Kl$zg8R#%$##x8,bVWvnOC?&W6%wy+,XW~?''##TR(nu%e#&#c#Fr9.,####IQ$g.-1?%gw'G7*OS+IQ%`<59l$"},
{1016.04,217.28,3.69,1.79,"PO9HA.Hc%mY#gY~Du#Eh7jY#6-&5##[*]fY#######4w)x5&*k@YH(###&H#M(]&Z#1&1o##d/3[6#N']:##$##0##SB4A,$CL:######R@'x']@g,E-)k##E.+Ua*|RZ&##{G$w,$T6)0,#kv*uP$###_Z%IU4Us595#3##'%*:W,xY$###cd*UG####'##"},
{1016.04,217.28,3.69,5.18,"ZP####M%*Cu#_@.###6v%]i)bG$A,#Mr9?&+###)[%#-'###a$+$##PQ%LQ$1RU###.8+v2)g6*-##vVUR7(###&v#qp6D>#{R,bG$'##WG#FVU1,#[d*p>#.z5(##MWU/l#pb#'##-44le/]#$0H&######QBHV,%<5#%##w7An>%w|.,Z$8H%Gu$Rk2kQ'"},
{1038.05,218.65,3.75,1.42,"9nL.,#(##Jm$is@;5#%M1@Z$#c#R#%kHD`P#iG#1?&Y>$###rnXlQ'&-&7x&}@YtY#[g0pI't?*I>#NFYVZ%_Z&qb#h@,UG#g,Gjv*D>#u5#$BY,v$n~1WH#Tn.rR#{AYdP#6#$(##io2_u%2V=######1A${vW*%$1T4d0#V,%+9#Z@Y0#####-##k-)pb#"},
{1038.05,218.65,3.75,4.38,"+-&8#$###'##;XG$##OG#S9#u]7x]#}o6ug%###Th#;U9?#$5/0p5%>u$+##2wV_G##I*qf#vJ3^8#xvV..%###AI#cXDHZ&je01,#}#&;#$:zV~?%n-+P5#Kf/YK&kvVp>$o5%fA%g6Pu6)Rc&###P,#y,&?TT(c#}#'xb#_w+uH%PuK/,#YG#J0&U[T###"},
{698.44,268.88,4.04,1.35,"c6*$##VH#k<?$##)6%3/%LlQ/##Xc&Gh*LlQ######GfHn?*Wn/Fp,Uc%J?%>?%b,7bA1pP$gE9Y18*V;B#$^>$C5#M@U###QH'+Z6Yw/'##AS.XN0yo5^5#F@UKc$/2?zH#;c%8##z?U)##7?&=@%4+I*Z#,p5zZ%#A-_H$r?UeG#v1>JI#+u#~##h?U&##"},
{698.44,268.88,4.04,5.23,"F-$X,%SY5|k#dS,}#')J*|l&uHM###jP#Gw%d82###[)-OH$j;)C,$231fY#_x.9l$m43*d'QnSB,$7l#gn$vf5P6'[.*WZ#_i*A,$%?$###NV(qb#A{*dx1fsET$(0?#Ai2n6)a=;eY#+##^g$KQ'######XV&|70)##;l$*I#(KN%##{5%###v%C######"},
{133.04,270.38,3.90,5.58,"9D+Nn/$##/,#W8/95#*##rm(.,####t-#R/2######S.$]m+UIJ_5%b##aW@=xX.,#*##(f*JZ&hP#@8'iJ0###@##xp1ZP#Q0/~l&Q&#mxXP{?5?&wU%YxXic&*B*o00'v'###vK$H81###g)'2wX%6#F]3>$#dyXYz+]98###i,=Yr93l$###Mf%XA2eY#"},
{666.04,273.54,3.56,1.60,"NZ%A~)Cy6M>#i{?TQ%}.-Sl$~eW95#}g7Su#+u####(gW:5#'80nY#[5EA?$UNAZP#@r7OH$]fW95#N^5(Q#ol'###-hW&##Y%/(##6~Ls>#2(;$##XM2Xv%1eW###z~-^R$P[+###*hW+##8d)0##E@P.H$xe1(##3n(v$'?eW&##Pm(`u#@w-###dfW0##"},
{130.06,289.08,3.27,1.08,"LQ&j?(s,%|G$y':b>$vQ%;p)r5&###lLO/S*jZ(###.]/N5#8f.6PE.,####Q0]Ov(fd*1~$UA2e?#X1]rH%z,'7##Q{@O5#F>Ghh7D>#$##t2]lJ*,d)2##_.,3s'./]'##=$)p##lw0:##,=Go#%&Q%bI#sfY50,s5&I?#Uu$L{'}I/(##37,]>#{k#B##"},
{130.06,289.08,3.27,4.43,"*6'/##qv+7Z#vL=1,#:I'0M)(c$c,#(<[Xd()@+8##2'3tI'5'7%##Ov)$?#&8[###o8-wq(f6*###|=[,J+ZP####Sy/=X75L73,#M$*&##{<[{5%[..nG#Ex-Dm%$;[@6'9#$%##']-=26cQ&K>#CR+###.TBK-(xY$###ze&vy.<_=;5#xl$(['fQ(2l#"},
{372.90,293.69,3.58,5.68,"RD)?o3J##G?'~.&5%-L&#OK6rG#ol'u'#ybO######u'#wbO{'11Q%XZ%;,E10,d[Ir&1|c&:0K~1;_Q#pI+N@,###u'#o2Bc~#B/2:[&y@.;.&)S+qN3lM=1dOJu$<Z#lb7=&2###A##o?&m##5?'^&#HU:oe#0e.'%#Af392)B&3@##P['E8%WI-###:,#"},
{346.95,298.03,3.73,4.63,"A]4$##2{>+##no59##1~RO#$4%X(##ZR+Ic#Km*###S%X'##.g7(##I+I'##pp8-$#~%XJ>#@&X;##Mn.O>#i%0###-&X$##f09###_kKF5#006'##>&Ind)9%X$##Gn)VK*581###o%X1##Dy7###AEC'##}/4@w'~/,N?'125?.+gc%9S(Fm(0,#>%X1##"},
{307.24,300.07,3.64,1.54,"zwW9,#,d)###Pm)5?#-wW?,#hITcG#-y5`6#R(;'##L07;,#WxWU5#(.,###cw-<n#xwWc#%0xW$Z#Jx1e@%/M<$##diBJ5#/xW%##lc&>5#Aw-{,$WR@`81WE?+d';A,_z0Wy6=5#JOG4##`wW*##}k#$##o@.cA&_FHrb#<s5Qy01x1_P#H'3bG#m^;###"},
{307.24,300.07,3.64,4.67,"l^;###Z93cG#TJ1_P#h<69g0iOHrb#Y..sJ&-u#$##c%X*##RXG4##k072,#EJ,Uq0VE?{Z'`[@U/15n-|,$lc&>5#1&X%##x{BK5#8V<$##T&2h@%8&X%Z#'&Xc#%bw-An#mv+###Z&XU5#V97<,#Ez:'##M96a6#P7TdG#6%X?,#Om):?#[?)###%&X9,#"},
{323.24,300.99,3.44,1.51,"^,N###`08###(A/>5#y@U>##B@U0,#DA0b##72>$##*B4%##M@U;,#h09###o.-O?${?U+##ABU],$9T4C##PWA9,#CL:$##e@Uw5#{&6####S-M@#,@U3H#>@U7,#.f1ge#A2?&##Th=H5#M@UK,#2T2/,#yI-eQ#$eN[R*8dQ>5#A%-1n&{]6%##8V=(##"},
{323.24,300.99,3.44,4.70,"$)>)##v&4###..++~'5@RI>#,oR.7*^R+gZ#RJ0###i7UO5#dh=I5#aq<$##N/1Me#A7UA5#,7U*H#H%-.I#+06###%8Up>#oy9$##[`@5,#x84?##]9UM#$.7U*##6A-}5$up:###a7U4,#|S2$##>r>$##'&0M##G7U/,#>8U=##}S2<5#P1;###0+F###"},
{313.82,309.05,4.01,1.58,"<S.ZG#yuO%##]wOI#$<i@2##05H(##zg:$##z+H$##D>####_95IH#PuO3,#JvO&Q#RtIn[#UHO%##9OEq>#VvO###,u####%;;AH#<wO`l&]vOOu#1;84^.kFG$##2HNV,#fvO###_#&###*|Bq,#TK/4H&r&,W'0EB1^n+pL4x>%*sC6##SvO$##ec'###"},
{393.93,337.59,4.23,4.61,"fWE###fQQ4,#6QQ###pq=i>#cOJ###'Q%,##A,$###TlP###5QQ(##zlLtY#=QQ=5#YV9f%%z=L###1Z%E##mP$$##?QQ###JQQkH$EsD3,#iSQmv&F^7^c#;bI$##V,%=##mP$'##FQQ%##%GL;,#Ji;rG$JQQH5#/:3=f(`DA###rP$'-#{k####tFJ$##"},
{387.29,347.43,3.89,1.27,"=6(4##VZ'2##95#I-#hdU>#$w;B>Z#`06KA&M..^$#2dUV,#g,&H##U7.+##yl&TB$EdUSG#E%U{Q$ge0f?$jZ(X%#2dU?,#mP$t##@^9&##7]3LT#2dU4##LdUNH#O]5%$#8Q&8&#2dU&##OG#%$#j_@$##TL;=.#3dUM##6dUr##*L9,$#ZP#a%#2dU###"},
{387.29,347.43,3.89,4.72,"eAY###lw.###$D<###KBY;5#*BY###],%$##dv*###vl'###CBY###z$+###rU;###FCY0l#RAY###m?(2c#?o2###o,&###jBY###v>%###ay7{Y#nIKMc$tAYM,$_d(xQ$/04###zY$$##SBY$##OG####cg7PH$YbJ1,#B8QVQ%on/X>#8(:###]P#,##"},
{335.70,353.49,3.80,4.67,"cL<###2##M)5hy9###C[$OL55&~###>5#T5#^./###[u%F5#w%~###K?%Hy+V96###K]M&8-0'~###oH(R5#Q97###3Q%7,#p&~###CH&<,#1g76,#W1X-?$3&~<5#(n)d-$F^9###a>$L5#)'~&##)c$###{x3vH$xZS&##A(~el$57,6,#fg9###iY#9,#"},
{290.64,367.28,4.20,5.74,"nD,2L995#.##b2*N|F.,#/##$M%9<D###;5#dI#cA3)##L#%o0+f~195#5##4]QO[*Nc&B$#|U-H~Q&M==##S(%EZQN#%%##]M(SS1###1##yIB:FF%Q%L##,w&w,AA(<'##8)4}[+x5&i5#fH#ac'######~/#H2A######&K#=PKeY#$##M|2c[,D>#[##"},
{865.82,476.16,4.20,4.47,"B##ZQ;6`B###x`,(YID>#$##U>7mP$######5R)######&##nY#xA%oZSCl$t[K_c&_>$+.'n^S######>##k]5######-##+H$gu#;~SIu$`,L|k#XG#v%&qZS######[&#^x4######k##<5#.,#]*0TiB;&2.,#X,#1qNbh?######CV%/v(######$$#"},
{840.43,505.48,3.78,1.41,"}k#eP#'N+'3A;vT:5#1##$H7$wTL#$L5$=%#v?)}-%kvT#########.:#NvT@h=###a##L^O#xTsR'Fl%lH#LQ%+/D}:=#########}$#WI-Pm*###K$#{82u{<W.'x~23Q$W5#G,;UwTpb#######r##:$)D>####9$#17,S,$T,$}m(DH'$##iP#;t5P6)"},
{447.53,507.96,3.58,1.92,"###4##.g7###>u$:##S+Hq[)wu'###Fm&);S%~.+##rw,1~A###&##q'8###B2@###Mh8o6%_6S###%##f[?5d).,#@p(~yP###>5#28/###}p8.,#]g47,#J;S^g:jY#^?#)O.26S|#$Ye(Y>#uu&r,&###R#$QG#ix.###LK&[;A#%)&##Rj2,QQ0,#Z##"},
{414.25,520.36,3.67,2.07,"(|1gZR/,#=##=W*yZREe+8u#C7%:&/]_R(c$###$##&a-PL;K.*=z;[m%1l#A^RY/31,#|6#bU3#s9d.(]l$9##q,LCB$zMB06'###Ve#&v'I,>`&4w5#pu%$%&+_R~>$(##+##;UL9-(###95####>6#V,%4##v5&U@#iZ(+##9&/_R$k-+###eH$/J)Ce/"},
{427.78,519.52,4.01,0.81,"eY#'-#~#&B,#8x12d#W>$A.%Ue/.##e5#EoF###.##t)+~L8/,#[-#[I-'H%%%O=Z#L#%XS+H]RrG#/R'@U.HZ&#z#^~R1I*+S.T,#0d);/*b'7:#$H>#ATNH[Ra43($(ix's>%Y2)baI3,#5:2iu%}>&A,#|23+7+j>%uu%Nd(>;'J=J~G#####_$cR-###"},
{325.70,524.06,3.61,1.92,"###(##&]2###eY#0##8>Jhu%HU:###d6'hb2vmXE##XA1C;$######tS0###u$,###`4E&Z#+nX###(Z$9'%.nX###)U.<^&######)7*###$d(###.O=&##hrXE>#M6'/##msX<5#|Z%c#$%##Hu$%Z$###3,#]P#;C3###LC*D,$x8.###{xEF>#<5####"},
{836.88,526.76,3.93,4.66,"|$##D>######c/#;)A######R~$.7,###$##z5#c6*######;w(0HQ###3c#P`,{GQ)##rA0oLQ/v(cG#67'L?'hZ'SR)~u$`bI1v(###FR%wIQ{k#1##RKQKHQ###jc%G/Fuc(.?$*HIY-&q/1######z>$(MQ######oQ&j$P###4v&aU-+u#v,#~IQ6Q$"},
{836.88,526.76,3.93,6.03,"###gH#j>%###&H$i6'5H&(-&|u$9Z%<v%@<AnY#$##5B'F-S###J'%Nc&###=%'*2SQu%U>#,fE'PK###XI)nD4p84-##l^8###cL#PL;###v94j*/Gq=N##l0Sm-S###)$#UJ+V.S###%#####5A#-iA###L$%n'/6OH###3A'[ZPD>####Q>#HaE######"},
{881.04,528.50,4.00,4.64,"p>%######1,#^~V######&##w(`###/,#:?#XS14##:8Zy,$yG%######1,#t&`######0##`(`@Q&D>#e##MK2OJ*k&`S5#*c$95#$##aP#S&_######Z$#G*`G1;.,##$#G(-&i<rM@%##$H%P5$.H#&R(8FH95#$##9C%;V;-@*jG#*V1T;<KQ&.[)yZ%"},
{911.69,527.27,3.87,5.33,"h.@o>%######zm>|01u#'im%Kk<J-(9c$$x('v%tP$?Z$?Z%oYA^,%/,#aG#>(W)Q%'##$^'v97L~-)R%_<9JZ%qH)dP#1Z$395xb#I>#'DVv(W[e0$##Xj7ey,c%WrG$?v&(R#eq>;5#G>#.,#$##(##Y(WLc%ub#Z>#snUKc$vG%V>#Aq8yP#JH'ic%d#&"},
{785.33,528.34,3.72,3.13,"B6(V>#P#%]Q$z$M0Q#0v(###R1[###d8W'##X>#`Z&x0[&##D>#(##i[&)[(zlOG5#oA,JZ%R0[###@,DDc#^#&<5#p4[d$)eY####v,#]%-(g7###VI#3{7J/[:5#cJ-QA&Km)RG#VfC5K4uG%###$##P?%n%0###*##iw(tQIOc&TG#`Q#,m$^80-~(mP$"},
{464.36,538.86,4.01,2.72,"I>#)N-VZ'###MT2Ar9.,#)-#nm+-~++c$ow$n5$9|C/,#WI(+[)0T$N*F=,#{ZR.6%'-&JK$x@/ER&ilJH@%(##`_:7n,fY#?$).##.HFf>$Z]R###,-%>$%ssC###N?;893$##H>#._)PZRD>####H).WQ'3,K###r,#,{*5[R.,#+$#Du7rH)f?(i@#w[R"},
{825.39,545.96,4.18,3.24,"dYLqb#NH&A5#|7W######M##3lB6#$###;##|K3+u#4,#oP#*FD###|I)O-(`9W######m#$O-T######r%#{K8$##/,#77#+m(###Sw'e9W.8W###iv'?9067W######5&#n?*$##K5#k%'.,####<[C$BQEjG###~(6pR?@tK######rK$|Z)###.##V6%"},
{973.95,544.83,4.09,2.77,"jQ'4?%qb#4,#$(ZUZ%>#$3&$G{;PD(K'ZlR&F,$K[$o3Da5$km)k,%Y>$@,#s'ZJv%3-(R$#p:7tW)D%Z1##FQ%q.%-S/###?o1|k####;,#I&ZO$$cA3`##K~-c:$D%Z'##(R'HI%wu'###b%/$c#9#$2,#w%Zy>#}T8/##xl&Km#D%Z###l5$rl%I#%###"},
{973.95,544.83,4.09,6.25,"$;:###9##($'=yX###](5J>#4l$###V5@V5$:5#RG#|5%tb#qA2###7-$mc%5yZ###9K-<Q$B%.'##7|Zh>$sb#'##87)L#%an0###dl#iZ$[wZ####e(w%'Ng9%##X}Zl.+[l&%##<'/xu%d$+###@5#($$GzZfG$|Y#Mn%ej@9~'h|ZFw)+Z$0H#~PHQ>#"},
{528.15,552.89,3.90,5.37,"OG####I##MZ&y-*###($$ux1}Q)C5#L-$$MN4u#y9395#WJBE>####y>#k,&kh9pb#;c#IZ%j@WA?#Lo3IT&9$'$n=ddVYQ#?H$L?(P>#U5$a^9.,#u5#|H(l@W,##/)9bc$7R*1H#oEW*Z$8n$G[+&##]-('g*uG%A5#T#%QEWpb#7H$eP#EU-W>$It0iG$"},
{282.58,558.76,4.08,0.32,"Vd+Vm#Nc&5##-QQYG####Io#Vo4D/*.,#R`*b>$~nJ###oI&QU8MZ#`$+-##'[SX5#-H&c%#%p6ON+B*F7$#$##J'KEZ&:5#xp6*c#le/&##j^SG,#i6*5##+x,AI#B^Sam+$##|Z&(x)O>M5m(###Pn'(6&4~G.,#ZG#Pl$>)*$m(yd%Ke/I$$U7.m##S[S"},
{536.48,567.04,3.68,4.65,"WD*W>$###nG$FS(rP$sb#VG#-u#&##X$(mG$1A0&##~#%kR%|J.###;.#;*Dd>I###O,#il%8[*###m$(]J'yd-q5#i#&#E',u#8@#l#;a98'YF1-$v7+^#%x-V###A5#qh$(y6J/#TQRcM%SG#He#2{Wyu'`sC$d$Id&hH(jwW###A,$.7#]`C)$#*wW0d#"},
{533.95,579.88,3.73,0.32,".6%D$)###W?%-.&,z8###s>#M,#8nHdP#-u#<%-or5?u#Kc$lG$ol$w5&a5$VH'L=5Qu%)##sq>&p-/Q#XV3781l5$yv#?3X###A%#(U8~P#q04F))M/3&##<0X8c$&##(8&CHQQy07##%[9###Y.#;4J###eC3U'*8g8###H0X'?$T,%u$#X`CYk1Q08,$#"},
{819.97,579.46,3.90,3.46,"9><HC9######SN3)M<###-##aFD95#7,#GQ$'jVpb#cP##?$]+8:sDOG#:##}fVzA3###@%#JmQWQ'Q5$u7#b)-|,&3m(B5#/GD>w-+u#@##7iVW>$###U##B54lo5###:##-S&/.,=#$F>#<O?Ru%~,$95#|gV#########?,@L5$###*##Oz6eY####8##"},
{327.64,583.86,3.73,1.26,"hI%'sA#Q#|,'*110Z%nN,5v(fZ'###hc9O5$CH&s#$h5%=5#8q01B1{k####tzZeY#R'1`G#8[*%##u|ZM>#:Q&Q,#gS0f5#Ri@yu$Fl%]##}wZ'##co4Y7#WI-E##=xZC,#s7-R5#^z=8##H`C>,#ZP#C&#LwZY##:p7q:#>u$z###xZF,#@-((##*:7]l$"},
{327.64,583.86,3.73,4.25,"lw/[c$'f1C,#rRY>5#.,#,%#RWD5V#<h='%####7;#{1?3,#<K5^##3S.s5#}RYH,#$m(K%#F^9(_#mRYR##{k#J%#A<Da#$A80b>#[>$#?#1VY^>#7.-3##mA01-#TTYQG#+u#)##D(5LJ-Pu$G>#z6*V#$Z'K<#$7-(###3:,3c$Ja;6?'@5#OG#4n%j*D"},
{556.13,581.41,3.97,1.07,"]#&'&#P#Q8##C942?#jI.Le#^x2%##Ed$@E/.,#7.$pU0oK60,#8g#.*E;5#>-AR6%o-+Q6%N^S3##$L6xR'_5%U8#4[SV5$j5%y-$Nx3U,$Y2=U>#&-'AW0q[SN%%5_=H/%N#%6U#`ZS&##lu$3-$:J0%###41_H(xY$k>#MP83R%w(>>,#7,#F7#z[SZP#"},
{327.64,583.86,3.73,1.26,"hI%'sA#Q#|,'*110Z%nN,5v(fZ'###hc9O5$CH&s#$h5%=5#8q01B1{k####tzZeY#R'1`G#8[*%##u|ZM>#:Q&Q,#gS0f5#Ri@yu$Fl%]##}wZ'##co4Y7#WI-E##=xZC,#s7-R5#^z=8##H`C>,#ZP#C&#LwZY##:p7q:#>u$z###xZF,#@-((##*:7]l$"},
{327.64,583.86,3.73,4.25,"lw/[c$'f1C,#rRY>5#.,#,%#RWD5V#<h='%####7;#{1?3,#<K5^##3S.s5#}RYH,#$m(K%#F^9(_#mRYR##{k#J%#A<Da#$A80b>#[>$#?#1VY^>#7.-3##mA01-#TTYQG#+u#)##D(5LJ-Pu$G>#z6*V#$Z'K<#$7-(###3:,3c$Ja;6?'@5#OG#4n%j*D"},
{831.27,583.98,3.80,3.41,";'V###I>#0,#IrT+C8######Lq-f{A###&##7N<ZP#)##Ju#p*@/,#=5#F?#o[>R+JZP#J##:'VJz;###=%#+dQqP$QG#|.#~g46?'###S,#+6H5/1L5$W##)(VL5$###1$#kk8U%/###U##/c#xG%&6$J$)L-JBH'{G$.,#S'V#########?j8mP$###&##"},
{558.65,591.25,3.65,4.40,"Ni6#7*fY#6##U(Z95#{m$>H${u'4##A+.?7,Ru%M,#70*:v(tXJ$##&l#O?#i%Z###[d)QA$Kw.T7$a'Zfc$TQ'|x(>K2UG#-cM1,#rb#P,#S'ZkY#6d)_,#w81y#%i*Z<#$5-(7c#]a3.$'GT5/,#C5#*(+d&Z<[*95#c7#)%(m9/1'2.,#g,&W@#^;<L>#"},
{536.54,606.20,3.60,0.21,"s6((&&qQ)###uwW4u####:$#2YK~8-'##6()+6&UpM95#rA)n]2Gn&fA3%##)xWt5$3l$N%#4r@h=0y^<I$#0,#?WS-H&###px4`G#eC54.*+zW<-&3-(6$%QJ,:U'PxWpP$F>#Ew'DB0)o0W>$###.7%a]4C(2WQ'+H%17)uS$o7/S],Q6)RZ#iZ(j-$^,M"},
{796.29,607.93,3.66,4.24,"L.%F829##@u$oj0q~2###oP#/59P6)###AT'C-&`>$vu$|;<>Z%95#_%#]iCP6R###T##Jq7ZLXOG####)4-OD5[Z'}d#5NX^>$###$$#/z:d@-###c###IQCsVBv)###sw+BX,`IX8l#06&6?$T,%###=5#/$$$I*###Gu$An%Cv)####Q$2%%c]6H5#1u#"},
{904.12,607.05,3.71,0.14,"*H%D>####_G#zWb.,####&##^5g/,#95####Um)tP#E,$%##TG#-u####Cc$:2gD>####uG#74g^P####4###]0;I%95#0,#G>##Z$###'Q$P1geY####66#M2ghY####,$#]J0>Q$.,#kZ#WG#zY$###lP#%mNeY####$H#P3g######q##.'6######7A#"},
{295.87,612.15,4.07,0.32,"*U$1_=###/u#ZR#=~U###7Q%$##Y`UP5$E,$7#$nW+l#&Lc$0[$mm*D>#1,#?Z$k,AOG####T&3aC1yG$Se*LQ'f5#*d%ApI%##gm##.,:5#,A,W(*<R+###0~U[#%###:A$.YJ.25/,#U3,###P$#D1<###zo0g~''K4$##R~UZu#@H'P%#]K6>=/iV@1$#"},
{295.87,612.15,4.07,4.49,"NA'{H*x##~e/dE1{k####rG$AS*%##E>#dP#95#&##7Z$}H)p#%###[q$JOHG'3###C7#Jh;0>I###-##0m$KQ'###$##|z,GbD###L}1@R)'Q%j5#VG55/1#[TI5#hQ%?I%gXJ7##hu&yz%|ZT(##4nKw8'=u#TS%.]T[P#O~Twc%8Z%C##Q[T0$#@o39.#"},
{352.14,618.47,3.43,2.37,"###]##waH/,#j5QO5#}6*mR#)dT5##6f2Iq%.,#h.%7YKNl$######QYCC,$IdT###Md&U%([eT(##G17Oc#a5%G##+iT'Q%.,####FD7>l$z:=###cc$Z^)seT###z,$T7%g06###H2'y`B;5#uP$#~)<5#D90###MQ$%Q#CoH*6'###2##-Q79K5}5#RR*"},
{534.09,629.82,3.77,5.46,"###*##}$*Tc&PG#Zc#W`;fG$>I+7l#^K,(7=OR*tU9###I]F###-##IU6D>#rI-2.&D*BWG#E7U3n'e@-0f%}R+)<U0U8V?####A?'yA..,#G7-lY#-#=qY#R7UH,#/)6du$Rv))7%)<U9H&MH#m>PjG$=5#3-$~L:e7*=5#*RFC$)3d%n>$8_-T,%tq%aH("},
{461.71,643.55,3.98,1.16,"0~%&2;zP#'%,r'2T,%h;),80[Z'###LG5_u%RZ&eu#`>$H>#8q.5w*Jm*###iLXSG#0g0}G$4I+9##cNX,l#ml'T##v7.K$%&{=jc$'/1}##}IXk5#U97&T#Vd+h.#xIXI5#NC8]5#KJ1_R#W1=C,#Uv)Hh&'ZP}v#*V=-V#sb#Y.#PIXN>#qL3vP#37,x>$"},
{461.71,643.55,3.98,4.21,"cc'd#$_E<wP#t-WSG#yb#L8#e;Bm^#Jq=kJ$%6&.E(o':mP#S%/'/$#04Z,#3.WL5#:m)z/#@g8yg#+.W,-#A7-(%#D)@>6$>e-vl$<c%y5#'2W6l#PR,E##.B0}#$=0WJ>#.7,&##sg0FR(B#$@5#kQ(Y#$I[<au%BH'###{U)?I+|_7Fl%]5$IH'Nn%F;;"},
{225.86,651.18,4.25,2.94,"-##)uBW>$###V?#J[S######w',|R./,#$##@c$S>#p,&/,#j-(d6AU-RRu#nz.4~SpA1m#$z_S.y5.c$|G#Do07,#,m(:,#t5&+##ATG-YFu@.9l#{|2O]SA[S`P#($%1[<cJ2gP#'d(7Q#$##*##i~%X08zY#G>#oc#D[S+c$%##>l#S|<eY#:,#jc&O,$"},
{885.03,663.77,4.03,6.10,"1##I?(&%#$<oH#$eY#`##G<o0,#TG#M5#:<o-u#95#5##C?ovA,nL;>#$bH(7G?sQ(=Q&&'3US.ke)>T43L82'4s,&XG#,q7g..0Q$k5$McJt;@nP#E$(Aj.UJ14Z#cE=p'1@(<?5#?#$4]'######d##v|FeY####@,#?L0|k#[P#(H#?T0pb#:u#7#$qH&"},
{598.78,670.03,3.62,2.55,"+c$rq'1I+&##gsFso/OG#M&#UJ1.7)pb#-<*2l#[WA###'14xd-S@#B3DAc#e?TcG#j#&KU#N'89o'bjDB@#&##,Y>-I)L5$Wv*&##b;8/m&YBT###l>$?v$PbD/##`#;`@-.,#O5#HU(x#S+u####@K):?&IQJD>#<###R%}SE_5%^##%%+_S&M$*R$#GHQ"},
{459.17,680.41,3.92,4.56,";<6w,&ZP####21Y###g$%g$)$?&)[#@r6y7/UL9T6%k6(~5$.uL###}k#<,#owY$##oR(~8&Wn/jc#cEUNR(h2?OZ$uf)RZ&7jD###j>#Y6'X{Y.[)&Q$HQ#T0/$M3<7C<l$=x1L-#5$H/c$'[)m5%`,#ih574@%j=###vd$R9$ZwRuR.###}Q('A(8g81##"},
{485.64,688.97,3.52,2.61,"T>#wX4Nc&###e944D8###y##)[)Od(oP$T8&2l#T_;95#a@*&~,wp)zsH3,#QwU{c'`,%UR#-]3qd$X4HMv#$##L{3(%,95#'@+,##7eJ}l&+xU###Mc$C7&SEB2##s??w81###F##QM+/tJ{k####er0;S,H)@###A##Za.?C995##$#fd=r,&Nl%WB#RwU"},
{541.30,700.94,3.63,5.42,"rb####-##t-*_v(dP#8~&#C4pG$hG#Pz(2:W%l#b#&4,#c<WaP#`P#]H'eG$hq0+%,el%6H%w8Wq-(>6%=E2xd+p2Wtc(`21###%##(d'rb#4^7'l#Mv%Pm(d7WUc#}B526$2I)-a2d:WJ>#|k####{b#,l#.@+###1##n?(f9W###*d%Sl$OK/:5#NH8&H%"},
{369.39,706.10,3.98,2.48,"I,$+U+du&$##6C8rS/###9/$NQ'_c%cG$A-;A,#Zo2?5#}sE^I-XI$uq>G5#37Vj>$1Q%X'$g]6Hd$=bHb8':##Yr9wT5]P#*R*%##z{8*-&L9V###*H$jH%y#L'##B$<(.+###4,#1=.MK6X>$###WK/rG$6TS####Z#+6$/;PA,$:$#%~-)%&]l&%_#wVA"},
{495.41,706.23,3.58,0.78,"{$&@T3:8$ZcQLhQ2-(I##>01tfQ######(##i%/.,####<,#&6#OK5Fh%KcQ:`@rb#dq&UcJDeQ###I5#p5$$;=######1##F58mg:CQ$N>#+U5M5$15<pl%NcQ###Bm'N-$5i@PG####5##AdQ:5#tb#0.#Ez<3##mkJe$$b4L###{?*g6#)A/k,%D>#7##"},
{356.24,708.67,3.76,4.71,"s$,-##w(>?,#fmU######T##EnU.##d&1JQ$R5$;##QoU]P#F-)###B*;iH&>nU###&##s6$]pUj$$wR.<##$H$gK#RmU###~c#eY#zp)Dx0x~2nY#8##oj2j`C@2&DI,P7#bP#y2)rK8###I8#>]4l>#al&WQ$6mK)##X%)qP$/v9%Q%T,#,I):U*bG$C##"},
{447.52,709.30,4.15,1.64,"(##gc%w^V.,#hi5N-)Q?'###Y^V###c&.5,#<&01,#Tn+:5#-##Du#]_V###eM>;5#p:5mP#1]V`c#HJ/Y5#C]2PQ#:05$##%-#W7-(2PNc&dR-bP#(a+o<@BmRa/'Jc%A])g_>VI%0Z%<##p%#v[VUI%(c$(v#`?L]])380`H(xn&cH(I.(Gp8~>#.,#(@#"},
{332.00,727.05,4.63,1.72,"W>$7I$GeQuP#nr2;J+qe0###K#B95#0<4QG#D>####Bo,lY#fZ&#c#,:S?5#&E>Q>#r6SG5#I8SN,#bz9hG#Z>$I##]+G.,#*&&pI,q8S{k#t.0iY#JvFTW8K6SQw#D&3/p&}k#`@#p99###-0#}bO2@'+u#WH#k6SJ@%ag6{>%5Y9-H&6R&XG#Vd%;c%$##"},
{885.96,764.49,3.84,1.56,"0X.2D=XZ#pcLJv<`Q([5#|l'dR%Df3.,####r.*)c$um*%l#N38/~KDE4`oRToRt&-6^06I*1z7fZ'6c#>c$vV;Jl$[-*6,#%C.:oRSn(%s:|mRD02xY$2{'zx5###Oc#3d%E@('-&;H&l5%6#$###PR#FoRyH*###%##nSEmP$###<##1:7/,#$##2Z#OT4"},
{894.74,773.89,4.04,1.27,"TzOQ%/]##Vn.^;3z,'###)##M-%d,%L6(?5#%##HL/-z995#T/T}T0Yk7+N=tvHrS2aP#~5#*k=A-(T6(8l#Rc&nR$B.TuP$U.T|N@)d(`;(Z_>R5$_Q%7%%2?NWm%p-+UH#tZ(K0$R-T~G#'@+###K,#K.@'-'###[?%1e+4n)W@'Ke.e#&BZ%3o$NPMjY#"},
{620.07,786.90,3.57,3.15,"67#a%X[P####ti,281######<y3######-##Z?(######/##<q+~?L.mNJ#$h*X7w-~P#f,#_@T######9##N$(+u####/##N@,L@%'$NW47`%X###mY#nr(#>J.,####6$#1$$xY$######.,#P>#yb#h*X^Q('##E>#:b0GH&9l$.,#j5#U##-H&######"},
{650.40,813.22,3.55,1.76,"y##ad)eY####$##En$i,&/,#.,#II#`}D^@-###B##8e@bK7LZ7g7.3l$f>$Xn+fH$c7.]G#$OFS>#0I%@=<0$(###+4.*FD'jVE@(~R(Vx1n?H0,#7-&JH&%fV###E5#C6%FC4cP#Gk;[5$uhV=YA|,&W-$nfV###$##gc#22:######Ml#|/,md(>e.%##"},
{895.44,844.52,3.54,6.05,"B,$*##GJ,TG#>]1/,#=l$mv%bH&#$'XG#FL6H>#Mu#VI'v@/Iv)Nc#{%1V##*5M&##(c$SY3WH(###E,#gzU###.##:}0qvUp~/g>#YR*|b#}{U$##SG#om$'-B=Q$Ec%<q6=5#[g#SwUnv*^5$###|.(W>$?%?ZP#7?$D,$at2JkDbG$$##=,#-Y1Lg9YG#"},
{998.52,850.42,3.87,1.49,"###c##8S.A,$^##iY=7&2###Cn$ip6(c$###@5#_?#jI.######i##E#O.,#oe(ga6}IWgP#{NWV&1o-*OQ#Td*A-#w=KE,####k##KbK###8d)%Q#{NW(x-6JWaP#QS)`)0qw/5##+KW3H$###B##bo3###95#$##hr.2T3<J0###76#j#=ol'###B+;TS)"},
{998.52,850.42,3.87,4.56,"VRH%~(j,&###/6#(H;%80###3N.Z]4D>#$##n82######H##t%WVu#X..N##u%*LW.O%WcP#K*W6x-Lv)@Q#e=J######'$#i{B=##Wd*W-#w6*WH#K*Wn/16%WiP#pw(UF4AkK######$$#?d*###B5#wQ#W>$###k.%-h8T/3###Z##^b;0~-6#$###]##"},
{656.71,878.20,3.37,1.33,"El#T,%S##gB77H&B$).##vf-]%/vG%9##l%@[P#iG$-z#`SX;5#3c#,S);w-q`E3,#o#%1o)1SX###+##/a+%f04,#eC'LjB8##LM;vG%###-?H>%-,u#(##pWX:m)###[##G02ja;/m(P##Kl#Jf4###s#%Ax&Rd+%##?6'|X/%*D%##zY$Jm#+WXPG####"},
{630.23,884.65,3.82,2.97,"*##0v$gG]v6+Lo2gG#`^MJ/)FK5######nI#ZP#######1-%L5#AG[vB]###*M=`'$`dSK-$Xe/L5$###A-#z,&95####m>#oG$6x?Dq<95#[bC]d$$m()##.],wu'###%##cv'{k####$##E5#Wl#<nUeY#.E1q,%?n.###`q3.,####*##~6)######=c#"},
{640.74,885.58,3.97,0.39,"6v%###~-#b./?r1xY$;##]u%;~*sl&R5$MT(/,#&##fH$$EW.,####BI#r[-8D=###n7&Tm(YxY.,#&##84-4e,c>$z6(n8C######'?$A?'4T1>#$$L1C#$;}Y-q;Pl%8v#+X*'yYl>%9$&###,##Yu%D>#+##1l#c^8.,#@]#3L8oA.F$)aB$vwY%##/['"},
{650.67,890.29,4.22,5.95,"###+##b[*eY#7m)*##wm*N@&DJ/9l$>,#Ko/dP#bH%~%-EZ&bl&]Q#s?*a>#6jEI5#mu&O69K?($##UR#a^WOG#4##'c4S}K*/-37%/n,A5#8aW/Q%WG#+e%k<7$D5Ad&A(8^P#9i&@wP=I+Y#%3,#An'bQ(.68hv+g5#Z$*&1&E^WhY#@u$pP$Ak03-(UH%"},
{962.71,901.31,4.05,1.71,"3lCD>#8'-cP#g[@>o1s1-Qm)4##q?(Wb5O@-6#$###?[$Z%.*I%N5$<k3pb#~<.ZkK~/'[?(7JQ4B330,G{6bG$-##iO6y@.###$##D{)Yl&^%.,u#_C$wiAQLQX?(=&,s9/Fw(`.*]KQ2Q%###T>#$d&OG#:5#$##N0/i,&*K'_S0Mf/;#$zc#*P@9&1;#$"},
{991.95,904.89,4.05,4.02,"(m(C-$z84v$$9z4S-(W,%7##xPHA-'3d&ae'mI*4)/>-GH-'^d++Q#_f3>Q#~aAh>$Ox02c#{nUT$%`A/RT'PA1)B$BpUu-'8w*Vc$Vd(id)3M:3$'MS+(n&RoUKn@:$)p?#}o5)<)%r@C$#3A/%v$926o47~$&DB2=%L#7+o-'t7N`G#gA-QK6ud&D>#l9$"},
{935.78,920.05,3.74,3.00,"+u#:5#9%#<_:Mv)~P#9##Oa5IS/PG#/##u8O;c%######VeR.,#%##&z(+dRJ[(ER'n>Ekh:#fR<-'b5%KQ#MB5T,%###r#%ku#NqPOb?}cR3H&U]%@fR+g30dR~>$iu%^C-Nn,M@-%##R6'ed(SF@c>$nB2.,#~H$wl&;6L^#&###%##elDhY#D>#'##}p6"},
{935.78,920.05,3.74,6.06,".##}y7`P#{k#%##LZGj,&.,#)6%VvO.,#Zu#3c#3h7D6&LL5%##Wc%GS-Z[,.H%+_,]6SS5$;8SV'4'?&<'&Wa<'3C0$$J;S###)I'&y5{k#b5%8v#%8SA-'T[J<M:{Q'87'WN,N6S/,#G,#'##]7ST,%###A##}/S4w,PG#T,#6F6`l&~P#VA#+1995####"},
{397.92,924.19,3.48,1.44,"6,#+u#h?#;J0IH'[l&K##=T+/S/nP$?##50FE>#8#$7:#^eXG>#fG#X/-iu&V4K3,#~H&(~&.eX###/##4X-k%0>5#Eg#MPI@,#='1+6'###OZIL$((c$%##ajXTR,###7##9917@JH#$3l#.l#=d)###gl%T%(n>%$##T-)O}/r1>###|Y$~?#BhXE>####"},
{383.39,928.37,3.84,0.74,"Hf,0Z%)##tl&>.&l/2fY#F5#{m,5Q$j5%iW(o#'###b[$6bZN>#OG#u##B~/B;=3l$u##lH&3pZL5$###ZM#bT2EH'q@)GX+######V##;[*JK3/,#w-$#m'LtZG{AVG#5H#H<(vnZ9l$kG#OG####I,$<c$@5#7,#3n,eY#,I#('5f$*6#$S.#eRXB5#u,&"},
{391.46,934.51,3.77,1.08,"Z>$bZ%f6)~G#kC=%##9u#6E)z70###6T#)1M8,#fY#)=-n?*TS.<R*:5#%c#fAWW>$###O'#*jBce0g6$n8&###Ow(azN}k#v$*95#D##}H*^hOBB6)##V,$mB)}@W95#0##/,#oEWK&3###/,####P5#sl'kl#m6*wu%M5$A###=EEl$7l$Q5$7U6$##cl$"},
{287.26,963.27,4.18,4.62,"{Z#TO;f/4###{@>905eY####T18######7##{G%######:c#I.*[~(XTZXG#:XZ1R(A7,Q,#8TZ######'##O%-######B,#xu'%##(YZm%-0SZ$##|e)nY2^5O######&@#_J0######.##{k####8q*^g8|Z)###IH#TJENc&###*##_I'c6)######U>#"},
{752.15,987.26,4.04,3.61,"#########'##ZP#######EI'5U8###$##yL(UV;ic':5#]g#zc(######P##;NB###8##Tm;-p5###y%#JqV+y2_5%1##Tb2Fw,=5####vG#,qV###,##MB&CbH###O$#.fI.,####+##~d<B,$$#####,4.`q=######R40g/4######`D&#########@d&"},
{706.05,992.20,3.87,4.51,",|C######S[#:W>I6'|b#pk7A.THl$TG#'8'K1<###$##<,#RuG{k####wl#Q$@<bHuY#f8.^/T[S/I>#]u$)+I###%##<,#d-T######RK&n2?06'(/&&zSU-TC,${Z$H^.FaG95####5##;/T######5##)%JY7*qH'w,%x-TQZ%7H%;-$zy:######B,#"},
{696.07,1000.51,4.07,4.48,"DH'######(##b/Qec'.,#W?$W6J+^4*6&ZR%c,QG>#:5#-$#gd,######>##i-QCu$1,#oN.<lJ;I)G#$,?:(-Q95####/$#Bw-######<##}0Q{.05,#sA*m/Q(4D6l#f.,=-QeY####C##%T3######~##J-QD,$Yl#Vt61YJ-Q%3J&n*>a,Q######B?#"},
{638.49,1010.43,3.84,4.75,"hU8######*##R0Po.0###?[#X/IVe0###8d#Q`>######*##4]1######-##ySS^Q(###})+C6LKQ'###]c99HO######AQ#TJ0######$##RVS+u####N-$+~H;n.###1J)Ij@I#%###R>#Ul%#########JnA2-(###8##.}073E$##Sm&-%AKQ'$##1u#"},
{368.36,1025.15,3.78,3.26,"######}Y#:^SOG#%##Qw(BAN;c%###e@'ND1qQ)###C5#3k3######Dp1Vz;Kv)L,##vNB$(/ZJJ6(?Q%3.?}Z'Pi5MH%z_S###mG#G;=###dv)'##l[S`P#gTIc4K>H&*e&f[&z_S;d',w)iY#]Q%rv)K6'/,#%##~|76H&2m&^#&cQ%y7+*-&j#&Y?$A/."},
{331.20,1028.87,3.56,0.45,"###pf/[$$HHP:v%sp9QS#c1;&ZIWv*/##S_&pR-Z8*xu%_b5.,#s9-,E2:8R*m(`$$s~E'W=^8RBT3|>$?{(Uw+b/GYl&I-#@,#-;R,m(~5$~m%5q0_J1{c'^J+8v(vb#gR*.?%mo*qQ)7,#lG#)&+=?&8#$oo*gu&$##eP#od&5],###kY#RG#TN,ZP#M5#"},
{1059.92,93.53,4.09,2.32,"/##7D(,z;###e0**KR|Z)m5#^-&ZO<'g6r5$###y5#^1=###&Q%r/#?HR&##VHRH$$ue1+o$|m,W5#'tC#8)###h##Oz:ge0vG%1##AKR/,#QJR%##qf0Iu#=*?$##]p0TsC[P#$##i.&THR######Jl8###{g2###1E+P5$Vd>5?'B##h7,>o%wJ3B[#:(;"},
{1083.99,113.43,4.76,1.13,"i,#f~17,#eY#Fl$[d*{?$jv)95#&##^D'lL;###'##tV'Vp9&B*;dT###2Q%$gT4:8k>$(2,BI+^('p&RO{76##X@&4(4g,&J//KJ/eZ$#bEHiTd~SC#$8[&He(HiT%80W>#CI&M%)T,%###0,#0z2].,/?&-Q#k+E_P#E>#5p.T:9###%##jB4OG####2##"},
{1083.99,113.43,4.76,5.14,"-d)nY#M#$;d$~K6Dl$4c#'0)O6($:)zA/?~*d09p##<$)Ah$9@,###>-#PV61[T>5#r.%QxE[H(g>#FY31^TNc&%##:c#JO/bm+###]##g/1u_Ta@,nZ&f7'7'.f258~KKd(OZ&dI+A?&)[%U,%[5$4,#wQ'r]-~p6###d5$/m#p]T{k####1,#O_T.,####"},
{798.20,115.26,4.43,4.75,"B6P######$##:AWSG#lZ&DQ$=[(4J-L04*l#l?Njl&pb#5##8[S######-##iAWP>#kQ'[>#M7,zu#dtCFu#bwS{b#ju%:[#3jF######W~&I@W7##x#&jo&7n-a6#fM7}EC3tH9,#N,$ur7wu'######%hLAL:.,#-##Fb1W,%0u#YQ#eAWdG$iY#8,#m}E"},
{252.59,118.12,4.79,4.66,"]#%######%##6V=######'##W$U###{*H###ZP####5:c###l[,######7##1iA######@##W9c###[tJ1##A,$###q:c$##[[,######'H#uh@######]5#O9c###_<C^,#>u$###~;c'##0Z%######X,#Ny8######O,#Q9c###kg94##0Z%###u:c%##"},
{782.01,117.83,4.30,4.65,"`}B######$##KOb,-'###G,#)90iK2S>#/o,Pq4yv*Tc${7.w?U######(##yKb$##Cu#Zm$~p5do2d'4k6&ktHa^7FQ&Ic#.-S######9##?Lb`P#d,%)?#n08Av$=IO5c#u[N$6%B.,^Q#Q;A######D%%[Kb3##pG$=&$Vx3Z@#>3=(%*YXHe>#'d'%K("},
{825.50,129.31,4.23,5.75,"Hi$C1<vY#~P#KTGo#'Gl#C9)g'7?l#7l$#f#F>#wG#EZ&)##cI>z,'###H##7#]{k####6$#y#JL.+95#,##(6%JQ%D>####'x~######M%#?x~######7'#s*@$?&###1##$w'ZP#######rw~###&##iw%%x~######m$#C^4######%##v-&95#######"},
{784.92,161.60,4.14,1.48,"###v%(L@-%##Gc$/[&pG$4d({'HF>#4##R6(Z,^1c$.,#iv$_5%9##pS2fG#U82<5#d>#EU/:FE7u#kQ#_{9B*^.,####MJ$8Q%yY#kn02,#B;5AR*2,#%Z#Rr95z5yY#gu$$+^cG$###en#6,#L>#OT/ZP#ZQ$JH'>5#Zu#[p8Rv(###o1%}&^1Z%&##$3%"},
{180.41,166.14,4.07,5.41,">YK9c#+c$|I#nR+p&&;y7%##+yVDZ$}>&###zD:6#$######cE=gl$/05%m#ze-$D)=7S3Z$6zVR-'@l$66$n>MD>####&##Et0E.-K[*$##A/*rI-T4;7z4?wV0,#Hl#S40*+F######F##8f'Pn/KH$[e/=>4hS2;##'I'ETH95####J-$o0-.,####(##"},
{217.71,184.78,4.75,1.32,".,#&##;>64p7,FFVZ&+B%_1QI0QK)=L5$$S$sG#)B,}>&######E$#]C7}>&~6&G6<{`Bc.-_1QV<CE,$-Q#}f3d,%X>$$#####O~$F2A###lH((9'`-QLm'6.Q?H&c#%We&Ep51J-.,#(##F6%+}:0S/&##H.*8(-O(<G>#g8PV7-0,#?5#a-'G=A######"},
{189.92,187.43,4.33,4.50,"B8O0A-L>#`9,tL53G6NK4$d$S5O%v#rl'iT$W-)aU0|Y$>$&cl&_P#^[&29OzS0Cu#zF;Ow+V7O|5$9R)<m$V6(8y*q^5n02######}V653>+?&1,#E9O0x.^5O.,#@Q$OE.]m+###4m$jAJ###4##eFE<$(j$)`%(K6O(-&}5O],%1,#mw$px5######J@$"},
{827.59,186.78,4.33,6.21,"95#$##WG#8xDImT###fG#YIEUy^######i?#|C8######)##I$'1u#iS+iY?xaEkY#Uz6J&/oz^$##$l#y>#yh>######<##%8,'Q$@oV16%fz<RG#Mo+]x)`x^###_G#>C%Jq;######9##1Z%;H#{{?7v&|<GRG#tY#-%%?x^######kA#DI+0u####9##"},
{981.37,196.75,3.90,2.54,"A*2BB0}>&H##{c$/+/t?V###+bB{Z%A7->##4&.A,$###i>#OZ&=[$L07Wc#x$*th%t?V###*AVsZ$^Q(p##Y83$##tG$n$&W#%t]$C4JJ5#6C9o($t?V6##<@VEl#I#%}$#Io3D?$uH)=6$Y,#uW+UB7.,#5M>n8#)::rd#A@VSl%###=%#)[(+0*j>%&##"},
{981.37,196.75,3.90,6.06,"xY$$##a_9/.(QJ1N,#*pYY>#AnY(##],%p[#&^7%##dv&5A+xc(1,#G<<y?&x%1###hsY;l#2oY###e6'>Z#3jA8l#;m)+##Vv((6%~x2)##CZ&/##[sY%Q$[nY###ne).R%MjF5,#{-'S].e>$6Z$;H&jY#.,#>,#*+>U,%,JW###K&*@7*hoY###^6#at?"},
{807.61,197.53,4.49,2.96,".,#F$#$vS$##%R)tY#tm,=,#8d)95#'##90.######=R$dg8oG#>C&ViC###EK2jR'-w,D5#jvSLu$'Q$)X0.@+W>#|JNVi<o%-@.&;,N2,#CM;^?&0/,m>$C{Sn;@Dl$)Q#Te*PzSm*GA5#,v&8K.[jG%##Dv&#%IwP$/,#j8&RwSE>####>6#<xS95####"},
{835.66,206.25,4.28,0.11,"'8,rP$QK.^?'pi~###Sc$2v$C%R######-##'.'{k#######5^6C,$0I&W7)&f~.,#Q,#FU'Cf~######h##cQ%I6'######pR,aP#JI)0-&gh~;5#qb#I##a&ZIZ#eY#5##(l#`Y2A,$%##=l$,##<.,###sh~$##?c%###IoTFn#|e2###Kw-ZX+}_AF##"},
{765.67,207.14,4.49,6.04,"(Q$Q&%IZR###EHMD-&'v&e[(3J/=Z$mU'jYIOG#'##&x$.~Ru?*88#IZR*##-[RV#$^Q([##'uK:%)eZ&*R%[ZRQ>#HH%t8H:@,i%#QlQ%[#vZR+##L5$M%#hcL/,#>5#&Q#N^R}&4Y>$~5#xY$W##p/4rK'fn0###$##{~$JU7.,####3##bw&;&2######"},
{143.42,214.28,4.67,5.23,"W6)*6$yc'nG#E;?D>####($#w%-fw*[v'U$(pH)D''<L6oc$4-(&##>K%D:9WNB###=##n&/smU+##41&(/FoI.3##xN,{|>###(##^U'Oy8>A0}Y#[7&)`>qrUCc%-I$-J)3|5z?(9VO*.*###sI$Xe.95#?5#/C*=$(ZP#_*-eB4:5#*##si*K(7go4N>#"},
{644.14,249.37,4.82,1.49,"s/5###9##^u$.yb###p..n,#Ce/&##ID`{,$qT7###]%*Lm%I97###J,#~u$8yb###)A/P,#;g8,##xzb0-%9g8$##I@(Ee(6'7.,#&##<,#8yb%##DR*Ll#i090##M&Q[B4~T6*##ye,m^6w84D>####$##Hyb###S>#rl$~V?ZP#j%%U5L(c$iY#Dd#awX"},
{756.03,271.67,4.81,3.28,"2o1&##_#&$##;g7###~.(vG$ep`###P~(4H%op`######f-#e~1###x>$?I*]GO###Tm&jI'{o[###yM6cl%Bs`###95#,$#CZ&###)d##3AB}H*##k-)vv'A+EV>#vHR=5#Qt`/,#4l$###.,#&##2J,_l&,e+jP#i6)mY#MB25R$@WC###/s`sY#;c%###"},
{756.03,271.67,4.81,6.19,"###k$#5qg###h&5G##M%Y<.$c$+)##(S-U-$1J/0,#N5$;##95#z$#'qg###Q+Kh5#_7ZVH#&e,1%&cM@I##OA(O/1j>%%##RG#&0#'qg###c3?Rd%hGP*##O-'k]'2NC$##/m#)K1SH(###@5##S#'qg###`:5r6%8rA###u#&kc#s:>###5H&:,#SH(2##"},
{675.75,282.14,4.62,4.54,"WB5######$##A9a6##z,'8##1z9un$sPP'##0T0WR'[?),##Q_<95#######P9aU,#(n-+##';;~/%,9a###83>_[&=6(%##cr@L5$###%##M9a_5#:/24##Ef26.$D9anG#.q98H$QZ&Rn'1);A,$######x9a/##%A/0,#JS0Q##|6Rm.(:@,$##'l#S+0"},
{660.36,283.84,4.63,4.54,"N07######)##voa4##~#&^##C;>Ux)J>O/##Yo-,L6_Q(0##Br?######%##+paD##4[*?###h9#K%Q~~&##H:7~%)|Z)0##x{A.,#######CpaX##281&##1:6no%noa###Ss@U[&*6'(##cV>pb####%##,pa^5#:/2+##N%-l-$w/`V>#lU9,?$g5%{Q%"},
{641.28,285.55,4.71,4.60,"FC;######?##=^b0u####i?#hr>pD?.,#tI*0$#7{@cP#J7-63D######%##L^bU#$tc(###J;:xuDy*I&##9^+cIW1Z%)##)|B#########Q^bD##A7-###qB63w$Yo`&##CM:d7)3-(-##P:9#########z^bE##(n-###Y8/0&&X6U###FV91m%5?''##"},
{252.00,316.29,4.33,6.00,"<]4######x##hmT(c$%##do%yC/_x4Y>#e|@pv#px5m##lZSN3F######d##smTL5$X5#MT(M/*('6L/$B_;Q.#NZR@$#/1;(tI######T##JoT6#$'##}o'Nx+@S0A$#6XCN##6.-?$#QFJb4L######s##,mT######wM)pv+eY#O##nnT+##(c$*$#o,R"},
{407.66,345.08,3.89,4.61,"4K_###g#&'##Qx3E5#0ZG:H$3K_;5#2-&'.$~&4###eu&+##VK_?5#t5&###518-e%+dRE5#IM_},%=-(c5#R^9###DH'<,#7K_###x>%I5#jy9C,#yII'%)3K_)###@'Sp(I97###nu&m,#zJ_###'##c>#5U9###Iw$GU5I..###Y,#HW3@H'###y5&p5#"},
{291.93,371.37,4.98,4.60,"9@,######tG#V+L######[n#V$U######fS&BmU###$##;9)+y6######K##emU6#$###[z)h<C,I*###@v:{mUpb#:5#NJ$)^7######~G#MpUJ[+###;Q#nRKvED###L$%#oUJ#%;5#f,#{70######<Z#@mS+u#%##h0*WaG`5%<##%`1?mU###0,#d6#"},
{760.28,371.99,4.66,2.11,"%5/G,$###0,#DNQeY#######Ch.uR.######F##7J0######b]18Z$95#5n&{D].,####Yw#K0]xY$###;##k5#Yl&/,####ll'%##*##4I9DB]%##&##f:)9B]fY#/,#J##[G#|Y$~P####8H&:5#j##FyXaz=$##3###.<w^:OG####;$#VG#vP$######"},
{556.92,456.21,4.33,6.27,"###SA-|k#95#_##~4B######=o'aw/###%##4I*OG####b,#####|0GC;###(p+I_P>::Y,$$NTJ95E,$'d#/U8######4%####b-#lJTMw.jv*bm#'s=QMT5IT`P#E#$=Q8$x1######M]####fH'F/'~h>v#%on0b>#OLT.v%2^8$##%N.Q@(h&5###I6#"},
{852.18,471.92,4.46,4.11,"9J#KgM1I+###[t/n[-######rS(6#$######]c$?u$######qq0W[@*.,w,%d+[/v(###c5#XHL######6##p-*=5#.,#,##-z:i5#ou&Zf'x%[######gV#bZS######m$#Nc&###'##Ll#981$##/Q#0J>|{D######dD$>%.######lI#+u####&##Xc$"},
{568.43,490.77,4.05,1.85,"$##Dl#>R+###D>#L##t|DH6'o99(##*w)s+1#D>$$#W'7E)%,##b,%=m)###+A0&##*C6;-%k@Y###L>#?`&OtK,##kRDO'']##<Q&X,%###Ue.SG#pn.)##NFYfY#;#$6##'WQT5$?_2j5$#~%Gl%.,#0,#t,$I,$L[*###0a-|Y$J6'$##6@;n,&6Z%###"},
{951.62,507.36,4.90,5.80,">m&Hl%###'##4>CQu%###C,#ni>PG#5-$%g1eY#=##]^/q-+g>#0OC95####G+6M7V95#+##F<Vz%+o[,@Q$g?(-i'f`B;#$aP#k.-kc$]6({bM#&.%v$8K+T7V:S$*WBSA#+w(0_'5kK$##R5$-c$sG#Q-'q}M###r5#)U+i6V4$#v6V[J#<?&>6#c|GL>#"},
{585.93,518.61,4.38,2.67,"eY#+f$Jm*###ZV>X&,pb#*&#^S1fR)95#_s.&##iNA(##Hr<9m)Q@#]3F1##^6T=Z#:-(}&#<q<by)RaD9.#%##)RGPw+U,%&m('##x=BO,$E9T###?-&TZ$rPH,##zY88~..,#9,#=M('lOD>####ZK+,Q%SHL###9H#-.(YpKj>%]##qm+y-%h%0T@#(jE"},
{434.89,529.00,4.08,4.54,"^h4P?'eY#2##7zX###rG#c6&Iv)I$#kK17/0X5$dR#H93mP$}<G$##1u#u5#SwX###jl$^y%#&17-#|(OO.(Ml%+I#Gi0O#%E<C###Pu$bZ$czXi,&sY#~Q#}U4}02*,8-Z$A,$VR#A]KgG$0R*R5$gG#)h.c#IP;=###~[#l%&iyXYl&###Ll%D:)SS1,##"},
{765.40,532.94,4.72,5.88,"t?'d#&95####RF;H?(###X>#MD7|k#MQ${/1fY#N##wg4uQ)7##`D7+u####Es2].UZP#(##x2UIJ+pR-=l#1[(Yi)niCPG####ud&C8*xY$$ZN38,V7(+J&0.U_e$a)C18#F@)JV']>P5,####@##>/,Au$6OH###Q-$Xy'K-U8$#V-U;]#Cu$L?#y3Hs>$"},
{572.88,547.39,4.16,1.05,"O[$GlFj5$U..T;2Vv*pV(,/1HZ&###Em:M[*F6(pY#k>$%c#613my.(n-$##R:WBu$)o/1##V@-5##i9W&##PG#WZ#]A11u#SNCqv$|Z)=$#77W`##vn1`&#w6+&&#/7W3##8Q&X##5o2ow%Qq=Z#$jY#EW(-WB7/$-06g'#Ku$mT#H,P,##2`<yP#]#&gR#"},
{572.88,547.39,4.16,4.22,"Oc&$7#6W<,Z#hcS,##3c$PK#0B5n'#]|GVw#1u#Nr'+2?Wu#mA3:J%ll'R##YdV1##_m+b%#@A11&#gdVU##`$+|##9EBgm$S81?#$E>#KQ#NgV%##PI,1##J&0+##@hVL#%xd-###JC2710]5$rY#vl'$c#jc8Xd*=Q&###0W(tJ3w:0X-*Y#$.80CR$/uE"},
{343.51,555.26,4.36,2.70,"{b#,0,(c$###Qy61o.###'J#t?*v?(D>#1Y2mG#_M>%##mFAl$+R%%hU;/##/ITs,%Ml%eK#`'9$T({<E3/$3##$7I68/=5#|6+)##l|>7u#WLT###1?%k#$/#E,##o$B4m(###K5#pN.V1=eY####CU*7Z%{vM###'6#[m'%2QZP#E?#jd+e@'H?(;8#=lO"},
{282.40,558.96,4.17,0.29,"i6*;@$-H&0##K,P,l####Yf#B]4>S*###ri+xP$][J###TJ(7C8sl#1I+/##>mS`>#j>%g%#T97S|-<V>2$####HDRQc&:5#WC8yY#oS02,#fpSF,#zc(4##}w,VI#MpS`m+###/I'V&+e#PF$)###p@'NQ&GSJ###N>#Ic$iM+z,'57%|@/Z6$EI,i##kmS"},
{784.35,571.85,4.62,0.15,"%j>$##S>#&$'y|A###]25kY#OG#A,#B)7XG####[##&04###?'6###wY#oZ&N&Z###2y*R-&H?(/##8qKk,%###*##Vx+4H&Iw.###&c#)-$<&Z###4-$4@&tMA###m{RDA,eY#C##jt=Bl$Zl&###E5#^6%~(Z&Q%'##ul#*@Ep26S(ZYZ%+##Zn%a&Z###"},
{784.35,571.85,4.62,2.70,"bv)$l####V~'a%Z<#$###1)&^OHt/&YJQYP7###Jw%p>HmZ'By4C,$###d5#H&Z'[&+u#Z-#EkF<F,8}J5##&##U_*yS3###381/,####&w$$&ZWw%X-*%R#hm*NW'[`D$##>l$_I'_5%'##.]3;u#rb#$?#.}ISR#Ny8}##El#^o%nU<###$['jG$###5##"},
{912.29,575.24,4.94,6.26,"t,%?n..,#)##ZhfU,%/,#*/#BpcG>#sP$#K#L#%eP#6d(LZ$_/0~^:###(##5if~#&###-&#<hf.m'D>#R$#,~+PZ&2u#y,$p07L5$###B5#9kf######L##PB[ZQ%ZP#3##.o,,7)D>#%##G%.######9l#?hf######x##zRXZZ%95#0##xb#58(LZ&95#"},
{952.39,594.41,4.71,5.90,"#.+Nc%###3##)D8U6)###,Q#VB/5Q%.u#n7(AATGQ%O5$D8%>l$7T0eY#(##^p/9vE9Q&(##QDTW0/%x1K,#Jl7x6B0p7&#####p?'J7'0v(FNC@7&M)9%v#w?Tkv#[?Te@#+~'@`1f?T0########Ge'R6)Pq;3,#qE=BZ$FATHw'MkKb,#?R+DQ#i@TNH#"},
{832.68,601.48,4.65,3.07,"C&*T5$O5$/,#/tYfY#.,#@##xINW>$###D##xG%N>#VJ/|G$B@,~G#LZ&J,#IoY[P#E>#L&#{oYCZ&###N$#=A/.,#y[)86%@~/-u#iY#B6#prYSH(###V$#f;Q}T8###O$#)#IZP#G5#^@$=~,OG#O>#;?%FdKL5$###g$#ToYW>$###$)#z07tP$95#2/#"},
{784.18,611.74,4.44,4.42,"sQ%Y-*###Ql$kJ-T,%###O%(^S0###B5#QJ(:5#$##Cv&(R)<x$gsH%##,l#ta0N)B###zQ%*l?*6'bP#2).=#$###&]'2vNKQ$p-+^6#+i@RfT2-(+##x04HiT7m)$##^O0U'-2sBrx%]fT^P####u,#Yq=X@,ZP#+##{XD1R<?o3###T-&Zf$ldT'l#/,#"},
{476.41,615.58,4.06,2.68,"zb#By)8Q&###@r>Hf/.,#=%##]2X6(###v`+p>$2`?%##@|9<I+IS$aiC1##7@U]l$(Q%W&#={@UB*_97N.#@5#4SH6[)sb#b$++##vkDn>$0CU###~u$*6$NuFE##_X:`R,95#k5##C+H6Tpb####nL-r>%JHI###2?#o?'H:M{k#P##RI+o&):$){##L@U"},
{420.25,620.98,4.68,0.28,"W,%jH#J-):#$<V>1Z#|k#.:%%&1;I)###UQ8%l#b*4#?%1t<-/0[v#i&5(##(nV]#$W>$@&#PD?468Q%/?%#F>#spKM5$B5#u~1Y>#>18W5$tpV3c#dl&4c#0y/@]%Y@Lf@.:5#{w(<n(@oVTQ'###R?%c6(1+>A,$7,#,m%D;(&~.F6$rd,*w$XC:B##:bE"},
{356.62,631.88,4.08,2.59,"mY#ry,~#&###}S1_]1###M@$j,&wc&N#%*O76?$x&4TG#XV7ZR,Se%4{?:,#PvS2H%*H%k'%881i$%2wSix)o5#(PBLB5=5#?R+'##qX>0Z$'yS###.6%fZ$J}D(##CpJpv*###^G#Bk2q(?pb####UU.CZ%BwS###q>#_I&cySpb#f?#+]0V@(?6(AL#KvS"},
{954.43,630.03,4.04,2.68,"aA2/K&=n.yQ#MhWS&*1v(>5#RfWr,&######U>#*]/######xn1RR#1~O1R%VeW<v$VK2y$$xeWST2###7##nY#b:5###$##If48##liWbK.;fW|P$SN=X6$r<;51:###(##n7,QR+###+##.q;PJ#-eWxH$X-HN-$0z;$##[z,=6(######ld'kT6######"},
{966.37,634.03,4.78,2.68,"wP$`W*Bv)(##4~'/.LL#%'Z#%c7lR-ZP####s$(M5$######Do3rh%;y7#$#yz;sU(ap90[#61V^['/H&2##Co0$J,######@dUsA$GS/Z%%mz>H?#RyPC~(w.VwQ((7)$d#l[*uT3###$##U.V{>#6[)R%%-)?f?#m/V.e',6GDH&K@,;,#<e(p-+######"},
{771.30,666.03,4.44,2.74,"s~+M:M~#&T5#OUJrB8###F5#9N66w(OG#######p8%mP$###3xRMh-8m)W5#|KXPI)/,##6#LwU,B*D>#.##.,#t~%3l$#Z#RJX([&3l$|-#4KXt,&###V%#Bs=R/.###?##DH'/H$###ZR$AQR$##R>#1{':[Q95#D5#c.&?}CuG%###C##/.+zI-###I##"},
{362.71,670.27,4.43,5.16,"&e&Cu$7{/$##>>3`./S?$:5#%&$.kI)3)Te/Q~-z0706#NN4$H%7,#W3W$##kC6jH)lR*_Z&,n$bXIl6#/f2OH#;vMZQ$X%/l>$E:%t.W###AcMgZ$X..$?$J82####$#%z7###;5#`6#{MBRI*J{&q-W*##c.WIH%+u#:$#SS1###(##B@&######^,#aQ("},
{545.99,679.07,4.28,0.44,"m,&&n#AA1]P#mbNZu#W>$W^$5V=&B,###6W&Z>$44-C6(I7&6x12R#n{C(##7SY6c#;c%i%#H`Bka-TD@;$#/,#rt18v(Nl$t.09,#lE?M$'0UYpG$#6&|$%f'2/J)j~P]-*C5#Kv%Hx+XwN.m':5#zQ%{[)Dp4F7+3,#Le%(S%4|C}>$Ml%$d#2)@A,#cn-"},
{914.23,739.42,4.14,6.07,"D/#W7Y<5#95#aN))]3###)##z6'B5#fY#dP####W,#k>%.,#@U*&GI5fM2g3B=Ydw/sZ(*-#WPGP>#L5$,#####7##AQ&OG#6S,L0)d9Y)$B#8Y1,#[v)T}/,OF<5#D>#C6#PG#%##hY#7H%l,&?u#j,%}@Iwu'###$##ns/$w,######{>#_5%######Hc#"},
{669.94,777.57,4.61,5.91,"8h#D~[%##TG#$s'u]7###=5#+m$<c#M?(:5####/$#7m)###-30u}G;#Hsx-Qb[P/2PQ'>I#R4Ba?%'$(-#####/-#<@,95#U7-D7%pkEy`[z~[K>#t,&Cv7e1:Gd(eY#Z##E>#?l#Bu$(Q$_Z&.,#Q,#fPHY$*######<y+i?(D>####I##U,%######jc#"},
{370.59,809.71,4.70,4.25,")%#q&Z95####h1')g7######>Z$|P$###^P#.,#xb#@5#tb#/L+7@K~`/f'8G+Z<n.*6$AZ$82=###1##%['######5##*v'bp1k#D{?C2g2r%Z?c#<S,vL$+aF###H##@I%######G##2H&_e-vp7*R%xn@EGN|b#D>#xM%E-)$##m>#PR&######j5#uG%"},
{413.22,817.96,4.73,5.81,"o_#8eZ{Y#Ol%?E'cJ2gY#/,#HZ#x6%},'######'I#xY$###wD-A)7DCX):4@kZ17+[l&&w#{M9XR&_5%-#####VQ#3l$wY#sQ(K&'CYDUiZdeZ3,#7l$VQ6GC9}>$ZP#r,#mP$/##D>#2m#o#'$##P,#h{;-e-######J8(6A/######J##[?)######G-#"},
{879.72,824.71,4.42,5.87,"###yg(]B7###gA0l=:+K3hP#HL4GK5N>#JS,X@#`f54Q#,g5###a##nX?.,#^'9Q5#+>6[m(CSWpb#K5#')1p[*5['SL/9kD-##k6(6/+D>#iH)3Q#O>AtZ$eeSJ#%`l%r`5V%)vT/48,E]38Q#D&.Tl%<c%q<Fq-%[%/ep$l+L###-Z#,XWfY#@,#-k/lRW"},
{884.95,841.30,4.13,0.22,"y$(1,#[$%gS1D<0-H&6##1?&'f)eS/Kl%&n&K>#cG#*.(&dK.,#$##SI$VB7-L6.,#j/*Sm($KZD>#>5#'},YR+$##&A)wOZ$##GI'#R(xY$Cf0qZ&Ep4qY#wOZ-d)zP$_d#{t1+}D:n-Z@&%##/v$hm+7H%U>#Lu$tZ%{.-.($yH*K[#8w-=)$4@TO>#>H&"},
{920.09,852.50,4.32,3.42,"###,##7Z%[P#yY$:-#;$)7##xu'7-$)6$u+?9##2@*K9$QjF%##J>#9d(.,#zq?QH#`p8^,#4x[aG#ju#;x@=I+`l#+t,)y[Q##GZ&t,&###th>yb#AC5N#$u}[dg80l#@c#Oo,|zM#].3u#)?$cG$.,#$##r,%J,$%R(%##6(*p<F^6'/,#l-&*'UL>#M>#"},
{637.17,866.54,4.49,6.18,"###V0%n)C###Vf-^h3Yy7J>#g//z2@2,#g'0OG#~8.An%(*@###QH#2GE###kT6,##Nk6lu%g8Y|k#A,#=9*@n+c$'`;9mT2###*##qJ,{k#3-(###us3YQ%]$R+u#_c$La18R(T$'xS,=&W###5##H$'+u#0o1*Q#hp6a$$S8Y/,#4c$#:I<$(^>#^f*-<Y"},
{880.25,866.25,4.68,2.19,"K5#*-%AZ^1U^|/+Tv'oW?ZP#PR)eY####7,##c#L5$:,#.,#MZ#%H2lU^qv*E%UKA$&T^N##=L8&?%.,#@,#'##<m(.,####??$AZ^AEDs?&C<>3.%@/1[R'DT*<A0###(##iP#C7,###0,#DJ0yH([>#%&'{>$%l#F.&a95Xm&_5%6,#)(6a>$eY#+##~m)"},
{673.39,876.05,4.39,2.11,"sb#:5#$l#9,#iG$46%OG#5#####0/&MH'S7,5##.'2.c#siB1,#y$(OG####U%.XU-D>#?###/Y0m%3c$yUKh-)J>#jJ&[jZ*Z$zFE###J,#*h:4f*95#Q.&@kZ{Y$###de%[+/yU5g~.M8-C'6Tm*0,#O@$d[+sG$$6%?_/S1(sS.9-'Pc%,V#-gZ6#$)##"},
{676.63,898.07,4.74,3.55,"9#$]5$aZ%<v'Hl%qc#]l&^H%du&Zv%<u#>G;4,#^$(>K$B8Y###$##A?&9m($_<I?#-D<SQ$[JZ@l#yP$pv:LR+$S%pk1HLZ1,#.,#-Q$#H%c'8%##~L/JZ%wOZ<|Bf,$v,$K&+zaY|-+YG#X5#L5$%##:5#*c#&l#Wv&[P#:~%AV=+f,/,#V7'gtF&Q$5,#"},
{980.30,898.84,4.19,1.03,"95#H$#skM?#$$m($7#Mg80),}G%|Y#t022HDaJ+vH(m?(r,&NQ'>$#3mTL5#hmTrG#Hy5_D'<%-Dc#bpTj/*D%+O#$o'6TG#P?'a5#.nT0,#TrT?[(t/3>-$29-L]*lnT`,$Ol%(##-)6u,&rb#-##PPB#-&Sq9mY#$o*l+3.m'du#*6FC0.###Q##{94nP$"},
{980.30,898.84,4.19,4.09,"C025l$###T##)vH]n+DZ%wu##f+NZ8RU7vY#sOARm'K#%3##nq74Q%KH')##[%UIu#@o,.'*]A1:[$l)U9I(5%UI>#Z6(t,#H:9I>#;%+v5${&Ur9,h[,Yc#t]5lM'r$U/H#D$UK5#s#':%#:R(oQ([T/[H'00/-.IqP$Nu#YU:=z*Yl&U.#F+JuP$###b$#"},
{663.31,910.92,4.13,4.49,"*$(/##4?&q#%CZ&e##f~1Dl#>#$Lm%H.-Qb7Jc$K#%D6#k2Y^,%N>#cu%f#%/x1>##7o29$$)/Y%##f5%sY3T~/;5#^_5{|Xvl&tY#r,&S>#;022,#[]1E5#64YK..]#%v##$G0:ROjK7uG#zc'ZG#e>$fc%(6&/,#z[&OZ%#C$t$,^[$aI-NC$IvP_P#;$("},
{426.44,915.84,4.45,0.06,"###*6$7O5Ux3/y6)##'o*D(*.NB######{.#[$&+u####5,#eG$I8&(QKYu%lYN;##'//NS$4o2###4,#)A&M>#.,#6,#Ol%|k#`#$X=/<8X2ZK%##vF0xJ2O`?###@c$-Q$xY$%##>c$BQ%,l#]7#%=XO8Xj]5j~#Q<X7-(XXGTG#a5$d5$tb#gY#~>#86&"},
{674.03,919.45,4.54,1.01,"Gm#Bo2qQ2B/20z)J.-%-%###xP#6S*j>%###-Z#oc'6#$%##v5%te$u.e@I+m,H'6$/~L~>#Z[)3f+nS1d5$+].uZ(iY#Wc#0Q%hRbH(e9##Sg9<J$6<C?$#-@+z,$pv*8$%H[+###1,#Wv#h..X5ESG#tJ&h-*<%(FZ&H-#s%.PQ&gG$Xc#4-(######I##"},
{674.03,919.45,4.54,4.20,"###S,#S6)###.c$CZ#[e-J6'r#'E$#g[,dR'%l#+/&r%/VxM1,#2d#fd,###gm*lc${$,o#$R4I,$#-q;0A$uMf1##P#%B,X^P#~c#6]/1$(}J2g5$r$)~A+[?FU>#o,F>Q$gSfV..=H%=@${k#%##Tu#3$(Fl%###Q5#_n*i5$###U'({R-i},k7/Fd#+g6"},
{432.23,937.41,4.43,0.04,"sn0I##dV-Ep7fGK###e>#]Z$x5&%##`#%,H%[P####;,#}5&*?&q-#X.d5[URJYc5#`);-Q$Fn.:5#:-%6v&6#$1,#_>#I$'lu&s3&N(dg#%:'`/.$j09W6##x.^#&/##qG$&##yu%<5#(l#mP$7,-4WBu&,nw04##V,%$L&1Q%#Q$[P#L5#$##rc%.,####"},
{717.29,997.29,4.31,4.51,"Qm*######+##q&T6?'###:I$xHG+^7.,#L@')%TZP####0##Ge/######I##<%TxY$1##xj0SRRsQ)vc#N[GH$T###$##zc#_&3######(##P(TX,%6,#I$%FnNzI-Q?%#J*O$T.,#$##>H#)d(#########*BM;c%###Q,#IpRO@-###/d#]$T######F##"},
{598.04,999.41,4.81,4.53,"^#JuG%$##K$?+.Q%Q%###2J@X{A`P#95#x##$##X>#EZ&####NW?d*###x>#mNWqQ)###[c#{EESG#;5#?,#$##tG#Bc%.,#cIW######A%#NIW######f&#KOI###E5#;6$###'##L6%{Y$]9WZP####i%$rKWOG####|w#Fz<###1,#>Z#######@l#(Q%"},
{671.84,1004.47,4.27,4.54,",%,######%##VDSmP$###/H#MSC'$(###Dl#QjC######0##Kw.######:##P@SfY#5,#TM)`AS(Z$Bl#4_1w?S######F##p./######*##zCS3v(0,#Av#,BSoS1:u#97%F?S######V##:.-######:##K@SnP$###W^&_?O.$(.,#Y:*f?S######G##"},
{316.13,1018.42,4.59,1.70,"D>#UI#NV>SG#IU+;o*{R.E/':L*n05?-'zI*k>$i-(Qw,.,#~#&=##'xRT#$'vR1##(p2_O16@+{5$.pQ~i8,l#1K.#~-J>#cl&0,#OxR4,#xzR<5#OL6ec#R)7au$5zR7l#~P#8l#i(6r6'd,$H>#_}B###Xa2@6(Pi7###}M)ukH)q60,#x,&?['oA0F?%"},
{582.92,1019.57,4.42,4.54,"######8#$=5#Qu%###<$(#J,$IV###RG#$B*&IV###$##O##@d*######e5#l3F95#4,#EoDytLeY####^NV,IV.,####r%$V82######)##UMV/v(###^$%=(Mg7/###pe*N+G.,####h5#).,#########'JV{k####M##HJVW>$###A$#oGO######H,#"},
{582.92,1019.57,4.42,5.65,"5_#%-S###.,#MO/q}M###'##~h)i7/.,#$##[##>$)QG####wA-=6(9##T?&H{Tr-+###q.#:?B:.N###Y##d0#TvT######K[).,#*##'Q$L&E^z=###Q##*^'MyT###(##PV%Or@######)###########t##qQ)######s$#`K7######4f#LK6######"},
{300.46,1021.49,4.31,0.97,"fG$(9#gB8%##Ro4'&#(K4G(&pm,P$#`x.{pN46&Dm$)h3sS1GH''B#oOK-##KvQW-$jI.o'#Uz<@q$}PP&J$f>#8)-[A2QG#3$'cQ#+vQ$##dzQi~-R7-T5#RA*cj0FaGl5%J>#=Z#[%-~%/G>####SX?7?&[A,e5%Y.)u|9ic%rd)+%*GJ-H>#;~(+6&I6("},
{259.55,1070.26,4.30,1.18,"Y,$LL*w7UB#$-u#Q[$az3t6)%]2yP#2%+<)+2?&EJ&bGML?%'##DF;OU3.,#bn-5X-Z;8*H%v7UB7'Z(=ov$DZ&9%#g6UtY#.##5:U#l#eY#@PJ+I<yY$Ru#T8U6C,.M>x##N5$d/$V6U$#####W-)9,#KQ'Lx/m5%&##L6'Pq1tN=r5&(#####/n&wu'###"},
{727.78,78.16,5.83,4.65,"P,$OG#######G_=#########yx`###rEF###A,$###7y`###2Z%######0##m{C######9##jx`###NkK*##>u$###Wy`###8Q$[H(###^G#j;A95####]##fx`###OjF=##uG%###iy`$##&e']@-$##1u#IFDD>####%##qx`###{^;%##%Q%###Ry`###"},
{996.69,88.12,4.99,3.73,"jdS:n-J>#vC$ye0/#:J;?J:%E>#9d#g|CbC3######,?%_U29E90eSmJ1`d&#gS=fS0Z%{5#X96&$%AQ&^]*###<##ZA1XZ&.80fI%^dSqQ%`dS@l$UH(#%#iWAMl$(c$;#####O5#07+'?%5J0.B$@o3#&#hcS%##pb#B(#$80$##95#FR$###$##]P#[W."},
{1007.75,98.33,5.36,2.17,"Ar2dc'(R#F^7M,#_m).1(e)C.,#KQ$4f+x&4L['9R(Oc$vu&o^%R:<w.#We04f*ndROR';80ldR0w+?9(]vI.R))l#V1$?dR'##W>$~'#KZRCv(hH)D$#odR&gR1v(cQ#LlE.p*yWB?;'K]4$##95#k?$37,n>$J#%D,#wB6s,$(v'n8+vv+g,#Ri=#'/.,#"},
{1007.75,98.33,5.36,3.72,"UZ&>R'D[(Nc$UvR}H)6u#c:$L]4A'+or@A;'.,#m,#gW?ix.194.,#?H$Co+U$I{vRQ7,ZK(.QE7yRvc(tZ#IR+(6$4$(Wo+<rA]5#GI*E:(PJ0>f*{vRQR'$wRN$)<$)Q$#S'6o>$U,%Q5#Ay5z)4YZ'uH#bn0h^%AL:2/#g5Q&##A,$e'#Wd+$##95#|Q$"},
{1076.47,99.41,5.43,0.88,"(##RC5OG####W-#UU:+##.,#DQ#t?*Y~$D,$/,#:5#&M%6l$-##-qQE>#@v(/V)nlQ0##)n,BpQeT6(c#J0*ul'@[$w|7tT195#Kn'/d'd?Gk~.Z&0bR&{<DRqQ4vO&c#,~&Q@(NLN9S.eu$fu&?6$QZNq~*8,#[2;O^8E>#r[&@PL;5#(##'C3%T2###6##"},
{1076.47,99.41,5.43,6.08,"j##[B6T9&`/4YD)e]6C7%pb#~l:>u$###$##X-(U>#Fl%(##Ov#+|8G9+z,'J_5DZ%*l=Td(cJP###F#$*J($o0-?#kI.K5#pU.$8+M5$$##`mJ6c$N[)ol;aW?U-(in(fLP%$'Ew#:3A9$(Ln+'##J#$Q.+NgMz,'###57%&s-3IPq$)Vm(+##_r,}WED>#"},
{14.37,130.58,5.42,1.25,"-S/######I##7T`###H>#SL#H|FF$#DWBgd####5$#M/3###N`C######U##?U`lY#W?&X*'TU9ou#_W`.B)###k,#ZU:###_'8######&##{Z`rm*#H%n5#Cq0]j2iT`U>####@[&Cy7.,#############3?%-u#######eG#,[%j>%######f[)W>$###"},
{798.14,171.44,5.50,1.30,"%##.LY+u####Jg-hJY######;=7^IY###Yu$UJ%GYM$##26&:5##8+q#'8c$N=AkQ(-##2m'TLYwb#1,#Xe(}6([Z&$Q#kNAuQ(nY#wG%bl$n1<Hl%;##JJ+|JYOG#$##02*g;=&n+Nw&gLYZl$oP$e>$K#$Y::xu'###mw#~KYr#'RG#r&$-h0,mP1Q%G#$"},
{49.71,179.37,5.69,3.65,"NE@*R*###Ei(@[(a<<=I)5ABQG#W[%'L0H18######Z1#E81xX:U,NQ.-^-$v(Ub%U6#$g##qA21I'pP$MR%&##8l#6m#'H%a97t-%x5O%R$4%UQ5$Gl%`%#hg8-Z$95#([%F>#M>#2,#qd*o.0xI#fd,)h$1uO###95#pC#Jm*######fT)A,$0,####Z47"},
{88.14,182.26,5.04,1.04,"###/(#H`C###^x4/(#=L:X.#8q;zw)<l$U3.F>#^L(:f/X80L5$q&#RlQ%##glQ`6#jp:Y%#0V<oM,AB65$####Vr*&@+###P5$9##dnQ###_fP<c$S_;2##=],>%'NFH&#####r,#g;>#########[~E+u#4&+ZP#~341_7nR(5?%E]/uT2###Y##9g4iY#"},
{115.81,182.17,5.67,2.06,"95#G9#v>Q5,#yL/Q:.#'6p[$0x(Gy0#B3yH&F5#1d%/y62,#_5%_##@@QW>#1?QM,#AD;_2*_v*>,#owN1B-95#>##wq<mY#4l$###vBQ>5#8CQyP$xi;,Q#)C-P?&)AQ+l#######D_0?d*######dmC###-x&D>#DQ:D>#|L%3A0fI&A,$Mu#fA14w%VH("},
{976.24,182.18,5.63,2.72,"y,D`7,Qc&W-#q7)u2-t-W(##%0Wzl$zd-nG#Jv)MZ#pm+Q#$8X>F@)|R.=6#kJ0DV%y-W:,#=0WE[%mw0F,#t%/C]%N'8Q>#q_%ZZJO-(fY#@931b07IV<v$V/WWd(WH(9R#;m'j~(x?*Ov&2~'>m(l-$<S+G7*u7)WI*NA,:#E&H%^P#ZQ#F818#$###z@%"},
{976.24,182.18,5.63,6.16,"###@##B5J(c$3d(NH&AM56R(0~,$x*bv'5C4zH&$I*}?#>EASc&4-$R&/?@),A/Bl#J:V+[(C7V3l#Un,Z233%-%##H:+i$Q7WB7l#kQ&qR%b./%##F<V([&b8V/,#lp0xu$Cq9aG#R`?k,%+2>SG#T5$@-$VH'@u#WZC*l#S7V###%M/C@(q[U$##fo(BRL"},
{65.66,214.10,4.95,4.12,"###,##B7&/&2pP$-l#Ev#,'5?;@###6##2+/<]4@##Ue+/$9###P&%pJ1D>#e1:ve,G,$,Q#@]V$&/###0M#,E?w.AJJ05o$yG$X9&hv+'##=3Asu%###m6#R_V$nOB?&CQ#``/F_VHH&%##>?&I$%Gv)O>#nXDPG####Y,#`y7Zu$bA*AR'3q4JZ%Er5Q5$"},
{855.31,239.52,5.91,6.08,"42@'##X>$V##;q<q[(###b##;l#:aA###bP#@w#zA4%##~P#<}JH5#j>%W##U1<F)+t99F##RuLr<;b.*9c:E9+go0L#C`'5utMK,$bu%AH#Dy6,6#@oLR$(W'NO+A@T-;q3Ox/w=65GI9(/}}MZP#8,#w$$F-)###Mx$9y4Bl#ru&o6#8#N#l#{5%SG#3u@"},
{871.56,237.59,5.32,3.38,"###%/%Oc&###D>#sG=rY#SQ&:R#}5N<H$2H&GS&Tm*<A.###3v(D0$t4L.##)6Nc]*^S+tc=ES+@e*o9NDE>r>HB,$9~+AH#h?)?##>9NX,%o9N|L8vA+5].kI((c9R6NbT-.B5*##-I)o0&######Yp'Nc&^,$M#%-A#O3E6,#q5%Q>#vW=OG####*l#{d&"},
{871.56,237.59,5.32,6.02,"b./jn)95#?##L#$StE%##p,%yR#'q;S>#9Z%q]*uG%######*80,((Z08=##WQNn'2Td'Il:]S)ko-_TNa{<znKQc%XH(-##vR-E##rb?@6'8UNu*A5%(&9/kw+aQ=C>K$1/oPNP,#AH'nT$pb####,8$LA1%Z#NZ&X-#cPNRG#jl%bP#FnJ5v(######&T&"},
{886.95,240.02,5.69,0.22,"?,#tq8%##&l#LA#gsGsY#8#$rz(VZ'######$[$N5$OG####C=E8:/gv'+j4)1.5:56(P}M;6(Pku&:@+A##Sf-bG$###$##yBMsC85f-/A+6L0*%Ae$P)h-i#P9##WR,r^$'8-Fl%###2##@Q$Jl%bm##GI2l#K6'PG#.%CWH(.,####29&E#$?u$###'##"},
{886.95,240.02,5.69,3.34,"###'##F#$J#%$##&9&+$(.,#E>#{mC$c#`Z&Pd#o+J-H$(Q%###2##3A-;c%*.,-q$f#P9##U$PsC.Vy/*%A~8-v%,|TL}L8###$##DK.L5$ed+D##6(PVc&6(P=`;DL.MU4,7(4j4g4G3:/OG####Q$%8#$######qh)VZ'%c#Z>$VJ#PaG&##2u#>,#({8"},
{249.42,293.75,5.35,1.13,"oc'bG$###-##gQT######@'#hQT;@#I97N'#*'59A%Kg9h##O]3{k####=##xRTCl$###('#yGMVD(`K7$.#y]7{H#h7/g%#7T..,####&##]_Q6L9###8##.|+8ST]#&S,#8C-wHS{k#z##L6$#########En#9@,######i%#AXH######A%#ZHT######"},
{759.70,306.14,5.04,3.20,"(d((##.6'###u^RTe/eY#+##2$59T^{k#4##t|0cL;######GQ&S>#c6*###_(9~$)^R*L>#{&Mv;?=h4zG$9Y^;c%%##8,#:Q&J5#Q6)###[sG$##3$&mP#BPM###n1.5d%QT^###;,#?Q#3l$(##&?%26'brC###~G#f6%`U;###6/(Km&QT^###4l#:H#"},
{759.70,306.14,5.04,5.84,"###TN*zMB###?@+[V%MPN###_o1^8(z,')##0H#3')r5&###0,#o:#b$V###:XD)i&EWC<,#'f1LJ$;n.%##xP$zZ#tc(>##:5#L$#.*V#I*5g79?#qZ;gg4*m(4@#zn1uP#>u$3##X>$=Q####ef'Lm;Bv)###78#i(V)m(6#$Z##?J0vG#+u####3,#vG#"},
{748.87,316.31,5.48,2.76,"<$);##3l$@##,fWpb#.,#@U#o|>jZ(###W3+3~+)/1###aZ#Cn.b>#J#%W##ABFFL7d$+e,#DjW>96D>#T##$p/YI-###2##0e.+##wl'z?#8)>JQ%R*C{e$;fW:#$=#$0D#sL<######`$#[?)###g?(6~%SiC###-&/^(&XdW######|D'Bv)###$##(L."},
{334.63,365.90,5.46,4.74,"2U7######@,#,PFh,&###j.&4&ZJ#%%##8~'EHQ###`Q'6,#{z?######f5#Z%Z###oG#z:,TIV###2&(tU36&Z###z,&D5#h:<######.##x'Z###K#$5c#lmSnP#|D756%K&ZH>##?%[>#&J.######,##2(Z###OG#####INEc#v;A###PSS]G#F6($##"},
{681.86,410.19,5.27,1.51,"+u#*##m#%1$'=H$;K.rH(cP#Bf.;d))##$A'_5%###-##p.'*6'9##Fg5o?'IN8iR:e:`7Z$A>`J924J*~'-$[)'##_x)=o.0d)###fx/Yl#jo4yG#W?`56HM:`sl%jA+V/?zd+YH%+91S#%VZ%eY#E5#N>#H[*###7[#?176w,sb#4##J_2>#$3u#F>#^c%"},
{459.18,438.24,5.72,1.49,"###`>#q5%.,#Q,#KV4@c%qb#EA)8K3G5#2c$gY#/,#x5$]P####F?#M&3###Hz3@f?/U`V,$LZ`8a=-x.e[%Nm*A,#?K1k5$###qY#v]4###$K3cZ#{Z`3bBWU`v6'rT.Yn>R[*Lx*0x1)##$##e5%:u####:5#eG$pI$u83<v(TZ$'@(W02_G#1U0]w+OG#"},
{547.95,495.23,5.79,1.22,"Dn)l,&1,#=5#M?':Q#DZ&q>#<R+L##tH(98@F>#9w&oR,+EZsR-2,#xY#iZ%z2BGc#q#'b##:AZ/Q#ce0x{';v(%_&#AZ_.){k####'Z$'v%kNCA5#]m*q>#ZDZX[%^f5K##3d':h$DAZM,$sb####eG$,##b['%?%G81$##/c3n%/yv+###on#in.+cA1&0"},
{978.22,564.54,5.56,2.86,"I>#0,#%I)QG####T##g@P(L7^h>4##X'-M}CEz</##^NA#R%9[)OG#qG$z>$+#MZG#]K+ShPw?*x>#xyL'BTYl&c>#eU/x@*E/-cc'###)##1DT{.)(~-a6$18,b;(`@T|G$###nw$VI,/,#'~+M5$/,#0,#8AT}Z$Nn/&##=-'x9'AcQ###/,#fd%6?':5#"},
{978.22,564.54,5.56,6.14,"FI,.,####C[%+vS###JS,]n%SH(###5&P?Z$[P####oZ%_#&,x1######Ow#?wSl>$I]-cM(2%,dH$C{SF.)###5,#W%*[Z'0',/7);c%f,#CpJBxSHR+m>#t]*C{S&QNL>#EQ%Yl$D-(95#C3;Lm'1x1$##{@(?YD[=I%##LeOCr?~P#0##NA/|k#$##$##"},
{408.65,588.37,5.20,1.90,"1L0F>#^C-6l$jLJcG$LR&tQ'6v658~DR(/,#J%#D9~3),pn0qnN-d%r[+5c#~A.qc%^R)>7*>0%s$VSn$Gn.3$#g9~@d$v~2Z06&/)c,%Iw(nGNCQ$mY#<Q$DS0$##|6#9U7:5#wP#Uw%H5OO^5yB3'l#RH$_}JG>####b##oA/D,$###bP#_G#VR$m%0zY$"},
{408.65,588.37,5.20,3.90,"Qm(kY#JEQ)H%/&S95#a@(i#%rn.-Q$b>$P?(####?#VI,I#%~?)t##},_kH%'PJf>$d}C5v$0x/6:3.,#&#####fh4I#%###yb#{r%g&_Fu#o%-UJ)3p6IH%6c$$23=6%2d)$##dq5ud(-H&%R)>|&JJ1LH#Je.VS%Cl$FK(Tu#El8R$']m*}L6;S.An%h@,"},
{420.43,620.99,4.66,0.25,"yG%rQ#L-)ub#kz>lG#fY#^C%U~0#R(###=v9ub#&N2o,%%,@.]2g$$7x2+##Q7Wu>$eY#L&#*i?7R;k-+<%#RG#A2QbG$N>#qx4[>#w]4i>$$:W[,$Ac%xP#v//RK&0-H&w,.,#L/)g%*$/UdQ(###QZ$Wm(0#Bpb#2,#($%vh()R*aQ$j[,07%r]5M,#7PF"},
{326.93,623.12,4.79,4.58,")0+$M51v(b>$k@L0Z%e%$4~,2Z%;,#:>5d@-3V<]G#nc&<H$ZjFqZ$nT7?##-][G>#J%*in%tI.2Z#^`[(m%(kCl5$Q951,#A;@&##&y1n%)F^[ZR(_Z'>R#_e-^M)W][*##ig6h,#](;/,#D-'mc'I6%vD5zB7]F6###YJ%{5%4w;P6)&##5m'$m$~l&_>#"},
{839.76,634.76,5.73,2.51,"X;#]ZS.,####+>/xL=Dc#+u#yw*oP$,~$3e,GJ.{k#A,#`/)ET,{H*/,#.##At~^Q(###C$#v~WdG$,##},%xw+Mu$$##q5$je0oP#fY#>-#ro~######8'#KvO7Z%.,#`##fZ$H9*pb####cf5&##.,#v:$I5N######B)#S@-o>#n6*wm$~P#%h$C]4fP#"},
{844.74,670.91,5.31,1.53,"a^]nY#mP$6-#Bh8mS*E&3)##DqZdc&{k#%##Z97$##_>$A,#m~](##%v'Hl2o96#w$3nU,v$D_]FZ$cZ'EQ$N=I###8,#Bu#oQ($##|5%fV-s#'9##ki3s,Lz~](##8?$/eFK}J######2?#t~/r5&###L6#,@)z,'-##pL2bu%+u#&##?28Kf4######hl#"},
{892.20,721.64,5.67,4.73,"m/%J0VzG%c5%(<+R^9###=5#$~*Pc&###7##+H%M>#YR,%##%M/Ot?lM50/V:3VuH'k$(d$)HM=1,####aH#xY$Q>#B-(/##xd*BIFu}6{ZJ~.VRw)hc&G>41z;5,#^#$hd$1,#`5$*H$:5#95####2/$2/VPQ'A,$(##mI?H,$C,$Vl#EH$$#####CH#}>&"},
{449.55,740.35,5.34,4.56,"|##@`;~#&###Gy$FEA_5%###/Q#t0&U^:###~P#Y{(_%/GQ%Fc$9^*PbL2,#aqTno0=~/b,#hf148%tYM7d(z>$='*FuJvy6j>%9##M&P<c$GmT5,#Xx/h8&%(9i>#<r?>I&-I)d8%wmTWv'G,$=,#3uCV#%53AD5#6S,28')/-Au#DEBd>#W$'5H$WpToP$"},
{877.67,748.27,5.06,3.02,"%##I(/@H'###lZ#wr@+u####SD8;m)###1##}z?###$##:##c_9NI@7e-^c$&_0y?JpcN&[&k:QrA2LQ'm5#.QM######4##Aj<-A+Ss6JM=4y0`:/;}GanIT6QfP#ic'RE-}*C######E##;5#/,#=x#s5Q0,#;5#c5#=8QI[+###&##$1,_&3######0##"},
{441.87,753.63,5.14,3.69,"?#$76$z70$##$:4V%.OG#'##9N52v(######KH&tY#vG%%##5,#*(*vy:###^A*&{RB,$Fl#&{R(_:###'6#_/0[@*k>%&##D>#]q%qWDmd*)@*i(,{6*rZ;$wR<U3###-W)$7)axR###&##`Q(l-%ER){|9)c$@7%+H%fH>iG$<|1|@*{n+vP#_~K,T,/6'"},
{657.53,782.35,5.49,6.15,"9Q%waC###1##Ee$wQT###*##5&$-/1|Y$>5#8,#g>#rm,###4J.0WTQf3CH$8s1ORT101Qm&0WT3y6O#%6H#(n)al$rQ)$##NB6uc$>#;PvL7^7Y5$|M9QVTjRTgY#Qu$_a/i.-(Z$95#8########5w#,FG.Z$<5#bl#IGGQ7+E>#.,#0S'e?)######S##"},
{631.95,787.58,5.81,3.09,"8?%Cn+,-%}5%`$$ZYL;5#A5#39,]%/I>#L>#QR+###-##_5$%n)G&EJK50,#>h.iKTrL<-c#$NT$19l>%-?##V:######:,#$&-S7'p8Fy3GIB1G%(*{<7MT^IT?#$Ec%DY3if1I#%###A##.,#&##d8$8YL###(##TH%QJTQ5$h5$p>%Ty.4Z#-$(###'##"},
{895.00,813.60,5.84,1.09,"wm)lc#=A1(##U6$S4=>C:06&v8'S^8$##el%g#$q?)'##?c$SR+e/%NK6J5#S'6@w$Xh<eN7Eo/^H'w7'YP@C(;9Q&+@#mgNzP$*x%x|EZP#Mt<>FC?m)k#$o?E]#Kb?'n-&@_QjsEL?$O{-$##5Q$p7Ir5&}u$bJ1$A(S;8gw':.-n5$f:45r*pZQDZ$4-'"},
{895.00,813.60,5.84,1.96,"&##iT*BR+###--&U14r5&*c#P/RK6&###p1,[P?{J+~^391R/Q%)K)PY=xG$>FC-x./R&tL,U7A,^1,A-?].,j-|1RTd+%c#cQ$d#&~a72u#d_1~x1l;7MZ${.*l->#~-(c#<K3Ah7A#$cI&D-&?l$mv)dl%-##`u$745]I-Qd*Bn'yw(H):[B7SG#&##aA("},
{895.00,813.60,5.84,5.91,">l#X@F_5%###h%%5XD,u#+Z$r7$*r@C5#:%,S$%PK6###>-&=K3_],P(6u#%UaB%e-[u#|j;=I%}D@)X6`'64#Dn]6hl%4[$[7.4l#|z1d$&<@OE-'E$&LP=V-(g7'/5>9g4o'8/,#4?%7A&P|D#m%0R):]%FRR7c#uw'QVRJl%9I#O(P`QRAm)<,#im(1H$"},
{372.30,827.97,5.79,3.23,"WH$(f0{Y$?5#)7#h|F######3'()R*######uQ'###%##$##rf/~mAOn.X,$n_0N&S{FG~Z%2)S}S2eu&m,#Ay395####'##|:2z~+Qu<C4HEU1K:07i?2)Sf$SO>#w#'J=1)J-pb#%##o>#.u####sn#8$S/,#/,#m5$7'SD>####cl%{L13,#95#m>$kY#"},
{358.47,829.12,5.05,3.11,"=7#uYNE>####b2*je1######Qo.#########,Q%###E5#SG#|p.K~NTcMul%msXcp7;$)v,#?nR######.##tH(###&##J>#[8-1)1fM?NqTtnXhG#D6(@?7_aFD>####.I#pH$+u#######$##/,#wP#JqXA,$###d#%Om?N,$95#j>$q5$Dc#W>$######"},
{871.89,839.18,5.28,2.34,")oBbJWzl'I~&@-$T_O&JW/H%M2;S[$(OG{b#/@'_m*######hJ%)LW,(9+6&de.{NWhd,'.&Sd(ml$6R&5L4Vm&+6'*##M&.zl%3F1{_@TG#TJSe?%6?'QZ#(^6###>##)T.^,%D>#%##5C42c$Xv%qx,X:90;=jP#a>$na3Qw..,#.##aw(Zl$j>%*##wP$"},
{613.05,843.56,5.35,3.66,"^J$U^:.,####4V5Gl%###<-#tc(######i'%bG$######De&eD)xOJtT-%l#{NWM$*(##TH#G`B######$:(######&##Lf1%%*Gv$*k;7eE(KWOG#:##1qLx6V###$##$x&,u####'##aA1{':ZH%pP$4m9z#'|R'Jv'AgNw].LQ&ub#8-%t-&.,####hl&"},
{892.70,861.02,5.11,5.55,"(##H@)v#%%Q%06%A~,46&_,$5[*AH$;m(jRA,n,vP$N7$xZL######S$&%.+mMAoY#/M7+m$=][@l$_#$gS@.n,un&R%=1^[###$##{u&^5$8(:###~U2IZ%Qb[85JoP#ou${w)Qb[/&0o>$7,#.,#vY#pb#vb####Rc#mZ((7%LJ07##N7-H5#z]1>5#o@+"},
{642.86,905.17,5.83,5.56,"$##Km):5#95#:5#c$'|u'&##|k#Nc#wf5gg3-@+aP#)[&Ea:(##vH(L#%###k]56p+Rw-~5#5B~dQ&CH$7<RVd*nc$=?5~C~w>#gH(:H%eY#sg77l#.U-#v&JG~UjF3Q#K%(~g,JG~Kr8<R)Z>#~P#l?#nl'Gv&ZP#j##uc'.S%4K4V##<m)_c$1=A+##'8,"},
{644.80,921.99,5.03,1.63,"6Z$[l#EkWW]`l@)(Z#+vCI#%J@*tG$5Q%>5#YI*OG####$##p-&G(HR{`)?&Uz8|{'4~W,##_S.sH&8S/%##(M8[>$PG#%##5H$@{LOw*8]1%W@%S%jl&t0,c7.I5#H94K#$QV:/,#wb#yb####1,#QH#[f4{k####Gc#wn-pb####+[%EQ&t6)D>#1##Rc%"},
{884.53,922.71,4.93,1.36,"Mu$S5#/v(%##W#$4&)~C;@5#A&*x?):v(0c#bP#9##7T_:5#H,$G##>]4'##zm(5i,JT_6,#(X_^e,h@.Z5#Dm)e##$T_%##{P#bP#)f1###nZ(S5#JY_2x.;T_G5#qA.]V0tQ)q##%T_,##a,#M#%a>$###$##Nu$F8'97,%Q%F5#g@+~f0OG#T##(T_3,#"},
{884.53,922.71,4.93,4.47,"@x^3,#OG#'$#Tn,<K00Z%I5#UA'QI,$##Nu$b>$###k5#X,%=x^,##]?)K$#Io.M;0Qx^K5#~}^M//yc(W5#T/2###+Z#cP#=x^&##f-*6$#<w-a5#;|^-&-_x^B5#um(C`,l84(##I,$I##Nx^E>#nY#R##S?(>l#aA*-I)x^;4,#>l#5x('$(&##m5%I,#"},
{416.18,947.44,5.24,3.58,"+u#&##D>#Ll#C,$gZ#'Q%T>#fu&D6$~#$1<82,#%$&C8$wM@5l$7##u5&+##$V<.%$4]3O##wJ~i5$_l#:TBmm+dH$5d7)L~95#1,#*-&{P$zq=G5#GM7gP#VP~k<D'Q$Kc#r8,mFY-M9e>$0,#xG$$l#H>#d5$kG$o~,###v]')|C.A+###z@*luH<#$=u#"},
{384.35,951.00,5.26,5.67,"&##Rd*###95#gG$5%+Y,%aG#/H&uP#_v)$cB(v'~P#Q@$?$P###jZ&(Z$OG#?`A#7)ml&s5#;]~tP$T,#8^F5.,&n%iQ6U^~jY####/c#Sl%,W?$##Kd%vQ(nb~j,N6,#TQ$ye*nb~=f0=c$/R)###$##4l#m5%###@##$I)`@%}98t>#|H*A,#-j;;l#|w/"},
{765.95,980.07,5.54,4.27,"Zl&######/##=-P######&9$//X######.`%z^<######X$#LQ&######&##RiR95####K?#JXY######`n#z]7###$##36#F>####.,#$##jSY.,####N%#VSY######6(#>x2######`,#(c$$##.,#+##hRY######w&#iRY######H'#4v($#####;##"},
{349.82,1002.31,5.09,5.24,"no)?A'R]Q95#BwAEI+9Q&###}C(X83######~$%W$*######Ge,iV+cmK>V;.~Q~,$^H&1T,S5B#H%.,#+##P_7=#$<5#E>#9V2?6'O2)&~QO05=@)?Z9[i=Ss?,8/7#$E##0ZL.,####{G#(N*^EB_?'7u#RQ$N365i=9R&I1;Lu$bG$^x$G{@######N$#"},
{643.15,1008.44,5.58,4.54,"hH)######'##pwW######Z$#6xW######e.#np:###)##.Q#5w-######4##kwW######))%xwW######Oi'/z;###$##V?$Cw-######'##~|W######V6#~|W95####)[#-D=######q>#Dd*######-##]xWOG####D8$1yWpb####5%%@NC###/,#I,#"},
{567.58,1016.84,5.91,4.74,"PG####86'###1_[###]u%###/`[######&##_l&'##iY#4,#######k#&;5#v][###FH&#%'q`[######lc#mn/$##z>%0,#cQ(###<5#zP#bcM+u#?5#nP2R^[bG$###k3*df3XG#zY$)##9m(######%##ryO}>&###r5#&tYll'###o>#cg60,#TG#=5#"},
{283.00,1022.30,5.86,0.48,"D>#N5#i[(ey3>Q&32*]#&Y~&hH)_*0ZP#-1(+u#&'+a7'20P=5#-6#(80dG$<I+G:$i,Po5%?.P:D/R?(b1&Mm(pm@a6(t/*cv#Y[&DA1bP#<7+{-%D}Aj6N%14CL089.u2@vP#ML.oR('Q%<1$1&0L5$'##M/*VG=`$*g$(E,#6.ASe.Kl%{>%*=?K,$8m$"},
{283.00,1022.30,5.86,3.95,"Y-)jQ&4l$fL0XsEF?&MQ%Xs0*7+H$&K8*g;;AH'###}&#{`D=7*RQ'J>#/{2nq78?H?8,DT+.%N6j<Kd))c#Ry7.,#J,#AZ%~06dl&}G$pM*G~/(n%]gOs;1VcO5Z$7'3XJ#2B4}6(M5$m5#g^4^`:dl&RC&I7-0E,~OIbo#l':e##,J/81$rl&B:)5B5mu#"},
{247.62,1031.20,4.93,5.91,"?d*###A,$0##?Ka###.,#o$#LKa9#$###z##gY#1m&RG#2,#E82$##A,$5##ZKa95####J$#JxZ$B10,#gG#%##Vs6$v&B,$'g7(##L5$.##RLa{k####N##;wKKe.C5#&e(^m(y/.:=3LL6MA04,#a5%$##0Oa>c%###&##Rr1-bF###Q>#3o&Kg_M[({P$"},
{254.67,1054.73,5.09,5.54,"D>#.##OG#$##oDB,##D>#~&#GmSbG$###N1#?//_x3}b#tc$ol'G,#g5%nY#PmSWG#E>#E0#joS{H*###q~#wM55]1$q0#h3E@,&##I-$SaA+pSHl%Q6#/i=F>7A6OjY#g>$&t0^$SId){b#(E2du&S##]o0OT'|d-p%$4YJV$&Iw+0z0HR+%G:o,&,Q$Y>$"},
{963.54,57.89,6.26,4.65,"W>$######ac%F_>######rP#+B`###$OF$##6#$###,1c###-H&###6##qu%h;B###,##^5$R0c###yFK(##(c$###91c###-6'###w>#0-&S^:###Q,#0$&F0c###7OG=##T,%###21c###;R$Q?(oY#1,#1B4fY#N>#iP#=&^###0V=*##I#%###i0c###"},
{1073.64,69.18,6.49,5.21,"C/%}#'|y/}>&?c7fY#qP$###fR?###Q?%###c9+###68C###V%)=/,6EA(##ZBTP5#4I+I##s9K######)##1H8<5#DC.###{6*j?$}%,OC3),N8##Ml%sh&{p7#l####$-#tW,Gw,E,$###+%+Zc%~x'E*>U,%@@)q5$zK*-u#d|:###c>#M,#RpN######"},
{880.88,92.94,6.65,3.93,"%##Mc$iu%MI+l#$WH'tb#2H%)l####$##c,%sP$###'##3c#tl%|+1bS1=Z$jnE<{8A,$9##7f0###'##Nl$VZ'###%##eB$X-)DM&M5IKU0S8Zsl%tb#|{'l]6######fB#7x2t##6#$>;#.,#&##CS'`PGFl%###$##|+1pm,C##ZP#j_#`$+()#d%0()#"},
{880.88,92.94,6.65,3.93,"%##Mc$iu%MI+l#$WH'tb#2H%)l####$##c,%sP$###'##3c#tl%|+1bS1=Z$jnE<{8A,$9##7f0###'##Nl$VZ'###%##eB$X-)DM&M5IKU0S8Zsl%tb#|{'l]6######fB#7x2t##6#$>;#.,#&##CS'`PGFl%###$##|+1pm,C##ZP#j_#`$+()#d%0()#"},
{1034.72,93.55,6.29,1.25,"###N&#sPP###eY#f1&>(<q%(yc&kS,E$(/r6GU.04BP5$c,%{k#M$#2QP2,#KOHf##KXER_)b?)A-$0UP4]+nY#GK.:K.Ql%}k#-##VSP$##oKO4Q%fFBEc#po,gR(jSP4l#95#~G#G_.fS0RG#'##p{1Y>$T.)qb#qE/-E7Ym'6v%l}?x_5%d';6&BJ-ic%"},
{1034.72,93.55,6.29,4.11,"j$&'&-c>$=#$asB?m(NQ%Qn'qE9(JMiw,K#$WP>,[(N5$;,#Xr757+RG#hu%(wQ3c#q@+vR&1E@4[$NxQwc%?wQ0,#Eu$Y##dw.e5%J?%,|=>wQ/E77?'b$%NNCer+~1=@%#WuQ&##pb#l&#xY#H%,7t?,o1(?%HW=?,#gD6.06m6'###EV&~OJ$#####0(#"},
{968.77,111.67,5.93,0.02,"v?'7?'###%##E8-Fw.###a##O%*&y4.u#1v#=5#^h&Ei?Ju$3$%wFEOG####k8,n%UD>#,##(&QsXIVG#ar+[c&:)0$d:h[M1,#RA-V8,ac'p3?~T/:94;,#l)ULIK|k#,8$ZQ&l)U%o-yv(###+##%V7j>%;7*O>#g;;.v%8M4Uf.Hc%--$3,#Y8)=8+%m("},
{745.84,120.07,6.73,3.90,"1c4O7-######P|V######(##@&0.,####P5#zQ(######-##K>E{Y$gY#0$#]_X.,####;@#Nh;eY#-##U7)J#%######7$%`/4$##95#+)$G^XFl%###Hh$)@@wu''##su$bZ'######0p((Q%Rc#C,$-W&W//)6%PG#%B%|sF95####Z_#8.-I##eY#(r$"},
{134.42,134.75,6.14,4.68,"o2C6##a-*`##j>P######_##:>M###u*F(##ZP####aL`###L3F|P#W>$w%#H3B}P$95#q##y]^95#i|E+##(c$###}L`$##We0]G#95#h/$my9Fc#bG$%$#8K`.,#|_?u##Qu%###eL`%##6Z%[Z&###Sc#0B5dP#.,#l$#1-S###R]4V##h,&###RL`###"},
{265.77,145.66,6.21,4.95,"2S%DjAE[*=5#M<-7.-######pOB#########e%P###zk@/,#FA*p@+FxRSG#&{R`#&<Z%M5#|QP######<##'yR###Wn*4,#%-'1,#&{R)x.-vR###@v$%k1|EG######I$#PxR###lY#-##&##J>#DT);f2Yl&###J##EO97m)######P6#gD<######1##"},
{668.30,154.01,6.37,1.51,"###)##wm%4[*g>#F/-.v%V%.Mw){?*X##_U8OG####l/#Yy8###E##(<:pb#9C1Ie=Z0~T,$84~uL7o[)tA*Cd*4,#S(-}m*###$##zh.%Q%^~0Oc#25~iPHQ0~,$&pA*#TDe@,XH%$p26c$###+##h?%I#%.,#.,#Y7#V;@tl'G>#RQ$83>0u#H>#+Q$_>$"},
{226.61,168.70,5.98,3.51,"###7##F%*s-+###Ig&419,u#RG#Gh+f[+^7+.,#s5#p$)c:8X>#<`Ti&1,.+,g-4#@-i<QI+C-(Au#nc%NAM###$##mm#QPI}-)<`TL,$Uf,>^TAD;###))&LV>######;U,######=##2f1M5$Dc$c#$/`Tlw0.,#####v9=6(#######f'######$##9h3"},
{771.22,174.25,5.73,1.69,"Ml$|l%@u$~#$A@HD>####t,$xMVuG%###>](0o,{K1ow+'W9kH'cG$|l']>#)O>~?)###5Z#^NVwQ)###D##of,MRJW>$(##bG#B,$n~)y?*@2?0Z%&##Vo'=JV8e,TG#rm#8$'%6B)l#u#'######1q#/IVSH(###^##CLV^v*.Z$.c#KKV$##>-&Iu#nIV"},
{160.10,177.69,6.28,5.61,"WL,5-%}}CnP$7rSbP#b#&###J11sY#S[*###)G4EZ&[P#####K0VS'vtH;o0)oSF5#8-',&*?|@=v${R.5##(HB:Z$xG%3,#RYEiG$]e%`oS,vLAZ%JV+vnSug/@.,=~-b5$9bD0c$}>%UH$y+5VL9]u$N6'6^(Y5JVL3E?':w'_6)1-'dS.GB+g,&>,#:[)"},
{385.32,183.37,6.65,1.50,"###gQ#Fl%###g5#C{53l$###1/(lA2J>#H>#cP#y5&(l#;5####H7%,]3###T:2+:F-g_9l#vl_ys@T~-07%O[+Ru$nT3g5$.##nu%Xy6###3T3Z-$vl_N=@2h_k$(y/-C7;`[*_d'je0K>#/##W$*'Z$###F>#`#%NS'C80F~-|Y$?##s/.[,$R-)######"},
{64.07,190.18,6.01,2.46,",##rP$M)-jI.QG#2,#nZ#(4G;6&.,#k##&}E1$&r>$H$&LJ0^$#pFHwm%G[+://FN8~6&'B1:wRt>%@-#yyR]6(lY#;r(AvRB##(A/w%#zuRuH)rR+-##lzRywRU,%2##Kv<f8+e&0JW3Uv)######O%#`,Q95####g,#W$O[5$B,$Y>#F].B5#)$%+]0pP$"},
{64.07,190.18,6.01,3.55,"xu&Fx/###0H#cz;a#&3,#Li)}$*t$*fu$7rS###AQ$Ep)ZPJ3o2iY#TG##]$C.DH(;<v(]A$(]E^mS~>$$e%I@+W7*>-'0R'KQ'k,#V$)W34w]4GS((>KYR%UnS26'I#%1&#Lp7+l#/,#MQ$###^)-nl&>/._@.[%$2I+sM(WjH###.,#Kz#xu'######Rf+"},
{36.36,221.55,6.57,4.81,"#########>'0W>$######beBFl%######x2U{k#######3:LhG$$##^P#a?C3aA######gt3%.U###Du#g42~m*lY#Y-'Vg-~P#*##Fu$Eq2{U::5#<5#Cq0(2UA['3v'~$%}S+At76^60Z$######7Z#VW>-Q$_,$Ve-p$*Q=1mA/%Q%%##,.%:*5tv)PQ'"},
{381.42,281.13,6.76,1.45,"5p6{k#######@Wm######$##MQR/##ty:###@1;.##ez>###xC=.,####$##ZWmQ5$###)##&.TYC4Ro4###o{B3[%HU:###Xz=;5####$##KWm@c$q>%0##Nr>KN:q%.Z,$m4J#@(wG%C5#tK7+u#######TWm:##dc''##3h87S(Ry7[G#'PHZ,%PG#3,#"},
{294.29,290.41,6.47,1.48,"sC<OG#######;rg9##zH*###R2?*-#)=HN,#6HO$##D>#=##W;@.,#######}qg2##gQ'R>#Sr@gH%az6dm&^5MsP$hY#)H#3)?C,$######bqgA,#Gv)'##'tCb~(>kJVG#~IPEc%pb#%##Fy5oP$######-sg*c#uc(###*E@wc$k4L2##gkHs#'+u#.##"},
{234.42,308.53,6.50,6.18,"~#&######e6#.7,######K@#Rd+######:yFZP#######Kz^xY$######j6$)ZP######s$#ny^~#&$##mG4PA*]2B-##Vz^9Q&######d>#>RT######2##hz^6#$(##U[&)T1S%/I##-5HUl%######%##&<A######$##+z^.,####Y@$M`?Nc&###$uA"},
{751.71,335.35,6.22,2.16,"/,#S,$.,#6,#K>#gG#:5#pe(V~W;5####F))d~W.,####i##B-(X#%###y6$IXEQG####4?7[y8######8aWu[W.,####F^'FR)U6%}Q)3A)g^WQu$[P#2_&*dHD>####Gg/mK+{k####ZQ$pb#5###$'y7@cx4$##.,#~b2|%/######Gn$lv'.,####y>#"},
{751.71,335.35,6.22,3.13,"Z0*g%YOG####Z/Cip:###;##uv;M..###<##|P#z82.,####n^:iY#P#%-##c'Y3l$###l%#8AQ+05###7%#`c$OU8######=&0C5#1Z%$##NMO:3DZ>$9##X~@<6SE>#.##EM7*$(###&##Z6):,#uG%'##W/0|H)WR(k5$/(Y;?'>Z$~l#>QP######NZ#"},
{462.77,349.33,6.54,4.47,"1I+######E##VRX######C%#rFK_5%/##u&,^>$l5%$[)Ad&Gf4######G##*SXDQ%.,#^$#uaC|UXA,$M,#Bp'x4JVZ'$##}T8######W##ESXv$&?d*t##&(5&I:AXH,##{UXtA1uG%X##Wp9######G##nQT0H#ET5m##y0:_-#VRXQ##VRX'##j>%]$#"},
{194.63,470.52,5.99,1.42,"###qc%+u####`,#|r6uG%###/p)pK6%##L>#*$$Qc&>5#H>#'##O%&TS1###dU2.T?2y`tY#;$af<=_R+U$%,S-TQ%)8.1Q$%##E[%wy8###795^m#P$ami;[z`4R'eS+sb6Ee,@n'Yv*'##%##7H%7Q%###OG#K>#Sf'~f3fA0sb#xG#9g0/d%iG$######"},
{304.37,527.41,6.51,1.36,")##k#&A,$###xY$D6#-]39c#Bv)_##Yw-T~?A5#Y%(v[,LVSM?%#.+6l$Gu#~(=Yc#lo5s5#CSZ|G#Qp7,s(Jv(b^'LSZrd'iG$;5#s>%'[&hD?P>#My6b>#aWZc?%<z;A##G7*0J$,UZ~v'&##xb#L5$/,#:v%Il$Py6###gF..B4M@+###Zo#&*DQC0'S,"},
{512.53,527.88,6.45,0.19,"$6$x>%95####UZ&M&(I#%&##EkJ2@(.,#Lh&6e-9z7$##f-;lp#WB5_5%###._87'(?C:T>#-8WLA)OR+`7#-o0c<WC03<.(1v<[6)RI(3c#<OB^,$dn-g~)@:W#d'`.+Y7&o~'Uo/&33m[T}x5$##9w(~x'S97######AA$[o*L-)/,#q#$P8$(*C&##&6&"},
{866.62,532.63,6.57,4.58,",v#D;@###x>$qE<4q<###sc#TlGwu')l#7:%Sm(jd)vm+F-%(~.OG####I-%39^######S?$99^(##7A/Y9&`Z':1+%9^dQ%?o3######*Z#Y8^######Y##Z:^F-&l&5c##i?(wI%G9^lY#Fl%###}##2EA-)@###A##L(,F~P~l%~l&R~%RQ%l#%c#Dw-+"},
{543.88,534.25,6.96,4.56,"+[':5#nG$al$$##`P#M7&>z4&-$QH'PK&33DVS*+u#h##[[+pp:*##}c'?T&lv(<21=1:R{0U6FZf1u]%T|D,I*###o=2u83NmTD,#y]3/v#6h;0%%tnT_Q$.nT4?%}n-p@%KB6{G#PpTWc$SO>:-&pQ()##Jg/z6+N>;V@+RmTYL3e-&Sy,S.-M|/aV?+##"},
{904.57,540.20,6.53,0.05,".,#AY?Gw*D:]S7&Zj;m6)j/.xd&{e+{k#)##1c#cB-9#$/u#&##W=]Cu$,u#t2(S:]L5$&##+S<q+HOG#W##FI+/v%^P#$d%+##Sq;######|4Eg]6###H##&:]0Z%$##p.#`e0E>#2l#7v#+##Zp7######`IQ{@/###(##W9]pb####C##d~-1Z%0,#R5#"},
{904.57,540.20,6.53,5.29,"g/(pb####/##`Z6AH'<5####Wt/R=GjB-sP$,m&gl&^/&+E>.K*######<,##$BD>####^R'#aZqQ)###C.$to17m)'##U?$}[+######wP377,######4bZeHO0Z%'##&%C~.'/o2fP#k,%@u$######t>42H&'##SG#bF8QG#F5#Q,$RK3H>#rb#Yu#;n*"},
{313.15,543.12,5.91,1.46,"D>#`##FS0vY#uG%Q##fI,kG7.,#av%yH)sUL+m$Oe,:v#jrA;e.Z5#BS0f5#S7X)##mm*63)zc(3&&08XxS)N==$p.YA2K>#k&47,#Ne.B5#|9XV#$h@.N,#HS-).$;8X,l#N6S'H#Tz;Q%%:-'`P#B%-&##d@>e-*;?'$##@T%/]28cB7v(KS03%,=n'WM2"},
{313.15,543.12,5.91,4.63,">e'023h./}R-GY?M?(nf%io416'$##~I=4R*`I,&##YZ&TG#NM=%S%L6Sk>#O8X*l#tn-f6$=w-J,#I:XP#$E~.6,#Rx26,#kA3@5#/Y>jx-=8Xc~(|c(mw%-%+ei(~7X'##*/1s>#Gn.T5#D$$5{>e6%#/-u6*2^H/,#17&Q7,Ql60Z%L##ZS1(c#OG#^##"},
{391.13,544.98,6.72,6.17,"mo#,w,95####6Q#NL*+u####|f4$V33,#&%'yQ)eH'%m#.6EfG61S-0p1)##&z01=/='75,#UAT@25(Q%P@$U3<VAT8##t;6~;=hG#-wM,^*S067$$54:#'/w@Tmy.&o,tw*06%uCT+9.A96Z>$L##o@-5q-ZP####k7#P(9lH)3,#c.%bf3gu%G-(,@%*XE"},
{357.37,554.04,6.57,2.76,"1J$E$(?u$###/m#K0395####WQ'%K.?Q%'e%Qc$(?&yY#H?E;'117'~n/Nm$Lp49;0]/3nG#=7U}%,>-(lz$n6)A-HR~.$W1;c%1##Xn+a`21L8r>#d]N;v&k:UV>#+o+H[&6y2U['0#55}H######KH%+%)@u$###83-(v'_^R{k#=I#Ae,}b5aV@1-#Q'7"},
{440.45,589.48,6.89,1.45,"1S*u5%/l#o7-bI(A['2?&w,%Y-*7v#?7+iQ9.?$hm(<[)&BMWOB.?%7@+A#$[aGC#$BT2y,$4wVL>#e^5RE*`J-B@(.yV>I'14HnY#L%,t$&CaE<u#C8/LZ$?|V&m'b&3/Z#9B-*7*0zVov(Kd)Gm)IZ$.R&Q[*.c$h,$I-%_{(gx49I&gG$ro#;dTf$'GQ&"},
{787.41,592.12,5.93,3.05,"5w,$##;@)(Q$;c%'##YM6Pa@<]4###hG#7%BW>$###^5#~p4RaC###J5#=7%zbN2,#kv%z_SQQ'X,#7nCd]S###(H#:f.bd+lHP######u$#v_S|Q&e6)^d$Af-<V'`[Sf#%Gl$'w%fn00,#0N>######O##UTQ-6%_5%,##$S)IK&q}M###lZ&^Z$qQ)###"},
{504.39,601.34,6.52,5.66,"ec$%R'1iUVu%DO;l,%`-*Im'Lm)G>#0##co256'tb#B,#@7+~%+]]'vdUa#$EgUjv&e:=(6#mp64PB2,#zb#<@'Q`@5,#NZ$j>%|(%xdUJy1G2@0@%%U52</~.-aG;0//,((V]0Ox2.##m42~>$Ul#^c%:C7D>#&##f'2-/.4K3V-'V8.0(*|x*{7TgH'Fx%"},
{522.28,604.54,6.48,0.10,"4$&`>$ZP#0[$T,$JT('m(<c$2;;}:8QG#ge(p$+Px0(##u-=|z%E'5KQ&%##f&/m{)6OGC5#rwTh^/|u'|@$,^5`{T/H&V])ql?OH'|36F6%7{;<?$;X;,o-UxTw%)q8/r7*jZ%_*4<_6uiA8.-;##1cCVx)=R+###<$#eK.?d%Y5$u[$BJ0_$%T&0>Z$'x/"},
{445.78,606.25,6.38,1.48,"N6%U?&[,%8u#o#'bZ#h$*t*2}Y#-R'|Z(^yKJI'OZ&am#cECZ&3@#$TI+8l#HmTG>#qv(o2)*e+rv'3pTcS)y7JQ-(g80O,$8g7A5#|@-,c#OqTlG$Mm)8Z#v01LH&ZpThu%8[T[,$UU4T/){Z(_>$f>$j5$u|/JI,AZ$wb#iU$xkM+o,5Z%Q6'&N<hc$C7)"},
{445.78,606.25,6.38,4.42,"I6&$f))~-Tg3-)88d(S'%|C<#v&/,#Il5;7,Wl$Sc$>-'hY#):8Ho%K%X)Z#8&XcZ%{A0<[$>w-T[#u'Xo5$HI+36$wo5jP#6x.<6'?QD<.*9&Xoq1V[*_.(g-)O*+}$XJ>#SS0s>#Mf4aG#rH#qRVN6&h#&#v&4)XA5#B$&&$'TV1EQ&$I%El$c,%v?(X#%"},
{489.48,614.86,6.16,2.73,".8).?%S6(/,#|Q%L]/.,####8I*Po/0Z$.%%'v&5l$V>#'mB<L3EJ(ZR,u?$*g3+N.~08~>#YwVqJ-^#&Yp#MS/i$KP6([206?'8##0J,@*2Lg84H#;zV1$&+zV(Z#(S*s?&5^4,/*e:/L[T2,#?Z%7$'*d&wG%{Y#)t5Ql%QgOCu$TR%[-(_,61{@?###q8"},
{843.55,615.70,6.27,2.82,"^$#n$W######DQ7f$W######i(OKQ'%##U5$9R&/-'.,#[G#.##Wu%QG#.,#mxQY,%D>#:##K*WOG####`##HS,7d&2v(9,#cP#o>%/,#7,#;RUJ>#D>#-&#|%W######_'#aA24##[]2j?&>c%95#L,#|Q(cmU###%##K7#o&W######9&#Jq<###-H$%)-"},
{369.46,626.49,6.95,2.67,"f#%B^(+aEP>#(Z#%M5L5$###S#%Nq3@H&1$%gu#&d'L-&l'4WH'%r&he0OS(F/0pq-d<D-?$URT*K-aS.?L()[%A}C[15nK09u#%Z#^d)+U.Ne/`>#`VT'm&~TTM>#;U-J/,Y16@I*;40GRT8$$y6*P5$4,#|k#z$)-C0p5%AHCZZ'1R%-f/@a.sQTvc#-U7"},
{857.25,637.44,6.42,2.12,"By#llJ######na,2v(######S0&nI*^l&$##xe$w6*e5%###|l#eu%.,####`~?95####?##]k[vb#NH&b,#%/*Hu$_n*uZ'O,$%l#.,#5l#<wV;5####Eh#Mf[###4,#XD#q@.###/Q#>J+vG%:5####CE1HU:######kZ6]2A1u####BB$37)d@,###)##"},
{283.46,639.08,6.86,0.12,"}k#-##M#$sZ(Pv%1M2Pc&U5$*{:z6OG>#r~'XR+Mq6<u#VrSf>#lQ&jv+95#=&1Gz&;$SM>#%xV-:._d+Do#Lg5?|VFe-]y+F%)mc'Jm)kG#z@/L,#SAJ`[*nyVH-%n]/{%,{[&(h1|;55wVUu%###B5#3I&ZP####/%$<7,I.(6c$J6$}S2XI$Rq9O,$wv*"},
{530.45,649.39,6.39,3.87,"&-%N[&q.NBI++a;Qv'a&-Xo3n?)%##Z?$}RV1%+J#$^d*/f/wx3qG#WWV3%&YSVX5$ZPG2-$M05:S,Pl%P5$V>#nm'Dr@95#S6)|1$pRVX7%0<BO7'Bz:p@'`6(ME?)?$md,}G#5E>qn-3l$YZ'{M'lT7+A$Dd*TT&lZ'<)-5L2hg1t$&8048A-,R(1l5$.*"},
{564.37,649.76,6.80,1.30,"=C4_R)n-)cv*zv*6v$fJ/&?%aQ(R[$cd*;c=e,$s-(@-'CBV@@VuP#f08ic#q|GVH#rbN_,#>@VS5#Z:9+j+aA.0n'z@VU.*UT5u[)3J+Bx*[h<9[%/r=Al#REVaZ%sL;W5#891|u$kAV:@(4H%}6)K>##%&{5%V5$#0,:m%d*-5d)Jo+Z>$2('Io0q@)[6("},
{450.58,664.38,6.78,4.51,".%)46&#R$_07ny2mP$Q?#&2;v]5$##dP#0(-nP#ll%mH%.lG8eRWu$QH&vd(12?,##)X7.^2vcRX##)04D1(:L4e(,KZQ*A*ndR:?'8-&=7&^4D]91{PCfl$PfRu#$@E@V>#WV<ec#LfRYQ&>J+|r@Y,$4d&4:.FtF<c%%m#9Y4jI,UH'cG#38&5aD;T,W6("},
{569.79,665.34,5.74,1.31,"+?&s#$}7.rY#_m+#?#+m'_46dP#$d%n[+khXCd&%v&)v#~jFVT5)?#r2C4##1SY9,#?@*:_'$S-]v%bTYKS)c5D.['+T2wb#Lx2_l##q9]G##WYs>$V[+U5#jK1=v$ATY[#%B#OS>#QL8_w&s5%6u#1J*;Z$/Z4UQ'$.)###e2&_e/v^3TQ'nQ('A)Am%$x+"},
{531.89,675.66,6.59,0.48,"tP$[P#~#$RB5~Z&tv%Wf3r[,L*Fov%{u'u;&r-+?D(b-)ec6B-&e-*~v*;H%l^:jR$m%Za5$d&ZYn&an/UA#Me.ZO-WC8}A,oK%nvTy6)[P#9W7v[,CN=<e)#(ZGA0td*/A'8%&`h<ef/nL6~A,7I*<y'Vn-PRP6#$4##X[&Ke-(I*%##tu$xZ$MJ1###C#$"},
{563.16,687.70,6.88,4.38,"###-##%6&Xl%F,$~G#Qc$Gy.Hl#3v&?t2+x.,U)NH'hI%Oc%v6+&##wP$>x(Je+/:-h'9__03nKy.-~O.h8//x-=5#B]D@H&v;?rY#&m(Y,#DEBx[$SAX'6$NAX|Q%$OB0w$a~0Im#'BXFl#;x)ol&rZ(###AK+'q:#`7o$)K5K>@GG7+b7(>R*uz+d3Gl>#"},
{418.17,692.53,6.59,0.05,"T?';5#D>#g$%H@'PK.vG%Y,$:1,y-Q/,#Nc$F$(t.Q5,#oq.$2'<x0MZ&4,#yd(3<,==I@5#j.QR)6},'[?$TC6c/QF>#2p)[%C9[)M`1+I&$z5[v$17IJ@)p-Q^&+F8,3e)r5%v-Bt%-~U9a@.Q##'4>Vx(d6*###/~#-f,c,%^P#Oe#2U7v#'`P#Y#$XM;"},
{418.17,692.53,6.59,3.60,"Fl$[C8c5%7,#'$#SWB^#&.,#fv$t~-9PFB,$w2+k`Ctn*6@+;J-Tz49d(5-%#B,I^6G3;>v&>AN~m*hB,D~-{$(q#'MU%plQjx3v$)gc%m#:m@.F,#XqQN_0rlQ0,#{@)O1+jd,c>#.m%5I*Lv(D]&r97%$;rb#4A%Bi;|a:Yl&$##@%%b&Ppb#>Q#kw)56'"},
{910.03,800.39,6.14,4.61,"'l#jc$js/vQNwd,Kl#Vz.=JS?IS<#$Z$'$J&cQ'Jn*O/1Ku#?Z$i43rLS1~,o.,=v;>r9f.,9JS}6&-%*xm%H%.*##A<4XZ&;5#43/|y7x~1JZ&IR&`A(VJS'lM~G#eG#weD37,###om#^I)ZZ'<#$P>#qZ%Ve0###$##|],}v,######OT+)c$###P##yH("},
{399.90,823.32,6.20,5.35,"L,#;/C7K4m5%ZV%:wU9u#;v'l}2VZ'###fG#Uc%9,#<#$%Q$TG#o[%0'N-y1pq--y*f3>m{:}{UTZ%Vc&yv$]9257)OG#(##.,#J5#GO5.m(u#'^H#IE1IyUuvU|b#X5$~?:;w,Om(###d,$?d'1I)uu&r,%;Q&###3##<12iZ(######/f'ZP#/,####D%*"},
{654.32,856.51,5.98,0.16,">l#$W=zb#V6(&##4Y4-e,:?'c%&u+?nZ';#$q`0<R+###)##7d'A%M~J14,#@K3DE9zu$+:1A.,2$%Nq*`(9[#MOG#Y>#a$%gd+{Y#~;3#m'/wPiY#5H$j&*LR*cP#5r0%nOcg9###PI$:y47#$$##I~%;H&RvP:5#@Z#c)+k19vQ&]D2GzPNc%>m#2zPCU9"},
{654.32,856.51,5.98,2.20,"$##6'-;$)###eR+hq4.,#7##LhQ<m(###g(,@*5~`9-#BUuE^c&I{0|q=a>$@bG[]2UG#%/%7$<ur<@c%~[%T1'cgQv6+)##1Z%0##.34z2B;;7>e-yS%BU7qm(,[>[H'gG##q9m&1.,##$##A+~P#[Z#G&16l#d>$.>5)~.%?&gZ%1L-X{:C7-;5#(##VU-"},
{997.20,865.43,6.33,4.57,"~{5'%'d7ByZ(LhQac%$n*:5#uv'9',[?)###`Z&=H%bG$T,#Hj4X?<x{@W,$&hQ:&)-J/`##IV3>t9Mw.4,#d,G=v(95#C##G@+ex,:_++=8XcQY#$CH&KM'P06??$v;8my1ncQSG#ZG#VA'NB(Tq<4.()&*x9/an0'##J&%w6'Jm*T,#;f-B%+8Q&$##6[%"},
{1008.29,890.58,6.22,0.60,"[R*Rw)z,&z[(*I(^y5-u#E,#'##]/0{#G######4,#~pN###D>#w^&e&2/J-Q?%QpN,?%Fm):<;z`CWw+.7'8l$7c#~pNy#&3I+2C#TV>1(0:m(A~'M+=emN<7MrB1r@+8w(CQ%|{-($N^5$gND#@$h&32C.&.+}S$JlN[$'Z[*-m%gS0=Z$%##pv$J*BrP$"},
{1008.29,890.58,6.22,5.30,".I$Hf2T^.Q<D2C#L_=aL0Td+.h&e&2|@-95#Gn)+6&@w(C@*'T$HlNBm&`m*A~'[4=amN:m(MpN}5%Qv)R?%]y5}k#G,#A[(9v%QA0=Z$OR*rB1f7+0n(PIM)jCUw+$.'jr:sA0OcF###(##&%%F*BfG$%##,*.UHNi>$CQ%6c#~pNnu%-c$4,#~pN######"},
{398.50,902.43,6.19,0.65,":5#s5$Od(eI,$##}e)GjB(c$LC1Op4B%+4d(P)9Qu%###8,#(Q$y4:vv+<#$$]0o/.^K2801v@.?5#NV-PW@$X>###4c#-[(d?(P7(CX;K$)--FH$)[#$-e)kT3{u&V8,_'Ks?'###Dj-|TT.,#0##Q&+jH)uT5D>#5-#=T0DVTHI,Q>#Wj.8'/gV60WT65?"},
{398.50,902.43,6.19,2.44,"5,#sT0L5$###p^7<U,A,$r5#C{Sj&4<5#j3.by+0BI2xSD`92@*D*7)A.u,%O=Ef%-I>#1y'c+3XlLn-*;.(h/'~zSk-+,##r#'V>#fB-~{;oU6*J-4%(2(4LI)vb>:n,a,$RL9Im(eY#??#L?'lY#0H$c.+{>%dP#3h+y7-w#'`c%r]*gN<h6*95#)##dC/"},
{413.97,920.07,6.26,2.15,"/l#mY#j>%###}>${c&###$##A5#$K,tb#~I+/H#U%.C,#8y5@,#7)<A,$###ww/wD9###h5#J8X{?'=5#n{Sim)G-%FQ;j:XM7+^q86v'`v&kz;z[,mY#VC.%=X(q:>#$P&)di,%=XDkGtI+]H(0,#n6&/h8w~->B2I6'bm'=8'6N<e5%@u#/w&*}BZ>$Gl#"},
{287.35,975.24,6.79,4.59,"gT(Jc%95#6v'Rm%h//mP$r#&6r4pH)###lG#gQ(######3H#<`Teu%fd'Rl%d(5/`27~T6l#<`TNJ.o?)5Q#uL;######e,#b[Tc5#nkBwu%;4IC5#G^K}M:8[T=5#a.(+m:u(>######t5#g[TG/(FA1.?##,G###B6#gJ0W#%###<##((3nQ':5#<5#@u#"},
{336.02,1014.36,6.39,0.82,"###~6#7&2######LD-%m(###zY#G;8a6)bQ'&~#>81Gl#Q?(#c#}y#,WBD>#&['xRBAS/D[*GIPEU5^v(9z/Dm'vJ'>uG{S/h$+yp#=K5Po(u-+LT&Lz1rKP#AP4y.6/,M'3]Z%5r)IYMUG#lGIFJ&a@.][#)]0V2*#:7>J(XZ#c`3hm*fY#O-)p''Be.B@&"},
{336.02,1014.36,6.39,2.46,"]6)^I#sn0CQ$IS-Bv'#B'H[OaP#*-%8{)nYOl-%T[+'C'9e.j-)L#$&M6i6)SU/f[OGp3}$*%16J[K,y.>&-B,$#?#;E228/###$##CV+]?)%R)U$([$?ve1?^O..,t/&K(4Z/)qg7jh(U(;######]6#xe1.,#O>#V]#-4I2~&Nv)fg#<V>x?#uaDAJ%5f3"},
{368.77,1016.94,5.96,6.04,"oB&yrCW5$&##%K$_p8######0<,;m)###'##Se)(c$###?u#a:2BQ&r7*%xD`e*td);10;16TnG}5&5Q%C##=9)780;#$cG#Gy/rITEc$:)33w)ZBNz;;i@,7@PR$)A?']-#`.+)7+86%Q-&O,#i%Ui>$ZP#yd$#&U######/2*M$U###'##{%$@$U%##E>#"},
{283.15,1022.29,5.90,0.43,"###F,#Rn)fB4a#&cV+4Z%ix)`?)D>795#00({k#<9+In&d@MH5#Z?#/f2M5$=[*eC%?ZOJc%2.PMq0Xc&(q&ov)*AD|?(8U.(A$u.*,804u#Z$)1.&@t>JlLty4+z18K-mg8xP#c)3v@,###s^$Cf3D>#2,#Cn&#wC=$(mQ'+##k8HPI+5l$f#&Ra?:l$`6%"},
{283.15,1022.29,5.90,3.94,",[)K?&Il%$V/O<D.6&N#$.a/7$(zZ&M%'b4Ei,&###`'#JPLTn,}k#SG#{h0<q5SkE=~*_&*/@Kr*AgQ'9l#7z:.,#R##{c'W95zb#f5$2W+3@+S[%V^O+i1LZOj,%6B1Ko#z986%*{k#~##::3N<=,?&w]%=I+UE.$jD]&$xo6c##NS0cU$T?'3q*D]4al#"},
{296.97,1049.44,6.20,4.80,"sf1v>$1T+JQ%5.H[l%fu%+v#HvGI/*w%/i['K8$s^3;w,cG$W{<DQ$^&-i.)i#Nt>$Pw)9n$~YDbA-A$(Tw(ZK0p-'>z*{OGi$NI#$mj>/]*&ZECQ%Z?Kz,%r[LIw-###E-%=%N}Q&S3:%S+)'N]5$-p31,#(4CCu#MRJ0c#iOGfY#K>#rJ*((4DM2Tn+d6("},
{933.71,98.45,7.78,5.71,"]7#2RO*c$###+2'd;B######k:%J>OX##du&/(#;4JS'#()@Oe'2)9wr,G}GXqQo?*=##iq2@j;*6'###<##qL'1T4###Kc#uP#16&mK(u2BS//yv*vJ/;z4V2-'z:{k#7##}:0@c%/,#Kv$V>#XZ%GZ&###&##|$%QC:/,#or=28->$))@#S28.x0e?'R1,"},
{820.57,114.71,7.25,5.62,"q%#VT4*l#ZP#GO,5&2O>#_l$L?7ie1###L###r#$3De##_5%z{'U&43,#[P#HEYxY$*##ro%kjE|Y$###g&#fI$;@,###&##{$I3l$jY#xZ#L(Z.,####Q%#a|A+c$###K##8Z$Z5$######AM?###C,#hW0Y%Z######='#V:6######+##'6%=c#v,&QG#"},
{745.93,120.16,6.67,3.88,"/Q5B@,######ajX######%##I04######D,#~-)######-##{,JeG$OG#G$#cgX######`$#I3C######o$&`5%######?$%yS3######fh#pfX######7C#vRF######Mc#=-(######9g';#$~G#$l#+E'*g5`P#####T#t|F######t(#*R*A##ZP#J2%"},
{136.73,154.71,7.29,5.06,"/x,###5J'Ec%f^S###sG$###e8J###Be)###bv'###z_S###LK3l%'cf2R>#>~S'H#;v(;##1(Q{Y$J>#&##jbA###2UKOG#Q@,r%'2J(*B-lZSEl#^P#Ug$.cL+Z$.,#E##F]S###Vx,4,#nv)Cl$bf'p02NJ0/Q$2l#R/&'2>$?%.,#c##K}B###WG#5##"},
{995.05,177.34,8.16,6.28,"A.*1u#/)92-&<x+nV4x-)}Z']P#i~(Un(N+HF>#L#$JL*b5PzY$A5#nh3wP$-C9vY#Q$'{C5Hl%(?#'B/Q>HJ5#'KJwA1CB5Xv*8,#tK+1-&g-S###L])]$'@95Q##8RFUM</l#VM0TyMT.SM7-###C$&?l#4q;###0~&Nw*CoRUv*gZ#aeP|#$?N=<m$)/S"},
{143.41,186.83,7.53,5.61,"V8)dc'Z-$FT483,s#&?(2^Q(ov>9v(L5$###u)+G6(_P#95#6k2l$*g,%sb#WV4}&*},MP.,)hR{G$bZ'DI'-<9fZ%I#%%##Q*A^P#x,$jn,kQM?Q%$n'WgR:-O<-'Ue(WJPo^4vZ'$-'cG#)m(/##/o'o|E7P9A_:U?'sd)iC.*=D1%*#R'c7*7Z%'l#fn+"},
{963.62,196.80,7.25,2.51,"lc9YyJKB6_##ay/c-<Nq=)##(x+z.+}$,,#####&m#l]1###@U3m<-LV?y##r%U4{)/&2_##NM>u9(_m+V##K$(qe+`?'mG#(p-`v</%-n,$_%Uww'@H'h?#/f0K(*fQ(y7(Hn.cl$fl&&0%By4|7+D$(Dv&w$UjY#95#T&$Q?(0l#hY#$<0{k#F6#,@+4w%"},
{792.26,199.08,7.39,1.41,"3z,kvR######VD0juR###0?%pJ$CcQ1,#YQ'E##&vRLc$6#$EwRJZ&###?H#8wR~>$`G#8c;qc&*R'{C0&dO&##J$)+%)%Q%MvR###$##K[#^yRG::`P#A@%#@&uBQ%02UZ&-##6o-~6)###=&2fY#b##M>@K93`NBg,#'8-$v#SvR?,#r?*vm%h@.###kG$"},
{554.84,306.60,7.68,3.06,"4[*######l5#lYO,##~>$PO0bc'S8'-T/@.Q###f;<{G$qC=mf6###U,#/?%]wZ1##Z-(b$$bJ2Ef&9KH]D?###Qz7k$&fwZ-o2###oC3&6%NxZ###?,#dv#;_;###En#o}G######V##]wZS..A##8jD4##byZ%Q%/,#p6#(V6YC<###Vf.$##cM@$##B]4"},
{533.02,321.19,8.65,4.78,"[FC+##$J.###CcJ,##5cK###{SS###PG####;94###W05###5dN###wc'<5#r4H$##BnK@?&$TS####6%ic$~D=1,#r%/%##)JR.,#L>#tY#2PD-~.;^(|W?aPM|r@>Q#'*3xJ,QdP/H&,##OK->u$###&##n0*oQS)##[,%#c#.SS%##tG$Y>#SRS###5H%"},
{238.45,333.93,8.93,6.12,"D>#######yv&,z;######^I#*RTW>$###gwB86%<R+5###oYY,%######Ku#BnY######H##bnYL5$%##B^,@m'?d*)##9ROr,&.,####$##q[U######)##dsYZ[,###ke'PT*04I###8sATu%######%##_K2#########Pn>(y6######W9#dnY######"},
{555.13,334.92,8.28,3.34,"######8,#G,$######(##MJ+#########oa?###:#$###PbA}Z)###$##Wc#0=H&##Gu$4d;:Q&_/'2L6+hT###[cAZ,%hV9yn0######*c#OgT5,#wQ(Ym$b:9p%%zyQB5K###M8+[-'@hTy5&###Gd#%y5ngT###`>#W5$l]S###K?#[z9gY#hY#.##VZI"},
{555.13,334.92,8.28,4.35,"xQ)######5##R5P'##c#&j.#e09~$#ucR_6#TcR>##/d)o&#0]3######[5#jdR|k#>[&Z'&G05<c$JeB{};ZcRgY#e5$j{*bR+######R>#~&G_K74##iG#`V'UcRpH#rx0p^)YcR'##Tn'T5$######%##>$$iZ(###&##S%$<h=###FH%]J#U<F)##l?)"},
{520.53,343.58,7.62,1.36,".##OL+z&3PT5EB&J3=1f+5'7/Q#$I%mpWwu'###,##uJ2###Dl$_]%cmW$##>1Wg&,&V=/##a@+Cn$cmW######5##po5###SH(h$#cmW(##gmWK##FiB0$#b$+w##kmW2#####2##1A0###qQ)}##/nW'Z#:[V7##k^:eI$~#&]##pmWG5####2##6[*###"},
{520.53,343.58,7.62,4.64,"}n1######)##+&Y[>#6w-$##_{@TZ#I%Y###w%Y0,#)c$$##*]2#########y%Y+##%S.%##v|G1##p%Y0,#k%Y###j,&(##(f1#########K%Y###Jv'GQ$jz>###`/Jxx/C%Y###BZ$j/(?@+######&##5(YQu%2##8Q$GL1k1>m.$T_9$//#sC:##Cg,"},
{492.17,346.37,7.03,1.53,":v'j#$a.X###:RM&Q$c83###GQ&B5#X.X#########:n-###eQ(###4/X###f.X###NL93##bm+###M/X,#####$##j/3###u5&###%0XnY#S.X###6_93%%gR-###p/Xa#$###)##Vp8###{k####IQ?tK68~/###x&,s0KT,%###}0XB@(###&##eU:###"},
{492.17,346.37,7.03,4.55,"'/1######0##I/Y+-&ZP#/##Az4-0K,d)###jF8Ii?6#$###]x4######9##l.YZ#$)R*@##(z8Wo&<.Y*##F/Y9Z$*6'H##O97######+##V.Y;##?%.7##T(=p##?.Y4##<.Y$##o#'r##ax4######.##|6W$##6w,)?#0p7*##1wTRu#<.Y###-c$<6#"},
{558.84,508.32,7.50,1.34,"y>#5?%W>$0,#95#s.%{[,nL5^G#S/-26&YBU'@'pl'fH#vHR(o-1c$[P#b>#%@U=,#p-*~c8+w)i0+E@U*JL-s6$8,Q'/*f2`S0###A5#n#%'DUQ6&X7.0Q#f'30g'~@US?&%AU/m$Zy7VR%pG$J>#Zc&lY#Ej-JK4t6*###Cy$d`By18^[*|[,yB/*v&lw)"},
{558.84,508.32,7.50,4.62,"HH%jI)J[(~(5#g1G?'[B$B6R2d(###FV+0]2mG$T,$<#$$##jL9Z7'$AUQQ$?AUG?&lp0>n(Yd+$Q#5EUou%'##6Q$wn.D>#,L1^-*p49.~+[@UM=<M%*_f+cR+D#6@@U6,#dG$_>#[%+Wl%|-$],MX6&l,&u,&ADU7l#s7+PS/v)8D>#8S$:Q&<5#I,#D6%"},
{439.71,524.54,7.94,4.52,"av'D7'lJ1Rl%o5I7Q%P6%5RNKJ.=,#{K/ci@<%*}Y$$l#;5#u'8|Z$xlPyP#?fWl,%&T-4p(Kf2.m$ZiWtm(5R(gu%a6(#Q$zd,F7*Mz6a$'BiWFeWJ?&]R%711*;VLFDi5#lZ(gG#'$'P@(@5##S+ID0&[),w(Z4Io#$/a7Uy4@81###b:&`Q(.,####~H#"},
{372.69,536.52,7.10,5.38,"EB*pC*/SSD5#B4;3f(|iD###ZV:.,####%##Gc%6#$###*##`7,Ar(YQS=-%GSSJy(j^;>m#`GG%D>###*##wc#9K5###$##eY#TD);_6p:;;g6Ip14o/XL0(_3;RSO5$Z]$?x+_@.KZ#^;1ZP#$##oJ*ln07#$Al#O*;(S+,25Yo2~l%eD(WF2m-+gH#~i4"},
{919.35,557.26,8.01,0.07,"2-#s{Zwl&oS0Oh&55I95#5,#Q@*%K,gY#d,$KZI.c$|5%Uu#Bn(ixZ######U|ZoT7###t##jaG+c$###|-$lwZ.,#*##H.$0'1Yy7######PzZ3l$###9##z6MGZ&###u>#VY@]]4eY#o,#%x/Jl%###YG#@zZ.,####*##|N?lR*'Z$~P#5&(&g3I-(}k#"},
{264.33,562.65,7.18,0.13,"_>$O5#eY#$##.Q%WK(du&'##1EB=&-.,#52(w-+]81'##r{S%?#:6%*R*###i082f$<aF?,#b~X<S(Jm)bK#4T2UaXIf/to,F.)1c$Cv(_#$Go3=,#os<{?(P_X'Z#6%)>7(K^.P/-p1/B~XW>$###F5#c-'W>$###:m#^I,GN7tc'J5#;R(Ng'wmU.,#jZ'"},
{440.84,588.61,7.21,1.46,"NS+Ou$K,$;n*@~)$%)q>%c#$vc(t?$i$*xd>fu$H?&YH'.CRD`;aZ%H/20l#]bL/u#u$*uu#jmUG>#G01W3,[x.Xm'xoUtI)~WC3,#n05o-&I,MRG#k?(KZ$qrUzH),&0}P#JB-^w-kpUv-([[*`-)~?&kH&u@-PG#dG#i?%#i(k]6$R%2c$Yo#~mU{l%;l$"},
{320.45,605.23,7.50,1.36,"LH%<-'G@$)YH:H%5Q$):+}6*Y>$G[$N[)[=>$c#1R&)['>?J+?Hs.*ET1CH$#_;###]BObQ$r>O#Q#zV=l|1pm*a](2?O_7+_>O?-%ax1(q(J^8j.%Q?O'6#eBOSo.+`@6##h@,Q%($@O$%(,I'{|CL#$Kd&^5$Kz,^97'##`e#8M;mJ/pb#97$aq=Fd'S?'"},
{320.45,605.23,7.50,4.34,"ge,57)n$%mJ17S-E>#qU%B:8MT4###Vu#Pg+iY#k5$q?%/OB%-O8@&zR-,y%L)A9?#K/O6T*N,O|##go3*^&xI,%`*C1;f~+nbKqh;D6'2(*Bz9#H:9YLp>#E.OQv#{p:-##Ux2-e#9?MC~'cc%2-OQ>#R7)P[)=z6nY#S.'0s2.7+I6'6u#0n&iD@AR',~*"},
{445.30,606.32,6.73,1.48,"&I&W6'Z>$bG#q#'M-$N$)'k5X,$IQ%kl&*)Si6'.c$cm#9OD(g6#l#N?'|P#-mSE>#d-'>N*Y%,x?'zoSQT,WmH.$'6T/]#%5_</,#H6'Du#|qSR#%_?(,Q#ML0b$)UpS&$&.uMFH%X91<8(9I*E>#VG#t5$Y<.{R.z>$J,$*1$8mS..)-c$aZ%RjAo5$@6&"},
{553.81,612.31,7.59,4.41,"F9,^c%QR&D,$}Y#NS&wZ'_?('I%y7*G&)cB6Yh1yY$L,#$J+]f3yP#X@R_c$)U2zS*I'7h%,2KVH6'h''_QL{J1###A2-*E=FL:p$$=YL(6#RXGw?$lIV3$%'JVBc$$(5mw&K^9/c#u9V~I&d$*^H&a6(@I&#~,t]4@U3nR)BLVAdSml&:.&U2:C)9b/2~>#"},
{369.95,627.22,6.93,2.72,"bP#I')3kI1,#6,#8<9mP$###Tl%lg0c?(ic$ql#]Q'TZ%{p5c>$p2)8g7f$',A-U`1d;?Y,$>wTln+bR*v&(2v%,*>[:3D90>#$&?#S7*&z0Pe/Q5#M{S4[(UxT@5#/f)}J-S{;%I)Y<,6wTp?$Vv)2c$nP#-u#s-(Sh2*H%SvGJl%8R%#T/M+2C-TJc#M06"},
{609.31,635.50,7.65,5.40,"uA%DZ&Hy#%e-*8&U-)y'%3:9oT.l1*=p40K0<3BV9(%e-?'(H@$ue+bMS&~.i;7LI)2g29[(~Z'ql#u&.&JSrQ)0##CZ%pR@>A)kW,iHSH6&9JSWe%E^9~S${K8###%##c0/OZ%D>####pv'L5$AK$SD<-W9of6B##Wc&Aa.Yl&###&##z~+2,####)##jd,"},
{858.83,635.73,7.71,1.85,"T,#hJV######[v;@<D######'[:yc(######@n%bQ(######9l#-I*######(qP)c$###`##^NVSl%:#$r##bg/WH'w>%XG#wb#.,#%##sy3e^:.,####C+0NIV###hP#bh$G%-%##CQ$Q].95#######F^J}>&gP#F>#`fFf7-+n,###a~$>['FT5###{#%"},
{564.93,648.74,7.18,1.31,"_J-7.)%e,$x/U@+xQ%#~,RH&Pc&l.&>7+nuAEZ$}H(al%`oU/-PEZ$:FGcl#;>M#6#,jEj,#~mUU>#AL7v+3uS.).'AnUiT2Ix2_d)Py1JJ*@NAt?&_1;6c#qrU-v&RC:c>#z]1gQ&6oU)%)A$'u-*bP#=R$o,%uP$+f)=I'A|+%I**]*P5$2q&U]2.7(WZ&"},
{881.92,651.84,7.83,1.33,"~($Un[######j7=O]5###'##Jq6.,####J##5R'95####-##|g1Ce/###+##}r[ac'###[##{&WW>$###.##K%+.,####(##:x2###)##KD3dn[#l#.,#+W''[LY$*95#F##dg71,#E>#4##6#$###<##ns[lw05,#DQ%`m9Pq<yG$qc&=@'E^9###WG#Wc#"},
{512.51,670.70,7.29,5.53,"q$$g4<z93e./z~+fvIFm)Ku$LZ&Md%<n)Z_=6#$R>#A'0S/0,](*7@HK5Tl$AfOZ^/|u'P6$up82q7E?&O?'w>%5o/W93F#$p>%i5#w'3TtEaWC%I'j$(*7@#@)cZHB6&=~Lu5$UL9>v#-eOD>####9##p{@nR-$##2Z#fE=b&2$&.1-&mgO9y&vz=DI$1eO"},
{336.68,680.15,7.90,1.48,"GR&7-';?$gC7_-%Kn);$&m,%LQ'8J%q$+IQ94Z$[?&9R'$BM['3E9.N)=?7+9NAP,$CE8aH%UdS4c#$18vN-r%-D&(MdSo~+b*F|#%'nRd.(E}H46$K4Ebl#+iSwx1S07[>#k/-wT0ndS+@'c-'.V:qm(ZR)jZ'n@(um*:Q$#T$z&50A(GZ&`~$|iDLc$yc'"},
{495.89,681.36,7.50,3.09,"'##%d%|6'OG#-##|h)nl'###]5#h{6-Q%.,#T>#^}E`>$,Q%fu%]r+sV<>u#}v'(|,+dLoG${uGh~+i@*-q3P?(_7*]v&qJQz7/7y%Xx2Ng&Jm)#c#LY5;g1KIQW,%88%sb:J;:RD4(6%g$J8w(tu&2u#z[$x,%We)kd&JZ%oL:P?&OH$S`,^$)516<5##]F"},
{495.89,681.36,7.50,5.54,"qD+<&2qG#bD;Gw#hcFB14Ux3S['n8P}H)g#&B,$>R$de,&q9bM)#aD*$$#J-%L*]9LQg8{G$}8Pg'00H&FQ$@%-u7-qH'p?)2H%I>#1J%Xz;&I)l?$5N;2>D/?Og@*#R'?+7S?'BGE]#%z{;EZ&###Y,#ro3D>####T,#k5Ltv+5,#e#$T}@T@+L%,w5%T'N"},
{399.10,689.20,7.25,5.51,"L?$#7'ANU9Q&0<>;6&F]2&@'vQ)<5#[>#m`@H?($##j5#S`8V8+zr*9IUpP#3KUL.&lWExl#]:9$r<*##nu%ZbAld,(##4I%/c$@3,AkG?K2M1;u%)ny6vz1<f.beN=v(jn%+Y9B%.'##Yg)###~G#7-%[]5.,#1##r<?`$*'n)b.->S,sA(/a-X97zb#Ko)"},
{564.84,687.22,7.17,4.38,"###+##<6'wl&H>#UG#]c$Jz3Ku#v5%:G2YA0d0+;Z%Uw%nu%;@,&##d>$Ix(en,&y+Up8U(1S7O.w*8t/@f/f%,H>#%:J'6&A);zb#7?'D##0N?#.%jSW($%NSWO.(.}DI.&Hn.mR$qSWAl#I.&w?)G$)###DJ(ikK,00&$&bM<+VWb-)bH%C[*M{.Bz<j5#"},
{342.33,699.54,6.99,1.54,"Oc#,e(Lc%zb#g,&3%$im+8,4nG$SH%Rv'E)R4Z$~#%f7$mK6YJ1?u#L9.p#%umVD5#5.+*a+X%-&8'fnVXC/:V6hd(3B/J,$.D=C5#H`;zY#3sVK[)Y%./Z#Pg/40-XnVol%U3EZl#:049%&-[)Kv$pe0eG#Kh)s/4KI)fY#ro#Y?R[@+'H%zl%sz9k>$KH%"},
{342.33,699.54,6.99,4.70,"`5$<?%~Q%r_;]R*u>%io#q$Uym)qb#}U)]x3N81rP#`v*4d$`K4L7&4sDp#$+/Wml%Q^/;9.Gn-*Z#W3W~m)(E:&c#?M=5,#P]/I,$ih7AI(;/W}0.q7-u7'D7+Q3+B.WC5#rx,+6%E81N#$].$aT5+Q$g,%%@(0DOyP$TH%:7,_Y3}>&pm#:H&|b#?Z#8e("},
{334.18,721.28,7.91,4.71,"95#mP#<#$16%0,#?$&|m'b)??$$h%/TS$uXJ1~'LH'|u#Mm)WS1###qb#Re$ex1RS(u^:yh4y7U1['I<..RN:&0>5#4:L~~/y(9E>#bG$^5#~E>#-$|mP4[&%8Ur~)=N=gn'en-'^&78UP,$(d$$Z$q5%D>#}]&xx4Jv&<$'sK3D(Nr#&hm&$@*PO/U7.`-$"},
{868.32,812.68,7.87,3.78,"}c'D~*ow%m7/w_7jm+K>#%##_T2QH'###j##yS3######y$#ac%iy+g;Tcq=,y*TYJzz3U#%_:T7w-A5#,[#uK8######~7$rZ##2Q%cGwG%{$*$8,W?<)10d7Tqb#+Z#cM)-r:######rZ%s>$#&/g]*no1RH'[d)j6)1P3#T,OB0DH'p%)z_''94.,#$##"},
{685.98,825.11,7.43,5.55,"LK%$)=S'+9l$Sv'Gw,d((zlO9f3###Q##BpJ-w,N#$SG#^[&Qq(fmSj^6Qd(7rS)*C/c$T-$];@######(g'|Y$######HT,Z6(qu$*BR=t:0nS8#$Y5$'X-.[J######+?#H7&95####wP#B,$gl#33;5h/1h8/,#.,#&&&c11######&##&J$ZP#######"},
{643.30,827.98,7.13,4.42,"1m#K|?=H%:#$Dm=bJ21,#tG#&{5Tc&###E##v5$8]3###2##e,%fR(cN0up9xIQz,&mQ%Bh)~M=j&)dn0Zv#QQ#(MQ&f2###Xu$}*5Ny0K6(aLQ_7)<?&MI%513Xx(EXFGl#>?&GK(SIQ)Q$]P#Dl#&K&n|D~x2L>#V5#ywDY?(G>#lA+}o/6#$/##958Y<B"},
{629.84,850.58,9.50,3.94,"Q.)qc'8$$019'n(_d+{P$kY#kA2G,$pb#~-#ll'######1R$My1n{16(Pd%P9o(aEA,QCp,%,'PsQ)B5#}$$=[*######e@(d~(6(PPj@bP#m^8uM9DX:LT'f5L#~)M#$J;)pJ)vB6###c,$%##9?&ag+%94Aq;;I*N@&]QF^/2bM,UK5Rv%7i2b>AOI,MQ#"},
{629.84,850.58,9.50,4.37,"~$&Yu%O5$###wS0{#'N#%*##Gz:D>####d##%Q%######.m#>d%Uw,dU+46Jq/)mYG|.)0-'BKJG7-4,#%m#~R)-e-###x5#nV4>?<%1O$[Hq~/X_0c[AWK3_.Ofd'AR)-0&+R&,.B<C:W>#0I$Xa8?A+dQ(7@*sT,Qy+=);FFA.v%(R(:V0Kc%=[$C.N0~,"},
{904.51,847.84,8.06,3.33,",u#:,#4Q%|&,.,#>6%dZ&d%-7##aJ06-%+$(j#&CR(I#%+##oH)hl%-.,LZ$j7V7$&mT*F<V*d'*9.F<VQZM9f0~d*f7,n,#.e*cH(c5%?#$F<Va8V_v(f[&GS(x)U$:V[Z&@9+$f.Yq:[>#>5#;5#Cc#'J/G?##03;C3~l&wI%Tq8vf2G[)DR*AH%w]2`g0"},
{904.51,847.84,8.06,5.72,"###%e&eK5;5#^n/+m%Ep2ez.]S.yu&):1P11/`9M-)%I(Rv$`m+%e'YR'|Z&|SXbH&'A(1sUH@+qw%#XVLUXzw.1##}w+=@&BI+###(w$]$)IXXLQIh6&FS(tn(DhM(?J6@(*c#7Q#Ln,uG$b6$VZ'F,#O#%sZ$e[,0##2R(-H$@%,7,#w6&@5#4H&47'Q,$"},
{388.16,872.52,7.59,4.67,"BS&Zz9pb####t_*Y08###&##;v$Mf-{k#&##@##hII3l$###|5${g3qz.o&4CyMAx176$zJ)Zg52^'3_=3##q,$Gb6.vMnP$Y>#Gk:t*=FJ0J$?2i:Gv&N@'F'/8T-AU6OQ%sb#8d#KyMif42,#.0-_g)C]4uV9l7,'d$z>;zx3ZP#T5#*|2######X@$TkI"},
{376.75,891.63,7.32,3.79,"Z#%u&'g;TIWBsv%[o1y~0jY#e'57Z%###4$#Jm*###+##W[$*[$m_Rj7T]P#8&(#q4]uIu5%m:T8-(&l#o?$YB6###+##2$%Zy'SkG</0;5#Ty3|I,n(8V39waG)-'eP#Pt2Y/.qI.###I,#`n,h6)1^.lq9JPH=Q&H>#bB,N-'4*:M5$ul${-$BxRk>%$##"},
{902.10,104.07,9.12,1.44,"{MBV6#D%Z)##(##S^,UB7###=u#h@)'(8_m)F>#I,$T>#N19YNE}$#E%Z)##1.(9?6q%ZcP#TUXfK2H?'Td&eu&?,#<?&Zd)h1>/$#`%Z-##xZ(?Q#]*Zmf2$&ZTl$^m&YD2F6(D@&-6''##4f3S##K%Z+#####8##{8/*[)c5%zl#>?&E@*;5#UJ)A,$###"},
{832.10,138.75,8.88,5.76,"2_#=#O######PVOB7-###%$#]lHll'###5$#nI#If4######`q/T7._5#+v&UaXQu%###N6#E^XeY####9##=Z#<m%y#'###|v,###@~$NHHq[X######2<(5xT######M##Z%+Y>#q?'Yu%$m(J,#6~*]^,iXJ######J_$1M=######i##lw-###$##W>#"},
{1046.18,160.45,8.74,1.69,"dl#VV1AC4.y6Zq5U?'T5$}956<4'80,##0?H7I$P6)AS&TYL[o,:[CGd(2Z$QxP]$)~H(X5#ds0T[LS{>,l#GzP^]5sv&l,%<7,M#$|1*pp4_uP$##gH'`~&a/2Eu#`u<:q7lwP1u#w,$#_0CQ&#l#sV){]2LM=###E$$;m&mx1a5%C##8S,nS/Tw.&##4~)"},
{1045.81,184.55,9.38,1.32,"W^36`;&6&1l#3OAAc%a,#p^5;o,+u#m^$1%U######~_%@$U=OFt,&@Q%3^'{=<@;4?::`Z#l)UFe-2f-y?':-'iG#9F>vl'NK6%##|l&cA&YA1q#$>ZA_M7K%UP#$Kw*u3.~e0t>#^RQp,$[S0sb#d5$]c$&S+E-)P##wx0faA}A1Y,%3/(zu&-o))'6)##"},
{680.39,264.00,8.56,4.44,"uOJD>####D##V-Qy7-###P##</'jD:WI-###9,#S8%[,Q###=-QD>####5##:-QIn(3H&Al#gg3y.'PK6Wu#Q>#,T$)#O###V.Q+u####0##buLz>%Lc%T^/qc')A%7-Q<g.H,#rJ%fFK###G+33'7###$##Hj/j,Q.,#k,$fZ#_1QC07U,$M##jU/}I/###"},
{680.39,264.00,8.56,5.44,"ic%l:-)A+hH&2##be,8F/.,#%##g>$j3.######z6).@'###$1#Xh8_[,###`$%tg79'3;5#i/0,$(Xr)^$)OG##Z#KQ9@c%D(#>g8,u####@2'C'7;R(l>$:>Fxx.'7'8[&S#%K;'kxShG$t##d%.######OL#}{C###$##iJ%XzS.,#&##%##540VZ'###"},
{56.46,272.76,9.06,0.28,"s^$^6R.,####]A*qh7w5&s,%hI-4##vP$3+:Fl$Jc%###N6Di,9*6'###(##`wI6#$###<v#;7RF5#/l#VS&WS.8f,Tn(Z91(%D######(##^:R######$y&(7RbG$:,#Ch(ke-1w*{d%2U4H8+#########-;R######F,#):QD>####O7%q#&Q#%###yw*"},
{467.69,329.16,8.20,4.33,"_}L######;%#t6Stc(###^g${?'Dw,*$(?S(gH)1Q#'~.g5#H6SkG$95#&%#BoQq8SI#%f##Vm=;6PL5$1##*g5ic%z,'l##x>QHc#=6(l$#D|CM_+ez>,$#a7Sjc&pb#N&#lh?H##`$+I$#g^;~##.A0J-#HkL>$#vy:P.#06S######5'#$x1n##Jx3w##"},
{870.20,514.25,8.65,1.53,"<##qM5L5$###0I%z]2n?)1,#>5#0D*cM?/,#c##DeSg,&###j@)*t4I+I>5#Y_Xj~/1f.ec$`m+{Z#[^XZ,$%##ue,2r@###+[)D,#ySE903:~Xz5$w~.0{0{w0n-$@]X8,####K##<5K~#&######nd#O07a-*<c#s[(:z6~>$%[#XJLN@-###)##sT,CGN"},
{870.20,514.25,8.65,4.59,"d0-&aF###1##>9POR,sb#I[#3x+po1}c(ju#>I#kf5######zOJA,$###b##ASY4##]7.L7$Po0;*1*SYVQ$7JBg'5(m(1##t:>###1##Wp1UTY+H$TR,wv#ue.b$%1VY+&/,QNfP#<w)m<2;c%###5$#,dLgg9/,#bP#3_(Xc&1,#:n&TC73l$###I,#|i7"},
{472.55,525.67,8.79,3.18,"OG#^m#q%,K$'}Y#~W+Z~.M5$HZ$m):SQ&3$(Zl$s|C)l#202tb#5K#/nStG$_&0UU(w6I:-'InSw6(@v%foK2@*%/*m-'~qSD>##.#vmSp]/B6(&##Dr(baCloS~@-(@#@k?zx-;fIdl%{%-=5#$H#RH'2.*###$##Kv#+%,8u#B,$6d#ST5+l#~l%)R&EQ&"},
{813.15,533.88,8.52,3.62,"x%CD>#######PR>>m)OG####p0&So1W(:.,#n>$`u%?L26$(<gO######*7'O)T.,####~,#86GX>$<u#AQ#*B2RG#N[$;z9tR-$##)##byKd$T1,#R>#B7?qtL######^-#&C21?&>5#J#$uP$YN6mP#2I=JB-Ca>=Z$gB3[V8{k####=,#@R&/u#######"},
{304.43,566.69,8.92,4.58,"x.+uP$K,$/I'-u#w#%W@'YtH&d$h-**8$&GLNJ+eY#h##2m(9aE1,#1n'KV0#:6cB*z`C(4=/7OU7*(|/t*Ec-*/,#WH?v$+Ft@lv%WjC7Z$J4@jI&37Pj$'T6P99,mD=/&'I$)cJ&B6PlP#kA*yK/nZ(Ol#dT'uXG3@(VQ&Z%-E:P:c$Ov%2A/?(0o#'{##"},
{735.60,581.34,9.41,0.47,"kG#?Z%E6(zY#%-$kV7#<>`?)<iU4T2|d'|n.YFF$##b+22aAeG$o,#,@+bG#Nm*O##~91Z04JfU9u#b%&~(9eiUbJD/?>Q02Ol%^m#E-)&##$D;J$%(-'5c#j7*3;0C&1Y5$50*o(PFNC%##bG$i,#3Z%k,$pl&DH$F6(#H$'##yw$?h=:5####zd$5a=fc'"},
{735.60,581.34,9.41,2.45,"3.UHq.+D/DrS#l#42)0(MA.UyY$,##$0)5HM%?%OG#%##^q:$$L1rRXZ'c6#993g?;?S.Ed'Q5$&m%;e(C`AK>#K#$w>$)p5Rn-5q*?kB[S-z`>#%+&##LJ*J#$L909Z$bI-lG#gH'$c#uQ)OG#}P#Q*;yw.|P$|k#L>#,K.(##[5$gG#>~/###`G#-Z#?d*"},
{980.68,592.04,8.53,4.68,"PG#Kx%?]2iY#;S-n$&)R)+8*Bq;D>#9##Y*4.,####6%#o<E4v'RB*u7*@v(n:58|@Fl#xB0fpV$%,B,#R?98S-#]*~^(noV0,#)l#B0%HQRb[++@)W?#RpV7qVWf3&##n>9w8,QqN^?'C?&D>#*##_u$M.+x#'}k#%##.:3+H%9@'###Wd'J>#U%'7#$P>#"},
{384.62,609.27,8.33,5.41,"PA%QB/fnHa?)#A*}4E#p4uG$*Q%%K*`o,Zo30?&56&ho-J%,>q.f30cjH)H$RdO_f'g]6/-$PB5Pm)>J*nv)~P#;u#`;4v,&}G%Zo(h&1pC:hM>Iv'wZ'm==H.*C06.Z${nI1,#<c$jp(lOH5d)A,$<##E19G7-E>#]>#M%C%K+t%/WQ&mgOte%#~-6i*hcO"},
{792.27,619.93,8.90,4.31,"'##_29e#&###dm%Xh<95####z05nP$$##$-$Yd*-u#Q>#n#%+##u&0|V0W>$m)15C9/I$i6)-BTA,$A##G_.uZ((Z$n'&yz9/##De+yf)T,%nw-(J.DS$kh:sDTp':=,#X9+|x*n@TC0*^7,9,#_'PCZ&###~I#Q@T<5#iY#7V$N?T###2,#D'#X?T######"},
{444.33,628.95,8.50,4.42,"h]2D>#3##e]-<-(<5#g-%&YAY#%;5#(_%3=H7@(O,$n$$s$*H;@5,#x-&/q-n_=)&')D=A;2NmOk-'9+5*)<1c$?Z#ypOiQ(ik>}~0G6'wP#B}?sw+/ZJal$dmO)s5QM<P-%Gl$$D'GlO1,#dA**[)g,%Um%H%%EmO3@)^>$C$'RoOpG$S$&6u#NT+EI*hH&"},
{324.48,643.93,8.92,4.61,"UR%=d(7Z$m#&6u#@.',S*]B5+J(Ww-?@$W>MO(8,u#rG#M'..T1`P#cG;qn-A;:p8*eV?1f-}.SGm'bh0[uHb5L_G#0j9+^/7f.K]%H.S:u#8F@{H&V5J]I(,/Sj<92q8c~'u)5RD/~FHL>#yZ'2q+`R,3n&87&af24n'*I(0~+V0STl$5J(rI*=y0h>$<e%"},
{233.64,650.35,10.01,1.42,"F]-;qP>F>y-&$]%ymP'w@DmP+G<ylPM?$*9IG%+C@*%D(`_5h&2-6%hZ&jR%hZ%L),d6Gx08]II6A,RH&*i3?e,a.'7n+.-%b/1C?'TQ'1##V#%Jd&zQEbT3*<@z>%yc%5s387*Al$@I(BH&}5$a~1.,#'##*##cB5}Y#n>%HQ#9T4'##'?%A6#I..###/,#"},
{233.64,650.35,10.01,4.78,"######E##^h=*##,H%B5#Np5Gu#W,%%##<f/.,#Q5#<%+,.,gG#e,%Pd(Bv).d%{23a;>SZ&v?D(K1W#%qQ%+-'6##f]0Mc%jx0d#&kv(Kv$u,&#o(EnE4S,*7J$QK>v%ip+p%)m~,1'4TR).f'_<CZ?(Qm%0Z#$40m,KPy5wr,d9RaK.c7R-;R2A0n6&MgL"},
{458.87,701.14,9.27,4.57,"=Q#Xo/|c(###B#$''%4^6r>%hg.{R*4-&>`;Fi:95#uP#P1/@~,[?'6~-6Z#Zy4@7&%[R`5$w[RXu$PI''dDrZRB,#F8,O+50Z%$##|~,;%*i&4AZ$v~R<-&]_RxXG-e)b[&JTH^)>]~-/Q$######0A&}R.H>#lu$8-;=$);p-({<8A)V362(/,o1$$&}~&"},
{248.88,722.11,9.67,1.92,"5p%(CX+V/^IVF`7*DXAH%K0M4r<f7+eY#RR#n-+QG#.,#l##dV>FZ$q['3AUy-(vL-=W>,%N(yQN~-2c$xQ$mI-E6&###2##_r2Jd)cl&E##<6&?$%8j9n'6`WBT6(2H$Hi2~Q'To-###(##2s=Xu%jP#I.'?6&Bw--##Gc%@5##g2$##Tu$###cu7######"},
{853.40,790.54,9.50,3.87,"1q/N.-Pu#ZP#~kEnZ(###+$#[U;######c.#(c$######~H$Gx):2:]oHw#&HDT'e-*Z#sd$8XF######l%%OG####^>#lQ&P-);13q38CC-yAT+n+56%N`+#-9sB8###Y#$mH#v$,b>#3c$'Q%{>%o%+YlAOQ$OO4S=Idd&Yv<evN4-(q,#Z30Xe0###1##"},
{669.13,840.69,9.02,5.61,"i6#VbGN5$###,.#]f30`18Q&:v'GQ&v]&#RM+C8<5#&##%M-An%Br7v@-[J+>1(:ZN'J+/I'%]N@y6`5$Rd&180######f@%lR-G5#3N0D^Nll&[#%`.DP3;YZN8#$<u#?z(113######6##yG@>AK8`7z:5|Y$zf),ZCuR*sw,iY#F>#@v$gn(#########"},
{419.77,886.77,8.21,5.58,"e.%l>Hc5%M>#'S$5^5pA,C,$'7)A6'jR%C{9h7/###%##1h){-'jd(>A-AF@<x(z>E<p6/['Q^Rm&3rb#vc#&W@######<%#tP$;$$]_Re[RRl%I6$G.OwO@+[R'l#wb#La/VaA######>##RH%88D1]R>I).,#g,#{&2Nq4WJ0###8,#|[(L].######%##"},
{659.28,892.73,8.41,0.58,"Cd(Bg)I81w-*-:4eu%]u$9>E'y0X,%yG#cJ/VG#c#%q#%D>#bOEeZ%Dw(6&-NxUaH'1%'}{ULv'?A(}{UrnSt,#1I)`7-.,#}e/-u#Y.'58/}{UpwU~v(k[&xI&}{UtdT{5%cQ%nm(Dv)P,#/,#C5#=~-]P#{>#v[-Z%+#H%hH'3J.I>#=-$i,&H>####D-#"},
{659.28,892.73,8.41,4.08,"###l##lR-RG#e#$pc&Yl%un/5J(.6'd>#z7/gQ%rv+gY#&##Cd*J##i-*g6&q7Uy5%c[&c<W+['oR'c<We8W-@%;q8Sh6bP#E%-%##V#$*x/c<W?,G1v&+]*g.'*<TN9WV?'BR*d/,uE:lH%~Z%mP$'##*H%bZ$km*(-%~u%-##My.Cq.$[)w~1?n*y$'b92"},
{659.28,892.73,8.41,5.63,"&##gQ%$7+gY#`P#/@%N(9Gh6oH%}6*~/,q3;;'1Fo2j7+q$')A/`n(hu&d>#W/W,['eI'W3WV?'e.'v{TA0WiH)-##Lg/>_/sn/b#$~Q&p#&W3Wp7Un,%d[&{S*W3WW>G2v&Ql%VQ$wv*MH%RG####o##J7-q?&Cd*K##=d))x/O.-$##Mu#y>%^Z%mP$'##"},
{399.41,939.58,8.61,0.89,"pu%s?&%00Pn.C(5|Y$a5#UW=v~,D>#l5#,~,Fu#bP#<H&wb#3|>56%hH&*8,1fVpl&s%''jVV[({$''jVa@Q8d'|Y#O%,5,#VT3QG#Qv%X80'jV_eVd6('~&/J''jVHfV9H%Hc%i[$/&15,#;5#0,#=$(*Z$?Z#/&1s[+kY#mH&-_9<5#@5#M5##p/?$).,#"},
{399.41,939.58,8.61,1.69,"@$$$03OG#&##XH%hR-;,#:%+%Q#+[)fu#|6+E5#,Z#(p395#{^4Y..ZG#{x'ioWav)o[%PsW<v&,&(dLLknWaG#46%%q5ZP#Rg.Ye/S-(-~'PsW_oW[c%GS'0e&S*W|ZLaQ'mG#W/)^7./,#Fc$~,%0e&E(9G@'QR+_>#Ew(sY#}7+O6$SQ'UG#$Q#r%.K#%"},
{399.41,939.58,8.61,3.54,"dG$A5#-l#Z10OG#)##h$'#C4fG${Y#/n&.&1)@+:5#'##RQ$kd,i>#46'<-%owV7-&qe'2ETnH(h@&gDR_eU8d(z5%b.+7Q$>%,&c#*n+^G#?|VKxVD7+gv%K~'?|VIxVw#&cA+yn*1x1,##3,#t>$#@*6#$R6#p^6Wh=/,#mS'E_8e7.Ac$CA/8Z$57+`x,"},
{399.41,939.58,8.61,5.75,"$##wZ''6&RG#;5#g6%;B2XU646%UH&*K/Q)9_R),8,7x0Z-'0&1kR*=#$*Q#p~WhH'c7'8aW{u&N~&8aWE]W>u$1##t^4Hq1TU68u#L>#ic%8aW.fTWc%d7''n')WRS[O=I(X>$q5#(%+u['pp2Gl%###X>#%n%So3%##?c$-Z#Je,3,#6m'1u#}l'F>#fG#"},
{281.06,1058.96,8.38,4.64,"D>####8##5v(_5%###}$#yuHc6*###)##jWV'Z$.,####jWV'7*###(##Gl$@<@O#$xQ$$h3+TVO-'J>#BX-Q1.ZL5&l#7C.g&4###&##z?$7FFtP#[x+cD22TVjl%h#%}@';3=)~*A&(?_9;<?/,#$##>,##@M|#%qC:o#$ASVZ5$F?'Aw#LOD7L/&x,Id'"},
{259.18,164.35,10.99,1.42,"9@W<-#[`DE##+##-',0p7###ru#_&3######j>#@H'######W@WkR#7>N3##)I'H54kvV@5#`gTiU6?c%U5#;[(%6&.,#$##;@W+-#$mSG##W$*8m#|CWid)qAW56&Xd''&)~$'A.,.,#$##.@We##E<E<##:5#iR#k;=qP$>o(Jn*Gd)nG$E-$TD?######"},
{259.18,164.35,10.99,4.79,"+u####P6#=i;*U5Ol%Q6$yd'PuGa5%sb#T##ly7###|0T###D>####hv%dJ06~($x,|/TaQ'62T:7+bH(Xl#p?S###=1T_G#eY####F?&,[(4H&J>#X'NlK3W-TF>#07'_t4R-T###s/T{>#######R5#Zl&######$$$?J/b$+###/##<T-BvR###);::##"},
{207.65,207.07,9.69,6.15,"Fb3}I/###q5#/D%h&5###&##@.$NYA_?([#%y#$+8ISx-;g7iD8m;>E,$|w$P-Gi%/D>#gI$^/0)R)Wu$R=1pZ'wY#t[&2oSC<;yjD78.ql$|oSpZ(###s6#]mS######_R$A/095####mZ$Ao3%##Er9{0)xlS######dC#*?Q######Cm#^w/######jm'"},
{92.82,256.02,11.10,1.34,"A7&U(4~[+dZ'9%Q&p2LZ&eI(sr=2v(_>#;/O:T/$?%_n(~(;T-)ad$k(:aZ%j6R~u$`,%)S#e-Ca2=-I)8w'P1/(^2%n*A$&w#'aZ$NG<.'5mvEIp7B?&~l#yj1{84v#&1S(1^1z-*U#%$R%$##R6<+)>T,%|I$U7R<c%###zN.x5R######;g(z5R######"},
{740.83,256.61,11.10,3.07,"<#$*##Wu%>l$:i:.,####R?#`[R######?:)4-(######5>?95#$##x5$HB5b5PG>#K>#pn&#o[95#?,#<51'.+wb#A7&=BP95#0##?e*h-*IeTfZ%t$+nY#it[rJ3C5#fl#Wp+9o[Ud&PZ%]c%Vo0Mf0FH&94HE6&Po2[,#)q[$m(###x##F8+_g:###l##"},
{397.36,285.28,11.71,1.42,"zw0W>$######N1h###6#$=##]3G*##?jA3-'>u$-##Ht;dPK7g7.u#D>#$##d1hI>#L5$;##.$RQZ$zn1i>$yH)De#~NC7@+'1:#l#D>#%##m2hku&.,#-##)ZGx~/@H'###`m*l6#@@X###5y5(c$######=2ha5%###%##,HKw#'Fl%###<7,t##RZR###"},
{912.21,934.14,10.52,4.45,"8~/######4##Jyc.,####h$#iK~r5&###'?#Jv'E>#@u#8u#5U9######G##Xyc.,####zR#NxZU,%(##2o(^~0(##wu&-d$7(<######P$#B{c.,####_.#r(Y{k####Xd$lC:.,#bP#Um%j_@######<K#:p`:5####Ci$5#GT$'rb#i7%Hh6|P$iY#;u#"},
{1052.78,124.78,13.37,3.88,"On,]J&L'7,Q#Bn)>23kR,&H$2%+0J*u5&l5$v$*6l$J,#=^4]W*^A-/(9#l#;-'@h.C:6}l&QlINK1|b#[u#BuBL5$6##je,k-<A6']>$M,$pq3[Z$Oi.-p4d{>`>$5J&@C/7lLYc#95#Ih#D>4######$##=AB.,#b$&7$(s]0Z#$+d$xd>Z-*j($@H'NM%"},
{1052.78,124.78,13.37,6.12,"m?$hi<sQ);c%5##yoS;v(P6)z]'SmS>c$=@,I-%'*2/nSOc&)~'J_/i~1D>#m5%$41I]2#Z$&7Nt/.A$)c$$=w+^z$rlSN##V>#<I'n-'D7,y,&B7'bf)Gp6r::@I&Jf1fh7]#%t1%{lSI>#m>%I5#U7(uJ/wc'O$'0S(fkB4H&7R#hmSwg0###|%#%mS.,#"},
{914.73,218.36,13.27,5.76,"feVkG$Wc$A3*ZKWZP####K%#YW@Y#$y#%?6'I6&Cw%RD=i?)ETIll'K##+0-{NWT,%###9m#@GK9[$@Q&*c#jY#l?#0|@kZ'][+OG#d%#WLW3JWOG#'##g>5FsCX>#`5$Jv%Fl%'##+A(JM4eY####z##R]2Ah<###&##)x'FWB###$##=I#UB2T,%/##7x("},
{708.09,252.92,13.41,3.00,".8.$##RG#VZ$5RO######H0'wv+######2c@>u$###'##nFEOB6?5#p>$'~)SnVB,$5,#=a,uJ1(v'<v$sqVfA1xG%I5#*BRVz;Rc$^J0;u#3sVIw.;5#{5#C</jmVIc$b,$2C,~mV1,#b>#R;>ZA0fc'x,#1nVA,$###~%#)x,1I+###B$#DQ$'~.$##=,#"},
{818.41,286.31,14.57,6.26,"<J.ZG#3KVm#%514<%,Ll$Tm(TYCg@.`P#du$jM0De/###Lu#X5$%K'IAXUG#]:5f1<|k#(##iAXNc&###@D3TA0+u####(k;###Wd#vLW.,#B4IMQ&(['nP#Q@X.,#$##{m%pm,###/##=A-######Z7>Ao3IJ1###.S('CXB@X###-##G,4tc(###Y##*I("},
{440.37,343.58,13.52,4.57,"lw0######2##0A~######%$#%}I###om)E6$######6@N&?&1U8D>####7##fA~3l$###*%##XC4@+WQ&_R$O#%SG#pA~5c$m09######j>#CD~m]5ZP#8##SyK;PKW>$-##d&2&l#1A~$##ln0######*##ZlO-?%46'xY#}A~Tc%KZ&CZ#`U;)##1A~/##"},
{800.25,389.39,13.48,1.44,"ZS$LwZ###*##a,5'::###d##c@++Z$.,#{##~P#TR$R_<4,#z0.LwZ###F#$X}Z/&2###M,#:?N{b#hP#c,%###V>#>/+4%-cm+D>#&##w]NswZ######bk03RM~#&$##m,$u$'A,$)##|-*Fl%###8##(fGkA2eY####wr-2x)H?(###;##ux,mP$###+##"},
{869.01,486.57,12.19,1.62,"1$(X>#eY#,-#hl%{h2%Q%)##8##oz.y:>###2##?AM:~/###a5%ZH#pw0I##5.(Gc9^?Q{b#vdJ'r7#vLRu#A,$ug.c@Q%##F>#y,$>/0PG#6$(AZ#}gOTq:A?QuP#)QEa4=4l$0##wCQgA3#l#(##vY#MH'95####dQ#rL;D>#$##*o&wrB/,#$##|:)f}K"},
{869.01,486.57,12.19,4.62,"A+2fYN###*##rD.p08OG#-##t?#_(=95####dG#T-)jY#$##`JRuQ)W>$q,#RwQ+48_GN@Q#d:PEV:S?(LZ#'&/A,$$##&6$WHR'##4u#LP9]aFju#y-IvE<e?R#c#kv'DH:(~.N,#eG$Dv#s$,###Z##0,B281###X5#$a4{k#L5#IH&Oh3ZP#<m#`5%6##"},
{1061.14,547.02,12.52,3.24,"ET5#########?sp######(##t-U###>Z$+l#E,$(##9C1;/0_g:#########'sp######Q##dT_###)##{G#hT-P5$.Z#*J,Zp9#########Osp######)##+0`F>#nY#R>#@I&E,$5u#=l#b/4#########4sp#########MmTqb####3##c?)95####'?#"},
{923.11,625.94,12.06,0.51,"[G4^C8######T?94f2######}?%[wI3Z%D>#)q.b~-eG$$##bZ8fv)]5$2,#K*Wjc'###Cl#US/uH(o5$oC8*:8&##.c#]^3Fx0/,#fe$_b>I&W######uD'jh?###(##Ji.I#%###(##>D7+u####@%#-(W^~1######yE,C%.pb####}$#|G%`P#D>#nP#"},
{975.33,635.87,12.87,0.59,"<u;<;8^u$3-'a7%xxH8/-#-'JS%#J.Jn(pb#2,#gY#uqR###xW-|jBeY#Al$qK/HW6~$)$J-YBR]#&TG#3[%_Z'E>#h,:#m&uqRP5$###ve&[e0;5#U5#VxOn^;###%##7=6E-)$##YI){5$qlR######j1%D&3$##.,#</&Vc&>5#D>#%[&wG%###|5&.##"},
{975.33,635.87,12.87,1.89,"###J5#;#$;^Re#%iQ'0u#o~R5v%r5&'##9E<.##OG#/##H81OG#@U*%(6O[(+%&:]RNZ&)l#759TB7###4##3R't5&&##PG#)Z$fL/c^R&.)%K1C&/Z;:f$'8]R*6'###J-#We)k,&?5#0,#2,#E(Qc83###p?']_RVJ1)##:E+|ZR###(##ez%R6).,####"},
{641.88,784.29,13.32,4.81,"I;8V,%###U,#A$(/T/ll&A$&w>$k%+{d,xb#]7);I)Wl%)Z#]58}U9yu%;#$=:,XqQ``AF#$XqQpGIUZ&DR%Rv)2Z$;v$SJ.IH@P6&]O>#B,~28A.'jdDi/OvmQtc'o-&O%B+e,:[*}P#~6'rI.a>#Be,AY69.-###C##sQ?c,%'?&q$$)L4x>#j%0BZ#B,$"},
{1030.03,862.46,13.40,3.37,"G[+######@##[vV###'##tV$8V=###=,#kD+N81###&##*|5S97######'##NUd###1##(o'#aE###t,#&c?)V;pb#$##`T%iS1#########JWd###$##_#$F0a###D##)X;2m(95#'##3;7g?)#########EWd######&##pUd######Bd$>u$######@~*"},
{260.96,1035.72,12.48,5.71,"Qu%######5##Tn[######M(#=o[#8-A,$m&#1U+~h<B,$>,#SS1######X##Yo[,c$###p'#dK[L|=)c$E$#Qi5~T3n>%4H#5J0######P$#)s[|H*###p$#W,3'aC###1##RU(%C6ql'x>${k#######:##Z~+xY$###H##Q:-@96###(##9f#L-R/,#UG#"},
{37.61,200.86,14.11,5.44,"############l-*######&##Jp/bG$Cy#.d)F,#6#$]2%pm,Z[,######F##x~ZpP$###j$#.eGC|B*Z#K##dB#OU:#N%%Q%-^8######~_%k^ZU[+4,#5K&^f*)fSqY#`P#a3*@B51##95#CZ&###$##2n<8H&.u#N?#l$I76%V6(jx$Ie.8g0pb#7##GI%"},
{707.94,254.85,14.09,3.01,"t812,#}Y$i,$sdP######fJ%/n,######IF9?u$###*##q<AEy7$##Lu#yI)7wUZP####=<*pS0v#'H5#t{U$n*,Q%|P#gyU9q:Gc$xR,K,$}{UN$*.,#(?#&+1SvU@5#3?%^B)hvURG#Cl#pq<U/1=c%y5#>wU+u####?%#`91v6+###i##ol${%1###2,#"},
{617.18,269.61,14.75,1.34,"z,'###J##uR)[jH(##:l#,2*7-(B,#<oPfR*######3+><d)D96&##G>#n#$`mPWZ&B?%wz*|,&<e&4oP1~+###B##2uJk>%F/0I>#95#b,#SnCTq:#H$3z2&-$GH;gGGyx5###n##FL45?'q(?:5####N7#{rDzb#_R'(qO/,#j$%Qb5JlP###?c#45<<R+"},
{617.18,269.61,14.75,4.77,"rbAXC<###f-$Pz4-9RQ82Iu####l)-~6RM5$.,#*##O-)O>#Zq1J3F###;##w4HhPL/6$tY8tP$[R'nwC3i=pb#$##?&.2,#,i;K#%###-##-;RX[*:c$D7&~d&-p,f7Rhc&1,#86%33C###t['vv*ZP####}&IdQ(###%##SA+IB-;YK_P#5##&C0Y1=###"},
{855.27,290.87,16.49,0.11,"^92[Z'z5%*U1T@EHw.1,#yu$_l6$x1###KZ#py.E>#G>##c#;&*V05(&0.,#_CXQ6)###5]'NEA{k####;Y2[M8######?Z#;~-uZ&OE9oG$kBX+u#'##,I$a6T###&##/K,f[,###$##J$&(c$###.|*lAX_vU###.##>rQ^'9###g##Dh/=c%###Y##.d("},
{273.59,342.42,16.54,5.05,"D>#######+Q#######>##o&.######I-#HdP######S?#8:8=R+######_?#;*E###(##wc9OQ'###Lh'u1V######Tu6jv+F95######l5#u1V###)##2%%L>H###e;,WR)######]fFOG#+v'######d>#m1V######'##~1V###s?$0u#<l$###-+/1v("},
{273.59,342.42,16.54,5.96,"T,%$##sG$qS+Ce/###$##6>4L5$###z##2AU.,####&$#olR9f3######qc$t@U{k####BE*xe-Yl&k##FBU1##ZP#3$#/@U}~0######&##szPvA4###Y5#'N(a?U%##t?)R9#WuQ&##ES/|l&.,####(##V.&Bv)######y&$;)A######P/#a?U######"},
{799.67,392.97,15.37,1.42,";B&URX###&##`=3s/5###j##S@+S#$3l$>$#lY#1y'GsCP,$sC-ZRX;5#{,%IXX7x2###-?#7#KuY#Cl$9Q%95#H,#})5ao42R*eY#G,#drU,SX.,####Ut07-II#%$##7H$8n)D>#@##<J.xY$###1d$q-E%//ZP####/q)c]-T,%###C##;y/L#%###9,#"},
{826.94,641.08,14.17,2.38,"xd#,~T######_g$j<G@5#.,#,g,8Q&@,#3m&QI)B,$###<Q#/g+K9MSu%S,#<`T(K4###&$#&6?6&2###'##Ey&=R+######ew/FZ#.@'f*.6[TT5$/H#eF1401ax2=R$/q6+W9LQ'$##2v$KQ'###J$#wxG1&136&?m$JZ8}k#Nc$ug)IsCKQ'###L,#+]."},
{888.45,743.92,14.79,3.04,"kY#;5#.S%bg97d%Oh50N1pZ(sn0G0+[R*[|/#w,$#####cs0|k#nl#_x0;d(N6&*)+;oPa?(.oPY.+UZ%=2,l083,#rY#UT+0,#`6%-~,$##Hu$Ud%BnBYV<xmPGI+UQ$m}9Ah53B4PG#-##w##%|>]P####&##;qP3Q$m>%nx${mP(##r5%`~#=mP######"},
{386.18,825.21,14.59,4.75,"tT/&?&###^##cZ&7L5|k#0?#;v&<A-XZ'j>#~965,#dP#jf*%{'m(;In+E>#Jy)4VQ:'57l#4VQ*3@zP$6m#8W:|m+76'7H#ey.:w&1uIo8+_K1(w&$W6XUQ=RQN?'il$]R@~7-Kg0*o1C,#6[*G,#m6(O,8qQ)###0###yIsP$zP$*o%%|?<5#_n*Y9,h,&"},
{242.77,958.99,17.68,5.93,"&@+######,##e(j######&$#K)j######6##@m%fY#/,####281######0##))jZP####$$#Fica@.###-##hQ#dA1.,####Ro4######1##s)j'6&###R###mCa`@######6B(K^8.,#%##?~/######(##0,jk6*###,##x9*-W@######|&'xC=######"},
{253.52,1009.75,15.13,5.84,"=6(######2##$^c######{%#M_c8m)###K$#M%&dR-######lw0######F##G^c,u####y%#T`c;i=###M##s;.{f6######L97######>##fac2v(###k##+G7(XD###&##p^(O;@######p#'######+##HN9Fl%###.##Op%^L;######Z~$BmS######"},
{710.93,163.07,17.95,3.35,".{:###a#$9-'VDSne1.,#SQ#Pe'r?S###8Z#bG#AF=>c%0,#|d)e/,O7+M5$IDSm?)qb#/##Ma>6w-###0##{c&e#&$##$c#u-*2A*rq0s(7|?SP#%%##9N*/^7eY####h'*{5&n#&###ux,95#Pc%GJ#5ASac'oy3%###]E,u#A'S###ag/###}CS###}H'"},
{1083.11,193.45,16.77,3.20,"TB7######,##ero######+$#VQS###?5#+I&###&##Mm''-'nU<######1##.so2u####'%#)IS9E<H>#Fw%g5#H2:Jd)~u%*L9######)##fto)v&PG#Y##WN;U#>Uu%.##D/(T'6/,#`,$K~0######2##rro`P#|#&Zn$5x0|6(Dy0I8+oB6zY$-##-z)"},
{748.14,225.36,17.46,2.93,"###8&&7|D95#=~.b5#iA1vH&QIT######tq'Qd*######U,A###vc$m2<OG#|U<=u#an*Vv&2KTuG%.,#4`(U^/g7/###PZC3,#qq3##K:5#d4HAV:rc'L?#.LTCZ&.,#1%#@r-nU<###2##=5#ql%?eA_B6a^936'au#-n?T2A######Kq&L7+OG####;Q$"},
{888.14,438.79,19.14,5.10,";Q%`K3^$&E?'eG#td*2g,ml';5#2,#?7>,WA###$##Uw$[^UiM)A]U8H&2u#qSGg*C.H%QQ$lR,Q,#B27p6*###1##q8.oJ2Y`U?%,<?&KD,3~Uj,%&m&LX.BR+P##lz:hZ%###*##`6()J-*{?'##TG#R7?VV>$@)UQ&i0)H>#*d&c//a,%8#$sb#wY#Bv'"},
{934.99,769.54,18.80,6.12,"&##/-%-(.biC<e(5%+^Z%*@LZZ&###]##}<Apb#)$#4L8wB3M>#CB(#wTGu$QdKoe,bf3[n%g06p#%M5$tH$###):#6vT'l#GQ&5j,LA1'.&X*:Ud(-@(lU0'T)C|@Du$H5#$##ji*5vT###k80^I)V6))I%0/,tM7#]/7A+J,#i8GX6).,####C555vT###"},
{846.47,775.37,17.86,3.51,"nl&S#$e9//|9TsDvY#b[*q[$n[-######s+6###%##]c${NAGw)7)1J3=JA.=mOH-$x/3B(&%N@'##vb#}_*###e##R/2u-*iY#Vc$X5<uB7ypOfR*$D6K],MAF/92-u#au#*##t<6vR..,#|k#U5#h:4(?&~w-~%&O^..N<$K-CpOpb#?Q$A$%ypO(c$Q5#"},
{607.83,817.61,18.24,3.91,"hR,K>#Ul$.00|9895#?5#zU)VZ'###+##1k3###$##l/.{n/>]0;y):-H*YA#APi?&nJ1iw#)29(-'-?$bD2$##jP#XN:lZ().*D^*T)6(vMI/O3W<7M1gs=S%Cl_>_#%lJ)e?'F-'m~-{6)######Kw%$QHWZ&]5$33-I@P<]1L,$M['(%B:.'8p4:~+qu%"},
{354.54,865.13,16.95,3.55,"~@+Ym'$C,LSO}6+B?$6g.Ya=*6'###2##;c@######(q$505mx,%UOQB0em(nQO)r2#J-+J$-h8OG#2##0;3/,####`-#'955l$Lv'3*.WrAf,><M:OL0:p/gfJ.d)5,#oI'j5$F6({G#fu&%##V5#x48=R+Wv(b#%:{.[SOWK4M5$'##z6A3u#Ol%X[&104"},
{528.78,279.13,21.64,1.63,"Co1$##J>#u#%HkBG$'Xx/1m%@5##6$[IJi6*dP#[P#[?%Gu${I-F5#6#$+##>STFZ%.$&#32sl'4c#IpJ*>H###K>#BGA_5%xI-|b#>u$$##mTT###;-'QH$-)=+##0WT}Q(###*##B1RmP$hl%NZ%_5%###PUTZP#1,#0,#kdNUG#EW3rH)$##&##wUT3l$"},
{702.68,808.97,20.85,5.78,"KZ%6f)3&1cH&3</#{76R*^Q%3S)FI%*m(5[&###Wm$6{6b~0}u$B1P`B7r,&g/P+11Y82^x%MtEV,$sb#Q{0O,$J>#g6$y}D*U2WM8K-%q4Cz/Pm@,Zv%oz7P6=Jx3###dl#c[$wbHfG$Gc$N</wq>D#$zG$zS%Jr>_#%R5$9z$m;B######7-$Q8OeY#oP#"},
{317.95,396.45,22.75,4.84,"eY####1##=-'^Q(######|@%:$)###T##F9K######0C&ytLD>#&##-c#~-(C5O######%n#1eY###hQ$'>3wu'###dLKg/1A$(QQ&Z5$.c#olP######_,#ZgY###5c#Al#(U7###,;NCc$v@-}k#:5#8##Q^5######'##ejYOG#$##%##L|6Xc&TN<I>#"},
{917.21,674.15,27.77,1.12,"a>#62<<-$:?'bc%p]0x'&-]3hl&.u#]1#+6R######R(#vWFi@%;CU(g2T#%W;U8}F?5#,-$:V;9?'v6##e*6,#<c$,3*696nZ&e6'_u7WbH5AU5Z%+$$;H8?]0G/1~-%KI)rl$F~-Pw+rb#~P#eG#<(2R@+QJ1###~c$%_+r,&YG#3]*pR+,H$n,%Sx-ZP#"}
};


//????,size????????????
vector<int> medianFilter(vector<int> &image, int size){
    int height = image[0];
    int width = image[1];
    vector<int> ret;
    ret.reserve(height*width+2);
    ret.push_back(height);
    ret.push_back(width);

    int bufr[2000];
    int bufg[2000];
    int bufb[2000];

    for(int y = 0; y < height; y++){
        for(int x = 0; x < width; x++){
            int cnt = 0;
            for(int i=max(y-size,0);i<=min(y+size,height-1);i++){
                for(int j=max(x-size,0);j<=min(x+size,width-1);j++){
                    int t = 2+j+i*width;
                    bufr[cnt] = (image[t]>>16)&255;
                    bufg[cnt] = (image[t]>>8)&255;
                    bufb[cnt] = image[t]&255;
                    cnt++;
                }
            }
            sort(bufr, bufr+cnt);
            sort(bufg, bufg+cnt);
            sort(bufb, bufb+cnt);
            ret.push_back((bufr[cnt/2]<<16) | (bufg[cnt/2]<<16) | bufb[cnt/2]);
        }
    }
    return ret;
}

//based on the code from: http://stackoverflow.com/questions/12380682/how-to-convert-rgb-to-hsl-in-c
void rgbToHsl(unsigned int rgb, double &h, double &s, double &l){
    double r = (rgb >> 16) & 255;
    double g = (rgb >> 8) & 255;
    double b = rgb & 255;
    r /= 255;
    g /= 255;
    b /= 255;
    double maxv = max(max(r, g), b);
    double minv = min(min(r, g), b);
    double d = maxv - minv;
    h = 0; s = 0; l = (maxv + minv) / 2;
    if (maxv != minv) {
        s = l > 0.5 ? d / (2 - maxv - minv) : d / (maxv + minv);
        if (maxv == r) { h = (g - b) / d + (g < b ? 6 : 0); }
        else if (maxv == g) { h = (b - r) / d + 2; }
        else if (maxv == b) { h = (r - g) / d + 4; }
        h /= 6;
    }
}

//????,s?????
vector<int> scale(vector<int> &image, int s){
    int height = image[0];
    int width = image[1];
    vector<int> ret;
    ret.reserve(height/s*width/s+2);
    ret.push_back(height/s);
    ret.push_back(width/s);

    for(int y = 0; y < height/s; y++){
        for(int x = 0; x < width/s; x++){
            int r = 0, g = 0, b = 0;
            double sum = 0;
            for(int i=0;i<s&&i+y*s<height;i++){
                for(int j=0;j<s&&j+x*s<width;j++){
                    int t = 2+(x*s+j)+(y*s+i)*width;
                    int tr = (image[t]>>16)&255;
                    int tg = (image[t]>>8)&255;
                    int tb = image[t]&255;
                    r+=tr;
                    g+=tg;
                    b+=tb;
                    sum++;
                }
            }
            r = (int)(r/sum+0.5);
            g = (int)(g/sum+0.5);
            b = (int)(b/sum+0.5);
            ret.push_back((r<<16)|(g<<8)|b);
        }
    }

    return ret;
}

struct siftPoint{
    double x, y;
    double sigma, angle;
    double vec[128];
};

vector<siftPoint> getSiftVec(VL::PgmBuffer &buffer){
    int    first          = 1 ;//-1
    int    octaves        = 3 ;//-1
    int    levels         = 4 ;//3
    //int    first          = -1 ;
    //int    octaves        = -1 ;
    //int    levels         = 3 ;
    float  threshold      = 0.04f / levels / 2.0f ;
    float  edgeThreshold  = 10.0f;
    float  magnif         = 3.0 ;
    int    nodescr        = 0 ;
    int    noorient       = 0 ;
    int    stableorder    = 0 ;
    int    savegss        = 0 ;
    int    verbose        = 0 ;
    int    binary         = 0 ;
    int    unnormalized   = 0 ;
    int    fp             = 0 ;

    vector<siftPoint> ret;

    // ---------------------------------------------------------------
    //                                            Gaussian scale space
    // ---------------------------------------------------------------    

    int         O      = octaves ;    
    int const   S      = levels ;
    int const   omin   = first ;
    float const sigman = .5 ;
    float const sigma0 = 1.6 * powf(2.0f, 1.0f / S) ;

    // optionally autoselect the number number of octaves
    // we downsample up to 8x8 patches
    if(O < 1) {
        O = max((int)(floor(log2(std::min(buffer.width,buffer.height))) - omin -3), 1) ;
    }

    // initialize scalespace
    VL::Sift sift(buffer.data, buffer.width, buffer.height, 
        sigman, sigma0,    O, S, omin, -1, S+1) ;


    // -------------------------------------------------------------
    //                                             Run SIFT detector
    // -------------------------------------------------------------    

    sift.detectKeypoints(threshold, edgeThreshold) ;

    // -------------------------------------------------------------
    //                  Run SIFT orientation detector and descriptor
    // -------------------------------------------------------------    

    /* set descriptor options */
    sift.setNormalizeDescriptor( ! unnormalized ) ;
    sift.setMagnification( magnif ) ;

    // -------------------------------------------------------------
    //            Run detector, compute orientations and descriptors
    // -------------------------------------------------------------
    for(VL::Sift::KeypointsConstIter iter = sift.keypointsBegin(); iter != sift.keypointsEnd(); ++iter) {

        // detect orientations
        VL::float_t angles [4] ;
        int nangles = sift.computeKeypointOrientations(angles, *iter) ;

        // compute descriptors
        for(int a = 0 ; a < nangles ; ++a) {
            siftPoint pt;
            pt.x = iter->x;
            pt.y = iter->y;
            pt.sigma = iter->sigma;
            pt.angle = angles[a];

            sift.computeKeypointDescriptor(pt.vec, *iter, angles[a]) ;
            ret.push_back(pt);
        } // next angle
    } // next keypoint

    return ret;
}

vector<siftPoint> getSiftVec(const char *fname){
    ifstream in(fname, ios::binary) ; 
    VL::PgmBuffer buffer ;
    extractPgm(in, buffer) ;
    return getSiftVec(buffer);
}

vector<siftPoint> getSiftVec(vector<int> &img){
    VL::PgmBuffer buffer ;
    extractVec(img, buffer) ;
    return getSiftVec(buffer);
}

//p1,p2???????
//q1,q2??????
//p1,q1,p2??, p1,q1??n?? p2?q2?m??,??????
//n???4,???????????
void locate_points(int n, double *p1, double *q1, int m, double *p2, double *q2) {
    //always n = 4 in this code
    double A[1000];
    double b[1000];
    double x[1000];
    for (int i = 0; i < n; i++) {
        A[i*2*8] = p1[i*2];
        A[i*2*8 + 1] = p1[i*2 + 1];
        A[i*2*8 + 2] = 1;
        A[i*2*8 + 3] = 0;
        A[i*2*8 + 4] = 0;
        A[i*2*8 + 5] = 0;
        A[i*2*8 + 6] = -q1[i*2] * p1[i*2];
        A[i*2*8 + 7] = -q1[i*2] * p1[i*2+1];
        A[i*2*8 + 8] = 0;
        A[i*2*8 + 9] = 0;
        A[i*2*8 + 10] = 0;
        A[i*2*8 + 11] = p1[i*2];
        A[i*2*8 + 12] = p1[i*2 + 1];
        A[i*2*8 + 13] = 1;
        A[i*2*8 + 14] = -q1[i*2 + 1] * p1[i*2];
        A[i*2*8 + 15] = -q1[i*2 + 1] * p1[i*2+1];
        b[i*2] = q1[i*2];
        b[i*2 + 1] = q1[i*2 + 1];
    }
    for (int i = 0; i < 2*n*8; i++) A[i]/=1e3;
    for (int i = 0; i < 2*n; i++) b[i]/=1e3;

    double u[1000], v[1000], w[1000];
    lsqr(2*n, 8, aprod, 0, A, b, v, w, x, NULL, 1e-18, 1e-18, 1e50, 2000, NULL);

    for (int i = 0; i < m; i++) {
        //p2[i*2]
        //p2[i*2+1]
        q2[i*2] = (x[0] * p2[i*2] + x[1] * p2[i*2+1] + x[2] * 1) / (x[6] * p2[i*2] + x[7] * p2[i*2 + 1] + 1);
        q2[i*2+1] = (x[3] * p2[i*2] + x[4] * p2[i*2+1] + x[5] * 1) / (x[6] * p2[i*2] + x[7] * p2[i*2 + 1] + 1);
    }


}

//typedef vector<pair<int, int> > matchresult
bool is_sim = false;
bool is_iss = false;

bool check_led(int height, int width, int x, int y, int r, int *red, int *green, int *blue) {
    int counter = 0;
    if (is_iss && r < 16) r = 10;
    
    for (int i = x - r; i <= x + r; i++){
        if (i < 0 || i >= width) continue;
        for (int j = y - r; j <= y + r; j++){
            if (j < 0 || j >= height) continue;
            if (is_iss) {
                if ((green[j * width + i] == 255 || blue[j * width + i] == 255) && red[j*width + i] < 180) counter++;
            }
            else {
                if (green[j * width + i] > (is_sim ? 25: 70) + red[j * width + i] || blue[j*width + i] > (is_sim ? 25: 70) + red[j * width + i]) {
                    counter++;
                }
            }
        }
    }
    
    if (is_iss) {
        if (r > 16) return counter > 60;
        else return counter > 30;
    }
    
    if (r > 16 && counter > (is_sim ? 10 : 50)) return true;
    if (r < 16 && counter > (is_sim ? 10 : 20)) return true;
    else return false;
}

struct point{
    double x, y;
};
//iss
double bpt1[] = {-1,-1,553,351,660,261,695,596,695,493,-1,-1,261,417,348,434,436,450,282,487,364,506,457,524,303,559,387,575,471,594,366,730,352,668,556,769,547,709,572,822,761,816,755,750};

//lab1-1
double bpt2[] = {597,252,496,223,661,256,699,588,701,491,-1,-1,265,419,351,435,442,448,287,489,374,504,459,520,308,559,395,574,478,590,374,725,360,667,568,764,558,706,585,820,773,812,767,750};

//lab1-3
double bpt3[] = {727,250,836,227,805,233,737,591,743,481,-1,-1,221,562,336,545,452,525,214,638,330,622,440,607,212,717,322,698,431,680,223,888,223,822,460,843,469,783,446,897,688,802,700,738};

//lab2-1
double bpt4[] = {488,154,560,135,538,147,489,383,494,311,-1,-1,151,350,227,343,302,332,146,404,220,397,293,388,142,455,213,448,284,438,142,575,143,532,299,553,304,510,288,590,449,532,460,489};

//lab2-2
double bpt5[] = {475,160,539,137,524,152,478,385,486,316,-436,-290,145,359,219,351,294,340,139,411,214,405,285,394,132,465,208,453,277,446,139,586,139,541,294,559,299,519,285,598,444,536,454,492};

//lab3-1
double bpt6[] = {374,128,392,21,419,132,415,370,414,310,409,428,118,295,181,299,245,300,115,358,179,361,242,362,114,422,178,425,241,426,120,590,122,532,271,596,271,537,267,656,420,599,420,542};

//lab3-2
double bpt7[] = {434,150,466,55,484,151,476,429,472,357,464,490,130,353,203,354,278,354,130,425,203,425,274,425,126,496,199,495,272,496,136,674,137,614,303,673,306,610,299,735,467,670,471,612};

//lab3-3
double bpt8[] = {313,169,368,57,339,171,383,397,361,344,375,461,146,334,197,336,244,336,155,404,206,404,252,403,162,475,213,469,261,469,192,658,183,593,307,645,298,580,313,705,409,628,400,571};

//lab1-2
double bpt9[] = {-1,-1,506,343,617,246,711,579,703,481,-1,-1,310,411,389,429,469,438,345,481,422,497,502,514,380,547,454,563,530,579,469,715,450,656,650,756,627,697,676,810,835,802,816,740};

//lab1-4
double bpt10[] = {-1,-1,807,367,875,238,877,591,871,487,-1,-1,314,563,439,550,558,532,330,643,449,624,566,610,342,719,462,700,577,682,400,894,388,828,653,849,645,785,-1,-1,893,809,893,746};
class RobonautEye{
private:
    void rec(vector<int> &_eyeImage, char *filename){
        vector<int> eyeImage = scale(_eyeImage, 2);
    }

    double sqr(double a){
        return a*a;
    }

    double dist(siftPoint &a, siftPoint &b){
        double sum = 0;
        for(int i = 0; i < 128;i++){
            sum += sqr(a.vec[i] - b.vec[i]);
        }
        return sum;
    }

    vector<pair<int, int> > mymatch(const char *filename, const char *filename2, vector<siftPoint> &pattern, vector<siftPoint> &pic){
        vector<int> match, backMatch;
        match.reserve(pattern.size());

        vector<double> score;
        backMatch.reserve(pic.size());
        for(int i = 0; i < pic.size(); i++)
            backMatch.push_back(-1);
        score.reserve(pic.size());
        for(int i = 0; i < pic.size(); i++)
            score.push_back(1e100);

        for(int i=0;i<pattern.size();i++){
            double mini = 1e100;
            double m2 = 1e100;
            double bx, by;
            int bi;
            for(int j=0;j<pic.size();j++){
                double d = dist(pattern[i], pic[j]);
                if(d < mini){
                    m2 = mini;
                    mini = d;
                    bx = pic[j].x;
                    by = pic[j].y;
                    bi = j;
                }else if(d < m2){
                    m2 = d;
                }
            }
            if(m2 > mini*1.5){
                if(backMatch[bi] != -1){ //????????????
                    if(score[bi] < mini){ //???????
                        match.push_back(-1);
                    }else{ //??????,??
                        score[bi] = mini;
                        match[backMatch[bi]] = -1;
                        backMatch[bi] = i;
                        match.push_back(bi);
                    }
                }else{
                    score[bi] = mini;
                    backMatch[bi] = i;
                    match.push_back(bi);
                }
            }else{
                match.push_back(-1);
            }
        }

//        FILE *fout = fopen(filename, "w");
//        FILE *fout2 = fopen(filename2, "w");
        vector<pair<int, int> > ret;
        for(int i = 0; i < match.size(); i++){
            if(match[i] != -1){
                ret.push_back(make_pair(i, match[i]));
//                fprintf(fout, "%f %f %d\n", pic[match[i]].x*2, pic[match[i]].y*2, i);
//                fprintf(fout2, "%f %f %d\n", pattern[i].x, pattern[i].y, i);
            }
        }
//        fclose(fout);
//        fclose(fout2);
        
        return ret;
    }

    pair<vector<point>,int> randmatch(vector<siftPoint> &pattern, vector<siftPoint> &pic, vector<pair<int, int> > &pairs, double *answer){

        double all_points[] = {
        349,    116,
        349,    116,
        394,    116,
        393,    355,
        394,    295,
        394,    414,
        97,    292,
        160,    293,
        223,    291,
        96,    355,
        161,    356,
        223,    355,
        96,    418,
        160,    418,
        224,    418,
        109,    584,
        108,    526,
        257,    585,
        259,    526,
        260,    643,
        408,    586,
        408,    526
        };


        vector<bool> used;
        used.reserve(pairs.size());
        for(int i=0;i<pairs.size();i++){
            used.push_back(false);
        }
        double p1[8], q1[8], bestp1[8], bestq1[8];
        double *p2 = new double[pairs.size()*2];
        double *q2 = new double[pairs.size()*2];
        double *q2ans = new double[pairs.size()*2];

        double anspic[44];

        for(int i=0;i<pairs.size();i++){
            int tp = pairs[i].first;
            int tq = pairs[i].second;
            p2[i*2] = pattern[tp].x;
            p2[i*2+1] = pattern[tp].y;
            q2ans[i*2] = pic[tq].x;
            q2ans[i*2+1] = pic[tq].y;
        }

        int best = 0;
        for(int iter = 0; iter<1000;iter++){
            int index[4];
            for(int i=0;i<4;i++){
                int t = 0;
                while(true){
                    t = rand()%pairs.size();
                    if(used[t])
                        continue;
                    index[i] = t;
                    used[t] = true;
                    break;
                }
                int tp = pairs[t].first;
                int tq = pairs[t].second;
                p1[i*2] = pattern[tp].x;
                p1[i*2+1] = pattern[tp].y;
                q1[i*2] = pic[tq].x;
                q1[i*2+1] = pic[tq].y;
            }
            for(int i=0;i<4;i++){
                used[index[i]] = false;
            }

            locate_points(4, p1, q1, pairs.size(), p2, q2);

            int ok = 0;
            for(int i = 0; i < pairs.size(); i++){
                if(sqr(q2ans[i*2]-q2[i*2])+sqr(q2ans[i*2+1]-q2[i*2+1]) < 36){
                    ok++;
                }
            }
            if(ok > best){
                best = ok;
                memcpy(bestp1, p1, sizeof(p1));
                memcpy(bestq1, q1, sizeof(q1));
            }

            if(best > 100) //??????
                break;
        }
        
        /*
    locate_points(4, bestp1, bestq1, 22, p2, q2);
    int idx = 0;
    for(int i = 0; i < pairs.size(); i++){
        if(sqr(q2ans[i*2]-q2[i*2])+sqr(q2ans[i*2+1]-q2[i*2+1]) < 36){
            p2[idx*2] = p2[i*2];
            p2[idx*2+1] = p2[i*2+1];
            q2[idx*2] = q2ans[i*2];
            q2[idx*2+1] = q2ans[i*2+1];
            idx++;
        }
    }
    locate_points(idx, p2, q2, 22, answer, anspic);
    delete[] p2;
    delete[] q2;
    delete[] q2ans;

        */
        
        
        delete[] p2;
        delete[] q2;
        delete[] q2ans;

        //fprintf(stderr, "best: %d\n", best);
        //fflush(stderr);
        locate_points(4, bestp1, bestq1, 22, answer, anspic);
        
        int adjust_n = 0;
        double adjust_p1[44];
        double adjust_q1[44];
        int adjust_m = 0;
        double adjust_p2[44];
        double adjust_q2[44];
        for (int i = 0; i < 22; i++){
            if (answer[i*2] <= -1){
                adjust_p2[adjust_m*2] = all_points[i*2];
                adjust_p2[adjust_m*2+1] = all_points[i*2+1];
                adjust_m++;
                continue;
            }
            adjust_p1[adjust_n*2] = all_points[i*2];
            adjust_p1[adjust_n*2+1] = all_points[i*2+1];
            adjust_q1[adjust_n*2] = anspic[i*2];
            adjust_q1[adjust_n*2+1] = anspic[i*2+1];
            adjust_n++;
        }
        if (adjust_n < 22) {
            locate_points(adjust_n, adjust_p1, adjust_q1, adjust_m, adjust_p2, adjust_q2);
            adjust_m = 0;
            for (int i = 0; i < 22; i++){
                if (answer[i*2] <= -1){
                    anspic[i*2] = adjust_q2[adjust_m*2];
                    anspic[i*2+1] = adjust_q2[adjust_m*2+1];
                    adjust_m++;
                }
            }
        }
        
        
        vector<point> ret;
        for(int i =0;i<22;i++){
            point p = {anspic[i*2]*scaled, anspic[i*2+1]*scaled}; //??????
            ret.push_back(p);
        }
        return make_pair(ret, best);
    }
    
    bool determine_line(vector<int> &x, vector<int> &y, double &param1, double &param2, double &param3) {
        fprintf(stderr, "x.size() = %d\n", x.size());
        if (x.size() < 80) return false;
        //??????????????
        /*
        for (int i = 0; i < 50; i++) {
            fprintf(stderr, "%d %d\n", x[i], y[i]);
        }*/
        //fprintf(stderr, "---\n");
        fflush(stderr);
        
        double best_err = 1e100;
        double best_x, best_y;
        for (int ca = 0; ca < 1000; ca++) {
            double r = rand() % x.size(); // % (x.size() / 3);
            double a = x[r];
            double b = y[r];
            double c = -1;
            r = rand() % x.size(); //% (x.size() / 3) + x.size() / 3;
            double d = x[r];
            double e = y[r];
            double f = -1;
            
            if (fabs(a*e-b*d) < 1e-13) continue;
            double cx = (c*e-f*b) / (a*e-b*d);
            double cy = (-c*d+f*a) / (a*e-b*d);
            
            double sum_err = 0;
            for (int i = 0; i < x.size(); i++) {
                //fprintf(stderr, "%f\n", sol[0] * x[i] + sol[1] * y[i] + 1);
                sum_err += sqrt(fabs(cx * x[i] + cy * y[i] + 1));
            }
            
            if (sum_err < best_err) {
                best_err = sum_err;
                best_x = cx;
                best_y = cy;
            }
        }
        
        param1 = best_x;
        param2 = best_y;
        param3 = 1;
                
        return true;
    }
    
    void cross(double a, double b, double c, double d, double e, double f, double &cx, double &cy){
            
            cx = -(c*e-f*b) / (a*e-b*d);
            cy = -(-c*d+f*a) / (a*e-b*d);
    }
    
    vector<int> rec_sim(int height, int width, int *red, int *green, int*blue) {
        
        double p1[] = {0, 0, 518, 0, 0, 748, 518, 748};
        double q1[] = {747,949,1801,839, 779,1881, 1509,1741};        
        double all_points[] = {
        349,    116,
        349,    116,
        394,    116,
        393,    355,
        394,    295,
        394,    414,
        97,    292,
        160,    293,
        223,    291,
        96,    355,
        161,    356,
        223,    355,
        96,    418,
        160,    418,
        224,    418,
        109,    584,
        108,    526,
        257,    585,
        259,    526,
        260,    643,
        408,    586,
        408,    526
        };

        //?????????
        bool *vis = new bool[width * height];
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                vis[i*width + j] = false;
            }
        }
        
        int dx[] = {-1,1,0,0};
        int dy[] = {0,0,-1,1};
        
        int best_size = 0;
        int best_point;
        
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                if (vis[i * width + j]) continue;
                if (red[i * width + j] > green[i * width + j] && red[i * width + j] > blue[i * width + j]) {
                    queue<int> q;
                    q.push(i*width +j);
                    vis[i*width +j] = true;
                    int size = 0;
                    while (!q.empty()) {
                        int x = q.front() % width;
                        int y = q.front() / width;
                        q.pop();
                        for (int k = 0; k < 4; k++) {
                            int nx = x + dx[k];
                            int ny = y + dy[k];
                            if (nx < 0 || nx >= 1600) continue;
                            if (ny < 0 || ny >= 1200) continue;
                            if (vis[ny*width+nx]) continue;
                            if (red[ny*width+nx] > green[ny*width+nx] && red[ny*width+nx] > blue[ny*width+nx]) {
                                q.push(ny*width+nx);
                                vis[ny*width+nx] = true;
                                size++;
                            }
                        }
                    }
                    if (size > best_size) {
                        best_size = size;
                        best_point = i*width +j;
                    }
                    
                }
            }
        }
        
        queue<int> q;
        int x = best_point % width;
        int y = best_point / width;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                vis[i*width + j] = false;
            }
        }
        
        fprintf(stderr, "red %d %d\n", x, y);
        fflush(stderr);
        for (int i = 0; i < 200; i++) {
            if (x - 50 >= 0) {
                if (red[(x - 50) + width * (y+i)] == 0) break;
                q.push((x - 50) + width * (y+i));
                vis[(x - 50) + width * (y+i)] = true;
            }
        }
        
        for (int i = 0; i < 200; i++) {
            if (x + 50 < 1600) {
                if (red[(x + 50) + width * (y+i)] == 0) break;
                q.push((x + 50) + width * (y+i));
                vis[(x + 50) + width * (y+i)] = true;
            }
        }
        
        

        int left_boundary[1200];
        int right_boundary[1200];
        int top_boundary[1600];
        int bottom_boundary[1600];
        
        
        
        for (int i = 0; i < 1200; i++) {
            left_boundary[i] = -1;
            right_boundary[i] = -1;
        }
        for (int i = 0; i < 1600; i++) {
            top_boundary[i] = -1;
            bottom_boundary[i] = -1;
        }
        
        //int sim_q1[1000];
        //int m = sizeof(p2)/sizeof(double)/2;
        double sim_q1[100];

            double corner1 = 1e100;
            double corner2 = -1e100;
            double corner3 = 1e100;
            double corner4 = -1e100;
        while (!q.empty()) {
            int x = q.front() % width;
            int y = q.front() / width;
            q.pop();
            for (int j = 0; j < 4; j++) {
                int nx = x + dx[j];
                int ny = y + dy[j];
                if (nx < 0 || nx >= 1600) continue;
                if (ny < 0 || ny >= 1200) continue;
                if (vis[ny*width+nx]) continue;
                if (red[ny*width+nx] == 0) continue;
                
                if (abs(red[ny*width+nx] - red[y* width +x]) < 3 && abs(green[ny*width+nx] - green[y*width + x]) < 3 && abs(blue[ny*width+nx] - blue[y * width + x]) < 3 ) {
                    q.push(ny*width+nx);
                    vis[ny*width+nx] = true;
                    if (top_boundary[nx] == -1 || top_boundary[nx] > ny) {
                        top_boundary[nx] = ny;
                    }
                    if (bottom_boundary[nx] == -1 || bottom_boundary[nx] < ny) {
                        bottom_boundary[nx] = ny;
                    }
                    if (left_boundary[ny] == -1 || left_boundary[ny] > nx) {
                        left_boundary[ny] = nx;
                    }
                    if (right_boundary[ny] == -1 || right_boundary[ny] < nx) {
                        right_boundary[ny] = nx;
                    }
                    if (0.5*nx+ny < corner1) {
                        corner1 = 0.5*nx+ny;
                        sim_q1[0] = nx;
                        sim_q1[1] = ny;
                    }
                    if (0.5*nx-ny > corner2) {
                        corner2 = 0.5*nx-ny;
                        sim_q1[2] = nx;
                        sim_q1[3] = ny;
                    }
                    if (0.5*nx-ny < corner3) {
                        corner3 = 0.5*nx-ny;
                        sim_q1[4] = nx;
                        sim_q1[5] = ny;
                    }
                    if (0.5*nx+ny > corner4) {
                        corner4 = 0.5*nx+ny;
                        sim_q1[6] = nx;
                        sim_q1[7] = ny;
                    }
                }
            }
        }
        
        vector<int> t_x, b_x, l_x, r_x;
        vector<int> t_y, b_y, l_y, r_y;
        
        fprintf(stderr, "sim_q1\n");
        for (int i = 0; i < 4; i++) {
            fprintf(stderr, "%f %f\n", sim_q1[i*2],sim_q1[i*2+1]);
        }
        fflush(stderr);
        for (int i = 0; i < 1600; i++) {
            if (top_boundary[i] != -1 && top_boundary[i] != 0 && i > sim_q1[0] && i < sim_q1[2]) {
                t_x.push_back(i);
                t_y.push_back(top_boundary[i]);
            }
            if (bottom_boundary[i] != -1 && bottom_boundary[i] != 1199 && i > sim_q1[4] && i < sim_q1[6]) {
                b_x.push_back(i);
                b_y.push_back(bottom_boundary[i]);
            }
        }
        
        fprintf(stderr, "sim_q1[1] = %f sim_q1[5] = %f\n", sim_q1[1], sim_q1[5]);
        fflush(stderr);
        
        for (int i = 0; i < 1200; i++) {
            if (left_boundary[i] != -1 && left_boundary[i] != 0 && i > sim_q1[1] && i < sim_q1[5]) {
                l_x.push_back(left_boundary[i]);
                l_y.push_back(i);
            }
            if (right_boundary[i] != -1 && right_boundary[i] != 1599 && i > sim_q1[3] && i < sim_q1[7]) {
                r_x.push_back(right_boundary[i]);
                r_y.push_back(i);
            }
        }
    
        fprintf(stderr, "---\n");
        for (int i = 0; i < b_x.size(); i+=20) {
            fprintf(stderr, "%d %d\n", b_x[i], b_y[i]);
        }
        fprintf(stderr, "---\n");
        
        //??????
        bool has_top, has_bottom, has_left, has_right;
        double t_a, t_b, t_c;
        double b_a, b_b, b_c;
        double l_a, l_b, l_c;
        double r_a, r_b, r_c;
        has_top = determine_line(t_x, t_y, t_a, t_b, t_c);
        has_bottom = determine_line(b_x, b_y, b_a, b_b, b_c);
        has_left = determine_line(l_x, l_y, l_a, l_b, l_c);
        has_right = determine_line(r_x, r_y, r_a, r_b, r_c);

        int n = 0;
        //double sim_q1[1000];
        double sim_p1[1000];
        
        if (has_top && has_left) {
            sim_p1[n] = p1[0];
            sim_p1[n + 1] = p1[1];
            cross(t_a, t_b, t_c, l_a, l_b, l_c, sim_q1[n], sim_q1[n + 1]);
            n += 2;
        }
        if (has_top && has_right) {
            sim_p1[n] = p1[2];
            sim_p1[n + 1] = p1[3];
            cross(t_a, t_b, t_c, r_a, r_b, r_c, sim_q1[n], sim_q1[n + 1]);
            n += 2;
        }
        if (has_bottom && has_left) {
            sim_p1[n] = p1[4];
            sim_p1[n + 1] = p1[5];
            cross(b_a, b_b, b_c, l_a, l_b, l_c, sim_q1[n], sim_q1[n + 1]);
            n += 2;    
        }
        if (has_bottom && has_right) {
            sim_p1[n] = p1[6];
            sim_p1[n + 1] = p1[7];
            cross(b_a, b_b, b_c, r_a, r_b, r_c, sim_q1[n], sim_q1[n + 1]);
            n += 2;
        }
        
        for (int i = 0; i < 4; i++) {
            fprintf(stderr, "%f %f\n", sim_q1[i*2],sim_q1[i*2+1]);
        }
        fprintf(stderr, "n=%d\n", n);
        fflush(stderr);
        /*
        if (n != 8) {
            exit(0);
        }
        */
        double q2[1000];
        locate_points(4, sim_p1, sim_q1, 22, all_points, q2);
                
        vector<int> ret;
        for(int i = 0; i < 22; i++){
            //printf("%f\t%f\n", q2[i*2],q2[i*2+1]);
            ret.push_back((int)(q2[i*2] + 0.5));
            ret.push_back((int)(q2[i*2+1] + 0.5));
        }
        
        delete vis;
        return ret;
    }

    vector<int> getPatternOrder(vector<vector<pair<int, int> > > &ms){
        vector<int> order;
        order.reserve(ms.size());
        for(int i = 0; i < ms.size(); i++){
            order.push_back(i);
        }
        for(int i = 0; i < order.size(); i++){
            for(int j = 0; j < order.size()-1; j++){
                if(ms[order[j]].size() < ms[order[j+1]].size()){
                    swap(order[j], order[j+1]);
                }
            }
        }
        return order;
    }

    vector<siftPoint> convertSiftPoint(siftPoint_save *sps, int len){
        vector<siftPoint> ret;
        ret.reserve(len);
        for(int i = 0; i < len; i++){
            siftPoint sp;
            sp.x = sps[i].x;
            sp.y = sps[i].y;
            sp.sigma = sps[i].sigma;
            sp.angle = sps[i].angle;
            for(int j = 0; j < 192; j++){
                if(sps[i].vec[j] == '~')
                    sps[i].vec[j] = '\\';
                sps[i].vec[j] -= 35;
            }
            for(int j = 0,k=0; j < 192; j+=3,k+=2){
                int num = ((sps[i].vec[j+2]*91)+sps[i].vec[j+1])*91+sps[i].vec[j];
                double v1 = num % 830;
                double v0 = num / 830;
                sp.vec[k] = v0/2000;
                sp.vec[k+1] = v1/2000;
            }
            ret.push_back(sp);
        }
        return ret;
    }
    
    inline bool is_red(int x, int y, int width, int height, int* r, int *g, int *b) {
        if (x < 0 || x >= width || y < 0 || y >= height) return false;
        int red = r[y*width+x];
        int green = g[y*width+x];
        int blue = b[y*width+x];
        if (red <= 40 && red > green + 20 && (is_iss || red > blue + 20)) return true;
        if (red <= 60 && red > green + 40 && (is_iss || red > blue + 40)) return true;
        if (red <= 90 && red > green + 60 && (is_iss || red > blue + 60)) return true;
        if (red > 120 && green < 70 && (is_iss || blue < 70)) return true;
        if (red > green + 110 && (is_iss || red > blue + 110)) return true;

        return false;
    }
    
    bool local_red_cover_search(int width, int height, int *red, int *green, int *blue, int sx, int sy, int &x, int &y) {
        //?????????true
        int y_dist = 130;
        int x_dist = 250;    
        fprintf(stderr, "local_red_cover_search %d %d\n", sx, sy);
        fflush(stderr);
        
        for (int j = sy + y_dist; j >= sy + 40; j--) {
            if (j < 0 || j >= height) continue;
            for (int dx = x_dist; dx >= 0 ; dx--) {
                int i = sx + dx;
                bool ok = true;
                for (int k = j; k >= j - 5; k--) {
                    if (!is_red(i, k, width, height, red, green, blue)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    x = i - 10;
                    y = j;
                    return false;
                }
                ok = true;
                i = sx - dx;
                for (int k = j; k >= j - 5; k--) {
                    if (!is_red(i, k, width, height, red, green, blue)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    x = i + 10;
                    y = j;
                    return false;
                }
            }
        }
        for (int dx = x_dist; dx >= 0 ; dx--) {
            for (int j = sy - y_dist; j < sy + 40; j++) {
                if (j < 0 || j >= height) continue;
            
                int i = sx + dx;
                bool ok = true;
                for (int k = j; k >= j - 5; k--) {
                    if (!is_red(i, k, width, height, red, green, blue)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    x = i - 10;
                    y = j;
                    return true;
                }
                i = sx - dx;
                ok = true;
                for (int k = j; k >= j - 5; k--) {
                    if (!is_red(i, k, width, height, red, green, blue)) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    x = i + 10;
                    y = j;
                    return true;
                }
            }
        }
        return true;
    }
    
    void local_black_button_search(int width, int height, int*red, int *green, int *blue, int &x, int &y) {
        bool *vis = new bool[width*height];
        memset(vis, 0, sizeof(bool) * width * height);
        
        int min_dist = 1000000;
        int best_x = x;
        int best_y = y;
        for (int i = x - 200; i <= x + 200; i++) {
            if (i < 0 || i >= width) continue;
            for (int j = y - 200; j <= y + 200; j++) {
                if (j < 0 || j >= height) continue;
                if (red[j * width + i] + green[j * width + i] + blue[j * width + i] == 0) {
                    if (abs(j - y) + abs(i - x) < min_dist) {
                        min_dist = abs(j - y) + abs(i - x);
                        best_x = i;
                        best_y = j;
                    }
                }
            }
        }
        queue<int> q;
        q.push(best_y * width + best_x);
        vis[best_y * width + best_x] = true;
        int dx[] = {-1, 1, 0, 0};
        int dy[] = {0, 0, -1, 1};
        int sum_x = 0;
        int sum_y = 0;
        int c = 0;
        while (!q.empty()) {
            int i = q.front() % width;
            int j = q.front() / width;
            q.pop();
            sum_x += i;
            sum_y += j;
            c++;
            for (int k = 0; k < 4; k++) {
                if (i + dx[k] < 0) continue;
                if (i + dx[k] >= width) continue;
                if (j + dy[k] < 0) continue;
                if (j + dy[k] >= height) continue;
                int nx = i + dx[k];
                int ny = j + dy[k];
                if (vis[ny * width + nx]) continue;
                vis[ny * width + nx] = true;
                if (red[ny * width + nx] + green[ny * width + nx] + blue[ny * width + nx] == 0) {
                    q.push(ny * width + nx);
                    
                }
            }
        }
        x = sum_x / c;
        y = sum_y / c;
        
        delete vis;
    }

    int scaled;
public:
    vector<string> recognizeObjects(vector<int> leftEyeImage, vector<int> rightEyeImage){
        int height = leftEyeImage[0];
        int width = leftEyeImage[1];
        vector<int> lretint;
        vector<int> rretint;
        fprintf(stderr, "begin %d\n", clock());
        int *lred = new int[height * width];
        int *lgreen = new int[height * width];
        int *lblue = new int[height * width];
        int *rred = new int[height * width];
        int *rgreen = new int[height * width];
        int *rblue = new int[height * width];
                
        for(int y = 0, t=2; y < height; y++){
            for(int x = 0; x < width; x++,t++){
                lred[x + y * width] = ((leftEyeImage[t] >> 16) & 255);
                lgreen[x + y * width] = ((leftEyeImage[t] >> 8) & 255);
                lblue[x + y * width] = (leftEyeImage[t] & 255);
            }
        }

        for(int y = 0, t=2; y < height; y++){
            for(int x = 0; x < width; x++,t++){
                rred[x + y * width] = ((rightEyeImage[t] >> 16) & 255);
                rgreen[x + y * width] = ((rightEyeImage[t] >> 8) & 255);
                rblue[x + y * width] = (rightEyeImage[t] & 255);
            }
        }
        
        int sim_counter = 0;
        int iss_counter = 0;
        for(int y = 0, t=2; y < height; y++){
            for(int x = 0; x < width; x++,t++){
                if (rred[x + y*width] == rgreen[x+y*width] && rred[x+y*width] == rblue[x+y*width]){
                    sim_counter++;
                }
                if (rblue[x+y*width] > rred[x+y*width] + 10 && rblue[x+y*width] > rgreen[x+y*width] + 10) {
                    iss_counter++;
                }
            }
        }
        if (sim_counter >= height * width / 2) {
            is_sim = true;
        }
        else if (iss_counter >= height * width / 2) {
            is_iss = true;
        }
        fprintf(stderr, "is_iss=%d\n", is_iss);
        
        if (is_sim) {
            fprintf(stderr, "is_sim\n");
            rretint = rec_sim(height, width, rred, rgreen, rblue);
            lretint = rec_sim(height, width, lred, lgreen, lblue);
        }
        else {
            fprintf(stderr, "#begin2 %d\n", clock());
            char *allPatternFiles[] = {"iss.pgm", "lab1-1.pgm", "lab1-3.pgm", "lab2-1.pgm", "lab2-2.pgm", "lab3-1.pgm", "lab3-2.pgm", "lab3-3.pgm", "lab1-2.pgm", "lab1-4.pgm"};
            siftPoint_save *spss[] = {sps1, sps2, sps3, sps4, sps5, sps6, sps7, sps8, sps9, sps10};
            int spsslen[] = {sizeof(sps1),sizeof(sps2), sizeof(sps3), sizeof(sps4), sizeof(sps5), sizeof(sps6), sizeof(sps7), sizeof(sps8), sizeof(sps9), sizeof(sps10)};
            double *bpts[] = {bpt1, bpt2, bpt3, bpt4, bpt5, bpt6, bpt7, bpt8, bpt9, bpt10};
            int N = 10;

            vector<vector<siftPoint> > vms;
            for(int i = 0; i < N; i++){
                //vector<siftPoint> vm = getSiftVec(allPatternFiles[i]);
                vector<siftPoint> vm = convertSiftPoint(spss[i], spsslen[i]/sizeof(siftPoint_save));
                vms.push_back(vm);
            }
            fprintf(stderr, "#pattern %d\n", clock());
            scaled = 1; //?????
            if(width > 1600){
                scaled = 2;
            }
            vector<int> scaleVector = scale(leftEyeImage, scaled);
            vector<siftPoint> vl = getSiftVec(scaleVector);

            scaleVector = scale(rightEyeImage, scaled);
            vector<siftPoint> vr = getSiftVec(scaleVector);

            fprintf(stderr, "#2pic %d\n", clock());

            pair<vector<point>, int> bestlp;
            pair<vector<point>, int> bestrp;
            bestlp.second = -1;
            bestrp.second = -1;

            //?????????
            vector<vector<pair<int, int> > > mls, mrs;
            for(int i = 0; i < vms.size(); i++){
                vector<siftPoint> &vm = vms[i];
                vector<pair<int, int> > ml = mymatch("lines_left.txt", "pl.txt", vm, vl);
                vector<pair<int, int> > mr = mymatch("lines_right.txt", "pr.txt", vm, vr);
                mls.push_back(ml);
                mrs.push_back(mr);
            }
            fprintf(stderr, "#match %d\n", clock());

            vector<int> order = getPatternOrder(mls);

            //for(int _i = 0; _i < order.size(); _i++){
            for(int _i = 0; _i < 2; _i++){
                int i = order[_i];
                vector<siftPoint> &vm = vms[i];
                double *bpt = bpts[i];
                vector<pair<int, int> > &ml = mls[i];

                int v1 = -1;
                if(ml.size() >= 4){
                    pair<vector<point>,int> tlp = randmatch(vm, vl, ml, bpt);
                    if(tlp.second > bestlp.second)
                        bestlp = tlp;
                    v1 = tlp.second;
                }

                fprintf(stderr, "%s %d/%d\n", allPatternFiles[i], v1, ml.size());
            }

            order = getPatternOrder(mrs);
            //for(int _i = 0; _i < order.size(); _i++){
            for(int _i = 0; _i < 2; _i++){
                int i = order[_i];
                vector<siftPoint> &vm = vms[i];
                double *bpt = bpts[i];
                vector<pair<int, int> > &mr = mrs[i];

                int v2 = -1;
                if(mr.size() >= 4){
                    pair<vector<point>,int> trp = randmatch(vm, vr, mr, bpt);
                    if(trp.second > bestrp.second)
                        bestrp = trp;
                    v2 = trp.second;
                }

                fprintf(stderr, "%s        %d/%d\n", allPatternFiles[i], v2, mr.size());
            }
            fprintf(stderr, "#match rnd %d\n", clock());

            
            //vector<pair<int, int> > mlr = mymatch("lines_left_right.txt", "plr.txt", vl, vr);
                
            vector<point> &lp = bestlp.first;
            vector<point> &rp = bestrp.first;

            /*if(mlr.size() >= 4){
                double pt[44];
                for(int i = 0; i < bestlp.first.size(); i++){
                    pt[i*2] = bestlp.first[i].x/2;
                    pt[i*2+1] = bestlp.first[i].y/2;
                }
                //scaled = 1;
                pair<vector<point>,int> tplr = randmatch(vl, vr, mlr, pt);
                fprintf(stderr, "l_r %d/%d\n", tplr.second, mlr.size());
                rp = tplr.first;
            }*/
            fflush(stderr);

            
            


            /*FILE *fout = fopen("tp.txt", "w");
            for(int i=0;i<vm.size();i++){
                double bx, by;
                bx = vm[i].x;
                by = vm[i].y;
                fprintf(fout, "%f %f %d\n", bx, by, i);
            }
            fclose(fout);*/
            
            //rec(leftEyeImage, "lines_left.txt");
            //rec(rightEyeImage, "lines_right.txt");

            for(int i = 0; i < 22; i++){
                char ch[100];
                if(lp[i].x >= width-1)
                    lp[i].x = width-1;
                if(lp[i].x <= 0)
                    lp[i].x = 0;
                if(rp[i].x >= width-1)
                    rp[i].x = width-1;
                if(rp[i].x <= 0)
                    rp[i].x = 0;
                if(lp[i].y >= height-1)
                    lp[i].y = height-1;
                if(lp[i].y <= 0)
                    lp[i].y = 0;
                if(rp[i].y >= height-1)
                    rp[i].y = height-1;
                if(rp[i].y <= 0)
                    rp[i].y = 0;
                    
                lretint.push_back( lp[i].x);
                lretint.push_back( lp[i].y);
                rretint.push_back( rp[i].x);
                rretint.push_back( rp[i].y);        
            }
        }

        vector<string> ret;
        char line[1000];
        
        int small_led_size = 15;
        if (is_sim) small_led_size = 15;
                
        
        fprintf(stderr, "0 1 pos %d %d\n",  rretint[0], rretint[1]);
        fflush(stderr);
        bool red_cover = true;
        if (abs(lretint[0] - lretint[4]) + abs(lretint[1] - lretint[5]) < 300) {
            red_cover &= local_red_cover_search(width, height, lred, lgreen, lblue, lretint[0], lretint[1], lretint[2], lretint[3]);
        }
        else if (abs(lretint[8] - lretint[4]) < 200) {
            red_cover &= local_red_cover_search(width, height, lred, lgreen, lblue, lretint[2], lretint[5], lretint[2], lretint[3]);
        }
        if (abs(rretint[0] - rretint[4]) + abs(rretint[1] - rretint[5]) < 300) {
            red_cover &= local_red_cover_search(width, height, rred, rgreen, rblue, rretint[0], rretint[1], rretint[2], rretint[3]);
        }
        else if (abs(rretint[8] - rretint[4]) < 200) {
            red_cover &= local_red_cover_search(width, height, rred, rgreen, rblue, rretint[2], rretint[5], rretint[2], rretint[3]);
        }
        
        if (is_sim) {
            //????????
            local_black_button_search(width, height, rred, rgreen, rblue, rretint[6], rretint[7]);
            local_black_button_search(width, height, lred, lgreen, lblue, lretint[6], lretint[7]);
        }
        
//        fprintf(stderr, "red_cover %d\n", red_cover);
//        fflush(stderr);
        
        ret.push_back("");
        ret.push_back("");
        ret.push_back("");
        /*
        if (check_led(height, width, lretint[4], lretint[5], small_led_size, lred, lgreen, lblue)
         || check_led(height, width, rretint[4], rretint[5], small_led_size, rred, rgreen, rblue) ) {
            sprintf(line, "UP,%d,%d,%d,%d", lretint[0], lretint[1], rretint[0], rretint[1]);
            ret.push_back(line);
            if (red_cover) sprintf(line, "UP,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            else sprintf(line, "DOWN,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[4], lretint[5], rretint[4], rretint[5]);
            ret.push_back(line);
        }
        else {
            sprintf(line, "DOWN,%d,%d,%d,%d", lretint[0], lretint[1], rretint[0], rretint[1]);
            ret.push_back(line);
            if (red_cover) sprintf(line, "UP,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            else sprintf(line, "DOWN,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[4], lretint[5], rretint[4], rretint[5]);
            ret.push_back(line);
        }
        */
        
        bool power_on = false;
        if (check_led(height, width, lretint[8], lretint[9], small_led_size, lred, lgreen, lblue)
         || check_led(height, width, rretint[8], rretint[9], small_led_size, rred, rgreen, rblue) ) {
            sprintf(line, "UP,%d,%d,%d,%d", lretint[6], lretint[7], rretint[6], rretint[7]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[8], lretint[9], rretint[8], rretint[9]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[10], lretint[11], rretint[10], rretint[11]);
            ret.push_back(line);
            power_on = true;
        }
        else if (check_led(height, width, lretint[10], lretint[11], small_led_size, lred, lgreen, lblue)
              || check_led(height, width, rretint[10], rretint[11], small_led_size, rred, rgreen, rblue) ) {
            sprintf(line, "DOWN,%d,%d,%d,%d", lretint[6], lretint[7], rretint[6], rretint[7]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[8], lretint[9], rretint[8], rretint[9]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[10], lretint[11], rretint[10], rretint[11]);
            ret.push_back(line);
            power_on = true;
        }
        else {
            sprintf(line, "CENTER,%d,%d,%d,%d", lretint[6], lretint[7], rretint[6], rretint[7]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[8], lretint[9], rretint[8], rretint[9]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[10], lretint[11], rretint[10], rretint[11]);
            ret.push_back(line);
        }
        
        for (int i = 0; i < 9; i++) {
            if (check_led(height, width, lretint[12+i*2], lretint[13+i*2], 19, lred, lgreen, lblue)
             || check_led(height, width, rretint[12+i*2], rretint[13+i*2], 19, rred, rgreen, rblue)) {
                sprintf(line, "ON,%d,%d,%d,%d", lretint[12+i*2], lretint[13+i*2], rretint[12+i*2], rretint[13+i*2]);
                ret.push_back(line);
                power_on = true;
            }
            else {
                sprintf(line, "OFF,%d,%d,%d,%d", lretint[12+i*2], lretint[13+i*2], rretint[12+i*2], rretint[13+i*2]);
                ret.push_back(line);
            }
        }
        
        if (check_led(height, width, lretint[32], lretint[33], small_led_size, lred, lgreen, lblue)
         || check_led(height, width, rretint[32], rretint[33], small_led_size, rred, rgreen, rblue)) {
            sprintf(line, "UP,%d,%d,%d,%d", lretint[30], lretint[31], rretint[30], rretint[31]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[32], lretint[33], rretint[32], rretint[33]);
            ret.push_back(line);
            power_on = true;
        }
        else {
            sprintf(line, "DOWN,%d,%d,%d,%d", lretint[30], lretint[31], rretint[30], rretint[31]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[32], lretint[33], rretint[32], rretint[33]);
            ret.push_back(line);
        }
        
        if (check_led(height, width, lretint[36], lretint[37], small_led_size, lred, lgreen, lblue)
         || check_led(height, width, rretint[36], rretint[37], small_led_size, rred, rgreen, rblue) ) {
            sprintf(line, "UP,%d,%d,%d,%d", lretint[34], lretint[35], rretint[34], rretint[35]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[36], lretint[37], rretint[36], rretint[37]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[38], lretint[39], rretint[38], rretint[39]);
            ret.push_back(line);
            power_on = true;
        }
        else if (check_led(height, width, lretint[38], lretint[39], small_led_size, lred, lgreen, lblue)
              || check_led(height, width, rretint[38], rretint[39], small_led_size, rred, rgreen, rblue) ) {
            sprintf(line, "DOWN,%d,%d,%d,%d", lretint[34], lretint[35], rretint[34], rretint[35]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[36], lretint[37], rretint[36], rretint[37]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[38], lretint[39], rretint[38], rretint[39]);
            ret.push_back(line);
            power_on = true;
        }
        else {
            sprintf(line, "CENTER,%d,%d,%d,%d", lretint[34], lretint[35], rretint[34], rretint[35]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[36], lretint[37], rretint[36], rretint[37]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[38], lretint[39], rretint[38], rretint[39]);
            ret.push_back(line);
        }
        
        if (check_led(height, width, lretint[42], lretint[43], small_led_size, lred, lgreen, lblue)
         || check_led(height, width, rretint[42], rretint[43], small_led_size, rred, rgreen, rblue)) {
            sprintf(line, "UP,%d,%d,%d,%d", lretint[40], lretint[41], rretint[40], rretint[41]);
            ret.push_back(line);
            sprintf(line, "ON,%d,%d,%d,%d", lretint[42], lretint[43], rretint[42], rretint[43]);
            ret.push_back(line);
            power_on = true;
        }
        else {
            sprintf(line, "DOWN,%d,%d,%d,%d", lretint[40], lretint[41], rretint[40], rretint[41]);
            ret.push_back(line);
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[42], lretint[43], rretint[42], rretint[43]);
            ret.push_back(line);
        }
        
        if (power_on || check_led(height, width, lretint[4], lretint[5], small_led_size, lred, lgreen, lblue)
         || check_led(height, width, rretint[4], rretint[5], small_led_size, rred, rgreen, rblue) ) {
            sprintf(line, "UP,%d,%d,%d,%d", lretint[0], lretint[1], rretint[0], rretint[1]);
            //ret.push_back(line);
            ret[0] = line;
            if (red_cover) sprintf(line, "UP,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            else sprintf(line, "DOWN,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            //ret.push_back(line);
            ret[1] = line;
            sprintf(line, "ON,%d,%d,%d,%d", lretint[4], lretint[5], rretint[4], rretint[5]);
            //ret.push_back(line);
            ret[2] = line;
        }
        else {
            sprintf(line, "DOWN,%d,%d,%d,%d", lretint[0], lretint[1], rretint[0], rretint[1]);
            //ret.push_back(line);
            ret[0] = line;
            if (red_cover) sprintf(line, "UP,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            else sprintf(line, "DOWN,%d,%d,%d,%d", lretint[2], lretint[3], rretint[2], rretint[3]);
            //ret.push_back(line);
            ret[1] = line;
            sprintf(line, "OFF,%d,%d,%d,%d", lretint[4], lretint[5], rretint[4], rretint[5]);
            //ret.push_back(line);
            ret[2] = line;
        }
        
        return ret;
    }
};

void checkAndConvert(){
    char *allPatternFiles[] = {"iss.pgm", "lab1-1.pgm", "lab1-3.pgm", "lab2-1.pgm", "lab2-2.pgm", "lab3-1.pgm", "lab3-2.pgm", "lab3-3.pgm", "lab1-2.pgm", "lab1-4.pgm"};
    int N = 10;

    int maxi = 0;
    int sum = 0;
    vector<vector<siftPoint> > vms;
    for(int i = 0; i < N; i++){
        vector<siftPoint> vm = getSiftVec(allPatternFiles[i]);
        printf("siftPoint_save sps%d[] = {\n", i+1);

        for(int j = 0;j < vm.size(); j++){
            siftPoint &sp = vm[j];
            printf("{%.2lf,%.2lf,%.2lf,%.2lf,\"", sp.x, sp.y, sp.sigma, sp.angle);
            for(int k=0;k<128;k++){
                sp.vec[k] *= 2000;
                if(sp.vec[k] >maxi)
                    maxi = sp.vec[k];
            }
            sum += 128;

            for(int k=0;k<128;k+=2){
                int num = (int)sp.vec[k] * 830 + (int)sp.vec[k+1];
                int t = num % 91 + 35;
                printf("%c", t=='\\'?'~':t);
                num /= 91;
                t = num % 91 + 35;
                printf("%c", t=='\\'?'~':t);
                num /= 91;
                t = num % 91 + 35;
                printf("%c", t=='\\'?'~':t);
            }
            printf("\"}");
            if(j != vm.size()-1)
                printf(",");
            printf("\n");
        }

        printf("};\n");

    }    
    printf("%d %d\n", sum, maxi);
}

