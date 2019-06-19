#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Original implementation in R code (R package "Bessel" v. 0.5-3) by        */
/*   Martin Maechler, Date: 23 Nov 2009, 13:39                               */
/*                                                                           */
/* Translated into C code by Kemal Dingic, Oct. 2011.                        */
/*                                                                           */
/* Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                 */
/*                                                                           */
/* Translated into C++ code by Peter Knaus, Mar. 2019.                       */
/*                                                                           */
/*---------------------------------------------------------------------------*/

double unur_bessel_k_nuasympt(double x, double nu, bool islog, bool expon_scaled){
  double M_LNPI = 1.14472988584940017414342735135;

  double z;
  double sz, t, t2, eta;
  double d, u1t,u2t,u3t,u4t;
  double res;


  z = x / nu;

  sz = hypot(1,z);
  t = 1. / sz;
  t2 = t*t;

  if (expon_scaled){
    eta = (1./(z + sz));
  } else {
    eta = sz;
  }

  eta += log(z) - log1p(sz);
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125.
                   + t2 * (-94121676.
                   + t2 * (349922430.
                   + t2 * (-446185740.
                   + t2 * 185910725.)))) / 39813120.;
                   d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;

  res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);

  if (islog){
    return res;
  } else {
    return exp(res);
  }
}
