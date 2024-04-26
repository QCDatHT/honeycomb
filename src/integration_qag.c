//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/integration.h>

static inline bool subinterval_too_small(double a1, double a2, double b2)
{
   const double tmp = (1 + _100_LOCAL_DBL_EPSILON_) * (fabs(a2) + _1000_LOCAL_DBL_MIN_);
   return fabs(a1) <= tmp && fabs(b2) <= tmp;
}

// Segregated prototypes, not exposed in the header
void gauss_kronrod_21(double (*fnc)(double, void *), void *p_fnc, double a, double b, double *result, double *abserr);
void gauss_kronrod_31(double (*fnc)(double, void *), void *p_fnc, double a, double b, double *result, double *abserr);
void gauss_kronrod_41(double (*fnc)(double, void *), void *p_fnc, double a, double b, double *result, double *abserr);
void gauss_kronrod_61(double (*fnc)(double, void *), void *p_fnc, double a, double b, double *result, double *abserr);

typedef void gauss_kronrod_rule(double (*fnc)(double, void *), void *p_fnc, double a, double b, double *result, double *abserr);

double qag_all_recursive_step(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, gauss_kronrod_rule *integration)
{
   double res = 0, err = 0;
   integration(fnc, p_fnc, a, b, &res, &err);
   const double center = (a + b) * 0.5;
   // check here also for zero integral! (Equivalent to just checking on the outside)
   // no depth check required, it is enough to check for the size of the interval
   // if (err <= epsabs || subinterval_too_small(a, center, b) || fabs(res) < _ZERO_THR_) return res;
   if (err <= epsabs || subinterval_too_small(a, center, b) ) return res;

   double res1 = qag_all_recursive_step(fnc, p_fnc, a, center, epsabs, integration);
   double res2 = qag_all_recursive_step(fnc, p_fnc, center, b, epsabs, integration);
   return res1 + res2;
}

void integration_qag(double (*fnc)(double, void *), void *p_fnc, double a, double b, double epsabs, double *result, tw3ev_integration_rule_e gk_type)
{
   gauss_kronrod_rule *integration = NULL;
   switch (gk_type) {
   case RULE_TW3EV_GK21:
      integration = gauss_kronrod_21;
      break;
   case RULE_TW3EV_GK31:
      integration = gauss_kronrod_31;
      break;
   case RULE_TW3EV_GK41:
      integration = gauss_kronrod_41;
      break;
   case RULE_TW3EV_GK61:
      integration = gauss_kronrod_61;
      break;
   default:
      integration = gauss_kronrod_21;
      break;
   }
   // double res = 0, err = 0;
   // integration(fnc, p_fnc, a, b, &res, &err);
   // // epsabs = max3(fabs((b - a) * 0.5 * res) * epsabs, 0.1 * EPS_TOLL, epsabs);
   // if (err < epsabs) {
   //    *result = res;
   //    return;
   // }
   epsabs = 1e-10;

   double temp = 0;
   const double width = (b - a);
   const size_t steps = 10;
// #pragma omp parallel for reduction(+ : temp)
   for (size_t il = 0; il < steps; il++)
      temp += qag_all_recursive_step(fnc, p_fnc, a + il * width / (double)steps, a + (il + 1) * width / (double)steps, epsabs, integration);
   // *result = qag_all_recursive_step(fnc, p_fnc, a, (a + b) * 0.5, epsabs, integration) + qag_all_recursive_step(fnc, p_fnc, (a + b) * 0.5, b, epsabs, integration);
   *result = temp;
}

 // Copyright (C) 2024 Simone Rodini; Lorenzo Rossi
 // This program is free software; you can redistribute it and/or modify
 // it under the terms of the GNU General Public License as published by
 // the Free Software Foundation; either version 2 of the License, or
 // (at your option) any later version.
 // 
 // This program is distributed in the hope that it will be useful,
 // but WITHOUT ANY WARRANTY; without even the implied warranty of
 // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 // GNU General Public License for more details.
 // 
 // You should have received a copy of the GNU General Public License along
 // with this program; if not, write to the Free Software Foundation, Inc.,
 // 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

