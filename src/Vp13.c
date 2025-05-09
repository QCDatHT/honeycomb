//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <honeycomb/kernels.h>
#include <honeycomb/kernels_common.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

extern stored_point_t **stored_points;

double subtr = 0;

static double Vp13_integrand_x2null_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;
   double FF = Fij(pa->x1 - v, pa->x2, pa->x3 + v, pa->ip, pa->jp, pa->N, pa->M, pa->c_fact);
   if (fabs(FF) < (_ZERO_THR_)) return 0;
   return FF * (SQ(pa->x1)) / (SQ(pa->x1 - v) * SQ(pa->x1 - v));
}
static double Vp13_integrand_1_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;

   double FF = Fij(pa->x1 - v, pa->x2, pa->x3 + v, pa->ip, pa->jp, pa->N, pa->M, pa->c_fact) - subtr;
   if (fabs(FF) < (_ZERO_THR_)) return 0;
   return FF * (3 * pa->x1 + pa->x3 - 2 * v) / (SQ(v - pa->x1));
}
static double Vp13_integrand_2_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;

   double FF = Fij(pa->x1 - v, pa->x2, pa->x3 + v, pa->ip, pa->jp, pa->N, pa->M, pa->c_fact) - subtr;
   if (fabs(FF) < (_ZERO_THR_)) return 0;
   return -FF * (3 * pa->x3 + pa->x1 + 2 * v) / (SQ(v + pa->x3));
}

double Vp13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW)
{

   integration_par_t *int_par = IW->int_par;

   int_par->ip = ip;
   int_par->jp = jp;
   if (i == ip && j == jp) subtr = 1;
   else subtr = 0;

   double res = 0;
   ;

   int_par->x1 = stored_points[i][j].x1;
   int_par->x2 = stored_points[i][j].x2;
   int_par->x3 = stored_points[i][j].x3;

   const double x1 = stored_points[i][j].x1;
   const double x2 = stored_points[i][j].x2;
   const double x3 = stored_points[i][j].x3;

   if (fabs(int_par->x1) < _ZERO_THR_) {
      if (i == ip && j == jp) return -2.0 / int_par->x3;
      else return 0;
   }
   if (fabs(int_par->x3) < _ZERO_THR_) {
      if (i == ip && j == jp) return 2.0 / int_par->x1;
      else return 0;
   }

   double vmin, vmax;
   vmin = max2(stored_points[i][j].x1 - stored_points[ip][jp].x1max, stored_points[ip][jp].x3min - stored_points[i][j].x3);
   vmax = min2(stored_points[i][j].x1 - stored_points[ip][jp].x1min, stored_points[ip][jp].x3max - stored_points[i][j].x3);
   if (vmin >= vmax) return 0;

   // integration(Vp13_integrand, (void *)int_par, vmin, vmax, EPS_TOLL, &res);
   if (fabs(int_par->x2) < _ZERO_THR_) {
      double res1 = 0;
      if (int_par->x1 > 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Vp13_integrand_x2null_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res += res1;
      } else if (int_par->x1 < 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Vp13_integrand_x2null_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res -= res1; //!
      }
      return res;
   } else {
      double res1 = 0;
      if (int_par->x1 > 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Vp13_integrand_1_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res += res1;
         if (i == ip && j == jp) res += ((upper - vmin) * (x1 + x3)) / ((upper - x1) * (vmin - x1)) - 2 * log(-upper + x1) + 2 * log(-vmin + x1);
      } else if (int_par->x1 < 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Vp13_integrand_1_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res -= res1; //!
         if (i == ip && j == jp) res += -(((lower - vmax) * (x1 + x3)) / ((lower - x1) * (-vmax + x1))) - 2 * log(lower - x1) + 2 * log(vmax - x1);
      }

      if (int_par->x3 > 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Vp13_integrand_2_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res += res1;
         if (i == ip && j == jp) res += ((lower - vmax) * (x1 + x3)) / ((lower + x3) * (vmax + x3)) + 2 * log(lower + x3) - 2 * log(vmax + x3);
      } else if (int_par->x3 < 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Vp13_integrand_2_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res -= res1; //!
         if (i == ip && j == jp) res += ((upper - vmin) * (x1 + x3)) / ((upper + x3) * (vmin + x3)) + 2 * log((upper + x3) / (vmin + x3));
      }
   }
   const double pref = x1 * x3 / CU(x2);
   return pref * res;
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
