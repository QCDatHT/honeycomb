//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/kernels.h>
#include <TW3EV_include/kernels_common.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

extern stored_point_t **stored_points;

static double Hhat31GG_integrand_1_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;
   double FF = Fij(pa->x1 + v, pa->x2, pa->x3 - v, pa->ip, pa->jp, pa->N, pa->M, pa->c_fact);
   if (fabs(FF) < (_ZERO_THR_)) return 0;
   return FF * (SQ(pa->x3) / (v * SQ(pa->x3 - v)));
}

static double Hhat31GG_integrand_2_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;
   double FF = Fij(pa->x1 + v, pa->x2, pa->x3 - v, pa->ip, pa->jp, pa->N, pa->M, pa->c_fact);
   if (fabs(FF) < (_ZERO_THR_)) return 0;
   return FF * (-SQ(pa->x1) / (v * SQ(v + pa->x1)));
}

static double Hhat31GG_integrand_sub1_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;
   double FFd = ((pa->x3 / (pa->x3 - v)) * Fij(pa->x1 + v, pa->x2, pa->x3 - v, (int32_t)pa->ip, (int32_t)pa->jp, pa->N, pa->M, pa->c_fact) - 1);
   if (fabs(FFd) < (_ZERO_THR_)) return 0;
   return -FFd * (pa->x3 / (-v * (pa->x3 - v)));
}

static double Hhat31GG_integrand_sub2_v_def_sign(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;
   double FFd = ((pa->x1 / (pa->x1 + v)) * Fij(pa->x1 + v, pa->x2, pa->x3 - v, (int32_t)pa->ip, (int32_t)pa->jp, pa->N, pa->M, pa->c_fact) - 1);
   if (fabs(FFd) < (_ZERO_THR_)) return 0;
   return -FFd * (pa->x1 / (v * (pa->x1 + v)));
}
double Hhat31GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW)
{

   integration_par_t *int_par = IW->int_par;

   int_par->ip = ip;
   int_par->jp = jp;

   double res=0;

   int_par->x1 = stored_points[i][j].x1;
   int_par->x2 = stored_points[i][j].x2;
   int_par->x3 = stored_points[i][j].x3;

   double vmin, vmax;
   vmin = max2(stored_points[i][j].x3 - stored_points[ip][jp].x3max, stored_points[ip][jp].x1min - stored_points[i][j].x1);
   vmax = min2(stored_points[i][j].x3 - stored_points[ip][jp].x3min, stored_points[ip][jp].x1max - stored_points[i][j].x1);
   if (vmin >= vmax) return 0;

   if (i == ip && j == jp) {
      double res1 = 0;
      if (int_par->x3 > 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Hhat31GG_integrand_sub1_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res += res1;
      } else if (int_par->x3 < 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Hhat31GG_integrand_sub1_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res -= res1; //!
      }

      if (int_par->x1 > 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Hhat31GG_integrand_sub2_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res += res1;
      } else if (int_par->x1 < 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Hhat31GG_integrand_sub2_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res -= res1; //!
      }

      res += (stored_points[i][j].x3 >= 0 ? log(1 - stored_points[i][j].x3 / vmin) : log(1 - stored_points[i][j].x3 / vmax));
      res += (stored_points[i][j].x1 >= 0 ? log(1 + stored_points[i][j].x1 / vmax) : log(1 + stored_points[i][j].x1 / vmin));

      if (fabs(stored_points[i][j].x1) < (_ZERO_THR_)) res += 1;
      if (fabs(stored_points[i][j].x3) < (_ZERO_THR_)) res += 1;
   } else {
      double res1 = 0;

      if (int_par->x3 > 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Hhat31GG_integrand_1_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res += res1;
      } else if (int_par->x3 < 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Hhat31GG_integrand_1_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res -= res1; //!
      }

      if (int_par->x1 > 0 && vmax > 0) {
         double lower = max2(vmin, 0);
         integration(Hhat31GG_integrand_2_v_def_sign, (void *)int_par, lower, vmax, EPS_TOLL, &res1);
         res += res1;
      } else if (int_par->x1 < 0 && vmin < 0) {
         double upper = min2(0, vmax);
         integration(Hhat31GG_integrand_2_v_def_sign, (void *)int_par, vmin, upper, EPS_TOLL, &res1);
         res -= res1; //!
      }
   }
   return res;
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

