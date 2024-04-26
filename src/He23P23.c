//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/kernels.h>
#include <TW3EV_include/kernels_common.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

extern stored_point_t **stored_points;

static double He23P23_integrand_v_def_sign(double v, void *p)
{
   // \theta(x3)\theta(-v)
   integration_par_t *pa = (integration_par_t *)p;
   double FF = (Fij(pa->x1, pa->x3 - v, pa->x2 + v, (int32_t)pa->ip, (int32_t)pa->jp, pa->N, pa->M, pa->c_fact));
   if (fabs(FF) < (_ZERO_THR_)) return 0;
   return FF * pa->x3 / (SQ(pa->x3 - v));
}

double He23P23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW)
{

   integration_par_t *pa = IW->int_par;

   pa->ip = ip;
   pa->jp = jp;

   double res=0;

   pa->x1 = stored_points[i][j].x1;
   pa->x2 = stored_points[i][j].x2;
   pa->x3 = stored_points[i][j].x3;

   double vmin, vmax;
   vmin = max2(stored_points[i][j].x3 - stored_points[ip][jp].x2max, stored_points[ip][jp].x3min - stored_points[i][j].x2);
   vmax = min2(stored_points[i][j].x3 - stored_points[ip][jp].x2min, stored_points[ip][jp].x3max - stored_points[i][j].x2);
   if (vmin >= vmax) return 0;

   ;

   // \delta(x3) contribution explicitly
   if (fabs(stored_points[i][j].x3) < (_ZERO_THR_)) return Fij(pa->x1, pa->x3, pa->x2, ip, jp, pa->N, pa->M, pa->c_fact);

   double res1 = 0;

   if (pa->x3 > 0 && vmin < 0) {
      double upper = min2(0, vmax);
      integration(He23P23_integrand_v_def_sign, (void *)pa, vmin, upper, EPS_TOLL, &res1);
      res += res1;
   } else if (pa->x3 < 0 && vmax > 0) {
      double lower = max2(vmin, 0);
      integration(He23P23_integrand_v_def_sign, (void *)pa, lower, vmax, EPS_TOLL, &res1);
      res -= res1; //!
   }
   // integration(He23P23_integrand, (void *)pa, vmin, vmax, EPS_TOLL, &res);
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

