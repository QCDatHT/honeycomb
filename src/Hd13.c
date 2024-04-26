//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/kernels.h>
#include <TW3EV_include/kernels_common.h>
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

extern stored_point_t **stored_points;

static double Hd13_integrand(double v, void *p)
{
   integration_par_t *pa = (integration_par_t *)p;
   double FF = Fij(pa->x1 - v, pa->x2, pa->x3 + v, pa->ip, pa->jp, pa->N, pa->M, pa->c_fact);
   return FF;
}

double Hd13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW)
{
   double th = -_Theta2_(stored_points[i][j].x1, stored_points[i][j].x3);
   if (fabs(th) < (_ZERO_THR_) || fabs(stored_points[i][j].x2) < (_ZERO_THR_)) return 0;
   integration_par_t *int_par = IW->int_par;

   int_par->ip = ip;
   int_par->jp = jp;

   double res = 0;

   int_par->x1 = stored_points[i][j].x1;
   int_par->x2 = stored_points[i][j].x2;
   int_par->x3 = stored_points[i][j].x3;

   double vmin, vmax;
   vmin = max2(stored_points[i][j].x1 - stored_points[ip][jp].x1max, stored_points[ip][jp].x3min - stored_points[i][j].x3);
   vmax = min2(stored_points[i][j].x1 - stored_points[ip][jp].x1min, stored_points[ip][jp].x3max - stored_points[i][j].x3);
   if (vmin >= vmax) return 0;

   integration(Hd13_integrand, (void *)int_par, vmin, vmax, EPS_TOLL, &res);

   return res * th * int_par->x1 * int_par->x3 / (CU(int_par->x2));
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

