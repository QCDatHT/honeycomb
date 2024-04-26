//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/model.h>
#include <TW3EV_include/ran2.h>

double model_zero_function(double x1, double x2, double x3, void *p) {
   (void)x1;
   (void)x2;
   (void)x3;
   (void)p;

   return 0.0; }

double Tu_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * cos(4.0 * x2);
}
double Td_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   double temp = (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
   return (2 - cos(3 * M_PI * temp)) * temp;
}
double Ts_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   return -0.3 * Td_test(x1, x2, x3, NULL);
}

double DTu_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   return (sin(x2 * M_PI) + 4 * (x1 * x1 - x3 * x3)) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}
double DTd_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(x2 * M_PI) * (2 - 2 * cos(3 * M_PI * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3))) / sqrt(r);
}
double DTs_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   return -0.3 * DTd_test(x1, x2, x3, NULL);
}

double Hu_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(x2 * M_PI) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) / (sqrt(r));
}
double Hd_test(double x1, double x2, double x3, void *p)
{
   (void)x1;
   (void)x2;
   (void)x3;
   (void)p;
   return 0;
}
double Hs_test(double x1, double x2, double x3, void *p)
{
   (void)x1;
   (void)x2;
   (void)x3;
   (void)p;
   return 0;
}

double Eu_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   return (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3);
}
double Ed_test(double x1, double x2, double x3, void *p)
{
   (void)x1;
   (void)x2;
   (void)x3;
   (void)p;
   return 0;
}
double Es_test(double x1, double x2, double x3, void *p)
{
   (void)x1;
   (void)x2;
   (void)x3;
   (void)p;
   return 0;
}

double TFp_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return sin(x1 - x3) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * sqrt(r);
}
double TFm_test(double x1, double x2, double x3, void *p)
{
   (void)p;
   double r = max3(fabs(x1), fabs(x2), fabs(x3));
   return cos(x1 - x3) * (1 - x1 * x1) * (1 - x2 * x2) * (1 - x3 * x3) * sqrt(r);
}

static bool check_symm_sign(double (*T)(double, double, double, void *), void *par, double s, char *distr, double xmin)
{

   static int64_t idum = -2;

   int32_t iBoot = 0;
   bool incorrect = false;

   while (iBoot < 1e+4) {
      double x1 = 2.0 * ran2(&idum) - 1.0;
      double x2 = 2.0 * ran2(&idum) - 1.0;
      double x3 = -x1 - x2;
      double r = max3(fabs(x1), fabs(x2), fabs(x3));
      if (r >= 1 || r <= xmin) continue;
      double val = T(x1, x2, x3, par);
      double val_sym = s * T(-x3, -x2, -x1, par);
      if (fabs(val - val_sym) >= 5 * _ZERO_THR_) tw3ev_log(TW3EV_WARNING, "Incorrect %s symmetry for (x1,x2,x3): (%.4e,%.4e,%.4e); violation: %.16e", distr, x1, x2, x3, fabs(val - val_sym));
      incorrect = incorrect || (fabs(val - val_sym) >= 5 * _ZERO_THR_);
      iBoot++;
   }
   return !incorrect;
}

static bool check_symm_gluon(double (*T)(double, double, double, void *), void *par, double s, char *distr, double xmin)
{
   static int64_t idum = -1;

   int32_t iBoot = 0;
   bool incorrect = false;

   while (iBoot < 1e+4) {
      double x1 = 2.0 * ran2(&idum) - 1.0;
      double x2 = 2.0 * ran2(&idum) - 1.0;
      double x3 = -x1 - x2;
      double r = max3(fabs(x1), fabs(x2), fabs(x3));
      if (r >= 1 || r <= xmin) continue;
      double val = T(x1, x2, x3, par);
      double val_sym = s * T(x3, x2, x1, par);
      if (fabs(val - val_sym) >= 5 * _ZERO_THR_) tw3ev_log(TW3EV_WARNING, "Incorrect %s symmetry for gluon-specific for (x1,x2,x3): (%.4e,%.4e,%.4e), violation: %.16e", distr, x1, x2, x3, fabs(val - val_sym));
      incorrect = incorrect || (fabs(val - val_sym) >= 5 * _ZERO_THR_);
      iBoot++;
   }
   return !incorrect;
}

bool check_T_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin)
{
   if (xmin < 0 || xmin >= 1) tw3ev_log(TW3EV_ERROR, "Invalid xmin (%.4e) in the check_T_symmetry for %s", xmin, distr);
   return check_symm_sign(T, par, 1, distr, xmin);
}
bool check_DT_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin)
{
   if (xmin < 0 || xmin >= 1) tw3ev_log(TW3EV_ERROR, "Invalid xmin (%.4e) in the check_DT_symmetry for %s", xmin, distr);
   return check_symm_sign(T, par, -1, distr, xmin);
}
bool check_TFp_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin)
{
   if (xmin < 0 || xmin >= 1) tw3ev_log(TW3EV_ERROR, "Invalid xmin (%.4e) in the check_TFp_symmetry for %s", xmin, distr);
   return check_symm_sign(T, par, 1, distr, xmin) && check_symm_gluon(T, par, -1, distr, xmin);
}
bool check_TFm_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin)
{
   if (xmin < 0 || xmin >= 1) tw3ev_log(TW3EV_ERROR, "Invalid xmin (%.4e) in the check_TFm_symmetry for %s", xmin, distr);
   return check_symm_sign(T, par, 1, distr, xmin) && check_symm_gluon(T, par, 1, distr, xmin);
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

