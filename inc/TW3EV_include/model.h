//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef MODEL_H
#define MODEL_H

#include <TW3EV_include/default.h>

double model_zero_function(double x1, double x2, double x3, void *p);

double Td_test(double x1, double x2, double x3, void *p);
double Tu_test(double x1, double x2, double x3, void *p);
double Ts_test(double x1, double x2, double x3, void *p);
double DTu_test(double x1, double x2, double x3, void *p);
double DTd_test(double x1, double x2, double x3, void *p);
double DTs_test(double x1, double x2, double x3, void *p);

double Hu_test(double x1, double x2, double x3, void *p);
double Hd_test(double x1, double x2, double x3, void *p);
double Hs_test(double x1, double x2, double x3, void *p);
double Eu_test(double x1, double x2, double x3, void *p);
double Ed_test(double x1, double x2, double x3, void *p);
double Es_test(double x1, double x2, double x3, void *p);
double TFp_test(double x1, double x2, double x3, void *p);
double TFm_test(double x1, double x2, double x3, void *p);

bool check_T_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin);
bool check_DT_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin);
bool check_TFp_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin);
bool check_TFm_symmetry(double (*T)(double, double, double, void *), void *par, char *distr, double xmin);

#endif
#ifdef __cplusplus
}
#endif
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

