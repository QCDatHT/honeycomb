//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef KERNELS_COMMON_H
#define KERNELS_COMMON_H

#include <TW3EV_include/kernels.h>

// Some forward declarations
double Fij(double x1, double x2, double x3, int32_t i, int32_t j, int32_t N, int32_t M, double c_fact);

// Fundamental function of the hexagon defined in grid space

// Nc - ns
double Hhat12(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hhat23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hplus12(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hplus23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
// 1/Nc - ns
double Hhat13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hplus13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hminus12(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hminus23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double He23P23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);

// Singlet
// q-q
double Hd13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
// q-g
double Vp13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Vm13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
// g-q
double Wp13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Wp13P23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Wm13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Wm13P23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double DW13(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double DW13P23(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);

// g-g
double Hhat12GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hhat23GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hhat31GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);

double Hplus12GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hplus13GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);

double Htildeplus12GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Htildeplus13GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);

double Hminus12GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);
double Hminus13GG(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW);

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

