//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef KERNELS_H
#define KERNELS_H

#include <TW3EV_include/default.h>
#include <TW3EV_include/integration.h>
#include <TW3EV_include/sparse_matrix.h>
#include <TW3EV_include/thpool.h>

typedef enum { IT_PLUS = 0, IT_MINUS = 1, IT_BOTH = 2 } interpolant_type_e;

typedef enum { LOG_GRID, PWR_GRID, IML_GRID, HYP_GRID } nonlinear_radial_grid_type_e;

typedef enum { LIN_GRID, COS_GRID } angular_grid_type_e;

typedef struct {
   int32_t ip, jp;
   // int32_t i, j; // For debug only
   double x1, x2, x3, c_fact;
   int32_t N, M;
} integration_par_t;

typedef struct {
   void *ws_integr;
   integration_par_t *int_par;
} integration_wrapper_t;

typedef struct {
   double x1, x2, x3;
   // int32_t i, j;
   double x1min, x2min, x3min;
   double x1max, x2max, x3max;
} stored_point_t;

void kernels_set_grid_type(interpolant_type_e F_t, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t, double pge);
void kernels_setup_for_computation(int32_t N, int32_t M, double c_fact, int32_t n_threads, interpolant_type_e F_t, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t,
                                   double pge, const char basefolder[]);
void kernels_free_setup(int32_t N, int32_t M);
sparse_mat_t *init_kernel(double (*Hk)(int32_t, int32_t, int32_t, int32_t, integration_wrapper_t *, int), int32_t N, int32_t M, double c_fact, int nf, printout_level_e pl);


double Fijk_plus(double ix, double jy, double i, double j, int32_t N);
double Fijk_minus(double ix, double jy, double i, double j, int32_t N);
double Fijk_both(double ix, double jy, double i, double j, int32_t N);

double H_NS(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_CO(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_J_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);

double H_qq_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_qg_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_gq_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_gg_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);

double H_qq_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_qg_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_gq_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);
double H_gg_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf);

double execution_time(void (*fnc)(void *), void *p);

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

