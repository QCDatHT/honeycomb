//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DIFF_EQ_H
#define DIFF_EQ_H

#include <TW3EV_include/default.h>
#include <TW3EV_include/kernels.h>
#include <TW3EV_include/sparse_matrix.h>

typedef struct {
   // Function pointer to computation of as = \alpha_s/4\pi at a given
   // << Log[\mu^2] >> point
   double (*as)(double, void *);
   // Parameters for as
   void *p_as;
   // Kernel non-singlet - do not own the memory
   sparse_mat_t *H_NS;
   sparse_mat_t *H_CO;
   // Kernel singlet - do not own the memory
   sparse_mat_t *H_qq_p, *H_qg_p, *H_gq_p, *H_gg_p;
   sparse_mat_t *H_qq_m, *H_qg_m, *H_gq_m, *H_gg_m;
   sparse_mat_t *H_J_p; // 4nf H^d_{13}

   double *nf;
   double bm_threshold, ch_threshold;
   bool chiral_odd_active, chiral_even_active;

   // Internal rk4th middle terms non-singlet - own the memory
   vector_t *k1, *k2, *k3, *k4, *temp;
   // Internal rk4th middle terms singlet - own the memory
   vector_t *k1_q, *k2_q, *k3_q, *k4_q, *temp_q;
   vector_t *k1_g, *k2_g, *k3_g, *k4_g, *temp_g;
   // Log[\mu^2] and step in Log[\mu^2]
   double t, dt, dth;
   // Solution at each step - do not own the memory
   vector_t *f_1_p, *f_2_p, *f_3_p, *f_4_p;
   vector_t *f_1_m, *f_2_m, *f_3_m, *f_4_m;
   vector_t *f_q_p, *f_q_m, *f_g_p, *f_g_m;
   vector_t *h_up, *h_dn, *h_st, *e_up, *e_dn, *e_st;
   vector_t **which_CO;
   int32_t how_many_CO;

} rk4th_internal_t;

void rescale_nf_kernels(rk4th_internal_t *rk, double curr_nf, double new_nf);

rk4th_internal_t *init_rk4th(vector_t *f0_1_p, vector_t *f0_1_m, vector_t *f0_2_p, vector_t *f0_2_m, vector_t *f0_3_p, vector_t *f0_3_m, vector_t *f0_4_p, vector_t *f0_4_m, vector_t *f0_q_p, vector_t *f0_g_p, vector_t *f0_q_m, vector_t *f0_g_m,
                             vector_t *h_up, vector_t *h_dn, vector_t *h_st, vector_t *e_up, vector_t *e_dn, vector_t *e_st, sparse_mat_t *H_NS, sparse_mat_t *H_CO, sparse_mat_t *H_J_p, sparse_mat_t *H_qq_p, sparse_mat_t *H_qg_p,
                             sparse_mat_t *H_gq_p, sparse_mat_t *H_gg_p, sparse_mat_t *H_qq_m, sparse_mat_t *H_qg_m, sparse_mat_t *H_gq_m, sparse_mat_t *H_gg_m, double (*as)(double, void *), void *p_as, double t0, double dt, double *nf, bool ce_act,
                             bool co_act, double ch_thr, double bt_thr, rk4th_internal_t *rk);

void free_rk4th(rk4th_internal_t **rk);
void rk4th_step(rk4th_internal_t *rk);
void reset_scale_rk4th(rk4th_internal_t *rk, double t, double dt);

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

