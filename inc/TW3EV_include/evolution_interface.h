//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef EVOLUTION_INTERFACE_H
#define EVOLUTION_INTERFACE_H


#include <TW3EV_include/default.h>
#include <TW3EV_include/diff_eq.h>
#include <TW3EV_include/integration.h>
#include <TW3EV_include/kernels.h>
#include <TW3EV_include/model.h>
#include <TW3EV_include/sparse_matrix.h>
#include <TW3EV_include/thpool.h>
#include <TW3EV_include/read_config.h>

typedef struct evolution_interface_t {
   printout_level_e pr_lev;
   interpolant_type_e F_t;
   tw3ev_integration_rule_e i_rule;

   nonlinear_radial_grid_type_e G_t;
   angular_grid_type_e AG_t;
   double grid_exponent;

   char *basefolder, *ker_subfolder, *res_subfolder;
   double xmin, c_fact;
   // Function pointer to computation of as = \alpha_s/4\pi at a given
   // << Log[\mu^2] >> point
   double (*as)(double, void *);
   // Parameters for as
   void *p_as;
   // Internal runge-kutta pointer - own the memory
   rk4th_internal_t *rk;
   // Initial, final and step scales
   double t0, tF, dt;
   int32_t nstep;
   // nf
   double nf;
   double bm_threshold, ch_threshold;

   // Select whether chiral-odd and even evolutins should be active
   bool chiral_odd_active, chiral_even_active;

   // Grid dimension
   int32_t N, M;
   // Non-singlet kernel - own the memory
   sparse_mat_t *H_NS;
   sparse_mat_t *H_CO;
   // Singlet kernels - own the memory
   sparse_mat_t *H_qq_p, *H_qg_p, *H_gq_p, *H_gg_p;
   sparse_mat_t *H_qq_m, *H_qg_m, *H_gq_m, *H_gg_m;
   sparse_mat_t *H_J_p; // 4nf H^d_{13}

   // Solution at each step - own the memory
   // evolution basis: u-d = f1; u+d-2s = f2; u+d+s = fq
   vector_t *f_1_p, *f_2_p, *f_3_p, *f_4_p;
   vector_t *f_1_m, *f_2_m, *f_3_m, *f_4_m;
   vector_t *f_q_p, *f_q_m, *f_g_p, *f_g_m;
   vector_t *h_up, *h_dn, *h_st, *e_up, *e_dn, *e_st;
   // Number of threads for computation
   int32_t n_thread;

} evolution_interface_t;

typedef struct {
   double (*Tu)(double, double, double, void *);
   double (*Td)(double, double, double, void *);
   double (*Ts)(double, double, double, void *);

   double (*DTu)(double, double, double, void *);
   double (*DTd)(double, double, double, void *);
   double (*DTs)(double, double, double, void *);

   double (*Hu)(double, double, double, void *);
   double (*Hd)(double, double, double, void *);
   double (*Hs)(double, double, double, void *);

   double (*Eu)(double, double, double, void *);
   double (*Ed)(double, double, double, void *);
   double (*Es)(double, double, double, void *);

   double (*TFp)(double, double, double, void *);
   double (*TFm)(double, double, double, void *);

   void *p_PDF;
} evolution_functions_t;

typedef struct {
   // Evolution basis
   vector_t *f_1_p, *f_2_p, *f_3_p, *f_4_p;
   vector_t *f_1_m, *f_2_m, *f_3_m, *f_4_m;
   vector_t *f_q_p, *f_q_m, *f_g_p, *f_g_m;

   // Definite-C-parity basis
   vector_t *up_p, *dn_p, *st_p, *ch_p, *bm_p, *gl_p;
   vector_t *up_m, *dn_m, *st_m, *ch_m, *bm_m, *gl_m;
   vector_t *e_up, *e_dn, *e_st;
   vector_t *h_up, *h_dn, *h_st;

   // Physical basis
   vector_t *T_up, *T_dn, *T_st, *T_ch, *T_bm, *Tp_gl;
   vector_t *DT_up, *DT_dn, *DT_st, *DT_ch, *DT_bm, *Tm_gl;

   // Info
   int32_t N, M;
   int nf;
   char *basefolder, *res_subfolder;
   double xmin, c_fact, grid_exponent, mu_2;
   interpolant_type_e F_t;
   nonlinear_radial_grid_type_e G_t;
   angular_grid_type_e AG_t;
   bool chiral_odd_active, chiral_even_active;

} saved_solution_t;

typedef enum { TU_f, DTU_f, TD_f, DTD_f, TS_f, DTS_f, TC_f, DTC_f, TB_f, DTB_f, TFP_f, TFM_f, HU_f, HD_f, HS_f, EU_f, ED_f, ES_f } function_type_e;
typedef enum { BASIS_PHYSICAL, BASIS_DEFINITE_C_PAR, BASIS_BOTH, BASIS_EVOLUTION } saved_basis_e;

typedef struct {
   // Generic info
   char *basefolder, *ker_subfolder, *res_subfolder;
   int32_t N, M, n_thread;

   // Physical information: initial scale (log(\mu^2)), final scale, step size
   // and heavy-quark thresholds
   double mu0_2, muF_2;
   int32_t nstep;
   double ch_thr_mu2, bm_thr_mu2;

   // Select whether chiral-odd evolution should be active or not
   bool chiral_odd_active, chiral_even_active;

   // Select which version. Probably to be removed and automatically compute the
   // average
   interpolant_type_e F_t;

   // For the type of grid
   double xmin;
   nonlinear_radial_grid_type_e G_t;
   angular_grid_type_e AG_t;
   double grid_exponent; // from grid to physical!

   // For test purposes
   char *model;

   tw3ev_integration_rule_e i_rule;
} input_parameters_t;

input_parameters_t *init_input_parameters(const char fn[], bool save_conf_out);
void supply_default_parameters(input_parameters_t *in_par);

// int read_configuration_file(input_parameters_t *dc, const char fn[]);
// int write_configuration_file(input_parameters_t *dc, const char fn[]);

evolution_functions_t *init_default_evolution_function(evolution_functions_t *);

evolution_interface_t *init_evolution_interface(double (*as)(double, void *), void *p_as, evolution_functions_t *ev, input_parameters_t *in_par, printout_level_e pl, evolution_interface_t *ei, bool start_from_previous);
void free_evolution_interface(evolution_interface_t **ei);
void reset_scales(evolution_interface_t *ei, double mu0_2, double muF_2, int32_t nstep);
saved_solution_t *execute_evolution(evolution_interface_t *ei);
saved_solution_t **execute_evolution_w_steps(evolution_interface_t *ei);

void save_model(saved_solution_t *sol, const char basename[], saved_basis_e SB);
saved_solution_t *load_model(input_parameters_t *in_par, const char basename[], int nf, saved_basis_e SB);

saved_solution_t *init_saved_solution(evolution_interface_t *ei);
saved_solution_t *init_saved_solution_for_tests(evolution_functions_t *ev, input_parameters_t *in_par, printout_level_e pl);
void free_saved_solution(saved_solution_t **);

double get_interpolated_value(vector_t *distr, double x1, double x2, double x3, saved_solution_t *sol);

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

