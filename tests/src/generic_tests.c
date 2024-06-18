//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
// Author: Lorenzo Rossi <mailto:lorenzo.rossi04@universitadipavia.it>
//
#ifdef __cplusplus
extern "C" {
#endif

#include <TW3EV_include/tw3ev.h>
#include <stdbool.h>

#define NTHREAD_LOC 10

typedef struct {
   double t_ch, t_bm; // charm and bottom thresholds
} as_parameters_t;

// this must be \alpha_s / 4\pi
double as(double t, void *p)
{
   (void)p;
   return 1.0 / (11.0 * (t + 3.));
}

double as_LO_GRV98(double t, void *p)
{
   as_parameters_t *ap = (as_parameters_t *)p;
   int32_t nf = 3;
   if (t > ap->t_ch && t <= ap->t_bm) nf = 4;
   else if (t > ap->t_bm) nf = 5;

   double lambda = 0;
   if (nf == 3) lambda = 0.204;
   if (nf == 4) lambda = 0.175;
   if (nf == 5) lambda = 0.132;

   double beta0 = 11. - 2.0 * nf / 3.0;

   double L2 = 2.0 * log(lambda);
   return 1.0 / (beta0 * (t - L2));
}

void init_default_evolution_function_local(evolution_functions_t *ev)
{
   ev->Tu = model_zero_function;
   ev->Td = model_zero_function;
   ev->Ts = model_zero_function;

   ev->DTu = model_zero_function;
   ev->DTd = model_zero_function;
   ev->DTs = model_zero_function;

   ev->Hu = Hu_test;
   ev->Hd = Hu_test;
   ev->Hs = model_zero_function;

   ev->Eu = Eu_test;
   ev->Ed = Eu_test;
   ev->Es = model_zero_function;

   ev->TFp = model_zero_function;
   ev->TFm = model_zero_function;

   ev->p_PDF = NULL;
   // return ev;
}

evolution_functions_t init_model_from_config(input_parameters_t *in_par)
{

   evolution_functions_t ev = {0};
   init_default_evolution_function_local(&ev);
   if (strcmp(in_par->model, "configA") == 0) {
      ev.Tu = Tu_test;
      ev.DTu = DTu_test;
      ev.TFp = TFp_test;
      ev.TFm = TFm_test;
   } else if (strcmp(in_par->model, "configB") == 0) {
      ev.Tu = Tu_test;
      ev.DTu = DTu_test;
      ev.Td = Td_test;
      ev.DTd = DTd_test;
      ev.TFp = TFp_test;
      ev.TFm = TFm_test;
   } else if (strcmp(in_par->model, "configC") == 0) {
      ev.Tu = Tu_test;
      ev.DTu = DTu_test;
      ev.TFp = TFp_test;
      ev.TFm = TFm_test;
      ev.Td = Td_test;
      ev.DTd = DTd_test;
      ev.Ts = Ts_test;
      ev.DTs = DTs_test;
   } else tw3ev_log(TW3EV_ERROR, "Undefined model for constructing evolution functions!");

   return ev;
}

evolution_interface_t *exec_given_config(const char *conf_name, evolution_interface_t *ei, const char *specname)
{
   input_parameters_t *in_par = init_input_parameters(conf_name, true);
   evolution_functions_t ev = init_model_from_config(in_par);
   double x = 3.14159;
   ev.p_PDF = &x;

   ei = init_evolution_interface(as, NULL, &ev, in_par, PL_ESSENTIAL, ei, false);

   saved_solution_t *solOut = execute_evolution(ei);

   save_model(solOut, specname, BASIS_BOTH);

   free_saved_solution(&solOut);

   free(in_par->basefolder);
   free(in_par->ker_subfolder);
   free(in_par->res_subfolder);
   free(in_par->model);
   free(in_par);
   return ei;
}

evolution_interface_t *compute_back_and_forth(const char *conf_name, evolution_interface_t *ei)
{
   input_parameters_t *in_par = init_input_parameters(conf_name, true);
   evolution_functions_t ev = init_model_from_config(in_par);

   ei = init_evolution_interface(as, NULL, &ev, in_par, PL_ESSENTIAL, ei, false);
   saved_solution_t *solIn = init_saved_solution(ei);
   save_model(solIn, "initial", BASIS_PHYSICAL);

   saved_solution_t *solOut = execute_evolution(ei);
   free_saved_solution(&solOut);
   tw3ev_log(TW3EV_WARNING, "Swapping scales");
   double temp = in_par->mu0_2;
   in_par->mu0_2 = in_par->muF_2;
   in_par->muF_2 = temp;
   ei = init_evolution_interface(as, NULL, &ev, in_par, PL_ESSENTIAL, ei, true);

   saved_solution_t *sol1 = execute_evolution(ei);
   save_model(sol1, "BandF", BASIS_PHYSICAL);

   free_saved_solution(&sol1);

   free(in_par->basefolder);
   free(in_par->ker_subfolder);
   free(in_par->res_subfolder);
   free(in_par->model);
   free(in_par);
   return ei;
}

void check_model_symmetries()
{
   evolution_functions_t *ev = NULL;
   ev = init_default_evolution_function(ev);

   bool a = check_T_symmetry(ev->Tu, NULL, "Tu", 0.01);

   a = a && check_T_symmetry(ev->Td, NULL, "Td", 0.01);
   a = a && check_T_symmetry(ev->Ts, NULL, "Ts", 0.01);

   a = a && check_T_symmetry(ev->Eu, NULL, "Eu", 0.01);
   a = a && check_T_symmetry(ev->Ed, NULL, "Ed", 0.01);
   a = a && check_T_symmetry(ev->Es, NULL, "Es", 0.01);

   a = a && check_DT_symmetry(ev->DTu, NULL, "DTu", 0.01);
   a = a && check_DT_symmetry(ev->DTd, NULL, "DTd", 0.01);
   a = a && check_DT_symmetry(ev->DTs, NULL, "DTs", 0.01);

   a = a && check_DT_symmetry(ev->Hu, NULL, "Hu", 0.01);
   a = a && check_DT_symmetry(ev->Hd, NULL, "Hd", 0.01);
   a = a && check_DT_symmetry(ev->Hs, NULL, "Hs", 0.01);

   a = a && check_TFp_symmetry(ev->TFp, NULL, "TFp", 0.01);
   a = a && check_TFm_symmetry(ev->TFm, NULL, "TFm", 0.01);

   if (a) tw3ev_log(TW3EV_INFO, "Ok");
   else tw3ev_log(TW3EV_WARNING, "Not ok");
}

int main(void)
{

   check_model_symmetries();

   evolution_interface_t *ei = NULL;
   tw3ev_log(TW3EV_WARNING, "Start first test...");
   ei = exec_given_config("../config_1.in", ei, "result");
   tw3ev_log(TW3EV_WARNING, "... end first test.");

   tw3ev_log(TW3EV_WARNING, "Start second test...");
   ei = exec_given_config("../config_2.in", ei, "result");
   tw3ev_log(TW3EV_WARNING, "... end second test.");

   tw3ev_log(TW3EV_WARNING, "Start third test...");
   ei = compute_back_and_forth("../config_3.in", ei);
   tw3ev_log(TW3EV_WARNING, "... end third test.");

   tw3ev_log(TW3EV_WARNING, "Start fourth test...");
   ei = exec_given_config("../config_4.in", ei, "result");
   tw3ev_log(TW3EV_WARNING, "... end fourth test.");

   tw3ev_log(TW3EV_WARNING, "Start fifth test...");
   ei = exec_given_config("../config_5.in", ei, "result");
   tw3ev_log(TW3EV_WARNING, "... end fifth test");

   free_evolution_interface(&ei);
   return 0;
}

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
