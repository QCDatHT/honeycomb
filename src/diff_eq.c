//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/diff_eq.h>

void rescale_nf_kernels(rk4th_internal_t *rk, double curr_nf, double new_nf)
{
   // Check for no need of rescaling
   if (fabs(curr_nf - new_nf) < _ZERO_THR_) return;
   tw3ev_log(TW3EV_INFO, "Rescaling kernel from %d to %d nf", (int)round(curr_nf), (int)round(new_nf));

   double nf_ratio = (new_nf / (curr_nf));
   double beta_diff = (2.0 / 3.0) * (new_nf - curr_nf);

   for (int32_t i = 0; i < rk->H_J_p->size; i++)
      rk->H_J_p->data[i] = nf_ratio * rk->H_J_p->data[i];

   for (int32_t i = 0; i < rk->H_qg_p->size; i++)
      rk->H_qg_p->data[i] = nf_ratio * rk->H_qg_p->data[i];

   for (int32_t i = 0; i < rk->H_qg_m->size; i++)
      rk->H_qg_m->data[i] = nf_ratio * rk->H_qg_m->data[i];



   for (int32_t i = 0; i < rk->H_gg_p->size; i++) {
      if (rk->H_gg_p->i_c[i] == rk->H_gg_p->i_r[i]) rk->H_gg_p->data[i] = rk->H_gg_p->data[i] + beta_diff;
   }
   

   for (int32_t i = 0; i < rk->H_gg_m->size; i++) {
      if (rk->H_gg_m->i_c[i] == rk->H_gg_m->i_r[i]) rk->H_gg_m->data[i] = rk->H_gg_m->data[i] + beta_diff;
   }
}

static void rk4th_step_parallel(rk4th_internal_t *rk, vector_t *f, sparse_mat_t *Hk);
static void rk4th_step_parallel_singlet(rk4th_internal_t *rk, vector_t *f_q, vector_t *f_g, sparse_mat_t *H_qq, sparse_mat_t *H_qg, sparse_mat_t *H_gq, sparse_mat_t *H_gg, sparse_mat_t *H_J_p);

void reset_scale_rk4th(rk4th_internal_t *rk, double t, double dt)
{
   rk->t = t;
   rk->dt = dt;
   rk->dth = 0.5 * dt;

   int nf = 3;
   nf += (t >= rk->ch_threshold ? 1 : 0);
   nf += (t >= rk->bm_threshold ? 1 : 0);

   if (nf != (int)round(*(rk->nf))) {
      rescale_nf_kernels(rk, *(rk->nf), nf);
      *(rk->nf) = nf;
   }
}

bool to_debug = false;

void (*rk4th_step_CE_fixed_nf)(rk4th_internal_t *rk);

static void rk4th_step_CE_nf_3(rk4th_internal_t *rk);
static void rk4th_step_CE_nf_4(rk4th_internal_t *rk);
static void rk4th_step_CE_nf_5(rk4th_internal_t *rk);

rk4th_internal_t *init_rk4th(vector_t *f0_1_p, vector_t *f0_1_m, vector_t *f0_2_p, vector_t *f0_2_m, vector_t *f0_3_p, vector_t *f0_3_m, vector_t *f0_4_p, vector_t *f0_4_m, vector_t *f0_q_p, vector_t *f0_g_p, vector_t *f0_q_m, vector_t *f0_g_m,
                             vector_t *h_up, vector_t *h_dn, vector_t *h_st, vector_t *e_up, vector_t *e_dn, vector_t *e_st, sparse_mat_t *H_NS, sparse_mat_t *H_CO, sparse_mat_t *H_J_p, sparse_mat_t *H_qq_p, sparse_mat_t *H_qg_p,
                             sparse_mat_t *H_gq_p, sparse_mat_t *H_gg_p, sparse_mat_t *H_qq_m, sparse_mat_t *H_qg_m, sparse_mat_t *H_gq_m, sparse_mat_t *H_gg_m, double (*as)(double, void *), void *p_as, double t0, double dt, double *nf, bool ce_act,
                             bool co_act, double ch_thr, double bt_thr,  rk4th_internal_t *rk)
{
   bool to_initialize = rk == NULL;
   int32_t nn = 0;
   if (ce_act) nn = f0_1_p->size;
   else if (co_act) {
      if (NULL != h_up) nn = h_up->size;
      else if (NULL != h_dn) nn = h_dn->size;
      else if (NULL != h_st) nn = h_st->size;
      else if (NULL != e_up) nn = e_up->size;
      else if (NULL != e_dn) nn = e_dn->size;
      else if (NULL != e_st) nn = e_st->size;
   } else {
      tw3ev_log(TW3EV_ERROR, "Trying to initialize rkth with no valid boundary conditions.");
      return NULL;
   }

   if (to_initialize) rk = (rk4th_internal_t *)calloc(1, sizeof(rk4th_internal_t));

   rk->nf = nf;
   rk->chiral_even_active = ce_act;
   rk->chiral_odd_active = co_act;

   rk->ch_threshold = ch_thr;
   rk->bm_threshold = bt_thr;

   rk->H_NS = H_NS;
   rk->H_CO = H_CO;

   rk->H_J_p = H_J_p;

   rk->H_qq_p = H_qq_p;
   rk->H_qg_p = H_qg_p;
   rk->H_gq_p = H_gq_p;
   rk->H_gg_p = H_gg_p;

   rk->H_qq_m = H_qq_m;
   rk->H_qg_m = H_qg_m;
   rk->H_gq_m = H_gq_m;
   rk->H_gg_m = H_gg_m;
   if (to_initialize) {
      rk->k1 = init_vector(nn);
      rk->k2 = init_vector(nn);
      rk->k3 = init_vector(nn);
      rk->k4 = init_vector(nn);
      rk->temp = init_vector(nn);
   }

   if (rk->chiral_even_active && NULL == rk->k1_q) {
      rk->k1_q = init_vector(nn);
      rk->k2_q = init_vector(nn);
      rk->k3_q = init_vector(nn);
      rk->k4_q = init_vector(nn);
      rk->temp_q = init_vector(nn);

      rk->k1_g = init_vector(nn);
      rk->k2_g = init_vector(nn);
      rk->k3_g = init_vector(nn);
      rk->k4_g = init_vector(nn);
      rk->temp_g = init_vector(nn);
   }

   rk->f_1_p = f0_1_p;
   rk->f_1_m = f0_1_m;

   rk->f_2_p = f0_2_p;
   rk->f_2_m = f0_2_m;

   rk->f_3_p = f0_3_p;
   rk->f_3_m = f0_3_m;

   rk->f_4_p = f0_4_p;
   rk->f_4_m = f0_4_m;

   rk->f_q_p = f0_q_p;
   rk->f_q_m = f0_q_m;

   rk->f_g_p = f0_g_p;
   rk->f_g_m = f0_g_m;

   rk->h_up = h_up;
   rk->h_dn = h_dn;
   rk->h_st = h_st;

   rk->e_up = e_up;
   rk->e_dn = e_dn;
   rk->e_st = e_st;

   rk->as = as;
   rk->p_as = p_as;

   // Takes care automatically to rescale the kernels by nf if needed
   reset_scale_rk4th(rk, t0, dt);

   if ((int)round(*(rk->nf)) == 3) rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_3;
   if ((int)round(*(rk->nf)) == 4) rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_4;
   if ((int)round(*(rk->nf)) == 5) rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_5;

   if (!to_initialize) free(rk->which_CO);

   rk->which_CO = (vector_t **)calloc(6, sizeof(vector_t *));
   rk->how_many_CO = 0;
   if (rk->h_up != NULL) {
      rk->which_CO[rk->how_many_CO] = rk->h_up;
      rk->how_many_CO += 1;
   }
   if (rk->h_dn != NULL) {
      rk->which_CO[rk->how_many_CO] = rk->h_dn;
      rk->how_many_CO += 1;
   }
   if (rk->h_st != NULL) {
      rk->which_CO[rk->how_many_CO] = rk->h_st;
      rk->how_many_CO += 1;
   }

   if (rk->e_up != NULL) {
      rk->which_CO[rk->how_many_CO] = rk->e_up;
      rk->how_many_CO += 1;
   }
   if (rk->e_dn != NULL) {
      rk->which_CO[rk->how_many_CO] = rk->e_dn;
      rk->how_many_CO += 1;
   }
   if (rk->e_st != NULL) {
      rk->which_CO[rk->how_many_CO] = rk->e_st;
      rk->how_many_CO += 1;
   }

   if (rk->how_many_CO == 0) {
      free(rk->which_CO);
      rk->which_CO = NULL;
   } else rk->which_CO = (vector_t **)realloc(rk->which_CO, rk->how_many_CO * sizeof(vector_t *));

   return rk;
}

void free_rk4th(rk4th_internal_t **rk_dp)
{
   if (NULL == rk_dp) return;
   rk4th_internal_t *rk = *rk_dp;
   if (rk == NULL) return;
   free_vector(&(rk->k1));
   free_vector(&(rk->k2));
   free_vector(&(rk->k3));
   free_vector(&(rk->k4));
   free_vector(&(rk->temp));

   free_vector(&(rk->k1_q));
   free_vector(&(rk->k2_q));
   free_vector(&(rk->k3_q));
   free_vector(&(rk->k4_q));
   free_vector(&(rk->temp_q));

   free_vector(&(rk->k1_g));
   free_vector(&(rk->k2_g));
   free_vector(&(rk->k3_g));
   free_vector(&(rk->k4_g));
   free_vector(&(rk->temp_g));

   if (NULL != rk->which_CO) free(rk->which_CO);

   free(rk);

   *rk_dp = NULL;
}

static void rk4th_step_CE_nf_3(rk4th_internal_t *rk)
{
   rk4th_step_parallel(rk, rk->f_1_p, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_2_p, rk->H_NS);

   rk4th_step_parallel(rk, rk->f_1_m, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_2_m, rk->H_NS);

   rk4th_step_parallel_singlet(rk, rk->f_q_p, rk->f_g_p, rk->H_qq_p, rk->H_qg_p, rk->H_gq_p, rk->H_gg_p, rk->H_J_p);
   rk4th_step_parallel_singlet(rk, rk->f_q_m, rk->f_g_m, rk->H_qq_m, rk->H_qg_m, rk->H_gq_m, rk->H_gg_m, NULL);
}
static void rk4th_step_CE_nf_4(rk4th_internal_t *rk)
{
   rk4th_step_parallel(rk, rk->f_1_p, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_2_p, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_3_p, rk->H_NS);

   rk4th_step_parallel(rk, rk->f_1_m, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_2_m, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_3_m, rk->H_NS);

   rk4th_step_parallel_singlet(rk, rk->f_q_p, rk->f_g_p, rk->H_qq_p, rk->H_qg_p, rk->H_gq_p, rk->H_gg_p, rk->H_J_p);
   rk4th_step_parallel_singlet(rk, rk->f_q_m, rk->f_g_m, rk->H_qq_m, rk->H_qg_m, rk->H_gq_m, rk->H_gg_m, NULL);
}
static void rk4th_step_CE_nf_5(rk4th_internal_t *rk)
{
   rk4th_step_parallel(rk, rk->f_1_p, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_2_p, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_3_p, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_4_p, rk->H_NS);

   rk4th_step_parallel(rk, rk->f_1_m, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_2_m, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_3_m, rk->H_NS);
   rk4th_step_parallel(rk, rk->f_4_m, rk->H_NS);

   rk4th_step_parallel_singlet(rk, rk->f_q_p, rk->f_g_p, rk->H_qq_p, rk->H_qg_p, rk->H_gq_p, rk->H_gg_p, rk->H_J_p);
   rk4th_step_parallel_singlet(rk, rk->f_q_m, rk->f_g_m, rk->H_qq_m, rk->H_qg_m, rk->H_gq_m, rk->H_gg_m, NULL);
}

void rk4th_step(rk4th_internal_t *rk)
{
   for (int8_t i = 0; i < rk->how_many_CO; i++)
      rk4th_step_parallel(rk, rk->which_CO[i], rk->H_CO);

   if (!rk->chiral_even_active) {
      rk->t += rk->dt;
      return;
   }

   if (rk->dt > 0) {
      bool cond1 = (rk->t < rk->ch_threshold && (rk->t + rk->dt) >= rk->ch_threshold);
      bool cond2 = (rk->t < rk->bm_threshold && (rk->t + rk->dt) >= rk->bm_threshold);

      if (cond1 || cond2) {
         double stagger_scale = 0;
         if (cond1) stagger_scale = rk->ch_threshold;
         if (cond2) stagger_scale = rk->bm_threshold;
         double stagger = stagger_scale - rk->t;

         double dt_loc = rk->dt;
         // form t to threshold
         rk->dt = stagger;
         rk->dth = 0.5 * rk->dt;
         rk4th_step_CE_fixed_nf(rk);
         rk->t += rk->dt;

         int nf = (int)round(*(rk->nf));
         rescale_nf_kernels(rk, nf, nf + 1);
         (*(rk->nf))++;

         if (cond1) {
            copy_vector(rk->f_q_p, rk->f_3_p);
            copy_vector(rk->f_q_m, rk->f_3_m);
            rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_4;
         }
         if (cond2) {
            copy_vector(rk->f_q_p, rk->f_4_p);
            copy_vector(rk->f_q_m, rk->f_4_m);
            rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_5;
         }

         // from threshold to t+dt
         rk->dt = dt_loc - stagger;
         rk->dth = 0.5 * rk->dt;
         rk4th_step_CE_fixed_nf(rk);
         rk->t += rk->dt;

         rk->dt = dt_loc;
         rk->dth = 0.5 * rk->dt;
      } else {
         rk4th_step_CE_fixed_nf(rk);
         rk->t += rk->dt;
      }
   } else {
      bool cond1 = (rk->t > rk->bm_threshold && (rk->t + rk->dt) <= rk->bm_threshold);
      bool cond2 = (rk->t > rk->ch_threshold && (rk->t + rk->dt) <= rk->ch_threshold);

      if (cond1 || cond2) {
         double stagger_scale = 0;
         if (cond1) stagger_scale = rk->bm_threshold;
         if (cond2) stagger_scale = rk->ch_threshold;
         double stagger = stagger_scale - rk->t;

         double dt_loc = rk->dt;
         // form t to threshold
         rk->dt = stagger;
         rk->dth = 0.5 * rk->dt;
         rk4th_step_CE_fixed_nf(rk);
         rk->t += rk->dt;

         int nf = (int)round(*(rk->nf));
         rescale_nf_kernels(rk, nf, nf - 1);
         (*(rk->nf))--;

         if (cond1) {
            rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_4;

            for (int i = 0; i < rk->f_4_p->size; i++) {
               rk->f_4_p->data[i] = (rk->f_q_p->data[i] - rk->f_4_p->data[i]) / 5.0;
               rk->f_4_m->data[i] = (rk->f_q_m->data[i] - rk->f_4_m->data[i]) / 5.0;
            }
         }
         if (cond2) {
            rk4th_step_CE_fixed_nf = rk4th_step_CE_nf_3;

            for (int i = 0; i < rk->f_3_p->size; i++) {
               rk->f_3_p->data[i] = (rk->f_q_p->data[i] - rk->f_3_p->data[i]) / 4.0;
               rk->f_3_m->data[i] = (rk->f_q_m->data[i] - rk->f_3_m->data[i]) / 4.0;
            }
         }

         // from threshold to t+dt
         rk->dt = dt_loc - stagger;
         rk->dth = 0.5 * rk->dt;
         rk4th_step_CE_fixed_nf(rk);
         rk->t += rk->dt;

         rk->dt = dt_loc;
         rk->dth = 0.5 * rk->dt;
      } else {
         rk4th_step_CE_fixed_nf(rk);
         rk->t += rk->dt;
      }
   }
}

static void weighted_sum_four_vectors_loc(vector_t *v1, vector_t *v2, vector_t *v3, vector_t *v4, vector_t *v_out)
{
   for (int i = 0; i < v_out->size; i++) 
      v_out->data[i] += (v1->data[i] + v2->data[i] * 2.0 + v3->data[i] + v4->data[i]) / 3.0;
   
}

void debug_print(vector_t *v, int32_t N)
{
   FILE *fp = fopen("debug_print_vector.dat", "w");
   for (int32_t a = 0; a < v->size; a++) {
      int32_t i, j;
      from_a_to_ij(a, N, &i, &j);
      fprintf(fp, "%d\t%d\t%.6e\n", i, j, v->data[a]);
   }
   fclose(fp);

   exit(0);
}

static void rk4th_step_parallel(rk4th_internal_t *rk, vector_t *f, sparse_mat_t *Hk)
{
   //  f += dt/6 (k1 + 2k2 + 2k3 + k4)
   // equivalent to f += (k1T + 2k2T + k3T + k4T)/3
   // with k1T = (dt/2)k1, k2T = (dt/2)k2, k3T = (dt)k3, k4T = (dt/2)k4
   // and k1T = -as(t)(dt/2)H.f
   //     k2T = -as(t+dt/2)(dt/2)H.(f+k1T)
   //     k3T = -as(t+dt/2)(dt)H.(f+k2T)
   //     k4T = -as(t+dt)(dt/2)H.(f+k3T)
   //  k1T = -as(t)(dt/2)H.f

   sparse_mat_vector_mul_parallel(Hk, f, rk->k1);
   scalar_vector_mul(-(rk->dth) * rk->as(rk->t, rk->p_as), rk->k1, rk->k1);
  
   sparse_mat_2_vector_mul_parallel(Hk, rk->k1, f, rk->k2);
   scalar_vector_mul(-(rk->dth) * rk->as(rk->t + rk->dth, rk->p_as), rk->k2, rk->k2);

   sparse_mat_2_vector_mul_parallel(Hk, rk->k2, f, rk->k3);
   scalar_vector_mul(-(rk->dt) * rk->as(rk->t + rk->dth, rk->p_as), rk->k3, rk->k3);

   sparse_mat_2_vector_mul_parallel(Hk, rk->k3, f, rk->k4);
   scalar_vector_mul(-(rk->dth) * rk->as(rk->t + rk->dt, rk->p_as), rk->k4, rk->k4);

   weighted_sum_four_vectors_loc(rk->k1, rk->k2, rk->k3, rk->k4, f);

   return;
}

static void rk4th_step_parallel_singlet(rk4th_internal_t *rk, vector_t *f_q, vector_t *f_g, sparse_mat_t *H_qq, sparse_mat_t *H_qg, sparse_mat_t *H_gq, sparse_mat_t *H_gg, sparse_mat_t *H_J_p)
{
   
   sparse_mat_vector_mul_parallel(H_qq, f_q, rk->k1_q);
   if (NULL != H_J_p) {
      sparse_mat_vector_mul_parallel(H_J_p, f_q, rk->temp_q);
      sum_vector(rk->k1_q, rk->temp_q, rk->k1_q);
   }

   sparse_mat_vector_mul_parallel(H_qg, f_g, rk->temp_q);
   sum_vector(rk->k1_q, rk->temp_q, rk->k1_q);

   sparse_mat_vector_mul_parallel(H_gg, f_g, rk->k1_g);

   sparse_mat_vector_mul_parallel(H_gq, f_q, rk->temp_g);
   sum_vector(rk->k1_g, rk->temp_g, rk->k1_g);

   scalar_vector_mul(-(rk->dth) * rk->as(rk->t, rk->p_as), rk->k1_q, rk->k1_q);
   scalar_vector_mul(-(rk->dth) * rk->as(rk->t, rk->p_as), rk->k1_g, rk->k1_g);


   sparse_mat_2_vector_mul_parallel(H_qq, rk->k1_q, f_q, rk->k2_q);
   if (NULL != H_J_p) {
      sparse_mat_2_vector_mul_parallel(H_J_p, rk->k1_q, f_q, rk->temp_q);
      sum_vector(rk->k2_q, rk->temp_q, rk->k2_q);
   }

   sparse_mat_2_vector_mul_parallel(H_qg, rk->k1_g, f_g, rk->temp_q);
   sum_vector(rk->k2_q, rk->temp_q, rk->k2_q);

   sparse_mat_2_vector_mul_parallel(H_gg, rk->k1_g, f_g, rk->k2_g);
   sparse_mat_2_vector_mul_parallel(H_gq, rk->k1_q, f_q, rk->temp_g);
   sum_vector(rk->k2_g, rk->temp_g, rk->k2_g);

   scalar_vector_mul(-(rk->dth) * rk->as(rk->t + rk->dth, rk->p_as), rk->k2_q, rk->k2_q);
   scalar_vector_mul(-(rk->dth) * rk->as(rk->t + rk->dth, rk->p_as), rk->k2_g, rk->k2_g);

   sparse_mat_2_vector_mul_parallel(H_qq, rk->k2_q, f_q, rk->k3_q);
   if (NULL != H_J_p) {
      sparse_mat_2_vector_mul_parallel(H_J_p, rk->k2_q, f_q, rk->temp_q);
      sum_vector(rk->k3_q, rk->temp_q, rk->k3_q);
   }

   sparse_mat_2_vector_mul_parallel(H_qg, rk->k2_g, f_g, rk->temp_q);
   sum_vector(rk->k3_q, rk->temp_q, rk->k3_q);

   sparse_mat_2_vector_mul_parallel(H_gg, rk->k2_g, f_g, rk->k3_g);
   sparse_mat_2_vector_mul_parallel(H_gq, rk->k2_q, f_q, rk->temp_g);
   sum_vector(rk->k3_g, rk->temp_g, rk->k3_g);

   scalar_vector_mul(-(rk->dt) * rk->as(rk->t + rk->dth, rk->p_as), rk->k3_q, rk->k3_q);
   scalar_vector_mul(-(rk->dt) * rk->as(rk->t + rk->dth, rk->p_as), rk->k3_g, rk->k3_g);

   sparse_mat_2_vector_mul_parallel(H_qq, rk->k3_q, f_q, rk->k4_q);
   if (NULL != H_J_p) {
      sparse_mat_2_vector_mul_parallel(H_J_p, rk->k3_q, f_q, rk->temp_q);
      sum_vector(rk->k4_q, rk->temp_q, rk->k4_q);
   }

   sparse_mat_2_vector_mul_parallel(H_qg, rk->k3_g, f_g, rk->temp_q);
   sum_vector(rk->k4_q, rk->temp_q, rk->k4_q);

   sparse_mat_2_vector_mul_parallel(H_gg, rk->k3_g, f_g, rk->k4_g);
   sparse_mat_2_vector_mul_parallel(H_gq, rk->k3_q, f_q, rk->temp_g);
   sum_vector(rk->k4_g, rk->temp_g, rk->k4_g);

   scalar_vector_mul(-(rk->dth) * rk->as(rk->t + rk->dt, rk->p_as), rk->k4_q, rk->k4_q);
   scalar_vector_mul(-(rk->dth) * rk->as(rk->t + rk->dt, rk->p_as), rk->k4_g, rk->k4_g);

   weighted_sum_four_vectors_loc(rk->k1_q, rk->k2_q, rk->k3_q, rk->k4_q, f_q);
   weighted_sum_four_vectors_loc(rk->k1_g, rk->k2_g, rk->k3_g, rk->k4_g, f_g);

   return;
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

