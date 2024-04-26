//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/kernels.h>
#include <TW3EV_include/kernels_common.h>

static const double Nc = 3.0;
static const double Ncm1 = 0.33333333333333;
static const double CF = 1.33333333333333;
static const double CA = 3.0;

stored_point_t **stored_points;
thpool_t *thp_kernels_h;
static int32_t n_threads_global;

static volatile int32_t done_tasks = 0;
static pthread_mutex_t lock_task = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t notify_task = PTHREAD_COND_INITIALIZER;

#define SIGNAL_DONE_TASK()                                                                                                                                                                                                                               \
   do {                                                                                                                                                                                                                                                  \
      pthread_mutex_lock(&lock_task);                                                                                                                                                                                                                    \
      done_tasks++;                                                                                                                                                                                                                                      \
      pthread_cond_signal(&notify_task);                                                                                                                                                                                                                 \
      pthread_mutex_unlock(&lock_task);                                                                                                                                                                                                                  \
   } while (0);

#define WAIT_ON_DONE_TASK(M)                                                                                                                                                                                                                             \
   {                                                                                                                                                                                                                                                     \
      do {                                                                                                                                                                                                                                               \
         pthread_mutex_lock(&lock_task);                                                                                                                                                                                                                 \
         while (done_tasks < (M))                                                                                                                                                                                                                        \
            pthread_cond_wait(&notify_task, &lock_task);                                                                                                                                                                                                 \
         pthread_mutex_unlock(&lock_task);                                                                                                                                                                                                               \
         break;                                                                                                                                                                                                                                          \
      } while (1);                                                                                                                                                                                                                                       \
      done_tasks = 0;                                                                                                                                                                                                                                    \
   }

static nonlinear_radial_grid_type_e G_type = LOG_GRID;
static angular_grid_type_e AG_type = LIN_GRID;
static double grid_exponent = 2;

void get_xmin_xmax2(int32_t i, int32_t j, stored_point_t *sp, int32_t N, int32_t M, double c_fact);

double (*Fijk)(double ix, double jy, double i, double j, int32_t N);
double Fijk_plus(double ix, double jy, double i, double j, int32_t N);
double Fijk_minus(double ix, double jy, double i, double j, int32_t N);
double Fijk_both(double ix, double jy, double i, double j, int32_t N);

static void get_grid_in_3D_space(int32_t i, int32_t j, int32_t N, int32_t M, double c_fact, double *x1, double *x2, double *x3);
void get_3D_in_grid_space_specific(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho);

void (*which_get_3D_in_grid_space_angle)(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho);
void (*which_get_3D_in_grid_space_radial)(int32_t M, double c_fact, double *rho);

void get_3D_in_grid_space_LIN_angle(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho);
void get_3D_in_grid_space_COS_angle(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho);
void get_3D_in_grid_space_LOG_radial(int32_t M, double c_fact, double *rho);
void get_3D_in_grid_space_PWR_radial(int32_t M, double c_fact, double *rho);
void get_3D_in_grid_space_IML_radial(int32_t M, double c_fact, double *rho);
void get_3D_in_grid_space_HYP_radial(int32_t M, double c_fact, double *rho);

void kernels_set_grid_type(interpolant_type_e F_t, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t, double pge)
{
   (void)F_t;
   G_type = G_t;
   AG_type = AG_t;
   grid_exponent = pge;
}

void get_grid_in_3D_space_double(double i, double j, int32_t N, int32_t M, double c_fact, double *x1, double *x2, double *x3);

void kernels_setup_for_computation(int32_t N, int32_t M, double c_fact, int32_t n_threads, interpolant_type_e F_t, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t, double pge, const char basefolder[])
{
   kernels_set_grid_type(F_t, G_t, AG_t, pge);

   if (n_threads <= 0) n_threads = 1;
   n_threads_global = n_threads;
   thp_kernels_h = thpool_create(n_threads, n_threads + 1);
   done_tasks = 0;
   FILE *fp = NULL;
   {
      char folder[BASEFOLDER_MAX_L] = "";
      strncpy(folder, basefolder, BASEFOLDER_MAX_L - 1);
      char grid_type_string[255] = "";
      switch (G_t) {
      case LOG_GRID: {
         which_get_3D_in_grid_space_radial = get_3D_in_grid_space_LOG_radial;
         strcat(grid_type_string, "_LOG");
         break;
      }
      case PWR_GRID: {
         which_get_3D_in_grid_space_radial = get_3D_in_grid_space_PWR_radial;
         strcat(grid_type_string, "_PWR");
         break;
      }
      case IML_GRID: {
         which_get_3D_in_grid_space_radial = get_3D_in_grid_space_IML_radial;
         strcat(grid_type_string, "_IML");
         break;
      }
      case HYP_GRID: {
         which_get_3D_in_grid_space_radial = get_3D_in_grid_space_HYP_radial;
         strcat(grid_type_string, "_HYP");
         break;
      }
      default: {
         tw3ev_log(TW3EV_ERROR, "Unsupported non-linear grid type. Supported "
                                "types: LOG_GRID, PWR_GRID, IML_GRID, HYP_GRID");
         break;
      }
      }

      switch (AG_t) {
      case LIN_GRID: {
         which_get_3D_in_grid_space_angle = get_3D_in_grid_space_LIN_angle;
         strcat(grid_type_string, "_LIN");
         break;
      }
      case COS_GRID: {
         which_get_3D_in_grid_space_angle = get_3D_in_grid_space_COS_angle;
         strcat(grid_type_string, "_COS");
         break;
      }
      default:
         break;
      }

      switch (F_t) {
      case IT_BOTH: {
         Fijk = Fijk_both;
         strcat(grid_type_string, "_BOTH");
         break;
      }
      case IT_PLUS: {
         Fijk = Fijk_plus;
         strcat(grid_type_string, "_PLUS");
         break;
      }
      case IT_MINUS: {
         Fijk = Fijk_minus;
         strcat(grid_type_string, "_MINUS");
         break;
      }
      default:
         break;
      }
      char *fName = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
      char N_as_string[255];
      char extension[6] = ".dat";
      sprintf(N_as_string, "_N%d_M%d", N, M);
      fName = clear_string(fName);
      strncpy(fName, "grid", BASEFOLDER_MAX_L);
      strcat(fName, grid_type_string);
      strcat(fName, N_as_string);
      strcat(fName, extension);
      fName = attach_path_to_string(fName, folder, "dat/grids/");
      fp = fopen(fName, "w");
      free(fName);
   }

   if (NULL != stored_points) kernels_free_setup(N, M);

   stored_points = (stored_point_t **)calloc((N + 1), sizeof(stored_point_t *));
   for (int32_t i = 0; i <= N; i++) {
      stored_points[i] = (stored_point_t *)calloc((M + 1), sizeof(stored_point_t));
      for (int32_t j = 0; j <= M; j++) {
         get_grid_in_3D_space(i, j, N, M, c_fact, &(stored_points[i][j].x1), &(stored_points[i][j].x2), &(stored_points[i][j].x3));

         get_xmin_xmax2(i, j, &(stored_points[i][j]), N, M, c_fact);
         if (NULL != fp) fprintf(fp, "%d\t%d\t%le\t%le\t%le\n", i, j, (stored_points[i][j].x1), (stored_points[i][j].x2), (stored_points[i][j].x3));
      }
   }
   if (NULL != fp) fclose(fp);
}

#define _UPDATE_XMAX_XMIN_()                                                                                                                                                                                                                             \
   do {                                                                                                                                                                                                                                                  \
      sp->x1min = min2(x1, sp->x1min);                                                                                                                                                                                                                   \
      sp->x2min = min2(x2, sp->x2min);                                                                                                                                                                                                                   \
      sp->x3min = min2(x3, sp->x3min);                                                                                                                                                                                                                   \
      sp->x1max = max2(x1, sp->x1max);                                                                                                                                                                                                                   \
      sp->x2max = max2(x2, sp->x2max);                                                                                                                                                                                                                   \
      sp->x3max = max2(x3, sp->x3max);                                                                                                                                                                                                                   \
   } while (0);

void get_xmin_xmax2(int32_t i, int32_t j, stored_point_t *sp, int32_t N, int32_t M, double c_fact)
{

   double x1 = 0, x2 = 0, x3 = 0;
   get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3);
   sp->x1min = x1;
   sp->x2min = x2;
   sp->x3min = x3;

   sp->x1max = x1;
   sp->x2max = x2;
   sp->x3max = x3;

   for (int li = -1; li <= 1; li++) {
      for (int lj = -1; lj <= 1; lj++) {
         if (j + lj < 0 || j + lj > M) continue;
         get_grid_in_3D_space(i + li, j + lj, N, M, c_fact, &x1, &x2, &x3);
         sp->x1min = min2(x1, sp->x1min);
         sp->x2min = min2(x2, sp->x2min);
         sp->x3min = min2(x3, sp->x3min);

         sp->x1max = max2(x1, sp->x1max);
         sp->x2max = max2(x2, sp->x2max);
         sp->x3max = max2(x3, sp->x3max);
      }
   }
}

void kernels_free_setup(int32_t N, int32_t M)
{
   (void)M;
   for (int32_t i = 0; i <= N; i++) {
      free(stored_points[i]);
   }
   free(stored_points);
   thpool_destroy(thp_kernels_h);
   stored_points = NULL;
}

double H_NS(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   (void)nf;
   double an_dim = (i == ip && j == jp) ? (-3 * CF) : 0.0;
   double Ncpart = Hhat12(i, j, ip, jp, IW) + Hhat23(i, j, ip, jp, IW) - 2.0 * Hplus12(i, j, ip, jp, IW);
   double Ncm1part = Hhat13(i, j, ip, jp, IW) - Hplus13(i, j, ip, jp, IW) - He23P23(i, j, ip, jp, IW) + 2 * Hminus12(i, j, ip, jp, IW);
   return Nc * Ncpart - Ncm1 * Ncm1part + an_dim;
}

double H_CO(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   (void)nf;
   double an_dim = (i == ip && j == jp) ? (-3 * CF) : 0.0;
   double Ncpart = Hhat12(i, j, ip, jp, IW) + Hhat23(i, j, ip, jp, IW) - 2.0 * Hplus12(i, j, ip, jp, IW) - 2.0 * Hplus23(i, j, ip, jp, IW);
   double Ncm1part = Hhat13(i, j, ip, jp, IW) + 2 * Hminus12(i, j, ip, jp, IW) + 2 * Hminus23(i, j, ip, jp, IW);
   return Nc * Ncpart - Ncm1 * Ncm1part + an_dim;
}

double H_J_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf) { return 4 * nf * Hd13(i, j, ip, jp, IW); }

double H_qq_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   (void)i;
   (void)j;
   (void)ip;
   (void)jp;
   (void)IW;
   (void)nf;

   tw3ev_log(TW3EV_ERROR, "The H_qq_p function is for debug-purpose only, it "
                          "should not be invoked directly.");
   return 0;
}

double H_qg_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf) { return nf * (Vp13(i, j, ip, jp, IW) - Vm13(i, j, ip, jp, IW)); }

double H_gq_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   (void)nf;
   double temp1 = Wp13(i, j, ip, jp, IW) + Wm13(i, j, ip, jp, IW) - 2.0 * DW13(i, j, ip, jp, IW);
   double temp2 = Wp13P23(i, j, ip, jp, IW) + Wm13P23(i, j, ip, jp, IW) - 2.0 * DW13P23(i, j, ip, jp, IW);
   return CA * (temp1 - temp2);
}

double H_gg_p(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   double temp1 = Hhat12GG(i, j, ip, jp, IW) + Hhat23GG(i, j, ip, jp, IW) + Hhat31GG(i, j, ip, jp, IW);
   double temp2 = -4.0 * (Hplus12GG(i, j, ip, jp, IW) + Hplus13GG(i, j, ip, jp, IW)) - 2.0 * (Htildeplus12GG(i, j, ip, jp, IW) + Htildeplus13GG(i, j, ip, jp, IW));
   double temp3 = 6.0 * (Hminus12GG(i, j, ip, jp, IW) + Hminus13GG(i, j, ip, jp, IW));
   double b0 = ((i == ip && j == jp) ? 11.0 - 2.0 * nf / 3.0 : 0);
   return CA * (temp1 + temp2 + temp3) - b0;
}

double H_qq_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   (void)i;
   (void)j;
   (void)ip;
   (void)jp;
   (void)IW;
   (void)nf;
   tw3ev_log(TW3EV_ERROR, "The H_qq_m function is for debug-purpose only, it "
                          "should not be invoked directly.");
   return 0;
}

double H_qg_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf) { return nf * (Vp13(i, j, ip, jp, IW) + Vm13(i, j, ip, jp, IW)); }

double H_gq_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   (void)nf;
   double temp1 = Wp13(i, j, ip, jp, IW) + Wm13(i, j, ip, jp, IW);
   double temp2 = Wp13P23(i, j, ip, jp, IW) + Wm13P23(i, j, ip, jp, IW);
   return -1.0 * ((Nc * Nc - 4.0) / Nc) * (temp1 + temp2);
}

double H_gg_m(int32_t i, int32_t j, int32_t ip, int32_t jp, integration_wrapper_t *IW, int nf)
{
   double temp1 = Hhat12GG(i, j, ip, jp, IW) + Hhat23GG(i, j, ip, jp, IW) + Hhat31GG(i, j, ip, jp, IW);
   double temp2 = -4.0 * (Hplus12GG(i, j, ip, jp, IW) + Hplus13GG(i, j, ip, jp, IW)) - 2.0 * (Htildeplus12GG(i, j, ip, jp, IW) + Htildeplus13GG(i, j, ip, jp, IW));
   double temp3 = -6.0 * (Hminus12GG(i, j, ip, jp, IW) + Hminus13GG(i, j, ip, jp, IW));
   double b0 = ((i == ip && j == jp) ? 11.0 - 2.0 * nf / 3.0 : 0);
   return CA * (temp1 + temp2 + temp3) - b0;
}

sparse_mat_t *init_kernel_single(double (*Hk)(int32_t, int32_t, int32_t, int32_t, integration_wrapper_t *, int), int32_t N, int32_t M, double c_fact, int32_t n_low, int32_t n_up, int nf)
{

   integration_wrapper_t *IW = (integration_wrapper_t *)calloc(1, sizeof(integration_wrapper_t));
   IW->ws_integr = NULL; 

   IW->int_par = (integration_par_t *)calloc(1, sizeof(integration_par_t));
   IW->int_par->N = N;
   IW->int_par->M = M;
   IW->int_par->c_fact = c_fact;

   int32_t NN = (N + 1) * M * M * (n_up - n_low);
   int32_t curr_size = 0;
   sparse_mat_t *sp_m = init_sparse_matrix(NN, NN, NN); 
   for (int32_t i = n_low; i < n_up; i++) {
      for (int32_t j = 0; j < M; j++) {
         for (int32_t ip = 0; ip <= N; ip++) {
            for (int32_t jp = j; jp < M; jp++) {
               double temp = Hk(i, j, ip, jp, IW, nf);
            
               if (fabs(temp) > 30 * _ZERO_THR_) {
                  sp_m->data[curr_size] = temp;
                  sp_m->i_r[curr_size] = from_ij_to_a(i, j, M);
                  sp_m->i_c[curr_size] = from_ij_to_a(ip, jp, M);
                  curr_size++;
                 
               }
            }
         }
      }
   }

   resize_sparse_matrix(sp_m, curr_size);
   free(IW->int_par);
   free(IW);

   return sp_m;
}

typedef struct {
   double (*Hk)(int32_t, int32_t, int32_t, int32_t, integration_wrapper_t *, int);
   int nf;
   int32_t N, M, n_low, n_up;
   sparse_mat_t *sp_m;
   double c_fact;

   sparse_mat_t *global_sp_m;
   int32_t offset;
} thread_container_t;

void init_kernel_thread_function(void *p)
{
   thread_container_t *tc = (thread_container_t *)p;
   tc->sp_m = init_kernel_single(tc->Hk, tc->N, tc->M, tc->c_fact, tc->n_low, tc->n_up, tc->nf);
   p = (void *)tc;
   SIGNAL_DONE_TASK();
}

static void local_copy_sp_memory(void *p)
{
   thread_container_t *tc = (thread_container_t *)p;
   memcpy(tc->global_sp_m->data + tc->offset, tc->sp_m->data, sizeof(double) * tc->sp_m->size);
   memcpy(tc->global_sp_m->i_c + tc->offset, tc->sp_m->i_c, sizeof(int32_t) * tc->sp_m->size);
   memcpy(tc->global_sp_m->i_r + tc->offset, tc->sp_m->i_r, sizeof(int32_t) * tc->sp_m->size);
   SIGNAL_DONE_TASK();
}

sparse_mat_t *init_kernel(double (*Hk)(int32_t, int32_t, int32_t, int32_t, integration_wrapper_t *, int), int32_t N, int32_t M, double c_fact, int nf, printout_level_e pl)
{
   struct timespec start, finish;
   double elapsed;
   if (pl == PL_ALL) clock_gettime(CLOCK_MONOTONIC, &start);

   int32_t n_thread = n_threads_global;

   int32_t n_min = 0;
   int32_t n_max = N;
   int32_t dn = (n_max - n_min) / (n_thread);

   thread_container_t tcs[n_thread];
   if (pl == PL_ALL) tw3ev_log(TW3EV_INFO, "# Subdivision of work for kernel initialization: ");
   for (int32_t i = 0; i < n_thread; i++) {
      (tcs[i]).nf = nf;
      (tcs[i]).Hk = Hk;
      (tcs[i]).N = N;
      (tcs[i]).M = M;
      (tcs[i]).c_fact = c_fact;
      (tcs[i]).n_low = n_min + i * dn;
      (tcs[i]).n_up = (i != n_thread - 1 ? n_min + (i + 1) * dn : n_max + 1);
      if (pl == PL_ALL) tw3ev_log(TW3EV_INFO, "# %d %d", (tcs[i]).n_low, (tcs[i]).n_up);
      thpool_add(thp_kernels_h, init_kernel_thread_function, (void *)&(tcs[i]));
   }
   WAIT_ON_DONE_TASK(n_thread);

   tcs[0].offset = 0;
   for (int32_t i = 1; i < n_thread; i++) {
      tcs[i].offset = tcs[i - 1].offset + tcs[i - 1].sp_m->size;
   }
   int32_t size_tot = tcs[n_thread - 1].offset + tcs[n_thread - 1].sp_m->size;
   int32_t NN = (N + 1) * (M + 1);
   sparse_mat_t *sp_m = init_sparse_matrix(size_tot, NN, NN);

   for (int32_t i = 0; i < n_thread; i++) {
      tcs[i].global_sp_m = sp_m;
      thpool_add(thp_kernels_h, local_copy_sp_memory, (void *)&(tcs[i]));
   }

   WAIT_ON_DONE_TASK(n_thread);

   for (int32_t i = 0; i < n_thread; i++)
      free_sparse_mat(&(tcs[i].sp_m));

   if (pl == PL_ALL) {
      clock_gettime(CLOCK_MONOTONIC, &finish);

      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      tw3ev_log(TW3EV_INFO, "Execution time for kernel initialization with N=%d M=%d: %.4e", N, M, elapsed);
   }

   generate_intervals_sparse_matrix(sp_m);

   return sp_m;
}


double Fijk_plus(double ix, double jy, double i, double j, int32_t N)
{
   (void)N;
   return max2(1.0 - max3(fabs(ix - i), fabs(jy - j), fabs(-ix - jy + i + j)), 0);
}

double Fijk_minus(double ix, double jy, double i, double j, int32_t N)
{
   (void)N;
   return max2(1.0 - max3(fabs(ix - i), fabs(jy - j), fabs(ix - jy - i + j)), 0);
}

double Fijk_both(double ix, double jy, double i, double j, int32_t N)
{
   (void)N;
   return (max2(1.0 - max3(fabs(ix - i), fabs(jy - j), fabs(-ix - jy + i + j)), 0) + max2(1.0 - max3(fabs(ix - i), fabs(jy - j), fabs(+ix - jy - i + j)), 0)) * 0.5;
}

void get_grid_in_3D_space_double(double i, double j, int32_t N, int32_t M, double c_fact, double *x1, double *x2, double *x3)
{
   if (i >= N + 1 || i < 0) i -= sign(i) * (N + 1);
   int32_t nT = (N + 1) / 6;
   int8_t sector = (int8_t)(i / ((double)nT));

   double rho = 1.0, phi = 0.0;
   switch (G_type) {
   case LOG_GRID: {
      rho = exp(c_fact * (((double)(j - M)) / ((double)M)));
      break;
   }
   case PWR_GRID: {
      rho = pow((((double)j / ((double)M)) + c_fact) / (1 + c_fact), grid_exponent);
      break;
   }
   case IML_GRID: {
      rho = 2.0 * exp(((double)(j - M)) / (c_fact * M)) / (1.0 + exp(((double)(j - M)) / (c_fact * M)));
      break;
   }
   case HYP_GRID: {
      rho = pow(cosh(((double)(j - M)) / ((double)M * c_fact)), -grid_exponent);
      break;
   }
   default:
      break;
   }
   switch (AG_type) {
   case LIN_GRID: {
      phi = i / ((double)nT) - (double)sector;
      break;
   }
   case COS_GRID: {
      phi = 0.5 - 0.5 * cos(M_PI * (i / ((double)nT) - (double)sector));
      break;
   }
   default:
      break;
   }

   switch (sector) {
   case 0: {
      *x1 = rho * (1 - phi);
      *x2 = rho * (phi);
      break;
   }
   case 1: {
      *x1 = rho * (-phi);
      *x2 = rho * (1);
      break;
   }
   case 2: {
      *x1 = rho * (-1);
      *x2 = rho * (1 - phi);
      break;
   }
   case 3: {
      *x1 = rho * (-1 + phi);
      *x2 = rho * (-phi);
      break;
   }
   case 4: {
      *x1 = rho * (phi);
      *x2 = rho * (-1);
      break;
   }
   case 5: {
      *x1 = rho * (1);
      *x2 = rho * (-1 + phi);
      break;
   }
   default:
      break;
   }
   *x3 = -(*x1) - (*x2);
}

static void get_grid_in_3D_space(int32_t i, int32_t j, int32_t N, int32_t M, double c_fact, double *x1, double *x2, double *x3)
{
   if (i >= N + 1 || i < 0) i -= isign(i) * (N + 1);

   int32_t nT = (N + 1) / 6;
   int8_t sector = i / nT;
   double rho = 1.0, phi = 0.0;
   switch (G_type) {
   case LOG_GRID: {
      rho = exp(c_fact * (((double)(j - M)) / ((double)M)));
      break;
   }
   case PWR_GRID: {
      rho = pow((((double)j / ((double)M)) + c_fact) / (1 + c_fact), grid_exponent);
      break;
   }
   case IML_GRID: {
      rho = 2.0 * exp(((double)(j - M)) / (c_fact * M)) / (1.0 + exp(((double)(j - M)) / (c_fact * M)));
      break;
   }
   case HYP_GRID: {
      rho = pow(cosh(((double)(j - M)) / ((double)M * c_fact)), -grid_exponent);
      break;
   }
   default:
      break;
   }
   switch (AG_type) {
   case LIN_GRID: {
      phi = ((double)(i % nT)) / ((double)nT);
      break;
   }
   case COS_GRID: {
      phi = 0.5 - 0.5 * cos(M_PI * ((double)(i % nT)) / ((double)nT));
      break;
   }
   default:
      break;
   }

   switch (sector) {
   case 0: {
      *x1 = rho * (1 - phi);
      *x2 = rho * (phi);
      break;
   }
   case 1: {
      *x1 = rho * (-phi);
      *x2 = rho * (1);
      break;
   }
   case 2: {
      *x1 = rho * (-1);
      *x2 = rho * (1 - phi);
      break;
   }
   case 3: {
      *x1 = rho * (-1 + phi);
      *x2 = rho * (-phi);
      break;
   }
   case 4: {
      *x1 = rho * (phi);
      *x2 = rho * (-1);
      break;
   }
   case 5: {
      *x1 = rho * (1);
      *x2 = rho * (-1 + phi);
      break;
   }
   default:
      break;
   }
   *x3 = -(*x1) - (*x2);
}


void get_3D_in_grid_space_specific(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho)
{
   *rho = max3(fabs(x1), fabs(x2), fabs(x3));
   which_get_3D_in_grid_space_angle(x1, x2, x3, N, M, c_fact, phi, rho);
   which_get_3D_in_grid_space_radial(M, c_fact, rho);
}

double Fij(double x1, double x2, double x3, int32_t i, int32_t j, int32_t N, int32_t M, double c_fact)
{
   double rho, phi;
   get_3D_in_grid_space_specific(x1, x2, x3, N, M, c_fact, &phi, &rho);
   if (i == 0 && phi > N && phi <= N + 1) phi -= N + 1;
   return Fijk(phi, rho, i, j, N);
}

void get_3D_in_grid_space_LIN_angle(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho)
{
   (void)M;
   (void)c_fact;
   int nT = (N + 1) / 6;
   if (x1 > 0 && x2 >= 0 && x3 < 0) *phi = (x2 / (*rho)) * nT;
   if (x1 <= 0 && x2 > 0 && x3 < 0) *phi = (1 - x1 / (*rho)) * nT;
   if (x1 < 0 && x2 > 0 && x3 >= 0) *phi = (3 - x2 / (*rho)) * nT;
   if (x1 < 0 && x2 <= 0 && x3 > 0) *phi = (3 - x2 / (*rho)) * nT;
   if (x1 >= 0 && x2 < 0 && x3 > 0) *phi = (4 + x1 / (*rho)) * nT;
   if (x1 > 0 && x2 < 0 && x3 <= 0) *phi = (6 + x2 / (*rho)) * nT;
}
void get_3D_in_grid_space_COS_angle(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho)
{
   (void)M;
   (void)c_fact;

   int nT = (N + 1) / 6;
   if (x1 > 0 && x2 >= 0 && x3 < 0) *phi = nT * ONE_O_PI * acos(1.0 - 2.0 * (x2 / (*rho)));
   else if (x1 <= 0 && x2 > 0 && x3 < 0) *phi = nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (-x1 / (*rho)));
   else if (x1 < 0 && x2 > 0 && x3 >= 0) *phi = 2 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (1 - x2 / (*rho)));
   else if (x1 < 0 && x2 <= 0 && x3 > 0) *phi = 3 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (-x2 / (*rho)));
   else if (x1 >= 0 && x2 < 0 && x3 > 0) *phi = 4 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (x1 / (*rho)));
   else if (x1 > 0 && x2 < 0 && x3 <= 0) *phi = 5 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (1 + x2 / (*rho)));
   else *phi = 0;
}

void get_3D_in_grid_space_LOG_radial(int32_t M, double c_fact, double *rho) { *rho = log((*rho)) * M / c_fact + M; }
void get_3D_in_grid_space_PWR_radial(int32_t M, double c_fact, double *rho) { *rho = (pow((*rho), 1.0 / grid_exponent) * (1 + c_fact) - c_fact) * M; }
void get_3D_in_grid_space_IML_radial(int32_t M, double c_fact, double *rho) { *rho = M * (1.0 + c_fact * log((*rho) / (2.0 - (*rho)))); }
void get_3D_in_grid_space_HYP_radial(int32_t M, double c_fact, double *rho) { *rho = M * (1.0 + c_fact * acosh(1.0 / pow((*rho), 1.0 / grid_exponent))); }

double execution_time(void (*fnc)(void *), void *p)
{
   clock_t start, end;

   start = clock();
   fnc(p);
   end = clock();
   return ((double)(end - start)) / CLOCKS_PER_SEC;
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
