//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/evolution_interface.h>
#include<float.h>

#define _QUARK_COMBINATION_TYPE_ 1
#define _GLUON_COMBINATION_TYPE_ 2

void get_3D_in_grid_space(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho, double grid_exp, nonlinear_radial_grid_type_e G_type, angular_grid_type_e AG_type);

static input_parameters_t _default_in_par = (input_parameters_t){
    .N = 20,
    .M = 25,
    .mu0_2 = 1,
    .muF_2 = 10,
    .nstep = 100,
    .n_thread = 1,
    .ch_thr_mu2 = 1.6129,
    .bm_thr_mu2 = 17.4724,
    .xmin = 0.01,
    .F_t = IT_BOTH,
    .G_t = HYP_GRID,
    .AG_t = LIN_GRID,
    .grid_exponent = 3,
    .chiral_even_active = true,
    .chiral_odd_active = true,
    .i_rule = RULE_TW3EV_GK21,
    .basefolder = NULL,
    .ker_subfolder = NULL,
    .res_subfolder = NULL,
    .model = NULL,
};

void supply_default_parameters(input_parameters_t *in_par)
{
   char *tmp_basefolder = in_par->basefolder;
   char *tmp_ker_subfolder = in_par->ker_subfolder;
   char *tmp_res_subfolder = in_par->res_subfolder;
   char *tmp_model = in_par->model;

   memcpy(in_par, &_default_in_par, sizeof(input_parameters_t));
   if (tmp_basefolder != NULL) in_par->basefolder = tmp_basefolder;
   else in_par->basefolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));

   strcpy(in_par->basefolder, "./");

   if (tmp_ker_subfolder != NULL) in_par->ker_subfolder = tmp_ker_subfolder;
   else in_par->ker_subfolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));

   strcpy(in_par->ker_subfolder, "");

   if (tmp_res_subfolder != NULL) in_par->res_subfolder = tmp_res_subfolder;
   else in_par->res_subfolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));

   strcpy(in_par->res_subfolder, "");

   if (tmp_basefolder != NULL) in_par->model = tmp_model;
   else {
      in_par->model = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   }
   strcpy(in_par->model, "");
}

static inline int local_parse_bool_from_string(char *buff)
{
   if (strcmp(buff, "false") == 0 || strcmp(buff, "False") == 0 || strcmp(buff, "f") == 0 || strcmp(buff, "F") == 0) return 0;
   if (strcmp(buff, "true") == 0 || strcmp(buff, "True") == 0 || strcmp(buff, "t") == 0 || strcmp(buff, "T") == 0) return 1;
   return -1;
}

static inline void tw3ev_write_configuration_file(input_parameters_t *in_par, const char *fname)
{
   FILE *fp = fopen(fname, "w");
   if (!fp) {
      tw3ev_log(TW3EV_WARNING, "Ubable to open file to dump configuration.");
      return;
   }
   fprintf(fp, "n\t%d\n", in_par->N);
   fprintf(fp, "M\t%d\n", in_par->M);
   fprintf(fp, "mu0_2\t%le\n", in_par->mu0_2);
   fprintf(fp, "mu0_2\t%le\n", in_par->muF_2);
   fprintf(fp, "nstep\t%d\n", in_par->nstep);
   fprintf(fp, "n_thread\t%d\n", in_par->n_thread);
   fprintf(fp, "charm_THR\t%le\n", in_par->ch_thr_mu2);
   fprintf(fp, "bottom_THR\t%le\n", in_par->bm_thr_mu2);
   fprintf(fp, "xmin\t%le\n", in_par->xmin);
   fprintf(fp, "grid_exponent\t%le\n", in_par->grid_exponent);

   switch (in_par->F_t) {
   case IT_MINUS:
      fprintf(fp, "interpolant_type\t%s\n", "IT_MINUS");
      break;
   case IT_PLUS:
      fprintf(fp, "interpolant_type\t%s\n", "IT_PLUS");
      break;
   case IT_BOTH:
      fprintf(fp, "interpolant_type\t%s\n", "IT_BOTH");
      break;
   default:
      break;
   }

   switch (in_par->G_t) {
   case LOG_GRID:
      fprintf(fp, "radial_grid_type\t%s\n", "LOG_GRID");
      break;
   case PWR_GRID:
      fprintf(fp, "radial_grid_type\t%s\n", "PWR_GRID");
      break;
   case IML_GRID:
      fprintf(fp, "radial_grid_type\t%s\n", "IML_GRID");
      break;
   case HYP_GRID:
      fprintf(fp, "radial_grid_type\t%s\n", "HYP_GRID");
      break;
   default:
      break;
   }

   switch (in_par->AG_t) {
   case LIN_GRID:
      fprintf(fp, "interpolant_type\t%s\n", "LIN_GRID");
      break;
   case COS_GRID:
      fprintf(fp, "interpolant_type\t%s\n", "COS_GRID");
      break;
   default:
      break;
   }

   switch (in_par->i_rule) {
   case RULE_TW3EV_GK21:
      fprintf(fp, "integration_rule\t%s\n", "GK21");
      break;
   case RULE_TW3EV_GK31:
      fprintf(fp, "integration_rule\t%s\n", "GK31");
      break;
   case RULE_TW3EV_GK41:
      fprintf(fp, "integration_rule\t%s\n", "GK41");
      break;
   case RULE_TW3EV_GK61:
      fprintf(fp, "integration_rule\t%s\n", "GK61");
      break;
   default:
      break;
   }

   fprintf(fp, "basefolder\t%s\n", in_par->basefolder);
   fprintf(fp, "kernel_subfolder\t%s\n", in_par->ker_subfolder);
   fprintf(fp, "result_subfolder\t%s\n", in_par->res_subfolder);
   fprintf(fp, "model\t%s\n", in_par->model);

   fclose(fp);

   return;
}

bool check_input_parameters_consistency(input_parameters_t *in_par)
{
   bool res = true;
   if (in_par->xmin >= 1 || in_par->xmin <= 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for xmin, resorting to default value of %.6e", _default_in_par.xmin);
      in_par->xmin = _default_in_par.xmin;
      res = res && false;
   }

   if (in_par->mu0_2 < 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for ch_thr_mu2, resorting to default value of %.6e", _default_in_par.mu0_2);
      in_par->mu0_2 = _default_in_par.mu0_2;
      res = res && false;
   }

   if (in_par->muF_2 < 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for muF_2, resorting to default value of %.6e", _default_in_par.muF_2);
      in_par->muF_2 = _default_in_par.muF_2;
      res = res && false;
   }

   if (in_par->ch_thr_mu2 < 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for ch_thr_mu2, resorting to default value of %.6e", _default_in_par.ch_thr_mu2);
      in_par->ch_thr_mu2 = _default_in_par.ch_thr_mu2;
      res = res && false;
   }

   if (in_par->bm_thr_mu2 <= in_par->ch_thr_mu2) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for ch_thr_mu2 and/or ch_thr_mu2, resorting to default values of %.6e  %.6e", _default_in_par.ch_thr_mu2, _default_in_par.bm_thr_mu2);
      in_par->bm_thr_mu2 = _default_in_par.bm_thr_mu2;
      in_par->ch_thr_mu2 = _default_in_par.ch_thr_mu2;
      res = res && false;
   }

   if (in_par->bm_thr_mu2 < 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for ch_thr_mu2, resorting to default value of %.6e", _default_in_par.bm_thr_mu2);
      in_par->bm_thr_mu2 = _default_in_par.bm_thr_mu2;
      res = res && false;
   }

   if (in_par->grid_exponent <= 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for grid_exponent, resorting to default value of %.6e", _default_in_par.grid_exponent);
      in_par->grid_exponent = _default_in_par.grid_exponent;
      res = res && false;
   }

   if (in_par->n_thread <= 0) {
      tw3ev_log(TW3EV_WARNING, "Invalid value for n_thread, resorting to default value of %d", _default_in_par.n_thread);
      in_par->grid_exponent = _default_in_par.n_thread;
      res = res && false;
   }

#ifdef _SC_NPROCESSORS_ONLN
   int64_t n_proc = sysconf(_SC_NPROCESSORS_ONLN);
   if (n_proc <= 0) tw3ev_log(TW3EV_WARNING, "In trying to detect number of online processors some error occurred. I will ignore it and continue with the provided number of threads: %d", in_par->n_thread);
   else {
      if (in_par->n_thread >= n_proc) {
         tw3ev_log(TW3EV_WARNING, "Incompatible number of threads (%d) compared to the detected number of online cores (%ld)", in_par->n_thread, n_proc);
         if (n_proc == 1) tw3ev_log(TW3EV_WARNING, "I will resort to use the only detected online core");
         else tw3ev_log(TW3EV_WARNING, "I will resort to use all but one active cores: %ld", n_proc - 1);

         in_par->n_thread = (n_proc == 1) ? 1 : n_proc - 1;
      }
   }
#endif

#ifdef USE_OMP
   char temp_env_var[255] = "";
   sprintf(temp_env_var, "%d", in_par->n_thread);
   setenv("OMP_NUM_THREADS", temp_env_var, 1);
   omp_set_num_threads(in_par->n_thread);
#endif

   size_t bf_len = strlen(in_par->basefolder);
   if (bf_len >= BASEFOLDER_MAX_L - 1) in_par->basefolder = (char *)realloc(in_par->basefolder, sizeof(char) * (BASEFOLDER_MAX_L + 1) * 2);
   if (in_par->basefolder[bf_len - 1] != '/') {
      in_par->basefolder[bf_len] = '/';
      in_par->basefolder[bf_len + 1] = '\0';
   }

   bf_len = strlen(in_par->ker_subfolder);
   if (bf_len >= BASEFOLDER_MAX_L - 1) in_par->ker_subfolder = (char *)realloc(in_par->ker_subfolder, sizeof(char) * (BASEFOLDER_MAX_L + 1) * 2);
   while (in_par->ker_subfolder[0] == '/') {
      for (size_t i = 0; i < bf_len - 1; i++)
         in_par->ker_subfolder[i] = in_par->ker_subfolder[i + 1];
      bf_len--;
   }
   if (bf_len != 0) {
      if (in_par->ker_subfolder[bf_len - 1] != '/') {
         in_par->ker_subfolder[bf_len] = '/';
         in_par->ker_subfolder[bf_len + 1] = '\0';
      }
   }

   bf_len = strlen(in_par->res_subfolder);
   if (bf_len >= BASEFOLDER_MAX_L - 1) in_par->res_subfolder = (char *)realloc(in_par->res_subfolder, sizeof(char) * (BASEFOLDER_MAX_L + 1) * 2);
   while (in_par->res_subfolder[0] == '/') {
      for (size_t i = 0; i < bf_len - 1; i++)
         in_par->res_subfolder[i] = in_par->res_subfolder[i + 1];
      bf_len--;
   }
   if (bf_len != 0) {
      if (in_par->res_subfolder[bf_len - 1] != '/') {
         in_par->res_subfolder[bf_len] = '/';
         in_par->res_subfolder[bf_len + 1] = '\0';
      }
   }

   struct stat st = {0};
   char local_buffer[BASEFOLDER_MAX_L] = "";

   strcpy(local_buffer, in_par->basefolder);
   if (stat(local_buffer, &st) == -1) {
      tw3ev_log(TW3EV_WARNING, "Creating directory <<%s>>", local_buffer);
      mkdir(local_buffer, 0777);
   }

   strcat(local_buffer, "/dat/");
   if (stat(local_buffer, &st) == -1) {
      tw3ev_log(TW3EV_WARNING, "Creating directory <<%s>>", local_buffer);
      mkdir(local_buffer, 0777);
   }

   strcat(local_buffer, "kernels/");
   if (stat(local_buffer, &st) == -1) {
      tw3ev_log(TW3EV_WARNING, "Creating directory <<%s>>", local_buffer);
      mkdir(local_buffer, 0777);
   }
   strcat(local_buffer, in_par->ker_subfolder);

   if (stat(local_buffer, &st) == -1) {
      tw3ev_log(TW3EV_WARNING, "Creating directory <<%s>>", local_buffer);
      mkdir(local_buffer, 0777);
   }

   strcpy(local_buffer, in_par->basefolder);
   strcat(local_buffer, "/dat/results/");
   if (stat(local_buffer, &st) == -1) {
      tw3ev_log(TW3EV_WARNING, "Creating directory <<%s>>", local_buffer);
      mkdir(local_buffer, 0777);
   }
   strcat(local_buffer, in_par->res_subfolder);

   if (stat(local_buffer, &st) == -1) {
      tw3ev_log(TW3EV_WARNING, "Creating directory <<%s>>", local_buffer);
      mkdir(local_buffer, 0777);
   }

   return res;
}

input_parameters_t *init_input_parameters(const char fn[], bool save_conf_out)
{
   input_parameters_t *in_par = (input_parameters_t *)calloc(1, sizeof(input_parameters_t));

   supply_default_parameters(in_par);

   tw3ev_conf_append_option("n", INTEGER);
   tw3ev_conf_append_option("M", INTEGER);
   tw3ev_conf_append_option("mu0_2", FLOAT);
   tw3ev_conf_append_option("muF_2", FLOAT);
   tw3ev_conf_append_option("nstep", INTEGER);
   tw3ev_conf_append_option("n_thread", INTEGER);
   tw3ev_conf_append_option("charm_THR", FLOAT);
   tw3ev_conf_append_option("bottom_THR", FLOAT);
   tw3ev_conf_append_option("xmin", EXPONENTIAL);
   tw3ev_conf_append_option("grid_exponent", FLOAT);
   tw3ev_conf_append_option("interpolant_type", STRING);
   tw3ev_conf_append_option("radial_grid_type", STRING);
   tw3ev_conf_append_option("angular_grid_type", STRING);
   tw3ev_conf_append_option("chiral_even_evolution", STRING);
   tw3ev_conf_append_option("chiral_odd_evolution", STRING);
   tw3ev_conf_append_option("basefolder", STRING);
   tw3ev_conf_append_option("model", STRING);
   tw3ev_conf_append_option("integration_rule", STRING);
   tw3ev_conf_append_option("kernel_subfolder", STRING);
   tw3ev_conf_append_option("result_subfolder", STRING);

   void *(par_local[20]);

   par_local[0] = (void *)&in_par->N;
   par_local[1] = (void *)&in_par->M;
   par_local[2] = (void *)&in_par->mu0_2;
   par_local[3] = (void *)&in_par->muF_2;
   par_local[4] = (void *)&in_par->nstep;
   par_local[5] = (void *)&in_par->n_thread;
   par_local[6] = (void *)&in_par->ch_thr_mu2;
   par_local[7] = (void *)&in_par->bm_thr_mu2;
   par_local[8] = (void *)&in_par->xmin;
   par_local[9] = (void *)&in_par->grid_exponent;

   char ipo_type[BASEFOLDER_MAX_L + 1] = "";
   char rad_type[BASEFOLDER_MAX_L + 1] = "";
   char ang_type[BASEFOLDER_MAX_L + 1] = "";
   char che_type[BASEFOLDER_MAX_L + 1] = "";
   char cho_type[BASEFOLDER_MAX_L + 1] = "";
   char iru_type[BASEFOLDER_MAX_L + 1] = "";

   par_local[10] = (void *)ipo_type;
   par_local[11] = (void *)rad_type;
   par_local[12] = (void *)ang_type;
   par_local[13] = (void *)che_type;
   par_local[14] = (void *)cho_type;

   par_local[15] = (void *)in_par->basefolder;

   par_local[16] = (void *)in_par->model;
   par_local[17] = (void *)iru_type;
   par_local[18] = (void *)in_par->ker_subfolder;
   par_local[19] = (void *)in_par->res_subfolder;

   tw3ev_read_configuration_file_sstream(par_local, fn);
   (void)check_input_parameters_consistency(in_par);
 
    int _b = -1;
   _b = local_parse_bool_from_string(che_type);
   if (_b >= 0) in_par->chiral_even_active = _b;
   else tw3ev_log(TW3EV_WARNING, "Unrecognized true/false statement for chiral even active.");

   _b = local_parse_bool_from_string(cho_type);
   if (_b >= 0) in_par->chiral_odd_active = _b;
   else tw3ev_log(TW3EV_WARNING, "Unrecognized true/false statement for chiral odd active.");

   if (strcmp(ipo_type, "IT_MINUS") == 0) in_par->F_t = IT_MINUS;
   else if (strcmp(ipo_type, "IT_PLUS") == 0) in_par->F_t = IT_PLUS;
   else if (strcmp(ipo_type, "IT_BOTH") == 0) in_par->F_t = IT_BOTH;
   else tw3ev_log(TW3EV_WARNING, "Unrecognized interpolant option.");


   if (strcmp(rad_type, "LOG_GRID") == 0) in_par->G_t = LOG_GRID;
   else if (strcmp(rad_type, "PWR_GRID") == 0) in_par->G_t = PWR_GRID;
   else if (strcmp(rad_type, "IML_GRID") == 0) in_par->G_t = IML_GRID;
   else if (strcmp(rad_type, "HYP_GRID") == 0) in_par->G_t = HYP_GRID;
   else tw3ev_log(TW3EV_WARNING, "Unrecognized radial grid option.");

   if (strcmp(ang_type, "LIN_GRID") == 0) in_par->AG_t = LIN_GRID;
   else if (strcmp(ang_type, "COS_GRID") == 0) in_par->AG_t = COS_GRID;
   else tw3ev_log(TW3EV_WARNING, "Unrecognized angular grid option.");

   if (strcmp(iru_type, "GK21") == 0 || strcmp(iru_type, "gk21") == 0) in_par->i_rule = RULE_TW3EV_GK21;
   else if (strcmp(iru_type, "GK31") == 0 || strcmp(iru_type, "gk31") == 0) in_par->i_rule = RULE_TW3EV_GK31;
   else if (strcmp(iru_type, "GK41") == 0 || strcmp(iru_type, "gk41") == 0) in_par->i_rule = RULE_TW3EV_GK41;
   else if (strcmp(iru_type, "GK61") == 0 || strcmp(iru_type, "gk61") == 0) in_par->i_rule = RULE_TW3EV_GK61;
   else tw3ev_log(TW3EV_WARNING, "Unrecognized integration rule option.");

   char cwd[BASEFOLDER_MAX_L + 1];
   char *check = getcwd(cwd, sizeof(cwd));
   if (check == NULL) {
      tw3ev_log(TW3EV_ERROR, "Unable to read current working directory to save `config.out\'. This should never happen.");
   }
   strcat(cwd, "/config.out");

   if (save_conf_out) tw3ev_write_configuration_file(in_par, cwd);

   tw3ev_reset_read_config();

   return in_par;
}

evolution_functions_t *init_default_evolution_function(evolution_functions_t *ev)
{
   if (ev == NULL) ev = calloc(1, sizeof(evolution_functions_t));
   ev->Tu = Tu_test;
   ev->Td = Td_test;
   ev->Ts = Ts_test;

   ev->DTu = DTu_test;
   ev->DTd = DTd_test;
   ev->DTs = DTs_test;

   ev->Hu = Hu_test;
   ev->Hd = Hd_test;
   ev->Hs = Hs_test;

   ev->Eu = Eu_test;
   ev->Ed = Ed_test;
   ev->Es = Es_test;

   ev->TFp = TFp_test;
   ev->TFm = TFm_test;

   ev->p_PDF = NULL;

   return ev;
}

static void get_grid_in_3D_space(int32_t i, int32_t j, int32_t N, int32_t M, double c_fact, double *x1, double *x2, double *x3, double grid_exponent, nonlinear_radial_grid_type_e G_type, angular_grid_type_e AG_type)
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
      *x1 = 0;
      *x2 = 0;
      break;
   }
   *x3 = -(*x1) - (*x2);
}

static sparse_mat_t *search_init_kernel(int32_t N, int32_t M, double c_fact, const char *folder, const char *subfolder, const char base_name[], const char grid_type_string[],
                                        double (*Hk)(int32_t, int32_t, int32_t, int32_t, integration_wrapper_t *, int), int nf, printout_level_e pl, tw3ev_integration_rule_e i_rule)
{

   sparse_mat_t *sp_m;
   char *fName = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   char N_as_string[255];
   char extension[6] = ".dat";
   sprintf(N_as_string, "_N%d_M%d", N + 1, M);
   fName = clear_string(fName);
   strncpy(fName, base_name, BASEFOLDER_MAX_L);
   strcat(fName, grid_type_string);
   strcat(fName, N_as_string);

   switch (i_rule) {
   case RULE_TW3EV_GK21: {
      strcat(fName, "_QAG_GK21");
      break;
   }
   case RULE_TW3EV_GK31: {
      strcat(fName, "_QAG_GK31");
      break;
   }
   case RULE_TW3EV_GK41: {
      strcat(fName, "_QAG_GK41");
      break;
   }
   case RULE_TW3EV_GK61: {
      strcat(fName, "_QAG_GK61");
      break;
   }
   default: {
      strcat(fName, "_QAG_GK21");
      break;
   }
   }

   strcat(fName, extension);
   char subpath[255] = "dat/kernels/";
   if (strlen(subfolder) > (256 - strlen(subpath))) {
      tw3ev_log(TW3EV_WARNING, "Ignoring subfolder request for saving kernels. Hardcoded limit of %d characters exceeded.", 256 - strlen(subpath));
   } else {
      strcat(subpath, subfolder);
   }
   fName = attach_path_to_string(fName, folder, subpath);
   if (pl != PL_NONE) tw3ev_log(TW3EV_INFO, "Looking for the saved kernel: %s", fName);
   FILE *fp = fopen(fName, "r");
   if (fp != NULL) {
      fclose(fp);
      if (pl != PL_NONE) tw3ev_log(TW3EV_INFO, "Found it!");
      sp_m = read_sparse_mat(fName, N, M);

   } else {
      if (pl != PL_NONE) tw3ev_log(TW3EV_INFO, "Not found it, I'm going to compute it and save it!");
      sp_m = init_kernel(Hk, N, M, c_fact, nf, pl);
      save_sparse_mat(sp_m, fName, N, M);
   }
   if (pl != PL_NONE) {
      tw3ev_log(TW3EV_INFO, "Initialization of the kernel done");
      tw3ev_log(TW3EV_INFO, "");
   }
   free(fName);
   return sp_m;
}

double zero_function(double a, double b, double c, void *p)
{
   (void)a;
   (void)b;
   (void)c;
   (void)p;
   return 0.0;
}

static void get_S_pm(int32_t N, int32_t M, double c_fact, double (*T)(double, double, double, void *), double (*DT)(double, double, double, void *), void *p_T, double s, vector_t **vt, double grid_exponent, nonlinear_radial_grid_type_e G_t,
                     angular_grid_type_e AG_t)
{
   if (T == NULL) {
      tw3ev_log(TW3EV_WARNING, "NULL function pointer <<%s>> for quark.", GET_VAR_NAME(T));
      T = zero_function;
   }
   if (DT == NULL) {
      tw3ev_log(TW3EV_WARNING, "NULL function pointer <<%s>>.", GET_VAR_NAME(DT));
      DT = zero_function;
   }
   int32_t dim = (N + 1) * (M + 1);
   if (*vt == NULL) *vt = init_vector(dim);

   for (int32_t i = 0; i <= N; i++) {
      for (int32_t j = 0; j <= M; j++) {
         int32_t a = from_ij_to_a(i, j, M);
         double x1, x2, x3;
         get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exponent, G_t, AG_t);

         double t123 = T(x1, x2, x3, p_T);
         double dt123 = DT(x1, x2, x3, p_T);
         double t321 = T(x3, x2, x1, p_T);
         double dt321 = DT(x3, x2, x1, p_T);
         (*vt)->data[a] = t123 - dt123 + s * t321 + s * dt321;
      }
   }
}

static void get_F_pm(int32_t N, int32_t M, double c_fact, double (*T)(double, double, double, void *), void *p_T, double s, vector_t **vt, double grid_exponent, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{
   if (T == NULL) {
      tw3ev_log(TW3EV_WARNING, "NULL function pointer <<%s>> for gluon.", GET_VAR_NAME(T));
      T = zero_function;
   }
   int32_t dim = (N + 1) * (M + 1);
   if ((*vt) == NULL) (*vt) = init_vector(dim);
   for (int32_t i = 0; i <= N; i++) {
      for (int32_t j = 0; j <= M; j++) {
         int32_t a = from_ij_to_a(i, j, M);
         double x1, x2, x3;
         get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exponent, G_t, AG_t);

         double t123 = T(x1, x2, x3, p_T);
         double t132 = T(x1, x3, x2, p_T);
         double t213 = T(x2, x1, x3, p_T);
         (*vt)->data[a] = t123 - s * (t132 - t213);
      }
   }
}

static void get_HE(int32_t N, int32_t M, double c_fact, double (*HE)(double, double, double, void *), void *p_T, vector_t **vt, double grid_exponent, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{
   if (HE == NULL) {
      tw3ev_log(TW3EV_WARNING, "NULL function pointer <<%s>> for chiral-odd.", GET_VAR_NAME(HE));
      HE = zero_function;
   }
   int32_t dim = (N + 1) * (M + 1);
   if ((*vt) == NULL) (*vt) = init_vector(dim);
   for (int32_t i = 0; i <= N; i++) {
      for (int32_t j = 0; j <= M; j++) {
         int32_t a = from_ij_to_a(i, j, M);
         double x1, x2, x3;
         get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exponent, G_t, AG_t);

         double t123 = HE(x1, x2, x3, p_T);
         (*vt)->data[a] = t123;
      }
   }
}

void reset_scales(evolution_interface_t *ei, double mu0_2, double muF_2, int32_t nstep)
{
   ei->t0 = log(mu0_2);
   ei->tF = log(muF_2);
   ei->nstep = nstep;
   ei->dt = (ei->tF - ei->t0) / ((double)nstep);

   if (ei->rk != NULL) {
      ei->rk->ch_threshold = ei->ch_threshold;
      ei->rk->bm_threshold = ei->bm_threshold;
      reset_scale_rk4th(ei->rk, ei->t0, ei->dt);
   }
}

static void set_boundary_conditions(evolution_interface_t *ei, evolution_functions_t *ev, int32_t N, int32_t M)
{
   if (ei->chiral_even_active) {
      vector_t *f_p_up = NULL, *f_m_up = NULL, *f_p_do = NULL, *f_m_do = NULL, *f_p_st = NULL, *f_m_st = NULL;

      get_S_pm(N, M, ei->c_fact, ev->Tu, ev->DTu, ev->p_PDF, +1.0, &f_p_up, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_S_pm(N, M, ei->c_fact, ev->Tu, ev->DTu, ev->p_PDF, -1.0, &f_m_up, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_S_pm(N, M, ei->c_fact, ev->Td, ev->DTd, ev->p_PDF, +1.0, &f_p_do, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_S_pm(N, M, ei->c_fact, ev->Td, ev->DTd, ev->p_PDF, -1.0, &f_m_do, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_S_pm(N, M, ei->c_fact, ev->Ts, ev->DTs, ev->p_PDF, +1.0, &f_p_st, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_S_pm(N, M, ei->c_fact, ev->Ts, ev->DTs, ev->p_PDF, -1.0, &f_m_st, ei->grid_exponent, ei->G_t, ei->AG_t);

      if (ei->f_1_p == NULL) ei->f_1_p = init_vector(f_p_up->size);
      dif_vector(f_p_up, f_p_do, ei->f_1_p);

      if (ei->f_1_m == NULL) ei->f_1_m = init_vector(f_m_up->size);
      dif_vector(f_m_up, f_m_do, ei->f_1_m);

      if (ei->f_2_p == NULL) ei->f_2_p = init_vector(f_p_up->size);
      scalar_vector_mul(-2.0, f_p_st, ei->f_2_p);
      sum_vector(ei->f_2_p, f_p_up, ei->f_2_p);
      sum_vector(ei->f_2_p, f_p_do, ei->f_2_p);

      if (ei->f_2_m == NULL) ei->f_2_m = init_vector(f_m_up->size);
      scalar_vector_mul(-2.0, f_m_st, ei->f_2_m);
      sum_vector(ei->f_2_m, f_m_up, ei->f_2_m);
      sum_vector(ei->f_2_m, f_m_do, ei->f_2_m);

      if (ei->f_3_p == NULL) ei->f_3_p = init_vector(f_p_up->size);
      else set_zero_vector(ei->f_3_p);

      if (ei->f_3_m == NULL) ei->f_3_m = init_vector(f_p_up->size);
      else set_zero_vector(ei->f_3_m);

      if (ei->f_4_p == NULL) ei->f_4_p = init_vector(f_p_up->size);
      else set_zero_vector(ei->f_4_p);

      if (ei->f_4_m == NULL) ei->f_4_m = init_vector(f_p_up->size);
      else set_zero_vector(ei->f_4_m);

      get_S_pm(N, M, ei->c_fact, ev->Tu, ev->DTu, ev->p_PDF, 1.0, &ei->f_q_p, ei->grid_exponent, ei->G_t, ei->AG_t);
      sum_vector(ei->f_q_p, f_p_do, ei->f_q_p);
      sum_vector(ei->f_q_p, f_p_st, ei->f_q_p);

      copy_vector(ei->f_q_p, ei->f_3_p);
      copy_vector(ei->f_q_p, ei->f_4_p);

      get_S_pm(N, M, ei->c_fact, ev->Tu, ev->DTu, ev->p_PDF, -1.0, &ei->f_q_m, ei->grid_exponent, ei->G_t, ei->AG_t);
      sum_vector(ei->f_q_m, f_m_do, ei->f_q_m);
      sum_vector(ei->f_q_m, f_m_st, ei->f_q_m);

      copy_vector(ei->f_q_m, ei->f_3_m);
      copy_vector(ei->f_q_m, ei->f_4_m);

      get_F_pm(N, M, ei->c_fact, ev->TFp, ev->p_PDF, 1.0, &ei->f_g_p, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_F_pm(N, M, ei->c_fact, ev->TFm, ev->p_PDF, -1.0, &ei->f_g_m, ei->grid_exponent, ei->G_t, ei->AG_t);

      free_vector(&(f_p_up));
      free_vector(&(f_m_up));
      free_vector(&(f_p_do));
      free_vector(&(f_m_do));
      free_vector(&(f_p_st));
      free_vector(&(f_m_st));
   }
   if (ei->chiral_odd_active) {
      get_HE(N, M, ei->c_fact, ev->Hu, ev->p_PDF, &ei->h_up, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_HE(N, M, ei->c_fact, ev->Hd, ev->p_PDF, &ei->h_dn, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_HE(N, M, ei->c_fact, ev->Hs, ev->p_PDF, &ei->h_st, ei->grid_exponent, ei->G_t, ei->AG_t);

      get_HE(N, M, ei->c_fact, ev->Eu, ev->p_PDF, &ei->e_up, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_HE(N, M, ei->c_fact, ev->Ed, ev->p_PDF, &ei->e_dn, ei->grid_exponent, ei->G_t, ei->AG_t);
      get_HE(N, M, ei->c_fact, ev->Es, ev->p_PDF, &ei->e_st, ei->grid_exponent, ei->G_t, ei->AG_t);
   }
}

static double get_c_fact_from_xmin(double xmin, double grid_exponent, nonlinear_radial_grid_type_e G_t)
{
   switch (G_t) {
   case LOG_GRID: {
      return -log(xmin);
      break;
   }
   case PWR_GRID: {
      return pow(xmin, 1.0 / grid_exponent) / (1.0 - pow(xmin, 1.0 / grid_exponent));
      break;
   }
   case IML_GRID: {
      return -1.0 / log(xmin / (2.0 - xmin));
      break;
   }
   case HYP_GRID: {
      return -1.0 / acosh(1.0 / pow(xmin, 1.0 / grid_exponent));
      break;
   }
   default: {
      tw3ev_log(TW3EV_ERROR, "Unsupported non-linear grid type. Supported "
                             "types: LOG_GRID, PWR_GRID, IML_GRID, HYP_GRID");
      break;
   }
   }
   return 0;
}

evolution_interface_t *init_evolution_interface(double (*as)(double, void *), void *p_as, evolution_functions_t *ev, input_parameters_t *in_par, printout_level_e pl, evolution_interface_t *ei, bool start_from_previous)
{
   bool to_initialize = NULL == ei;
   (void)check_input_parameters_consistency(in_par);
   int32_t N = 6 * in_par->N - 1;
   int32_t M = in_par->M;

   double c_fact = get_c_fact_from_xmin(fabs(in_par->xmin), in_par->grid_exponent, in_par->G_t);

   if (!to_initialize) {

      to_initialize = to_initialize || (M != ei->M);
      to_initialize = to_initialize || (N != ei->N);
      to_initialize = to_initialize || (in_par->F_t != ei->F_t);
      to_initialize = to_initialize || (in_par->G_t != ei->G_t);
      to_initialize = to_initialize || (in_par->AG_t != ei->AG_t);
      to_initialize = to_initialize || (in_par->grid_exponent != ei->grid_exponent);
      to_initialize = to_initialize || (in_par->i_rule != ei->i_rule);
      to_initialize = to_initialize || (fabs(ei->xmin - in_par->xmin) > _ZERO_THR_);
      to_initialize = to_initialize || (fabs(ei->c_fact - c_fact) > _ZERO_THR_);
      to_initialize = to_initialize || (strcmp(ei->ker_subfolder, in_par->ker_subfolder) != 0);
      to_initialize = to_initialize || (strcmp(ei->res_subfolder, in_par->res_subfolder) != 0);

      if (to_initialize) free_evolution_interface(&ei);
   }

   if (to_initialize) start_from_previous = false;

   if (to_initialize) ei = (evolution_interface_t *)calloc(1, sizeof(evolution_interface_t));

   ei->xmin = fabs(in_par->xmin);
   ei->c_fact = c_fact;

   ei->pr_lev = pl;
   ei->F_t = in_par->F_t;
   ei->G_t = in_par->G_t;
   ei->AG_t = in_par->AG_t;
   ei->grid_exponent = in_par->grid_exponent;
   if (to_initialize) {
      ei->basefolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
      ei->ker_subfolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
      ei->res_subfolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   }

   strncpy(ei->basefolder, in_par->basefolder, BASEFOLDER_MAX_L);
   strncpy(ei->ker_subfolder, in_par->ker_subfolder, BASEFOLDER_MAX_L);
   strncpy(ei->res_subfolder, in_par->res_subfolder, BASEFOLDER_MAX_L);

   ei->N = N;
   ei->M = M;
   ei->as = as;
   ei->p_as = p_as;

   if (to_initialize) ei->nf = 3;
   ei->ch_threshold = log(in_par->ch_thr_mu2);
   ei->bm_threshold = log(in_par->bm_thr_mu2);
   reset_scales(ei, in_par->mu0_2, in_par->muF_2, in_par->nstep);
   ei->n_thread = in_par->n_thread;

#ifndef _SKIP_KERNELS_
   bool old_ce_act = false, old_co_act = false;
   if (!to_initialize) {
      old_ce_act = ei->chiral_even_active;
      old_co_act = ei->chiral_odd_active;
   }
#endif
   ei->chiral_even_active = in_par->chiral_even_active;
   ei->chiral_odd_active = in_par->chiral_odd_active;

   ei->i_rule = in_par->i_rule;

   switch (ei->i_rule) {
   case RULE_TW3EV_GK21: {
      set_integration_routine(integration_rule_gk21);
      break;
   }
   case RULE_TW3EV_GK31: {
      set_integration_routine(integration_rule_gk31);
      break;
   }
   case RULE_TW3EV_GK41: {
      set_integration_routine(integration_rule_gk41);
      break;
   }
   case RULE_TW3EV_GK61: {
      set_integration_routine(integration_rule_gk61);
      break;
   }
   default: {
      set_integration_routine(integration_rule_gk21);
      break;
   }
   }

   char grid_type_string[255] = "";
   switch (ei->G_t) {
   case LOG_GRID: {
      strcat(grid_type_string, "_LOG");
      break;
   }
   case PWR_GRID: {
      strcat(grid_type_string, "_PWR");
      break;
   }
   case IML_GRID: {
      strcat(grid_type_string, "_IML");
      break;
   }
   case HYP_GRID: {
      strcat(grid_type_string, "_HYP");
      break;
   }
   default: {
      tw3ev_log(TW3EV_ERROR, "Unsupported non-linear grid type. Supported "
                             "types: LOG_GRID, PWR_GRID, IML_GRID, HYP_GRID");
      break;
   }
   }

   switch (ei->AG_t) {
   case LIN_GRID: {
      strcat(grid_type_string, "_LIN");
      break;
   }
   case COS_GRID: {
      strcat(grid_type_string, "_COS");
      break;
   }
   default:
      break;
   }

   switch (ei->F_t) {
   case IT_BOTH: {
      strcat(grid_type_string, "_BOTH");
      break;
   }
   case IT_PLUS: {
      strcat(grid_type_string, "_PLUS");
      break;
   }
   case IT_MINUS: {
      strcat(grid_type_string, "_MINUS");
      break;
   }
   default:
      break;
   }

   if (ei->pr_lev == PL_ALL)
      tw3ev_log(TW3EV_INFO,
                "Initialization with xmin: %le, c_fact: %lf, N: %d, M: %d, "
                "n_thread: %d",
                ei->xmin, ei->c_fact, ei->N, ei->M, ei->n_thread);

   struct timespec start, finish;
   double elapsed;
   if (ei->pr_lev != PL_NONE) clock_gettime(CLOCK_MONOTONIC, &start);

   kernels_setup_for_computation(N, M, ei->c_fact, in_par->n_thread, in_par->F_t, in_par->G_t, in_par->AG_t, in_par->grid_exponent, in_par->basefolder);

   if ((to_initialize && ei->chiral_odd_active) || (!to_initialize && !old_co_act && ei->chiral_odd_active))
      ei->H_CO = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_CO", grid_type_string, H_CO, ei->nf, ei->pr_lev, in_par->i_rule);
   if (!to_initialize && old_co_act && !ei->chiral_odd_active) free_sparse_mat(&ei->H_CO);

   if ((to_initialize && ei->chiral_even_active) || (!to_initialize && !old_ce_act && ei->chiral_even_active)) {
      ei->H_NS = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_NS", grid_type_string, H_NS, ei->nf, ei->pr_lev, in_par->i_rule);

      ei->H_qq_p = ei->H_NS;
      ei->H_J_p = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_J_p_nf3", grid_type_string, H_J_p, 3, ei->pr_lev, in_par->i_rule);

      ei->H_qg_p = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_qg_p_nf3", grid_type_string, H_qg_p, 3, ei->pr_lev, in_par->i_rule);
      ei->H_gq_p = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_gq_p", grid_type_string, H_gq_p, ei->nf, ei->pr_lev, in_par->i_rule);
      ei->H_gg_p = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_gg_p_nf3", grid_type_string, H_gg_p, 3, ei->pr_lev, in_par->i_rule);
      ei->H_qq_m = ei->H_NS;
      ei->H_qg_m = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_qg_m_nf3", grid_type_string, H_qg_m, 3, ei->pr_lev, in_par->i_rule);
      ei->H_gq_m = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_gq_m", grid_type_string, H_gq_m, ei->nf, ei->pr_lev, in_par->i_rule);
      ei->H_gg_m = search_init_kernel(N, M, ei->c_fact, ei->basefolder, ei->ker_subfolder, "saved_kernel_H_gg_m_nf3", grid_type_string, H_gg_m, 3, ei->pr_lev, in_par->i_rule);
   }

   if (ei->pr_lev != PL_NONE) {
      clock_gettime(CLOCK_MONOTONIC, &finish);

      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      tw3ev_log(TW3EV_INFO, "Initialization time: %.5e s, with N=%d and M=%d", elapsed, N + 1, M);
   }

   if (to_initialize || !start_from_previous) set_boundary_conditions(ei, ev, N, M);

   ei->rk = init_rk4th(ei->f_1_p, ei->f_1_m, ei->f_2_p, ei->f_2_m, ei->f_3_p, ei->f_3_m, ei->f_4_p, ei->f_4_m, ei->f_q_p, ei->f_g_p, ei->f_q_m, ei->f_g_m, ei->h_up, ei->h_dn, ei->h_st, ei->e_up, ei->e_dn, ei->e_st, ei->H_NS, ei->H_CO, ei->H_J_p,
                       ei->H_qq_p, ei->H_qg_p, ei->H_gq_p, ei->H_gg_p, ei->H_qq_m, ei->H_qg_m, ei->H_gq_m, ei->H_gg_m, as, p_as, ei->t0, ei->dt, &(ei->nf), ei->chiral_even_active, ei->chiral_odd_active, ei->ch_threshold, ei->bm_threshold, ei->rk);

   kernels_free_setup(N, M);

   return ei;
}

void free_evolution_interface(evolution_interface_t **ei_dp)
{

   if (NULL == ei_dp) return;
   evolution_interface_t *ei = *ei_dp;
   if (ei == NULL) return;
   free(ei->basefolder);
   free(ei->ker_subfolder);
   free(ei->res_subfolder);
   free_rk4th(&(ei->rk));

   free_sparse_mat(&(ei->H_CO));
   free_sparse_mat(&(ei->H_NS));
   free_sparse_mat(&(ei->H_J_p));
   free_sparse_mat(&(ei->H_qg_p));
   free_sparse_mat(&(ei->H_gq_p));
   free_sparse_mat(&(ei->H_gg_p));
   free_sparse_mat(&(ei->H_qg_m));
   free_sparse_mat(&(ei->H_gq_m));
   free_sparse_mat(&(ei->H_gg_m));

   free_vector(&(ei->f_1_p));
   free_vector(&(ei->f_2_p));
   free_vector(&(ei->f_3_p));
   free_vector(&(ei->f_4_p));
   free_vector(&(ei->f_1_m));
   free_vector(&(ei->f_2_m));
   free_vector(&(ei->f_3_m));
   free_vector(&(ei->f_4_m));
   free_vector(&(ei->f_q_p));
   free_vector(&(ei->f_q_m));
   free_vector(&(ei->f_g_p));
   free_vector(&(ei->f_g_m));
   free_vector(&(ei->h_up));
   free_vector(&(ei->h_dn));
   free_vector(&(ei->h_st));
   free_vector(&(ei->e_up));
   free_vector(&(ei->e_dn));
   free_vector(&(ei->e_st));

   free(ei);
   *ei_dp = NULL;
   return;
}

saved_solution_t *execute_evolution(evolution_interface_t *ei)
{
   struct timespec start, finish;
   double elapsed;
   if (ei->pr_lev != PL_NONE) clock_gettime(CLOCK_MONOTONIC, &start);

   for (int32_t i = 0; i < ei->nstep; i++)
      rk4th_step(ei->rk);

   if (ei->pr_lev != PL_NONE) {
      clock_gettime(CLOCK_MONOTONIC, &finish);

      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      tw3ev_log(TW3EV_INFO, "Evolution time: %.5e s", elapsed);
   }

   return init_saved_solution(ei);
}

saved_solution_t **execute_evolution_w_steps(evolution_interface_t *ei)
{
   struct timespec start, finish;
   double elapsed;

   if (ei->pr_lev != PL_NONE) clock_gettime(CLOCK_MONOTONIC, &start);

   saved_solution_t **sol_arr = (saved_solution_t **)calloc(ei->nstep, sizeof(saved_solution_t *));

   for (int32_t i = 0; i < ei->nstep; i++) {
      rk4th_step(ei->rk);
      sol_arr[i] = init_saved_solution(ei);
   }

   if (ei->pr_lev != PL_NONE) {
      clock_gettime(CLOCK_MONOTONIC, &finish);

      elapsed = (finish.tv_sec - start.tv_sec);
      elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
      tw3ev_log(TW3EV_INFO, "Evolution time: %.5e s", elapsed);
   }

   return sol_arr;
}

static char *generate_file_name(int32_t N, int32_t M, int nf, const char basefolder[], const char basename[], const char specname[])
{
   char *fName = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   char N_as_string[255];
   char extension[6] = ".dat";
   sprintf(N_as_string, "_N%d_M%d", N + 1, M);

   char nf_as_string[16];
   sprintf(nf_as_string, "%d", nf);

   fName = clear_string(fName);
   strncpy(fName, basename, BASEFOLDER_MAX_L);
   strcat(fName, N_as_string);
   strcat(fName, specname);
   strcat(fName, nf_as_string);
   strcat(fName, extension);
   fName = attach_path_to_string(fName, basefolder, "");

   return fName;
}

static void single_save(int32_t N, int32_t M, double c_fact, int nf, const char basefolder[], const char basename[], const char specname[], vector_t *toSave, double grid_exponent, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{
   char *fName = generate_file_name(N, M, nf, basefolder, basename, specname);

   FILE *fp = fopen(fName, "w");
   int32_t i, j;
   for (int32_t a = 0; a < toSave->size; a++) {
      from_a_to_ij(a, M, &i, &j);
      double x1, x2, x3;
      get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exponent, G_t, AG_t);

      fprintf(fp, "%d\t%d\t%.16e\t%.16e\t%.16e\t%.16e\n", i, j, x1, x2, x3, toSave->data[a]);
   }
   fclose(fp);

   free(fName);
}

static void single_load(int32_t N, int32_t M, double c_fact, int nf, const char basefolder[], const char basename[], const char specname[], vector_t *toSave)
{
   (void)c_fact;
   char *fName = generate_file_name(N, M, nf, basefolder, basename, specname);

   FILE *fp = fopen(fName, "r");
   if (!fp) tw3ev_log(TW3EV_ERROR, "Unable to locate the file <<%s>>. Please ensure that the file is correctly placed in the dat folder.", fName);

   char buffer[BASEFOLDER_MAX_L + 1];
   int32_t i, j;
   while (NULL != fgets(buffer, BASEFOLDER_MAX_L, fp)) {

      double x1, x2, x3, temp;
      sscanf(buffer, "%d\t%d\t%le\t%le\t%le\t%le\n", &i, &j, &x1, &x2, &x3, &temp);

      if (i > N || j > M) {
         tw3ev_log(TW3EV_ERROR, "Out-of-range indexe(s) %d %d in file <<%s>>.", i, j, fName);
      }
      int32_t a = from_ij_to_a(i, j, M);

      toSave->data[a] = temp;
   }

   fclose(fp);
   free(fName);
}

static int32_t get_permuted_a_index(double x1, double x2, double x3, double c_fact, int32_t N, int32_t M, double grid_exp, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{
   int32_t i_per, j_per;
   double dummy_i, dummy_j;
   get_3D_in_grid_space(x1, x2, x3, N, M, c_fact, &dummy_i, &dummy_j, grid_exp, G_t, AG_t);

   i_per = round(dummy_i);
   j_per = round(dummy_j);
   double norm_i = (i_per != 0 ? (double)abs(i_per) : 1);
   double norm_j = (j_per != 0 ? (double)abs(j_per) : 1);

   if (fabs(i_per - dummy_i) > 10 * _ZERO_THR_ * norm_i) tw3ev_log(TW3EV_ERROR, "Something went wrong with the grids. Unable to get permuted indexes. [i]");
   if (fabs(j_per - dummy_j) > 10 * _ZERO_THR_ * norm_j) tw3ev_log(TW3EV_ERROR, "Something went wrong with the grids. Unable to get permuted indexes. [j]");
   if (i_per < 0) i_per += N + 1;
   if (i_per > N) i_per -= N + 1;

   return from_ij_to_a(i_per, j_per, M);
}

static void single_save_physical(int32_t N, int32_t M, double c_fact, int nf, const char basefolder[], const char basename[], const char specname[], vector_t *T_in, vector_t *DT_in, int dist_type, double grid_exponent,
                                 nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{
   char specname_T[255] = "";
   char specname_DT[255] = "";

   if (_QUARK_COMBINATION_TYPE_ == dist_type) {
      strcpy(specname_T, "_T");
      strcpy(specname_DT, "_DT");
   } else if (_GLUON_COMBINATION_TYPE_ == dist_type) {
      strcpy(specname_T, "_Tp");
      strcpy(specname_DT, "_Tm");
   }
   strcat(specname_T, specname);
   strcat(specname_DT, specname);

   char *fName_T = generate_file_name(N, M, nf, basefolder, basename, specname_T);
   char *fName_DT = generate_file_name(N, M, nf, basefolder, basename, specname_DT);

   FILE *fp_T = fopen(fName_T, "w");
   FILE *fp_DT = fopen(fName_DT, "w");
   int32_t i, j;

   for (int32_t a = 0; a < T_in->size; a++) {
      from_a_to_ij(a, M, &i, &j);
      double x1, x2, x3;
      get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exponent, G_t, AG_t);

      fprintf(fp_T, "%d\t%d\t%.16e\t%.16e\t%.16e\t%.16e\n", i, j, x1, x2, x3, T_in->data[a]);
      fprintf(fp_DT, "%d\t%d\t%.16e\t%.16e\t%.16e\t%.16e\n", i, j, x1, x2, x3, DT_in->data[a]);
   }

   fclose(fp_T);
   fclose(fp_DT);

   free(fName_T);
   free(fName_DT);
}

static void single_load_physical(int32_t N, int32_t M, double c_fact, int nf, const char basefolder[], const char basename[], const char specname[], vector_t *T_in, vector_t *DT_in, vector_t *v_pl, vector_t *v_mn, int dist_type, double grid_exp,
                                 nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{

   char specname_T[255] = "";
   char specname_DT[255] = "";

   if (_QUARK_COMBINATION_TYPE_ == dist_type) {
      strcpy(specname_T, "_T");
      strcpy(specname_DT, "_DT");
   } else if (_GLUON_COMBINATION_TYPE_ == dist_type) {

      strcpy(specname_T, "_Tp");
      strcpy(specname_DT, "_Tm");
   }
   strcat(specname_T, specname);
   strcat(specname_DT, specname);

   char *fName_T = generate_file_name(N, M, nf, basefolder, basename, specname_T);
   char *fName_DT = generate_file_name(N, M, nf, basefolder, basename, specname_DT);

   FILE *fp_T = fopen(fName_T, "r");
   if (!fp_T) tw3ev_log(TW3EV_ERROR, "Unable to locate the file <<%s>>. Please ensure that the file is correctly placed in the dat folder.", fName_T);

   char buffer[BASEFOLDER_MAX_L + 1];
   while (NULL != fgets(buffer, BASEFOLDER_MAX_L, fp_T)) {

      int32_t i, j;
      double x1, x2, x3, temp;
      sscanf(buffer, "%d\t%d\t%le\t%le\t%le\t%le\n", &i, &j, &x1, &x2, &x3, &temp);

      if (i > N || j > M) {
         tw3ev_log(TW3EV_ERROR, "Out-of-range indexe(s) %d %d in file <<%s>>.", i, j, fName_T);
      }
      int32_t a = from_ij_to_a(i, j, M);

      T_in->data[a] = temp;
   }
   fclose(fp_T);

   FILE *fp_DT = fopen(fName_DT, "r");
   if (!fp_DT) tw3ev_log(TW3EV_ERROR, "Unable to locate the file <<%s>>. Please ensure that the file is correctly placed in the dat folder.", fName_DT);

   while (NULL != fgets(buffer, BASEFOLDER_MAX_L, fp_DT)) {
      int32_t i, j;

      double x1, x2, x3, temp;
      sscanf(buffer, "%d\t%d\t%le\t%le\t%le\t%le\n", &i, &j, &x1, &x2, &x3, &temp);

      if (i > N || j > M) {
         tw3ev_log(TW3EV_ERROR, "Out-of-range indexe(s) %d %d in file <<%s>>.", i, j, fName_DT);
      }
      int32_t a = from_ij_to_a(i, j, M);

      DT_in->data[a] = temp;
   }
   fclose(fp_DT);

   for (int32_t a = 0; a < v_pl->size; a++) {
      int32_t i, j;
      int32_t a_321;
      int32_t a_213, a_132;

      from_a_to_ij(a, M, &i, &j);
      double x1, x2, x3;
      get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exp, G_t, AG_t);

      a_321 = get_permuted_a_index(x3, x2, x1, c_fact, N, M, grid_exp, G_t, AG_t);
      a_213 = get_permuted_a_index(x2, x1, x3, c_fact, N, M, grid_exp, G_t, AG_t);
      a_132 = get_permuted_a_index(x1, x3, x2, c_fact, N, M, grid_exp, G_t, AG_t);

      if (_QUARK_COMBINATION_TYPE_ == dist_type) {
         v_pl->data[a] = T_in->data[a] - DT_in->data[a] + T_in->data[a_321] + DT_in->data[a_321];
         v_mn->data[a] = T_in->data[a] - DT_in->data[a] - T_in->data[a_321] - DT_in->data[a_321];
      } else if (_GLUON_COMBINATION_TYPE_ == dist_type) {
         v_pl->data[a] = T_in->data[a] + T_in->data[a_132] - T_in->data[a_213];
         v_mn->data[a] = DT_in->data[a] - DT_in->data[a_132] + DT_in->data[a_213];
      }
   }

   free(fName_T);
   free(fName_DT);
}

static void generate_physical_combinations(int32_t N, int32_t M, double c_fact, vector_t *v_pl, vector_t *v_mn, int dist_type, vector_t *T_out, vector_t *DT_out, double grid_exp, nonlinear_radial_grid_type_e G_t, angular_grid_type_e AG_t)
{
   int32_t i, j, a_321;

   for (int32_t a = 0; a < v_pl->size; a++) {
      from_a_to_ij(a, M, &i, &j);
      double x1, x2, x3;
      get_grid_in_3D_space(i, j, N, M, c_fact, &x1, &x2, &x3, grid_exp, G_t, AG_t);

      a_321 = get_permuted_a_index(x3, x2, x1, c_fact, N, M, grid_exp, G_t, AG_t);

      double res_T = 0, res_DT = 0;
      if (_QUARK_COMBINATION_TYPE_ == dist_type) {
         res_T = 0.25 * (v_pl->data[a] + v_pl->data[a_321] + v_mn->data[a] - v_mn->data[a_321]);
         res_DT = -0.25 * (v_pl->data[a] - v_pl->data[a_321] + v_mn->data[a] + v_mn->data[a_321]);
      } else if (_GLUON_COMBINATION_TYPE_ == dist_type) {
         res_T = 0.5 * (v_pl->data[a] - v_pl->data[a_321]);  
         res_DT = 0.5 * (v_mn->data[a] + v_mn->data[a_321]); 
      }

      T_out->data[a] = res_T;
      DT_out->data[a] = res_DT;
   }
}

saved_solution_t *init_saved_solution(evolution_interface_t *ei)
{
   uint32_t dim = (ei->N + 1) * (ei->M + 1);
   saved_solution_t *sol = (saved_solution_t *)calloc(1, sizeof(saved_solution_t));
   sol->N = ei->N;
   sol->M = ei->M;
   sol->xmin = ei->xmin;
   sol->c_fact = ei->c_fact;
   sol->grid_exponent = ei->grid_exponent;
   sol->mu_2 = exp(ei->rk->t);
   sol->F_t = ei->F_t;
   sol->G_t = ei->G_t;
   sol->AG_t = ei->AG_t;
   sol->chiral_even_active = ei->chiral_even_active;
   sol->chiral_odd_active = ei->chiral_odd_active;

   sol->nf = (int)(ei->nf);
   sol->basefolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   sol->res_subfolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   strncpy(sol->basefolder, ei->basefolder, BASEFOLDER_MAX_L);
   strncpy(sol->res_subfolder, ei->res_subfolder, BASEFOLDER_MAX_L);

   if (sol->chiral_even_active) {
      sol->f_1_p = init_vector(dim);
      sol->f_2_p = init_vector(dim);
      sol->f_3_p = init_vector(dim);
      sol->f_4_p = init_vector(dim);

      sol->f_1_m = init_vector(dim);
      sol->f_2_m = init_vector(dim);
      sol->f_3_m = init_vector(dim);
      sol->f_4_m = init_vector(dim);

      sol->f_q_p = init_vector(dim);
      sol->f_q_m = init_vector(dim);
      sol->f_g_p = init_vector(dim);
      sol->f_g_m = init_vector(dim);

      sol->up_p = init_vector(dim);
      sol->dn_p = init_vector(dim);
      sol->st_p = init_vector(dim);
      sol->ch_p = init_vector(dim);
      sol->bm_p = init_vector(dim);
      sol->gl_p = sol->f_g_p;

      sol->up_m = init_vector(dim);
      sol->dn_m = init_vector(dim);
      sol->st_m = init_vector(dim);
      sol->ch_m = init_vector(dim);
      sol->bm_m = init_vector(dim);
      sol->gl_m = sol->f_g_m;

      sol->T_up = init_vector(dim);
      sol->T_dn = init_vector(dim);
      sol->T_st = init_vector(dim);
      sol->T_ch = init_vector(dim);
      sol->T_bm = init_vector(dim);
      sol->Tp_gl = init_vector(dim);

      sol->DT_up = init_vector(dim);
      sol->DT_dn = init_vector(dim);
      sol->DT_st = init_vector(dim);
      sol->DT_ch = init_vector(dim);
      sol->DT_bm = init_vector(dim);
      sol->Tm_gl = init_vector(dim);
   }
   if (sol->chiral_odd_active) {
      sol->e_up = init_vector(dim);
      sol->e_dn = init_vector(dim);
      sol->e_st = init_vector(dim);

      sol->h_up = init_vector(dim);
      sol->h_dn = init_vector(dim);
      sol->h_st = init_vector(dim);
   }

   if (sol->chiral_even_active) {
      memcpy(sol->gl_p->data, ei->f_g_p->data, sizeof(double) * dim);
      memcpy(sol->gl_m->data, ei->f_g_m->data, sizeof(double) * dim);
      memcpy(sol->f_1_p->data, ei->f_1_p->data, sizeof(double) * dim);
      memcpy(sol->f_2_p->data, ei->f_2_p->data, sizeof(double) * dim);
      memcpy(sol->f_3_p->data, ei->f_3_p->data, sizeof(double) * dim);
      memcpy(sol->f_4_p->data, ei->f_4_p->data, sizeof(double) * dim);
      memcpy(sol->f_1_m->data, ei->f_1_m->data, sizeof(double) * dim);
      memcpy(sol->f_2_m->data, ei->f_2_m->data, sizeof(double) * dim);
      memcpy(sol->f_3_m->data, ei->f_3_m->data, sizeof(double) * dim);
      memcpy(sol->f_4_m->data, ei->f_4_m->data, sizeof(double) * dim);
      memcpy(sol->f_q_p->data, ei->f_q_p->data, sizeof(double) * dim);
      memcpy(sol->f_q_m->data, ei->f_q_m->data, sizeof(double) * dim);
   }
   if (sol->chiral_odd_active) {
      memcpy(sol->e_up->data, ei->e_up->data, sizeof(double) * dim);
      memcpy(sol->e_dn->data, ei->e_dn->data, sizeof(double) * dim);
      memcpy(sol->e_st->data, ei->e_st->data, sizeof(double) * dim);
      memcpy(sol->h_up->data, ei->h_up->data, sizeof(double) * dim);
      memcpy(sol->h_dn->data, ei->h_dn->data, sizeof(double) * dim);
      memcpy(sol->h_st->data, ei->h_st->data, sizeof(double) * dim);
   }

   if (!sol->chiral_even_active) return sol;

   switch ((int)ei->nf) {
   case 3: {
      for (uint32_t i = 0; i < dim; i++) {
         sol->up_p->data[i] = (3.0 * ei->f_1_p->data[i] + ei->f_2_p->data[i] + 2.0 * ei->f_q_p->data[i]) / 6.0;

         sol->dn_p->data[i] = (-3.0 * ei->f_1_p->data[i] + ei->f_2_p->data[i] + 2.0 * ei->f_q_p->data[i]) / 6.0;
         sol->st_p->data[i] = (-2.0 * ei->f_2_p->data[i] + 2.0 * ei->f_q_p->data[i]) / 6.0;

         sol->ch_p->data[i] = ei->f_3_p->data[i];
         sol->bm_p->data[i] = ei->f_4_p->data[i];

         sol->up_m->data[i] = (3.0 * ei->f_1_m->data[i] + ei->f_2_m->data[i] + 2.0 * ei->f_q_m->data[i]) / 6.0;

         sol->dn_m->data[i] = (-3.0 * ei->f_1_m->data[i] + ei->f_2_m->data[i] + 2.0 * ei->f_q_m->data[i]) / 6.0;
         sol->st_m->data[i] = (-2.0 * ei->f_2_m->data[i] + 2.0 * ei->f_q_m->data[i]) / 6.0;
         sol->ch_m->data[i] = ei->f_3_m->data[i];
         sol->bm_m->data[i] = ei->f_4_m->data[i];
      }
      break;
   }
   case 4: {
      for (uint32_t i = 0; i < dim; i++) {
         sol->up_p->data[i] = (6.0 * ei->f_1_p->data[i] + 2.0 * ei->f_2_p->data[i] + ei->f_3_p->data[i] + 3.0 * ei->f_q_p->data[i]) / 12.0;
         sol->dn_p->data[i] = (-6.0 * ei->f_1_p->data[i] + 2.0 * ei->f_2_p->data[i] + ei->f_3_p->data[i] + 3.0 * ei->f_q_p->data[i]) / 12.0;
         sol->st_p->data[i] = (-4.0 * ei->f_2_p->data[i] + ei->f_3_p->data[i] + 3.0 * ei->f_q_p->data[i]) / 12.0;
         sol->ch_p->data[i] = (-ei->f_3_p->data[i] + ei->f_q_p->data[i]) / 4.0;
         sol->bm_p->data[i] = ei->f_4_p->data[i];

         sol->up_m->data[i] = (6.0 * ei->f_1_m->data[i] + 2.0 * ei->f_2_m->data[i] + ei->f_3_m->data[i] + 3.0 * ei->f_q_m->data[i]) / 12.0;
         sol->dn_m->data[i] = (-6.0 * ei->f_1_m->data[i] + 2.0 * ei->f_2_m->data[i] + ei->f_3_m->data[i] + 3.0 * ei->f_q_m->data[i]) / 12.0;
         sol->st_m->data[i] = (-4.0 * ei->f_2_m->data[i] + ei->f_3_m->data[i] + 3.0 * ei->f_q_m->data[i]) / 12.0;
         sol->ch_m->data[i] = (-ei->f_3_m->data[i] + ei->f_q_m->data[i]) / 4.0;
         sol->bm_m->data[i] = ei->f_4_m->data[i];
      }
      break;
   }
   case 5: {

      for (uint32_t i = 0; i < dim; i++) {
         sol->up_p->data[i] = (30.0 * ei->f_1_p->data[i] + 10.0 * ei->f_2_p->data[i] + 5.0 * ei->f_3_p->data[i] + 3.0 * ei->f_4_p->data[i] + 12.0 * ei->f_q_p->data[i]) / 60.0;
         sol->dn_p->data[i] = (-30.0 * ei->f_1_p->data[i] + 10.0 * ei->f_2_p->data[i] + 5.0 * ei->f_3_p->data[i] + 3.0 * ei->f_4_p->data[i] + 12.0 * ei->f_q_p->data[i]) / 60.0;
         sol->st_p->data[i] = (-20.0 * ei->f_2_p->data[i] + 5.0 * ei->f_3_p->data[i] + 3.0 * ei->f_4_p->data[i] + 12.0 * ei->f_q_p->data[i]) / 60.0;
         sol->ch_p->data[i] = (-5.0 * ei->f_3_p->data[i] + ei->f_4_p->data[i] + 4.0 * ei->f_q_p->data[i]) / 20.0;
         sol->bm_p->data[i] = (-ei->f_4_p->data[i] + ei->f_q_p->data[i]) / 5.0;

         sol->up_m->data[i] = (30.0 * ei->f_1_m->data[i] + 10.0 * ei->f_2_m->data[i] + 5.0 * ei->f_3_m->data[i] + 3.0 * ei->f_4_m->data[i] + 12.0 * ei->f_q_m->data[i]) / 60.0;
         sol->dn_m->data[i] = (-30.0 * ei->f_1_m->data[i] + 10.0 * ei->f_2_m->data[i] + 5.0 * ei->f_3_m->data[i] + 3.0 * ei->f_4_m->data[i] + 12.0 * ei->f_q_m->data[i]) / 60.0;
         sol->st_m->data[i] = (-20.0 * ei->f_2_m->data[i] + 5.0 * ei->f_3_m->data[i] + 3.0 * ei->f_4_m->data[i] + 12.0 * ei->f_q_m->data[i]) / 60.0;
         sol->ch_m->data[i] = (-5.0 * ei->f_3_m->data[i] + ei->f_4_m->data[i] + 4.0 * ei->f_q_m->data[i]) / 20.0;
         sol->bm_m->data[i] = (-ei->f_4_m->data[i] + ei->f_q_m->data[i]) / 5.0;
      }

      break;
   }
   default:
      break;
   }

   if (sol->chiral_even_active) {
      generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->up_p, sol->up_m, _QUARK_COMBINATION_TYPE_, sol->T_up, sol->DT_up, sol->grid_exponent, sol->G_t, sol->AG_t);
      generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->dn_p, sol->dn_m, _QUARK_COMBINATION_TYPE_, sol->T_dn, sol->DT_dn, sol->grid_exponent, sol->G_t, sol->AG_t);
      generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->st_p, sol->st_m, _QUARK_COMBINATION_TYPE_, sol->T_st, sol->DT_st, sol->grid_exponent, sol->G_t, sol->AG_t);
      generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->ch_p, sol->ch_m, _QUARK_COMBINATION_TYPE_, sol->T_ch, sol->DT_ch, sol->grid_exponent, sol->G_t, sol->AG_t);
      generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->bm_p, sol->bm_m, _QUARK_COMBINATION_TYPE_, sol->T_bm, sol->DT_bm, sol->grid_exponent, sol->G_t, sol->AG_t);

      generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->gl_p, sol->gl_m, _GLUON_COMBINATION_TYPE_, sol->Tp_gl, sol->Tm_gl, sol->grid_exponent, sol->G_t, sol->AG_t);
   }

   return sol;
}

saved_solution_t *load_model(input_parameters_t *in_par, const char basename_[], int nf, saved_basis_e SB)
{

   if (BASIS_BOTH == SB) {
      tw3ev_log(TW3EV_WARNING, "Loading both basis cannot be done (the loaded one is converted into the other). I will try to load physical basis.");
      SB = BASIS_PHYSICAL;
   }
   if (BASIS_EVOLUTION == SB) {
      tw3ev_log(TW3EV_WARNING, "Loading evolution basis is for debug purposes only. I will try to load physical basis.");
      SB = BASIS_PHYSICAL;
   }
   saved_solution_t *sol = (saved_solution_t *)calloc(1, sizeof(saved_solution_t));
   sol->N = 6 * in_par->N - 1;
   sol->M = in_par->M;
   int32_t dim = (sol->N + 1) * (sol->M + 1);
   sol->xmin = in_par->xmin;
   sol->grid_exponent = in_par->grid_exponent;
   sol->c_fact = get_c_fact_from_xmin(in_par->xmin, in_par->grid_exponent, in_par->G_t);
   sol->F_t = in_par->F_t;
   sol->G_t = in_par->G_t;
   sol->AG_t = in_par->AG_t;

   sol->chiral_even_active = in_par->chiral_even_active;
   sol->chiral_odd_active = in_par->chiral_odd_active;

   sol->nf = nf;
   sol->basefolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   strncpy(sol->basefolder, in_par->basefolder, BASEFOLDER_MAX_L);
   sol->res_subfolder = (char *)calloc(BASEFOLDER_MAX_L + 1, sizeof(char));
   strncpy(sol->res_subfolder, in_par->res_subfolder, BASEFOLDER_MAX_L);

   char subfolder[BASEFOLDER_MAX_L] = "/dat/results/";
   strcat(subfolder, sol->res_subfolder);

   char basename[BASEFOLDER_MAX_L] = "";
   strcpy(basename, subfolder);
   strcat(basename, basename_);

   if (sol->chiral_even_active) {
      sol->f_1_p = init_vector(dim);
      sol->f_2_p = init_vector(dim);
      sol->f_3_p = init_vector(dim);
      sol->f_4_p = init_vector(dim);

      sol->f_1_m = init_vector(dim);
      sol->f_2_m = init_vector(dim);
      sol->f_3_m = init_vector(dim);
      sol->f_4_m = init_vector(dim);

      sol->f_q_p = init_vector(dim);
      sol->f_q_m = init_vector(dim);
      sol->f_g_p = init_vector(dim);
      sol->f_g_m = init_vector(dim);

      sol->up_p = init_vector(dim);
      sol->dn_p = init_vector(dim);
      sol->st_p = init_vector(dim);
      sol->ch_p = init_vector(dim);
      sol->bm_p = init_vector(dim);
      sol->gl_p = sol->f_g_p;

      sol->up_m = init_vector(dim);
      sol->dn_m = init_vector(dim);
      sol->st_m = init_vector(dim);
      sol->ch_m = init_vector(dim);
      sol->bm_m = init_vector(dim);
      sol->gl_m = sol->f_g_m;

      sol->T_up = init_vector(dim);
      sol->T_dn = init_vector(dim);
      sol->T_st = init_vector(dim);
      sol->T_ch = init_vector(dim);
      sol->T_bm = init_vector(dim);
      sol->Tp_gl = init_vector(dim);

      sol->DT_up = init_vector(dim);
      sol->DT_dn = init_vector(dim);
      sol->DT_st = init_vector(dim);
      sol->DT_ch = init_vector(dim);
      sol->DT_bm = init_vector(dim);
      sol->Tm_gl = init_vector(dim);
   }
   if (sol->chiral_odd_active) {
      sol->e_up = init_vector(dim);
      sol->e_dn = init_vector(dim);
      sol->e_st = init_vector(dim);

      sol->h_up = init_vector(dim);
      sol->h_dn = init_vector(dim);
      sol->h_st = init_vector(dim);
   }
   if (sol->chiral_even_active) {
      if (BASIS_DEFINITE_C_PAR == SB) {
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_up_nf", sol->up_p);
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_dn_nf", sol->dn_p);
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_st_nf", sol->st_p);
         if (nf > 3) single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_ch_nf", sol->ch_p);
         if (nf > 4) single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_bm_nf", sol->bm_p);
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_gl_nf", sol->gl_p);

         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_up_nf", sol->up_m);
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_dn_nf", sol->dn_m);
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_st_nf", sol->st_m);
         if (nf > 3) single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_ch_nf", sol->ch_m);
         if (nf > 4) single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_bm_nf", sol->bm_m);
         single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_gl_nf", sol->gl_m);

         generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->up_p, sol->up_m, _QUARK_COMBINATION_TYPE_, sol->T_up, sol->DT_up, sol->grid_exponent, sol->G_t, sol->AG_t);
         generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->dn_p, sol->dn_m, _QUARK_COMBINATION_TYPE_, sol->T_dn, sol->DT_dn, sol->grid_exponent, sol->G_t, sol->AG_t);
         generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->st_p, sol->st_m, _QUARK_COMBINATION_TYPE_, sol->T_st, sol->DT_st, sol->grid_exponent, sol->G_t, sol->AG_t);
         generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->ch_p, sol->ch_m, _QUARK_COMBINATION_TYPE_, sol->T_ch, sol->DT_ch, sol->grid_exponent, sol->G_t, sol->AG_t);
         generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->bm_p, sol->bm_m, _QUARK_COMBINATION_TYPE_, sol->T_bm, sol->DT_bm, sol->grid_exponent, sol->G_t, sol->AG_t);

         generate_physical_combinations(sol->N, sol->M, sol->c_fact, sol->gl_p, sol->gl_m, _GLUON_COMBINATION_TYPE_, sol->Tp_gl, sol->Tm_gl, sol->grid_exponent, sol->G_t, sol->AG_t);
      }
      if (BASIS_PHYSICAL == SB) {

         single_load_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_up_nf", sol->T_up, sol->DT_up, sol->up_p, sol->up_m, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_load_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_dn_nf", sol->T_dn, sol->DT_dn, sol->dn_p, sol->dn_m, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_load_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_st_nf", sol->T_st, sol->DT_st, sol->st_p, sol->st_m, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_load_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_ch_nf", sol->T_ch, sol->DT_ch, sol->ch_p, sol->ch_m, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_load_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_bm_nf", sol->T_bm, sol->DT_bm, sol->bm_p, sol->bm_m, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);

         single_load_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_gl_nf", sol->Tp_gl, sol->Tm_gl, sol->gl_p, sol->gl_m, _GLUON_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
      }
      for (size_t a = 0; a < (size_t)sol->up_p->size; a++) {
         sol->f_1_p->data[a] = sol->up_p->data[a] - sol->dn_p->data[a];
         sol->f_2_p->data[a] = sol->up_p->data[a] + sol->dn_p->data[a] - 2.0 * sol->st_p->data[a];
         sol->f_3_p->data[a] = sol->up_p->data[a] + sol->dn_p->data[a] + sol->st_p->data[a] - 3.0 * sol->ch_p->data[a];
         sol->f_4_p->data[a] = sol->up_p->data[a] + sol->dn_p->data[a] + sol->st_p->data[a] + sol->ch_p->data[a] - 4.0 * sol->bm_p->data[a];
         sol->f_q_p->data[a] = sol->up_p->data[a] + sol->dn_p->data[a] + sol->st_p->data[a] + sol->ch_p->data[a] + sol->bm_p->data[a];

         sol->f_1_m->data[a] = sol->up_m->data[a] - sol->dn_m->data[a];
         sol->f_2_m->data[a] = sol->up_m->data[a] + sol->dn_m->data[a] - 2.0 * sol->st_m->data[a];
         sol->f_3_m->data[a] = sol->up_m->data[a] + sol->dn_m->data[a] + sol->st_m->data[a] - 3.0 * sol->ch_m->data[a];
         sol->f_4_m->data[a] = sol->up_m->data[a] + sol->dn_m->data[a] + sol->st_m->data[a] + sol->ch_m->data[a] - 4.0 * sol->bm_m->data[a];
         sol->f_q_m->data[a] = sol->up_m->data[a] + sol->dn_m->data[a] + sol->st_m->data[a] + sol->ch_m->data[a] + sol->bm_m->data[a];
      }
   }
   if (sol->chiral_odd_active) {
      single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_E_up_nf", sol->e_up);
      single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_E_dn_nf", sol->e_dn);
      single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_E_st_nf", sol->e_st);
      single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_H_up_nf", sol->h_up);
      single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_H_dn_nf", sol->h_dn);
      single_load(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_H_st_nf", sol->h_st);
   }

   return sol;
}

void free_saved_solution(saved_solution_t **sol_dp)
{
   if (NULL == sol_dp) return;
   saved_solution_t *sol = *sol_dp;
   if (NULL == sol) return;
   if (sol->chiral_even_active) {
      free_vector(&(sol->up_p));
      free_vector(&(sol->dn_p));
      free_vector(&(sol->st_p));
      free_vector(&(sol->ch_p));
      free_vector(&(sol->bm_p));

      free_vector(&(sol->up_m));
      free_vector(&(sol->dn_m));
      free_vector(&(sol->st_m));
      free_vector(&(sol->ch_m));
      free_vector(&(sol->bm_m));

      free_vector(&(sol->f_1_p));
      free_vector(&(sol->f_2_p));
      free_vector(&(sol->f_3_p));
      free_vector(&(sol->f_4_p));
      free_vector(&(sol->f_1_m));
      free_vector(&(sol->f_2_m));
      free_vector(&(sol->f_3_m));
      free_vector(&(sol->f_4_m));
      free_vector(&(sol->f_q_p));
      free_vector(&(sol->f_q_m));
      free_vector(&(sol->f_g_p));
      free_vector(&(sol->f_g_m));

      free_vector(&(sol->T_up));
      free_vector(&(sol->T_dn));
      free_vector(&(sol->T_st));
      free_vector(&(sol->T_ch));
      free_vector(&(sol->T_bm));
      free_vector(&(sol->Tp_gl));

      free_vector(&(sol->DT_up));
      free_vector(&(sol->DT_dn));
      free_vector(&(sol->DT_st));
      free_vector(&(sol->DT_ch));
      free_vector(&(sol->DT_bm));
      free_vector(&(sol->Tm_gl));
   }
   if (sol->chiral_odd_active) {
      free_vector(&(sol->e_up));
      free_vector(&(sol->e_dn));
      free_vector(&(sol->e_st));
      free_vector(&(sol->h_up));
      free_vector(&(sol->h_dn));
      free_vector(&(sol->h_st));
   }

   free(sol->basefolder);
   free(sol->res_subfolder);
   free(sol);
   *sol_dp = NULL;
}

void save_model(saved_solution_t *sol, const char basename_[], saved_basis_e SB)
{
   if (NULL == sol) {
      tw3ev_log(TW3EV_WARNING, "Trying to save solution, but NULL ptr as input has been provided. I will ignore this call.");
      return;
   }

   char subfolder[BASEFOLDER_MAX_L] = "/dat/results/";
   strcat(subfolder, sol->res_subfolder);

   char buffer[2 * BASEFOLDER_MAX_L];
   strcpy(buffer, basename_);
   strcat(buffer, "_card.txt");
   attach_path_to_string(buffer, sol->basefolder, subfolder);

   char basename[BASEFOLDER_MAX_L] = "";
   strcpy(basename, subfolder);
   strcat(basename, basename_);


   FILE *fp = fopen(buffer, "w");
   fprintf(fp, "mu_2:\t%.6e\n", sol->mu_2);
   fprintf(fp, "grid_exponent:\t%.6e\n", sol->grid_exponent);
   fprintf(fp, "xmin:\t%.6e\n", sol->xmin);
   switch (sol->G_t) {
   case LOG_GRID:
      fprintf(fp, "radial_grid_type:\t%s\n", "LOG_GRID");
      break;
   case PWR_GRID:
      fprintf(fp, "radial_grid_type:\t%s\n", "PWR_GRID");
      break;
   case IML_GRID:
      fprintf(fp, "radial_grid_type:\t%s\n", "IML_GRID");
      break;
   case HYP_GRID:
      fprintf(fp, "radial_grid_type:\t%s\n", "HYP_GRID");
      break;
   default:
      break;
   }

   switch (sol->AG_t) {
   case LIN_GRID:
      fprintf(fp, "interpolant_type:\t%s\n", "LIN_GRID");
      break;
   case COS_GRID:
      fprintf(fp, "interpolant_type:\t%s\n", "COS_GRID");
      break;
   default:
      break;
   }

   fclose(fp);

   if (sol->chiral_even_active) {

      if (BASIS_EVOLUTION == SB) {
         tw3ev_log(TW3EV_WARNING, "Saving evolution basis should be done solely for debug purposes! I save also physical basis along side it.");
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_1_p_nf", sol->f_1_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_2_p_nf", sol->f_2_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_3_p_nf", sol->f_3_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_4_p_nf", sol->f_4_p, sol->grid_exponent, sol->G_t, sol->AG_t);

         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_1_m_nf", sol->f_1_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_2_m_nf", sol->f_2_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_3_m_nf", sol->f_3_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_4_m_nf", sol->f_4_m, sol->grid_exponent, sol->G_t, sol->AG_t);

         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_q_p_nf", sol->f_q_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_q_m_nf", sol->f_q_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_g_p_nf", sol->f_g_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_f_g_m_nf", sol->f_g_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         SB = BASIS_PHYSICAL;
      }

      if (BASIS_DEFINITE_C_PAR == SB || BASIS_BOTH == SB) {
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_up_nf", sol->up_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_dn_nf", sol->dn_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_st_nf", sol->st_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         if (sol->nf > 3) single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_ch_nf", sol->ch_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         if (sol->nf > 4) single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_bm_nf", sol->bm_p, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_P_gl_nf", sol->gl_p, sol->grid_exponent, sol->G_t, sol->AG_t);

         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_up_nf", sol->up_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_dn_nf", sol->dn_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_st_nf", sol->st_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         if (sol->nf > 3) single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_ch_nf", sol->ch_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         if (sol->nf > 4) single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_bm_nf", sol->bm_m, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_M_gl_nf", sol->gl_m, sol->grid_exponent, sol->G_t, sol->AG_t);
      }

      if (BASIS_PHYSICAL == SB || BASIS_BOTH == SB) {
         single_save_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_up_nf", sol->T_up, sol->DT_up, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_dn_nf", sol->T_dn, sol->DT_dn, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         single_save_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_st_nf", sol->T_st, sol->DT_st, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         if (sol->nf > 3) single_save_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_ch_nf", sol->T_ch, sol->DT_ch, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
         if (sol->nf > 4) single_save_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_bm_nf", sol->T_bm, sol->DT_bm, _QUARK_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);

         single_save_physical(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_gl_nf", sol->Tp_gl, sol->Tm_gl, _GLUON_COMBINATION_TYPE_, sol->grid_exponent, sol->G_t, sol->AG_t);
      }
   }
   if (sol->chiral_odd_active) {
      single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_E_up_nf", sol->e_up, sol->grid_exponent, sol->G_t, sol->AG_t);
      single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_E_dn_nf", sol->e_dn, sol->grid_exponent, sol->G_t, sol->AG_t);
      single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_E_st_nf", sol->e_st, sol->grid_exponent, sol->G_t, sol->AG_t);
      single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_H_up_nf", sol->h_up, sol->grid_exponent, sol->G_t, sol->AG_t);
      single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_H_dn_nf", sol->h_dn, sol->grid_exponent, sol->G_t, sol->AG_t);
      single_save(sol->N, sol->M, sol->c_fact, sol->nf, sol->basefolder, basename, "_H_st_nf", sol->h_st, sol->grid_exponent, sol->G_t, sol->AG_t);
   }
}

double get_interpolated_value_TYPE(vector_t *distr, double x1, double x2, double x3, saved_solution_t *sol, double (*Fijk_t)(double, double, double, double, int32_t))
{

   static int32_t i = 0, j = 0;
   static double ix = 0, jy = 0;

   if (fabs(x1) <= sol->xmin && fabs(x2) <= sol->xmin && fabs(x3) <= sol->xmin) return DBL_MAX;

   get_3D_in_grid_space(x1, x2, x3, sol->N, sol->M, sol->c_fact, &ix, &jy, sol->grid_exponent, sol->G_t, sol->AG_t);


   i = (int32_t)round(ix);
   j = (int32_t)round(jy);

   if (i == sol->N + 1) {
      i = 0;
      ix -= sol->N + 1;
   }

   int32_t ip1 = i + 1, im1 = i - 1;
   if (im1 < 0) im1 += sol->N + 1;
   if (ip1 > sol->N) ip1 -= sol->N + 1;


   double res = 0;
   res += distr->data[from_ij_to_a(i, j, sol->M)] * Fijk_t(ix, jy, i, j, sol->N);
   res += distr->data[from_ij_to_a(ip1, j, sol->M)] * Fijk_t(ix, jy, i + 1, j, sol->N);
   res += distr->data[from_ij_to_a(im1, j, sol->M)] * Fijk_t(ix, jy, i - 1, j, sol->N);
   if (j < sol->M) {
      res += distr->data[from_ij_to_a(i, j + 1, sol->M)] * Fijk_t(ix, jy, i, j + 1, sol->N);
      if (Fijk_t == Fijk_minus) res += distr->data[from_ij_to_a(ip1, j + 1, sol->M)] * Fijk_t(ix, jy, i + 1, j + 1, sol->N);
      if (Fijk_t == Fijk_plus) res += distr->data[from_ij_to_a(im1, j + 1, sol->M)] * Fijk_t(ix, jy, i - 1, j + 1, sol->N);
   }
   if (j > 0) {
      res += distr->data[from_ij_to_a(i, j - 1, sol->M)] * Fijk_t(ix, jy, i, j - 1, sol->N);
      if (Fijk_t == Fijk_plus) res += distr->data[from_ij_to_a(ip1, j - 1, sol->M)] * Fijk_t(ix, jy, i + 1, j - 1, sol->N);
      if (Fijk_t == Fijk_minus) res += distr->data[from_ij_to_a(im1, j - 1, sol->M)] * Fijk_t(ix, jy, i - 1, j - 1, sol->N);
   }

   return res;
}

double get_interpolated_value(vector_t *distr, double x1, double x2, double x3, saved_solution_t *sol)
{
   switch (sol->F_t) {
   case IT_PLUS:
      return get_interpolated_value_TYPE(distr, x1, x2, x3, sol, Fijk_plus);
   case IT_MINUS:
      return get_interpolated_value_TYPE(distr, x1, x2, x3, sol, Fijk_minus);
   case IT_BOTH: {
      return (get_interpolated_value_TYPE(distr, x1, x2, x3, sol, Fijk_plus) + get_interpolated_value_TYPE(distr, x1, x2, x3, sol, Fijk_minus)) * 0.5;
   }
   }
   return DBL_MAX;
}

void get_3D_in_grid_space(double x1, double x2, double x3, int32_t N, int32_t M, double c_fact, double *phi, double *rho, double grid_exp, nonlinear_radial_grid_type_e G_type, angular_grid_type_e AG_type)
{
   *rho = max3(fabs(x1), fabs(x2), fabs(x3));
   int nT = (N + 1) / 6;

   switch (AG_type) {
   case LIN_GRID: {
      if (x1 > 0 && x2 >= 0 && x3 < 0) *phi = (x2 / (*rho)) * nT;
      if (x1 <= 0 && x2 > 0 && x3 < 0) *phi = (1 - x1 / (*rho)) * nT;
      if (x1 < 0 && x2 > 0 && x3 >= 0) *phi = (3 - x2 / (*rho)) * nT;
      if (x1 < 0 && x2 <= 0 && x3 > 0) *phi = (3 - x2 / (*rho)) * nT;
      if (x1 >= 0 && x2 < 0 && x3 > 0) *phi = (4 + x1 / (*rho)) * nT;
      if (x1 > 0 && x2 < 0 && x3 <= 0) *phi = (6 + x2 / (*rho)) * nT;
      break;
   }
   case COS_GRID: {
      if (x1 > 0 && x2 >= 0 && x3 < 0) *phi = nT * ONE_O_PI * acos(1.0 - 2.0 * (x2 / (*rho)));
      if (x1 <= 0 && x2 > 0 && x3 < 0) *phi = nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (-x1 / (*rho)));
      if (x1 < 0 && x2 > 0 && x3 >= 0) *phi = 2 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (1 - x2 / (*rho)));
      if (x1 < 0 && x2 <= 0 && x3 > 0) *phi = 3 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (-x2 / (*rho)));
      if (x1 >= 0 && x2 < 0 && x3 > 0) *phi = 4 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (x1 / (*rho)));
      if (x1 > 0 && x2 < 0 && x3 <= 0) *phi = 5 * nT + nT * ONE_O_PI * acos(1.0 - 2.0 * (1 + x2 / (*rho)));
      break;
   }
   default:
      break;
   }

   switch (G_type) {
   case LOG_GRID: {
      *rho = log((*rho)) * M / c_fact + M;
      break;
   }
   case PWR_GRID: {
      *rho = (pow((*rho), 1.0 / grid_exp) * (1 + c_fact) - c_fact) * M;
      break;
   }
   case IML_GRID: {
      *rho = M * (1.0 + c_fact * log((*rho) / (2.0 - (*rho))));
      break;
   }
   case HYP_GRID: {
      *rho = M * (1.0 + c_fact * acosh(1.0 / pow((*rho), 1.0 / grid_exp)));
      break;
   }
   default:
      break;
   }
}

#undef _QUARK_COMBINATION_TYPE_
#undef _GLUON_COMBINATION_TYPE_

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
