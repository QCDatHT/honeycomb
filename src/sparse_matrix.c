//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/sparse_matrix.h>
#include <TW3EV_include/thpool.h>

sparse_mat_t *init_sparse_matrix(int32_t size, int32_t n, int32_t m)
{
   sparse_mat_t *sp_m = (sparse_mat_t *)calloc(1, sizeof(sparse_mat_t));
   sp_m->size = size;
   sp_m->data = (double *)calloc(size, sizeof(double));
   sp_m->i_c = (int32_t *)calloc(size, sizeof(int32_t));
   sp_m->i_r = (int32_t *)calloc(size, sizeof(int32_t));
   sp_m->rows = n;
   sp_m->cols = m;
   sp_m->interval = (index_triplet_t *)calloc(n, sizeof(index_triplet_t));
   sp_m->n_intervals = n;
   return sp_m;
}

sparse_mat_t *init_sparse_matrix_non_zero(int32_t size, int32_t n, int32_t m)
{
   sparse_mat_t *sp_m = (sparse_mat_t *)calloc(1, sizeof(sparse_mat_t));
   sp_m->size = size;
   sp_m->data = (double *)malloc(size * sizeof(double));
   sp_m->i_c = (int32_t *)malloc(size * sizeof(int32_t));
   sp_m->i_r = (int32_t *)malloc(size * sizeof(int32_t));
   sp_m->rows = n;
   sp_m->cols = m;
   return sp_m;
}

void free_sparse_mat(sparse_mat_t **sp_m)
{
   if (sp_m == NULL || *sp_m == NULL) return;

   free((*sp_m)->i_c);
   free((*sp_m)->i_r);
   free((*sp_m)->data);
   free((*sp_m)->interval);
   free(*sp_m);
   *sp_m = NULL;
}

sparse_mat_t *from_matrix_to_sparse(double **in_mat, int32_t n, int32_t m)
{
   int32_t size = 0;
   for (int32_t i = 0; i < n; i++) {
      for (int32_t j = 0; j < m; j++) {
         if (fabs(in_mat[i][j]) > _50_LOCAL_DBL_EPSILON_) {
            size++;
         }
      }
   }
   sparse_mat_t *sp_m = init_sparse_matrix(size, n, m);
   int32_t index = 0;
   for (int32_t i = 0; i < n; i++) {
      for (int32_t j = 0; j < m; j++) {
         if (fabs(in_mat[i][j]) > _50_LOCAL_DBL_EPSILON_) {
            sp_m->i_r[index] = i;
            sp_m->i_c[index] = j;
            sp_m->data[index] = in_mat[i][j];
            index++;
         }
      }
   }

   return sp_m;
}

void save_sparse_mat(sparse_mat_t *sp_m, char fileName[], int32_t N, int32_t M)
{
   (void)N;
   if (sp_m == NULL) return;
   FILE *fp = fopen(fileName, "w");
   for (int32_t k = 0; k < sp_m->size; k++) {
      int32_t i, j, ip, jp;
      from_a_to_ij(sp_m->i_r[k], M, &i, &j);
      from_a_to_ij(sp_m->i_c[k], M, &ip, &jp);
      fprintf(fp, "%d\t%d\t%d\t%d\t%.16e\n", i, j, ip, jp, sp_m->data[k]);
   }
   fclose(fp);
}

void generate_intervals_sparse_matrix(sparse_mat_t *sp_m)
{

   sp_m->interval = realloc(sp_m->interval, sp_m->size * sizeof(index_triplet_t));
   if (sp_m->size == 0) {
      sp_m->n_intervals = 0;
      return;
   }

   int32_t start = 0;
   int32_t i0 = sp_m->i_r[0];
   int32_t curr = 0;
   sp_m->interval[curr].i = i0;
   sp_m->interval[curr].a = start;
   for (int32_t j = 0; j < sp_m->size; j++) {
      if (sp_m->i_r[j] == i0) continue;
      sp_m->interval[curr].b = j - 1;
      start = j;
      i0 = sp_m->i_r[j];
      curr++;
      sp_m->interval[curr].i = i0;
      sp_m->interval[curr].a = start;
   }
   sp_m->interval[curr].b = sp_m->size - 1;
   sp_m->n_intervals = curr + 1;
   sp_m->interval = realloc(sp_m->interval, sp_m->n_intervals * sizeof(index_triplet_t));
}

sparse_mat_t *read_sparse_mat(char fileName[], int32_t N, int32_t M)
{
   int32_t NN = (N + 1) * (M + 1);
   FILE *fp = fopen(fileName, "r");

   fseek(fp, 0L, SEEK_END);
   uint64_t file_size = ftell(fp);
   rewind(fp);

   char *global_buffer = (char *)calloc(file_size + 1, sizeof(char));
   size_t n_char = fread(global_buffer, 1, file_size, fp);
   fclose(fp);
   if (n_char != file_size) tw3ev_log(TW3EV_ERROR, "Panicking, read %ld bytes instead of expected %ld for %s.", n_char, file_size, fileName);

   char delim[2] = "\n";

   char *line = strtok(global_buffer, delim);
   int32_t size = NN;
   sparse_mat_t *sp_m = init_sparse_matrix(size, NN, NN);
   int32_t i = 0;
   while (NULL != line) {
      int32_t i0, j0, ip0, jp0;
      sscanf(line, "%d\t%d\t%d\t%d\t%le\n", &i0, &j0, &ip0, &jp0, &(sp_m->data[i]));

      sp_m->i_r[i] = from_ij_to_a(i0, j0, M);
      sp_m->i_c[i] = from_ij_to_a(ip0, jp0, M);

      i++;

      if (i == sp_m->size) {
         size += NN;
         resize_sparse_matrix(sp_m, size);
      }
      line = strtok(NULL, delim);
   }
   resize_sparse_matrix(sp_m, i);

   free(global_buffer);

   generate_intervals_sparse_matrix(sp_m);

   return sp_m;
}

sparse_mat_t *read_sparse_mat_old(char fileName[], int32_t N, int32_t M)
{
   int32_t NN = (N + 1) * (M + 1);
   char buffer[255];
   FILE *fp = fopen(fileName, "r");

   int32_t size = NN;
   sparse_mat_t *sp_m = init_sparse_matrix(size, NN, NN);
   int32_t i = 0;

   while (fgets(buffer, 255, fp) != NULL) {
      int32_t i0, j0, ip0, jp0;
      sscanf(buffer, "%d\t%d\t%d\t%d\t%le\n", &i0, &j0, &ip0, &jp0, &(sp_m->data[i]));

      sp_m->i_r[i] = from_ij_to_a(i0, j0, M);
      sp_m->i_c[i] = from_ij_to_a(ip0, jp0, M);

      i++;

      if (i == sp_m->size) {
         size += NN;
         resize_sparse_matrix(sp_m, size);
      }
   }
   fclose(fp);
   resize_sparse_matrix(sp_m, i);
   return sp_m;
}

void sparse_mat_vector_mul_parallel(sparse_mat_t *sp_m, vector_t *v_in, vector_t *v_out)
{
   memset(v_out->data, 0, sizeof(double) * v_out->size);
   const int32_t limit = sp_m->n_intervals;
   double *out = v_out->data;
   double *inp = v_in->data;
   double *data = sp_m->data;
   int32_t j = 0;
   index_triplet_t itr;
   int32_t out_index;

#ifdef USE_OMP
#pragma omp parallel for private(j, itr, out_index)
#endif
   for (int32_t i = 0; i < limit; i++) {
      itr = sp_m->interval[i];
      out_index = itr.i;
      for (j = itr.a; j <= itr.b; j++)
         out[out_index] += data[j] * inp[sp_m->i_c[j]];
   }
}

void sparse_mat_2_vector_mul_parallel(sparse_mat_t *sp_m, vector_t *v_in_1, vector_t *v_in_2, vector_t *v_out)
{
   memset(v_out->data, 0, sizeof(double) * v_out->size);
   const int32_t limit = sp_m->n_intervals;
   double *out = v_out->data;
   double *inp_1 = v_in_1->data;
   double *inp_2 = v_in_2->data;
   double *data = sp_m->data;
   int32_t *i_c = sp_m->i_c;

   int32_t j = 0;
   index_triplet_t itr;
   int32_t out_index;
#ifdef USE_OMP
#pragma omp parallel for private(j, itr, out_index)
#endif
   for (int32_t i = 0; i < limit; i++) {
      itr = sp_m->interval[i];
      out_index = itr.i;
      for (j = itr.a; j <= itr.b; j++) {
         out[out_index] += data[j] * (inp_1[i_c[j]] + inp_2[i_c[j]]);
      }
   }
}

void resize_sparse_matrix(sparse_mat_t *sp_m, int32_t new_size)
{
   sp_m->size = new_size;
   sp_m->data = (double *)realloc(sp_m->data, sizeof(double) * sp_m->size);
   sp_m->i_c = (int32_t *)realloc(sp_m->i_c, sizeof(int32_t) * sp_m->size);
   sp_m->i_r = (int32_t *)realloc(sp_m->i_r, sizeof(int32_t) * sp_m->size);
   sp_m->interval = (index_triplet_t *)realloc(sp_m->interval, sizeof(index_triplet_t) * sp_m->size);
}

vector_t *init_vector(int32_t n)
{
   vector_t *vt = (vector_t *)calloc(1, sizeof(vector_t));
   vt->size = n;
   vt->data = (double *)calloc(n, sizeof(double));

   return vt;
}

void copy_vector(vector_t *vin, vector_t *vout)
{
   if (vout == NULL || vin == NULL) return;
   if (vin->size != vout->size) return;
   memcpy(vout->data, vin->data, sizeof(double) * vin->size);
}

void sum_vector(vector_t *v1, vector_t *v2, vector_t *vout)
{
   if (vout == NULL) return;
   if (v1->size != v2->size || v1->size != vout->size || v2->size != vout->size) return;

#ifdef USE_OMP
#pragma omp parallel for
#endif
   for (int32_t i = 0; i < v1->size; i++)
      vout->data[i] = v1->data[i] + v2->data[i];
}

void dif_vector(vector_t *v1, vector_t *v2, vector_t *vout)
{
   if (vout == NULL) return;
   if (v1->size != v2->size || v1->size != vout->size || v2->size != vout->size) return;
#ifdef USE_OMP
#pragma omp parallel for
#endif
   for (int32_t i = 0; i < v1->size; i++)
      vout->data[i] = v1->data[i] - v2->data[i];
}

void weighted_sum_five_vectors(vector_t *v1, double c1, vector_t *v2, double c2, vector_t *v3, double c3, vector_t *v4, double c4, vector_t *v5, double c5, vector_t *v_out)
{
#ifdef USE_OMP
#pragma omp parallel for
#endif
   for (int i = 0; i < v_out->size; i++) {
      v_out->data[i] = v1->data[i] * c1 + v2->data[i] * c2 + v3->data[i] * c3 + v4->data[i] * c4 + v5->data[i] * c5;
   }
}

void scalar_vector_mul(double c, vector_t *v1, vector_t *vout)
{
#ifdef USE_OMP
#pragma omp parallel for
#endif
   for (int32_t i = 0; i < v1->size; i++)
      vout->data[i] = c * v1->data[i];
}

void free_vector(vector_t **v)
{
   if (v == NULL || *v == NULL) return;
   free((*v)->data);
   free(*v);
   *v = NULL;
}

void set_zero_vector(vector_t *v)
{
   if (v == NULL) return;

   memset(v->data, 0, sizeof(double) * v->size);
   return;
}

void scalar_sparse_mat_mul(double x, sparse_mat_t *sp_m)
{
#ifdef USE_OMP
#pragma omp parallel for
#endif
   for (int i = 0; i < sp_m->size; i++)
      sp_m->data[i] *= x;
}
void scalar_sparse_mat_add(double x, sparse_mat_t *sp_m)
{
#ifdef USE_OMP
#pragma omp parallel for
#endif
   for (int i = 0; i < sp_m->size; i++)
      sp_m->data[i] += x;
}

int32_t from_ij_to_a(int32_t i, int32_t j, int32_t M) { return j + (M + 1) * i; }
void from_a_to_ij(int32_t a, int32_t M, int32_t *i, int32_t *j)
{
   *j = a % (M + 1);
   *i = a / (M + 1);
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
