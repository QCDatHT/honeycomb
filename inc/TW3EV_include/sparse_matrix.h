//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <TW3EV_include/default.h>
#include <TW3EV_include/thpool.h>

typedef struct {
   int32_t i, a, b;
} index_triplet_t;

typedef struct {
   int32_t *i_r, *i_c;
   double *data;
   int32_t size, rows, cols;
   index_triplet_t *interval;
   int32_t n_intervals;
} sparse_mat_t;

typedef struct {
   int32_t size;
   double *data;
} vector_t;

// The save-read routines in tensor notation are the one that should be used
// The others are for more generic matrices. They can be used but
// depends on the specific choice of the flattening function for the tensor
// indexes.


sparse_mat_t *from_matrix_to_sparse(double **in_mat, int32_t n, int32_t m);
sparse_mat_t *init_sparse_matrix(int32_t size, int32_t n, int32_t m);
sparse_mat_t *init_sparse_matrix_non_zero(int32_t size, int32_t n, int32_t m);
void free_sparse_mat(sparse_mat_t **);
// void save_sparse_mat(sparse_mat_t *sp_m, char fileName[]);
void save_sparse_mat(sparse_mat_t *sp_m, char fileName[], int32_t N, int32_t M);
// sparse_mat_t *read_sparse_mat(char fileName[], int32_t n, int32_t m);
sparse_mat_t *read_sparse_mat(char fileName[], int32_t N, int32_t M);
void sparse_mat_vector_mul(sparse_mat_t *sp_m, vector_t *v, vector_t *v_out);
void sparse_mat_2_vector_mul(sparse_mat_t *sp_m, vector_t *v1, vector_t *v2, vector_t *v_out);
void sparse_mat_vector_mul_parallel(sparse_mat_t *sp_m, vector_t *v, vector_t *v_out);
void sparse_mat_2_vector_mul_parallel(sparse_mat_t *sp_m, vector_t *v1, vector_t *v2, vector_t *v_out);
void resize_sparse_matrix(sparse_mat_t *sp_m, int32_t new_size);

void generate_intervals_sparse_matrix(sparse_mat_t *sp_m);

void scalar_sparse_mat_mul(double x, sparse_mat_t *sp_m);
void scalar_sparse_mat_add(double x, sparse_mat_t *sp_m);

// vector algebra
vector_t *init_vector(int32_t n);
void free_vector(vector_t **v);
void set_zero_vector(vector_t *v);
// The vector sum can work in-place, i.e. either v1 or v2 or both can be equal
// to vout
void copy_vector(vector_t *vin, vector_t *vout);
void sum_vector(vector_t *v1, vector_t *v2, vector_t *vout);
void dif_vector(vector_t *v1, vector_t *v2, vector_t *vout);
void weighted_sum_five_vectors(vector_t *v1, double c1, vector_t *v2, double c2, vector_t *v3, double c3, vector_t *v4, double c4, vector_t *v5, double c5, vector_t *v_out);
void scalar_vector_mul(double c, vector_t *v1, vector_t *vout);

// some util functions
int32_t from_ij_to_a(int32_t i, int32_t j, int32_t M);
void from_a_to_ij(int32_t a, int32_t M, int32_t *i, int32_t *j);

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

