//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef THPOOL_H
#define THPOOL_H

#include <TW3EV_include/default.h>

#define MAX_THREADS 64
#define MAX_QUEUE 65536

typedef struct thpool_t thpool_t;

thpool_t *thpool_create(int32_t thread_count, int32_t queue_size);

int32_t thpool_add(thpool_t *pool, void (*routine)(void *), void *arg);

int32_t thpool_destroy(thpool_t *pool);

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

