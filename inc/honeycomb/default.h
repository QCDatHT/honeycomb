//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DEFAULT_H
#define DEFAULT_H

#define EPS_TOLL 1.0e-12
#define _ZERO_THR_ 2.2204460492503131e-15
#define _LOCAL_DBL_EPSILON_ 2.2204460492503131e-16
#define _LOCAL_DBL_MIN_ 2.2250738585072014e-308
#define _1000_LOCAL_DBL_MIN_ 2.2250738585072014e-305
#define _100_LOCAL_DBL_EPSILON_ 2.2204460492503131e-14
#define _50_LOCAL_DBL_EPSILON_ 1.1102230246251565e-14

#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef USE_OMP
#include <omp.h>
#endif

#define WSL_SIZE 500
#define SQRT2 1.4142135623730950488016887242097
#define SQRT3 1.7320508075688772935274463415058
#define SQRT6 2.4494897427831780981972840747059
#define INV_SQRT2 0.7071067811865475
#define INV_SQRT6 0.4082482904638631
#define THREE_O_PI 0.9549296585513720146133025802351
#define PI_O_THREE 1.0471975511965977461542144610932

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define ONE_O_PI 0.318309886183790671537767526745

#define DA_INCREMENT 1024

#define BASEFOLDER_MAX_L 1023
#define FNAME_MAX_L 254

inline double max2(double a, double b) { return a >= b ? a : b; }
inline double min2(double a, double b) { return a < b ? a : b; }
inline double max3(double a, double b, double c) { return a >= b ? (a >= c ? a : c) : (b >= c ? b : c); }

#define Sin(x) sin((x))
#define Cos(x) cos((x))
#define Log(x) log((x))
#define Power(x, a) pow((x), (a))
#define Sqrt(x) sqrt((x))

#define _theta_(x) ((x) > 0 ? 1 : 0)
#define _Theta2_(x, y) (_theta_((x)) * _theta_((y)) - _theta_(-(x)) * _theta_(-(y)))

#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))

#define SWAP(a, b, T)                                                                                                                                                                                                                                    \
   {                                                                                                                                                                                                                                                     \
      T temp = a;                                                                                                                                                                                                                                        \
      a = b;                                                                                                                                                                                                                                             \
      b = temp;                                                                                                                                                                                                                                          \
   }

// Dynamic array implementation
#define DA_APPEND(str, item)                                                                                                                                                                                                                             \
   {                                                                                                                                                                                                                                                     \
      (str)->content[(str)->count] = (item);                                                                                                                                                                                                             \
      (str)->count++;                                                                                                                                                                                                                                    \
      if ((str)->count >= (str)->size) {                                                                                                                                                                                                                 \
         (str)->size += DA_INCREMENT;                                                                                                                                                                                                                    \
         (str)->content = realloc((str)->content, (str)->size * sizeof((str)->content[0]));                                                                                                                                                              \
      }                                                                                                                                                                                                                                                  \
   }

#define STR_APPEND_NULL(str) DA_APPEND(str, '\0')

#define DA_APPEND_MANY(str, items, how_many)                                                                                                                                                                                                             \
   {                                                                                                                                                                                                                                                     \
      bool _to_realloc_ = false;                                                                                                                                                                                                                         \
      while ((str)->count + (how_many) >= (str)->size) {                                                                                                                                                                                                 \
         (str)->size += DA_INCREMENT;                                                                                                                                                                                                                    \
         _to_realloc_ = true;                                                                                                                                                                                                                            \
      }                                                                                                                                                                                                                                                  \
      if (_to_realloc_) {                                                                                                                                                                                                                                \
         (str)->content = realloc((str)->content, (str)->size * sizeof((str)->content[0]));                                                                                                                                                              \
         printf("REALLOC\n");                                                                                                                                                                                                                            \
      }                                                                                                                                                                                                                                                  \
      memcpy((str)->content + (str)->count, items, how_many);                                                                                                                                                                                            \
      (str)->count += (how_many);                                                                                                                                                                                                                        \
   }

// for (uint64_t j = 0; j < (how_many); j++) {
//    (str)->content[(str)->count + j] = items[j];
// }
// Does not contemplate the possibility of calling it with already initialized content pointer!
#define DA_INIT(str, Tstr, T)                                                                                                                                                                                                                            \
   if (NULL == str) {                                                                                                                                                                                                                                    \
      (str) = (Tstr *)calloc(1, sizeof(Tstr));                                                                                                                                                                                                           \
      (str)->content = (T *)calloc(DA_INCREMENT, sizeof(T));                                                                                                                                                                                             \
      (str)->count = 0;                                                                                                                                                                                                                                  \
      (str)->size = DA_INCREMENT;                                                                                                                                                                                                                        \
   } else {                                                                                                                                                                                                                                              \
      (str)->content = (T *)calloc(DA_INCREMENT, sizeof(T));                                                                                                                                                                                             \
      (str)->count = 0;                                                                                                                                                                                                                                  \
      (str)->size = DA_INCREMENT;                                                                                                                                                                                                                        \
   }

// Magic to get a string that contains the name of the variable
#define GET_VAR_NAME(var) #var

typedef enum { TW3EV_INFO = 0, TW3EV_WARNING = 1, TW3EV_ERROR = 2 } tw3ev_log_level_e;

typedef enum { PL_NONE, PL_ESSENTIAL, PL_ALL } printout_level_e;

char *attach_path_to_string(char *s, const char *baseFolder, const char *subdir);
char *clear_string(char *s);
void tw3ev_log(tw3ev_log_level_e level, char *msg, ...);
void set_local_error_handling(void (*f)());


void tw3ev_quicksort(int32_t *arr, int32_t low, int32_t high);

inline double sign(const double x) { return (0.0 < x) - (x < 0.0); }
inline int32_t isign(const int32_t x) { return (0 < x) - (x < 0); }

void *free_ptr(void *p);

// inline functions are not linked, i.e. all translation units must know the implementation at compile time

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

