//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/default.h>
#include <pthread.h>

pthread_mutex_t log_lock = PTHREAD_MUTEX_INITIALIZER;

void local_error_handling();

void (*tw3ev_error_handling)() = local_error_handling;

void local_error_handling()
{
   fprintf(stderr, "DEFAULT ERROR HANDLER INVOKED: EXITING THE PROGRAM!\n");
   exit(1);
}

void set_local_error_handling(void (*f)()) { tw3ev_error_handling = f; }

void tw3ev_log(tw3ev_log_level_e level, char *msg, ...)
{
   static char const messages[3][255] = {"[INFO]", "[WARNING]", "[ERROR]"};
   static char const padding[3][255] = {"   ", "", "  "};
   static char const colors[3][255] = {"\033[0;32m", "\033[0;33m", "\033[0;31m"};
   static char buffer[1024] = "";
   va_list argptr;
   va_start(argptr, msg);
   vsprintf(buffer, msg, argptr);
   va_end(argptr);

   pthread_mutex_lock(&log_lock);

   fprintf(stderr, "%s", colors[level]);
   fprintf(stderr, "%s%s ", messages[level], padding[level]);
   fprintf(stderr, "%s", "\033[0m"); 
   fprintf(stderr, "%s\n", buffer);
   if (level == TW3EV_ERROR) tw3ev_error_handling();
   pthread_mutex_unlock(&log_lock);
}

char *attach_path_to_string(char *s, const char *baseFolder, const char *subdir)
{
   char buffer[BASEFOLDER_MAX_L + 1] = {0};
   strncpy(buffer, s, BASEFOLDER_MAX_L);
   strcpy(s, baseFolder);

   if ('\0' == s[0]) {
      s[0] = '.';
      s[1] = '/';
      s[2] = '\0';
   } else {
      size_t len = strlen(s);
      int64_t i = len - 1;
      bool skipping = false;
      for (; i >= 0; i--) {
         if ('/' == s[i] && !skipping) {
            skipping = true;
            continue;
         }
         if ('/' == s[i] && skipping) {
            for (int64_t j = i; j < (int64_t)len; j++) {
               s[j] = s[j + 1]; 
            }
            len--;   
            continue;
         }
         if ('/' != s[i]) skipping = false;
      }
      if ('/' != s[len - 1]) s[len++] = '/';
   }
   strcat(s, subdir);
   strcat(s, buffer);
   return s;
}

char *clear_string(char *s)
{
   strcpy(s, "");
   return s;
}

void *free_ptr(void *p)
{
   if (NULL == p) return NULL;
   free(p);
   return NULL;
}

static void partition(int32_t *arr, int32_t low, int32_t high, int32_t *lt_out, int32_t *gt_out)
{
   int32_t pivot = arr[(high + low) / 2];
   int32_t lt = low, eq = low, gt = high;
   while (eq <= gt) {
      if (arr[eq] < pivot) {
         SWAP(arr[eq], arr[lt], int32_t);
         lt++;
         eq++;
      } else if (arr[eq] > pivot) {
         SWAP(arr[eq], arr[gt], int32_t);
         gt--;
      } else eq += 1;
   }
   *lt_out = lt;
   *gt_out = gt;
   return;
}

void tw3ev_quicksort(int32_t *arr, int32_t low, int32_t high)
{
   if (low >= 0 && low < high) {
      int32_t lt, gt;
      partition(arr, low, high, &lt, &gt);
      tw3ev_quicksort(arr, low, lt - 1);
      tw3ev_quicksort(arr, gt + 1, high);
   }
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

