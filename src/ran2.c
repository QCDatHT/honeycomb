//
// Author: Simone Rodini <mailto:simone.rodini@desy.de>
//

#include <TW3EV_include/ran2.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1 + IMM1 / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

double ran2(int64_t *idum)
{
   int j;
   int64_t k;
   static int64_t idum2 = 123456789L;
   static int64_t iy = 0;
   static int64_t iv[NTAB];
   double temp;

   if (*idum <= 0) {
      if (-(*idum) < 1)
         *idum = 1;
      else
         *idum = -(*idum);
      idum2 = (*idum);
      for (j = NTAB + 7; j >= 0; j--) {
         k = (*idum) / IQ1;
         *idum = IA1 * (*idum - k * IQ1) - k * IR1;
         if (*idum < 0)
            *idum += IM1;
         if (j < NTAB)
            iv[j] = *idum;
      }
      iy = iv[0];
   }
   k = (*idum) / IQ1;
   *idum = IA1 * (*idum - k * IQ1) - k * IR1;
   if (*idum < 0)
      *idum += IM1;
   k = idum2 / IQ2;
   idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
   if (idum2 < 0)
      idum2 += IM2;
   j = iy / NDIV;
   iy = iv[j] - idum2;
   iv[j] = *idum;
   if (iy < 1)
      iy += IMM1;
   if ((temp = AM * iy) > RNMX)
      return RNMX;
   else
      return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

double ran2N(int64_t *idum)
{
   static int iset = 0;
   static double gset;
   double fac, rsq, v1, v2, deviate;
   if (*idum < 0)
      iset = 0; 
   if (iset == 0) {
      do {
         v1 = 2.0 * ran2(idum) - 1.0;
         v2 = 2.0 * ran2(idum) - 1.0;
         rsq = v1 * v1 + v2 * v2;
      } while (rsq >= 1.0 || rsq == 0.0); 
      fac = sqrt(-2.0 * log(rsq) / rsq);

      gset = v1 * fac;
      iset = 1; 
      deviate = v2 * fac;
   } else {           
      iset = 0;      
      deviate = gset; 
   }

   return deviate;
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

