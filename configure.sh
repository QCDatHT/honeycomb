#!/usr/bin/env bash

#
# Author: Simone Rodini <mailto:simone.rodini@desy.de>
#

PREFERRED_COMP=${1:-""}
LIBOMP_PATH=${LIBOMP_PATH:=""}

set -eu

# List of Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'

NC='\033[0m' # No Color

if [ -z "${PREFERRED_COMP}" ] || [ -z "$(which ${PREFERRED_COMP})" ]; then
   if [ -z "${PREFERRED_COMP}" ]; then
      echo -e "${BLUE}Empty compiler option${NC}"
   elif [ -z "$(which ${PREFERRED_COMP})" ]; then
      echo -e "${ORANGE}Asked compiler \"${PREFERRED_COMP}\" not found in the system${NC}"
   fi
   COMPILERS=("gcc-13" "gcc-12" "gcc-11" "gcc" "clang-16" "clang-15" "clang-14" "clang" "cc" "zig")
   echo -e "List of compilers that I am going to search: ${GREEN}${COMPILERS[*]}${NC}"
   PREFERRED_COMP=""
   for value in "${COMPILERS[@]}"; do
      if [ ! -z "$(which ${value})" ]; then
         echo -e "First supported compiler found: ${GREEN}${value}${NC}"
         PREFERRED_COMP="${value}"
         break
      fi
   done
fi
if [ -z "${PREFERRED_COMP}" ]; then
   echo -e "${RED}Panic:${NC} no compiler found"
   exit 1
else
   if [ "${PREFERRED_COMP}" = "zig" ]; then
      PREFERRED_COMP="zig cc"
   fi
   echo -e "I will build using: ${GREEN} ${PREFERRED_COMP}${NC}"
fi

if [ ! -d "./dat" ]; then mkdir ./dat; fi
if [ ! -d "./dat/kernels" ]; then mkdir ./dat/kernels; fi
if [ ! -d "./obj" ]; then mkdir ./obj; fi
if [ ! -d "./lib" ]; then mkdir ./lib; fi
if [ ! -d "./tests/bin" ]; then mkdir ./tests/bin; fi

if [ -f ./bin/libtw3ev.a ]; then
   rm ./bin/libtw3ev.a
fi

for file in ./obj/*.o; do
   if [ -f $file ]; then
      echo -e "${RED}Removing object file ${file}${NC}"
      rm $file
   fi
done

# get number of logical processors for default compilation

NFORMAKE=1
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
   NFORMAKE=$(nproc)
elif [[ "$OSTYPE" == "darwin"* ]]; then
   NFORMAKE=$(sysctl -n hw.ncpu)
else
   echo "Unsupported OS, aborting"
   exit 1
fi

JFLAG="-j"
if [ $NFORMAKE -eq 1 ]; then
   JFLAG="${JFLAG}${NFORMAKE}"
else
   TMP=$(($NFORMAKE - 1))
   JFLAG="${JFLAG}${TMP}"
fi

# Build variables for Makefile, prepend them
echo "" >Makefile
echo -e "SRCDIR=$(pwd)/src" >>Makefile
echo -e "TSTDIR=$(pwd)/tests" >>Makefile
echo -e "TSTBINDIR=$(pwd)/tests/bin" >>Makefile
echo -e "OBJDIR=$(pwd)/obj" >>Makefile
echo -e "INCDIR=$(pwd)/inc" >>Makefile
echo -e "LIBDIR=$(pwd)/lib\n" >>Makefile

# -ffast-math #=> unsafe
CFLAGS_COMMON='CFLAGS=-O3 -std=gnu17  -fPIC -I$(INCDIR) -march=native -Wall -Wextra -Werror -fno-math-errno -fno-trapping-math -fno-rounding-math -fno-signaling-nans -fexcess-precision=fast -fcx-limited-range -fomit-frame-pointer -funroll-loops'

if [ -f "temp.c" ]; then
   rm "temp.c"
fi
touch "temp.c"
echo '#include<stdio.h> ' >>temp.c
echo 'int main(){ ' >>temp.c
echo 'double res=0; ' >>temp.c
echo '#pragma omp parallel for ' >>temp.c
echo 'for(int i=0; i<1000; i++) ' >>temp.c
echo 'res += i; ' >>temp.c
echo 'return 0; ' >>temp.c
echo '}' >>temp.c

echo -e "${ORANGE}Checking for libomp support...${NC}"
set +e

tmp_FLAG=""
if [ ! -z "${LIBOMP_PATH}" ]; then
   export PATH="${LIBOMP_PATH}:${PATH}"
   tmp_FLAG="-L${LIBOMP_PATH}"
fi
$PREFERRED_COMP temp.c ${tmp_FLAG} -fopenmp -o compilation_check_tw3ev.out > /dev/null 2>&1

if [ $? -ne 0 ]; then
   echo -e "${RED}...with the selected compiler and library path the compilation fails, disabling OMP support${NC}"
else
   echo -e "${GREEN}...test compilation successful, enabling OMP support${NC}" 
   CFLAGS_COMMON="${CFLAGS_COMMON} -fopenmp -DUSE_OMP"
fi
set -e

rm -rf "temp.c"
rm -rf "compilation_check_tw3ev.out"


echo -e "CC=${PREFERRED_COMP}" >>Makefile
echo -e "${CFLAGS_COMMON}" >>Makefile

ARCHIVER="ar"
if [ "${PREFERRED_COMP}" = "zig cc" ]; then
   ARCHIVER="zig ar"
elif [ $(echo "$PREFERRED_COMP" | grep "gcc") ]; then
   oldIFS=$IFS
   IFS='-' read -ra TEMP_ARR <<<"$PREFERRED_COMP"

   length="${#TEMP_ARR[@]}"
   if [ "$length" -eq 1 ]; then
      ARCHIVER="${TEMP_ARR[0]}-ar"
   else
      ARCHIVER="${TEMP_ARR[0]}-ar-${TEMP_ARR[1]}"
   fi
   if [ -z "$(which ${ARCHIVER})" ]; then
      echo -e "${RED}${ARCHIVER} not found, defaulting to \`ar'${NC}"
      ARCHIVER="ar"
   fi
   IFS=$oldIFS
fi

echo -e "ARCHIV=${ARCHIVER}" >>Makefile

templ_Makefile=$(echo '
LDFLAGS = -ltw3ev -lm -lpthread 
LIBRARY=libtw3ev
CSRC=$(wildcard $(SRCDIR)/*.c)
HEADERS=$(wildcard $(INCDIR)/TW3EV_include/*.h)
OBJS=$(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(CSRC))

TCSRC=$(wildcard $(TSTDIR)/src/*.c)
TBINS=$(patsubst $(TSTDIR)/src/%.c, %, $(TCSRC))


all: $(LIBDIR)/$(LIBRARY) 
test: $(TBINS)

$(LIBDIR)/$(LIBRARY): $(OBJS)
	$(ARCHIV) r $(LIBDIR)/$(LIBRARY).a $(OBJS) 
	

$(OBJDIR)/%.o : $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@ 

$(TBINS): % : $(TSTDIR)/src/%.c
	$(CC) $(CFLAGS)  $^ -L$(LIBDIR)')

if [ ! -z "${LIBOMP_PATH}" ]; then
   templ_Makefile="${templ_Makefile} -L${LIBOMP_PATH}"
fi
templ_Makefile="${templ_Makefile} $(echo ' $(LDFLAGS) -o $(TSTBINDIR)/$@.out

clean:
	rm -rf  $(OBJDIR)/* $(LIBDIR)/*.a && clear')"

echo -e "${templ_Makefile}" >>Makefile

echo -e "#!/usr/bin/env bash" >compile.sh

if [ ! -z "${LIBOMP_PATH}" ]; then
   echo -e "export PATH=\"${LIBOMP_PATH}:\${PATH}\"" >>compile.sh
fi

echo "make ${JFLAG}" >>compile.sh
echo "make test ${JFLAG}" >>compile.sh

# echo -e "#!/usr/bin/env bash" >run_tests.sh
# echo -e "cd ./tests/bin && ./t_multiple.out && cd ../.." >>run_tests.sh

set +eu

# Copyright (C) 2024 Simone Rodini; Lorenzo Rossi
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.