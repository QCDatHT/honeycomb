
SRCDIR=/mnt/d/Documents/Code/tw3evolution/src
TSTDIR=/mnt/d/Documents/Code/tw3evolution/tests
TSTBINDIR=/mnt/d/Documents/Code/tw3evolution/tests/bin
OBJDIR=/mnt/d/Documents/Code/tw3evolution/obj
INCDIR=/mnt/d/Documents/Code/tw3evolution/inc
LIBDIR=/mnt/d/Documents/Code/tw3evolution/lib

CC=gcc-13
CFLAGS=-O3 -std=gnu17  -fPIC -I$(INCDIR) -march=native -Wall -Wextra -Werror -fno-math-errno -fno-trapping-math -ffinite-math-only -fno-rounding-math -fno-signaling-nans -fexcess-precision=fast -fcx-limited-range -fomit-frame-pointer -funroll-loops -fopenmp -DUSE_OMP
ARCHIV=gcc-ar-13

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
	$(CC) $(CFLAGS)  $^ -L$(LIBDIR)  $(LDFLAGS) -o $(TSTBINDIR)/$@.out

clean:
	rm -rf  $(OBJDIR)/* $(LIBDIR)/*.a && clear
