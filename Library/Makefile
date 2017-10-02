SHELL = sh
CC    = gcc-4.6
CXX   = g++-4.6


GFLAGS   = -std=c99 -O2 -g -march=native -Wall -Wextra -pedantic -fPIC -Wno-overlength-strings -pipe -mavx -DUSE_FLOAT=1 -DUSE_OPENMP=1 -DUSE_OPENCL=1 -fopenmp
INCLUDES = -I../../H2Lib/Library -I./
BFLAGS   = -std=c++98 -O2 -g -march=native -Wall -Wextra -pedantic -fPIC -Wno-overlength-strings -pipe -mavx -DUSE_FLOAT=1 -DUSE_OPENMP=1 -DUSE_OPENCL=1 -fopenmp
LDFLAGS  = ../../H2Lib/libh2.a  -shared

LIBRARY   = libMT.so
BENCH     = benchmark
SOURCES   = fastaddevalgca.c gcopencl.c greencross.c laplace2d.c laplace3dgca.c ocl_system.c oclintegralinfos.c oclworkpkgs.c singquadgca.c
CLSOURCES = clbasic.cl clfastaddevalgca.cl clgcidxinfo.cl clgeom.cl clgreencross.cl clintegralinfo.cl clkernels.cl clsingquad.cl
OBJS      = $(SOURCES:.c=.o)
CLOBJS    = $(CLSOURCES:.cl=.o)

lib: $(OBJS)
	$(CC) -o $(LIBRARY) $(OBJS) $(GFLAGS) $(LDFLAGS)

$(OBJS): %.o: %.c
	$(CC) -c $< $(INCLUDES) $(GFLAGS) -o $@

bench:
	$(CXX) -o $(BENCH) ../Tests/benchmark.cpp fastaddevalgca.o gcopencl.o greencross.o laplace2d.o laplace3dgca.o ocl_system.o oclintegralinfos.o oclworkpkgs.o singquadgca.o $(INCLUDES) $(BFLAGS) -lOpenCL -lblas -lcairo -llapack ../../H2Lib/libh2.a
clean:
	-rm -f ${LIBRARY} ${OBJS} ${BENCH}