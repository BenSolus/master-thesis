SHELL = sh
CC    = gcc
CXX   = g++


GFLAGS   = -std=c99 -O0 -g -march=native -Wall -Wextra -pedantic -fPIC -Wno-overlength-strings -pipe -DUSE_OPENMP=1 -DUSE_OPENCL=1 -fopenmp
INCLUDES = -I../../H2Lib/Library -I./
BFLAGS   = -std=c++98 -O0 -g -march=native -Wall -Wextra -pedantic -fPIC -Wno-overlength-strings -pipe -DUSE_OPENMP=1 -DUSE_OPENCL=1 -fopenmp
LDFLAGS  = ../../H2Lib/libh2.a  -shared

LIBRARY   = libMT.so
BENCH     = benchmark
SOURCES   = gcopencl.c greencross.c laplace2d.c laplace3d.c ocl_system.c
CLSOURCES = clbasic.cl clgcidxinfo.cl clgeom.cl clgreencross.cl clquad.cl
OBJS      = $(SOURCES:.c=.o)
CLOBJS    = $(CLSOURCES:.cl=.o)

lib: $(OBJS)
	$(CC) -o $(LIBRARY) $(OBJS) $(GFLAGS) $(LDFLAGS)

$(OBJS): %.o: %.c
	$(CC) -c $< $(INCLUDES) $(GFLAGS) -o $@

bench:
	$(CXX) -o $(BENCH) ../Tests/benchmark.cpp gcopencl.o greencross.o laplace2d.o laplace3d.o ocl_system.o $(INCLUDES) $(BFLAGS) -lOpenCL -lblas -lcairo -llapack ../../H2Lib/libh2.a
clean:
	-rm -f ${LIBRARY} ${OBJS} ${BENCH}