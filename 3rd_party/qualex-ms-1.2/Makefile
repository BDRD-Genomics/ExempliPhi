# Makefile of QUALEX-MS solver for GNU make

# Please substitute correct values for LAPACKLIB and BLASLIB
# to use it on your environment

CFLAGS = -DNDEBUG -O3 -funroll-all-loops -Wall
LINKFLAGS = -s

CC = gcc
CXX = g++
LIBS = -lgfortran
# LAPACKLIB = -L/export/anaconda3/envs/exempliphi/lib -llapack
LAPACKLIB = -llapack
#BLASLIB = -lf77blas -latlas
BLASLIB = -lblas

.SUFFIXES: .o .cc .c

OBJS_c = bool_vector.cc graph.cc greedy_clique.cc main.cc preproc_clique.cc qualex.cc refiner.cc eigen.c mdv.c
OBJS_o = bool_vector.o graph.o greedy_clique.o main.o preproc_clique.o qualex.o refiner.o eigen.o mdv.o

qualex-ms: $(OBJS_c) $(OBJS_o)
	$(CXX) $(LINKFLAGS) $(CFLAGS) -o $@ $(OBJS_o) $(LAPACKLIB) $(BLASLIB) $(LIBS) -Wl,-rpath=/share/apps/anaconda3/envs/amd-pipeline/lib

clean:
	rm *.o qualex-ms

.c.o:
	$(CC) $(CFLAGS) -c $<

.cc.o:
	$(CXX) $(CFLAGS) -c $<
