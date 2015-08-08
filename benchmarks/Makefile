CLANG_PREFIX=/usr
CXX = $(CLANG_PREFIX)/bin/clang++
CC = $(CLANG_PREFIX)/bin/clang
G = g++
TERRA = ../../terra/terra

FLAGS = -O3
LINUX = -pthread

#FLAGS = -g -O3 -I/usr/local/atlasgccavx/include
#FLAGS = -O3 -I/opt/intel/composer_xe_2013.0.088/mkl/include
#LIBS = -framework Accelerate
# LIBS = -L/usr/local/atlasgccavx/lib -latlas -lcblas
#LIBS = -L/opt/intel/composer_xe_2013.0.088/mkl/lib -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

#for linux, flag -pthread is necessary
send: 
	scp convolution.c mnovaes@tucson:~/cnn/benchmarks

all: install run 

install:
	$(CC) $(FLAGS) -o conv convolution.c ../bin/my_numconv.o #../bin/my_naivenumconv.o

run:
	./conv

clean:
	rm ../bin/my_conv.o
	# rm *o

server-linux: 
	$(G) $(FLAGS) $(LINUX) -o conv convolution.c ../bin/my_numconv.o

mkl:
	icc -w -vec-report0 -I/opt/intel/system_studio_2015.3.055/mkl/include convolution.c ../bin/my_numconv.o ../bin/my_naivenumconv.o -L"/opt/intel/system_studio_2015.3.055/mkl/lib/intel64" -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L"/opt/intel/system_studio_2015.3.055/compiler/lib/intel64" -liomp5 -lm -ldl -lpthread  -o conv.out

%:	%.cpp
	$(CXX) $(FLAGS) $(LIBS) -o $@ $<

%:	%.c
	$(CC) $(FLAGS) $(LIBS) -o $@ $<
