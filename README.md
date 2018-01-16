% Developers are Prof. Haavard Rue and Dr. Alexander Litvinenko
% CEMSE, KAUST
Interface between INLA package and sparse matrix solver PARDISO

to compile pardiso examples
~/PARDISO> gcc mypardiso_sym.c -g -o mypardiso_sym.exe -L/home/litvina/PARDISO libpardiso500-GNU461-X86-64.so -llapack -lblas  -lgfortran -fopenmp -lpthread -lm
