% Developers are Prof. Haavard Rue and Dr. Alexander Litvinenko
% CEMSE, KAUST
Interface between INLA package and sparse matrix solver PARDISO

to compile pardiso examples
~/PARDISO> gcc mypardiso_sym.c -g -o mypardiso_sym.exe -L/home/litvina/PARDISO libpardiso500-GNU461-X86-64.so -llapack -lblas  -lgfortran -fopenmp -lpthread -lm



cd pardiso_6/

vi f150_p1.txt

vi many_runs.sh

vi many_runs7.sh

cd /project/

cd k1051

cd pardiso_6/


vi f150_p1.txt

module unload PrgEnv-cray/6.0.4

module load PrgEnv-intel/6.0.4

vi many_runs7.sh

export OMP_NUM_THREADS=32
./mypardiso_unsym_olaf.exe -nx=150 -ny=150 -nz=1 -nt=150 -mm=20 > f150_p32.txt

export OMP_NUM_THREADS=32

./mypardiso_unsym_olaf.exe -nx=150 -ny=150 -nz=1 -nt=150 -mm=20 > f150_p32.txt

vi f150_p32.txt

./mypardiso_unsym_olaf.exe -nx=150 -ny=150 -nz=1 -nt=150 -mm=20 

vi mypardiso_unsym_olaf.c

vi many_runs7.sh

vi shaheen-pardiso.sh

vi output.txt

vi error.txt

sbatch shaheen-pardiso.sh
