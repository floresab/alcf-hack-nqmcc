gfortran base install => ./install.sh
cd build 
make

mpirun -n n_proc ./build/bin/qmc < inputs/be8.in

c13 for a real world example
