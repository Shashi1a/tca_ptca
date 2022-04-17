gfortran -I/usr/local/Cellar/open-mpi/4.1.2/include -Wl,-flat_namespace -Wl,-commons,use_dylibs -I/usr/local/Cellar/open-mpi/4.1.2/lib \
-L/usr/local/Cellar/open-mpi/4.1.2/lib -L/usr/local/opt/libevent/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -o a.out ptca_hubbard_centered_new.f90 -fallow-argument-mismatch -lblas -llapack


## to compile use this on mac
## mpirun --use-hwthread-cpus  a.out