## compiler
FC = gfortran

# compiler flags
# standard
CFLAGS = -std=f2008ts

# warning flags
CFLAGS += -Wall 

## lib flags
LFLAGS = -Wl,  -flat_namespace -Wl,-commons, use_dylibs -I/usr/local/Cellar/open-mpi/4.1.2/lib -L/usr/local/Cellar/open-mpi/4.1.2/lib \
		-L/usr/local/opt/libevent/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lblas -llapack -fallow-argument-mismatch

## source files
SRCS = ptca_hubbard_centered_new
OBJS = $(SRCS:=.o)

## executable
MAIN = main 

## compile project
all :$(MAIN)
	@echo Model compiled

$(MAIN) : $(OBJS)
	$(FC) $(CFLAGS) -o $(MAIN) %(OBJS) $(LFLAGS)

.SUFFIXES :.o .f90

.f90.o :
	$(FC) $(CFLAGS) $(LFLAGS) -c $<

clean :
	$(RM) *.o *.so *.mod $(MATH)