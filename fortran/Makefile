# ######## GNU ########
F90 = gfortran
CC  = gcc
FFLAGS = -O3 -Wall
#FFLAGS = -g -Wall
# -ffpe-trap=zero,invalid,underflow -fbacktrace
CFLAGS = -O3 -Wall
LDFLAGS = -O3

# ######## INTEL ########
# F90 = ifort
# FFLAGS = -warn all -O3 -xSSE4.2 -ipo
# #FFLAGS = -xSSE4.2 -ipo 
# #FFLAGS = -O3 -no-vec
# LDFLAGS = -O3

# ######## PGI ########
# F90 = pgfortran
# FFLAGS = -O4 -fast -Minform=warn -Minfo=all -Mpreprocess 
# ##-Mipa=inline -Mvect=simd:256
# ##FFLAGS = -g -Minform=warn -Minfo -Mpreprocess -Ktrap=fp
# LDFLAGS = -O4
# ##LDFLAGS = -g

# FFTW flags
FFLAGS +=  -I$(FFTW_INC_DIR)
LDFLAGS += -L$(FFTW_LIB_DIR) -lfftw3_threads -lfftw3

SRCDIR = .
OBJ = \
	m_Precision.o \
	wallclock.o \
	m_Timing.o \
	m_Monitoring.o \
	fftw_test_3d_multithread.o

fftw_test_3d_multithread: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LDFLAGS) 

clean:
	rm -f *.o *.mod fftw_test_3d_multithread


# fortran rule
%.o:    $(SRCDIR)/%.f90
	$(F90) -cpp $(FFLAGS) -c $<

# C rule
%.o:    $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $<

