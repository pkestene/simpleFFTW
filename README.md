# simpleFFTW

A simple example of use for FFTW (serial, multithread and MPI) to compare the performance of a serial and parallel execution of the fast Fourier transform on a given platform.

We use both Complex-to-Complex as well as Real-to-Complex transformation for testing.

In all the examples, we assume environment variables FFTW_INC_DIR and FFTW_LIB_DIR
to specify the location of fftw headers and libraries.

```bash
# example config (if FFTW is installed by your package manager)
export FFTW_INC_DIR=/usr/include
export FFTW_LIB_DIR=/usr/lib/x86_64-linux-gnu
```

## C FFTW API

How to build ?

### serial

```bash
gcc -Wall -o fftw_test_3d_serial_double fftw_test_3d_serial.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3
```

### multithread

```bash
gcc -o fftw_test_3d_mt_double fftw_test_3d_multithread.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3_threads -lfftw3
```

### MPI

```bash
mpicc -o fftw_test_3d_mpi_double fftw_test_3d_mpi.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3_mpi -lfftw3
```

## Fortran FFTW API

Use the Makefile to build the multithreaded version.

```bash
make
```
