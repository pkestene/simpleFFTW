/*
 *
 * gcc -Wall -DUSE_FLOAT -o fftw_test_3d_serial_float  fftw_test_3d_serial.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3f
 *
 * gcc -Wall             -o fftw_test_3d_serial_double fftw_test_3d_serial.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3
 *
 */

#include <math.h>
#include <sys/time.h> // for gettimeofday
#include <stdlib.h> // for atoi


#include <fftw3.h>

#define SQR(x) ((x)*(x))

// fftw wrapper for single / double precision
#if defined(USE_FLOAT)
typedef float FFTW_REAL;
typedef fftwf_complex FFTW_COMPLEX;
typedef fftwf_plan    FFTW_PLAN;
#define FFTW_WISDOM_FILENAME         ("fftwf_wisdom.txt")
#define FFTW_PLAN_DFT_R2C_2D         fftwf_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D         fftwf_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_R2C_3D         fftwf_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D         fftwf_plan_dft_c2r_3d
#define FFTW_PLAN_DFT_3D             fftwf_plan_dft_3d
#define FFTW_DESTROY_PLAN            fftwf_destroy_plan
#define FFTW_EXECUTE                 fftwf_execute
#define FFTW_EXPORT_WISDOM_TO_FILE   fftwf_export_wisdom_to_file
#define FFTW_IMPORT_WISDOM_FROM_FILE fftwf_import_wisdom_from_file
#define FFTW_PLAN_WITH_NTHREADS      fftwf_plan_with_nthreads
#define FFTW_INIT_THREADS            fftwf_init_threads
#define FFTW_CLEANUP                 fftwf_cleanup
#define FFTW_CLEANUP_THREADS         fftwf_cleanup_threads
#define FFTW_PRINT_PLAN              fftwf_print_plan
#define FFTW_FLOPS                   fftwf_flops
#define FFTW_MALLOC                  fftwf_malloc
#define FFTW_FREE                    fftwf_free
#else
typedef double FFTW_REAL;
typedef fftw_complex FFTW_COMPLEX;
typedef fftw_plan    FFTW_PLAN;
#define FFTW_WISDOM_FILENAME         ("fftw_wisdom.txt")
#define FFTW_PLAN_DFT_R2C_2D         fftw_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D         fftw_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_R2C_3D         fftw_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D         fftw_plan_dft_c2r_3d
#define FFTW_PLAN_DFT_3D             fftw_plan_dft_3d
#define FFTW_DESTROY_PLAN            fftw_destroy_plan
#define FFTW_EXECUTE                 fftw_execute
#define FFTW_EXPORT_WISDOM_TO_FILE   fftw_export_wisdom_to_file
#define FFTW_IMPORT_WISDOM_FROM_FILE fftw_import_wisdom_from_file
#define FFTW_PLAN_WITH_NTHREADS      fftw_plan_with_nthreads
#define FFTW_INIT_THREADS            fftw_init_threads
#define FFTW_CLEANUP                 fftw_cleanup
#define FFTW_CLEANUP_THREADS         fftw_cleanup_threads
#define FFTW_PRINT_PLAN              fftw_print_plan
#define FFTW_FLOPS                   fftw_flops
#define FFTW_MALLOC                  fftw_malloc
#define FFTW_FREE                    fftw_free
#endif


#ifdef NO_FFTW_EXHAUSTIVE
#define MY_FFTW_FLAGS (FFTW_ESTIMATE)
#else
#define MY_FFTW_FLAGS (FFTW_EXHAUSTIVE)
#endif // NO_FFTW_EXHAUSTIVE

// FFTW_PATIENT, FFTW_MEASURE, FFTW_ESTIMATE
#define O_METHOD  FFTW_ESTIMATE
#define O_METHOD_STR "FFTW_ESTIMATE"


int main(int argc, char **argv){

  ptrdiff_t NX = 100, NY = 100, NZ = 100;

  // time measurement
  struct timeval tv_start, tv_stop;
  double deltaT;
  int N_ITER;

  //int k0=5;


  if (argc < 2) {
    printf("\n"
	   "Test/bench for serial FFTW speed\n"
	   "\n"
	   "\n"
	   "fftw_test_3d_serial_double <l> <n> [<i>]\n"
	   " l : number of FFT points\n"
	   " i : how many FFT calls to make and average their exec times\n\n");
    return -1;
  }
  NX = atoi(argv[1]);
  NY = NX;
  NZ = NX;
  printf("geometry : %ld %ld %ld\n",NX,NY,NZ);

  N_ITER = atoi(argv[2]);

  /*
   * 1. complex to complex transform
   */
  {
    FFTW_PLAN plan;
    FFTW_COMPLEX *data;
    ptrdiff_t i, j, k;

    data = (FFTW_COMPLEX *) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * NX*NY*NZ*sizeof(FFTW_COMPLEX));
    
    /* create plan for forward DFT */
    gettimeofday(&tv_start, NULL);
    plan = FFTW_PLAN_DFT_3D(NX, NY, NZ, data, data,
			    FFTW_FORWARD, FFTW_ESTIMATE);
    gettimeofday(&tv_stop, NULL);
    deltaT =  tv_stop.tv_sec  - tv_start.tv_sec +
      1e-6 * (tv_stop.tv_usec - tv_start.tv_usec);
    printf("Plan init time (c2c)  : %f sec\n", deltaT);
    
    /* initialize data */
    for (i = 0; i < NX; ++i) {
      for (j = 0; j < NY; ++j) {
	for (k = 0; k < NZ; ++k) {
	  /*data[(i*NY + j)*NZ + k][0] = 0;
	    data[(i*NY + j)*NZ + k][1] = k;*/
	  data[(i*NY + j)*NZ + k][0] = (i*NY + j)*NZ + k;
	  data[(i*NY + j)*NZ + k][1] = 0;
	  /*data[(i*NY + j)*NZ + k][0] = cos(2*M_PI*k0*(i)/NX);
	    data[(i*NY + j)*NZ + k][1] = sin(2*M_PI*k0*(i)/NX);*/
	  /*if (i<3 and j<3 and k<3) {
	    data[(i*NY + j)*NZ + k][0] = 1;
	    data[(i*NY + j)*NZ + k][1] = 0;
	    } else {
	    data[(i*NY + j)*NZ + k][0] = 0;
	    data[(i*NY + j)*NZ + k][1] = 0;
	    }*/
	}
      }
    }
    
    /* compute transforms, in-place, as many times as desired */
    gettimeofday(&tv_start, NULL);
    for (i=0; i<N_ITER; i++)
      FFTW_EXECUTE(plan);
    gettimeofday(&tv_stop, NULL);
    deltaT =  tv_stop.tv_sec  - tv_start.tv_sec +
      1e-6 * (tv_stop.tv_usec - tv_start.tv_usec);
    deltaT = deltaT / N_ITER;
    
    // print perf statistics
    double add, mul, fma, flops;
    FFTW_FLOPS(plan, &add, &mul, &fma);
    flops = add + mul + 2*fma;
    printf("Plan exec time (c2c)  : %f sec (average over %d iteration(s))\n", deltaT, N_ITER);
    printf("Theoretical FLOP count: %f\n", flops);
    printf("Plan MFLOPS estimate  : %f\n", 1e-6*flops/deltaT);
    
    FFTW_DESTROY_PLAN(plan);
    FFTW_FREE(data);

  } // end complex to complex transform

  printf("\n");
  
  /*
   * 2. complex to real transform
   */
  {
    FFTW_PLAN plan;
    FFTW_REAL *data;
    FFTW_COMPLEX *dataComplex;

    ptrdiff_t i, j, k;

    int NZ_r2c = 2*(NZ/2+1);

    // in-place transform
    data        = (FFTW_REAL *) FFTW_MALLOC(sizeof(FFTW_COMPLEX) * NX*NY*NZ_r2c*sizeof(FFTW_REAL));
    dataComplex = (FFTW_COMPLEX *) data;


    /* create plan for forward DFT */
    gettimeofday(&tv_start, NULL);
    plan = FFTW_PLAN_DFT_R2C_3D(NX, NY, NZ, data, dataComplex,
				FFTW_ESTIMATE);
    gettimeofday(&tv_stop, NULL);
    deltaT =  tv_stop.tv_sec  - tv_start.tv_sec +
      1e-6 * (tv_stop.tv_usec - tv_start.tv_usec);
    printf("Plan init time (r2c)  : %f sec\n", deltaT);
    
    /* initialize data */
    for (i = 0; i < NX; ++i) {
      for (j = 0; j < NY; ++j) {
	for (k = 0; k < NZ; ++k) {
	  data[(i*NY + j)*NZ_r2c + k] = (i * NY + j) * NZ + k;
	  /*data[(i*NY + j)*NZ_r2c + k] = cos(2*M_PI*k0*(i)/NX);*/
	  /*if (i<3 and j<3 and k<3) {
	    data[(i*NY + j)*NZ_r2c + k] = 1;
	    } else {
	    data[(i*NY + j)*NZ_r2c + k] = 0;
	    }*/
	}
      }
    }
    
    /* compute transforms, in-place, as many times as desired */
    gettimeofday(&tv_start, NULL);
    for (i=0; i<N_ITER; i++)
      FFTW_EXECUTE(plan);
    gettimeofday(&tv_stop, NULL);
    deltaT =  tv_stop.tv_sec  - tv_start.tv_sec +
      1e-6 * (tv_stop.tv_usec - tv_start.tv_usec);
    deltaT = deltaT / N_ITER;
    
    // print perf statistics
    double add, mul, fma, flops;
    FFTW_FLOPS(plan, &add, &mul, &fma);
    flops = add + mul + 2*fma;
    printf("Plan exec time (r2c)  : %f sec (average over %d iteration(s))\n", deltaT, N_ITER);
    printf("Theoretical FLOP count: %f\n", flops);
    printf("Plan MFLOPS estimate  : %f\n", 1e-6*flops/deltaT);
    
    FFTW_DESTROY_PLAN(plan);
    FFTW_FREE(data);

  } // end complex to complex transform

  


  return 0;
}
