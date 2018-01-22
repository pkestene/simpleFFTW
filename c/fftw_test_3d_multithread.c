/*
 * gcc  -DUSE_FLOAT -o fftw_test_3d_mt_float  fftw_test_3d_multithread.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3f_threads -lfftw3f
 *
 * gcc              -o fftw_test_3d_mt_double fftw_test_3d_multithread.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3_threads -lfftw3
 */

#include <fftw3.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

/* USE_FLOAT allows to switch computations from double to float type */
#if defined(USE_FLOAT)
typedef float FFTW_REAL;
typedef fftwf_complex FFTW_COMPLEX;
typedef fftwf_plan    FFTW_PLAN;
#define FFTW_PLAN_DFT_R2C_2D fftwf_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D fftwf_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_R2C_3D fftwf_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D fftwf_plan_dft_c2r_3d
#define FFTW_DESTROY_PLAN fftwf_destroy_plan
#define FFTW_EXECUTE fftwf_execute
#define FFTW_EXPORT_WISDOM_TO_FILE fftwf_export_wisdom_to_file
#define FFTW_IMPORT_WISDOM_FROM_FILE fftwf_import_wisdom_from_file
#define FFTW_WISDOM_FILENAME ("fftwf_wisdom.txt")
#define FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define FFTW_INIT_THREADS fftwf_init_threads
#define FFTW_CLEANUP fftwf_cleanup
#define FFTW_CLEANUP_THREADS fftwf_cleanup_threads
#define FFTW_PRINT_PLAN fftwf_print_plan
#define FFTW_FLOPS fftwf_flops
#else
typedef double FFTW_REAL;
typedef fftw_complex FFTW_COMPLEX;
typedef fftw_plan    FFTW_PLAN;
#define FFTW_PLAN_DFT_R2C_2D fftw_plan_dft_r2c_2d
#define FFTW_PLAN_DFT_C2R_2D fftw_plan_dft_c2r_2d
#define FFTW_PLAN_DFT_R2C_3D fftw_plan_dft_r2c_3d
#define FFTW_PLAN_DFT_C2R_3D fftw_plan_dft_c2r_3d
#define FFTW_DESTROY_PLAN fftw_destroy_plan
#define FFTW_EXECUTE fftw_execute
#define FFTW_EXPORT_WISDOM_TO_FILE fftw_export_wisdom_to_file
#define FFTW_IMPORT_WISDOM_FROM_FILE fftw_import_wisdom_from_file
#define FFTW_WISDOM_FILENAME ("fftw_wisdom.txt")
#define FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define FFTW_INIT_THREADS fftw_init_threads
#define FFTW_CLEANUP fftw_cleanup
#define FFTW_CLEANUP_THREADS fftw_cleanup_threads
#define FFTW_PRINT_PLAN fftw_print_plan
#define FFTW_FLOPS fftw_flops
#endif

// FFTW_PATIENT, FFTW_MEASURE, FFTW_ESTIMATE
//#define O_METHOD  FFTW_PATIENT
//#define O_METHOD_STR "FFTW_PATIENT"
//#define O_METHOD  FFTW_MEASURE
//#define O_METHOD_STR "FFTW_MEASURE"
#define O_METHOD  FFTW_ESTIMATE
#define O_METHOD_STR "FFTW_ESTIMATE"

int main(int argc, char *argv[]) {

  FFTW_PLAN plan;
  double add, mul, fma, flops, deltaT;
  FFTW_REAL* in;
  FFTW_COMPLEX* out;

  struct timeval tv_start, tv_stop;

  int N_THREADS;
  int FFT_LEN;
  int N_ITER;
  int i;
  int dim_x, dim_y, dim_z;
  int dim_xyz, dim_xyz2;

  if (argc < 3) {
    printf("\n"
	   "Test/bench for multithread FFTW speed\n"
	   "Note: you must have compiled FFTW with './configure --enable-threads'\n"
	   "\n"
	   "fftw_test_3d_mt_double <l> <n> [<i>]\n"
	   " l : number of FFT points\n"
	   " n : number of SMP threads\n"
	   " i : how many FFT calls to make and average their exec times\n\n");
    return -1;
  }
  FFT_LEN = atoi(argv[1]);
  dim_x = FFT_LEN;
  dim_y = FFT_LEN;
  dim_z = FFT_LEN;
  dim_xyz  = dim_x * dim_y *  dim_z;
  dim_xyz2 = dim_x * dim_y * (dim_z/2+1);


  N_THREADS = atoi(argv[2]);
  N_ITER = 1;
  if (argc == 4) {
    N_ITER = atoi(argv[3]);
  }
  if (FFT_LEN < 8 || FFT_LEN > 1024) {
    printf("FFT len not between 8 and 1024\n"); return -1;
  }
  if (N_THREADS < 1) N_THREADS=1;
  if (N_THREADS > 16) N_THREADS=16;

  printf("Using FFT_LEN=%d, N_THREADS=%d and N_ITER=%d\n",
	 FFT_LEN,
	 N_THREADS,
	 N_ITER);

  FFTW_INIT_THREADS();
  FFTW_PLAN_WITH_NTHREADS(N_THREADS);

  in  = (FFTW_REAL *)    malloc(dim_xyz *sizeof(FFTW_REAL   ));
  out = (FFTW_COMPLEX *) malloc(dim_xyz2*sizeof(FFTW_COMPLEX));
  for (i = 0; i<dim_xyz; i++) {
    in[i] = drand48();
  }

  /* initialize the plan: this can take a long time! */
  //   printf("Allocated in=%p out=%p, FFT len=%d, threads=%d\n", in, out, FFT_LEN, N_THREADS);
  printf("Creating real-to-complex 3D DFT plan with " O_METHOD_STR "...\n");
  gettimeofday(&tv_start, NULL);
  plan = FFTW_PLAN_DFT_R2C_3D(dim_x, dim_y, dim_z, 
			      in, out, O_METHOD);
  gettimeofday(&tv_stop, NULL);
  deltaT =  tv_stop.tv_sec  - tv_start.tv_sec + 
    1e-6 * (tv_stop.tv_usec - tv_start.tv_usec);
  printf("Plan init time        : %f sec\n", deltaT);

#ifdef I_AM_AN_FFTW_GURU
  FFTW_PRINT_PLAN(plan);
  printf("\n");
#endif

  /* now run one or more FFT interations */
  gettimeofday(&tv_start, NULL);
  for (i=0; i<N_ITER; i++)
    FFTW_EXECUTE(plan);
  gettimeofday(&tv_stop, NULL);
  deltaT =  tv_stop.tv_sec  - tv_start.tv_sec +
    1e-6 * (tv_stop.tv_usec - tv_start.tv_usec);
  deltaT = deltaT / N_ITER;

  FFTW_FLOPS(plan, &add, &mul, &fma);
  flops = add + mul + 2*fma;
  printf("Plan exec time        : %f sec (average over %d iteration(s))\n", deltaT, N_ITER);
  printf("Theoretical FLOP count: %f\n", flops);
  printf("Plan MFLOPS estimate  : %f\n", 1e-6*flops/deltaT);

  free(in);
  free(out);

  FFTW_CLEANUP_THREADS();

  return 0;
}
