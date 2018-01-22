/*
 * mpicc -DUSE_FLOAT -o fftw_test_3d_mpi_float  fftw_test_3d_mpi.c -I$FFTW_INC_DIR -L$FFTW_LIB_DIR -lfftw3f_mpi -lfftw3f
 *
 * mpicc             -o fftw_test_3d_mpi_double fftw_test_3d_mpi.c -I$FFTW_INC_DIR  -L$FFTW_LIB_DIR -lfftw3_mpi -lfftw3
 *
 * Run poincare frontal:
 * mpirun --mca btl sm,self --mca mtl ^psm -np 16 ./fftw_test_3d_mpi_double 256 16
 *
 */


#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include <fftw3-mpi.h>

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
#define FFTW_MPI_PLAN_DFT_R2C_3D     fftwf_mpi_plan_dft_r2c_3d
#define FFTW_DESTROY_PLAN            fftwf_destroy_plan
#define FFTW_EXECUTE                 fftwf_execute
#define FFTW_EXPORT_WISDOM_TO_FILE   fftwf_export_wisdom_to_file
#define FFTW_IMPORT_WISDOM_FROM_FILE fftwf_import_wisdom_from_file
#define FFTW_PLAN_WITH_NTHREADS      fftwf_plan_with_nthreads
#define FFTW_INIT_THREADS            fftwf_init_threads
#define FFTW_MPI_INIT                fftwf_mpi_init
#define FFTW_CLEANUP                 fftwf_cleanup
#define FFTW_CLEANUP_THREADS         fftwf_cleanup_threads
#define FFTW_PRINT_PLAN              fftwf_print_plan
#define FFTW_FLOPS                   fftwf_flops
#define FFTW_MALLOC                  fftwf_malloc
#define FFTW_FREE                    fftwf_free
#define FFTW_MPI_LOCAL_SIZE_3D       fftwf_mpi_local_size_3d
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
#define FFTW_MPI_PLAN_DFT_R2C_3D     fftw_mpi_plan_dft_r2c_3d
#define FFTW_DESTROY_PLAN            fftw_destroy_plan
#define FFTW_EXECUTE                 fftw_execute
#define FFTW_EXPORT_WISDOM_TO_FILE   fftw_export_wisdom_to_file
#define FFTW_IMPORT_WISDOM_FROM_FILE fftw_import_wisdom_from_file
#define FFTW_PLAN_WITH_NTHREADS      fftw_plan_with_nthreads
#define FFTW_INIT_THREADS            fftw_init_threads
#define FFTW_MPI_INIT                fftw_mpi_init
#define FFTW_CLEANUP                 fftw_cleanup
#define FFTW_CLEANUP_THREADS         fftw_cleanup_threads
#define FFTW_PRINT_PLAN              fftw_print_plan
#define FFTW_FLOPS                   fftw_flops
#define FFTW_MALLOC                  fftw_malloc
#define FFTW_FREE                    fftw_free
#define FFTW_MPI_LOCAL_SIZE_3D       fftw_mpi_local_size_3d
#endif


#ifdef NO_FFTW_EXHAUSTIVE
#define MY_FFTW_FLAGS (FFTW_ESTIMATE)
#else
#define MY_FFTW_FLAGS (FFTW_EXHAUSTIVE)
#endif // NO_FFTW_EXHAUSTIVE

// FFTW_PATIENT, FFTW_MEASURE, FFTW_ESTIMATE
#define O_METHOD  FFTW_ESTIMATE
#define O_METHOD_STR "FFTW_ESTIMATE"

int main(int argc, char *argv[]) 
{

  int NX,NY,NZ, NZ_r2c;
  int N_ITER;

  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int nbMpiProc;
  MPI_Comm_size(MPI_COMM_WORLD, &nbMpiProc);

  if (rank == 0) {
    // parse command line
    if (argc < 3) {
      printf("\n"
	     "Test/bench for mpi FFTW (C2C and R2C, in-place)\n"
	     "\n"
	     "\n"
	     "fftw_test_3d_mpi_double <l> [<i>]\n"
	     " l : number of FFT points\n"
	     " i : how many FFT calls to make and average their exec times\n\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    NX = atoi(argv[1]);

    printf("geometry : %d %d %d\n",NX,NX,NX);
    
    N_ITER = atoi(argv[2]);    
  }

  // broadcast NX and N_ITER
  MPI_Bcast(&NX    ,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&N_ITER,1,MPI_INT,0,MPI_COMM_WORLD);
  NY = NX;
  NZ = NX;
  NZ_r2c = 2*(NZ/2+1);
  
  // initialize FFTW in MPI
  FFTW_MPI_INIT();

  /*
   * 2. real to complex transform
   */
  {
    FFTW_REAL   * in; // local data
    FFTW_COMPLEX* out;
    
    ptrdiff_t  i, j, k; // for    row major order
    
    // fftw resources
    FFTW_PLAN  plan;
    ptrdiff_t alloc_local, local_n0, local_0_start;
    
    // time measurement
    double tStart, tStop, deltaT;

    // get local data size and allocate
    alloc_local = FFTW_MPI_LOCAL_SIZE_3D(NX, NY, NZ/2+1, MPI_COMM_WORLD,
					 &local_n0, &local_0_start);
    in  = (FFTW_REAL    *) FFTW_MALLOC( sizeof(FFTW_REAL) * 2 * alloc_local);
    out = (FFTW_COMPLEX *) in;
    printf("[rank %d] alloc_local %ld local_n0 %ld local_0_start %ld\n", 
	   rank, alloc_local, local_n0, local_0_start);


    for (i = 0; i < local_n0; i++) {
      for (j = 0; j < NY; j++) {
	for (k = 0; k < NZ; k++) {
	  in[(i*NY + j)*NZ_r2c + k] = drand48();
	}
      }
    }
    
    /* initialize the plan with global sizes: this can take a long time! */
    if (rank==0) printf("Creating MPI real-to-complex 3D inplace DFT plan with " O_METHOD_STR "...\n");
    tStart = MPI_Wtime();
    plan = FFTW_MPI_PLAN_DFT_R2C_3D(NX, NY, NZ, 
				    in, out, 
				    MPI_COMM_WORLD, O_METHOD);
    tStop = MPI_Wtime();
    deltaT =  tStop - tStart;
    if (rank==0) printf("Plan init time        : %f sec\n", deltaT);
    
#ifdef I_AM_AN_FFTW_GURU
    if (rank==0) {
      FFTW_PRINT_PLAN(plan);
      printf("\n");
    }
#endif

    /* now run one or more FFT interations */
    tStart = MPI_Wtime();
    for (i=0; i<N_ITER; i++)
      FFTW_EXECUTE(plan);
    tStop = MPI_Wtime();
    deltaT =  (tStop - tStart)/N_ITER;
    
    /* collect and print statistics */
    double add, mul, fma, flops;
    
    FFTW_FLOPS(plan, &add, &mul, &fma);
    flops = add + mul + 2*fma;
    
    if (rank==0) {
      printf("Plan exec time (r2c)  : %f sec (average over %d iteration(s))\n", deltaT, N_ITER);
      printf("Theoretical FLOP count: %f\n", flops);
      printf("Plan MFLOPS estimate  : %f\n", 1e-6*flops/deltaT);
    }
    
    FFTW_DESTROY_PLAN(plan);
    FFTW_FREE(in);

  }

  MPI_Finalize();

  return 0;
}
