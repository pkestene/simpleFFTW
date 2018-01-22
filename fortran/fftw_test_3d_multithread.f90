!!!! -*- Mode: F90 -*- !!!!


! FFTW_PATIENT, FFTW_MEASURE, FFTW_ESTIMATE
! #define O_METHOD  FFTW_PATIENT
! #define O_METHOD_STR "FFTW_PATIENT"
! #define O_METHOD  FFTW_MEASURE
! #define O_METHOD_STR "FFTW_MEASURE"
#define O_METHOD  FFTW_ESTIMATE
#define O_METHOD_STR 'FFTW_ESTIMATE'

program fftw_test_3d_multithread

  use Precision
  use Monitoring

  include 'fftw3.f'

  integer*8          :: plan
  integer            :: iRet
  real   (fp_kind)   :: add, mul, fma, flops, deltaT
  character(LEN=100) :: option

  integer(int_kind)  :: N_THREADS
  integer(int_kind)  :: FFT_LEN
  integer(int_kind)  :: N_ITER
  integer(int_kind)  :: i
  integer(int_kind)  :: dim_x, dim_y, dim_z

  real(fp_kind)   , dimension(:,:,:), allocatable :: in
  complex(fp_kind), dimension(:,:,:), allocatable :: out


  narg = COMMAND_ARGUMENT_COUNT()
  if (narg < 2) then
     write(*,*) ""
     write(*,*) "Test/bench for multithread FFTW speed"
     write(*,*) "Note: you must have compiled FFTW with './configure --enable-threads'"
     write(*,*) ""
     write(*,*) "fftw_test_3d_mt_double <l> <n> [<i>]"
     write(*,*) " l : number of FFT points"
     write(*,*) " n : number of SMP threads"
     write(*,*) " i : how many FFT calls to make and average their exec times"
     stop
  end if 

  ! retrieve 1st input parameter : FFT size
  call GET_COMMAND_ARGUMENT(1,option)
  read(option,*)  FFT_LEN

  dim_x    = FFT_LEN
  dim_y    = FFT_LEN
  dim_z    = FFT_LEN

  ! retrieve 2nd input parameter : number of threads
  call GET_COMMAND_ARGUMENT(2,option)
  read(option,*) N_THREADS

  ! retrieve 3rd input parameter : number of iteration
  N_ITER = 1
  if (narg>=3) then
     call GET_COMMAND_ARGUMENT(3,option)
     read(option,*) N_ITER
  end if

  if (FFT_LEN < 8 .or. FFT_LEN > 1024) then
     write(*,*) "FFT len not between 8 and 1024"
     stop
  end if
  
  if (N_THREADS < 1) N_THREADS=1
  if (N_THREADS > 16) N_THREADS=16

  write(*,*) "Using ", FFT_LEN, " as fft length"
  write(*,*) "Using ", N_THREADS, " threads"
  write(*,*) "Using ", N_ITER, " iteration"

  ! allocate memory and initialize
  allocate (in  (     dim_x     , dim_y , dim_z) )
  allocate (out (     dim_x/2+1 , dim_y , dim_z) )

  ! 
  DO k = 1,dim_z
     DO j = 1,dim_y
        DO i = 1,dim_x
           in(i,j,k) = 1.d0 * (i+j+k)
        END DO
     END DO
  END DO

  ! initialize FFTW
  call dfftw_init_threads(iRet)
  write(*,*) 'fftw init thread status : ',iRet

  call dfftw_plan_with_nthreads(N_THREADS)

  !   initialize the plan: this can take a long time !
  
  write(*,*) 'Creating real-to-complex 3D DFT plan with ' // O_METHOD_STR // ' ...'
  call timerStart(fft_plan_timer)
  call dfftw_plan_dft_r2c_3d(plan, dim_x, dim_y, dim_z, in, out, FFTW_ESTIMATE)
  call timerStop(fft_plan_timer)

  write(*,*) 'Plan init time        : ', fft_plan_timer%elapsed , ' sec'

#ifdef I_AM_AN_FFTW_GURU
   dfftw_print_plan(plan)
#endif

   !  now run one or more FFT interations 
   call timerStart(fft_exec_timer)

   DO i = 1,N_ITER
      call dfftw_execute_dft_r2c(plan, in, out)
   end DO

   call timerStop (fft_exec_timer)
   deltaT = fft_exec_timer%elapsed / N_ITER

   
   call dfftw_flops(plan, add, mul, fma)
   flops = add + mul + 2*fma
   write(*,*) 'Plan exec time        : ',deltaT,' sec (average over ', N_ITER,' iteration(s))'
   write(*,*) 'Theoretical FLOP count: ', flops
   write(*,*) 'Plan MFLOPS estimate  : ', 1e-6*flops/deltaT

   call dfftw_destroy_plan(plan)

   deallocate(in)
   deallocate(out)
   
   call dfftw_cleanup_threads();
   
end program fftw_test_3d_multithread

