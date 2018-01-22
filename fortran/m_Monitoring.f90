module Monitoring

  use Timing
  use Precision

  ! declare a Timer type
  type Timer
     real(fp_kind) :: elapsed=0.0
     real(fp_kind) :: start=0.0
     real(fp_kind) :: stop=0.0
  end type Timer
  
  type(Timer) :: fft_plan_timer
  type(Timer) :: fft_exec_timer


  contains

    ! start timer
    subroutine timerStart(t)
      implicit none
      type(Timer), intent(inout) :: t

      !call cpu_time(t%start)
      t%start = wallclock()

    end subroutine timerStart

    ! stop timer and accumulate timings
    subroutine timerStop(t)
      implicit none
      type(Timer), intent(inout) :: t

      !call cpu_time(t%stop)
      t%stop = wallclock()

      t%elapsed = t%elapsed + t%stop - t%start

    end subroutine timerStop

end module Monitoring
