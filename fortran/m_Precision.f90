!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!< REAL and INTEGER precision
!! Note that REAL32 and INT32 are not available in gfortran-4.4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Precision

  integer, parameter :: sp_kind   = kind(1.0)  !< single precision
  integer, parameter :: dp_kind   = kind(1.d0) !< double precision
  integer, parameter :: int_kind  = kind(1)

  ! default: use double precison
  integer, parameter :: fp_kind = dp_kind

  ! TODO: use iso_fortran_env module (Fortran 2003)

  contains
    
    function isBigEndian()
      integer, parameter :: short = selected_int_kind(4)
      integer( short ) :: source = 1_short
      logical :: isBigEndian
      isBigEndian = .false.
      if ( iachar( transfer( source, 'a' ) ) == 0 ) isBigEndian = .true.
      return
    end function isBigEndian

    function useDoublePrecision()
      real(fp_kind) :: tmp
      logical       :: useDoublePrecision
      useDoublePrecision = .false.
      if (kind(tmp) == dp_kind ) useDoublePrecision=.true.
    end function useDoublePrecision

end module Precision
