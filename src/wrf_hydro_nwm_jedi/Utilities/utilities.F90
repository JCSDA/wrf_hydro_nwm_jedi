!> General utilities module
module wrf_hydro_nwm_jedi_util_mod

use netcdf
use datetime_mod

implicit none

public


!> Indices type
type :: indices
   integer :: ind_x, ind_y, ind_z
end type indices


contains

  
!> A netcdf error handler.  
!> Check the error flag from a NetCDF function call, and print
!>  appropriate/optional error or success message.
subroutine error_handler(status, failure, success)
  integer,                    intent(in) :: status   !< the returned code
  character(len=*), optional, intent(in) :: failure  !< mesage for failure
  character(len=*), optional, intent(in) :: success  !< message for success
  
  if (status .ne. NF90_NOERR) then
     if (present(failure)) then
        write(*,'(/," ***** ", A)') failure
     endif
     write(*,'(" ***** ",A,/)') nf90_strerror(status)
     write(*,*) 'FATAL ERROR: Stopped'
     stop 123
  endif
  
  if (present(success)) write(*,'(A)') success
end subroutine error_handler


!> Date time equality function. Shorthand.
logical function datetime_eq(dt1, dt2)
  implicit none
  type(datetime), intent(in) :: dt1, dt2

  if( (dt1 <= dt2) .and. (dt1 >= dt2)) then
     datetime_eq = .TRUE.
  else
     datetime_eq = .FALSE.
  end if
end function datetime_eq

  
end module wrf_hydro_nwm_jedi_util_mod
