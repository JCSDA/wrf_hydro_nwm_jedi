!> General utilities module
module wrf_hydro_nwm_jedi_util_mod

use netcdf
use datetime_mod

use fckit_configuration_module, only: fckit_configuration
use duration_mod, only: duration, duration_to_string

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


! ------------------------------------------------------------------------------
!> Generate filename (based on oops/qg)
!!
!! The configuration \p f_conf is expected to provide the following
!! - "datadir" : the directory the filenames should be prefixed with
!! - "exp" : experiment name
!! - "type" : one of "fc", "an", "incr", "ens"
!! - "member" : required only if "type == ens"

function genfilename (f_conf,length,vdate,domain_type)
   type(fckit_configuration),  intent(in) :: f_conf
   integer,                    intent(in) :: length
   type(datetime),             intent(in) :: vdate
   character(len=3), optional, intent(in) :: domain_type
 
   character(len=length)                  :: genfilename
   character(len=length) :: fdbdir, expver, typ, validitydate, referencedate, sstep, &
        & prefix, mmb
   type(datetime) :: rdate
   type(duration) :: step
   integer lenfn
   character(len=:), allocatable :: str
 
   call f_conf%get_or_die("datadir", str)
   fdbdir = str
   call f_conf%get_or_die("exp", str)
   expver = str
   call f_conf%get_or_die("type", str)
   typ = str
  
   if (typ=="ens") then
      call f_conf%get_or_die("member", str)
      mmb = str
      lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
      prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(domain_type) // "." // TRIM(typ) // "." // TRIM(mmb)
   else
      lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
      prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(domain_type) // "." // TRIM(typ)
   endif
 
   if (typ=="fc" .or. typ=="ens") then
      call f_conf%get_or_die("date", str)
      referencedate = str
      call datetime_to_string(vdate, validitydate)
      call datetime_create(TRIM(referencedate), rdate)
      call datetime_diff(vdate, rdate, step)
      call duration_to_string(step, sstep)
      lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
      genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep)
   endif
 
   if (typ=="an" .or. typ=="incr") then
      call datetime_to_string(vdate, validitydate)
      lenfn = lenfn + 1 + LEN_TRIM(validitydate)
      genfilename = TRIM(prefix) // "." // TRIM(validitydate)
   endif
 
   if (lenfn>length) &
        & call abor1_ftn("fields:genfilename: filename too long")
 
    if ( allocated(str) ) deallocate(str)
 
 end function genfilename

end module wrf_hydro_nwm_jedi_util_mod
