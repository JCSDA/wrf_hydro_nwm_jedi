module wrf_hydro_jedi_geometry_mod

use fckit_configuration_module, only: fckit_configuration
use netcdf

implicit none
private

public :: wrf_hydro_jedi_geometry

!------------------------------------------------------------------------------

type :: wrf_hydro_jedi_geometry
   integer :: ix, iy
   integer :: xstart, ystart
   integer :: xend, yend
   real :: dx, dy
contains
  procedure :: init   => wrf_hydro_jedi_geometry_init
  procedure :: clone  => wrf_hydro_jedi_geometry_clone
  procedure :: delete => wrf_hydro_jedi_geometry_delete
  procedure :: get_info => wrf_hydro_jedi_geometry_get_info
end type

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine wrf_hydro_jedi_geometry_init(self, f_conf)
  class(wrf_hydro_jedi_geometry),   intent(out) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  character(512) :: wrfinput_flnm
  integer :: ierr,ncid
  
  !call abor1_ftn("ERROR: wrf_hydro_jedi_geometry_init() needs to be implemented.")

  wrfinput_flnm = '/home/alex/Downloads/JEDI/src/build/wrf_hydro_jedi/test/testinput/wrfinput.nc'

  self%dx = 1
  ierr = 0
  write(*,*) 'Reading from NETCDF file the Geometry'

  ierr = nf90_open(wrfinput_flnm, NF90_NOWRITE, ncid)

  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn("ERROR: wrfinput file needed by Geometry not found")
  end if

  ierr = nf90_get_att(ncid, NF90_GLOBAL, "DX", self%dx)
  
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "DY", self%dy)

  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn("ERROR: dx and/or dy need by Geometry not found")
  end if

  write(*,*) "Dx and dy from Fortran", self%dx, self%dy
  
end subroutine

!------------------------------------------------------------------------------

subroutine wrf_hydro_jedi_geometry_clone(self, other)
  class(wrf_hydro_jedi_geometry),  intent(out) :: self
  class(wrf_hydro_jedi_geometry),   intent(in) :: other
      
  !call abor1_ftn("ERROR: wrf_hydro_jedi_geometry_clone() needs to be implemented.")
end subroutine
  
!------------------------------------------------------------------------------

subroutine wrf_hydro_jedi_geometry_delete(self)
  class(wrf_hydro_jedi_geometry),  intent(out) :: self
       
  !call abor1_ftn("ERROR: wrf_hydro_jedi_geometry_delete() needs to be implemented.")
end subroutine

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

subroutine wrf_hydro_jedi_geometry_get_info(self,dx_,dy_)
  class(wrf_hydro_jedi_geometry),  intent(in) :: self
  real, intent(out) :: dx_,dy_

  dx_ = self%dx
  dy_ = self%dy

end subroutine wrf_hydro_jedi_geometry_get_info

!------------------------------------------------------------------------------

end module
