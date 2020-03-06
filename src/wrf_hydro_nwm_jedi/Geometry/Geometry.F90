module wrf_hydro_nwm_jedi_geometry_mod

use fckit_configuration_module, only: fckit_configuration
use netcdf

implicit none
private

public :: wrf_hydro_nwm_jedi_geometry

!------------------------------------------------------------------------------

type :: wrf_hydro_nwm_jedi_geometry
   integer :: ix, iy !unsused
   integer :: xstart, ystart !unsused
   integer :: xend, yend !unsused
   integer :: npx, npy, npz
   real :: dx, dy
 contains
   procedure :: init   => wrf_hydro_nwm_jedi_geometry_init
   procedure :: clone  => wrf_hydro_nwm_jedi_geometry_clone
   procedure :: delete => wrf_hydro_nwm_jedi_geometry_delete
   procedure :: get_info => wrf_hydro_nwm_jedi_geometry_get_info
end type wrf_hydro_nwm_jedi_geometry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_init(self, f_conf)
  class(wrf_hydro_nwm_jedi_geometry),   intent(out) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  character(512) :: wrfinput_flnm
  character(len=:), allocatable :: str
  integer :: ierr,ncid

  call f_conf%get_or_die("input_file",str)
  wrfinput_flnm = str
  write(*,*) 'Input file for Geometry: ',wrfinput_flnm
  deallocate(str)

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
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "west_east", self%npx)
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "south_north", self%npy)
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "soil_layers_stag", self%npz)
  
  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn("ERROR: dx and/or dy need by Geometry not found")
  end if
  
end subroutine

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: other

  self%dx = other%dx
  self%dy = other%dy
  self%npx = other%npx
  self%npy = other%npy
  self%npz = other%npz
  
end subroutine
  
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_delete(self)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self

end subroutine

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_get_info(self,dx_,dy_,npx_,npy_,npz_)
  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self
  real, intent(out) :: dx_,dy_
  integer, intent(out) :: npx_,npy_,npz_

  dx_ = self%dx
  dy_ = self%dy
  npx_ = self%npx
  npy_ = self%npy
  npz_ = self%npz  

end subroutine wrf_hydro_nwm_jedi_geometry_get_info

!------------------------------------------------------------------------------

end module
