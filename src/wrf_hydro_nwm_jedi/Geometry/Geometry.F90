module wrf_hydro_nwm_jedi_geometry_mod

use fckit_configuration_module, only: fckit_configuration
use netcdf

implicit none
private

public :: wrf_hydro_nwm_jedi_geometry, indices, error_handler

!------------------------------------------------------------------------------
! General:
! The dims are generically named: 1=x, 2=y, 3=z

type wrf_hydro_nwm_lsm_geometry
   integer :: xstart, ystart, zstart
   integer :: xend, yend, zend ! same as xdim_len and ydim_len? redundant?
   ! why not zstart, zend?
   integer :: xdim_len, ydim_len, zdim_len
   integer :: n_snow_layers
   real    :: dx, dy
   real, allocatable :: lat(:,:), lon(:,:)
end type wrf_hydro_nwm_lsm_geometry

type wrf_hydro_nwm_stream_geometry
   integer :: xstart, xend
   integer :: xdim_len
   real    :: dx
   real, allocatable :: lat(:), lon(:)
end type wrf_hydro_nwm_stream_geometry

type :: wrf_hydro_nwm_jedi_geometry
   ! integer :: ix, iy  ! unused
   type(wrf_hydro_nwm_lsm_geometry) :: lsm
   type(wrf_hydro_nwm_stream_geometry) :: stream
 contains
   procedure :: init   => wrf_hydro_nwm_jedi_geometry_init
   procedure :: clone  => wrf_hydro_nwm_jedi_geometry_clone
   procedure :: delete => wrf_hydro_nwm_jedi_geometry_delete
   procedure :: get_info => wrf_hydro_nwm_jedi_geometry_get_info
   procedure :: get_nn => wrf_hydro_nwm_jedi_geometry_get_nn
end type wrf_hydro_nwm_jedi_geometry

! This should go into Utilities
type :: indices
   integer :: ind_1, ind_2, ind_3
end type indices


!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------


subroutine wrf_hydro_nwm_jedi_geometry_init(self, f_conf)
  class(wrf_hydro_nwm_jedi_geometry),   intent(out) :: self
  type(fckit_configuration),  intent(in) :: f_conf

  ! LSM init required
  call wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf)
  ! Streamflow optional TODO JLM how to signal?  
end subroutine wrf_hydro_nwm_jedi_geometry_init


subroutine get_geom_dim_len(ncid, meta_name, dim_len)
  integer,          intent(in)  :: ncid
  character(len=*), intent(in)  :: meta_name
  integer,          intent(out) :: dim_len

  integer :: ierr, dimid
  character(len=512) :: dim_name

  ! The metadata name / attribute reveals the actual variable name.
  ! This provides a buffer against changes to the model files.
  ierr = nf90_get_att(ncid, NF90_GLOBAL, meta_name, dim_name)
  call error_handler(ierr, &
       "STOP:  Problem finding '"//trim(meta_name)//"' attribute")
  ierr = nf90_inq_dimid(ncid, trim(dim_name), dimid)
  call error_handler(ierr, &
       "STOP:  Problem finding '"//trim(dim_name)//"' dimension")
  ierr = nf90_inquire_dimension(ncid, dimid, len=dim_len)
  call error_handler(ierr, &
       "STOP:  Problem getting '"//trim(dim_name)//"' dimension length")
end subroutine get_geom_dim_len


subroutine wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf)
  ! LSM component is not optional
  class(wrf_hydro_nwm_jedi_geometry),   intent(inout) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  character(512) :: geometry_flnm
  character(len=:), allocatable :: str
  integer :: ierr, ncid, dimid

  call f_conf%get_or_die("input_file",str)
  geometry_flnm = str
  deallocate(str)
  ierr = 0
  write(*,*) "Reading geometry file: "//trim(geometry_flnm)//"'"
  ierr = nf90_open(geometry_flnm, NF90_NOWRITE, ncid)
  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn( "ERROR: geometry file not found")
  end if

  ! Hard-coded values
  self%lsm%xstart = 1
  self%lsm%ystart = 1
  self%lsm%zstart = 1
  self%lsm%n_snow_layers = 3  ! This is max, hard-coded.

  ! Metadata / atts
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "lsm_dx", self%lsm%dx)
  call error_handler(ierr, "STOP:  Problem finding DX attribute")
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "lsm_dy", self%lsm%dy)
  call error_handler(ierr, "STOP:  Problem finding DY attribute")  

  ! Dimension data
  call get_geom_dim_len(ncid, "lsm_lon_dim_name", self%lsm%xdim_len)
  self%lsm%xend = self%lsm%xdim_len
  call get_geom_dim_len(ncid, "lsm_lat_dim_name", self%lsm%ydim_len)
  self%lsm%yend = self%lsm%ydim_len
  call get_geom_dim_len(ncid, "lsm_z_dim_name", self%lsm%zdim_len)
  self%lsm%zend = self%lsm%zdim_len

  ! Positon data
  allocate( &
       self%lsm%lat(self%lsm%xdim_len, self%lsm%ydim_len), &
       self%lsm%lon(self%lsm%xdim_len, self%lsm%ydim_len))
  call get_2d("XLAT", ncid, self%lsm%lat, self%lsm%xdim_len, self%lsm%ydim_len)
  call get_2d("XLONG", ncid, self%lsm%lon, self%lsm%xdim_len, self%lsm%ydim_len)
end subroutine wrf_hydro_nwm_jedi_lsm_geometry_init


subroutine wrf_hydro_nwm_jedi_geometry_get_nn(self, lat, lon, ind)!dim1_idx, dim2_idx)
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: self
  real, intent(in) :: lat, lon
  real,dimension(2) :: minimum
  real, allocatable :: diff_lat(:,:), diff_lon(:,:), l2_norm(:,:)
  type(indices), intent(out) :: ind

  allocate(diff_lat, source=self%lsm%lat)
  allocate(diff_lon, source=self%lsm%lon)
  allocate(l2_norm, source=self%lsm%lon)
  
  diff_lat = diff_lat - lat
  diff_lon = diff_lon - lon

  l2_norm = sqrt( diff_lon**2 + diff_lat**2 )

  minimum = minloc(l2_norm)

  ind%ind_1 = minimum(1); ind%ind_2 = minimum(2)

  deallocate(l2_norm)
  deallocate(diff_lat)
  deallocate(diff_lon)  
end subroutine wrf_hydro_nwm_jedi_geometry_get_nn


subroutine wrf_hydro_nwm_jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: other

  self%lsm%dx = other%lsm%dx
  self%lsm%dy = other%lsm%dy
  self%lsm%xdim_len = other%lsm%xdim_len
  self%lsm%ydim_len = other%lsm%ydim_len
  self%lsm%zdim_len = other%lsm%zdim_len
  self%lsm%n_snow_layers = other%lsm%n_snow_layers
  self%lsm%xstart = other%lsm%xstart
  self%lsm%ystart = other%lsm%ystart
  self%lsm%xend = other%lsm%xend
  self%lsm%yend = other%lsm%yend
  allocate(self%lsm%lat, source = other%lsm%lat) 
  allocate(self%lsm%lon, source = other%lsm%lon)  
end subroutine
  

subroutine wrf_hydro_nwm_jedi_geometry_delete(self)
  class(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: self

  deallocate(self%lsm%lat)
  deallocate(self%lsm%lon)
end subroutine


subroutine wrf_hydro_nwm_jedi_geometry_get_info( &
     self, &
     dx, dy, &
     xdim_len, ydim_len, zdim_len)
  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self
  real, intent(out) :: dx, dy
  integer, intent(out) :: xdim_len, ydim_len, zdim_len

  dx = self%lsm%dx
  dy = self%lsm%dy
  xdim_len = self%lsm%xdim_len
  ydim_len = self%lsm%ydim_len
  zdim_len = self%lsm%zdim_len 
end subroutine wrf_hydro_nwm_jedi_geometry_get_info


!Subroutine copied from module_wrfinputfile.F
subroutine get_2d(name, ncid, array, idim, jdim)
  ! From the NetCDF unit <ncid>, read the variable named <name> with  
  ! dimensions <idim> and <jdim>, filling the pre-dimensioned array <array>
  implicit none
  character(len=*),           intent(in)  :: name
  integer,                    intent(in)  :: ncid
  integer,                    intent(in)  :: idim
  integer,                    intent(in)  :: jdim
  real, dimension(idim,jdim), intent(out) :: array
  ! Local:
  integer                                 :: ierr
  integer                                 :: varid

  ierr = nf90_inq_varid(ncid,  name,  varid)
  ! If the variable is "XLAT", and "XLAT" is not found, look for "XLAT_M"
  ! If the variable is "XLAT_M", and "XLAT_M" is not found, look for "XLAT"
  ! If the variable is "XLONG", and "XLONG" is not found, look for "XLONG_M"
  ! If the variable is "XLONG_M", and "XLONG_M" is not found, look for "XLONG"
  if (name == "XLAT") then
     if (ierr /= 0) then
        ierr = nf90_inq_varid(ncid,  "XLAT_M",  varid)
     endif
  else if (name == "XLAT_M") then
     if (ierr /= 0) then
        ierr = nf90_inq_varid(ncid,  "XLAT",  varid)
     endif
  else  if (name == "XLONG") then
     if (ierr /= 0) then
        ierr = nf90_inq_varid(ncid,  "XLONG_M",  varid)
     endif
  else if (name == "XLONG_M") then
     if (ierr /= 0) then
        ierr = nf90_inq_varid(ncid,  "XLONG",  varid)
     endif
  endif
  call error_handler(ierr, &
       "STOP:  READ_WRFINPUT: Problem finding variable '"//trim(name)//"' in wrfinput file.")

  ierr = nf90_get_var(ncid, varid, array)
  call error_handler(ierr, &
       "STOP:  READ_WRFINPUT:  Problem retrieving variable '"//trim(name)//"' from wrfinput file.")
end subroutine get_2d


! Move this to utilities.
subroutine error_handler(status, failure, success)
  !
  ! Check the error flag from a NetCDF function call, and print appropriate
  ! error message.
  !
  implicit none
  integer,                    intent(in) :: status
  character(len=*), optional, intent(in) :: failure
  character(len=*), optional, intent(in) :: success
  
  if (status .ne. NF90_NOERR) then
     if (present(failure)) then
        write(*,'(/," ***** ", A)') failure
     endif
     write(*,'(" ***** ",A,/)') nf90_strerror(status)
     stop 'FATAL ERROR: In module_wrfinputfile.F -- Stopped'
  endif
  
  if (present(success)) then
     write(*,'(A)') success
  endif 
end subroutine error_handler

end module
