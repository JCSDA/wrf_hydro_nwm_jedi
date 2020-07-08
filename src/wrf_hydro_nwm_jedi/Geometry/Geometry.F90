module wrf_hydro_nwm_jedi_geometry_mod

use fckit_configuration_module, only: fckit_configuration
use netcdf

implicit none
private

public :: wrf_hydro_nwm_jedi_geometry, indices, error_handler

!------------------------------------------------------------------------------
! Data types
! General:
!   The dims are generically named: 1=x, 2=y, 3=z

type wrf_hydro_nwm_lsm_geometry
   integer :: xstart, ystart, zstart
   integer :: xend, yend, zend ! same as xdim_len and ydim_len? redundant?
   integer :: xdim_len, ydim_len, zdim_len
   integer :: n_snow_layers
   real    :: dx, dy
   real, allocatable :: lat(:,:), lon(:,:)
end type wrf_hydro_nwm_lsm_geometry


type wrf_hydro_nwm_stream_geometry
   integer :: xstart, xend
   integer :: xdim_len
   real, allocatable :: lat(:), lon(:), dx(:)
end type wrf_hydro_nwm_stream_geometry


type :: wrf_hydro_nwm_jedi_geometry
   ! integer :: ix, iy  ! unused
   type(wrf_hydro_nwm_lsm_geometry) :: lsm
   type(wrf_hydro_nwm_stream_geometry) :: stream
 contains
   procedure :: init   => wrf_hydro_nwm_jedi_geometry_init
   procedure :: clone  => wrf_hydro_nwm_jedi_geometry_clone
   procedure :: delete => wrf_hydro_nwm_jedi_geometry_delete
   procedure :: get_lsm_info => wrf_hydro_nwm_jedi_geometry_get_lsm_info
   procedure :: get_lsm_nn => wrf_hydro_nwm_jedi_geometry_get_lsm_nn
end type wrf_hydro_nwm_jedi_geometry


! This should go into Utilities
type :: indices
   integer :: ind_1, ind_2, ind_3
end type indices


!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!> Initialize the geometry object from the preprocessed geometry input file.
subroutine wrf_hydro_nwm_jedi_geometry_init(self, f_conf)
  class(wrf_hydro_nwm_jedi_geometry),   intent(out) :: self  !< the geom object
  type(fckit_configuration),  intent(in) :: f_conf  !< the yaml file

  character(512) :: geometry_flnm
  character(len=:), allocatable :: str
  integer :: ierr, ncid

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

  ! LSM init required
  call wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf, ncid)

  ! Streamflow optional TODO JLM how to signal/handle?
  call wrf_hydro_nwm_jedi_stream_geometry_init(self, f_conf, ncid)

end subroutine wrf_hydro_nwm_jedi_geometry_init


!> Clone the geometry object
subroutine wrf_hydro_nwm_jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self  !< geom object
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: other  !< the new geom object

  ! LSM
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

  ! Stream
end subroutine wrf_hydro_nwm_jedi_geometry_clone


!> Destroy the geometry object
subroutine wrf_hydro_nwm_jedi_geometry_delete(self)
  class(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: self  !< geom object

  deallocate(self%lsm%lat)
  deallocate(self%lsm%lon)
end subroutine wrf_hydro_nwm_jedi_geometry_delete


!-----------------------------------------------------------------------------
! LSM Section

!> Init the LSM geometry object
subroutine wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf, ncid)
  ! LSM component is not optional
  class(wrf_hydro_nwm_jedi_geometry), intent(inout) :: self  !< the full geom object
  type(fckit_configuration), intent(in) :: f_conf  !< the yaml file
  integer, intent(in) :: ncid
  
  integer :: ierr, dimid

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
  call get_geom_dim_len(ncid, "lsm_xdim_name", self%lsm%xdim_len)
  self%lsm%xend = self%lsm%xdim_len
  call get_geom_dim_len(ncid, "lsm_ydim_name", self%lsm%ydim_len)
  self%lsm%yend = self%lsm%ydim_len
  call get_geom_dim_len(ncid, "lsm_zdim_name", self%lsm%zdim_len)
  self%lsm%zend = self%lsm%zdim_len

  ! Positon data
  allocate( &
       self%lsm%lat(self%lsm%xdim_len, self%lsm%ydim_len), &
       self%lsm%lon(self%lsm%xdim_len, self%lsm%ydim_len))
  ! Could this be polymorphic?
  call get_geom_data_2d("lsm_lat_name", ncid, self%lsm%lat)
  call get_geom_data_2d("lsm_lon_name", ncid, self%lsm%lon)
end subroutine wrf_hydro_nwm_jedi_lsm_geometry_init


!> Get the LSM nearest neighbor from a lat and lon pair.
! @todo change indices to have dimension?
subroutine wrf_hydro_nwm_jedi_geometry_get_lsm_nn( &
     self, lat, lon, ind)
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: self  !< geom object
  real, intent(in) :: lat, lon  !< the lat, lon of interest
  type(indices), intent(out) :: ind  !< the index pair of the nearest neighbor

  real,dimension(2) :: minimum
  real, allocatable :: diff_lat(:,:), diff_lon(:,:), l2_norm(:,:)

  allocate(diff_lat, source=self%lsm%lat)
  allocate(diff_lon, source=self%lsm%lon)
  allocate(l2_norm, source=self%lsm%lon)
  
  diff_lat = diff_lat - lat
  diff_lon = diff_lon - lon
  l2_norm = sqrt( diff_lon**2 + diff_lat**2 )
  minimum = minloc(l2_norm)
  ind%ind_1 = minimum(1)
  ind%ind_2 = minimum(2)

  deallocate(l2_norm)
  deallocate(diff_lat)
  deallocate(diff_lon)  
end subroutine wrf_hydro_nwm_jedi_geometry_get_lsm_nn


!> Get lsm info
subroutine wrf_hydro_nwm_jedi_geometry_get_lsm_info( &
     self, &
     dx, dy, &
     xdim_len, ydim_len, zdim_len)
  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self  !< geom object
  real, intent(out) :: dx, dy
  integer, intent(out) :: xdim_len, ydim_len, zdim_len

  dx = self%lsm%dx
  dy = self%lsm%dy
  xdim_len = self%lsm%xdim_len
  ydim_len = self%lsm%ydim_len
  zdim_len = self%lsm%zdim_len 
end subroutine wrf_hydro_nwm_jedi_geometry_get_lsm_info


!-----------------------------------------------------------------------------
! Streamflow section

!> Init the streamflow geometry object
subroutine wrf_hydro_nwm_jedi_stream_geometry_init(self, f_conf, ncid)
  ! stream component IS optional
  class(wrf_hydro_nwm_jedi_geometry), intent(inout) :: self  !< the full geom object
  type(fckit_configuration), intent(in) :: f_conf  !< the yaml file
  integer, intent(in) :: ncid !< the geom file netcdf id

  integer :: ierr, dimid

  ! Hard-coded values
  self%stream%xstart = 1

  ! Metadata / atts: None

  ! Dimension data
  call get_geom_dim_len(ncid, "stream_xdim_name", self%stream%xdim_len)
  self%stream%xend = self%stream%xdim_len

  ! Positon data
  allocate( &
       self%stream%lat(self%stream%xdim_len), &
       self%stream%lon(self%stream%xdim_len), &
       self%stream%dx(self%stream%xdim_len))

  call get_geom_data_1d("stream_lon_name", ncid, self%stream%lat)
  call get_geom_data_1d("stream_lat_name", ncid, self%stream%lon)
  call get_geom_data_1d("stream_dx_name",  ncid, self%stream%dx)
end subroutine wrf_hydro_nwm_jedi_stream_geometry_init


! -----------------------------------------------------------------------------
! Helper functions for preprocessed geometry file operations

!> Helper function for reading geometry file dimension information
subroutine get_geom_dim_len(ncid, meta_name, dim_len)
  integer,          intent(in)  :: ncid  !< the netcdf file id
  character(len=*), intent(in)  :: meta_name  !< the meta name in file attrs
  integer,          intent(out) :: dim_len  !< the dim length

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


!> Helper function for reading geometry file data (1D)
! @todo make this polymorphic with 1,2,3D versions
subroutine get_geom_data_1d(meta_name, ncid, vector)
  character(len=*),   intent(in)  :: meta_name
  integer,            intent(in)  :: ncid
  real, dimension(:), intent(out) :: vector

  integer :: ierr, varid
  character(len=512) :: var_name

  ! The metadata name / attribute reveals the actual variable name.
  ! This provides a buffer against changes to the model files.
  ierr = nf90_get_att(ncid, NF90_GLOBAL, meta_name, var_name)
  call error_handler(ierr, &
       "STOP:  Problem finding '"//trim(meta_name)//"' attribute")
  ierr = nf90_inq_varid(ncid, var_name, varid)
  call error_handler(ierr, &
       "STOP:  READ_WRFINPUT: Problem finding variable '"// &
       trim(var_name)//"' in geometry file")
  ierr = nf90_get_var(ncid, varid, vector)
  call error_handler(ierr, &
       "STOP:  READ_WRFINPUT: Problem getting variable '"// &
       trim(var_name)//"' in geometry file")
end subroutine get_geom_data_1d


!> Helper function for reading geometry file data (2D)
! @todo make this polymorphic
subroutine get_geom_data_2d(meta_name, ncid, array)
  character(len=*),      intent(in)  :: meta_name
  integer,               intent(in)  :: ncid
  real, dimension(:, :), intent(out) :: array

  integer :: ierr, varid
  character(len=512) :: var_name

  ! The metadata name / attribute reveals the actual variable name.
  ! This provides a buffer against changes to the model files.
  ierr = nf90_get_att(ncid, NF90_GLOBAL, meta_name, var_name)
  call error_handler(ierr, &
       "STOP:  Problem finding '"//trim(meta_name)//"' attribute")
  ierr = nf90_inq_varid(ncid, var_name, varid)
  call error_handler(ierr, &
       "STOP:  READ_WRFINPUT: Problem finding variable '"// &
       trim(var_name)//"' in geometry file")
  ierr = nf90_get_var(ncid, varid, array)
  call error_handler(ierr, &
       "STOP:  READ_WRFINPUT: Problem getting variable '"// &
       trim(var_name)//"' in geometry file")
end subroutine get_geom_data_2d


! -----------------------------------------------------------------------------
! Utilities to be moved to a separate file.

!> A netcdf error handler.
subroutine error_handler(status, failure, success)
  !
  ! Check the error flag from a NetCDF function call, and print appropriate
  ! error message.
  !
  integer,                    intent(in) :: status   !< the returned code
  character(len=*), optional, intent(in) :: failure  !< mesage for failure
  character(len=*), optional, intent(in) :: success  !< message for success
  
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
