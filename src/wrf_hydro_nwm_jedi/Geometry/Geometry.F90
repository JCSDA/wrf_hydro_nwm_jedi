
!> Geometry (location information) implementation for wrf_hydro_nwm - jedi integration.
module wrf_hydro_nwm_jedi_geometry_mod

use fckit_configuration_module, only: fckit_configuration
use wrf_hydro_nwm_jedi_util_mod, only: error_handler, indices
use netcdf
use fckit_log_module,            only: fckit_log

implicit none
private

! For doxygen purposes, this public statement is just a summary for the code reader?
public :: wrf_hydro_nwm_jedi_geometry, get_lsm_nn, get_geoval_levels


! General:
!   The dims are generically named: 1=x, 2=y, 3=z

!> Geometry for LSM (land surface model) fields (2D and 3D)
type, private :: wrf_hydro_nwm_lsm_geometry
   integer :: xstart, ystart, zstart
   integer :: xend, yend, zend !< @todo same as xdim_len and ydim_len? redundant?
   integer :: xdim_len, ydim_len, zdim_len
   integer :: n_snow_layers
   real    :: dx, dy
   real, allocatable :: lat(:,:), lon(:,:), sfc_elev(:,:)
end type wrf_hydro_nwm_lsm_geometry


!> Geometry for stream channel/reach network (1D)
type, private :: wrf_hydro_nwm_stream_geometry
   integer :: xstart, xend
   integer :: xdim_len
   real, allocatable :: lat(:), lon(:), dx(:)
end type wrf_hydro_nwm_stream_geometry


!> Fortran geometry object (contains specific geometry components).
type, public :: wrf_hydro_nwm_jedi_geometry
   type(wrf_hydro_nwm_lsm_geometry), allocatable :: lsm        !< LSM geom component (optional)
   type(wrf_hydro_nwm_stream_geometry), allocatable :: stream  !< stream geom component (optional)
 contains
   procedure :: init   => wrf_hydro_nwm_jedi_geometry_init    !< init/create
   procedure :: clone  => wrf_hydro_nwm_jedi_geometry_clone   !< copy
   procedure :: delete => wrf_hydro_nwm_jedi_geometry_delete  !< delete
   procedure :: get_lsm_info => get_lsm_info  !< get lsm info
   procedure :: get_lsm_nn => get_lsm_nn  !< get lsm nearest neighbor
   procedure :: get_geoval_levels => get_geoval_levels
   ! procedure :: get_stream_info => wrf_hydro_nwm_jedi_geometry_get_stream_info
   ! procedure :: get_stream_nn => wrf_hydro_nwm_jedi_geometry_get_stream_nn
   procedure :: lsm_active => lsm_active  !< query if the lsm component is active/allocated
   procedure :: stream_active => stream_active  !< query if the stream component is active/allocated
end type wrf_hydro_nwm_jedi_geometry

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------


!> Initialize the geometry object from the preprocessed geometry input file.
subroutine wrf_hydro_nwm_jedi_geometry_init(self, f_conf)
  class(wrf_hydro_nwm_jedi_geometry), intent(out) :: self  !< the geom object
  type(fckit_configuration), intent(in) :: f_conf  !< the yaml file

  character(512) :: geometry_flnm
  character(len=:), allocatable :: str
  integer :: ierr, ncid

  call f_conf%get_or_die("input_file",str)
  geometry_flnm = str
  deallocate(str)
  ierr = 0
  write(*,*) "Reading geometry file: "//trim(geometry_flnm)//"'"
  ierr = nf90_open(geometry_flnm, NF90_NOWRITE, ncid)
  call error_handler(ierr, "STOP: geometry file not found")

  call wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf, ncid)
  call wrf_hydro_nwm_jedi_stream_geometry_init(self, f_conf, ncid)

  ierr = nf90_close(ncid)
  call error_handler(ierr, "STOP: geometry file can not be closed")
end subroutine wrf_hydro_nwm_jedi_geometry_init


!> Clone the geometry object
subroutine wrf_hydro_nwm_jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self  !< geom object
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: other  !< the new geom object

  ! LSM
  if (other%lsm_active()) then
     allocate(self%lsm)
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
  end if

  ! Stream
  if (other%stream_active()) then
     allocate(self%stream)
     self%stream%dx = other%stream%dx
     self%stream%xdim_len = other%stream%xdim_len
     self%stream%xstart = other%stream%xstart
     self%stream%xend = other%stream%xend
     allocate(self%stream%lat, source = other%stream%lat)
     allocate(self%stream%lon, source = other%stream%lon)
  end if
end subroutine wrf_hydro_nwm_jedi_geometry_clone


!> Destroy the geometry object
subroutine wrf_hydro_nwm_jedi_geometry_delete(self)
  class(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: self  !< geom object

  ! TODO JLM: this should be recursive?
  if (self%lsm_active()) deallocate(self%lsm)
  if (self%stream_active()) deallocate(self%stream)
end subroutine wrf_hydro_nwm_jedi_geometry_delete


!-----------------------------------------------------------------------------
! LSM Section


!> Check if the lsm geometry object is allocated/active.
function lsm_active(self) result(is_active)
  class(wrf_hydro_nwm_jedi_geometry), intent(in) :: self  !< the full geom object
  logical :: is_active

  is_active = allocated(self%lsm)
end function lsm_active


!> Init the LSM geometry object
subroutine wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf, ncid)
  ! LSM component is not optional
  class(wrf_hydro_nwm_jedi_geometry), intent(inout) :: self  !< the full geom object
  type(fckit_configuration), intent(in) :: f_conf  !< the yaml file
  integer, intent(in) :: ncid
  
  integer :: ierr, dimid
  character(len=512) :: src_file_name

  ! Check LSM active
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "lsm_src_file", src_file_name)
  call error_handler(ierr, "STOP:  lsm_src_file attribute not defined in geometry file")
  if (ierr .ne. NF90_NOERR) then
     return
  end if

  allocate(self%lsm)

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
       self%lsm%lon(self%lsm%xdim_len, self%lsm%ydim_len), &
       self%lsm%sfc_elev(self%lsm%xdim_len, self%lsm%ydim_len))
  ! Could this be polymorphic?
  call get_geom_data_2d("lsm_lat_name", ncid, self%lsm%lat)
  call get_geom_data_2d("lsm_lon_name", ncid, self%lsm%lon)
  call get_geom_data_2d("lsm_sfc_elev_name", ncid, self%lsm%sfc_elev)
end subroutine wrf_hydro_nwm_jedi_lsm_geometry_init


!> Get the LSM nearest neighbor from a lat and lon pair.
! @todo change indices to have dimension?
subroutine get_lsm_nn( &
     self, lat, lon, ind)
  class(wrf_hydro_nwm_jedi_geometry), intent(in) :: self  !< geom object
  real, intent(in) :: lat  !< the lat of interest
  real, intent(in) :: lon  !< the lon of interest
  type(indices), intent(out) :: ind  !< the index pair of the nearest neighbor

  real,dimension(2) :: minimum
  real, allocatable :: diff_lat(:,:), diff_lon(:,:), l2_norm(:,:)

  ! TODO JLM: I dont love this...
  if (.not. self%lsm_active()) return

  allocate(diff_lat, source=self%lsm%lat)
  allocate(diff_lon, source=self%lsm%lon)
  allocate(l2_norm, source=self%lsm%lon)

  diff_lat = diff_lat - lat
  diff_lon = diff_lon - lon
  l2_norm = sqrt( diff_lon**2 + diff_lat**2 )
  minimum = minloc(l2_norm)
  ind%ind_x = minimum(1)
  ind%ind_y = minimum(2)

  deallocate(l2_norm)
  deallocate(diff_lat)
  deallocate(diff_lon)
end subroutine get_lsm_nn


!> Get lsm geometry info
subroutine get_lsm_info( &
     self, &
     dx, dy, &
     xdim_len, ydim_len, zdim_len)
  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self  !< geom object
  real, intent(out) :: dx  !< the x resolution
  real, intent(out) :: dy  !< the y resolution
  integer, intent(out) :: xdim_len  !< x dimension size
  integer, intent(out) :: ydim_len  !< y dimension size
  integer, intent(out) :: zdim_len  !< z dimension size

  !> @todo JLM: I dont love this...
  if (.not. self%lsm_active()) return

  dx = self%lsm%dx
  dy = self%lsm%dy
  xdim_len = self%lsm%xdim_len
  ydim_len = self%lsm%ydim_len
  zdim_len = self%lsm%zdim_len
end subroutine get_lsm_info

!> Get lsm geoval levels
subroutine get_geoval_levels(self, vars, nvars, nlevels)

  use, intrinsic :: iso_c_binding, only: c_size_t
  use oops_variables_mod,          only: oops_variables

  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self  !< geom object
  type(oops_variables),      intent(in) :: vars
  integer(c_size_t),         intent(in) :: nvars
  integer(c_size_t),      intent(inout) :: nlevels(nvars)
  character(len=*),           parameter :: myname = &
                                  & "geometry_mod:get_geoval_levels"
  ! local variables
  integer :: ivar

  call fckit_log%debug(myname // ' : start')

  nlevels = 0
  do ivar = 1, nvars

    select case (vars%variable(ivar))

      case ( "SNICE", "SNLIQ" )

        nlevels(ivar) = 3
        call fckit_log%debug("Found "//trim(myname)//":" &
                                   & //trim(vars%variable(ivar)))

      case ( "SNEQV", "SNOWH", "swe", "snow_depth", "LAI")

        nlevels(ivar) = 1
        call fckit_log%debug("Found "//trim(myname)//":" &
                                   & //trim(vars%variable(ivar)))

      case default
        call abor1_ftn(trim(myname)//":"//trim(vars%variable(ivar)) &
                                 & //" not found")
    end select
  end do

  call fckit_log%debug(myname // ' : end')

end subroutine get_geoval_levels


!-----------------------------------------------------------------------------
! Streamflow section


!> Check if the streamflow geometry object is allocated/active.
function stream_active(self) result(is_active)
  class(wrf_hydro_nwm_jedi_geometry), intent(in) :: self  !< the full geom object
  logical :: is_active

  is_active = allocated(self%stream)
end function stream_active


!> Init the streamflow geometry object
subroutine wrf_hydro_nwm_jedi_stream_geometry_init(self, f_conf, ncid)
  ! stream component IS optional
  class(wrf_hydro_nwm_jedi_geometry), intent(inout) :: self  !< the full geom object
  type(fckit_configuration), intent(in) :: f_conf  !< the yaml file
  integer, intent(in) :: ncid !< the geom file netcdf id

  integer :: ierr, dimid
  character(len=512) :: src_file_name

  ! Check LSM active
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "stream_src_file", src_file_name)
  call error_handler(ierr, "STOP:  stream_src_file attribute not defined in geometry file")
  if (ierr .ne. NF90_NOERR) then
     return
  end if

  allocate(self%stream)
  
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


!> Get the stream nearest neighbor from a lat and lon pair.
! @todo change indices to have dimension?
subroutine wrf_hydro_nwm_jedi_geometry_get_stream_nn( &
     self, lat, lon, ind)
  class(wrf_hydro_nwm_jedi_geometry), intent(in) :: self  !< geom object
  real, intent(in) :: lat, lon  !< the lat, lon of interest
  type(indices), intent(out) :: ind  !< the index pair of the nearest neighbor

  real, dimension(1) :: minimum
  real, allocatable  :: diff_lat(:), diff_lon(:), l2_norm(:)

  if (.not. self%stream_active()) return

  allocate(diff_lat, source=self%stream%lat)
  allocate(diff_lon, source=self%stream%lon)
  allocate(l2_norm,  source=self%stream%lon)

  diff_lat = diff_lat - lat
  diff_lon = diff_lon - lon
  l2_norm = sqrt( diff_lon**2 + diff_lat**2 )
  minimum = minloc(l2_norm)
  ind%ind_x = minimum(1)

  deallocate(l2_norm)
  deallocate(diff_lat)
  deallocate(diff_lon)  
end subroutine wrf_hydro_nwm_jedi_geometry_get_stream_nn


! !> Get stream info
! subroutine wrf_hydro_nwm_jedi_geometry_get_stream_info( &
!      self, &
!      dx, dy, &
!      xdim_len, ydim_len, zdim_len)
!   class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self  !< geom object
!   real, intent(out) :: dx, dy
!   integer, intent(out) :: xdim_len, ydim_len, zdim_len

!   dx = self%stream%dx
!   dy = self%stream%dy
!   xdim_len = self%stream%xdim_len
!   ydim_len = self%stream%ydim_len
!   zdim_len = self%stream%zdim_len 
! end subroutine wrf_hydro_nwm_jedi_geometry_get_lsm_info


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


end module wrf_hydro_nwm_jedi_geometry_mod
