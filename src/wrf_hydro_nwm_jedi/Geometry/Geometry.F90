module wrf_hydro_nwm_jedi_geometry_mod

  use fckit_configuration_module, only: fckit_configuration
use netcdf

implicit none
private

public :: wrf_hydro_nwm_jedi_geometry, indices, error_handler

!------------------------------------------------------------------------------

type wrf_hydro_nwm_lsm_geometry
   integer :: xstart, ystart
   integer :: xend, yend ! same as xdim_len and ydim_len
   integer :: xdim_len, ydim_len
   integer :: npz, snow_layers  ! TODO JLM: I dont like npz
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
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_init(self, f_conf)
  class(wrf_hydro_nwm_jedi_geometry),   intent(out) :: self
  type(fckit_configuration),  intent(in) :: f_conf

  ! LSM init required
  call wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf)
  ! Streamflow optional TODO JLM how to signal?
  
end subroutine wrf_hydro_nwm_jedi_geometry_init

  
subroutine wrf_hydro_nwm_jedi_lsm_geometry_init(self, f_conf)
  ! LSM component is not optional
  class(wrf_hydro_nwm_jedi_geometry),   intent(inout) :: self
  type(fckit_configuration),  intent(in) :: f_conf
  character(512) :: wrfinput_flnm
  character(len=:), allocatable :: str
  integer :: ierr, ncid, dimid

  call f_conf%get_or_die("input_file",str)
  wrfinput_flnm = str
  write(*,*) 'Input file for Geometry: ', wrfinput_flnm
  deallocate(str)

  self%lsm%xstart = 1
  self%lsm%ystart = 1
  self%lsm%snow_layers = 3  ! This is max, hard-coded.

  ierr = 0
  write(*,*) 'Reading from NETCDF file the Geometry'
  ierr = nf90_open(wrfinput_flnm, NF90_NOWRITE, ncid)
  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn("ERROR: wrfinput file required by Geometry not found.")
  end if
  
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "DX", self%lsm%dx)
  call error_handler(ierr, "STOP:  Problem finding DX attribute")
  ierr = nf90_get_att(ncid, NF90_GLOBAL, "DY", self%lsm%dy)
  call error_handler(ierr, "STOP:  Problem finding DY attribute")  

  ierr = nf90_inq_dimid(ncid, "west_east", dimid)
  call error_handler(ierr, "STOP:  Problem finding west_east dimension")
  ierr = nf90_inquire_dimension(ncid, dimid, len=self%lsm%xdim_len)
  call error_handler(ierr, "STOP:  Problem getting west_east dimension length")
  self%lsm%xend = self%lsm%xdim_len

  ierr = nf90_inq_dimid(ncid, "south_north", dimid)
  call error_handler(ierr, "STOP:  Problem finding south_north dimension")
  ierr = nf90_inquire_dimension(ncid, dimid, len=self%lsm%ydim_len)
  call error_handler(ierr, "STOP:  Problem getting south_north dimension length")
  self%lsm%yend = self%lsm%ydim_len

  ! This is unused currently... 
  ierr = nf90_inq_dimid(ncid, "soil_layers_stag", self%lsm%npz)
  call error_handler(ierr, "STOP:  Problem finding soil_layers_stag dimension")
  
  allocate( &
       self%lsm%lat(self%lsm%xdim_len, self%lsm%ydim_len), &
       self%lsm%lon(self%lsm%xdim_len, self%lsm%ydim_len))
  
  call get_2d("XLAT", ncid, self%lsm%lat, self%lsm%xdim_len, self%lsm%ydim_len)
  call get_2d("XLONG", ncid, self%lsm%lon, self%lsm%xdim_len, self%lsm%ydim_len)
  
  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn("ERROR: dx and/or dy need by Geometry not found")
  end if
  
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

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: other

  self%lsm%dx = other%lsm%dx
  self%lsm%dy = other%lsm%dy
  self%lsm%xdim_len = other%lsm%xdim_len
  self%lsm%ydim_len = other%lsm%ydim_len
  self%lsm%npz = other%lsm%npz
  self%lsm%snow_layers = other%lsm%snow_layers
  self%lsm%xstart = other%lsm%xstart
  self%lsm%ystart = other%lsm%ystart
  self%lsm%xend = other%lsm%xend
  self%lsm%yend = other%lsm%yend
  allocate(self%lsm%lat, source = other%lsm%lat) 
  allocate(self%lsm%lon, source = other%lsm%lon)
  
end subroutine
  
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_delete(self)
  class(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: self

  deallocate(self%lsm%lat)
  deallocate(self%lsm%lon)

end subroutine

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_get_info(self,dx_,dy_,xdim_len,ydim_len,npz_)
  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self
  real, intent(out) :: dx_,dy_
  integer, intent(out) :: xdim_len,ydim_len,npz_

  dx_ = self%lsm%dx
  dy_ = self%lsm%dy
  xdim_len = self%lsm%xdim_len
  ydim_len = self%lsm%ydim_len
  npz_ = self%lsm%npz  

end subroutine wrf_hydro_nwm_jedi_geometry_get_info

!------------------------------------------------------------------------------

!Subroutine copied from module_wrfinputfile.F
subroutine get_2d(name, ncid, array, idim, jdim)
  !
  ! From the NetCDF unit <ncid>, read the variable named <name> with  
  ! dimensions <idim> and <jdim>, filling the pre-dimensioned array <array>
  !
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
