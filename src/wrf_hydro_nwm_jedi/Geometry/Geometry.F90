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
   integer :: idim, jdim, npz
   real    :: dx, dy
   real, allocatable :: xlat(:,:), xlong(:,:)
 contains
   procedure :: init   => wrf_hydro_nwm_jedi_geometry_init
   procedure :: clone  => wrf_hydro_nwm_jedi_geometry_clone
   procedure :: delete => wrf_hydro_nwm_jedi_geometry_delete
   procedure :: get_info => wrf_hydro_nwm_jedi_geometry_get_info
   procedure :: coo_to_grid => wrf_hydro_nwm_jedi_geometry_coo_to_grid
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
  integer :: ierr,ncid,dimid

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

  ierr = nf90_inq_dimid(ncid, "west_east", dimid)
  call error_handler(ierr, "STOP:  Problem finding west_east dimension")
  ierr = nf90_inquire_dimension(ncid, dimid, len=self%idim)
  call error_handler(ierr, "STOP:  Problem finding west_east dimension")

  ierr = nf90_inq_dimid(ncid, "south_north", dimid)
  call error_handler(ierr, "STOP:  Problem finding south_north dimension")
  ierr = nf90_inquire_dimension(ncid, dimid, len=self%jdim)
  call error_handler(ierr, "STOP:  Problem finding south_north dimension")
  
  ierr = nf90_inq_dimid(ncid,"soil_layers_stag", self%npz)

  allocate(self%xlat(self%idim, self%jdim), self%xlong(self%idim, self%jdim))
  call get_2d("XLAT", ncid, self%xlat, self%idim, self%jdim)
  call get_2d("XLONG", ncid, self%xlong, self%idim, self%jdim)
  
  if(ierr /= 0) then
     write(*,*) ierr
     call abor1_ftn("ERROR: dx and/or dy need by Geometry not found")
  end if

  ! For testing purposes
  block 
    integer:: x,y
    real :: lat, long

    lat = self%xlat(10,10)
    long = self%xlong(12,12)

    write(*,*) "lat long requested: ",lat,long
    
    call wrf_hydro_nwm_jedi_geometry_coo_to_grid(self, lat, long, x, y)

    write(*,*) "X found: ",x," Y found: ", y

    write(*,*) "Value in those positions: lat:",self%xlat(y,x)," long:", self%xlong(y,x)

  end block
  
end subroutine wrf_hydro_nwm_jedi_geometry_init

subroutine wrf_hydro_nwm_jedi_geometry_coo_to_grid(self, lat, long, x, y)
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: self
  real, intent(in) :: lat, long
  real,dimension(2) :: minimum
  real, allocatable :: diff_lat(:,:), diff_long(:,:), l2_norm(:,:)
  integer, intent(out) :: x, y
  integer :: i, j

  allocate(diff_lat, source=self%xlat)
  allocate(diff_long, source=self%xlong)
  allocate(l2_norm, source=self%xlong)
  
  diff_lat = diff_lat - lat
  diff_long = diff_long - long

  l2_norm = sqrt( diff_long**2 + diff_lat**2 )

  minimum = minloc(l2_norm)

  x = minimum(2); y = minimum(1)

  deallocate(l2_norm)
  deallocate(diff_lat)
  deallocate(diff_long)
  
end subroutine wrf_hydro_nwm_jedi_geometry_coo_to_grid

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm_jedi_geometry),  intent(out) :: self
  class(wrf_hydro_nwm_jedi_geometry),   intent(in) :: other

  self%dx = other%dx
  self%dy = other%dy
  self%idim = other%idim
  self%jdim = other%jdim
  self%npz = other%npz
  allocate(self%xlat, source = other%xlat) 
  allocate(self%xlong, source = other%xlong)
  
end subroutine
  
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_delete(self)
  class(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: self

  deallocate(self%xlat)
  deallocate(self%xlong)

end subroutine

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_get_info(self,dx_,dy_,idim_,jdim_,npz_)
  class(wrf_hydro_nwm_jedi_geometry),  intent(in) :: self
  real, intent(out) :: dx_,dy_
  integer, intent(out) :: idim_,jdim_,npz_

  dx_ = self%dx
  dy_ = self%dy
  idim_ = self%idim
  jdim_ = self%jdim
  npz_ = self%npz  

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
  call error_handler(ierr, "STOP:  READ_WRFINPUT: Problem finding variable '"//trim(name)//"' in wrfinput file.")

  ierr = nf90_get_var(ncid, varid, array)
  call error_handler(ierr, "STOP:  READ_WRFINPUT:  Problem retrieving variable '"//trim(name)//"' from wrfinput file.")

end subroutine get_2d

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
