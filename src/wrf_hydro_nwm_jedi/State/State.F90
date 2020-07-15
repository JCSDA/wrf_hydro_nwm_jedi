! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_state_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use fckit_mpi_module
use oops_variables_mod

use wrf_hydro_nwm_jedi_field_mod, only: wrf_hydro_nwm_jedi_fields,checksame
! use fv3jedi_constants_mod,       only: rad2deg, constoz
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry, error_handler
! use fv3jedi_increment_utils_mod, only: fv3jedi_increment
! use fv3jedi_interpolation_mod,   only: field2field_interp
! use fv3jedi_kinds_mod,           only: kind_real
use iso_c_binding, only : c_float
!use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
!use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
use wrf_hydro_nwm_jedi_state_utils_mod, only: wrf_hydro_nwm_jedi_state
!use fv3jedi_getvalues_mod,       only: getvalues
use netcdf
!use wind_vt_mod, only: a2d

!use mpp_domains_mod,             only: east, north, center

implicit none
private

public :: wrf_hydro_nwm_jedi_state, create, delete, zeros, copy, axpy,&
     create_from_other,&!add_incr, &
     read_file, get_mean_stddev, &! write_file, gpnorm, rms, &
     change_resol, state_print !getvalues, analytic_IC, state_print

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, vars)

  implicit none
  type(wrf_hydro_nwm_jedi_state),  intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry),   intent(in)    :: geom
  type(oops_variables), intent(in)    :: vars

  self%nf = vars%nvars()
  call self%fields_obj%create(geom,vars)
  
end subroutine create

! ------------------------------------------------------------------------------

subroutine create_from_other(self, other)

  ! Passed variables
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self  !< Fields
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: other !< Other fields

  ! Create new state from other state
  self = other

  ! Initialize all arrays to zero
  !call zeros(self)

end subroutine create_from_other

! ------------------------------------------------------------------------------

subroutine delete(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self

  call self%fields_obj%deallocate_field()

end subroutine delete

! ------------------------------------------------------------------------------

function get_n_fields(self) result(nf)
  implicit none
  integer :: nf
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self

  nf = self%nf

end function get_n_fields

! ------------------------------------------------------------------------------

subroutine zeros(self)

implicit none
type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
integer :: var

! do var = 1, self%nf
!   self%fields(var)%array = 0.d0
! enddo

end subroutine zeros

! ------------------------------------------------------------------------------

subroutine copy(self,rhs)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state), intent(in)    :: rhs

  integer :: var

  call checksame(self%fields_obj,rhs%fields_obj,"wrf_hydro_nwm_jedi_state_mod.copy")
  
  ! write(*,*) "Copying field self%nf= ",self%nf

  ! Deep copy
  self%fields_obj = rhs%fields_obj

  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init

end subroutine copy

! ------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)

implicit none
type(wrf_hydro_nwm_jedi_state),  intent(inout) :: self
real(kind=c_float), intent(in)    :: zz
type(wrf_hydro_nwm_jedi_state),  intent(in)    :: rhs

! integer :: var

! call checksame(self%fields,rhs%fields,"wrf_hydro_nwm_jedi_state_mod.axpy")

! do var = 1, self%nf
!   self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
! enddo

end subroutine axpy

! ------------------------------------------------------------------------------

subroutine change_resol(self,geom,rhs,geom_rhs)

implicit none
type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
type(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: geom
type(wrf_hydro_nwm_jedi_state), intent(in)    :: rhs
type(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: geom_rhs

! ! Interpolation
! integer :: var
! type(field2field_interp) :: interp
! logical :: integer_interp = .false.

! call checksame(self%fields,rhs%fields,"wrf_hydro_nwm_jedi_state_mod.change_resol")

! if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1) == 0) then

!   call copy(self, rhs)

! else

!   ! Check if integer interp needed
!   do var = 1, self%nf
!     if (rhs%fields(var)%integerfield) integer_interp = .true.
!   enddo

!   call interp%create(geom%interp_method, integer_interp, geom_rhs, geom)
!   call interp%apply(self%nf, geom_rhs, rhs%fields, geom, self%fields)
!   call interp%delete()

!   self%calendar_type = rhs%calendar_type
!   self%date_init = rhs%date_init

! endif

end subroutine change_resol

! ------------------------------------------------------------------------------

subroutine read_file(geom, self, c_conf, vdate)
  use string_utils

  implicit none

  type(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: geom     !< Geometry
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self     !< State
  type(c_ptr),         intent(in)    :: c_conf   !< Configuration
  type(datetime),      intent(inout) :: vdate    !< DateTime
  
  character(len=10) :: filetype
  character(len=255) :: filename
  integer :: flipvert
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str
  
  integer :: ixfull, jxfull, var
  
  ! Fortran configuration
  ! ---------------------
  f_conf = fckit_configuration(c_conf)
  
  call f_conf%get_or_die("model_filename",str)
  filename = str
  deallocate(str)
  
  ! ixfull = geom%xend-geom%xstart+1
  ! jxfull = geom%yend-geom%ystart+1

  call self%fields_obj%read_fields_from_file(filename,geom%xstart,geom%xend,geom%ystart,geom%yend)

end subroutine read_file

! ------------------------------------------------------------------------------

subroutine state_print(self,string)
  use iso_c_binding, only : c_null_char
 implicit none
 type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
 character(len=1,kind=c_char) :: string(8192)
 character(len=8192) :: tmp_str
 integer :: s_len,i
! call self%fields_obj%print_single_field("swe",tmp_str)
 call self%fields_obj%print_all_fields(tmp_str)

 s_len = len_trim(tmp_str)

 do i = 1, s_len
    string(i:i) = tmp_str(i:i)
 end do
 string(s_len+1) = c_null_char

end subroutine state_print

! ------------------------------------------------------------------------------

subroutine get_mean_stddev(self,nf,pstat)

 implicit none
 type(wrf_hydro_nwm_jedi_state), intent(in) :: self
 integer :: nf
 real(kind=c_float), intent(inout) ::  pstat(3,nf)

 ! pstat(1,1) = 0.0
 ! pstat(2,1) = 0.0
 ! pstat(3,1) = 0.0

 ! call self%fields(1)%mean_stddev(pstat(1,1),pstat(2,1),pstat(3,1))
 
end subroutine get_mean_stddev

! ------------------------------------------------------------------------------

subroutine gpnorm(self, nf, pstat)

 implicit none
 type(wrf_hydro_nwm_jedi_state),  intent(in)    :: self
 integer,              intent(in)    :: nf
 real(kind=c_float), intent(inout) :: pstat(3, nf)

! if (nf .ne. self%nf) then
!   call abor1_ftn("wrf_hydro_nwm_jedi_state: gpnorm | nf passed in does not match expeted nf")
! endif

! call fields_gpnorm(nf, self%fields, pstat, self%f_comm)

end subroutine gpnorm

! ------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(wrf_hydro_nwm_jedi_state),  intent(in)  :: self
real(kind=c_float), intent(out) :: prms

! call fields_rms(self%nf, self%fields, prms, self%f_comm)

end subroutine rms

! ------------------------------------------------------------------------------

subroutine get_from_restart_3d(restart_filename_remember, parallel_xstart, parallel_xend, subwindow_xstart, ixfull, jxfull, name, array, return_error)
  implicit none
    character(len=*), intent(in) :: restart_filename_remember
    integer,                            intent(in) :: parallel_xstart
    integer,                            intent(in) :: parallel_xend
    integer,                            intent(in) :: subwindow_xstart
    integer,                            intent(in) :: ixfull
    integer,                            intent(in) :: jxfull
    character(len=*),                   intent(in)  :: name
    real(kind=c_float),             dimension(:,:,:), intent(out) :: array
    integer,          optional,         intent(out) :: return_error
    integer :: ierr
    integer :: ncid
    integer :: varid
    integer, dimension(4) :: nstart
    integer, dimension(4) :: ncount
! #ifdef _PARALLEL_
!     ierr = nf90_open_par(trim(restart_filename_remember), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
! #else
    ierr = nf90_open(trim(restart_filename_remember), NF90_NOWRITE, ncid)
    !#endif
    call error_handler(ierr, "GET_FROM_RESTART: Problem opening restart file '"//trim(restart_filename_remember)//"'")
    nstart = (/parallel_xstart-subwindow_xstart+1,1, 1, 1/)
    ncount = (/parallel_xend-parallel_xstart+1, size(array,2), size(array,3), 1/)
    if (present(return_error)) then
       ierr = nf90_inq_varid(ncid, name, varid)
       if (ierr == NF90_NOERR) then
          return_error = 0
          call error_handler(ierr, "Problem finding variable in restart file '"//trim(name)//"'")
          ierr = nf90_get_var(ncid, varid, array, start=nstart(1:4))
          call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
       else
          return_error = 1
          write(*,'("Did not find optional variable ''",A,"'' in restart file ''", A, "''")') trim(name), trim(restart_filename_remember)
       endif
    else
       ierr = nf90_inq_varid(ncid, name, varid)
       call error_handler(ierr, "Problem finding required variable in restart file: '"//trim(name)//"'")
       ierr = nf90_get_var(ncid, varid, array, start=nstart(1:4))
       call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
    endif
    ierr = nf90_close(ncid)
    call error_handler(ierr, "Problem closing restart file")
    
  end subroutine get_from_restart_3d

  subroutine get_from_restart_2d_float(restart_filename_remember, parallel_xstart, parallel_xend, subwindow_xstart, ixfull, jxfull, name, array, return_error)
    implicit none
    character(len=*), intent(in) :: restart_filename_remember
    integer,                            intent(in) :: parallel_xstart
    integer,                            intent(in) :: parallel_xend
    integer,                            intent(in) :: subwindow_xstart
    integer,                            intent(in) :: ixfull
    integer,                            intent(in) :: jxfull
    character(len=*),                   intent(in)  :: name
    real,             dimension(parallel_xstart:parallel_xend,jxfull), intent(out) :: array
    integer,          optional,         intent(out) :: return_error
    integer :: ierr
    integer :: ncid
    integer :: varid
    integer, dimension(4) :: nstart
    integer, dimension(4) :: ncount
    integer :: rank
! #ifdef _PARALLEL_
!     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
!     if (ierr /= MPI_SUCCESS) stop "MPI_COMM_RANK"
!     ierr = nf90_open_par(trim(restart_filename_remember), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
! #else
    rank = 0
    ierr = nf90_open(trim(restart_filename_remember), NF90_NOWRITE, ncid)
!#endif
    call error_handler(ierr, "GET_FROM_RESTART: Problem opening restart file '"//trim(restart_filename_remember)//"'")
    nstart = (/ parallel_xstart-subwindow_xstart+1, 1,  1, -99999 /)
    ncount = (/ parallel_xend-parallel_xstart+1,   size(array,2),  1, -99999 /)
    if (present(return_error)) then
       ierr = nf90_inq_varid(ncid, name, varid)
       if (ierr == NF90_NOERR) then
          return_error = 0
          call error_handler(ierr, "Problem finding variable in restart file '"//trim(name)//"'")
          ierr = nf90_get_var(ncid, varid, array, start=nstart(1:3))
          call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
       else
          return_error = 1
          if (rank == 0) write(*,'("Did not find optional variable ''",A,"'' in restart file ''", A, "''")') trim(name), trim(restart_filename_remember)
       endif
    else
       ierr = nf90_inq_varid(ncid, name, varid)
       call error_handler(ierr, "Problem finding required variable in restart file: '"//trim(name)//"'")
       ierr = nf90_get_var(ncid, varid, array, start=nstart(1:3))
       call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
    endif
    ierr = nf90_close(ncid)
    call error_handler(ierr, "Problem closing restart file")
    
  end subroutine get_from_restart_2d_float


end module wrf_hydro_nwm_jedi_state_mod
