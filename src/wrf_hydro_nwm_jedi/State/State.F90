! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State (a fields/variable manager) for wrf_hydro_nwm - jedi.
module wrf_hydro_nwm_jedi_state_mod

use iso_c_binding, only: c_char, c_float, c_ptr, c_null_char
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use fckit_mpi_module
use oops_variables_mod

use wrf_hydro_nwm_jedi_fields_mod,   only: wrf_hydro_nwm_jedi_fields, checksame
! use fv3jedi_constants_mod,       only: rad2deg, constoz
use wrf_hydro_nwm_jedi_geometry_mod,only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_util_mod,    only: error_handler, indices
! use fv3jedi_increment_utils_mod, only: fv3jedi_increment
! use fv3jedi_interpolation_mod,   only: field2field_interp
! use fv3jedi_kinds_mod,           only: kind_real
!use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
!use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
!use fv3jedi_getvalues_mod,       only: getvalues
use netcdf
!use wind_vt_mod, only: a2d

!use mpp_domains_mod,             only: east, north, center

implicit none
private

public :: &
     wrf_hydro_nwm_jedi_state, &
     create, &
     delete, &
     zeros, &
     ones,  &
     copy,  &
     create_from_other, &
     add_incr, &
     read_state_from_file,  &
     write_state_to_file, &
     change_resol, &
     state_print, &
     axpy
     ! get_mean_stddev, &
     ! rms !&
     ! gpnorm, &
     ! getvalues, &
     ! analytic_IC,


!> Fortran mirror of C class
type, public :: wrf_hydro_nwm_jedi_state
  !Local copies of grid for convenience
!  integer :: isc, iec, jsc, jec
  integer :: npx, npy, npz
  integer :: ntiles, ntile
  logical :: hydrostatic = .true.
  integer :: calendar_type, date_init(6) !Read/write for GFS
  integer :: nf
  logical :: have_agrid
  logical :: have_dgrid
  type(fckit_mpi_comm) :: f_comm
  type(wrf_hydro_nwm_jedi_fields), allocatable :: fields_obj
end type wrf_hydro_nwm_jedi_state


contains


!> init method for state. It calls all fields
subroutine create(self, geom, vars)
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in)    :: geom
  type(oops_variables),              intent(in)    :: vars
  ! type(datetime),                    intent(inout) :: time

  self%nf = vars%nvars()
  allocate(self%fields_obj)
  call self%fields_obj%create(geom, vars)
  call self%fields_obj%zero()
end subroutine create


subroutine create_from_other(self, other)
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self  !< Self State object
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: other !< Other State object

  self = other
end subroutine create_from_other


subroutine delete(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self

  call self%fields_obj%deallocate_field()
end subroutine delete


function get_n_fields(self) result(nf)
  implicit none
  integer :: nf
  type(wrf_hydro_nwm_jedi_state), intent(in) :: self

  nf = self%nf
end function get_n_fields


subroutine zeros(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  integer :: var

  call self%fields_obj%zero()
end subroutine zeros


subroutine ones(self)
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  call self%fields_obj%one()
end subroutine ones


subroutine copy(self, rhs)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state), intent(in   ) :: rhs

  integer :: var

  call checksame( &
       self%fields_obj, rhs%fields_obj, &
       "wrf_hydro_nwm_jedi_state_mod.copy")

  ! Deep copy
  self%fields_obj = rhs%fields_obj
  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init

  ! write(*,*) "Checking for deep copy in copy"
  ! call rhs%fields_obj%print_all_fields()
  ! call self%fields_obj%print_all_fields()
end subroutine copy


subroutine add_incr(self, rhs)
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state),    intent(in)    :: rhs

  call self%fields_obj%add_increment(rhs%fields_obj)
end subroutine add_incr


subroutine change_resol(self, geom, rhs, geom_rhs)
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(inout) :: geom
  type(wrf_hydro_nwm_jedi_state),    intent(in)    :: rhs
  type(wrf_hydro_nwm_jedi_geometry), intent(inout) :: geom_rhs

  ! Interpolation
  integer :: var
  !type(field2field_interp) :: interp
  !logical :: integer_interp = .false.

  ! call checksame(self%fields_obj,rhs%fields_obj,"wrf_hydro_nwm_jedi_state_mod.change_resol")

  ! if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1) == 0) then

  !   call copy(self, rhs)

  ! else

  !   ! Check if integer interp needed
  !   do var = 1, self%nf
  !     if (rhs%fields_obj(var)%integerfield) integer_interp = .true.
  !   enddo

  !   call interp%create(geom%interp_method, integer_interp, geom_rhs, geom)
  !   call interp%apply(self%nf, geom_rhs, rhs%fields_obj, geom, self%fields_obj)
  !   call interp%delete()

  !   self%calendar_type = rhs%calendar_type
  !   self%date_init = rhs%date_init

  ! endif
end subroutine change_resol


!> Read the state from files (via read_fields_from_file)
subroutine read_state_from_file(self, c_conf, f_dt)
  use string_utils
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self   !< State
  type(c_ptr),                       intent(in)    :: c_conf !< Configuration
  type(datetime),                    intent(inout) :: f_dt   !< DateTime

  character(len=10) :: filetype
  character(len=255) :: filename_lsm, filename_hydro
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str
  character(len=30) :: fstring
  ! integer :: flipvert
  ! integer :: ixfull, jxfull, var

  f_conf = fckit_configuration(c_conf)

  call f_conf%get_or_die("filename_lsm", str)
  filename_lsm = str
  deallocate(str)

  call f_conf%get_or_die("filename_hydro", str)
  filename_hydro = str
  deallocate(str)

  call self%fields_obj%read_fields_from_file(filename_lsm, filename_hydro, f_dt)
  call datetime_to_string(f_dt, fstring)
  write(*,*) 'read_state_from_file f_dt to string: '//fstring
end subroutine read_state_from_file


subroutine write_state_to_file(self, c_conf, f_dt)
  use string_utils
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self   !< State
  type(c_ptr),                    intent(in)    :: c_conf !< Configuration
  type(datetime),                 intent(inout) :: f_dt

  character(len=10) :: filetype
  character(len=255) :: filename_lsm, filename_hydro
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str
  character(len=30) :: fstring
  ! integer :: flipvert
  ! integer :: ixfull, jxfull, var
 
  f_conf = fckit_configuration(c_conf)

  call f_conf%get_or_die("filename_lsm", str)
  filename_lsm = str
  deallocate(str)

  call f_conf%get_or_die("filename_hydro", str)
  filename_hydro = str
  deallocate(str)

  call self%fields_obj%write_fields_to_file(filename_lsm, filename_hydro, f_dt)
  call datetime_to_string(f_dt, fstring)
end subroutine write_state_to_file

!> Print state
subroutine state_print(self, string)
  use iso_c_binding, only : c_null_char, c_new_line
  implicit none
  type(wrf_hydro_nwm_jedi_state),           intent(in)  :: self         !< State
  character(len=1, kind=c_char),  optional, intent(out) :: string(8192) !< The output string

  character(len=1, kind=c_char) :: local_string(8192)

  if(present(string)) then
     call self%fields_obj%print_all_fields(string)
  else
     write(*,*) c_new_line
     !                      "Print State (C++) ------------------------ ";
     write(*,*) c_new_line//"Print State (Fortran) -------------------- "
     call self%fields_obj%print_all_fields()
     write(*,*) c_new_line//"End Print State (Fortran) ---------------- "
     write(*,*) c_new_line//c_new_line
  end if
end subroutine state_print


subroutine axpy(self, scalar, other_in)
  use iso_c_binding, only: c_float
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  real(kind=c_float),             intent(   in) :: scalar
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: other_in

  type(wrf_hydro_nwm_jedi_state) :: other

  call create_from_other(other, other_in)
  call other%fields_obj%scalar_mult(scalar)  ! = scalar * other
  call self%fields_obj%add_increment(other%fields_obj)  ! = self + (scalar*other)
end subroutine axpy


! subroutine get_mean_stddev(self, mean, stddev)
!   implicit none
!   type(wrf_hydro_nwm_jedi_state), intent(in ) :: self
!   real(kind=c_float),             intent(out) :: mean
!   real(kind=c_float),             intent(out) :: stddev

!   call self%fields_obj%mean_stddev())
! end subroutine get_mean_stddev


! subroutine rms(self, prms)
!   implicit none
!   type(wrf_hydro_nwm_jedi_state), intent(in)  :: self
!   real(kind=c_float),             intent(out) :: prms

!   ! call fields_rms(self%nf, self%fields_obj, prms, self%f_comm)
! end subroutine rms


! subroutine gpnorm(self, nf, pstat)
!   implicit none
!   type(wrf_hydro_nwm_jedi_state),  intent(in)    :: self
!   integer,                         intent(in)    :: nf
!   real(kind=c_float),              intent(inout) :: pstat(3, nf)

!   ! if (nf .ne. self%nf) then
!   !   call abor1_ftn("wrf_hydro_nwm_jedi_state: gpnorm | nf passed in does not match expeted nf")
!   ! endif

!   ! call fields_gpnorm(nf, self%fields_obj, pstat, self%f_comm)
! end subroutine gpnorm


end module wrf_hydro_nwm_jedi_state_mod
