! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> State (a fields/variable manager) for wrf_hydro_nwm - jedi.
module wrf_hydro_nwm_jedi_state_mod

use iso_c_binding
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
use iso_c_binding, only : c_float
!use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
!use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
!use fv3jedi_getvalues_mod,       only: getvalues
use netcdf
!use wind_vt_mod, only: a2d

!use mpp_domains_mod,             only: east, north, center

implicit none
private

public :: wrf_hydro_nwm_jedi_state, create, delete, zeros, copy, axpy,&
     create_from_other, add_incr, &
     read_state_from_file, get_mean_stddev, write_state_to_file, &!gpnorm, rms, &
     change_resol, state_print !getvalues, analytic_IC, state_print


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
  type(wrf_hydro_nwm_jedi_fields),allocatable :: fields_obj
end type wrf_hydro_nwm_jedi_state


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


!> The init method for state. It calles all fields
subroutine create(self, geom, vars)
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in)    :: geom
  type(oops_variables),              intent(in)    :: vars

  self%nf = vars%nvars()
  allocate(self%fields_obj)
  call self%fields_obj%create(geom, vars)
  
end subroutine create

! ------------------------------------------------------------------------------

subroutine create_from_other(self, other)

  ! Passed variables
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self  !< Fields
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: other !< Other fields

  ! Create new state from other state
  write(*,*) "Calling create from other"
  self = other

  ! Initialize all arrays to zero
  !call zeros(self)

end subroutine create_from_other

subroutine delete(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self

  call self%fields_obj%deallocate_field()
end subroutine delete


function get_n_fields(self) result(nf)
  implicit none
  integer :: nf
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self

  nf = self%nf

end function get_n_fields


subroutine zeros(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  integer :: var
  
  ! do var = 1, self%nf
  !   self%fields(var)%array = 0.d0
  ! enddo
end subroutine zeros


subroutine copy(self, rhs)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state), intent(in)    :: rhs

  integer :: var

  call checksame( &
       self%fields_obj,rhs%fields_obj, &
       "wrf_hydro_nwm_jedi_state_mod.copy")
  
  ! Deep copy
  self%fields_obj = rhs%fields_obj
  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init

  ! write(*,*) "Checking for deep copy in copy"
  ! call rhs%fields_obj%print_all_fields()
  ! call self%fields_obj%print_all_fields()
  
end subroutine copy


subroutine axpy(self, zz, rhs)
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


!> Read the state from files (via read_fields_from_file)
subroutine read_state_from_file(self, c_conf)
  use string_utils
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self   !< State
  type(c_ptr),                       intent(in)    :: c_conf !< Configuration
  
  character(len=10) :: filetype
  character(len=255) :: filename_lsm, filename_hydro
  integer :: flipvert
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str
  integer :: ixfull, jxfull, var

  f_conf = fckit_configuration(c_conf)

  call f_conf%get_or_die("filename_lsm", str)
  filename_lsm = str
  deallocate(str)

  call f_conf%get_or_die("filename_hydro", str)
  filename_hydro = str
  deallocate(str)
  
  call self%fields_obj%read_fields_from_file(filename_lsm, filename_hydro)
end subroutine read_state_from_file

subroutine write_state_to_file(self, c_conf, vdate)
  use string_utils
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self   !< State
  type(c_ptr),                       intent(in)    :: c_conf !< Configuration
  type(datetime), intent(inout) :: vdate
  
  character(len=10) :: filetype
  character(len=255) :: filename_lsm, filename_hydro
  character(len=255) :: validitydate
  integer :: flipvert
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str
  integer :: ixfull, jxfull, var

  f_conf = fckit_configuration(c_conf)

  call datetime_to_string(vdate, validitydate)
  
  filename_lsm = "lsm."//validitydate
  filename_hydro = "hydro."//validitydate

  ! call f_conf%get_or_die("filename_lsm", str)
  ! filename_lsm = str
  ! deallocate(str)

  ! call f_conf%get_or_die("filename_hydro", str)
  ! filename_hydro = str
  ! deallocate(str)
  write(*,*) "State written in ",filename_lsm,filename_hydro
  !  call self%fields_obj%read_fields_from_file(filename_lsm, filename_hydro)
  !call self%fields_obj%write_state_to_file(trim(filename_lsm),trim(filename_hydro))
  
end subroutine write_state_to_file

  function genfilename (c_conf,length,vdate)

    use iso_c_binding
    use datetime_mod
    use duration_mod

    type(c_ptr),    intent(in) :: c_conf  !< Configuration
    integer,        intent(in) :: length
    type(datetime), intent(in) :: vdate

    type(fckit_configuration)     :: f_conf
    character(len=:), allocatable :: str
    character(len=length)         :: genfilename
    character(len=length)         :: fdbdir, expver, typ, validitydate, referencedate, sstep
    character(len=length)         :: prefix, mmb
    type(datetime)                :: rdate
    type(duration)                :: step
    integer                       :: lenfn

    f_conf = fckit_configuration(c_conf)

    ! here we should query the length and then allocate "string".
    ! But Fortran 90 does not allow variable-length allocatable strings.
    call f_conf%get_or_die("datadir", str)
    fdbdir = trim(str)
    deallocate(str)
    call f_conf%get_or_die("exp", str)
    expver = trim(str)
    deallocate(str)
    call f_conf%get_or_die("type", str)
    typ = trim(str)
    deallocate(str)

    if (typ=="ens") then
      call f_conf%get_or_die("member", str)
      mmb = trim(str)
      deallocate(str)

      lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
      prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
    else
      lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
      prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
    endif

    if (typ=="fc" .or. typ=="ens") then
      call f_conf%get_or_die("date", str)
      referencedate = trim(str)
      deallocate(str)
      call datetime_to_string(vdate, validitydate)
      call datetime_create(TRIM(referencedate), rdate)
      call datetime_diff(vdate, rdate, step)
      call duration_to_string(step, sstep)
      lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
      genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep) // ".nc"
    endif

    if (typ=="an") then
      call datetime_to_string(vdate, validitydate)
      lenfn = lenfn + 1 + LEN_TRIM(validitydate)
      genfilename = TRIM(prefix) // "." // TRIM(validitydate) // ".nc"
    endif

    if (lenfn>length) &
      & call abor1_ftn("sw_state_mod:genfilename: filename too long")

  end function genfilename

!> Print the state
subroutine state_print(self, string)
  use iso_c_binding, only : c_null_char
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in)  :: self  !< State
  character(len=1, kind=c_char),  intent(out) :: string(8192) !< The output string

  character(len=8192) :: tmp_str
  integer :: s_len, i

  ! call self%fields_obj%print_single_field("swe",tmp_str)
  call self%fields_obj%print_all_fields(tmp_str)
  s_len = len_trim(tmp_str)
  do i = 1, s_len
     string(i:i) = tmp_str(i:i)
  end do
  string(s_len+1) = c_null_char
end subroutine state_print


subroutine get_mean_stddev(self, nf, pstat)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in)    :: self
  real(kind=c_float),             intent(inout) :: pstat(3, nf)

  integer :: nf

  ! pstat(1,1) = 0.0
  ! pstat(2,1) = 0.0
  ! pstat(3,1) = 0.0
  ! call self%fields(1)%mean_stddev(pstat(1,1), pstat(2,1), pstat(3,1))
end subroutine get_mean_stddev


subroutine gpnorm(self, nf, pstat)
  implicit none
  type(wrf_hydro_nwm_jedi_state),  intent(in)    :: self
  integer,                         intent(in)    :: nf
  real(kind=c_float),              intent(inout) :: pstat(3, nf)

  ! if (nf .ne. self%nf) then
  !   call abor1_ftn("wrf_hydro_nwm_jedi_state: gpnorm | nf passed in does not match expeted nf")
  ! endif

  ! call fields_gpnorm(nf, self%fields, pstat, self%f_comm)
end subroutine gpnorm


subroutine rms(self, prms)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in)  :: self
  real(kind=c_float),             intent(out) :: prms

  ! call fields_rms(self%nf, self%fields, prms, self%f_comm)
end subroutine rms


end module wrf_hydro_nwm_jedi_state_mod
