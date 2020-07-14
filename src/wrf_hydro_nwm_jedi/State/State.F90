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

public :: wrf_hydro_nwm_jedi_state, create, delete, zeros, copy, axpy, &
     !add_incr, &
     read_state_from_file, get_mean_stddev, &! write_file, gpnorm, rms, &
     change_resol, state_print !getvalues, analytic_IC, state_print


! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

  
subroutine create(self, geom, vars)
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in)    :: geom
  type(oops_variables),              intent(in)    :: vars

  write(*,*) "state create"
  self%nf = vars%nvars()
  call self%fields_obj%create(geom, vars)
end subroutine create

! ------------------------------------------------------------------------------

! subroutine create_from_other(self, other)

!   ! Passed variables
!   type(wrf_hydro_nwm_jedi_state), intent(inout) :: self  !< Fields
!   type(wrf_hydro_nwm_jedi_state), intent(   in) :: other !< Other fields

!   ! Create new state from other state
!   self = shallow_water_state_type(other%get_geometry())

!   ! Initialize all arrays to zero
!   !call zeros(self)

! end subroutine create_from_other

subroutine delete(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self

  call self%fields_obj%deallocate_field()
end subroutine delete


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
  
  ! write(*,*) "Copying field self%nf= ",self%nf
  ! Deep copy
  self%fields_obj = rhs%fields_obj
  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init
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


!subroutine add_incr(geom,self,rhs)

!implicit none
! type(wrf_hydro_nwm_jedi_geometry),      intent(inout) :: geom
! type(wrf_hydro_nwm_jedi_state),     intent(inout) :: self
! type(fv3jedi_increment), intent(in)    :: rhs

! integer :: var,i,j,k
! logical :: found_neg
! type(fv3jedi_field), pointer :: field_pointer

! real(kind=kind_real), allocatable :: rhs_ud(:,:,:), rhs_vd(:,:,:)
! real(kind=kind_real), allocatable :: rhs_delp(:,:,:)

! real(kind=kind_real), pointer :: self_ua(:,:,:)
! real(kind=kind_real), pointer :: self_va(:,:,:)
! real(kind=kind_real), pointer :: self_ud(:,:,:)
! real(kind=kind_real), pointer :: self_vd(:,:,:)
! real(kind=kind_real), pointer :: self_t (:,:,:)
! real(kind=kind_real), pointer :: self_pt(:,:,:)
! real(kind=kind_real), pointer :: self_delp(:,:,:)
! real(kind=kind_real), pointer :: self_ps(:,:,:)

! real(kind=kind_real), pointer :: rhs_ua(:,:,:)
! real(kind=kind_real), pointer :: rhs_va(:,:,:)
! real(kind=kind_real), pointer :: rhs_t (:,:,:)
! real(kind=kind_real), pointer :: rhs_pt(:,:,:)
! real(kind=kind_real), pointer :: rhs_ps(:,:,:)

! ! Handle special cases, e.g. D grid winds to A grid winds
! ! -------------------------------------------------------

! ! Get D-Grid winds if necessary
! if (rhs%have_agrid) then !A-Grid in increment
!   if (.not.rhs%have_dgrid) then !D-Grid not in increment
!     if (self%have_dgrid) then !D-grid in state
!       allocate(rhs_ud(rhs%isc:rhs%iec  ,rhs%jsc:rhs%jec+1,1:rhs%npz))
!       allocate(rhs_vd(rhs%isc:rhs%iec+1,rhs%jsc:rhs%jec  ,1:rhs%npz))
!       call pointer_field_array(rhs%fields, 'ua', rhs_ua)
!       call pointer_field_array(rhs%fields, 'va', rhs_va)
!       call a2d(geom, rhs_ua, rhs_va, rhs_ud, rhs_vd) !Linear
!     endif
!   endif
! endif


! ! Convert ps to delp if necessary
! ! -------------------------------
! if (has_field(rhs%fields, 'ps')) then !ps in increment
!   if (.not.has_field(rhs%fields, 'delp')) then !delp not in increment
!     if (has_field(self%fields, 'delp')) then !delp in state
!       allocate(rhs_delp(rhs%isc:rhs%iec,rhs%jsc:rhs%jec,1:rhs%npz))
!       call pointer_field_array(rhs%fields, 'ps', rhs_ps)
!       do k = 1,rhs%npz
!         rhs_delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*rhs_ps(:,:,1) !TLM
!       enddo
!     endif
!   endif
! endif

! !Fields to add determined from increment
! do var = 1,rhs%nf

!   !Winds are a special case
!   if (rhs%fields(var)%fv3jedi_name == 'ua') then

!     if (has_field(self%fields, 'ua')) then
!       call pointer_field_array(rhs%fields,  'ua', rhs_ua)
!       call pointer_field_array(self%fields, 'ua', self_ua)
!       self_ua = self_ua + rhs_ua
!     endif
!     if (has_field(self%fields, 'ud') .and. .not.has_field(rhs%fields, 'ud')) then
!       call pointer_field_array(self%fields, 'ud', self_ud)
!       self_ud = self_ud + rhs_ud
!     endif

!   elseif (rhs%fields(var)%fv3jedi_name == 'va') then

!     if (has_field(self%fields, 'va')) then
!       call pointer_field_array(rhs%fields,  'va', rhs_va)
!       call pointer_field_array(self%fields, 'va', self_va)
!       self_va = self_va + rhs_va
!     endif
!     if (has_field(self%fields, 'vd') .and. .not.has_field(rhs%fields, 'vd')) then
!       call pointer_field_array(self%fields, 'vd', self_vd)
!       self_vd = self_vd + rhs_vd
!     endif

!   elseif (rhs%fields(var)%fv3jedi_name == 't') then

!     if (has_field(self%fields, 't')) then
!       call pointer_field_array(rhs%fields,  't', rhs_t)
!       call pointer_field_array(self%fields, 't', self_t)
!       self_t = self_t + rhs_t
!     endif

!     if (has_field(self%fields, 'pt')) then
!       call pointer_field_array(rhs%fields,  'pt', rhs_pt)
!       call pointer_field_array(self%fields, 'pt', self_pt)
!       self_pt = self_pt + rhs_pt
!     endif

!   elseif (rhs%fields(var)%fv3jedi_name == 'ps') then

!     if (has_field(self%fields, 'ps')) then
!       call pointer_field_array(rhs%fields,  'ps', rhs_ps)
!       call pointer_field_array(self%fields, 'ps', self_ps)
!       self_ps = self_ps + rhs_ps
!     endif

!     if (has_field(self%fields, 'delp') .and. .not. has_field(rhs%fields, 'delp')) then
!       call pointer_field_array(self%fields, 'delp', self_delp)
!       self_delp = self_delp + rhs_delp
!     endif

!   else

!     !Get pointer to state
!     call pointer_field(self%fields, rhs%fields(var)%fv3jedi_name, field_pointer)

!     !Add increment to state
!     field_pointer%array = field_pointer%array + rhs%fields(var)%array

!     !Nullify pointer
!     nullify(field_pointer)

!   endif

! enddo

! if (allocated(rhs_ud)) deallocate(rhs_ud)
! if (allocated(rhs_vd)) deallocate(rhs_vd)
! if (allocated(rhs_delp)) deallocate(rhs_delp)

! !Check for negative tracers and increase to 0.0
! do var = 1,self%nf

!   if (self%fields(var)%tracer) then

!     found_neg = .false.

!     do k = 1,self%fields(var)%npz
!       do j = geom%jsc,geom%jec
!         do i = geom%isc,geom%iec
!           if (self%fields(var)%array(i,j,k) < 0.0_kind_real) then
!             found_neg = .true.
!             self%fields(var)%array(i,j,k) = 0.0_kind_real
!           endif
!         enddo
!       enddo
!     enddo

!     !Print message warning about negative tracer removal
!     if (found_neg .and. self%f_comm%rank() == 0) print*, &
!              'wrf_hydro_nwm_jedi_state_mod.add_incr: Removed negative values for '&
!              //trim(self%fields(var)%fv3jedi_name)

!   endif

! enddo

!end subroutine add_incr


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


!subroutine read_state_from_file(self, c_conf)
subroutine read_state_from_file(self, geom, c_conf, vdate)
  use string_utils
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self   !< State
  type(wrf_hydro_nwm_jedi_geometry), intent(inout) :: geom   !< Geometry
  type(c_ptr),                       intent(in)    :: c_conf !< Configuration
  type(datetime),                    intent(inout) :: vdate  !< DateTime
  
  character(len=10) :: filetype
  character(len=255) :: filename_lsm, filename_hydro
  integer :: flipvert
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str
  integer :: ixfull, jxfull, var

  write(*,*) "foobar"
  
  f_conf = fckit_configuration(c_conf)

  call f_conf%get_or_die("filename_lsm", str)
  filename_lsm = str
  deallocate(str)

  call f_conf%get_or_die("filename_hydro", str)
  filename_hydro = str
  deallocate(str)
  
  call self%fields_obj%read_fields_from_file(filename_lsm, filename_hydro)
end subroutine read_state_from_file


! subroutine write_file(geom, self, c_conf, vdate)
!   use fv3jedi_io_latlon_mod
!   implicit none
!   type(wrf_hydro_nwm_jedi_geometry),  intent(inout) :: geom     !< Geometry
!   type(wrf_hydro_nwm_jedi_state), intent(inout) :: self     !< State
!   type(c_ptr),         intent(in)    :: c_conf   !< Configuration
!   type(datetime),      intent(inout) :: vdate    !< DateTime

!   ! type(fv3jedi_io_gfs)  :: gfs
!   ! type(fv3jedi_io_geos) :: geos

!   character(len=10) :: filetype
!   integer :: flipvert
!   type(fckit_configuration) :: f_conf
!   character(len=:), allocatable :: str


!   ! Fortran configuration
!   ! ---------------------
!   f_conf = fckit_configuration(c_conf)


!   call f_conf%get_or_die("filetype",str)
!   filetype = str
!   deallocate(str)

!   if (trim(filetype) == 'gfs') then

!     flipvert = 0
!     if (f_conf%has("flip_vertically")) then
!       call f_conf%get_or_die("flip_vertically",flipvert)
!     endif
!     if (flipvert==1) call flip_array_vertical(self%nf, self%fields)

!     call gfs%setup(f_conf)
!     call gfs%write_all(geom, self%fields, vdate, self%calendar_type, self%date_init)

!     if (flipvert==1) call flip_array_vertical(self%nf, self%fields)

!   elseif (trim(filetype) == 'geos') then

!     call geos%setup(geom, self%fields, vdate, 'write', f_conf)
!     call geos%write_all(geom, self%fields, vdate)
!     call geos%delete()

!   else

!      call abor1_ftn("wrf_hydro_nwm_jedi_state_mod.write: restart type not supported")

!   endif

! end subroutine write_file


subroutine state_print(self, string)
  use iso_c_binding, only : c_null_char
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  
  character(len=1, kind=c_char) :: string(8192)
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
