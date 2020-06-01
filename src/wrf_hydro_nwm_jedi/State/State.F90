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

use wrf_hydro_nwm_jedi_field_mod, only: wrf_hydro_nwm_jedi_fields
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
     !add_incr, &
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
  
  call self%fields_obj%create(geom,vars)
  
end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
! integer :: var

! do var = 1, self%nf
!   call self%fields(var)%deallocate_field()
! enddo
! deallocate(self%fields)

end subroutine delete

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

! call checksame(self%fields,rhs%fields,"wrf_hydro_nwm_jedi_state_mod.copy")

! write(*,*) "Copying field self%nf= ",self%nf

! do var = 1, self%nf
!   self%fields(var) = rhs%fields(var)
! enddo

! self%calendar_type = rhs%calendar_type
! self%date_init = rhs%date_init

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

!SNLIQ(Time, south_north, snow_layers, west_east)

ixfull = geom%xend-geom%xstart+1
jxfull = geom%yend-geom%ystart+1

!call get_from_restart_3d(filename, geom%xstart, geom%xend, geom%xstart, ixfull, jxfull, "SNLIQ"   , self%fields(1)%array )

!write(*,*) "Print first column from State ",self%fields(1)%array(:,1,1)

!  do var = 1, self%nf
!        call get_from_restart_2d_float(filename, geom%xstart, geom%xend, geom%xstart, ixfull, jxfull, self%fields(var)%wrf_hydro_nwm_name, self%fields(var)%array)

!       ! write(*,*) "Print first column from State ",self%fields(var)%array(:,1,1)
! enddo

end subroutine read_file

! ------------------------------------------------------------------------------

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

! ------------------------------------------------------------------------------

subroutine state_print(self)

 implicit none
 type(wrf_hydro_nwm_jedi_state), intent(in) :: self

 ! call fields_print(self%nf, self%fields, "SNEQV", self%f_comm)
 call self%fields_obj%print_all_fields()

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
