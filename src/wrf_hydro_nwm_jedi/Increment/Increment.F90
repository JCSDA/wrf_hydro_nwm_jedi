! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module wrf_hydro_nwm_jedi_increment_mod

  use fckit_configuration_module, only: fckit_configuration
  use datetime_mod
  use fckit_mpi_module
  use oops_variables_mod

  use wrf_hydro_nwm_jedi_field_mod,    only: wrf_hydro_nwm_jedi_fields,checksame
  use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
  use wrf_hydro_nwm_jedi_util_mod,     only: error_handler
  use iso_c_binding, only : c_float
  use wrf_hydro_nwm_jedi_state_mod, only: wrf_hydro_nwm_jedi_state

  !GetValues
  use ufo_locs_mod,          only: ufo_locs
  use ufo_geovals_mod,       only: ufo_geovals

  implicit none

  private
  public :: create, diff_incr!, create_from_other, delete, zeros, random, copy, &
          ! self_add, self_schur, self_sub, self_mul, axpy_inc, axpy_state, &
          ! dot_prod, &
          ! read_file, write_file, &
          ! gpnorm, rms, &
          ! change_resol, &
          ! getpoint, setpoint, &
          ! ug_coord, increment_to_ug, increment_from_ug, dirac, jnormgrad, &
          ! increment_print

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

  ! Instantiate a shallow water increment (state) object from the geometry
    ! ----------------------------------------------------------------------
    
!  self = shallow_water_state_type(geom)

  ! Initialize all arrays to zero
!  call zeros(self)

  end subroutine create

! ------------------------------------------------------------------------------

! subroutine create_from_other(self, other)

!   ! Passed variables
!   type(shallow_water_state_type), intent(inout) :: self  !< Fields
!   type(shallow_water_state_type), intent(   in) :: other !< Other fields

!   ! Create new state from other state
!   self = shallow_water_state_type(other%get_geometry())

!   ! Initialize all arrays to zero
!   call zeros(self)

! end subroutine create_from_other

! ! ------------------------------------------------------------------------------

! subroutine delete(self)

!   type(shallow_water_state_type), intent(inout) :: self

! end subroutine delete

! ! ------------------------------------------------------------------------------

! subroutine zeros(self)

!   type(shallow_water_state_type), intent(inout) :: self

!   type(shallow_water_geometry_type) :: geom
!   real(r8kind), pointer             :: u(:,:), v(:,:), h(:,:)
!   integer                           :: i, j

!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   geom = self%get_geometry()

!   do j=geom%get_yps(), geom%get_ype()
!     do i=geom%get_xps(), geom%get_xpe()
!       u(i,j) = 0.0_r8kind
!       v(i,j) = 0.0_r8kind
!       h(i,j) = 0.0_r8kind
!     end do
!   end do

! end subroutine zeros

! ! ------------------------------------------------------------------------------

! subroutine random(self)

!   type(shallow_water_state_type), intent(inout) :: self

!   integer, parameter                :: rseed = 7
!   type(shallow_water_geometry_type) :: geom
!   real(r8kind), pointer             :: u(:,:), v(:,:), h(:,:)
!   integer                           :: xps, xpe, yps, ype

!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   geom = self%get_geometry()

!   xps = geom%get_xps()
!   xpe = geom%get_xpe()
!   yps = geom%get_yps()
!   ype = geom%get_ype()

!   call normal_distribution(u(xps:xpe, yps:ype), 0.0_r8kind, 1.0_r8kind, rseed)
!   call normal_distribution(v(xps:xpe, yps:ype), 0.0_r8kind, 1.0_r8kind, rseed)
!   call normal_distribution(h(xps:xpe, yps:ype), 0.0_r8kind, 1.0_r8kind, rseed)

! end subroutine random

! ! ------------------------------------------------------------------------------

! subroutine copy(self, rhs)

!   implicit none

!   type(shallow_water_state_type), intent(inout) :: self
!   type(shallow_water_state_type), intent(   in) :: rhs

!   self = rhs

! end subroutine copy

! ! ------------------------------------------------------------------------------

! subroutine self_add(self, rhs)

!   type(shallow_water_state_type), intent(inout) :: self
!   type(shallow_water_state_type), intent(   in) :: rhs

!   type(shallow_water_geometry_type) :: geom_self, geom_rhs
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: rhs_u(:,:), rhs_v(:,:), rhs_h(:,:)
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_rhs = self%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_rhs%get_nx()     .and. &
!            geom_self%get_ny() == geom_rhs%get_ny()     .and. &
!            geom_self%get_xmax() == geom_rhs%get_xmax() .and. &
!            geom_self%get_ymax() == geom_rhs%get_ymax())

!   if (check) then
!     call self%get_u_ptr(self_u)
!     call self%get_v_ptr(self_v)
!     call self%get_h_ptr(self_h)
!     call rhs%get_u_ptr(rhs_u)
!     call rhs%get_v_ptr(rhs_v)
!     call rhs%get_h_ptr(rhs_h)

!     do j=geom_self%get_yps(), geom_self%get_ype()
!        do i=geom_self%get_xps(), geom_self%get_xpe()
!           self_u(i,j) = self_u(i,j) + rhs_u(i,j)
!           self_v(i,j) = self_v(i,j) + rhs_v(i,j)
!           self_h(i,j) = self_h(i,j) + rhs_h(i,j)
!        end do
!     end do
!   else
!      call abor1_ftn("sw increment: self_add not implemented for mismatched resolutions")
!   endif

! end subroutine self_add

! ! ------------------------------------------------------------------------------

! subroutine self_schur(self, rhs)

!   type(shallow_water_state_type), intent(inout) :: self
!   type(shallow_water_state_type), intent(   in) :: rhs

!   type(shallow_water_geometry_type) :: geom_self, geom_rhs
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: rhs_u(:,:), rhs_v(:,:), rhs_h(:,:)
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_rhs = self%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_rhs%get_nx()     .and. &
!            geom_self%get_ny() == geom_rhs%get_ny()     .and. &
!            geom_self%get_xmax() == geom_rhs%get_xmax() .and. &
!            geom_self%get_ymax() == geom_rhs%get_ymax())

!   if (check) then
!     call self%get_u_ptr(self_u)
!     call self%get_v_ptr(self_v)
!     call self%get_h_ptr(self_h)
!     call rhs%get_u_ptr(rhs_u)
!     call rhs%get_v_ptr(rhs_v)
!     call rhs%get_h_ptr(rhs_h)

!     do j=geom_self%get_yps(), geom_self%get_ype()
!        do i=geom_self%get_xps(), geom_self%get_xpe()
!           self_u(i,j) = self_u(i,j) * rhs_u(i,j)
!           self_v(i,j) = self_v(i,j) * rhs_v(i,j)
!           self_h(i,j) = self_h(i,j) * rhs_h(i,j)
!        end do
!     end do
!   else
!      call abor1_ftn("sw increment:  self_schur not implemented for mismatched resolutions")
!   endif

! end subroutine self_schur

! ! ------------------------------------------------------------------------------

! subroutine self_sub(self, rhs)

!   type(shallow_water_state_type), intent(inout) :: self
!   type(shallow_water_state_type), intent(   in) :: rhs

!   type(shallow_water_geometry_type) :: geom_self, geom_rhs
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: rhs_u(:,:), rhs_v(:,:), rhs_h(:,:)
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_rhs = self%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_rhs%get_nx()     .and. &
!            geom_self%get_ny() == geom_rhs%get_ny()     .and. &
!            geom_self%get_xmax() == geom_rhs%get_xmax() .and. &
!            geom_self%get_ymax() == geom_rhs%get_ymax())

!   if (check) then
!     call self%get_u_ptr(self_u)
!     call self%get_v_ptr(self_v)
!     call self%get_h_ptr(self_h)
!     call rhs%get_u_ptr(rhs_u)
!     call rhs%get_v_ptr(rhs_v)
!     call rhs%get_h_ptr(rhs_h)

!     do j=geom_self%get_yps(), geom_self%get_ype()
!        do i=geom_self%get_xps(), geom_self%get_xpe()
!           self_u(i,j) = self_u(i,j) - rhs_u(i,j)
!           self_v(i,j) = self_v(i,j) - rhs_v(i,j)
!           self_h(i,j) = self_h(i,j) - rhs_h(i,j)
!        end do
!     end do
!   else
!      call abor1_ftn("sw increment:  self_sub not implemented for mismatched resolutions")
!   endif

! end subroutine self_sub

! ! ------------------------------------------------------------------------------

! subroutine self_mul(self, zz)

!   implicit none

!   type(shallow_water_state_type), intent(inout) :: self
!   real(kind=r8kind),              intent(   in) :: zz


!   type(shallow_water_geometry_type) :: geom
!   real(r8kind), pointer             :: u(:,:), v(:,:), h(:,:)
!   integer                           :: i, j

!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   geom = self%get_geometry()

!   do j=geom%get_yps(), geom%get_ype()
!     do i=geom%get_xps(), geom%get_xpe()
!       u(i,j) = zz * u(i,j)
!       v(i,j) = zz * v(i,j)
!       h(i,j) = zz * h(i,j)
!     end do
!   end do

! end subroutine self_mul

! ! ------------------------------------------------------------------------------

! subroutine axpy_inc(self, zz, rhs)

!   type(shallow_water_state_type), intent(inout) :: self
!   real(kind=r8kind),              intent(   in) :: zz
!   type(shallow_water_state_type), intent(   in) :: rhs

!   type(shallow_water_geometry_type) :: geom_self, geom_rhs
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: rhs_u(:,:), rhs_v(:,:), rhs_h(:,:)
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_rhs = self%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_rhs%get_nx()     .and. &
!            geom_self%get_ny() == geom_rhs%get_ny()     .and. &
!            geom_self%get_xmax() == geom_rhs%get_xmax() .and. &
!            geom_self%get_ymax() == geom_rhs%get_ymax())

!   if (check) then
!     call self%get_u_ptr(self_u)
!     call self%get_v_ptr(self_v)
!     call self%get_h_ptr(self_h)
!     call rhs%get_u_ptr(rhs_u)
!     call rhs%get_v_ptr(rhs_v)
!     call rhs%get_h_ptr(rhs_h)

!     do j=geom_self%get_yps(), geom_self%get_ype()
!        do i=geom_self%get_xps(), geom_self%get_xpe()
!           self_u(i,j) = self_u(i,j) + zz * rhs_u(i,j)
!           self_v(i,j) = self_v(i,j) + zz * rhs_v(i,j)
!           self_h(i,j) = self_h(i,j) + zz * rhs_h(i,j)
!        end do
!     end do
!   else
!      call abor1_ftn("sw increment:  axpy_inc not implemented for mismatched resolutions")
!   endif

! end subroutine axpy_inc

! ! ------------------------------------------------------------------------------

! subroutine axpy_state(self, zz, rhs)

!   type(shallow_water_state_type), intent(inout) :: self
!   real(kind=r8kind),              intent(   in) :: zz
!   type(shallow_water_state_type), intent(   in) :: rhs

!   type(shallow_water_geometry_type) :: geom_self, geom_rhs
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: rhs_u(:,:), rhs_v(:,:), rhs_h(:,:)
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_rhs = self%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_rhs%get_nx()     .and. &
!            geom_self%get_ny() == geom_rhs%get_ny()     .and. &
!            geom_self%get_xmax() == geom_rhs%get_xmax() .and. &
!            geom_self%get_ymax() == geom_rhs%get_ymax())

!   if (check) then
!     call self%get_u_ptr(self_u)
!     call self%get_v_ptr(self_v)
!     call self%get_h_ptr(self_h)
!     call rhs%get_u_ptr(rhs_u)
!     call rhs%get_v_ptr(rhs_v)
!     call rhs%get_h_ptr(rhs_h)

!     do j=geom_self%get_yps(), geom_self%get_ype()
!        do i=geom_self%get_xps(), geom_self%get_xpe()
!           self_u(i,j) = self_u(i,j) + zz * rhs_u(i,j)
!           self_v(i,j) = self_v(i,j) + zz * rhs_v(i,j)
!           self_h(i,j) = self_h(i,j) + zz * rhs_h(i,j)
!        end do
!     end do
!   else
!      call abor1_ftn("sw increment:  axpy_state not implemented for mismatched resolutions")
!   endif

! end subroutine axpy_state

! ! ------------------------------------------------------------------------------

! subroutine dot_prod(self, other, zprod)

!   type(shallow_water_state_type), intent(   in) :: self
!   type(shallow_water_state_type), intent(   in) :: other
!   real(kind=r8kind),              intent(inout) :: zprod

!   type(fckit_mpi_comm)              :: f_comm
!   type(shallow_water_geometry_type) :: geom_self, geom_other
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: other_u(:,:), other_v(:,:), other_h(:,:)
!   real(kind=r8kind)                 :: zp
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_other = other%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_other%get_nx()     .and. &
!            geom_self%get_ny() == geom_other%get_ny()     .and. &
!            geom_self%get_xmax() == geom_other%get_xmax() .and. &
!            geom_self%get_ymax() == geom_other%get_ymax())

!   if (check) then
!     f_comm = fckit_mpi_comm()

!     call self%get_u_ptr(self_u)
!     call self%get_v_ptr(self_v)
!     call self%get_h_ptr(self_h)
!     call other%get_u_ptr(other_u)
!     call other%get_v_ptr(other_v)
!     call other%get_h_ptr(other_h)

!     zp = 0.0_r8kind
!     do j=geom_self%get_yps(), geom_self%get_ype()
!        do i=geom_self%get_xps(), geom_self%get_xpe()
!           zp = zp + self_u(i,j) * other_u(i,j)
!           zp = zp + self_v(i,j) * other_v(i,j)
!           zp = zp + self_h(i,j) * other_h(i,j)
!        end do
!     end do

!     !Get global dot product
!     call f_comm%allreduce(zp, zprod, fckit_mpi_sum())

!     !For debugging print result:
!     if (f_comm%rank() == 0) print*, "Dot product test result: ", zprod

!   else
!      call abor1_ftn("sw increment:  dot_prod not implemented for mismatched resolutions")
!   endif

! end subroutine dot_prod

! ! ------------------------------------------------------------------------------

subroutine diff_incr(self, x1, x2)

  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: x1
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: x2

  type(wrf_hydro_nwm_jedi_geometry) :: geom_self, geom_x1, geom_x2
  real, pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
  real, pointer             :: x1_u(:,:), x1_v(:,:), x1_h(:,:)
  real, pointer             :: x2_u(:,:), x2_v(:,:), x2_h(:,:)
  integer                           :: i, j
  logical                           :: check

  ! ! Get geometries
  ! geom_self = self%get_geometry()
  ! geom_x1 = x1%get_geometry()
  ! geom_x2 = x2%get_geometry()

  ! ! Check for matching resolution
  ! check = (geom_self%get_nx() == geom_x1%get_nx()     .and. &
  !          geom_self%get_ny() == geom_x1%get_ny()     .and. &
  !          geom_self%get_xmax() == geom_x1%get_xmax() .and. &
  !          geom_self%get_ymax() == geom_x1%get_ymax() .and. &
  !          geom_self%get_nx() == geom_x2%get_nx()     .and. &
  !          geom_self%get_ny() == geom_x2%get_ny()     .and. &
  !          geom_self%get_xmax() == geom_x2%get_xmax() .and. &
  !          geom_self%get_ymax() == geom_x2%get_ymax())

  ! if (check) then
  !   call self%get_u_ptr(self_u)
  !   call self%get_v_ptr(self_v)
  !   call self%get_h_ptr(self_h)
  !   call x1%get_u_ptr(x1_u)
  !   call x1%get_v_ptr(x1_v)
  !   call x1%get_h_ptr(x1_h)
  !   call x2%get_u_ptr(x2_u)
  !   call x2%get_v_ptr(x2_v)
  !   call x2%get_h_ptr(x2_h)

  !   do j=geom_self%get_yps(), geom_self%get_ype()
  !      do i=geom_self%get_xps(), geom_self%get_xpe()
  !         self_u(i,j) = x1_u(i,j) - x2_u(i,j)
  !         self_v(i,j) = x1_v(i,j) - x2_v(i,j)
  !         self_h(i,j) = x1_h(i,j) - x2_h(i,j)
  !      end do
  !   end do
  ! else
  !    call abor1_ftn("sw increment:  self_sub not implemented for mismatched resolutions")
  ! endif

end subroutine diff_incr

! ! ------------------------------------------------------------------------------

! subroutine change_resol(self, rhs)

!   type(shallow_water_state_type), intent(inout) :: self
!   type(shallow_water_state_type), intent(   in) :: rhs

!   type(shallow_water_geometry_type) :: geom_self, geom_rhs
!   real(r8kind), pointer             :: self_u(:,:), self_v(:,:), self_h(:,:)
!   real(r8kind), pointer             :: rhs_u(:,:), rhs_v(:,:), rhs_h(:,:)
!   integer                           :: i, j
!   logical                           :: check

!   ! Get geometries
!   geom_self = self%get_geometry()
!   geom_rhs = self%get_geometry()

!   ! Check for matching resolution
!   check = (geom_self%get_nx() == geom_rhs%get_nx()     .and. &
!            geom_self%get_ny() == geom_rhs%get_ny()     .and. &
!            geom_self%get_xmax() == geom_rhs%get_xmax() .and. &
!            geom_self%get_ymax() == geom_rhs%get_ymax())

!   if (check) then
!     call copy(self, rhs)
!   else
!      call abor1_ftn("sw increment:  change_resol not implemented yet")
!   endif

! end subroutine change_resol

! ! ------------------------------------------------------------------------------

! subroutine read_file(geom, self, c_conf, vdate)

!   implicit none

!   type(shallow_water_geometry_type), intent(inout) :: geom     !< Geometry
!   type(shallow_water_state_type),    intent(inout) :: self     !< Increment
!   type(c_ptr),                       intent(   in) :: c_conf   !< Configuration
!   type(datetime),                    intent(inout) :: vdate    !< DateTime

!   type(fckit_configuration)     :: f_conf
!   character(len=:), allocatable :: str
!   character(len=200)            :: sdate
!   character(len=255)            :: filename

!   f_conf = fckit_configuration(c_conf)

!   ! Get the input file name
!   call f_conf%get_or_die("filename", str)
!   filename = trim(str)
!   deallocate(str)

!   ! Read the state from the file
!   call self%read(trim(filename))

!   ! Set the valid time to the initialization date
!   call f_conf%get_or_die("date", str)
!   sdate = trim(str)
!   deallocate(str)
!   call datetime_set(sdate, vdate)

! end subroutine read_file

! ! ------------------------------------------------------------------------------

! subroutine write_file(geom, self, c_conf, vdate)

!   implicit none

!   type(shallow_water_geometry_type), intent(inout) :: geom     !< Geometry
!   type(shallow_water_state_type),    intent(   in) :: self     !< Increment
!   type(c_ptr),                       intent(   in) :: c_conf   !< Configuration
!   type(datetime),                    intent(inout) :: vdate    !< DateTime

!   integer, parameter               :: max_string_length=800 ! Yuk!
!   character(len=max_string_length) :: filename
!   type(fckit_mpi_comm)             :: f_comm

!   ! Get MPI communicator
!   f_comm = fckit_mpi_comm()

!   ! Generate an output file name
!   filename = genfilename(c_conf, max_string_length, vdate)

!   ! Write the increment to the file
!   if (f_comm%rank() == 0) then
!     write(*,*) 'sw_increment_mod:write_file: writing ' // trim(filename)
!   end if
!   call self%write(trim(filename))

! end subroutine write_file

! ! ------------------------------------------------------------------------------

! subroutine increment_print(self)

!   implicit none
!   type(shallow_water_state_type), intent(in) :: self

!   type(fckit_mpi_comm)              :: f_comm
!   type(shallow_water_geometry_type) :: geom
!   integer                           :: xps, xpe, yps, ype, i
!   real(r8kind)                      :: n
!   real(r8kind)                      :: tmp(3), pstat(3)
!   real(r8kind), pointer             :: field(:,:)
!   character(len=1)                  :: fields(3)

!   ! Get MPI communicator
!   f_comm = fckit_mpi_comm()

!   ! Set fields
!   fields(1) = "U"
!   fields(2) = "V"
!   fields(3) = "H"

!   ! Get local indices
!   geom = self%get_geometry()
!   xps = geom%get_xps()
!   xpe = geom%get_xpe()
!   yps = geom%get_yps()
!   ype = geom%get_ype()

!   ! Get total global size of fields
!   n = real(geom%get_nx() * geom%get_ny(), r8kind)
  
!   ! Print header
!   if (f_comm%rank() == 0) then
!     write(*,"(A70)") "----------------------------------------------------------------------"
!     write(*,"(A17)") "|     Increment print" 
!     write(*,"(A70)") "----------------------------------------------------------------------"
!   endif

!   ! Print field stats
!   do i = 1, size(fields)
!     select case(fields(i))
!       case("U")
!         call self%get_u_ptr(field)
!       case("V")
!         call self%get_v_ptr(field)
!       case("H")
!         call self%get_h_ptr(field)
!       case DEFAULT
!         call abor1_ftn("sw_increment: increment_print, invalid field" // fields(i))         
!     end select

!     tmp(1) = minval(field(xps:xpe, yps:ype))
!     tmp(2) = maxval(field(xps:xpe, yps:ype))
!     tmp(3) =    sum(field(xps:xpe, yps:ype)**2)
!     call f_comm%allreduce(tmp(1), pstat(1), fckit_mpi_min())
!     call f_comm%allreduce(tmp(2), pstat(2), fckit_mpi_max())
!     call f_comm%allreduce(tmp(3), pstat(3), fckit_mpi_sum())
!     pstat(3) = sqrt(pstat(3) / n)
!     if (f_comm%rank() == 0) write(*,"(A10,A6,ES14.7,A6,ES14.7,A6,ES14.7)") &
!                                      fields(i),                            &
!                                      "| Min=",real(pstat(1),4),            &
!                                      ", Max=",real(pstat(2),4),            &
!                                      ", RMS=",real(pstat(3),4)
!   end do

!   if (f_comm%rank() == 0) &
!     write(*,"(A70)") "---------------------------------------------------------------------"

! end subroutine increment_print

! ! ------------------------------------------------------------------------------

! subroutine gpnorm(self, nf, pstat)

!   type(shallow_water_state_type), intent(   in) :: self
!   integer,                        intent(   in) :: nf
!   real(kind=r8kind),              intent(inout) :: pstat(3, nf)

!   type(fckit_mpi_comm)              :: f_comm
!   type(shallow_water_geometry_type) :: geom
!   integer                           :: xps, xpe, yps, ype, i
!   real(r8kind)                      :: n
!   real(r8kind)                      :: tmp(3)
!   real(r8kind), pointer             :: field(:,:)
!   character(len=1)                  :: fields(3)

!   ! Get MPI communicator
!   f_comm = fckit_mpi_comm()

!   ! Set fields
!   fields(1) = "U"
!   fields(2) = "V"
!   fields(3) = "H"

!   ! Get local indices
!   geom = self%get_geometry()
!   xps = geom%get_xps()
!   xpe = geom%get_xpe()
!   yps = geom%get_yps()
!   ype = geom%get_ype()

!   ! Get total global size of fields
!   n = real(geom%get_nx() * geom%get_ny(), r8kind)
  
!   ! Calculate field stats
!   do i = 1, size(fields)
!     select case(fields(i))
!       case("U")
!         call self%get_u_ptr(field)
!       case("V")
!         call self%get_v_ptr(field)
!       case("H")
!         call self%get_h_ptr(field)
!       case DEFAULT
!         call abor1_ftn("sw_state: , gpnorm invalid field" // fields(i))
!     end select

!     tmp(1) = minval(field(xps:xpe, yps:ype))
!     tmp(2) = maxval(field(xps:xpe, yps:ype))
!     tmp(3) =    sum(field(xps:xpe, yps:ype)**2)
!     call f_comm%allreduce(tmp(1), pstat(1, i), fckit_mpi_min())
!     call f_comm%allreduce(tmp(2), pstat(2, i), fckit_mpi_max())
!     call f_comm%allreduce(tmp(3), pstat(3, i), fckit_mpi_sum())
!     pstat(3, i) = sqrt(pstat(3, i) / n)

!   end do

! end subroutine gpnorm

! ! ------------------------------------------------------------------------------

! subroutine rms(self, prms)

!   type(shallow_water_state_type), intent( in) :: self
!   real(kind=r8kind),              intent(out) :: prms

!   type(fckit_mpi_comm)              :: f_comm
!   type(shallow_water_geometry_type) :: geom
!   integer                           :: xps, xpe, yps, ype
!   integer                           :: i, j, n
!   real(r8kind), pointer             :: u(:,:), v(:,:), h(:,:)
!   real(r8kind)                      :: sum_sqr

!   ! Get MPI communicator
!   f_comm = fckit_mpi_comm()

!   ! Get local indices
!   geom = self%get_geometry()
!   xps = geom%get_xps()
!   xpe = geom%get_xpe()
!   yps = geom%get_yps()
!   ype = geom%get_ype()

!   ! Get total global number of field elements
!   n = geom%get_nx() * geom%get_ny() * 3

!   ! Get pointers to fields
!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   ! Calculate local sum of squares
!   sum_sqr = 0.0_r8kind
!   do j = yps, ype
!     do i = xps, xpe
!       sum_sqr = sum_sqr + u(i,j)**2
!       sum_sqr = sum_sqr + v(i,j)**2
!       sum_sqr = sum_sqr + h(i,j)**2
!     enddo
!   enddo

!   ! Get global sum of squares
!   call f_comm%allreduce(sum_sqr, prms, fckit_mpi_sum())

!   ! Calculate rms
!   prms = sqrt(prms / real(n, r8kind))

! end subroutine rms

! ! ------------------------------------------------------------------------------

! subroutine dirac(self, c_conf, geom)

!   type(shallow_water_state_type),    intent(inout) :: self
!   type(c_ptr),                       intent(   in) :: c_conf
!   type(shallow_water_geometry_type), intent(   in) :: geom

!   call abor1_ftn("sw increment:  dirac not implemented yet")

! end subroutine dirac

! ! ------------------------------------------------------------------------------

! subroutine getpoint(self, geoiter, values)

!   type(shallow_water_state_type), intent(   in) :: self
!   type(sw_geom_iter),             intent(   in) :: geoiter
!   real(r8kind),                   intent(inout) :: values(3)

!   ! Pointers to the fields
!   real(r8kind), pointer :: u(:,:), v(:,:), h(:,:)

!   ! Get pointers to fields
!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   values(1) = u(geoiter%ilon, geoiter%ilat)
!   values(2) = v(geoiter%ilon, geoiter%ilat)
!   values(3) = h(geoiter%ilon, geoiter%ilat)

! end subroutine getpoint

! ! ------------------------------------------------------------------------------

! subroutine setpoint(self, geoiter, values)

!   ! Passed variables
!   type(shallow_water_state_type), intent(inout) :: self
!   type(sw_geom_iter),             intent(   in) :: geoiter
!   real(r8kind),                   intent(   in) :: values(3)

!   ! Pointers to the fields
!   real(r8kind), pointer :: u(:,:), v(:,:), h(:,:)

!   ! Get pointers to fields
!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   ! Set values
!   u(geoiter%ilon, geoiter%ilat) = values(1)
!   v(geoiter%ilon, geoiter%ilat) = values(2)
!   h(geoiter%ilon, geoiter%ilat) = values(3)

! end subroutine setpoint

! ! ------------------------------------------------------------------------------

! subroutine ug_size(self, ug)

!   type(shallow_water_state_type), intent(   in) :: self
!   type(unstructured_grid),        intent(inout) :: ug

!   type(shallow_water_geometry_type) :: geom

!   geom = self%get_geometry()

!   ! Set number of grids
!   if (ug%colocated==1) then
!      ! Colocatd
!      ug%ngrid = 1
!   else
!      ! Not colocated
!      ug%ngrid = 1
!      call abor1_ftn("fv3jedi_increment_mod.ug_size: Uncolocated grids not coded yet")
!   end if

!   ! Allocate grid instances
!   if (.not.allocated(ug%grid)) allocate(ug%grid(ug%ngrid))

!   ! Set local number of points
!   ug%grid(1)%nmga = geom%get_npx() * geom%get_npy()

!   ! Set number of levels
!   ug%grid(1)%nl0 = 1

!   ! Set number of variables
!   ug%grid(1)%nv = 3

!   ! Set number of timeslots
!   ug%grid(1)%nts = ug%nts

! end subroutine ug_size

! ! ------------------------------------------------------------------------------

! subroutine ug_coord(self, ug, geom)

!   type(shallow_water_state_type),    intent(   in) :: self
!   type(unstructured_grid),           intent(inout) :: ug
!   type(shallow_water_geometry_type), intent(   in) :: geom

!   integer :: n, i, j
!   real(r8kind) :: dx, dy

!   ! Define size
!   call ug_size(self, ug)

!   ! Allocate unstructured grid coordinates
!   call allocate_unstructured_grid_coord(ug)

!   ! Get grid spacing
!   dx = geom%get_dx()
!   dy = geom%get_dy()

!   if (ug%colocated==1) then
!     n = 0
!     do j = geom%get_yps(), geom%get_ype()
!       do i = geom%get_xps(), geom%get_xpe()
!         n = n + 1
!         ug%grid(1)%lon(n) = dx * (i - 1)
!         ug%grid(1)%lat(n) = dy * (j - 1)
!         ug%grid(1)%area(n) = 1.0_r8kind
!         ug%grid(1)%vunit(n, 1) = 1.0_r8kind ! 'fake' sigma coordinates
!         ug%grid(1)%lmask(n, 1) = .true.
!       enddo
!     enddo
!   endif

! end subroutine ug_coord

! ! ------------------------------------------------------------------------------

! subroutine increment_to_ug(self, ug, its)

!   type(shallow_water_state_type), intent(   in) :: self
!   type(unstructured_grid),        intent(inout) :: ug
!   integer,                        intent(   in) :: its

!   type(shallow_water_geometry_type) :: geom
!   integer                           :: n, i, j
!   real(r8kind)                      :: dx, dy
!   real(r8kind), pointer             :: u(:,:), v(:,:), h(:,:)

!   ! Get the geometry
!   geom = self%get_geometry()

!   ! Define size
!   call ug_size(self, ug)

!   ! Allocate unstructured grid increment
!   call allocate_unstructured_grid_field(ug)

!   ! Get grid spacing
!   dx = geom%get_dx()
!   dy = geom%get_dy()

!   ! Get pointers to fields
!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   ! Copy increment
!   ug%grid(1)%fld(:,:,:,its) = 0.0_r8kind

!   if (ug%colocated==1) then

!     ! U
!     n = 0
!     do j = geom%get_yps(), geom%get_ype()
!       do i = geom%get_xps(), geom%get_xpe()
!         n = n + 1
!         ug%grid(1)%fld(n, 1, 1, its) = u(i, j)
!       enddo
!     enddo

!     ! V
!     n = 0
!     do j = geom%get_yps(), geom%get_ype()
!       do i = geom%get_xps(), geom%get_xpe()
!         n = n + 1
!         ug%grid(1)%fld(n, 1, 2, its) = v(i, j)
!       enddo
!     enddo

!     ! H
!     n = 0
!     do j = geom%get_yps(), geom%get_ype()
!       do i = geom%get_xps(), geom%get_xpe()
!         n = n + 1
!         ug%grid(1)%fld(n, 1, 3, its) = h(i, j)
!       enddo
!     enddo
!   endif

! end subroutine increment_to_ug

! ! -----------------------------------------------------------------------------

! subroutine increment_from_ug(self, ug, its)

!   type(shallow_water_state_type), intent(inout) :: self
!   type(unstructured_grid),        intent(   in) :: ug
!   integer,                        intent(   in) :: its

!   type(shallow_water_geometry_type) :: geom
!   integer                           :: n, i, j
!   real(r8kind)                      :: dx, dy
!   real(r8kind), pointer             :: u(:,:), v(:,:), h(:,:)

!   ! Get the geometry
!   geom = self%get_geometry()

!   ! Get grid spacing
!   dx = geom%get_dx()
!   dy = geom%get_dy()

!   ! Get pointers to fields
!   call self%get_u_ptr(u)
!   call self%get_v_ptr(v)
!   call self%get_h_ptr(h)

!   ! Copy increment
!   if (ug%colocated==1) then

!       n = 0
!       do j = geom%get_yps(), geom%get_ype()
!         do i = geom%get_xps(), geom%get_xps()
!           n = n + 1
!           u(i,j) = ug%grid(1)%fld(n,1,1,its)
!         enddo
!       enddo

!       n = 0
!       do j = geom%get_yps(), geom%get_ype()
!         do i = geom%get_xps(), geom%get_xps()
!           n = n + 1
!           v(i,j) = ug%grid(1)%fld(n,1,1,its)
!         enddo
!       enddo

!       n = 0
!       do j = geom%get_yps(), geom%get_ype()
!         do i = geom%get_xps(), geom%get_xps()
!           n = n + 1
!           h(i,j) = ug%grid(1)%fld(n,1,1,its)
!         enddo
!       enddo

!   endif

! end subroutine increment_from_ug

! ! ------------------------------------------------------------------------------

! subroutine jnormgrad(self, geom, ref, c_conf)

!   implicit none

!   type(shallow_water_state_type),    intent(inout) :: self
!   type(shallow_water_geometry_type), intent(   in) :: geom
!   type(shallow_water_state_type),    intent(   in) :: ref !To linearize around if nl
!   type(c_ptr),                       intent(   in) :: c_conf

!   call abor1_ftn("sw increment:  jnormgrad not implemented yet")

! end subroutine jnormgrad

! ! ------------------------------------------------------------------------------

! function genfilename (c_conf,length,vdate)

!   use iso_c_binding
!   use datetime_mod
!   use duration_mod

!   type(c_ptr),    intent(in) :: c_conf  !< Configuration
!   integer,        intent(in) :: length
!   type(datetime), intent(in) :: vdate

!   type(fckit_configuration)     :: f_conf
!   character(len=:), allocatable :: str
!   character(len=length)         :: genfilename
!   character(len=length)         :: fdbdir, expver, typ, validitydate, referencedate, sstep
!   character(len=length)         :: prefix, mmb
!   type(datetime)                :: rdate
!   type(duration)                :: step
!   integer                       :: lenfn

!   f_conf = fckit_configuration(c_conf)

!   ! here we should query the length and then allocate "string".
!   ! But Fortran 90 does not allow variable-length allocatable strings.
!   call f_conf%get_or_die("datadir", str)
!   fdbdir = trim(str)
!   deallocate(str)
!   call f_conf%get_or_die("exp", str)
!   expver = trim(str)
!   deallocate(str)
!   call f_conf%get_or_die("type", str)
!   typ = trim(str)
!   deallocate(str)

!   if (typ=="ens") then
!     call f_conf%get_or_die("member", str)
!     mmb = trim(str)
!     deallocate(str)

!     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ) + 1 + LEN_TRIM(mmb)
!     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ) // "." // TRIM(mmb)
!   else
!     lenfn = LEN_TRIM(fdbdir) + 1 + LEN_TRIM(expver) + 1 + LEN_TRIM(typ)
!     prefix = TRIM(fdbdir) // "/" // TRIM(expver) // "." // TRIM(typ)
!   endif

!   if (typ=="fc" .or. typ=="ens") then
!     call f_conf%get_or_die("date", str)
!     referencedate = trim(str)
!     deallocate(str)
!     call datetime_to_string(vdate, validitydate)
!     call datetime_create(TRIM(referencedate), rdate)
!     call datetime_diff(vdate, rdate, step)
!     call duration_to_string(step, sstep)
!     lenfn = lenfn + 1 + LEN_TRIM(referencedate) + 1 + LEN_TRIM(sstep)
!     genfilename = TRIM(prefix) // "." // TRIM(referencedate) // "." // TRIM(sstep) // ".nc"
!   endif

!   if (typ=="an") then
!     call datetime_to_string(vdate, validitydate)
!     lenfn = lenfn + 1 + LEN_TRIM(validitydate)
!     genfilename = TRIM(prefix) // "." // TRIM(validitydate) // ".nc"
!   endif

!   if (lenfn>length) &
!     & call abor1_ftn("sw_increment_mod:genfilename: filename too long")

! end function genfilename

! ------------------------------------------------------------------------------

end module wrf_hydro_nwm_jedi_increment_mod
