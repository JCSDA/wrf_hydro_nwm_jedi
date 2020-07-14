! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module wrf_hydro_nwm_jedi_state_interface_mod

!use fv3jedi_kinds_mod
use datetime_mod
use duration_mod
use iso_c_binding, only: c_int, c_float, c_ptr, c_char
use oops_variables_mod
use fckit_configuration_module, only: fckit_configuration

use wrf_hydro_nwm_jedi_state_mod
use wrf_hydro_nwm_jedi_state_utils_mod, only: wrf_hydro_nwm_jedi_state_registry
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_mod_c, only: wrf_hydro_nwm_jedi_geometry_registry
!use fv3jedi_increment_utils_mod, only: fv3jedi_increment, fv3jedi_increment_registry

!GetValues
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
!use fv3jedi_getvalues_traj_mod, only: fv3jedi_getvalues_traj, fv3jedi_getvalues_traj_registry

private
public :: wrf_hydro_nwm_jedi_state_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

  
subroutine wrf_hydro_nwm_jedi_state_create_c(c_key_self, c_key_geom, c_vars) &
     bind(c,name='wrf_hydro_nwm_jedi_state_create_f90')
  implicit none
  integer(c_int), intent(inout)  :: c_key_self
  integer(c_int), intent(in)     :: c_key_geom !< Geometry
  type(c_ptr), value, intent(in) :: c_vars     !< List of variables
  
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  type(oops_variables)         :: vars
  
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_state_registry%init()
  call wrf_hydro_nwm_jedi_state_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  
  vars = oops_variables(c_vars)
  call create(self, geom, vars)
end subroutine wrf_hydro_nwm_jedi_state_create_c


! subroutine wrf_hydro_nwm_jedi_state_create_from_other_c(c_key_self, c_key_other) bind(c,name='wrf_hydro_nwm_jedi_state_create_from_other_f90')

!   implicit none
!   integer(c_int),intent(inout) :: c_key_self  !< Fields
!   integer(c_int),intent(   in) :: c_key_other !< Other fields

!   ! Local variables
!   type(wrf_hydro_nwm_jedi_state), pointer :: self
!   type(wrf_hydro_nwm_jedi_state), pointer :: other

!   ! Interface
!   call wrf_hydro_nwm_jedi_state_registry%get(c_key_other, other)
!   call wrf_hydro_nwm_jedi_state_registry%init()
!   call wrf_hydro_nwm_jedi_state_registry%add(c_key_self)
!   call wrf_hydro_nwm_jedi_state_registry%get(c_key_self, self)

!   ! Call Fortran
!   call create_from_other(self, other)
! end subroutine wrf_hydro_nwm_jedi_state_create_from_other_c


subroutine wrf_hydro_nwm_jedi_state_delete_c(c_key_self) &
     bind(c,name='wrf_hydro_nwm_jedi_state_delete_f90')

  implicit none
  integer(c_int), intent(inout) :: c_key_self
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  
  call delete(self)
  
  call wrf_hydro_nwm_jedi_state_registry%remove(c_key_self)
end subroutine wrf_hydro_nwm_jedi_state_delete_c


subroutine wrf_hydro_nwm_jedi_state_zero_c(c_key_self) &
     bind(c,name='wrf_hydro_nwm_jedi_state_zero_f90')
implicit none
integer(c_int), intent(in) :: c_key_self
! type(wrf_hydro_nwm_jedi_state), pointer :: self

! call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
! call zeros(self)
end subroutine wrf_hydro_nwm_jedi_state_zero_c


subroutine wrf_hydro_nwm_jedi_state_copy_c(c_key_self,c_key_rhs) &
     bind(c,name='wrf_hydro_nwm_jedi_state_copy_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs
  
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: rhs
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_rhs,rhs)
  
  call copy(self, rhs)
end subroutine wrf_hydro_nwm_jedi_state_copy_c


subroutine wrf_hydro_nwm_jedi_state_axpy_c(c_key_self,c_zz,c_key_rhs) &
     bind(c,name='wrf_hydro_nwm_jedi_state_axpy_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  real(c_float), intent(in) :: c_zz
  integer(c_int), intent(in) :: c_key_rhs
  
  ! type(wrf_hydro_nwm_jedi_state), pointer :: self
  ! type(wrf_hydro_nwm_jedi_state), pointer :: rhs
  ! real(kind=kind_real) :: zz
  
  ! call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  ! call wrf_hydro_nwm_jedi_state_registry%get(c_key_rhs,rhs)
  ! zz = c_zz
  
  ! call axpy(self,zz,rhs)
end subroutine wrf_hydro_nwm_jedi_state_axpy_c


! subroutine wrf_hydro_nwm_jedi_state_add_incr_c(c_key_geom,c_key_self,c_key_rhs) bind(c,name='wrf_hydro_nwm_jedi_state_add_incr_f90')

! implicit none
! integer(c_int), intent(in) :: c_key_geom
! integer(c_int), intent(in) :: c_key_self
! integer(c_int), intent(in) :: c_key_rhs
! type(fv3jedi_geom), pointer :: geom
! type(wrf_hydro_nwm_jedi_state), pointer :: self
! type(fv3jedi_increment), pointer :: rhs

! call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
! call fv3jedi_increment_registry%get(c_key_rhs,rhs)
! call fv3jedi_geom_registry%get(c_key_geom, geom)

! call add_incr(geom,self,rhs)

!end subroutine wrf_hydro_nwm_jedi_state_add_incr_c

! ------------------------------------------------------------------------------

!subroutine wrf_hydro_nwm_jedi_state_change_resol_c(c_key_state,c_key_geom,c_key_rhs,c_key_geom_rhs) bind(c,name='wrf_hydro_nwm_jedi_state_change_resol_f90')

! implicit none
! integer(c_int), intent(in) :: c_key_state
! integer(c_int), intent(in) :: c_key_geom
! integer(c_int), intent(in) :: c_key_rhs
! integer(c_int), intent(in) :: c_key_geom_rhs

! type(wrf_hydro_nwm_jedi_state), pointer :: state, rhs
! type(fv3jedi_geom),  pointer :: geom, geom_rhs

! call wrf_hydro_nwm_jedi_state_registry%get(c_key_state,state)
! call fv3jedi_geom_registry%get(c_key_geom, geom)
! call wrf_hydro_nwm_jedi_state_registry%get(c_key_rhs,rhs)
! call fv3jedi_geom_registry%get(c_key_geom_rhs, geom_rhs)

! call change_resol(state,geom,rhs,geom_rhs)

!end subroutine wrf_hydro_nwm_jedi_state_change_resol_c

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_state_read_file_c(c_key_geom, c_key_state, c_conf, c_dt) &
     bind(c,name='wrf_hydro_nwm_jedi_state_read_file_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_state  !< State
  type(c_ptr), intent(in)    :: c_conf !< Configuration
  type(c_ptr), intent(inout) :: c_dt   !< DateTime
  integer(c_int), intent(in) :: c_key_geom  !< Geometry
  
  type(wrf_hydro_nwm_jedi_state), pointer :: state
  type(datetime) :: fdate
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  
  write(*,*) "Key_geom from read_state_from_file", c_key_geom
  
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_state, state)
  call c_f_datetime(c_dt, fdate)
  call read_state_from_file(state, geom, c_conf, fdate)
  ! call read_state_from_file(state, c_conf)
end subroutine wrf_hydro_nwm_jedi_state_read_file_c


subroutine wrf_hydro_nwm_jedi_state_analytic_init_c(c_key_state, c_key_geom, c_conf, c_dt) &
     bind(c,name='wrf_hydro_nwm_jedi_state_analytic_init_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_state  !< State
  integer(c_int), intent(in) :: c_key_geom  !< Geometry
  type(c_ptr), intent(in)    :: c_conf !< Configuration
  type(c_ptr), intent(inout) :: c_dt   !< DateTime
  
  type(wrf_hydro_nwm_jedi_state), pointer :: state
  type(wrf_hydro_nwm_jedi_geometry), pointer :: geom
  type(datetime) :: fdate
  
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_state,state)
  call c_f_datetime(c_dt, fdate)
  !call analytic_IC(state, geom, c_conf, fdate)
end subroutine wrf_hydro_nwm_jedi_state_analytic_init_c


subroutine wrf_hydro_nwm_jedi_state_write_file_c(c_key_geom, c_key_state, c_conf, c_dt) &
     bind(c,name='wrf_hydro_nwm_jedi_state_write_file_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_state  !< State
  type(c_ptr), intent(in) :: c_conf !< Configuration
  type(c_ptr), intent(in) :: c_dt   !< DateTime
  integer(c_int), intent(in) :: c_key_geom  !< Geometry
  
  type(wrf_hydro_nwm_jedi_state), pointer :: state
  type(datetime) :: fdate
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  
  !call fv3jedi_geom_registry%get(c_key_geom, geom)
  !call wrf_hydro_nwm_jedi_state_registry%get(c_key_state,state)
  !call c_f_datetime(c_dt, fdate)
  !call write_file(geom, state, c_conf, fdate)
end subroutine wrf_hydro_nwm_jedi_state_write_file_c


subroutine wrf_hydro_nwm_jedi_state_gpnorm_c(c_key_state, kf, pstat) &
     bind(c,name='wrf_hydro_nwm_jedi_state_gpnorm_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_state
  integer(c_int), intent(in) :: kf
  real(c_float), intent(inout) :: pstat(3*kf)
  
  ! type(wrf_hydro_nwm_jedi_state), pointer :: state
  ! real(kind=kind_real) :: zstat(3, kf)
  ! integer :: jj, js, jf
  
  ! call wrf_hydro_nwm_jedi_state_registry%get(c_key_state,state)
  
  ! call gpnorm(state, kf, zstat)
  ! jj=0
  ! do jf = 1, kf
  !   do js = 1, 3
  !     jj=jj+1
  !     pstat(jj) = zstat(js,jf)
  !   enddo
  ! enddo
end subroutine wrf_hydro_nwm_jedi_state_gpnorm_c


subroutine wrf_hydro_nwm_jedi_state_print_c(c_key_self,string) &
     bind(c,name='wrf_hydro_nwm_jedi_state_print_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  character(len=1,kind=c_char) :: string(8192)
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  
  call state_print(self,string)
end subroutine wrf_hydro_nwm_jedi_state_print_c


subroutine wrf_hydro_nwm_jedi_state_get_mean_stddev_c(c_key_self,nf,pstat) &
     bind(c,name='wrf_hydro_nwm_jedi_state_get_mean_stddev_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  integer(c_int), intent(in) :: nf
  real(c_float), intent(inout) :: pstat(3,nf)
  
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  
  call get_mean_stddev(self,nf,pstat)
end subroutine wrf_hydro_nwm_jedi_state_get_mean_stddev_c


subroutine wrf_hydro_nwm_jedi_state_rms_c(c_key_state, prms) &
     bind(c,name='wrf_hydro_nwm_jedi_state_rms_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_state
  real(c_float), intent(inout) :: prms

  !type(wrf_hydro_nwm_jedi_state), pointer :: state
  !real(kind=c_float) :: zz
  
  ! call wrf_hydro_nwm_jedi_state_registry%get(c_key_state,state)
  
  ! call rms(state, zz)
  
  ! prms = zz
end subroutine wrf_hydro_nwm_jedi_state_rms_c


! subroutine wrf_hydro_nwm_jedi_state_getvalues_notraj_c(c_key_geom, c_key_state,c_key_loc,c_vars,c_key_gom) bind(c,name='wrf_hydro_nwm_jedi_state_getvalues_notraj_f90')

! implicit none
! integer(c_int), intent(in) :: c_key_state !< State to be interpolated
! integer(c_int), intent(in) :: c_key_loc   !< List of requested locations
! type(c_ptr), value, intent(in) :: c_vars  !< List of requested variables
! integer(c_int), intent(in) :: c_key_gom   !< Interpolated values
! integer(c_int), intent(in) :: c_key_geom  !< Geometry
! type(wrf_hydro_nwm_jedi_state), pointer :: state
! type(ufo_locs),  pointer :: locs
! type(ufo_geovals),  pointer :: gom
! type(oops_variables) :: vars
! type(fv3jedi_geom),  pointer :: geom

! vars = oops_variables(c_vars)

! call fv3jedi_geom_registry%get(c_key_geom, geom)
! call wrf_hydro_nwm_jedi_state_registry%get(c_key_state, state)
! call ufo_locs_registry%get(c_key_loc, locs)
! call ufo_geovals_registry%get(c_key_gom, gom)

! call getvalues(geom, state, locs, vars, gom)

! end subroutine wrf_hydro_nwm_jedi_state_getvalues_notraj_c

! ------------------------------------------------------------------------------

! subroutine wrf_hydro_nwm_jedi_state_getvalues_c(c_key_geom, c_key_state,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='wrf_hydro_nwm_jedi_state_getvalues_f90')

! implicit none
! integer(c_int), intent(in)     :: c_key_state  !< State to be interpolated
! integer(c_int), intent(in)     :: c_key_loc  !< List of requested locations
! type(c_ptr), value, intent(in) :: c_vars     !< List of requested variables
! integer(c_int), intent(in)     :: c_key_gom  !< Interpolated values
! integer(c_int), intent(in)     :: c_key_traj !< Trajectory for interpolation/transforms
! integer(c_int), intent(in)     :: c_key_geom  !< Geometry

! type(wrf_hydro_nwm_jedi_state), pointer :: state
! type(ufo_locs),  pointer :: locs
! type(ufo_geovals),  pointer :: gom
! type(oops_variables) :: vars
! type(fv3jedi_getvalues_traj), pointer :: traj
! type(fv3jedi_geom),  pointer :: geom

! vars = oops_variables(c_vars)

! call wrf_hydro_nwm_jedi_state_registry%get(c_key_state, state)
! call fv3jedi_geom_registry%get(c_key_geom, geom)
! call ufo_locs_registry%get(c_key_loc, locs)
! call ufo_geovals_registry%get(c_key_gom, gom)
! call fv3jedi_getvalues_traj_registry%get(c_key_traj, traj)

! call getvalues(geom, state, locs, vars, gom, traj)

! end subroutine wrf_hydro_nwm_jedi_state_getvalues_c


subroutine wrf_hydro_nwm_jedi_state_sizes_c(c_key_self,nx,ny,nf) &
     bind(c,name='wrf_hydro_nwm_jedi_state_sizes_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(inout) :: nx,ny,nf
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  
  ! call wrf_hydro_nwm_jedi_state_registry%get(c_key_self,self)
  
  ! nf = self%nf
  ! nx = self%npx
  ! ny = self%npy
end subroutine wrf_hydro_nwm_jedi_state_sizes_c


end module wrf_hydro_nwm_jedi_state_interface_mod
