! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module wrf_hydro_nwm_jedi_covariance_interface_mod

use iso_c_binding
use wrf_hydro_nwm_jedi_covariance_mod
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_mod_c, only : wrf_hydro_nwm_jedi_geometry_registry
use wrf_hydro_nwm_jedi_increment_registry_mod, only: wrf_hydro_nwm_jedi_increment_registry
use wrf_hydro_nwm_jedi_state_mod
use wrf_hydro_nwm_jedi_state_utils_mod, only: wrf_hydro_nwm_jedi_state_registry
use wrf_hydro_nwm_jedi_state_mod, only: wrf_hydro_nwm_jedi_state

private
public :: wrf_hydro_nwm_jedi_covar_registry
! ------------------------------------------------------------------------------

#define LISTED_TYPE wrf_hydro_nwm_jedi_covar

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: wrf_hydro_nwm_jedi_covar_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
  subroutine c_wrf_hydro_nwm_jedi_b_setup_f90(c_key_self, c_conf, c_vars) &
       & bind (c,name='wrf_hydro_nwm_jedi_b_setup_f90')

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure
  type(c_ptr),    intent(   in) :: c_conf      !< Configuration
  ! integer(c_int), intent(   in) :: c_key_geom  !< Geometry
  type(c_ptr), value, intent(in) :: c_vars !< List of variables

  type(wrf_hydro_nwm_jedi_covar),  pointer :: self
  ! type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  type(oops_variables) :: vars
  
  ! call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_covar_registry%init()
  call wrf_hydro_nwm_jedi_covar_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self, self)

  vars = oops_variables(c_vars)
  
  call wrf_hydro_nwm_jedi_covar_setup(self, c_conf, vars)

end subroutine c_wrf_hydro_nwm_jedi_b_setup_f90

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_b_delete(c_key_self) bind (c,name='wrf_hydro_nwm_jedi_b_delete_f90')

  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure

  type(wrf_hydro_nwm_jedi_covar), pointer :: self

  call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self,self)
  call wrf_hydro_nwm_jedi_covar_delete(self)
  call wrf_hydro_nwm_jedi_covar_registry%remove(c_key_self)

end subroutine wrf_hydro_nwm_jedi_b_delete

! ------------------------------------------------------------------------------

!> Multiply increment vector by covariance

subroutine c_wrf_hydro_nwm_jedi_b_mult(c_key_self, c_key_in, c_key_out) bind(c,name='wrf_hydro_nwm_jedi_b_mult_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_in
  integer(c_int), intent(in) :: c_key_out

  type(wrf_hydro_nwm_jedi_covar), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: xin
  type(wrf_hydro_nwm_jedi_state), pointer :: xout

  call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_in, xin)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_out, xout)

  call wrf_hydro_nwm_jedi_covar_mult(self, xin, xout)

end subroutine c_wrf_hydro_nwm_jedi_b_mult

! ------------------------------------------------------------------------------

!> Generate randomized increment

! subroutine c_wrf_hydro_nwm_jedi_b_randomize(c_key_self, c_key_out) bind(c,name='wrf_hydro_nwm_jedi_b_randomize_f90')

!   implicit none
!   integer(c_int), intent(in) :: c_key_self
!   integer(c_int), intent(in) :: c_key_out

!   type(wrf_hydro_nwm_jedi_covar),                 pointer :: self
!   type(wrf_hydro_nwm_jedi_state), pointer :: xout
!   type(wrf_hydro_nwm_jedi_state)          :: tmp

!   call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self,self)
!   call wrf_hydro_nwm_jedi_increment_registry%get(c_key_out,xout)

!   ! Make a randomized increment
!   call copy(tmp, xout)
!   call random(tmp)

!   ! Apply B to the random increment
!   call sw_covar_mult(self, tmp, xout)

! end subroutine c_wrf_hydro_nwm_jedi_b_randomize

! ------------------------------------------------------------------------------

end module wrf_hydro_nwm_jedi_covariance_interface_mod
