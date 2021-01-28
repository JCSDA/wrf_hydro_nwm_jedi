! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
module wrf_hydro_nwm_jedi_covariance_interface_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use oops_variables_mod

use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_mod_c, only : wrf_hydro_nwm_jedi_geometry_registry
use wrf_hydro_nwm_jedi_state_mod, only: wrf_hydro_nwm_jedi_state, create_from_other
use wrf_hydro_nwm_jedi_state_utils_mod, only: wrf_hydro_nwm_jedi_state_registry
use wrf_hydro_nwm_jedi_increment_registry_mod, only: wrf_hydro_nwm_jedi_increment_registry
use wrf_hydro_nwm_jedi_increment_mod, only: random_normal
use wrf_hydro_nwm_jedi_covariance_mod, only: &
     wrf_hydro_nwm_jedi_covar, &
     wrf_hydro_nwm_jedi_covar_setup, &
     wrf_hydro_nwm_jedi_covar_delete, &
     wrf_hydro_nwm_jedi_covar_mult

private
public :: wrf_hydro_nwm_jedi_covar_registry

! ------------------------------------------------------------------------------

#define LISTED_TYPE wrf_hydro_nwm_jedi_covar

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: wrf_hydro_nwm_jedi_covar_registry


contains


!> Linked list implementation
#include "oops/util/linkedList_c.f"

subroutine c_wrf_hydro_nwm_jedi_b_setup_f90(c_key_self, c_key_bkg, c_conf, c_vars) &
     & bind (c, name='wrf_hydro_nwm_jedi_b_setup_f90')
  implicit none
  integer(c_int),     intent(inout) :: c_key_self  !< Background error covariance structure
  integer(c_int),     intent(in   ) :: c_key_bkg   !< Background
  type(c_ptr),        intent(   in) :: c_conf      !< Configuration
  type(c_ptr), value, intent(   in) :: c_vars      !< List of variables

  type(wrf_hydro_nwm_jedi_covar), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: bkg
  type(oops_variables) :: vars
  type(fckit_configuration) :: f_conf

  ! integer(c_int), intent(   in) :: c_key_geom  !< Geometry
  ! type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  ! call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)

  call wrf_hydro_nwm_jedi_covar_registry%init()
  call wrf_hydro_nwm_jedi_covar_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_bkg, bkg)
  vars = oops_variables(c_vars)
  f_conf = fckit_configuration(c_conf)
  call wrf_hydro_nwm_jedi_covar_setup(self, bkg, f_conf, vars)
end subroutine c_wrf_hydro_nwm_jedi_b_setup_f90


subroutine wrf_hydro_nwm_jedi_b_delete(c_key_self) &
     bind (c, name='wrf_hydro_nwm_jedi_b_delete_f90')
  implicit none
  integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure

  type(wrf_hydro_nwm_jedi_covar), pointer :: self

  call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self,self)
  call wrf_hydro_nwm_jedi_covar_delete(self)
  call wrf_hydro_nwm_jedi_covar_registry%remove(c_key_self)
end subroutine wrf_hydro_nwm_jedi_b_delete


!> Multiply increment vector by covariance
subroutine c_wrf_hydro_nwm_jedi_b_mult(c_key_self, c_key_in, c_key_out) &
     bind(c, name='wrf_hydro_nwm_jedi_b_mult_f90')
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


!> Generate randomized increment
subroutine c_wrf_hydro_nwm_jedi_b_randomize(c_key_self, c_key_out) &
     bind(c, name='wrf_hydro_nwm_jedi_b_randomize_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_out

  type(wrf_hydro_nwm_jedi_covar), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: xout
  type(wrf_hydro_nwm_jedi_state)          :: tmp

  call wrf_hydro_nwm_jedi_covar_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_out, xout)

  ! Make a randomized increment
  call create_from_other(tmp, xout)
  call random_normal(tmp)

  ! Apply B to the random increment
  call wrf_hydro_nwm_jedi_covar_mult(self, tmp, xout)
end subroutine c_wrf_hydro_nwm_jedi_b_randomize


end module wrf_hydro_nwm_jedi_covariance_interface_mod
