! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_covariance_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use oops_variables_mod

use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_state_mod, only: wrf_hydro_nwm_jedi_state
use wrf_hydro_nwm_jedi_fields_mod, only: base_field

implicit none

private
public :: &
     wrf_hydro_nwm_jedi_covar, &
     wrf_hydro_nwm_jedi_covar_setup, &
     wrf_hydro_nwm_jedi_covar_delete, &
     wrf_hydro_nwm_jedi_covar_mult

!> Fortran derived type to hold configuration data for the background/model covariance
type :: wrf_hydro_nwm_jedi_covar
  type(oops_variables) :: vars
  real :: normfactor
end type wrf_hydro_nwm_jedi_covar


contains


!> Setup for the model's error covariance matrices (B and Q_i)
!> This routine queries the configuration for the parameters that define the
!> covariance matrix, and stores the relevant values in the
!> error covariance structure.
subroutine wrf_hydro_nwm_jedi_covar_setup(self, bkg, f_conf, vars)
  implicit none
  type(wrf_hydro_nwm_jedi_covar),  intent(inout) :: self    !< Covariance structure
  type(wrf_hydro_nwm_jedi_state),  intent(in   ) :: bkg     !< Background
  type(fckit_configuration),       intent(in   ) :: f_conf  !< Configuration
  type(oops_variables),            intent(in   ) :: vars    ! Variables
  ! type(wrf_hydro_nwm_jedi_geometry), intent(   in) :: geom    !< Geometry

  self%vars = vars
  ! Get field normalization factors from configuration
  call f_conf%get_or_die("normfactor", self%normfactor)
end subroutine wrf_hydro_nwm_jedi_covar_setup


subroutine wrf_hydro_nwm_jedi_covar_delete(self)
  implicit none
  type(wrf_hydro_nwm_jedi_covar), intent(inout) :: self  !< Covariance structure
end subroutine wrf_hydro_nwm_jedi_covar_delete


subroutine wrf_hydro_nwm_jedi_covar_mult(self, xin, xout)
  implicit none
  type(wrf_hydro_nwm_jedi_covar), intent(   in) :: self
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: xin
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: xout

  class(base_field), pointer :: in_f
  class(base_field), pointer :: out_f

  integer :: i

  ! Currently applying a simple scalar multiplication to all variables
  ! listed in the yaml file  
  do i = 1, self%vars%nvars()
     call xin%fields_obj%search_field( &
              trim(self%vars%variable(i)), in_f, .true.)
     call xout%fields_obj%search_field( &
              trim(self%vars%variable(i)), out_f, .true.)
     write(*,*) "norm factor"
     write(*,*) self%normfactor
     call out_f%apply_cov(in_f, self%normfactor)
  end do
end subroutine wrf_hydro_nwm_jedi_covar_mult


end module wrf_hydro_nwm_jedi_covariance_mod
