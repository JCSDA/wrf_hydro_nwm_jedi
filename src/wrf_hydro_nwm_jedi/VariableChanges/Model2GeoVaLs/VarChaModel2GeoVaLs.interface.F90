! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module wrf_hydro_nwm_jedi_varchamodel2geovals_interface_mod

    use iso_c_binding
    
    use datetime_mod
    
    use fckit_configuration_module, only: fckit_configuration
    
    use wrf_hydro_nwm_jedi_geometry_mod,            only: wrf_hydro_nwm_jedi_geometry
    use wrf_hydro_nwm_jedi_geometry_mod_c,          only: wrf_hydro_nwm_jedi_geometry_registry
    use wrf_hydro_nwm_jedi_increment_registry_mod,  only: wrf_hydro_nwm_jedi_increment_registry
    use wrf_hydro_nwm_jedi_state_mod,               only: wrf_hydro_nwm_jedi_state, copy
    use wrf_hydro_nwm_jedi_state_interface_mod,     only: wrf_hydro_nwm_jedi_state_registry
    use wrf_hydro_nwm_jedi_varchamodel2geovals_mod, only: wrf_hydro_nwm_jedi_varchamodel2geovals
     
    implicit none
    
    private
    public :: wrf_hydro_nwm_jedi_varchamodel2geovals_registry
    
    ! --------------------------------------------------------------------------------------------------
    
#define LISTED_TYPE wrf_hydro_nwm_jedi_varchamodel2geovals
!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"
    
    !> Global registry
    type(registry_t) :: wrf_hydro_nwm_jedi_varchamodel2geovals_registry
    
    ! --------------------------------------------------------------------------------------------------
    
    contains

      
    ! --------------------------------------------------------------------------------------------------
    
!> Linked list implementation
#include "oops/util/linkedList_c.f"
    
    ! --------------------------------------------------------------------------------------------------
    
    subroutine c_wrf_hydro_nwm_jedi_varchamodel2geovals_create(c_key_self, c_key_geom, c_conf) &
               bind (c, name='wrf_hydro_nwm_jedi_varchamodel2geovals_create_f90')
    
    implicit none
    integer(c_int), intent(inout) :: c_key_self
    integer(c_int), intent(in)    :: c_key_geom
    type(c_ptr),    intent(in)    :: c_conf
    
    type(wrf_hydro_nwm_jedi_varchamodel2geovals), pointer :: self
    type(wrf_hydro_nwm_jedi_geometry),            pointer :: geom
    type(fckit_configuration)                             :: conf
    
    ! Linked list
    ! -----------
    call wrf_hydro_nwm_jedi_varchamodel2geovals_registry%init()
    call wrf_hydro_nwm_jedi_varchamodel2geovals_registry%add(c_key_self)
    call wrf_hydro_nwm_jedi_varchamodel2geovals_registry%get(c_key_self, self)
    call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom,geom)
    
    ! APIs
    ! ----
    conf = fckit_configuration(c_conf)
    
    ! Implementation
    ! --------------
    call self%create(geom, conf)
    
    end subroutine c_wrf_hydro_nwm_jedi_varchamodel2geovals_create
    
    ! --------------------------------------------------------------------------------------------------
    
    subroutine c_wrf_hydro_nwm_jedi_varchamodel2geovals_delete(c_key_self) &
               bind (c, name='wrf_hydro_nwm_jedi_varchamodel2geovals_delete_f90')
    
    implicit none
    integer(c_int), intent(inout) :: c_key_self  !< Change variable structure
    
    type(wrf_hydro_nwm_jedi_varchamodel2geovals), pointer :: self
    
    ! Linked list
    ! -----------
    call wrf_hydro_nwm_jedi_varchamodel2geovals_registry%get(c_key_self,self)
    
    ! Implementation
    ! --------------
    call self%delete()
    
    ! Linked list
    ! -----------
    call wrf_hydro_nwm_jedi_varchamodel2geovals_registry%remove(c_key_self)
    
    end subroutine c_wrf_hydro_nwm_jedi_varchamodel2geovals_delete
    
    ! --------------------------------------------------------------------------------------------------
    
    subroutine c_wrf_hydro_nwm_jedi_varchamodel2geovals_changevar(c_key_self, c_key_geom, c_key_xm, c_key_xg) &
               bind (c, name='wrf_hydro_nwm_jedi_varchamodel2geovals_changevar_f90')
    
    implicit none
    integer(c_int), intent(in) :: c_key_self
    integer(c_int), intent(in) :: c_key_geom
    integer(c_int), intent(in) :: c_key_xm
    integer(c_int), intent(in) :: c_key_xg
    
    type(wrf_hydro_nwm_jedi_varchamodel2geovals), pointer :: self
    type(wrf_hydro_nwm_jedi_geometry),            pointer :: geom
    type(wrf_hydro_nwm_jedi_state),               pointer :: xm
    type(wrf_hydro_nwm_jedi_state),               pointer :: xg
    
    ! Linked list
    ! -----------
    call wrf_hydro_nwm_jedi_varchamodel2geovals_registry%get(c_key_self,self)
    call wrf_hydro_nwm_jedi_state_registry%get(c_key_xm,xm)
    call wrf_hydro_nwm_jedi_state_registry%get(c_key_xg,xg)
    call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom,geom)
    
    ! Implementation
    ! --------------
    call self%changevar(geom, xm, xg)
    
    end subroutine c_wrf_hydro_nwm_jedi_varchamodel2geovals_changevar
    
    !-------------------------------------------------------------------------------
!> C++ interface for linear change of variables from geovals to model
!!
!! Only the identity operator is need for the linear variables.
!! \throws abor1_ftn aborts if the field name cannot be in the "getval_name*"
!! section of the variable metadata
    subroutine c_wrf_hydro_nwm_jedi_model2geovals_linear_changevarAD(c_key_geom, c_key_dxin, c_key_dxout) &
    bind(c,name='wrf_hydro_nwm_jedi_model2geovals_linear_changevarAD_f90')
    
    integer(c_int), intent(in) :: c_key_geom, c_key_dxin, c_key_dxout
    
    type(wrf_hydro_nwm_jedi_geometry),            pointer :: geom
    type(wrf_hydro_nwm_jedi_state),           pointer :: dxin, dxout
    
    call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
    call wrf_hydro_nwm_jedi_increment_registry%get(c_key_dxin, dxin)
    call wrf_hydro_nwm_jedi_increment_registry%get(c_key_dxout, dxout)
  
    ! Implementation
    ! --------------
    call copy(dxout, dxin)
  
    end subroutine c_wrf_hydro_nwm_jedi_model2geovals_linear_changevarAD
    ! --------------------------------------------------------------------------------------------------
    
    end module wrf_hydro_nwm_jedi_varchamodel2geovals_interface_mod
    
