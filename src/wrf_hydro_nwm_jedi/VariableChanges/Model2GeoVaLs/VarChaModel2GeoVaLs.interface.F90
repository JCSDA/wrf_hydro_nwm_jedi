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
    use wrf_hydro_nwm_jedi_geometry_mod_c,          only : wrf_hydro_nwm_jedi_geometry_registry
    use wrf_hydro_nwm_jedi_state_mod,               only: wrf_hydro_nwm_jedi_state
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
    
    ! --------------------------------------------------------------------------------------------------
    
    end module wrf_hydro_nwm_jedi_varchamodel2geovals_interface_mod
    