! (C) Copyright 2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_varchamodel2geovals_mod

    use iso_c_binding
    
    use fckit_configuration_module, only: fckit_configuration
    use fckit_log_module,           only: fckit_log
    
    use datetime_mod
    
    use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
    use wrf_hydro_nwm_jedi_state_mod,    only: wrf_hydro_nwm_jedi_state, copy
    use wrf_hydro_nwm_jedi_fields_mod,   only: base_field
      
    implicit none
    
    private
    public :: wrf_hydro_nwm_jedi_varchamodel2geovals
    
    type :: wrf_hydro_nwm_jedi_varchamodel2geovals
      integer :: isc, iec, jsc, jec, npz
      character(len=10) :: tropprs_method
      contains
        procedure, public :: create
        procedure, public :: delete
        procedure, public :: changevar
    end type wrf_hydro_nwm_jedi_varchamodel2geovals
    
    ! --------------------------------------------------------------------------------------------------
    
    contains
    
    ! --------------------------------------------------------------------------------------------------

    subroutine create(self, geom, conf)

      class(wrf_hydro_nwm_jedi_varchamodel2geovals), intent(inout) :: self
      type(wrf_hydro_nwm_jedi_geometry),             intent(in)    :: geom
      type(fckit_configuration),                     intent(in)    :: conf

    end subroutine create

    ! --------------------------------------------------------------------------------------------------

    subroutine delete(self)
      
      class(wrf_hydro_nwm_jedi_varchamodel2geovals), intent(inout) :: self
    
    end subroutine delete

    ! --------------------------------------------------------------------------------------------------

    subroutine changevar(self, geom, xm, xg)
      
      class(wrf_hydro_nwm_jedi_varchamodel2geovals), intent(inout) :: self
      type(wrf_hydro_nwm_jedi_geometry),             intent(in)    :: geom
      type(fv3jedi_state),                           intent(in)    :: xm
      type(fv3jedi_state),                           intent(inout) :: xg

      call copy(xg%fields, xm%fields)

    end subroutine changevar
