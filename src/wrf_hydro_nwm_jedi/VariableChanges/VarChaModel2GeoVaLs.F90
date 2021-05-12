! (C) Copyright 2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_vc_model2geovals_mod

    use iso_c_binding
    
    use fckit_configuration_module, only: fckit_configuration
    use fckit_log_module,           only: fckit_log
    
    use datetime_mod
    
    use fv3jedi_constants_mod, only: constoz, grav
    use fv3jedi_geom_mod,      only: fv3jedi_geom
    use fv3jedi_fieldfail_mod, only: field_fail
    use fv3jedi_field_mod,     only: copy_subset, field_clen
    use fv3jedi_kinds_mod,     only: kind_real
    use fv3jedi_state_mod,     only: fv3jedi_state
      
    implicit none
    
    private
    public :: fv3jedi_vc_model2geovals
    
    type :: fv3jedi_vc_model2geovals
      integer :: isc, iec, jsc, jec, npz
      character(len=10) :: tropprs_method
      contains
        procedure, public :: create
        procedure, public :: delete
        procedure, public :: changevar
    end type fv3jedi_vc_model2geovals
    
    ! --------------------------------------------------------------------------------------------------
    
    contains
    
    ! --------------------------------------------------------------------------------------------------