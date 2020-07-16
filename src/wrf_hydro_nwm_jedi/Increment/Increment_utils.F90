! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Registry for increment to work around circular make dependencies
module wrf_hydro_nwm_jedi_increment_registry_mod

use wrf_hydro_nwm_jedi_state_mod, only: wrf_hydro_nwm_jedi_state

implicit none

private
public :: wrf_hydro_nwm_jedi_increment_registry

#define LISTED_TYPE wrf_hydro_nwm_jedi_state

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: wrf_hydro_nwm_jedi_increment_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------


!> Linked list implementation
#include "oops/util/linkedList_c.f"


end module wrf_hydro_nwm_jedi_increment_registry_mod
