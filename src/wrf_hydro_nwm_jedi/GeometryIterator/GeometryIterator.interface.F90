! (C) Copyright 2019-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_geometry_iter_mod_c

use iso_c_binding

use wrf_hydro_nwm_jedi_geometry_iter_mod
use wrf_hydro_nwm_jedi_geometry_mod_c, only: wrf_hydro_nwm_jedi_geometry_registry
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry

implicit none
public

! Setup the C/Fortran interface registry
#define LISTED_TYPE wrf_hydro_nwm_jedi_geometry_iter

#include "oops/util/linkedList_i.f"

type(registry_t) :: wrf_hydro_nwm_jedi_geometry_iter_registry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! Setup the C/Fortran interface registry
#include "oops/util/linkedList_c.f"

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_nwm_jedi_geom_iter_setup(c_key_self, c_key_geom, c_iindex, c_jindex) bind(c,name='wrf_hydro_nwm_jedi_geom_iter_setup_f90')
  integer(c_int), intent(inout) :: c_key_self !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_geom !< Geometry
  integer(c_int), intent(   in) :: c_iindex    !< Index
  integer(c_int), intent(   in) :: c_jindex    !< Index
  
  ! Local variables
  type(wrf_hydro_nwm_jedi_geometry_iter), pointer :: self
  type(wrf_hydro_nwm_jedi_geometry),      pointer :: geom
  
  call wrf_hydro_nwm_jedi_geometry_iter_registry%init()
  call wrf_hydro_nwm_jedi_geometry_iter_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_self, self)  
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)

  ! Call Fortran
  call self%init(geom, c_iindex, c_jindex)
  
end subroutine c_wrf_hydro_nwm_jedi_geom_iter_setup

! ------------------------------------------------------------------------------

end module
