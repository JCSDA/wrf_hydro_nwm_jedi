! (C) Copyright 2019-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_jedi_geometry_mod_c

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use wrf_hydro_jedi_geometry_mod, only: wrf_hydro_jedi_geometry

implicit none
private

! Setup the C/Fortran interface registry
#define LISTED_TYPE wrf_hydro_jedi_geometry
#include "oops/util/linkedList_i.f"
type(registry_t) :: wrf_hydro_jedi_geometry_registry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! Setup the C/Fortran interface registry
#include "oops/util/linkedList_c.f"

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_jedi_geometry_setup(c_key_self, c_conf) bind(c,name='wrf_hydro_jedi_geometry_setup_f90')
  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf
  
  type(wrf_hydro_jedi_geometry), pointer :: self
  
  call wrf_hydro_jedi_geometry_registry%init()
  call wrf_hydro_jedi_geometry_registry%add(c_key_self)
  call wrf_hydro_jedi_geometry_registry%get(c_key_self, self)  
  call self%init(fckit_configuration(c_conf))
  
end subroutine c_wrf_hydro_jedi_geometry_setup

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_jedi_geometry_clone(c_key_self, c_key_other) bind(c,name='wrf_hydro_jedi_geometry_clone_f90')

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in)    :: c_key_other

  type(wrf_hydro_jedi_geometry), pointer :: self, other
  
  call wrf_hydro_jedi_geometry_registry%add(c_key_self)
  call wrf_hydro_jedi_geometry_registry%get(c_key_self , self )
  call wrf_hydro_jedi_geometry_registry%get(c_key_other, other)
  call self%clone(other)

end subroutine c_wrf_hydro_jedi_geometry_clone

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_jedi_geometry_delete(c_key_self) bind(c,name='wrf_hydro_jedi_geometry_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(wrf_hydro_jedi_geometry), pointer :: self
  
  call wrf_hydro_jedi_geometry_registry%get(c_key_self , self )
  call self%delete()
  call wrf_hydro_jedi_geometry_registry%remove(c_key_self)

end subroutine c_wrf_hydro_jedi_geometry_delete

!------------------------------------------------------------------------------

subroutine c_wrf_hydro_jedi_geometry_info(c_key_self,dx,dy) bind(c,name='wrf_hydro_jedi_geometry_info_f90')
  integer(c_int), value, intent(in) :: c_key_self
  real(c_float), intent(out) :: dx,dy

  type(wrf_hydro_jedi_geometry), pointer :: self
  
  call wrf_hydro_jedi_geometry_registry%get(c_key_self , self)
  call self%get_info(dx,dy)

end subroutine c_wrf_hydro_jedi_geometry_info

end module
