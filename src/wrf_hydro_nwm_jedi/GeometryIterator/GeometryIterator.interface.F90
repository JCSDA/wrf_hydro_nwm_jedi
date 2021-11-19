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
!> C++ interface for wrf_hydro_nwm_jedi_geometry_iter_mod::wrf_hydro_nwm_jedi_geometry_iter::setup()
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
!> C++ interface for wrf_hydro_nwm_jedi_geometry_iter_mod::wrf_hydro_nwm_jedi_geometry_iter::clone()
subroutine c_wrf_hydro_nwm_jedi_geom_iter_clone(c_key_self, c_key_other) bind(c, name='wrf_hydro_nwm_jedi_geom_iter_clone_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator

  ! Local variables
  type(wrf_hydro_nwm_jedi_geometry_iter), pointer :: self, other

  ! Interface
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_other, other)
  call wrf_hydro_nwm_jedi_geometry_iter_registry%init()
  call wrf_hydro_nwm_jedi_geometry_iter_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%clone(other)

end subroutine c_wrf_hydro_nwm_jedi_geom_iter_clone

! ------------------------------------------------------------------------------
!> C++ interface for wrf_hydro_nwm_jedi_geometry_iter_mod::wrf_hydro_nwm_jedi_geometry_iter::delete()
subroutine c_wrf_hydro_nwm_jedi_geom_iter_delete(c_key_self) bind(c, name='wrf_hydro_nwm_jedi_geom_iter_delete_f90')
  integer(c_int), intent(inout) :: c_key_self !< Geometry iterator

  ! Clear interface
  call wrf_hydro_nwm_jedi_geometry_iter_registry%remove(c_key_self)

end subroutine c_wrf_hydro_nwm_jedi_geom_iter_delete

! ------------------------------------------------------------------------------
!> C++ interface for wrf_hydro_nwm_jedi_geometry_iter_mod::wrf_hydro_nwm_jedi_geometry_iter::equals()
subroutine c_wrf_hydro_nwm_jedi_geom_iter_equals(c_key_self, c_key_other, c_equals) bind(c, name='wrf_hydro_nwm_jedi_geom_iter_equals_f90')
  integer(c_int), intent(inout) :: c_key_self  !< Geometry iterator
  integer(c_int), intent(   in) :: c_key_other !< Other geometry iterator
  integer(c_int), intent(inout) :: c_equals    !< Equality flag

  ! Local variables
  type(wrf_hydro_nwm_jedi_geometry_iter),pointer :: self,other

  ! Interface
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_other, other)

  ! Call Fortran
  call self%equals(other, c_equals)

end subroutine c_wrf_hydro_nwm_jedi_geom_iter_equals

! ------------------------------------------------------------------------------
!> C++ interface for wrf_hydro_nwm_jedi_geometry_iter_mod::wrf_hydro_nwm_jedi_geometry_iter::current()
subroutine c_wrf_hydro_nwm_jedi_geom_iter_current(c_key_self, c_lon, c_lat) bind(c, name='wrf_hydro_nwm_jedi_geom_iter_current_f90')
  integer(c_int), intent(   in) :: c_key_self !< Geometry iterator
  real(c_double), intent(inout) :: c_lat      !< Latitude
  real(c_double), intent(inout) :: c_lon      !< Longitude

  ! Local variables
  type(wrf_hydro_nwm_jedi_geometry_iter), pointer :: self

  ! Interface
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_self, self)

  ! Call Fortran
  call self%current(c_lon, c_lat)

end subroutine c_wrf_hydro_nwm_jedi_geom_iter_current

end module
