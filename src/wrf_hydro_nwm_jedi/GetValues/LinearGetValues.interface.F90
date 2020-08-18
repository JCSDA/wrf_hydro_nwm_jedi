! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module wrf_hydro_nwm_jedi_lineargetvalues_interface_mod

! Intrinsic
use iso_c_binding

! oops dependencies
use datetime_mod
use duration_mod
use oops_variables_mod

! ufo dependencies
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry

! self dependency
use wrf_hydro_nwm_jedi_getvalues_mod, only: wrf_hydro_nwm_jedi_getvalues
use wrf_hydro_nwm_jedi_getvalues_interface_mod, only: wrf_hydro_nwm_jedi_getvalues_registry
! wrf_hydro_nwm_jedi dependencies
use wrf_hydro_nwm_jedi_geometry_mod_c, only: wrf_hydro_nwm_jedi_geometry_registry
use wrf_hydro_nwm_jedi_geometry_mod,            only: wrf_hydro_nwm_jedi_geometry
!use wrf_hydro_nwm_jedi_kinds_mod,           only: kind_real
use wrf_hydro_nwm_jedi_state_utils_mod, only: wrf_hydro_nwm_jedi_state_registry
use wrf_hydro_nwm_jedi_increment_registry_mod, only: wrf_hydro_nwm_jedi_increment_registry
use wrf_hydro_nwm_jedi_state_mod,           only: wrf_hydro_nwm_jedi_state

implicit none
private
!public :: wrf_hydro_nwm_jedi_getvalues_registry

! --------------------------------------------------------------------------------------------------

!> Linked list interface
! #define LISTED_TYPE wrf_hydro_nwm_jedi_getvalues
! #include "oops/util/linkedList_i.f"
! type(registry_t) :: wrf_hydro_nwm_jedi_getvalues_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
!#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

! subroutine wrf_hydro_nwm_jedi_getvalues_create_c(c_key_self, c_key_geom, c_key_locs) &
!            bind (c, name='wrf_hydro_nwm_jedi_getvalues_create_f90')

! integer(c_int),     intent(inout) :: c_key_self      !< Key to self
! integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
! integer(c_int),     intent(in)    :: c_key_locs      !< Key to observation locations

! type(wrf_hydro_nwm_jedi_getvalues), pointer :: self
! type(wrf_hydro_nwm_jedi_geometry),      pointer :: geom
! type(ufo_locs),          pointer :: locs

! ! Create object
! call wrf_hydro_nwm_jedi_getvalues_registry%init()
! call wrf_hydro_nwm_jedi_getvalues_registry%add(c_key_self)
! call wrf_hydro_nwm_jedi_getvalues_registry%get(c_key_self, self)

! ! Others
! call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
! call ufo_locs_registry%get(c_key_locs, locs)

! ! Call method
! call self%create(geom, locs)

! end subroutine wrf_hydro_nwm_jedi_getvalues_create_c

  ! --------------------------------------------------------------------------------------------------

 !  subroutine wrf_hydro_nwm_jedi_lineargetvalues_set_trajectory_c(c_key_self, c_key_geom, c_key_state, c_t1, &
!        c_t2, c_key_locs, c_key_geovals) &
!        bind (c,name='wrf_hydro_nwm_jedi_lineargetvalues_set_trajectory_f90')

! integer(c_int), intent(in) :: c_key_self
! integer(c_int), intent(in) :: c_key_geom
! integer(c_int), intent(in) :: c_key_state
! type(c_ptr),    intent(in) :: c_t1
! type(c_ptr),    intent(in) :: c_t2
! integer(c_int), intent(in) :: c_key_locs
! integer(c_int), intent(in) :: c_key_geovals

! type(wrf_hydro_nwm_jedi_getvalues), pointer :: self
! type(wrf_hydro_nwm_jedi_geometry),      pointer :: geom
! type(wrf_hydro_nwm_jedi_state),     pointer :: state
! type(datetime)                         :: t1
! type(datetime)                         :: t2
! type(ufo_locs),                pointer :: locs
! type(ufo_geovals),             pointer :: geovals

! ! Get objects
! call sw_lineargetvalues_registry%get(c_key_self, self)
! call sw_geom_registry%get(c_key_geom, geom)
! call sw_state_registry%get(c_key_state, state)
! call c_f_datetime(c_t1, t1)
! call c_f_datetime(c_t2, t2)
! call ufo_locs_registry%get(c_key_locs, locs)
! call ufo_geovals_registry%get(c_key_geovals, geovals)

! ! Call method
! call self%set_trajectory(geom, state, t1, t2, locs, geovals)

! end subroutine sw_lineargetvalues_set_trajectory_c

! subroutine wrf_hydro_nwm_jedi_getvalues_delete_c(c_key_self) bind (c, name='wrf_hydro_nwm_jedi_getvalues_delete_f90')

! integer(c_int), intent(inout) :: c_key_self !< Key to self

! type(wrf_hydro_nwm_jedi_getvalues), pointer :: self

! ! Get object
! call wrf_hydro_nwm_jedi_getvalues_registry%get(c_key_self, self)

! ! Call method
! call self%delete()

! ! Remove object
! call wrf_hydro_nwm_jedi_getvalues_registry%remove(c_key_self)

! end subroutine wrf_hydro_nwm_jedi_getvalues_delete_c

! --------------------------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_c( &
     c_key_self, c_key_geom, c_key_state, &
     c_t1, c_t2, &
     c_key_locs, c_key_geovals) &
     bind (c, name='wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_f90')

  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_int), intent(in) :: c_key_state
  type(c_ptr),    intent(in) :: c_t1
  type(c_ptr),    intent(in) :: c_t2
  integer(c_int), intent(in) :: c_key_locs
  integer(c_int), intent(in) :: c_key_geovals

  type(wrf_hydro_nwm_jedi_getvalues), pointer :: self
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  type(wrf_hydro_nwm_jedi_state),     pointer :: state
  type(datetime)                   :: t1
  type(datetime)                   :: t2
  type(ufo_locs),          pointer :: locs
  type(ufo_geovals),       pointer :: geovals

  ! Get objects
  call wrf_hydro_nwm_jedi_getvalues_registry%get(c_key_self, self)
  write(*,*) "After getvalues registry"
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  write(*,*) "After geometry registry"
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_state, state)
  write(*,*) "After increment registry"
  call c_f_datetime(c_t1, t1)
  call c_f_datetime(c_t2, t2)
  call ufo_locs_registry%get(c_key_locs, locs)
  write(*,*) "After ufo_locs registry"
  call ufo_geovals_registry%get(c_key_geovals, geovals)
  write(*,*) "After ufo_geovals registry"
  ! Call method
  call self%fill_geovals(geom, state%fields_obj, t1, t2, locs, geovals)
end subroutine wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_c


subroutine wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_ad_c( &
     c_key_self, c_key_geom, c_key_state, &
     c_t1, c_t2, &
     c_key_locs, c_key_geovals) &
     bind (c, name='wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_ad_f90')

  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_geom
  integer(c_int), intent(in) :: c_key_state
  type(c_ptr),    intent(in) :: c_t1
  type(c_ptr),    intent(in) :: c_t2
  integer(c_int), intent(in) :: c_key_locs
  integer(c_int), intent(in) :: c_key_geovals

  type(wrf_hydro_nwm_jedi_getvalues), pointer :: self
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  type(wrf_hydro_nwm_jedi_state),     pointer :: state
  type(datetime)                   :: t1
  type(datetime)                   :: t2
  type(ufo_locs),          pointer :: locs
  type(ufo_geovals),       pointer :: geovals

  ! Get objects
  call wrf_hydro_nwm_jedi_getvalues_registry%get(c_key_self, self)
  write(*,*) "After getvalues registry"
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  write(*,*) "After geometry registry"
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_state, state)
  write(*,*) "After increment registry"
  call c_f_datetime(c_t1, t1)
  call c_f_datetime(c_t2, t2)
  call ufo_locs_registry%get(c_key_locs, locs)
  write(*,*) "After ufo_locs registry"
  call ufo_geovals_registry%get(c_key_geovals, geovals)
  write(*,*) "After ufo_geovals registry"
  ! Call method
  call self%fill_geovals_ad(geom, state%fields_obj, t1, t2, locs, geovals)
end subroutine wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_ad_c


end module wrf_hydro_nwm_jedi_lineargetvalues_interface_mod
