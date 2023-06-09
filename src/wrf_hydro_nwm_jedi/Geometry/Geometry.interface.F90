! (C) Copyright 2019-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_geometry_mod_c

use iso_c_binding

use atlas_module, only: atlas_fieldset, atlas_functionspace_pointcloud
use fckit_configuration_module, only: fckit_configuration
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry, wrf_hydro_nwm_jedi_geometry_set_lonlat, wrf_hydro_nwm_jedi_geometry_fill_extra_fields

use oops_variables_mod,          only: oops_variables

implicit none
public

! Setup the C/Fortran interface registry
#define LISTED_TYPE wrf_hydro_nwm_jedi_geometry
#include "oops/util/linkedList_i.f"
type(registry_t) :: wrf_hydro_nwm_jedi_geometry_registry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

! Setup the C/Fortran interface registry
#include "oops/util/linkedList_c.f"

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_nwm_jedi_geometry_setup(c_key_self, c_conf) bind(c,name='wrf_hydro_nwm_jedi_geometry_setup_f90')
  integer(c_int), intent(inout) :: c_key_self
  type(c_ptr),       intent(in) :: c_conf
  
  type(wrf_hydro_nwm_jedi_geometry), pointer :: self
  
  call wrf_hydro_nwm_jedi_geometry_registry%init()
  call wrf_hydro_nwm_jedi_geometry_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self, self)  
  call self%init(fckit_configuration(c_conf))
  
end subroutine c_wrf_hydro_nwm_jedi_geometry_setup

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_set_lonlat_c(c_key_self, c_afieldset, c_include_halo) &
  & bind(c,name='wrf_hydro_nwm_jedi_geometry_set_lonlat_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr), intent(in), value :: c_afieldset
  logical(c_bool), intent(in) :: c_include_halo

  type(wrf_hydro_nwm_jedi_geometry), pointer :: self
  type(atlas_fieldset) :: afieldset
  logical :: include_halo

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self, self)
  afieldset = atlas_fieldset(c_afieldset)
  include_halo = c_include_halo
  call wrf_hydro_nwm_jedi_geometry_set_lonlat(self, afieldset,include_halo)

end subroutine wrf_hydro_nwm_jedi_geometry_set_lonlat_c

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_set_functionspace_pointer_c(c_key_self, c_afunctionspace, c_afunctionspace_incl_halo) &
 & bind(c,name='wrf_hydro_nwm_jedi_geometry_set_functionspace_pointer_f90')
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr), intent(in), value :: c_afunctionspace
  type(c_ptr), intent(in), value :: c_afunctionspace_incl_halo

  type(wrf_hydro_nwm_jedi_geometry), pointer :: self

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self, self)
  self%lsm%afunctionspace = atlas_functionspace_pointcloud(c_afunctionspace)
  self%lsm%afunctionspace_incl_halo = atlas_functionspace_pointcloud(c_afunctionspace_incl_halo)

end subroutine wrf_hydro_nwm_jedi_geometry_set_functionspace_pointer_c

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_geometry_fill_extra_fields_c(c_key_self, c_afieldset) &
 & bind(c,name='wrf_hydro_nwm_jedi_geometry_fill_extra_fields_f90')
  integer(c_int),intent(in) :: c_key_self     !< Geometry
  type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

  type(wrf_hydro_nwm_jedi_geometry), pointer :: self
  type(atlas_fieldset) :: afieldset

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self, self)
  afieldset = atlas_fieldset(c_afieldset)
  call wrf_hydro_nwm_jedi_geometry_fill_extra_fields(self, afieldset)

end subroutine wrf_hydro_nwm_jedi_geometry_fill_extra_fields_c

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_nwm_jedi_geometry_clone(c_key_self, c_key_other) bind(c,name='wrf_hydro_nwm_jedi_geometry_clone_f90')

  integer(c_int), intent(inout) :: c_key_self
  integer(c_int), intent(in)    :: c_key_other

  type(wrf_hydro_nwm_jedi_geometry), pointer :: self, other
  
  call wrf_hydro_nwm_jedi_geometry_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self , self )
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_other, other)
  call self%clone(other)

end subroutine c_wrf_hydro_nwm_jedi_geometry_clone

! -----------------------------------------------------------------------------

subroutine c_wrf_hydro_nwm_jedi_geometry_delete(c_key_self) bind(c,name='wrf_hydro_nwm_jedi_geometry_delete_f90')
  integer(c_int), intent(inout) :: c_key_self

  type(wrf_hydro_nwm_jedi_geometry), pointer :: self
  
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self , self )
  call self%delete()
  call wrf_hydro_nwm_jedi_geometry_registry%remove(c_key_self)

end subroutine c_wrf_hydro_nwm_jedi_geometry_delete

!------------------------------------------------------------------------------

! subroutine c_wrf_hydro_nwm_jedi_geometry_get_nn(c_key_self,lat,long,dim1_idx,dim2_idx) bind(c,name='wrf_hydro_nwm_jedi_geometry_coo_to_grid_f90')
!   integer(c_int), value, intent(in) :: c_key_self
!   real(c_float), intent(in) :: lat,long
!   integer(c_int), intent(out) :: dim1_idx, dim2_idx

!   type(wrf_hydro_nwm_jedi_geometry), pointer :: self

!   call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self , self)
!   call self%get_nn(lat,long,dim1_idx,dim2_idx)

! end subroutine c_wrf_hydro_nwm_jedi_geometry_get_nn

!------------------------------------------------------------------------------

subroutine c_wrf_hydro_nwm_jedi_geometry_info(c_key_self,dx,dy,npx,npy,npz) bind(c,name='wrf_hydro_nwm_jedi_geometry_info_f90')
  integer(c_int), value, intent(in) :: c_key_self
  real(c_float), intent(out) :: dx,dy
  integer(c_int), intent(out) :: npx,npy,npz

  type(wrf_hydro_nwm_jedi_geometry), pointer :: self

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self , self)
  call self%get_lsm_info(dx, dy, npx, npy, npz)

end subroutine c_wrf_hydro_nwm_jedi_geometry_info

!------------------------------------------------------------------------------

subroutine c_wrf_hydro_nwm_jedi_geoval_levels(c_key_self, c_vars, c_nvars, c_nlevels) bind(c,name='wrf_hydro_nwm_jedi_geoval_levels_f90')

  integer(c_int), intent(in)       :: c_key_self
  type(c_ptr), value, intent(in)   :: c_vars
  integer(c_size_t), intent(in)    :: c_nvars            !< size of Variables
  integer(c_size_t), intent(inout) :: c_nlevels(c_nvars)

! local variables
  type(wrf_hydro_nwm_jedi_geometry), pointer :: self
  type(oops_variables)                       :: vars

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_self, self)
  vars = oops_variables(c_vars)
  call self%get_geoval_levels(vars, c_nvars, c_nlevels)

end subroutine c_wrf_hydro_nwm_jedi_geoval_levels

end module
