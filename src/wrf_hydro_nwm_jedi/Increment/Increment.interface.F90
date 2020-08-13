! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module wrf_hydro_nwm_jedi_increment_interface_mod

use datetime_mod
use duration_mod
use iso_c_binding
use oops_variables_mod,         only: oops_variables
use fckit_configuration_module, only: fckit_configuration

use datetime_mod
use duration_mod
use iso_c_binding, only: c_int, c_float, c_ptr, c_char
use oops_variables_mod
use fckit_configuration_module, only: fckit_configuration

use wrf_hydro_nwm_jedi_state_mod
use wrf_hydro_nwm_jedi_state_interface_mod, only: wrf_hydro_nwm_jedi_state_registry
use wrf_hydro_nwm_jedi_increment_mod, only: diff_incr
use wrf_hydro_nwm_jedi_increment_registry_mod, only: wrf_hydro_nwm_jedi_increment_registry
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_mod_c, only: wrf_hydro_nwm_jedi_geometry_registry

!use unstructured_grid_mod,     only: unstructured_grid, unstructured_grid_registry

!GetValues
use ufo_locs_mod
use ufo_locs_mod_c,         only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c,      only: ufo_geovals_registry

implicit none

private

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_increment_create_c(c_key_self, c_key_geom, c_vars) bind(c, name='wrf_hydro_nwm_jedi_increment_create_f90')

  implicit none
  integer(c_int),     intent(inout) :: c_key_self
  integer(c_int),     intent(   in) :: c_key_geom !< Geometry
  type(c_ptr), value, intent(   in) :: c_vars     !< List of variables

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom
  type(oops_variables)         :: vars

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_increment_registry%init()
  call wrf_hydro_nwm_jedi_increment_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)

  vars = oops_variables(c_vars)
  write(*,*) "Invoking create from Increment ",trim(vars%variable(1))
  call create(self, geom, vars)

end subroutine wrf_hydro_nwm_jedi_increment_create_c

! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_increment_create_from_other_c(c_key_self, c_key_other) bind(c,name='wrf_hydro_nwm_jedi_increment_create_from_other_f90')

  implicit none
  integer(c_int),intent(inout) :: c_key_self  !< Fields
  integer(c_int),intent(   in) :: c_key_other !< Other fields

  ! Local variables
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: other

  ! Interface
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_other, other)
  call wrf_hydro_nwm_jedi_increment_registry%init()
  call wrf_hydro_nwm_jedi_increment_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)

  ! Call Fortran
  call create_from_other(self, other)

end subroutine wrf_hydro_nwm_jedi_increment_create_from_other_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_delete_c(c_key_self) bind(c, name='sw_increment_delete_f90')

!   implicit none

!   integer(c_int), intent(inout) :: c_key_self

!   type(shallow_water_state_type), pointer :: self

!   call sw_increment_registry%get(c_key_self, self)

!   call delete(self)

!   call sw_increment_registry%remove(c_key_self)

! end subroutine sw_increment_delete_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_zero_c(c_key_self) bind(c, name='sw_increment_zero_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self

!   type(shallow_water_state_type), pointer :: self

!   call sw_increment_registry%get(c_key_self, self)
!   call zeros(self)

! end subroutine sw_increment_zero_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_dirac_c(c_key_self, c_conf, c_key_geom) bind(c, name='sw_increment_dirac_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   type(c_ptr),    intent(in) :: c_conf     !< Configuration
!   integer(c_int), intent(in) :: c_key_geom !< Geometry

!   type(shallow_water_state_type),    pointer :: self
!   type(shallow_water_geometry_type), pointer :: geom

!   call sw_geom_registry%get(c_key_geom, geom)
!   call sw_increment_registry%get(c_key_self, self)
!   call dirac(self, c_conf, geom)

! end subroutine sw_increment_dirac_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_random_c(c_key_self) bind(c, name='sw_increment_random_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self

!   type(shallow_water_state_type), pointer :: self

!   call sw_increment_registry%get(c_key_self, self)
!   call random(self)

! end subroutine sw_increment_random_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_ug_coord_c(c_key_inc, c_key_ug, c_key_geom) bind (c,name='sw_increment_ug_coord_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_inc
!   integer(c_int), intent(in) :: c_key_ug
!   integer(c_int), intent(in) :: c_key_geom !< Geometry

!   type(shallow_water_state_type),     pointer :: inc
!   type(unstructured_grid),            pointer :: ug
!   type(shallow_water_geometry_type),  pointer :: geom

!   call sw_increment_registry%get(c_key_inc,inc)
!   call unstructured_grid_registry%get(c_key_ug,ug)
!   call sw_geom_registry%get(c_key_geom, geom)

!   call ug_coord(inc, ug, geom)

! end subroutine sw_increment_ug_coord_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_increment_to_ug_c(c_key_inc, c_key_ug, c_its) bind (c,name='sw_increment_increment_to_ug_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_inc
!   integer(c_int), intent(in) :: c_key_ug
!   integer(c_int), intent(in) :: c_its

!   type(shallow_water_state_type), pointer :: inc
!   type(unstructured_grid),        pointer :: ug
!   integer                                 :: its

!   its = c_its+1

!   call sw_increment_registry%get(c_key_inc,inc)
!   call unstructured_grid_registry%get(c_key_ug,ug)

!   call increment_to_ug(inc, ug, its)

! end subroutine sw_increment_increment_to_ug_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_increment_from_ug_c(c_key_inc, c_key_ug, c_its) bind (c,name='sw_increment_increment_from_ug_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_inc
!   integer(c_int), intent(in) :: c_key_ug
!   integer(c_int), intent(in) :: c_its

!   type(shallow_water_state_type), pointer :: inc
!   type(unstructured_grid),        pointer :: ug
!   integer                                 :: its

!   its = c_its+1

!   call sw_increment_registry%get(c_key_inc,inc)
!   call unstructured_grid_registry%get(c_key_ug,ug)

!   call increment_from_ug(inc, ug, its)

! end subroutine sw_increment_increment_from_ug_c

! ! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_increment_copy_c(c_key_self, c_key_rhs) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_copy_f90')

  implicit none

  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: rhs

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self,self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_rhs,rhs)
  
  call copy(self, rhs)

end subroutine wrf_hydro_nwm_jedi_increment_copy_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_self_add_c(c_key_self, c_key_rhs) bind(c, name='sw_increment_self_add_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: self
!   type(shallow_water_state_type), pointer :: rhs

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_increment_registry%get(c_key_rhs, rhs)

!   call self_add(self, rhs)

! end subroutine sw_increment_self_add_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_self_schur_c(c_key_self, c_key_rhs) bind(c, name='sw_increment_self_schur_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: self
!   type(shallow_water_state_type), pointer :: rhs

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_increment_registry%get(c_key_rhs, rhs)

!   call self_schur(self, rhs)

! end subroutine sw_increment_self_schur_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_self_sub_c(c_key_self, c_key_rhs) bind(c, name='sw_increment_self_sub_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: self
!   type(shallow_water_state_type), pointer :: rhs

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_increment_registry%get(c_key_rhs, rhs)

!   call self_sub(self, rhs)

! end subroutine sw_increment_self_sub_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_self_mul_c(c_key_self, c_zz) bind(c, name='sw_increment_self_mul_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   real(c_double), intent(in) :: c_zz

!   type(shallow_water_state_type), pointer :: self
!   real(kind=r8kind)                       :: zz

!   call sw_increment_registry%get(c_key_self, self)
!   zz = c_zz

!   call self_mul(self, zz)

! end subroutine sw_increment_self_mul_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_axpy_inc_c(c_key_self, c_zz, c_key_rhs) bind(c, name='sw_increment_axpy_inc_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   real(c_double), intent(in) :: c_zz
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: self
!   type(shallow_water_state_type), pointer :: rhs
!   real(kind=r8kind)                       :: zz

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_increment_registry%get(c_key_rhs, rhs)
!   zz = c_zz

!   call axpy_inc(self, zz, rhs)

! end subroutine sw_increment_axpy_inc_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_axpy_state_c(c_key_self, c_zz, c_key_rhs) bind(c, name='sw_increment_axpy_state_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   real(c_double), intent(in) :: c_zz
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: self
!   type(shallow_water_state_type), pointer :: rhs
!   real(kind=r8kind)                       :: zz

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_state_registry%get(c_key_rhs, rhs)
!   zz = c_zz

!   call axpy_state(self, zz, rhs)

! end subroutine sw_increment_axpy_state_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_dot_prod_c(c_key_inc1, c_key_inc2, c_prod) bind(c, name='sw_increment_dot_prod_f90')

!   implicit none

!   integer(c_int), intent(   in) :: c_key_inc1, c_key_inc2
!   real(c_double), intent(inout) :: c_prod

!   real(kind=r8kind)                       :: zz
!   type(shallow_water_state_type), pointer :: inc1, inc2

!   call sw_increment_registry%get(c_key_inc1, inc1)
!   call sw_increment_registry%get(c_key_inc2, inc2)

!   call dot_prod(inc1, inc2, zz)

!   c_prod = zz

! end subroutine sw_increment_dot_prod_c

! ! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_increment_diff_incr_c(c_key_lhs, c_key_x1, c_key_x2) bind(c, name='wrf_hydro_nwm_jedi_increment_diff_incr_f90')

  implicit none

  integer(c_int), intent(in) :: c_key_lhs
  integer(c_int), intent(in) :: c_key_x1
  integer(c_int), intent(in) :: c_key_x2

  type(wrf_hydro_nwm_jedi_state), pointer :: lhs
  type(wrf_hydro_nwm_jedi_state), pointer :: x1
  type(wrf_hydro_nwm_jedi_state), pointer :: x2

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_lhs, lhs)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_x1, x1)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_x2, x2)

  call diff_incr(lhs, x1, x2)

end subroutine wrf_hydro_nwm_jedi_increment_diff_incr_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_change_resol_c(c_key_inc, c_key_rhs) bind(c, name='sw_increment_change_resol_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_inc
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: inc, rhs

!   call sw_increment_registry%get(c_key_inc, inc)
!   call sw_increment_registry%get(c_key_rhs, rhs)

!   call change_resol(inc, rhs)

! end subroutine sw_increment_change_resol_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_read_file_c(c_key_geom, c_key_inc, c_conf, c_dt) bind(c, name='sw_increment_read_file_f90')

!   implicit none

!   integer(c_int), intent(   in) :: c_key_inc  !< Increment
!   type(c_ptr),    intent(   in) :: c_conf     !< Configuration
!   type(c_ptr),    intent(inout) :: c_dt       !< DateTime
!   integer(c_int), intent(   in) :: c_key_geom !< Geometry

!   type(shallow_water_state_type),    pointer :: inc
!   type(datetime)                             :: fdate
!   type(shallow_water_geometry_type), pointer :: geom

!   call sw_geom_registry%get(c_key_geom, geom)
!   call sw_increment_registry%get(c_key_inc, inc)
!   call c_f_datetime(c_dt, fdate)
!   call read_file(geom, inc, c_conf, fdate)

! end subroutine sw_increment_read_file_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_write_file_c(c_key_geom, c_key_inc, c_conf, c_dt) bind(c, name='sw_increment_write_file_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_inc  !< Increment
!   type(c_ptr),    intent(in) :: c_conf     !< Configuration
!   type(c_ptr),    intent(in) :: c_dt       !< DateTime
!   integer(c_int), intent(in) :: c_key_geom !< Geometry

!   type(shallow_water_state_type),    pointer :: inc
!   type(datetime)                             :: fdate
!   type(shallow_water_geometry_type), pointer :: geom

!   call sw_geom_registry%get(c_key_geom, geom)
!   call sw_increment_registry%get(c_key_inc, inc)
!   call c_f_datetime(c_dt, fdate)
!   call write_file(geom, inc, c_conf, fdate)

! end subroutine sw_increment_write_file_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_gpnorm_c(c_key_inc, kf, pstat) bind(c, name='sw_increment_gpnorm_f90')

!   implicit none

!   integer(c_int), intent(   in) :: c_key_inc
!   integer(c_int), intent(   in) :: kf
!   real(c_double), intent(inout) :: pstat(3*kf)

!   type(shallow_water_state_type), pointer :: inc
!   real(kind=r8kind)                       :: zstat(3, kf)
!   integer                                 :: jj, js, jf

!   call sw_increment_registry%get(c_key_inc, inc)

!   call gpnorm(inc, kf, zstat)
!   jj = 0
!   do jf = 1, kf
!     do js = 1, 3
!       jj = jj + 1
!       pstat(jj) = zstat(js, jf)
!     enddo
!   enddo

! end subroutine sw_increment_gpnorm_c

! ! ------------------------------------------------------------------------------

subroutine wrf_hydro_nwm_jedi_increment_print_c(c_key_self, string) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_print_f90')

  implicit none

  integer(c_int), intent(in) :: c_key_self
  character(len=1,kind=c_char) :: string(8192)
  type(wrf_hydro_nwm_jedi_state), pointer :: self
  
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self,self)
  
  call state_print(self,string)
  
end subroutine wrf_hydro_nwm_jedi_increment_print_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_rms_c(c_key_inc, prms) bind(c, name='sw_increment_rms_f90')

!   implicit none

!   integer(c_int), intent(   in) :: c_key_inc
!   real(c_double), intent(inout) :: prms

!   type(shallow_water_state_type), pointer :: inc
!   real(kind=r8kind)                       :: zz

!   call sw_increment_registry%get(c_key_inc, inc)

!   call rms(inc, zz)

!   prms = zz

! end subroutine sw_increment_rms_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_getpoint_c(c_key_self, c_key_geoiter, c_values) bind(c, name='sw_increment_getpoint_f90')

!   implicit none

!   integer(c_int), intent(   in) :: c_key_self       !< Increment
!   integer(c_int), intent(   in) :: c_key_geoiter    !< GeometryIterator
!   real(c_double), intent(inout) :: c_values(3)      !< Values

!   type(shallow_water_state_type), pointer :: self
!   type(sw_geom_iter),             pointer :: geoiter

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_geom_iter_registry%get(c_key_geoiter, geoiter)

!   call getpoint(self, geoiter, c_values)

! end subroutine sw_increment_getpoint_c

! ! ------------------------------------------------------------------------------

! subroutine qg_fields_setpoint_c(c_key_self, c_key_geoiter, c_values) bind(c, name='sw_increment_setpoint_f90')

!   implicit none

!   ! Passed variables
!   integer(c_int), intent(in) :: c_key_self        !< Increment
!   integer(c_int), intent(in) :: c_key_geoiter     !< Geometry iterator
!   real(c_double), intent(in) :: c_values(3)       !< Values

!   ! Local variables
!   type(shallow_water_state_type), pointer :: self
!   type(sw_geom_iter),             pointer :: geoiter

!   ! Interface
!   call sw_increment_registry%get(c_key_self, self)
!   call sw_geom_iter_registry%get(c_key_geoiter, geoiter)

!   ! Call Fortran
!   call setpoint(self, geoiter, c_values)

! end subroutine qg_fields_setpoint_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_sizes_c(c_key_self, nx, ny, nv) bind(c, name='sw_increment_sizes_f90')

!   implicit none

!   integer(c_int), intent(   in) :: c_key_self
!   integer(c_int), intent(inout) :: nx, ny, nv

!   type(shallow_water_state_type), pointer :: self
!   type(shallow_water_geometry_type)       :: geom

!   call sw_increment_registry%get(c_key_self, self)

!   geom = self%get_geometry()

!   nv = 3
!   nx = geom%get_nx()
!   ny = geom%get_ny()

! end subroutine sw_increment_sizes_c

! ! ------------------------------------------------------------------------------

! subroutine sw_increment_jnormgrad_c(c_key_self, c_key_geom, c_key_state, c_conf) bind(c, name='sw_increment_jnormgrad_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_self
!   integer(c_int), intent(in) :: c_key_geom
!   integer(c_int), intent(in) :: c_key_state
!   type(c_ptr),    intent(in) :: c_conf

!   type(shallow_water_state_type),    pointer :: self
!   type(shallow_water_geometry_type), pointer :: geom
!   type(shallow_water_state_type),    pointer :: state

!   call sw_increment_registry%get(c_key_self, self)
!   call sw_geom_registry%get(c_key_geom, geom)
!   call sw_state_registry%get(c_key_state, state)

!   call jnormgrad(self, geom, state, c_conf)

! end subroutine sw_increment_jnormgrad_c

! ------------------------------------------------------------------------------

end module wrf_hydro_nwm_jedi_increment_interface_mod
