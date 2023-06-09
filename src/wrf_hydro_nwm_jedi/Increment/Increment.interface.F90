! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

module wrf_hydro_nwm_jedi_increment_interface_mod

use atlas_module, only: atlas_fieldset
use datetime_mod
use duration_mod
use iso_c_binding
use oops_variables_mod,         only: oops_variables
use fckit_configuration_module, only: fckit_configuration

use wrf_hydro_nwm_jedi_increment_mod, only: &
     increment_create, &
     increment_print, &
     diff_incr, &
     self_mult, &
     axpy_inc, &
     dot_prod, &
     schur_prod, &
     random_normal, &
     dirac, &
     to_fieldset_inc, &
     from_fieldset_inc, &
     to_fieldset_ad_inc, &
     getpoint, &
     setpoint

use wrf_hydro_nwm_jedi_increment_registry_mod, only: &
     wrf_hydro_nwm_jedi_increment_registry
use wrf_hydro_nwm_jedi_state_mod
use wrf_hydro_nwm_jedi_state_interface_mod, only: wrf_hydro_nwm_jedi_state_registry
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_mod_c, only: wrf_hydro_nwm_jedi_geometry_registry
use wrf_hydro_nwm_jedi_geometry_iter_mod_c, only: wrf_hydro_nwm_jedi_geometry_iter_registry
use wrf_hydro_nwm_jedi_geometry_iter_mod, only: wrf_hydro_nwm_jedi_geometry_iter

!GetValues
use ufo_locations_mod
use ufo_geovals_mod
use ufo_geovals_mod_c,      only: ufo_geovals_registry

implicit none

private


contains


subroutine wrf_hydro_nwm_jedi_increment_create_c( &
     c_key_self, c_key_geom, c_vars) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_create_f90')
  implicit none
  integer(c_int),     intent(inout) :: c_key_self
  integer(c_int),     intent(   in) :: c_key_geom  !< Geometry
  type(c_ptr), value, intent(   in) :: c_vars      !< List of variables

  type(wrf_hydro_nwm_jedi_state),    pointer :: self
  type(wrf_hydro_nwm_jedi_geometry), pointer :: geom
  type(oops_variables)                       :: vars

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_increment_registry%init()
  call wrf_hydro_nwm_jedi_increment_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)

  vars = oops_variables(c_vars)
  call increment_create(self, geom, vars)
end subroutine wrf_hydro_nwm_jedi_increment_create_c


subroutine wrf_hydro_nwm_jedi_increment_create_from_other_c( &
     c_key_self, c_key_other) &
     bind(c,name='wrf_hydro_nwm_jedi_increment_create_from_other_f90')
  implicit none
  integer(c_int),intent(inout) :: c_key_self
  integer(c_int),intent(   in) :: c_key_other

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: other

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_other, other)
  call wrf_hydro_nwm_jedi_increment_registry%init()
  call wrf_hydro_nwm_jedi_increment_registry%add(c_key_self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)

  call create_from_other(self, other)
end subroutine wrf_hydro_nwm_jedi_increment_create_from_other_c


subroutine wrf_hydro_nwm_jedi_increment_read_file_c( &
     c_key_geom, c_key_inc, c_conf, c_dt) &
     bind(c,name='wrf_hydro_nwm_jedi_increment_read_file_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_geom !< Geometry
  integer(c_int), intent(in) :: c_key_inc  !< Increment
  type(c_ptr), intent(in)    :: c_conf     !< Configuration
  type(c_ptr), intent(inout) :: c_dt       !< DateTime

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(datetime) :: f_dt
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom

  write(*,*) "Key_geom from read_increment_from_file", c_key_geom

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, self)
  call c_f_datetime(c_dt, f_dt)
  call read_state_from_file(self, c_conf, f_dt)
  call f_c_datetime(f_dt, c_dt)
end subroutine wrf_hydro_nwm_jedi_increment_read_file_c


subroutine wrf_hydro_nwm_jedi_increment_write_file_c( &
     c_key_geom, c_key_inc, c_conf, c_dt) &
     bind(c,name='wrf_hydro_nwm_jedi_increment_write_file_f90')

  implicit none
  integer(c_int), intent(in) :: c_key_geom !< Geometry
  integer(c_int), intent(in) :: c_key_inc  !< Increment
  type(c_ptr), intent(in)    :: c_conf     !< Configuration
  type(c_ptr), intent(inout) :: c_dt       !< DateTime

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(datetime) :: fdate
  type(wrf_hydro_nwm_jedi_geometry),  pointer :: geom

  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, self)
  call c_f_datetime(c_dt, fdate)
  call write_state_to_file(self, c_conf, fdate)
end subroutine wrf_hydro_nwm_jedi_increment_write_file_c


subroutine wrf_hydro_nwm_jedi_increment_random_c(c_key_self) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_random_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call random_normal(self)
  ! call increment_print(self)
end subroutine wrf_hydro_nwm_jedi_increment_random_c


subroutine wrf_hydro_nwm_jedi_increment_copy_c(c_key_self, c_key_rhs) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_copy_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: rhs

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_rhs, rhs)

  call copy(self, rhs) !Implemented in State
end subroutine wrf_hydro_nwm_jedi_increment_copy_c


subroutine wrf_hydro_nwm_jedi_increment_print_c(c_key_self, string) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_print_f90')
  implicit none
  integer(c_int),               intent( in) :: c_key_self
  character(len=1,kind=c_char), intent(out) :: string(8192)

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call increment_print(self, string=string)
end subroutine wrf_hydro_nwm_jedi_increment_print_c


subroutine wrf_hydro_nwm_jedi_increment_add_c(c_key_self, c_key_rhs) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_add_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: rhs

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_rhs, rhs)

  call add_incr(self, rhs) !Implemented in State
end subroutine wrf_hydro_nwm_jedi_increment_add_c


subroutine wrf_hydro_nwm_jedi_increment_sub_c(c_key_self, c_key_rhs) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_sub_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: rhs

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_rhs, rhs)
  call diff_incr(self, self, rhs)
end subroutine wrf_hydro_nwm_jedi_increment_sub_c


subroutine wrf_hydro_nwm_jedi_increment_mul_c(c_key_self, c_zz) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_mul_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  real(c_float),  intent(in) :: c_zz

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self,self)
  call self_mult(self, c_zz)
end subroutine wrf_hydro_nwm_jedi_increment_mul_c


subroutine wrf_hydro_nwm_jedi_increment_axpy_c( &
     c_key_self, c_aa, c_key_yy) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_axpy_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_aa
  integer(c_int), intent(in) :: c_key_yy

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: yy
  real(c_float) :: scalar
  scalar = c_aa

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_yy, yy)

  call axpy_inc(self, scalar, yy)
end subroutine wrf_hydro_nwm_jedi_increment_axpy_c


subroutine wrf_hydro_nwm_jedi_increment_accumul_c( &
     c_key_self, c_aa, c_key_yy) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_accumul_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  real(c_double), intent(in) :: c_aa
  integer(c_int), intent(in) :: c_key_yy

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: yy
  real(c_float) ::  aa_float

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_state_registry%get(c_key_yy, yy)
  aa_float = c_aa
  call axpy_inc(self, aa_float, yy)
end subroutine wrf_hydro_nwm_jedi_increment_accumul_c


subroutine wrf_hydro_nwm_jedi_increment_dot_prod_c( &
     c_key_inc1, c_key_inc2, c_prod &
     ) bind(c, name='wrf_hydro_nwm_jedi_increment_dot_prod_f90')

  implicit none

  integer(c_int), intent(   in) :: c_key_inc1, c_key_inc2
  real(c_double), intent(inout) :: c_prod

  type(wrf_hydro_nwm_jedi_state), pointer :: inc1, inc2

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc1, inc1)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc2, inc2)

  c_prod = 0.d0
  c_prod = dot_prod(inc1, inc2)
end subroutine wrf_hydro_nwm_jedi_increment_dot_prod_c


subroutine wrf_hydro_nwm_jedi_increment_diff_incr_c( &
     c_key_lhs, c_key_x1, c_key_x2) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_diff_incr_f90')
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


subroutine wrf_hydro_nwm_jedi_increment_zero_c(c_key_self) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_zero_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call zeros(self)
end subroutine wrf_hydro_nwm_jedi_increment_zero_c


subroutine wrf_hydro_nwm_jedi_increment_ones_c(c_key_self) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_ones_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call ones(self)
end subroutine wrf_hydro_nwm_jedi_increment_ones_c


subroutine wrf_hydro_nwm_jedi_increment_dirac_c(c_key_self, c_conf) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_dirac_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  type(c_ptr),    intent(in) :: c_conf

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(fckit_configuration) :: f_conf
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  f_conf = fckit_configuration(c_conf)
  call dirac(self, f_conf)
end subroutine wrf_hydro_nwm_jedi_increment_dirac_c


function wrf_hydro_nwm_jedi_increment_rms_c(c_key_inc) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_rms_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_inc  !> Increment registry key
  real(c_double) :: wrf_hydro_nwm_jedi_increment_rms_c  !> return value

  type(wrf_hydro_nwm_jedi_state), pointer :: increment
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, increment)
  wrf_hydro_nwm_jedi_increment_rms_c = increment%fields_obj%rms()
end function wrf_hydro_nwm_jedi_increment_rms_c


subroutine wrf_hydro_nwm_jedi_increment_schur_c(c_key_self, c_key_rhs) &
     bind(c, name='wrf_hydro_nwm_jedi_increment_schur_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_self
  integer(c_int), intent(in) :: c_key_rhs

  type(wrf_hydro_nwm_jedi_state), pointer :: self
  type(wrf_hydro_nwm_jedi_state), pointer :: rhs

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_self, self)
  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_rhs, rhs)
  call schur_prod(self, rhs)
end subroutine wrf_hydro_nwm_jedi_increment_schur_c


! subroutine sw_increment_delete_c(c_key_self) bind(c, name='sw_increment_delete_f90')
!   implicit none
!   integer(c_int), intent(inout) :: c_key_self
!   type(shallow_water_state_type), pointer :: self
!   call sw_increment_registry%get(c_key_self, self)
!   call delete(self)
!   call sw_increment_registry%remove(c_key_self)
! end subroutine sw_increment_delete_c

! subroutine sw_increment_change_resol_c(c_key_inc, c_key_rhs) bind(c, name='sw_increment_change_resol_f90')

!   implicit none

!   integer(c_int), intent(in) :: c_key_inc
!   integer(c_int), intent(in) :: c_key_rhs

!   type(shallow_water_state_type), pointer :: inc, rhs

!   call sw_increment_registry%get(c_key_inc, inc)
!   call sw_increment_registry%get(c_key_rhs, rhs)

!   call change_resol(inc, rhs)
! end subroutine sw_increment_change_resol_c


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


! subroutine wrf_hydro_nwm_jedi_increment_setpoint_c(c_key_self, c_key_geoiter, c_values) bind(c, name='sw_increment_setpoint_f90')
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
! end subroutine wrf_hydro_nwm_jedi_increment_setpoint_c


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

subroutine wrf_hydro_nwm_jedi_increment_to_fieldset_c(c_key_inc, c_key_geom, c_vars, c_afieldset) bind (c, name='wrf_hydro_nwm_jedi_increment_to_fieldset_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_inc
  integer(c_int), intent(in) :: c_key_geom
  type(c_ptr), value, intent(in) :: c_vars
  type(c_ptr), intent(in), value :: c_afieldset

  type(wrf_hydro_nwm_jedi_state), pointer :: inc
  type(wrf_hydro_nwm_jedi_geometry), pointer :: geom
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, inc)
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call to_fieldset_inc(inc, geom, vars, afieldset)

end subroutine wrf_hydro_nwm_jedi_increment_to_fieldset_c

subroutine wrf_hydro_nwm_jedi_increment_to_fieldset_ad_c(c_key_inc, c_key_geom, c_vars, c_afieldset) bind (c, name='wrf_hydro_nwm_jedi_increment_to_fieldset_ad_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_inc
  integer(c_int), intent(in) :: c_key_geom
  type(c_ptr), value, intent(in) :: c_vars
  type(c_ptr), intent(in), value :: c_afieldset

  type(wrf_hydro_nwm_jedi_state), pointer :: inc
  type(wrf_hydro_nwm_jedi_geometry), pointer :: geom
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, inc)
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call to_fieldset_ad_inc(inc, geom, vars, afieldset)

end subroutine wrf_hydro_nwm_jedi_increment_to_fieldset_ad_c

subroutine wrf_hydro_nwm_jedi_increment_from_fieldset_c(c_key_inc, c_key_geom, c_vars, c_afieldset) bind (c, name='wrf_hydro_nwm_jedi_increment_from_fieldset_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_inc
  integer(c_int), intent(in) :: c_key_geom
  type(c_ptr), value, intent(in) :: c_vars
  type(c_ptr), intent(in), value :: c_afieldset

  type(wrf_hydro_nwm_jedi_state), pointer :: inc
  type(wrf_hydro_nwm_jedi_geometry), pointer :: geom
  type(oops_variables) :: vars
  type(atlas_fieldset) :: afieldset

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, inc)
  call wrf_hydro_nwm_jedi_geometry_registry%get(c_key_geom, geom)
  vars = oops_variables(c_vars)
  afieldset = atlas_fieldset(c_afieldset)

  call from_fieldset_inc(inc, geom, vars, afieldset)

end subroutine wrf_hydro_nwm_jedi_increment_from_fieldset_c

!> C++ interface for wrf_hydro_nwm_jedi_increment_mod::wrf_hydro_nwm_jedi_increment::getpoint()
subroutine c_wrf_hydro_nwm_jedi_increment_getpoint(c_key_inc,c_key_iter,values, values_len) bind(c, name='wrf_hydro_nwm_jedi_increment_getpoint_f90')
  implicit none
  integer(c_int), intent(in) :: c_key_inc
  integer(c_int), intent(in) :: c_key_iter
  integer(c_int), intent(in) :: values_len
  real(c_double), intent(inout) :: values(values_len)

  type(wrf_hydro_nwm_jedi_state),      pointer :: inc
  type(wrf_hydro_nwm_jedi_geometry_iter),      pointer :: iter

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc, inc)
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_iter, iter)

  call getpoint(inc, iter, values_len, values)
end subroutine c_wrf_hydro_nwm_jedi_increment_getpoint

!> C++ interface for wrf_hydro_nwm_jedi_increment_mod::wrf_hydro_nwm_jedi_increment::setpoint()
subroutine c_wrf_hydro_nwm_jedi_increment_setpoint(c_key_inc,c_key_iter,values, values_len) bind(c,name='wrf_hydro_nwm_jedi_increment_setpoint_f90')
  integer(c_int), intent(inout) :: c_key_inc
  integer(c_int), intent(in) :: c_key_iter
  integer(c_int), intent(in) :: values_len
  real(c_double),  intent(in) :: values(values_len)

  type(wrf_hydro_nwm_jedi_state),      pointer :: inc
  type(wrf_hydro_nwm_jedi_geometry_iter), pointer :: iter

  call wrf_hydro_nwm_jedi_increment_registry%get(c_key_inc,inc)
  call wrf_hydro_nwm_jedi_geometry_iter_registry%get(c_key_iter,iter)

  call setpoint(inc, iter, values_len, values)

end subroutine c_wrf_hydro_nwm_jedi_increment_setpoint


end module wrf_hydro_nwm_jedi_increment_interface_mod