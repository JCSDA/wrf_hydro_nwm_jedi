! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module wrf_hydro_nwm_jedi_increment_mod

use atlas_module, only: atlas_fieldset
use iso_c_binding, only: c_char, c_new_line
use fckit_configuration_module, only: fckit_configuration
use datetime_mod
use fckit_mpi_module
use oops_variables_mod

use wrf_hydro_nwm_jedi_fields_mod,    only: wrf_hydro_nwm_jedi_fields, checksame
use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_iter_mod, only: wrf_hydro_nwm_jedi_geometry_iter
use wrf_hydro_nwm_jedi_util_mod,     only: error_handler
use wrf_hydro_nwm_jedi_constants_mod, only: &
     zero_c_double, zero_c_float, one_c_double, one_c_float

use iso_c_binding, only : c_float, c_double
use wrf_hydro_nwm_jedi_state_mod, only: &
     wrf_hydro_nwm_jedi_state, create_from_other, state_print

!GetValues
use ufo_locations_mod
use ufo_geovals_mod,       only: ufo_geovals

implicit none

private
public :: &
     increment_create, &
     increment_print, &
     diff_incr, &
     self_mult, &
     axpy_inc, &
     dot_prod, &
     schur_prod, &
     random_normal, &
     zeros, &
     ones, &
     dirac, &
     set_atlas_inc, &
     to_atlas_inc, &
     from_atlas_inc, &
     to_atlas_ad_inc, &
     getpoint, &
     setpoint
! create_from_other, &
! delete, &
! copy, &
! self_add, &
! self_schur, &
! self_sub, &
! self_mul, &
! read_file, &
! write_file, &
! gpnorm, &
! rms, &
! change_resol, &
! getpoint, &
! setpoint, &
! ug_coord, &
! increment_to_ug, &
! increment_from_ug, &
! dirac, &
! jnormgrad, &


contains

  
subroutine increment_create(self, geom, vars)
  implicit none
  type(wrf_hydro_nwm_jedi_state),    intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in)    :: geom
  type(oops_variables),              intent(in)    :: vars

  self%nf = vars%nvars()
  allocate(self%fields_obj)
  call self%fields_obj%create(geom, vars) !Initializes to zero by default
  ! Initialize all arrays to zero
  ! call self%fields_obj%zeros()
end subroutine increment_create


subroutine increment_print(self, string)
  implicit none
  type(wrf_hydro_nwm_jedi_state),           intent(inout) :: self         !< increment is a state
  character(len=1, kind=c_char),  optional, intent(out  ) :: string(8192) !< The output string

  integer :: ii

  if(present(string)) then
     call self%fields_obj%print_all_fields(string)
  else
     write(*,*) c_new_line
     !                      "Print Increment (C++) -------------------- ";
     write(*,*) c_new_line//"Print Increment (Fortran) ---------------- "
     call self%fields_obj%print_all_fields()
     write(*,*) c_new_line//"End Print Increment (Fortran) ------------ "
     write(*,*) c_new_line//c_new_line
  endif
end subroutine increment_print


subroutine random_normal(self)
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  integer, parameter :: rseed = 1 ! constant for reproducability of tests
    ! NOTE from soca: random seeds are not quite working the way expected,
    !  it is only set the first time normal_distribution() is called with a seed
  integer :: ff

  call self%fields_obj%set_random_normal(rseed)
end subroutine random_normal


subroutine self_mult(self, zz)
  use iso_c_binding, only: c_float
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  real(c_float),                  intent(   in) :: zz

  call self%fields_obj%scalar_mult(zz)
end subroutine self_mult


!> Method for linear transform axpy == a*x+y == scalar*other + self
!> @todo implement the resolution check requested by oops?
subroutine axpy_inc(self, scalar, other_in)
  use iso_c_binding, only: c_float
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  real(kind=c_float),             intent(   in) :: scalar
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: other_in

  type(wrf_hydro_nwm_jedi_state) :: other

  call create_from_other(other, other_in)
  call other%fields_obj%scalar_mult(scalar)  ! = scalar * other
  call self%fields_obj%add_increment(other%fields_obj)  ! = self + (scalar*other)
end subroutine axpy_inc


function dot_prod(self, other) result(the_result)
  use iso_c_binding, only: c_double
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in) :: self
  type(wrf_hydro_nwm_jedi_state), intent(in) :: other
  real(c_double)                             :: the_result

  the_result = zero_c_double
  the_result = self%fields_obj%dot_prod(other%fields_obj)
end function dot_prod


subroutine diff_incr(self, x1, x2)
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: x1
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: x2

  call self%fields_obj%difference(x1%fields_obj, x2%fields_obj)
end subroutine diff_incr


! subroutine delete(self)
!   type(shallow_water_state_type), intent(inout) :: self
! end subroutine delete


subroutine zeros(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  call self%fields_obj%zero()
end subroutine zeros


subroutine ones(self)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  call self%fields_obj%one()
end subroutine ones


subroutine dirac(self, f_conf)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(fckit_configuration) :: f_conf
  call self%fields_obj%dirac(f_conf)
end subroutine dirac


subroutine schur_prod(self, rhs)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_state), intent(   in) :: rhs
  call self%fields_obj%schur_prod(rhs%fields_obj)
end subroutine schur_prod

subroutine set_atlas_inc(self, geom, vars, afieldset, opt_include_halo)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(inout) :: afieldset
  logical :: opt_include_halo
  call self%fields_obj%set_atlas(geom, vars, afieldset, opt_include_halo)
end subroutine set_atlas_inc

subroutine to_atlas_inc(self, geom, vars, afieldset, opt_include_halo)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(inout) :: afieldset
  logical :: opt_include_halo
  call self%fields_obj%to_atlas(geom, vars, afieldset, opt_include_halo)
end subroutine to_atlas_inc

subroutine from_atlas_inc(self, vars, afieldset)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(inout) :: self
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(in) :: afieldset
  call self%fields_obj%from_atlas(vars, afieldset)
end subroutine from_atlas_inc

subroutine to_atlas_ad_inc(self, geom, vars, afieldset)
  implicit none
  type(wrf_hydro_nwm_jedi_state), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(inout) :: afieldset
  call self%fields_obj%to_atlas_ad(geom, vars, afieldset)
end subroutine to_atlas_ad_inc

! ------------------------------------------------------------------------------
!> Get the values at a specific grid point

subroutine getpoint(self, geoiter, values_len, values)
  use iso_c_binding, only: c_double, c_int 
  implicit none
  type(wrf_hydro_nwm_jedi_state),      intent(in)  :: self
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter !< iterator pointing to desired gridpoint
  !> return values for every field in a vertical column
  integer(c_int),   intent(in) :: values_len
  real(c_double), intent(inout) :: values(values_len)
  call self%fields_obj%get_point(geoiter, values_len, values)
end subroutine getpoint

subroutine setpoint(self, geoiter, values_len, values)
  use iso_c_binding, only: c_double, c_int 
  implicit none
  type(wrf_hydro_nwm_jedi_state),      intent(inout)  :: self
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter !< iterator pointing to desired gridpoint
  !> return values for every field in a vertical column
  integer(c_int),   intent(in) :: values_len
  real(c_double), intent(in) :: values(values_len)
  call self%fields_obj%set_point(geoiter, values_len, values)
end subroutine setpoint

end module wrf_hydro_nwm_jedi_increment_mod