! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fields (model states) implementation for wrf_hydro_nwm - jedi integration
module wrf_hydro_nwm_jedi_fields_mod

use atlas_module, only: atlas_field, atlas_fieldset, atlas_real
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding, only: &
       c_int, c_float, c_double, c_char, c_null_char
use datetime_mod
use fckit_mpi_module
use random_mod

use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_geometry_iter_mod, only: wrf_hydro_nwm_jedi_geometry_iter
use wrf_hydro_nwm_jedi_util_mod, only: error_handler, indices, datetime_eq
use wrf_hydro_nwm_jedi_constants_mod, only: &
     unopened_ncid, zero_c_double, zero_c_float, one_c_float
use oops_variables_mod
use netcdf
! use mpp_domains_mod,   only: east, north, center

implicit none

private

public :: &
     wrf_hydro_nwm_jedi_fields, &
     checksame, &
     base_field !, &
! has_field, &
! long_name_to_wrf_hydro_name, &
! pointer_field, &
! pointer_field_array, &
! copy_field_array, &
! allocate_copy_field_array, &
! fields_rms, &
! fields_gpnorm, &
! fields_print, &
! checksame, &
! flip_array_vertical, &
! copy_subset, &
! mean_stddev

!-----------------------------------------------------------------------------
! Abstract definitions


!> Abstract field definition
! @todo should this be public?
type, abstract, public :: base_field
   character(len=32) :: short_name = "null"   !< Short name (to match file name)
   character(len=20) :: wrf_hydro_nwm_name = "null" !< Common name
   character(len=64) :: long_name = "null"    !< More descriptive name
   character(len=32) :: units = "null"        !< Units for the field
   integer :: ncid_index                      !< Index for restart file (1=lsm, 2=hydro)
   logical :: tracer = .false.                !< Tracer or not?
   ! class(cov_obj) :: covariance
 contains
   procedure (print_field_interface),     pass(self), deferred :: print_field
   procedure (print_field_dims_interface),pass(self), deferred :: print_field_dims
   procedure (read_file_interface),       pass(self), deferred :: read_file
   procedure (read_file_interface),       pass(self), deferred :: write_file
   procedure (get_value_field_interface), pass(self), deferred :: get_value
   procedure (set_value_field_interface), pass(self), deferred :: set_value
   procedure (apply_covariance_mult_interface), pass(self), deferred :: apply_cov
   procedure (difference_interface), deferred :: diff
   procedure (zero_interface), deferred :: zero
   procedure (one_interface), deferred :: one
   procedure (dirac_interface), deferred :: dirac
   procedure (set_random_normal_interface), deferred :: set_random_normal
   procedure (dot_prod_interface), deferred :: dot_prod
   procedure (schur_prod_interface), deferred :: schur_prod
   procedure (rms_interface), deferred :: rms
   procedure (add_incr_interface), pass(self), deferred :: add_incr
   procedure (scalar_mult_interface), pass(self), deferred :: scalar_mult
   procedure (set_atlas_interface), deferred :: set_atlas
   procedure (to_atlas_interface), deferred :: to_atlas
   procedure (from_atlas_interface), deferred :: from_atlas
   procedure (get_point_interface), deferred :: get_point
   procedure (set_point_interface), deferred :: set_point
   ! OVERLOADED OPERATORS
   generic, public :: operator(-) => diff
   ! generic, public :: operator(+) => add_incr
end type base_field


abstract interface
   subroutine print_field_interface(self, string, print_array)
     import base_field
     class(base_field),          intent(in ) :: self
     character(len=*), optional, intent(out) :: string
     logical,          optional, intent(in ) :: print_array
   end subroutine print_field_interface

   subroutine print_field_dims_interface(self)
     import base_field
     class(base_field), intent(in) :: self
   end subroutine print_field_dims_interface

   subroutine read_file_interface(self, ncid_vector)
     import base_field
     class(base_field),     intent(inout) :: self
     integer, dimension(2), intent(in) :: ncid_vector
   end subroutine read_file_interface

   subroutine write_file_interface(self, ncid_vector)
     import base_field
     class(base_field),     intent(inout) :: self
     integer, dimension(2), intent(in) :: ncid_vector
   end subroutine write_file_interface

   pure function get_value_field_interface(self, ind) result(val)
     use iso_c_binding, only : c_float
     import base_field, indices
     class(base_field), intent(in) :: self
     type(indices), intent(in) :: ind
     real(kind=c_float) :: val
   end function get_value_field_interface

   subroutine set_value_field_interface(self, ind, val)
     use iso_c_binding, only : c_float
     import base_field, indices
     class(base_field),  intent(inout) :: self
     type(indices),      intent(in)    :: ind
     real(kind=c_float), intent(in)    :: val
   end subroutine set_value_field_interface

   subroutine apply_covariance_mult_interface(self, in_f, scalar)
     use iso_c_binding, only : c_float
     import base_field
     class(base_field), intent(inout) :: self
     class(base_field), intent(in)    :: in_f
     real(kind=c_float),intent(in)    :: scalar
   end subroutine apply_covariance_mult_interface

   function difference_interface(a, b) result(val)
     import base_field
     class(base_field), intent(in) :: a
     class(base_field), intent(in) :: b
     class(base_field), allocatable :: val
   end function difference_interface

   subroutine add_incr_interface(self, inc)
     import base_field
     class(base_field), intent(inout) :: self
     class(base_field), intent(in) :: inc
   end subroutine add_incr_interface

   subroutine scalar_mult_interface(self, scalar)
     use iso_c_binding, only: c_float
     import base_field
     class(base_field), intent(inout) :: self
     real(c_float),     intent(in   ) :: scalar
   end subroutine scalar_mult_interface

   subroutine zero_interface(self)
     import base_field
     class(base_field), intent(inout) :: self
   end subroutine zero_interface

   subroutine one_interface(self)
     import base_field
     class(base_field), intent(inout) :: self
   end subroutine one_interface

   subroutine dirac_interface(self, ix, iy, iz)
     import base_field
     class(base_field), intent(inout) :: self
     integer, intent(in) :: ix
     integer, intent(in) :: iy
     integer, intent(in) :: iz
   end subroutine dirac_interface

   subroutine set_random_normal_interface(self, seed)
     import base_field
     class(base_field), intent(inout) :: self
     integer,           intent(in   ) :: seed
   end subroutine set_random_normal_interface

   function dot_prod_interface(self, other) result(dot_prod)
     use iso_c_binding, only: c_double
     import base_field
     class(base_field), intent(in) :: self
     class(base_field), intent(in) :: other
     real(c_double)                :: dot_prod
   end function dot_prod_interface

   subroutine schur_prod_interface(self, other)
     import base_field
     class(base_field), intent(inout) :: self
     class(base_field), intent(in) :: other
   end subroutine schur_prod_interface

   function rms_interface(self) result(rms)
     use iso_c_binding, only: c_float
     import base_field
     class(base_field), intent(in) :: self
     real(c_float)                 :: rms
   end function rms_interface

   subroutine set_atlas_interface(self, geom, afieldset)
     use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
     use oops_variables_mod
     use atlas_module, only: atlas_fieldset
     import base_field
     class(base_field), intent(in) :: self
     type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
     type(atlas_fieldset), intent(inout) :: afieldset
   end subroutine set_atlas_interface

   subroutine to_atlas_interface(self, geom, afieldset)
     use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
     use oops_variables_mod
     use atlas_module, only: atlas_fieldset
     import base_field
     class(base_field), intent(in) :: self
     type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
     type(atlas_fieldset), intent(inout) :: afieldset
   end subroutine to_atlas_interface

   subroutine from_atlas_interface(self, afieldset)
     use oops_variables_mod
     use atlas_module, only: atlas_fieldset
     import base_field
     class(base_field), intent(inout) :: self
     type(atlas_fieldset), intent(in) :: afieldset
   end subroutine from_atlas_interface

   subroutine get_point_interface(self, geoiter, values_len, values)
    use iso_c_binding, only : c_double, c_int
    use wrf_hydro_nwm_jedi_geometry_iter_mod, only: wrf_hydro_nwm_jedi_geometry_iter
    import base_field
    class(base_field), intent(in) :: self
    integer(c_int),    intent(in) :: values_len
    real(c_double),  intent(inout) :: values(values_len)
    type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter
  end subroutine get_point_interface   

  subroutine set_point_interface(self, geoiter, values_len, values)
    use iso_c_binding, only : c_double, c_int
    use wrf_hydro_nwm_jedi_geometry_iter_mod, only: wrf_hydro_nwm_jedi_geometry_iter
    import base_field
    class(base_field), intent(inout) :: self
    integer(c_int),       intent(in) :: values_len
    real(c_double),     intent(in) :: values(values_len)
    type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter
  end subroutine set_point_interface 

end interface


!> An elementary field. The function is twofold:
! 1) It allows one to create an array of base class object
! 2) It makes it possible to organize the element data structures
! other than simple arrays (e.g. binary search tree, hash tables, etc...)
type, private :: elem_field
   class(base_field), allocatable :: field
end type elem_field


!> Fortran mirror of C class
type, public :: wrf_hydro_nwm_jedi_fields
   integer :: nf
   ! integer :: xdim_len, ydim_len, dim3_len !refer to geometry and depth
   ! real(kind=c_float), allocatable :: array(:,:,:)
   ! logical :: integerfield = .false.
   type(elem_field), allocatable :: fields(:)
 contains
   procedure :: create
   procedure :: search_field
   procedure :: read_fields_from_file
   procedure :: write_fields_to_file
   procedure :: print_single_field
   procedure :: print_all_fields
   procedure :: add_increment
   procedure :: difference
   procedure :: scalar_mult
   procedure :: zero
   procedure :: one
   procedure :: dirac
   procedure :: set_random_normal
   procedure :: rms
   procedure :: dot_prod
   procedure :: schur_prod
   procedure :: set_atlas
   procedure :: to_atlas
   procedure :: from_atlas
   procedure :: get_point
   procedure :: set_point
   ! procedure :: allocate_field
   ! procedure :: equals
   ! procedure :: copy => field_copy
   ! generic :: assignment(=) => equals
   procedure :: deallocate_field
   ! procedure :: mean_stddev
end type wrf_hydro_nwm_jedi_fields


!-----------------------------------------------------------------------------
! Field implementations

type, private, extends(base_field) :: field_1d
   integer :: xdim_len
   real(c_float), allocatable :: array(:)
 contains
   procedure, pass(self) :: print_field => print_field_1d
   procedure, pass(self) :: print_field_dims => print_field_dims_1d
   procedure, pass(self) :: read_file => read_file_1d
   procedure, pass(self) :: write_file => write_file_1d
   procedure, pass(self) :: get_value => get_value_1d
   procedure, pass(self) :: set_value => set_value_1d
   procedure, pass(self) :: fill => fill_field_1d
   procedure, pass(self) :: mean_stddev => mean_stddev_1d
   procedure, pass(self) :: apply_cov => apply_cov_1d
   procedure :: diff => diff_1d
   procedure :: add_incr => add_incr_1d
   procedure :: zero => zero_1d
   procedure :: one => one_1d
   procedure :: dirac => dirac_1d
   procedure :: set_random_normal => set_random_normal_1d
   procedure :: dot_prod => dot_prod_1d
   procedure, pass(self) :: schur_prod => schur_prod_1d
   procedure, pass(self) :: rms => rms_1d
   procedure :: scalar_mult => scalar_mult_1d
   procedure, pass(self) :: set_atlas => set_atlas_1d
   procedure, pass(self) :: to_atlas => to_atlas_1d
   procedure, pass(self) :: from_atlas => from_atlas_1d
   procedure, pass(self) :: get_point => get_point_1d
   procedure, pass(self) :: set_point => set_point_1d
   ! Destructor
   ! final :: destroy_field_1d
end type field_1d


type, private, extends(base_field) :: field_2d
   integer :: xdim_len, ydim_len
   real(c_float), allocatable :: array(:, :)
 contains
   procedure, pass(self) :: print_field => print_field_2d
   procedure, pass(self) :: print_field_dims => print_field_dims_2d
   procedure, pass(self) :: read_file => read_file_2d
   procedure, pass(self) :: write_file => write_file_2d
   procedure, pass(self) :: get_value => get_value_2d
   procedure, pass(self) :: set_value => set_value_2d
   procedure, pass(self) :: fill => fill_field_2d
   procedure, pass(self) :: mean_stddev => mean_stddev_2d
   procedure, pass(self) :: rms => rms_2d
   procedure, pass(self) :: apply_cov => apply_cov_2d
   procedure :: diff => diff_2d
   procedure :: add_incr => add_incr_2d
   procedure :: scalar_mult => scalar_mult_2d
   procedure :: zero => zero_2d
   procedure :: one => one_2d
   procedure :: dirac => dirac_2d
   procedure :: set_random_normal => set_random_normal_2d
   procedure :: dot_prod => dot_prod_2d
   procedure, pass(self) :: schur_prod => schur_prod_2d
   procedure, pass(self) :: set_atlas => set_atlas_2d
   procedure, pass(self) :: to_atlas => to_atlas_2d
   procedure, pass(self) :: from_atlas => from_atlas_2d
   procedure, pass(self) :: get_point => get_point_2d
   procedure, pass(self) :: set_point => set_point_2d
   ! Destructor
   ! final :: destroy_field_2d
end type field_2d


type, private, extends(base_field) :: field_3d
   integer :: xdim_len, ydim_len, zdim_len
   real(c_float), allocatable :: array(:, :, :)
 contains
   procedure, pass(self) :: print_field => print_field_3d
   procedure, pass(self) :: print_field_dims => print_field_dims_3d
   procedure, pass(self) :: read_file => read_file_3d
   procedure, pass(self) :: write_file => write_file_3d
   procedure, pass(self) :: get_value => get_value_3d
   procedure, pass(self) :: set_value => set_value_3d
   procedure, pass(self) :: fill => fill_field_3d
   procedure, pass(self) :: mean_stddev => mean_stddev_3d
   procedure, pass(self) :: rms => rms_3d
   procedure, pass(self) :: apply_cov => apply_cov_3d
   procedure :: diff => diff_3d
   procedure :: add_incr => add_incr_3d
   procedure :: scalar_mult => scalar_mult_3d
   procedure :: zero => zero_3d
   procedure :: one => one_3d
   procedure :: dirac => dirac_3d
   procedure :: set_random_normal => set_random_normal_3d
   procedure :: dot_prod => dot_prod_3d
   procedure, pass(self) :: schur_prod => schur_prod_3d
   procedure, pass(self) :: set_atlas => set_atlas_3d
   procedure, pass(self) :: to_atlas => to_atlas_3d
   procedure, pass(self) :: from_atlas => from_atlas_3d
   procedure, pass(self) :: get_point => get_point_3d
   procedure, pass(self) :: set_point => set_point_3d
   ! Destructor
   ! final :: destroy_field_3d
end type field_3d


contains


!-----------------------------------------------------------------------------
! Create

subroutine create(self, geom, vars)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in)    :: geom
  type(oops_variables),              intent(in)    :: vars

  class(field_1d), allocatable :: tmp_1d_field
  class(field_2d), allocatable :: tmp_2d_field
  class(field_3d), allocatable :: tmp_3d_field

  integer :: var, vcount, nlev

  ! Total fields
  ! ------------
  self%nf = vars%nvars()

  ! Allocate fields structure
  ! -------------------------
  allocate(self%fields(self%nf))

  vcount = 0
  do var = 1, self%nf
     select case (trim(vars%variable(var)))

     case("qlink1")
        vcount = vcount + 1
        allocate(tmp_1d_field)
        call tmp_1d_field%fill( &
             xdim_len=geom%stream%xdim_len, &
             short_name=vars%variable(var), &
             long_name='streamflow', &
             wrf_hydro_nwm_name='qlink1', &
             units='cms', &
             ncid_index=2)
        allocate(self%fields(vcount)%field, source=tmp_1d_field)
        deallocate(tmp_1d_field)

     case("SNEQV")
        vcount = vcount + 1
        allocate(tmp_2d_field)
        call tmp_2d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             short_name=vars%variable(var), &
             long_name='snow_water_equivalent', &
             wrf_hydro_nwm_name='SNEQV', &
             units='mm', &
             ncid_index=1)
        allocate(self%fields(vcount)%field, source=tmp_2d_field)
        deallocate(tmp_2d_field)

      case("swe")
         vcount = vcount + 1
         allocate(tmp_2d_field)
         call tmp_2d_field%fill( &
              xdim_len=geom%lsm%xdim_len, &
              ydim_len=geom%lsm%ydim_len, &
              short_name=vars%variable(var), &
              long_name='snow_water_equivalent', &
              wrf_hydro_nwm_name='SNEQV', &
              units='mm', &
              ncid_index=1)
         allocate(self%fields(vcount)%field, source=tmp_2d_field)
         deallocate(tmp_2d_field)

     case("SNOWH")
        vcount = vcount + 1
        allocate(tmp_2d_field)
        call tmp_2d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             short_name=vars%variable(var), &
             long_name='snow_depth', &
             wrf_hydro_nwm_name='SNOWH', &
             units='m', &
             ncid_index=1)
        allocate(self%fields(vcount)%field, source=tmp_2d_field)
        deallocate(tmp_2d_field)

      case("snow_depth")
         vcount = vcount + 1
         allocate(tmp_2d_field)
         call tmp_2d_field%fill( &
              xdim_len=geom%lsm%xdim_len, &
              ydim_len=geom%lsm%ydim_len, &
              short_name=vars%variable(var), &
              long_name='snow_depth', &
              wrf_hydro_nwm_name='SNOWH', &
              units='m', &
              ncid_index=1)
         allocate(self%fields(vcount)%field, source=tmp_2d_field)
         deallocate(tmp_2d_field)

     case("BULK_SNICE")
        vcount = vcount + 1
        allocate(tmp_2d_field)
        call tmp_2d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             short_name=vars%variable(var), &
             long_name='bulk_snow_ice', &
             wrf_hydro_nwm_name='BULK_SNICE', &
             units='mm', &  ! TODO JLM CHECK
             ncid_index=1)
        allocate(self%fields(vcount)%field, source=tmp_2d_field)
        deallocate(tmp_2d_field)

     case("BULK_SNLIQ")
        vcount = vcount + 1
        allocate(tmp_2d_field)
        call tmp_2d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             short_name=vars%variable(var), &
             long_name='bulk_snow_liquid', &
             wrf_hydro_nwm_name='BULK_SNLIQ', &
             units='mm', &  ! TODO JLM CHECK
             ncid_index=1)
        allocate(self%fields(vcount)%field, source=tmp_2d_field)
        deallocate(tmp_2d_field)

     case("BULK_SNOW_T")
        vcount = vcount + 1
        allocate(tmp_2d_field)
        call tmp_2d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             short_name=vars%variable(var), &
             long_name='bulk_mass_weighted_snow_temp', &
             wrf_hydro_nwm_name='BULK_SNOW_T', &
             units='degK', &  ! TODO JLM CHECK
             ncid_index=1)
        allocate(self%fields(vcount)%field, source=tmp_2d_field)
        deallocate(tmp_2d_field)

     case("LAI")
        vcount = vcount + 1
        allocate(tmp_2d_field)
        call tmp_2d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             short_name=vars%variable(var),&
             long_name='leaf_area', &
             wrf_hydro_nwm_name='LAI', units='m^2m^-2', ncid_index=1)
        allocate(self%fields(vcount)%field, source=tmp_2d_field)
        deallocate(tmp_2d_field)

     case("SNLIQ")
        vcount = vcount + 1
        allocate(tmp_3d_field)
        !MIDDLE DIMENSION HARDCODED
        call tmp_3d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             zdim_len=3, &
             short_name=vars%variable(var),&
             long_name='snow_liquid', &
             wrf_hydro_nwm_name='SNLIQ', &
             units='liter', &
             ncid_index=1) !> @todo: unit invented
        allocate(self%fields(vcount)%field, source=tmp_3d_field)
        deallocate(tmp_3d_field)

     case("SNICE")
        vcount = vcount + 1
        allocate(tmp_3d_field)
        !MIDDLE DIMENSION HARDCODED
        call tmp_3d_field%fill( &
             xdim_len=geom%lsm%xdim_len, &
             ydim_len=geom%lsm%ydim_len, &
             zdim_len=3, &
             short_name=vars%variable(var),&
             long_name='snow_ice', &
             wrf_hydro_nwm_name='SNICE', &
             units='liter', &
             ncid_index=1) !> @todo: unit invented
        allocate(self%fields(vcount)%field, source=tmp_3d_field)
        deallocate(tmp_3d_field)

     case default
        call abor1_ftn("Create: unknown variable "//trim(vars%variable(var)))
     end select
  end do
end subroutine create


! --------------------------------------------------------------------------------------------------
! Deallocatae method

subroutine deallocate_field(self)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(inout) :: self

  ! write(*,*) "Deallocating fields"
  if(allocated(self%fields)) deallocate(self%fields)
  ! self%lalloc = .false.
end subroutine deallocate_field


!-----------------------------------------------------------------------------
! Fill field method

subroutine fill_field_1d(self, &
     xdim_len, &
     short_name, long_name, wrf_hydro_nwm_name, &
     units, tracer, ncid_index)
  implicit none
  class(field_1d),   intent(inout) :: self
  integer,           intent(in)    :: xdim_len
  character(len=*),  intent(in)    :: short_name
  character(len=*),  intent(in)    :: long_name
  character(len=*),  intent(in)    :: wrf_hydro_nwm_name
  character(len=*),  intent(in)    :: units
  logical, optional, intent(in)    :: tracer
  integer,           intent(in)    :: ncid_index

  self%xdim_len = xdim_len
  if(.not.allocated(self%array)) then
     allocate( self%array(xdim_len) )
     self%array = zero_c_float
  else
     call abor1_ftn("Fields.F90.allocate_field: Field already allocated")
  end if
  self%short_name   = trim(short_name)
  self%long_name    = trim(long_name)
  self%wrf_hydro_nwm_name = trim(wrf_hydro_nwm_name)
  self%units        = trim(units)
  self%ncid_index   = ncid_index
end subroutine fill_field_1d


subroutine fill_field_2d(self, &
     xdim_len, ydim_len, &
     short_name, long_name, wrf_hydro_nwm_name, &
     units, tracer, ncid_index)
  implicit none
  class(field_2d),   intent(inout) :: self
  integer,           intent(in)    :: xdim_len, ydim_len
  character(len=*),  intent(in)    :: short_name
  character(len=*),  intent(in)    :: long_name
  character(len=*),  intent(in)    :: wrf_hydro_nwm_name
  character(len=*),  intent(in)    :: units
  logical, optional, intent(in)    :: tracer
  integer,           intent(in)    :: ncid_index

  self%xdim_len = xdim_len
  self%ydim_len = ydim_len
  if(.not.allocated(self%array)) then
     allocate( self%array(xdim_len, ydim_len) )
     self%array = zero_c_float
  else
     call abor1_ftn("Fields.F90.allocate_field: Field already allocated")
  end if
  self%short_name   = trim(short_name)
  self%long_name    = trim(long_name)
  self%wrf_hydro_nwm_name = trim(wrf_hydro_nwm_name)
  self%units        = trim(units)
  self%ncid_index   = ncid_index
end subroutine fill_field_2d


subroutine fill_field_3d(self, &
     xdim_len, ydim_len, zdim_len, &
     short_name, long_name, wrf_hydro_nwm_name, &
     units, tracer, ncid_index)
  implicit none
  class(field_3d),  intent(inout) :: self
  integer,          intent(in)    :: xdim_len, ydim_len, zdim_len
  character(len=*), intent(in)    :: short_name
  character(len=*), intent(in)    :: long_name
  character(len=*), intent(in)    :: wrf_hydro_nwm_name
  character(len=*), intent(in)    :: units
  logical, optional,intent(in)    :: tracer
  integer,          intent(in)    :: ncid_index

  self%xdim_len = xdim_len
  self%ydim_len = ydim_len
  self%zdim_len = zdim_len
  if(.not.allocated(self%array)) then
     allocate( self%array(xdim_len, ydim_len, zdim_len) )
     self%array = zero_c_float
  else
     call abor1_ftn("Fields.F90.allocate_field: Field already allocated")
  end if
  self%short_name   = trim(short_name)
  self%long_name    = trim(long_name)
  self%wrf_hydro_nwm_name = trim(wrf_hydro_nwm_name)
  self%units        = trim(units)
  self%ncid_index   = ncid_index
end subroutine fill_field_3d


! -----------------------------------------------------------------------------
! search field
! This subroutine unifies the long_name_to_wrf_hydro_name and pointer_field routines

subroutine search_field(self, long_name, field_pointer, pass_wrf_hydro_name)
  class(wrf_hydro_nwm_jedi_fields), target, intent(in)  :: self
  character(len=*),                         intent(in)  :: long_name
  class(base_field),               pointer, intent(out) :: field_pointer
  logical,optional :: pass_wrf_hydro_name

  integer :: n
  character(len=255) :: wrf_hydro_nwm_name

  field_pointer => null()

  if(.not.present(pass_wrf_hydro_name)) then
     !Mapping between wrf_hydro variable name and obs name
     select case (long_name)
     case ("swe")
        wrf_hydro_nwm_name = "SNEQV"
     case ("snow_depth")
        wrf_hydro_nwm_name = "SNOWH"
     case ("leaf_area")
        wrf_hydro_nwm_name = "LAI"
     case default
        wrf_hydro_nwm_name = "null"
     end select
  else
     wrf_hydro_nwm_name = long_name
  end if
  !linear search
  do n = 1, self%nf
     if (trim(wrf_hydro_nwm_name) == trim(self%fields(n)%field%wrf_hydro_nwm_name)) then
        field_pointer => self%fields(n)%field
        return
     endif
  enddo

  call abor1_ftn( &
       "wrf_hydro_nwm_jedi_fields_mod.long_name_to_wrf_hydro_nwm_jedi_name" &
       //"long_name "//trim(long_name)//" not found in fields.")
end subroutine search_field


!-----------------------------------------------------------------------------
! Get value method

pure function get_value_1d(self, ind) result(val)
  class(field_1d), intent(in) :: self
  type(indices), intent(in) :: ind

  real(c_float) :: val
  val = self%array(ind%ind_x)
end function get_value_1d


pure function get_value_2d(self, ind) result(val)
  class(field_2d), intent(in) :: self
  type(indices), intent(in) :: ind

  real(c_float) :: val
  val = self%array(ind%ind_x, ind%ind_y)
end function get_value_2d


pure function get_value_3d(self, ind) result(val)
  class(field_3d), intent(in) :: self
  type(indices), intent(in) :: ind

  real(c_float) :: val
  val = self%array(ind%ind_x, ind%ind_y, ind%ind_z)
end function get_value_3d


!-----------------------------------------------------------------------------
! Set value method

subroutine set_value_1d(self, ind, val)
  class(field_1d), intent(inout) :: self
  type(indices),   intent(in)    :: ind
  real(c_float),   intent(in)    :: val

  self%array(ind%ind_x) = val
end subroutine set_value_1d


subroutine set_value_2d(self, ind, val)
  class(field_2d), intent(inout) :: self
  type(indices),   intent(in)    :: ind
  real(c_float),   intent(in)    :: val

  self%array(ind%ind_x, ind%ind_y) = val
end subroutine set_value_2d


subroutine set_value_3d(self, ind, val)
  class(field_3d), intent(inout) :: self
  type(indices),   intent(in)    :: ind
  real(c_float),   intent(in)    :: val

  self%array(ind%ind_x, ind%ind_y, ind%ind_z) = val
end subroutine set_value_3d


!-----------------------------------------------------------------------------
! Equality

subroutine checksame(self, other, method)
  implicit none
  type(wrf_hydro_nwm_jedi_fields), intent(in) :: self
  type(wrf_hydro_nwm_jedi_fields), intent(in) :: other
  character(len=*),                intent(in) :: method
  integer :: var

  if (self%nf .ne. other%nf) then
!   write(*,*) self%nf, other,%nf
!   call abor1_ftn(trim(method)//"(checksame): Different number of fields")
  endif

  do var = 1, self%nf
     if (self%fields(var)%field%wrf_hydro_nwm_name .ne. other%fields(var)%field%wrf_hydro_nwm_name) then
        write(*,*) self%fields(var)%field%wrf_hydro_nwm_name, other%fields(var)%field%wrf_hydro_nwm_name
!        call abor1_ftn(trim(method)//"(checksame): field "//trim(self%fields(var)%field%wrf_hydro_nwm_name)//&
!             " not in the equivalent position in the right hand side")
     endif
  enddo
end subroutine checksame


!-----------------------------------------------------------------------------
! Diff methods

!> This difference subroutine sets self from two passed fields.
subroutine difference(self, f1, f2)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  class(wrf_hydro_nwm_jedi_fields),  intent(in)    :: f1
  class(wrf_hydro_nwm_jedi_fields),  intent(in)    :: f2

  integer :: ff

  do ff = 1, self%nf
     ! - operator is overloaded from the base_field
     self%fields(ff)%field = f1%fields(ff)%field - f2%fields(ff)%field
  end do
end subroutine difference


function diff_1d(a, b) result(val)
  class(field_1d),   intent(in)  :: a
  class(base_field), intent(in)  :: b
  class(base_field), allocatable :: val

  class(field_1d), allocatable :: tmp

  select type(b)
  type is (field_1d)
     allocate(tmp)
     tmp = a
     tmp%array = a%array - b%array
     tmp%xdim_len = a%xdim_len
     allocate(val, source = tmp)
     deallocate(tmp)
  class default
     call abor1_ftn("diff_1d: b not a 1d_field")
  end select
end function diff_1d


function diff_2d(a, b) result(val)
  class(field_2d),   intent(in)  :: a
  class(base_field), intent(in)  :: b
  class(base_field), allocatable :: val

  class(field_2d), allocatable :: tmp

  select type(b)
  type is (field_2d)
     allocate(tmp)
     tmp = a
     tmp%array = a%array - b%array
     tmp%xdim_len = a%xdim_len
     tmp%ydim_len = a%ydim_len
     allocate(val,source=tmp)
     deallocate(tmp)
  class default
     call abor1_ftn("diff_2d: b not a 2d_field")
  end select
end function diff_2d


function diff_3d(a, b) result(val)
  class(field_3d),   intent(in)  :: a
  class(base_field), intent(in)  :: b
  class(base_field), allocatable :: val

  class(field_3d), allocatable :: tmp

  select type(b)
  type is (field_3d)
     allocate(tmp)
     tmp = a
     tmp%array = a%array - b%array
     tmp%xdim_len = a%xdim_len
     tmp%ydim_len = a%ydim_len
     tmp%zdim_len = a%zdim_len
     allocate(val,source=tmp)
     deallocate(tmp)
  class default
     call abor1_ftn("diff_3d: b not a 3d_field")
  end select
end function diff_3d


!-----------------------------------------------------------------------------
! Increment addition

subroutine add_increment(self, inc)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  class(wrf_hydro_nwm_jedi_fields),  intent(in) :: inc

  integer :: ff

  ! write(*,*) "Increment invoked in Fields.F90"
  do ff = 1, self%nf
     call self%fields(ff)%field%add_incr(inc%fields(ff)%field)
  end do
end subroutine add_increment


subroutine add_incr_1d(self, inc)
  class(field_1d),   intent(inout) :: self
  class(base_field), intent(in)    :: inc

  select type(inc)
  type is (field_1d)
     self%array = self%array + inc%array
  class default
     call abor1_ftn("add_incr_1d: inc not a 1d_field")
  end select
end subroutine add_incr_1d


subroutine add_incr_2d(self, inc)
  class(field_2d), intent(inout) :: self
  class(base_field), intent(in) :: inc

  select type(inc)
  type is (field_2d)
     self%array = self%array + inc%array
  class default
     call abor1_ftn("add_incr_2d: inc not a 2d_field")
  end select
end subroutine add_incr_2d


subroutine add_incr_3d(self, inc)
  class(field_3d), intent(inout) :: self
  class(base_field), intent(in) :: inc

  select type(inc)
  type is (field_3d)
     self%array = self%array + inc%array
  class default
     call abor1_ftn("add_incr_3d: inc not a 3d_field")
  end select
end subroutine add_incr_3d


!-----------------------------------------------------------------------------
! Scalar multiplication

subroutine scalar_mult(self, scalar)
  use iso_c_binding, only: c_float
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  real(c_float),                     intent(in   ) :: scalar

  integer :: ff

  do ff = 1, self%nf
     ! write(*,*) ff, scalar
     call self%fields(ff)%field%scalar_mult(scalar)
  end do
end subroutine scalar_mult


subroutine scalar_mult_1d(self, scalar)
  class(field_1d), intent(inout) :: self
  real(c_float),   intent(in   ) :: scalar
  self%array = self%array * scalar
end subroutine scalar_mult_1d


subroutine scalar_mult_2d(self, scalar)
  class(field_2d), intent(inout) :: self
  real(c_float),   intent(in   ) :: scalar
  self%array = self%array * scalar
end subroutine scalar_mult_2d


subroutine scalar_mult_3d(self, scalar)
  class(field_3d), intent(inout) :: self
  real(c_float),   intent(in   ) :: scalar
  self%array = self%array * scalar
end subroutine scalar_mult_3d


!-----------------------------------------------------------------------------
! Zero fields

subroutine zero(self)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self

  integer :: f

  do f = 1, self%nf
     call self%fields(f)%field%zero()
  end do
end subroutine zero


subroutine zero_1d(self)
  class(field_1d), intent(inout) :: self

  self%array = zero_c_float
end subroutine zero_1d


subroutine zero_2d(self)
  class(field_2d), intent(inout) :: self

  self%array = zero_c_float
end subroutine zero_2d


subroutine zero_3d(self)
  class(field_3d), intent(inout) :: self

  self%array = zero_c_float
end subroutine zero_3d


!-----------------------------------------------------------------------------
! One (identity) fields

subroutine one(self)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self

  integer :: f

  do f = 1, self%nf
     call self%fields(f)%field%one()
  end do
end subroutine one


subroutine one_1d(self)
  class(field_1d), intent(inout) :: self

  self%array = one_c_float
end subroutine one_1d


subroutine one_2d(self)
  class(field_2d), intent(inout) :: self

  self%array = one_c_float
end subroutine one_2d


subroutine one_3d(self)
  class(field_3d), intent(inout) :: self

  self%array = one_c_float
end subroutine one_3d


!-----------------------------------------------------------------------------
! Dirac field

subroutine dirac(self, f_conf)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  type(fckit_configuration) :: f_conf
  integer :: ndir, jdir, ff
  integer,allocatable :: ix(:), iy(:), iz(:)
  logical :: found
  character(len=32), allocatable :: field(:)
  character(len=:), allocatable :: str_array(:)

  call f_conf%get_or_die("ndir", ndir)

  if (ndir>0) then
    allocate(ix(ndir))
    allocate(iy(ndir))
    allocate(iz(ndir))
    allocate(field(ndir))
    call f_conf%get_or_die("ix", ix)
    call f_conf%get_or_die("iy", iy)
    call f_conf%get_or_die("iz", iz)
    call f_conf%get_or_die("field",str_array)
    field = str_array
    deallocate(str_array)

    do jdir=1,ndir
      found = .false.
      do ff=1,self%nf
         if (trim(self%fields(ff)%field%long_name)==trim(field(jdir))) then
           call self%fields(ff)%field%dirac(ix(jdir), iy(jdir), iz(jdir))
           found = .true.
         end if
      end do
      if (.not.found) call abor1_ftn("dirac: field "//trim(field(jdir))//" not found")
    enddo
  end if
end subroutine dirac


subroutine dirac_1d(self, ix, iy, iz)
  class(field_1d), intent(inout) :: self
  integer, intent(in) :: ix
  integer, intent(in) :: iy
  integer, intent(in) :: iz

  ! TODO: check ix
  self%array(ix) = one_c_float
end subroutine dirac_1d


subroutine dirac_2d(self, ix, iy, iz)
  class(field_2d), intent(inout) :: self
  integer, intent(in) :: ix
  integer, intent(in) :: iy
  integer, intent(in) :: iz

  ! TODO: check ix, iy
  self%array(ix, iy) = one_c_float
end subroutine dirac_2d


subroutine dirac_3d(self, ix, iy, iz)
  class(field_3d), intent(inout) :: self
  integer, intent(in) :: ix
  integer, intent(in) :: iy
  integer, intent(in) :: iz

  ! TODO: check ix, iy, iz
  self%array(ix, iy, iz) = one_c_float
end subroutine dirac_3d


!-----------------------------------------------------------------------------
! set_random_normal

subroutine set_random_normal(self, seed)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  integer,                           intent(in)    :: seed

  integer :: ff

  do ff = 1, self%nf
     call self%fields(ff)%field%set_random_normal(seed)
  end do
end subroutine set_random_normal


subroutine set_random_normal_1d(self, seed)
  class(field_1d), intent(inout) :: self
  integer,         intent(in)    :: seed

  call normal_distribution( &
       self%array, &
       zero_c_float, one_c_float, seed)
end subroutine set_random_normal_1d


subroutine set_random_normal_2d(self, seed)
  class(field_2d), intent(inout) :: self
  integer,         intent(in)    :: seed

  integer :: jj
  do jj = 1, size(self%array, 2)
     call normal_distribution( &
          self%array(:, jj), &
          zero_c_float, one_c_float, seed)
  end do
end subroutine set_random_normal_2d


subroutine set_random_normal_3d(self, seed)
  class(field_3d), intent(inout) :: self
  integer,         intent(in)    :: seed

  integer :: jj, kk
  do jj = 1, size(self%array, 2)
     do kk = 1, size(self%array, 3)
        call normal_distribution( &
             self%array(:, jj, kk), &
             zero_c_float, one_c_float, seed)
     end do
  end do
end subroutine set_random_normal_3d


! -----------------------------------------------------------------------------
! Dot product

function dot_prod(self, other)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(in) :: self
  class(wrf_hydro_nwm_jedi_fields),  intent(in) :: other
  real(c_double)                                :: dot_prod

  integer :: ff

  dot_prod = zero_c_double
  do ff = 1, self%nf
     dot_prod = dot_prod + self%fields(ff)%field%dot_prod(other%fields(ff)%field)
  end do
end function dot_prod


function dot_prod_1d(self, other)
  implicit none
  class(field_1d),   intent(in) :: self
  class(base_field), intent(in) :: other
  real(c_double)                :: dot_prod_1d

  integer :: ii
  dot_prod_1d = zero_c_double

  select type(other)
  type is (field_1d)
     do ii = 1, size(self%array, 1)
        dot_prod_1d = dot_prod_1d + self%array(ii) * other%array(ii)
     end do
  class default
     call abor1_ftn("dot_prod_1d: other is not a 1d_field")
  end select
end function dot_prod_1d


function dot_prod_2d(self, other)
  implicit none
  class(field_2d),   intent(in) :: self         !> This field
  class(base_field), intent(in) :: other        !> The other field
  real(c_double)                :: dot_prod_2d  !> The return value

  integer :: ii, jj
  dot_prod_2d = zero_c_double

  select type(other)
  type is (field_2d)
     do ii = 1, size(self%array, 1)
        do jj = 1, size(self%array, 2)
          dot_prod_2d = dot_prod_2d + self%array(ii, jj) * other%array(ii, jj)
        end do
     end do
  class default
     call abor1_ftn("dot_prod_2d: other is not a 2d_field")
  end select

  ! write(*,*) "Dot product in 2d_field ", dot_prod_2d
end function dot_prod_2d


function dot_prod_3d(self, other)
  implicit none
  class(field_3d),   intent(in) :: self
  class(base_field), intent(in) :: other
  real(c_double)                :: dot_prod_3d

  integer :: ii, jj, kk
  dot_prod_3d = zero_c_double

  select type(other)
  type is (field_3d)
     do ii = 1, size(self%array, 1)
        do jj = 1, size(self%array, 2)
           do kk = 1, size(self%array, 3)
              dot_prod_3d = dot_prod_3d + self%array(ii, jj, kk) * other%array(ii, jj, kk)
           end do
        end do
     end do
  class default
     call abor1_ftn("dot_prod_3d: other is not a 3d_field")
  end select
end function dot_prod_3d


!-----------------------------------------------------------------------------
! Mean and standard deviation

subroutine mean_stddev_1d(self, mean, stddev)
  implicit none
  class(field_1d), target, intent(in) :: self
  real(c_float), intent(inout) :: mean, stddev

  real(c_double) :: tmp
  integer :: nn, ii

  nn = size(self%array, 1)

  tmp = zero_c_double
  tmp = sum(self%array(:))
  tmp = tmp / nn
  mean = real(tmp, kind=c_float)

  !Computing stddev
  tmp = zero_c_double
  do ii = 1, size(self%array, 1)
        tmp = tmp + (self%array(ii) - mean)**2
  end do
  tmp = sqrt(tmp / nn)
  stddev = real(tmp, kind=c_float)
end subroutine mean_stddev_1d


subroutine mean_stddev_2d(self, mean, stddev)
  implicit none
  class(field_2d), target, intent(in) :: self
  real(c_float), intent(inout) :: mean, stddev

  real(c_double) :: tmp
  integer :: nn, ii, jj

  nn = size(self%array, 1) * size(self%array, 2)

  ! Compute mean
  tmp = zero_c_double
  tmp = sum(self%array(:, :))
  tmp = tmp / nn
  mean = real(tmp, kind=c_float)

  !Computing stddev
  tmp = zero_c_double
  do ii = 1,size(self%array, 1)
     do jj = 1,size(self%array, 2)
        tmp = tmp + (self%array(ii, jj) - mean)**2
     end do
  end do
  tmp = sqrt(tmp / nn)
  stddev = real(tmp, kind=c_float)
end subroutine mean_stddev_2d


subroutine mean_stddev_3d(self, mean, stddev, zlayer)
  implicit none
  class(field_3d), target, intent(in) :: self
  real(c_float),           intent(inout) :: mean,stddev
  integer, optional,       intent(in) :: zlayer
  real(c_double) :: tmp
  integer :: nn, ii, jj, kk

  ! TODO: make optional code for z layers.
  ! nn = size(self%array,1) * size(self%array,3)
  nn = size(self%array,1) * size(self%array,2) * size(self%array,3)

  ! Compute mean
  tmp = zero_c_double
  !tmp = sum(self%array(:,zlayer,:))
  tmp = sum(self%array(:,:,:))
  tmp = tmp / nn
  mean = real(tmp, kind=c_float)

  !Computing stddev
  tmp = zero_c_double
  ! do i=1,size(self%array,1)
  !    do j=1,size(self%array,2)
  !       tmp = tmp + (self%array(i,zlayer,j) - mean)**2
  !    end do
  ! end do
  do ii=1,size(self%array, 1)
     do jj=1,size(self%array, 2)
        do kk=1,size(self%array, 3)
           tmp = tmp + (self%array(ii, jj, kk) - mean)**2
        end do
     end do
  end do

  tmp = sqrt(tmp / nn)
  stddev = real(tmp, kind=c_float)
end subroutine mean_stddev_3d


!-----------------------------------------------------------------------------
! RMS implementations

function rms(self)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(in) :: self
  real(c_float)                                :: rms

  integer :: ff

  ! @todo: jlm this seems bonkers to add all the fields together if there is
  ! more than one
  rms = zero_c_float
  do ff = 1, size(self%fields)
     rms = rms + self%fields(ff)%field%rms()
  end do
end function rms


!> RMS method (1D)
function rms_1d(self) result(rms)
  implicit none
  class(field_1d), intent(in) :: self
  real(kind=c_float) :: rms
  !type(fckit_mpi_comm), intent(in)    :: f_comm

  real(kind=c_double) :: dot_prod_w_self, rms_double

  !do i = 1, self%xdim_len
  !   zz = zz + self%array(i)**2
  !   ii = ii + 1
  !enddo
  dot_prod_w_self = self%dot_prod(self)
  rms_double = sqrt(dot_prod_w_self / (self%xdim_len))
  rms = real(rms_double, c_float)

  ! !Get global values
  ! call f_comm%allreduce(zz,rms,fckit_mpi_sum())
  ! call f_comm%allreduce(ii,iisum,fckit_mpi_sum())
end function rms_1d


function rms_2d(self) result(rms)
  implicit none
  class(field_2d), intent(in) :: self
  real(kind=c_float) :: rms
  !type(fckit_mpi_comm), intent(in)    :: f_comm

  real(kind=c_double) :: dot_prod_w_self, rms_double

  ! do j = 1, self%ydim_len
  !    do i = 1, self%xdim_len
  !       zz = zz + self%array(i, j)**2
  !       ii = ii + 1
  !    enddo
  ! enddo
  dot_prod_w_self = self%dot_prod(self)
  rms_double = sqrt(dot_prod_w_self / (self%xdim_len * self%ydim_len))
  rms = real(rms_double, c_float)

  ! !Get global values
  ! call f_comm%allreduce(zz,rms,fckit_mpi_sum())
  ! call f_comm%allreduce(ii,iisum,fckit_mpi_sum())
end function rms_2d


function rms_3d(self) result(rms)
  implicit none
  class(field_3d),  intent(in) :: self
  ! integer, intent(in) :: zlayer
  real(kind=c_float) :: rms
  ! ?optionally: specify a layer?
  !type(fckit_mpi_comm), intent(in)    :: f_comm

  real(kind=c_double) :: dot_prod_w_self, rms_double

  dot_prod_w_self = self%dot_prod(self)
  rms_double = sqrt(dot_prod_w_self / (self%xdim_len * self%ydim_len * self%zdim_len))
  rms = real(rms_double, c_float)

  ! !Get global values
  ! call f_comm%allreduce(zz,rms,fckit_mpi_sum())
  ! call f_comm%allreduce(ii,iisum,fckit_mpi_sum())
end function rms_3d


! -----------------------------------------------------------------------------
! Schur product

subroutine schur_prod(self, other)
  implicit none
  class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
  class(wrf_hydro_nwm_jedi_fields),  intent(in   ) :: other
  integer :: ff

  do ff = 1, self%nf
     call self%fields(ff)%field%schur_prod( other%fields(ff)%field )
  end do
end subroutine schur_prod


subroutine schur_prod_1d(self, other)
  implicit none
  class(field_1d),   intent(inout) :: self
  class(base_field), intent(in) :: other

  select type(other)
  type is (field_1d)
     self%array = self%array * other%array
  class default
     call abor1_ftn("schur_prod_1d: other is not a 1d_field")
  end select
end subroutine schur_prod_1d


subroutine schur_prod_2d(self, other)
  implicit none
  class(field_2d),   intent(inout) :: self         !> This field
  class(base_field), intent(in) :: other        !> The other field

  select type(other)
  type is (field_2d)
     self%array = self%array * other%array
  class default
     call abor1_ftn("schur_prod_2d: other is not a 2d_field")
  end select
end subroutine schur_prod_2d


subroutine schur_prod_3d(self, other)
  implicit none
  class(field_3d),   intent(inout) :: self
  class(base_field), intent(in) :: other

  select type(other)
  type is (field_3d)
     self%array = self%array * other%array
  class default
     call abor1_ftn("schur_prod_3d: other is not a 3d_field")
  end select
end subroutine schur_prod_3d


!-----------------------------------------------------------------------------
! Apply covariance

subroutine apply_cov_1d(self, in_f, scalar)
  class(field_1d),   intent(inout) :: self
  class(base_field), intent(in)    :: in_f
  real(kind=c_float),intent(in)    :: scalar

  select type(in_f)
  type is (field_1d)
     self%array = in_f%array * scalar
  class default
     call abor1_ftn("apply_cov_1d: in_f not a 1d_field")
  end select
end subroutine apply_cov_1d


subroutine apply_cov_2d(self, in_f, scalar)
  class(field_2d),   intent(inout) :: self
  class(base_field), intent(in)    :: in_f
  real(kind=c_float),intent(in)    :: scalar

  select type(in_f)
  type is (field_2d)
     self%array = in_f%array * scalar
  class default
     call abor1_ftn("apply_cov_2d: in_f not a 2d_field")
  end select
end subroutine apply_cov_2d


subroutine apply_cov_3d(self, in_f, scalar)
  class(field_3d),   intent(inout) :: self
  class(base_field), intent(in)    :: in_f
  real(kind=c_float),intent(in)    :: scalar

  integer :: layer

  select type(in_f)
  type is (field_3d)
     do layer = 1, size(self%array, 2) !zlayer is dim2
        self%array(:,layer,:) = in_f%array(:,layer,:) * scalar
     end do
  class default
     call abor1_ftn("apply_cov_3d: in_f not a 3d_field")
  end select
end subroutine apply_cov_3d


! --------------------------------------------------------------------------------------------------
! Print dimension method

subroutine print_field_dims_1d(self)
  use iso_c_binding, only : c_new_line, c_float
  implicit none
  class(field_1d), intent(in) :: self

  write(*,*) 'Printing size', self%xdim_len
end subroutine print_field_dims_1d


subroutine print_field_dims_2d(self)
  implicit none
  class(field_2d), intent(in) :: self

  write(*,*) 'Printing size ', self%xdim_len, self%ydim_len
end subroutine print_field_dims_2d


subroutine print_field_dims_3d(self)
  implicit none
  class(field_3d), intent(in) :: self

  write(*,*) 'Printing size ', self%xdim_len, self%zdim_len, self%ydim_len
end subroutine print_field_dims_3d


!-----------------------------------------------------------------------------
! Print

subroutine print_single_field(self, long_name, string)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(inout) :: self
  character(len=*),                 intent(in)    :: long_name
  character(len=*), optional,       intent(out)   :: string

  class(base_field), pointer :: field_pointer
  call self%search_field(long_name, field_pointer)
  if(associated(field_pointer)) then
     call field_pointer%print_field(string)
  end if
end subroutine print_single_field


!> Print all fields
subroutine print_all_fields(self, string)
  use iso_c_binding, only : c_new_line
  implicit none
  class(wrf_hydro_nwm_jedi_fields),        intent(in ) :: self
  character(len=1, kind=c_char), optional, intent(out) :: string(8192) !< The output string
  integer :: ff, ii, s_len
  character(len=8192) :: tmp_str, agg_str

  tmp_str = c_null_char
  agg_str = ''
  if(present(string)) then
     do ff = 1, self%nf
        call self%fields(ff)%field%print_field(tmp_str)
        ! this aggregates the string over the passed fields
        agg_str = trim(agg_str) // trim(tmp_str)
     end do
     ! To manually/visually check the aggregation of the string.
     ! write(*,*) "[asdf "//trim(agg_str)//"asdf]"
     s_len = len_trim(agg_str)
     do ii = 1, s_len
        string(ii:ii) = agg_str(ii:ii)
     end do
     string(s_len+1) = c_null_char
  else
     do ff = 1, self%nf
        call self%fields(ff)%field%print_field()
     end do
  end if
end subroutine print_all_fields


!> Print method for a 1d field
subroutine print_field_1d(self, string, print_array)
  use iso_c_binding, only : c_new_line, c_float
  implicit none
  class(field_1d),            intent(in ) :: self
  character(len=*), optional, intent(out) :: string
  logical,          optional, intent(in ) :: print_array

  logical :: print_array_
  real(kind=c_float) :: mean, stddev
  character(len=255) :: float_str1, float_str2, float_str3
  character(len=255) :: int_str1
  character(len=8192) :: message

  if(.not.(present(print_array))) print_array_ = .FALSE.

  call self%mean_stddev(mean, stddev)
  write(float_str1, *) mean
  write(float_str2, *) stddev
  write(float_str3, *) self%rms()
  write(int_str1, *) self%xdim_len

  write(message,*) &
       c_new_line//self%wrf_hydro_nwm_name // &
       c_new_line//self%long_name //  &
       c_new_line//'1D Field shape: ( ' // &
       trim(int_str1) // ' )' // &
       c_new_line//'Mean: '//trim(float_str1) // &
       c_new_line//'Std Dev: '//trim(float_str2) // &
       c_new_line//'RMS: '//trim(float_str3) // &
       c_new_line

  if(present(string)) then
     write(string, *) trim(message)
  else
     write(*, *) trim(message)
     if(print_array_) write(*,*) self%array
  end if
end subroutine print_field_1d


subroutine print_field_2d(self, string, print_array)
  use iso_c_binding, only : c_new_line, c_float
  implicit none
  class(field_2d),            intent(in ) :: self
  character(len=*), optional, intent(out) :: string
  logical,          optional, intent(in ) :: print_array

  logical :: print_array_
  real(kind=c_float) :: mean, stddev
  character(len=255) :: float_str1, float_str2, float_str3
  character(len=255) :: int_str1, int_str2
  character(len=8192) :: message

  if(.not.(present(print_array))) print_array_ = .FALSE.

  call self%mean_stddev(mean, stddev)
  write(float_str1,*) mean
  write(float_str2,*) stddev
  write(float_str3,*) self%rms()

  write(int_str1, *) self%xdim_len
  write(int_str2, *) self%ydim_len

  write(message, *) &
       c_new_line//self%wrf_hydro_nwm_name // &
       c_new_line//self%long_name //  &
       c_new_line//'2D Field shape: ( ' // &
       trim(int_str1) // ', ' // &
       trim(int_str2) // ' )' // &
       c_new_line//'Mean: '//trim(float_str1) // &
       c_new_line//'Std Dev: '//trim(float_str2) // &
       c_new_line//'RMS: '//trim(float_str3) // &
       c_new_line

  if(present(string)) then
     write(string, *) trim(message)
  else
     write(*,*) trim(message)
     if(print_array_) write(*,*) self%array
  end if
end subroutine print_field_2d


subroutine print_field_3d(self, string, print_array)
  use iso_c_binding, only : c_new_line
  implicit none
  class(field_3d),            intent(in ) :: self
  character(len=*), optional, intent(out) :: string
  logical,          optional, intent(in ) :: print_array

  logical :: print_array_
  ! integer :: z_layer
  real(kind=c_float) :: mean, stddev
  character(len=255) :: float_str1, float_str2, float_str3
  character(len=255) :: int_str1, int_str2, int_str3, zlayer_str
  character(len=8192) :: message

  if(.not.(present(print_array))) print_array_ = .FALSE.

  !do z_layer = 1, self%ydim_len  !z is the 2nd dimension...

  call self%mean_stddev(mean, stddev)
  ! , z_layer)

  float_str1 = ' '
  float_str2 = ' '
  float_str3 = ' '
  zlayer_str = ' '
  write(float_str1, *) mean
  write(float_str2, *) stddev
  write(float_str3, *) self%rms()

  write(int_str1, *) self%xdim_len
  write(int_str2, *) self%ydim_len
  write(int_str3, *) self%zdim_len
  ! write(zlayer_str, *) z_layer

  write(message, *) &
          c_new_line//self%wrf_hydro_nwm_name // &
          c_new_line//self%long_name //  &
          c_new_line//'3D Field shape: ( ' // &
          trim(int_str1) // ', ' // &
          trim(int_str2) // ', ' // &
          trim(int_str3) // ' )' // &
          c_new_line//'Mean: '//trim(float_str1) // &
          c_new_line//'Std Dev: '//trim(float_str2) // &
          c_new_line//'RMS: '//trim(float_str3) // &
          c_new_line

  if(present(string)) then
     write(string,*) trim(message)
  else
     write(*,*) trim(message)
     if(print_array_) write(*,*) self%array
  end if
end subroutine print_field_3d


!-----------------------------------------------------------------------------
! File related subroutines/functions

subroutine read_fields_from_file(self, filename_lsm, filename_hydro, f_dt)
  class(wrf_hydro_nwm_jedi_fields), target, intent(inout) :: self
  character(len=*),                         intent(in   ) :: filename_lsm, filename_hydro
  type(datetime),                           intent(out  ) :: f_dt

  integer :: ncid_lsm, ncid_hydro, n, ierr
  integer, dimension(2) :: ncid_vector
  logical :: read_file_lsm, read_file_hydro
  type(datetime)  :: file_time_lsm, file_time_hydro

  ! open files here and pass ncids
  ! check if lsm geom is defined: open the RESTART file
  ncid_lsm = open_get_restart_ncid(self, filename_lsm=filename_lsm)
  ! check if any hydro variables are defined: open the HYDRO_RST file
  ncid_hydro = open_get_restart_ncid(self, filename_hydro=filename_hydro)

  ! verify that the files are at the same time if both are being used.
  file_time_lsm = get_lsm_file_time(ncid_lsm)
  !file_time_hydro = get_lsm_file_time(ncid_hydro)

  !if (datetime_eq(file_time_lsm, file_time_hydro)) then
  f_dt = file_time_lsm
  !else
  !   call abor1_ftn("read_fields_from_file: lsm and hydryo restart files have different times")
  !end if

  ncid_vector = (/ ncid_lsm, ncid_hydro /)
  do n = 1, self%nf
     call self%fields(n)%field%read_file(ncid_vector)
  enddo

  ! close nc files
  call close_restart_ncid(ncid_lsm, filename_lsm)
  call close_restart_ncid(ncid_hydro, filename_hydro)
end subroutine read_fields_from_file


subroutine write_fields_to_file(self, filename_lsm, filename_hydro, f_dt)
  class(wrf_hydro_nwm_jedi_fields), target, intent(inout) :: self
  character(len=*),                         intent(in   ) :: filename_lsm, filename_hydro
  type(datetime),                           intent(out  ) :: f_dt

  integer :: ncid_lsm, ncid_hydro, n, ierr
  integer, dimension(2) :: ncid_vector
  logical :: write_file_lsm, write_file_hydro
  type(datetime)  :: file_time_lsm, file_time_hydro

  ! open files here and pass ncids
  ! check if lsm geom is defined: open the RESTART file
  ncid_lsm = open_get_restart_ncid(self, filename_lsm=filename_lsm, mode=nf90_write)
  ! check if any hydro variables are defined: open the HYDRO_RST file
  ncid_hydro = open_get_restart_ncid(self, filename_hydro=filename_hydro, mode=nf90_write)

  ! verify that the files are at the same time if both are being used.
  file_time_lsm = get_lsm_file_time(ncid_lsm)
  !file_time_hydro = get_lsm_file_time(ncid_hydro)

  !if (datetime_eq(file_time_lsm, file_time_hydro)) then
  f_dt = file_time_lsm
  !else
  !   call abor1_ftn("read_fields_from_file: lsm and hydryo restart files have different times")
  !end if

  ncid_vector = (/ ncid_lsm, ncid_hydro /)
  do n = 1, self%nf
     call self%fields(n)%field%write_file(ncid_vector)
  enddo

  ! close nc files
  call close_restart_ncid(ncid_lsm, filename_lsm)
  call close_restart_ncid(ncid_hydro, filename_hydro)
end subroutine write_fields_to_file


!> Helper function get get ncids for the two restart files (RESTART and
!> HYDRO_RST. Contains harded coded information RESTART is index 1 and
!> RESTART is index 2. If a file is not returned, the unopened_ncid
!> constant is returned.
function open_get_restart_ncid(self, filename_lsm, filename_hydro, mode) result(ncid)
  class(wrf_hydro_nwm_jedi_fields), target, intent(inout) :: self
  character(len=*), optional, intent(in) :: filename_lsm, filename_hydro
  integer, optional, intent(in) :: mode

  integer :: ncid, mode_local
  logical :: read_file = .false.
  character(len=256) :: filename
  integer :: ierr, n, ncid_index

  if (present(filename_lsm) .and. present(filename_hydro)) then
     write(*,*) "FATAL ERROR: get_restart_ncid: both optional file arguments not allowed."
     stop 124
  else if (.not.(present(filename_lsm)) .and. .not.(present(filename_hydro))) then
     write(*,*) &
          "FATAL ERROR: get_restart_ncid: at least one optional file argument required."
     stop 125
  endif

  ! Hard coded association. A dictionary/hashtable could be used (code in utilties).
  if (present(filename_lsm)) then
     filename = filename_lsm
     ncid_index = 1
  else
     filename = filename_hydro
     ncid_index = 2
  end if

  ! Can this be removed?
  if (.not. present(mode)) then
     mode_local = nf90_nowrite
  else
     mode_local = mode
  end if

  do n = 1, self%nf
     if (self%fields(n)%field%ncid_index .eq. ncid_index) read_file = .true.
  enddo
  write(*, *) filename
  write(*, *) read_file

  if (read_file) then
     ierr = nf90_open(filename, mode_local, ncid)
     call error_handler(ierr, "STOP: file can not be opened: "//trim(filename))
     write(*, *) filename
  else
     ncid = unopened_ncid
  end if
end function open_get_restart_ncid


subroutine close_restart_ncid(ncid, filename)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: filename

  integer :: ierr

  if (ncid .eq. unopened_ncid) return

  ierr = nf90_close(ncid)
  call error_handler(ierr, "STOP: file can not be closed: "//trim(filename))
end subroutine close_restart_ncid


! Get file time ----------------------
function get_lsm_file_time(ncid) result(f_dt)
  integer,        intent(in) :: ncid
  type(datetime)             :: f_dt

  !character*19, dimension(1) :: time_char
  character*19 :: time_char
  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, "Times", varid)
  call error_handler(ierr, "Problem finding variable 'Times' in LSM restart file")
  ierr = nf90_get_var(ncid, varid, time_char)
  call error_handler(ierr, "Problem finding value of 'Times' in LSM restart file")
  ! Create a DateTime from an ISO8601 string, e.g.: 2020-09-02T21:08:08Z
  ! time_char is of the form: 2017-01-01_00:00:00
  write(*,*) "time_char: "//time_char
  write(*,*) "ISO08601 : "//time_char(1:10)//'T'//time_char(12:19)//'Z'
  call datetime_create(time_char(1:10)//'T'//time_char(12:19)//'Z', f_dt)
end function get_lsm_file_time


function get_hydro_file_time(ncid) result(f_dt)
  integer,        intent(in) :: ncid
  type(datetime)             :: f_dt

  character*19, dimension(1) :: time_char
  integer :: ierr
  integer :: varid

  ! Restart_Time = "2017-01-01_00:00:00" ;

  !ierr = nf90_inq_varid(ncid, "Times", varid)
  !call error_handler(ierr, &
  !     "Problem finding variable in restart file '"//trim(name)//"'")
  !ierr = nf90_get_var(ncid, varid, time_char)
  !call error_handler(ierr, &
  !     "Problem finding variable in restart file: '"//trim(name)//"'")
  !write(*,*) "time_char: "//time_char
end function get_hydro_file_time


! read/write 1-d ----------------------
subroutine read_file_1d(self, ncid_vector)
  class(field_1d),       intent(inout) :: self
  integer, dimension(2), intent(in)    :: ncid_vector

  call get_from_restart_1d_float( &
       ncid_vector(self%ncid_index), &
       self%wrf_hydro_nwm_name, &
       self%array)
end subroutine read_file_1d


subroutine write_file_1d(self, ncid_vector)
  class(field_1d),       intent(inout) :: self
  integer, dimension(2), intent(in)    :: ncid_vector

  call write_to_restart_1d_float( &
       ncid_vector(self%ncid_index), &
       self%wrf_hydro_nwm_name, &
       self%array)
end subroutine write_file_1d


subroutine get_from_restart_1d_float(ncid, name, array)
  implicit none
  integer,                            intent(in) :: ncid
  character(len=*),                   intent(in)  :: name
  real, dimension(:),                 intent(out) :: array

  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, name, varid)
  call error_handler(ierr, &
       "Problem finding variable in restart file '"//trim(name)//"'")
  ierr = nf90_get_var(ncid, varid, array)
  call error_handler(ierr, &
       "Problem finding variable in restart file: '"//trim(name)//"'")
end subroutine get_from_restart_1d_float


subroutine write_to_restart_1d_float(ncid, name, array)
  implicit none
  integer,                            intent(in) :: ncid
  character(len=*),                   intent(in)  :: name
  real, dimension(:),                 intent(out) :: array

  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, name, varid)
  call error_handler(ierr, &
       "Problem finding variable in restart file '"//trim(name)//"'")
  ierr = nf90_put_var(ncid, varid, array)
  call error_handler(ierr, &
       "Problem writing variable in restart file: '"//trim(name)//"'")
end subroutine write_to_restart_1d_float


! read/write 2-d ----------------------
subroutine read_file_2d(self, ncid_vector)
  class(field_2d),       intent(inout) :: self
  integer, dimension(2), intent(in)    :: ncid_vector

  call get_from_restart_2d_float( &
       ncid_vector(self%ncid_index), &
       self%wrf_hydro_nwm_name, &
       self%array)
end subroutine read_file_2d


subroutine write_file_2d(self, ncid_vector)
  class(field_2d),       intent(inout) :: self
  integer, dimension(2), intent(in)    :: ncid_vector

  call write_to_restart_2d_float( &
       ncid_vector(self%ncid_index), &
       self%wrf_hydro_nwm_name, &
       self%array)
end subroutine write_file_2d


subroutine get_from_restart_2d_float(ncid, name, array)
  implicit none
  integer,                            intent(in)  :: ncid
  character(len=*),                   intent(in)  :: name
  real, dimension(:, :),              intent(out) :: array

  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, name, varid)
  call error_handler(ierr, &
       "Problem finding variable in restart file '"//trim(name)//"'")
  ierr = nf90_get_var(ncid, varid, array)
  call error_handler(ierr, &
       "Problem finding variable in restart file: '"//trim(name)//"'")
end subroutine get_from_restart_2d_float


subroutine write_to_restart_2d_float(ncid, name, array)
  implicit none
  integer,                            intent(in)  :: ncid
  character(len=*),                   intent(in)  :: name
  real, dimension(:, :),              intent(out) :: array

  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, name, varid)
  call error_handler(ierr, &
       "Problem finding variable in restart file '"//trim(name)//"'")
  ierr = nf90_put_var(ncid, varid, array)
  call error_handler(ierr, &
       "Problem writing variable in restart file: '"//trim(name)//"'")
end subroutine write_to_restart_2d_float


! read/write 3-d ----------------------
subroutine read_file_3d(self, ncid_vector)
  class(field_3d),       intent(inout) :: self
  integer, dimension(2), intent(in)    :: ncid_vector

  call get_from_restart_3d_float( &
       ncid_vector(self%ncid_index), &
       self%wrf_hydro_nwm_name, &
       self%array)
end subroutine read_file_3d


subroutine write_file_3d(self, ncid_vector)
  class(field_3d),       intent(inout) :: self
  integer, dimension(2), intent(in)    :: ncid_vector

  call get_from_restart_3d_float( &
       ncid_vector(self%ncid_index), &
       self%wrf_hydro_nwm_name, &
       self%array)
end subroutine write_file_3d


subroutine get_from_restart_3d_float(ncid, name, array)
  implicit none
  integer,                            intent(in)  :: ncid
  character(len=*),                   intent(in)  :: name
  real, dimension(:, :, :),           intent(out) :: array

  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, name, varid)
  call error_handler(ierr, &
       "Problem finding variable in restart file '"//trim(name)//"'")
  ierr = nf90_get_var(ncid, varid, array)
  call error_handler(ierr, &
       "Problem finding variable in restart file: '"//trim(name)//"'")
end subroutine get_from_restart_3d_float


subroutine write_to_restart_3d_float(ncid, name, array)
  implicit none
  integer,                            intent(in)  :: ncid
  character(len=*),                   intent(in)  :: name
  real, dimension(:, :, :),           intent(out) :: array

  integer :: ierr
  integer :: varid

  ierr = nf90_inq_varid(ncid, name, varid)
  call error_handler(ierr, &
       "Problem finding variable in restart file '"//trim(name)//"'")
  ierr = nf90_put_var(ncid, varid, array)
  call error_handler(ierr, &
       "Problem finding variable in restart file: '"//trim(name)//"'")
end subroutine write_to_restart_3d_float


! --------------------------------------------------------------------------------------------------
! subroutine allocate_field(self,xdim_len,ydim_len,dim3_len,short_name,long_name,&
!                           wrf_hydro_nwm_name,units,tracer,integerfield)

! implicit none
! class(wrf_hydro_nwm_jedi_field), target,  intent(inout) :: self
! integer,                       intent(in)    :: xdim_len,ydim_len,dim3_len
! character(len=*),              intent(in)    :: short_name
! character(len=*),              intent(in)    :: long_name
! character(len=*),              intent(in)    :: wrf_hydro_nwm_name
! character(len=*),              intent(in)    :: units
! logical, optional,             intent(in)    :: tracer
! logical, optional,             intent(in)    :: integerfield

! self%xdim_len = xdim_len
! self%ydim_len = ydim_len
! self%dim3_len = dim3_len

! if(.not.allocated(self%array)) then
!    allocate( self%array(xdim_len, ydim_len, dim3_len) )
!    self%lalloc = .true.
! else
!    call abor1_ftn("Fields.F90.allocate_field: Field already allocated")
! end if

! if(dim3_len == 1) self%is_2D = .true.

! self%short_name   = trim(short_name)
! self%long_name    = trim(long_name)
! self%wrf_hydro_nwm_name = trim(wrf_hydro_nwm_name)
! self%units        = trim(units)

! ! self%staggerloc   = staggerloc

! ! if (present(tracer)) then
! !   self%tracer = tracer
! ! else
! !   self%tracer = .false.
! ! endif

! ! ! Fields that are e.g. types
! ! if (present(integerfield)) then
! !   self%integerfield = integerfield
! ! else
! !   self%integerfield = .false.
! ! endif

! end subroutine allocate_field



! --------------------------------------------------------------------------------------------------

! subroutine field_copy(self, rhs)
!   implicit none
!   ! class(wrf_hydro_nwm_jedi_field), intent(inout) :: self
!   ! type(wrf_hydro_nwm_jedi_field),  intent(in)    :: rhs

! end subroutine field_copy

! --------------------------------------------------------------------------------------------------

! subroutine equals(self,rhs)

! implicit none
! ! class(wrf_hydro_nwm_jedi_field), intent(inout) :: self
! ! type (wrf_hydro_nwm_jedi_field), intent(in)    :: rhs

! ! call self%allocate_field( rhs%isc,rhs%iec,rhs%jsc,rhs%jec,rhs%npz, &
! !                           short_name=rhs%short_name, &
! !                           long_name=rhs%long_name, &
! !                           fv3jedi_name=rhs%fv3jedi_name, &
! !                           units=rhs%units, &
! !                           staggerloc=rhs%staggerloc, &
! !                           tracer = rhs%tracer)

! ! self%array = rhs%array

! end subroutine equals

! --------------------------------------------------------------------------------------------------

! logical function has_field(fields, fv3jedi_name)

! ! type(wrf_hydro_nwm_jedi_field), target,  intent(in)  :: fields(:)
! ! character(len=*),             intent(in)  :: fv3jedi_name

! ! integer :: var

! ! has_field = .false.
! ! do var = 1, size(fields)
! !   if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
! !     has_field = .true.
! !     exit
! !   endif
! ! enddo

! end function has_field

! --------------------------------------------------------------------------------------------------

! subroutine allocate_copy_field_array(fields, fv3jedi_name, field_array)

! ! type(wrf_hydro_nwm_jedi_field),               intent(in)  :: fields(:)
! ! character(len=*),                  intent(in)  :: fv3jedi_name
! ! real(kind=c_float), allocatable, intent(out) :: field_array(:,:,:)

! ! integer :: var
! ! logical :: found

! ! if(allocated(field_array)) deallocate(field_array)

! ! found = .false.
! ! do var = 1,size(fields)
! !   if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
! !     allocate(field_array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,fields(var)%npz))
! !     field_array = fields(var)%array
! !     found = .true.
! !     exit
! !   endif
! ! enddo

! ! if (.not.found) call abor1_ftn("fv3jedi_field.allocate_copy_field_array: field "&
! !                                 //trim(fv3jedi_name)//" not found in fields")

! end subroutine allocate_copy_field_array


! subroutine fields_print(nf, fields, name, f_comm)
  ! implicit none
  ! integer,              intent(in)    :: nf
  ! type(wrf_hydro_nwm_jedi_field),  intent(in)    :: fields(nf)
  ! character(len=*),     intent(in)    :: name
  ! type(fckit_mpi_comm), intent(in)    :: f_comm

  ! integer :: i

  ! do i = 1, nf
  !    write(*,*) fields(i)%wrf_hydro_nwm_name
  !    write(*,*) "From Fields.F90: first column:", fields(i)%array(:,1,1)
  !    write(*,*) "------"
  ! end do

! integer :: var
! real(kind=kind_real) :: tmp(3), pstat(3), gs3, gs3g
! character(len=34) :: printname

! integer :: ngrid, sngrid

! ngrid = (fields(1)%iec-fields(1)%isc+1)*(fields(1)%iec-fields(1)%isc+1)
! call f_comm%allreduce(ngrid,sngrid,fckit_mpi_sum())
! sngrid = nint(sqrt(real(sngrid,kind_real)/6.0_kind_real))

! printname = "|     "//trim(name)//" print"

! if (f_comm%rank() == 0) then
!   write(*,"(A70)") "----------------------------------------------------------------------"
!   write(*,"(A34)") printname
!   write(*,"(A70)") "----------------------------------------------------------------------"
!   write(*,"(A70)") " "
!   write(*,"(A27,I5)") "    Cube sphere face size: ", sngrid
!   write(*,"(A70)") " "
! endif

! do var = 1,nf

!   gs3 = real((fields(var)%iec-fields(var)%isc+1)*(fields(var)%jec-fields(var)%jsc+1)*fields(var)%npz, kind_real)
!   call f_comm%allreduce(gs3,gs3g,fckit_mpi_sum())

!   tmp(1) = minval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
!   tmp(2) = maxval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
!   tmp(3) =    sum(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz)**2)

!   call f_comm%allreduce(tmp(1),pstat(1),fckit_mpi_min())
!   call f_comm%allreduce(tmp(2),pstat(2),fckit_mpi_max())
!   call f_comm%allreduce(tmp(3),pstat(3),fckit_mpi_sum())
!   pstat(3) = sqrt(pstat(3)/gs3g)

!   if (f_comm%rank() == 0) write(*,"(A10,A6,ES14.7,A6,ES14.7,A6,ES14.7)") &
!                                    trim(fields(var)%short_name),&
!                                    "| Min=",real(pstat(1),4),&
!                                    ", Max=",real(pstat(2),4),&
!                                    ", RMS=",real(pstat(3),4)

! enddo

! if (f_comm%rank() == 0) &
!   write(*,"(A70)") "---------------------------------------------------------------------"
!end subroutine fields_print

! --------------------------------------------------------------------------------------------------



! subroutine flip_array_vertical(nf,fields)
! implicit none
! integer,             intent(in)    :: nf
! type(wrf_hydro_nwm_jedi_field), intent(inout) :: fields(nf)

! integer :: n, lev_in, lev_out
! real(kind=kind_real), allocatable :: array_tmp(:,:,:)

! do n = 1, nf

!   if (fields(n)%npz > 1) then

!     if (fields(n)%staggerloc == center) then
!       allocate(array_tmp(fields(n)%isc:fields(n)%iec,fields(n)%jsc:fields(n)%jec,1:fields(n)%npz))
!     elseif (fields(n)%staggerloc == north) then
!       allocate(array_tmp(fields(n)%isc:fields(n)%iec,fields(n)%jsc:fields(n)%jec+1,1:fields(n)%npz))
!     elseif (fields(n)%staggerloc == east) then
!       allocate(array_tmp(fields(n)%isc:fields(n)%iec+1,fields(n)%jsc:fields(n)%jec,1:fields(n)%npz))
!     endif

!     do lev_in = 1,fields(n)%npz

!       lev_out = fields(n)%npz-lev_in+1
!       array_tmp(:,:,lev_out) = fields(n)%array(:,:,lev_in)

!     enddo

!     fields(n)%array = array_tmp

!     deallocate(array_tmp)

!   endif

! enddo
!end subroutine flip_array_vertical


! subroutine copy_subset(rhs,lhs,not_copied)
! implicit none
! type(wrf_hydro_nwm_jedi_field),        intent(in)    :: rhs(:)
! type(wrf_hydro_nwm_jedi_field),        intent(inout) :: lhs(:)
! character(len=10), allocatable, optional, intent(out)   :: not_copied(:)

! integer :: var
! character(len=10) :: not_copied_(10000)
! integer :: num_not_copied

! ! Loop over fields and copy if existing in both
! num_not_copied = 0
! do var = 1, size(lhs)
!   if (has_field(rhs, lhs(var)%fv3jedi_name )) then
!     call copy_field_array(rhs, lhs(var)%fv3jedi_name, lhs(var)%array)
!   else
!     num_not_copied = num_not_copied + 1
!     not_copied_(num_not_copied) = lhs(var)%fv3jedi_name
!   endif
! enddo

! ! Send back list of variables not retrivable from rhs
! if (present(not_copied) .and. num_not_copied > 0) then
!   allocate(not_copied(num_not_copied))
!   not_copied(1:num_not_copied) = not_copied_(1:num_not_copied)
! endif
!end subroutine copy_subset

!-----------------------------------------------------------------------------

! subroutine fields_gpnorm(nf, fields, pstat, f_comm)

! ! implicit none
! ! integer,              intent(in)    :: nf
! ! type(wrf_hydro_nwm_jedi_field),  intent(in)    :: fields(nf)
! ! real(kind=c_float), intent(inout) :: pstat(3, nf)
! ! type(fckit_mpi_comm), intent(in)    :: f_comm

! ! integer :: var
! ! real(kind=kind_real) :: tmp(3),  gs3, gs3g

! ! do var = 1,nf

! !   gs3 = real((fields(var)%iec-fields(var)%isc+1)*(fields(var)%jec-fields(var)%jsc+1)*fields(var)%npz, kind_real)
! !   call f_comm%allreduce(gs3,gs3g,fckit_mpi_sum())

! !   tmp(1) = minval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
! !   tmp(2) = maxval(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz))
! !   tmp(3) =    sum(fields(var)%array(fields(var)%isc:fields(var)%iec,fields(var)%jsc:fields(var)%jec,1:fields(var)%npz)**2)

! !   call f_comm%allreduce(tmp(1),pstat(1,var),fckit_mpi_min())
! !   call f_comm%allreduce(tmp(2),pstat(2,var),fckit_mpi_max())
! !   call f_comm%allreduce(tmp(3),pstat(3,var),fckit_mpi_sum())
! !   pstat(3,var) = sqrt(pstat(3,var)/gs3g)

! ! enddo
! end subroutine fields_gpnorm


! subroutine long_name_to_wrf_hydro_name(fields, long_name, wrf_hydro_nwm_name)

!   type(elem_field), intent(in)  :: fields(:)
!   character(len=*),    intent(in)  :: long_name
!   character(len=*),    intent(out) :: wrf_hydro_nwm_name

!   integer :: n

  ! select case (long_name)
  ! case ("swe")
  !    wrf_hydro_nwm_name = "SNEQV"
  ! case ("snow_depth")
  !    wrf_hydro_nwm_name = "SNOWH"
  ! case ("leaf_area")
  !    wrf_hydro_nwm_name = "LAI"
  ! case default
  !    wrf_hydro_nwm_name = "null"
  ! end select

! do n = 1, size(fields)
!   if (trim(long_name) == trim(fields(n)%long_name)) then
!     wrf_hydro_nwm_name = trim(fields(n)%wrf_hydro_nwm_name)
!     return
!   endif
! enddo

! ! Try with increment_of_ prepended to long_name
! do n = 1, size(fields)
!    if ("increment_of_"//trim(long_name) == trim(fields(n)%long_name)) then
!       wrf_hydro_nwm_name = trim(fields(n)%wrf_hydro_nwm_name)
!       return
!    endif
! enddo

! call abor1_ftn("wrf_hydro_nwm_jedi_fields_mod.long_name_to_wrf_hydro_nwm_jedi_name long_name "//trim(long_name)//&
!      " not found in fields.")
! end subroutine long_name_to_wrf_hydro_name


! subroutine copy_field_array(fields, fv3jedi_name, field_array)
! ! type(wrf_hydro_nwm_jedi_field),  intent(in)  :: fields(:)
! ! character(len=*),     intent(in)  :: fv3jedi_name
! ! real(kind=c_float), intent(out) :: field_array(:,:,:)

! ! integer :: var
! ! logical :: found

! ! found = .false.
! ! do var = 1,size(fields)
! !   if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
! !     field_array = fields(var)%array
! !     found = .true.
! !     exit
! !   endif
! ! enddo

! ! if (.not.found) call abor1_ftn("fv3jedi_field.copy_field_array: field "&
! !                                 //trim(fv3jedi_name)//" not found in fields")
! end subroutine copy_field_array


! subroutine pointer_field(fields, wrf_hydro_nwm_name, field_pointer)
! ! type(wrf_hydro_nwm_jedi_field), target,  intent(in)  :: fields(:)
! ! character(len=*),             intent(in)  :: wrf_hydro_nwm_name
! ! type(wrf_hydro_nwm_jedi_field), pointer, intent(out) :: field_pointer

! ! integer :: var
! ! logical :: found

! ! if(associated(field_pointer)) nullify(field_pointer)

! ! found = .false.
! ! do var = 1,size(fields)
! !   if ( trim(fields(var)%wrf_hydro_nwm_name) == trim(wrf_hydro_nwm_name)) then
! !     field_pointer => fields(var)
! !     found = .true.
! !     exit
! !   endif
! ! enddo

! ! if (.not.found) call abor1_ftn("wrf_hydro_field.pointer_field: field "&
! !                                 //trim(wrf_hydro_nwm_name)//" not found in fields")
! end subroutine pointer_field


! subroutine pointer_field_array(fields, wrf_hydro_nwm_jedi_name, array_pointer)
! ! type(wrf_hydro_nwm_jedi_field), target,   intent(in)    :: fields(:)
! ! character(len=*),              intent(in)    :: wrf_hydro_nwm_jedi_name
! ! real(kind=c_float), pointer, intent(out)   :: array_pointer(:,:,:)

! ! integer :: var
! ! logical :: found

! ! if(associated(array_pointer)) nullify(array_pointer)

! ! found = .false.
! ! do var = 1,size(fields)
! !   if ( trim(fields(var)%fv3jedi_name) == trim(fv3jedi_name)) then
! !     array_pointer => fields(var)%array
! !     found = .true.
! !     exit
! !   endif
! ! enddo

! ! if (.not.found) call abor1_ftn("fv3jedi_field.pointer_field_array: field "&
! !                                 //trim(fv3jedi_name)//" not found in fields")
! end subroutine pointer_field_array

! -----------------------------------------------------------------------------
! set_atlas

subroutine set_atlas(self, geom, vars, afieldset)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(inout) :: afieldset
  integer :: jvar, ff
  logical :: found

  do jvar=1,vars%nvars()
     found = .false.
     do ff=1,self%nf
        if (trim(self%fields(ff)%field%short_name)==trim(vars%variable(jvar))) then
           call self%fields(ff)%field%set_atlas(geom, afieldset)
           found = .true.
        end if
    end do
    if (.not.found) call abor1_ftn("set_atlas: field "//trim(vars%variable(jvar))//" not found")
  enddo
end subroutine set_atlas

subroutine set_atlas_1d(self, geom, afieldset)
  implicit none
  class(field_1d), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(atlas_fieldset), intent(inout) :: afieldset

  call abor1_ftn("set_atlas_1d: no set_atlas interface for 1d fields")
end subroutine set_atlas_1d

subroutine set_atlas_2d(self, geom, afieldset)
  implicit none
  class(field_2d), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(atlas_fieldset), intent(inout) :: afieldset

  type(atlas_field) :: afield

  if (.not.afieldset%has_field(self%long_name)) then
    afield = geom%lsm%afunctionspace%create_field(name=self%long_name,kind=atlas_real(c_double),levels=0)
    call afieldset%add(afield)
    call afield%final()
  end if
end subroutine set_atlas_2d

subroutine set_atlas_3d(self, geom, afieldset)
  implicit none
  class(field_3d), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(atlas_fieldset), intent(inout) :: afieldset

  type(atlas_field) :: afield

  if (.not.afieldset%has_field(self%long_name)) then
    afield = geom%lsm%afunctionspace%create_field(name=self%long_name,kind=atlas_real(c_double),levels=self%zdim_len)
    call afieldset%add(afield)
    call afield%final()
  end if
end subroutine set_atlas_3d

! -----------------------------------------------------------------------------
! to_atlas

subroutine to_atlas(self, geom, vars, afieldset)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(inout) :: afieldset
  integer :: jvar, ff
  logical :: found

  do jvar=1,vars%nvars()
     found = .false.
     do ff=1,self%nf
        if (trim(self%fields(ff)%field%short_name)==trim(vars%variable(jvar))) then
           call self%fields(ff)%field%to_atlas(geom, afieldset)
           found = .true.
        end if
    end do
    if (.not.found) call abor1_ftn("to_atlas: field "//trim(vars%variable(jvar))//" not found")
  enddo
end subroutine to_atlas

subroutine to_atlas_1d(self, geom, afieldset)
  implicit none
  class(field_1d), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(atlas_fieldset), intent(inout) :: afieldset

  call abor1_ftn("to_atlas_1d: no to_atlas interface for 1d fields")

end subroutine to_atlas_1d

subroutine to_atlas_2d(self, geom, afieldset)
  implicit none
  class(field_2d), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(atlas_fieldset), intent(inout) :: afieldset

  integer :: ix, iy, inode
  real(c_double), pointer :: ptr(:)
  type(atlas_field) :: afield

  if (afieldset%has_field(self%long_name)) then
    afield = afieldset%field(trim(self%long_name))
  else
    afield = geom%lsm%afunctionspace%create_field(name=self%long_name,kind=atlas_real(c_double),levels=0)
    call afieldset%add(afield)
  end if
  call afield%data(ptr)
  inode = 0
  do iy=1, self%ydim_len
    do ix=1, self%xdim_len
       inode = inode+1
       ptr(inode) = real(self%array(ix, iy), c_double)
    enddo
  enddo
  call afield%final()
end subroutine to_atlas_2d

subroutine to_atlas_3d(self, geom, afieldset)
  implicit none
  class(field_3d), intent(in) :: self
  type(wrf_hydro_nwm_jedi_geometry), intent(in) :: geom
  type(atlas_fieldset), intent(inout) :: afieldset

  integer :: ix, iy, iz, inode
  real(c_double), pointer :: ptr(:,:)
  type(atlas_field) :: afield

  if (afieldset%has_field(self%long_name)) then
    afield = afieldset%field(trim(self%long_name))
  else
    afield = geom%lsm%afunctionspace%create_field(name=self%long_name,kind=atlas_real(c_double),levels=self%zdim_len)
    call afieldset%add(afield)
  end if
  call afield%data(ptr)
  do iz=1, self%zdim_len
    inode = 0
    do iy=1, self%ydim_len
      do ix=1, self%xdim_len
         inode = inode+1
         ptr(iz, inode) = real(self%array(ix, iy, iz), c_double)
      enddo
    enddo
  enddo
  call afield%final()
end subroutine to_atlas_3d

! -----------------------------------------------------------------------------
! from_atlas

subroutine from_atlas(self, vars, afieldset)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(inout) :: self
  type(oops_variables), intent(in) :: vars
  type(atlas_fieldset), intent(in) :: afieldset
  integer :: jvar, ff
  logical :: found

  do jvar=1,vars%nvars()
     found = .false.
     do ff=1,self%nf
        if (trim(self%fields(ff)%field%short_name)==trim(vars%variable(jvar))) then
           call self%fields(ff)%field%from_atlas(afieldset)
           found = .true.
        end if
    end do
    if (.not.found) call abor1_ftn("from_atlas: field "//trim(vars%variable(jvar))//" not found")
  enddo
end subroutine from_atlas

subroutine from_atlas_1d(self, afieldset)
  implicit none
  class(field_1d), intent(inout) :: self
  type(atlas_fieldset), intent(in) :: afieldset

  call abor1_ftn("from_atlas_1d: no from_atlas interface for 1d fields")

end subroutine from_atlas_1d

subroutine from_atlas_2d(self, afieldset)
  implicit none
  class(field_2d), intent(inout) :: self
  type(atlas_fieldset), intent(in) :: afieldset

  integer :: ix, iy, inode
  real(c_double), pointer :: ptr(:)
  type(atlas_field) :: afield

  afield = afieldset%field(trim(self%long_name))
  call afield%data(ptr)
  inode = 0
  do iy=1, self%ydim_len
    do ix=1, self%xdim_len
       inode = inode+1
       self%array(ix, iy) = real(ptr(inode), c_float)
    enddo
  enddo
  call afield%final()
end subroutine from_atlas_2d

subroutine from_atlas_3d(self, afieldset)
  implicit none
  class(field_3d), intent(inout) :: self
  type(atlas_fieldset), intent(in) :: afieldset

  integer :: ix, iy, iz, inode
  real(c_double), pointer :: ptr(:,:)
  type(atlas_field) :: afield

  afield = afieldset%field(trim(self%long_name))
  call afield%data(ptr)
  do iz=1, self%zdim_len
    inode = 0
    do iy=1, self%ydim_len
      do ix=1, self%xdim_len
         inode = inode+1
         self%array(ix, iy, iz) = real(ptr(iz, inode), c_float)
      enddo
    enddo
  enddo
  call afield%final()
end subroutine from_atlas_3d

! -----------------------------------------------------------------------------
! get_point

subroutine get_point(self, geoiter, values_len, values)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(in) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(inout) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter
  integer :: ff
  do ff=1,self%nf
       call self%fields(ff)%field%get_point(geoiter, values_len, values)
  end do
end subroutine get_point

subroutine get_point_1d(self, geoiter, values_len, values)
  implicit none
  class(field_1d), intent(in) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(inout) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter

  call abor1_ftn("get_point_1d: no get_point interface for 1d fields")
end subroutine get_point_1d

subroutine get_point_2d(self, geoiter, values_len, values)
  implicit none
  class(field_2d), intent(in) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(inout) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter

  values = self%array(geoiter%iind, geoiter%jind)
end subroutine get_point_2d

subroutine get_point_3d(self, geoiter, values_len, values)
  implicit none
  class(field_3d), intent(in) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(inout) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter

  ! call abor1_ftn("get_point_3d: no get_point interface for 3d fields")
  values = self%array(geoiter%iind, geoiter%jind, :) 
end subroutine get_point_3d

! -----------------------------------------------------------------------------
! set_point

subroutine set_point(self, geoiter, values_len, values)
  implicit none
  class(wrf_hydro_nwm_jedi_fields), intent(inout) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(in) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter
  integer :: ff
  do ff=1,self%nf
       call self%fields(ff)%field%set_point(geoiter, values_len, values)
  end do
end subroutine set_point

subroutine set_point_1d(self, geoiter, values_len, values)
  implicit none
  class(field_1d), intent(inout) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(in) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter

  call abor1_ftn("set_point_1d: no set_point interface for 1d fields")
end subroutine set_point_1d

subroutine set_point_2d(self, geoiter, values_len, values)
  implicit none
  class(field_2d), intent(inout) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(in) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter

  self%array(geoiter%iind, geoiter%jind) = values(1)
end subroutine set_point_2d

subroutine set_point_3d(self, geoiter, values_len, values)
  implicit none
  class(field_3d), intent(inout) :: self
  integer(c_int),   intent(in) :: values_len
  real(kind=c_double),  intent(in) :: values(values_len)
  type(wrf_hydro_nwm_jedi_geometry_iter),  intent(in) :: geoiter

  ! call abor1_ftn("set_point_3d: no set_point interface for 3d fields")
  self%array(geoiter%iind, geoiter%jind, :) = values
end subroutine set_point_3d

end module wrf_hydro_nwm_jedi_fields_mod