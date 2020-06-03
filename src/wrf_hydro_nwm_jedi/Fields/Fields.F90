! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_field_mod

  use iso_c_binding, only: c_int, c_float,c_double
  use fckit_mpi_module
  use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry, error_handler, indices
  use oops_variables_mod
  use netcdf
!use fv3jedi_kinds_mod, only: kind_real
!use mpp_domains_mod,   only: east, north, center

  implicit none

  public

  type, abstract :: base_field
     character(len=32) :: short_name = "null"   !Short name (to match file name)
     character(len=10) :: wrf_hydro_nwm_name = "null" !Common name
     character(len=64) :: long_name = "null"    !More descriptive name
     character(len=32) :: units = "null"        !Units for the field
     logical :: tracer = .false.
   contains
     procedure (print_field_interface), pass(self), deferred :: print_field
     procedure (read_file_interface), pass(self), deferred :: read_file
     procedure (get_value_field_interface), pass(self), deferred :: get_value
  end type base_field

  abstract interface
     subroutine print_field_interface(self)
       import base_field
       class(base_field), intent(in) :: self
     end subroutine print_field_interface

     subroutine read_file_interface(self,filename,xstart,xend,ystart,yend)
       import base_field
       class(base_field), intent(inout) :: self
       character(len=*), intent(in) :: filename
       integer, intent(in) :: xstart, xend, ystart, yend
     end subroutine read_file_interface

     function get_value_field_interface(self,ind) result(val)
       use iso_c_binding, only : c_float
       import base_field, indices
       class(base_field), intent(in) :: self
       type(indices), intent(in) :: ind
       real(kind=c_float) :: val
     end function get_value_field_interface
  end interface

  type, extends(base_field) :: field_2d
     integer :: dim1_len, dim2_len
     real(c_float), allocatable :: array(:,:)
   contains
     procedure, pass(self) :: print_field => print_field_2d
     procedure, pass(self) :: read_file => read_file_2d
     procedure, pass(self) :: get_value => get_value_2d
     procedure, pass(self) :: fill => fill_field_2d
     !procedure, pass(self) :: mean_stddev
     procedure, pass(self) :: rms => field_rms_2d
     !    Destructor
     ! final :: destroy_field_2d
  end type field_2d

  type, extends(base_field) :: field_3d
     integer :: dim1_len, dim2_len, dim3_len
     real(c_float), allocatable :: array(:,:,:)
   contains
     procedure, pass(self) :: print_field => print_field_3d
     procedure, pass(self) :: read_file => read_file_3d
     procedure, pass(self) :: get_value => get_value_3d
     procedure, pass(self) :: fill => fill_field_3d
     procedure, pass(self) :: rms => field_rms_3d
     !    Destructor
     ! final :: destroy_field_3d
  end type field_3d

  !The function of this element is twofold:
  ! 1) It allows one to create an array of base class object
  ! 2) It makes it possible to organize the element data structures
  ! other than simple arrays (e.g. binary search tree, hash tables, etc...)
  type elem_field
     class(base_field), allocatable :: field
  end type elem_field
  
public :: wrf_hydro_nwm_jedi_fields, checksame!, &
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
!     mean_stddev

!Field type mirror of C class
type :: wrf_hydro_nwm_jedi_fields
   integer :: nf
 ! integer :: dim1_len, dim2_len, dim3_len !refer to geometry and depth
 ! real(kind=c_float), allocatable :: array(:,:,:)
   ! logical :: integerfield = .false.
   type(elem_field), allocatable :: fields(:)
contains
  procedure :: create
  procedure :: search_field
  procedure :: read_fields_from_file
  procedure :: print_all_fields
  ! procedure :: allocate_field
  ! procedure :: equals
  ! procedure :: copy => field_copy
  ! generic :: assignment(=) => equals
  ! procedure :: deallocate_field
  ! procedure :: mean_stddev
end type wrf_hydro_nwm_jedi_fields

! --------------------------------------------------------------------------------------------------

contains

  subroutine create(self, geom, vars)
    implicit none
    class(wrf_hydro_nwm_jedi_fields),  intent(inout) :: self
    type(wrf_hydro_nwm_jedi_geometry),   intent(in)    :: geom
    type(oops_variables), intent(in)    :: vars

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
       case("SNEQV")
          vcount=vcount+1;

          allocate(tmp_2d_field)
          call tmp_2d_field%fill(geom%dim1_len, geom%dim2_len, &
               short_name = vars%variable(var), &
               long_name = 'snow_water_equivalent', &
               wrf_hydro_nwm_name = 'SNEQV', units = 'mm')

          allocate(self%fields(vcount)%field, source = tmp_2d_field)
          deallocate(tmp_2d_field)

       case("SNOWH")
          vcount=vcount+1;

          allocate(tmp_2d_field)
          call tmp_2d_field%fill(geom%dim1_len, geom%dim2_len, &
               short_name = vars%variable(var), &
               long_name = 'snow_depth', &
               wrf_hydro_nwm_name = 'SNOWH', units = 'm')
          
          allocate(self%fields(vcount)%field, source = tmp_2d_field)
          deallocate(tmp_2d_field)

       case("LAI")
          vcount=vcount+1;
          
          call tmp_2d_field%fill(geom%dim1_len, geom%dim2_len, &
               short_name = vars%variable(var),&
               long_name = 'leaf_area', &
               wrf_hydro_nwm_name = 'LAI', units = 'm^2m^-2')

          allocate(self%fields(vcount)%field, source = tmp_2d_field)
          deallocate(tmp_2d_field)

       case("SNLIQ")
          vcount=vcount+1;
          !MIDDLE DIMENSION HARDCODED
          call tmp_3d_field%fill(geom%dim1_len, 3, geom%dim2_len, &
               short_name = vars%variable(var),&
               long_name = 'snow_liquid', &
               wrf_hydro_nwm_name = 'SNLIQ', units = 'liter') !unit invented

          allocate(self%fields(vcount)%field, source = tmp_3d_field)
          deallocate(tmp_3d_field)

       case default
          call abor1_ftn("Create: unknown variable "//trim(vars%variable(var)))          
       end select
    end do

  end subroutine create

  subroutine fill_field_2d(self,dim1_len,dim2_len,short_name,long_name,&
       wrf_hydro_nwm_name,units,tracer)
    implicit none

    class(field_2d),intent(inout) :: self
    integer, intent(in)    :: dim1_len, dim2_len
    character(len=*),              intent(in)    :: short_name
    character(len=*),              intent(in)    :: long_name
    character(len=*),              intent(in)    :: wrf_hydro_nwm_name
    character(len=*),              intent(in)    :: units
    logical, optional,             intent(in)    :: tracer

    self%dim1_len = dim1_len
    self%dim2_len = dim2_len

    if(.not.allocated(self%array)) then
       allocate( self%array(dim1_len, dim2_len) )
    else
       call abor1_ftn("Fields.F90.allocate_field: Field already allocated")
    end if

    self%short_name   = trim(short_name)
    self%long_name    = trim(long_name)
    self%wrf_hydro_nwm_name = trim(wrf_hydro_nwm_name)
    self%units        = trim(units)
    
  end subroutine fill_field_2d

  subroutine fill_field_3d(self,dim1_len,dim2_len,dim3_len,short_name,long_name,&
       wrf_hydro_nwm_name,units,tracer)
    implicit none
    
    class(field_3d),intent(inout) :: self
    integer, intent(in)    :: dim1_len, dim2_len, dim3_len
    character(len=*),              intent(in)    :: short_name
    character(len=*),              intent(in)    :: long_name
    character(len=*),              intent(in)    :: wrf_hydro_nwm_name
    character(len=*),              intent(in)    :: units
    logical, optional,             intent(in)    :: tracer

    self%dim1_len = dim1_len
    self%dim2_len = dim2_len
    self%dim3_len = dim3_len
    
    if(.not.allocated(self%array)) then
       allocate( self%array(dim1_len, dim2_len, dim3_len) )
    else
       call abor1_ftn("Fields.F90.allocate_field: Field already allocated")
    end if

    self%short_name   = trim(short_name)
    self%long_name    = trim(long_name)
    self%wrf_hydro_nwm_name = trim(wrf_hydro_nwm_name)
    self%units        = trim(units)
    
  end subroutine fill_field_3d

  subroutine read_file_2d(self,filename,xstart,xend,ystart,yend)
    class(field_2d),intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: xstart, xend, ystart, yend
    integer :: ixfull, jxfull

    ixfull = xend-xstart+1
    jxfull = yend-ystart+1

    call get_from_restart_2d_float(filename, xstart, xend, xstart, ixfull, jxfull, self%wrf_hydro_nwm_name, self%array)
    
  end subroutine read_file_2d

  subroutine read_file_3d(self,filename,xstart,xend,ystart,yend)
    class(field_3d),intent(inout) :: self
    character(len=*), intent(in) :: filename
    integer, intent(in) :: xstart, xend, ystart, yend
    integer :: ixfull, jxfull

    ixfull = xend-xstart+1
    jxfull = yend-ystart+1

    call get_from_restart_3d(filename, xstart, xend, xstart, ixfull, jxfull, self%wrf_hydro_nwm_name, self%array)
    
  end subroutine read_file_3d

! --------------------------------------------------------------------------------------------------
  
  function get_value_2d(self,ind) result(val)
    class(field_2d), intent(in) :: self
    type(indices), intent(in) :: ind
    real(c_float) :: val

    val = self%array(ind%ind_1,ind%ind_2)
    
  end function get_value_2d

  ! --------------------------------------------------------------------------------------------------

  function get_value_3d(self,ind) result(val)
    class(field_3d), intent(in) :: self
    type(indices), intent(in) :: ind
    real(c_float) :: val

    val = self%array(ind%ind_1,ind%ind_2,ind%ind_3)
    
  end function get_value_3d

  ! --------------------------------------------------------------------------------------------------
  
! subroutine allocate_field(self,dim1_len,dim2_len,dim3_len,short_name,long_name,&
!                           wrf_hydro_nwm_name,units,tracer,integerfield)

! implicit none
! class(wrf_hydro_nwm_jedi_field), target,  intent(inout) :: self
! integer,                       intent(in)    :: dim1_len,dim2_len,dim3_len
! character(len=*),              intent(in)    :: short_name
! character(len=*),              intent(in)    :: long_name
! character(len=*),              intent(in)    :: wrf_hydro_nwm_name
! character(len=*),              intent(in)    :: units
! logical, optional,             intent(in)    :: tracer
! logical, optional,             intent(in)    :: integerfield

! self%dim1_len = dim1_len
! self%dim2_len = dim2_len
! self%dim3_len = dim3_len

! if(.not.allocated(self%array)) then
!    allocate( self%array(dim1_len, dim2_len, dim3_len) )
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

! subroutine deallocate_field(self)

! ! implicit none
! ! class(wrf_hydro_nwm_jedi_field), intent(inout) :: self

! ! if(self%lalloc) deallocate(self%array)
! ! self%lalloc = .false.

! end subroutine deallocate_field

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

  ! --------------------------------------------------------------------------------------------------
  ! This subroutine unifies the long_name_to_wrf_hydro_name and pointer_field routines
  subroutine search_field(self,long_name,field_pointer)
    class(wrf_hydro_nwm_jedi_fields),  target, intent(inout) :: self
    character(len=*),    intent(in)  :: long_name
    class(base_field), pointer, intent(out) :: field_pointer

    integer :: n
    character(len=255) :: wrf_hydro_nwm_name
    
    field_pointer => null()

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
    
    !linear search
    do n = 1, size(self%fields)
       if (trim(wrf_hydro_nwm_name) == trim(self%fields(n)%field%wrf_hydro_nwm_name)) then
          field_pointer => self%fields(n)%field
          return
       endif
    enddo
    
    call abor1_ftn("wrf_hydro_nwm_jedi_field_mod.long_name_to_wrf_hydro_nwm_jedi_name long_name "//trim(long_name)//&
         " not found in fields.")
    
  end subroutine search_field

  ! --------------------------------------------------------------------------------------------------

  subroutine read_fields_from_file(self,filename,xstart,xend,ystart,yend)
    class(wrf_hydro_nwm_jedi_fields),  target, intent(inout) :: self
    character(len=*),    intent(in)  :: filename
    integer, intent(in) :: xstart, xend, ystart, yend
    
    integer :: n

    do n = 1, size(self%fields)
       call self%fields(n)%field%read_file(filename,xstart,xend,ystart,yend)
    enddo
    
  end subroutine read_fields_from_file

  ! --------------------------------------------------------------------------------------------------
  
  subroutine print_all_fields(self)
    implicit none
    class(wrf_hydro_nwm_jedi_fields),  intent(in) :: self
    integer :: f

    do f = 1,size(self%fields)
       call self%fields(f)%field%print_field()
    end do
    
  end subroutine print_all_fields

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

! call abor1_ftn("wrf_hydro_nwm_jedi_field_mod.long_name_to_wrf_hydro_nwm_jedi_name long_name "//trim(long_name)//&
!      " not found in fields.")

! end subroutine long_name_to_wrf_hydro_name

! --------------------------------------------------------------------------------------------------

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

! --------------------------------------------------------------------------------------------------

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

! --------------------------------------------------------------------------------------------------

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

! --------------------------------------------------------------------------------------------------

! subroutine mean_stddev(self,mean,stddev,rms)
!   implicit none
!   ! class(wrf_hydro_nwm_jedi_field), target,  intent(in) :: self
!   ! real(c_float),intent(inout) :: mean,stddev,rms
!   ! real(c_double) :: tmp
!   ! integer :: n, i, j

!   ! n = size(self%array,1)*size(self%array,2)

!   ! tmp = 0.d0
!   ! tmp = sum(self%array(:,:,1))
!   ! tmp = tmp/n

!   ! mean = real(tmp,kind=c_float)

!   ! tmp = 0.d0
  
!   ! !Computing stddev
!   ! do i=1,size(self%array,1)
!   !    do j=1,size(self%array,2)
!   !       tmp = tmp + (self%array(i,j,1) - mean)**2
!   !    end do
!   ! end do
!   ! tmp = tmp / n
!   ! stddev = real(tmp,kind=c_float)

!   ! tmp = 0.d0

!   ! do i = 1,size(self%array,1)
!   !    do j = 1,size(self%array,2)
!   !       tmp = tmp + self%array(i,j,1)**2
!   !    enddo
!   ! enddo

!   ! tmp = sqrt(tmp/real(n,kind=c_double))

!   ! rms = real(tmp,kind=c_float)
  
! end subroutine mean_stddev

! --------------------------------------------------------------------------------------------------

function field_rms_2d(self) result(rms)
  implicit none
  class(field_2d),  intent(in) :: self
  real(kind=c_float) :: rms
  !type(fckit_mpi_comm), intent(in)    :: f_comm
  
  integer :: i, j, ii
  real(kind=c_float) :: zz
  
  zz = 0.0
  ii = 0
  
  do j = 1,self%dim2_len
     do i = 1,self%dim1_len
        zz = zz + self%array(i,j)**2
        ii = ii + 1 !unnecessary
     enddo
  enddo

! !Get global values
! call f_comm%allreduce(zz,rms,fckit_mpi_sum())
! call f_comm%allreduce(ii,iisum,fckit_mpi_sum())

  rms = sqrt(zz/real(ii,c_float))

end function field_rms_2d

function field_rms_3d(self) result(rms)
  implicit none
  class(field_3d),  intent(in) :: self
  real(kind=c_float) :: rms
  !type(fckit_mpi_comm), intent(in)    :: f_comm
  
  integer :: i, j, k, ii
  real(kind=c_float) :: zz
  
  zz = 0.0
  ii = 0

  do k = 1,self%dim2_len
     do j = 1, 3!HARDCODED
        do i = 1,self%dim1_len
           zz = zz + self%array(i,j,k)**2
           ii = ii + 1 !unnecessary
        enddo
     enddo
  end do
! !Get global values
! call f_comm%allreduce(zz,rms,fckit_mpi_sum())
! call f_comm%allreduce(ii,iisum,fckit_mpi_sum())

  rms = sqrt(zz/real(ii,c_float))

end function field_rms_3d

! --------------------------------------------------------------------------------------------------

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

! --------------------------------------------------------------------------------------------------

subroutine print_field_2d(self)
  implicit none
  class(field_2d),intent(in) :: self

  write(*,*) 'Printing ', self%long_name
  
end subroutine print_field_2d

subroutine print_field_3d(self)
  implicit none
  class(field_3d),intent(in) :: self

  write(*,"(A70)") 'Printing ', self%long_name
  
end subroutine print_field_3d

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

subroutine checksame(self,other,method)

implicit none
type(wrf_hydro_nwm_jedi_fields), intent(in) :: self
type(wrf_hydro_nwm_jedi_fields), intent(in) :: other
character(len=*),    intent(in) :: method

integer :: var

if (size(self%fields) .ne. size(other%fields)) then
  call abor1_ftn(trim(method)//"(checksame): Different number of fields")
endif

do var = 1,size(self%fields)
  if (self%fields(var)%field%wrf_hydro_nwm_name .ne. other%fields(var)%field%wrf_hydro_nwm_name) then
     call abor1_ftn(trim(method)//"(checksame): field "//trim(self%fields(var)%field%wrf_hydro_nwm_name)//&
                     " not in the equivalent position in the right hand side")
  endif
enddo

end subroutine checksame

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

! ------------------------------------------------------------------------------

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


subroutine get_from_restart_3d(restart_filename_remember, parallel_xstart, parallel_xend, subwindow_xstart, ixfull, jxfull, name, array, return_error)
  implicit none
    character(len=*), intent(in) :: restart_filename_remember
    integer,                            intent(in) :: parallel_xstart
    integer,                            intent(in) :: parallel_xend
    integer,                            intent(in) :: subwindow_xstart
    integer,                            intent(in) :: ixfull
    integer,                            intent(in) :: jxfull
    character(len=*),                   intent(in)  :: name
    real(kind=c_float),             dimension(:,:,:), intent(out) :: array
    integer,          optional,         intent(out) :: return_error
    integer :: ierr
    integer :: ncid
    integer :: varid
    integer, dimension(4) :: nstart
    integer, dimension(4) :: ncount
! #ifdef _PARALLEL_
!     ierr = nf90_open_par(trim(restart_filename_remember), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
! #else
    ierr = nf90_open(trim(restart_filename_remember), NF90_NOWRITE, ncid)
    !#endif
    call error_handler(ierr, "GET_FROM_RESTART: Problem opening restart file '"//trim(restart_filename_remember)//"'")
    nstart = (/parallel_xstart-subwindow_xstart+1,1, 1, 1/)
    ncount = (/parallel_xend-parallel_xstart+1, size(array,2), size(array,3), 1/)
    if (present(return_error)) then
       ierr = nf90_inq_varid(ncid, name, varid)
       if (ierr == NF90_NOERR) then
          return_error = 0
          call error_handler(ierr, "Problem finding variable in restart file '"//trim(name)//"'")
          ierr = nf90_get_var(ncid, varid, array, start=nstart(1:4))
          call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
       else
          return_error = 1
          write(*,'("Did not find optional variable ''",A,"'' in restart file ''", A, "''")') trim(name), trim(restart_filename_remember)
       endif
    else
       ierr = nf90_inq_varid(ncid, name, varid)
       call error_handler(ierr, "Problem finding required variable in restart file: '"//trim(name)//"'")
       ierr = nf90_get_var(ncid, varid, array, start=nstart(1:4))
       call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
    endif
    ierr = nf90_close(ncid)
    call error_handler(ierr, "Problem closing restart file")
    
  end subroutine get_from_restart_3d

  subroutine get_from_restart_2d_float(restart_filename_remember, parallel_xstart, parallel_xend, subwindow_xstart, ixfull, jxfull, name, array, return_error)
    implicit none
    character(len=*), intent(in) :: restart_filename_remember
    integer,                            intent(in) :: parallel_xstart
    integer,                            intent(in) :: parallel_xend
    integer,                            intent(in) :: subwindow_xstart
    integer,                            intent(in) :: ixfull
    integer,                            intent(in) :: jxfull
    character(len=*),                   intent(in)  :: name
    real,             dimension(parallel_xstart:parallel_xend,jxfull), intent(out) :: array
    integer,          optional,         intent(out) :: return_error
    integer :: ierr
    integer :: ncid
    integer :: varid
    integer, dimension(4) :: nstart
    integer, dimension(4) :: ncount
    integer :: rank
! #ifdef _PARALLEL_
!     call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
!     if (ierr /= MPI_SUCCESS) stop "MPI_COMM_RANK"
!     ierr = nf90_open_par(trim(restart_filename_remember), NF90_NOWRITE, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)
! #else
    rank = 0
    ierr = nf90_open(trim(restart_filename_remember), NF90_NOWRITE, ncid)
!#endif
    call error_handler(ierr, "GET_FROM_RESTART: Problem opening restart file '"//trim(restart_filename_remember)//"'")
    nstart = (/ parallel_xstart-subwindow_xstart+1, 1,  1, -99999 /)
    ncount = (/ parallel_xend-parallel_xstart+1,   size(array,2),  1, -99999 /)
    if (present(return_error)) then
       ierr = nf90_inq_varid(ncid, name, varid)
       if (ierr == NF90_NOERR) then
          return_error = 0
          call error_handler(ierr, "Problem finding variable in restart file '"//trim(name)//"'")
          ierr = nf90_get_var(ncid, varid, array, start=nstart(1:3))
          call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
       else
          return_error = 1
          if (rank == 0) write(*,'("Did not find optional variable ''",A,"'' in restart file ''", A, "''")') trim(name), trim(restart_filename_remember)
       endif
    else
       ierr = nf90_inq_varid(ncid, name, varid)
       call error_handler(ierr, "Problem finding required variable in restart file: '"//trim(name)//"'")
       ierr = nf90_get_var(ncid, varid, array, start=nstart(1:3))
       call error_handler(ierr, "Problem finding variable in restart file: '"//trim(name)//"'")
    endif
    ierr = nf90_close(ncid)
    call error_handler(ierr, "Problem closing restart file")
    
  end subroutine get_from_restart_2d_float

! --------------------------------------------------------------------------------------------------

end module wrf_hydro_nwm_jedi_field_mod
