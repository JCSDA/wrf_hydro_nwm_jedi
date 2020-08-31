! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Get values module
module wrf_hydro_nwm_jedi_getvalues_mod

! fckit
use fckit_mpi_module,              only: fckit_mpi_comm

! oops
use datetime_mod,                   only: datetime
!use type_bump,                      only: bump_type
!use unstructured_interpolation_mod, only: unstrc_interp

! ufo
use ufo_locs_mod,                   only: ufo_locs, ufo_locs_time_mask
use ufo_geovals_mod,                only: ufo_geovals

! fv3jedi uses
!use fv3jedi_constants_mod,          only: rad2deg
!use fv3jedi_bump_mod,               only: bump_init, bump_apply
use wrf_hydro_nwm_jedi_fields_mod,  only: wrf_hydro_nwm_jedi_fields, base_field
use wrf_hydro_nwm_jedi_geometry_mod,only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_util_mod,    only: indices
!use fv3jedi_interpolation_mod,      only: unsinterp_integer_apply, unsinterp_nearest_apply
!use fv3jedi_kinds_mod,              only: kind_real
use iso_c_binding, only: c_float
use wrf_hydro_nwm_jedi_state_mod,        only: wrf_hydro_nwm_jedi_state

! --------------------------------------------------------------------------------------------------

implicit none
private
public :: wrf_hydro_nwm_jedi_getvalues, wrf_hydro_nwm_jedi_getvalues_base

type, abstract, public :: wrf_hydro_nwm_jedi_getvalues_base
   integer              :: isc, iec, jsc, jec, npz, ngrid
   character(len=2048)  :: interp_method
   !type(bump_type)      :: bump
   integer              :: nnear = 4
   !type(unstrc_interp)  :: unsinterp
   type(fckit_mpi_comm) :: comm
 contains
   procedure, public :: create
   procedure, public :: delete
   procedure, public :: fill_geovals
   procedure, public :: fill_geovals_ad
   generic, public :: set_trajectory => fill_geovals
   generic, public :: fill_geovals_tl => fill_geovals
end type wrf_hydro_nwm_jedi_getvalues_base


type, public, extends(wrf_hydro_nwm_jedi_getvalues_base) :: wrf_hydro_nwm_jedi_getvalues
end type wrf_hydro_nwm_jedi_getvalues


! --------------------------------------------------------------------------------------------------
contains
! --------------------------------------------------------------------------------------------------

  
subroutine create(self, geom, locs)
  class(wrf_hydro_nwm_jedi_getvalues_base), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry),        intent(in)    :: geom
  type(ufo_locs),                           intent(in)    :: locs

  ! Geometry
  ! --------
  ! self%isc = geom%isc
  ! self%iec = geom%iec
  ! self%jsc = geom%jsc
  ! self%jec = geom%jec
  ! self%npz = geom%npz
  ! self%ngrid = geom%ngrid
  ! self%comm = geom%f_comm
  
  ! Create interpolation weights
  ! ----------------------------
  ! self%interp_method = trim(geom%interp_method)
  ! if (trim(self%interp_method) == 'bump') then
  !   call bump_init(geom, locs%nlocs, locs%lat, locs%lon, self%bump, 55555)
  ! endif
  
  ! Always create unstructured interpolation as it is used for special case fields, e.g. integers
  ! call self%unsinterp%create( geom%f_comm, self%nnear, 'barycent', &
  !                             self%ngrid, rad2deg*geom%lat_us, rad2deg*geom%lon_us, &
  !                             locs%nlocs, locs%lat, locs%lon )
end subroutine create


subroutine delete(self)
  class(wrf_hydro_nwm_jedi_getvalues_base), intent(inout) :: self

  ! if (trim(self%interp_method) == 'bump') call self%bump%dealloc

  ! call self%unsinterp%delete()
end subroutine delete


!> Set the geovals from the increment/state
subroutine fill_geovals(self, geom, fields_obj, t1, t2, locs, geovals)
  class(wrf_hydro_nwm_jedi_getvalues_base), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry),        intent(in)    :: geom
  type(wrf_hydro_nwm_jedi_fields),          intent(in)    :: fields_obj
  type(datetime),                           intent(in)    :: t1
  type(datetime),                           intent(in)    :: t2
  type(ufo_locs),                           intent(in)    :: locs
  type(ufo_geovals),                        intent(inout) :: geovals
  
  integer :: gv, ii
  ! integer :: idx_1, idx_2
  class(base_field), pointer :: field
  ! character(len=10) :: wrf_hydro_nwm_name
  logical, allocatable :: time_mask(:)
  real(kind=c_float) :: lat, long
  ! real(kind=c_float), allocatable :: field_us(:)
  real(kind=c_float), allocatable :: geovals_all(:,:), geovals_tmp(:)
  type(indices) :: ind
  ! Get mask for locations in this time window
  ! ------------------------------------------
  call ufo_locs_time_mask(locs, t1, t2, time_mask)


  ! Allocate geovals
  ! ----------------
  if (.not. geovals%linit) then
     do gv = 1, geovals%nvar
        geovals%geovals(gv)%nval = 1!fields(gv)%dim3_len
        allocate(geovals%geovals(gv)%vals(geovals%geovals(gv)%nval, geovals%geovals(gv)%nlocs))
        geovals%geovals(gv)%vals = 0.0
     enddo
  endif
  geovals%linit = .true.

  self%ngrid = 1

  ! ! Loop over GeoVaLs
  ! ! -----------------
  ! allocate(field_us(self%ngrid))
  ! allocate(geovals_all(locs%nlocs, self%npz+1))
  allocate(geovals_all(locs%nlocs, 1+1))
  allocate(geovals_tmp(locs%nlocs))

  do gv = 1, geovals%nvar

     !   ! Get GeoVaLs field
     !   ! -----------------
     !   call long_name_to_wrf_hydro_name(fields, trim(geovals%variables(gv)), wrf_hydro_nwm_name)

     !  call pointer_field(fields, wrf_hydro_nwm_name, field)

     call fields_obj%search_field(trim(geovals%variables(gv)),field)

     !   ! Interpolation
     !   ! -------------
     geovals_all = 0.0
     geovals_tmp = 0.0
     ! Work on geovals_tmp

     do ii = 1, locs%nlocs
        lat = locs%lat(ii)
        long = locs%lon(ii)
        call geom%get_lsm_nn(lat, long, ind)
        geovals_tmp(ii) = field%get_value(ind)!field%array(idx_1, idx_2, 1)
     end do

     geovals_all(1:locs%nlocs, 1) = geovals_tmp(1:locs%nlocs)   

     ! Can optionally interpolate real valued magnitude fields with bump
     ! -----------------------------------------------------------------
     ! if ( trim(self%interp_method) == 'bump' .and. &
     !      .not.field%integerfield .and. trim(field%space)=='magnitude' ) then

     !   call bump_apply(field%npz, geom, field%array, locs%nlocs, geovals_all(:,1:field%npz), self%bump)

     ! else ! Otherwise use unstructured interpolation

     !   do jlev = 1, field%npz
     !     n = 0
     !     do jj = field%jsc, field%jec
     !       do ji = field%isc, field%iec
     !         n = n + 1
     !         field_us(n) = field%array(ji, jj, jlev)
     !       enddo
     !     enddo

     !     ! Conditions for integer and directional fields
     !     ! ---------------------------------------------
     !     if (.not. field%integerfield .and. trim(field%space)=='magnitude') then
     !       call self%unsinterp%apply(field_us, geovals_tmp)
     !     elseif (field%integerfield) then
     !       call unsinterp_integer_apply(self%unsinterp, field_us, geovals_tmp)
     !     elseif (trim(field%space)=='direction') then
     !       call unsinterp_nearest_apply(self%unsinterp, field_us, geovals_tmp)
     !     else
     !       call abor1_ftn("fv3jedi_getvalues_mod.fill_geovals: interpolation for this kind of "// &
     !                      "field is not supported. FieldName: "// trim(field%fv3jedi_name))
     !     endif
     !     geovals_all(1:locs%nlocs, jlev) = geovals_tmp(1:locs%nlocs)
     !   enddo

     ! endif

     ! Fill GeoVaLs relevant to this window
     ! ------------------------------------
     do ii = 1,locs%nlocs
        if (time_mask(ii)) geovals%geovals(gv)%vals(1:1, ii) = geovals_all(ii, 1:1)
     enddo
  enddo

  write(*,*) "End of fill_geovals"

  ! deallocate(field_us)
  deallocate(geovals_all)
  deallocate(geovals_tmp)
  deallocate(time_mask)
end subroutine fill_geovals


!> Set the increment/state from the geovals
subroutine fill_geovals_ad(self, geom, fields_obj, t1, t2, locs, geovals)
  class(wrf_hydro_nwm_jedi_getvalues_base), intent(inout) :: self
  type(wrf_hydro_nwm_jedi_geometry),        intent(in)    :: geom
  type(wrf_hydro_nwm_jedi_fields),          intent(inout) :: fields_obj
  type(datetime),                           intent(in)    :: t1
  type(datetime),                           intent(in)    :: t2
  type(ufo_locs),                           intent(in)    :: locs
  type(ufo_geovals),                        intent(in)    :: geovals

  integer :: gv, ii
  class(base_field), pointer :: field
  ! character(len=10) :: wrf_hydro_nwm_name
  logical, allocatable :: time_mask(:)
  real(kind=c_float) :: lat, long
  ! real(kind=c_float), allocatable :: field_us(:)
  ! real(kind=c_float), allocatable :: geovals_all(:,:), geovals_tmp(:)
  type(indices) :: ind
  ! Get mask for locations in this time window
  ! ------------------------------------------
  call ufo_locs_time_mask(locs, t1, t2, time_mask)

  self%ngrid = 1

  do gv = 1, geovals%nvar

     call fields_obj%search_field(trim(geovals%variables(gv)), field)
     do ii = 1, locs%nlocs
        if (time_mask(ii)) then
           lat = locs%lat(ii)
           long = locs%lon(ii)
           call geom%get_lsm_nn(lat, long, ind)
           call field%set_value(ind, real(geovals%geovals(gv)%vals(1, ii), kind=c_float))
        endif
     end do

     ! Can optionally interpolate real valued magnitude fields with bump
     ! -----------------------------------------------------------------
     ! if ( trim(self%interp_method) == 'bump' .and. &
     !      .not.field%integerfield .and. trim(field%space)=='magnitude' ) then

     !   call bump_apply(field%npz, geom, field%array, locs%nlocs, geovals_all(:,1:field%npz), self%bump)

     ! else ! Otherwise use unstructured interpolation

     !   do jlev = 1, field%npz
     !     n = 0
     !     do jj = field%jsc, field%jec
     !       do ji = field%isc, field%iec
     !         n = n + 1
     !         field_us(n) = field%array(ji, jj, jlev)
     !       enddo
     !     enddo

     !     ! Conditions for integer and directional fields
     !     ! ---------------------------------------------
     !     if (.not. field%integerfield .and. trim(field%space)=='magnitude') then
     !       call self%unsinterp%apply(field_us, geovals_tmp)
     !     elseif (field%integerfield) then
     !       call unsinterp_integer_apply(self%unsinterp, field_us, geovals_tmp)
     !     elseif (trim(field%space)=='direction') then
     !       call unsinterp_nearest_apply(self%unsinterp, field_us, geovals_tmp)
     !     else
     !       call abor1_ftn("fv3jedi_getvalues_mod.fill_geovals: interpolation for this kind of "// &
     !                      "field is not supported. FieldName: "// trim(field%fv3jedi_name))
     !     endif
     !     geovals_all(1:locs%nlocs, jlev) = geovals_tmp(1:locs%nlocs)
     !   enddo

     ! endif

     ! Fill GeoVaLs relevant to this window
     ! ------------------------------------
     ! do n = 1,locs%nlocs
     !    if (time_mask(n)) geovals%geovals(gv)%vals(1:1, n) = geovals_all(n, 1:1)
     ! enddo
  enddo

  write(*,*) "End of fill_geovals_ad"

  ! deallocate(field_us)
  ! deallocate(geovals_all)
  ! deallocate(geovals_tmp)
  deallocate(time_mask)
end subroutine fill_geovals_ad


end module wrf_hydro_nwm_jedi_getvalues_mod
