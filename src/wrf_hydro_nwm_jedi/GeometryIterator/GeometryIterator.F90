! (C) Copyright 2019-2021 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Geometry iterator
module wrf_hydro_nwm_jedi_geometry_iter_mod

    use kinds, only : kind_real
    use wrf_hydro_nwm_jedi_geometry_mod, only: wrf_hydro_nwm_jedi_geometry
    
    implicit none
    private
    
    
    ! ------------------------------------------------------------------------------
    ! ------------------------------------------------------------------------------
    
    ! ------------------------------------------------------------------------------
    !> Geometry iterator
    !!
    !! When initialized, the iterator points to the first valid local grid cell.
    !! Calls to wrf_hydro_nwm_jedi_geometry_iter::next() moves the iterator forward, and calls to
    !! wrf_hydro_nwm_jedi_geometry_iter::current() retrieves the lat/lon of the current grid cell.
    !! The iterator is mainly used by soca_increment_mod::soca_increment::getpoint()
    !! and soca_increment_mod::soca_increment::setpoint()
    type, public :: wrf_hydro_nwm_jedi_geometry_iter
      type(wrf_hydro_nwm_jedi_geometry), pointer :: geom => null() !< Geometry
    
      integer :: iind = 1  !< i index of current grid point
      integer :: jind = 1  !< j index of current grid point
    
    contains
    
      !> \copybrief wrf_hydro_nwm_jedi_geom_iter_setup \see wrf_hydro_nwm_jedi_geom_iter_setup
      procedure :: init => wrf_hydro_nwm_jedi_geometry_iter_init
    
      !> \copybrief wrf_hydro_nwm_jedi_geom_iter_clone \see wrf_hydro_nwm_jedi_geom_iter_clone
      procedure :: clone => wrf_hydro_nwm_jedi_geometry_iter_clone
    
      !> \copybrief wrf_hydro_nwm_jedi_geom_iter_equals \see wrf_hydro_nwm_jedi_geom_iter_equals
      procedure :: equals => wrf_hydro_nwm_jedi_geometry_iter_equals
    
      !> \copybrief wrf_hydro_nwm_jedi_geom_iter_current \see wrf_hydro_nwm_jedi_geom_iter_current
      procedure :: current => wrf_hydro_nwm_jedi_geometry_iter_current

      !> \copybrief wrf_hydro_nwm_jedi_geom_iter_next \see wrf_hydro_nwm_jedi_geom_iter_next
      procedure :: orog => wrf_hydro_nwm_jedi_geometry_iter_orography      
    
      !> \copybrief wrf_hydro_nwm_jedi_geom_iter_next \see wrf_hydro_nwm_jedi_geom_iter_next
      procedure :: next => wrf_hydro_nwm_jedi_geom_iter_next
    
    end type wrf_hydro_nwm_jedi_geometry_iter
        
    ! ------------------------------------------------------------------------------
    contains
    ! ------------------------------------------------------------------------------
    
    
    ! ------------------------------------------------------------------------------
    !> Setup for the geometry iterator
    
    subroutine wrf_hydro_nwm_jedi_geometry_iter_init(self, geom, iind, jind)
      class(wrf_hydro_nwm_jedi_geometry_iter),    intent(inout) :: self
      type(wrf_hydro_nwm_jedi_geometry), pointer, intent(   in) :: geom !< Pointer to geometry
      integer,                                    intent(   in) :: iind, jind  !< starting index
    
      ! Associate geometry
      self%geom => geom
    
      ! Define iind/jind for local tile
      self%iind = iind
      self%jind = jind
    
    end subroutine wrf_hydro_nwm_jedi_geometry_iter_init

    ! ------------------------------------------------------------------------------
!> Clone for the geometry iterator from \p other to \p self

subroutine wrf_hydro_nwm_jedi_geometry_iter_clone(self, other)
    class(wrf_hydro_nwm_jedi_geometry_iter), intent(inout) :: self
    type(wrf_hydro_nwm_jedi_geometry_iter),  intent(   in) :: other !< Other geometry iterator to clone from
  
    ! Associate geometry
    self%geom => other%geom
  
    ! Copy iind/jind
    self%iind = other%iind
    self%jind = other%jind
  
  end subroutine wrf_hydro_nwm_jedi_geometry_iter_clone

! ------------------------------------------------------------------------------
!> Check for the geometry iterator equality (pointing to same i/j location)
!!
!! \relates soca_geom_iter_mod::soca_geom_iter
subroutine wrf_hydro_nwm_jedi_geometry_iter_equals(self, other, equals)
    class(wrf_hydro_nwm_jedi_geometry_iter), intent( in) :: self
    type(wrf_hydro_nwm_jedi_geometry_iter),  intent( in) :: other  !< Other geometry iterator
    integer,               intent(out) :: equals !< Equality flag
  
    ! Initialization
    equals = 0
  
    ! Check equality
    if (associated(self%geom, other%geom) .and. (self%iind==other%iind) &
        .and. (self%jind==other%jind)) equals = 1
  
  end subroutine wrf_hydro_nwm_jedi_geometry_iter_equals  

! ------------------------------------------------------------------------------
!> Get geometry iterator current lat/lon
!!
!! \throws abor1_ftn aborts if iterator is out of bounds
!! \relates wrf_hydro_nwm_jedi_geom_iter_mod::wrf_hydro_nwm_jedi_geom_iter
  subroutine wrf_hydro_nwm_jedi_geometry_iter_current(self, lon, lat)
    class(wrf_hydro_nwm_jedi_geometry_iter), intent( in) :: self
    real(kind_real),                         intent(out) :: lat  !< Latitude
    real(kind_real),                         intent(out) :: lon  !< Longitude
  
    ! Check iind/jind
    if (self%iind == -1 .AND. self%jind == -1) then
      ! special case of {-1,-1} means end of the grid
      lat = self%geom%lsm%lat(self%geom%lsm%xdim_len,self%geom%lsm%ydim_len)
      lon = self%geom%lsm%lon(self%geom%lsm%xdim_len,self%geom%lsm%ydim_len)
    elseif (self%iind < 1 .OR. self%iind > self%geom%lsm%xdim_len .OR. &
            self%jind < 1 .OR. self%jind > self%geom%lsm%ydim_len) then
      ! outside of the grid
      call abor1_ftn('wrf_hydro_nwm_jedi_geometry_iter_current: iterator out of bounds')
    else
      ! inside of the grid
      lat = self%geom%lsm%lat(self%iind,self%jind)
      lon = self%geom%lsm%lon(self%iind,self%jind)
    endif
  
  end subroutine wrf_hydro_nwm_jedi_geometry_iter_current  

! ------------------------------------------------------------------------------
  !> Get geometry iterator current surface elevation
  subroutine wrf_hydro_nwm_jedi_geometry_iter_orography(self, oro)

    ! Passed variables
    class(wrf_hydro_nwm_jedi_geometry_iter), intent( in) :: self
    real(kind_real),    intent(out) :: oro  !< Orography

    ! Check iindex/jindex
    if (self%iind == -1 .AND. self%jind == -1) then
      ! special case of {-1,-1} means end of the grid
      oro = self%geom%lsm%sfc_elev(self%geom%lsm%xdim_len,self%geom%lsm%ydim_len)
    elseif (self%iind < 1 .OR. self%iind > self%geom%lsm%xdim_len .OR. &
        self%jind < 1 .OR. self%jind > self%geom%lsm%ydim_len) then
  ! outside of the grid
      call abor1_ftn('wrf_hydro_nwm_jedi_geom_iter_orography: iterator out of bounds')
    else
      ! inside of the grid
      oro = self%geom%lsm%sfc_elev(self%iind,self%jind)
    endif

  end subroutine wrf_hydro_nwm_jedi_geometry_iter_orography

! ------------------------------------------------------------------------------
!> Update geometry iterator to next point
!!
!! \todo skip over masked points
!! \relates _geom_iter_mod::soca_geom_iter
  subroutine wrf_hydro_nwm_jedi_geom_iter_next(self)
    class(wrf_hydro_nwm_jedi_geometry_iter), intent(inout) :: self
    integer :: iind, jind
  
    iind = self%iind
    jind = self%jind
  

      ! increment by 1
      if (iind.lt.self%geom%lsm%xdim_len) then
        iind = iind + 1
      elseif (iind.eq.self%geom%lsm%xdim_len) then
        iind = 1
        jind = jind + 1
      end if
  
    if (jind > self%geom%lsm%ydim_len) then
        iind=-1
        jind=-1
    end if
  
    self%iind = iind
    self%jind = jind
  
  end subroutine wrf_hydro_nwm_jedi_geom_iter_next


end module wrf_hydro_nwm_jedi_geometry_iter_mod