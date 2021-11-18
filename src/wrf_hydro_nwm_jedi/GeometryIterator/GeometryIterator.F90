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
    
      !> \copybrief soca_geom_iter_setup \see soca_geom_iter_setup
      procedure :: init => wrf_hydro_nwm_jedi_geometry_iter_init
    
      !> \copybrief soca_geom_iter_clone \see soca_geom_iter_clone
      procedure :: clone => wrf_hydro_nwm_jedi_geometry_iter_clone
    
      !> \copybrief soca_geom_iter_equals \see soca_geom_iter_equals
!      procedure :: equals => soca_geom_iter_equals
    
      !> \copybrief soca_geom_iter_current \see soca_geom_iter_current
!      procedure :: current => soca_geom_iter_current
    
      !> \copybrief soca_geom_iter_next \see soca_geom_iter_next
!      procedure :: next => soca_geom_iter_next
    
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

end module wrf_hydro_nwm_jedi_geometry_iter_mod