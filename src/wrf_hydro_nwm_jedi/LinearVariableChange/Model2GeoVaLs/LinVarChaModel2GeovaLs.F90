! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module wrf_hydro_nwm_jedi_lvc_model2geovals_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,           only: fckit_log

use datetime_mod

use wrf_hydro_nwm_jedi_geometry_mod,      only: wrf_hydro_nwm_jedi_geometry
use wrf_hydro_nwm_jedi_state_mod,         only: wrf_hydro_nwm_jedi_state

implicit none

private
public :: wrf_hydro_nwm_jedi_lvc_model2geovals

type :: wrf_hydro_nwm_jedi_lvc_model2geovals
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: multiply
    procedure, public :: multiplyadjoint
end type wrf_hydro_nwm_jedi_lvc_model2geovals

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, dummyconf)

class(wrf_hydro_nwm_jedi_lvc_model2geovals), intent(inout) :: self
type(wrf_hydro_nwm_jedi_geometry),           intent(in)    :: geom
type(wrf_hydro_nwm_jedi_state),              intent(in)    :: bg
type(wrf_hydro_nwm_jedi_state),              intent(in)    :: fg
type(fckit_configuration),                   intent(in)    :: dummyconf


end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(wrf_hydro_nwm_jedi_lvc_model2geovals), intent(inout) :: self

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, geom, dxm, dxg)

class(wrf_hydro_nwm_jedi_lvc_model2geovals), intent(inout) :: self
type(wrf_hydro_nwm_jedi_geometry),           intent(inout) :: geom
type(wrf_hydro_nwm_jedi_state),          intent(in)    :: dxm
type(wrf_hydro_nwm_jedi_state),          intent(inout) :: dxg

    call abor1_ftn("wrf_hydro_nwm_jedi_lvc_model2geovals_mod.multiply not implemented.")

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiplyadjoint(self, geom, dxg, dxm)

class(wrf_hydro_nwm_jedi_lvc_model2geovals), intent(inout) :: self
type(wrf_hydro_nwm_jedi_geometry),               intent(inout) :: geom
type(wrf_hydro_nwm_jedi_state),          intent(in)    :: dxg
type(wrf_hydro_nwm_jedi_state),          intent(inout) :: dxm

    call abor1_ftn("wrf_hydro_nwm_jedi_lvc_model2geovals_mod.multiplyadjoint not implemented.")

end subroutine multiplyadjoint

! --------------------------------------------------------------------------------------------------

end module wrf_hydro_nwm_jedi_lvc_model2geovals_mod
