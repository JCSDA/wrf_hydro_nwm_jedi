module wrf_hydro_nwm-jedi_geometry_mod

use fckit_configuration_module, only: fckit_configuration

implicit none
private

public :: wrf_hydro_nwm-jedi_geometry

!------------------------------------------------------------------------------

type :: wrf_hydro_nwm-jedi_geometry
contains
  procedure :: init   => wrf_hydro_nwm-jedi_geometry_init
  procedure :: clone  => wrf_hydro_nwm-jedi_geometry_clone
  procedure :: delete => wrf_hydro_nwm-jedi_geometry_delete
end type

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm-jedi_geometry_init(self, f_conf)
  class(wrf_hydro_nwm-jedi_geometry),   intent(out) :: self
  type(fckit_configuration),  intent(in) :: f_conf

  call abor1_ftn("ERROR: wrf_hydro_nwm-jedi_geometry_init() needs to be implemented.")
end subroutine

!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm-jedi_geometry_clone(self, other)
  class(wrf_hydro_nwm-jedi_geometry),  intent(out) :: self
  class(wrf_hydro_nwm-jedi_geometry),   intent(in) :: other
      
  call abor1_ftn("ERROR: wrf_hydro_nwm-jedi_geometry_clone() needs to be implemented.")
end subroutine
  
!------------------------------------------------------------------------------

subroutine wrf_hydro_nwm-jedi_geometry_delete(self)
  class(wrf_hydro_nwm-jedi_geometry),  intent(out) :: self
       
  call abor1_ftn("ERROR: wrf_hydro_nwm-jedi_geometry_delete() needs to be implemented.")
end subroutine

!------------------------------------------------------------------------------

end module