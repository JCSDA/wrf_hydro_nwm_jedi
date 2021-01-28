!> General utilities module
module wrf_hydro_nwm_jedi_constants_mod

use iso_c_binding, only: c_int, c_float, c_double
use netcdf

implicit none

public

! named numeric constants
real(c_float),  parameter :: zero_c_float = 0.d0
real(c_double), parameter :: zero_c_double = 0.d0
real(c_float),  parameter :: one_c_float = 1.d0
real(c_double), parameter :: one_c_double = 1.d0

! conventions
integer,        parameter :: unopened_ncid = -9999

contains

end module wrf_hydro_nwm_jedi_constants_mod
