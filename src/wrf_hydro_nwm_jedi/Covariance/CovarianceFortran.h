/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_COVARIANCE_COVARIANCEFORTRAN_H_
#define WRF_HYDRO_NWM_JEDI_COVARIANCE_COVARIANCEFORTRAN_H_

#include "wrf_hydro_nwm_jedi/Increment/Increment.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace wrf_hydro_nwm_jedi {

  typedef int F90bmat;

extern "C" {

  void wrf_hydro_nwm_jedi_b_setup_f90(F90bmat &,
                                      const F90state &,
                                      const eckit::Configuration * const *,
                                      const oops::Variables &);
  void wrf_hydro_nwm_jedi_b_delete_f90(F90bmat &);
  /* void wrf_hydro_nwm_jedi_b_linearize_f90(const F90bmat &, */
  /*                         const eckit::Configuration * const *); */
  void wrf_hydro_nwm_jedi_b_mult_f90(const F90bmat &,
                                     const F90inc &,
                                     const F90inc &);
  void wrf_hydro_nwm_jedi_b_randomize_f90(const F90bmat &,
                                          const F90inc &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
#endif  // WRF_HYDRO_NWM_JEDI_COVARIANCE_COVARIANCEFORTRAN_H_
