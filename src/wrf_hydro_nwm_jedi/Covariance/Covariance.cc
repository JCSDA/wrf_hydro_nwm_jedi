/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "wrf_hydro_nwm_jedi/Covariance/Covariance.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/State/State.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {

// ----------------------------------------------------------------------------

  Covariance::Covariance(const Geometry & geom,
                         const oops::Variables & vars,
                         const eckit::Configuration & conf,
                         const State & x1, const State & x2) {
    util::abor1_cpp("Covariance::Covariance() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  Covariance::~Covariance() {
    util::abor1_cpp("Covariance::~Covariance() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Covariance::multiply(const Increment & dxin, Increment & dxout) const {
    util::abor1_cpp("Covariance::multiply() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Covariance::inverseMultiply(const Increment &, Increment &) const {
    util::abor1_cpp("Covariance::inverseMultiply() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Covariance::randomize(Increment &) const {
    util::abor1_cpp("Covariance::randomize() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Covariance::print(std::ostream & os) const {
    util::abor1_cpp("Covariance::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the Covariance here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
