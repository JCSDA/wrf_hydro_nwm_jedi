/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm-jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm-jedi/ModelAux/ModelAuxCovariance.h"

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm-jedi {

// ----------------------------------------------------------------------------

  ModelAuxCovariance::ModelAuxCovariance(const eckit::Configuration & conf,
                                         const Geometry & geom)
  : conf_(conf) {
    util::abor1_cpp(
      "ModelAuxCovariance::ModelAuxCovariance() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  ModelAuxCovariance::~ModelAuxCovariance() {
    util::abor1_cpp(
      "ModelAuxCovariance::~ModelAuxCovariance() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::linearize(const ModelAuxControl &,
                                     const Geometry &) {
    util::abor1_cpp(
      "ModelAuxCovariance::linearize() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::multiply(const ModelAuxIncrement &,
                                    ModelAuxIncrement &) {
    util::abor1_cpp(
      "ModelAuxCovariance::multiply() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::inverseMultiply(const ModelAuxIncrement &,
                                           ModelAuxIncrement &) const {
    util::abor1_cpp(
      "ModelAuxCovariance::inverseMultiply() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::randomize(ModelAuxIncrement &) const {
    util::abor1_cpp(
      "ModelAuxCovariance::randomize() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void ModelAuxCovariance::print(std::ostream & os) const {
    util::abor1_cpp("ModelAuxCovariance::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the ModelAuxCovariance here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm-jedi
