/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxControl.h"

#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {

// ----------------------------------------------------------------------------

  ModelAuxControl::ModelAuxControl(const Geometry &,
                                   const eckit::Configuration &) {
    util::abor1_cpp(
      "ModelAuxControl::ModelAuxControl() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  ModelAuxControl::ModelAuxControl(const Geometry &,
                                   const ModelAuxControl &) {
    util::abor1_cpp(
      "ModelAuxControl::ModelAuxControl() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  ModelAuxControl::ModelAuxControl(const ModelAuxControl &, const bool) {
    util::abor1_cpp(
      "ModelAuxControl::ModelAuxControl() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  ModelAuxControl::~ModelAuxControl() {
    util::abor1_cpp(
      "ModelAuxControl::~ModelAuxControl() needs to be implemented.",
      __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void ModelAuxControl::print(std::ostream & os) const {
    util::abor1_cpp("ModelAuxControl::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the ModelAuxControl here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
