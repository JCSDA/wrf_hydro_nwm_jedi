/*
 * (C) Copyright 2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "wrf_hydro_nwm_jedi/LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLsFortran.h"
#include "wrf_hydro_nwm_jedi/Traits.h"


// Forward declarations
namespace eckit {
  class Configuration;
}

namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class State;
  class Increment;

// -------------------------------------------------------------------------------------------------

class LinVarChaModel2GeoVaLs: public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "wrf_hydro_nwm_jedi::LinVarChaModel2GeoVaLs";}
  explicit LinVarChaModel2GeoVaLs(const State &, const State &, const Geometry &,
                                  const eckit::LocalConfiguration &);
  ~LinVarChaModel2GeoVaLs();
  void multiply(const Increment &, Increment &) const override;
  void multiplyInverse(const Increment &, Increment &) const override;
  void multiplyAD(const Increment &, Increment &) const override;
  void multiplyInverseAD(const Increment &, Increment &) const override;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90lvc_M2G keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -------------------------------------------------------------------------------------------------
}  // namespace wrf_hydro_nwm_jedi
