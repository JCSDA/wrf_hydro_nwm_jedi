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
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Traits.h"
#include "VarChaModel2GeoVaLsFortran.h"
#include "oops/base/VariableChangeBase.h"


// Forward declarations
namespace eckit {
  class Configuration;
}

namespace wrf_hydro_nwm_jedi {

// -------------------------------------------------------------------------------------------------

class VarChaModel2GeoVaLs: public oops::VariableChangeBase<Traits>, 
                           private util::ObjectCounter<VarChaModel2GeoVaLs> {
 public:
  static const std::string classname() {return "wrf_hydro_nwm_jedi::VarChaModel2GeoVaLs";}

  explicit VarChaModel2GeoVaLs(const Geometry &, const eckit::Configuration &);
  ~VarChaModel2GeoVaLs();
  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  F90vc_M2G keyFtnConfig_;
  std::shared_ptr<const Geometry> geom_;
  void print(std::ostream &) const override;
};
// -------------------------------------------------------------------------------------------------
}  // namespace wrf_hydro_nwm_jedi
