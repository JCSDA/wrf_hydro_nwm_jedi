/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_MODELAUX_MODELAUXCOVARIANCE_H_
#define WRF_HYDRO_NWM_JEDI_MODELAUX_MODELAUXCOVARIANCE_H_

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class ModelAuxControl;
  class ModelAuxIncrement;
}

// -----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

class ModelAuxCovariance : public util::Printable,
                           private util::ObjectCounter<ModelAuxCovariance> {
 public:
  static const std::string classname() {return "wrf_hydro_nwm-jedi::ModelAuxCovariance";}

/// Constructor, destructor
  ModelAuxCovariance(const eckit::Configuration & conf, const Geometry &): conf_(conf) {}
  ~ModelAuxCovariance() {}

/// Linear algebra operators
  void linearize(const ModelAuxControl &, const Geometry &) {}
  void multiply(const ModelAuxIncrement &, ModelAuxIncrement &) const {}
  void inverseMultiply(const ModelAuxIncrement &, ModelAuxIncrement &) const {}
  void randomize(ModelAuxIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_MODELAUX_MODELAUXCOVARIANCE_H_
