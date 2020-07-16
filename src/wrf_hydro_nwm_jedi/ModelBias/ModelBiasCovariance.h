/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_MODELBIASCOVARIANCE_H_
#define WRF_HYDRO_NWM_JEDI_MODELBIASCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace wrf_hydro_nwm_jedi {
  class ModelBias;
  class ModelBiasIncrement;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasCovariance : public util::Printable,
  private boost::noncopyable,
  private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname()
                                 {return "wrf_hydro_nwm_jedi::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration & conf,
		      const Geometry &): conf_(conf) {}
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const Geometry &) {}
  void multiply(const ModelBiasIncrement &,
                ModelBiasIncrement) const {}
  void inverseMultiply(const ModelBiasIncrement &,
                ModelBiasIncrement) const {}
  void randomize(ModelBiasIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi

#endif  // JEDI_SRC_MODELBIAS_MODELBIASCOVARIANCESW_H_
