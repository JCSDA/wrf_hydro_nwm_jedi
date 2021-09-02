/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_MODEL_MODEL_H_
#define WRF_HYDRO_NWM_JEDI_MODEL_MODEL_H_

#include <memory>
#include <ostream>

#include "oops/interface/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"

#include "wrf_hydro_nwm_jedi/Traits.h"

// forward declarations
namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class ModelAuxControl;
  struct Traits;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // Model class
  class Model : public oops::interface::ModelBase<Traits>,
                private util::ObjectCounter<Model>
  {
   public:
    Model(const Geometry &, const eckit::Configuration &);
    ~Model();

    void initialize(State &) const;
    void step(State &, const ModelAuxControl &) const;
    void finalize(State &) const;

    const util::Duration & timeResolution() const { return tstep_;}
    const oops::Variables & variables() const { return vars_; }

   private:
    void print(std::ostream &) const;
    int keyConfig_;
    util::Duration tstep_;
    const std::unique_ptr<Geometry> geom_;
    const oops::Variables vars_;
  };

}  // namespace wrf_hydro_nwm_jedi
#endif  // WRF_HYDRO_NWM_JEDI_MODEL_MODEL_H_
