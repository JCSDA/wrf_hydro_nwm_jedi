/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Model/Model.h"
#include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxControl.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/Traits.h"

#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {

  static oops::interface::ModelMaker<Traits, Model> modelmaker_("WRF_HYDRO_NWM_JEDI");

  Model::Model(const Geometry & geom,
               const eckit::Configuration & conf)
    : keyConfig_(0),
      tstep_("PT1H"),
      geom_(new Geometry(geom)),
      vars_(conf, "model variables") {
    std::cout << "Creating model" << std::endl;

    // adapted from soca
    // Log::trace() << "Model::Model" << std::endl;
    // Log::trace() << "Model vars: " << vars_ << std::endl;
    // tstep_ = util::Duration(conf.getString("tstep"));
    // const eckit::Configuration * configc = &conf;
    // soca_setup_f90(&configc, geom_->toFortran(), keyConfig_);
    // NOTE: tstep_ should be set to the actual timestep
  }

  Model::~Model() {
    // util::abor1_cpp("Model::~Model() needs to be implemented.",
    //             __FILE__, __LINE__);
  }

  void Model::initialize(State & xx) const {
    // util::abor1_cpp("Model::initialize() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }

  void Model::step(State & xx, const ModelAuxControl & xx_bias) const {
    // util::abor1_cpp("Model::step() needs to be implemented.",
    //                 __FILE__, __LINE__);
    xx.validTime() += tstep_;
  }

  void Model::finalize(State & xx) const {
    // util::abor1_cpp("Model::finalize() needs to be implemented.",
    //                  __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Model::print(std::ostream & os) const {
    // util::abor1_cpp("Model::print() needs to be implemented.",
    //                 __FILE__, __LINE__);
    os << "Geometry: "
       << "(TODO, print diagnostic info about the geometry here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
