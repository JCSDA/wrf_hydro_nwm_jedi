/*
 * (C) Copyright 2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/Traits.h"
#include "wrf_hydro_nwm_jedi/VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLs.h"
#include "wrf_hydro_nwm_jedi/VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLsFortran.h"

namespace wrf_hydro_nwm_jedi {

// -------------------------------------------------------------------------------------------------
static VariableChangeMaker<VarChaModel2GeoVaLs>
       makerVarChaModel2GeoVaLs_("Model2GeoVaLs");
static VariableChangeMaker<VarChaModel2GeoVaLs> makerVarChaDefault_("default");
// -------------------------------------------------------------------------------------------------

VarChaModel2GeoVaLs::VarChaModel2GeoVaLs(const Geometry & geom, const eckit::Configuration & conf) :
  geom_(new Geometry(geom)) 
{
  oops::Log::trace() << "VarChaModel2GeoVaLs::VarChaModel2GeoVaLs start" << std::endl;
  const eckit::Configuration * configc = &conf;
  wrf_hydro_nwm_jedi_varchamodel2geovals_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << "VarChaModel2GeoVaLs::VarChaModel2GeoVaLs done" << std::endl;
}

VarChaModel2GeoVaLs::~VarChaModel2GeoVaLs() {
  wrf_hydro_nwm_jedi_varchamodel2geovals_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}

void VarChaModel2GeoVaLs::changeVar(const State & xin,
                                         State & xout) const {
  oops::Log::trace() << classname() << "VarChaModel2GeoVaLs::changeVar starting" << std::endl;
  wrf_hydro_nwm_jedi_varchamodel2geovals_changevar_f90(keyFtnConfig_,
                                                       geom_->toFortran(),
                                                       xin.toFortran(),
                                                       xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << "VarChaModel2GeoVaLs::changeVar done" << std::endl;
  }

  void VarChaModel2GeoVaLs::changeVarInverse(const State & xg, State & xm) const {
    util::abor1_cpp("VarChaModel2GeoVaLs::changeVarInverse not implemented");
  }

  void VarChaModel2GeoVaLs::print(std::ostream & os) const {
    os << classname() << "VarChaModel2GeoVaLs";
  }

}  // namespace wrf_hydro_nwm_jedi
