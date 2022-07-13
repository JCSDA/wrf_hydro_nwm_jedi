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
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/Traits.h"

#include "wrf_hydro_nwm_jedi/LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h"

namespace wrf_hydro_nwm_jedi {

// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaModel2GeoVaLs>
  makerLinVarChaModel2GeoVaLs_("Model2GeoVaLs");

static LinearVariableChangeMaker<LinVarChaModel2GeoVaLs>
         makerLinVarChaModel2GeoDef_("default");

// -------------------------------------------------------------------------------------------------

LinVarChaModel2GeoVaLs::LinVarChaModel2GeoVaLs(const State & bg, const State & fg,
                                    const Geometry & resol, const eckit::LocalConfiguration & conf)
  : LinearVariableChangeBase(), geom_(new Geometry(resol))
{
  util::Timer timer(classname(), "LinVarChaModel2GeoVaLs");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  wrf_hydro_nwm_jedi_lvc_model2geovals_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                       fg.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}

// -----------------------------------------------------------------------------

LinVarChaModel2GeoVaLs::~LinVarChaModel2GeoVaLs() {
  util::Timer timer(classname(), "~LinVarChaModel2GeoVaLs");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  wrf_hydro_nwm_jedi_lvc_model2geovals_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::changeVarTL(const Increment &dxin,
                                            Increment &dxout) const {
//  util::abor1_cpp("LinVarChaModel2GeoVaLs::multiply not implemented");
  wrf_hydro_nwm_jedi_lvc_model2geovals_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                                          dxin.toFortran(),
                                                          dxout.toFortran());
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::changeVarInverseTL(const Increment &,
                                                Increment &) const {
  util::abor1_cpp("LinVarChaModel2GeoVaLs::multiplyInverse not implemented");
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::changeVarAD(const Increment &dxin,
                                              Increment &dxout) const {
//  util::abor1_cpp("LinVarChaModel2GeoVaLs::multiplyAD not implemented");
  wrf_hydro_nwm_jedi_lvc_model2geovals_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                          dxin.toFortran(),
                                                          dxout.toFortran());
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::changeVarInverseAD(const Increment &,
                                                  Increment &) const {
  util::abor1_cpp("LinVarChaModel2GeoVaLs::multiplyInverseAD not implemented");
}
// -------------------------------------------------------------------------------------------------
void LinVarChaModel2GeoVaLs::print(std::ostream & os) const {
  os << classname() << " variable change";
}

}  // namespace wrf_hydro_nwm_jedi
