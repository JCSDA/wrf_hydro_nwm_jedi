/*
 * (C) Copyright 2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/interface/LinearVariableChange.h"
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
static oops::LinearVariableChangeMaker<Traits,
       oops::LinearVariableChange<Traits, LinVarChaModel2GeoVaLs> >
  makerLinVarChaModel2GeoVaLs_("Model2GeoVaLs");

static oops::LinearVariableChangeMaker<Traits,
       oops::LinearVariableChange<Traits, LinVarChaModel2GeoVaLs> >
  makerLinVarChaModel2GeoDef_("default");

// -------------------------------------------------------------------------------------------------

LinVarChaModel2GeoVaLs::LinVarChaModel2GeoVaLs(const State & bg, const State &fg,
                                        const Geometry &geom,
                                        const eckit::Configuration &conf)
  : geom_(new Geometry(geom)) {
}

// -----------------------------------------------------------------------------

LinVarChaModel2GeoVaLs::~LinVarChaModel2GeoVaLs() {
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::multiply(const Increment &dxin,
                                            Increment &dxout) const {
  wrf_hydro_nwm_jedi_model2geovals_linear_changevar_Ident_f90(geom_->toFortran(),
                                                          dxin.toFortran(),
                                                          dxout.toFortran());
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::multiplyInverse(const Increment &,
                                                Increment &) const {
  util::abor1_cpp("LinVarChaModel2GeoVaLs::multiplyInverse not implemented");
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::multiplyAD(const Increment &dxin,
                                              Increment &dxout) const {
  wrf_hydro_nwm_jedi_model2geovals_linear_changevar_Ident_f90(geom_->toFortran(),
                                                          dxin.toFortran(),
                                                          dxout.toFortran());
}

// -----------------------------------------------------------------------------

void LinVarChaModel2GeoVaLs::multiplyInverseAD(const Increment &,
                                                  Increment &) const {
  util::abor1_cpp("LinVarChaModel2GeoVaLs::multiplyInverseAD not implemented");
}

}  // namespace wrf_hydro_nwm_jedi
