/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

// #include "Fortran.h"
#include "wrf_hydro_nwm_jedi/Traits.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/LinearModel/TlmId.h"

namespace wrf_hydro_nwm_jedi {

// -----------------------------------------------------------------------------
  static oops::LinearModelMaker<Traits, TlmId>
  makerIdTLM_("WRF_HYDRO_NWM_IdTLM");
// -----------------------------------------------------------------------------
  TlmId::TlmId(const Geometry & resol,
	       const eckit::Configuration & tlConf)
    : keyConfig_(0), tstep_(), resol_(resol), linvars_(tlConf)
  {
    tstep_ = util::Duration(tlConf.getString("tstep"));
    oops::Log::trace() << "TlmId created" << std::endl;
  }
// -----------------------------------------------------------------------------
  TlmId::~TlmId() {
    oops::Log::trace() << "TlmIdSW destructed" << std::endl;
  }
  // -----------------------------------------------------------------------------
  void TlmId::setTrajectory(const State &, State &,
			    const ModelAuxControl &) {std::cout << "In trajectory" << std::endl;}
// -----------------------------------------------------------------------------
  void TlmId::initializeTL(Increment & dx) const {
    oops::Log::debug() << "TlmId::initializTL" << std::endl;
  }
// -----------------------------------------------------------------------------
  void TlmId::stepTL(Increment & dx,
		     const ModelAuxIncrement &) const {
    std::cout << "StepTL" << std::endl;
    dx.updateTime(tstep_);
  }
// -----------------------------------------------------------------------------
  void TlmId::finalizeTL(Increment & dx) const {
    oops::Log::debug() << "TlmId::finalizeTL" << std::endl;
  }
// -----------------------------------------------------------------------------
void TlmId::initializeAD(Increment & dx) const {
  oops::Log::debug() << "TlmId::initializAD" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::stepAD(Increment & dx,
		   ModelAuxIncrement &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeAD(Increment & dx) const {
  oops::Log::debug() << "TlmId::finalizeAD" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::print(std::ostream & os) const {
  os << "IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace sw
