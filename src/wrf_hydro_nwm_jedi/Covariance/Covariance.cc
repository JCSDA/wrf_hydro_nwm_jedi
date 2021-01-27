/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "wrf_hydro_nwm_jedi/Covariance/Covariance.h"
#include "wrf_hydro_nwm_jedi/Covariance/CovarianceFortran.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/State/State.h"

#include "eckit/config/Configuration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

namespace wrf_hydro_nwm_jedi {

// ----------------------------------------------------------------------------

  Covariance::Covariance(const Geometry & geom,     // unused but required by the interface.
                         const oops::Variables & vars,
                         const eckit::Configuration & conf,
                         const State & bkg,
                         const State & traj) {
    // time_ = util::DateTime(conf.getString("date"));
    const eckit::Configuration * configc = &conf;
    vars_ = oops::Variables(conf, "analysis variables");
    wrf_hydro_nwm_jedi_b_setup_f90(
        keyFtnConfig_, bkg.toFortran(), &configc, vars_);
    oops::Log::trace() << "Covariance created" << std::endl;
  }

// ----------------------------------------------------------------------------

  Covariance::~Covariance() {
    //  util::abor1_cpp("Covariance::~Covariance() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Covariance::multiply(const Increment & dxin, Increment & dxout) const {
    wrf_hydro_nwm_jedi_b_mult_f90(keyFtnConfig_, dxin.toFortran(),
                                  dxout.toFortran());
  }

// ----------------------------------------------------------------------------

  void Covariance::inverseMultiply(const Increment & dxin,
                                   Increment & dxout) const {
    dxout = dxin;
  }

// ----------------------------------------------------------------------------

  void Covariance::randomize(Increment & dx) const {
    wrf_hydro_nwm_jedi_b_randomize_f90(keyFtnConfig_, dx.toFortran());
  }

// ----------------------------------------------------------------------------

  void Covariance::print(std::ostream & os) const {
    util::abor1_cpp("Covariance::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the Covariance here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
