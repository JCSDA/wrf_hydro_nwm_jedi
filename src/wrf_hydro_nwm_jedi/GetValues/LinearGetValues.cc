/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "wrf_hydro_nwm_jedi/GetValues/LinearGetValues.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace wrf_hydro_nwm_jedi {
  
  // ----------------------------------------------------------------------------
  
  LinearGetValues::LinearGetValues() {
    util::abor1_cpp("LinearGetValues::LinearGetValues() needs to be implemented.",
                    __FILE__, __LINE__);
  }
  
// ----------------------------------------------------------------------------

  LinearGetValues::LinearGetValues(const Geometry &, const ufo::Locations &) {
    util::abor1_cpp("LinearGetValues::LinearGetValues(const Geometry &, const ufo::Locations &) needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  LinearGetValues::~LinearGetValues() {
    util::abor1_cpp("LinearGetValues::~LinearGetValues() needs to be implemented.",
		    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------
  void LinearGetValues::setTrajectory(const State & state, const util::DateTime & t1,
				      const util::DateTime & t2,
				      ufo::GeoVaLs & geovals)
  {
    util::abor1_cpp("LinearGetValues::setTrajectory() needs to be implemented.", __FILE__, __LINE__);
  }
  
  // ----------------------------------------------------------------------------

  void LinearGetValues::fillGeoVaLsTL(const Increment & inc,
				      const util::DateTime & t1,
				      const util::DateTime & t2,
				      ufo::GeoVaLs & geovals) const {
    oops::Log::trace() << "LinearGetValuesSW::fillGeovalsTL starting"
		       << std::endl;
    
    // const util::DateTime * t1p = &t1;
    // const util::DateTime * t2p = &t2;
    
  // sw_lineargetvalues_fill_geovals_tl_f90(keyLinearGetValues_,
  //                                        geom_->toFortran(),
  //                                        inc.toFortran(),
  //                                        &t1p, &t2p,
  //                                        locs_.toFortran(),
  //                                        geovals.toFortran());
  // oops::Log::trace() << "LinearGetValuesSW::fillGeovalsTL done" << std::endl;
  }
  
  // -------------------------------------------------------------------------------------------------
  
  void LinearGetValues::fillGeoVaLsAD(Increment & inc,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
				      const ufo::GeoVaLs & geovals) const {
    oops::Log::trace() << "LinearGetValuesSW::fillGeovalsAD starting"
                     << std::endl;
    
    // const util::DateTime * t1p = &t1;
    // const util::DateTime * t2p = &t2;
    
    // sw_lineargetvalues_fill_geovals_ad_f90(keyLinearGetValues_,
    //                                        geom_->toFortran(),
    //                                        inc.toFortran(),
    //                                        &t1p, &t2p,
    //                                        locs_.toFortran(),
    //                                        geovals.toFortran());
    
    // oops::Log::trace() << "LinearGetValuesSW::fillGeovalsAD done"
    //                    << std::endl;
  }
  
  // -------------------------------------------------------------------------------------------------

  void LinearGetValues::print(std::ostream & os) const {
    util::abor1_cpp("LinearGetValues::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the LinearGetValues here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
