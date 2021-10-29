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

  LinearGetValues::LinearGetValues(
      const Geometry & geom,
      const ufo::Locations & locs,
      const eckit::Configuration &) :
    locs_(locs), geom_(new Geometry(geom)) {
    wrf_hydro_nwm_jedi_getvalues_create_f90(
        keyGetValues_, geom.toFortran(), locs_);
    oops::Log::trace() << "LinearGetValues::LinearGetValues done" << std::endl;
  }


  LinearGetValues::~LinearGetValues() {
    std::cout << "LinearGetValues destructor invoked but not implemented" << std::endl;
    // util::abor1_cpp("LinearGetValues::~LinearGetValues() needs to be implemented.",
    //              __FILE__, __LINE__);
  }


  void LinearGetValues::setTrajectory(const State & state, const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      ufo::GeoVaLs & geovals)
  {
    // wrf_hydro_nwm_jedi_lineargetvalues_set_trajectory_f90(keyLinearGetValues_,
    //                                                    geom_->toFortran(),
    //                                                    state.toFortran(),
    //                                                    &t1p,
    //                                                    &t2p,
    //                                                    locs_.toFortran(),
    //                                                    geovals.toFortran());
    wrf_hydro_nwm_jedi_getvalues_fill_geovals_f90(keyGetValues_,
                                                  state.geometry()->toFortran(),
                                                  state.toFortran(),
                                                  t1, t2,
                                                  locs_,
                                                  geovals.toFortran());
    // util::abor1_cpp(
    //     "LinearGetValues::setTrajectory() needs to be implemented.", __FILE__, __LINE__);
  }


  void LinearGetValues::fillGeoVaLsTL(const Increment & inc,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      ufo::GeoVaLs & geovals) const {
    oops::Log::trace() << "LinearGetValues::fillGeovalsTL starting"
                       << std::endl;

    std::cout << "Before invoking lineargetvalues" << std::endl;

    wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_f90(keyGetValues_,
                                                        geom_->toFortran(),
                                                        inc.toFortran(),
                                                        t1, t2,
                                                        locs_,
                                                        geovals.toFortran());
    oops::Log::trace() << "LinearGetValues::fillGeovalsTL done" << std::endl;
  }


  void LinearGetValues::fillGeoVaLsAD(Increment & inc,
                                      const util::DateTime & t1,
                                      const util::DateTime & t2,
                                      const ufo::GeoVaLs & geovals) const {
    oops::Log::trace() << "LinearGetValues::fillGeovalsAD starting"
                     << std::endl;

    // sw_lineargetvalues_fill_geovals_ad_f90(keyLinearGetValues_,
    //                                        geom_->toFortran(),
    //                                        inc.toFortran(),
    //                                        &t1p, &t2p,
    //                                        locs_.toFortran(),
    //                                        geovals.toFortran());

    wrf_hydro_nwm_jedi_lineargetvalues_fill_geovals_ad_f90(keyGetValues_,
                                                           geom_->toFortran(),
                                                           inc.toFortran(),
                                                           t1, t2,
                                                           locs_,
                                                           geovals.toFortran());

    oops::Log::trace() << "LinearGetValues::fillGeovalsAD done" << std::endl;
  }


  void LinearGetValues::print(std::ostream & os) const {
    util::abor1_cpp("LinearGetValues::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "(TODO, print diagnostic info about the LinearGetValues here)"
       << std::endl;
  }


}  // namespace wrf_hydro_nwm_jedi
