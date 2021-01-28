/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_GETVALUES_LINEARGETVALUES_H_
#define WRF_HYDRO_NWM_JEDI_GETVALUES_LINEARGETVALUES_H_

#include <memory>
#include <ostream>

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

#include "ufo/Locations.h"

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/GetValues/GetValues.interface.h"
#include "wrf_hydro_nwm_jedi/GetValues/LinearGetValues.interface.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h"

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class State;
  class Increment;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // LinearGetValues class
  class LinearGetValues : public util::Printable {
   public:
    // LinearGetValues();
    LinearGetValues(const Geometry &, const ufo::Locations &);
    virtual ~LinearGetValues();
    void setTrajectory(const State & state,
                       const util::DateTime & t1,
                       const util::DateTime & t2,
                       ufo::GeoVaLs & geovals);
    void fillGeoVaLsTL(const Increment & inc,
                       const util::DateTime & t1,
                       const util::DateTime & t2,
                       ufo::GeoVaLs & geovals) const;
    void fillGeoVaLsAD(Increment & inc,
                       const util::DateTime & t1,
                       const util::DateTime & t2,
                       const ufo::GeoVaLs & geovals) const;
   private:
    void print(std::ostream & os) const;
    F90getvalues keyGetValues_;
    ufo::Locations locs_;
    std::shared_ptr<const Geometry> geom_;
    std::unique_ptr<VarChaModel2GeoVaLs> model2geovals_;
  };

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_GETVALUES_LINEARGETVALUES_H_
