/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/GetValues/GetValues.interface.h"
#include "wrf_hydro_nwm_jedi/State/State.h"

// -------------------------------------------------------------------------------------------------

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace oops {
  class Variables;
}

namespace wrf_hydro_nwm_jedi {
  class State;
  class Geometry;
  class VarChaModel2GeoVaLs;

// -------------------------------------------------------------------------------------------------

  class GetValues : public util::Printable, private util::ObjectCounter<GetValues> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm_jedi::GetValues";}

    GetValues(const Geometry &, const ufo::Locations & locs,
              const eckit::Configuration & config);
    virtual ~GetValues();

    void fillGeoVaLs(const State &, const util::DateTime &, const util::DateTime &,
                     ufo::GeoVaLs &) const;

   private:
    void print(std::ostream &) const;
    F90getvalues keyGetValues_;
    ufo::Locations locs_;
    std::shared_ptr<const Geometry> geom_;
    std::unique_ptr<VarChaModel2GeoVaLs> model2geovals_;
  };

  // -------------------------------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
