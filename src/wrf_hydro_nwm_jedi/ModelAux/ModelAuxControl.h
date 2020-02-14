/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM-JEDI_MODELAUX_MODELAUXCONTROL_H_
#define WRF_HYDRO_NWM-JEDI_MODELAUX_MODELAUXCONTROL_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace wrf_hydro_nwm-jedi {
  class Geometry;
}

//-----------------------------------------------------------------------------

namespace wrf_hydro_nwm-jedi {

  // ModelAuxControl class
  class ModelAuxControl : public util::Printable,
                          private util::ObjectCounter<ModelAuxControl> {
   public:
    static const std::string classname() {return "wrf_hydro_nwm-jedi::ModelAuxControl";}

    ModelAuxControl(const Geometry &, const eckit::Configuration &);
    ModelAuxControl(const Geometry &, const ModelAuxControl &);
    ModelAuxControl(const ModelAuxControl &, const bool);
    ~ModelAuxControl();

   private:
    void print(std::ostream & os) const;
  };
}  // namespace wrf_hydro_nwm-jedi
#endif  // WRF_HYDRO_NWM-JEDI_MODELAUX_MODELAUXCONTROL_H_
