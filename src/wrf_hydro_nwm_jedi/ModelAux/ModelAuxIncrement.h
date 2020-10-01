/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_MODELAUX_INCREMENT_H_
#define WRF_HYDRO_NWM_JEDI_MODELAUX_INCREMENT_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class ModelAuxControl;
}

//-----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // ModelAuxIncrement class
  class ModelAuxIncrement :
    public util::Printable,
    private util::ObjectCounter<ModelAuxIncrement>,
    public util::Serializable {
    
   public:
    static const std::string classname() {return "wrf_hydro_nwm_jedi::ModelAuxIncrement";}

    ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &);
    ModelAuxIncrement(const ModelAuxIncrement &, const bool);
    ModelAuxIncrement(const Geometry &, const eckit::Configuration &);
    ~ModelAuxIncrement();

    // Linear algebra operators
    void diff(const ModelAuxControl &, const ModelAuxControl &) {}
    void zero();
    ModelAuxIncrement & operator*=(const double);
    ModelAuxIncrement & operator+=(const ModelAuxIncrement &);
    ModelAuxIncrement & operator-=(const ModelAuxIncrement &);
    double norm() const;
    void axpy(const double, const ModelAuxIncrement &);
    double dot_product_with(const ModelAuxIncrement &) const;

    /// Serialize and deserialize
    size_t serialSize() const {return 0;}
    void serialize(std::vector<double> &) const {}
    void deserialize(const std::vector<double> &, size_t &) {}
    
   private:
    void print(std::ostream &) const;
  };
}  // namespace wrf_hydro_nwm_jedi
#endif  // WRF_HYDRO_NWM-JEDI_MODELAUX_INCREMENT_H_
