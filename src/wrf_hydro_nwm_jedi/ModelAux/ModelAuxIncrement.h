/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_MODELAUX_MODELAUXINCREMENT_H_
#define WRF_HYDRO_NWM_JEDI_MODELAUX_MODELAUXINCREMENT_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

// forward declarations
namespace eckit {
  class Configuration;
}
namespace wrf_hydro_nwm_jedi {
  class ModelAuxControl;
  class ModelAuxCovariance;
  class Geometry;
}

namespace wrf_hydro_nwm_jedi {
  class ModelAuxIncrement :
    public util::Printable,
    // private util::ObjectCounter<ModelAuxIncrement>,
    public util::Serializable {
   public:
    static const std::string classname() {return "wrf_hydro_nwm_jedi::ModelAuxIncrement";}

    ModelAuxIncrement(const Geometry &, const eckit::Configuration &) {}
    ModelAuxIncrement(const ModelAuxIncrement &, const bool) {}
    ModelAuxIncrement(const ModelAuxIncrement &, const eckit::Configuration &) {}
    ~ModelAuxIncrement() {}

    // Linear algebra operators
    void diff(const ModelAuxControl &, const ModelAuxControl &) {}
    void zero() {}
    ModelAuxIncrement & operator=(const ModelAuxIncrement &) {return *this;}
    ModelAuxIncrement & operator+=(const ModelAuxIncrement &) {return *this;}
    ModelAuxIncrement & operator-=(const ModelAuxIncrement &) {return *this;}
    ModelAuxIncrement & operator*=(const double) {return *this;}

    void axpy(const double, const ModelAuxIncrement &) {}
    // JLM TODO: zero dot product?
    double dot_product_with(const ModelAuxIncrement &) const {return 0.0;}

    void read(const eckit::Configuration&) {}
    void write(const eckit::Configuration&) {}

    double norm() const {return 0.0;}

    /// Serialize and deserialize
    size_t serialSize() const  override {return 0;}
    void serialize(std::vector<double> &) const override {}
    void deserialize(const std::vector<double> &, size_t &) override {}

   private:
    explicit ModelAuxIncrement(const ModelAuxCovariance &);
    void print(std::ostream & os) const override {}
  };
}  // namespace wrf_hydro_nwm_jedi
#endif  // WRF_HYDRO_NWM_JEDI_MODELAUX_MODELAUXINCREMENT_H_
