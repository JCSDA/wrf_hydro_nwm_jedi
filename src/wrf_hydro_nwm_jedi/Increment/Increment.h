/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef WRF_HYDRO_NWM_JEDI_INCREMENT_INCREMENT_H_
#define WRF_HYDRO_NWM_JEDI_INCREMENT_INCREMENT_H_

#include <memory>
#include <ostream>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "atlas/field.h"

#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/IncrementFortran.h"

// forward declarations
namespace oops {
  class GridPoint;
  class UnstructuredGrid;
  class Variables;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace wrf_hydro_nwm_jedi {
  class Fields;
  class Geometry;
  class GetValuesTraj;
  class State;
  typedef int F90inc;
}


// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm_jedi {

  // Increment class
  class Increment :
    public util::Printable,
    public util::Serializable
    // public oops::GeneralizedDepartures,
    // public util::Printable,
    // private util::ObjectCounter<IncrementSW> {
  {
   public:
    // Constructor, destructor
    Increment(const Geometry &,
              const oops::Variables &,
              const util::DateTime &);
    Increment(const Geometry &, Increment &);
    Increment(const Increment &, const bool);
    ~Increment();

    void read(const eckit::Configuration &);
    void write(const eckit::Configuration &);
    double norm() const;

    void diff(const State &, const State &);
    void zero();
    void zero(const util::DateTime &);
    void ones();
    Increment & operator =(const Increment &);
    Increment & operator-=(const Increment &);
    Increment & operator+=(const Increment &);
    Increment & operator*=(const double &);
    void axpy(const double &, const Increment &, const bool check = true);
    double dot_product_with(const Increment &) const;
    void schur_product_with(const Increment &);
    void random();
    void dirac(const eckit::Configuration &);

    // ATLAS
    void setAtlas(atlas::FieldSet *) const;
    void toAtlas(atlas::FieldSet *) const;
    void fromAtlas(atlas::FieldSet *);

    // Interpolate increment to observation location
    void getValuesTL(const ufo::Locations &,
                     const oops::Variables &,
                     ufo::GeoVaLs &,
                     const GetValuesTraj &) const;

    void getValuesAD(const ufo::Locations &,
                     const oops::Variables &,
                     const ufo::GeoVaLs &,
                     const GetValuesTraj &);

    // time manipulation
    const util::DateTime & validTime() const;
    util::DateTime & validTime();
    void updateTime(const util::Duration &);

    // unstructured grid conversions
    void ug_coord(oops::UnstructuredGrid &) const;
    void field_to_ug(oops::UnstructuredGrid &, const int &) const;
    void field_from_ug(const oops::UnstructuredGrid &, const int &);

    // Iterator access
    /* oops::GridPoint getPoint(const GeometryIterator &) const; */
    void setPoint(const oops::GridPoint &, const GeometryIterator &);

    /// Serialize and deserialize
    size_t serialSize() const {return 0;}
    void serialize(std::vector<double> &) const {}
    void deserialize(const std::vector<double> &, size_t &) {}

    /// Other / utils
    void accumul(const double &, const State &);

    std::shared_ptr<const Geometry> geometry() const;

    F90inc & toFortran() {return keyInc_;}
    const F90inc & toFortran() const {return keyInc_;}

   private:
    oops::Variables vars_;
    util::DateTime time_;
    F90inc keyInc_;
    void print(std::ostream &) const;
    std::unique_ptr<Fields> fields_;
  };

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_INCREMENT_INCREMENT_H_
