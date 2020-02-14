/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef WRF_HYDRO_NWM-JEDI_INCREMENT_INCREMENT_H_
#define WRF_HYDRO_NWM-JEDI_INCREMENT_INCREMENT_H_

#include <memory>
#include <ostream>

#include <boost/shared_ptr.hpp>

#include "oops/util/DateTime.h"
#include "oops/util/Printable.h"

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
namespace wrf_hydro_nwm-jedi {
  class Fields;
  class Geometry;
  class GetValuesTraj;
  class State;
}

// ----------------------------------------------------------------------------

namespace wrf_hydro_nwm-jedi {

  // Increment class
  class Increment : public util::Printable {
   public:
    // Constructor, destructor
    Increment(const Geometry &, const oops::Variables &,
              const util::DateTime &);
    Increment(const Increment &, const bool);
    ~Increment();

    void read(const eckit::Configuration &);
    double norm() const;
    void random();

    Increment & operator =(const Increment &);
    Increment & operator-=(const Increment &);
    Increment & operator+=(const Increment &);
    Increment & operator*=(const double &);
    void axpy(const double &, const Increment &, const bool check = true);
    double dot_product_with(const Increment &) const;
    void zero();
    void diff(const State &, const State &);
    void schur_product_with(const Increment &);

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
    oops::GridPoint getPoint(const GeometryIterator &) const;
    void setPoint(const oops::GridPoint &, const GeometryIterator &);


    boost::shared_ptr<const Geometry> geometry() const;

   private:
    void print(std::ostream &) const;
    std::unique_ptr<Fields> fields_;
  };
}  // namespace wrf_hydro_nwm-jedi

#endif  // WRF_HYDRO_NWM-JEDI_INCREMENT_INCREMENT_H_
