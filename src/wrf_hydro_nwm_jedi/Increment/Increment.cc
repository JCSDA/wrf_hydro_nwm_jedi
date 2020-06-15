/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "wrf_hydro_nwm_jedi/Fields/Fields.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"

//#include "oops/base/GridPoint.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

namespace wrf_hydro_nwm_jedi {

// ----------------------------------------------------------------------------

  Increment::Increment(const Geometry & geom,
                       const oops::Variables & vars,
                       const util::DateTime & vt) {
    // util::abor1_cpp("Increment::Increment() needs to be implemented.",
    //                 __FILE__, __LINE__);
    wrf_hydro_nwm_jedi_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  }

// ----------------------------------------------------------------------------

  Increment::Increment(const Increment & other, const bool copy) {
    // util::abor1_cpp("Increment::Increment() needs to be implemented.",
    //                  __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  Increment::~Increment() {
    // util::abor1_cpp("Increment::~Increment() needs to be implemented.",
    //                 __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator =(const Increment &) {
    // util::abor1_cpp("Increment::operator= needs to be implemented.",
    //                 __FILE__, __LINE__);
    // Needed by test
    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator -=(const Increment &) {
    // util::abor1_cpp("Increment::operator-= needs to be implemented.",
    //                 __FILE__, __LINE__);
    // Needed by test
    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator +=(const Increment &) {
    // util::abor1_cpp("Increment::operator+= needs to be implemented.",
    //                 __FILE__, __LINE__);
    // Needed by test
    return *this;
  }

// ----------------------------------------------------------------------------

  Increment & Increment::operator *=(const double &) {
    // util::abor1_cpp("Increment::operator*= needs to be implemented.",
    //                 __FILE__, __LINE__);
    // Needed by test
    return *this;
  }

// ----------------------------------------------------------------------------

  void Increment::axpy(const double &, const Increment &, const bool check) {
    // util::abor1_cpp("Increment::axpy() needs to be implemented.",
    //                 __FILE__, __LINE__);
    // Needed by test
  }

// ----------------------------------------------------------------------------

  double Increment::dot_product_with(const Increment &) const {
    // util::abor1_cpp("Increment::dot_product_with() needs to be implemented.",
    //                 __FILE__, __LINE__);
    // Needed by test
    return 0.0;
  }

// ----------------------------------------------------------------------------

  void Increment::zero() {
    util::abor1_cpp("Increment::zero() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::diff(const State & x1, const State & x2) {
    util::abor1_cpp("Increment::diff() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::schur_product_with(const Increment & ) {
    util::abor1_cpp("Increment::schur_product_with() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::random() {
    // util::abor1_cpp("Increment::random() needs to be implemented.",
    //                 __FILE__, __LINE__);
    //Needed by the test
  }

// ----------------------------------------------------------------------------

  double Increment::norm() const {
    // util::abor1_cpp("Increment::norm() needs to be implemented.",
    //                 __FILE__, __LINE__);
    // return fields_->norm();
    //Needed by the test
    return 0.0;
  }

// ----------------------------------------------------------------------------

  void Increment::getValuesTL(const ufo::Locations &,
                              const oops::Variables &,
                              ufo::GeoVaLs &,
                              const GetValuesTraj &) const {
    util::abor1_cpp("Increment::getValuesTL() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::getValuesAD(const ufo::Locations &,
                              const oops::Variables &,
                              const ufo::GeoVaLs &,
                              const GetValuesTraj &) {
    util::abor1_cpp("Increment::getValuesAD() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::print(std::ostream & os) const {
    util::abor1_cpp("Increment::print() needs to be implemented.",
                    __FILE__, __LINE__);
    os << "Increment: "
       << "(TODO, print diagnostic info about the increment here)"
       << std::endl;
  }

// ----------------------------------------------------------------------------

  const util::DateTime & Increment::validTime() const {
    return fields_->time();
  }

// ----------------------------------------------------------------------------

  util::DateTime & Increment::validTime() {
    return fields_->time();
  }

// ----------------------------------------------------------------------------

  void Increment::updateTime(const util::Duration & dt) {
    fields_->time() += dt;
  }

// ----------------------------------------------------------------------------

  boost::shared_ptr<const Geometry> Increment::geometry() const {
    return fields_->geometry();
  }

// ----------------------------------------------------------------------------

  void Increment::read(const eckit::Configuration & conf) {
    util::abor1_cpp("Increment::read() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::ug_coord(oops::UnstructuredGrid &) const {
    util::abor1_cpp("Increment::ug_coord() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::field_to_ug(oops::UnstructuredGrid &, const int &) const {
    util::abor1_cpp("Increment::field_to_ug() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  void Increment::field_from_ug(const oops::UnstructuredGrid &, const int &) {
    util::abor1_cpp("Increment::field_from_ug() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

  // oops::GridPoint Increment::getPoint(const GeometryIterator &) const {
  //   util::abor1_cpp("Increment::getPoint() needs to be implemented.",
  //                   __FILE__, __LINE__);

  //   oops::Variables vars;
  //   std::vector<double> vals;
  //   std::vector<int> varlens;
  //   return oops::GridPoint(vars, vals, varlens);
  // }

// ----------------------------------------------------------------------------

  void Increment::setPoint(const oops::GridPoint &, const GeometryIterator &) {
    util::abor1_cpp("Increment::setPoint() needs to be implemented.",
                    __FILE__, __LINE__);
  }

// ----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
