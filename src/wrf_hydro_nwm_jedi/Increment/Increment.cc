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
#include "wrf_hydro_nwm_jedi/State/State.h"

//#include "oops/base/GridPoint.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "eckit/exception/Exceptions.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

using oops::Log;


namespace wrf_hydro_nwm_jedi {

  Increment::Increment(const Geometry & geom,
                       const oops::Variables & vars,
                       const util::DateTime & time)
    : fields_(new Fields(geom, vars)),
      vars_(vars),
      time_(time)
  {
    Log::trace() << "Increment::Increment 1 constructor " << std::endl;
    wrf_hydro_nwm_jedi_increment_create_f90(
	keyInc_,
	fields_->geometry()->toFortran(),
	vars_);
    Log::trace() << "Increment::Increment 1 constructor END " << std::endl;
  }


  Increment::Increment(const Increment & other,
		       const bool copy)
    : fields_(new Fields(*other.fields_->geometry(), other.vars_)),
      vars_(other.vars_),
      time_(other.time_) {
    Log::trace() << "Increment::Increment 2 constructor " << std::endl;
    wrf_hydro_nwm_jedi_increment_create_from_other_f90(keyInc_, other.keyInc_);
    if(copy) {
      wrf_hydro_nwm_jedi_increment_copy_f90(keyInc_, other.keyInc_);
    } else {
      wrf_hydro_nwm_jedi_increment_zero_f90(keyInc_);
    }
    // this->print(std::cout);
    Log::trace() << "Increment::Increment 2 constructor END" << std::endl;
  }


  Increment::Increment(const Geometry &,
		       Increment & other)
    : time_(other.time_),
      vars_(other.vars_),
      fields_(new Fields(*other.fields_->geometry(), other.vars_))
  {
    Log::trace() << "Increment::Increment 3 constructed from other" << std::endl;
    // wrf_hydro_nwm_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
    wrf_hydro_nwm_jedi_increment_create_from_other_f90(keyInc_, other.keyInc_);
    //wrf_hydro_nwm_increment_change_resol_f90(toFortran(), other.keyFlds_);
    Log::trace() << "Increment::Increment 3 constructed from other END" << std::endl;
  }

  Increment::~Increment() {
    //util::abor1_cpp("Increment::~Increment() needs to be implemented.",
    //                __FILE__, __LINE__);
  }


  // ----------------------------------------------------------------------------

  Increment & Increment::operator =(const Increment &rhs) {
    wrf_hydro_nwm_jedi_increment_copy_f90(keyInc_, rhs.keyInc_);
    time_ = rhs.time_;
    return *this;
  }


  Increment & Increment::operator -=(const Increment &dx) {
    ASSERT(this->validTime() == dx.validTime());
    wrf_hydro_nwm_jedi_increment_sub_f90(keyInc_, dx.keyInc_);
    return *this;
  }


  Increment & Increment::operator +=(const Increment &dx) {
    ASSERT(this->validTime() == dx.validTime());
    wrf_hydro_nwm_jedi_increment_add_f90(keyInc_, dx.keyInc_);
    return *this;
  }


  Increment & Increment::operator *=(const double &zz) {
    wrf_hydro_nwm_jedi_increment_mul_f90(keyInc_, static_cast<float>(zz));
    return *this;
  }


  void Increment::axpy(const double &scalar,
		       const Increment &other,
		       const bool check) {
    ASSERT(this->validTime() == other.validTime());
    wrf_hydro_nwm_jedi_increment_axpy_f90(
	keyInc_, static_cast<float>(scalar), other.keyInc_);
  }


  double Increment::dot_product_with(const Increment & other) const {
    double result = 0.0;
    wrf_hydro_nwm_jedi_increment_dot_prod_f90(keyInc_, other.keyInc_, result);
    std::cout << "dot product result: " << result << std::endl;
    return result;
  }


  void Increment::zero() {
    wrf_hydro_nwm_jedi_increment_zero_f90(keyInc_);
  }


  void Increment::zero(const util::DateTime & time) {
    wrf_hydro_nwm_jedi_increment_zero_f90(keyInc_);
    time_ = time;
  }


  void Increment::ones() {
    wrf_hydro_nwm_jedi_increment_ones_f90(keyInc_);
  }


  void Increment::diff(const State & x1, const State & x2) {
    ASSERT(this->validTime() == x1.validTime());
    ASSERT(this->validTime() == x2.validTime());
    wrf_hydro_nwm_jedi_increment_diff_incr_f90(
	keyInc_, x1.toFortran(), x2.toFortran());
  }


  void Increment::schur_product_with(const Increment & dx) {
    wrf_hydro_nwm_jedi_increment_schur_f90(toFortran(), dx.toFortran());
  }


  void Increment::random() {
    wrf_hydro_nwm_jedi_increment_random_f90(toFortran());
  }


  double Increment::norm() const {
    double norm = 0.0;
    norm = wrf_hydro_nwm_jedi_increment_rms_f90(toFortran());
    return norm;
  }


  void Increment::getValuesTL(const ufo::Locations &,
                              const oops::Variables &,
                              ufo::GeoVaLs &,
                              const GetValuesTraj &) const {
    util::abor1_cpp("Increment::getValuesTL() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  void Increment::getValuesAD(const ufo::Locations &,
                              const oops::Variables &,
                              const ufo::GeoVaLs &,
                              const GetValuesTraj &) {
    util::abor1_cpp("Increment::getValuesAD() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  void Increment::print(std::ostream & os) const {
    char *string = new char[8192];
    wrf_hydro_nwm_jedi_increment_print_f90(keyInc_, string);
    os << "Print Increment (C++) -------------------- ";
    os << string;
    os << "End Print Increment (C++) ---------------- " << std::endl  << std::endl;
    delete[] string;
  }


  // const util::DateTime & time() const {return time_;}
  // util::DateTime & time() {return time_;}


  const util::DateTime & Increment::validTime() const { return time_;  }
  util::DateTime & Increment::validTime() { return time_; }


  void Increment::updateTime(const util::Duration & dt) {
    // fields_->time() += dt;}
    time_ += dt;}


  std::shared_ptr<const Geometry> Increment::geometry() const {return fields_->geometry();}


  void Increment::read(const eckit::Configuration & conf) {
    util::abor1_cpp("Increment::read() needs to be implemented.",
                    __FILE__, __LINE__);
  }

 void Increment::write(const eckit::Configuration & conf) {
    util::abor1_cpp("Increment::write() needs to be implemented.",
                    __FILE__, __LINE__);
  }

  void Increment::ug_coord(oops::UnstructuredGrid &) const {
    util::abor1_cpp("Increment::ug_coord() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  void Increment::field_to_ug(oops::UnstructuredGrid &, const int &) const {
    util::abor1_cpp("Increment::field_to_ug() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  void Increment::field_from_ug(const oops::UnstructuredGrid &, const int &) {
    util::abor1_cpp("Increment::field_from_ug() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  // oops::GridPoint Increment::getPoint(const GeometryIterator &) const {
  //   util::abor1_cpp("Increment::getPoint() needs to be implemented.",
  //                   __FILE__, __LINE__);

  //   oops::Variables vars;
  //   std::vector<double> vals;
  //   std::vector<int> varlens;
  //   return oops::GridPoint(vars, vals, varlens);
  // }


  void Increment::setPoint(const oops::GridPoint &, const GeometryIterator &) {
    util::abor1_cpp("Increment::setPoint() needs to be implemented.",
                    __FILE__, __LINE__);
  }


  void Increment::accumul(const double & zz, const State & xx) {
    oops::Log::trace() << "Increment accumul starting" << std::endl;
    wrf_hydro_nwm_jedi_increment_accumul_f90(keyInc_, zz, xx.toFortran());
    oops::Log::trace() << "Increment accumul END" << std::endl;
  }


}  // namespace wrf_hydro_nwm_jedi
