/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_LINEARMODEL_H_
#define WRF_HYDRO_NWM_JEDI_LINEARMODEL_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/LinearModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "wrf_hydro_nwm_jedi/Traits.h"

// Forward declarations
namespace wrf_hydro_nwm_jedi {
  class Geometry;
  class State;
  class Increment;
  struct Traits;
}

namespace eckit {
  class Configuration;
}

namespace wrf_hydro_nwm_jedi {
  
// -----------------------------------------------------------------------------
/// SW linear identity model definition.
/*!
 *  SW linear identity model definition and configuration parameters.
 */

  class TlmId: public oops::LinearModelBase<Traits>,
    private util::ObjectCounter<TlmId> {
  public:
      static const std::string classname() {return "wrf_hydro_nwm_jedi::TlmId";}

      TlmId(const Geometry &, const eckit::Configuration &);
      ~TlmId();

/// Model trajectory computation
      void setTrajectory(const State &, State &,
			 const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(Increment &) const override;
  void stepTL(Increment &, const ModelBiasIncrement &)
               const override;
  void finalizeTL(Increment &) const override;

  void initializeAD(Increment &) const override;
  void stepAD(Increment &, ModelBiasIncrement &)
                const override;
  void finalizeAD(Increment &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const Geometry & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;

// Data
  int keyConfig_;
  util::Duration tstep_;
  const Geometry resol_;
  const oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace sw
#endif  // JEDI_SRC_TLMID_TLMIDSW_H_
