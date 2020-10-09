/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "wrf_hydro_nwm_jedi/Traits.h"

#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/HofX.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateModelFactory<wrf_hydro_nwm_jedi::Traits>();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::HofX<wrf_hydro_nwm_jedi::Traits, ufo::ObsTraits> hofx;
  return run.execute(hofx);
}
