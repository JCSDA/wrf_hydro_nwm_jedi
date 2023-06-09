/*
 * (C) Copyright 2019-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/runs/LocalEnsembleDA.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsErrorFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/instantiateObsLocFactory.h"
#include "ufo/ObsTraits.h"
#include "wrf_hydro_nwm_jedi/ObsLocalization/instantiateObsLocFactory.h"
#include "wrf_hydro_nwm_jedi/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsLocFactory<wrf_hydro_nwm_jedi::Traits>();
  ufo::instantiateObsErrorFactory();
  ufo::instantiateObsFilterFactory();
  wrf_hydro_nwm_jedi::instantiateObsLocFactory();
  oops::LocalEnsembleDA<wrf_hydro_nwm_jedi::Traits, ufo::ObsTraits> letkf;
  return run.execute(letkf);
}
