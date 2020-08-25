/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/runs/MakeObs.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"
#include "wrf_hydro_nwm_jedi/Traits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::MakeObs<wrf_hydro_nwm_jedi::Traits, ufo::ObsTraits> mkobs;
  return run.execute(mkobs);
}
