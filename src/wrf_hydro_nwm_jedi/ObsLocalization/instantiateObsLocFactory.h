/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/ObsLocalizationBase.h"
#include "ufo/instantiateObsLocFactory.h"
#include "ufo/ObsTraits.h"
#include "wrf_hydro_nwm_jedi/ObsLocalization/ObsLocVerticalBrasnett.h"
#include "wrf_hydro_nwm_jedi/Traits.h"

namespace wrf_hydro_nwm_jedi {
void instantiateObsLocFactory() {
  ufo::instantiateObsLocFactory<wrf_hydro_nwm_jedi::Traits>();
  static oops::ObsLocalizationMaker<wrf_hydro_nwm_jedi::Traits, ufo::ObsTraits,
                                    wrf_hydro_nwm_jedi::ObsLocVerticalBrasnett>
         makerVerticalBrasnett_("Vertical Brasnett");
}

}  // namespace wrf_hydro_nwm_jedi
