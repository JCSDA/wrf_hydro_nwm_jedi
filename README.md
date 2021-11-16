# WRF-Hydro/NWM JEDI Implementation

Develop Branch: [![ Build Status](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoibkViSGpiTjh4TDhMVnhnMU1SbFlYOEQvVTNYb1E0cnh3Qi9FeHNYRFVkY3hLYzhBOFJmT3ZaUE9oWmxkMjB6K0ZsREpORUZSUEVRdE91NEtVeGZNWHBBPSIsIml2UGFyYW1ldGVyU3BlYyI6ImpJNmpwNjRCQzNyeGJBaFYiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://console.aws.amazon.com/codesuite/codebuild/469205354006/projects/wrf_hydro_nwm_jedi-internal-clang)


----
##  Overview

Primarily aimed at snow stated updating of NoahMP within WRF-Hydro using JEDI.


## Description of Sub Directories

### preprocess/
Currently preprocess focuses on gathering geometry information
from disparate WRF-Hydro/NWM input files and putting them into
a single file.
