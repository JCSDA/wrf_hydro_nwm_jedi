import pathlib
import xarray as xr

# Create a restart that has slightly different SWE
in_dir = pathlib.Path(
    '/Users/james/jedi/bundle-nwm/wrf_hydro_nwm_jedi/test/Data/wrf_hydro_nwm_files')
in_file = in_dir / 'RESTART.2017010100_DOMAIN1'
out_file = in_dir / 'RESTART.2017010100_sneqv_plus5_DOMAIN1'
ds = xr.open_dataset(in_file)
ds['SNEQV'] = ds['SNEQV'] + 5.0
ds.to_netcdf(out_file)
