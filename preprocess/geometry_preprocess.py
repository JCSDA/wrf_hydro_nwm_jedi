import pathlib
import wrfhydropy
import xarray as xr

# Testing purposes:
domain_dir = pathlib.Path(
    '/Users/jamesmcc/Downloads/croton_NY')
hrldas_patch_filename = 'hrldas_namelist_patches.json'
hydro_patch_filename = 'hrldas_namelist_patches.json'
config = 'nwm_ana'

# This relies on wrfhydropy, make arguments to simply specify the individual files by hand.
# def geometry_preprocess(
#     domain_dir: pathlib.Path,
#     config: str,
#     hrldas_patch_filename = 'hrldas_namelist_patches.json',
#     hydro_patch_filename = 'hydro_namelist_patches.json'):

hrldas_patch_file = domain_dir / hrldas_patch_filename
hydro_patch_file = domain_dir / hydro_patch_filename

hrldas_patches = wrfhydropy.namelist.JSONNamelist(hrldas_patch_file)
hydro_patches = wrfhydropy.namelist.JSONNamelist(hydro_patch_file)

# HRLDAS
hrldas_namelist = hrldas_patches.get_config(config)
wrfinput_file = pathlib.Path(
    hrldas_namelist['noahlsm_offline']['hrldas_setup_file'])
if not wrfinput_file.exists():
    wrfinput_file = domain_dir / wrfinput_file
if not wrfinput_file.exists():
    raise FileExistsError('wrfinput file not found: ' + wrfinput_file)
  
wrfinput_ds = xr.open_dataset(wrfinput_file)

hrldas_drop_list = [
    'CANWAT', 'DZS', 'HGT', 'ISLTYP', 'IVGTYP',
    'LAI', 'MAPFAC_MX', 'MAPFAC_MY', 'SEAICE', 'SHDMAX',
    'SHDMIN', 'SMOIS', 'SNOW', 'TMN', 'TSK',
     'TSLB', 'XLAND']
hrldas_geom = wrfinput_ds.drop_vars(hrldas_drop_list).squeeze('Time')

hrldas_keep_attrs = ['DX', 'DY']
hrldas_attrs = {
    key: val for key,val in wrfinput_ds.attrs.items() if key in hrldas_keep_attrs }

hrldas_attrs['lsm_lat_dim_name'] = 'south_north'
hrldas_attrs['lsm_lon_dim_name'] = 'west_east'
hrldas_attrs['lsm_z_dim_name'] = 'soil_layers_stag'

hrldas_attrs['lsm_lat_name'] = 'XLAT'
hrldas_attrs['lsm_lon_name'] = 'XLON'
hrldas_attrs['lsm_z_name'] = 'ZS'

hrldas_geom.attrs = hrldas_attrs
