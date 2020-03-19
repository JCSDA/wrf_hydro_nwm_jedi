import pathlib
import wrfhydropy

# Testing purposes:
domain_dir = pathlib.Path(
    '/Users/james/Downloads/croton/example_case')
#    '/Users/james/Downloads/croton_NY')
hrldas_patch_filename = 'hrldas_namelist_patches.json'
hydro_patch_filename = 'hrldas_namelist_patches.json'
config = 'nwm_ana'

# def geometry_preprocess(
#     domain_dir: pathlib.Path,
#     config: str,
#     hrldas_patch_filename = 'hrldas_namelist_patches.json',
#     hydro_patch_filename = 'hydro_namelist_patches.json'):

hrldas_patch_file = domain_dir / hrldas_patch_filename
hydro_patch_file = domain_dir / hydro_patch_filename

hrldas_patches = wrfhydropy.namelist.JSONNamelist(hrldas_patch_file)
hydro_patches = wrfhydropy.namelist.JSONNamelist(hydro_patch_file)

hrldas_namelist = hrldas_patches.get_config(config)


adsfj

    

    
