Preprocessing
=============


## Geometry 

```
(368) jamesmcc@relleno[549]:~/jedi/repos/wrf_hydro_nwm_jedi/preprocess>  ./geometry_preprocess.py \
  --domain_dir /Users/jamesmcc/Downloads/croton_NY \
  --config nwm_ana \
  --hrldas_patch_filename hrldas_namelist_patches.json \
  --hydro_patch_filename hrldas_namelist_patches.json 

Writing geometry file for config 'nwm_ana' to:
/Users/jamesmcc/Downloads/croton_NY/NWM/DOMAIN/geometry_nwm_ana.nc

(368) jamesmcc@relleno[550]:~/jedi/repos/wrf_hydro_nwm_jedi/preprocess> ncdump -h /Users/jamesmcc/Downloads/croton_NY/NWM/DO
MAIN/geometry_nwm_ana.nc
netcdf geometry_nwm_ana {
dimensions:
	Time = UNLIMITED ; // (0 currently)
	south_north = 16 ;
	west_east = 15 ;
	soil_layers_stag = 4 ;
variables:
	float XLAT(south_north, west_east) ;
		XLAT:_FillValue = NaNf ;
		XLAT:FieldType = 104 ;
		XLAT:MemoryOrder = "XY " ;
		XLAT:units = "degrees latitude" ;
		XLAT:description = "Latitude on mass grid" ;
		XLAT:stagger = "M" ;
		XLAT:sr_x = 1 ;
		XLAT:sr_y = 1 ;
	float XLONG(south_north, west_east) ;
		XLONG:_FillValue = NaNf ;
		XLONG:FieldType = 104 ;
		XLONG:MemoryOrder = "XY " ;
		XLONG:units = "degrees longitude" ;
		XLONG:description = "Longitude on mass grid" ;
		XLONG:stagger = "M" ;
		XLONG:sr_x = 1 ;
		XLONG:sr_y = 1 ;
	float ZS(soil_layers_stag) ;
		ZS:_FillValue = -1.e+36f ;
		ZS:units = "m" ;

// global attributes:
		:DX = 1000.f ;
		:DY = 1000.f ;
		:lsm_lat_dim_name = "south_north" ;
		:lsm_lon_dim_name = "west_east" ;
		:lsm_z_dim_name = "soil_layers_stag" ;
		:lsm_lat_name = "XLAT" ;
		:lsm_lon_name = "XLON" ;
		:lsm_z_name = "ZS" ;
}
```
