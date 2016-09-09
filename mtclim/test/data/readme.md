# MTCLIM test data sets

### Daymet Data
1.  Seattle, WA location

        Latitude: 47.6519  Longitude: -122.3046
        X & Y on Lambert Conformal Conic: -1583161.66 759257.95
        Tile: 12269
        Elevation: 30 meters
        1980-2013
2.  Telluride, CO location

        Latitude: 37.9146  Longitude: -107.8456
        X & Y on Lambert Conformal Conic: -659717.37 -455261.18
        Tile: 11377
        Elevation: 3252 meters
        1980-2013

### VIC Stehekin Test Dataset

1.  `vic42_daily` - raw binary VIC forcings at a daily timestep.  The forcings are described below:

        FORCE_FORMAT	BINARY	# BINARY or ASCII
        FORCE_ENDIAN	LITTLE	# LITTLE (PC/Linux) or BIG (SUN)
        N_TYPES		4	# Number of variables (columns)
        FORCE_TYPE	PREC	UNSIGNED	40
        FORCE_TYPE	TMAX	SIGNED	100
        FORCE_TYPE	TMIN	SIGNED	100
        FORCE_TYPE	WIND	SIGNED	100
        FORCE_DT	24	# Forcing time step length (hours)
        FORCEYEAR	1949	# Year of first forcing record
        FORCEMONTH	01	# Month of first forcing record
        FORCEDAY	01	# Day of first forcing record
        FORCEHOUR	00	# Hour of first forcing record
        GRID_DECIMAL	4	# Number of digits after decimal point in forcing file names
        WIND_H          10.0    # height of wind speed measurement (m)
        MEASURE_H       2.0     # height of humidity measurement (m)
        ALMA_INPUT	FALSE	# TRUE = ALMA-compliant input variable units; FALSE = standard VIC units

2.  `vic42_disagg` - VIC/MTCLIM disaggregated forcings developed from `vic42_daily` using VIC.4.2.b and `OUTPUTFORCE = TRUE`.  The outputs were developed using default settings.
