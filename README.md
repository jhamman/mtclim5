MT-CLIM
======
**Mountain Climate Simulator**

## Summary

MT-CLIM is a computer program that uses observations of daily maximum temperature, minimum temperature, and precipitation from one location (the "base") to estimate the temperature, precipitation, radiation, and humidity at another location (the "site"). The base and the site can be at different elevations, and can have different slopes and aspects. Better results are obtained when the base and site are relatively close to one another (at the bottom and the top of a valley, for example).

Temperature estimates at the site are based on the base temperatures and a user-supplied temperature lapse rate. Separate lapse rates can be supplied for daily maximum and minimum temperature. Precipitation estimates at the site are based on the daily record of preciupitation from the base, and a user-specified ratio of annual total precipitation between the site and the base.

The estimation of radiation and humidity are more complex, since these parameters are not assumed to be measured at the base station. Humidity estimates are based on the observation that daily minimum temperature is usually very close to dewpoint temperature. The MT-CLIM algorithm includes corrections to this assumption for arid climates. Radiation estimates are based on the observation that the diurnal temperature range (from minimum temperature to maximum temperature) is closely related to the daily average atmospheric transmittance. In conjunction with information about the latitude, elevation, slope, and aspect of the site, this relationship can be used to estimate daily total radiation with a typical error range of +/- 15%. Some of the limitations of MT-CLIM are the use of a single base station for observations and the need for the user to specify the temperature and precipitation relationships with elevation.

## Development History

The development of MT-CLIM was originally developed by Peter Thornton at the Numerical Terradynamic Simulation Group (NTSG) at the University of Montana.  MT-CLIM has been used as a run-time meteorological forcing disaggregator inside the [Variable Infiltration Capacity (VIC) model](https://github.com/UW-Hydro/VIC).

## Version 5

Development of MT-CLIM version 5.0 began with the version of MT-CLIM version 4.3 embedded in the [Variable Infiltration Capacity (VIC) model](https://github.com/UW-Hydro/VIC) version 4.2.  From that point, the C source code was extracted from VIC and converted to [Python](https://www.python.org/).  Version 5.0 will include the following new features:

- [Python](https://www.python.org/)/[Pandas](http://pandas.pydata.org/) interface.
- Algorithmic flexibility (resampling methods, model choices, and parameter space).
- Subaily estimation or temporal disaggregation of all variables.
- Improved solar geometry functions.
- Improved estimation of tranmissivity by calculating the mean cos(solar zenith angle) and integrating that quantity when disaggregating shortwave to sub-daily time steps and calculating transmissivity.
- ASCII text and netCDF input/output options.

## Publications

* Bohn, T. J., B. Livneh, J. W. Oyler, S. W. Running, B. Nijssen, and D. P. Lettenmaier, 2013a: Global evaluation of MT-CLIM and related algorithms for forcing of ecological and hydrological models, Agr. Forest. Meteorol., 176, 38-49, doi:10.1016/j.agrformet.2013.03.003.

* Bristow, K.L., and G.S. Campbell, 1984. On the relationship between incoming solar radiation and daily maximum and minimum temperature. Agricultural and Forest Meteorology, 31:159-166.

* Running, S.W., R.R. Nemani, and R.D. Hungerford, 1987. Extrapolation of synoptic meteorological data in mountainous terrain and its use for simulating forest evaporation and photosynthesis. Canadian Journal of Forest Research, 17:472-483.

* Glassy, J.M., and S.W. Running, 1994. Validating diurnal climatology of the MT-CLIM model across a climatic gradient in Oregon. Ecological Applications, 4(2):248-257.

* Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for estimating surface humidity from daily minimum temperature. Agricultural and Forest Meteorology, 85:87-98.

* Thornton, P.E., and S.W. Running, 1999. An improved algorithm for estimating incident daily solar radiation from measurements of temperature, humidity, and precipitation. Agricultural and Forest Meteorology, 93:211-228.
