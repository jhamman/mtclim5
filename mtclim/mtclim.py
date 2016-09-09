'''
mtclim: Mountain Climate Simulator
'''

# Mountain Climate Simulator, meteorological forcing disaggregator
# Copyright (C) 2015  Joe Hamman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
from warnings import warn
from statsmodels.tools.eval_measures import rmse
from .physics import svp, calc_pet, atm_pres
from .share import default_parameters, default_options, constants


class MtClim(object):
    '''The University of Montana Mountain Climate Simulator'''

    def __init__(self, run=False, data=None, parameters=None, options=None):
        '''
        Initialize MtClim object.

        Parameters
        ----------
        data : pandas.DataFrame, optional
            Input data: pandas DataFrame with at least `tmax`, `tmin`, and
            `prcp`. Timestep freq should be `D`.
        run : bool, optional.
            Compute all values imediately (default: False)
        parameters : dict-like, optional
            Dictionary of parameters to use apart from the default parameters.
        options : dict-like, optional
            Dictionary of options to use apart from the default options.
        '''

        # Set model parameters
        self.params = default_parameters
        if parameters is not None:
            self.params.update(parameters)

        # Set model options
        self.options = default_options
        if options is not None:
            self.options.update(options)

        # Initialize data attributes

        if data is not None:
            self.data = data
        else:
            self.data = pd.DataFrame()

        if run:
            self.init()
            self.run()

    def init(self):
        self.tinystepspday = 86400 / constants['SRADDT']
        self.tiny_radfract = np.zeros(shape=(366, self.tinystepspday),
                                      dtype=np.float64)

    def run(self):
        self.calc_tair()
        self.calc_prcp()
        self.snowpack()
        self.calc_srad_humidity_iterative()
        self.calc_longwave()

    def resample(self):
        return

    @property
    def data(self):
        '''The objects DataFrame'''
        return self._data

    @data.setter
    def data(self, df):
        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                'data must be a Pandas DataFrame instance, got %s' % type(df))
        if not all([v in df for v in ['tmax', 'tmin', 'prcp']]):
            raise ValueError('data must include tmax, tmin, and prcp')

        self._data = df.resample('D').mean()
        self.ndays = len(self._data)

    def __repr__(self):
        r = 'MtClim object\nparameters: {0}\noptions: {1}\ndata: {2}'.format(
            self.params, self.options, self.data.head())
        return r

    def __str__(self):
        return 'MtClim object'

    def calc_tair(self):
        '''
        Calculates daily air temperatures.
        '''
        # calculate elevation difference in meters
        dz = (self.params['site_elev'] -
              self.params['base_elev'])

        # apply lapse rate corrections to tmax and tmin
        self.data['s_tmax'] = self.data['tmax'] + \
            (dz * self.params['tmax_lr'])
        self.data['s_tmin'] = self.data['tmin'] + \
            (dz * self.params['tmin_lr'])

        # Since tmax lapse rate usually has a larger absolute value than tmin
        # lapse rate, it is possible at high elevation sites for these
        # corrections to result in tmin > tmax. Check for that occurrence and
        # force tmin = corrected tmax - 0.5 deg C.
        self.data['s_tmin'].where(self.data['s_tmin'] > self.data['s_tmax'],
                                  other=self.data['s_tmax'] - 0.5,
                                  inplace=True)

        # derived temperatures
        tmean = self.data[['s_tmax', 's_tmin']].mean(axis=1)
        self.data['s_tday'] = ((self.data['s_tmax'] - tmean) *
                               self.params['TDAYCOEF']) + tmean

    def calc_prcp(self):
        '''
        Calculates daily total precipitation
        '''
        if (self.params['site_isoh'] is not None and
                self.params['base_isoh'] is not None):
            factor = self.params['site_isoh'] / self.params['base_isoh']
            self.data['s_prcp'] = self.data['prcp'] * factor
        else:
            self.data['s_prcp'] = self.data['prcp']

    def snowpack(self):
        '''
        estimates the accumulation and melt of snow for radiation algorithm
        corrections
        '''
        # initialize SWE array
        self.data['s_swe'] = 0.

        # first pass
        self._simple_snowpack(0.)

        # use the first pass to set the initial snowpack conditions for the
        # first day of data
        start_yday = self.data.index.dayofyear[0]
        prev_yday = (self.data.index[0] - pd.Timedelta(1, unit='D')).dayofyear
        count = 0
        swe_sum = 0.
        for i in range(self.ndays):
            if (self.data.index.dayofyear[i] == start_yday or
                    self.data.index.dayofyear[i] == prev_yday):
                count += 1
                swe_sum += self.data['s_swe'][i]

        # Proceed with correction if there are valid days to reinitialize
        # the snowpack estiamtes. Otherwise use the first-pass estimate.
        if count:
            snowpack = swe_sum / count
            self._simple_snowpack(snowpack)

    def _simple_snowpack(self, snowpack):
        '''
        '''
        for i in range(self.ndays):
            if (self.data['s_tmin'][i] <= self.params['SNOW_TCRIT']):
                snowpack += self.data['s_prcp'][i]
            else:
                snowpack -= self.params['SNOW_TRATE'] * \
                    (self.data['s_tmin'][i] - self.params['SNOW_TCRIT'])
            snowpack = np.maximum(snowpack, 0.)
            self.data['s_swe'][i] = snowpack

    def calc_srad_humidity_iterative(self, tol=0.01, win_type='boxcar'):
        '''
        Iterative estimation of shortwave radiation and humidity

        TODO: simplify
        '''
        ndays = self.ndays

        daylength = np.zeros(366)
        window = np.zeros(ndays + 90)

        ttmax0 = np.zeros(366)
        flat_potrad = np.zeros(366)
        slope_potrad = np.zeros(366)
        t_fmax = np.zeros(ndays)
        self.data['s_tfmax'] = 0.

        # calculate diurnal temperature range for transmittance calculations
        self.data['tmax'] = np.maximum(self.data['tmax'], self.data['tmin'])
        dtr = self.data['tmax'] - self.data['tmin']

        # smooth dtr array: After Bristow and Campbell, 1984
        # use 30-day antecedent smoothing window
        sm_dtr = pd.rolling_window(dtr, window=30, freq='D',
                                   win_type=win_type).fillna(method='bfill')
        if self.ndays <= 30:
            warn('Timeseries is shorter than rolling mean window, filling '
                 'missing values with unsmoothed data.')
            sm_dtr.fillna(dtr, inplace=True)

        # calculate the annual total precip
        sum_prcp = self.data['s_prcp'].values.sum()

        ann_prcp = (sum_prcp / self.ndays) * 365.25
        if (ann_prcp == 0.):
            ann_prcp = 1.0

        # Generate the effective annual precip, based on a 3-month
        # moving-window. Requires some special case handling for the
        # beginning of the record and for short records.

        # check if there are at least 90 days in this input file, if not,
        # use a simple total scaled to effective annual precip
        if (ndays < 90):
            sum_prcp = self.data['s_prcp'].values.sum()

            effann_prcp = (sum_prcp / self.ndays) * 365.25

            # if the effective annual precip for this period
            # is less than 8 cm, set the effective annual precip to 8 cm
            # to reflect an arid condition, while avoiding possible
            # division-by-zero errors and very large ratios (PET/Pann)
            effann_prcp = np.maximum(effann_prcp, 8.)

            parray = effann_prcp

        else:
            # Check if the yeardays at beginning and the end of this input file
            # match up. If so, use parts of the three months at the end
            # of the input file to generate effective annual precip for
            # the first 3-months. Otherwise, duplicate the first 90 days
            # of the record.
            start_yday = self.data.index.dayofyear[0]
            end_yday = self.data.index.dayofyear[ndays - 1]
            if (start_yday != 1):
                if end_yday == start_yday - 1:
                    isloop = True
            else:
                if end_yday == 365 or end_yday == 366:
                    isloop = True

            # fill the first 90 days of window
            for i in range(90):
                if (isloop):
                    window[i] = self.data['s_prcp'][ndays - 90 + i]
                else:
                    window[i] = self.data['s_prcp'][i]
            # fill the rest of the window array
            window[90:] = self.data['s_prcp']

            # for each day, calculate the effective annual precip from
            # scaled 90-day total
            for i in range(self.ndays):
                sum_prcp = 0.
                for j in range(90):
                    sum_prcp += window[i + j]
                    sum_prcp = (sum_prcp / 90.) * 365.25

            # if the effective annual precip for this 90-day period
            # is less than 8 cm, set the effective annual precip to 8 cm
            # to reflect an arid condition, while avoiding possible
            # division-by-zero errors and very large ratios (PET/Pann)

        # start of the main radiation algorithm

        # before starting the iterative algorithm between humidity and
        # radiation, calculate all the variables that don't depend on
        # humidity so they only get done once.

        trans1 = self._calc_trans()

        # STEP (3) build 366-day array of ttmax0, potential rad, and daylength
        # precalculate the transcendentals
        # check for (+/-) 90 degrees latitude, throws off daylength calc
        lat = np.clip(self.params['site_lat'] * constants['RADPERDEG'],
                      -np.pi / 2., np.pi / 2.)
        coslat = np.cos(lat)
        sinlat = np.sin(lat)
        cosslp = np.cos(self.params['site_slope'] * constants['RADPERDEG'])
        sinslp = np.sin(self.params['site_slope'] * constants['RADPERDEG'])
        cosasp = np.cos(self.params['site_aspect'] *
                        constants['RADPERDEG'])
        sinasp = np.sin(self.params['site_aspect'] *
                        constants['RADPERDEG'])
        # cosine of zenith angle for east and west horizons
        coszeh = np.cos(np.pi / 2. - (self.params['site_east_horiz'] *
                                      constants['RADPERDEG']))
        coszwh = np.cos(np.pi / 2. - (self.params['site_west_horiz'] *
                                      constants['RADPERDEG']))

        # sub-daily time and angular increment information
        dt = constants['SRADDT']  # set timestep
        dh = dt / constants['SECPERRAD']  # calculate hour-angle step

        # begin loop through yeardays
        for i in range(365):
            # calculate cos and sin of declination
            decl = constants['MINDECL'] * np.cos((i + constants['DAYSOFF']) *
                                                 constants['RADPERDAY'])
            cosdecl = np.cos(decl)
            sindecl = np.sin(decl)

            # do some precalculations for beam-slope geometry (bsg)
            bsg1 = -sinslp * sinasp * cosdecl
            bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl
            bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl

            # calculate daylength as a function of lat and decl
            cosegeom = coslat * cosdecl
            sinegeom = sinlat * sindecl
            coshss = np.clip(-sinegeom / cosegeom, -1, 1)

            hss = np.cos(coshss)  # hour angle at sunset (radians)
            # daylength (seconds)
            daylength[i] = np.maximum(2.0 * hss * constants['SECPERRAD'],
                                      86400)

            # solar constant as a function of yearday (W/m^2)
            sc = 1368.0 + 45.5 * np.sin((2.0 * np.pi * i / 365.25) + 1.7)

            # extraterrestrial radiation perpendicular to beam, total over
            # the timestep (J)
            dir_beam_topa = sc * dt

            sum_trans = 0.
            sum_flat_potrad = 0.
            sum_slope_potrad = 0.

            # begin sub-daily hour-angle loop, from -hss to hss
            for h in np.arange(-hss, hss, dh):
                # precalculate cos and sin of hour angle
                cosh = np.cos(h)
                sinh = np.sin(h)

                # calculate cosine of solar zenith angle
                cza = cosegeom * cosh + sinegeom

                # calculate cosine of beam-slope angle
                cbsa = sinh * bsg1 + cosh * bsg2 + bsg3

                # check if sun is above a flat horizon
                if (cza > 0.):
                    # when sun is above the ideal (flat) horizon, do all the
                    # flat-surface calculations to determine daily total
                    # transmittance, and save flat-surface potential radiation
                    # for later calculations of diffuse radiation

                    # potential radiation for this time period, flat surface,
                    # top of atmosphere
                    dir_flat_topa = dir_beam_topa * cza

                    # determine optical air mass
                    am = 1.0 / (cza + 0.0000001)
                    if (am > 2.9):
                        ami = int((np.cos(cza) / constants['RADPERDEG'])) - 69
                        if (ami < 0):
                            ami = 0
                        if (ami > 20):
                            ami = 20
                        am = constants['OPTAM'][ami]

                    # correct instantaneous transmittance for this optical
                    # air mass
                    trans2 = np.power(trans1, am)

                    # instantaneous transmittance is weighted by potential
                    # radiation for flat surface at top of atmosphere to get
                    # daily total transmittance
                    sum_trans += trans2 * dir_flat_topa

                    # keep track of total potential radiation on a flat
                    # surface for ideal horizons
                    sum_flat_potrad += dir_flat_topa

                    # keep track of whether this time step contributes to
                    # component 1 (direct on slope)
                    if ((h < 0. and cza > coszeh and cbsa > 0.) or
                            (h >= 0. and cza > coszwh and cbsa > 0.)):
                        # sun between east and west horizons, and direct on
                        # slope. this period contributes to component 1
                        sum_slope_potrad += dir_beam_topa * cbsa
                else:
                    dir_flat_topa = -1

                tinystep = np.clip(((12 * 3600 + h * constants['SECPERRAD']) /
                                    constants['SRADDT']),
                                   0, self.tinystepspday - 1)

                if dir_flat_topa > 0:
                    self.tiny_radfract[i, tinystep] = dir_flat_topa
                else:
                    self.tiny_radfract[i, tinystep] = 0

            if daylength[i] and sum_flat_potrad > 0:
                self.tiny_radfract[i] /= sum_flat_potrad

            # calculate maximum daily total transmittance and daylight average
            # flux density for a flat surface and the slope
            if daylength[i]:
                ttmax0[i] = sum_trans / sum_flat_potrad
                flat_potrad[i] = sum_flat_potrad / daylength[i]
                slope_potrad[i] = sum_slope_potrad / daylength[i]
            else:
                ttmax0[i] = 0.
                flat_potrad[i] = 0.
                slope_potrad[i] = 0.

        # force yearday 366 = yearday 365
        ttmax0[365] = ttmax0[364]
        flat_potrad[365] = flat_potrad[364]
        slope_potrad[365] = slope_potrad[364]
        daylength[365] = daylength[364]

        self.tiny_radfract[365] = self.tiny_radfract[364]

        # STEP (4)  calculate the sky proportion for diffuse radiation

        # uses the product of spherical cap defined by average horizon angle
        # and the great-circle truncation of a hemisphere. this factor does not
        # vary by yearday.
        avg_horizon = (self.params['site_east_horiz'] +
                       self.params['site_west_horiz']) / 2.0
        horizon_scalar = 1.0 - np.sin(avg_horizon * constants['RADPERDEG'])
        if (self.params['site_slope'] > avg_horizon):
            slope_excess = self.params['site_slope'] - avg_horizon
        else:
            slope_excess = 0.
        if (2.0 * avg_horizon > 180.):
            slope_scalar = 0.
        else:
            slope_scalar = np.clip(
                1.0 - (slope_excess / (180.0 - 2.0 * avg_horizon)), 0, None)

        sky_prop = horizon_scalar * slope_scalar

        # b parameter, and t_fmax not varying with Tdew, so these can be
        # calculated once, outside the iteration between radiation and humidity
        # estimates. Requires storing t_fmax in an array.
        # b parameter from 30-day average of DTR
        b = self.params['B0'] + self.params['B1'] * \
            np.exp(-self.params['B2'] * sm_dtr)

        # proportion of daily maximum transmittance
        t_fmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, self.params['C']))

        # correct for precipitation if this is a rain day
        inds = np.nonzero(self.data['prcp'] >
                          self.options['SW_PREC_THRESH'])[0]
        t_fmax[inds] *= self.params['RAIN_SCALAR']
        self.data['s_tfmax'] = t_fmax

        # Initial values of vapor pressure, etc
        if 'tdew' in self.data:
            # Observed Tdew supplied
            tdew = self.data['tdew']
        else:
            # Estimate Tdew
            tdew = self.data['s_tmin']
        if 's_hum' in self.data:
            # Observed vapor pressure supplied
            pva = self.data['s_hum']
        else:
            # convert dewpoint to vapor pressure
            pva = svp(tdew)

        # Other values needed for srad_humidity calculation
        pa = atm_pres(self.params['site_elev'])
        yday = self.data.index.dayofyear - 1
        self.data['s_dayl'] = daylength[yday]
        tdew_save = tdew
        pva_save = pva

        # Initial estimates of solar radiation, cloud fraction, etc.
        tdew, pva, pet = self._compute_srad_humidity_onetime(
            tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength,
            parray, pa, dtr)

        # estimate annual PET
        sum_pet = pet.values.sum()
        ann_pet = (sum_pet / self.ndays) * 365.25

        # Reset humidity terms if no iteration desired
        if (('tdew' in self.data) or ('s_hum' in self.data) or
                (self.options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
                 ann_pet / ann_prcp >= 2.5)):
                tdew = tdew_save[:]
                pva = pva_save[:]

        # Set up srad-humidity iterations
        if (self.options['VP_ITER'].upper() == 'VP_ITER_ALWAYS' or
            (self.options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
             ann_pet / ann_prcp >= 2.5) or
                self.options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
            if (self.options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
                max_iter = 100
            else:
                max_iter = 2
        else:
            max_iter = 1

        # srad-humidity iterations
        # iter_i = 1
        rmse_tdew = tol + 1
        # while (rmse_tdew > tol and iter_i < max_iter):
        #     tdew_save = tdew[:]
        #
        #     tdew = self._compute_srad_humidity_onetime(
        #         tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop,
        #         daylength, parray, pa, dtr)
        #
        #     rmse_tdew = rmse(tdew, tdew_save)
        #     iter_i += 1

        def f(tdew, *args):
            rmse_tdew = rmse(self._compute_srad_humidity_onetime(tdew, *args),
                             tdew)
            return rmse_tdew

        res = minimize_scalar(f, tdew, args=(pva, ttmax0, flat_potrad,
                                             slope_potrad, sky_prop, daylength,
                                             parray, pa, dtr),
                              tol=rmse_tdew, options={'maxiter': max_iter})
        tdew = res.x

        pva = svp(tdew)

        # save humidity in output data structure
        if 's_hum' not in self.data:
            self.data['s_hum'] = pva
        # output humidity as vapor pressure deficit (Pa)
        # calculate saturated VP at tday
        pvs = svp(self.data['s_tday'])
        vpd = pvs - pva
        self.data['s_vpd'] = np.maximum(vpd, 0.)

    def _compute_srad_humidity_onetime(self, tdew, pva, ttmax0, flat_potrad,
                                       slope_potrad, sky_prop, daylength,
                                       parray, pa, dtr):
        '''
        Initial estimates of solar radiation, cloud fraction, etc.

        Parameters
        ----------
        tdew : pandas.Series
            description
        pva : pandas.Series
            description
        ttmax0 : pandas.Series
            description
        flat_potrad : pandas.Series
            description
        slope_potrad : pandas.Series
            description
        sky_prop : pandas.Series
            description
        daylength : pandas.Series
            description
        parray : pandas.Series
            description
        pa : pandas.Series
            description
        dtr : pandas.Series
            description

        Returns
        ----------
        tdew : pandas.Series
            description
        pva : pandas.Series
            description
        pet : pandas.Series
            description
        '''

        yday = self.data.index.dayofyear - 1

        # Compute SW radiation
        t_tmax = ttmax0[yday] + (self.params['ABASE'] * pva)

        # this is mainly for the case of observed VP supplied, for
        # which t_tmax sometimes ends up being negative (when potential
        # radiation is low and VP is high)
        t_tmax = np.minimum(t_tmax, 0.0001)

        self.data['s_ttmax'] = t_tmax

        # final daily total transmittance
        t_final = t_tmax * self.data['s_tfmax']

        # estimate fraction of radiation that is diffuse, on an
        # instantaneous basis, from relationship with daily total
        # transmittance in Jones (Plants and Microclimate, 1992)
        # Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
        # Fig 6.14, p. 122.
        pdif = np.clip(-1.25 * t_final + 1.25, 0., 1.)

        # estimate fraction of radiation that is direct, on an
        # instantaneous basis
        pdir = 1.0 - pdif

        # the daily total radiation is estimated as the sum of the
        # following two components:
        # 1. The direct radiation arriving during the part of
        # the day when there is direct beam on the slope.
        # 2. The diffuse radiation arriving over the entire daylength
        # (when sun is above ideal horizon).

        # component 1
        srad1 = slope_potrad[yday] * t_final * pdir

        # component 2 (diffuse)
        # includes the effect of surface albedo in raising the diffuse
        # radiation for obstructed horizons
        srad2 = (flat_potrad[yday] * t_final * pdif) * \
            (sky_prop + self.params['DIF_ALB'] * (1.0 - sky_prop))

        # snow pack influence on radiation
        sc = np.zeros_like(self.data['s_swe'])
        if (self.options['MTCLIM_SWE_CORR']):
            inds = np.nonzero(self.data['s_swe'] > 0. * daylength[yday] > 0.)
            # convert to W/m2 and check for zero daylength
            # snow correction in J/m2/day
            sc[inds] = (1.32 + 0.096 * self.data['s_swe'][inds]) *\
                1.0e6 / daylength[yday][inds]

            # set a maximum correction of 100 W/m2
            sc = np.maximum(sc, 100.)  # JJH - this is fishy
            # this could also be sc[inds] = np.maximum(sc[inds], 100.) optimize

        # save daily radiation
        # save cloud transmittance when rad is an input
        if 's_swrad' in self.data:
            potrad = (srad1 + srad2 + sc) * daylength[yday] / t_final / 86400
            self.data['s_tfmax'] = 1.0
            inds = np.nonzero((potrad > 0.) * (self.data['s_swrad'] > 0.) *
                              (daylength[yday] > 0))[0]
            # both of these are 24hr mean rad. here
            self.data['s_tfmax'][inds] = (self.data['s_swrad'][inds] /
                                          (potrad[inds] * t_tmax[inds]))
            self.data['s_tfmax'] = np.maximum(self.data['s_tfmax'], 1.)
        else:
            self.data['s_swrad'] = srad1 + srad2 + sc

        if (self.options['LW_CLOUD'].upper() == 'LW_CLOUD_DEARDORFF'):
            self.data['s_tskc'] = (1. - self.data['s_tfmax'])
        else:
            self.data['s_tskc'] = np.sqrt((1. - self.data['s_tfmax']) / 0.65)
        self.data['s_fdir'] = pdir

        # Compute PET using SW radiation estimate, and update Tdew, pva **
        tmink = self.data['s_tmin'] + constants['KELVIN']
        pet = calc_pet(self.data['s_swrad'], self.data['s_tday'], pa,
                       self.data['s_dayl'])

        # calculate ratio (PET/effann_prcp) and correct the dewpoint
        ratio = pet / parray
        self.data['s_ppratio'] = ratio * 365.25
        tdewk = tmink * (-0.127 + 1.121 *
                         (1.003 - 1.444 * ratio + 12.312 *
                          np.power(ratio, 2) - 32.766 * np.power(ratio, 3)) +
                         0.0006 * dtr)
        tdew = tdewk - constants['KELVIN']

        return tdew

    def calc_longwave(self):
        '''This routine estimates long wave radiation based on the fractional
        cloud cover (self.data['s_tskc']), the current air temperature (C), and
        the atmospheric vapor pressure (Pa).
        '''

        # See Bras, R. F. , "Hydrology, an introduction to hydrologic science",
        # Addison-Wesley, 1990, p. 42-45

        # convert to Kelvin
        air_temp = self.data['s_tday'] + constants['KELVIN']
        # convert to mbar
        vp = self.data['s_vp'] / 100

        if (self.options['LW_TYPE'].upper() == 'TVA'):
            # TVA (1972) - see Bras 2.35
            emissivity_clear = 0.740 + 0.0049 * vp
        elif (self.options['LW_TYPE'].upper() == 'ANDERSON'):
            # Anderson (1964)
            emissivity_clear = 0.68 + 0.036 * np.power(vp, 0.5)
        elif (self.options['LW_TYPE'].upper() == 'BRUTSAERT'):
            # Brutsaert (1975)
            x = vp / air_temp
            emissivity_clear = 1.24 * np.power(x, 0.14285714)
        elif (self.options['LW_TYPE'].upper() == 'SATTERLUND'):
            # Satterlund (1979)
            emissivity_clear = 1.08 * \
                (1 - np.exp(-1 * np.power(vp, (air_temp / 2016))))
        elif (self.options['LW_TYPE'].upper() == 'IDSO'):
            # Idso (1981)
            emissivity_clear = 0.7 + 5.95e-5 * vp * np.exp(1500 / air_temp)
        elif (self.options['LW_TYPE'].upper() == 'PRATA'):
            # Prata (1996)
            x = 46.5 * vp / air_temp
            emissivity_clear = 1 - (1 + x) * \
                np.exp(-1 * np.power((1.2 + 3 * x), 0.5))
        else:
            raise ValueError('Unknown LW_TYPE {0}'.format(
                self.options['LW_TYPE']))

        tskc = self.data['s_tskc']

        if (self.options['LW_CLOUD'].upper() == 'LW_CLOUD_DEARDORFF'):
            # Assume emissivity of clouds is 1.0, and that total emissivity is
            # weighted average of cloud emission plus clear-sky emission,
            # weighted by fraction of sky occupied by each
            # (method of Deardorff, 1978)
            emissivity = tskc * 1.0 + (1 - tskc) * emissivity_clear
        else:
            # see Bras 2.43
            emissivity = (1.0 + (0.17 * tskc * tskc)) * emissivity_clear

        self.data['s_lwrad'] = (emissivity * constants['STEFAN_B'] *
                                np.power(air_temp, 4))

    def _calc_trans(self):
        # STEP (1) calculate pressure ratio (site/reference) = f(elevation)
        pratio = np.power((1.0 - (constants['LR_STD'] *
                                  self.params['site_elev']) /
                           constants['T_STD']),
                          (constants['G_STD'] / (constants['LR_STD'] *
                           (constants['R'] / constants['MA']))))

        # STEP (2) correct initial transmittance for elevation
        return np.power(self.params['TBASE'], pratio)
