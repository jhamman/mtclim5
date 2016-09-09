'''
MTCLIM disaggregation routines.
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


# def disagg_precip():
#     # initialize */
#     for (rec = 0; rec < global_param.nrecs; rec++):
#         for (j = 0; j < NF; j++):
#             atmos[rec].prec[j] = 0
#                 for (idx = 0; idx < Ndays; idx++):
#     if (local_forcing_data[PREC][idx] > 0):
#       # Note: prcp_peaktime and prcp_duration must be in units of hours for this to work
#       h_offset = global_param.starthour - hour_offset_int
#             rec = (int)((idx*24-h_offset))/global_param.dt
#             tpeak = soil_con->prcp_peaktime[dmy[rec].month-1]
#             if (tpeak != NODATA_PD):
#         duration = soil_con->prcp_duration[dmy[rec].month-1]
#                 if (duration < 1):
#           duration = 1
#                   }
#         Pmax = 2*local_forcing_data[PREC][idx]/duration
#                 tstart = tpeak-0.5*duration
#                 tend = tpeak+0.5*duration
#                 hstart = (int)(tstart/options.SNOW_STEP)
#                 hpeak = (int)(tpeak/options.SNOW_STEP)
#                 hend = (int)(tend/options.SNOW_STEP)
#                 rec_start = (int)((idx*24+hstart-h_offset)/global_param.dt)
#                 j_start = (int)((idx*24+hstart-h_offset-rec_start*global_param.dt)/options.SNOW_STEP)
#                 hdiff = idx*24+hstart - (rec_start*global_param.dt+j_start*options.SNOW_STEP+h_offset)
#                 rec = rec_start
#                 j = j_start
#                 for (h=hstart; h<=hend; h+=options.SNOW_STEP, j+=options.SNOW_STEP):
#           hour = h-hdiff
#                     if (j==NF):
#             j=0
#                         rec++
#                       }
#           if (rec >= global_param.nrecs): break; }
#           if (hour+options.SNOW_STEP <= hpeak):
#             if (hour == hstart):
#               atmos[rec].prec[j] += Pmax/duration*(hour+options.SNOW_STEP-tstart)*(hour+options.SNOW_STEP-tstart)
#                           }
#             else:
#               atmos[rec].prec[j] += Pmax/duration*((hour+options.SNOW_STEP-tstart)*(hour+options.SNOW_STEP-tstart)-(hour-tstart)*(hour-tstart))
#                           }
#           }
#           else if (hour > hpeak):
#             if (hour == hend):
#               atmos[rec].prec[j] += Pmax/duration*( (duration*(tend-tpeak)-(tend-tpeak)*(tend-tpeak)) - (duration*(hour-tpeak)-(hour-tpeak)*(hour-tpeak)) )
#                           }
#             else:
#               atmos[rec].prec[j] += Pmax/duration*( (duration*(hour+options.SNOW_STEP-tpeak)-(hour+options.SNOW_STEP-tpeak)*(hour+options.SNOW_STEP-tpeak)) - (duration*(hour-tpeak)-(hour-tpeak)*(hour-tpeak)) )
#                           }
#           }
#           else if (hour == hpeak):
#             if (hour == hstart):
#               atmos[rec].prec[j] += 0.25*Pmax*duration
#                           }
#             else:
#               atmos[rec].prec[j] += Pmax/duration*((tpeak-tstart)*(tpeak-tstart)-(hour-tstart)*(hour-tstart))
#                           }
#             if (hour == hend):
#               atmos[rec].prec[j] += 0.25*Pmax*duration
#                           }
#             else:
#               atmos[rec].prec[j] += Pmax/duration*( duration*(hour+options.SNOW_STEP-tpeak)-(hour+options.SNOW_STEP-tpeak)*(hour+options.SNOW_STEP-tpeak) )

def set_max_min_hour(hourlyrad, ndays, tmaxhour, tminhour):
    '''
    This function estimates the times of minimum and maximum temperature for
    each day of the simulation, based on the hourly cycle of incoming solar
    radiation.'''

    for i in range(ndays):
        risehour = -999
        sethour = -999
        for hour in range(12):
            if (hourlyrad[i * 24 + hour] > 0 and
                    (i * 24 + hour == 0 or hourlyrad[i * 24 + hour - 1] <= 0)):
                risehour = hour
        for hour in range(12, 24):
            if (hourlyrad[i * 24 + hour] <= 0 and
                    hourlyrad[i * 24 + hour - 1] > 0):
                sethour = hour
        if i == ndays - 1 and sethour == -999:
            sethour = 23
        if risehour >= 0 and sethour >= 0:
            tmaxhour[i] = 0.67 * (sethour - risehour) + risehour
            tminhour[i] = risehour - 1
        else:
            # arbitrarily set the min and max times to 2am and 2pm
            tminhour[i] = 2
            tmaxhour[i] = 14


def disagg_tair(tair, method='hermite', **kwargs):

    return disagg_tair


def disagg_thermal():
    pass


def disagg_wind():
    pass
