import datetime
import pandas as pd
import numpy as np

examples = {'Seattle': './data/daymet_seattle.csv',
            'Telluride': './data/daymet_telluride.csv'}


def parse_daymet_dates(year, yday):
    return np.array([datetime.datetime(int(y), 1, 1) +
                     datetime.timedelta(days=int(d) - 1)
                     for y, d in zip(year, yday)])


def load_daymet_example(fname=examples['Seattle']):
    renames = {}
    units = {}
    df = pd.read_csv(fname, skiprows=6,
                     index_col='datetime',
                     parse_dates={'datetime': ['year', 'yday']},
                     date_parser=parse_daymet_dates)
    for c in df.columns:
        n, u = c.split(' ', 1)
        renames[c] = n
        units[n] = u

    df = df.rename(columns=renames)

    return df
