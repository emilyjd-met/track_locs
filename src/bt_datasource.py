'''
Dictionaries of input information. These are set up as
value 0: Start date
value 1: End date
value 2: String representing subdirectories
value 4: Directory path
value 5: Pattern of filename
value 6: Optional value for number of offset days for icedrift

The dictionaries should be ordered so that the "top priority" files are
listed first, with lowest number. These will be preferentially selected.
'''

import os
from datetime import date, datetime, timedelta
import collections
from netCDF4 import Dataset

valid_vers = {'conc': ['v3'],
              'drift': ['cont'],
              'type': ['v3']}

today = date.today()
dm3 = date.today() - timedelta(days=3)
dm14 = date.today() - timedelta(days=14)
dm16 = date.today() - timedelta(days=16)
dm20 = date.today() - timedelta(days=20)

# ===================== ICE CONCENTRATION =====================
ic = {}

# V3 (SSMI) from thredds
ic['v3'] = {1: (date(1978, 10, 25), date(2020, 12, 31), 'ym',
                'https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/conc_450a_files/',
                'ice_conc_{a:}_ease2-250_cdr-v3p0_{d:%Y%m%d}1200.nc'
            ),
            2: (date(2021, 1, 1), dm3, 'ym',
                'https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/conc_cra_files/',
                'ice_conc_{a:}_ease2-250_icdr-v3p0_{d:%Y%m%d}1200.nc'
            )
}
ic['v3']['ref'] = "OSI SAF Global sea ice concentration climate data record 1978-2020 (v3.0, 2022), OSI-450-a, doi:10.15770/EUM_SAF_OSI_0013. EUMETSAT Ocean and Sea Ice Satellite Application Facility. OSI SAF Global sea ice concentration interim climate data record (v3.0, 2022), OSI-430-a, doi:10.15770/EUM_SAF_OSI_0014. EUMETSAT Ocean and Sea Ice Satellite Application Facility. All data extracted from OSI SAF thredds server."

# ===================== ICE DRIFT =====================

idr = {}

idr['cont'] = {1: (date(1991, 1, 1), date(2020, 12, 31), 'ym',
                  'https://thredds.met.no/thredds/dodsC/osisaf/met.no/reprocessed/ice/drift_455m_files/merged/',
                  'ice_drift_{a:}_ease2-750_cdr-v1p0_24h-{d:%Y%m%d}1200.nc',
                 1
                  ),
#               2: (date(2013, 3, 4), dm3, 'ym',
# Avoid summer season data gaps
               2: (date(2017, 5, 27), dm3, 'ym',
                  'https://thredds.met.no/thredds/dodsC/osisaf/met.no/ice/drift_lr/merged/',
                  'ice_drift_{a:}_polstere-625_multi-oi_{m:%Y%m%d}1200-{d:%Y%m%d}1200.nc',
                  2
               )
           }
idr['cont']['ref'] = "OSI SAF Global Low Resolution Sea Ice Drift data record 1991-2020 (v1, 2022), OSI-455, doi: 10.15770/EUM_SAF_OSI_0012. EUMETSAT Ocean and Sea Ice Satellite Application Facility. OSI SAF Global Low Resolution Sea Ice Drift, OSI-405-c, doi: 10.15770/EUM_SAF_OSI_NRT_2007. EUMETSAT Ocean and Sea Ice Satellite Application Facility. All data extracted from OSI SAF thredds server {}."

# ===================== ICE TYPE=====================

it = {}

it['v3'] = {1: (date(1978, 10, 25), date(2020, 12, 31), 'ym',
                    'https://thredds.met.no/thredds/dodsC/c3s/cdr_ice_type_v3p0_files/',
                    'ice_type_{a:}_ease2-250_cdr-v3p0_{d:%Y%m%d}1200.nc'
                    ),
                2: (date(2021, 1, 1), dm16, 'ym',
                    'https://thredds.met.no/thredds/dodsC/c3s/icdr_ice_type_v3p0_files/',
                    'ice_type_{a:}_ease2-250_icdr-v3p0_{d:%Y%m%d}1200.nc'
                   )
               }


# ==========================================================

def check_ice_file(fname):
    '''Check the validity of a filename for a file or a THREDDS server file'''

    goodfile = False
    # Testing THREDDS server files
    if "thredds" in fname:
        try:
            # Testing...
            with Dataset(fname) as dataset:
                la = dataset['lat']
            goodfile = True
        except:
            pass
    # Brief check for other files
    else:
        if os.path.exists(fname):
            goodfile = True

    return goodfile


def find_ice_file(pdate, hemi='nh', mode='conc', version='v3'):
    '''Find an ice file for a particular date and mode'''

    # Checking the input is not in datetime format, and convert if it is
    if isinstance(pdate, datetime):
        pdate = pdate.date()

    # Checking inputs
    valid_modes = ['conc', 'drift', 'type']
    valid_hemis = ['nh', 'sh']
    if hemi not in valid_hemis:
        raise ValueError("Variable hemi must be in {}: {}".format(valid_hemis, hemi))
    if mode not in valid_modes:
        raise ValueError("Variable mode must be in {}: {}".format(valid_modes, mode))

    if version not in valid_vers[mode]:
        raise ValueError("Variable version must be in {}: {}".format(valid_vers[mode], version))

    # Fetching the dictionary of input information, excluding the references
    exclude = ['ref']
    if mode == 'conc':
        inp_dict = {k: ic[version][k] for k in set(list(ic[version].keys())) - set(exclude)}
    elif mode == 'drift':
        inp_dict = {k: idr[version][k] for k in set(list(idr[version].keys())) - set(exclude)}
    elif mode == 'type':
        inp_dict = {k: it[version][k] for k in set(list(it[version].keys())) - set(exclude)}

    # Locating a file by date
    dfile = None
    for _, val in collections.OrderedDict(sorted(inp_dict.items())).items():
        # For each entry in the dictionary, check if the required date
        # falls in the date range specified (if there are more than one
        # file specified by entries in the dictionary, the first entry
        # from the ordered dictionary will be returned)
        sdate = val[0]
        edate = val[1]

        # Set a minus date according to the dictionary entry if there is one
        timestep = 1
        if len(val) >= 6 and val[5] is not None:
            timestep = val[5]
        mdate = pdate - timedelta(days=timestep)

        if pdate >= sdate and pdate <= edate:
            # Find the full filepath for the given date
            if val[2] == 'ym':
                subdir = '{:%Y/%m}'.format(pdate)
            elif val[2] == 'ymm':
                subdir = '{:%Y/%m}'.format(mdate)
            elif val[2] == 'ymd':
                subdir = '{:%Y/%m/%d}'.format(pdate)
            elif val[2] is None:
                subdir = ''
            filename = val[4].format(a=hemi, d=pdate, m=mdate)
            filepath = os.path.join(val[3], subdir, filename)

            # Check if this file exists (including if it is on THREDDS)
            if check_ice_file(filepath):
                dfile = filepath
                break

    # The file is returned if present, else None is returned
    if mode == 'conc':
        return dfile
    elif mode == 'drift':
        return dfile, timestep
    elif mode == 'type':
        return dfile
