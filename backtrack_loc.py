'''
Calculate back trajectories for ice stations (mobile with the ice) or
moorings (fixed with respect to the ice) and write to an .nc file.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Created by Emily Down at Met Norway based on code by Thomas Lavergne,
funded by the SUDARCO project.

Please credit as "SUDARCO project/Met Norway"
'''


#TODO:
#Thomas suggests: I think this constant change between lat/lon and dx/dy is most useful when we change grid/projection (e.g. changing from CDR to NRT). But if we stay for several days in one projection, we could possibly stay with the dx/dy and not go through lat/lon. (we can still compute lat/lon, but do not need to go back).

import sys
import os
import re
import argparse
from datetime import datetime, date, timedelta
import warnings
from dateutil.relativedelta import relativedelta
import numpy as np
from netCDF4 import Dataset, date2num
import pyresample as pr
try:
    from pyresample.area_config import check_and_wrap
except ImportError as e:
    from pyresample.utils import check_and_wrap
import cartopy.crs as ccrs
from pyresample import geometry, _spatial_mp
import json

# We use some tools and routines developed for the OSI-450/SICCI SIC CDR
import bt_utils as utils
from bt_datasource import valid_vers, check_ice_file, find_ice_file, ic, idr

warnings.filterwarnings("ignore", message="Possible more than 10 neighbours within 150000 m for some data points")
warnings.filterwarnings("ignore", message="You will likely lose important projection information when converting to a PROJ string from another format. See: https://proj.org/faq.html#what-is-the-best-format-for-describing-coordinate-reference-systems")

# Dummy class to overwrite type of proj4_string (from unicode to str),
# because of bug in (older) pyproj
#from pyresample import geometry, _spatial_mp
#class myAreaDefinition(geometry.AreaDefinition):
#    def __init__(self, area_def):
#        for a in vars(area_def):
#            setattr(self, a, getattr(area_def, a))
#    @property
#    def proj4_string(self):
#        '''Returns projection definition as Proj.4 string'''
#        return str(super(myAreaDefinition, self).proj4_string)
#
#    def Proj(self):
#        return _spatial_mp.Proj(self.proj4_string)


# ################################################################
# INPUT ARGUMENTS
# ################################################################

def parse_args():

    p = argparse.ArgumentParser(description='')
    p.add_argument('-i', '--bid', required=False, default=None,
                   help="Name of the ice-station or mooring. Can be "
                   "a comma-separated list")
    p.add_argument('-e', '--enddate', required=False, default=None,
                   help="Date as YYYYmmddHHMM. This is the fixed date at "
                   "which the lat/lon is given and from which the back- "
                   "or forward-tracking should be made. For a mooring, "
                   "this is the end date of the timeseries required. "
                   "Can be a comma-separated list corresponding to each "
                   "ID, matching the length of that list")
    p.add_argument('-fe', '--firstenddate', required=False, default=None,
                   help="This should not be given for an ice station. "
                   "For a mooring, this gives the start time from which "
                   "trajectories are required. In YYYYmmddHHMM format. "
                   "Can be a comma-separated list corresponding to each "
                   "ID, matching the length of that list")
    p.add_argument('-la', '--elat', required=False, default=None,
                   help="Latitude at the enddate. This can be a comma-"
                   "separated list, then must match the length of the "
                   "ID list. Should be limited to positions in one "
                   "hemisphere only")
    p.add_argument('-lo', '--elon', required=False, default=None,
                   help="Longitude at the enddate. This can be a comma-"
                   "separated list, then must match the length of the "
                   "ID list")
    p.add_argument('-jf', '--json_file', required=False, default=None,
                   help="Path to a json file containing an input dataset, "
                   "with lons, lats and enddates")
    p.add_argument('-jd', '--json_data', required=False, default=None,
                   help="Name of the dataset within the json file. Required "
                   "if json file is used")
    p.add_argument('-o', '--outname', required=False, default='.',
                   help="Full output filepath, or path to directory in which "
                   "a default filename can be created. Default is '.'")
    p.add_argument('-p', '--period', required=False, default='-4y',
                   help="Maximum period over which to track the "
                   "trajectory. If with a d at the end, this is days. If "
                   "with a y, this is years")
    p.add_argument('-a', '--forwardtrack', default=False,
                   action='store_true',
                   help="Forward track instead of backtrack")
    p.add_argument('-n', '--projname', required=False, default=None,
                   help="Optional project name, otherwise the first "
                   "ID  or the json_data name will be used")
    p.add_argument('-ot', '--ow_threshold', required=False,
                   default='15',
                   help="Open water threshold to look for ow_repeat "
                   "days of lower concentration than this for closing "
                   "off a trajectory")
    p.add_argument('-otm', '--ow_threshold_mid', required=False,
                   default=None,
                   help="Intermediate open water threshold to flag. This "
                   "should be the higher concentration threshold.")
    p.add_argument('-or', '--ow_repeat', required=False, default='3',
                   help="Days of open water classified concentration "
                   "before a trajectory is closed off")
    p.add_argument('-v', '--iceconc_v', required=False, default='v3',
                   choices=valid_vers['conc'],
                   help="Version of the ice concentration file to use. "
                   "Options are {}".format(valid_vers['conc']))
    p.add_argument('-w', '--icedrift_v', required=False, default='cont',
                   choices=valid_vers['drift'],
                   help="Version of the ice drift file to use. "
                   "Options are {}".format(valid_vers['drift']))
    p.add_argument('-q', '--reproc', action='store_true', default=False,
                   help="Add YYYY/mm subdirs to the file creation")
    p.add_argument('-k', '--maxskip', required=False, default=4,
                   help="Maximum number of days to skip back if ice "
                   "concentration, drift or type file is missing")
    p.add_argument('-lmf', '--landmask_fp', required=False, default=None,
                   help='Filepath to extra landmask to be used for '
                   'checking where the trajectory ends')
    p.add_argument('-lmv', '--landmask_var', required=False,
                   default='coastmask_250',
                   help='Name of the landmask variable to be used for '
                   'checking where the trajectory ends')
    p.add_argument('-lmx', '--landmask_val', required=False, default='2',
                   help='Name of the landmask value to be used for '
                   'checking where the trajectory ends')
    p.add_argument('-z', '--verbose', action='store_true', default=False,
                   help="Print extra information while running the code")
    p.add_argument('-f', '--force', action='store_true', default=False,
                   help="Force overwrite of backtracking file if it "
                   "already exists")

    args = p.parse_args()

    # Checking a correct set of inputs is defined
    if args.bid or args.enddate or args.elat or args.elon:
        if not (args.bid and args.enddate and args.elat and args.elon):
            raise ValueError("Arguments -i/--bid, -e/--enddate, -la/--elat "
                             "and -lo/--elon must be used together. If one is "
                             "defined, they all should be")
    if args.json_file or args.json_data:
        if not (args.json_file and args.json_data):
            raise ValueError("Arguments "
                             "must be used together. json_data is the name "
                             "of the dataset to be used within the json file")
    if args.bid and args.json_file:
        raise ValueError("Either the dataset should be defined directly "
                         "using -i/--bid, -e/--enddate, -la/--elat and "
                         "-lo/--elon, or json inputs should be used with "
                         "-jf/--json_file and -jd/--json_data. Do not use "
                         "both")
    if not args.bid and not args.json_file:
        raise ValueError("The dataset should either be defined directly "
                         "using -i/--bid, -e/--enddate, -la/--elat and "
                         "-lo/--elon, or json inputs should be used with "
                         "-jf/--json_file and -jd/--json_data")


    return args


def print_percentage(step, total_steps):

    percentage = int((step / total_steps) * 100)
    sys.stdout.write('\r')
    #sys.stdout.write(f'{percentage:.2f}%')
    sys.stdout.write(f'Tracking... {percentage}%')
    sys.stdout.flush()


def create_bt_name(outname, firstenddate, enddate, bid, json_data, projname,
                   forwardtrack=False, reproc=True):
    '''Find the name of the backtracking file'''

    if firstenddate is not None:
        datestr = '2d_{:%Y%m%d}-{:%Y%m%d}'.format(firstenddate[0], enddate[0])
    else:
        datestr = '1d_{:%Y%m%d}'.format(enddate[0])

    # If no project name is specified, use the first ID
    if projname is None:
        if json_data is not None:
            projname = json_data
        else:
            projname = bid[0].replace(" ", "")

    if forwardtrack:
        outf = 'forwardtrack_loc_{}_{}.nc'.format(projname, datestr)
    else:
        outf = 'backtrack_loc_{}_{}.nc'.format(projname, datestr)

    if reproc:
        outname = os.path.join(outname, '{:%Y}'.format(enddate),
                               '{:%m}'.format(enddate))
        os.makedirs(outname, exist_ok=True)

    outname = os.path.join(outname, outf)

    return outname, projname


def stripdate(listdirs):
    '''Strip month and year subdirectories from a list of directory names'''

    pattm = '^[0-9]{2}$'
    patty = '^[0-9]{4}$'
    cleandirs = []

    for dr in listdirs:
        cleaning = True
        while cleaning:
            if (re.match(pattm, os.path.basename(dr)) or
                re.match(patty, os.path.basename(dr))):
                dr = os.path.dirname(dr)
            else:
                cleandirs.append(dr)
                cleaning = False
                break

    return cleandirs


def begindateandper(enddate, period, forwardtrack=False):
    '''Given the enddate and period, find the beginning date'''

    if period.endswith("d") or period.endswith("w"):
        per = int(period[:-1])
        if period.endswith("w"):
            per = int(per) * 7
        if forwardtrack:
            begindate = enddate + timedelta(days=int(per))
        else:
            begindate = enddate - timedelta(days=int(per))
    elif period.endswith("y"):
        if forwardtrack:
            begindate = enddate + relativedelta(years=int(period[:-1]))
            per = (begindate - enddate).days
        else:
            begindate = enddate - relativedelta(years=int(period[:-1]))
            per = (enddate - begindate).days
    else:
        raise ValueError("Period should be specified in terms of days, weeks or years")

    # Checking for forward tracking that the "begindate" is not after
    # today - 16 days (OSI SAF's NRT data delay.)
    shiftdays = 16
    shiftper = timedelta(days=shiftdays)
    if forwardtrack:
        if begindate > datetime.today() - shiftper:
            if enddate < datetime.today() - shiftper:
                begindate = datetime.today() - shiftper
                per = per - (begindate - datetime.today()).days - shiftdays
                print("WARNING: forward tracking is selected and the period "
                      "gives a date after today, setting the tracking period "
                      "up to today.")
            else:
                print("The date to forward track from is too late in time. "
                      "This cannot be later than today - 16 days, as current "
                      "ice concentration data is needed.")
    return begindate, per


def convtodt(dt):
    '''Check a time is in python datetime format and convert if not. If Nan,
    just return NaN'''

    # If this is not a number, return immediately
    if np.isnan(dt):
        return dt

    # If this is already python datetime, return immediately
    if isinstance(dt, datetime):
        return dt

    # Convert datetime64
    if isinstance(dt, np.datetime64):
        unix_epoch = np.datetime64(0, 's')
        one_second = np.timedelta64(1, 's')
        seconds_since_epoch = (dt - unix_epoch) / one_second
        return datetime.utcfromtimestamp(seconds_since_epoch)

    raise ValueError("Unknown type of time {}".format(dt))


# ################################################################
# LANDMASK FOR TRIMMING TRAJECTORIES
# ################################################################

def read_extra_landmask(landmask_fp, landmask_var, landmask_val):

    with Dataset(landmask_fp, 'r') as dataset:

        lats = dataset['lat'][:]
        lons = dataset['lon'][:]
        ex_lmask = np.zeros_like(lats, dtype=int)
        lmdata = dataset[landmask_var][:]
        ex_lmask[lmdata >= landmask_val] = 1

        ex_lmask_swath = pr.geometry.SwathDefinition(lons=lons, lats=lats)

    return ex_lmask, ex_lmask_swath


# ################################################################
# LOCATE AND LOAD ICE DRIFT PRODUCT
# ################################################################

def retrieve_icedrift_file(dt, hemi='nh', icedrift_v='cdr-v1', maxskip=4,
                           timestep=1, forwardtrack=False, verbose=False):
    '''Locate the ice drift file. Note that dt is the 'end-date' '''

    # TODO - with forward tracking, should this look forward in time? Or
    # can it still look back?
    drift_skip = 0
    fdt = dt
    fname = None
    drift_found = False
    # Look for the file up to maxskip days back in time (backtracking) or
    # maxskip days forward in time (forward tracking)
    for try_dt in range(0, maxskip):
        if try_dt > 0:
            if forwardtrack:
                fdt = fdt + timedelta(days=1)
            else:
                fdt = fdt - timedelta(days=1)
            drift_skip += 1
            if verbose:
                print("\nWARNING: for ice drift try again with {}".format(fdt))
        fname, dtimestep = find_ice_file(fdt, hemi=hemi, mode='drift',
                                         version=icedrift_v)

        if fname is not None:
            if check_ice_file(fname):
            #if os.path.exists(fname):
                # Read the icedrift file here
                try:
                    lats1, lons1, bdXs, bdYs = read_icedrift_file(fname, fdt,
                                                    dtimestep=dtimestep,
                                                    timestep=timestep)
                    drift_found = True
                    break
                except:
                    continue

    if not drift_found:
        raise ValueError("Missing sea-ice drift file for version {} "
                         "date {} hemisphere {}".format(icedrift_v, dt, hemi))
    elif drift_skip > 0:
        print("\tGAP WARNING: No ice drift data available for {}, using data "
              "for {} instead (skip {} days)".format(
                  datetime.strftime(dt, '%Y%m%d'),
                  datetime.strftime(fdt, '%Y%m%d'), drift_skip))

    driftdata = {}
    driftdata['icedrift_file'] = fname
    driftdata['dirname'] = os.path.dirname(fname)
    driftdata['lats1'] = lats1
    driftdata['lons1'] = lons1
    driftdata['bdXs'] = bdXs
    driftdata['bdYs'] = bdYs
    driftdata['dtimestep'] = dtimestep
    driftdata['drift_skip'] = drift_skip

    return driftdata


def check_vars_in_dataset(varlist, dataset):
    '''Snippet to check if all required variables are in a NetCDF dataset'''
    allvars = True
    for var in varlist:
        if not var in dataset.variables:
            allvars = False

    return allvars


def read_icedrift_file(driftfile, pdate, dtimestep=1, timestep=1):
    '''Read lats, lons, and drifts from the icedrift file'''

    with Dataset(driftfile, 'r') as dataset:

        # Checking which file type
        reqvarmet = ['lat1', 'lon1', 'dX', 'dY']
        reqvarnsidc = ['latitude', 'longitude', 'u', 'v']

        # Read in MET final format file
        if check_vars_in_dataset(reqvarmet, dataset):
            lats1 = dataset.variables['lat1'][0, :].compressed()
            lons1 = dataset.variables['lon1'][0, :].compressed()
            dXs = dataset.variables['dX'][0, :].compressed()
            dYs = dataset.variables['dY'][0, :].compressed()
            units = dataset.variables['dX'].units
            if (lons1.size != lats1.size or dXs.size != lats1.size
                or dYs.size != lats1.size):
                raise ValueError("Fields lat1, lon1, dX, and dY do not have "
                                 "the same mask!")

            # Note that the lons1 and lats1 we read are the
            # 'end-positions' of the drift vectors, and are not on a
            # regular grid. This is because we do backtracking.

            # Since we want to go back in time, return the opposite of the
            # icedrift vector.
            if units in ['km', 'kilometer', 'kilometers']:
                scale = -1. * 1000.
            elif units in ['m', 'meter', 'meters']:
                scale = -1.
            dXs = scale * dXs
            dYs = scale * dYs

            # The y-axis is flipped in OSI SAF data
            # No, not anymore!!
            #dYs = -1. * dYs

        # Read in NSDIC ice drift file
        elif check_vars_in_dataset(reqvarnsidc, dataset):

            # The lats and lons given for NSIDC data are at the start
            # of the drift vector
            lat = dataset.variables['latitude'][:]
            lon = dataset.variables['longitude'][:]
            u = dataset.variables['u'][:]
            v = dataset.variables['v'][:]
            # The drift u and v are indexed by (time, y, x), so first
            # determine and extract the correct timestep of the drift
            if isinstance(pdate, datetime):
                dayindex = (pdate - datetime(1970, 1, 1)).days
            else:
                dayindex = (pdate - date(1970, 1, 1)).days
            vartime = dataset.variables['time'][:]
            index = np.where(vartime == dayindex)[0][0]
            dx = u[index, :, :]
            dy = v[index, :, :]
            # Exit gracefully if all data is masked otherwise it throws an
            # error
            if dx.mask.all():
                return None

            # Lats and lons are for the whole field, not just the drift
            # positions, so extract the relevant lat/lons
            lat.mask = dx.mask
            lon.mask = dx.mask
            # Compress the arrays
            lats0 = lat.compressed()
            lons0 = lon.compressed()
            dXs = dx.compressed()
            dYs = dy.compressed()

            # Convert the drift units from cm/s to m/day
            dXs = dXs * (24 * 3600.) / (100.)
            dYs = dYs * (24 * 3600.) / (100.)
            # Find the data CRS of the dataset (hardcoded as EASE for
            # NSIDC
            if lats0[0] > 0:
                data_ccrs = ccrs.LambertAzimuthalEqualArea(
                    central_longitude=0, central_latitude=90,
                    false_easting=0, false_northing=0)
            else:
                data_ccrs = ccrs.LambertAzimuthalEqualArea(
                    central_longitude=0, central_latitude=-90,
                    false_easting=0, false_northing=0)
            # Need to modify the lats and lons since we want the end
            # position for these rather than the start position, since
            # we are doing backtracking
            lats1 = np.empty_like(lats0)
            lons1 = np.empty_like(lons0)
            pc = ccrs.PlateCarree()
            for i in range(len(lats0)):
                x0, y0 = data_ccrs.transform_point(lons0[i], lats0[i], pc)
                x1 = x0 + dXs[i]
                y1 = y0 + dYs[i]
                lons1[i], lats1[i] = pc.transform_point(x1, y1, data_ccrs)
            # The negative must be returned since the motion is tracked
            # backwards
            scale = -1.
            dXs = scale * dXs
            dYs = scale * dYs
        # No further ice drift file types implemented yet
        else:
            raise ValueError("File {} appears to be neither a MET final "
                             "format ice drift file nor an NSDIC file"
                             "".format(driftfile))

    # If the ice drift needs to be scaled to the correct timestep, do
    # this here
    if dtimestep != timestep:
        dXs = dXs * timestep / dtimestep
        dYs = dYs * timestep / dtimestep

    return lats1, lons1, dXs, dYs


# ################################################################
# LOCATE AND LOAD ICE CONCENTRATION PRODUCT
# ################################################################

def find_iceconc_file(dt, hemi='nh', iceconc_v='cdr-v2', maxskip=4,
                      forwardtrack=False, verbose=False):
    '''Locate an ice concentration file for the correct version and date'''

    # Look for the file up to maxskip days back in time
    # TODO - with forward tracking it is still OK for it to look backward
    # in time if it doesn't find the file it wants?
    conc_skip = 0
    fdt = dt
    fname = None
    for try_dt in range(0, maxskip):
        if try_dt > 0:
            if forwardtrack:
                fdt = fdt - timedelta(days=1)
            else:
                fdt = fdt - timedelta(days=1)
            conc_skip += 1
            if verbose:
                print("\nWARNING: for iceconc try again with {}".format(fdt))
        fname = find_ice_file(fdt.date(), hemi=hemi, version=iceconc_v)
        if fname is not None:
            if check_ice_file(fname):
#                if os.path.exists(fname):
                break

    if fname is None:
        raise ValueError("Missing sea-ice concentration file for version {} "
                         "date {} hemisphere {}".format(iceconc_v, dt, hemi))
    elif conc_skip > 0:
        print("GAP WARNING: No ice conc data available for {}, using data "
              "fo {} instead (skip {} days)".format(
                  datetime.strftime(dt, '%Y%m%d'),
                  datetime.strftime(fdt, '%Y%m%d'), conc_skip))

    return fname, conc_skip


def iceconc_areadef_lmask(dt, hemi='nh', iceconc_v='cdr-v2', maxskip=4,
                          forwardtrack=False, verbose=False):
    '''Find the area definition and the landmask from the ice concentration'''

    concf, _ = find_iceconc_file(dt, hemi=hemi, iceconc_v=iceconc_v,
                                 maxskip=maxskip, forwardtrack=forwardtrack,
                                 verbose=verbose)
    iceconc_area_defs = utils.load_grid_defs_from_OSISAF_ncCF_file(concf,
                            verbose=verbose)
    iceconc_area_defs.area_id = hemi

    with Dataset(concf) as _:
        sflags = _.variables['status_flag'][0, :]
        if iceconc_v == 'cdr-v1':
            iceconc_lmask = (sflags == 100)
        else:
            iceconc_lmask = (sflags & 1 == 1)

    return iceconc_area_defs, iceconc_lmask


def read_iceconc_file(f):
    '''Read ice concentration variable from the file'''

    with Dataset(f, 'r') as _:
        iceconc = _.variables['ice_conc'][0, :]

    return iceconc


def get_ice_conc(dt, area_def, lmask, iceconc_v='cdr-v2', maxskip=4,
                 forwardtrack=False, verbose=False):
    '''Retrieve the ice concentration, interpolated over land'''

    f, conc_skip = find_iceconc_file(dt, hemi=area_def.area_id,
                                     iceconc_v=iceconc_v, maxskip=maxskip,
                                     forwardtrack=forwardtrack, verbose=verbose)
    dat = read_iceconc_file(f)

    # Interpolate sea-ice concentration over land (in case a trajectory
    # drifts onto land)
    if lmask.shape != dat.shape or lmask.shape != area_def.shape:
        raise ValueError("Error with provided land mask to perform land "
                         "interpolation of ice conc")
    sic_lons, sic_lats = area_def.get_lonlats()
    src = pr.geometry.SwathDefinition(lons=sic_lons[~lmask],
                                      lats=sic_lats[~lmask])
    trg = pr.geometry.SwathDefinition(lons=sic_lons[lmask],
                                      lats=sic_lats[lmask])
    landsic = pr.kd_tree.resample_nearest(src, dat[~lmask], trg,
                                          radius_of_influence=150000,
                                          fill_value=None)
    dat[lmask] = landsic

    return dat, conc_skip, f


def extract_iceconc(dt, lons, lats, area_def, lmask, iceconc_v='cdr-v2',
                    maxskip=4, forwardtrack=False, verbose=False):
    '''Retrieve the ice concentration, resampled to the area definition
    provided'''

    trg = pr.geometry.SwathDefinition(lons=np.array(lons), lats=np.array(lats))
    dat, conc_skip, concfile = get_ice_conc(dt, area_def, lmask,
                                            iceconc_v=iceconc_v,
                                            maxskip=maxskip,
                                            forwardtrack=forwardtrack,
                                            verbose=verbose)
    resampled_iceconc = pr.kd_tree.resample_nearest(area_def, dat, trg,
                                                    radius_of_influence=20000,
                                                    fill_value=None)
    concdir = os.path.dirname(concfile)

    return resampled_iceconc, conc_skip, concdir


# ################################################################
# BACKTRACKING AND ICE TYPE CLASSIFICATION
# ################################################################

def backtrace_traj_step(dt, lons, lats, area_def, forwardtrack=False,
                        icedrift_v='cdr-v1', maxskip=4, timestep=1,
                        verbose=False):
    '''Backtrack a single step by adjusting the input lats/lons according
    to the ice drift'''

    # TODO - do I actually need to pass forwardtrack here?
    
    trg = pr.geometry.SwathDefinition(lons=lons, lats=lats)

    # Locate and read the ice drift file
    try:
        driftdata = retrieve_icedrift_file(dt, hemi=area_def.area_id,
                                           icedrift_v=icedrift_v,
                                           maxskip=maxskip,
                                           timestep=timestep,
                                           forwardtrack=forwardtrack,
                                           verbose=verbose)
        lats1 = driftdata['lats1']
        lons1 = driftdata['lons1']
        bdXs = driftdata['bdXs']
        bdYs = driftdata['bdYs']
        drift_skip = driftdata['drift_skip']

    except Exception as ex:
        raise ValueError("Drift file not locatable/readable for {:%Y-%m-%d}. "
                         "({})".format(dt, ex))
    # If lats1 is an empty array, this causes a crash
    if lats1.size == 0:
        raise ValueError("No latitudes read from file {} \n".format(f))

    # NOTE: Thomas suggests to not change back from dx,dy to lat,lon all
    # the time (see TODO note at top of code)
    
    # Extract the dXs and dYs to the locations of the input (lons, lats)
    # TODO - Is this close enough, since they are vectors? No rotation
    # accounted for.
    src = pr.geometry.SwathDefinition(lons=lons1, lats=lats1)
    resampled_bdrift = pr.kd_tree.resample_gauss(src, np.dstack((bdXs, bdYs)),
                                                 trg,
                                                 radius_of_influence=150000,
                                                 sigmas=[62500, 62500],
                                                 neighbours=10)

    # Move the input lons and lats by the (bdXs, bdYs) vectors
    # This requires knowledge of the pyproj object associated to the ice
    # drift file
    pobj = _spatial_mp.Proj_MP(area_def.proj_dict)

    # Get x/y position of input lons/lats:
    x0s, y0s = pobj(lons, lats)
    # Move the points by the resampled bdXs, bdYs
    if forwardtrack:
        x1s = x0s - resampled_bdrift[:, 0]
        y1s = y0s - resampled_bdrift[:, 1]
    else:
        x1s = x0s + resampled_bdrift[:, 0]
        y1s = y0s + resampled_bdrift[:, 1]
    # Transform back to lons/lats
    new_lons, new_lats = pobj(x1s, y1s, inverse=True)
    return new_lons, new_lats, drift_skip


def count_repeats_2d(a, repeats=3):
    '''Count the number of repeats of a condition along the 0-axis (time)
    of the 2d array.
    Adapted from the 1d example at http://stackoverflow.com/a/24343375 '''

    b = np.concatenate((a[None, 0, :], a[:-1, :] != a[1:, :],
                        np.ones((1, a.shape[1])).astype('bool')), axis=0)
    reps = np.zeros(a.shape[1]).astype('bool')
    for i in range(a.shape[1]):
        reps[i] = (np.diff(np.where(b[:, i])[0])[::2] >= repeats).any()

    return reps


def check_ow(conc, ow_threshold=15, ow_repeat=3):
    '''Check the conc is not less than the open water threshold for more
    than ow_repeat steps - if so, set these steps to NaN'''

    # The filled(150.) is because we do not want masked values to
    # interfere with the checking
    fconc = conc.filled(150)
    for traj in range(len(fconc[0, :])):
        if np.all(fconc[-ow_repeat:, traj] < ow_threshold):
            fconc[-ow_repeat:, traj] = np.nan

    return conc


# ################################################################
# MAIN PROCESSOR
# ################################################################

def backtrack_traj(enddates, endlons, endlats, begindates, periods,
                   forwardtrack=False,
                   hemi='nh', iceconc_v='cdr-v2', icedrift_v='cdr-v1',
                   maxskip=4, ow_threshold=15, ow_threshold_mid=None,
                   ow_repeat=3, ex_lmask=None, ex_lmask_swath=None,
                   verbose=False):
    '''Backtrack or forward track until the time limit or a gap in the ice'''

    concdirs = []
    driftdirs = []

    maxper = np.nanmax(np.abs(periods))
    
    # Setting up the area definitions
    iceconc_area_def, iceconc_lmask = iceconc_areadef_lmask(enddates[0],
                                                            hemi=hemi,
                                                            iceconc_v=iceconc_v,
                                                            maxskip=maxskip,
                                                    forwardtrack=forwardtrack,
                                                            verbose=verbose)

    # Finding the initial icedrift file (just use timestep of 1 - this gets
    # overwritten by thedrift timestep in a moment)
    driftdata = retrieve_icedrift_file(enddates[0].date(), hemi=hemi,
                                       icedrift_v=icedrift_v,
                                       maxskip=maxskip,
                                       timestep=1, verbose=verbose)
    if driftdata['dirname'] not in driftdirs:
        driftdirs.append(driftdata['dirname'])
    drift_file = driftdata['icedrift_file']
    max_drift_skip = driftdata['drift_skip']
    adef = utils.load_grid_defs_from_OSISAF_ncCF_file(drift_file,
                                                      verbose=verbose)
    area_def = adef
    area_def.area_id = hemi

    # Setting up the timestep required for the backtracking. At the moment
    # this is always 1 day
    timestep = 1
    # Initialize tracking
    arrlen = len(endlons)
    # The status array keeps track of trajectory status. Status 0 is that
    # this trajectory has not begun yet. Status 1 is a trajectory in progress.
    # Status 2 is a completed trajectory.
    status = np.full((arrlen), 0)
    count_ow = np.full((arrlen), 0)
    count_ow_mid = np.full((arrlen), 0)
    count_land = np.full((arrlen), 0)
    # The trajstep array keeps track of which timestep each trajectory is on
    # for filling in the array
    trajstep = np.full((arrlen), -1)
    lons = np.full((arrlen, maxper + 2), np.nan, dtype='float')
    lats = np.full((arrlen, maxper + 2), np.nan, dtype='float')
    conc = np.full((arrlen, maxper + 2), 150., dtype='float')
    flag = np.full((arrlen, maxper + 2), 0, dtype='int')
    lons[:, 0] = endlons
    lats[:, 0] = endlats

    if forwardtrack:
        tpoint = np.nanmin(enddates)
    else:
        tpoint = np.nanmax(enddates)
    valid_date = True

    loopcount = 0
    while valid_date and not np.all(status == 2):

        print_percentage(loopcount, maxper)
        
        # Check if any new trajectories are happening at this timepoint,
        # and if so the concentration needs to be initialised
        if np.any(enddates == tpoint):
            status[enddates == tpoint] = 1
            trajstep[enddates == tpoint] = 0
            # The lats/lons are already initialised
            concinit, max_conc_skip, concdir = extract_iceconc(tpoint,
                                                endlons[enddates == tpoint],
                                                endlats[enddates == tpoint],
                                                iceconc_area_def,
                                                iceconc_lmask,
                                                iceconc_v=iceconc_v,
                                                maxskip=maxskip,
                                                forwardtrack=forwardtrack,
                                                verbose=verbose)
            conc[enddates == tpoint, 0] = concinit
            if concdir not in concdirs:
                concdirs.append(concdir)

        # Do tracking at the points where the status is > 1
        trp = status >= 1
        lons[trp, trajstep[trp] + 1], lats[trp, trajstep[trp] + 1], \
            drift_skip = backtrace_traj_step(tpoint,
                                             lons[trp, trajstep[trp]],
                                             lats[trp, trajstep[trp]],
                                             area_def,
                                             forwardtrack=forwardtrack,
                                             icedrift_v=icedrift_v,
                                             maxskip=maxskip,
                                             timestep=timestep,
                                             verbose=verbose)

        if drift_skip > max_drift_skip:
            max_drift_skip = drift_skip

        # If an extra landmask is provided for trimming trajectories that
        # go over land, check this here and close any of these trajectories
        if ex_lmask is not None:

            swath_newll = pr.geometry.SwathDefinition(
                lons=lons[trp, trajstep[trp] + 1],
                lats=lats[trp, trajstep[trp] + 1])
            land_points = np.full((arrlen), 0)
            land_points[trp] = pr.kd_tree.resample_nearest(ex_lmask_swath,
                                                      ex_lmask, swath_newll,
                                                    radius_of_influence=150000,
                                                      fill_value=None)

            count_land[land_points == 1] += 1
            count_land[land_points == 0] = 0
            # TODO - For now use the same number of repeats for open water
            # as land
            close = count_land >= ow_repeat
            status[close] = 2            
#            status[land_points] = 2

        # Go one timestep in time (to the time of the new position)
        # before finding the ice concentration at this time
        if forwardtrack:
            tpoint = tpoint + timedelta(days=timestep)
        else:
            tpoint = tpoint - timedelta(days=timestep)
        new_conc, conc_skip, concdir = extract_iceconc(tpoint,
                                            lons[trp, trajstep[trp] + 1],
                                            lats[trp, trajstep[trp] + 1],
                                            iceconc_area_def,
                                            iceconc_lmask,
                                            iceconc_v=iceconc_v,
                                            maxskip=maxskip,
                                            forwardtrack=forwardtrack,
                                            verbose=verbose)
        tmp_conc = np.full_like(status, 150., dtype=float)
        tmp_conc[trp] = new_conc


        # If there is a mid open water threshold, check this and flag
        # as necessary
        # Flags are 0 = nominal, 1 = concentration below "mid" concentration
        # and -1 for a closed trajectory
        if ow_threshold_mid is not None:
            # Flag any arrays that are over concentration between the mid
            # and standard ice concentrations for more than ow_repeat
            # successive steps
            ow_mid = np.logical_and(
                np.logical_or(np.isnan(tmp_conc),
                              tmp_conc <= ow_threshold_mid),
#                np.logical_or(tmp_conc > ow_threshold,
#                              tmp_conc <= ow_threshold_mid)),
                status == 1)
            count_ow_mid[ow_mid] += 1
            count_ow_mid[~ow_mid] = 0
            mid = count_ow_mid >= ow_repeat
            # Flag all that row after this trajectory is closed
            for i in range(arrlen):
                if mid[i]:
                    flag[i, trajstep[i] - ow_repeat:] = 1

        # Close any arrays that are over open water according to the standard
        # concentration threshold for more than ow_repeat successive steps
        ow = np.logical_and(np.logical_or(np.isnan(tmp_conc),
                                          tmp_conc <= ow_threshold),
                            status == 1)
        count_ow[ow] += 1
        count_ow[~ow] = 0
        close = count_ow >= ow_repeat
        status[close] = 2
        # When the trajectory is closed, we don't want these last lons/lats
        for i in range(arrlen):
            if close[i]:
                lons[i, trajstep[i] - ow_repeat: trajstep[i]] = np.nan
                lats[i, trajstep[i] - ow_repeat: trajstep[i]] = np.nan
                flag[i, trajstep[i] - ow_repeat:] = -1


#        if sum(close) > 0:
#            for i in range(ow_repeat):
#                lons[close, trajstep[close] + 1 - i] = np.nan
#                lats[close, trajstep[close] + 1 - i] = np.nan
        #lons[close, trajstep[close] + 1 - ow_repeat: trajstep[close] + 1] = np.nan
        #lats[close, trajstep[close] + 1 - ow_repeat: trajstep[close] + 1] = np.nan
#            flag[close, trajstep[close] + 1 - ow_repeat: ] = -1

        # Increment the trajectory index
        trajstep[status == 1] += 1
        # Check if any trajectories should be closed because this exceeds
        # the maximum period
        status[trajstep >= periods] = 2

        if concdir not in concdirs:
            concdirs.append(concdir)
        if conc_skip > max_conc_skip:
            max_conc_skip = conc_skip
#        if verbose:
#            if forwardtrack:
#                print('Forward tracked to {}'.format(tpoint))
#            else:
#                print('Backtracked to {}'.format(tpoint))

        # Check if the latest date is still in range
        if forwardtrack:
            if tpoint >= max(begindates):
                valid_date = False
        else:
            if tpoint <= min(begindates):
                valid_date = False

        
        loopcount += 1
    # End of while loop

    # Trim the lon/lats to the actual timesteps needed for valid trajectories
    maxts = np.max(trajstep) + 1
    lons = lons[:, :maxts]
    lats = lats[:, :maxts]
    flag = flag[:, :maxts]

    # Find the beginning dates
    newbegindates = []
    for i, edatei in enumerate(enddates):
        try:
            dur = int(np.where(flag[i, :] == -1)[0][0])
        except:
            dur = len(flag[i, :])
        newbegindates.append(edatei - timedelta(days=dur))

    max_skips = {'drift': max_drift_skip,
                 'conc': max_conc_skip}
    bt_out = {}
    bt_out['traj_lons'] = lons
    bt_out['traj_lats'] = lats
    bt_out['flag'] = flag
    bt_out['max_timestep'] = maxts
    bt_out['max_skips'] = max_skips
    bt_out['concdir'] = concdirs
    bt_out['driftdir'] = driftdirs
    bt_out['enddates'] = enddates
    bt_out['begindates'] = newbegindates

    return bt_out


# ################################################################
# MAIN CODE
# ################################################################


def backtrack_loc(bid=None, enddate=None, firstenddate=None, elon=None,
                  elat=None, json_file=None, json_data=None, outname='.',
                  period='-4y', forwardtrack=False, projname=None,
                  ow_threshold=15., ow_threshold_mid=None, ow_repeat=3,
                  iceconc_v='v3', icedrift_v='cont', reproc=False, maxskip=4,
                  landmask_fp=None, landmask_var='coastmask_250',
                  landmask_val=2, verbose=False, force=False):

    # If there is an extra landmask, read this into a swath definition
    if landmask_fp is not None:
        ex_lmask, ex_lmask_swath = read_extra_landmask(landmask_fp,
                                                       landmask_var,
                                                       landmask_val)
    else:
        ex_lmask = None
        ex_lmask_swath = None

    # If the inputs are supplied by json file, read this into arrays
    if json_file:
        with open(json_file, 'r') as jfile:
            jdata = json.load(jfile)[json_data]
        bid = []
        enddate = []
        elon = []
        elat = []
        mode = '1d'
        if 'firstenddate' in jdata[list(jdata.keys())[0]]:
            firstenddate = []
            mode = '2d'
        for jkey in jdata:
            bid.append(jkey)
            enddate.append(jdata[jkey]['enddate'])
            elon.append(jdata[jkey]['lon'])
            elat.append(jdata[jkey]['lat'])
            if mode == '2d':
                firstenddate.append(jdata[jkey]['firstenddate'])
        # Dateify the dates
        enddate = [datetime.strptime(str(x), "%Y%m%d%H%M") for x in enddate]
        enddate = np.array(enddate)
        if mode == '2d':
            firstenddate = [datetime.strptime(str(x), "%Y%m%d%H%M")
                            for x in firstenddate]
            firstenddate = np.array(firstenddate)
        # Converting to numpy arrays
        bid = np.array(bid)
        elon = np.array(elon)
        elat = np.array(elat)

    # Otherwise convert the variables as needed
    else:

        # Converting the dates to datetime format. If firstenddate is specified,
        # then this is 2-d mode (fixed moorings) else this is 1-d mode (ice
        # stations)
        try:
            enddate = [datetime.strptime(x, "%Y%m%d%H%M")
                       for x in enddate.split(',')]
        except:
            raise ValueError("Please specify enddate as a comma-separated list "
                             "of values in format YYYYmmddHHMM: {}".format(
                                 enddate))
        if firstenddate is not None:
            try:
                firstenddate = [datetime.strptime(x, "%Y%m%d%H%M")
                                for x in firstenddate.split(',')]
            except:
                raise ValueError("Please specify firstenddate as a "
                                 "comma-separated list of values in format "
                                 "YYYYmmddHHMM: {}".format(firstenddate))
            mode = '2d'
        else:
            mode = '1d'

        # Converting the endpoint lat/lon and buoy IDs to a list
        elon = np.array([float(x) for x in elon.split(',')])
        elat = np.array([float(x) for x in elat.split(',')])
        bid = np.array([x for x in bid.strip().split(',')])

    # Checking if the lats/lons need wrapping
    elon, elat = check_and_wrap(elon, elat)

    # Checking which hemisphere we are working in
    if np.all(elat >= 0):
        hemi = 'nh'
    elif np.all(elat <= 0):
        hemi = 'sh'
    else:
        raise ValueError("Please do not specify locations in both hemispheres")

    # Make the enddates, endlats, endlons and begindates into 1-d arrays
    if mode == '1d':

        print("{}: Running backtracking in 1-D mode".format(bid[0]))

        elat_1d = np.array(elat)
        elon_1d = np.array(elon)
        enddate_1d = [e.replace(hour=12, minute=0) for e in enddate]
        begindate_1d = []
        per_1d = []
        for ed in enddate_1d:
            bd, p = begindateandper(ed, period, forwardtrack=forwardtrack)
            begindate_1d.append(bd)
            per_1d.append(p)
        enddate_1d = np.array(enddate_1d, dtype='datetime64[s]')
        begindate_1d = np.array(begindate_1d, dtype='datetime64[s]')
        per_1d = np.array(per_1d)
        # For the 1-d case, the timespan is 1 since there is only one
        # date in this dimension
        maxtimespan = 1

    elif mode == '2d':

        print("{}: Running backtracking in 2-D mode".format(bid[0]))

        # Find the timespans from the enddate and the firstenddate
        timespans = []
        for i in range(len(bid)):
            # Make a list of enddates from the enddate and the firstenddate
            timespans.append((enddate[i].replace(hour=12, minute=0)
                - firstenddate[i].replace(hour=12, minute=0)).days + 1)
        if np.all(np.array(timespans) <= 0):
            stride = -1
        elif np.all(np.array(timespans) >= 0):
            stride = 1
        else:
            raise ValueError("Timespans are both negative and positive")
        maxtimespan = np.max(np.abs(timespans))

        # A 2-d array of size (bidlen, maxtimespan) is made for the enddates,
        # lats, lons, and begindates, and then this is reshaped to a 1-d array
        enddatearr = np.full((len(bid), maxtimespan), np.nan,
                             dtype='datetime64[s]')
        elatarr = np.full((len(bid), maxtimespan), np.nan)
        elonarr = np.full((len(bid), maxtimespan), np.nan)
        begindatearr = np.full((len(bid), maxtimespan), np.nan,
                               dtype='datetime64[s]')
        perarr = np.full((len(bid), maxtimespan), -1)
        for i in range(len(bid)):
            for j in range(0, timespans[i], stride):
                nenddate = enddate[i].replace(hour=12, minute=0) - timedelta(days=j)
                bd, p = begindateandper(nenddate, period,
                                        forwardtrack=forwardtrack)
                enddatearr[i, j] = np.datetime64(nenddate)
                begindatearr[i, j] = np.datetime64(bd)
                perarr[i, j] = p
                elatarr[i, j] = elat[i]
                elonarr[i, j] = elon[i]
        enddate_1d = enddatearr.ravel()
        elat_1d = elatarr.ravel()
        elon_1d = elonarr.ravel()
        begindate_1d = begindatearr.ravel()
        per_1d = perarr.ravel()

    # Need these 1d arrays of time to be in python datetime format, this
    # is ugly code - and we still need numpy arrays, not lists
    enddate_1d = np.array([convtodt(x) for x in enddate_1d])
    begindate_1d = np.array([convtodt(x) for x in begindate_1d])

    # Find the backtracking name if it is not specified already
    if not outname.endswith('.nc'):
        outname, projname = create_bt_name(outname, firstenddate, enddate, bid,
                                           json_data, projname,
                                           forwardtrack=forwardtrack,
                                           reproc=reproc)

    # Check if the file already exists and then only run if force is set
    if os.path.isfile(outname):
        if force:
            print("File {} exists and force is set, overwriting."
                  "".format(outname))
        else:
            print("Force is not set and the file {} exists, exiting"
                  "".format(outname))
            sys.exit(0)

    btoutput = backtrack_traj(enddate_1d, elon_1d, elat_1d, begindate_1d,
                              per_1d,
                              forwardtrack=forwardtrack,
                              hemi=hemi,
                              iceconc_v=iceconc_v,
                              icedrift_v=icedrift_v,
                              maxskip=maxskip,
                              ow_threshold=ow_threshold,
                              ow_threshold_mid=ow_threshold_mid,
                              ow_repeat=ow_repeat,
                              ex_lmask=ex_lmask,
                              ex_lmask_swath=ex_lmask_swath,
                              verbose=verbose)

    # Working out the time coverage
    maxdate = np.nanmax(btoutput['enddates'])
    mindate = np.nanmin(btoutput['begindates'])
    duration = maxdate - mindate

    # Formatting the arrays
    dimdata = (len(bid), maxtimespan, btoutput['max_timestep'])
    lona = np.full(dimdata, np.nan, dtype=float)
    lata = np.full(dimdata, np.nan, dtype=float)
    flaga = np.full(dimdata, -1, dtype=int)

    lona = btoutput['traj_lons'].reshape(dimdata)
    lata = btoutput['traj_lats'].reshape(dimdata)
    flagarr = btoutput['flag'].reshape(dimdata)

    # Masking NaNs in the lat/lon arrays
    lonarr = np.ma.masked_invalid(lona)
    latarr = np.ma.masked_invalid(lata)

    # Converting the ID names to char arrays with max length 16
    bidpad = [idstr.ljust(16) if len(idstr) <= 16 else idstr[:16]
              for idstr in bid]
    idarr = np.array([[char for char in idstr] for idstr in bidpad], dtype='S1')

    # Skip information
    max_skips_conc = 0
    max_skips_drift = 0
    allconcdir = []
    alldriftdir = []
    if btoutput['max_skips']['conc'] > max_skips_conc:
        max_skips_conc = btoutput['max_skips']['conc']
    if btoutput['max_skips']['drift'] > max_skips_drift:
        max_skips_drift = btoutput['max_skips']['drift']

    # Conc and drift directories
    allconcdir.extend(btoutput['concdir'])
    alldriftdir.extend(btoutput['driftdir'])
    concdirs = ','.join(set(stripdate(allconcdir)))
    driftdirs = ','.join(set(stripdate(alldriftdir)))

    # Checking flag values
    flagset = np.unique(flagarr)
    if ow_threshold_mid is not None:
        allowed = np.array([-1, 0, 1])
    else:
        allowed = np.array([-1, 0])
    for fl in flagset:
        if fl not in allowed:
            print("ERROR. The set of allowed status flags is {} and {} was found!".format(allowed, fl))

    # Writing out the NetCDF file
    with Dataset(outname, 'w', format='NETCDF4') as outf:

        # Setting up the dimensions
        dimidstrlen = outf.createDimension('idstrlen', 16)
        dimid = outf.createDimension('id', len(bid))
        dimtimestep = outf.createDimension('timestep', btoutput['max_timestep'])
        dimtrajenddate = outf.createDimension('trajenddate', maxtimespan)
        lldims = ('id', 'trajenddate', 'timestep')

        # ID number
        idnumvar = outf.createVariable('idnum', 'i4', dimensions=('id'))
        idnumvar.standard_name = 'platform_id'
        idnumvar.long_name = 'Ice station ID number'
        idnumvar[:] = [x for x in range(len(bid))]

        # ID name
        idnamevar = outf.createVariable('idname', 'S1', dimensions=('id',
                                                                    'idstrlen'))
        idnamevar.standard_name = 'platform_name'
        idnamevar.long_name = 'Ice station name'
        idnamevar[:] = idarr

        # Timestep
        timestepvar = outf.createVariable(varname='timestep',
                                          dimensions=('timestep'),
                                          datatype='int32')
        timestepvar.long_name = 'Timestep of days backwards in time'
        timestepvar.units = 'day'
        timestepvar[:] = range(btoutput['max_timestep'])

        # Trajectory end date
        timeunits = 'seconds since 1970-01-01 00:00'
        enddatevar = outf.createVariable(varname='trajenddate',
                                         dimensions=('id', 'trajenddate'),
                                         datatype='float64')
        enddatevar.standard_name = 'time'
        enddatevar.long_name = 'End datetime of the trajectory'
        enddatevar.units = timeunits
        enddatevar[:] = date2num(enddate_1d.reshape((len(bid), maxtimespan)),
                                 units=timeunits, calendar='standard')

        # Longitude
        lonvar = outf.createVariable('longitude', 'f4', dimensions=lldims)
        lonvar.standard_name = 'longitude'
        lonvar.long_name = 'Longitude'
        lonvar.units = 'degrees_east'
        lonvar[:] = lonarr

        # Latitude
        latvar = outf.createVariable('latitude', 'f4', dimensions=lldims)
        latvar.standard_name = 'latitude'
        latvar.long_name = 'Latitude'
        latvar.units = 'degrees_north'
        latvar[:] = latarr

        # Status flag
        flag = outf.createVariable('status_flag', 'i4', dimensions=lldims)
        flag.standard_name = 'status_flag'
        flag.long_name = 'Status flag'
        if ow_threshold_mid is not None:
            flag.flag_values = "-1, 0, 1"
            flag.flag_meanings = "No data, Nominal trajectory above sea-ice concentration {}% (sic_threshold), Nominal trajectory between sea-ice concentration {}% (sic_threshold) and {}% (sic_threshold_flag)".format(ow_threshold, ow_threshold, ow_threshold_mid)
        else:
            flag.flag_values = "-1, 0"
            flag.flag_meanings = "No data, Nominal trajectory above sea-ice concentration {}% (sic_threshold)".format(ow_threshold)
        flag[:] = flagarr

        # Global metadata
        outf.title = "Backtracked trajectories of on-ice positions"
        outf.summary = "Backtracked trajectories of on-ice positions. With each daily timestep backwards in time, the position is shifted according to the sea-ice drift from that date and location, and compared to the sea-ice concentration from that date and location. Trajectories are ended when the calculated position is in a region of sea-ice concentration less than the threshold limit for a certain consecutive number of days."
        #outf.comment
        outf.creator_name = "MET Norway"
        outf.creator_url = "https://www.met.no"
        outf.geospatial_lat_min = np.nanmin(latarr)
        outf.geospatial_lat_max = np.nanmax(latarr)
        outf.geospatial_lon_min = np.nanmin(lonarr)
        outf.geospatial_lon_max = np.nanmax(lonarr)
        outf.geospatial_vertical_min = 0
        outf.geospatial_vertical_max = 0
        outf.history = "Created {:%Y%m%dT%H%M%S}".format(datetime.now())
        outf.date_created = "{:%Y%m%dT%H%M%S}".format(datetime.now())
        outf.institution = "MET Norway"
        outf.contributor_name = "Emily Down, Thomas Lavergne, Signe Aaboe"
        outf.contributor_role = "PrincipalInvestigator,author,author"
        outf.license = "Data and products are licensed under Creative Commons 4.0 BY International (http://creativecommons.org/licenses/by/4.0/)"
        outf.project = "SUDARCO"
        outf.standard_name_vocabulary = "CF Standard Name Table (Version 83, 17 October 2023)"
        outf.keywords = "GCMDSK:Earth Science > Cryosphere > Sea Ice, GCMDSK:Earth Science > Oceans > Sea Ice, GCMDPROV: Government Agencies-non-US > Norway > NO/MET > Norwegian Meteorological Institute"
        outf.keywords_vocabulary = "GCMDSK:GCMD Science Keywords:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/sciencekeywords, GCMDPROV:GCMD Providers:https://gcmd.earthdata.nasa.gov/kms/concepts/concept_scheme/providers"
        outf.source = "Sea-ice Drift: {} \n Sea-ice Concentration: {}".format(idr[icedrift_v]['ref'].format("{:%B %Y}".format(datetime.now())), ic[iceconc_v]['ref'].format("{:%B %Y}".format(datetime.now())))        
        outf.time_coverage_start = "{:%Y%m%dT%H%M%S}".format(mindate)
        outf.time_coverage_end = "{:%Y%m%dT%H%M%S}".format(maxdate)
        outf.time_coverage_duration = "{} days".format(duration.days)
        outf.dataset_name = projname
        outf.sic_threshold = 'Trajectories closed when below {}% ice concentration'.format(ow_threshold)
        if ow_threshold_mid is not None:
            outf.sic_threshold_flag = 'Trajectories flagged (status_flag = 1) when below {}% ice concentration'.format(ow_threshold_mid)            
        outf.ow_repeat = "Trajectories closed after {} days below threshold ice concentration (sic_threshold)".format(ow_repeat)
        if ex_lmask is not None:
            outf.trimmed_over_land = 'True (trimmed over land and coastal regions by a landmask)' 
        else:
            outf.trimmed_over_land = 'False (not trimmed over land and coastal regions by a landmask)'
        perlen = period[:-1]
        if period.endswith('d'):
            perunit = 'days'
        elif period.endswith('w'):
            perunit = 'weeks'
        elif period.endswith('y'):
            perunit = 'years'
        outf.max_backtracking_period = '{} {} (maximum time period over which positions are tracked)'.format(perlen, perunit)
        outf.max_gap_days_drift = str(max_skips_drift)
        outf.max_gap_days_conc = str(max_skips_conc)
        #outf.conc_dir = concdirs
        #outf.drift_dir = driftdirs
        outf.Conventions = 'CF-1.7,ACDD-1.3'
    if verbose:
        print("\nOutput file ready in {}\n".format(outname))

    return outname
    

def main():

    args = parse_args()
    bid = args.bid
    enddate = args.enddate
    firstenddate = args.firstenddate
    elon = args.elon
    elat = args.elat
    json_file = args.json_file
    json_data = args.json_data
    outname = args.outname
    period = args.period
    forwardtrack = args.forwardtrack
    projname = args.projname
    ow_threshold = float(args.ow_threshold)
    if args.ow_threshold_mid is not None:
        ow_threshold_mid = float(args.ow_threshold_mid)
    else:
        ow_threshold_mid = None
    ow_repeat = int(args.ow_repeat)
    iceconc_v = args.iceconc_v
    icedrift_v = args.icedrift_v
    reproc = args.reproc
    maxskip = int(args.maxskip)
    landmask_fp = args.landmask_fp
    landmask_var = args.landmask_var
    landmask_val = int(args.landmask_val)
    verbose = args.verbose
    force = args.force


    bt_out = backtrack_loc(bid=bid, enddate=enddate, firstenddate=firstenddate,
                           elon=elon, elat=elat, json_file=json_file,
                           json_data=json_data, outname=outname, period=period,
                           forwardtrack=forwardtrack, projname=projname,
                           ow_threshold=ow_threshold,
                           ow_threshold_mid=ow_threshold_mid,
                           ow_repeat=ow_repeat, iceconc_v=iceconc_v,
                           icedrift_v=icedrift_v, reproc=reproc,
                           maxskip=maxskip, landmask_fp=landmask_fp,
                           landmark_var=landmask_var,
                           landmask_val=landmask_val, verbose=verbose,
                           force=force)

if __name__ == '__main__':

    main()
