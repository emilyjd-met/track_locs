'''Validation functions'''

import os
import sys
import math
import random
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import xarray as xr
import pandas as pd
from backtrack_loc import backtrack_loc
from plot_trajectories import plot_traj


# Function to read data from a trajectory file into an array
def read_trajfile(trajfile):

    trajdata = {}

    with Dataset(trajfile, 'r') as dataset:

        # Dimensions
        trajdata['dim_id'] = dataset.dimensions['id'].size
        trajdata['dim_timestep'] = dataset.dimensions['timestep'].size
        trajdata['dim_trajenddate'] = dataset.dimensions['trajenddate'].size

        # Getting the IDs
        trajdata['idnames'] = []
        for i in range(trajdata['dim_id']):
            trajdata['idnames'].append(''.join([x.decode('utf-8') for x in
                                                list(dataset.variables['idname'][i, :])]))

        # Getting the data
        trajdata['idnum'] = dataset.variables['idnum'][:]
        trajdata['timestep'] = dataset.variables['timestep'][:]
        trajdata['trajenddate'] = dataset.variables['trajenddate'][:]
        trajdata['latitude'] = dataset.variables['latitude'][:]
        trajdata['longitude'] = dataset.variables['longitude'][:]
        trajdata['status_flag'] = dataset.variables['status_flag'][:]

        # Dictionary of ids
        trajdata['iddict'] = dict(zip(trajdata['idnames'], trajdata['idnum']))

    return trajdata


def validate_traj(newtrajfile, knowntrajfile):
    '''Function to validate a trajectory against a saved file (i.e. a
    check that nothing has changed while debugging'''

    # Read in the data
    knowntraj = read_trajfile(knowntrajfile)
    newtraj = read_trajfile(newtrajfile)

    latdiffs = []
    londiffs = []

    # Assembling the lat/lon differences
    for idn in knowntraj['idnames']:
        # For each named station in the known output, check that this named station exists in the new output
        if idn in newtraj['idnames']:
            print("Comparing data for station {}".format(idn))
            knownidnum = knowntraj['iddict'][idn]
            newidnum = newtraj['iddict'][idn]

            # Now cyling over the trajectory end dates
            for i, ted in enumerate(knowntraj['trajenddate'][knownidnum, :]):
                if ted in newtraj['trajenddate'][newidnum, :]:
                    j = np.where(newtraj['trajenddate'][newidnum, :] == ted)[0][0]

                    # Selecting the minimum of the timesteps from the two datasets 
                    tsdim = min(knowntraj['dim_timestep'], newtraj['dim_timestep'])

                    # Finding the lat/lon data
                    knownlats = knowntraj['latitude'][knownidnum, i, :tsdim-1]
                    knownlons = knowntraj['longitude'][knownidnum, i, :tsdim-1]
                    newlats = newtraj['latitude'][knownidnum, j, :tsdim-1]
                    newlons = newtraj['longitude'][knownidnum, j, :tsdim -1]

                    # Differencing these datasets
                    latdf = knownlats - newlats
                    londf = knownlats - newlats

                    # Removing the masked values
                    latdf = latdf[latdf.mask == False]
                    londf = londf[londf.mask == False]

                    # Appending to the overall arrays
                    latdiffs.extend(list(latdf))
                    londiffs.extend(list(londf))

    # Finding some statistics
    latdiffs = np.array(latdiffs)
    londiffs = np.array(londiffs)
    print("")
    print("LATITUDE: Mean {}; Std. Dev {}".format(np.nanmean(latdiffs), np.nanstd(latdiffs)))
    print("LONGITUDE: Mean {}; Std. Dev {}".format(np.nanmean(londiffs), np.nanstd(londiffs)))


def compute_distance(lat0, lon0, lat1, lon1):

    rearth = 6371.
    deg_to_rad = math.pi / 180.
    lat0 *= deg_to_rad
    lon0 *= deg_to_rad
    lat1 *= deg_to_rad
    lon1 *= deg_to_rad
    dlat = 0.5 * (lat1 - lat0)
    dlon = 0.5 * (lon1 - lon0)
    a = (math.sin(dlat) * math.sin(dlat)) + (
        math.cos(lat0) * math.cos(lat1) * math.sin(dlon) * math.sin(dlon))
    dist = rearth * 2. * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return dist


def forward_backward_test(tname, tdate, tlon, tlat, tper, outdir='.'):
    '''A test suggested by Thomas to run forwards tracking followed by
    backwards tracking to then check the difference'''

    # Forward tracking
    fw_out = backtrack_loc(bid=tname, enddate=tdate, elon=tlon, elat=tlat,
                           period=tper, outname=outdir, forwardtrack=True,
                           force=True)

    # Finding the end lat/lon and date from the forwards tracking
    with Dataset(fw_out) as ds:
        elon = ds['longitude'][0, 0, -1]
        elat = ds['latitude'][0, 0, -1]
        print("End longitude, latitude : {:.2f}, {:.2f} ".format(elon, elat))
        edate = datetime.strftime(
            datetime.strptime(tdate, '%Y%m%d%H%M')
            + timedelta(days=int(tper[0: -1])), '%Y%m%d%H%M')

    # Run the backwards tracking on the end lat/lon and time
    # Backward tracking
    bk_out = backtrack_loc(bid=tname, enddate=edate, elon=str(elon),
                           elat=str(elat), period=tper, outname=outdir,
                           force=True)

    # Finding the equivalent start lat/lon after the forwards-and-backtracking
    with Dataset(bk_out) as ds:
        blon = ds['longitude'][0, 0, -1]
        blat = ds['latitude'][0, 0, -1]
    print("Beginning longitude, latitude : {:.2f}, {:.2f} ".format(blon, blat))

    # Plotting
    plot_traj(fw_out, region='ease-nh', colmode='timestep', output=outdir)
    plot_traj(bk_out, region='ease-nh', colmode='timestep', output=outdir)

    # Printing results
    print("Start of the forward trajectory: {:.2f}, {:.2f}".format(float(tlat),
                                                                   float(tlon)))
    print("End of the backtrack trajectory: {:.2f}, {:.2f}".format(float(blat),
                                                                   float(blon)))
    print("Distance between these: {:.2f} km".format(
        compute_distance(float(tlat), float(tlon), float(blat), float(blon))))


def nearest(items, pivot):
    '''Function to return nearest date to pandas datetime column
        From https://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date
    '''
    return min(items, key=lambda x: abs(x - pivot))


def pos_at_interval(buoyfile, period, sdate='first'):
    '''Accepts the path to a NetCDF buoy dataset and a time interval and returns the start and end lat/lons
        buoyfile: Path to NetCDF file
        period: Can be in 'd' (days) or 'y' (years) and should be an integer
        sdate: Either 'first' for the first datetime in the buoyfile, or a specified datetime
    '''
    ds = xr.open_dataset(buoyfile)
    df = ds.to_dataframe()
    #print(ds)
    #print(df)

    # Finding the initial datetime
    if sdate == 'first':
        firstdatetime = df.index[0]
        firstdatetime_py = firstdatetime.to_pydatetime()
        startdatetime_exact = firstdatetime_py.replace(hour=12, minute=0, second=0, microsecond=0)
    else:
        startdatetime_exact = sdate
    startdatetime_real = nearest(df.index, startdatetime_exact)

    # And the initial lat/lon
    startlat = df[df.index == startdatetime_real]['lat'][0]
    startlon = df[df.index == startdatetime_real]['lon'][0]

    # Finding the enddate
    period_unit = period[-1]
    period_val = int(period[0: -1])
    if period_unit == 'w':
        period_val *= 7
        period_unit = 'd'
    if period_unit == 'd':
        enddatetime_exact = startdatetime_exact + timedelta(days=period_val)
    elif period_unit == 'y':
        enddatetime_exact = startdatetime_exact + relativedelta(years=period_val)
    else:
        raise ValueError("Unrecognised time unit {}".format(period_unit))
    enddatetime_real = nearest(df.index, enddatetime_exact)

    # Checking that the end date is reasonably close
    if abs((enddatetime_real - enddatetime_exact).days) > 2:
        raise ValueError("No end datetime found within 2 days of requested")

    # Find the final lat/lon
    endlat = df[df.index == enddatetime_real]['lat'][0]
    endlon = df[df.index == enddatetime_real]['lon'][0]

    buoydict = {}
    buoydict['startlat'] = startlat
    buoydict['startlon'] = startlon
    buoydict['startdatetime_exact'] = startdatetime_exact
    buoydict['startdatetime_real'] = startdatetime_real
    buoydict['endlat'] = endlat
    buoydict['endlon'] = endlon
    buoydict['enddatetime_exact'] = enddatetime_exact
    buoydict['enddatetime_real'] = enddatetime_real
    buoydict['buoyname'] = '{}_{}'.format(ds.network, ds.id)
    buoydict['period'] = period

    return buoydict


def per_to_days(period, startdate=None):

    period_unit = period[-1]
    period_val = int(period[0: -1])
    if period_unit == 'd':
        period_days = period_val
    elif period_unit == 'w':
        period_days = period_val * 7
    elif period_unit == 'y':
        enddate = startdate + relativedelta(years=period_val)
        period_days = (enddate - startdate).days

    return(period_days)


def compare_buoy_traj(buoyfile, interval, workdir='./workdir'):
    '''Compare the forward tracking to a single buoy'''

    # Make the working directory if it doesn't exist
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    # First retrieve the start and end lat/lons from the buoy file for the required interval
    try:
        bd = pos_at_interval(buoyfile, interval)
    except:
        offset = "buoy fail"
        return offset
        #raise ValueError("No lat/lons possible from buoy file {} with interval {}".format(buoyfile, interval))

    #print(bd)

    # Run forward tracking from the start lat/lon and date
    #per = (bd['enddatetime_exact'] - bd['startdatetime_exact']).days
    # Need to get date as YYYYmmddHHMM format **********
    try:
        bdse = bd['startdatetime_exact']
        if isinstance(bdse, str):
            try:
                dt = datetime.strptime(bdse, '%Y%m%d%H%M')
            except ValueError:
                dt = datetime.strptime(bdse, '%Y-%m-%d %H:%M:%S')
        elif isinstance(bdse, datetime):
            dt = bdse
        else:
            raise ValueError('bdse must be a string or a datetime object')
        refdate = datetime.strftime(bdse, '%Y%m%d%H%M')
    except:
        offset = "refdate fail"
        return offset

    try:
        ft_bm = backtrack_loc(bid=bd['buoyname'], enddate=refdate,
                              elon=str(bd['startlon']),
                              elat=str(bd['startlat']),
                              period=bd['period'], outname=workdir,
                              forwardtrack=True,
                              force=True)
    except:
        offset = "tracking fail"
        return offset

    # Retrieve the end lat/lons
    with Dataset(ft_bm) as ds:

        # Checking that the trajectory ran as long as the buoy trajectory did
        trajp = ds.time_coverage_duration.replace(" days", "d").replace(" years", "y")
        traj_len = per_to_days(trajp, startdate=bd['startdatetime_exact'])
        buoy_traj_len = per_to_days(bd['period'], startdate=bd['startdatetime_exact'])
        if abs(traj_len - buoy_traj_len) > 2:
            #raise ValueError("This trajectory did not run as long as required, passing")
            offset = "traj len fail"
            return offset

        try:
            elon = ds['longitude'][0, 0, -1]
            elat = ds['latitude'][0, 0, -1]
        except:
            offset = "get lat/lon fail"
            return offset

    # And calculate the distance
    try:
        offset = compute_distance(bd['endlat'], bd['endlon'], elat, elon)
    except:
        offset = "distance compute fail"
        return offset

    return offset


def filter_buoy_loc_time(bfile):
    '''Some rough filtering for the goodness of the buoy files'''

    goodbuoy = True
    with Dataset(bfile, 'r') as dataset:

        if float(dataset.geospatial_lat_max) < 78:
            goodbuoy = False
        if float(dataset.geospatial_lon_min) < 160.:
            goodbuoy = False
        if float(dataset.geospatial_lon_max) > 280.:
            goodbuoy = False
        # For this rough filtering we want the whole trajectory within
        # this period
        tcs = datetime.strptime(dataset.time_coverage_start, '%Y-%m-%dT%H:%M:%SZ')
        tce = datetime.strptime(dataset.time_coverage_end, '%Y-%m-%dT%H:%M:%SZ')
        if tcs < datetime(1991, 10, 1):
            goodbuoy = False
        if tce > datetime(2020, 12, 31):
            goodbuoy = False

    return goodbuoy


# Randomly sample files from a directory
# From https://stackoverflow.com/questions/36949029/most-efficient-way-to-randomly-sample-directory-in-python
def randomfiles(mydir, samplenum):

    files = sorted(os.listdir(mydir))
    file_count = len(files)
    #print("file_count = ", file_count)
    random_nums = random.sample(range(0, file_count), samplenum)
    random_files = [os.path.join(mydir, files[x]) for x in random_nums]
    #print(random_files)

    return random_files


def run_buoy_comparison(buoydir, samplenum, period, workdir='./workdir'):

    rfiles = randomfiles(buoydir, samplenum)

    # Do I want some rough filtering?

    buoy_fail = 0
    tracking_fail = 0
    traj_len_fail = 0
    distance_compute_fail = 0
    offsets = []

    for i, rf in enumerate(rfiles):
#        try:
        print(i, rf)
        offset = compare_buoy_traj(rf, period, workdir=workdir)
        if offset == 'buoy fail':
            buoy_fail += 1
        elif offset == 'tracking fail':
            tracking_fail += 1
        elif offset == 'traj len fail':
            traj_len_fail += 1
        elif offset == 'distance compute fail':
            distance_compute_fail += 1
        else:
            offsets.append(offset)
#        except:
#            pass

    print("====================================")
    print("OFFSETS: ", offsets)

    #print("{}/{} files succeeded, average distance {}".format(len(offsets), samplenum, np.array(offsets).mean()))
    print("{}/{} files succeeded".format(len(offsets), samplenum))
    print("{}/{} files failed due to problem reading required length trajectory from buoy file".format(
        buoy_fail, samplenum))
    print("{}/{} files failed due to problem computing trajectory".format(
        tracking_fail, samplenum))
    print("{}/{} files failed due to problem reading required length trajectory from ice drift file".format(
        traj_len_fail, samplenum))
    print("{}/{} files failed due to problem computing distance".format(
        distance_compute_fail, samplenum))
    print("Average distance {:.2f} km".format(np.array(offsets).mean()))
