'''Plot tracked trajectories
'''

import os.path
import json
import argparse
from argparse import RawDescriptionHelpFormatter
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.ops import cascaded_union
from pyresample import utils, AreaDefinition
from netCDF4 import Dataset
import cmocean

from grid_info import region_params, valid_regions
from bg_params import bg_plot_setup

col_land = 'Tan' # grey
col_sea = '#04629a'
col_iceconc50 = '#9b9b9b'
col_iceconc100 = 'white'
col_sea_trans = (0.01568627450980392, 0.3843137254901961, 0.6039215686274509, 0.9)

titlefont = 12
labelfont = 8
#copyfontsize = 6
copyfont = 6
deffont = 10
plt.rc('font', size=deffont) # controls default text sizes

mie_file = os.path.join(os.path.dirname(__file__), 'min_ice_extent.json')
min_sie_ver = 'cdr-v2p2'

# Defining the matplotlib colours by https://sashamaps.net/docs/resources/20-colors/
# 100%
qblack = '#000000'
qgrey = '#a9a9a9'
qwhite = '#ffffff'
qyellow = '#ffe119'
qblue = '#4363d8'
# 99.99%
qmaroon = '#800000'
qorange = '#f58231'
qnavy = '#000075'
qlavender = '#dcbeff'
# 99%
qbrown = '#9A6324'
qred = '#e6194B'
qpink = '#fabed4'
qbeige = '#fffac8'
qteal = '#469990'
qgreen = '#3cb44b'
qcyan = '#42d4f4'
qmint = '#aaffc3'
qmagenta = '#f032e6'
# 95%
qapricot = '#ffd8b1'
qolive = '#808000'
qlime = '#bfef45'
qpurple = '#911eb4'
cmap_list20 = [qmaroon, qnavy, qorange, qgreen, qmagenta, qyellow, qteal, qbrown, qcyan, qred, qlavender, qmint, qbeige, qpink, qolive, qpurple, qlime, qapricot, qgrey, qblue]
#cmap_list20 = [qblue, qmaroon, qnavy, qorange, qgreen, qmagenta, qteal, qbrown, qcyan, qred, qyellow, qlavender, qmint, qbeige, qpink, qolive, qpurple, qlime, qapricot, qgrey]
cmap_20 = colors.LinearSegmentedColormap.from_list('cmap20', cmap_list20, 20)

monthcols = {1: (qmaroon, 'Jan'),
             2: (qred, 'Feb'),
             3: (qorange, 'Mar'),
             4: (qyellow, 'Apr'),
             5: (qlime, 'May'),
             6: (qgreen, 'Jun'),
             7: (qteal, 'Jul'),
             8: (qcyan, 'Aug'),
             9: (qblue, 'Sep'),
             10: (qnavy, 'Oct'),
             11: (qpurple, 'Nov'),
             12: (qmagenta, 'Dec')}
#mocolmap = colors.ListedColormap([monthcols[x][0] for x in monthcols.keys()])
molabs = [monthcols[x][1] for x in monthcols.keys()]


def colbar_mo_discrete():

    #cmap_orig = cm.get_cmap('cmo.thermal', 12)
    cmap_orig = cm.get_cmap('viridis', 12)
    colors_i = np.linspace(0.0, 1.0, 12)
    c12 = cmap_orig(colors_i)
    cmap_mo = colors.LinearSegmentedColormap.from_list('listmon', c12, 12)

    return cmap_mo


def parse_args():

    valid_colmode = ['id', 'enddate', 'timestep', 'myi', 'month']

    p = argparse.ArgumentParser("traj_plot",
                                formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('trajinp',
                   help="NetCDF file with the trajectories")
    p.add_argument('-r', '--region', required=False, choices=valid_regions,
                   default='ease-nh-wide',
                   help="Region to plot")
    p.add_argument('-o', '--output', required=False, default='.',
                   help="Either an output directory or a full filepath")
    p.add_argument('-l', '--label_traj', action='store_true', default=False,
                   help="Label trajectories with station ID")
    p.add_argument('-l2', '--label_traj', action='store_true', default=False,
                   help="Label trajectories directly with station ID (no key)")
    p.add_argument('-c', '--colmode', required=False, default='timestep',
                   choices=valid_colmode,
                   help="Colouring mode for trajectory, default is 'timestep', "
                   "valid choices are {}".format(valid_colmode))
    p.add_argument('-t', '--title', required=False, default='',
                   help="Title for the plot. Leaving this blank sets the "
                   "title to the dataset title. Setting to 'notitle' leaves "
                   " this blank")
    p.add_argument('-s', '--stride', required=False, default='1d',
                   help="Stride of trajectory enddates to plot, e.g. 1m, "
                   "7d, default is 1d")
    p.add_argument('-b', '--bathymetry', action='store_true', default=False,
                   help="Plot bathymetry lines")
    p.add_argument('-g', '--latlongrid', action='store_true', default=False,
                   help="Plot a lat/lon grid?")
    p.add_argument('-g2', '--latlongrid2', action='store_true', default=False,
                   help="Plot a lat/lon grid with labels?")
    p.add_argument('-bgf', '--bg_file', required=False, default=None,
                   help="Filepath to the background file to plot")
    p.add_argument('-bgv', '--bg_var', required=False, default=None,
                   help="Variable name of the background to plot")
    p.add_argument('-bgi', '--bg_ind', required=False, default=None,
                   help="Index of the frame of background file to plot")
    p.add_argument('-bgd', '--bg_date', required=False, default=None,
                   help="Date label of the background frame to display")
    p.add_argument('-lw', '--linew', required=False, default=1.0,
                   help="Line width for the trajectories (default=1)")
    p.add_argument('-icf', '--ignore_mid_conc', required=False, default=False,
                   help="Don't flag mid concentration part of the trajectory "
                   "with transparency even if flags are set")
    p.add_argument('-j', '--plot_jpg', action='store_true', default=False,
                   help="Plot a jpeg instead of a png")
    p.add_argument('-v', '--verbose', action='store_true', default=False,
                   help="Print extra information while running the code")

    args = p.parse_args()

    return args


def read_traj(ncfile):

    trajdata = {}

    # Reading from the NetCDF file
    with Dataset(ncfile, 'r') as dataset:

        # Finding the dimensions
        dims = dataset.variables['longitude'][:].shape
        trajdata['dimid'], trajdata['dimreferencedate'], trajdata['dimtimestep'] = dims

        # idnum and idname are size (dimid)
        trajdata['idnum'] = dataset.variables['idnum'][:]
        trajdata['idname'] = [''.join([x.decode('utf8') for x in y]).strip()
                              for y in dataset.variables['idname'][:]]
        # timestep is size (dimstimestep)
        trajdata['timestep'] = dataset.variables['timestep'][:]
        # referencedate is dimensions (dimid, dimreferencedate)
        trajdata['referencedatesec'] = dataset.variables['referencedate'][:, :]
        trajdata['timeunits'] = dataset.variables['referencedate'].units
        ted = []
        for i in range(trajdata['dimid']):
            ted.append([datetime.fromtimestamp(t)
                        for t in trajdata['referencedatesec'][i, :]])
        trajdata['referencedate'] = np.array(ted)
        trajdata['lat'] = dataset.variables['latitude'][:]
        trajdata['lon'] = dataset.variables['longitude'][:]
        trajdata['flag'] = dataset.variables['status_flag'][:]
        trajdata['dataset_name'] = dataset.dataset_name

    return trajdata


def filter_by_stride(trajdata, stride):
    '''Filter the data by trajectory end date according to the stride'''

    strideunits = stride[-1:]
    if strideunits not in ['m', 'd']:
        raise ValueError("Unknown stride units: {}".format(strideunits))
    strideval = int(stride[0:-1])

    aenddates = []
    alons = []
    alats = []
    aflags = []
    for i in range(len(trajdata['idname'])):
        nenddates = []
        nlats = []
        nlons = []
        nflags = []
        nextstride = trajdata['referencedate'][i, 0]
        for j, ed in enumerate(trajdata['referencedate'][i, :]):
            if ed == nextstride:
                if strideunits == 'm':
                    nextstride = nextstride - relativedelta(months=strideval)
                elif strideunits == 'd':
                    nextstride = nextstride - timedelta(days=strideval)
                nenddates.append(ed)
                nlons.append(trajdata['lon'][i, j, :])
                nlats.append(trajdata['lat'][i, j, :])
                nflags.append(trajdata['flag'][i, j, :])
        aenddates.append(nenddates)
        alons.append(nlons)
        alats.append(nlats)
        aflags.append(nflags)

    trajdata['referencedate'] = np.array(aenddates)
    trajdata['lon'] = np.array(alons)
    trajdata['lat'] = np.array(alats)
    trajdata['flag'] = np.array(aflags)

    return trajdata


def find_mie_date(zdate, mie_list):
    '''Find the minimum sea ice extent date for a given date given the list'''

    if isinstance(zdate, datetime):
        zdate = zdate.date()

    try:
        mie_date = datetime.strptime(mie_list[str(zdate.year)], '%Y%m%d').date()
        if mie_date <= zdate:
            return datetime.strptime('{}{}'.format(mie_date.strftime('%Y%m%d'),
                                                   '12'), '%Y%m%d%H')
        mie_date = datetime.strptime(mie_list[str(zdate.year - 1)],
                                     '%Y%m%d').date()
        return datetime.strptime('{}{}'.format(mie_date.strftime('%Y%m%d'),
                                               '12'), '%Y%m%d%H')
    except:
        raise ValueError("Problem with retrieving the minimum ice extent date "
                         "for date {}".format(zdate))


def nc_read(ncfile, var, bg_ind=None, skip=None):

    ncdata = {}

    if isinstance(var, list):
        v = var[0]
    else:
        v = var

    # Reading from the NetCDF file
    with Dataset(ncfile, 'r') as dataset:
        try:
            ncdata['grid_mapping'] = dataset.variables[v].__dict__['grid_mapping']
        except:
            ncdata['grid_mapping'] = 'crs'
        gridmap_vars = dataset.variables[ncdata['grid_mapping']].__dict__
        if 'proj4_string' in gridmap_vars:
            proj4var = 'proj4_string'
        elif 'proj4' in gridmap_vars:
            proj4var = 'proj4'
        elif 'proj4text' in gridmap_vars:
            proj4var = 'proj4text'
        else:
            raise ValueError("No proj4 variable found in grid mapping "
                             "variable {}".format(ncdata['grid_mapping']))
        ncdata['proj4_string'] = dataset.variables[ncdata['grid_mapping']].__dict__[proj4var]
        try:
            ncdata['proj_dict'] = utils._proj4.proj4_str_to_dict(ncdata['proj4_string'])
        except:
            ncdata['proj_dict'] = utils.proj4.proj4_str_to_dict(ncdata['proj4_string'])

        if 'xc' in dataset.variables:
            xvar = 'xc'
        elif 'x' in dataset.variables:
            xvar = 'x'
        else:
            raise ValueError("No x-variable found")
        if 'yc' in dataset.variables:
            yvar = 'yc'
        elif 'y' in dataset.variables:
            yvar = 'y'
        else:
            raise ValueError("No y-variable found")
        ncdata['xc'] = dataset[xvar][:]
        ncdata['yc'] = dataset[yvar][:]

        try:
            ncdata['fv'] = dataset.variables[v].__dict__['_FillValue']
        except:
            pass

        if isinstance(var, list):
            varlist = var
        else:
            varlist = [var]
        for item in varlist:
            vardata = dataset[item][:]

            # Extra step for ice age - set the very old ice 11-19 years old
            # as 10 years for 10+ years (very few pixels)
            # No - go for 5+ years
            if item == 'age_of_sea_ice':
                vold = np.logical_and(vardata > 5, vardata <= 19)
                vardata[vold] = 5
                vardata[vardata >= 20] = 0
                vardata[vardata < 1] = 0
                #vmask = np.zeros_like(vardata, dtype=bool)
                #vardata.mask = vmask
                #vardata.mask[vardata == 0] = 1

            # Extra step for ice type so ambig displays nicely between FYI
            # and MYI
            vardata[vardata == 3] = 100
            vardata[vardata == 4] = 3
            vardata[vardata == 100] = 4

            if len(vardata.shape) == 3:
                # Either use bg_ind if provided, or default to 0
                if bg_ind is None:
                    bg_ind = 0
                if skip:
                    vardata = vardata[bg_ind, ::skip, ::skip]
                else:
                    vardata = vardata[bg_ind, :, :]
            else:
                if skip:
                    vardata = vardata[::skip, ::skip]
                else:
                    vardata = vardata[:, :]

            # Landmask test
            #vd = np.zeros_like(vardata, dtype=int)
            #vd[vardata >= 2] = 1
            #vardata = vd

            # This is for reducing flags to the lowest integers
            if var in ['statusflag', 'status_flag', 'flag']:
                vardata = np.asarray(vardata, float) + 100
                uniques = np.unique(vardata)
                uniques = uniques[np.logical_not(np.isnan(uniques))]
                for newval, origval in enumerate(uniques):
                    vardata[vardata == origval] = newval
                ncdata['sf_labs'] = [str(int(u - 100)) for u in uniques]

            # NOTE: Be very careful with the fill value here. Trying to
            # use ncdata['fv'] as the fill_value for a status array such
            # as 'flag' means that flag values of 0 (i.e. nominal) are
            # masked out
            if 'fv' in ncdata.keys() and item not in ['statusflag',
                                                      'status_flag', 'flag']:
                ncdata[item] = ma.array(vardata, fill_value=ncdata['fv'])
            else:
                ncdata[item] = ma.array(vardata)

            if var in ['statusflag', 'status_flag', 'flag']:
                ncdata[item].mask = None
            else:
                ncdata[item].mask = ncdata[item].data == ncdata[item].fill_value

        if 'lon' in dataset.variables:
            lonvar = 'lon'
        elif 'longitude' in dataset.variables:
            lonvar = 'longitude'
        else:
            raise ValueError("No lon-variable found")
        if 'lat' in dataset.variables:
            latvar = 'lat'
        elif 'latitude' in dataset.variables:
            latvar = 'latitude'
        else:
            raise ValueError("No lat-variable found")
        if skip:
            ncdata['lon'] = dataset.variables[lonvar][::skip, ::skip]
            ncdata['lat'] = dataset.variables[latvar][::skip, ::skip]
        else:
            ncdata['lon'] = dataset.variables[lonvar][:]
            ncdata['lat'] = dataset.variables[latvar][:]

        # Try fetching time info
        try:
            try:
                d0 = datetime.strptime(dataset.start_date,'%Y-%m-%d %H:%M:%S')
                d1 = datetime.strptime(dataset.stop_date,'%Y-%m-%d %H:%M:%S')
            except:
                try:
                    d0 = datetime.strptime(dataset.start_date_and_time,'%Y-%m-%dT%H:%M:%SZ')
                    d1 = datetime.strptime(dataset.end_date_and_time,'%Y-%m-%dT%H:%M:%SZ')
                except:
                    d0 = datetime.strptime(dataset.time_coverage_start,'%Y-%m-%dT%H:%M:%SZ')
                    d1 = datetime.strptime(dataset.time_coverage_end,'%Y-%m-%dT%H:%M:%SZ')
            d0_00 = datetime.combine(d0.date(),time(0))
            d1_00 = datetime.combine(d1.date(),time(0))
            ncdata['sdate'] = d0_00
            ncdata['edate'] = d1_00
            ncdata['tspan_hours'] = (d1_00 - d0_00).total_seconds() / (60.*60.)
        except:
            pass

        try:
            ncdata['time'] = dataset.variables['time'][0]
            ncdata['time_bnds0'] = dataset.variables['time_bnds'][0][0]
            ncdata['time_bnds1'] = dataset.variables['time_bnds'][0][1]
        except:
            pass

        # Having read the grid info, if the grid is just called 'crs', set this
        if ncdata['grid_mapping'] == 'crs':
            gmn = dataset.variables['crs'].__dict__['grid_mapping_name']
            if 'lambert' in gmn or 'laea' in gmn:
                ncdata['grid_mapping'] = 'LambertAzimuthalEqualArea'
            if 'stereographic' in gmn or 'polstere' in gmn:
                ncdata['grid_mapping'] = 'Polar_Stereographic_Grid'

    # Grid spacing
    sorted_xc = np.sort(ncdata['xc'])
    sorted_yc = np.sort(ncdata['yc'])
    smallest_xc = sorted_xc[0]
    second_smallest_xc = sorted_xc[1]
    smallest_yc = sorted_yc[0]
    second_smallest_yc = sorted_yc[1]
    ncdata['ax'] = second_smallest_xc - smallest_xc
    ncdata['ay'] = second_smallest_yc - smallest_yc

    # Area definitions and extents
    if abs(float(ncdata['xc'][0])) < 10000.:
        sf = 1000.
    else:
        sf = 1.
    ncdata['area_extent'] = (float(ncdata['xc'][0] * sf),
                             float(ncdata['yc'][0] * sf),
                             float(ncdata['xc'][-1] * sf),
                             float(ncdata['yc'][-1] * sf))
    ncdata['mpl_extent'] = (ncdata['area_extent'][0],
                            ncdata['area_extent'][2],
                            ncdata['area_extent'][3],
                            ncdata['area_extent'][1])
    # Shifting by half a pixel for divergences and convergences, which
    # are calculated between pixels
    ncdata['mpl_extent_divs'] = (ncdata['area_extent'][0]
                                 + (0.5 * sf * ncdata['ax']),
                                 ncdata['area_extent'][2]
                                 + (0.5 * sf * ncdata['ax']),
                                 ncdata['area_extent'][3]
                                 + (0.5 * sf * ncdata['ay']),
                                 ncdata['area_extent'][1]
                                 + (0.5 * sf * ncdata['ay']))
    ncdata['area_def'] = AreaDefinition('data', 'data', 'data',
                                        ncdata['proj_dict'],
                                        ncdata['xc'].shape[0],
                                        ncdata['yc'].shape[0],
                                        ncdata['area_extent'])
    ncdata['data_crs'] = ncdata['area_def'].to_cartopy_crs()

    if ncdata['grid_mapping'] in ['Polar_Stereographic_Grid',
                                  'projection_stere']:
        data_globe = ccrs.Globe(semimajor_axis=ncdata['proj_dict']['a'],
                                semiminor_axis=ncdata['proj_dict']['b'])
        if ncdata['lat'][0, 0] > 0:
            ncdata['data_ccrs'] = ccrs.NorthPolarStereo(central_longitude=-45.0,
                                                        globe=data_globe)
            ncdata['hemi'] = 'nh'
        else:
            ncdata['data_ccrs'] = ccrs.SouthPolarStereo(central_longitude=0.0,
                                                        globe=data_globe)
            ncdata['hemi'] = 'sh'
    elif ncdata['grid_mapping'] in ['LambertAzimuthalEqualArea',
                                    'Lambert_Azimuthal_Equal_Area',
                                    'Lambert_Azimuthal_Grid',
                                    'projection_laea']:
        if ncdata['lat'][0, 0] > 0:
            ncdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=90,
                false_easting=0, false_northing=0)
            ncdata['hemi'] = 'nh'
        else:
            ncdata['data_ccrs'] = ccrs.LambertAzimuthalEqualArea(
                central_longitude=0, central_latitude=-90,
                false_easting=0, false_northing=0)
            ncdata['hemi'] = 'sh'
    else:
        raise ValueError("Unrecognised grid mapping {}".format(ncdata['grid_mapping']))


    return ncdata


def plot_traj(trajinp, region='ease-nh-wide', output='.', label_traj=False,
              label_traj2=False, colmode='timestep', title='', stride='1d',
              bathymetry=False, latlongrid=False, latlongrid2=False,
              bg_file=None, bg_var=None, bg_ind=None, bg_date=None, linew=1.0,
              ignore_mid_conc=False, plot_jpg=False, verbose=False):

    # Finding the region parameters
    rp = region_params(region)
    if rp['lllat'] > 0.:
        hemi = 'nh'
    else:
        hemi = 'sh'
    if region.startswith('ease'):
        gridtype = 'ease'
    elif region.startswith('polstere30'):
        gridtype = 'polstere30'
    elif region.startswith('pol'):
        gridtype = 'polstere'
    else:
        raise ValueError("Unrecognised grid type for region {} - "
                         "expected to start with 'ease' or 'pol'"
                         "".format(region))
    # Define grid based on region
    if gridtype == 'polstere':
        if hemi == 'nh':
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': 90.,
                                 'lat_ts' : 70.,
                                 'lon_0': -45.0,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.NorthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)
        else:
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': -90.,
                                 'lat_ts' : -70.,
                                 'lon_0': 0.,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.SouthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)

    elif gridtype == 'polstere30':
        if hemi == 'nh':
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': 90.,
                                 'lat_ts' : 70.,
                                 'lon_0': 30,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.NorthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)
        else:
            plot_proj4_params = {'proj': 'stere',
                                 'lat_0': -90.,
                                 'lat_ts' : -70.,
                                 'lon_0': 0.,
                                 'a': 6378273,
                                 'b': 6356889.44891}
            plot_globe = ccrs.Globe(semimajor_axis=plot_proj4_params['a'],
                                    semiminor_axis=plot_proj4_params['b'])
            plot_crs = ccrs.SouthPolarStereo(
                central_longitude=plot_proj4_params['lon_0'], globe=plot_globe)

    elif gridtype == 'ease':
        if hemi == 'nh':
            plot_crs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                      central_latitude=90,
                                                      false_easting=0,
                                                      false_northing=0)
        else:
            plot_crs = ccrs.LambertAzimuthalEqualArea(central_longitude=0,
                                                      central_latitude=-90,
                                                      false_easting=0,
                                                      false_northing=0)
    else:
        raise ValueError("Unrecognised region {}".format(region))

    # Reading in background file if given and setting up background
    # quantities for the plot
    if bg_file is not None:
        bg_nc = nc_read(bg_file, bg_var, bg_ind=bg_ind)
        bg_dict = bg_plot_setup(bgvar=bg_var, bgnc=bg_nc)
    else:
        bg_dict = None

    pc = ccrs.PlateCarree()

    # Finding the plot corners
    llx, lly = plot_crs.transform_point(rp['lllon'], rp['lllat'], src_crs=pc)
    urx, ury = plot_crs.transform_point(rp['urlon'], rp['urlat'], src_crs=pc)
#    llx, lly = plot_crs.transform_point(min(rp['lllon'], rp['urlon']), min(rp['lllat'], rp['urlat']), src_crs=pc)
#    urx, ury = plot_crs.transform_point(max(rp['lllon'], rp['urlon']), max(rp['lllat'], rp['urlat']), src_crs=pc)

    # The projection keyword determines how the plot will look
    fig = plt.figure(figsize=(6, 6))
#    fig = plt.figure(figsize=(14, 12))
    ax = plt.axes(projection=plot_crs)

    # Setting the region
# TODO: Maybe switch to using PlateCarree for extents (change in background
# settings file)
#    ax.set_extent([rp['lllon'], rp['urlon'], rp['lllat'], rp['urlat']], crs=pc)
    ax.set_extent([llx, urx, lly, ury], crs=plot_crs)

    # Plotting the background if this is set
#    print("plot_crs = ", plot_crs.__dict__)
#    print("bg_nc['mpl_extent'] = ", bg_nc['mpl_extent'])
#    print("bg_nc['data_ccrs'] = ", bg_nc['data_ccrs'].__dict__)
    if bg_var is not None:
#        bg = ax.imshow(bg_nc[bg_var], extent=(-2000000, 2000000, 2000000, -2000000),
        bg = ax.imshow(bg_nc[bg_var], extent=bg_nc['mpl_extent'],
                       transform=bg_nc['data_ccrs'],
                       cmap=bg_dict['cmap'], vmax=bg_dict['max'],
                       vmin=bg_dict['min'], norm=bg_dict['norm'],
                       interpolation='none', alpha=bg_dict['alpha'])
    #if colbar:
        if bg_date is not None:
            bg_cb_lab = '{} at {}'.format(bg_dict['colbar_label'], bg_date)
        else:
            bg_cb_lab = bg_dict['colbar_label']
        if bg_dict['colbar_type'] == 'discrete':
            cb = plt.colorbar(bg, ticks=bg_dict['cmap_lvl'],
#                              location='right',
#                              orientation='vertical', pad=0.0, shrink=0.4
                              orientation='horizontal', pad=0.0, shrink=0.3
            )
            cb.ax.set_xticklabels(bg_dict['cmap_labs'])
            cb.set_label(bg_cb_lab, fontsize=labelfont)
        else:
            #fig = plt.gcf()
            ## [left, bottom, width, height] # [left, bottom, width, height]
            # cax = fig.add_axes([0.2, 0.02, 0.6, 0.04])
            #plt.colorbar(bg, cax=cax, ticks=bg_dict['cmap_lvl'],
            #             orientation='horizontal', pad=0.05, shrink=0.4
            #).set_label(bg_dict['colbar_label'], fontsize=labelfont)
            #cax.tick_params(axis='x', which='major', labelsize=labelfont)
            plt.colorbar(bg, ticks=bg_dict['cmap_lvl'],
                         location='right',
#                         orientation='horizontal', pad=0.05, shrink=0.65
                         orientation='vertical', pad=0.05, shrink=0.4
            ).set_label(bg_dict['colbar_label'])

    # If there is no background, set a pale sea colour
    else:
        ax.set_facecolor(col_sea_trans)


    # Reading the trajectories
    try:
        trajdata = read_traj(trajinp)
    except:
        raise ValueError("{} is not able to be read as a trajectory file"
                         "".format(trajinp))

    # If a stride for end trajectories is set, then it is necessary to
    # filter the trajectories
    if stride != '1d':
        trajdata = filter_by_stride(trajdata, stride)

    # Checking if there is valid data
    if (trajdata['dimid'] == 0 or trajdata['dimreferencedate'] == 0
        or trajdata['dimtimestep'] == 0):
        raise ValueError("No valid trajectory data to plot, exiting.")

    # Setting up the colour map by end date, if this is selected
    if colmode == 'enddate':
        edrangei = []
        minreferencedatei = []
        maxreferencedatei = []
        for i in range(trajdata['dimid']):
            maxtedi = np.nanmax(trajdata['referencedate'][i, :])
            mintedi = np.nanmin(trajdata['referencedate'][i, :])
            maxreferencedatei.append(maxtedi)
            minreferencedatei.append(mintedi)
            edrangei.append((maxtedi - mintedi).days)
        edrange = max(edrangei)
        edcolmap = cm.viridis_r
        edcb = ax.imshow(np.array([[0.0, 0.0], [0.0, 0.0]]), cmap=edcolmap,
                         vmin=0, vmax=edrange)

    # Setting up the colour map by timestep, if this is selected
    if colmode == 'timestep':
        tsrange = trajdata['timestep'][-1] - trajdata['timestep'][0]
        tscolmap = cm.inferno
        tscb = ax.imshow(np.array([[0.0, 0.0], [0.0, 0.0]]), cmap=tscolmap,
                         vmin=0, vmax=tsrange)

    # Setting up the colour map by month if selected
    if colmode == 'month':
        mocolmap = colbar_mo_discrete()
        mocb = ax.imshow(np.array([[0.0, 0.0], [0.0, 0.0]]), cmap=mocolmap,
                         vmin=0, vmax=12)

    # Setting up the colour information for colour mode myi
    with open(mie_file) as mfile:
        miedict = json.load(mfile)
    mie_list = miedict["mie"][min_sie_ver][hemi]
    if hemi == 'nh':
        mie_list.update(miedict["mie_special"][min_sie_ver][hemi])

    # Setting a default plotting colour
    qcol = 'black'

    # Loop over each ID
    for i in range(trajdata['dimid']):
        if verbose:
            print("Plotting for id {}".format(trajdata['idname'][i]))

        # Setting the colour if per ID
        if colmode == 'id':

            if i < 20:
                qcol = cmap_list20[i]
            else:
                qcol = cmap_list20[i % 20]

        # Loop over each trajectory enddate
        # Don't use the dimension here in case this has been filtered to
        # be a shorter array
        # Reverse the array so the latest enddate is plotted last
        for j in reversed(range(len(trajdata['referencedate'][0, :]))):
            if verbose:
                print("Plotting for end date {}"
                      "".format(trajdata['referencedate'][i][j]))

            # Setting the colour if per enddate
            if colmode == 'enddate':
                if edrange > 0:
                    edindex = (trajdata['referencedate'][i][j]
                               - minreferencedatei[i]).days
                    qcol = edcolmap(edindex / edrange)
                else:
                    qcol = 'black'

            # Setting the colour if myi colouring required
            if colmode == 'myi':
                mie_date = find_mie_date(trajdata['referencedate'][i][j],
                                         mie_list)
                myindex = (trajdata['referencedate'][i][j] - mie_date).days

            # Plot the trajectory
            for k in range(trajdata['dimtimestep']):

                # Filter out masked and NaN lon/lats
                if (not np.ma.is_masked(trajdata['lon'][i, j, k])
                    and not np.ma.is_masked(trajdata['lat'][i, j, k])
                    and not np.isnan(trajdata['lon'][i, j, k])
                    and not np.isnan(trajdata['lat'][i, j, k])):

                    #print("i, j, k, lon, lat = ", i, j, k, trajdata['lon'][i, j, k], trajdata['lat'][i, j, k])

                    # Setting the colour if per timestep or per month
                    if colmode == 'timestep':
                        qcol = tscolmap(trajdata['timestep'][k] / tsrange)
                    if colmode == 'month':
                        mdate = trajdata['referencedate'][i][j] - timedelta(days=int(trajdata['timestep'][k]))
                        qcol = mocolmap(mdate.month - 1)

                    # Setting the colour if per myi
                    if colmode == 'myi':
                        if k <= myindex:
                            qcol = qmaroon
                        else:
                            qcol = qnavy

                    # Setting transparency if between the concentration
                    # thresholds
                    lalpha = 1
                    lstyle = 'solid'
                    if not ignore_mid_conc and trajdata['flag'][i, j, k] == 1:
                        lalpha = 1
                        lstyle = 'dotted'

                    # Just a dot is plotted for the first step
                    if k == 0:
                        x1, y1 = plot_crs.transform_point(
                            trajdata['lon'][i, j, k],
                            trajdata['lat'][i, j, k],
                            src_crs=pc)
#                        ax.plot(x1, y1, '-o', markersize=2, color=qcol,
#                                alpha=lalpha, linestyle=lstyle)
                        ax.plot(x1, y1, 'o', markersize=2, color=qcol,
                                alpha=lalpha, linestyle=lstyle)

#                        # Plot the ID if this option is selected and it's the
#                        # first enddate
#                        xlaboffset = 300000
#                        ylaboffset = 10000
#                        if label_traj and j == 0:
#                            label = '{}'.format(trajdata['idname'][i])
#                            if ((x1 > llx) and (x1 < urx) and (y1 > lly)
#                                and (y1 < ury)):
#                                ax.text(x1 - xlaboffset, y1 - ylaboffset,
#                                        label, fontsize=labelfont)
                        # Plot the ID if this option is selected and it's the
                        # first enddate
                        xlaboffset = 80000
                        ylaboffset = 50000
                        if label_traj and j == 0:
                            label = '{}'.format(str(i + 1))
                            ax.text(x1 - xlaboffset, y1 - ylaboffset,
                                    label, fontsize=labelfont)
                        if label_traj2 and j == 0:
                            label = '{}'.format(trajdata['idname'][i])
                            ax.text(x1 - xlaboffset, y1 - ylaboffset,
                                    label, fontsize=labelfont)

                    # For subsequent points, the lines are joined up
                    elif k > 0:
                        if (not np.isnan(trajdata['lon'][i, j, k]) and
                            not np.isnan(trajdata['lat'][i, j, k])):
                            x0 = x1
                            y0 = y1
                            x1, y1 = plot_crs.transform_point(
                                trajdata['lon'][i, j, k],
                                trajdata['lat'][i, j, k],
                                src_crs=pc)
#                            ax.plot([x0, x1], [y0, y1], '-', color=qcol,
#                                    linewidth=linew, alpha=lalpha,
#                                    linestyle=lstyle)
                            ax.plot([x0, x1], [y0, y1], color=qcol,
                                    linewidth=linew, alpha=lalpha,
                                    linestyle=lstyle)

    # Add a text box if labels
    if label_traj:
        leglist = []
        for i in range(trajdata['dimid']):
            leglist.append('{} - {}'.format(str(i + 1), trajdata['idname'][i]))
        textleg = '\n'.join(leglist)
        props = dict(boxstyle='round', facecolor='white', alpha=1)
        ax.text(0.05, 0.95, textleg, transform=ax.transAxes,
                fontsize=labelfont,
                verticalalignment='top', bbox=props)

    # Add colourbars
    props = dict(boxstyle='round', facecolor='white', alpha=1)
    if colmode == 'id':
        handles = []
        for i, idn in enumerate(trajdata['idname']):
            if i < 20:
                handles.append(mlines.Line2D([], [], color=cmap_list20[i],
                                             label=idn))
            else:
                handles.append(mlines.Line2D([], [], color=cmap_list20[i % 20],
                                             label=idn))
        ax.legend(handles=handles, fontsize=labelfont, loc=2,
                  facecolor='white', edgecolor='black', framealpha=1,
                  fancybox=True)
    elif colmode == 'myi':
        handles = []
        handles.append(mlines.Line2D([], [], color=qmaroon, label='MYI'))
        handles.append(mlines.Line2D([], [], color=qnavy,
                                     label='May not be MYI'))
        ax.legend(handles=handles, fontsize=labelfont, loc=2)
    elif colmode == 'enddate' and edrange > 0:
        if edrange <= 20:
            cmap_lvl = np.arange(0, edrange)
        elif edrange > 20 and edrange <= 200:
            cmap_lvl = np.arange(0, edrange, 10)
        elif edrange > 200 and edrange <= 2000:
            cmap_lvl = np.arange(0, edrange, 100)
        else:
            cmap_lvl = np.arange(0, edrange, 1000)
        plt.colorbar(edcb, ticks=cmap_lvl,
                     orientation='horizontal', pad=0.1, shrink=0.4
    ).set_label('End date of trajectory relative to the earliest date (days)')
    elif colmode == 'timestep':
        if tsrange <= 20:
            cmap_lvl = np.arange(0, tsrange)
        elif tsrange > 20 and tsrange <= 60:
            cmap_lvl = np.arange(0, tsrange, 10)
        elif tsrange > 60 and tsrange <= 300:
            cmap_lvl = np.arange(0, tsrange, 50)
        elif tsrange > 300 and tsrange <= 2000:
            cmap_lvl = np.arange(0, tsrange, 100)
        else:
            cmap_lvl = np.arange(0, tsrange, 1000)
        plt.colorbar(tscb, ticks=cmap_lvl,
                     orientation='horizontal', pad=0.1, shrink=0.4
        ).set_label('Elapsed time in trajectory (days)')
    elif colmode == 'month':
        tickpos = [x + 0.5 for x in range(12)]
        mcb = plt.colorbar(mocb, ticks=tickpos, orientation='horizontal',
                     pad=0.1, shrink=0.4).set_ticklabels(molabs)#.set_label('Month')


    # Bathymetry
    if bathymetry:
        bathym = cfeature.NaturalEarthFeature(name='bathymetry_J_1000',
        #    bathym = cfeature.NaturalEarthFeature(name='bathymetry_K_200',
                                              scale='10m', category='physical')
        bathym = cascaded_union(list(bathym.geometries()))
        ax.add_geometries(bathym, facecolor='none', edgecolor='darkgrey',
                          linestyle='solid', linewidth=0.5,
                          crs=ccrs.PlateCarree())
#    ax.add_feature(bathym, facecolor='none', edgecolor='darkgrey',
#                   linewidth=0.5)
#    ax.add_wms(
#        wms="http://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv?",
#        layers=["GEBCO_LATEST"],
#    )

    # Coastlines
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land',
                                                '50m', edgecolor='black',
                                                facecolor=col_land,
                                                linewidth=0.5))
    ice_shelves = cfeature.NaturalEarthFeature(
        category='physical',
        name='antarctic_ice_shelves_polys',
        scale='10m',
        facecolor=col_land,
        edgecolor=col_land)
    ax.add_feature(ice_shelves)

    # Gridlines
    if latlongrid:
        llats = [0.1 * x for x in range(3600)]
        for llon in [-80. + 5. * x for x in range(33)]:
            llons = [llon] * 3600
            xyes = [plot_crs.transform_point(lalo[0], lalo[1], src_crs=pc)
                    for lalo in zip(llats, llons)]
            xes = []
            yes = []
            for item in xyes:
                xes.append(item[0])
                yes.append(item[1])
            ax.plot(xes, yes, color='black', linewidth=0.2, linestyle='--')
        llats = [x * 20. for x in range(18)]
        for llat in llats:
            x0, y0 = plot_crs.transform_point(llat, 80., src_crs=pc)
            x1, y1 = plot_crs.transform_point(llat, -80., src_crs=pc)
            ax.plot([x0, x1], [y0, y1], color='black', linewidth=0.2,
                    linestyle='--')

    # Alternative gridlines with labels
    if latlongrid2:
        xticks = np.arange(-180, 181, 5)
        yticks = np.arange(-90, 91, 2)
        gl = ax.gridlines(crs=pc, draw_labels=True,
                          linewidth=1, color='darkgray', alpha=0.6, linestyle=':')
        gl.xlocator = mticker.FixedLocator(xticks)
        gl.ylocator = mticker.FixedLocator(yticks)

    # Title
    if title == '':
        title = trajdata['dataset_name']
        #title = '{}, {}'.format(trajdata['idname'][0],
        #                       datetime.strftime(trajdata['referencedate'][0][0],
        #                                         '%Y%m%d'))
    elif title == 'notitle':
        title = ''
    plt.title(title, fontsize=titlefont)

    # If no output is given, the input name is taken replacing .nc with .png,
    # if only a directory is given, then the name is taken replacing .nc with
    # .png and it is put into the defined directory, else if a file is given
    # then this is used
    if not os.path.isfile(output):
        if colmode in ['id', 'myi']:
            cmn = colmode
        elif colmode == 'timestep':
            cmn = 'ts'
        elif colmode == 'month':
            cmn = 'mo'
        elif colmode == 'enddate':
            cmn = 'ed'

        o = '{}_{}.png'.format(trajinp.replace('.nc', ''), cmn)
        if os.path.isdir(output):
            output = os.path.join(output, os.path.basename(o))
    plt.tight_layout()
    if plot_jpg:
        output = output.replace('.png', '.jpg')
        plt.savefig(output, bbox_inches='tight', format='jpeg', dpi=300)
    else:
        plt.savefig(output, bbox_inches='tight')
    plt.show()
    plt.close()
    print("Figure is in {}".format(output))


def main():

    args = parse_args()
    trajinp = args.trajinp
    region = args.region
    output = args.output
    label_traj = args.label_traj
    label_traj2 = args.label_traj2
    colmode = args.colmode
    title = args.title
    stride = args.stride
    bathymetry = args.bathymetry
    latlongrid = args.latlongrid
    latlongrid2 = args.latlongrid2
    bg_file = args.bg_file
    bg_var = args.bg_var
    bg_ind = args.bg_ind
    if bg_ind is not None:
        bg_ind = int(bg_ind)
    bg_date = args.bg_date
    linew = args.linew
    ignore_mid_conc = args.ignore_mid_conc
    plot_jpg = args.plot_jpg
    verbose = args.verbose

    plot_traj(trajinp, region=region, output=output, label_traj=label_traj,
              label_traj2=label_traj2, colmode=colmode, title=title,
              stride=stride, bathymetry=bathymetry, latlongrid=latlongrid,
              latlongrid2 = args.latlongrid2, bg_file=bg_file, bg_var=bg_var,
              bg_ind=bg_ind, bg_date=bg_date, linew=linew,
              ignore_mid_conc=ignore_mid_conc, plot_jpg=plot_jpg,
              verbose=verbose)


if __name__ == '__main__':

    main()
