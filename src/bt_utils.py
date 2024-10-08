'''
More general useful functions for use with the tracking code.
'''
import os
import sys
import re
from math import fabs
from datetime import datetime
import numpy as np
import pyresample as pr
from netCDF4 import Dataset


def load_grid_defs_from_OSISAF_ncCF_file(ncCF_file, verbose=False):
    '''
    Loads grid definitions from an OSI SAF netcdf file

    Searches for a variable with attribute 'grid_mapping'. Use the
    value of that attribute to load proj string information, and then
    use xc and yc for defining the extent of the pyresample area_def object.
    '''

    def get_area_id(area, projection, resolution):
        reso = int(float(resolution)*10)
        area_id = '{}{}'.format(area, reso)
        pcs_id = '{}_{}-{}'.format(area, projection, reso)
        return area_id, pcs_id

    with Dataset(ncCF_file, 'r') as gf:
        # proj_dict
        crs_name = None
        for varn in list(gf.variables.keys()):
            # First try to get the CRS name from the global variables
            try:
                a_var = gf.variables[varn]
                crs_name = a_var.grid_mapping
                break
            except AttributeError:
                pass
        if crs_name is None:
            raise ValueError("ERROR: did not find any variable with "
                             "'grid_mapping' attribute")
        else:
            if verbose:
                print("Info: read grid_mapping information from {} (found in "
                      "{})".format(crs_name, varn))

        try:
            proj_str = gf.variables[crs_name].proj4_string
        except:
            proj_str = gf.variables[crs_name].proj4text
        proj_dict = dict([s[1:].split('=') for s in proj_str.split()
                          if '=' in s])

        # xcs
        try:
            xcs = gf.variables['xc']
        except:
            xcs = gf.variables['x']
        if xcs.units == 'km':
            xcs = xcs[:] * 1000.
        elif xcs.units in ['m', 'meters']:
            xcs = xcs[:]
        else:
            raise NotImplementedError("Do not know how to handle 'xc' with "
                                      "units {}".format(xcs.units))
        # ycs
        try:
            ycs = gf.variables['yc']
        except:
            ycs = gf.variables['y']
        if ycs.units == 'km':
            ycs = ycs[:] * 1000.
        elif ycs.units in ['m', 'meters']:
            ycs = ycs[:]
        else:
            raise NotImplementedError("Do not know how to handle 'yc' with "
                                      "units {}".format(ycs.units))

        # Shape
        shape = a_var.shape[::-1]
        if len(shape) == 3:
            shape = (shape[0], shape[1])

        # Grid spacing
        xgs = fabs(xcs[1] - xcs[0])
        ygs = fabs(ycs[1] - ycs[0])

        # Extent
        extent = [xcs.min(), ycs.min(), xcs.max(), ycs.max()]
        extent[0] = extent[0] - (0.5 * xgs)
        extent[1] = extent[1] - (0.5 * ygs)
        extent[2] = extent[2] + (0.5 * xgs)
        extent[3] = extent[3] + (0.5 * ygs)

        # Try and guess an area_id from global attributes
        area_id, pcs_id = crs_name, crs_name
        try:
            area_id, pcs_id = get_area_id(gf.area, gf.projection,
                                          gf.resolution.split(' ')[0])
        #except AttributeError as ae:
        #    print("Missing global attribute {} to create an area ID".format(ae))
        except Exception as ex:
            if verbose:
                print("Got {}".format(ex))

            area_def = pr.geometry.AreaDefinition(area_id, crs_name, pcs_id,
                                                  proj_dict,
                                                  shape[0], shape[1],
                                                  extent)

        return area_def


def print_percentage(step, total_steps):
    '''
    Print the completion percentage.
    '''
    percentage = int((step / total_steps) * 100)
    sys.stdout.write('\r')
    #sys.stdout.write(f'{percentage:.2f}%')
    sys.stdout.write(f'Tracking... {percentage}%')
    sys.stdout.flush()


def strip_date(listdirs):
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
