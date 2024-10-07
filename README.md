The file req_bt_locs.yml can be used to create a conda environment which has the necessary python dependencies to run the code.

The main code is backtrack_loc.py which traces trajectories of sea ice parcels back in time to where the sea ice is found to be below 
a certain concentration, in order to find a probable origin location for this sea ice. This code uses sea ice drift and concentration 
files from the OSI SAF project which it accesses via a thredds server. The output is a NetCDF file with the sea ice position for each date.

This code can be used in two modes:
* Ice station mode. In this mode, the trajectory of one or more packets of sea ice are traced back in time to either the time limit set,
  or to when the sea-ice concentration falls below a threshold. The inputs to this are a location and date for each packet of ice at the
  end of the trajectory.
* Moorings mode. In this mode, one or more stationary locations are set, together with start and end dates, and the tracking is performed
  for the location for each date in the date range between the start and end dates. This mode is to allow information to be gathered about
  how sea-ice origin for sea ice drifting past a fixed location might change over a time period.

In addition it is possible to track the sea ice forward by a fixed time period.

As a supplement to this main tracking code, there is a plot_trajectories.py code which can be used to prepare PDF plots of the trajectories
from the NetCDF files.

This code was prepared at MET Norway with funding from the SUDARCO project. 
It is licensed under the GNU Affero General Public License v3.0.
Please credit use with "<Data/plot> prepared with use of code funded by the SUDARCO project and prepared at Met Norway".
