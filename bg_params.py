import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import numpy as np
import cmocean


def colbar_discrete_greys():

   # Discrete colorbar
    cmap = cm.Greys
    cmaplist = [cmap(i) for i in range(cmap.N)]
    collen = len(cmaplist)
    col0 = cmaplist[0]
    col1 = cmaplist[int(collen / 6) - 1]
    col2 = cmaplist[int(collen / 3) - 1]
    # OK, a hack to make values up to 2.1 white... (knowing that the
    # full range is 1)
    whitenum = 0
    for i in range(collen):
        if i <= whitenum:
            cmaplist[i] = col0 # Leave colourmap at white
        elif (i > whitenum) and (i < (collen - 1)):
            cmaplist[i] = col1
        else:
            cmaplist[i] = col2
    # Create the new map
    discrete_greys = matplotlib.colors.LinearSegmentedColormap.from_list(
        'discrete_greys', cmaplist, cmap.N)
    return discrete_greys


def colbar_sic_discrete():

    # SIC discrete colorbar with bottom 5% as white
#    cmap_conc = cm.get_cmap('gray', 24)
    cmap_conc = cm.get_cmap(cmocean.cm.ice, 24)
    colors_i = np.linspace(0.0, 1.0, 24)
    c24 = cmap_conc(colors_i)
    c20 = c24[4:]
    # Below 5% conc set to white
    c20[0] = colors.colorConverter.to_rgba_array('white')
    # This is only needed for 10% conc
    cmap_conc20 = colors.LinearSegmentedColormap.from_list('listconc20',
                                                           c20, 20)
    print("colbar = ", cmap_conc20)
    print("type(colbar) = ", type(cmap_conc20))
    return cmap_conc20


def colbar_simask_discrete():

    # Sea ice mask discrete colorbar with close ice as white
    cmap_conc = cm.get_cmap('PuBu_r', 3)
    colors_i = np.linspace(0.0, 1.0, 3)
    c3 = cmap_conc(colors_i)
    # Set high value set to white
    c3[-1] = colors.colorConverter.to_rgba_array('white')
    cmap_mask3 = colors.LinearSegmentedColormap.from_list('listmask3',
                                                           c3, 3)
    return cmap_mask3


def colbar_simask_discrete2():
    # Sea ice mask discrete colorbar with close ice as white, matching
    # new ice conc quicklooks

    cmap_sim = colors.LinearSegmentedColormap.from_list('listsim',
                                [col_sea, col_iceconc50, col_iceconc100], 3)

    return cmap_sim


def colbar_sit_discrete():

    # SIT discrete colorbar
#    cmap_sit = cm.get_cmap('Blues_r', 24)
    cmap_sit = cm.get_cmap(cmocean.cm.ice, 24)
    colors_i = np.linspace(0.0, 1.0, 24)
    csit28 = cmap_sit(colors_i)
    csit20 = csit28[8:]
    # Below 5% sit set to white
    cmap_sit20 = colors.LinearSegmentedColormap.from_list('listsit20',
                                                           csit20, 20)
    return cmap_sit20


def colbar_sity_discrete():

    # Sea-ice type discrete colorbar
    cmap_sity= cm.get_cmap(cmocean.cm.ice_r, 6)
    colors_i = np.linspace(0.0, 1.0, 6)
#    # Reverse ambig and MYI colours
#    am = colors_i[-1]
#    my = colors_i[-2]
#    colors_i[-1] = my
#    colors_i[-2] = am
    csityl = cmap_sity(colors_i)
    csityr = csityl[0:4]
    csity = csityr[::-1]
    cmap_sity = colors.LinearSegmentedColormap.from_list('listsity4',
                                                           csity, 4)
    return cmap_sity

def colbar_sity_like_age():

    # Sea ice age discrete colorbar
    cmap_sia = cm.get_cmap(cmocean.cm.ice, 10)
    colors_i = np.linspace(0.0, 1.0, 10)
    csia_l = cmap_sia(colors_i)
    csia = csia_l[4:]

    cmap_sity= cm.get_cmap(cmocean.cm.ice_r, 4)
    colors_i = np.linspace(0.0, 1.0, 4)
    csity = cmap_sity(colors_i)
    csity[0] = colors.colorConverter.to_rgba_array('white')
    csity[1] = csia[1]
    csity[2] = csia[2]
    csity[3] = csia[3]
    cmap_sity = colors.LinearSegmentedColormap.from_list('listsity4',
                                                         csity, 4)
    return cmap_sity



def colbar_ffsf_discrete(numcol):

    # Final format status flags discrete colorbar
    cmap_ffsf= cm.get_cmap(cm.Spectral, numcol)
    colors_i = np.linspace(0.0, 1.0, numcol)
    cffsf = cmap_ffsf(colors_i)
    cmap_ffsf = colors.LinearSegmentedColormap.from_list('listffsf',
                                                           cffsf, numcol)
    return cmap_ffsf


def colbar_type_status_flags():

    # Sea-ice type flags discrete colorbar
    cmap_sitysf = cm.get_cmap(cmocean.cm.deep, 4)
    colors_i = np.linspace(0.0, 1.0, 4)
    csitysf = cmap_sitysf(colors_i)
    cmap_sitysf = colors.LinearSegmentedColormap.from_list('listsitysf4',
                                                           csitysf, 4)
#    bounds = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
#              1000]
#    bounds = [x - 0.5 for x in bounds]
#    norm = colors.BoundaryNorm(bounds, cmap_sitysf.N)
    return cmap_sitysf#, bounds, norm


#def colbar_ice_age():
#
#    # Sea ice age discrete colorbar
#    cmap_sia = cm.get_cmap(cmocean.cm.ice, 15)
#    colors_i = np.linspace(0.0, 1.0, 15)
#    csia_l = cmap_sia(colors_i)
#    csia = csia_l[4:]
#    # Set last segment to white
#    csia[-1] = colors.colorConverter.to_rgba_array('white')
#    cmap_sia = colors.LinearSegmentedColormap.from_list('listsia',
#                                                        csia, 11)
#    return cmap_sia


# UP TO 10+ YEARS
#def colbar_ice_age():
#
#    # Sea ice age discrete colorbar
#    cmap_sia = cm.get_cmap(cmocean.cm.ice, 15)
#    colors_i = np.linspace(0.0, 1.0, 15)
#    csia_l = cmap_sia(colors_i)
#    csia = csia_l[4:]
#    # Set last segment to white
#    csia[0] = colors.colorConverter.to_rgba_array('white')
#    cmap_sia = colors.LinearSegmentedColormap.from_list('listsia',
#                                                        csia, 11)
#    return cmap_sia

def colbar_ice_age():

    # Sea ice age discrete colorbar
    cmap_sia = cm.get_cmap(cmocean.cm.ice, 10)
    colors_i = np.linspace(0.0, 1.0, 10)
    csia_l = cmap_sia(colors_i)
    csia = csia_l[4:]
    # Set last segment to white
    csia[0] = colors.colorConverter.to_rgba_array('white')
    cmap_sia = colors.LinearSegmentedColormap.from_list('listsia',
                                                        csia, 6)
    return cmap_sia


def bg_plot_setup(bgvar=None, bgnc=None):
    '''Set up background quantities for plotting'''

    bg = {}
    bg['alpha'] = 1
    bg['norm'] = None
    bg['bounds'] = None
    bg['colbar_type'] = 'cont'

    if bgvar in ['coastmask_250']:
        bg['lbl']  = 'Landmask'
        bg['min'] = 0
        bg['max'] = 6
        bg['cmap_lvl']  = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
        bg['cmap_labs']  = ['0', '1', '2', '3', '4', '5']
        bg['cmap'] = colbar_ffsf_discrete(6)
        bg['colbar_label'] = ""
        bg['colbar_type'] = 'discrete'

#    # Sea-ice age up to 10+ years
#    if bgvar in ['age_of_sea_ice']:
#        bg['lbl']  = 'Age of Sea Ice'
#        bg['min'] = -0.5
#        bg['max'] = 10.5
#        bg['cmap_lvl']  = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#        bg['cmap_labs']  = ['', '1', '2', '3', '4', '5', '6', '7', '8', '9', '1#0+']
#        bg['cmap'] = colbar_ice_age()
#        bg['colbar_label'] = "Age of Sea Ice (years)"
#        bg['colbar_type'] = 'discrete'

    # Sea-ice age up to 5+ years
    if bgvar in ['age_of_sea_ice']:
        bg['lbl']  = 'Age of Sea Ice'
        bg['min'] = -0.5
        bg['max'] = 5.5
        bg['cmap_lvl']  = [0, 1, 2, 3, 4, 5]
        bg['cmap_labs']  = ['', '1', '2', '3', '4', '5+']
        bg['cmap'] = colbar_ice_age()
        bg['colbar_label'] = "Age of Sea Ice (years)"
        bg['colbar_type'] = 'discrete'

        
    # Average drift in background
    if bgvar in ['avdX', 'avdY']:
        # TODO - Make drift period flexible
        bg['lbl']  = 'Average {} component over 24 hours [km]'.format(bgvar)
        bg['lim'] = 30
        bg['max'] = bg['lim']
        bg['min'] = -bg['lim']
        bg['cmap_lvl']  = [bg['min'], 0, bg['max']]
        bg['cmap'] = cm.RdBu
        bg['colbar_label'] = "metres"

    # TODO - WRITE A CODE WHICH CREATES A TEMP NC FILE WITH DISPLACEMENT
    # AS A VARIABLE
    elif bgvar == 'displacement':
       bg['lbl']  = 'Displacement over 24 hours [km]'
       #bg['fld']  = pow((pow(u,2) + pow(v,2)), 0.5) / 1000
       bg['max'] = 30
       bg['min']  = 0
       bg['cmap'] = cmocean.cm.thermal
       bg['cmap_lvl'] = np.arange(bg['min'], bg['max'], 5)
       bg['cmap_fmt'] = '%2.0f'
       bg['colbar_label'] = bg['lbl']

#    # Original status flag array
#    elif bgvar in ['statusflag', 'status_flag', 'flag']:
#        bg['lbl']  = 'Status flags'
#        bg['max'] = np.nanmax(bgnc[bgvar])
#        bg['min'] = 0
##        bg['cmap'] = cm.viridis
#        bg['cmap'] = cm.Spectral
#        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'])
#        bg['colbar_label'] = ""
#        bg['colbar_type'] = 'discrete'

    # Final format status flag array
    elif bgvar in ['statusflag', 'status_flag', 'flag']:
        bg['lbl']  = 'Status flags'
        bg['max'] = np.nanmax(bgnc[bgvar]) + 1
        bg['min'] = 0
#        bg['cmap'] = cm.viridis
#        bg['cmap'] = cm.Spectral
        bg['cmap'] = colbar_ffsf_discrete(int(bg['max']))
#        bg['cmap_lvl'] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 11.5, 12.5,
#                          13.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5, 30.5]
#        bg['cmap_labs'] = ['missing_input_data', 'over_land', 'no_ice',
#                           'close_to_coast_or_edge', 'summer_period',
#                           'masked_outside_icedrift', 'processing_failed',
#                           'too_low_correlation', 'not_enough_neighbours',
#                           'filtered_by_neighbours', 'smaller_pattern',
#                           'corrected_by_neighbours', 'interpolated',
#                           'gapfilled_wind_parameter', 'filled_by_wind',
#                           'blended_with_wind', 'nominal_quality']
        #bg['cbar_ticks'] = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
        #                    9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5]
#        bg['cmap_labs'] = ['0', '1', '2', '3', '4', '5', '10', '11', '12',
#                           '13', '20', '21', '22', '23', '24', '25', '30']
        bg['cmap_labs'] = bgnc['sf_labs']
        bg['cmap_lvl'] = np.arange(bg['min'], bg['max'])
        bg['cmap_lvl'] = [x + 0.5 for x in bg['cmap_lvl']]
        bg['colbar_label'] = ""
        bg['colbar_type'] = 'discrete'

#    # Status flags for ice type
#    elif bgvar in ['statusflag', 'status_flag']:
#        bg['lbl']  = 'Status flags'
#        bg['max'] = 4
#        bg['min'] = 0
#        bg['cmap'] = colbar_type_status_flags()
##        bg['cmap'], bg['bounds'], bg['norm'] = colbar_type_status_flags()
##        bg['cmap_lvl'] = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
#        bg['cmap_lvl'] = [0.5, 1.5, 2.5, 3.5]
#        bg['cmap_labs'] = ['Land', 'Climatology', 'Drift Mask', 'Unclassified coast']
##        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'], 0.1)
#        bg['colbar_type'] = 'discrete'
#        bg['colbar_label'] = "Status flags"

    elif bgvar in ['uncertainty']:
        bg['lbl']  = 'Uncertainty'
        bg['max'] = 1
        bg['min'] = 0
        bg['cmap'] = cmocean.cm.dense
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'], 0.1)
        bg['colbar_label'] = "Uncertainty"

    elif bgvar in ['uncert_dX_and_dY']:
        bg['lbl']  = 'Uncertainty'
        bg['max'] = ceil(np.nanmax(bgnc[bgvar]))
        bg['min'] = 0
        bg['cmap'] = cmocean.cm.dense
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'], 1)
        bg['colbar_label'] = "Uncertainty"

#    elif bgvar == 'avavdiv':
#        avavdiv = dataset.variables['avavdiv'][:, 0:-1, 0:-1]
#        avavdiv_arr = ma.array(avavdiv[0,:,:], fill_value=fv)
#        deform_pow = 6
#        bg['lbl'] = ('Divergence from monthly average drift '
#                          '[10^{} x 1/s]'.format(deform_pow))
#        bg['fld'] = avavdiv_arr * (10.**deform_pow)
#        bg['lim'] = 2
#        bg['max'] = bg['lim']
#        bg['min'] = -bg['lim']
#        bg['cmap'] = cm.seismic_r
#        bg['cmap_lvl']  = [bg['min'], 0, bg['max]
#        bg['colbar_label'] = ""

#    elif bgvar == 'avdivs':
#        avdivs = dataset.variables['avdivs'][:, 0:-1, 0:-1]
#        avdivs_arr = ma.array(avdivs[0,:,:], fill_value=fv)
#        deform_pow = 6
#        bg['lbl'] = ('Average divergence from individual drifts '
#                          '[10^{} x 1/s]'.format(deform_pow))
#        bg['fld'] = avdivs_arr * (10.**deform_pow)
#        bg['lim'] = 3
#        bg['max'] = bg['lim']
#        bg['min'] = 0
#        bg['cmap_lvl'] = [bg['min'], bg['max']]
#        bg['cmap'] = cm.Blues
#        bg['colbar_label'] = ""

#    elif bgvar == 'avconvs':
#        avconvs = dataset.variables['avconvs'][:, 0:-1, 0:-1]
#        avconvs_arr = ma.array(avconvs[0,:,:], fill_value=fv)
#        deform_pow = 6
#        bg['lbl'] = ('Average convergence from individual drifts '
#                          '[10^{} x 1/s]'.format(deform_pow))
#        bg['fld'] = avconvs_arr * (10.**deform_pow)
#        bg['lim'] = -3
#        bg['max'] = 0
#        bg['min'] = bg['lim']
#        bg['cmap_lvl']  = [bg['min'], bg['max']]
#        bg['cmap'] = cm.Reds_r
#        bg['colbar_label'] = ""

#    elif bgvar == 'divsconvs':
#        avconvs = dataset.variables['avconvs'][:, 0:-1, 0:-1]
#        avconvs_arr = ma.array(avconvs[0,:,:], fill_value=fv)
#        avdivs = dataset.variables['avdivs'][:, 0:-1, 0:-1]
#        avdivs_arr = ma.array(avdivs[0,:,:], fill_value=fv)
#        deform_pow = 6
#        bg['lbl'] = ('Average convergence/divergence from '
#                          'individual drifts [10^{} x 1/s]'.format(deform_pow))
#        bg['divs'] = avdivs_arr * (10.**deform_pow)
#        bg['convs'] = avconvs_arr * (10.**deform_pow)
#        bg['lim'] = 3
#        bg['max'] = bg['lim']
#        bg['min'] = -bg['lim']
#        bg['cmap_lvl']  = [bg['min'], 0, bg['max']]
#        bg['cmap'] = cm.seismic_r
#        bg['colbar_label'] = "metres"

#    elif bgvar == 'avdivcube':
#        avdivcube = dataset.variables['avdivcube'][:, 0:-1, 0:-1]
#        avdivcube_arr = ma.array(avdivcube[0,:,:], fill_value=fv)
#        deform_pow = 6
#        bg['lbl'] = ('Average net convergence/divergence from '
#                          'individual drifts [10^{} x 1/s]'.format(deform_pow))
#        bg['fld'] = avdivcube_arr * (10.**deform_pow)
#        bg['lim'] = 3
#        bg['max'] = bg['lim
#        bg['min'] = -bg['lim
#        bg['cmap_lvl']  = [bg['min'], bg['max']]
#        bg['cmap'] = cm.seismic_r
#        bg['colbar_label'] = "metres"

#    elif bgvar in ['middleicetype']:
#        bg['lbl'] = 'Sea-ice type'
#        bg['min'] = 2
#        bg['max'] = 3
#        bg['cmap_lvl'] = [2, 3]
#        bg['cmap'] = cmocean.cm.thermal#discrete_greys
#        bg['colbar_label'] = "Sea-ice type"

    elif bgvar in ['icetype', 'ice_type']:
        bg['lbl'] = 'Sea-ice type'
        bg['min'] = 0.5
        bg['max'] = 4.5
        bg['colbar_type'] = 'discrete'
        bg['cmap_lvl'] = [1, 2, 3, 4]
        bg['cmap_labs'] = ['Open Water', 'FYI', 'FYI/MYI', 'MYI']
#        bg['cmap'] = colbar_sity_discrete()
        bg['cmap'] = colbar_sity_like_age()
#        bg['cmap'] = cmocean.cm.ice_r#discrete_greys
        bg['colbar_label'] = "Sea-ice type"

    elif bgvar in ['ice_edge']:
        bg['lbl'] = 'Sea-ice edge'
        bg['min'] = 0.5
        bg['max'] = 3.5
        bg['cmap'] = colbar_simask_discrete()
        bg['cmap_lvl'] = [1, 2, 3]
        bg['cmap_labs'] = ['Open Water', 'Open Ice', 'Close Ice']
        bg['colbar_label'] = "Sea-ice edge"
        bg['colbar_type'] = 'discrete'

#    elif bgvar == 'avice':
#        avice = dataset2.variables['avicetype'][:, :, :]
#        avice_arr = ma.array(avice[0,:,:], fill_value=fv)
#        bg['lbl'] = dataset2.variables['avicetype'].__dict__['long_name']
#        bg['fld'] = avice_arr
#        bg['max'] = 3
#        bg['min'] = 2
#        bg['cmap_lvl'] = [2, 3]
#        bg['cmap'] = discrete_greys
#        bg['colbar_label'] = ""

#    elif bgvar == 'Vorticity':
#        deform_pow = 6
#        bg['lbl'] = 'Vorticity [10^{} x 1/s]'.format(deform_pow)
#        bg['fld'] = deform_curl * (10.**deform_pow)
#        bg['max'] = +3
#        if bg['lim'] is not None:
#            bg['max'] = float(bg['lim'])
#        bg['min'] = -bg['max']
#        bg['cmap'] = cm.PiYG_r
#        bg['cmap_lvl']  = [-3,-1.5,0,1.5,3]

    elif bgvar == 'sit':
        bg['lbl'] = 'Sea-ice thickness [m]'
        bg['min'] = 0
        bg['max'] = 2.5
        bg['cmap'] = colbar_sit_discrete()
        bg['cmap_lvl'] = np.arange(bg['min'], bg['max'], 0.5)
        bg['cmap_fmt'] = '%2.0f'
        bg['colbar_label'] = "Sea-ice thickness (m)"

    elif bgvar == 'sit_anom':
        bg['lbl'] = 'Sea-ice thickness anomaly [m]'
        bg['min'] = -1.2
        bg['max'] = 1.2
        bg['cmap'] = cmocean.cm.balance_r
#        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'],0.5)
        bg['cmap_lvl'] = [-1, -0.5, 0, 0.5, 1]
        bg['cmap_fmt'] = '%2.0f'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                         vcenter=0,
                                                         vmax=bg['max'])
        bg['colbar_label'] = "Sea-ice thickness anomaly (m)"

#    elif bgvar == 'conc':
    elif bgvar == 'ice_conc':
#        conc = dataset2.variables['av_ice_conc'][:, :, :]
#        conc_arr = ma.array(conc[0,:,:], fill_value=fv)
#        conc_arr[conc_arr < 5.] = fv
#        conc_arr.mask[conc_arr < 5.] = True
#        bg['fld'] = conc_arr
        bg['lbl'] = 'Sea-ice concentration [%]'
        bg['min'] = 0
        bg['max'] = 100
        bg['cmap'] = colbar_sic_discrete()
#        bg['cmap'] = cmocean.cm.gray
#        bg['cmap'] = cm.PuBu_r
#        bg['cmap'] = cm.Blues_r
        bg['cmap_lvl']  = np.arange(bg['min'],bg['max'],10)
        bg['cmap_fmt']  = '%2i'
        bg['colbar_label'] = "Sea-ice concentration (%)"

    elif bgvar == 'ice_mask':
        bg['lbl'] = 'Sea-ice concentration'
        bg['min'] = 0.5
        bg['max'] = 3.5
        bg['cmap'] = colbar_simask_discrete2()
        bg['cmap_lvl'] = [1, 2, 3]
        bg['cmap_labs'] = ['Open Water', 'Open Ice', 'Close Ice']
        bg['colbar_label'] = "Sea-ice concentration"
        bg['colbar_type'] = 'discrete'

    elif bgvar == 'conc_anom':
        bg['lbl'] = 'Sea-ice concentration anomaly [%]'
        bg['min'] = -100
        bg['max'] = 100
        bg['cmap'] = cmocean.cm.balance_r
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'],20)
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Sea-ice concentration anomaly (%)"

    elif bgvar == 'sst':
        bg['lbl'] = 'Sea surface temperature [K]'
        bg['min'] = 270
        bg['max'] = 290
        bg['cmap'] = cmocean.cm.speed
        bg['cmap_lvl'] = np.arange(bg['min'],bg['max'],5)
        bg['cmap_fmt'] = '%2i'
        bg['colbar_label'] = "Sea surface temperature (K)"

    elif bgvar == 'sst_anom':
        bg['lbl'] = 'Sea surface temperature anomaly [K]'
        bg['min'] = -6
        bg['max'] = 6
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-4, -2, 0, 2, 4]
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Sea surface temperature anomaly (K)"

    elif bgvar in ['uwind_avg', 'vwind_avg']:
        bg['lbl'] = 'Wind component {}'.format(bgvar)
        bg['min'] = -12
        bg['max'] = 12
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-12, -6, 0, 6, 12]
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Wind component {}".format(bgvar)

    elif bgvar in ['driftX', 'driftY', 'dX', 'dY']:
        bg['lbl'] = '{} component over 24 hours [km]'.format(bgvar)
        bg['min'] = -30
        bg['max'] = 30
        #bg['cmap'] = cmocean.cm.balance_r
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [-30, -15, 0, 15, 30]
        bg['cmap_fmt'] = '%2i'
        bg['norm'] = matplotlib.colors.TwoSlopeNorm(vmin=bg['min'],
                                                    vcenter=0,
                                                    vmax=bg['max'])
        bg['colbar_label'] = "Ice drift component {}".format(bgvar)

    elif bgvar in ['t0', 't1']:
        if bgvar == 't0':
            bg['lbl'] = 'Start time [hours UTC]'
        elif bgvar == 't1':
            bg['lbl'] = 'Stop time [hours UTC]'
#        bg['min'] = -24
#        bg['max'] = 24
        bg['min'] = -12
        bg['max'] = 12
        bg['cmap'] = cmocean.cm.balance
#        bg['cmap_lvl'] = [-24, -18, -12, -6, 0, 6, 12, 18, 24]
        bg['cmap_lvl'] = [-12, -9, -6, -3, 0, 3, 6, 9, 12]
        bg['cmap_fmt'] = '%2i'
        bg['colbar_label'] = bg['lbl']

    elif bgvar in ['tdiff']:
        bg['lbl'] = 'Time span [hours]'
        bg['min'] = 21
        bg['max'] = 27
        bg['cmap'] = cmocean.cm.balance
        bg['cmap_lvl'] = [21, 22, 23, 24, 25, 26, 27]
        bg['cmap_fmt'] = '%2i'
        bg['colbar_label'] = bg['lbl']

    return bg
