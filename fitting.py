# Import required packages.
from timeit import default_timer
from wisp_analysis import *
from mpfit import *

'''The steps for adding additional emission lines to the model are the following:
1.) Define the vacuum rest wavelength for your emission line in the list below.
2.) Define the index number that will represent the model amplitude for that line.
3.) Set the value of 'first_node_index' to one greater than the last line index.
4.) Add the line to the three locations where redshifted centroids are calculated.
NOTE: Ensure that each model parameter index coded below is only used once.'''

# Define the emission line vacuum wavelengths.
la_1216_vac = 1215.670
n5_1238_vac = 1238.821
n5_1242_vac = 1242.804
c4_1548_vac = 1548.203
c4_1550_vac = 1550.777
h2_1640_vac = 1640.420
o3_1660_vac = 1660.8092
o3_1666_vac = 1666.1497
s3_1883_vac = 1882.707
s3_1892_vac = 1892.030
c3_1907_vac = 1906.680
c3_1909_vac = 1908.730
m2_2796_vac = 2796.352
m2_2803_vac = 2803.531
o2_3727_vac = 3727.092
o2_3730_vac = 3729.875
hg_4342_vac = 4341.684
o3_4363_vac = 4364.436
h2_4686_vac = 4687.020
hb_4863_vac = 4862.683
o3_4959_vac = 4960.295
o3_5007_vac = 5008.240
o1_6300_vac = 6302.046
o1_6363_vac = 6365.536
n2_6550_vac = 6549.850
ha_6565_vac = 6564.610
n2_6585_vac = 6585.280
s2_6716_vac = 6718.290
s2_6731_vac = 6732.670
s3_9069_vac = 9071.100
s3_9532_vac = 9533.200
he10830_vac = 10832.86

# Define the index number for each model parameter that will be fit by MPFIT.
# This allows us to avoid hard-coding indices everywhere else within the code.
# Strongly suggest not changing the first five parameters and keeping dz's grouped.
# Otherwise, the code must be modified to initialize values and fix params for mpfit.
nnodes_idx  = 0
fitreg_idx  = 1
t_wave_idx  = 2
z_idx       = 3
fwhm_grism_idx = 4 # fwhm_red
fwhm_ratio_idx = 5 # fwhm_blue = (ratio_fwhm * fwhm_red)

la_1216_dzx = 6  # dz
c4_1548_dzx = 7  # dz
uv_line_dzx = 8  # dz
m2_2796_dzx = 9  # dz
o2_3727_dzx = 10 # dz
o3_5007_dzx = 11 # dz
s3_he_dzx   = 12 # dz

la_1216_idx = 13
n5_1238_idx = 14
n5_1242_idx = 15 # ratio
c4_1548_idx = 16
c4_1550_idx = 17 # ratio
h2_1640_idx = 18
o3_1660_idx = 19
o3_1666_idx = 20 # ratio
s3_1883_idx = 21
s3_1892_idx = 22 # ratio
c3_1907_idx = 23
c3_1909_idx = 24 # ratio
m2_2796_idx = 25
m2_2803_idx = 26 # ratio
o2_3727_idx = 27
o2_3730_idx = 28 # ratio
hg_4342_idx = 29
o3_4363_idx = 30
h2_4686_idx = 31
hb_4863_idx = 32
o3_4959_idx = 33
o3_5007_idx = 34 # ratio
o1_6300_idx = 35
o1_6363_idx = 36 # ratio
n2_6550_idx = 37
ha_6565_idx = 38
n2_6585_idx = 39 # ratio
s2_6716_idx = 40
s2_6731_idx = 41 # ratio
s3_9069_idx = 42
s3_9532_idx = 43 # ratio
he10830_idx = 44
# lyman alpha wing.
la_wing_amp_idx = 45
la_wing_sig_idx = 46 # ratio

################################################################################

first_line_index = 6  # This should be the first emission line parameter index.
first_node_index = 47 # This must be one larger than last line parameter index.
ratio_indices = [15, 17, 20, 22, 24, 26, 28, 34, 36, 39, 41, 43, 46] # A list of all parameter indices that are ratios below.

def get_ratio_indices():
    return ratio_indices

def get_fitpar_indices():
    return first_line_index, first_node_index

################################################################################

def emissionline_model(pars, x):
    # The emission line model includes 'pars' for the continuum and emission line
    # parameters. The input 'x' is the wavelength array, which can include gaps.
    # mpfit.py allows for 'pegged', 'fixed', and 'limited' parameters with 'limits'.
    # See the mpfit.py code for a detailed description, but it is useful to note
    # that 'pegged' parameters are simply those 'fixed' at their upper or lower limit.

    # fixed parameters
    nnodes          = int(pars[nnodes_idx]) # calculated from len(node_wave) and passed here.
    fit_region      = pars[fitreg_idx]
    transition_wave = pars[t_wave_idx]
    z               = pars[z_idx]
    ratio_fwhm      = pars[fwhm_ratio_idx]
    fwhm_red        = pars[fwhm_grism_idx]
    fwhm_blue       = (ratio_fwhm * fwhm_red)
    sigma_blue      = (fwhm_blue / 2.3548)
    sigma_red       = (fwhm_red  / 2.3548)

    # redshift wiggle room
    la_1216_dz      = pars[la_1216_dzx]
    c4_1548_dz      = pars[c4_1548_dzx]
    uv_line_dz      = pars[uv_line_dzx]
    m2_2796_dz      = pars[m2_2796_dzx]
    o2_3727_dz      = pars[o2_3727_dzx]
    o3_5007_dz      = pars[o3_5007_dzx]
    s3_he_dz        = pars[s3_he_dzx]

    # define the emission line amplitude model parameters.
    # those set as ratios are defined relative to the other line in the doublet.
    # the physically allowable range for each ratio is set later in the code.
    # the following ratios are fixed:
    # 1666 = 2.46 * 1660; 5007 = 3 * 4959; 6300 = 3 * 6363; 6585 = 3 * 6550; 9532 = 2.48 * 9069
    la_1216_amp = pars[la_1216_idx]
    h2_1640_amp = pars[h2_1640_idx]
    n5_1238_amp = pars[n5_1238_idx]
    n5_1242_amp = n5_1238_amp / pars[n5_1242_idx]
    c4_1548_amp = pars[c4_1548_idx]
    c4_1550_amp = c4_1548_amp / pars[c4_1550_idx]
    o3_1660_amp = pars[o3_1660_idx]
    o3_1666_amp = o3_1660_amp / pars[o3_1666_idx]
    s3_1883_amp = pars[s3_1883_idx]
    s3_1892_amp = s3_1883_amp / pars[s3_1892_idx]
    c3_1907_amp = pars[c3_1907_idx]
    c3_1909_amp = c3_1907_amp / pars[c3_1909_idx]
    m2_2796_amp = pars[m2_2796_idx]
    m2_2803_amp = m2_2796_amp / pars[m2_2803_idx]
    o2_3727_amp = pars[o2_3727_idx]
    o2_3730_amp = o2_3727_amp / pars[o2_3730_idx]
    hg_4342_amp = pars[hg_4342_idx]
    o3_4363_amp = pars[o3_4363_idx]
    h2_4686_amp = pars[h2_4686_idx]
    hb_4863_amp = pars[hb_4863_idx]
    o3_4959_amp = pars[o3_4959_idx]
    o3_5007_amp = o3_4959_amp / pars[o3_5007_idx]
    o1_6300_amp = pars[o1_6300_idx]
    o1_6363_amp = o1_6300_amp / pars[o1_6363_idx]
    ha_6565_amp = pars[ha_6565_idx]
    n2_6550_amp = pars[n2_6550_idx]
    n2_6585_amp = n2_6550_amp / pars[n2_6585_idx]
    s2_6716_amp = pars[s2_6716_idx]
    s2_6731_amp = s2_6716_amp / pars[s2_6731_idx]
    s3_9069_amp = pars[s3_9069_idx]
    s3_9532_amp = s3_9069_amp / pars[s3_9532_idx]
    he10830_amp = pars[he10830_idx]
    # lyman alpha wing.
    la_wing_amp = pars[la_wing_amp_idx]
    la_wing_sig = pars[la_wing_sig_idx]

    # define the observed wavelengths
    la_1216_obs = la_1216_vac * (1 + z + la_1216_dz)
    n5_1238_obs = n5_1238_vac * (1 + z + la_1216_dz)
    n5_1242_obs = n5_1242_vac * (1 + z + la_1216_dz)
    c4_1548_obs = c4_1548_vac * (1 + z + c4_1548_dz)
    c4_1550_obs = c4_1550_vac * (1 + z + c4_1548_dz)
    h2_1640_obs = h2_1640_vac * (1 + z + uv_line_dz)
    o3_1660_obs = o3_1660_vac * (1 + z + uv_line_dz)
    o3_1666_obs = o3_1666_vac * (1 + z + uv_line_dz)
    s3_1883_obs = s3_1883_vac * (1 + z + uv_line_dz)
    s3_1892_obs = s3_1892_vac * (1 + z + uv_line_dz)
    c3_1907_obs = c3_1907_vac * (1 + z + uv_line_dz)
    c3_1909_obs = c3_1909_vac * (1 + z + uv_line_dz)
    m2_2796_obs = m2_2796_vac * (1 + z + m2_2796_dz)
    m2_2803_obs = m2_2803_vac * (1 + z + m2_2796_dz)
    o2_3727_obs = o2_3727_vac * (1 + z + o2_3727_dz)
    o2_3730_obs = o2_3730_vac * (1 + z + o2_3727_dz)
    hg_4342_obs = hg_4342_vac * (1 + z + o3_5007_dz)
    o3_4363_obs = o3_4363_vac * (1 + z + o3_5007_dz)
    h2_4686_obs = h2_4686_vac * (1 + z + o3_5007_dz)
    hb_4863_obs = hb_4863_vac * (1 + z + o3_5007_dz)
    o3_4959_obs = o3_4959_vac * (1 + z + o3_5007_dz)
    o3_5007_obs = o3_5007_vac * (1 + z + o3_5007_dz)
    o1_6300_obs = o1_6300_vac * (1 + z)
    o1_6363_obs = o1_6363_vac * (1 + z)
    n2_6550_obs = n2_6550_vac * (1 + z)
    ha_6565_obs = ha_6565_vac * (1 + z)
    n2_6585_obs = n2_6585_vac * (1 + z)
    s2_6716_obs = s2_6716_vac * (1 + z)
    s2_6731_obs = s2_6731_vac * (1 + z)
    s3_9069_obs = s3_9069_vac * (1 + z + s3_he_dz)
    s3_9532_obs = s3_9532_vac * (1 + z + s3_he_dz)
    he10830_obs = he10830_vac * (1 + z + s3_he_dz)

    # initialize the continuum and emission line models as lists of zeros.
    cont_model = x * 0.
    line_model = x * 0 # 0 or 0.?

    # initialize the x-values (wave) that will be used as nodes for the spline fit.
    spline_x_values = pars[first_node_index + nnodes : first_node_index + 2 * nnodes]

    # initialize the y-values (flux) that will be used as nodes for the spline fit.
    spline_y_values = pars[first_node_index : first_node_index + nnodes]

    # split the spline model into short wavelength (blue) and long wavelength (red) halves.
    spline_x_values_blue = spline_x_values[np.where(spline_x_values < transition_wave)]
    spline_y_values_blue = spline_y_values[np.where(spline_x_values < transition_wave)]
    spline_x_values_red  = spline_x_values[np.where(spline_x_values >= transition_wave)]
    spline_y_values_red  = spline_y_values[np.where(spline_x_values >= transition_wave)]

    # fitting a cubic spline requires > 4 nodes so append more points if required.
    if len(spline_x_values_blue) > 0:
        i = 0
        while np.size(spline_x_values_blue) < 4:
            i = i + 1
            spline_x_values_blue = np.append(spline_x_values_blue, transition_wave + (1000 * i))
            spline_y_values_blue = np.append(spline_y_values_blue, spline_y_values_blue[-1])

        continuum_spline_model_blue = interpolate.splrep(spline_x_values_blue, spline_y_values_blue, s=0, k=3)

        w = np.where((x > la_1216_obs - fit_region) & (x < transition_wave))
        if np.size(w) > 0:
            cont_model[w] = interpolate.splev(x[w], continuum_spline_model_blue, der=0)

    # fitting a cubic spline requires > 4 nodes so append more points if required.
    # if len(spline_x_values_red) > 0:
    #     i = 0
    #     while np.size(spline_x_values_red) < 4:
    #         i = i + 1
    #         spline_x_values_red = np.append(spline_x_values_red, transition_wave - (1000 * i))
    #         spline_y_values_red = np.append(spline_y_values_red, spline_y_values_red[0])
    #
    #     continuum_spline_model_red = interpolate.splrep(spline_x_values_red, spline_y_values_red, s=0, k=3)
    #
    #     w = np.where((x >= transition_wave) & (x < he10830_obs  + fit_region))
    #     if np.size(w) > 0:
    #         cont_model[w] = interpolate.splev(x[w], continuum_spline_model_red, der=0)
    if len(spline_x_values_red) > 0:
        # Allow for linear fits, which are sometimes best for the grism data.
        # if len(spline_x_values_red) == 1:
        #     continuum_spline_model_red = ([spline_x_values_red, spline_x_values_red], [0.3, 0.3])
        if len(spline_x_values_red) == 2:
            continuum_spline_model_red = interpolate.splrep(spline_x_values_red, spline_y_values_red, s=0, k=1)
        if len(spline_x_values_red) == 3:
            continuum_spline_model_red = interpolate.splrep(spline_x_values_red, spline_y_values_red, s=0, k=2)
        if len(spline_x_values_red) > 3:
            continuum_spline_model_red = interpolate.splrep(spline_x_values_red, spline_y_values_red, s=0, k=3)

        # print 'continuum_spline_model_red =', continuum_spline_model_red
        # print 'type(continuum_spline_model_red) =', type(continuum_spline_model_red)

        w = np.where((x >= transition_wave) & (x < he10830_obs  + fit_region))

        if np.size(w) > 0:
            cont_model[w] = interpolate.splev(x[w], continuum_spline_model_red, der=0)

    ############################################################################
    ############################################################################

    '''
    The line_model that will be fit to the data consists of a flux continuum
    and emission lines. The emission lines are added to the model one at a time
    using the add_emission_line_to_model() function defined below. The function
    extracts a wavelength array centered on the line with a width controlled by
    the 'fit_region' parameter defined in the 'default.config' file. There are
    several groups of emission lines such as [O III]+HB, and HA+[N II]+[S II]
    that have special constraints and are added using specialized code below.
    Closely-spaced doublet lines are added with the add_doublet_lines_to_model()
    function so that both components are fit over identical wavelength intervals.
    '''

    def add_emission_line_to_model(centroid, amplitude):
        if (centroid < transition_wave): sigma = sigma_blue
        if (centroid > transition_wave): sigma = sigma_red
        w = np.where((x > (centroid - fit_region)) & (x < (centroid + fit_region)))
        if (np.size(w) > 0):
            line_model[w] = line_model[w] + \
            (amplitude * gaussian(x[w], centroid, sigma))

    def add_doublet_lines_to_model(centroid_1, amplitude_1, centroid_2, amplitude_2):
        if (centroid_1 < transition_wave): sigma = sigma_blue
        if (centroid_1 > transition_wave): sigma = sigma_red
        w = np.where((x > (centroid_1 - fit_region)) & (x < (centroid_2 + fit_region)))
        if (np.size(w) > 0):
            line_model[w] = line_model[w] + \
            (amplitude_1 * gaussian(x[w], centroid_1, sigma)) + \
            (amplitude_2 * gaussian(x[w], centroid_2, sigma))

    ############################################################################

    #add_emission_line_to_model(la_1216_obs, la_1216_amp)
    add_doublet_lines_to_model(n5_1238_obs, n5_1238_amp, n5_1242_obs, n5_1242_amp)
    add_doublet_lines_to_model(c4_1548_obs, c4_1548_amp, c4_1550_obs, c4_1550_amp)
    add_emission_line_to_model(h2_1640_obs, h2_1640_amp)
    add_doublet_lines_to_model(o3_1660_obs, o3_1660_amp, o3_1666_obs, o3_1666_amp)
    add_doublet_lines_to_model(s3_1883_obs, s3_1883_amp, s3_1892_obs, s3_1892_amp)
    add_doublet_lines_to_model(c3_1907_obs, c3_1907_amp, c3_1909_obs, c3_1909_amp)
    add_doublet_lines_to_model(m2_2796_obs, m2_2796_amp, m2_2803_obs, m2_2803_amp)
    add_doublet_lines_to_model(o2_3727_obs, o2_3727_amp, o2_3730_obs, o2_3730_amp)
    add_emission_line_to_model(hg_4342_obs, hg_4342_amp)
    add_emission_line_to_model(o3_4363_obs, o3_4363_amp)
    add_emission_line_to_model(h2_4686_obs, h2_4686_amp)
    add_doublet_lines_to_model(o1_6300_obs, o1_6300_amp, o1_6363_obs, o1_6363_amp)
    add_emission_line_to_model(s3_9069_obs, s3_9069_amp)
    add_emission_line_to_model(s3_9532_obs, s3_9532_amp)
    add_emission_line_to_model(he10830_obs, he10830_amp)

    ############################################################################

    # lyman alpha is highly asymmetric so add a red wing component.
    # main.
    if (la_1216_obs < transition_wave): sigma = sigma_blue
    if (la_1216_obs > transition_wave): sigma = sigma_red
    w = np.where((x > (la_1216_obs - fit_region)) & (x < (la_1216_obs + fit_region)))
    if (np.size(w) > 0):
        line_model[w] = line_model[w] + \
        (la_1216_amp * gaussian(x[w], la_1216_obs, sigma))

    # wing.
    w = np.where((x > (la_1216_obs)) & (x < (la_1216_obs + fit_region)))
    if (np.size(w) > 0):
        line_model[w] = line_model[w] + \
        (la_wing_amp * gaussian(x[w], la_1216_obs, la_wing_sig*sigma))

    ############################################################################

    # The [O III] and Hb lines are fit simultaneously in a single wavelength range.
    if (hb_4863_obs < transition_wave): sigma = sigma_blue
    if (hb_4863_obs > transition_wave): sigma = sigma_red
    w = np.where((x > (hb_4863_obs - fit_region)) & (x < (o3_5007_obs + fit_region)))
    if np.size(w) > 0:
        line_model[w] = line_model[w] + \
        o3_5007_amp * gaussian(x[w], o3_5007_obs, sigma) + \
        o3_4959_amp * gaussian(x[w], o3_4959_obs, sigma) + \
        hb_4863_amp * gaussian(x[w], hb_4863_obs, sigma)

    ############################################################################

    # The Ha, [N II], and [S II] lines are fit simultaneously in a single wavelength range.
    if (ha_6565_obs < transition_wave): sigma = sigma_blue
    if (ha_6565_obs > transition_wave): sigma = sigma_red
    w = np.where((x > (n2_6550_obs - fit_region)) & (x < (s2_6731_obs + fit_region)))
    if np.size(w) > 0:
        line_model[w] = line_model[w] + \
        ha_6565_amp * gaussian(x[w], ha_6565_obs, sigma) + \
        n2_6550_amp * gaussian(x[w], n2_6550_obs, sigma) + \
        n2_6585_amp * gaussian(x[w], n2_6585_obs, sigma) + \
        s2_6716_amp * gaussian(x[w], s2_6716_obs, sigma) + \
        s2_6731_amp * gaussian(x[w], s2_6731_obs, sigma)

    ############################################################################

    model = cont_model + line_model

    return model

    ############################################################################
    ############################################################################


def model_resid(pars, fjac = None, lam = None, flux = None, err = None):
    model = emissionline_model(pars, lam)
    resid = (flux - model) / err
    status = 0
    return [status, resid]


def fit_obj(input_list):
    print('\nRunning fit_obj...\n')
    lam_spec    = input_list[0]
    flux_spec   = input_list[1]
    error_spec  = input_list[2]
    config_pars = input_list[3]
    z_in        = input_list[4] # the best-guess redshift or the current line-list redshift
    fwhm_guess  = input_list[5] # red fwhm guess
    beam_name   = input_list[6] # object id.
    ab_image_max= input_list[7] # MDR 2022/06/30
    fast_fit    = input_list[8] # MDR 2022/06/30
    line_mask   = config_pars['line_mask']
    fit_region  = config_pars['fit_region']
    disp_red    = config_pars['dispersion_red']
    disp_blue   = config_pars['dispersion_blue']
    node_wave   = config_pars['node_wave']
    dz          = config_pars['delta_z'] # this is the maximum shift from the line list redshift/best estimate.
    transition_wave = config_pars['transition_wave']
    nnodes      = len(node_wave)

    # the fluxes are divided by a scale factor for easier handling by mpfit.
    scl        = 1.0e-18
    flux_spec  = (flux_spec  / scl)
    error_spec = (error_spec / scl)

    # calculate the approximate line locations based on the current guess of 'z_in'.
    la_1216_obs = la_1216_vac * (1 + z_in)
    n5_1238_obs = n5_1238_vac * (1 + z_in)
    n5_1242_obs = n5_1242_vac * (1 + z_in)
    c4_1548_obs = c4_1548_vac * (1 + z_in)
    c4_1550_obs = c4_1550_vac * (1 + z_in)
    h2_1640_obs = h2_1640_vac * (1 + z_in)
    o3_1660_obs = o3_1660_vac * (1 + z_in)
    o3_1666_obs = o3_1666_vac * (1 + z_in)
    s3_1883_obs = s3_1883_vac * (1 + z_in)
    s3_1892_obs = s3_1892_vac * (1 + z_in)
    c3_1907_obs = c3_1907_vac * (1 + z_in)
    c3_1909_obs = c3_1909_vac * (1 + z_in)
    m2_2796_obs = m2_2796_vac * (1 + z_in)
    m2_2803_obs = m2_2803_vac * (1 + z_in)
    o2_3727_obs = o2_3727_vac * (1 + z_in)
    o2_3730_obs = o2_3730_vac * (1 + z_in)
    hg_4342_obs = hg_4342_vac * (1 + z_in)
    o3_4363_obs = o3_4363_vac * (1 + z_in)
    h2_4686_obs = h2_4686_vac * (1 + z_in)
    hb_4863_obs = hb_4863_vac * (1 + z_in)
    o3_4959_obs = o3_4959_vac * (1 + z_in)
    o3_5007_obs = o3_5007_vac * (1 + z_in)
    o1_6300_obs = o1_6300_vac * (1 + z_in)
    o1_6363_obs = o1_6363_vac * (1 + z_in)
    n2_6550_obs = n2_6550_vac * (1 + z_in)
    ha_6565_obs = ha_6565_vac * (1 + z_in)
    n2_6585_obs = n2_6585_vac * (1 + z_in)
    s2_6716_obs = s2_6716_vac * (1 + z_in)
    s2_6731_obs = s2_6731_vac * (1 + z_in)
    s3_9069_obs = s3_9069_vac * (1 + z_in)
    s3_9532_obs = s3_9532_vac * (1 + z_in)
    he10830_obs = he10830_vac * (1 + z_in)

    ############################################################################
    ############################################################################
    '''
    First, we fit only the continuum with the emission lines masked out so we
    have an accurate first guess when fitting the entire spectrum with emission
    lines. It may be possible to skip this step, but doing this for each object
    ensures a robust fit. The mask_emission_lines() function below is used to
    mask all of the emission lines that are being fit. The code will mask from
    'blue_line' to 'red_line'. This can be used to mask a single emission line,
    or to mask groups of lines. For example, to mask out HA + [N II] + [S II]
    we can call 'mask_emission_lines(n2_6550_obs, s2_6731_obs)'.
    '''

    def mask_emission_lines(blue_line, red_line):
        w = np.where((lam_spec > (blue_line - line_mask)) & (lam_spec < (red_line + line_mask)))
        mask_spec[w] = 1.0

    ############################################################################

    # First set the mask equal to zero everywhere and mask regions outside the bluest and reddest line.
    mask_spec = (lam_spec * 0.0)
    w = np.where((lam_spec < la_1216_obs - fit_region) | (lam_spec > he10830_obs + fit_region))
    mask_spec[w] = 1.0

    # only mask the strong lines so there is sufficient continuum.
    # testing showed that commenting out mask_emission_lines(o1_6300_obs, o1_6363_obs)
    # made the code run much more slowly and inexplicably changed the emission
    # line S/N estimates by very large factors for one test object, so change with care.
    mask_emission_lines(la_1216_obs, la_1216_obs)
    #mask_emission_lines(n5_1238_obs, n5_1242_obs)
    mask_emission_lines(c4_1548_obs, c4_1550_obs)
    mask_emission_lines(h2_1640_obs, h2_1640_obs)
    mask_emission_lines(o3_1660_obs, o3_1666_obs)
    mask_emission_lines(s3_1883_obs, s3_1892_obs)
    mask_emission_lines(c3_1907_obs, c3_1909_obs)
    mask_emission_lines(m2_2796_obs, m2_2803_obs)
    mask_emission_lines(o2_3727_obs, o2_3730_obs)
    mask_emission_lines(hg_4342_obs, o3_4363_obs)
    #mask_emission_lines(h2_4686_obs, h2_4686_obs)
    mask_emission_lines(hb_4863_obs, o3_5007_obs)
    mask_emission_lines(o1_6300_obs, o1_6363_obs)
    mask_emission_lines(n2_6550_obs, s2_6731_obs)
    mask_emission_lines(s3_9069_obs, s3_9069_obs)
    mask_emission_lines(s3_9532_obs, s3_9532_obs)
    mask_emission_lines(he10830_obs, he10830_obs)

    w = np.where(mask_spec == 0.0)
    cont_guess = np.median(flux_spec[w])

    ############################################################################
    ############################################################################

    '''
    Next, we must pass the model parameter guesses into the 'pguess' list, which
    is provided to mpfit(). In the first pass we're only interested in obtaining
    an accurate estimate of the continuum with the emission lines masked, so we
    can set the amplitude guesses to zero.
    '''

    # initialize a list with a length equal to the number of model parameters.
    pguess_cont = np.zeros(first_node_index + 2 * nnodes)

    # initialize all non-continuum values to zero (amplitudes and dz's).
    pguess_cont[0 : first_node_index] = 0

    # define the fixed parameters.
    pguess_cont[nnodes_idx]     = nnodes
    pguess_cont[fitreg_idx]     = fit_region
    pguess_cont[t_wave_idx]     = transition_wave
    pguess_cont[z_idx]          = z_in
    pguess_cont[fwhm_ratio_idx] = 0.058 # 0.027 # blue/red fwhm ratio
    pguess_cont[fwhm_grism_idx] = fwhm_guess # fwhm_red.

    # pass the continuum guess y-values into the pguess_cont list.
    pguess_cont[first_node_index : first_node_index + nnodes] = cont_guess

    # pass the continuum guess x-values into the pguess_cont list.
    pguess_cont[first_node_index + nnodes : first_node_index + 2 * nnodes] = node_wave

    # the lines are not fit but set the ratio parameters to unity to avoid division by zero.
    for idx in ratio_indices:
        pguess_cont[idx] = 1.0

    # determine the number of parameters for the model.
    num_of_params = len(pguess_cont)

    # initialize the parinfo_cont list with the arguments required by mpfit.
    parinfo_cont = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]} for i in range(num_of_params)]

    # fill the parinfo_cont list for mpfit with the values defined above.
    for i in range(num_of_params): parinfo_cont[i]['value'] = pguess_cont[i]

    # the amplitudes are zero and we want to fix them for the continuum fit.
    for j in range(first_node_index): parinfo_cont[j]['fixed'] = 1

    # the continuum x-values for the nodes should be fixed.
    for i in range(first_node_index + nnodes, first_node_index + 2 * nnodes): parinfo_cont[i]['fixed'] = 1

    # initialize the wavelength, flux, and error function arguments.
    fa_cont = {'lam':lam_spec[w], 'flux':flux_spec[w], 'err':error_spec[w]}

    # call mpfit for the continuum-only fit model.
    out = mpfit(model_resid, pguess_cont, functkw = fa_cont, parinfo = parinfo_cont, quiet = True)

    ############################################################################
    ############################################################################

    '''
    Next, we want to repeat the spectral fit using the best fit model continuum
    from the first run of mpfit. We only attempt the second iteration if the
    first run on the continuum was successful (i.e. mpfit returns 'status' > 0).
    First, we use the estimate_emission_line_amplitudes() function to estimate
    the amplitude of each emission line and add them to the model guesses.
    '''

    if out.status > 0:

        def estimate_emission_line_amplitudes(centroid):
            if ((centroid > np.min(lam_spec)) & (centroid < np.max(lam_spec))):
                if (centroid < transition_wave): disp = disp_blue
                if (centroid > transition_wave): disp = disp_red
                wline = np.where((lam_spec > (centroid - 5.0 * disp)) & (lam_spec < (centroid + 5.0 * disp)))
                if (np.size(wline) > 0):
                    peakval = np.max(flux_spec[wline])
                    amplitude_guess = (peakval - emissionline_model(out.params, np.array([centroid])))[0] # Peak minus continuum.
                    if (amplitude_guess < 0.0):
                        amplitude_guess = 0.0
                else: amplitude_guess = 1.0
            else: amplitude_guess = 1.0 # Set to unity if line is not being fit.

            return amplitude_guess

        ########################################################################

        '''
        The mpfit() documentation provides a detailed description of the parinfo
        array. There is a very brief description of each below for completeness:

        'value':   > the parameter value.
        'fixed':   > set 0 = False or 1 = True.
        'limited': > set 0 = False or 1 = True. (specify for [lower_bound, upper_bound])
        'limits':  > set lower_bound and/or upper_bound limits, must match 'limited'.
        'tied':    > fixes one parameter in terms of another as a string expression.
        '''

        # initialize the parinfo array for mpfit with zeros for all parameters.
        parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.], 'tied':''} for i in range(num_of_params)]

        # define the fixed parameters.
        parinfo[nnodes_idx]['value'] = nnodes
        parinfo[nnodes_idx]['fixed'] = 1
        parinfo[fitreg_idx]['value'] = fit_region
        parinfo[fitreg_idx]['fixed'] = 1
        parinfo[t_wave_idx]['value'] = transition_wave
        parinfo[t_wave_idx]['fixed'] = 1

        # allow fwhm_ratio to vary because extended sources in the grism will have broadened profiles.
        parinfo[fwhm_ratio_idx]['value'] = 0.058 # 0.027 # ratio of the MUSE/WFC3 data dispersions.
        parinfo[fwhm_ratio_idx]['fixed'] = 0
        parinfo[fwhm_ratio_idx]['limited'][0] = 1
        parinfo[fwhm_ratio_idx]['limited'][1] = 1
        parinfo[fwhm_ratio_idx]['limits'][0]  = (0.5 * 0.058) # currently 2x. grism bigger = ratio lower.
        parinfo[fwhm_ratio_idx]['limits'][1]  = (1.0 * 0.058) # currently 1x.

        parinfo[z_idx]['value'] = z_in
        parinfo[z_idx]['limited'][0] = 1
        parinfo[z_idx]['limited'][1] = 1
        parinfo[z_idx]['limits'][0]  = z_in - dz
        parinfo[z_idx]['limits'][1]  = z_in + dz

        parinfo[fwhm_grism_idx]['value'] = fwhm_guess
        parinfo[fwhm_grism_idx]['limited'][0] = 1
        parinfo[fwhm_grism_idx]['limited'][1] = 1
        parinfo[fwhm_grism_idx]['limits'][0]  = config_pars['min_fwhm_scl'] * fwhm_guess
        parinfo[fwhm_grism_idx]['limits'][1]  = config_pars['max_fwhm_scl'] * fwhm_guess

        ########################################################################

        # set the 'dz' parameter values for each line with an allowed shift.
        # the values were already initialized to zero in 'parinfo =' above.

        loop_dzx = ['la_1216_dzx', 'c4_1548_dzx', 'uv_line_dzx', \
        'm2_2796_dzx', 'o2_3727_dzx', 'o3_5007_dzx', 's3_he_dzx']

        for dzx in loop_dzx:
            parinfo[eval(dzx)]['limited'][0] = 1
            parinfo[eval(dzx)]['limited'][1] = 1
            parinfo[eval(dzx)]['limits'][0]  = -1.0*dz
            parinfo[eval(dzx)]['limits'][1]  = dz

        ########################################################################

        # determine when the 'dz' parameter for each line can vary relative to Ha.

        parinfo[la_1216_dzx]['fixed'] = 0 # lya is highly asymmetric so always let it vary.
        parinfo[c4_1548_dzx]['fixed'] = 0 # the far uv lines are always in MUSE so allow them to vary.
        parinfo[uv_line_dzx]['fixed'] = 0 # the far uv lines are always in MUSE so allow them to vary.
        parinfo[s3_he_dzx]['fixed']   = 1 # previous version said always fix dz_siii because in grism.

        if (m2_2796_obs >= transition_wave):
            parinfo[m2_2796_dzx]['fixed'] = 1 # fix dz to zero if all lines are in the grism.
        else:
            parinfo[m2_2796_dzx]['fixed'] = 0 # allow dz to vary if all lines are not in the grism.

        if (o2_3727_obs >= transition_wave):
            parinfo[o2_3727_dzx]['fixed'] = 1 # fix dz to zero if all lines are in the grism.
        else:
            parinfo[o2_3727_dzx]['fixed'] = 0 # allow dz to vary if all lines are not in the grism.

        if (hb_4863_obs >= transition_wave):
            parinfo[o3_5007_dzx]['fixed'] = 1 # fix dz to zero if all lines are in the grism.
        else:
            parinfo[o3_5007_dzx]['fixed'] = 0 # allow dz to vary if all lines are not in the grism.

        ########################################################################

        # fill the parinfo array for mpfit with the amplitude guesses.
        loop_lines = [\
        'la_1216', 'n5_1238', 'c4_1548', 'h2_1640', 'o3_1660', \
        's3_1883', 'c3_1907', 'm2_2796', 'o2_3727', 'hg_4342', \
        'o3_4363', 'h2_4686', 'hb_4863', 'o3_4959', 'o1_6300', \
        'ha_6565', 'n2_6550', 's2_6716', 's3_9069', 'he10830']

        # fix the lower limits to zero because the gaussian amplitudes must be positive.
        for line in loop_lines:
            parinfo[eval(line+'_idx')]['value'] = estimate_emission_line_amplitudes(eval(line+'_obs'))
            parinfo[eval(line+'_idx')]['limited'][0] = 1
            parinfo[eval(line+'_idx')]['limits'][0]  = 0

        # fill the parinfo array for mpfit with the amplitude guesses.
        # guess 6585 is 10% of Ha based on literature for consistency.
        parinfo[n2_6585_idx]['value'] = (estimate_emission_line_amplitudes(ha_6565_obs) / 10.0)

        ########################################################################

        # fix the upper and lower limits for ratios to their physically allowed values.
        loop_ratios = ['n5_1242', 'c4_1550', 's3_1892', 'c3_1909', 'm2_2803', 'o2_3730', 's2_6731']

        # allow the ratios to vary and set the upper and lower limits.
        for line in loop_ratios:
            parinfo[eval(line+'_idx')]['fixed']      = 0 # allowed to vary.
            parinfo[eval(line+'_idx')]['limited'][0] = 1 # lower limit = True.
            parinfo[eval(line+'_idx')]['limited'][1] = 1 # upper limit = True.

        parinfo[n5_1242_idx]['value']     = 1.95 # based on a grid of cloudy models.
        parinfo[n5_1242_idx]['limits'][0] = 1.15
        parinfo[n5_1242_idx]['limits'][1] = 2.00

        parinfo[c4_1550_idx]['value']     = 1.98 # based on a grid of cloudy models.
        parinfo[c4_1550_idx]['limits'][0] = 1.22
        parinfo[c4_1550_idx]['limits'][1] = 2.00

        parinfo[s3_1892_idx]['value']     = 0.50 # based on a grid of cloudy models.
        parinfo[s3_1892_idx]['limits'][0] = 0.00005
        parinfo[s3_1892_idx]['limits'][1] = 1.69

        parinfo[c3_1909_idx]['value']     = 0.92 # based on a grid of cloudy models.
        parinfo[c3_1909_idx]['limits'][0] = 0.00015
        parinfo[c3_1909_idx]['limits'][1] = 1.54

        parinfo[m2_2803_idx]['value']     = 1.98 # based on a grid of cloudy models.
        parinfo[m2_2803_idx]['limits'][0] = 1.26
        parinfo[m2_2803_idx]['limits'][1] = 2.00

        parinfo[o2_3730_idx]['value']     = 1.0 / 1.50
        parinfo[o2_3730_idx]['limits'][0] = 1.0 / 1.50 # 3730/3727 lower limit ratio is reverse of osterbrock plot so 1/1.50 ~ 0.66
        parinfo[o2_3730_idx]['limits'][1] = 1.0 / 0.38 # 3730/3727 upper limit ratio is reverse of osterbrock plot so 1/0.38 ~ 2.63

        parinfo[s2_6731_idx]['value']     = 1.40
        parinfo[s2_6731_idx]['limits'][0] = 0.42 # 6716/6731 = 0.44 lower limit for T = 10,000 K (0.42 in high T ~ 20,000 limit)
        parinfo[s2_6731_idx]['limits'][1] = 1.47 # 6716/6731 = 1.44 upper limit for T = 10,000 K (1.47 in low T ~ 5,000 limit)

        ########################################################################

        # lyman alpha wing.
        parinfo[la_wing_amp_idx]['value'] = estimate_emission_line_amplitudes(la_1216_obs) / 5.0 # just a guess.
        parinfo[la_wing_amp_idx]['limited'][0] = 1
        parinfo[la_wing_amp_idx]['limits'][0]  = 0
        parinfo[la_wing_sig_idx]['value'] = 5.0 # guess 5x the default line width guess.
        parinfo[la_wing_sig_idx]['fixed'] = 0
        parinfo[la_wing_sig_idx]['limited'][0] = 1
        parinfo[la_wing_sig_idx]['limited'][1] = 1
        parinfo[la_wing_sig_idx]['limits'][0]  = 1.000
        parinfo[la_wing_sig_idx]['limits'][1]  = 100.0

        ########################################################################

        # define emission lines that are fixed relative to their stronger doublet for mpfit.
        # the amplitudes of the parent lines are fixed to be positive so does not need to be specified.

        parinfo[o3_1666_idx]['value'] = (1.0 / 2.46) # blue line / red line
        parinfo[o3_1666_idx]['fixed'] = 1

        parinfo[o3_5007_idx]['value'] = (1.0 / 3.0) # blue line / red line
        parinfo[o3_5007_idx]['fixed'] = 1

        parinfo[o1_6363_idx]['value'] = (3.0 / 1.0) # blue line / red line
        parinfo[o1_6363_idx]['fixed'] = 1

        parinfo[n2_6585_idx]['value'] = (1.0 / 3.0) # blue line / red line
        parinfo[n2_6585_idx]['fixed'] = 1

        '''
        If Ha and [N II] are in the grism then lock the flux ratio because the
        lines are not resolved. As noted in earlier code versions, "do not even
        think you should change this for grism data". Indeed, testing showed that
        allowing this to vary often puts all of the flux into [N II], which is
        unphysical. Following the model of Henry et al. 2021, ApJ, 919, 143, we
        fix the flux of [N II] 6585 to 10% of Ha. This corresponds to 13% of the
        total flux including both [N II] lines. This overall contribution is in
        agreement with Erb et al. 2006, ApJ, 644, 813 (Table 2) for moderate mass
        galaxies. The range is ~13 - 36%, increasing with host galaxy stellar mass.
        '''

        # 6550 is 1/30 of total, 6585 is 3/30 of total, thus [N II] = 4/30 = 0.13 of the total flux.
        if (ha_6565_obs > transition_wave):
            parinfo[n2_6550_idx]['tied']  = 'p['+str(ha_6565_idx)+'] / 30.0'

        # if the [S III] lines are split between MUSE and WFC3 then scale the fixed height ratio by the difference in fwhm.
        # while the below scaling is correct, the gap between the MUSE and WFC3 (~668 Ang) is larger than the separation between
        # the [S III] lines (~463 Ang), thus it is not possible to have 9069 in MUSE and 9532 in WFC3 so the scaling is ommitted.
        # if ((s3_9532_obs > transition_wave) & (s3_9069_obs <= transition_wave)):
        #     parinfo[s3_9532_idx]['tied']  = 'p['+str(s3_9069_idx)+'] * 2.48 * p['+str(fwhm_ratio_idx)+']'
        # else:
        #     parinfo[s3_9532_idx]['tied']  = 'p['+str(s3_9069_idx)+'] * 2.48'
        parinfo[s3_9532_idx]['value'] = (1.0 / 2.48) # blue line / red line
        parinfo[s3_9532_idx]['fixed'] = 1

        ########################################################################

        # pass the continuum node values from the first pass fit to the current array.
        for i in range(first_node_index, first_node_index + 2 * nnodes): parinfo[i]['value'] = out.params[i]

        # the continuum x-values for the nodes should be fixed.
        for i in range(first_node_index+nnodes, first_node_index + 2 * nnodes): parinfo[i]['fixed'] = 1

        # define the wavelength range to fit.
        w = np.where((lam_spec > (la_1216_obs - fit_region)) & (lam_spec < (he10830_obs + fit_region)))

        # initialize the wavelength, flux, and error function arguments.
        fa = {'lam':lam_spec[w], 'flux':flux_spec[w], 'err':error_spec[w]}

        pguess = [] # define the pguess list.

        # pass the parinfo information from above into the pguess array for mpfit.
        for i in range(num_of_params): pguess.append(parinfo[i]['value'])

        start_time  = default_timer() # time the fit for performance analyses.

        if fast_fit == True:
            print 'Running fast line fit...'

            # for a fast fit fix the dispersion ratio to its initial guess.
            parinfo[fwhm_ratio_idx]['fixed'] = 1

            # for a fast fit fix all varying line ratios to their initial guess.
            for line in loop_ratios:
                parinfo[eval(line+'_idx')]['fixed'] = 1

            # for a fast fit fix all dz's to their initial guess of zero.
            for dzx in loop_dzx:
                parinfo[eval(dzx)]['fixed'] = 1

        else:
            print 'Running full line fit...'

        out = mpfit(model_resid, pguess, functkw = fa, parinfo = parinfo, quiet = True, maxiter = 100)

        end_time = default_timer() # time the fit for performance analyses.

        fit_time = np.format_float_positional(end_time - start_time, precision = 2)

        print '\033[92m' + 'MPFIT completed', out.niter, 'iterations in', fit_time, 'seconds.' + '\033[0m'

        chisq = out.fnorm / (len(w[0]) - num_of_params)

        ########################################################################

        # evaluate continuum spline.
        modelpars_nolines = cp.deepcopy(out.params)

        # this is needed to evaluate the continuum for fluxes below.
        modelpars_nolines[first_line_index:first_node_index] = 0.

        # the lines are not fit but set the ratio parameters to unity to avoid division by zero.
        for idx in ratio_indices:
            modelpars_nolines[idx] = 1.0

        # re-evaluate the observed wavelengths based on the refined z estimate.
        z_out = out.params[z_idx]

        la_1216_obs = la_1216_vac * (1 + z_out)
        n5_1238_obs = n5_1238_vac * (1 + z_out)
        n5_1242_obs = n5_1242_vac * (1 + z_out)
        c4_1548_obs = c4_1548_vac * (1 + z_out)
        c4_1550_obs = c4_1550_vac * (1 + z_out)
        h2_1640_obs = h2_1640_vac * (1 + z_out)
        o3_1660_obs = o3_1660_vac * (1 + z_out)
        o3_1666_obs = o3_1666_vac * (1 + z_out)
        s3_1883_obs = s3_1883_vac * (1 + z_out)
        s3_1892_obs = s3_1892_vac * (1 + z_out)
        c3_1907_obs = c3_1907_vac * (1 + z_out)
        c3_1909_obs = c3_1909_vac * (1 + z_out)
        m2_2796_obs = m2_2796_vac * (1 + z_out)
        m2_2803_obs = m2_2803_vac * (1 + z_out)
        o2_3727_obs = o2_3727_vac * (1 + z_out)
        o2_3730_obs = o2_3730_vac * (1 + z_out)
        hg_4342_obs = hg_4342_vac * (1 + z_out)
        o3_4363_obs = o3_4363_vac * (1 + z_out)
        h2_4686_obs = h2_4686_vac * (1 + z_out)
        hb_4863_obs = hb_4863_vac * (1 + z_out)
        o3_4959_obs = o3_4959_vac * (1 + z_out)
        o3_5007_obs = o3_5007_vac * (1 + z_out)
        o1_6300_obs = o1_6300_vac * (1 + z_out)
        o1_6363_obs = o1_6363_vac * (1 + z_out)
        n2_6550_obs = n2_6550_vac * (1 + z_out)
        ha_6565_obs = ha_6565_vac * (1 + z_out)
        n2_6585_obs = n2_6585_vac * (1 + z_out)
        s2_6716_obs = s2_6716_vac * (1 + z_out)
        s2_6731_obs = s2_6731_vac * (1 + z_out)
        s3_9069_obs = s3_9069_vac * (1 + z_out)
        s3_9532_obs = s3_9532_vac * (1 + z_out)
        he10830_obs = he10830_vac * (1 + z_out)

        ########################################################################
        ########################################################################

        '''
        The blocks of code below calculate the flux, error, and equivalent
        width for each emission line. The below function determines these for
        each line and adds them to the output. There are groups of lines such
        as HA+[N II] that have special constraints and are handled separately.
        '''

        def calculate_emission_line_flux(centroid, amplitude_index, rest_wave):
            with np.errstate(invalid = 'ignore', divide = 'ignore'):
                if ((centroid > np.min(lam_spec)) & (centroid < np.max(lam_spec))):
                    if (centroid < transition_wave):
                        sig     = ((out.params[fwhm_grism_idx] * out.params[fwhm_ratio_idx]) / 2.3548) # red fwhm * blue/red fwhm ratio.
                        sig_err = ((out.perror[fwhm_grism_idx] * out.params[fwhm_ratio_idx]) / 2.3548)
                    if (centroid > transition_wave):
                        sig     = (out.params[fwhm_grism_idx] / 2.3548) # red fwhm
                        sig_err = (out.perror[fwhm_grism_idx] / 2.3548)
                    covar_term  = (2.0 * (out.covar[fwhm_grism_idx][amplitude_index] / (out.params[fwhm_grism_idx] * out.params[amplitude_index])))
                    # The emission line flux is 'line_flux' = (sqrt(2*pi) * (height) * (sigma)).
                    line_flux = (np.sqrt(2.0 * math.pi) * (out.params[amplitude_index]) * (sig))
                    if (line_flux > 0.0):
                        line_err = line_flux * (np.sqrt(((out.perror[amplitude_index] / out.params[amplitude_index]) ** 2.0) + ((sig_err / sig) ** 2.0) + covar_term))
                    else:
                        w = np.where((lam_spec > (centroid - fwhm_guess)) & (lam_spec < (centroid + fwhm_guess)))
                        line_err = np.sqrt(np.sum(error_spec[w] ** 2.0))
                    line_cont   = emissionline_model(modelpars_nolines, rest_wave * np.array([1.0 + out.params[z_idx]]))[0] # redshift.
                    line_ew_obs = (line_flux / line_cont)
                else:
                    line_flux   = (-1 / scl)
                    line_err    = (-1 / scl)
                    line_ew_obs =  -1

            return line_flux, line_err, line_ew_obs

        ############################################################################

        # calculate the emission line fluxes and return them to measure_z_interactive().
        h2_1640_flux, h2_1640_err, h2_1640_ew_obs = calculate_emission_line_flux(h2_1640_obs, h2_1640_idx, h2_1640_vac)
        hg_4342_flux, hg_4342_err, hg_4342_ew_obs = calculate_emission_line_flux(hg_4342_obs, hg_4342_idx, hg_4342_vac)
        o3_4363_flux, o3_4363_err, o3_4363_ew_obs = calculate_emission_line_flux(o3_4363_obs, o3_4363_idx, o3_4363_vac)
        h2_4686_flux, h2_4686_err, h2_4686_ew_obs = calculate_emission_line_flux(h2_4686_obs, h2_4686_idx, h2_4686_vac)
        hb_4863_flux, hb_4863_err, hb_4863_ew_obs = calculate_emission_line_flux(hb_4863_obs, hb_4863_idx, hb_4863_vac)
        n2_6550_flux, n2_6550_err, n2_6550_ew_obs = calculate_emission_line_flux(n2_6550_obs, n2_6550_idx, n2_6550_vac) # tied to 1/3 of 6585.
        ha_6565_flux, ha_6565_err, ha_6565_ew_obs = calculate_emission_line_flux(ha_6565_obs, ha_6565_idx, ha_6565_vac)
        n2_6585_flux, n2_6585_err, n2_6585_ew_obs = calculate_emission_line_flux(n2_6585_obs, n2_6585_idx, n2_6585_vac)
        he10830_flux, he10830_err, he10830_ew_obs = calculate_emission_line_flux(he10830_obs, he10830_idx, he10830_vac)

        '''
        Calculate the flux, error, and equivalent width values for halpha and
        [N II] separately, which is required because [N II] is 'tied' to halpha
        in the grism, and mpfit does not return uncertainties for 'tied' params.
        '''

        if (ha_6565_flux > 0.0): # if HA is detected.
            if (n2_6585_flux > 0.0): # if [N II] is also detected.
                if (ha_6565_obs <= transition_wave): # if in muse.
                    ha_6550_6565_6585_flux   = (n2_6550_flux + ha_6565_flux + n2_6585_flux)
                    ha_6550_6565_6585_err    = (n2_6550_err + ha_6565_err + n2_6585_err)
                    ha_6550_6565_6585_ew_obs = (n2_6550_ew_obs + ha_6565_ew_obs + n2_6585_ew_obs)
                if (ha_6565_obs > transition_wave): # if in wfc3.
                    ha_6550_6565_6585_flux   = (n2_6550_flux + ha_6565_flux + n2_6585_flux) # could also hard-code as 1.133 * ha_6565_flux.
                    # [N II] / Ha = 0.133 is hard-coded for grism and errors are not calculated for 'tied' params.
                    ha_6550_6565_6585_err    = (1.133 * ha_6565_err)
                    ha_6550_6565_6585_ew_obs = (n2_6550_ew_obs + ha_6565_ew_obs + n2_6585_ew_obs)

            if (n2_6585_flux <= 0.0): # if [N II] is not detected.
                ha_6550_6565_6585_flux   = ha_6565_flux
                ha_6550_6565_6585_err    = ha_6565_err
                ha_6550_6565_6585_ew_obs = n2_6585_ew_obs
        else: # if the lines are not detected then return -1.
            ha_6550_6565_6585_flux   = ha_6565_flux
            ha_6550_6565_6585_err    = ha_6565_err
            ha_6550_6565_6585_ew_obs = ha_6565_ew_obs

        ############################################################################

        def calculate_doublet_line_flux(centroid_1, amplitude_index_1, rest_wave_1, centroid_2, amplitude_index_2, rest_wave_2):
            with np.errstate(invalid = 'ignore', divide = 'ignore'):
                if ((centroid_1 > np.min(lam_spec)) & (centroid_2 < np.max(lam_spec))):
                    if (centroid_2 < transition_wave):
                        sig     = ((out.params[fwhm_grism_idx] * out.params[fwhm_ratio_idx]) / 2.3548) # red fwhm * blue/red fwhm ratio.
                        sig_err = ((out.perror[fwhm_grism_idx] * out.params[fwhm_ratio_idx]) / 2.3548)
                    if (centroid_1 > transition_wave):
                        sig     = (out.params[fwhm_grism_idx] / 2.3548) # red fwhm
                        sig_err = (out.perror[fwhm_grism_idx] / 2.3548)
                    covar_term  = (2.0 * (out.covar[fwhm_grism_idx][amplitude_index_1] / (out.params[fwhm_grism_idx] * out.params[amplitude_index_1])))
                    # The emission line flux is 'line_flux' = (sqrt(2*pi) * (height) * (sigma)).
                    line_flux_1 = (np.sqrt(2.0 * math.pi) * (out.params[amplitude_index_1]) * (sig))
                    line_flux_2 = (line_flux_1 / out.params[amplitude_index_2])
                    line_flux   = (line_flux_1 + line_flux_2)
                    if (line_flux > 0.0):
                        line_err_1 = line_flux_1 * (np.sqrt(((out.perror[amplitude_index_1] / out.params[amplitude_index_1]) ** 2.0) + ((sig_err / sig) ** 2.0) + covar_term))
                        # line_err_2 = line_flux_2 * (np.sqrt(((out.perror[amplitude_index_2] / out.params[amplitude_index_2]) ** 2.0) + ((sig_err / sig) ** 2.0) + covar_term)) # Same fractional error?
                        line_err   = (line_err_1 * (1.0 + 1.0 / out.params[amplitude_index_2])) # check this.
                        line_err_2 = (line_err - line_err_1) # mpfit doesn't return the required parameter uncertainties for ratios so the error is the total minus the first line.
                    else:
                        w_1 = np.where((lam_spec > (centroid_1 - fwhm_guess)) & (lam_spec < (centroid_1 + fwhm_guess)))
                        line_err_1 = np.sqrt(np.sum(error_spec[w_1] ** 2.0))
                        w_2 = np.where((lam_spec > (centroid_2 - fwhm_guess)) & (lam_spec < (centroid_2 + fwhm_guess)))
                        line_err_2 = np.sqrt(np.sum(error_spec[w_2] ** 2.0))
                        w = np.where((lam_spec > (centroid_1 - fwhm_guess)) & (lam_spec < (centroid_2 + fwhm_guess)))
                        line_err = np.sqrt(np.sum(error_spec[w] ** 2.0))
                    line_cont_1   = emissionline_model(modelpars_nolines, rest_wave_1 * np.array([1.0 + out.params[z_idx]]))[0] # redshift.
                    line_ew_obs_1 = (line_flux_1 / line_cont_1)
                    line_cont_2   = emissionline_model(modelpars_nolines, rest_wave_2 * np.array([1.0 + out.params[z_idx]]))[0] # redshift.
                    line_ew_obs_2 = (line_flux_2 / line_cont_2)
                    line_cont     = emissionline_model(modelpars_nolines, ((rest_wave_1 + rest_wave_2) / 2.0) * np.array([1.0 + out.params[z_idx]]))[0] # redshift.
                    line_ew_obs   = (line_flux / line_cont)
                else:
                    line_flux_1   = (-1 / scl)
                    line_err_1    = (-1 / scl)
                    line_ew_obs_1 =  -1
                    line_flux_2   = (-1 / scl)
                    line_err_2    = (-1 / scl)
                    line_ew_obs_2 =  -1
                    line_flux     = (-1 / scl)
                    line_err      = (-1 / scl)
                    line_ew_obs   =  -1

            return line_flux, line_err, line_ew_obs, line_flux_1, line_err_1, line_ew_obs_1, line_flux_2, line_err_2, line_ew_obs_2

        ############################################################################

        '''
        Calculate the emission line fluxes and return them to measure_z_interactive().
        lyman alpha has a half-gaussian wing and the function returns gaussian params
        so divide the wing amplitude by two to divide the flux in half.
        '''
        out.params[la_wing_amp_idx] = (out.params[la_wing_amp_idx] / 2.0)

        la_1216_wing_flux, la_1216_wing_err, la_1216_wing_ew_obs, la_1216_flux, la_1216_err, la_1216_ew_obs, la_wing_flux, la_wing_err, la_wing_ew_obs = \
        calculate_doublet_line_flux(la_1216_obs, la_1216_idx, la_1216_vac, la_1216_obs, la_wing_amp_idx, la_1216_vac)

        n5_1238_1242_flux, n5_1238_1242_err, n5_1238_1242_ew_obs, n5_1238_flux, n5_1238_err, n5_1238_ew_obs, n5_1242_flux, n5_1242_err, n5_1242_ew_obs = \
        calculate_doublet_line_flux(n5_1238_obs, n5_1238_idx, n5_1238_vac, n5_1242_obs, n5_1242_idx, n5_1242_vac)

        c4_1548_1550_flux, c4_1548_1550_err, c4_1548_1550_ew_obs, c4_1548_flux, c4_1548_err, c4_1548_ew_obs, c4_1550_flux, c4_1550_err, c4_1550_ew_obs = \
        calculate_doublet_line_flux(c4_1548_obs, c4_1548_idx, c4_1548_vac, c4_1550_obs, c4_1550_idx, c4_1550_vac)

        o3_1660_1666_flux, o3_1660_1666_err, o3_1660_1666_ew_obs, o3_1660_flux, o3_1660_err, o3_1660_ew_obs, o3_1666_flux, o3_1666_err, o3_1666_ew_obs = \
        calculate_doublet_line_flux(o3_1660_obs, o3_1660_idx, o3_1660_vac, o3_1666_obs, o3_1666_idx, o3_1666_vac)

        s3_1883_1892_flux, s3_1883_1892_err, s3_1883_1892_ew_obs, s3_1883_flux, s3_1883_err, s3_1883_ew_obs, s3_1892_flux, s3_1892_err, s3_1892_ew_obs = \
        calculate_doublet_line_flux(s3_1883_obs, s3_1883_idx, s3_1883_vac, s3_1892_obs, s3_1892_idx, s3_1892_vac)

        c3_1907_1909_flux, c3_1907_1909_err, c3_1907_1909_ew_obs, c3_1907_flux, c3_1907_err, c3_1907_ew_obs, c3_1909_flux, c3_1909_err, c3_1909_ew_obs = \
        calculate_doublet_line_flux(c3_1907_obs, c3_1907_idx, c3_1907_vac, c3_1909_obs, c3_1909_idx, c3_1909_vac)

        m2_2796_2803_flux, m2_2796_2803_err, m2_2796_2803_ew_obs, m2_2796_flux, m2_2796_err, m2_2796_ew_obs, m2_2803_flux, m2_2803_err, m2_2803_ew_obs = \
        calculate_doublet_line_flux(m2_2796_obs, m2_2796_idx, m2_2796_vac, m2_2803_obs, m2_2803_idx, m2_2803_vac)

        o2_3727_3730_flux, o2_3727_3730_err, o2_3727_3730_ew_obs, o2_3727_flux, o2_3727_err, o2_3727_ew_obs, o2_3730_flux, o2_3730_err, o2_3730_ew_obs = \
        calculate_doublet_line_flux(o2_3727_obs, o2_3727_idx, o2_3727_vac, o2_3730_obs, o2_3730_idx, o2_3730_vac)

        o3_4959_5007_flux, o3_4959_5007_err, o3_4959_5007_ew_obs, o3_4959_flux, o3_4959_err, o3_4959_ew_obs, o3_5007_flux, o3_5007_err, o3_5007_ew_obs = \
        calculate_doublet_line_flux(o3_4959_obs, o3_4959_idx, o3_4959_vac, o3_5007_obs, o3_5007_idx, o3_5007_vac)

        o1_6300_6363_flux, o1_6300_6363_err, o1_6300_6363_ew_obs, o1_6300_flux, o1_6300_err, o1_6300_ew_obs, o1_6363_flux, o1_6363_err, o1_6363_ew_obs = \
        calculate_doublet_line_flux(o1_6300_obs, o1_6300_idx, o1_6300_vac, o1_6363_obs, o1_6363_idx, o1_6363_vac)

        s2_6716_6731_flux, s2_6716_6731_err, s2_6716_6731_ew_obs, s2_6716_flux, s2_6716_err, s2_6716_ew_obs, s2_6731_flux, s2_6731_err, s2_6731_ew_obs = \
        calculate_doublet_line_flux(s2_6716_obs, s2_6716_idx, s2_6716_vac, s2_6731_obs, s2_6731_idx, s2_6731_vac)

        s3_9069_9532_flux, s3_9069_9532_err, s3_9069_9532_ew_obs, s3_9069_flux, s3_9069_err, s3_9069_ew_obs, s3_9532_flux, s3_9532_err, s3_9532_ew_obs = \
        calculate_doublet_line_flux(s3_9069_obs, s3_9069_idx, s3_9069_vac, s3_9532_obs, s3_9532_idx, s3_9532_vac)

        ############################################################################
        ############################################################################

        fit_results = {}

        fit_results['fit_parameters'] = out.params
        fit_results['fit_status']     = out.status
        fit_results['chisq']          = chisq
        fit_results['scl_factor']     = scl
        fit_results['fwhm_muse']      = out.params[fwhm_grism_idx] * out.params[fwhm_ratio_idx]
        fit_results['fwhm_muse_error']= out.perror[fwhm_grism_idx] * out.perror[fwhm_ratio_idx] # correct or add in quadrature?
        fit_results['fwhm_g141']      = out.params[fwhm_grism_idx]
        fit_results['fwhm_g141_error']= out.perror[fwhm_grism_idx]
        fit_results['redshift']       = out.params[z_idx]
        fit_results['redshift_error'] = out.perror[z_idx]
        fit_results['la_1216_dz']     = out.params[la_1216_dzx]
        fit_results['c4_1548_dz']     = out.params[c4_1548_dzx]
        fit_results['uv_line_dz']     = out.params[uv_line_dzx]
        fit_results['m2_2796_dz']     = out.params[m2_2796_dzx]
        fit_results['o2_3727_dz']     = out.params[o2_3727_dzx]
        fit_results['o3_5007_dz']     = out.params[o3_5007_dzx]
        fit_results['s3_he_dz']       = out.params[s3_he_dzx]

        # calculate the flux, error, and equivalent width for each emission line.
        result_lines = [\
        'la_1216', 'la_wing', 'la_1216_wing', \
        'n5_1238', 'n5_1242', 'n5_1238_1242', \
        'c4_1548', 'c4_1550', 'c4_1548_1550', \
        'h2_1640', \
        'o3_1660', 'o3_1666', 'o3_1660_1666', \
        's3_1883', 's3_1892', 's3_1883_1892', \
        'c3_1907', 'c3_1909', 'c3_1907_1909', \
        'm2_2796', 'm2_2803', 'm2_2796_2803', \
        'o2_3727', 'o2_3730', 'o2_3727_3730', \
        'hg_4342', \
        'o3_4363', \
        'h2_4686', \
        'hb_4863', \
        'o3_4959', 'o3_5007', 'o3_4959_5007', \
        'o1_6300', 'o1_6363', 'o1_6300_6363', \
        'n2_6550', 'ha_6565', 'n2_6585', 'ha_6550_6565_6585', \
        's2_6716', 's2_6731', 's2_6716_6731', \
        's3_9069', 's3_9532', 's3_9069_9532', 'he10830']

        for line in result_lines:
            fit_results[line+'_flux']   = eval(line+'_flux') * scl
            fit_results[line+'_error']  = eval(line+'_err')  * scl
            fit_results[line+'_ew_obs'] = eval(line+'_ew_obs')

            # print the line fluxes and uncertainties to the terminal.
            # if (fit_results[line+'_flux'] > 0.0):
            #     perc_error = 100.0 * (fit_results[line+'_error'] / fit_results[line+'_flux'])
            #     if (len(line+'_flux') < 13):
            #         print line+'_flux =', '{:.2e}'.format(fit_results[line+'_flux'], 2), ' (+/- '+str(np.around(perc_error, 1))+' %)'

    else:
        fit_results = {}
        fit_results['fit_status'] = out.status

    # for j in range(len(out.params)):
    #     print 'out.params['+str(j)+']', out.params[j]

    return fit_results
