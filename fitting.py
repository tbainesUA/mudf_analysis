from wisp_analysis import *
from mpfit import *

def emissionline_model(pars, x):
    # The emission line model includes 'pars' for the continuum and emission line
    # parameters. The input 'x' is the wavelength array, which can include gaps.

    # mpfit.py allows for 'pegged', 'fixed', and 'limited' parameters with 'limits'.
    # See the mpfit.py code for a detailed description of how each is defined, but
    # it's useful to note that 'pegged' parameters are simply thos 'fixed' at their
    # upper or lower limit, which is confusing terminology.

    # pegged parameters:
    nnodes          = int(pars[0])   #### calculated from len(node_wave) and passed here.
    fit_region      = pars[1]  ##### could be possibly gotten rid of/hard coded.  ### not particularly useful now w/global cont.
    transition_wave = pars[2]
    z            = pars[3]
    oiii_5007_dz = pars[4]
    oii_3727_dz  = pars[5]
    siii_hei_dz  = pars[6]
    hg_4342_dz   = (oiii_5007_dz + oii_3727_dz) / 2.0
    ratio_fwhm   = pars[7]
    fwhm_red     = pars[8]
    fwhm_blue    = (ratio_fwhm * fwhm_red)

    sigma_blue = (fwhm_blue / 2.35)
    sigma_red  = (fwhm_red  / 2.35)

    ##### line amplitudes
    ### current sims include ha, hb, hg, oiii, oii, nii, sii
    ha_6565_amp   = pars[9]
    sii_6716_amp  = pars[10]
    sii_ratio     = pars[11]    ####  ~1.4 in the LDL.
    sii_6731_amp  = sii_6716_amp / sii_ratio
    nii_ha_frac   = pars[12]
    nii_6585_amp  = nii_ha_frac * ha_6565_amp
    nii_6550_amp  = nii_ha_frac * ha_6565_amp /3.
    oiii_5007_amp = pars[13]
    oiii_4959_amp = oiii_5007_amp / 3.
    hb_4863_amp   = pars[14]
    oii_3727_amp  = pars[15]
    oii_3730_amp  = oii_3727_amp * 1.4
    hg_4342_amp   = pars[16]
    siii_9069_amp = pars[17]
    siii_9532_amp =  2.48 * siii_9069_amp
    hei_10830_amp = pars[18]
    lya_1216_amp  = pars[19] # M.D.R 01/06/2021
    lya_1216_dz   = pars[20] # M.D.R 01/06/2021
    #print('pars = ', pars) # M.D.R 01/07/2021

    #if nnodes >= 2:
        ### y-values of spline nodes
    #cont_node_values = pars[19:19+nnodes] # M.D.R 01/06/2021
    cont_node_values = pars[21:21+nnodes] # M.D.R 01/06/2021
    #print('\ncont_node_values = ', cont_node_values) # M.D.R 01/07/2021

    #### pass the x-values (lam) for the spline nodes.  they'll be fixed.
    #clam = pars[19+nnodes:19 + 2* nnodes] # M.D.R 01/06/2021
    clam = pars[21+nnodes:21 + 2* nnodes] # M.D.R 01/06/2021


    #### define the observed wavelengths:
    #lya_1216_obs = 1215.67 * (1+z) # M.D.R 01/06/2021
    #lya_1216_obs = 1215.67 * (1+z+dz_lya) # M.D.R 01/06/2021  # MDR 2022/05/25 - removed and switched all to lya_1216_obs.
    #lya_obs = 1215.67 * (1+z) # M.D.R 01/06/2021
    lya_1216_obs  = 1215.670 * (1+z+lya_1216_dz) # M.D.R 01/06/2021
    oii_3727_obs  = 3727.092 * (1+z+oii_3727_dz) #3727.0 * (1+z+dz_oii) - M.D.R 01/06/2021
    oii_3730_obs  = 3729.875 * (1+z+oii_3727_dz) #3730.0 * (1+z+dz_oii) - M.D.R 01/06/2021
    hg_4342_obs   = 4341.684 * (1+z+hg_4342_dz) #4342.0 * (1+z + dz_hg) - M.D.R 01/06/2021
    hb_4863_obs   = 4862.683 * (1+z+oiii_5007_dz) #4863.0 * (1+z + dz_oiii) - M.D.R 01/06/2021
    oiii_4959_obs = 4960.295 * (1+z+oiii_5007_dz)# 4960.0 * (1+z+dz_oiii) - M.D.R 01/06/2021
    oiii_5007_obs = 5008.240 * (1+z+oiii_5007_dz) #5008.0 * (1+z+dz_oiii) - M.D.R 01/06/2021
    nii_6550_obs  = 6549.850 * (1+z) #6550.0 * (1+z) - M.D.R 01/06/2021
    ha_6565_obs   = 6564.610 * (1+z) #6565.0 * (1+z) - M.D.R 01/06/2021
    nii_6585_obs  = 6585.280 * (1+z) #6585.0 * (1+z) - M.D.R 01/06/2021
    sii_6716_obs  = 6718.290 * (1+z) #6718.0 * (1+z) - M.D.R 01/06/2021
    sii_6731_obs  = 6732.670 * (1+z) #6733.0 * (1+z) - M.D.R 01/06/2021
    siii_9069_obs = 9071.100 * (1+z+siii_hei_dz)   #### looked the vac waves up on nist
    siii_9532_obs = 9533.200 * (1+z+siii_hei_dz)
    hei_10830_obs = 10832.86 * (1+z+siii_hei_dz)

    #### if the siii lines have different resolution
    if ((siii_9532_obs > transition_wave) & (siii_9069_obs <= transition_wave)):
        siii_9532_amp = 2.48 * siii_9069_amp * sigma_blue/sigma_red

    cont_model = x * 0.
    line_model = x * 0

    ##### evaluate the continuum:  blue side
    #w1=np.where((x > oii_3727_obs - fit_region) & (x < transition_wave)) # M.D.R 01/06/2021
    # w1=np.where((x > lya_1216_obs - fit_region) & (x < transition_wave)) # M.D.R 01/06/2021
    w1=np.where((x > lya_1216_obs - fit_region) & (x < transition_wave)) # MDR 2022/05/25
    w2 = np.where(clam < transition_wave)
    clam_blue = clam[w2]
    node_values_blue = cont_node_values[w2]
    ### add junk nodes on if there are too few for a cubic spline.
    #### this gets around scipy's stupidity on cubic splines.
    if len(clam_blue) > 0:
        i=0
        while np.size(clam_blue) < 4:
            i= i +1
            clam_blue = np.append(clam_blue, transition_wave + 1000*i)
            node_values_blue = np.append(node_values_blue, node_values_blue[-1])

        cont_spline_rep_blue = interpolate.splrep(clam_blue, node_values_blue, s=0, k=3)
        if np.size(w1) > 0:
            cont_model[w1] = interpolate.splev(x[w1], cont_spline_rep_blue, der=0)

    ### evaluate the continuum, red side
    w1 = np.where((x >= transition_wave) & (x < hei_10830_obs  + (fit_region)))
    w2 = np.where(clam >= transition_wave)
    clam_red = clam[w2]
    node_values_red = cont_node_values[w2]

    ### add junk nodes on if there are too few for a cubic spline.
    #### this gets around scipy's stupidity on cubic splines.
    if len(clam_red) > 0:
        i=0
        while np.size(clam_red) < 4:
            i = i  +1
            clam_red = np.append(transition_wave - 1000 * i, clam_red)
            node_values_red = np.append(node_values_red[0], node_values_red)

        cont_spline_rep_red = interpolate.splrep(clam_red, node_values_red, s=0, k=3)
        if np.size(w1) > 0:
            cont_model[w1] = interpolate.splev(x[w1], cont_spline_rep_red, der=0)

    # The blocks of code below add each emission line to the line_model, which
    # consists of the continuum + emission lines. The lines are added manually
    # in each block which is unsustainable for adding many new lines. The below
    # function adds these lines one at a time, with a few lines that have special
    # constrains still added in individual blocks below.

    ############################################################################
    ############################################################################

    def add_emission_line_to_model(centroid, amplitude):

        if (centroid < transition_wave): sigma = sigma_blue
        if (centroid > transition_wave): sigma = sigma_red
        w = np.where((x > (centroid - fit_region)) & (x < (centroid + fit_region)))
        if (np.size(w) > 0): line_model[w] = line_model[w] + (amplitude * gaussian(x[w], centroid, sigma))

    ############################################################################

    add_emission_line_to_model(lya_1216_obs, lya_1216_amp)
    add_emission_line_to_model(oii_3727_obs, oii_3727_amp)
    add_emission_line_to_model(oii_3730_obs, oii_3730_amp)
    add_emission_line_to_model(hg_4342_obs, hg_4342_amp)
    add_emission_line_to_model(siii_9069_obs, siii_9069_amp)
    add_emission_line_to_model(siii_9532_obs, siii_9532_amp)
    add_emission_line_to_model(hei_10830_obs, hei_10830_amp)

    ############################################################################

    # The [O III] and HB lines are fit simultaneously in a single wavelength range.
    if (hb_4863_obs < transition_wave): sigma = sigma_blue
    if (hb_4863_obs > transition_wave): sigma = sigma_red
    # Incase hb has sigma_blue but oiii has sigma_red.
    if (oiii_5007_obs < transition_wave): sigma_oiii = sigma_blue
    if (oiii_5007_obs > transition_wave): sigma_oiii = sigma_red

    w = np.where((x > (hb_4863_obs - fit_region)) & (x < (oiii_5007_obs + fit_region)))

    if np.size(w) > 0:
        line_model[w] = line_model[w] + \
        oiii_5007_amp * gaussian(x[w], oiii_5007_obs, sigma_oiii) + \
        oiii_4959_amp * gaussian(x[w], oiii_4959_obs, sigma_oiii)+ \
        hb_4863_amp * gaussian(x[w], hb_4863_obs, sigma)

    ############################################################################

    # The Ha, [N II], and [S II] lines are fit simultaneously in a single wavelength range.
    if (ha_6565_obs < transition_wave): sigma = sigma_blue
    if (ha_6565_obs > transition_wave): sigma = sigma_red

    w = np.where((x > (nii_6550_obs - fit_region)) & (x < (sii_6731_obs + fit_region)))

    if np.size(w) > 0:
        line_model[w] = line_model[w] + \
        ha_6565_amp  * gaussian(x[w], ha_6565_obs,  sigma) + \
        nii_6550_amp * gaussian(x[w], nii_6550_obs, sigma) + \
        nii_6585_amp * gaussian(x[w], nii_6585_obs, sigma) + \
        sii_6716_amp * gaussian(x[w], sii_6716_obs, sigma) + \
        sii_6731_amp * gaussian(x[w], sii_6731_obs, sigma)

    ############################################################################

    model = cont_model + line_model

    return model

    ############################################################################
    ############################################################################
    #
    # #### start adding lines:
    # ### M.D.R 01/06/2021 ###
    # if lya_1216_obs > transition_wave:
    #    sigma = sigma_red
    # else:
    #     sigma = sigma_blue
    # w=np.where( (x > lya_1216_obs - fit_region) & (x < lya_1216_obs + fit_region))
    # if np.size(w) > 0:
    #     line_lya = lya_1216_amp  * gaussian(x[w], lya_1216_obs, sigma)
    #     line_model[w] = line_model[w] + line_lya
    # ### M.D.R 01/06/2021 ###
    #
    #
    #   #### region oii
    # if oii_3727_obs > transition_wave:
    #    sigma = sigma_red
    # else:
    #     sigma = sigma_blue
    # w=np.where((x > oii_3727_obs - fit_region) & (x < oii_3727_obs + fit_region ))
    # if np.size(w) > 0:
    #     line_oii = oii_3727_amp * gaussian(x[w], oii_3727_obs, sigma) + \
    #                oii_3730_amp * gaussian(x[w], oii_3730_obs, sigma)
    #     line_model[w] = line_model[w] + line_oii
    #
    #
    # if hg_4342_obs > transition_wave:
    #    sigma = sigma_red
    # else:
    #     sigma = sigma_blue
    # w=np.where( (x > hg_4342_obs - fit_region) & (x < hg_4342_obs + fit_region))
    # if np.size(w) > 0:
    #     line_hg = hg_4342_amp  * gaussian(x[w], hg_4342_obs, sigma)
    #     line_model[w] = line_model[w] + line_hg
    #
    #
    # ### region around oiii/hb
    # if hb_4863_obs > transition_wave:
    #    sigma = sigma_red
    # else:
    #     sigma = sigma_blue
    #
    #  ### in case hb has sigma_blue but oiii has sigma_red.
    # if oiii_5007_obs >  transition_wave:
    #     sigma_oiii = sigma_red
    # else:
    #     sigma_oiii = sigma_blue
    # w=np.where( (x > hb_4863_obs - fit_region)  & (x < oiii_5007_obs + fit_region))
    # if np.size(w) > 0:
    #
    #     lines_oiii = oiii_5007_amp * gaussian(x[w], oiii_5007_obs, sigma_oiii) + \
    #     oiii_5007_amp/3. * gaussian(x[w], oiii_4959_obs, sigma_oiii)+ \
    #     hb_4863_amp * gaussian(x[w], hb_4863_obs, sigma)
    #
    #     line_model[w] = line_model[w] + lines_oiii
    #
    # if ha_6565_obs > transition_wave:
    #    sigma = sigma_red
    # else:
    #     sigma = sigma_blue
    # if sii_6716_obs >  transition_wave:
    #     sigma_sii = sigma_red
    # else:
    #     sigma_sii = sigma_blue
    # w=np.where( (x >= ha_6565_obs - fit_region) & (x <= sii_6731_obs + fit_region))
    # if np.size(w) > 0:
    #     lines_ha = ha_6565_amp * gaussian(x[w], ha_6565_obs, sigma) + \
    #          nii_6550_amp * gaussian(x[w], nii_6550_obs, sigma) + \
    #          nii_6585_amp * gaussian(x[w], nii_6585_obs, sigma) + \
    #          sii_6716_amp * gaussian(x[w], sii_6716_obs, sigma_sii) + \
    #          sii_6731_amp * gaussian(x[w], sii_6731_obs, sigma_sii)
    #     line_model[w] = line_model[w] + lines_ha
    #
 #
 #
 #    ##### SIII 9069
 #    if siii_9069_obs > transition_wave:
 #       sigma = sigma_red
 #    else:
 #        sigma = sigma_blue
 #    w=np.where( (x > siii_9069_obs - fit_region) & (x < siii_9069_obs + fit_region))
 #    if np.size(w) > 0:
 #        line_siii_9069 = siii_9069_amp  * gaussian(x[w], siii_9069_obs, sigma)
 #        line_model[w] = line_model[w] +  line_siii_9069
 #
 #
 # ##### SIII 9532
 #    if siii_9532_obs > transition_wave:
 #       sigma = sigma_red
 #    else:
 #        sigma = sigma_blue
 #    w=np.where( (x > siii_9532_obs - fit_region) & (x < siii_9532_obs + fit_region))
 #    if np.size(w) > 0:
 #        line_siii_9532 = siii_9532_amp  * gaussian(x[w], siii_9532_obs, sigma)
 #        line_model[w] = line_model[w] + line_siii_9532
 #
 #
 # ##### hei_10830
 #    if hei_10830_obs > transition_wave:
 #       sigma = sigma_red
 #    else:
 #        sigma = sigma_blue
 #    w=np.where( (x > hei_10830_obs - fit_region) & (x < hei_10830_obs + fit_region))
 #    if np.size(w) > 0:
 #        line_hei_10830 = hei_10830_amp  * gaussian(x[w], hei_10830_obs, sigma)
 #        line_model[w] = line_model[w] + line_hei_10830

# ha moved from here.
# oiii moved from here.


def model_resid(pars, fjac=None, lam = None, flux = None, err = None):
    model = emissionline_model(pars, lam)
    resid = (flux - model) / err
    status = 0
    return [status, resid]


def fit_obj(input_list):
    lam_spec = input_list[0]
    flux_spec = input_list[1]
    error_spec = input_list[2]
    config_pars = input_list[3]
    z_in = input_list[4]   ### this is the best-guess redshift, or the current line-list redshift
    fwhm_guess = input_list[5]    #### red fwhm guess
    beam_name = input_list[6] ### not used here.

    line_mask = config_pars['line_mask']
    fit_region= config_pars['fit_region']
    disp_red  = config_pars['dispersion_red']   ### the dispersion in setting a window for finding the peak line pixels.
    disp_blue = config_pars['dispersion_blue']
    node_wave = config_pars['node_wave']
    nnodes = len(node_wave)
    transition_wave = config_pars['transition_wave']

    #### set dz as the range of wavelengths over which we limit z, can change.
    dz = config_pars['delta_z']  ### this is the maximum shift from the line list redshift/best estimate.
    scl = 1e-18    ### the fitter works best when numbers are not too large/small.
    ### guess the wavelengths of the lines.
    ### in real data, this will come from a peak pixel estimate of the line center
    ### which we'll get from the process of finding the line.

    #### these only need to be approximate for masking the lines.
    # I'm overly cautious and so updated these to exact values. # MDR 2022/05/25
    lya_1216_obs  = 1215.670 * (1+z_in) # M.D.R 01/06/2021
    oii_3727_obs  = 3727.092 * (1+z_in)
    hg_4342_obs   = 4341.684 * (1+z_in)
    hb_4863_obs   = 4862.683 * (1+z_in)
    oiii_5007_obs = 5008.240 * (1+z_in)
    ha_6565_obs   = 6564.610 * (1+z_in)
    sii_obs       = 6725.480 * (1+z_in)
    sii_6716_obs  = 6718.290 * (1+z_in)
    sii_6731_obs  = 6732.670 * (1+z_in)
    siii_9069_obs = 9071.100 * (1+z_in)
    siii_9532_obs = 9533.200 * (1+z_in)
    hei_10830_obs = 10832.86 * (1+z_in)

    flux_spec = flux_spec / scl
    error_spec = error_spec/scl


    ############################################################
    #### step 1, fit only the continuum, ignoring the regions nearest the lines.
    ### the goal here is to get the best guesses on the continuum to use when fitting the whole thing.
    ### I am no longer sure of the importance of this step; I should test working without it.
    #############################################################

    #### find the windows that we estimate to cover the continuum adjacent to the lines.
    ### this doesn't have to be perfect, because we're only using it to get gueses.
  #  w=np.where( ((lam_spec > ha_obs -fit_region) & (lam_spec < ha_obs - line_mask)) |
  #              ((lam_spec > sii_obs + line_mask) & (lam_spec < sii_obs  + fit_region)) |
  #         ((lam_spec > hb_obs - fit_region) & (lam_spec < hb_obs - line_mask))
  #        | ((lam_spec > oiii_obs + line_mask) &  (lam_spec < oiii_obs + fit_region))
  #        | ((lam_spec > hg_obs - fit_region) & (lam_spec < hg_obs - line_mask))
  #        | ((lam_spec > hg_obs + line_mask) & (lam_spec < hg_obs + fit_region))
  #        | ((lam_spec > oii_obs - fit_region) & (lam_spec < oii_obs - line_mask))
  #        | ((lam_spec > oii_obs + line_mask) & (lam_spec < oii_obs + fit_region))
  #        | ((lam_spec > he1_obs - fit_region) & (lam_spec < he1_obs - line_mask))
  #        | ((lam_spec > he1_obs + line_mask) & (lam_spec < he1_obs + fit_region))
  #        | ((lam_spec > siii_9069_obs - fit_region) & (lam_spec < siii_9069_obs - line_mask))
  #        | ((lam_spec > siii_9069_obs + line_mask) & (lam_spec < siii_9069_obs + fit_region))
  #        | ((lam_spec > siii_9532_obs - fit_region) & (lam_spec < siii_9532_obs - line_mask))
  #        | ((lam_spec > siii_9532_obs + line_mask) & (lam_spec < siii_9532_obs + fit_region)))

    ### low/high regions to ignore for sure.
    mask_spec = lam_spec * 0.
    w=np.where((lam_spec < lya_1216_obs - fit_region) | (lam_spec > hei_10830_obs + fit_region)) # M.D.R 01/06/2021
#    w=np.where((lam_spec < oii_obs - fit_region) | (lam_spec > he1_obs + fit_region)) # M.D.R 01/06/2021 - original code
    mask_spec[w] = 1.
    w=np.where( (lam_spec > ha_6565_obs - line_mask) & (lam_spec < sii_obs + line_mask ))
    mask_spec[w] = 1.
    w=np.where( (lam_spec > hb_4863_obs - line_mask) & (lam_spec < oiii_5007_obs + line_mask))
    mask_spec[w] = 1.
    w=np.where( (lam_spec > hg_4342_obs - line_mask) & (lam_spec < hg_4342_obs+ line_mask))
    mask_spec[w] = 1.
    w=np.where( (lam_spec > oii_3727_obs - line_mask) & (lam_spec < oii_3727_obs + line_mask))
    mask_spec[w] = 1.
    w=np.where( (lam_spec > hei_10830_obs- line_mask) & (lam_spec < hei_10830_obs + line_mask))
    mask_spec[w] = 1.
    w=np.where( (lam_spec > siii_9069_obs - line_mask) & (lam_spec < siii_9069_obs + line_mask))
    mask_spec[w] = 1.
    w=np.where( (lam_spec > siii_9532_obs - line_mask) & (lam_spec < siii_9532_obs + line_mask))
    mask_spec[w] = 1.

    w=np.where(mask_spec == 0.)

    cont_guess = np.median(flux_spec[w])
    ####### define model constraints
    # z_in, fwhm, nnodes, line amps, spline y-guesses
    #pguess = np.zeros(19 + 2 * nnodes) # M.D.R 01/06/2021
    pguess = np.zeros(21 + 2 * nnodes) # M.D.R 01/06/2021

    ### pegged parameters
    pguess[0] = nnodes
    pguess[1] = fit_region
    pguess[2] = transition_wave
    pguess[3] = z_in
    pguess[4] = 0   ### oiii_5007_dz
    pguess[5] = 0   ### oii_3727_dz
    pguess[6] = 0   ### siii_hei_dz
    pguess[7] = 0.027    ### blue/red fwhm ratio. # pguess[7] = 0.5    ### blue/red fwhm ratio. - M.D.R. - 10/15/2020
    pguess[8] = fwhm_guess    ### fwhm_red.
    pguess[11] = 1.4  ### sii ratio.  ### since we divide by this number it shouldn't be zero.

    pguess[20] = 0   ### lya_1216_dz # M.D.R 01/06/2021
    #pguess[19:19 +nnodes] = cont_guess # M.D.R 01/06/2021
    pguess[21:21 +nnodes] = cont_guess # M.D.R 01/06/2021
    #pguess[19+nnodes: 19 + 2 * nnodes] = node_wave # M.D.R 01/06/2021
    pguess[21+nnodes: 21 + 2 * nnodes] = node_wave # M.D.R 01/06/2021

    ### the rest of the pguess are amplitudes and they are zeros.

    ### define the parinfo and fill in the values.
    npars= len(pguess)
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
              for i in range(npars)]
    for i in range(npars): parinfo[i]['value'] = pguess[i]

    ### apply limits and fix pegged parameters
    parinfo[0]['fixed'] = 1   # nnodes
    parinfo[1]['fixed'] = 1   # fit_region
    parinfo[2]['fixed'] = 1   # transition wave
    parinfo[3]['fixed'] = 1   # z  ### fix the z for the first iteration
    parinfo[4]['fixed'] = 1   # dz   These don't matter for the continuum only fit.
    parinfo[5]['fixed'] = 1   # dz
    parinfo[6]['fixed'] = 1   # dz
    parinfo[7]['fixed'] = 1  ### ratio of blue to red fwhm. pegged at 0.5

    #### parameters fixed for the continuum fitting
    parinfo[6]['fixed'] = 1  # z

    ### should be fixed to zero because of how parinfo was defined.
    parinfo[9]['fixed'] = 1    #ha
    parinfo[10]['fixed'] = 1   #sii
    parinfo[11]['fixed'] = 1   #sii_ratio
    parinfo[12]['fixed'] = 1   #nii_frac
    parinfo[13]['fixed'] = 1   #oiii
    parinfo[14]['fixed'] = 1   #hb
    parinfo[15]['fixed'] = 1   #oii
    parinfo[16]['fixed'] = 1   #hg
    parinfo[17]['fixed'] = 1   # siii_9069
    parinfo[18]['fixed'] = 1   # he1 10830
    parinfo[19]['fixed'] = 1   # lya 1216 # M.D.R 01/06/2021
    parinfo[20]['fixed'] = 1   # lya_1216_dz # M.D.R 01/06/2021

    ### also fix the node x-values.
    #for i in range(19+nnodes, 19+2*nnodes): parinfo[i]['fixed'] = 1  # M.D.R 01/06/2021
    for i in range(21+nnodes, 21+2*nnodes): parinfo[i]['fixed'] = 1  # M.D.R 01/06/2021

    ### set up the arguments that go to the fitters
    fa = {'lam':lam_spec[w], 'flux':flux_spec[w], 'err':error_spec[w]}

    out = mpfit(model_resid, pguess, functkw=fa, parinfo = parinfo, quiet=True)
    #### if out.status is zero then the spectrum was junk.  do not go on.
    if out.status > 0:
        #fit = emissionline_model(out.params, lam_spec[w])

        ############################################################
        ###### step 2, repeat the fit, now using the continuum guesses
        #############################################################
        #### first select the new wavelength windows for fitting:
        #w=np.where((lam_spec < oii_obs - fit_region) | (lam_spec > he1_10830 + fit_region))
        ### evaulate continuum fit again, this time indluding the wavlengths under the line:
        #fit = emissionline_model(out.params, lam_spec[w])

        # overplot the continuum only fit from step 1, this time evaluting it under the emission line as well.
        #plt.plot(lam_spec[w], fit, color = 'green')

         #####################  2a ############ improve peak guesses
         #### adjust redshift guess according to peak of lines.
         #### simplified!  keeping zguess instead of re-guessing.


        # The blocks of code below determine the best guess amplitude for each
        # emission line. The amplitudes are added manually to the model, which
        # is unsustainable for adding many new lines. The below function
        # determines the amplitudes one at a time and adds them to the input.

        ############################################################################
        ############################################################################

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
            #print('Running estimate_emission_line_amplitudes()...')

            return amplitude_guess

        ############################################################################

        lya_1216_amp_guess  = estimate_emission_line_amplitudes(lya_1216_obs)
        oii_3727_amp_guess  = estimate_emission_line_amplitudes(oii_3727_obs)
        hg_4342_amp_guess   = estimate_emission_line_amplitudes(hg_4342_obs)
        hb_4863_amp_guess   = estimate_emission_line_amplitudes(hb_4863_obs)
        oiii_5007_amp_guess = estimate_emission_line_amplitudes(oiii_5007_obs)
        ha_6565_amp_guess   = estimate_emission_line_amplitudes(ha_6565_obs)
        sii_6716_amp_guess  = estimate_emission_line_amplitudes(sii_6716_obs)
        sii_6731_amp_guess  = estimate_emission_line_amplitudes(sii_6731_obs)
        siii_9069_amp_guess = estimate_emission_line_amplitudes(siii_9069_obs)
        hei_10830_amp_guess = estimate_emission_line_amplitudes(hei_10830_obs)

        ############################################################################
        ############################################################################
        #
        #
        # ###### calculate guesses for amplitudes and redshifts of strongest lines.
        # if ((ha_6565_obs > np.min(lam_spec)) & (ha_6565_obs < np.max(lam_spec) ) ) :
        #     if ha_6565_obs > transition_wave:
        #        disp    = disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > ha_6565_obs - 5* disp) & (lam_spec < ha_6565_obs + 5*disp))
        #     if np.size(wline) >0:
        #         peakval = np.max(flux_spec[wline])
        #         ha_6565_amp_guess = peakval - emissionline_model(out.params, np.array([ha_6565_obs]))  ### flam - flam_cont
        #         ha_6565_amp_guess = ha_6565_amp_guess[0] ### ensure it is scalar
        #         if ha_6565_amp_guess < 0:
        #             ha_6565_amp_guess = 0.
        #     else:
        #         ha_6565_amp_guess =1.
        # else:
        #     ha_6565_amp_guess = 1.    ### doesn't matter not being fit.
        #
        #
        # if ((oiii_5007_obs > np.min(lam_spec) ) & (oiii_5007_obs < np.max(lam_spec)) ) :
        #      if oiii_5007_obs >  transition_wave:
        #          disp = disp_red
        #      else:
        #          disp = disp_blue
        #      wline = np.where((lam_spec > oiii_5007_obs - 5 * disp) & (lam_spec < oiii_5007_obs + 5*disp))
        #      ### a gap in the spectrum can still break the code.
        #      if np.size(wline) > 0:
        #          peakval = np.max(flux_spec[wline])
        #          oiii_5007_amp_guess=  peakval - emissionline_model(out.params, np.array([oiii_5007_obs])) ### flam - flam_cont
        #          oiii_5007_amp_guess = oiii_5007_amp_guess[0]
        #          if oiii_5007_amp_guess < 0:
        #              oiii_5007_amp_guess = 0.
        #      else:
        #          oiii_5007_amp_guess = 1.
        # else:
        #      oiii_5007_amp_guess = 1.
        #
        #
        # if ((hg_4342_obs > np.min(lam_spec)) & ( hg_4342_obs < np.max(lam_spec))):
        #     if hg_4342_obs > transition_wave:
        #         disp = disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > hg_4342_obs - 5*disp) & (lam_spec < hg_4342_obs + 5*disp))
        #     if np.size(wline) > 0 :
        #         peakval = np.max(flux_spec[wline])
        #         hg_4342_amp_guess = peakval - emissionline_model(out.params, np.array([hg_4342_obs]))  ### flam - flam_cont
        #         hg_4342_amp_guess = hg_4342_amp_guess[0]
        #         if hg_4342_amp_guess < 0:
        #             hg_4342_amp_guess = 0.
        #     else:
        #         hg_4342_amp_guess = 1
        # else:
        #     hg_4342_amp_guess = 1
        #
        #
        # if ((hb_4863_obs > np.min(lam_spec)) & ( hb_4863_obs < np.max(lam_spec))):
        #     if hb_4863_obs > transition_wave:
        #         disp = disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > hb_4863_obs - 5*disp) & (lam_spec < hb_4863_obs + 5*disp))
        #     if np.size(wline) > 0:
        #         peakval = np.max(flux_spec[wline])
        #         hb_4863_amp_guess = peakval - emissionline_model(out.params, np.array([hb_4863_obs]))  ### flam - flam_cont
        #         hb_4863_amp_guess = hb_4863_amp_guess[0]
        #         if hb_4863_amp_guess < 0:
        #             hb_4863_amp_guess = 0.
        #     else:
        #         hb_4863_amp_guess=1.
        # else:
        #     hb_4863_amp_guess = 1
        #
        #
        # if ((oii_3727_obs > np.min(lam_spec)) & ( oii_3727_obs < np.max(lam_spec) ) ):
        #     if oii_3727_obs > transition_wave:
        #         disp =  disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > oii_3727_obs - 5*disp) & (lam_spec < oii_3727_obs + 5*disp))
        #     if np.size(wline) > 0 :
        #         peakval = np.max(flux_spec[wline])
        #         oii_3727_amp_guess = peakval - emissionline_model(out.params, np.array([oii_3727_obs]))  ### flam - flam_cont
        #         oii_3727_amp_guess = oii_3727_amp_guess[0]
        #         if oii_3727_amp_guess < 0:
        #             oii_3727_amp_guess = 0.
        #     else:
        #         oii_3727_amp_guess = 1.
        # else:
        #     oii_3727_amp_guess = 1
        #
        #
        # if ((sii_6716_obs > np.min(lam_spec)) & ( sii_6716_obs < np.max(lam_spec) ) ):
        #     if sii_6716_obs > transition_wave:
        #         disp =  disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > sii_6716_obs - 5*disp) & (lam_spec < sii_6716_obs + 5*disp))
        #     if np.size(wline) > 0:
        #         peakval = np.max(flux_spec[wline])
        #         sii_6716_amp_guess = peakval - emissionline_model(out.params, np.array([sii_6716_obs]))  ### flam - flam_cont
        #         sii_6716_amp_guess = sii_6716_amp_guess[0]
        #         if sii_6716_amp_guess < 0:
        #             sii_6716_amp_guess = 0.
        #     else:
        #         sii_6716_amp_guess = 1
        # else:
        #     sii_6716_amp_guess = 1
        #
        #
        # if ((siii_9069_obs > np.min(lam_spec)) & ( siii_9069_obs < np.max(lam_spec) ) ):
        #     if siii_9069_obs > transition_wave:
        #         disp =  disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > siii_9069_obs - 5*disp) & (lam_spec < siii_9069_obs + 5*disp))
        #     if np.size(wline) > 0:
        #         peakval = np.max(flux_spec[wline])
        #         siii_9069_amp_guess = peakval - emissionline_model(out.params, np.array([siii_9069_obs]))  ### flam - flam_cont
        #         siii_9069_amp_guess = siii_9069_amp_guess[0]
        #         if siii_9069_amp_guess < 0:
        #             siii_9069_amp_guess = 0.
        #     else:
        #         siii_9069_amp_guess = 1.
        # else:
        #     siii_9069_amp_guess = 1
        #
        #
        # if ((hei_10830_obs > np.min(lam_spec)) & ( hei_10830_obs < np.max(lam_spec) ) ):
        #     if hei_10830_obs > transition_wave:
        #         disp =  disp_red
        #     else:
        #         disp = disp_blue
        #     wline = np.where((lam_spec > hei_10830_obs - 5*disp) & (lam_spec < hei_10830_obs + 5*disp))
        #     if np.size(wline) >0:
        #         peakval = np.max(flux_spec[wline])
        #         hei_10830_amp_guess = peakval - emissionline_model(out.params, np.array([hei_10830_obs]))  ### flam - flam_cont
        #         hei_10830_amp_guess = hei_10830_amp_guess[0]
        #         if hei_10830_amp_guess < 0:
        #             hei_10830_amp_guess = 0.
        #     else:
        #         hei_10830_amp_guess = 1
        # else:
        #     hei_10830_amp_guess = 1
        #
# ### M.D.R 01/06/2021 ###
#         if ((lya_1216_obs > np.min(lam_spec)) & ( lya_1216_obs < np.max(lam_spec))):
#             if lya_1216_obs > transition_wave:
#                 disp = disp_red
#             else:
#                 disp = disp_blue
#             wline = np.where((lam_spec > lya_1216_obs - 5*disp) & (lam_spec < lya_1216_obs + 5*disp))
#             if np.size(wline) > 0 :
#                 peakval = np.max(flux_spec[wline])
#                 lya_1216_amp_guess = peakval - emissionline_model(out.params, np.array([lya_1216_obs]))  ### flam - flam_cont
#                 lya_1216_amp_guess = lya_1216_amp_guess[0]
#                 if lya_1216_amp_guess < 0:
#                     lya_1216_amp_guess = 0.
#             else:
#                 lya_1216_amp_guess = 1
#         else:
#             lya_1216_amp_guess = 1
# ### M.D.R 01/06/2021 ###



        ######### 2b fill in the parinfo array with this information:
        #pguess2 = [redshift_guess, fwhm_guess, nnodes,  ha_6565_amp_guess, 0.1 * ha_6565_amp_guess, 0.1 * ha_6565_amp_guess, 0.1 * ha_6565_amp_guess,
        #          oiii_5007_amp_guess, 0.3 * oiii_5007_amp_guess, hg_4342_amp_guess, oii_3727_amp_guess, out.params[11:11+nnodes*4], fit_region, line_mask ]
#
        parinfo2 = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
                  for i in range(npars)]

        #for i in range(npars): parinfo2[i]['value'] = pguess2[i]

        #### define and fix pegged parameters
        parinfo2[0]['value'] = nnodes
        parinfo2[1]['value']  = fit_region
        parinfo2[2]['value'] = transition_wave
        parinfo2[3]['value'] = z_in
        parinfo2[4]['value'] = 0 ### dz
        parinfo2[5]['value'] = 0 ### dz
        parinfo2[6]['value'] = 0  ### dz

        parinfo2[0]['fixed'] = 1   # nnodes
        parinfo2[1]['fixed'] = 1   # fit_region
        parinfo2[2]['fixed'] = 1   # transition_wave
        parinfo2[3]['fixed'] = 0   # z_in
        if z_in > 1.55:
            parinfo2[4]['fixed'] = 1   # fix oiii_5007_dz
        else:
            parinfo2[4]['fixed'] = 0   # vary oiii_5007_dz
        if z_in > 2.4:
            parinfo2[5]['fixed'] = 0   # fix oii_3727_dz  # M.D.R 01/06/2021 - was 1, try 0?
        else:
            parinfo2[5]['fixed'] = 0   # vary oii_3727_dz
        parinfo2[6]['fixed'] = 1   # always fix dz_siii

        parinfo2[7]['fixed'] = 1   ### blue to red fwhm ratio.  fixed at 0.5
        parinfo2[7]['value'] = 0.027 # parinfo2[7]['value'] = 0.5 # - M.D.R. - 10/15/2020
        parinfo2[8]['value'] = fwhm_guess   ### this guess is either the red fwhm or 2x the blue fwhm.
        parinfo2[9]['value'] = ha_6565_amp_guess
        parinfo2[10]['value'] = sii_6716_amp_guess
        parinfo2[11]['value'] = 1.4 ### will deal with this carefully below
        parinfo2[12]['value'] = 0.1  ### this parameter is a fraction of ha amplitude, unlike others.   ## de
        parinfo2[13]['value'] = oiii_5007_amp_guess
        parinfo2[14]['value'] = hb_4863_amp_guess
        parinfo2[15]['value'] = oii_3727_amp_guess
        parinfo2[16]['value'] = hg_4342_amp_guess
        parinfo2[17]['value'] = siii_9069_amp_guess
        parinfo2[18]['value'] = hei_10830_amp_guess
        parinfo2[19]['value'] = lya_1216_amp_guess # M.D.R 01/06/2021
        parinfo2[20]['fixed'] = 0   # vary lya_1216_dz # M.D.R 01/06/2021

        #for i in range(19, 19 + 2 * nnodes):  parinfo2[i]['value'] = out.params[i]  #### set x and y values of nodes. # M.D.R 01/06/2021
        for i in range(21, 21 + 2 * nnodes):  parinfo2[i]['value'] = out.params[i]  #### set x and y values of nodes. # M.D.R 01/06/2021
        #for i in range(19+nnodes, 19+2*nnodes): parinfo2[i]['fixed'] = 1   ### fix x values of nodes. # M.D.R 01/06/2021
        for i in range(21+nnodes, 21+2*nnodes): parinfo2[i]['fixed'] = 1   ### fix x values of nodes. # M.D.R 01/06/2021

        #z
        parinfo2[3]['limited'][0] = 1
        parinfo2[3]['limited'][1] = 1
        parinfo2[3]['limits'][0]  = z_in - dz
        parinfo2[3]['limits'][1]  = z_in + dz

        #oiii_5007_dz
        parinfo2[4]['limited'][0] = 1
        parinfo2[4]['limited'][1] = 1
        parinfo2[4]['limits'][0]  = -1*dz
        parinfo2[4]['limits'][1]  = dz

        ## oii_3727_dz
        parinfo2[5]['limited'][0] = 1
        parinfo2[5]['limited'][1] = 1
        parinfo2[5]['limits'][0]  = -1*dz
        parinfo2[5]['limits'][1]  = dz

        #dz_siii/he1
        parinfo2[6]['limited'][0] = 1
        parinfo2[6]['limited'][1] = 1
        parinfo2[6]['limits'][0]  = -1* dz
        parinfo2[6]['limits'][1]  = dz

        #fwhm
        parinfo2[8]['limited'][0] = 1
        parinfo2[8]['limits'][0]  = config_pars['min_fwhm_scl'] * fwhm_guess
        parinfo2[8]['limited'][1] = 1
        parinfo2[8]['limits'][1]  = config_pars['max_fwhm_scl'] * fwhm_guess

       ### amplitudes have to be positive

    #ha
        parinfo2[9]['limited'][0] = 1
        parinfo2[9]['limits'][0]  = 0

    #sii
        parinfo2[10]['limited'][0] = 1
        parinfo2[10]['limits'][0]  = 0


    #### these are no longer features in the configuration file, for wisp at least, but
    ### but they can be changed by a quick code edit

    # sii   & nii ratios are fixed to the values given above.  do not even think you should change this for grism data.
        parinfo2[11]['fixed'] = 1
        parinfo2[12]['fixed'] = 1

    #  oiii
        parinfo2[13]['limited'][0] = 1
        parinfo2[13]['limits'][0]  = 0

    # hb
        parinfo2[14]['limited'][0] = 1
        parinfo2[14]['limits'][0]  = 0

    # oii
        parinfo2[15]['limited'][0] = 1
        parinfo2[15]['limits'][0]  = 0

    # hg
        parinfo2[16]['limited'][0] = 1
        parinfo2[16]['limits'][0]  = 0

    # siii 9069
        parinfo2[17]['limited'][0] = 1
        parinfo2[17]['limits'][0]  = 0

    # he 1
        parinfo2[18]['limited'][0] = 1
        parinfo2[18]['limits'][0]  = 0
### M.D.R 01/06/2021 ###
    # lya
        parinfo2[19]['limited'][0] = 1
        parinfo2[19]['limits'][0]  = 0

    ## lya_1216_dz
        parinfo2[20]['limited'][0] = 1
        parinfo2[20]['limited'][1] = 1
        parinfo2[20]['limits'][0]  = -1*dz
        parinfo2[20]['limits'][1]  = dz
### M.D.R 01/06/2021 ###


        # step 2c, do the fit:
        #w=np.where((lam_spec > oii_obs - fit_region) & (lam_spec < he1_obs + fit_region)) # M.D.R 01/07/2021
        w=np.where((lam_spec > lya_1216_obs - fit_region) & (lam_spec < hei_10830_obs + fit_region)) # M.D.R 01/07/2021
        fa2 = {'lam':lam_spec[w], 'flux':flux_spec[w], 'err':error_spec[w]}

        pguess2 = []
        for i in range(npars): pguess2.append(parinfo2[i]['value'])
        #print parinfo2

        out = mpfit(model_resid, pguess2, functkw=fa2, parinfo = parinfo2, quiet=True)
        chisq  = out.fnorm / (len(w[0]) - npars)

    ### evaluate continuum spline.
        modelpars_nolines = cp.deepcopy(out.params)
        #modelpars_nolines[9:19] = 0. # M.D.R 01/06/2021
        modelpars_nolines[9:21] = 0. # M.D.R 01/06/2021
        modelpars_nolines[11] = 1 ### sii ratio-- division so do this to avoid divide by zero.

        #####re-evaluate observed wavelengths based on the refined fit to z;
        z_out         = out.params[3]
        # lya_obs       = 1216 * (1+z_out) # M.D.R 01/06/2021
        # oii_obs       = 3727 * (1+z_out)
        # hg_obs        = 4342 * (1+z_out)
        # hb_obs        = 4863 * (1+z_out)
        # oiii_obs      = 5008 * (1+z_out)
        # ha_obs        = 6565 * (1+z_out)
        # sii_obs       = 6725 * (1+z_out)
        # sii_6716_obs  = 6718 * (1+z_out)
        # sii_6731_obs  = 6733 * (1+z_out)
        # siii_9069_obs = 9071 * (1+z_out)
        # siii_9532_obs = 9533 * (1+z_out)
        # he1_obs       = 10833* (1+z_out)
        # I'm overly cautious and so updated these to exact values. # MDR 2022/05/25
        lya_1216_obs  = 1215.670 * (1+z_out) # M.D.R 01/06/2021
        oii_3727_obs  = 3727.092 * (1+z_out)
        hg_4342_obs   = 4341.684 * (1+z_out)
        hb_4863_obs   = 4862.683 * (1+z_out)
        oiii_5007_obs = 5008.240 * (1+z_out)
        ha_6565_obs   = 6564.610 * (1+z_out)
        sii_obs       = 6725.480 * (1+z_out)
        sii_6716_obs  = 6718.290 * (1+z_out)
        sii_6731_obs  = 6732.670 * (1+z_out)
        siii_9069_obs = 9071.100 * (1+z_out)
        siii_9532_obs = 9533.200 * (1+z_out)
        hei_10830_obs = 10832.86 * (1+z_out)


        ##will use this below for all the direct fitting.
        if config_pars['fitfluxes'] != True:
            cont_model = emissionline_model(modelpars_nolines, lam_spec)
            flam_cont_subbed = flux_spec - cont_model

        #####################################################
        ##### STEP 3  #### evaluate fluxes and return outputs.
        #### divide by zeros are handled, so we can suppress the warnings.
        with np.errstate(invalid='ignore', divide='ignore'):
            #if zguess_ha > 0:   ### if ha is in the bandpass:
            if ((ha_6565_obs > np.min(lam_spec)) & (ha_6565_obs < np.max(lam_spec) ) ) :
                if ha_6565_obs > transition_wave:
                    sig = out.params[8] / 2.35
                    sig_err = out.perror[8] / 2.35
                    covar_term = 2  * out.covar[8][9] / (out.params[8] * out.params[9])
                else:
                    sig = (out.params[8] * out.params[7])/2.35   #### fwhm_blue = out.params[7] * out.params[8] ; 7 is the ratio and is fixed.
                    sig_err = (out.perror[8] * out.params[7])/2.35
                    covar_term = 2  * out.covar[8][9] / (out.params[8] * out.params[9])    #### here the covar and the param would be multiplied by out.params[7] so they cancel.

                if config_pars['fitfluxes'] == True:
                    ha_flux = np.sqrt(2 * math.pi) *  out.params[9] * sig

                    hanii_flux = (1 + out.params[12]) * ha_flux # remember, nii flux is out.params[6] * ha_flux
                    if ha_flux > 0:
                        ha_err =  ha_flux * np.sqrt( (out.perror[9]/out.params[9])**2 + (sig_err/sig)**2 + covar_term)
                        hanii_err = (1 + out.params[12]) * ha_err
                    else :
                        w=np.where((lam_spec > ha_6565_obs - 2 * 2.35 * sig) & (lam_spec < ha_6565_obs + 2 * 2.35 * sig))   ### this is actually a rather poor approximation
                        hanii_err = np.sqrt(np.sum(error_spec[w]**2))

                else:
                    w=np.where((lam_spec > ha_6565_obs - 2 * 2.35 * sig) & (lam_spec < sii_obs + 2 * 2.35 * sig))
                    hanii_flux = integrate.trapz(flam_cont_subbed[w], lam_spec[w])
                    hanii_err = np.sqrt(np.sum(error_spec[w]**2))

                ha_cont = emissionline_model(modelpars_nolines, 6564 * np.array([1+out.params[3]] ))
                hanii_ew_obs = hanii_flux/ha_cont[0]

            else:
                hanii_flux =-1./scl
                hanii_err = -1./scl
                hanii_ew_obs = -1


            if ((oiii_5007_obs > np.min(lam_spec)) & (oiii_5007_obs < np.max(lam_spec) ) ) : ### if oiii 5007 is in the bandpass:
                if oiii_5007_obs >  transition_wave:
                    sig = out.params[8] / 2.35
                    sig_err = out.perror[8] / 2.35
                    covar_term = 2  * out.covar[8][13] / (out.params[8] * out.params[13])
                else:
                    sig = (out.params[7] * out.params[8] )/ 2.35
                    sig_err = (out.perror[8]  * out.params[7])/ 2.35
                    covar_term = 2  * out.covar[8][13] / (out.params[8] * out.params[13])

                if config_pars['fitfluxes'] == True:
                    oiii_flux = 1.3 * np.sqrt(2 * math.pi) *  out.params[13] * sig    #### to include both lines here.
                    if oiii_flux > 0:
                        oiii_err = oiii_flux * np.sqrt( (out.perror[13]/out.params[13])**2 + (sig_err/sig)**2 + covar_term)
                    else:
                        w=np.where((lam_spec > oiii_5007_obs - 2 * 2.35 * sig) & (lam_spec < oiii_5007_obs + 2 * 2.35 * sig))
                        oiii_err = np.sqrt(np.sum(error_spec[w]**2))
                else:
                     w=np.where((lam_spec > hb_4863_obs - 2 * 2.35 * sig) & (lam_spec < oiii_5007_obs + 2 * 2.35 * sig))
                     oiii_flux = integrate.trapz(flam_cont_subbed[w], lam_spec[w])
                     oiii_err = np.sqrt(np.sum(error_spec[w]**2))

                oiii_cont = emissionline_model(modelpars_nolines, 5008 * np.array( [1+out.params[3]] ) )
                oiii_ew_obs = oiii_flux / oiii_cont[0]

            else :
                 oiii_flux = -1./scl
                 oiii_ew_obs = -1.
                 oiii_err = -1./scl


            if ((sii_6716_obs > np.min(lam_spec)) & ( sii_6731_obs < np.max(lam_spec))):  ### sii covered.
                if sii_6716_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][10] / (out.params[8] * out.params[10])
                else:
                    sig = (out.params[8]*out.params[7])/2.35
                    sig_err = (out.params[7]* out.perror[8])/2.35
                    covar_term = 2  * out.covar[8][10] / (out.params[8] * out.params[10])

                sii_6716_flux = np.sqrt(2 * math.pi) *  out.params[10] * sig
                sii_6731_flux = sii_6716_flux / out.params[11]
                sii_flux = sii_6716_flux + sii_6731_flux
                if sii_flux > 0:
                    sii_6716_err = sii_6716_flux * np.sqrt((out.perror[10]/out.params[10])**2  + (sig_err/sig)**2 + covar_term)
                    sii_err = sii_6716_err * (1 + 1/out.params[11])   ### nominally this scale factor is 1.7.
                else:
                    w=np.where((lam_spec > sii_6716_obs - fwhm_guess) & (lam_spec < sii_6731_obs + fwhm_guess))
                    sii_err = np.sqrt(np.sum(error_spec[w]**2))
                sii_cont = emissionline_model(modelpars_nolines, 6725 *np.array([1+out.params[3]]))
                sii_ew_obs = sii_flux/sii_cont[0]
            else:
              sii_flux = -1/scl
              sii_err = -1 /scl
              sii_ew_obs = -1


            if ((hb_4863_obs > np.min(lam_spec)) & ( hb_4863_obs < np.max(lam_spec))):  ### hb covered.
                if hb_4863_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][14] / (out.params[8] * out.params[14])
                else:
                    sig = (out.params[7] * out.params[8])/2.35
                    sig_err = (out.params[7]* out.perror[8])/2.35
                    covar_term = 2  * out.covar[8][14] / (out.params[8] * out.params[14])

                hb_flux = np.sqrt(2 * math.pi) * out.params[14]  * sig
                if hb_flux > 0:
                    hb_err =  hb_flux * np.sqrt( (out.perror[14]/out.params[14])**2 + (sig_err/sig)**2  +  covar_term)
                else :
                    w=np.where((lam_spec > hb_4863_obs - fwhm_guess) & (lam_spec < hb_4863_obs + fwhm_guess))
                    hb_err = np.sqrt(np.sum(error_spec[w]**2))
                hb_cont = emissionline_model(modelpars_nolines, 4862 * np.array( [1 + out.params[3]] ))
                hb_ew_obs = hb_flux/hb_cont[0]
            else:
                hb_flux =-1/scl
                hb_err = -1/scl
                hb_ew_obs = -1


            if ((hg_4342_obs > np.min(lam_spec)) & ( hg_4342_obs < np.max(lam_spec))):  ### hg covered.
                if hg_4342_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][16] / (out.params[8] * out.params[16])
                else:
                    sig = (out.params[7] * out.params[8])/2.35
                    sig_err = (out.perror[8] * out.params[7]) /2.35
                    covar_term = 2  * out.covar[8][16] / (out.params[8] * out.params[16])

                hg_flux = np.sqrt(2 * math.pi) * out.params[16]  * sig
                if hg_flux > 0:
                    hg_err =  hg_flux * np.sqrt( (out.perror[16]/out.params[16])**2 + (sig_err/sig)**2  +  covar_term)
                else :
                    w=np.where((lam_spec > hg_4342_obs - fwhm_guess) & (lam_spec < hg_4342_obs + fwhm_guess))
                    hg_err = np.sqrt(np.sum(error_spec[w]**2))
                hg_cont = emissionline_model(modelpars_nolines, 4342 * np.array( [1 + out.params[3]] ))
                hg_ew_obs = hg_flux/hg_cont[0]
            else:
                hg_flux =-1/scl
                hg_err = -1/scl
                hg_ew_obs = -1


            if ((oii_3727_obs > np.min(lam_spec)) & ( oii_3727_obs < np.max(lam_spec))):  ### oii covered.
                if oii_3727_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][15] / (out.params[8] * out.params[15])
                else:
                    sig = (out.params[7] * out.params[8])/2.35
                    sig_err = (out.perror[8] * out.params[7]) /2.35
                    covar_term = 2  * out.covar[8][15] / (out.params[8] * out.params[15])

                oii_flux = 2.4 * np.sqrt(2 * math.pi) * out.params[15]  * sig    #### here the oii lines are also fit as a doublet, but the ratio is hard-coded.
                if oii_flux > 0:
                    oii_err =  oii_flux * np.sqrt( (out.perror[15]/out.params[15])**2 + (sig_err/sig)**2  +  covar_term)
                else:
                    w=np.where((lam_spec > oii_3727_obs - fwhm_guess) & (lam_spec < oii_3727_obs + fwhm_guess))
                    oii_err = np.sqrt(np.sum(error_spec[w]**2))
                oii_cont = emissionline_model(modelpars_nolines, 3727 * np.array( [1 + out.params[3]] ))
                oii_ew_obs = oii_flux/oii_cont[0]
            else:
                oii_flux =-1  /scl
                oii_err = -1/scl
                oii_ew_obs = -1

            if ((siii_9069_obs > np.min(lam_spec)) & ( siii_9069_obs < np.max(lam_spec))):  ### siii_9069 covered.
                if siii_9069_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][17] / (out.params[8] * out.params[17])
                else:
                    sig = (out.params[7] * out.params[8]) /2.35
                    sig_err = (out.perror[8] * out.params[7])/2.35
                    covar_term = 2  * out.covar[8][17] / (out.params[8] * out.params[17])

                siii_9069_flux = np.sqrt(2 * math.pi) * out.params[17]  * sig
                if siii_9069_flux > 0:
                    siii_9069_err =  siii_9069_flux * np.sqrt( (out.perror[17]/out.params[17])**2 + (sig_err/sig)**2  +  covar_term)
                else :
                    w=np.where((lam_spec > siii_9069_obs - fwhm_guess) & (lam_spec < siii_9069_obs + fwhm_guess))
                    siii_9069_err = np.sqrt(np.sum(error_spec[w]**2))
                siii_9069_cont = emissionline_model(modelpars_nolines, 9069 * np.array( [1 + out.params[3]] ))
                siii_9069_ew_obs = siii_9069_flux/ siii_9069_cont[0]
            else:
                siii_9069_flux =-1/scl
                siii_9069_err = -1/scl
                siii_9069_ew_obs = -1


            if ((siii_9532_obs > np.min(lam_spec)) & ( siii_9532_obs < np.max(lam_spec))):  ### siii_9069 covered.
                #### this block of code basically repeats the above block (rather than taking its results)
                #### and multiplies by 2.48.
                #### this implementation was confusing.
                if siii_9069_obs > transition_wave:    #### calculating the siii_9069_flux first and then scaling.  hence, use sig for 9069.
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][17] / (out.params[8] * out.params[17])
                else:
                    sig = (out.params[7] * out.params[8])/2.35
                    sig_err = (out.perror[8] * out.params[7])/2.35
                    covar_term = 2  * out.covar[8][17] / (out.params[8] * out.params[17])

                siii_9532_flux = 2.48* np.sqrt(2 * math.pi) * out.params[17]  * sig
                if siii_9532_flux > 0:
                    siii_9532_err =  siii_9532_flux * np.sqrt( (out.perror[17]/out.params[17])**2 + (sig_err/sig)**2  +  covar_term)
                else :
                    w=np.where((lam_spec > siii_9532_obs - fwhm_guess) & (lam_spec < siii_9532_obs + fwhm_guess))
                    siii_9532_err = np.sqrt(np.sum(error_spec[w]**2))
                siii_9532_cont = emissionline_model(modelpars_nolines, 9532 * np.array( [1 + out.params[3]] ))
                siii_9532_ew_obs = siii_9532_flux/ siii_9532_cont[0]

            else:
                siii_9532_flux =-1/scl
                siii_9532_err = -1/scl
                siii_9532_ew_obs = -1


            if ( (hei_10830_obs > np.min(lam_spec)) & (hei_10830_obs < np.max(lam_spec))):  ### he1 covered.
                if hei_10830_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][18] / (out.params[8] * out.params[18])
                else:
                    sig = (out.params[7] * out.params[8])/2.35
                    sig_err = (out.perror[8] * out.params[7])/2.35
                    covar_term = 2  * out.covar[8][18] / (out.params[8] * out.params[18])

                hei_flux =  np.sqrt(2 * math.pi) * out.params[18]  * sig
                if hei_flux > 0:
                    hei_err =  hei_flux * np.sqrt( (out.perror[18]/out.params[18])**2 + (sig_err/sig)**2  +  covar_term)
                else:
                    w=np.where((lam_spec > hei_10830_obs - fwhm_guess) & (lam_spec < hei_10830_obs + fwhm_guess))
                    hei_err = np.sqrt(np.sum(error_spec[w]**2))
                hei_cont = emissionline_model(modelpars_nolines, 10830. * np.array([1 + out.params[3]]))
                hei_ew_obs = hei_flux/hei_cont[0]
            else:
                hei_flux =-1/scl
                hei_err = -1/scl
                hei_ew_obs = -1

### M.D.R 01/06/2021 ###
            if ( (lya_1216_obs > np.min(lam_spec)) & (lya_1216_obs < np.max(lam_spec))):  ### lya covered.
                if lya_1216_obs > transition_wave:
                    sig = out.params[8]/2.35
                    sig_err = out.perror[8]/2.35
                    covar_term = 2  * out.covar[8][19] / (out.params[8] * out.params[19])
                else:
                    sig = (out.params[7] * out.params[8])/2.35
                    sig_err = (out.perror[8] * out.params[7])/2.35
                    covar_term = 2  * out.covar[8][19] / (out.params[8] * out.params[19])

                lya_flux =  np.sqrt(2 * math.pi) * out.params[19]  * sig
                if lya_flux > 0:
                    lya_err =  lya_flux * np.sqrt( (out.perror[19]/out.params[19])**2 + (sig_err/sig)**2  +  covar_term)
                else:
                    w=np.where((lam_spec > lya_1216_obs - fwhm_guess) & (lam_spec < lya_1216_obs + fwhm_guess))
                    lya_err = np.sqrt(np.sum(error_spec[w]**2))
                lya_cont = emissionline_model(modelpars_nolines, 1216. * np.array([1 + out.params[3]]))
                lya_ew_obs = lya_flux/lya_cont[0]
            else:
                lya_flux =-1/scl
                lya_err = -1/scl
                lya_ew_obs = -1
### M.D.R 01/06/2021 ###


        fit_results = {}
        fit_results['redshift'] = out.params[3]
        fit_results['redshift_err'] = out.perror[3]
        fit_results['oiii_5007_dz'] =  out.params[4]
        fit_results['oii_3727_dz'] =  out.params[5]
        fit_results['siii_hei_dz'] =  out.params[6]
        fit_results['fwhm_g141'] = out.params[8]
        fit_results['fwhm_g141_err'] = out.perror[8]
        fit_results['hanii_flux'] = hanii_flux * scl
        fit_results['hanii_error'] = hanii_err*scl
        fit_results['hanii_ew_obs'] = hanii_ew_obs
        fit_results['oiii_flux'] = oiii_flux * scl
        fit_results['oiii_error'] = oiii_err *scl
        fit_results['oiii_ew_obs'] = oiii_ew_obs
        fit_results['chisq'] = chisq
        fit_results['fit_status'] =  out.status
        fit_results['fit_scale_factor'] = scl
        fit_results['fit_parameters'] = out.params
        fit_results['sii_flux']  = sii_flux * scl
        fit_results['sii_error'] = sii_err* scl
        fit_results['sii_ew_obs'] = sii_ew_obs
        fit_results['hb_flux']  = hb_flux * scl
        fit_results['hb_error'] = hb_err* scl
        fit_results['hb_ew_obs'] = hb_ew_obs
        fit_results['hg_flux']  = hg_flux * scl
        fit_results['hg_error'] = hg_err * scl
        fit_results['hg_ew_obs'] = hg_ew_obs
        fit_results['oii_flux'] = oii_flux * scl
        fit_results['oii_error'] = oii_err*  scl
        fit_results['oii_ew_obs'] = oii_ew_obs
        fit_results['siii_9069_flux'] = siii_9069_flux *scl
        fit_results['siii_9069_error'] = siii_9069_err * scl
        fit_results['siii_9069_ew_obs'] = siii_9069_ew_obs
        fit_results['siii_9532_flux'] = siii_9532_flux * scl
        fit_results['siii_9532_error'] = siii_9532_err * scl
        fit_results['siii_9532_ew_obs'] = siii_9532_ew_obs
        fit_results['hei_flux'] = hei_flux  * scl
        fit_results['hei_error'] = hei_err  * scl
        fit_results['hei_ew_obs']  = hei_ew_obs
        fit_results['lya_flux'] = lya_flux  * scl # M.D.R 01/06/2021
        fit_results['lya_error'] = lya_err  * scl # M.D.R 01/06/2021
        fit_results['lya_ew_obs']  = lya_ew_obs # M.D.R 01/06/2021
        fit_results['lya_1216_dz'] =  out.params[20] # M.D.R 01/06/2021

    else:
        fit_results = {}
        fit_results['fit_status'] = out.status

    #return [out.params[0], out.params[1], ha_flux*scl, ha_err *scl, ha_ew_obs,  oiii_flux*scl, oiii_err * scl,  oiii_ew_obs, chisq, out.status, scl,  out.params]
    return fit_results

# # The below is not used and can be removed and/or commented out? # MDR 2022/05/25
# def fitandplot(path, par, beam_number, z_nate, fwhm_nate, config):
#       #### path = location of .dat spectra
#       #### par  = integer par number
#       ### beam_number = integer beam number
#       ### z_nate = z from line list or z_estimate
#       ### fwhm_nate = fwhm of line in G141, or 2 * FWHM of line in G102.
#       ### config = dictionary of configuration file
#
#       bluespec = path + 'MU1_' + str(beam_number) + '_G102.1D.dat' # bluespec = path + 'Par' + str(par) + '_G102_BEAM_' + str(beam_number) + 'A.dat' # - M.D.R. - 10/15/2020
#       redspec  = path + 'MU1_' + str(beam_number) + '_G141.1D.dat' # redspec  = path + 'Par' + str(par) + '_G141_BEAM_' + str(beam_number) + 'A.dat' # - M.D.R. - 10/15/2020
#
#       #config = read_config('refit.config')
#       tab_blue = asciitable.read(bluespec, names = ['lambda',  'flux', 'ferror', 'contam', 'junk'])
#       tab_red =  asciitable.read(redspec,  names = ['lambda',  'flux', 'ferror', 'contam', 'junk'])
#       spec = trim_spec(tab_blue, tab_red, config)
#       lam =spec[0]
#       flux = spec[1]
#       error = spec[2]
#       contam = spec[3]
#
#       input_list = [spec[0], spec[1]-spec[3], spec[2], config, z_nate,  fwhm_nate,  'notused']
#       out = fit_obj(input_list)  ### dumping the inputs into a list was used in a previously parallelized version.
#
#       ### define windows for model plotting.
#       ha_6565_obs   = 6565 * (1+z_nate)
#       sii_obs       = 6725 * (1+z_nate)
#       sii_6731_obs  = 6731 * (1+z_nate)
#       hb_4863_obs   = 4863 * (1+z_nate)
#       oiii_5007_obs = 5008 * (1+z_nate)
#       hg_4342_obs   = 4342 * (1+z_nate)
#       oii_3727_obs  = 3727 * (1+z_nate)
#       siii_9069_obs = 9071 * (1+z_nate)
#       siii_9532_obs = 9533 * (1+z_nate)
#       hei_10830_obs = 10833* (1+z_nate)
#       lya_1216_obs  = 1216 * (1+z_nate) # M.D.R 01/06/2021
#
#       lam_spec =lam
#       fit_region = config['fit_region']
#       flux_decontam = spec[1] - spec[3]
#
#       plt.plot(spec[0], flux_decontam, ls='steps-mid', color = 'k')
#       #w=np.where((lam_spec >= oii_obs - fit_region) & (lam_spec<hei_obs + fit_region)) # M.D.R 01/07/2021
#       w=np.where((lam_spec >= lya_1216_obs - fit_region) & (lam_spec<hei_10830_obs + fit_region)) # M.D.R 01/07/2021
#       model_fit = emissionline_model(out['fit_parameters'], lam[w]) * out['fit_scale_factor']
#       plt.plot(lam[w], model_fit, color = 'r')
#       fitpars = out['fit_parameters']
#       fitpars_nolines = cp.deepcopy(fitpars)
#       #fitpars_nolines[9:19] = 0. # M.D.R 01/06/2021
#       fitpars_nolines[9:20] = 0. # M.D.R 01/06/2021
#       fitpars_nolines[11] = 1.4  ### can't kill this one or divide by zero.
#       fitpars_nolines[12] = 0.1
#
#       cont_fit = emissionline_model(fitpars_nolines, lam[w]) * out['fit_scale_factor']
#       plt.plot(lam[w], cont_fit, color = 'blue',ls = '--', lw=2)
#      # print model_resid(out['fit_parameters'], lam = spec[0], flux = flux, err=error
#
#       plt.plot([ha_6565_obs, ha_6565_obs], [-1e-18, 6e-18], color='r')
#       plt.plot([hb_4863_obs, hb_4863_obs], [-1e-18, 6e-18], color = 'g')
#       plt.plot([oiii_5007_obs, oiii_5007_obs], [-1e-18, 6e-18], color = 'g')
#       plt.plot([hg_4342_obs, hg_4342_obs], [-1e-18, 6e-18], color= 'cyan')
#       plt.plot([oii_3727_obs, oii_3727_obs], [-1e-18, 6e-18], color = 'b')
#       plt.plot([siii_9069_obs, siii_9069_obs], [-1e-18, 6e-18], color = 'magenta')
#       plt.plot([siii_9532_obs, siii_9532_obs], [-1e-18, 6e-18], color = 'magenta')
#       plt.plot([hei_10830_obs, hei_10830_obs], [-1e-18, 6e-18], color= 'orange')
#       plt.plot([lya_1216_obs, lya_1216_obs], [-1e-18, 6e-18], color= 'orange') # M.D.R 01/06/2021
#       plt.xlim(config['lambda_min']-100, config['lambda_max']+100)
#
#       #w=np.where( (lam_spec > 9500)  & (lam_spec < 15500)) # M.D.R 01/07/2021
#       w=np.where( (lam_spec > 11500)  & (lam_spec < 15500)) # M.D.R 01/07/2021
#       ymax = np.max(model_fit)
#       #plt.ylim(-1e-18, 1.5 * ymax) # M.D.R 01/06/2021
#       plt.ylim(-1e-18, 1.2 * ymax) # M.D.R 01/06/2021
#
#       #fit_region = config['fit_region']
#       ### overplot the nodes
#       #dx = fit_region / 2
#       #clam1 = oii_obs - dx
#       #dx = (hg_obs - oii_obs)/2.
#       #clam2 = hg_obs - dx
#       #dx = (hb_obs - hg_obs)/2.
#       #clam3 = hb_obs - dx
#       #dx = (ha_obs - oiii_obs)/3.
#       #clam4 = ha_obs - 2 * dx
#       #clam5 = ha_obs -  dx
#       #dx = (siii_9069_obs - sii_6731_obs)/4
#       #clam5 = siii_9069_obs - 3 * dx
#       #clam6 = siii_9069_obs - 2 * dx
#       #clam7 = siii_9069_obs - dx
#       #clam8 = (siii_9069_obs + siii_9532_obs) / 2.
#       #clam9 = (siii_9532_obs + hei_obs) / 2.
#       #dx = fit_region/ 2
#       #clam10   =  hei_obs  + dx
#       #clam = np.array([clam1, clam2, clam3, clam4, clam5, clam6, clam7, clam8, clam9, clam10])
#
#       clam = np.array(config['node_wave'])
#       plt.plot(clam, clam * 0 + 0.5 * ymax, 'ko')
#       plt.title('BEAM ' + str(beam_number))
#
#       if not os.path.exists('fitfigs/'):
#         os.makedirs('fitfigs/')
#       plt.savefig('fitfigs/fig_beam' + str(beam_number) + '.pdf')
#
#       if config['showfig'] == True:
#           plt.show()
#       plt.close()
#       return out
