#!/usr/bin/env python
##########################################################################
##########################################################################
# get_z_interactive.py By Nathaniel R. Ross, UCLA, nross@astro.ucla.edu
# usage: python get_z_interactive.py line_list
# Reads in list of emission lines from the WISP survey HST/WFC3/IR Grisms and
# plots the spectra, iterates through the detected emission lines, allowing the
# user to make line identifications or reject spurious lines quickly and
# efficiently.
#
# Version 1.0 updates the program to look for wavelength-calibrated 2d grism
# stamps first. Also, will now look for default line list name and make
# linelistfile an optional parameter. Furthermore, I have added the option
# to save the incomplete line list to file as you progress, but default is NOT to
# do this.
# Optional arguments for go() are 1. linelistfile (String) - path to line list file
#                                 2. save_temp (Bool) - Save progress in
#                                    linelist/Par???lines_with_redshifts.incomplete.dat
#                                 3. recover_temp (Bool) - Recover progress from previous
#                                    session using the .incomplete.dat files
#
# ** Major change in version 1.0:
# We now use xpa instead of IRAF to display the 2d grism stamps and full 2d
# direct images (instead of cutout stamp). Reason for this change was the desire
# to show the grism stamps that have a wavelength solution applied. When using
# the IRAF display command, ds9 would fail to recognize this coordinate system.
# The XPA system can be downloaded here: http://hea-www.harvard.edu/RD/xpa/index.html
#
##########################################################################
import os
import shutil
from glob import glob
import tarfile
import distutils
import numpy as np # MDR 2022/05/17
import fileinput
import scipy
import pylab as plt
from scipy.interpolate import spline
from astropy.table import Table
from distutils.sysconfig import *  # question-- what is this for?
import sys
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
import matplotlib.patheffects as PathEffects
# Explicitly import readline to make the text entry process easier on OSX
import readline
# SQLLite database support for data persistence
from .WISPLFDatabaseManager import WISPLFDatabaseManager as WDBM

from wisp_analysis import *

verbose = True # MDR 2022/05/17

if (verbose == True):
    print('\nCompiling measure_z_interactive...\n') # MDR 2022/05/17
    print('The current working director is...\n') # MDR 2022/05/17
    print((os.getcwd())) # MDR 2022/05/17

#######################################################################

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

# Make all catalog header and data write commands loops over the 'flux_strings' variable.

'''
Note that the code expects that the strings in 'flux_strings' to match those in
'fitresults' with +'_flux' appended. This way the results can be looped over.
All emission line labels and parameters are defined at the top of the code,
except for the contamination flags, which must be defined within a function so
they reset for each object.
'''

supported_lines = [\
la_1216_vac, n5_1238_vac, c4_1548_vac, h2_1640_vac, o3_1660_vac, s3_1883_vac, \
c3_1907_vac, m2_2796_vac, o2_3727_vac, hg_4342_vac, o3_4363_vac, h2_4686_vac, \
hb_4863_vac, o3_5007_vac, o1_6300_vac, ha_6565_vac, s2_6716_vac, s3_9069_vac, \
s3_9532_vac, he10830_vac]

# These lines are close to their doublets so are not plotted in ax1.
# However, their fluxes and fitresults are still saved to the catalog.
supported_lines_extra = [\
n5_1242_vac, c4_1550_vac, \
o3_1666_vac, s3_1892_vac, c3_1909_vac, \
m2_2803_vac, o2_3730_vac, o3_4959_vac, o1_6363_vac, \
n2_6550_vac, n2_6585_vac, s2_6731_vac]

supported_lines_strings = [\
r'Ly$\alpha$', 'N V', 'C IV', 'He II', 'O III]', 'Si III]', \
'C III]', 'Mg II', '[O II]', r'H$\gamma$', '[O III]', 'He II', \
r'H$\beta$', '[O III]', '[O I]', r'H$\alpha$', '[S II]', '[S III]', \
'[S III]', 'He I']

flux_strings = [\
'la_1216_wing', 'n5_1238_1242', 'c4_1548_1550', 'h2_1640', 'o3_1660_1666', 's3_1883_1892', \
'c3_1907_1909', 'm2_2796_2803', 'o2_3727_3730', 'hg_4342', 'o3_4363', 'h2_4686', \
'hb_4863', 'o3_4959_5007', 'o1_6300_6363', 'ha_6550_6565_6585', 's2_6716_6731', 's3_9069', \
's3_9532', 'he10830']

contam_flags_string = \
'la_1216, n5_1238, n5_1242, c4_1548, c4_1550, h2_1640, \
o3_1660, o3_1666, s3_1883, s3_1892, c3_1907, c3_1909, \
m2_2796, m2_2803, o2_3727, o2_3730, hg_4342, o3_4363, \
h2_4686, hb_4863, o3_4959, o3_5007, o1_6300, o1_6363, \
n2_6550, ha_6565, n2_6585, s2_6716, s2_6731, s3_9069, \
s3_9532, he10830, cont'

# colors to help split up terminal output
# helpmsg = light blue
# obj = purple/magenta
# heading = green
# accept = green
# interim = cyan
# endc = clear color
setcolors = {'helpmsg':'\033[94m', \
             'obj':'\033[95m', \
             'heading':'\033[92m', \
             'accept':'\033[92m', \
             'working':'\033[0m', \
             'interim':'\033[36m', \
             'random':'\033[31m', \
             'endc':'\033[0m'}
#######################################################################


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def getzeroorders(zeroorderpath, g='G141', magcut=23.5):  # MB: changed from 23.0
    if (verbose == True): print('Running getzeroorders...\n') # MDR 2022/05/17
    """
    Changing to return a table
    """
    zop = open(zeroorderpath, 'r')
    zox = []
    zoy = []
    zoid = []
    zmag = []
    for line in zop:
        if len(line) > 60:
            linesplit = line.split()
            zox.append(float(linesplit[1][0:-1]))
            zoy.append(float(linesplit[2][0:-3]))
            zoid.append(int(linesplit[-2].split('{')[-1]))
            zmag.append(float(linesplit[-1][1:-2]))  # get mag from reg file

    zop.close()
    zoid = np.array(zoid)
    zoy = np.array(zoy)
    zox = np.array(zox)
    zmag = np.array(zmag)
    cond = (zmag <= magcut)
    t = Table([zox[cond], zoy[cond], zoid[cond]], names=('x', 'y', 'objid'))
    return t


def getfirstorders(firstorderpath):
    if (verbose == True): print('Running getfirstorders...\n') # MDR 2022/05/17
    """
    Changing to return a table
    """
    fop = open(firstorderpath, 'r')
    fox = []
    foy = []
    folen = []
    fowid = []
#    foid=[] nothing is done with 1st order IDs anymore
    for line in fop:
        # if line[0]!='#':
        linesplit = line.split()
        fox.append(float(linesplit[1][0:-1]))  # [0:-1] strips off the comma.
        foy.append(float(linesplit[2][0:-1]))
        folen.append(float(linesplit[3][0:-1]))
        # python is weird.
        fowid.append(float(linesplit[-1].split('{')[-1].split('}')[0]))
    print(foid)
    t = Table([fox, foy, folen, fowid], names=(
        'x', 'y', 'len', 'width', 'objid'))
    return t


def get_remaining_objects(full_obj_list, objid_done):
    if (verbose == True): print('Running get_remaining_objects...') # MDR 2022/05/17
    # array of flags
    wdone = np.in1d(np.array(full_obj_list), objid_done)
    remaining = np.copy(full_obj_list)
    mask = np.ones(full_obj_list.shape, dtype=bool)
    mask[wdone] = False
    remaining = remaining[mask]
    return remaining


def print_help_message():
    """
    Just putting this here to keep it out of the way.
    """
    msg = setcolors['helpmsg'] + "Available Options:\n"
    msg += setcolors['heading'] + "\tOBJECT SPECIFIC OPTIONS:\n"
    msg += setcolors['helpmsg'] + "\ta = accept object fit\n" \
        "\tac = accept object fit, noting contamination\n"  \
        "\tr = reject object\n" \
        "\tc = add comment\n" \
        "\tuser = toggle between previously saved fits\n" \
        "\tcontam = specify contamination to line flux and/or continuum\n" \
        "\treset = reset interactive options back to default for this object\n" \
        "\ts = print the (in progress) object summary\n\n"
    msg += setcolors['heading'] + "\tEMISSION LINE SPECIFIC OPTIONS:\n"
    msg += setcolors['helpmsg'] + "\tz = enter a different z guess\n" \
        "\tw = enter a different emission line wavelength guess\n" \
        "\tha,  or hb, o31, o32, o2, s2, s31, s32 = change redshift guess\n" \
        "\tdz = change the allowable redshift difference between lines\n" \
        "\tn = skip to next brightest line found in this object\n\n"
    msg += setcolors['heading'] + "\tSPECTRUM SPECIFIC OPTIONS:\n"
    msg += setcolors['helpmsg'] + "\tfw = change the fwhm guess in pixels\n" \
        "\tt = change transition wavelength\n" \
        "\tm1, m2, or m3 = mask up to three discontinuous wavelength regions\n" \
        "\tnodes = change the wavelengths for the continuum spline nodes\n" \
        "\tbluecut = change the blue cutoff of the G102 spec\n" \
        "\tredcut  = change the red cutoff of the G141 spec\n\n"
    msg += setcolors['heading'] + "\tDS9 SPECIFIC OPTIONS:\n"
    msg += setcolors['helpmsg'] + "\tlin = linear z-scale\n" \
        "\tlog = logarithmic\n" \
        "\tzs102 = z1,z2 comma-separated range for G102 zscale\n" \
        "\tzs141 = z1,z2 comma-separated range for G141 zscale\n" \
        "\tdc = recenter direct images\n" \
        "\treload = reload direct images\n" \
        "\tdr = reload direct image reg files\n\n"
    msg += setcolors['heading'] + "\tSOFTWARE SPECIFIC OPTIONS: \n"
    msg += setcolors['helpmsg'] + "\th = print this message\n" \
        "\tq = quit\n" + setcolors['endc']

    print(msg)


def check_masked_lines(fitresults, snr_meas_array, spdata):
    if (verbose == True): print('Running check_masked_lines...\n') # MDR 2022/05/17
    """ """
    spec_lam = spdata[0]

    z = fitresults['redshift']
    fwhm = fitresults['fwhm_g141']

    for i,(wave,line) in enumerate(zip(supported_lines,flux_strings)):
        waveobs = wave * (1. + z)
        w = (spdata[1].mask[(spec_lam >= (waveobs-fwhm/2.)) & (spec_lam <= (waveobs+fwhm/2.))]).any()
        if w:
            for dtype in ['flux', 'error', 'ew_obs']:
                fitresults['%s_%s'%(line,dtype)] = -1.0
            snr_meas_array[i] = 0.0

    return fitresults, snr_meas_array


def comment_out_obj(par, obj, catalogname):
    if (verbose == True): print('Running comment_out_obj...\n') # MDR 2022/05/17
    # if a row already exists for this object, comment it out
    # objstr = '{:<8d}'.format(par) + '{:<6d}'.format(obj)
    objstr = '{:<6d}'.format(obj)
    if os.path.exists(catalogname):
        for line in fileinput.input(catalogname, inplace=True):
            if objstr in line:
                print("#%s" % line, end=' ')
            else:
                print('%s' % line, end=' ')


def print_prompt(prompt, prompt_type='obj'):
    print((setcolors[prompt_type] + prompt + setcolors['endc']))


def write_object_summary(par, obj, fitresults, snr_meas_array, contamflags, summary_type='accept'):
    if (verbose == True): print('Running write_object_summary...\n') # MDR 2022/05/17
    """ """
    # string names for output
    linefluxes = np.array([fitresults['%s_flux'%fs] for fs in flux_strings])

    # initial message
    msg = setcolors[summary_type] + '#'*72
    msg += '\n## Par{} Obj {}:\n##   Fit Redshift: z = {:.4f}\n'.format(par, obj, fitresults['redshift'])

    # lines with S/N > 3
    good_snr = np.where(snr_meas_array > 3)
    msg += '##   Lines fit with S/N > 3:\n'
    for gsnr in good_snr[0]:
        # msg += '##\t%s: Flux = %.3e    S/N = %.2f\n'%(supported_lines_strings[gsnr],
        #                         linefluxes[gsnr], snr_meas_array[gsnr])
        msg += '##\t%s: Flux = %.3e    S/N = %.2f\n'%(flux_strings[gsnr],
                                linefluxes[gsnr], snr_meas_array[gsnr])

    cfout = ['%s:%i'%(cf,contamflags[cf]) for cf in contamflags if contamflags[cf] > 0]
    msg += '##   Contamination flags set:\n##\t' + ', '.join(cfout) + '\n'
    msg += '#'*72 + setcolors['endc']
    print(msg)


def make_tarfile(outdir):
    if (verbose == True): print('Running make_tarfile...\n') # MDR 2022/05/17
    """ """
    # copy defauly.config into output directory to keep a copy
    shutil.copy('default.config', outdir)
    with tarfile.open('%s.tar.gz'%outdir, 'w:gz') as tar:
        tar.add(outdir, arcname=os.path.basename(outdir))


def plot_object(zguess, zfit, spdata, config_pars, snr_meas_array, snr_tot_others, full_fitmodel, full_contmodel, current_lam, lamlines_found, index_of_strongest_line, contmodel, plottitle,outdir, zset=None):
    if (verbose == True): print('\nRunning plot_object...\n') # MDR 2022/05/17
    """
    # save the figure for everything, junk objects and all
    # previous figures are overwritten
    """
    # the expected wavelengths of emission lines given the zguess
    lamobs = (1 + zguess) * np.array(supported_lines)
    lamobs_extra = (1 + zguess) * np.array(supported_lines_extra)

    # define the filename that will be used for figures.
    plotfilename = os.path.join(outdir, 'figs', '%s_fit.png' % plottitle)

    # initialize the wavelength, flux, uncertainty, contamination, and zero order arrays.
    spec_lam = spdata[0]
    spec_val = spdata[1]
    spec_unc = spdata[2]
    spec_con = spdata[3]
    spec_zer = spdata[4]

    # apply the mask to the wavelength array
    masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)

    # limit the S/N array to MUSE so we can zoom-in on those lines.
    snr_muse = snr_meas_array[np.argwhere(lamobs < config_pars['mask_region1'][0])]

    # replace the occasional NaN value with zero to prevent gridspec crashes.
    snr_muse = np.nan_to_num(snr_muse)

    # determine the number of emission lines with S/N > 3.
    num_gridspec_plots = len(np.where(snr_muse >= 3.0)[0])

    # the maximum number of zoom-in plots is 7 for visibility.
    if (num_gridspec_plots > 7): num_gridspec_plots = 7

    #if (verbose == True): print 'There are '+str(num_gridspec_plots)+' lines with S/N > 3 in MUSE.\n'

    # generate the plot grid.
    plt.ion()
    fig = plt.figure(1, figsize=(11, 8), dpi=75)
    plt.clf()
    # gs = gridspec.GridSpec(3, 4)
    gs = gridspec.GridSpec(3, num_gridspec_plots+1)
    ax1 = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2:, 0:1])

    if (num_gridspec_plots == 1):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        all_ax = [ax1, ax2, ax3]
    elif (num_gridspec_plots == 2):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        all_ax = [ax1, ax2, ax3, ax4]
    elif (num_gridspec_plots == 3):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        all_ax = [ax1, ax2, ax3, ax4, ax5]
    elif (num_gridspec_plots == 4):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    elif (num_gridspec_plots == 5):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        ax7 = fig.add_subplot(gs[2:, 5:6])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
    elif (num_gridspec_plots == 6):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        ax7 = fig.add_subplot(gs[2:, 5:6])
        ax8 = fig.add_subplot(gs[2:, 6:7])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    elif (num_gridspec_plots >= 7):
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        ax7 = fig.add_subplot(gs[2:, 5:6])
        ax8 = fig.add_subplot(gs[2:, 6:7])
        ax9 = fig.add_subplot(gs[2:, 7:8])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    else:
        all_ax = [ax1, ax2]

    xmin = np.ma.min(spec_lam) - 10.0 # 200.0 - M.D.R. - 10/22/2020
    xmax = np.ma.max(spec_lam) + 10.0 # 200.0 - M.D.R. - 10/22/2020
    ymin = -0.2 * np.ma.max(spec_val)
    ymax = 1.2 * np.ma.max(spec_val)

    # the line widths for the data and overlaid fit.
    lw_data = 2.0
    lw_fits = 1.75

    ax1.plot(spec_lam, spec_val, 'k', spec_lam, spec_con, 'hotpink', ls='steps', lw=lw_data)
    ax1.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)

    # transforms for plotting in data and axes coordinates
    ax1trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax2trans = mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)

    # contamination model
    ax1.fill_between(spec_lam, spec_con, -1, color='#ff69b4', alpha=0.1, step='pre')

    # plot observed wavelengths of all the possible lines.
    for li, lstring, sn_meas in zip(lamobs, supported_lines_strings, snr_meas_array):
        if (li > xmin + 100) & (li < xmax - 100):
            for ax in all_ax: #[ax1, ax2, ax3, ax4]: # for ax in [ax1, ax2]: # - M.D.R. - 10/22/2020
                if (ax != ax2): # skip s/n plot.
                    ax.axvline(x=li, color='b')
                    # add blue lines for close doublets without text label.
                    # for li in lamobs_extra:
                    #     if (li > xmin + 1) & (li < xmax - 1):
                    #         ax.axvline(x=li, color='b', linestyle='--')
                    #
            stringplot = lstring + '  (' + str(round(sn_meas, 1)) + ') '

            # use data coordinates for x-axis and axes coords for y-axis
            ax1.text(li, 0.96, stringplot, rotation='vertical', ha='right', fontsize='16', transform=ax1trans)

    ax1.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
    ax1.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)

    # plot 0th orders
    w = np.where(spec_zer == 3)
    spec_zero_bad = spec_zer * 0 - 1
    spec_zero_bad[w] = 1.
    # mild zeroth orders
    w = np.where(spec_zer == 2)
    spec_zero_mild = spec_zer * 0 - 1
    spec_zero_mild[w] = 1.
    for ax in [ax1, ax2]:
        # use data coordinates for x-axis and axes coords for y-axis
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        if np.any(spec_zero_bad[spec_zero_bad != -1]):
            ax.fill_between(spec_lam, 0, 1, where=spec_zero_bad == 1,
                            color='red', alpha=0.3, transform=trans,
                            label='Major 0th order contam')
        if np.any(spec_zero_mild[spec_zero_mild != -1]):
            ax.fill_between(spec_lam, 0, 1, where=spec_zero_mild == 1,
                            color='orange', alpha=0.3, transform=trans,
                            label='Minor 0th order contam')

    # plot any masked regions
    for mr in ['mask_region1', 'mask_region2', 'mask_region3']:
        if (config_pars[mr][0] != 0.) & (config_pars[mr][1] != 0.):
            for ax in all_ax:#[ax1, ax2]:
                trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
                handles, labels = ax.get_legend_handles_labels()
                if 'masked regions' in labels:
                    maskedlabel = None
                else:
                    maskedlabel = 'masked regions'
                ax.fill_between(config_pars[mr], 0, 1, color='grey', alpha=0.3, transform=trans, label=maskedlabel)
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 0:
        ax1.legend(bbox_to_anchor=[1.0, 1.08], loc='upper right')

    # find values of spec_lam nearest to the nodes
    nodelam = config_pars['node_wave']
    nl_arr = []
    cont_node = []
    for nl in nodelam:
        w = np.argmin(np.abs(spec_lam - nl))
        nl_arr.append(spec_lam[w])
        cont_node.append(full_contmodel[w])
    ax1.plot(nl_arr, cont_node, 'ko', ms=9)

    # repeat for line_candidates
    lf_lam = []
    lf_cont = []
    for lf in lamlines_found:
        w = np.argmin(np.abs(spec_lam - lf))
        lf_lam.append(spec_lam[w])
        lf_cont.append(full_contmodel[w])
    ax1.plot(lf_lam, lf_cont, 'bo', ms=9)

    # indicate "current" line
#    current_lam = lamlines_found[index_of_strongest_line]
    current_cont = contmodel[np.argmin(np.abs(np.ma.compressed(masked_spec_lam) - current_lam))]
    ax1.plot(current_lam, current_cont, 'ro', ms=10)

    ax1.set_xlabel('Observed Wavelength ($\AA$)', size='xx-large')
    #ax1.set_ylabel(r'F$_\lambda$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$', size='xx-large')
    ax1.set_ylabel('Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)', size='xx-large')
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    ax1.set_title(plottitle, size='xx-large')

    # second panel for s/n
    s2n = (spec_val - full_contmodel) / spec_unc
    s2n_lam = spec_lam
    mask = np.logical_and(s2n > -10000., s2n < 10000.)
    s2n = s2n[mask]
    s2n_lam = s2n_lam[mask]
    ax2.plot(s2n_lam, s2n, 'k-', linestyle='steps', lw=lw_data)
    ymin = -0.2 * s2n.max()
    ymax = 1.2 * s2n.max() #ymax = 1.5 * s2n.max() # M.D.R 01/07/2021
    ax2.axhline(y=config_pars['n_sigma_above_cont'], c='r')
    # for li in lamobs:
    #     ax2.axvline(x=li, color='b')
    ax2.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)
    #ax2.set_xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
    #ax2.set_ylabel(r'S/N', size='xx-large')
    ax2.set_title(r'S/N', size='xx-large')
    ax2.set_xlim([xmin, xmax])
    ax2.set_ylim(ymin, ymax)
    # fig = plt.gcf() a

    # - M.D.R. - 10/22/2020
    # find max S/N in list of detected lines.
    # pull index location, determine lambda of that line.
    # set ymin/max based on centroid +/-i range based on dispersion.
    #print('snr_meas_array =', snr_meas_array)  # M.D.R 01/27/2021
    #print(np.max(snr_meas_array)) # M.D.R 01/27/2021
    #print(np.argmax(snr_meas_array)) # M.D.R 01/27/2021
    #lamobs_maxsnr = lamobs[np.argmax(snr_meas_array)]
    lamobs_blue = lamobs[np.argwhere(lamobs < config_pars['mask_region1'][0])]
    #lamobs_red = lamobs[np.argwhere(lamobs > config_pars['mask_region1'][1])]

    # DEFINE FUNCTION FOR ADDING EMISSION LINE LABELS TO ZOOM-IN PLOTS.
    def add_line_labels(ax, axtrans, xmin, xmax):
        #if (verbose == True): print 'Calling add_line_labels()...\n'
        for li, lstring, sn_meas in zip(lamobs, supported_lines_strings, snr_meas_array):
            if (li > xmin + 1) & (li < xmax - 1):
                ax.axvline(x=li, color='b')
                stringplot = lstring + '  (' + str(round(sn_meas, 1)) + ') '
                # use data coordinates for x-axis and axes coords for y-axis
                ax.text(li, 0.92, stringplot, ha='right', fontsize='14', transform=axtrans)
                #ax.set_title(stringplot, fontsize='14')
        return 'Added line label to', str(ax)

    def add_extra_lines(ax, axtrans, xmin, xmax):
        #if (verbose == True): print 'Calling add_extra_lines()...\n'
        for li in lamobs_extra:
            if (li > xmin + 1) & (li < xmax - 1):
                ax.axvline(x=li, color='b')#, linestyle='--')
                #stringplot = lstring + '  (' + str(round(sn_meas, 1)) + ') '
                # use data coordinates for x-axis and axes coords for y-axis
                #ax.text(li, 0.92, stringplot, ha='right', fontsize='14', transform=axtrans)
                #ax.set_title(stringplot, fontsize='14')
        return 'Added line label to', str(ax)

    # array of the indices for lines with s/n > 3 from strongest to weakest.
    snr_muse_idx = np.flip(np.argsort(snr_muse, axis=0))

    # loop over lines with s/n > 3 and create zoom-in plots.
    for x in range(num_gridspec_plots):
        #print 'x =', x
        if (snr_muse[snr_muse_idx[x]][0] >= 3.0):

            idx_i = snr_muse_idx[x][0]
            snr_i = snr_muse[snr_muse_idx[x]][0][0]
            lam_i = lamobs_blue[snr_muse_idx[x]][0][0]
            # print 'idx_i = ', idx_i
            # print 'snr_i = ', snr_i
            # print 'lam_i = ', lam_i

            # the wavelength array index corresponding to the detected line wavelength.
            snr_muse_idy = min(enumerate(spec_lam), key = lambda x: abs(lam_i - x[1]))[0]
            # half-width of zoom-in plots in angstroms
            half_xrange_ang = 30
            # half-width of zoom-in plots in pixels
            half_xrange_pix = int((half_xrange_ang / config_pars['dispersion_blue']))

            xmin = lam_i-half_xrange_ang
            xmax = lam_i+half_xrange_ang

            # update to a sigma rejected min?
            ymin = 0.9*np.min(spec_val[snr_muse_idy-half_xrange_pix:snr_muse_idy+half_xrange_pix])
            ymax = 1.2*np.max(spec_val[snr_muse_idy-half_xrange_pix:snr_muse_idy+half_xrange_pix])

            if (x == 0):
                ax3trans = mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
                ax3.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax3.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax3.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax3.set_xlim([xmin, xmax])
                ax3.set_ylim([ymin, ymax])
                add_line_labels(ax3, ax3trans, xmin, xmax)
                add_extra_lines(ax3, ax3trans, xmin, xmax)
            elif (x == 1):
                ax4trans = mtransforms.blended_transform_factory(ax4.transData, ax4.transAxes)
                ax4.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax4.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax4.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax4.set_xlim([xmin, xmax])
                ax4.set_ylim([ymin, ymax])
                add_line_labels(ax4, ax4trans, xmin, xmax)
                add_extra_lines(ax4, ax4trans, xmin, xmax)
            elif (x == 2):
                ax5trans = mtransforms.blended_transform_factory(ax5.transData, ax5.transAxes)
                ax5.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax5.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax5.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax5.set_xlim([xmin, xmax])
                ax5.set_ylim([ymin, ymax])
                add_line_labels(ax5, ax5trans, xmin, xmax)
                add_extra_lines(ax5, ax5trans, xmin, xmax)
            elif (x == 3):
                ax6trans = mtransforms.blended_transform_factory(ax6.transData, ax6.transAxes)
                ax6.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax6.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax6.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax6.set_xlim([xmin, xmax])
                ax6.set_ylim([ymin, ymax])
                add_line_labels(ax6, ax6trans, xmin, xmax)
                add_extra_lines(ax6, ax6trans, xmin, xmax)
            elif (x == 4):
                ax7trans = mtransforms.blended_transform_factory(ax7.transData, ax7.transAxes)
                ax7.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax7.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax7.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax7.set_xlim([xmin, xmax])
                ax7.set_ylim([ymin, ymax])
                add_line_labels(ax7, ax7trans, xmin, xmax)
                add_extra_lines(ax7, ax7trans, xmin, xmax)
            elif (x == 5):
                ax8trans = mtransforms.blended_transform_factory(ax8.transData, ax8.transAxes)
                ax8.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax8.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax8.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax8.set_xlim([xmin, xmax])
                ax8.set_ylim([ymin, ymax])
                add_line_labels(ax8, ax8trans, xmin, xmax)
                add_extra_lines(ax8, ax8trans, xmin, xmax)
            elif (x == 6):
                ax9trans = mtransforms.blended_transform_factory(ax9.transData, ax9.transAxes)
                ax9.plot(spec_lam, spec_val, 'k', lw=lw_data)
                ax9.plot(spec_lam, full_fitmodel, color='r', lw=lw_fits)
                ax9.plot(spec_lam, full_contmodel, color='b', linestyle='--', lw=lw_fits)
                ax9.set_xlim([xmin, xmax])
                ax9.set_ylim([ymin, ymax])
                add_line_labels(ax9, ax9trans, xmin, xmax)
                add_extra_lines(ax9, ax9trans, xmin, xmax)
            else:
                print('Failed to set XY limits.')

            # try:
            #     print 'eval(ax+str(x+3)) = ', eval('ax'+str(x+3))
            #     eval('ax'+str(x+3)).set_xlim([lam_i-15, lam_i+15])
            #     eval('ax'+str(x+3)).set_ylim([-0.2*1.2*np.max(spec_val[lam_i-15:lam_i+15]), 1.2*np.max(spec_val[lam_i-15:lam_i+15])])
            # except:
            #     print '\nFailed to set XY limits.'
            #     pass

    if zset is None:
        addtext = 'In progress (z = {:.4f}'.format(zfit) +  ', ' + 'wSNR = {:.2f})'.format(snr_tot_others)
        addtextcolor = 'orange'
    elif zset == 0:
        addtext = 'Rejected'
        addtextcolor = 'red'
    elif zset == 1:
        addtext = 'Accepted (z = {:.4f}'.format(zfit) +  ', ' + 'wSNR = {:.2f})'.format(snr_tot_others)
        addtextcolor = 'green'

    fig.text(0.3, 0.96, addtext, ha='right', va='bottom', color=addtextcolor,
             fontsize=18, fontweight=500,
             path_effects=[PathEffects.withStroke(linewidth=0.5,foreground="k")])

    plt.tight_layout() # MDR 2022/05/19 - Added to fill screen when fitting.
    fig.savefig(plotfilename)
    plt.draw()

    #if (verbose == True): print('Exiting plot_object...\n') # MDR 2022/05/17


def inspect_object(user, par, obj, objinfo, lamlines_found, ston_found, g102zeros, g141zeros, linelistoutfile, commentsfile, remaining, allobjects, show_dispersed=True, stored_fits = False, path_to_wisp_data = ' '):
    if (verbose == True): print('Running inspect_object...\n') # MDR 2022/05/17
    """
    An attempt to move all object-specific tasks
    """
    # set up and filenames
    outdir = 'Par%s_output_%s'%(par,user)
    if path_to_wisp_data ==  ' ':
        specnameg102 = 'MU1_%i_G102.1D.dat' % (obj) # specnameg102 = 'Par%i_G102_BEAM_%iA.dat' % (par, obj) - M.D.R. - 10/08/2020 ... From e.g. Par406_G141_BEAM_422A.dat
        specnameg141 = 'MU1_%i_G141.1D.dat' % (obj) # specnameg141 = 'Par%i_G141_BEAM_%iA.dat' % (par, obj) - M.D.R. - 10/08/2020 ... To e.g. MU1_295_G141.1D.dat
    else :
        specnameg102 = path_to_wisp_data + '/Par' + str(par) + '/Spectra/MU1_%i_G102.1D.dat' % (obj)
        specnameg141 = path_to_wisp_data + '/Par' + str(par) + '/Spectra/MU1_%i_G141.1D.dat' % (obj)
    #plottitle = 'Par%i_%i' % (par, obj) # plottitle = 'Par%i_BEAM_%i' % (par, obj) - M.D.R. - 10/08/2020
    plottitle = 'mudf_%i' % (obj)
    fitdatafilename = os.path.join(outdir, 'fitdata/%s_fitspec' % plottitle)
    # read in 1D spectrum
    availgrism = ''
    if os.path.exists(specnameg102):
        availgrism += 'g102'
        tab_blue = asciitable.read(
            specnameg102, names=['lambda', 'flux', 'ferror', 'contam', 'zero'])
    else:
        tab_blue = None
    # check for g141, too. there are maybe 4 G102-only fields
    if os.path.exists(specnameg141):
        availgrism += 'g141'
        tab_red = asciitable.read(
            specnameg141, names=['lambda', 'flux', 'ferror', 'contam', 'zero'])
    else:
        tab_red = None

    availgrism = 'both' if availgrism == 'g102g141' else availgrism

    # display the object
    if g102zeros is not None:
        #show2dNEW('G102', par, obj, g102firsts, g102zeros, 'linear')
        show2dNEW('G102', par, obj, g102zeros, user, 'linear', path_to_wisp_data = path_to_wisp_data)
    if g141zeros is not None:
        show2dNEW('G141', par, obj, g141zeros, user, 'linear', path_to_wisp_data = path_to_wisp_data)
    # pan full images to the new object
    showDirectNEW(obj, par, path_to_wisp_data = path_to_wisp_data)
    if show_dispersed:
        showDispersed(obj, par, path_to_wisp_data = path_to_wisp_data)

    # define parameters for this object
    ra = objinfo['ra']
    dec = objinfo['dec']
    a_image = objinfo['a_image']
    b_image = objinfo['b_image']
    jmag = objinfo['jmag']
    jerr = objinfo['jerr']
    hmag = objinfo['hmag']
    herr = objinfo['herr']

    # start with a fresh set of config pars
    config_pars = read_config('default.config', availgrism=availgrism)

    # Data have been successfully loaded for this object. If it has been inspected
    # previously, the original results will have been stored in the SQLLite database
    # and a retrieval obtion should be offered.
#    databaseManager = WDBM(dbFileNamePrefix='Par{}'.format(par))
    databaseManager = WDBM(dbFileNamePrefix=os.path.join(outdir,'Par{}'.format(par))) # databaseManager = WDBM(dbFileNamePrefix=os.path.join(outdir,'Par{}'.format(par))) - M.D.R. - 10/08/2020
    mostRecentObjectData = databaseManager.getMostRecentObject()
    # if mostRecentObjectData is not None:
    #     print('Most recent object in database: Par {}, Obj {}, Date {}'.format(*mostRecentObjectData))
    # else:
    #     print('Database is empty.')
    catalogueEntryData = databaseManager.loadCatalogueEntry(parNumber=par, objectId=obj)

    ### # MDR 2022/05/17
    # GET REFIT OPTION WORKING.
    print(('\nSearching for previous fits to object {}...\n'.format(str(obj))))
    #print('file = ', outdir + '/done_%s'%user)

    if os.path.exists(outdir + '/done_%s'%user): # needed for the first run before file exists
        with open(outdir + '/done_%s'%user) as f:
            for index, line in enumerate(f):
                #print("Line {}: {}".format(index, line.strip()))
                if (int(line.strip()) == obj):
                    catalogueEntryData = 1
                    print('Match found...\n')
                    break

        if (catalogueEntryData != 1):
            print('No match found...\n')

    print('Done searching for previous fits...')
    ### # MDR 2022/05/17

    rejectPrevFit = True

    # if (verbose == True):
    #     print('databaseManager = ', databaseManager) # MDR 2022/05/17
    #     print('mostRecentObjectData = ', mostRecentObjectData) # MDR 2022/05/17
    #     print('catalogueEntryData = ', catalogueEntryData) # MDR 2022/05/17
    #     print('rejectPrevFit = ', rejectPrevFit) # MDR 2022/05/17

    # if catalogueEntryData is not None :
    if catalogueEntryData == 1:
        #nonFitResults, fitResults = catalogueEntryData
        #(par_db, obj_db, ra_db, dec_db, jmag_db, hmag_db, a_image_db, b_image_db, contamflag_db, entrytime_db) = nonFitResults
#        print('Found previous fit results for Pointing {}, Object {}.\nEnter "y" to accept the earlier fit.'.format(par_db, obj_db))
        print(('You have already fit Obj {}. Refit? [y/N]').format(obj))
        rejectPrevFit = input('> ').strip().lower() == 'y'
        if rejectPrevFit:
            comment_out_obj(par, obj, linelistoutfile)

#        print('Accepting previous fit.' if acceptPrevFit else 'Re-fitting this object.')
#    else :
#        print('No previous fit results found. Fitting this object now.')

    # get line, fwhm, z estimate
    # choose the lamline that has the highest S/N estimate
    s = np.argsort(ston_found)
    # reverse s/n order

    ston_found = ston_found[s[::-1]]
    lamlines_found = lamlines_found[s[::-1]]
    index_of_strongest_line = 0
    lamline = lamlines_found[index_of_strongest_line]

    # MDR 2022/06/10 - Changed the first guess line to [O II] for MUSE.
    #zguess = lamline / ha_6565_vac - 1
    #zguess = lamline / ((o2_3727_vac + o2_3730_vac) / 2.0) - 1.0
    zguess = (lamline / o2_3730_vac) - 1.0
    # fwhm is defined for the red side, regardless of where line is
    fwhm_guess = 2.35 * a_image * config_pars['dispersion_red']

    if stored_fits != False:
        first_stored_fit = stored_fits[0]
        users = [path.split('/')[-3].split('_')[-1] for path in stored_fits]
        fileObject = open(first_stored_fit,'r')
        alldata = pickle.load(fileObject)
        config_pars = alldata[10]
        fitresults_old = alldata[8]
        zguess = fitresults_old['redshift']
        fwhm_guess = fitresults_old['fwhm_g141']
        print('using stored fit from: '  + users[0])
        print('available stored fits: ')
        print(users)
        ### also need to figure out what else to add?
            ### config pars for nodes can also be entered here.

### replace this with printouts from pickle files

    # print object info to screen
    if rejectPrevFit:
        print(' ')
        print_prompt("=" * 72)
        print_prompt("Par%i Obj %i:" % (par, obj))
        print_prompt("Initial redshift guess: z = %f" % (zguess))
        print_prompt("\nWhat would you like to do with this object?\nSee the README for options, or type 'h' to print them all to the screen.")

    comment = ''

    # set this to loop over flags created at top of code.
    contamflags = {\
    'la_1216':0, 'n5_1238':0, 'n5_1242':0, 'c4_1548':0, 'c4_1550':0, 'h2_1640':0, \
    'o3_1660':0, 'o3_1666':0, 's3_1883':0, 's3_1892':0, 'c3_1907':0, 'c3_1909':0, \
    'm2_2796':0, 'm2_2803':0, 'o2_3727':0, 'o2_3730':0, 'hg_4342':0, 'o3_4363':0, \
    'h2_4686':0, 'hb_4863':0, 'o3_4959':0, 'o3_5007':0, 'o1_6300':0, 'o1_6363':0, \
    'n2_6550':0, 'ha_6565':0, 'n2_6585':0, 's2_6716':0, 's2_6731':0, 's3_9069':0, \
    's3_9532':0, 'he10830':0, 'cont':0} # MDR 2022/07/22

    # Skip if previous fit is to be accepted
    done = 0 if rejectPrevFit else 1
    fast_fit = False # MDR 2022/06/30 - move to configuration file?

    while (done == 0):
        # sticking with the while loop to determine whether user is finished with object
        # get spectrum for obj. do this every time because sometimes we
        # re-read with a mask or a different transition wavelength
        spdata = trim_spec(tab_blue, tab_red, config_pars, mask_zeros=True, return_masks=True)

        spec_lam = spdata[0]
        spec_val = spdata[1]
        spec_unc = spdata[2]
        spec_con = spdata[3]
        spec_zer = spdata[4]
        mask_flg = spdata[5]

        # Determine the largest extent of the object so broadening of the lines can be accounted for in the fitting. MDR 2022/06/30
        ab_image_max = np.max([objinfo['a_image'][0], objinfo['b_image'][0]])

        # apply the mask to the wavelength array
        masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)
        # compress the masked arrays for fitting
        fit_inputs = [np.ma.compressed(masked_spec_lam), np.ma.compressed(spec_val), np.ma.compressed(spec_unc), config_pars, zguess, fwhm_guess, str(obj), ab_image_max, fast_fit]
        # parsing the input to facilitate parallel processing when fitting is done in batch mode.
        fitresults = fit_obj(fit_inputs)
        zfit = fitresults['redshift']
        fitpars = fitresults['fit_parameters']
        fitpars_nolines = cp.deepcopy(fitpars)

        ############################################################################
        '''
        The model parameters for the emission line amplitudes must be set to zero
        for the continuum fit, while the values for the line ratios must be non-zero
        to avoid division by zero. The 'get_fitpar_indices()' and 'get_ratio_indices()'
        functions below are defined in fitting.py to send these model parameter indices
        back to measure_z_interactive().
        '''
        first_line_index, first_node_index = get_fitpar_indices()

        fitpars_nolines[first_line_index:first_node_index] = 0.

        for idx in get_ratio_indices():
            fitpars_nolines[idx] = 1.0
        ############################################################################

        fitmodel = emissionline_model(fitpars, np.ma.compressed(masked_spec_lam)) * fitresults['scl_factor']
        contmodel = emissionline_model(fitpars_nolines, np.ma.compressed(masked_spec_lam)) * fitresults['scl_factor']

        # the fitting is done on compressed arrays, so we need to
        # create masked versions of the fit and continuum models
        full_fitmodel = np.zeros(spec_lam.shape, dtype=float)
        full_contmodel = np.zeros(spec_lam.shape, dtype=float)
        full_fitmodel[np.ma.nonzero(spec_val)] = fitmodel
        full_contmodel[np.ma.nonzero(spec_val)] = contmodel
        full_fitmodel = np.ma.masked_where(np.ma.getmask(spec_val), full_fitmodel)
        full_contmodel = np.ma.masked_where(np.ma.getmask(spec_val), full_contmodel)

        # loop over the lines specified in 'flux_strings' and save the results to an array.
        snr_meas_array = []
        for line in flux_strings:
            snr_meas_array.append(fitresults[line+'_flux'] / fitresults[line+'_error'])
        snr_meas_array = np.array(snr_meas_array)

        signal_lines = []
        for line in flux_strings:
            signal_lines.append(fitresults[line+'_flux'])
        signal_lines = np.array(signal_lines)

        err_lines = []
        for line in flux_strings:
            err_lines.append(fitresults[line+'_error'])
        err_lines = np.array(err_lines)

        # calculate a weighted s/n ratio for all lines and print on the plot.
        #signal_lines = np.array([fitresults['o2_total_flux'], fitresults['hg_4342_flux'], fitresults['hb_4863_flux'],fitresults['ha_total_flux'], fitresults['s2_total_flux']])
        #err_lines   = np.array([fitresults['o2_total_error'], fitresults['hg_4342_error'], fitresults['hb_4863_error'], fitresults['ha_total_error'], fitresults['s2_total_error']])

        w = np.where(signal_lines > 0)

        # MDR 2022/06/10 - updated the definition of weight SNR to be sum of SNRs that are  > 3 weighted by their contribution to the total flux.
        total_flux = 0
        for line in flux_strings:
            snr_line = (fitresults[line+'_flux'] / fitresults[line+'_error'])
            if (snr_line >= 3.0):
                total_flux = (total_flux + fitresults[line+'_flux'])
        # MDR 2022/06/10
        snr_tot_others = []
        for line in flux_strings:
            snr_line = (fitresults[line+'_flux'] / fitresults[line+'_error'])
            if (snr_line >= 3.0):
                snr_tot_weight = (fitresults[line+'_flux'] / total_flux)
                if (np.isfinite(snr_line * snr_tot_weight)):
                    snr_tot_others.append(snr_line * snr_tot_weight)
                else:
                    snr_tot_others.append(0.0)
        # MDR 2022/06/10
        snr_tot_others = np.sum(snr_tot_others)

        # plot the whole darn thing
        plot_object(zguess, fitresults['redshift'],
                    spdata, config_pars, snr_meas_array, snr_tot_others, full_fitmodel,
                    full_contmodel, lamline, lamlines_found,
                    index_of_strongest_line, contmodel, plottitle, outdir)
#        print "    Guess Redshift: z = %f" % (zguess)
        print_prompt("    Fit Redshift:   z = %f\n" % (zfit))
        # input
        option = input("> ")

        # checking user's input. keeping this format the same as before
        # any time done is set to 1, the object is considered fit

        # reject object
        if option.strip().lower() == 'r':
           zset = 0
           if len(comment) > 0:
               comment = 'rejected'  + ', ' + comment
           else :
                comment = 'rejected'
           done = 1

        # # accept object
        # elif option.strip().lower() == 'a':
        #     done = 1
        #     zset = 1
        #     flagcont = 1
        # accept object MDR 2022/06/30
        elif option.strip().lower() == 'a':
            if fast_fit == False:
                done = 1
                zset = 1
                flagcont = 1
            elif fast_fit == True:
                print('\nWARNING: Still using fast fit mode, type full for refined fit.')

        ### MDR 2022/06/30
        elif option.strip().lower() == 'full':
            fast_fit = False
        elif option.strip().lower() == 'fast':
            fast_fit = True
        ### MDR 2022/06/30

        # accept object and note contamination
        elif option.strip().lower() == 'ac':
            done = 1
            zset = 1
            flagcont = 2
            # add to contamination flags
            for k,v in contamflags.items():
                contamflags[k] = contamflags[k] | 1

        # change redshift guess
        elif option.strip().lower() == 'z':
            print_prompt("The current redshift guess is: %f\nEnter Redshift Guess:" % zguess)
            try:
                zguess = float(input("> "))
            except ValueError:
                print_prompt('Invalid Entry.')

        # change wavelength guess
        elif option.strip().lower() == 'w':
            print_prompt("The current emission line wavelength is: %f\nEnter Wavelength Guess in Angstroms:" % lamline)
            # save previous line guess (if not Ha)
            old_rest_wave = lamline / (1. + zguess)
            try:
                newwave = float(input("> "))
            except ValueError:
                print_prompt('Invalid Entry.')
            else:
                zguess = newwave / old_rest_wave - 1.
                lamline = newwave

        # change the fwhm guess
        elif option.strip().lower() == 'fw':
            print_prompt("Enter a Guess for FWHM in pixels")
            print_prompt("The current fwhm_fit is:  " + str(fitresults['fwhm_g141'] / config_pars['dispersion_red']) + " and 2*A_image is: " + str(2 * a_image[0]))
            try:
                fwhm_guess = config_pars['dispersion_red'] * float(input("> "))
            except ValueError:
                print_prompt('Invalid Entry.')

        # change delta_z range
        elif option.strip().lower() == 'dz':
            print_prompt("Enter a new delta_z limit")
            print_prompt("The current delta_z is:  " + str(config_pars['delta_z']))
            try:
                config_pars['delta_z'] = float(input("> "))
            except ValueError:
                print_prompt('Invalid Entry.')

        # mask out 1, 2, or 3 regions of the spectrum
        elif option.strip().lower() == 'm1':
            print_prompt("Enter wavelength window to mask out:  blue, red:")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]),
                            float(maskstr.split(",")[1])]
            except (IndexError,ValueError):
                print_prompt('Invalid entry. Enter wavelengths separated by commas')
            else:
                config_pars['mask_region1'] = maskwave

        elif option.strip().lower() == 'm2':
            print_prompt("Enter wavelength window to mask out:  blue, red:")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]),
                            float(maskstr.split(",")[1])]
            except (IndexError,ValueError):
                print_prompt('Invalid entry. Enter wavelengths separated by commas')
            else:
                config_pars['mask_region2'] = maskwave

        elif option.strip().lower() == 'm3':
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]),
                            float(maskstr.split(",")[1])]
            except (IndexError,ValueError):
                print_prompt('Invalid entry. Enter wavelengths separated by commas')
            else:
                config_pars['mask_region3'] = maskwave


        # change the transition wavelength between the grisms
        elif option.strip().lower() == 't':
            print_prompt("The current transition wavelength is: " + str(config_pars['transition_wave']) + "\nEnter the wavelength for the G102 to G141 transition:")
            try:
                config_pars['transition_wave'] = float(input("> "))
            except ValueError:
                print_prompt('Invalid entry. Enter wavelength of grism transition.')


        # change the nodes used for the continuum spline
        elif option.strip().lower() == 'nodes':
            strnw = ','.join(str(nw) for nw in config_pars['node_wave'])
            print_prompt("Enter Wavelengths for Continuum Spline: w1, w2, w3, w4, ....")
            print_prompt("current node wavelengths are: %s" % strnw)
            nodestr = input("> ")
            nodesplit = nodestr.split(',')
            node_arr = []
            try:
                for nodelam in nodesplit:
                    node_arr.append(float(nodelam))
            except ValueError:
                print_prompt('Invalid entry. Enter wavelengths separated by commas')
            else:
                node_arr = np.array(node_arr)
                # sort by wavelength
                node_arr = np.sort(node_arr)
                config_pars['node_wave'] = node_arr


        elif option.strip().lower() == 'user':

            if stored_fits !=  False:
                print_prompt("Enter name of user for toggling between stored fits")
                user_input = input("> ")
                try:
                    w = users.index(user_input)
                except ValueError:
                    print_prompt('Invalid entry. Enter a valid user name.')

                different_stored_fit = stored_fits[w]
                fileObject = open(different_stored_fit,'r')
                alldata = pickle.load(fileObject)
                config_pars = alldata[10]
                fitresults_old = alldata[8]
                zguess = fitresults_old['redshift']
                fwhm_guess = fitresults_old['fwhm_g141']
            else:
                print('there are no stored fits!')


        # reset all options
        elif option.strip().lower() == 'reset':
            print_prompt("Reset configuration parameters, fwhm guess, and zguess to default values")
            config_pars = read_config('default.config', availgrism=availgrism)
            fwhm_guess = 2.35 * a_image * config_pars['dispersion_red']
            # reset strongest line, too
            index_of_strongest_line = 0
            lamline = lamlines_found[index_of_strongest_line]
            #zguess = lamline / ha_6565_vac - 1
            #zguess = lamline / ((o2_3727_vac + o2_3730_vac) / 2.0) - 1.0
            zguess = (lamline / o2_3730_vac) - 1.0
            # reset contamflags
            for k, v in contamflags.items():
                    contamflags[k] = contamflags[k] & 0
            ### if use stored = true, this should set us back to using the pickle file values

        # change the blue cutoff of G102 (or whichever grism is present?)
        elif option.strip().lower() == 'bluecut':
            print_prompt("The current blue cutoff is: " + str(config_pars['lambda_min']) + "\nChange the blue cutoff of G102:")
            try:
                config_pars['lambda_min'] = float(input("> "))
            except ValueError:
                print_prompt('Invalid entry. Enter wavelength of blue cutoff.')


        # change the red cutoff of G141
        elif option.strip().lower() == 'redcut':
            print_prompt("The current red cutoff is: " + str(config_pars['lambda_max']) + "\nChange the red cutoff of G141:")
            try:
                config_pars['lambda_max'] = float(input("> "))
            except ValueError:
                print_prompt('Invalid entry. Enter wavelength of red cutoff.')


        # change to next brightest line
        elif option.strip().lower() == 'n':
            nlines_found_cwt = np.size(lamlines_found)
            index_of_strongest_line = index_of_strongest_line + 1
            if index_of_strongest_line < (nlines_found_cwt):
                lamline = lamlines_found[index_of_strongest_line]
                #zguess = lamline / 6564.610 - 1.0
                #zguess = lamline / ((o2_3727_vac + o2_3730_vac) / 2.0) - 1.0
                zguess = (lamline / o2_3730_vac) - 1.0
            else:
                print_prompt('There are no other automatically identified peaks. Select another option.')
                # stay at current line
                index_of_strongest_line -= 1


        # change to another line
        elif option.strip().lower() == 'lya':
            zguess = (lamline / la_1216_vac) - 1.0
        elif option.strip().lower() == 'c4':
            zguess = (lamline / c4_1548_vac) - 1.0
        elif option.strip().lower() == 'o2':
            zguess = (lamline / o2_3730_vac) - 1.0
        elif option.strip().lower() == 'hb':
            zguess = (lamline / hb_4863_vac) - 1.0
        elif option.strip().lower() == 'o31':
            zguess = (lamline / o3_4959_vac) - 1.0
        elif option.strip().lower() == 'o32':
            zguess = (lamline / o3_5007_vac) - 1.0
        elif option.strip().lower() == 'ha':
            zguess = (lamline / ha_6565_vac) - 1.0
        elif option.strip().lower() == 's2':
            zguess = (lamline / s2_6716_vac) - 1.0
        elif option.strip().lower() == 's31':
            zguess = (lamline / s3_9069_vac) - 1.0
        elif option.strip().lower() == 's32':
            zguess = (lamline / s3_9532_vac) - 1.0

        # note contamination
        elif option.strip().lower() == 'contam':
            print_prompt("Specify contamination.\nEnter a comma-separated list of identifiers choosing from:\n"+contam_flags_string)
            cf = input("> ")
            cflags = [thing.strip().lower() for thing in cf.split(',')]
            # continuum contamination sets bit 1 for all lines and the continuum itself
            if 'c' in cflags:
                for k, v in contamflags.items():
                    contamflags[k] = contamflags[k] | 2
            cflaglines = [thing for thing in cflags if thing != 'c']
            # specific line contamination sets bit 2 for all lines
            for contamflag in cflaglines:
                try:
                    contamflags[contamflag] = contamflags[contamflag] | 4
                except KeyError:
                    print_prompt('{} not known. Skipping'.format(contamflag))
            # sqlite3 database support - automatically creates and initializes DB if required
            #databaseManager.setFlags(par, obj, [(flagName, flagValue) for flagName, flagValue in contamflags.iteritems()])

        # add a comment
        elif option.strip().lower() == 'c':
            print_prompt("Enter your comment here:")
            if len(comment) > 0:
                comment = input("> ") + ', ' + comment
            else :
                comment =  input("> ")

            # sqlite3 database support - automatically creates and initializes DB if required
           # databaseManager.saveAnnotation((par, obj, comment.decode('utf-8')))

        # set or unset one or more flags
        elif option.strip().lower() == 'flag':
            print_prompt('Enter a comma-separated list of flag, value pairs e.g. CONTAM, 1, CONTIN, 2:')
            print_prompt('Valid flags are {}'.format(WISPLFDatabaseManager.WISPLFDatabaseManager.validMutableFlags))
            flagList = input("> ")
            # sqlite3 database support - automatically creates and initializes DB if required
            #databaseManager.setFlagsFromString(par, obj, flagList.decode('utf-8'))

        # write object summary
        elif option.strip().lower() == 's':
            write_object_summary(par, obj, fitresults, snr_meas_array,
                                 contamflags, summary_type='working')

        # print help message
        elif option.strip().lower() == 'h':
            print_help_message()

        ### image/display options ###
        # change 2d stamp scaling to linear
        elif option.strip().lower() == 'lin':
            if g102zeros is not None:
                show2dNEW('G102', par, obj, g102zeros, user, 'linear', path_to_wisp_data = path_to_wisp_data)
            if g141zeros is not None:
                show2dNEW('G141', par, obj, g141zeros, user, 'linear', path_to_wisp_data = path_to_wisp_data)

        # change 2d stamp scaling to log
        elif option.strip().lower() == 'log':
            if g102zeros is not None:
                show2dNEW('G102', par, obj, g102zeros, user, 'log', path_to_wisp_data = path_to_wisp_data )
            if g141zeros is not None:
                show2dNEW('G141', par, obj, g141zeros, user, 'log', path_to_wisp_data = path_to_wisp_data)

        # change g102 2d stamp scaling to zscale
        elif option.strip().lower() == 'zs102':
            print_prompt("Enter comma-separated range for G102 zscale: z1,z2")
            zscale = input("> ")
            zs = zscale.split(',')
            try:
                z1 = float(zs[0])
                z2 = float(zs[1])
            except ValueError:
                print_prompt('Invalid entry.')
            else:
                if g102zeros is not None:
                    show2dNEW('G102', par, obj, g102zeros, user, 'linear',
                              zran1=z1, zran2=z2, path_to_wisp_data = path_to_wisp_data)

        # change g141 2d stamp scaling to zscale
        elif option.strip().lower() == 'zs141':
            print_prompt("Enter comma-separated range for G141 zscale: z1,z2")
            zscale = input("> ")
            zs = zscale.split(',')
            try:
                z1 = float(zs[0])
                z2 = float(zs[1])
            except ValueError:
                print_prompt('Invalid entry.')
            else:
                if g141zeros is not None:
                    show2dNEW('G141', par, obj, g141zeros, user, 'linear',
                              zran1=z1, zran2=z2, path_to_wisp_data = path_to_wisp_Data )

        # recenter full images
        elif option.strip().lower() == 'dc':
            showDirectNEW(obj, path_to_wisp_data = path_two_wisp_data )
            if show_dispersed:  # MB
                showDispersed(obj, path_to_wisp_data = path_to_wisp_data)

        # reload full iamges
        elif option.strip().lower() == 'reload':
            showDirectNEW(obj, load_image=True, path_to_wisp_data = path_to_wisp_data)
            if show_dispersed:
                showDispersed(obj, load_image=True, path_to_wisp_data = path_to_wisp_data)

        # reload direct image region files
        elif option.strip().lower() == 'dr':
            reloadReg()

        # new options dealing with iterating objects
        # can't actually go back or choose another object now,
        # but allow for them now just in case
        elif option.strip().lower() == 'b':
            print_prompt('Please either reject or accept this object first.')
        elif 'obj' in option:
            print_prompt('Please either reject or accept this object first.')
        # print remaining objects that have not yet been inspected
        elif option.strip().lower() == 'left':
            print_prompt('Remaining objects:')
            print(remaining)
        # print all objects in line list
        elif option.strip().lower() == 'list':
            print_prompt('All objects:')
            print(allobjects)

        # quit this object
        elif option.strip().lower() == 'q':
            print_prompt('Quitting Obj %i. Nothing saved to file' % (obj))
            print_prompt('-' * 72)
            return 0

        # catch-all for everything else
        else:
            print_prompt("Invalid entry.  Try again.")

    # only re-save data if the previous fit was discarded
    if rejectPrevFit :
        # plot the whole darn thing
        plot_object(zguess, fitresults['redshift'],
                    spdata, config_pars, snr_meas_array, snr_tot_others, full_fitmodel,
                    full_contmodel, lamline, lamlines_found,
                    index_of_strongest_line, contmodel, plottitle,
                    outdir, zset=zset)

        # write to file if object was accepted
        if zset == 1:
            if np.any(spdata[1].mask):
                fitresults, snr_meas_array = check_masked_lines(fitresults, snr_meas_array, spdata)

            # write object summary
            write_object_summary(par, obj, fitresults, snr_meas_array,
                                 contamflags)

            # sqlite3 database support - automatically creates and initializes DB if required
            #databaseManager.saveCatalogueEntry(databaseManager.layoutCatalogueData(par, obj, ra[0], dec[0], a_image[0],
            #                                                                       b_image[0], jmag[0], hmag[0], fitresults, flagcont))

            writeToCatalog(linelistoutfile, par, obj, ra, dec, a_image,
                           b_image, jmag, hmag, snr_tot_others, fitresults, contamflags)

            writeFitdata(fitdatafilename, spec_lam, spec_val, spec_unc,
                         spec_con, spec_zer, full_fitmodel, full_contmodel, mask_flg)

            fitspec_pickle = open(fitdatafilename + '.pickle', 'wb')
            output_meta_data = [par, obj, ra, dec, a_image, b_image,
                                jmag, hmag, fitresults, flagcont, config_pars]
            pickle.dump(output_meta_data, fitspec_pickle)
            fitspec_pickle.close()
        # else:
        #     # done == 1, but zset == 0 => rejected
        #     databaseManager.saveCatalogueEntry(databaseManager.layoutCatalogueData(par, obj, ra[0], dec[0], a_image[0],
        #                                                                           b_image[0], jmag[0], hmag[0],
        #                                                                           None,
        #                                                                           None))
        #     databaseManager.setFlags(par, obj, [('REJECT', 1)])
        #     databaseManager.saveAnnotation((par, obj, 'REJECTED'))

        # write comments to file
        # if we go back to the previous objects, duplicate comments will still be
        # written
        writeComments(commentsfile, par, obj, comment)

        # write object to done file, incase process gets interrupted
        if not os.path.exists(outdir + '/done_%s'%user):
            f = open(outdir + '/done_%s'%user, 'w')
        else:
            f = open(outdir + '/done_%s'%user, 'a')
        f.write('%i\n' % obj)
        f.close()


def check_input_objid(objlist, objid, nextup):
    if (verbose == True): print('\nRunning check_input_objid...\n') # MDR 2022/05/17

    fulllist = ', '.join(['%i' % o for o in objlist])
    if objid not in objlist:
        print_prompt('Obj %i is not in the line list'%(objid),
                     prompt_type='interim')
        print_prompt('Full line list: \n%s'%(fulllist), prompt_type='interim')
    while objid not in objlist:
        if nextup == 0:
            o = input("Please try again, or hit enter or 'q' to quit. > ")
        else:
            o = input('Please try again, or hit enter to continue with Obj %s: > '%nextup)
        if o.strip().lower() == 'q':
            return False
        elif isFloat(o.strip()):
            objid = int(o)
        elif 'obj' in o:
            objid = int(re.search('\d+', o).group())
        elif o.strip() == '':
            return False
    return objid


def measure_z_interactive(linelistfile=" ", path_to_wisp_data = ' ', show_dispersed=True, path_to_stored_fits = ' ', print_colors=True):
    if (verbose == True): print('\nRunning measure_z_interactive...\n') # MDR 2022/05/17
    # turn off color printing to terminal if required
    if print_colors is False:
        global setcolors
        for k,v in setcolors.items():
            setcolors[k] = '\033[0m'

    if path_to_wisp_data == ' ':
        ### running from the Spectra directory
        path_to_wisp_data = '../../'


    #if path_to_stored_fits == ' ':
    #    use_stored_fits  = False
    #elif os.path.exists(path_to_stored_fits) :
    #    use_stored_fits = True
    #    print 'looking for stored fit data'
    #else:
    #    use_stored_fits = False
    #    print 'not using stored fit data'


    #### STEP 0:   set ds9 window to tile mode ################################
    ###########################################################################
    # not the best way to do this, but matching the method in guis.py
    cmd = 'xpaset -p ds9 tile'
    os.system(cmd)
    if show_dispersed:
        cmd = 'xpaset -p ds9 tile grid layout 2 3'
    else:
        cmd = 'xpaset -p ds9 tile grid layout 2 2'
    os.system(cmd)

    #### STEP 1:   get linelist ###############################################
    ###########################################################################
    if linelistfile == " ":
        files = glob('linelist/Par*lines.dat')
        if len(files) == 0:
            print_prompt('No line list file found', prompt_type='interim')
            return 0
        else:
            linelistfile = files[0]
    if not os.path.exists(linelistfile):
        print_prompt("Invalid path to line list file: %s" % (linelistfile), prompt_type='interim')
        return 0
    else:
        print_prompt('Found line list file: %s' % (linelistfile), prompt_type='interim')

    #### STEP 1b:   read the list of candidate lines  ####################
    ###########################################################################

    llin = asciitable.read(linelistfile, names=[
                           'parnos', 'grism', 'objid', 'wavelen', 'npix', 'ston'])
    parnos = llin['parnos']
    grism = llin['grism']
    objid = llin['objid']
    wavelen = llin['wavelen']
    npix = llin['npix']
    ston = llin['ston']
    objid_unique = np.unique(objid)
    par = '0' #par = parnos[0] # MDR 2022/05/17


    #### STEP 2:  set user name and output directory #########################
    ###########################################################################
    if (verbose == True):
        print('')
        print('Searching for data files in...\n') # MDR 2022/05/17
        print((os.getcwd()))
        print('')

    tmp = glob(path_to_wisp_data + 'Par' + str(par) + '/Spectra/*.dat')[0] # MDR 2022/05/17
    print_prompt('You are about to inspect emission lines identified in {}'.format(par), prompt_type='interim')
    print_prompt('Please enter your name or desired username', prompt_type='interim')
    while True:
        user = input('> ')
        if len(user) == 0:
            print_prompt('Username required', prompt_type='interim')
            continue
        else:
            break
    user = user.strip().lower()
    # create output directory
    outdir = 'Par%s_output_%s'%(par,user)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # if (verbose == True):
    #     print('tmp =', tmp) # MDR 2022/05/17
    #     print('user =', user) # MDR 2022/05/17
    #     print('outdir =', outdir) # MDR 2022/05/17

    #### STEP 3: define filenames and check for partially complete work #####
    #########################################################################
    if (verbose == True): print('\nCreating figs and fitdata directories...\n') # MDR 2022/05/17
    if not os.path.exists(os.path.join(outdir,'figs')):
        os.makedirs(os.path.join(outdir,'figs'))
    if not os.path.exists(os.path.join(outdir,'fitdata')):
        os.makedirs(os.path.join(outdir,'fitdata'))

    parts = os.path.splitext(os.path.basename(linelistfile))
    linelistoutfile = os.path.join(outdir,'%s_catalog_%s.dat'%(parts[0],user))
    commentsfile = os.path.join(outdir,'%s_comments_%s.dat'%(parts[0],user))
    # the file that will be used to determine which objects are "done"
    donefile = outdir+ '/done_%s'%user

    # if (verbose == True):
    #     print('parts =', parts) # MDR 2022/05/17
    #     print('linelistoutfile =', linelistoutfile) # MDR 2022/05/17
    #     print('commentsfile =', commentsfile) # MDR 2022/05/17
    #     print('donefile =', donefile) # MDR 2022/05/17

    if os.path.isfile(linelistoutfile):
        print_prompt('\nOutput file: \n  %s \nalready exists\n' % linelistoutfile, prompt_type='interim')
        ask = input('Append? [Y/n] ')
        if ask.lower() == 'n':
            os.unlink(linelistoutfile)
            os.unlink(commentsfile)
            # starting over, no objects have been done
            os.unlink(donefile)
            objid_done = np.array([])

            # sqlite3 database support - automatically creates and initializes DB if required
            # If the field has been previously examined, but those results are to be discarded,
            # then reset the database tables.
            # All Par numbers in the putative line list file should be the same, so the zeroth
            # element corresponds to the current field ID.
            #databaseManager = WDBM(dbFileNamePrefix=os.path.join(outdir,'Par{}'.format(parnos[0])))
            #databaseManager.resetDatabaseTables()
        else:
            # an object may be written to the comment file before it has
            # actually been inspected, so use donefile for a list
            # of the "done" objects
            objid_done = np.atleast_1d(np.genfromtxt(donefile, dtype=int))
    else:
        if os.path.exists(donefile):
            os.unlink(donefile)
        objid_done = np.array([], dtype=int)

    #### STEP 4: create trace.reg files ############################
    #########################################################################
    if (verbose == True): print('Creating trace.reg files...\n') # MDR 2022/05/17

    trace102 = open(path_to_wisp_data + '/Par'  + str(par) + '/Spectra/G102_trace.reg', 'w')
    trace102.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace102.write('wcs;\n')
    # sensitivity drops below 25% of max at wave < 8250 and wave > 11540
    # so box should be 3290 angstroms wide and be centered at 9895.
    trace102.write('box(9895,0,3290,1,1.62844e-12)\n')
    trace102.close()
    trace141 = open(path_to_wisp_data + '/Par' + str(par) + '/Spectra/G141_trace.reg', 'w')
    trace141.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    trace141.write('wcs;\n')
    # sensitivity drops below 25% of max at wave < 10917 and wave > 16904
    # so box should be 5897 angstroms wide and be centered at 13910.5
    trace141.write('box(13910.5,0,5897,1,0)\n')
    trace141.close()

    #### STEP 5:  Get zero and first order positions; unpack them ###########
    #########################################################################
    g102zeroordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G102_0th.reg'
    g102firstordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G102_1st.reg'
    g141zeroordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G141_0th.reg'
    g141firstordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G141_1st.reg'

    #### STEP 6:  Get object information from SExtractor catalog ############
    #########################################################################
    # a_image will give an order of magnitude estimate on the FWHM of the line
    #   this determines the initial guess and sets an upper limit on how broad
    #   it can be.
    # ra/dec, b_image, jmag, jerr, hmag, herr will be carried forward into
    #   the output linelist.
    # find all available cats
    secats = glob(os.getcwd()+'/DATA/DIRECT_GRISM/fin_f*.cat') # MDR 2022/05/17
    secats.sort()
    cat = asciitable.read(secats[0])

    if (verbose == True):
        print('I found the following catalogs...\n') # MDR 2022/05/17
        print(secats) # MDR 2022/05/17
        print('\nThe catalog was read in as...\n') # MDR 2022/05/17
        print(cat)
        print('') # MDR 2022/05/17

    #beam = cat['col2'] - M.D.R. - 10/08/2020
    #a_image = cat['col5'] - M.D.R. - 10/08/2020
    #b_image = cat['col6'] - M.D.R. - 10/08/2020
    #ra = cat['col8'] - M.D.R. - 10/08/2020
    #dec = cat['col9'] - M.D.R. - 10/08/2020

    # Edited to call specific columns. - M.D.R. - 10/08/2020
    beam = cat['NUMBER']
    a_image = cat['A_IMAGE']
    b_image = cat['B_IMAGE']
    ra = cat['ALPHA_J2000']
    dec = cat['DELTA_J2000']

    # which filter is this?
    # if os.path.basename(secats[0]) == 'fin_F110.cat': # MDR 2022/05/17
    #     jmag = cat['col13'] # MDR 2022/05/17
    #     jerr = cat['col14'] # MDR 2022/05/17
    #     # in case there is a F110-only field # MDR 2022/05/17
    #     hmag = np.ones(ra.shape, dtype=float) * 99. # MDR 2022/05/17
    #     herr = np.ones(ra.shape, dtype=float) * 99. # MDR 2022/05/17
    # else:
    jmag = np.ones(ra.shape, dtype=float) * 99.
    jerr = np.ones(ra.shape, dtype=float) * 99.
    hmag = cat['MAG_AUTO'] # hmag = cat['col13'] - M.D.R. - 10/08/2020
    herr = cat['MAGERR_AUTO'] # herr = cat['col14'] - M.D.R. - 10/08/2020
    # read in second file if there are two
    # if len(secats) == 2: # MDR 2022/05/17
    #     cat2 = asciitable.read(secats[1]) # MDR 2022/05/17
    #     # second catalog should be hband # MDR 2022/05/17
    #     hmag = cat2['col13'] # MDR 2022/05/17
    #     herr = cat2['col14'] # MDR 2022/05/17
    objtable = Table([beam, ra, dec, a_image, b_image, jmag, jerr, hmag, herr],
                     names=('obj', 'ra', 'dec', 'a_image', 'b_image', 'jmag',
                            'jerr', 'hmag', 'herr'))

    #### STEP 7:  Set up initial ds9 display ################################
    #########################################################################
    if os.path.exists(g102zeroordreg):
        g102zeroarr = getzeroorders(g102zeroordreg, g='G102')
        # nothing is done with the first orders anymore
        # g102firstarr=getfirstorders(g102firstordreg)
        show2dNEW('G102', parnos[0], objid_unique[0], g102zeroarr, user, 'linear', path_to_wisp_data = path_to_wisp_data)
    else:
        g102zeroarr = None
        g102firstarr = None
    if os.path.exists(g141zeroordreg):
        g141zeroarr = getzeroorders(g141zeroordreg, g='G102')
        # g141firstarr=getfirstorders(g141firstordreg)
        show2dNEW('G141', parnos[0], objid_unique[0], g141zeroarr, user, 'linear', path_to_wisp_data = path_to_wisp_data)
    else:
        g141zeroarr = None
        g141firstarr = None

    showDirectNEW(objid_unique[0], parnos[0], load_image=True, path_to_wisp_data = path_to_wisp_data)
    if show_dispersed:  # MB
        showDispersed(objid_unique[0], parnos[0], load_image=True, path_to_wisp_data  = path_to_wisp_data)

    #### STEP 8:  Loop through objects ############
    #########################################################################
    if (verbose == True): print('\nStarting loop through objects...\n') # MDR 2022/05/17

    remaining_objects = get_remaining_objects(objid_unique, objid_done)
    allobjects = [unique_obj for unique_obj in objid_unique]

    # if (verbose == True):
    #     print('remaining_objects =', remaining_objects) # MDR 2022/05/17
    #     print('allobjects =', allobjects) # MDR 2022/05/17

    print_prompt('\nAs you loop through the objects, you can choose from the following\noptions at any time:\n\txxx = skip to object xxx\n\tb = revisit the previous object\n\tleft = list all remaining objects that need review\n\tlist = list all objects in line list\n\tany other key = continue with the next object\n\th = help/list interactive commands\n\tq = quit\n', prompt_type='interim')

    while remaining_objects.shape[0] > 0:
        if path_to_stored_fits == ' ':
            use_stored_fits  = False
        elif os.path.exists(path_to_stored_fits) :
            use_stored_fits = True
            print('looking for stored fit data')
        else:
            use_stored_fits = False
            print('not using stored fit data')

        ndone = len(np.unique(objid_done))
        progress = float(ndone) / float(len(objid_unique)) * 100.
        print_prompt("\nProgress: %.1f percent" % (progress),prompt_type='interim')

        # do some things as long as there are still objects to inspect
        next_obj = remaining_objects[0]
        print_prompt('Next up: Obj %i' % (next_obj),prompt_type='interim')
        o = input(
            "Enter 'xxx' to skip to Obj xxx or hit any key to continue. > ")


        if o.strip().lower() == 'left':
            #remaining_list = ', '.join(['%i'%i for i in remaining_objects])
            print_prompt('Remaining objects:', prompt_type='interim')
            print(remaining_objects)
            o = input('> ')

        if o.strip().lower() == 'list':
            print_prompt('All objects:', prompt_type='interim')
            print(allobjects)
            o = input('> ')

        if o.strip().lower() == 'b':
            previous_obj = int(objid_done[-1])
            # need to figure out what object came before this one
            #w = np.where(objid_unique == remaining_objects[0])
            # if on first object, this will roll around to previous object
            next_obj = previous_obj
            print_prompt("Going back to previous object: Obj %i" % (next_obj),prompt_type='interim')

        elif o.strip().lower() == 'q':
            print_prompt("Quitting; saved through previously completed object.",prompt_type='interim')
            return 0

        elif o.strip().lower() == 'random':
            print_prompt('Hi Vihang!', prompt_type='random')
            next_obj_idx = np.random.randint(remaining_objects.shape[0])
            next_obj = remaining_objects[next_obj_idx]

        elif isFloat(o.strip()):
            next_obj = int(re.search('\d+', o).group())
            # confirm that requested object is in line list
            next_obj = check_input_objid(objid_unique, next_obj, remaining_objects[0])
            next_obj = next_obj if next_obj else remaining_objects[0]
        elif 'obj' in o:
            next_obj = int(re.search('\d+', o).group())
            # confirm that requested object is in line list
            next_obj = check_input_objid(objid_unique, next_obj, remaining_objects[0])
            next_obj = next_obj if next_obj else remaining_objects[0]


        # pass the information for this object
        wlinelist = np.where(objid == next_obj)
        lamlines_found = wavelen[wlinelist]
        ston_found = ston[wlinelist]
        wcatalog = np.where(objtable['obj'] == next_obj)
        objinfo = objtable[wcatalog]
        #inspect_object(user, parnos[0], next_obj, objinfo, lamlines_found,
        #               ston_found, g102zeroarr, g141zeroarr, linelistoutfile,
        #               commentsfile, remaining_objects, allobjects,
        #               show_dispersed=show_dispersed)


        if (use_stored_fits == True):
            ### get pickle files:
            inpickles = []
            path_pickle1 = path_to_stored_fits + '/Par0_output_a/fitdata/Par0_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle1 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_mbagley/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle2 = path_to_stored_fits + '/Par'  + str(parnos[0]) +    '_output_marc/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle3 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_ben/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle4 = path_to_stored_fits + '/Par'  + str(parnos[0]) +     '_output_claudia/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle5 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_vihang/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle6 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_ivano/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle7 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_mbeck/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle8 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_karlenoid/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle9 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_mjr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle10 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_sophia/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle11 = '/Volumes/Thunderbay/wisps/mzr_refit/Par'  + str(parnos[0]) + '_output_marc-mzr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle12 = '/Volumes/Thunderbay/wisps/mzr_refit/Par'  + str(parnos[0]) + '_output_alaina-mzr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'



            ### put new fits first
            # if os.path.exists(path_pickle11):
            #     inpickles.append(path_pickle11)
            # if os.path.exists(path_pickle12):
            #     inpickles.append(path_pickle12)
            if os.path.exists(path_pickle1):
                inpickles.append(path_pickle1)
            # if os.path.exists(path_pickle2):
            #     inpickles.append(path_pickle2)
            # if os.path.exists(path_pickle3):
            #     inpickles.append(path_pickle3)
            # if os.path.exists(path_pickle4):
            #     inpickles.append(path_pickle4)
            # if os.path.exists(path_pickle5):
            #     inpickles.append(path_pickle5)
            # if os.path.exists(path_pickle6):
            #     inpickles.append(path_pickle6)
            # if os.path.exists(path_pickle7):
            #     inpickles.append(path_pickle7)
            # if os.path.exists(path_pickle8):
            #     inpickles.append(path_pickle8)
            # if os.path.exists(path_pickle9):
            #     inpickles.append(path_pickle9)
            # if os.path.exists(path_pickle10):
            #     inpickles.append(path_pickle10)


            if len(inpickles) == 0:
                use_stored_fits = False


        if use_stored_fits == True:
             inspect_object(user, parnos[0], next_obj, objinfo,
                                lamlines_found, ston_found, g102zeroarr,
                                g141zeroarr, linelistoutfile, commentsfile,
                                remaining_objects, allobjects,
                                 show_dispersed=show_dispersed, stored_fits = inpickles, path_to_wisp_data = path_to_wisp_data)
        else:
            inspect_object(user, parnos[0], next_obj, objinfo,
                            lamlines_found, ston_found, g102zeroarr,
                            g141zeroarr, linelistoutfile, commentsfile,
                            remaining_objects, allobjects,
                            show_dispersed=show_dispersed, stored_fits = False, path_to_wisp_data = path_to_wisp_data)

        objid_done = np.append(objid_done, next_obj)
        remaining_objects = get_remaining_objects(objid_unique, objid_done)

    # outside the while loop, field is done
    redo = ' '
    while redo != 'q':
        print_prompt("You've finished this field.\nEnter an object ID below to revisit a particular object.\nOtherwise enter 'q' to quit the field.", prompt_type='interim')
        redo = input("> ").strip().lower()
        if redo != 'q':
            try:
                next_obj = int(re.search('\d+', redo).group())
            except (ValueError,AttributeError):
                print_prompt("Invalid entry. Enter an object ID or enter 'q' to quit", prompt_type='interim')
            else:
                next_obj = check_input_objid(objid_unique, next_obj, 0)
                if next_obj:
                    # pass the information for this object
                    wlinelist = np.where(objid == next_obj)
                    lamlines_found = wavelen[wlinelist]
                    ston_found = ston[wlinelist]
                    wcatalog = np.where(objtable['obj'] == next_obj)
                    objinfo = objtable[wcatalog]


                    if (use_stored_fits == True):
                       ### get pickle files:
                       inpickles = []
                       path_pickle1 = path_to_stored_fits + '/Par0_output_a/fitdata/Par0_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle1 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_mbagley/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle2 = path_to_stored_fits + '/Par'  + str(parnos[0]) +    '_output_marc/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle3 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_claudia/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle4 = path_to_stored_fits + '/Par'  + str(parnos[0]) +     '_output_ben/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle5 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_vihang/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle6 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_ivano/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle7 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_mbeck/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle8 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_karlenoid/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle9 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_mjr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle10 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_sophia/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle11 = '/Volumes/Thunderbay/wisps/mzr_refit/Par'  + str(parnos[0]) + '_output_marc-mzr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
                       # path_pickle12 = '/Volumes/Thunderbay/wisps/mzr_refit/Par'  + str(parnos[0]) + '_output_alaina-mzr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'


                       ### put new fits first
                       # if os.path.exists(path_pickle11):
                       #     inpickles.append(path_pickle11)
                       # if os.path.exists(path_pickle12):
                       #     inpickles.append(path_pickle12)
                       if os.path.exists(path_pickle1):
                          inpickles.append(path_pickle1)
                       # if os.path.exists(path_pickle2):
                       #    inpickles.append(path_pickle2)
                       # if os.path.exists(path_pickle3):
                       #    inpickles.append(path_pickle3)
                       # if os.path.exists(path_pickle4):
                       #    inpickles.append(path_pickle4)
                       # if os.path.exists(path_pickle5):
                       #    inpickles.append(path_pickle5)
                       # if os.path.exists(path_pickle6):
                       #    inpickles.append(path_pickle6)
                       # if os.path.exists(path_pickle7):
                       #    inpickles.append(path_pickle7)
                       # if os.path.exists(path_pickle8):
                       #    inpickles.append(path_pickle8)
                       # if os.path.exists(path_pickle9):
                       #    inpickles.append(path_pickle9)
                       # if os.path.exists(path_pickle10):
                       #    inpickles.append(path_pickle10)


                       if len(inpickles) == 0:
                           use_stored_fits = False


                    if (use_stored_fits ==True):
                       inspect_object(user, parnos[0], next_obj, objinfo,
                                lamlines_found, ston_found, g102zeroarr,
                                g141zeroarr, linelistoutfile, commentsfile,
                                remaining_objects, allobjects,
                                 show_dispersed=show_dispersed, stored_fits = inpickles, path_to_wisp_data = path_to_wisp_data)
                    else:
                        inspect_object(user, parnos[0], next_obj, objinfo,
                            lamlines_found, ston_found, g102zeroarr,
                            g141zeroarr, linelistoutfile, commentsfile,
                            remaining_objects, allobjects,
                            show_dispersed=show_dispersed, stored_fits = False, path_to_wisp_data = path_to_wisp_data)
                else:
                    break

    make_tarfile(outdir)
    print_prompt('A tarfile of your outputs has been created: %s.tar.gz'%outdir, prompt_type='interim')


   # Clean up temp files
    if os.path.exists('./tempcoo.dat') == 1:
        os.unlink('./tempcoo.dat')
    if os.path.exists('./temp_zero_coords.coo') == 1:
        os.unlink('./temp_zero_coords.coo')
    if os.path.exists('./temp110.fits') == 1:
        os.unlink('./temp110.fits')
    if os.path.exists('./temp140.fits') == 1:
        os.unlink('./temp140.fits')
    if os.path.exists('./temp160.fits') == 1:
        os.unlink('./temp160.fits')
    if os.path.exists('./temp_zero_coords.reg') == 1:
        os.unlink('./temp_zero_coords.reg')
    if os.path.exists('G102_trace.reg') == True:
        os.unlink('G102_trace.reg')
    if os.path.exists('G141_trace.reg') == True:
        os.unlink('G141_trace.reg')


# parnos, objid are scalar not array.
def writeToCatalog(catalogname, parnos, objid, ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, snr_tot_others, fitresults, contamflags):
    if not os.path.exists(catalogname):
        cat = open(catalogname, 'w')
        cat.write('#1 objid \n')
        cat.write('#2 redshift \n')
        cat.write('#3 redshift_error \n')
        cat.write('#4 ra_obj \n')
        cat.write('#5 dec_obj \n')
        cat.write('#6 f140w_mag \n')
        cat.write('#7 a_image_obj \n')
        cat.write('#8 b_image_obj \n')
        cat.write('#9 snr_tot_others \n')
        cat.write('#10 chisq \n')
        cat.write('#11 fwhm_muse \n')
        cat.write('#12 fwhm_muse_error \n')
        cat.write('#13 fwhm_g141 \n')
        cat.write('#14 fwhm_g141_error \n')
        cat.write('#15 la_1216_dz \n')
        cat.write('#16 c4_1548_dz \n')
        cat.write('#17 uv_line_dz \n')
        cat.write('#18 m2_2796_dz \n')
        cat.write('#19 o2_3727_dz \n')
        cat.write('#20 o3_5007_dz \n')
        cat.write('#21 s3_he_dz \n')

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

        results_idx = 22

        for line in result_lines:
            cat.write('#'+str(results_idx+0)+' '+line+'_flux \n')
            cat.write('#'+str(results_idx+1)+' '+line+'_error \n')
            cat.write('#'+str(results_idx+2)+' '+line+'_ew_obs \n')
            cat.write('#'+str(results_idx+3)+' '+line+'_contam \n')
            results_idx = results_idx + 4

        cat.close()
    # does not leave space before RA?
    outstr = '{:<6d}'.format(objid) +\
        '{:>8.5f}'.format(fitresults['redshift']) + \
        '{:>8.5f}'.format(fitresults['redshift_error']) +\
        '{:>12.6f}'.format(ra_obj[0]) + \
        '{:>12.6f}'.format(dec_obj[0]) + \
        '{:>8.2f}'.format(hmag_obj[0]) + \
        '{:>8.3f}'.format(a_image_obj[0]) + \
        '{:>8.3f}'.format(b_image_obj[0]) + \
        '{:>10.2f}'.format(snr_tot_others) +\
        '{:>10.2f}'.format(fitresults['chisq']) + \
        '{:>10.2f}'.format(fitresults['fwhm_muse']) + \
        '{:>10.2f}'.format(fitresults['fwhm_muse_error'])  +  \
        '{:>10.2f}'.format(fitresults['fwhm_g141']) + \
        '{:>10.2f}'.format(fitresults['fwhm_g141_error'])  +  \
        '{:>10.5f}'.format(fitresults['la_1216_dz'])   + \
        '{:>10.5f}'.format(fitresults['c4_1548_dz'])   + \
        '{:>10.5f}'.format(fitresults['uv_line_dz'])   + \
        '{:>10.5f}'.format(fitresults['m2_2796_dz'])   + \
        '{:>10.5f}'.format(fitresults['o2_3727_dz'])   + \
        '{:>10.5f}'.format(fitresults['o3_5007_dz'])  + \
        '{:>10.5f}'.format(fitresults['s3_he_dz']) +\
        '{:>13.2e}'.format(fitresults['la_1216_flux']) + \
        '{:>13.2e}'.format(fitresults['la_1216_error']) + \
        '{:>13.2e}'.format(fitresults['la_1216_ew_obs']) +\
        '{:>6d}'.format(contamflags['la_1216']) +\
        '{:>13.2e}'.format(fitresults['la_wing_flux']) + \
        '{:>13.2e}'.format(fitresults['la_wing_error']) + \
        '{:>13.2e}'.format(fitresults['la_wing_ew_obs']) +\
        '{:>6d}'.format(contamflags['la_1216']) +\
        '{:>13.2e}'.format(fitresults['la_1216_wing_flux']) + \
        '{:>13.2e}'.format(fitresults['la_1216_wing_error']) + \
        '{:>13.2e}'.format(fitresults['la_1216_wing_ew_obs']) +\
        '{:>6d}'.format(contamflags['la_1216']) +\
        '{:>13.2e}'.format(fitresults['n5_1238_flux'])  + \
        '{:>13.2e}'.format(fitresults['n5_1238_error']) +  \
        '{:>13.2e}'.format(fitresults['n5_1238_ew_obs']) +\
        '{:>6d}'.format(contamflags['n5_1238']) +\
        '{:>13.2e}'.format(fitresults['n5_1242_flux'])  + \
        '{:>13.2e}'.format(fitresults['n5_1242_error']) +  \
        '{:>13.2e}'.format(fitresults['n5_1242_ew_obs']) +\
        '{:>6d}'.format(contamflags['n5_1242']) +\
        '{:>13.2e}'.format(fitresults['n5_1238_1242_flux'])  + \
        '{:>13.2e}'.format(fitresults['n5_1238_1242_error']) +  \
        '{:>13.2e}'.format(fitresults['n5_1238_1242_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['n5_1238'], contamflags['n5_1242']])) +\
        '{:>13.2e}'.format(fitresults['c4_1548_flux'])  + \
        '{:>13.2e}'.format(fitresults['c4_1548_error']) +  \
        '{:>13.2e}'.format(fitresults['c4_1548_ew_obs']) +\
        '{:>6d}'.format(contamflags['c4_1548']) +\
        '{:>13.2e}'.format(fitresults['c4_1550_flux'])  + \
        '{:>13.2e}'.format(fitresults['c4_1550_error']) +  \
        '{:>13.2e}'.format(fitresults['c4_1550_ew_obs']) +\
        '{:>6d}'.format(contamflags['c4_1550']) +\
        '{:>13.2e}'.format(fitresults['c4_1548_1550_flux'])  + \
        '{:>13.2e}'.format(fitresults['c4_1548_1550_error']) +  \
        '{:>13.2e}'.format(fitresults['c4_1548_1550_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['c4_1548'], contamflags['c4_1550']])) +\
        '{:>13.2e}'.format(fitresults['h2_1640_flux'])  + \
        '{:>13.2e}'.format(fitresults['h2_1640_error']) +  \
        '{:>13.2e}'.format(fitresults['h2_1640_ew_obs']) +\
        '{:>6d}'.format(contamflags['h2_1640']) +\
        '{:>13.2e}'.format(fitresults['o3_1660_flux'])  + \
        '{:>13.2e}'.format(fitresults['o3_1660_error']) +  \
        '{:>13.2e}'.format(fitresults['o3_1660_ew_obs']) +\
        '{:>6d}'.format(contamflags['o3_1660']) +\
        '{:>13.2e}'.format(fitresults['o3_1666_flux'])  + \
        '{:>13.2e}'.format(fitresults['o3_1666_error']) +  \
        '{:>13.2e}'.format(fitresults['o3_1666_ew_obs']) +\
        '{:>6d}'.format(contamflags['o3_1666']) +\
        '{:>13.2e}'.format(fitresults['o3_1660_1666_flux'])  + \
        '{:>13.2e}'.format(fitresults['o3_1660_1666_error']) +  \
        '{:>13.2e}'.format(fitresults['o3_1660_1666_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['o3_1660'], contamflags['o3_1666']])) +\
        '{:>13.2e}'.format(fitresults['s3_1883_flux'])  + \
        '{:>13.2e}'.format(fitresults['s3_1883_error']) +  \
        '{:>13.2e}'.format(fitresults['s3_1883_ew_obs']) +\
        '{:>6d}'.format(contamflags['s3_1883']) +\
        '{:>13.2e}'.format(fitresults['s3_1892_flux'])  + \
        '{:>13.2e}'.format(fitresults['s3_1892_error']) +  \
        '{:>13.2e}'.format(fitresults['s3_1892_ew_obs']) +\
        '{:>6d}'.format(contamflags['s3_1892']) +\
        '{:>13.2e}'.format(fitresults['s3_1883_1892_flux'])  + \
        '{:>13.2e}'.format(fitresults['s3_1883_1892_error']) +  \
        '{:>13.2e}'.format(fitresults['s3_1883_1892_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['s3_1883'], contamflags['s3_1892']])) +\
        '{:>13.2e}'.format(fitresults['c3_1907_flux'])  + \
        '{:>13.2e}'.format(fitresults['c3_1907_error']) +  \
        '{:>13.2e}'.format(fitresults['c3_1907_ew_obs']) +\
        '{:>6d}'.format(contamflags['c3_1907']) +\
        '{:>13.2e}'.format(fitresults['c3_1909_flux'])  + \
        '{:>13.2e}'.format(fitresults['c3_1909_error']) +  \
        '{:>13.2e}'.format(fitresults['c3_1909_ew_obs']) +\
        '{:>6d}'.format(contamflags['c3_1909']) +\
        '{:>13.2e}'.format(fitresults['c3_1907_1909_flux'])  + \
        '{:>13.2e}'.format(fitresults['c3_1907_1909_error']) +  \
        '{:>13.2e}'.format(fitresults['c3_1907_1909_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['c3_1907'], contamflags['c3_1909']])) +\
        '{:>13.2e}'.format(fitresults['m2_2796_flux'])  + \
        '{:>13.2e}'.format(fitresults['m2_2796_error']) +  \
        '{:>13.2e}'.format(fitresults['m2_2796_ew_obs']) +\
        '{:>6d}'.format(contamflags['m2_2796']) +\
        '{:>13.2e}'.format(fitresults['m2_2803_flux'])  + \
        '{:>13.2e}'.format(fitresults['m2_2803_error']) +  \
        '{:>13.2e}'.format(fitresults['m2_2803_ew_obs']) +\
        '{:>6d}'.format(contamflags['m2_2803']) +\
        '{:>13.2e}'.format(fitresults['m2_2796_2803_flux'])  + \
        '{:>13.2e}'.format(fitresults['m2_2796_2803_error']) +  \
        '{:>13.2e}'.format(fitresults['m2_2796_2803_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['m2_2796'], contamflags['m2_2803']])) +\
        '{:>13.2e}'.format(fitresults['o2_3727_flux'])  + \
        '{:>13.2e}'.format(fitresults['o2_3727_error']) +  \
        '{:>13.2e}'.format(fitresults['o2_3727_ew_obs']) +\
        '{:>6d}'.format(contamflags['o2_3727']) +\
        '{:>13.2e}'.format(fitresults['o2_3730_flux'])  + \
        '{:>13.2e}'.format(fitresults['o2_3730_error']) +  \
        '{:>13.2e}'.format(fitresults['o2_3730_ew_obs']) +\
        '{:>6d}'.format(contamflags['o2_3730']) +\
        '{:>13.2e}'.format(fitresults['o2_3727_3730_flux'])  + \
        '{:>13.2e}'.format(fitresults['o2_3727_3730_error']) +  \
        '{:>13.2e}'.format(fitresults['o2_3727_3730_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['o2_3727'], contamflags['o2_3730']])) +\
        '{:>13.2e}'.format(fitresults['hg_4342_flux']) +\
        '{:>13.2e}'.format(fitresults['hg_4342_error']) + \
        '{:>13.2e}'.format(fitresults['hg_4342_ew_obs']) +\
        '{:>6d}'.format(contamflags['hg_4342']) +\
        '{:>13.2e}'.format(fitresults['o3_4363_flux']) +\
        '{:>13.2e}'.format(fitresults['o3_4363_error']) + \
        '{:>13.2e}'.format(fitresults['o3_4363_ew_obs']) +\
        '{:>6d}'.format(contamflags['o3_4363']) +\
        '{:>13.2e}'.format(fitresults['h2_4686_flux']) +\
        '{:>13.2e}'.format(fitresults['h2_4686_error']) + \
        '{:>13.2e}'.format(fitresults['h2_4686_ew_obs']) +\
        '{:>6d}'.format(contamflags['h2_4686']) +\
        '{:>13.2e}'.format(fitresults['hb_4863_flux']) + \
        '{:>13.2e}'.format(fitresults['hb_4863_error']) + \
        '{:>13.2e}'.format(fitresults['hb_4863_ew_obs']) +\
        '{:>6d}'.format(contamflags['hb_4863']) +\
        '{:>13.2e}'.format(fitresults['o3_4959_flux']) + \
        '{:>13.2e}'.format(fitresults['o3_4959_error']) + \
        '{:>13.2e}'.format(fitresults['o3_4959_ew_obs']) +\
        '{:>6d}'.format(contamflags['o3_4959']) +\
        '{:>13.2e}'.format(fitresults['o3_5007_flux']) + \
        '{:>13.2e}'.format(fitresults['o3_5007_error']) + \
        '{:>13.2e}'.format(fitresults['o3_5007_ew_obs']) +\
        '{:>6d}'.format(contamflags['o3_5007']) +\
        '{:>13.2e}'.format(fitresults['o3_4959_5007_flux'])  + \
        '{:>13.2e}'.format(fitresults['o3_4959_5007_error']) +  \
        '{:>13.2e}'.format(fitresults['o3_4959_5007_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['o3_4959'], contamflags['o3_5007']])) +\
        '{:>13.2e}'.format(fitresults['o1_6300_flux'])  + \
        '{:>13.2e}'.format(fitresults['o1_6300_error']) +  \
        '{:>13.2e}'.format(fitresults['o1_6300_ew_obs']) +\
        '{:>6d}'.format(contamflags['o1_6300']) +\
        '{:>13.2e}'.format(fitresults['o1_6363_flux'])  + \
        '{:>13.2e}'.format(fitresults['o1_6363_error']) +  \
        '{:>13.2e}'.format(fitresults['o1_6363_ew_obs']) +\
        '{:>6d}'.format(contamflags['o1_6363']) +\
        '{:>13.2e}'.format(fitresults['o1_6300_6363_flux'])  + \
        '{:>13.2e}'.format(fitresults['o1_6300_6363_error']) +  \
        '{:>13.2e}'.format(fitresults['o1_6300_6363_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['o1_6300'], contamflags['o1_6363']])) +\
        '{:>13.2e}'.format(fitresults['n2_6550_flux'])  + \
        '{:>13.2e}'.format(fitresults['n2_6550_error']) +  \
        '{:>13.2e}'.format(fitresults['n2_6550_ew_obs']) +\
        '{:>6d}'.format(contamflags['n2_6550']) +\
        '{:>13.2e}'.format(fitresults['ha_6565_flux'])  + \
        '{:>13.2e}'.format(fitresults['ha_6565_error']) +  \
        '{:>13.2e}'.format(fitresults['ha_6565_ew_obs']) +\
        '{:>6d}'.format(contamflags['ha_6565']) +\
        '{:>13.2e}'.format(fitresults['n2_6585_flux'])  + \
        '{:>13.2e}'.format(fitresults['n2_6585_error']) +  \
        '{:>13.2e}'.format(fitresults['n2_6585_ew_obs']) +\
        '{:>6d}'.format(contamflags['n2_6585']) +\
        '{:>13.2e}'.format(fitresults['ha_6550_6565_6585_flux']) + \
        '{:>13.2e}'.format(fitresults['ha_6550_6565_6585_error']) + \
        '{:>13.2e}'.format(fitresults['ha_6550_6565_6585_ew_obs']) + \
        '{:>6d}'.format(np.max([contamflags['n2_6550'], contamflags['ha_6565'], contamflags['n2_6585']])) +\
        '{:>13.2e}'.format(fitresults['s2_6716_flux'])  + \
        '{:>13.2e}'.format(fitresults['s2_6716_error']) +  \
        '{:>13.2e}'.format(fitresults['s2_6716_ew_obs']) +\
        '{:>6d}'.format(contamflags['s2_6716']) +\
        '{:>13.2e}'.format(fitresults['s2_6731_flux'])  + \
        '{:>13.2e}'.format(fitresults['s2_6731_error']) +  \
        '{:>13.2e}'.format(fitresults['s2_6731_ew_obs']) +\
        '{:>6d}'.format(contamflags['s2_6731']) +\
        '{:>13.2e}'.format(fitresults['s2_6716_6731_flux']) + \
        '{:>13.2e}'.format(fitresults['s2_6716_6731_error']) + \
        '{:>13.2e}'.format(fitresults['s2_6716_6731_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['s2_6716'], contamflags['s2_6731']])) +\
        '{:>13.2e}'.format(fitresults['s3_9069_flux']) + \
        '{:>13.2e}'.format(fitresults['s3_9069_error']) + \
        '{:>13.2e}'.format(fitresults['s3_9069_ew_obs']) +\
        '{:>6d}'.format(contamflags['s3_9069']) +\
        '{:>13.2e}'.format(fitresults['s3_9532_flux']) + \
        '{:>13.2e}'.format(fitresults['s3_9532_error']) + \
        '{:>13.2e}'.format(fitresults['s3_9532_ew_obs']) +\
        '{:>6d}'.format(contamflags['s3_9532']) +\
        '{:>13.2e}'.format(fitresults['s3_9069_9532_flux']) + \
        '{:>13.2e}'.format(fitresults['s3_9069_9532_error']) + \
        '{:>13.2e}'.format(fitresults['s3_9069_9532_ew_obs']) +\
        '{:>6d}'.format(np.max([contamflags['s3_9069'], contamflags['s3_9532']])) +\
        '{:>13.2e}'.format(fitresults['he10830_flux']) + \
        '{:>13.2e}'.format(fitresults['he10830_error']) + \
        '{:>13.2e}'.format(fitresults['he10830_ew_obs']) +\
        '{:>6d}'.format(contamflags['he10830']) + '\n'

    """
    # if a row already exists for this object, comment it out
    objstr = '{:<8d}'.format(parnos) + '{:<6d}'.format(objid)
    for line in fileinput.input(catalogname, inplace=True):
        if objstr in line:
            print "#%s" % line,
        else:
            print '%s' % line,
#    """

    cat = open(catalogname, 'a')
    cat.write(outstr)
    cat.close()


def writeFitdata(filename, lam, flux, eflux, contam, zero, fit, continuum, masks):
    if (verbose == True): print('\nRunning writeFitdata...\n') # MDR 2022/05/17
    """ """
    fitspec_file = filename + '.dat'
    t = Table([lam, flux, eflux, contam, zero, fit, continuum, masks],
              names=('Lam', 'Flam', 'Flam_err', 'Contam', 'Zero', 'Fitmodel', 'Contmodel', 'Masked'))
    t['Lam'].format = '{:.1f}'
    t['Flam'].format = '{:.5e}'
    t['Flam_err'].format = '{:.5e}'
    t['Contam'].format = '{:.5e}'
    t['Zero'].format = '{:.0f}'
    t['Fitmodel'].format = '{:.5e}'
    t['Contmodel'].format = '{:.5e}'
    t['Masked'].format = '{:.0f}'
    asciitable.write(t, fitspec_file, fill_values=[(asciitable.masked, '--')],
                     overwrite=True, format='fixed_width_two_line',
                     position_char='#', delimiter_pad=' ')


def writeComments(filename, parnos, objid, comment):
    if (verbose == True): print('Running writeComments...\n') # MDR 2022/05/17
    if os.path.exists(filename) == False:
        cat = open(filename, 'w')
    else:
        cat = open(filename, 'a')

    outstr = '{:<8d}'.format(parnos) + \
        '{:<6d}'.format(objid) +\
        comment + '\n'

    cat.write(outstr)
    cat.close()


def UpdateCatalog(linelistoutfile):
    if (verbose == True): print('Running UpdateCatalog...\n') # MDR 2022/05/17
    allDirectoryFiles = os.listdir('./fitdata/')
    objid_list = []
    for obj in allDirectoryFiles:
        x = obj.split('_')[2]
        objid_list.append(int(x))
        Parno = obj.split('_')[0]   # this is inefficient, but I don't care.
    objid_list = np.sort(np.unique(np.array(objid_list)))
    for obj in objid_list:
        print(obj)
        inpickle = './fitdata/' + Parno + \
            '_BEAM_' + str(obj) + '_fitspec.pickle'
        fileObject = open(inpickle, 'r')
        alldata = pickle.load(fileObject)
        # definition from above
       #                      0          1                 2      3        4           5            6         7         8           9         10
       # output_meta_data = [parnos[0], objid_unique[i], ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont, config_pars]
        parnos = alldata[0]
        objid_unique = alldata[1]
        ra_obj = alldata[2]
        dec_obj = alldata[3]
        a_image_obj = alldata[4]
        b_image_obj = alldata[5]
        jmag_obj = alldata[6]
        hmag_obj = alldata[7]
        fitresults = alldata[8]
        flagcont = alldata[9]
        # config_pars = alldata[10] ## not used here.

        if (verbose == True): print('Writing to catalog...\n') # MDR 2022/05/17
        WriteToCatalog(linelistoutfile, parnos, objid_unique, ra_obj, dec_obj,
                       a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont)
