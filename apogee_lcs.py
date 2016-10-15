import matplotlib
matplotlib.use('agg')
from keblat import *
import matplotlib.pyplot as plt
import sys, itertools
import emcee
import lmfit
if lmfit.__version__[2] != '9':
    print "Version >= 0.9 of lmfit required..."
    sys.exit()
from lmfit import minimize, Parameters, report_fit
import scipy.optimize
from helper_funcs import *

from mpl_toolkits.axes_grid1 import make_axes_locatable

#kic = int(float(sys.argv[1]))
#period, tpe, esinw, ecosw, rsum, rrat, b, frat, q1, q2, q3, q4 = np.array(sys.argv[2:-2], dtype=float)
#nwalkers, niter = int(sys.argv[-2]), int(sys.argv[-1])

clobber_lc=False #overwrite LC only fits?
clobber_sed=False #overwrite SED only fits?
kiclist, perlist, pdeplist, sdeplist, morphlist = np.loadtxt('data/kebproperties_0216_full.dat',
                                          usecols=(0, 1, 3, 4, 8), unpack=True, delimiter=';')

#goodlist = (morphlist<0.6) & (pdeplist>0.1) & (sdeplist>0.01) & (perlist > 1.)
goodlist = (perlist>0)

excludelist = get_excludelist(fname='data/sed_flag_file_0328')

kic = int(sys.argv[1])
keblat = Keblat(preload=False)
keblat.loadiso2()
keblat.loadsed(sedfile='data/kepsedall_0216_full_test.dat')
keblat.loadvkeb(filename='data/kebproperties_0216_full.dat')
goodlist_ind = np.where(kiclist[goodlist].astype(int) == kic)[0]

goodv = keblat.kiclookup(kic, target=keblat.vkeb[:, 0])
keblat.loadlc(kic, keblat.vkeb[goodv, [1, 2, 5, 6, 7]])
magsobs, emagsobs, extinction, glat, z0 = keblat.getmags(kic)
ebv = extinction[0]
if np.isnan(z0):
    z0 = keblat.zsun
#print "Loading SED data, excluding ", excludelist[kic]
try:
    exclusive = excludelist[kic]
except:
    exclusive = []
keblat.isoprep(magsobs, emagsobs, extinction, glat, z0, exclude=exclusive)#exclude=[]) #excludelist[goodlist_ind])#'gmag','rmag','imag','zmag'])
period, tpe = keblat.vkeb[goodv, 1], keblat.vkeb[goodv, 2]
ecosw, esinw = keblat.vkeb[goodv, -2], keblat.vkeb[goodv, -1]
frat = (keblat.vkeb[goodv, 4]/keblat.vkeb[goodv, 3])

ebv_arr, ebv_sig, ebv_dist_bounds, ebv_bounds = None, None, None, None#get3dmap(kic)

lcbounds = np.array([(.2, 24.), (0.1, 1e4), (1e-4, 10.), (keblat.period-5., keblat.period+5.),
                   (keblat.tpe-10., keblat.tpe+10.), (0., 0.99), (0., 3.), (1e-6, 10.), (0.,1.),
                   (0.,1.), (0.,1.), (0.,1.), (-14., -3.)])

if keblat.swidth < 0.:
    print "No secondary eclipses detected. Exiting."
    sys.exit()
def make_sed_plots(kic, mlpars, prefix, suffix='', lc_constraints=None, savefig=True):
    isoerr = mlpars[-1]
    if isoerr < 0:
        isoerr = np.exp(isoerr)
    magsmod = keblat.isofit(mlpars)
    #magsmod, r, T, logg = isofit_single(mlpars[:5])
    if np.any(np.isinf(magsmod)):
        print "Input isopars give -inf magsmod: ", mlpars, magsmod
        return magsmod
    plt.figure()
    plt.subplot(311)
    plt.errorbar(keblat.maglams, keblat.magsobs, np.sqrt(keblat.emagsobs**2 + isoerr**2), fmt='k.')
    plt.plot(keblat.maglams, magsmod, 'r.')
    plt.ylabel('magnitude')

    plt.subplot(312)
    plt.errorbar(keblat.maglams, keblat.magsobs-magsmod, np.sqrt(keblat.emagsobs**2 + isoerr**2), fmt='k.')
    plt.xlabel('wavelength (angstrom)')
    plt.ylabel('data-model')
    plt.suptitle('KIC '+str(kic) +' SED (only)')

    plt.subplot(313)
    if lc_constraints is not None:
        plt.plot(lc_constraints, 'kx')
    lc_inputs = np.array([(keblat.r1 +keblat.r2)/(mlpars[0]+mlpars[1])**(1./3.), keblat.r2/keblat.r1, keblat.frat])
    plt.plot(lc_inputs, 'rx')
    plt.xlim((-1, 3))
    if savefig:
        plt.savefig(prefix + suffix+'.png')
    return magsmod


def make_lc_plots(kic, lcpars, prefix, suffix='', savefig=True, polyorder=2):
    keblat.updatephase(lcpars[4], lcpars[3])

    lcmod, lcpol = keblat.lcfit(lcpars[:13], keblat.jd[keblat.clip], keblat.quarter[keblat.clip],
                                keblat.flux[keblat.clip], keblat.fluxerr[keblat.clip],
                                keblat.crowd[keblat.clip], polyorder=polyorder)
    # phase = ((keblat.jd[keblat.clip]-lcpars[4]) % lcpars[3])/lcpars[3]
    # phase[phase<-np.clip(keblat.pwidth*3., 0., 0.2)]+=1.
    # phase[phase>np.clip(keblat.sep+keblat.swidth*3., keblat.sep, 1.0)]-=1.
    lcres = keblat.flux[keblat.clip] - lcmod*lcpol

    lcerr=0
    if len(lcpars)==14:
        lcerr=lcpars[-1]
        if lcerr<0:
            lcerr=np.exp(lcerr)

    # fig = plt.figure(figsize=(16,16))
    # ax = fig.add_subplot(111)
    # ax.errorbar(keblat.jd, keblat.flux, keblat.fluxerr, fmt='k.')
    # #ax.plot(keblat.jd[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'k.')
    #
    # ax.plot(keblat.jd[keblat.clip], lcmod*lcpol, 'r.')
    # ax.plot(keblat.jd[keblat.clip], lcpol, 'g.')
    # ax.set_ylabel('Kepler Flux')
    #
    # divider = make_axes_locatable(ax)
    # axb = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    # axb.errorbar(keblat.jd[keblat.clip], keblat.flux[keblat.clip]-lcmod*lcpol,
    #          keblat.fluxerr[keblat.clip], fmt='k.')
    # axb.set_ylabel('Data - Model')
    # axb.set_xlabel('BJD-24554833')
    #axb.set_yticklabels(axb.yaxis.get_majorticklabels()[1:])

    pe = (keblat.phase[keblat.clip] >= -1.2*keblat.pwidth) * (keblat.phase[keblat.clip] <= 1.2*keblat.pwidth)
    se = (keblat.phase[keblat.clip] >= -1.2*keblat.swidth+keblat.sep) * (keblat.phase[keblat.clip] <= 1.2*keblat.swidth+keblat.sep)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')

    ax1.errorbar(keblat.phase[keblat.clip][pe], keblat.flux[keblat.clip][pe]/lcpol[pe],
                 keblat.dflux[keblat.clip][pe], fmt='k.', ecolor='gray')
    ax1.plot(keblat.phase[keblat.clip][pe], lcmod[pe], 'r.')
    ax2.errorbar(keblat.phase[keblat.clip][se], keblat.flux[keblat.clip][se]/lcpol[se],
                 keblat.dflux[keblat.clip][se], fmt='k.', ecolor='gray')
    ax2.plot(keblat.phase[keblat.clip][se], lcmod[se], 'r.')
    ax3.errorbar(keblat.phase[keblat.clip][pe], lcres[pe],
                 keblat.dflux[keblat.clip][pe], fmt='k.', ecolor='gray')
    ax4.errorbar(keblat.phase[keblat.clip][se], lcres[se],
                 keblat.dflux[keblat.clip][se], fmt='k.', ecolor='gray')


    # ax1.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)

    ax1.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    # ax1.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))


    # ax2.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    # ax2.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))
    ax2.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))

    # ax3.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    # ax4.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)


    # ax4.set_ylim((np.min(lcres), np.max(lcres)))
    fig.suptitle('KIC '+str(kic)+' LC (only)')
    plt.subplots_adjust(hspace=0)
    if savefig:
        plt.savefig(prefix + suffix+'.png')
    return True


def make_sedlc_plots(kic, allpars, prefix, suffix='', savefig=True, polyorder=2):
    if allpars[-1] < 0:
        allpars[-1] = np.exp(allpars[-1])
    if allpars[4] < 10:
        allpars[4] = np.exp(allpars[4])
    if allpars[16] < 0:
        allpars[16] = np.exp(allpars[16])
    residuals = keblat.lnlike(allpars, lc_constraints=None, qua=np.unique(keblat.quarter),
                              polyorder=polyorder, residual=True)
    lcpars = [keblat.pars['m1'] + keblat.pars['m2'], keblat.r2+keblat.r1, keblat.r2/keblat.r1,
    		keblat.pars['period'], keblat.pars['tpe'], keblat.pars['esinw'], keblat.pars['ecosw'],
    		keblat.pars['b'], keblat.frat, keblat.pars['q1'], keblat.pars['q2'], keblat.pars['q3'],
    		keblat.pars['q4']]
    lcmod, lcpol = keblat.lcfit(lcpars, keblat.jd[keblat.clip], keblat.quarter[keblat.clip],
    				keblat.flux[keblat.clip], keblat.fluxerr[keblat.clip],
    				keblat.crowd[keblat.clip], polyorder=2)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(keblat.maglams, keblat.magsobs, keblat.emagsobs, fmt='k.')
    ax.plot(keblat.maglams, residuals[:len(keblat.maglams)] * \
             np.sqrt(keblat.emagsobs**2 + keblat.pars['isoerr']**2) + keblat.magsobs, 'r.')
    for ii in range(len(keblat.maglams)):
        ax.text(keblat.maglams[ii], keblat.magsobs[ii], keblat.ipname[ii].replace('mag', ''))
    ax.set_ylabel('Magnitude')
    #ax.set_yticklabels(ax.yaxis.get_majorticklabels()[:-1])

    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    ax2.errorbar(keblat.maglams, residuals[:len(keblat.maglams)] * \
             np.sqrt(keblat.emagsobs**2 + keblat.pars['isoerr']**2),
                 np.sqrt(keblat.emagsobs**2 + keblat.pars['isoerr']**2), fmt='k.')
    ax2.set_xlabel('Wavelength (Angstrom)')
    ax2.set_ylabel('Data - Model')
    ax2.set_ylim((-0.3, 0.3))
    #ax2.set_yticklabels(ax2.yaxis.get_majorticklabels()[:-1])
    plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_yticklabels()[-1], visible=False)

    plt.suptitle('KIC '+str(kic)+' SED (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_SED.png')

    phase = ((keblat.jd[keblat.clip]-lcpars[4]) % lcpars[3])/lcpars[3]
    phase[phase<-np.clip(keblat.pwidth*3., 0., 0.2)]+=1.
    phase[phase>np.clip(keblat.sep+keblat.swidth*3., keblat.sep, 1.0)]-=1.

    lcres = residuals[-np.sum(keblat.clip):] * \
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2)

    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(121)
    ax.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    ax.errorbar(phase, keblat.flux[keblat.clip]/lcpol,
                 keblat.fluxerr[keblat.clip], fmt='k.', ecolor='gray')
    ax.plot(phase, lcmod, 'r.')
    ax.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    ax.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))
    ax.set_ylabel('Kepler Flux')

    divider = make_axes_locatable(ax)
    axb = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    axb.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    axb.errorbar(phase, lcres,
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2), fmt='k.', ecolor='gray')

    axb.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    axb.set_ylim((np.min(lcres), np.max(lcres)))
    axb.set_ylabel('Data - Model')
    axb.set_xlabel('Phase (Primary Eclipse)')
    #axb.set_yticklabels(axb.yaxis.get_majorticklabels()[1:])

    ax2 = fig.add_subplot(122)
    ax2.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    ax2.errorbar(phase, keblat.flux[keblat.clip]/lcpol,
                 keblat.fluxerr[keblat.clip], fmt='k.', ecolor='gray')
    ax2.plot(phase, lcmod, 'r.')
    ax2.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))
    ax2.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))

    divider2 = make_axes_locatable(ax2)
    ax2b = divider2.append_axes("bottom", size=2.0, pad=0, sharex=ax2)
    ax2b.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    ax2b.errorbar(phase, lcres,
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2), fmt='k.', ecolor='gray')

    ax2b.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))
    ax2b.set_ylim((np.min(lcres), np.max(lcres)))
    ax2b.set_xlabel('Phase (Secondary Eclipse)')

    #ax2b.set_yticklabels(ax2b.yaxis.get_majorticklabels()[1:])


    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.suptitle('KIC '+str(kic)+' LC (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_LC.png')
    return True


def make_lcrv_plots(kic, allpars, prefix, suffix='', savefig=True, polyorder=2):
    residuals = lnlike_lcrv(allpars, qua=np.unique(keblat.quarter), polyorder=polyorder,
                            residual=True)
    lcpars = keblat.getpars(partype='lc')[:13]
    lcmod, lcpol = keblat.lcfit(lcpars, keblat.jd[keblat.clip], keblat.quarter[keblat.clip],
    				keblat.flux[keblat.clip], keblat.fluxerr[keblat.clip],
    				keblat.crowd[keblat.clip], polyorder=2)

    phase = ((keblat.jd[keblat.clip]-lcpars[4]) % lcpars[3])/lcpars[3]
    phase[phase<-np.clip(keblat.pwidth*3., 0., 0.2)]+=1.
    phase[phase>np.clip(keblat.sep+keblat.swidth*3., keblat.sep, 1.0)]-=1.

    lcres = keblat.flux[keblat.clip] - lcmod*lcpol

    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(121)
    ax.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    ax.errorbar(phase, keblat.flux[keblat.clip]/lcpol,
                 keblat.fluxerr[keblat.clip], fmt='k.', ecolor='gray')
    ax.plot(phase, lcmod, 'r.')
    ax.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    ax.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))
    ax.set_ylabel('Kepler Flux')

    divider = make_axes_locatable(ax)
    axb = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    axb.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    axb.errorbar(phase, lcres,
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2), fmt='k.', ecolor='gray')

    axb.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    axb.set_ylim((np.min(lcres), np.max(lcres)))
    axb.set_ylabel('Data - Model')
    axb.set_xlabel('Phase (Primary Eclipse)')
    #axb.set_yticklabels(axb.yaxis.get_majorticklabels()[1:])

    ax2 = fig.add_subplot(122)
    ax2.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    ax2.errorbar(phase, keblat.flux[keblat.clip]/lcpol,
                 keblat.fluxerr[keblat.clip], fmt='k.', ecolor='gray')
    ax2.plot(phase, lcmod, 'r.')
    ax2.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))
    ax2.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))

    divider2 = make_axes_locatable(ax2)
    ax2b = divider2.append_axes("bottom", size=2.0, pad=0, sharex=ax2)
    ax2b.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    ax2b.errorbar(phase, lcres,
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2), fmt='k.', ecolor='gray')

    ax2b.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))
    ax2b.set_ylim((np.min(lcres), np.max(lcres)))
    ax2b.set_xlabel('Phase (Secondary Eclipse)')

    #ax2b.set_yticklabels(ax2b.yaxis.get_majorticklabels()[1:])


    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.suptitle('KIC '+str(kic)+' LC (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_LC.png')

    # rvpars = keblat.getpars(partype='rv')
    rvpars = np.array([keblat.pars['msum'], keblat.pars['mrat'], keblat.pars['period'], keblat.pars['tpe'],
                       keblat.pars['esinw'], keblat.pars['ecosw'], keblat.pars['inc'], keblat.pars['k0'], keblat.pars['rverr']])

    rv_fit = keblat.rvfit(rvpars, keblat.rv_t)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rvphase = (keblat.rv_t - keblat.pars['tpe'])%keblat.pars['period']/keblat.pars['period']
    ax.errorbar(rvphase[~keblat.bad1], keblat.rv1_obs[~keblat.bad1], keblat.rv1_err_obs[~keblat.bad1], fmt='b*')
    ax.errorbar(rvphase[~keblat.bad2], keblat.rv2_obs[~keblat.bad2], keblat.rv2_err_obs[~keblat.bad2], fmt='r*')
    rvt = np.linspace(0, 1, 100)*keblat.pars['period']+keblat.pars['tpe']
    rvmod = keblat.rvfit(rvpars, rvt)
    ax.plot(np.linspace(0, 1, 100), rvmod[0], 'b-')
    ax.plot(np.linspace(0, 1, 100), rvmod[1], 'r-')
    ax.set_ylabel('RV [m/s]')
    #ax.set_yticklabels(ax.yaxis.get_majorticklabels()[:-1])

    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    ax2.errorbar(rvphase[~keblat.bad1], (keblat.rv1_obs-rv_fit[0])[~keblat.bad1], np.sqrt(keblat.rv1_err_obs**2+rvpars[-1]**2)[~keblat.bad1], fmt='b.')
    ax2.errorbar(rvphase[~keblat.bad2], (keblat.rv2_obs-rv_fit[1])[~keblat.bad2], np.sqrt(keblat.rv2_err_obs**2+rvpars[-1]**2)[~keblat.bad2], fmt='r.')

    ax2.set_xlabel('Phase')
    ax2.set_ylabel('Data - Model')
    #ax2.set_yticklabels(ax2.yaxis.get_majorticklabels()[:-1])
    plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_yticklabels()[-1], visible=False)

    plt.suptitle('KIC '+str(kic)+' RV (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_RV.png')

    return True

def make_sedlcrv_plots(kic, allpars, prefix, suffix='', savefig=True, polyorder=2):
    residuals = keblat.lnlike_all(allpars, lc_constraints=None, qua=np.unique(keblat.quarter),
                              polyorder=polyorder, residual=True)
    lcpars = keblat.getpars(partype='lc')[:13]
    lcmod, lcpol = keblat.lcfit(lcpars, keblat.jd[keblat.clip], keblat.quarter[keblat.clip],
    				keblat.flux[keblat.clip], keblat.fluxerr[keblat.clip],
    				keblat.crowd[keblat.clip], polyorder=2)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(keblat.maglams, keblat.magsobs, keblat.emagsobs, fmt='k.')
    ax.plot(keblat.maglams, residuals[:len(keblat.maglams)] * \
             np.sqrt(keblat.emagsobs**2 + keblat.pars['isoerr']**2) + keblat.magsobs, 'r.')
    for ii in range(len(keblat.maglams)):
        ax.text(keblat.maglams[ii], keblat.magsobs[ii], keblat.ipname[ii].replace('mag', ''))
    ax.set_ylabel('Magnitude')
    #ax.set_yticklabels(ax.yaxis.get_majorticklabels()[:-1])

    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    ax2.errorbar(keblat.maglams, residuals[:len(keblat.maglams)] * \
             np.sqrt(keblat.emagsobs**2 + keblat.pars['isoerr']**2),
                 np.sqrt(keblat.emagsobs**2 + keblat.pars['isoerr']**2), fmt='k.')
    ax2.set_xlabel('Wavelength (Angstrom)')
    ax2.set_ylabel('Data - Model')
    ax2.set_ylim((-0.3, 0.3))
    #ax2.set_yticklabels(ax2.yaxis.get_majorticklabels()[:-1])
    plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_yticklabels()[-1], visible=False)

    plt.suptitle('KIC '+str(kic)+' SED (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_SED.png')

    phase = ((keblat.jd[keblat.clip]-lcpars[4]) % lcpars[3])/lcpars[3]
    phase[phase<-np.clip(keblat.pwidth*3., 0., 0.2)]+=1.
    phase[phase>np.clip(keblat.sep+keblat.swidth*3., keblat.sep, 1.0)]-=1.

    lcres = keblat.flux[keblat.clip] - lcmod*lcpol

    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(121)
    ax.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    ax.errorbar(phase, keblat.flux[keblat.clip]/lcpol,
                 keblat.fluxerr[keblat.clip], fmt='k.', ecolor='gray')
    ax.plot(phase, lcmod, 'r.')
    ax.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    ax.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))
    ax.set_ylabel('Kepler Flux')

    divider = make_axes_locatable(ax)
    axb = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    axb.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    axb.errorbar(phase, lcres,
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2), fmt='k.', ecolor='gray')

    axb.set_xlim((-1.2*keblat.pwidth, 1.2*keblat.pwidth))
    axb.set_ylim((np.min(lcres), np.max(lcres)))
    axb.set_ylabel('Data - Model')
    axb.set_xlabel('Phase (Primary Eclipse)')
    #axb.set_yticklabels(axb.yaxis.get_majorticklabels()[1:])

    ax2 = fig.add_subplot(122)
    ax2.plot(keblat.phase[keblat.clip], keblat.flux[keblat.clip]/lcpol, 'g.', alpha=0.4)
    ax2.errorbar(phase, keblat.flux[keblat.clip]/lcpol,
                 keblat.fluxerr[keblat.clip], fmt='k.', ecolor='gray')
    ax2.plot(phase, lcmod, 'r.')
    ax2.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))
    ax2.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))

    divider2 = make_axes_locatable(ax2)
    ax2b = divider2.append_axes("bottom", size=2.0, pad=0, sharex=ax2)
    ax2b.plot(keblat.phase[keblat.clip], lcres, 'g.', alpha=0.4)
    ax2b.errorbar(phase, lcres,
                 np.sqrt(keblat.fluxerr[keblat.clip]**2 + keblat.pars['lcerr']**2), fmt='k.', ecolor='gray')

    ax2b.set_xlim((-1.2*keblat.swidth+keblat.sep, 1.2*keblat.swidth+keblat.sep))
    ax2b.set_ylim((np.min(lcres), np.max(lcres)))
    ax2b.set_xlabel('Phase (Secondary Eclipse)')

    #ax2b.set_yticklabels(ax2b.yaxis.get_majorticklabels()[1:])


    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.suptitle('KIC '+str(kic)+' LC (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_LC.png')

    rvpars = keblat.getpars(partype='rv')
    rv_fit = keblat.rvfit(rvpars, keblat.rv_t)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    rvphase = (keblat.rv_t - keblat.pars['tpe'])%keblat.pars['period']/keblat.pars['period']
    ax.errorbar(rvphase, keblat.rv1_obs, keblat.rv1_err_obs, fmt='b*')
    ax.errorbar(rvphase, keblat.rv2_obs, keblat.rv2_err_obs, fmt='r*')
    rvt = np.linspace(0, 1, 100)*keblat.pars['period']+keblat.pars['tpe']
    rvmod = keblat.rvfit(rvpars, rvt)
    ax.plot(np.linspace(0, 1, 100), rvmod[0], 'b-')
    ax.plot(np.linspace(0, 1, 100), rvmod[1], 'r-')
    ax.set_ylabel('RV [cm/s]')
    #ax.set_yticklabels(ax.yaxis.get_majorticklabels()[:-1])

    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
    ax2.errorbar(rvphase, (keblat.rv1_obs-rv_fit[0]), np.sqrt(keblat.rv1_err_obs**2+rvpars[-1]**2), fmt='b.')
    ax2.errorbar(rvphase, (keblat.rv2_obs-rv_fit[1]), np.sqrt(keblat.rv2_err_obs**2+rvpars[-1]**2), fmt='r.')

    ax2.set_xlabel('Phase')
    ax2.set_ylabel('Data - Model')
    #ax2.set_yticklabels(ax2.yaxis.get_majorticklabels()[:-1])
    plt.setp(ax.get_xticklabels(), visible=False)
    #plt.setp(ax2.get_yticklabels()[-1], visible=False)

    plt.suptitle('KIC '+str(kic)+' RV (simultaneous)')
    if savefig:
        plt.savefig(prefix+suffix+'_RV.png')

    return True

def get_pars2vals(fisopars, partype='lc', crow=False):
    if partype == 'lc':
        parnames = ['msum', 'rsum', 'rrat', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'frat',
                    'q1', 'q2', 'q3', 'q4']
        if crow:
            parnames += ['cr' + str(ii) for ii in np.unique(keblat.quarter)]
    elif partype == 'sed':
        parnames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'isoerr']
    elif partype == 'rv':
        parnames = ['msum', 'mrat', 'period', 'tpe', 'esinw', 'ecosw', 'inc', 'k0', 'rverr']
    elif partype == 'lcsed':
        parnames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'period', 'tpe',
                         'esinw', 'ecosw', 'b', 'q1', 'q2', 'q3', 'q4', 'lcerr', 'isoerr']
    elif partype == 'lcrv':
        parnames = ['msum', 'mrat', 'rsum', 'rrat', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'frat', 'q1',
                'q2', 'q3', 'q4', 'lcerr', 'k0', 'rverr']
    else:
        print "Partype options are 'lc', 'sed', 'rv', 'lcsed', 'lcrv'. Try again."
        return

    parvals = np.zeros(len(parnames))
    if isinstance(fisopars, lmfit.parameter.Parameters):
        for j in range(len(parnames)):
            parvals[j] = fisopars[parnames[j]].value
    elif isinstance(fisopars, dict):
        for j in range(len(parnames)):
            parvals[j] = fisopars[parnames[j]]
    else:
        for j in range(len(parnames)):
            parvals[j] = fisopars[j]
    return parvals

def rez(fit_params, polyorder=0):
    guess = get_pars2vals(fit_params, partype='lc')
    if len(guess)>13:
        crowd_fits = keblat.broadcast_crowd(keblat.quarter, guess[-len(np.unique(keblat.quarter)):])
    else:
        crowd_fits = keblat.crowd
    lcmod, lcpol = keblat.lcfit(guess[:13], keblat.jd[keblat.clip],
                                keblat.quarter[keblat.clip], keblat.flux[keblat.clip],
                                keblat.fluxerr[keblat.clip], crowd_fits[keblat.clip],
                                polyorder=polyorder)
    if np.any(np.isinf(lcmod)):
        return np.ones(np.sum(keblat.clip))*1e10#(1.-keblat.flux[keblat.clip])/keblat.dflux[keblat.clip]
    return (keblat.flux[keblat.clip] - lcmod*lcpol)/keblat.fluxerr[keblat.clip]

def lnlike_rv(fisopars, residual=True):
    pars = get_pars2vals(fisopars, partype='rv')
    rv1, rv2 = keblat.rvfit(pars, keblat.rv_t)
    res = np.concatenate(((rv1[~keblat.bad1] - keblat.rv1_obs[~keblat.bad1]) /
                          np.sqrt(keblat.rv1_err_obs[~keblat.bad1]**2 + pars[-1]**2),
                          (rv2[~keblat.bad2] - keblat.rv2_obs[~keblat.bad2]) /
                          np.sqrt(keblat.rv2_err_obs[~keblat.bad2]**2 + pars[-1]**2)))
    if residual:
        if np.any(np.isinf(res)) or np.sum(np.isnan(res)) > 0.05 * len(res):
            return np.ones(np.sum(~keblat.bad1)+np.sum(~keblat.bad2)) * 1e20
        return res
    return np.sum(res**2)

def lnlike_lcrv(fisopars, qua=[1], polyorder=2, residual=True):
    pars = get_pars2vals(fisopars, partype='lcrv')
    res = keblat.lnlike_lcrv(pars, qua=qua, polyorder=polyorder, residual=residual)
    if residual:
        if np.any(np.isinf(res)) or np.sum(np.isnan(res)) > 0.05 * len(res):
            return np.ones(np.sum(~keblat.bad1)+np.sum(~keblat.bad2) + np.sum(keblat.clip)) * 1e20
        return res
    return np.sum(res**2)

def lnlike_lmfit(fisopars, lc_constraints=None, ebv_arr=None, qua=[1], polyorder=2, residual=False):
    allpars = get_pars2vals(fisopars, partype='lcsed')

    if ebv_arr is not None:
        allpars[5] = np.interp(allpars[4], ebv_dist, ebv_arr)
    #print "sum of clips = ", keblat.clip.sum()
    res = keblat.lnlike(allpars, lc_constraints=lc_constraints, qua=qua, polyorder=polyorder,
                         residual=residual)
    if np.any(np.isinf(res)):
        print "Inf res"
        extra = 0 if lc_constraints is None else len(lc_constraints)
        return np.ones(len(keblat.magsobs)+len(keblat.flux[keblat.clip]) + extra)*1e20
    bads = np.isnan(res)
    if np.sum(bads) > 0.05*len(res):
        print "Seems to be a lot of Nans..."
        extra = 0 if lc_constraints is None else len(lc_constraints)
        return np.ones(len(keblat.magsobs)+len(keblat.flux[keblat.clip]) + extra)*1e20
    return res

def ew_search_lmfit(ew_trials, pars0, argpars, fit_ecosw=True, polyorder=0):
    fit_ew = Parameters()
    fit_ew.add('esinw', value=pars0[5], min=-0.9, max=0.9, vary=True)
    fit_ew.add('ecosw', value=pars0[6], min=-0.9, max=0.9, vary=fit_ecosw)
    chisq = 1e18
    for ii in range(len(ew_trials)):
        fit_ew['esinw'].value=ew_trials[ii][0]
        fit_ew['ecosw'].value=ew_trials[ii][1]
        result = minimize(ew_search, fit_ew, kws={'pars0':pars0, 'argpars':argpars, 'polyorder':polyorder})
        if result.redchi < chisq or ii==0:
            chisq = result.redchi * 1.0
            ew_best = result.params['esinw'].value, result.params['ecosw'].value
            print "Better redchi: ", chisq, result.redchi, ew_best
    return ew_best

def ew_search(ew, pars0=None, argpars=None, polyorder=0, retmod=False):
    pars = np.array(pars0).copy()
    #esinw, ecosw = ew
    try:
        pars[5], pars[6] = ew['esinw'].value, ew['ecosw'].value
    except:
        pars[5], pars[6] = ew[0], ew[1]
    mod, poly = keblat.lcfit(pars, keblat.jd, keblat.quarter, keblat.flux, keblat.fluxerr, keblat.crowd, polyorder=polyorder)
    # _phasesort = np.argsort(keblat.phase)
    # _phase = keblat.phase[_phasesort]
    # _flux = keblat.flux[_phasesort]
    # _dflux = keblat.dflux[_phasesort]
    # Nbins = int(1./min(keblat.pwidth, keblat.swidth))*4
    # bins = np.linspace(_phase[0], _phase[-1], Nbins)
    # digitized = np.digitize(_phase, bins)
    # x_means = [_phase[digitized == i].mean() for i in range(1, len(bins))]
    # y_means = [_flux[digitized == i].mean() for i in range(1, len(bins))]
    # y_means /= y_means
    #
    if np.any(np.isinf(mod)):
        return np.ones_like(keblat.jd+1)*1e10
    if retmod:
        return mod, poly
    return np.append((keblat.flux - mod*poly)/keblat.fluxerr, tse_residuals((pars[5],pars[6]), *argpars))


def get_age_trials(mass):
    if np.log10(mass) <= 0.02:
        return np.log10(np.exp(np.linspace(np.log(3e6), np.log(7e9), 5)))
    else:
        agelim = -2.7 * np.log10(mass) + 9.9
        return np.log10(np.exp(np.linspace(np.log(3e6), np.log(10**agelim), 5)))

def sed_chisq(sedpars0, lc_constraints):
    # ll,_ = keblat.ilnlike(sedpars0, lc_constraints=lc_constraints)
    # return ll/(-0.5)
    residuals = keblat.ilnlike(sedpars0, lc_constraints=lc_constraints, residual=True)
    if np.all(residuals)==1e12:
        return residuals
    residuals[:-3] *= np.sqrt(keblat.emagsobs**2+np.exp(sedpars0[-1])**2)
    lc_inputs = np.array([(keblat.r1+keblat.r2)/(sedpars0[0]+sedpars0[1])**(1./3.),
                          keblat.r2/keblat.r1, keblat.frat])

    residuals[-3:] = np.log(lc_constraints) - np.log(lc_inputs)
    #residuals = np.append(residuals, lc_inputs)
    #residuals = np.append(residuals, lc_inputs)
    return residuals


def opt_sed(sedpars0, lc_constraints, ebv_dist, ebv_arr, fit_ebv=True, ret_lcpars=True, vary_z0=True):
    m1, m2, z0, age, dist, ebv, h0, isoerr = sedpars0
    fit_params2 = Parameters()
    fit_params2.add('m1', value=m1, min=0.1, max=12.)
    fit_params2.add('m2', value=m2, min=0.1, max=12.)
    fit_params2.add('z0', value=z0, min=0.001, max=0.06, vary=vary_z0)
    fit_params2.add('age', value=age, min=6.0, max=10.1, vary=True)
    fit_params2.add('dist', value=dist, min=10., max=15000.)
    kws = {'lc_constraints': lc_constraints, 'residual': True, 'ebv_arr': ebv_arr, 'ebv_dist': ebv_dist}
    if fit_ebv:
        fit_params2.add('ebv', value=ebv, min=0.0, max=1.0)#ebv[0], vary=False)
        kws['ebv_arr'] = None
        kws['ebv_dist'] = None
    fit_params2.add('h0', value=h0, vary=False)
    fit_params2.add('isoerr', value=isoerr, min=-8, max=0.,
                    vary=False)
#    fit_params2.add('mbound', expr='m1-m2>0')

#    if isoerr>0.1:
#    fit_params2['isoerr'].vary=True

    fit_kws={'maxfev':2000*(len(fit_params2)+1)}
    if len(keblat.magsobs)<6:
        fit_params2['ebv'].vary=False

    print "=========================================================================="
    print "======================== Starting SED ONLY fit... ========================"
    print "=========================================================================="

    result2 = minimize(keblat.ilnlike, fit_params2, kws=kws, iter_cb=MinimizeStopper(30), **fit_kws)
    isores = keblat.ilnlike(result2.params, lc_constraints=lc_constraints, residual=True)
    redchi2 = np.sum(isores**2) / len(isores)
    print redchi2, result2.redchi
    report_fit(result2)
    fit_params2 = result2.params.copy()

    # fit_params2['z0'].vary=True
    # result2 = minimize(keblat.ilnlike, fit_params2, kws=kws, **fit_kws)
    # isores = keblat.ilnlike(result2.params, lc_constraints=lc_constraints, residual=True)
    # current_redchi2 = np.sum(isores**2) / len(isores)
    # if current_redchi2 < redchi2:
    #     print current_redchi2, result2.redchi
    #     report_fit(result2)
    #     print result2.message
    #     fit_params2 = result2.params.copy()
    #     redchi2 = current_redchi2*1.

    niter=0
    while (redchi2>1.) and (niter<3):
        result2 = minimize(keblat.ilnlike, fit_params2, kws=kws, iter_cb=MinimizeStopper(30), **fit_kws)
        #redchi2_0 = keblat.ilnlike(result2.params, lc_constraints=lc_constraints)
        isores = keblat.ilnlike(result2.params, lc_constraints=lc_constraints, residual=True)
        current_redchi2 = np.sum(isores**2) / len(isores)
        print "Iteration: ", niter, current_redchi2, result2.redchi, result2.nfev
        if current_redchi2 < redchi2:
            print "Saving the following results:"
            report_fit(result2)
            redchi2 = current_redchi2*1.0
            fit_params2 = result2.params.copy()
        niter+=1

    #print result2.params, fit_params2
    if ret_lcpars:
        return fit_params2, keblat.r1, keblat.r2, keblat.frat
    return fit_params2


def opt_lc(lcpars0, jd, phase, flux, dflux, crowd, clip, set_upperb=2., fit_crowd=False, fit_se=False,
           vary_msum=True):
    msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4 = lcpars0

    fit_params = Parameters()
    fit_params.add('esinw', value=esinw, min=-.999, max=0.999, vary=False)
    fit_params.add('ecosw', value=ecosw, min=-.999, max=0.999, vary=False)#ecosw-0.05, max=ecosw+0.05, vary=False)
    fit_params.add('rsum', value=rsum, min=0.1, max=10000., vary=False)
    fit_params.add('rrat', value=rrat, min=1e-4, max=1e3, vary=False)
    fit_params.add('b', value=b, min=0., max=set_upperb, vary=False)
    fit_params.add('frat', value=frat, min=1e-6, max=1e2, vary=False)
    fit_params.add('msum', value=msum, min=0.2, max=24., vary=False)

    fit_params.add('period', value=period, min=period-0.005, max=period+0.005, vary=False)
    fit_params.add('tpe', value=tpe, min=tpe-10., max=tpe+10., vary=False)
    fit_params.add('q1', value=q1, min=0., max=1., vary=False)
    fit_params.add('q2', value=q2, min=0., max=1., vary=False)
    fit_params.add('q3', value=q3, min=0., max=1., vary=False)
    fit_params.add('q4', value=q4, min=0., max=1., vary=False)

    fit_params['rrat'].vary=True
    fit_params['rsum'].vary=True
    fit_params['b'].vary=True
    fit_params['frat'].vary=True
    fit_params['esinw'].vary=True
    fit_params['ecosw'].vary=True
    fit_params['tpe'].vary=True

    fit_kws={'maxfev':100*(len(fit_params)+1)}

    if fit_crowd:
        print "Fitting crowding parameters..."
        for ii in fit_params.keys():
            fit_params[ii].vary=True
        for ii in range(len(np.unique(keblat.quarter))):
            fit_params.add('cr'+str(np.unique(keblat.quarter)[ii]), value=keblat._crowdsap[ii], min=0.1, max=1.0)
        result0 = minimize(rez, fit_params, kws={'polyorder':2}, iter_cb=MinimizeStopper(10), **fit_kws)
        report_fit(result0)

        for ii in result0.params.keys():
            result0.params[ii].vary=True
        niter=0
        redchi2 = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=2))**2) / np.sum(keblat.clip)
        guess=get_pars2vals(result0.params, partype='lc')
        while (redchi2>1.) and (niter<5):
            result0 = minimize(rez, result0.params, kws={'polyorder':2}, iter_cb=MinimizeStopper(10), **fit_kws)
            current_chi = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=2))**2) / np.sum(keblat.clip)
            if current_chi < redchi2:
                redchi2=current_chi*1.0
                report_fit(result0)
                guess=get_pars2vals(result0.params, partype='lc')
            niter+=1

        return guess

    print "=========================================================================="
    print "==================== Starting LIGHTCURVE ONLY fit... ====================="
    print "=========================================================================="
    if fit_se:
        if ecosw == 0.:
            ecosw=1e-5
        if esinw == 0.:
            esinw=1e-5
        fit_params['esinw'].min=esinw-0.1*abs(esinw)
        fit_params['ecosw'].min=ecosw-0.05*abs(ecosw)
        fit_params['esinw'].max=esinw+0.1*abs(esinw)
        fit_params['ecosw'].max=ecosw+0.05*abs(ecosw)
        
        result0 = minimize(rez, fit_params, kws={'polyorder': 0}, iter_cb=MinimizeStopper(10), **fit_kws)
        report_fit(result0)
        fit_params = result0.params
    
    result0 = minimize(rez, fit_params, kws={'polyorder': 1}, iter_cb=MinimizeStopper(10), **fit_kws)
    report_fit(result0)


#    redchi2 = np.sum((result0.residual)**2) / (len(result0.residual)-result0.nfev)
    #guess = get_lcvals(result0.params)
    fit_params = result0.params
    redchi2 = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=1))**2) / np.sum(keblat.clip)

    fit_params['msum'].vary=vary_msum
    fit_params['tpe'].vary=True
    fit_params['period'].vary=True
    fit_params['b'].vary=True
    fit_params['frat'].vary=True
    fit_params['esinw'].vary=True
    fit_params['ecosw'].vary=True
    fit_params['rsum'].vary=True
    fit_params['rrat'].vary=True
    fit_params['q1'].vary=True
    fit_params['q2'].vary=True
    fit_params['q3'].vary=True
    fit_params['q4'].vary=True

    result0 = minimize(rez, fit_params, kws={'polyorder': 1}, iter_cb=MinimizeStopper(10), **fit_kws)
#    current_redchi = np.sum((result0.residual)**2) / (len(result0.residual)-result0.nfev)
    current_redchi = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=1))**2) / np.sum(keblat.clip)

    if current_redchi < redchi2:
        redchi2 = current_redchi * 1.
        #guess = get_lcvals(result0.params)
        fit_params = result0.params
        print "polyorder = 1: ", current_redchi, result0.redchi
        report_fit(result0)

    result0 = minimize(rez, fit_params, kws={'polyorder': 2}, iter_cb=MinimizeStopper(10), **fit_kws)
#    current_redchi = np.sum((result0.residual)**2) / (len(result0.residual)-result0.nfev)
    current_redchi = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=2))**2) / np.sum(keblat.clip)
    if current_redchi < redchi2:
        redchi2 = current_redchi * 1.
        #guess = get_lcvals(result0.params)
        fit_params = result0.params
        print "polyorder = 2: ", current_redchi, result0.redchi
        report_fit(result0)


#     fit_params['rsum'].vary=True
#     fit_params['rrat'].vary=True
#     niter=0
#
#     while (redchi2>1.) and (niter<5):
#         result0 = minimize(rez, fit_params, kws={'polyorder': 2}, iter_cb=MinimizeStopper(10), **fit_kws)
#         current_redchi = np.sum((rez(get_lcvals(result0.params), polyorder=2))**2) / np.sum(keblat.clip)
#         print "Iteration: ", niter, redchi2, current_redchi, result0.redchi, result0.nfev#, get_lcvals(result0.params)
# #        current_redchi = np.sum((result0.residual)**2) / (len(result0.residual)-result0.nfev)
#         if current_redchi < redchi2:
#             print "Saving the following results:"
#             report_fit(result0)
#             redchi2 = current_redchi * 1.
#             #guess = get_lcvals(result0.params)
#             fit_params = result0.params
#         niter+=1
    guess = get_pars2vals(fit_params, partype='lc')
    return guess


def opt_sedlc(fit_params2, guess, ebv_dist, ebv_arr, jd, phase, flux, dflux, crowd, clip, mciso=None, fit_ebv=True,
              set_upperb = None, init_porder=1, init_varyb=False, init_varyew=False, init_varyza=False, lc_constraints=False):
    isonames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'isoerr']

    fit_params = Parameters()
    """
    for ii in range(len(keblat.pars)):
        fit_params.add(keblat.pars.keys()[ii], value=initvals[ii], min=keblat.parbounds.values()[ii][0],
                        max=keblat.parbounds.values()[ii][1])
    fit_params.pop('lcerr')
    fit_params['isoerr'].vary=False
    fit_params['h0'].vary=False
    kws = {'qua': np.unique(keblat.quarter), 'polyorder': 2, 'residual': True, 'ebv_arr': None}

    if not fit_ebv:
        fit_params.pop['ebv']
        kws = {'qua': np.unique(keblat.quarter), 'polyorder': 2, 'residual': True, 'ebv_arr': ebv_arr}
    """
    fit_params = fit_params2
    if mciso is not None:
        for ii in range(len(isonames)):
            fit_params[isonames[ii]].value = mciso[ii]
    fit_params['ebv'].max = 1.0
    fit_params['age'].vary=True
    fit_params['m1'].vary=True
    fit_params['m2'].vary=True
    fit_params['z0'].vary=True
    fit_params['ebv'].vary=True
    fit_params['dist'].vary=True
    fit_params['isoerr'].vary=False

    kws = {'lc_constraints': lc_constraints,
           'qua': np.unique(keblat.quarter), 'polyorder': init_porder, 'residual': True, 'ebv_arr': None}

    fit_params.add('esinw', value=guess[5], min=-.999, max=0.999, vary=False)
    fit_params.add('ecosw', value=guess[6], min=guess[6]-0.05, max=guess[6]+0.05, vary=False)
    fit_params.add('tpe', value=guess[4], min=tpe-10., max=tpe+10., vary=False)
    fit_params.add('period', value=guess[3], min=period-0.005, max=period+0.005, vary=False)
    fit_params.add('q1', value=guess[-4], min=0., max=1., vary=False)
    fit_params.add('q2', value=guess[-3], min=0., max=1., vary=False)
    fit_params.add('q3', value=guess[-2], min=0., max=1., vary=False)
    fit_params.add('q4', value=guess[-1], min=0., max=1., vary=False)
    fit_params.add('lcerr', value=1e-6, min=0., max=1e-3, vary=False)
    if set_upperb is None:
        fit_params.add('b', value=guess[7], min=0., vary=False)
    else:
        fit_params.add('b', value=guess[7], min=0., max=set_upperb, vary=False)

    fit_params['esinw'].vary=False
    fit_params['ecosw'].vary=False
    if init_varyew:
        fit_params['esinw'].vary=True
    if init_varyb:
        fit_params['b'].vary=True
    if init_varyza:
        fit_params['z0'].vary=False
        fit_params['age'].vary=False

    fit_kws={'maxfev':2000*(len(fit_params)+1)}

    print "=========================================================================="
    print "================= Starting SED + LC simultaneous fit... =================="
    print "=========================================================================="

    print fit_params
    result3 = minimize(lnlike_lmfit, fit_params, kws=kws, iter_cb=MinimizeStopper(60), **fit_kws)


    fit_params = result3.params.copy()
    _allres = lnlike_lmfit(result3.params, lc_constraints=lc_constraints, qua=np.unique(keblat.quarter), polyorder=2, residual=True)
    redchi2 = np.sum(_allres**2) / len(_allres)
    print redchi2, result3.redchi
    report_fit(result3)
    print result3.message

    fit_params['age'].vary=True
    fit_params['m1'].vary=True
    fit_params['m2'].vary=True
    fit_params['z0'].vary=True
    fit_params['dist'].vary=True
    fit_params['esinw'].vary=True
    fit_params['ecosw'].vary=True
    fit_params['b'].vary=True
    fit_params['q1'].vary=True
    fit_params['q2'].vary=True
    fit_params['q3'].vary=True
    fit_params['q4'].vary=True

    result3 = minimize(lnlike_lmfit, fit_params, kws=kws, iter_cb=MinimizeStopper(60), **fit_kws)
    _allres = lnlike_lmfit(result3.params, lc_constraints=lc_constraints, qua=np.unique(keblat.quarter), polyorder=2, residual=True)
    current_redchi2 = np.sum(_allres**2) / len(_allres)
    if current_redchi2 < redchi2:
        print "The following results are saved:", current_redchi2, result3.redchi
        report_fit(result3)
        print result3.message
        fit_params = result3.params.copy()
        redchi2 = current_redchi2

    kws['polyorder'] = 2

    fit_params['tpe'].vary=True
    fit_params['period'].vary=True

    niter=0
    while (niter<10): #(redchi2>1.) and (niter<10):
        result3 = minimize(lnlike_lmfit, fit_params, kws=kws, iter_cb=MinimizeStopper(60), **fit_kws)
        _allres = lnlike_lmfit(result3.params, lc_constraints=lc_constraints, qua=np.unique(keblat.quarter), polyorder=2, residual=True)
        current_redchi2 = np.sum(_allres**2) / len(_allres)
        print "Iteration: ", niter, current_redchi2, result3.redchi

        if current_redchi2 < redchi2:
            print "The following results are saved:"
            report_fit(result3)
            fit_params = result3.params.copy()
            redchi2 = current_redchi2
        niter+=1
    #print "logL of best allpars = ", keblat.lnlike(allpars, lc_constraints=None, qua=np.unique(keblat.quarter), polyorder=2)
    allpars = get_pars2vals(fit_params, partype='lcsed')
    #print fit_params
    return allpars

def opt_rv(**kwargs):
    fit_pars = Parameters()
    for name, val in kwargs.items():
        print name, val
        fit_pars.add(name, value=val, min=keblat.parbounds[name][0], max=keblat.parbounds[name][1])
    fit_pars['rverr'].vary=False
    fit_pars['period'].vary=False
    fit_pars['tpe'].vary=False
    # fit_pars['m1'].min=fit_pars['m1'].value*0.5
    # fit_pars['m2'].min=fit_pars['m2'].value*0.5
    # fit_pars['m1'].max=fit_pars['m1'].value*1.5
    # fit_pars['m2'].max=fit_pars['m2'].value*1.5
    # fit_pars['k0'].min=fit_pars['k0'].value - abs(fit_pars['k0'].value)
    # fit_pars['k0'].min=fit_pars['k0'].value + abs(fit_pars['k0'].value)

    fit_kws = {'maxfev': 2000 * (len(fit_pars) + 1)}

    print "=========================================================================="
    print "========================= Starting RV ONLY fit... ========================"
    print "=========================================================================="
    print fit_pars
    result3 = minimize(lnlike_rv, fit_pars, iter_cb=MinimizeStopper(60), **fit_kws)
    fit_params = result3.params.copy()
    _allres = lnlike_rv(result3.params)
    redchi2 = np.sum(_allres ** 2) / len(_allres)
    print redchi2, result3.redchi
    report_fit(result3)
    print result3.message
    niter = 0
    while (niter < 10):  # (redchi2>1.) and (niter<10):
        result3 = minimize(lnlike_rv, fit_params,iter_cb=MinimizeStopper(60), **fit_kws)
        _allres = lnlike_rv(result3.params)
        current_redchi2 = np.sum(_allres ** 2) / len(_allres)
        print "Iteration: ", niter, current_redchi2, result3.redchi
        if current_redchi2 < redchi2:
            print "The following results are saved:"
            report_fit(result3)
            fit_params = result3.params.copy()
            redchi2 = current_redchi2
        niter += 1
    # print "logL of best allpars = ", keblat.lnlike(allpars, lc_constraints=None, qua=np.unique(keblat.quarter), polyorder=2)
    rvpars = get_pars2vals(fit_params, partype='rv')
    # print fit_params
    return rvpars

def opt_lcrv(**kwargs):
    fit_pars = Parameters()
    for name, val in kwargs.items():
        if name in keblat.parbounds.keys():
            fit_pars.add(name, value=val, min=keblat.parbounds[name][0], max=keblat.parbounds[name][1])
        else:
            fit_pars.add(name, value=val)
    kws = {'qua': np.unique(keblat.quarter), 'polyorder': 2, 'residual': True}
    fit_kws = {'maxfev': 2000 * (len(fit_pars) + 1)}


    print "=========================================================================="
    print "================= Starting LC + RV simultaneous fit... ==================="
    print "=========================================================================="
    result3 = minimize(lnlike_lcrv, fit_pars, kws=kws, iter_cb=MinimizeStopper(60), **fit_kws)
    fit_params = result3.params.copy()
    _allres = lnlike_lcrv(result3.params, qua=np.unique(keblat.quarter), polyorder=2, residual=True)
    redchi2 = np.sum(_allres ** 2) / len(_allres)
    print redchi2, result3.redchi
    report_fit(result3)
    print result3.message
    niter = 0
    while (niter < 10):  # (redchi2>1.) and (niter<10):
        result3 = minimize(lnlike_lcrv, fit_params, kws=kws, iter_cb=MinimizeStopper(60), **fit_kws)
        _allres = lnlike_lcrv(result3.params, qua=np.unique(keblat.quarter), polyorder=2, residual=True)
        current_redchi2 = np.sum(_allres ** 2) / len(_allres)
        print "Iteration: ", niter, current_redchi2, result3.redchi
        if current_redchi2 < redchi2:
            print "The following results are saved:"
            report_fit(result3)
            fit_params = result3.params.copy()
            redchi2 = current_redchi2
        niter += 1
    # print "logL of best allpars = ", keblat.lnlike(allpars, lc_constraints=None, qua=np.unique(keblat.quarter), polyorder=2)
    lcrvpars = get_pars2vals(fit_params, partype='lcrv')
    # print fit_params
    return lcrvpars

def estimate_rsum(rsum, period, eclipse_widths, msum=1.0):
    if msum is None:
        msum=rsum
    res = abs(rsum/compute_a(period, msum, unit_au=False) - eclipse_widths)
    return res

def compute_tse(e, w, period, tpe_obs):
    t0 = tpe_obs - sudarsky(np.pi/2. - w, e, period)
    tse = t0 + sudarsky(-np.pi/2. - w, e, period)
    return tse

def estimate_ew(ew, *args):
    period, tpe_obs, tse_obs = args[0], args[1], args[2]
    e, w = ew
    if (e<0) or (e>1):
        return 1e8
    tse = compute_tse(e, w, period, tpe_obs)
    return ((tse % period - tse_obs % period)/0.001)**2

def tse_residuals(ew, *args):
    esinw, ecosw = ew
    e = np.sqrt(esinw**2 + ecosw**2)
    w = np.arctan2(esinw, ecosw)
    if e>1:
        return 1e3
    period, tpe, tse0 = args[0], args[1], args[2]
    tse = tpe - sudarsky(np.pi/2.-w, e, period) + sudarsky(-np.pi/2.-w, e, period)
    return (tse % period - tse0 % period) / 0.01#**2

def flatbottom(x, y, sep, swidth):
    check = (x<sep+swidth/3.) * (x>sep-swidth/3.)
    grady = np.gradient(y)
    grady_m = np.polyfit(x[check], grady[check], 1)[0]
    if abs(grady_m)<0.1:
        return 0.01
    elif abs(grady_m)>10.:
        return 0.4
    else:
        return 0.1

def guess_rrat(sdep, pdep):
    if (pdep>0.2):
        val = sdep/pdep*1.4
        if val>1.:
            return 0.95
        elif val<0.5:
            return 0.7
        return val
    else:
        return np.clip(sdep/pdep, 0.1, 0.95)

def check_lcresiduals(x, y, ey):
    degrees = [2, 5, 9, 13]
    bic = np.zeros(len(degrees))
    for i in range(len(degrees)):
        z = np.poly1d(np.polyfit(x, y, degrees[i]))
        bic[i] = np.sum(((z(x) - y)/ey)**2) + degrees[i]*np.log(len(y))
    bic_slope = np.median(np.diff(bic))/bic[0]
    if bic_slope < -0.1:
        return bic, bic_slope, True
    return bic, bic_slope, False
        
def ilnprob(isopars, lc_constraints=None):
    lp = keblat.ilnprior(isopars)
    if np.isinf(lp):
        return -np.inf, str((0, 0, 0))
    ll, blobs = keblat.ilnlike(isopars, lc_constraints=lc_constraints)
    if (np.isnan(ll) or np.isinf(ll)):
        return -np.inf, str((0, 0, 0))
    return lp + ll, blobs
    

def lnprob(allpars, lc_constraints=None):
    lp = keblat.lnprior(allpars)
    if np.isinf(lp):
        return -np.inf, str((-np.inf, -np.inf, -np.inf))
    ll = keblat.lnlike(allpars, lc_constraints=lc_constraints, qua=np.unique(keblat.quarter))
    if (np.isnan(ll) or np.isinf(ll)):
        return -np.inf, str((-np.inf, -np.inf, -np.inf, -np.inf))
    return lp + ll, str((keblat.r1, keblat.r2, keblat.frat))

def mix_lnprior(allpars):
    m1, m2, z0, age, dist, ebv, h0, period, tpe, esinw, ecosw, \
        b, q1, q2, q3, q4, lnlcerr = allpars[:17]
    Pb, Yb, lnisoerr = allpars[17:]
    e = np.sqrt(esinw**2 + ecosw**2)
    pars2check = np.array([m1, m2, z0, age, dist, ebv, h0, \
        period, tpe, e, b, q1, q2, q3, q4, lnlcerr, Pb, Yb, lnisoerr])
    bounds = np.array([(.1, 12.), (.1, 12.), (0.001, 0.06), (6., 10.1),
                       (10., 15000.), (0.0, 1.0), (119-20., 119+20.), #(10., 15000.), (0.0, 1.0), (119-20., 119+20.),
                        (5., 3000.), (0., 1e8), (0., 0.99), (0., 10.), (0.,1.),
                        (0.,1.), (0.,1.), (0.,1.), (-14, -4.5),
                       (0., 1.), (np.min(keblat.magsobs-2.), np.max(keblat.magsobs)+2.), (-8, 0.)])

    pcheck = np.all((pars2check >= bounds[:,0]) & \
                    (pars2check <= bounds[:,1]))

    if pcheck:
        return 0.0 + age*np.log(10.) + np.log(np.log(10.))
    else:
        return -np.inf

def mix_lnlike2(allpars, polyorder=2):
    m1, m2, z0, age, dist, ebv, h0, period, tpe, esinw, \
        ecosw, b, q1, q2, q3, q4, lcerr = allpars[:17]
    #dist = np.exp(dist)
    Pb, Yb = allpars[17:-1]
    isoerr = np.exp(allpars[-1])
    lcerr = np.exp(lcerr)
    ldcoeffs = np.array([q1, q2, q3, q4])
    keblat.updatepars(m1=m1, m2=m2, z0=z0, age=age, dist=dist, ebv=ebv,
                    h0=h0, period=period, tpe=tpe, esinw=esinw,
                    ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4,
                    lcerr=lcerr, isoerr=isoerr)
    isopars = [m1, m2, z0, age, dist, ebv, h0, isoerr]
    magsmod = keblat.isofit(isopars)


    isores = (magsmod - keblat.magsobs) / keblat.emagsobs

    Lin = 1./(np.sqrt(TWOPI)*keblat.emagsobs) * np.exp(-0.5 * isores**2)
    Lout = 1./np.sqrt(TWOPI * (isoerr**2 + keblat.emagsobs**2)) * \
        np.exp(-0.5 * (keblat.magsobs-Yb)**2 / (isoerr**2+keblat.emagsobs**2))
    lnll = np.sum(np.log((1.-Pb) * Lin + Pb * Lout))
    print lnll
    lcpars = np.concatenate((np.array([m1+m2, keblat.r1+keblat.r2,
                                       keblat.r2/keblat.r1, period,
                                       tpe, esinw, ecosw, b, keblat.frat]),
                                       ldcoeffs))

    clip = keblat.clip
    lcmod, lcpol = keblat.lcfit(lcpars, keblat.jd[clip],
                              keblat.quarter[clip], keblat.flux[clip],
                            keblat.fluxerr[clip], keblat.crowd[clip],
                            polyorder=polyorder)

    lcres = (lcmod*lcpol - keblat.flux[clip]) / np.sqrt(keblat.fluxerr[clip]**2 + lcerr**2)

    if np.any(np.isinf(lcmod)):
        return -np.inf

    lnll += -0.5 * (np.sum(lcres**2) + np.sum(np.log((keblat.fluxerr[clip]**2 + lcerr**2))))
    lnll += -0.5 * ((isopars[5] - keblat.ebv)/(keblat.debv))**2
    #lnll += -0.5 * ((isopars[2] - keblat.z0)/(0.2 * np.log(10) * keblat.z0))**2

    return lnll

def mix_lnlike(allpars, polyorder=2, split=False):
    m1, m2, z0, age, dist, ebv, h0, period, tpe, esinw, \
        ecosw, b, q1, q2, q3, q4, lcerr = allpars[:17]
    #dist = np.exp(dist)
    Pb, Yb = allpars[17:-1]
    isoerr = np.exp(allpars[-1])
    lcerr = np.exp(lcerr)

    ldcoeffs = np.array([q1, q2, q3, q4])
    keblat.updatepars(m1=m1, m2=m2, z0=z0, age=age, dist=dist, ebv=ebv,
                    h0=h0, period=period, tpe=tpe, esinw=esinw,
                    ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4,
                    lcerr=lcerr, isoerr=isoerr)
    isopars = [m1, m2, z0, age, dist, ebv, h0, isoerr]
    magsmod = keblat.isofit(isopars)
    if np.any(np.isinf(magsmod)) or np.isinf(keblat.r1) or np.isinf(keblat.r2):
        return -np.inf #/ np.sqrt(self.emagsobs**2 + isoerr**2)

    isores = (magsmod - keblat.magsobs) / keblat.emagsobs

    # lnll = (1.-Pb)/(np.sqrt(TWOPI)*keblat.emagsobs) * np.exp(-0.5 * isores**2) + \
    #     Pb/np.sqrt(TWOPI * (isoerr**2 + keblat.emagsobs**2)) * \
    #     np.exp(-0.5 * (keblat.magsobs-Yb)**2 / (isoerr**2+keblat.emagsobs**2))
    # lnll = np.sum(np.log(lnll))

    # in case of numerical instabilities with small log sums...
    Lin = -0.5 * isores**2 + np.log((1.-Pb)/(np.sqrt(TWOPI)*keblat.emagsobs))
    Lout = -0.5 * (keblat.magsobs-Yb)**2 / (isoerr**2 + keblat.emagsobs**2) + \
                    np.log(Pb/np.sqrt(TWOPI * (isoerr**2 + keblat.emagsobs**2)))

    lnll = np.logaddexp(Lin, Lout)
    lnll = np.sum(lnll)

    #now the light curve fitting part
    lcpars = np.concatenate((np.array([m1+m2, keblat.r1+keblat.r2,
                                       keblat.r2/keblat.r1, period,
                                       tpe, esinw, ecosw, b, keblat.frat]),
                                       ldcoeffs))

    clip = keblat.clip
    lcmod, lcpol = keblat.lcfit(lcpars, keblat.jd[clip],
                              keblat.quarter[clip], keblat.flux[clip],
                            keblat.fluxerr[clip], keblat.crowd[clip],
                            polyorder=polyorder)

    lcres = (lcmod*lcpol - keblat.flux[clip]) / np.sqrt(keblat.fluxerr[clip]**2 + lcerr**2)

    if np.any(np.isinf(lcmod)):
        return -np.inf

    lnll += -0.5 * (np.sum(lcres**2) + np.sum(np.log((keblat.fluxerr[clip]**2 + lcerr**2))))
    lnll += -0.5 * ((isopars[5] - keblat.ebv)/(keblat.debv))**2
    #lnll += -0.5 * ((isopars[2] - keblat.z0)/(0.2 * np.log(10) * keblat.z0))**2

    if split:
        return lnll, -0.5 * (np.sum(lcres**2) + np.sum(np.log((keblat.fluxerr[clip]**2 + lcerr**2)))) + \
                -0.5 * ((isopars[5] - keblat.ebv)/(keblat.debv))**2
    return lnll

def mix_lnprob(allpars, polyorder=2):
    lp = mix_lnprior(allpars)
    if np.isinf(lp):
        return -np.inf, str((0,0,0))
    ll = mix_lnlike(allpars, polyorder=polyorder)
    if np.isinf(ll) or np.isnan(ll):
        return -np.inf, str((0,0,0))
    return lp+ll, str((keblat.r1, keblat.r2, keblat.frat))


def k_lnprior(lcpars):
    msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4, lcerr = lcpars
    e = np.sqrt(esinw**2 + ecosw**2)
    pars2check = np.array([msum, rsum, rrat, period, tpe, e, b, frat, q1, q2, q3, q4, lcerr])

    pcheck = np.all((pars2check >= lcbounds[:,0]) & (pars2check <= lcbounds[:,1]))

    if pcheck:
        return 0.0
    else:
        return -np.inf

def k_lnlike(lcpars, polyorder=2):
    lcmod, lcpol = keblat.lcfit(lcpars[:-1], keblat.jd[keblat.clip],
                                  keblat.quarter[keblat.clip], keblat.flux[keblat.clip],
                                keblat.fluxerr[keblat.clip], keblat.crowd[keblat.clip],
                                polyorder=polyorder)
    lcerr = np.exp(lcpars[-1])
    lcres = (lcmod*lcpol - keblat.flux[keblat.clip]) / np.sqrt(keblat.fluxerr[keblat.clip]**2 + lcerr**2)

    if np.any(np.isinf(lcmod)):
        return -np.inf

    chisq = np.sum(lcres**2) + np.sum(np.log((keblat.fluxerr[keblat.clip]**2 + lcerr**2)))
    #chisq += ((lcpars[3] - keblat.period)/(0.003*keblat.period))**2
    #chisq += ((lcpars[4] - keblat.tpe)/(0.03*keblat.tpe))**2
    return -0.5 * chisq

def k_lnprob(lcpars, polyorder=2):
    lp = k_lnprior(lcpars)
    if np.isinf(lp):
        return -np.inf, str((0,0,0))
    ll = k_lnlike(lcpars, polyorder=polyorder)
    if np.isnan(ll) or np.isinf(ll):
        return -np.inf, str((0,0,0))
    return lp + ll, str((0,0,0))

def plot_mc(filename, header, footer, nwalkers, ndim, niter, burnin=40000, plot=True, posteriors=False, huber_truths=[],
            isonames=None, iso_extras=False):
    iwalker = np.arange(nwalkers)
    data = np.loadtxt(filename)
    if data.shape[0]/nwalkers < niter/20:
        print "MC file not complete... returning the last ball of walkers (1, nwalkers, ndim)"
        return None, None, None, None, None, False
    afrac = np.empty((data.shape[0]/nwalkers, nwalkers))
    logli, r1, temp1, logg1 = afrac*0., afrac*0., afrac*0., afrac*0.
    params = np.empty((data.shape[0]/nwalkers, nwalkers, len(isonames)))
    strays = []
    for jj in iwalker:
        afrac[:,jj] = data[jj::nwalkers,2]
        logli[:,jj] = data[jj::nwalkers,3]
        r1[:,jj] = data[jj::nwalkers,4]
        temp1[:,jj] = data[jj::nwalkers,5]
        logg1[:,jj] = data[jj::nwalkers,6]
        if len(afrac[:,jj][(afrac[:,jj]<0.1)])>=0.66*len(afrac[:,jj]):
            strays.append(jj)

        for ii in range(len(isonames)):
            params[:, jj, ii] = data[jj::nwalkers, ii+7]

    if plot:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        print "Making plots now."
        fig = plt.figure(figsize=(16, 16))
        for ii in range(len(isonames)):
            ax = fig.add_subplot(int(len(isonames)/2)+1, 2, ii+1)
            ax.plot(params[:, :, ii])
            ax.plot([burnin/10, burnin/10], plt.ylim(), 'y-', lw=2.0)
            ax.set_xlabel('N/10 iteration')
            ax.set_ylabel(isonames[ii])
            divider = make_axes_locatable(ax)
            axhist = divider.append_axes("right", size=1.2, pad=0.1, sharey=ax)
            axhist.hist(params[:,:,ii], 100, histtype='step', alpha=0.6, normed=True,
                        orientation='horizontal')
            axhist.hist(params[:,:,ii].ravel(), 100, histtype='step', color='k',
                        normed=True, orientation='horizontal')
            plt.setp(axhist.get_yticklabels(), visible=False)
        plt.savefig(header+footer+'_parameters.png')

    r1[np.isinf(r1)] = 0.
    temp1[np.isinf(temp1)] = 0.
    logg1[np.isinf(logg1)] = 0.

    mostlike = np.where(logli == np.max(logli))
    mlpars = params[:,:,:][mostlike][0]
    print "Max likelihood out of all samples: ", logli[:,:][mostlike]
    for kk in range(len(isonames)):
        print("""{0} = {1}""".format(str(isonames[kk]), mlpars[kk]))
    if burnin/10>=params.shape[0]:
        print "Burn-in shorter than length of MCMC run, adjusting..."
        burnin = params.shape[0]*3/4*10
    afrac, logli = afrac[burnin/10:,:], logli[burnin/10:,:]
    r1, temp1, logg1 = r1[burnin/10:,:], temp1[burnin/10:,:], logg1[burnin/10:,:]
    params = params[burnin/10:,:,:]

    print "bad/stray walkers =", strays, len(strays)

    keep = iwalker[~np.in1d(iwalker, strays)]
    if len(strays)>=0.33*nwalkers:
        keep = iwalker

    bfpars = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                 zip(*np.percentile(params[:,keep,:].reshape((-1, ndim)),
                                    [16, 50, 84], axis=0)))
    print "MCMC result: "
    print "Accep. Frac = ", np.mean(afrac[:, keep])
    for kk in range(len(isonames)):
        print("""{0} = {1[0]} +{1[1]} -{1[2]}""".format(str(isonames[kk]),
              bfpars[kk]))

    # mostlike = np.where(logli[:,keep] == np.max(logli[:,keep]))
    # mlpars = params[:,keep,:][mostlike][0]
    # print "Max likelihood: ", logli[:,keep][mostlike]
    # for kk in range(len(isonames)):
    #     print("""{0} = {1}""".format(str(isonames[kk]), mlpars[kk]))

    if plot:
        if len(mlpars) == 8:
            iso_pars = mlpars.copy()
            iso_pars[-1] = np.exp(iso_pars[-1])
            make_sed_plots(kic, iso_pars, header+footer, suffix='', savefig=True)
        elif len(mlpars) == 14:
            make_lc_plots(kic, mlpars, header+footer, suffix='', savefig=True, polyorder=2)
        elif len(mlpars) == 18:
            make_sedlc_plots(kic, mlpars, header+footer, suffix='', savefig=True, polyorder=2)
        elif len(mlpars) == 20:
            mlpars_sedlc = mlpars[:-3]
            mlpars_sedlc = np.append(mlpars_sedlc, mlpars[-1])
            make_sedlc_plots(kic, mlpars_sedlc, None, header+footer, suffix='', savefig=True, polyorder=2)
        else:
            print "Not recognized mcmc run"
    if posteriors:
        import corner
        plt.figure(figsize=(14, 14))
        samples = np.concatenate((params[:, :, :], ((r1 + temp1) / (params[:,:,0]+params[:,:,0])**(1./3.))[:,:,np.newaxis],
                              (temp1/r1)[:,:,np.newaxis], logg1[:,:,np.newaxis]), axis=2)
        if iso_extras:
            isonames = isonames + ['rsum', 'rrat', 'frat']
        if len(huber_truths) != len(isonames):
            huber_truths = [None] * len(isonames)
        post_inds = np.arange(len(isonames))
        post_inds = np.delete(post_inds, np.where(np.std(samples, axis=(0,1)) == 0)[0])
        try:
            corner.corner(samples[:, keep,:][:,:,post_inds].reshape((-1, len(post_inds))),
                            labels=np.array(isonames)[post_inds], truths=np.array(huber_truths)[post_inds], truth_color='red')
            plt.savefig(header+footer+'_posteriors.png')
        except Exception, e:
            print str(e)
    return params, r1, temp1, logg1, mlpars, True

def make_p0_ball(p_init, ndim, nwalkers, period_scale=1e-7, mass_scale=1e-4, age_scale=1e-5):
    p0_scale = np.ones(ndim)*1e-4
    p0_scale[[0,1]] = mass_scale
    p0_scale[3] = age_scale
    p0_scale[[7,8]] = period_scale
    p0 = [p_init + p0_scale * p_init * np.random.randn(ndim) for ii in range(nwalkers)]
    p0 = np.array(p0)
    p0[:,6] = 119.
    #p0[:,16] = 1e-4
    p0[:,12:16] = np.clip(p0[:,12:16], 0., 1.0)
    return p0

prefix = 'kics/'+str(kic)+'/'
check_dir_exists(prefix)
keblat.start_errf(prefix+'lcfit.err')

rvdata = np.loadtxt('data/{0}.rv'.format(kic), delimiter=';')

# uncomment the code segment below if want to fit RV
# //load rv data
# //make init guess for masses + K offset
# m1, m2, k0 = keblat.rvprep(rvdata[:,0], rvdata[:,1], rvdata[:,3], rvdata[:,2], rvdata[:,4])
# //run light-curve opt first
# //make sure keblat.pars are updated...
# lcmod, lcpol = keblat.lcfit(opt_lcpars, keblat.jd[keblat.clip].....)
# //update the bounds to make them stricter
# keblat.updatebounds('period', 'tpe', 'esinw', 'ecosw')
# rvpars = [m1+m2, m2/m1, opt_lcpars[3], opt_lcpars[4], opt_lcpars[5], opt_lcpars[6], keblat.pars['inc'], k0, 0]
# //optimize rvparameters using opt_lc + init rv guesses
# opt_rvpars = opt_rv(msum=m1+m2, mrat=m2/m1, period=opt_lcpars[3], tpe=opt_lcpars[4], esinw=opt_lcpars[5],
#                       ecosw=opt_lcpars[6], inc=keblat.pars['inc'], k0=k0, rverr=0)
# //fix msum from rv fit to lc fit
# opt_lcpars[0] = opt_rvpars[0]
# lcpars2 = opt_lc(opt_lcpars, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, keblat.clip, set_upperb = 2.0, vary_msum=False)
# //then optimize both simultaneously
# opt_lcrvpars = opt_lcrv(msum=opt_rvpars[0], mrat=opt_rvpars[1],
#                         rsum=lcpars2[1], rrat=lcpars2[2], period=lcpars2[3],
#                         tpe=lcpars2[4], esinw=lcpars2[5], ecosw=lcpars2[6],
#                         b=lcpars2[7], frat=lcpars2[8], q1=lcpars2[-4],
#                         q2=lcpars2[-3], q3=lcpars2[-2], q4=lcpars2[-1],
#                         lcerr=0.0, k0=opt_rvpars[-2], rverr=0.)

q1, q2, q3, q4 = 0.01, 0.01, 0.01, 0.01
age, h0, dist = 9.2, 119., 850.
chunks = identify_gaps(keblat.cadnum, retbounds_inds=True)
chunks = np.delete(chunks, np.where(np.diff(chunks)<2)[0])
lcchi2_threshold = 3/np.nanmedian(np.array([np.nanmedian(abs(keblat.flux[chunks[ii]:chunks[ii+1]] -
                                                          np.nanmedian(keblat.flux[chunks[ii]:chunks[ii+1]])))
                                          for ii in range(len(chunks)-1)]))
if not os.path.isfile(prefix+'lcpars.lmfit') or clobber_lc:

    # make initial guesses for rsum and f2/f1, assuming main sequence equal mass binary
    rsum = scipy.optimize.fmin_l_bfgs_b(estimate_rsum, 1.0,
                                        args=(period, 2*(keblat.pwidth+keblat.swidth)),
                                        bounds=[(1e-3, 1e3)], approx_grad=True)[0][0]

    # ew = scipy.optimize.fmin(tse_residuals, np.array([1e-3, ecosw]),
    #                          args=(period, tpe, tpe+keblat.sep*period))
    b = flatbottom(keblat.phase[keblat.clip], keblat.flux[keblat.clip], keblat.sep, keblat.swidth)
    # if sdeplist[goodlist][goodlist_ind] < 0.02 and pdeplist[goodlist][goodlist_ind] < 0.04:
    #     b = 1.0
    rrat = guess_rrat(sdeplist[goodlist][goodlist_ind], pdeplist[goodlist][goodlist_ind])
    frat = rrat**(2.5)
    if rsum > 10:
        msum = 2.0
    else:
        msum = rsum
    ew_trials = [[esinw, ecosw], [-esinw, ecosw], [-0.521, ecosw], [-0.332, ecosw], [-0.142, ecosw], [0.521, ecosw], [0.332, ecosw], [0.142, ecosw]]
    lcpars0 = np.array([msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4])
    ew = ew_search_lmfit(ew_trials, lcpars0, (period, tpe, tpe+keblat.sep*period), fit_ecosw=False, polyorder=1)

    b_trials = [0.01, 0.1, 0.4]
    rrat_trials = [0.3, 0.7, 0.95]


    b_trials = [b] + [float(jj) for jj in np.array(b_trials)[~np.in1d(b_trials, b)]]
    rrat_trials = [rrat] + [float(jj) for jj in np.array(rrat_trials)[~np.in1d(rrat_trials, rrat)]]
    lc_search_counts=0
    bestlcchi2 = 1e25

    ###################################################################################
    ########################### LC ONLY OPTIMIZATION FIRST ############################
    ###################################################################################

    # for i_b, i_rrat, i_ew in list(itertools.product(b_trials, rrat_trials, ew_trials)):
    for i_b, i_rrat in list(itertools.product(b_trials, rrat_trials)):
        lcpars0 = np.array([rsum, rsum, i_rrat, period, tpe, ew[0], ew[1], i_b, i_rrat**(2.5),
                            q1, q2, q3, q4])
        upper_b = 2.*i_b if i_b==0.01 else 3.0

        opt_lcpars0 = opt_lc(lcpars0, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, \
                            keblat.clip, set_upperb=upper_b, fit_se=False)

        lcchi2 = np.sum(rez(opt_lcpars0, polyorder=2)**2)/np.sum(keblat.clip)
        if (lcchi2 < bestlcchi2) or (lc_search_counts < 1):
            print "Saving from this run:", lcchi2, bestlcchi2, lc_search_counts
            bestlcchi2 = lcchi2*1.0
            opt_lcpars = opt_lcpars0
        lc_search_counts+=1

        if (bestlcchi2 < 15) and opt_lcpars[2]<=1.0:
            print "These init b, rrat, esinw, ecosw lcpars are: ", i_b, i_rrat, ew
            break

        # opt_lcpars0 = opt_lc(lcpars0, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, \
        #                     keblat.clip, set_upperb=upper_b, prefix=prefix)
        # lcchi2 = np.sum(rez(opt_lcpars0, polyorder=2)**2)/np.sum(keblat.clip)
        # if lcchi2 < bestlcchi2:
        #     bestlcchi2 = lcchi2*1.0
        #     opt_lcpars = opt_lcpars0 * 1.0
        #     make_lc_plots(kic, opt_lcpars0, prefix, polyorder=2, suffix='lc_opt')

    try:
        make_lc_plots(kic, opt_lcpars, prefix, polyorder=2, suffix='lc_opt')
    except Exception, e:
        print str(e)


    if bestlcchi2 < lcchi2_threshold:
        print "Saving lmfit lcpars..."
        np.savetxt(prefix+'lcpars.lmfit', opt_lcpars)
    else:
        print("Bestlcchi2 = {0}, exiting.".format(bestlcchi2))
        #sys.exit()
else:
    print "Loading lcpars lmfit"
    opt_lcpars = np.loadtxt(prefix+'lcpars.lmfit')

lcmod, lcpol = keblat.lcfit(opt_lcpars, keblat.jd, keblat.quarter, keblat.flux, keblat.fluxerr, keblat.crowd, polyorder=0)
#print blah
m1, m2, k0 = keblat.rvprep(rvdata[:,0], rvdata[:,3]*1e3, rvdata[:,1]*1e3, rvdata[:,4]*1e3, rvdata[:,2]*1e3)
keblat.updatepars(m1=m1, m2=m2, msum=m1+m2, mrat=m2/m1)
keblat.updatebounds('period', 'tpe', 'msum', 'mrat')#'esinw', 'ecosw')
rvpars = [m1+m2, m2/m1, opt_lcpars[3], opt_lcpars[4], opt_lcpars[5], opt_lcpars[6], keblat.pars['inc'], k0, 0]
# //optimize rvparameters using opt_lc + init rv guesses
opt_rvpars = opt_rv(msum=m1+m2, mrat=m2/m1, period=opt_lcpars[3], tpe=opt_lcpars[4], esinw=opt_lcpars[5],
                       ecosw=opt_lcpars[6], inc=keblat.pars['inc'], k0=k0, rverr=0)
# //fix msum from rv fit to lc fit
opt_lcpars[0] = opt_rvpars[0]
lcpars2 = opt_lc(opt_lcpars, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, keblat.clip, set_upperb = 2.0, vary_msum=False)
# //then optimize both simultaneously
opt_lcrvpars = opt_lcrv(msum=opt_rvpars[0], mrat=opt_rvpars[1],
                        rsum=lcpars2[1], rrat=lcpars2[2], period=lcpars2[3],
                        tpe=lcpars2[4], esinw=lcpars2[5], ecosw=lcpars2[6],
                        b=lcpars2[7], frat=lcpars2[8], q1=lcpars2[-4],
                        q2=lcpars2[-3], q3=lcpars2[-2], q4=lcpars2[-1],
                        lcerr=0.0, k0=opt_rvpars[-2], rverr=0.)

make_lcrv_plots(kic, opt_lcrvpars, prefix, suffix='lcrv_opt', savefig=True, polyorder=2)
np.savetxt(prefix+'lcrv.lmfit', opt_lcrvpars)

print blah
opt_lcpars = np.append(opt_lcpars, np.log(np.median(abs(np.diff(keblat.flux)))))
if np.isinf(k_lnprior(opt_lcpars)):
    if (np.sqrt(opt_lcpars[5]**2+opt_lcpars[6]**2) > 1.):
        opt_lcpars[5] = ew[0]
        opt_lcpars[6] = ew[1]
        opt_lcpars[3] = period
        opt_lcpars[4] = tpe
    print "Some stuff seems to be out of bounds... somehow... changing them to mean bound values"

###################################################################################
################################## LC ONLY MCMC ###################################
###################################################################################
ndim = len(opt_lcpars)
nwalkers = 64
niter = 20000
header = prefix+'lc_'#postml_'
footer = str(nwalkers)+'x'+str(niter/1000)+'k'
mcfile = header+footer+'.mcmc'

#opt_lcpars[-1] = np.log(np.median(abs(np.diff(keblat.flux))))
p0_scale = np.ones(ndim)*1e-4
p0_scale[3] = 1e-7
p0_scale[4] = 1e-6
p0_scale[[2, 8]] = 1e-5
p0 = [opt_lcpars + p0_scale * opt_lcpars * np.random.randn(ndim) for ii in range(nwalkers)]
p0 = np.array(p0)
isonames = ['msum', 'rsum', 'rrat', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'frat', 'q1', 'q2', 'q3', 'q4', 'lcerr']
success = False
if os.path.isfile(mcfile):
    params, r1, temp1, logg1, mlpars, success = plot_mc(mcfile, header, footer, nwalkers, ndim, niter,
                                                        burnin=niter*3/4, plot=True, posteriors=True,
                                        huber_truths=[None]*(len(isonames)+3), isonames=isonames, iso_extras=False)

if not success or not os.path.isfile(mcfile) or clobber_lc:
    if not os.path.isfile(mcfile) or clobber_lc:
        print "MCMC file does not exist, creating..."
        outf = open(mcfile, "w")
        outf.close()
    else:
        if not success:
            print "MCMC file not complete, appending..."
            #niter = niter-mlpars
            outf = open(mcfile, "a")
            outf.close()
    start_time = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, k_lnprob, threads=4)

    print("Running "+str(niter/1000)+"k MCMC chain")

    for res in sampler.sample(p0, iterations=niter, storechain=False):
        if sampler.iterations % 10 == 0:
            position = res[0]
            outf = open(mcfile, "a")
            for k in range(position.shape[0]):
                blobs = np.array(res[3][k][1:-1].split(","), dtype=float)
                outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(sampler.iterations,
                           k, sampler.acceptance_fraction[k], str(res[1][k]),
                           str(blobs[0]), str(blobs[1]), str(blobs[2]),
                            " ".join([str(ii) for ii in position[k]])))
            outf.close()
        if sampler.iterations % 10000 == 0:
            print "Time Elapsed since niter = ", sampler.iterations, time.time()-start_time

    print "Total Time Elapsed for MCMC Run = ", time.time()-start_time

    print "Tot. Acceptance Fraction = ", np.mean(sampler.acceptance_fraction)
    try:
        print "Tot autocorr time = ", np.mean(sampler.acor)
    except:
        print "Could not compute autocorr time..."

    params, r1, temp1, logg1, mlpars, success = plot_mc(mcfile, header, footer, nwalkers, ndim, niter,
                                                        burnin=niter*3/4, plot=True, posteriors=True,
                                        huber_truths=[None]*(len(isonames)+3), isonames=isonames, iso_extras=False)

###################################################################################
############################# SED ONLY OPTIMIZATION ###############################
###################################################################################

opt_lcpars = mlpars
upper_b = 1.4
if opt_lcpars[-6]*2. > upper_b:
    upper_b = opt_lcpars[-6]*2.

if keblat.magsobs[-1] > 12.:
    dist = 1250.
if keblat.magsobs[-1] > 14.:
    dist = 1500.

lc_constraints = np.array([opt_lcpars[1]/opt_lcpars[0]**(1./3.), opt_lcpars[2], opt_lcpars[8]])
print "LC Constraints (rsum/msum^(1/3), rrat, frat) are: ", lc_constraints

if os.path.isfile(prefix+'isopars.lmfit') and not clobber_sed:
    print "isopars lmfit file already exists, loading..."
    opt_isopars = np.loadtxt(prefix+'isopars.lmfit')
    fit_isopars = opt_sed(opt_isopars, lc_constraints, ebv_dist, ebv_arr, fit_ebv=True, ret_lcpars=False)
else:
    #ebv_trials = [ebv] + list(np.linspace(0.01, 0.1, 4))
    msum_trials = [opt_lcpars[0]] + [1.0, 1.5, 2.0, 2.5, 4.0]
    mrat_trials = [opt_lcpars[2]**(1./0.8)] + [0.3, 0.5, 0.9]
    iso_bestchi2 = 1e25
    iso_counter=0
    for i_msum, i_mrat in list(itertools.product(msum_trials, mrat_trials)):
        m1 = i_msum/(1.+i_mrat)
        m2 = i_msum/(1.+1./i_mrat)
        ### HERE ###
        for i_age in get_age_trials(m1):
            isopars0 = [m1, m2, keblat.z0, i_age, dist, ebv, h0, np.log(0.05)]
            fit_isopars0 = opt_sed(isopars0, lc_constraints, ebv_dist, ebv_arr, fit_ebv=True, ret_lcpars=False)
            isores = keblat.ilnlike(fit_isopars0, lc_constraints=lc_constraints, residual=True)
            iso_redchi2 = np.sum(isores**2) / len(isores)
            if (iso_counter<1) or (iso_redchi2 < iso_bestchi2):
                fit_isopars = fit_isopars0.copy()
                iso_bestchi2 = iso_redchi2*1.0
            iso_counter+=1
        if iso_bestchi2 <= 1.:
            break
    #
    # msum_trials = [opt_lcpars[0]] + [1.0, 1.5, 2.0, 2.5]
    # mrat_trials = opt_lcpars[2] + [0.3, 0.5, 0.9]
    # iso_bestchi2 = 1e25
    # iso_counter=0
    # for i_msum, i_mrat in list(itertools.product(msum_trials, mrat_trials)):
    #     m1 = i_msum/(1.+i_mrat)
    #     m2 = i_msum/(1.+1./i_mrat)
    #     ### HERE ###
    #     for i_age in get_age_trials(m1):
    #         isopars0 = [m1, m2, keblat.z0, i_age, dist, ebv, h0, np.log(0.05)]
    #         fit_isopars0 = minimize(sed_chisq, isopars0, method='nelder-mead', args=(lc_constraints,), options={'xtol':1e-7, 'disp':True})
    #         isores = keblat.ilnlike(fit_isopars0.x, lc_constraints=lc_constraints, residual=True)
    #         iso_redchi2 = np.sum(isores**2) / len(isores)
    #         if (iso_counter<1) or (iso_redchi2 < iso_bestchi2):
    #             fit_isopars = fit_isopars0.copy()
    #             iso_bestchi2 = iso_redchi2*1.0
    #         iso_counter+=1
    #     if iso_bestchi2 <= 1.:
    #         break

    opt_isopars = keblat.ilnlike(fit_isopars, retpars=True)
    np.savetxt(prefix+'isopars.lmfit', opt_isopars)

if np.isinf(ilnprob(opt_isopars, lc_constraints = lc_constraints)[0]):
    print "Somehow isopars are out of bounds? You should double check this, DIana; using initial iso parameters instead"
    opt_isopars = isopars0
###################################################################################
################################## SED ONLY MCMC ##################################
###################################################################################
ndim = len(opt_isopars)
nwalkers = 32
niter = 20000
p0_scale = np.ones(ndim)*1e-4
if (opt_isopars[3] > 10.) or (opt_isopars[0] > 1.):
    p0_scale[3] = 1e-5
p0 = [opt_isopars + p0_scale * opt_isopars * np.random.randn(ndim) for ii in range(nwalkers)]
p0 = np.array(p0)
p0[:,6] = 119.

header = prefix+'sed_'#postml_'
footer = str(nwalkers)+'x'+str(niter/1000)+'k'
mcfile = header+footer+'.mcmc'
isonames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'isoerr']

success = False
if os.path.isfile(mcfile):
    params, r1, temp1, logg1, mlpars, success = plot_mc(mcfile, header, footer, nwalkers, ndim, niter,
                                                        burnin=niter*3/4, plot=True, posteriors=True,
                                        huber_truths=[None]*(len(isonames)) + list(lc_constraints), isonames=isonames, iso_extras=True)

if not success or not os.path.isfile(mcfile) or clobber_sed:
    if not os.path.isfile(mcfile) or clobber_sed:
        print "MCMC file does not exist, creating..."
        outf = open(mcfile, "w")
        outf.close()
    else:
        if not success:
            print "MCMC file not complete, appending..."
            #niter = niter-mlpars
            outf = open(mcfile, "a")
            outf.close()
    start_time = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ilnprob, args=(lc_constraints,), threads=4)

    print("Running "+str(niter/1000)+"k MCMC chain")

    for res in sampler.sample(p0, iterations=niter, storechain=False):
        if sampler.iterations % 10 == 0:
            position = res[0]
            outf = open(mcfile, "a")
            for k in range(position.shape[0]):
                blobs = np.array(res[3][k][1:-1].split(","), dtype=float)
                outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(sampler.iterations,
                           k, sampler.acceptance_fraction[k], str(res[1][k]),
                           str(blobs[0]), str(blobs[1]), str(blobs[2]),
                            " ".join([str(ii) for ii in position[k]])))
            outf.close()
        if sampler.iterations % 10000 == 0:
            print "Time Elapsed since niter = ", sampler.iterations, time.time()-start_time

    print "Total Time Elapsed for MCMC Run = ", time.time()-start_time

    print "Tot. Acceptance Fraction = ", np.mean(sampler.acceptance_fraction)
    try:
        print "Tot autocorr time = ", np.mean(sampler.acor)
    except:
        print "Could not compute autocorr time..."

    params, r1, temp1, logg1, mlpars, success = plot_mc(mcfile, header, footer, nwalkers, ndim, niter,
                                                        burnin=niter*3/4, plot=True, posteriors=True,
                                        huber_truths=[None]*(len(isonames)) + list(lc_constraints), isonames=isonames, iso_extras=True)

###################################################################################
#################################### SEDLC OPT ####################################
###################################################################################
opt_allpars0 = np.concatenate((np.array(mlpars)[:-1], opt_lcpars[3:8], opt_lcpars[9:-1], np.array([1e-5, 0.005])))
from scipy.optimize import least_squares

#bounds = ((0.1, 0.1, 0.001, 6.0, 10., 0., 118., keblat.period-0.1, keblat.tpe-5., opt_lcpars[5]-0.05, opt_lcpars[6]-0.05,  opt_lcpars[7]-0.1, 0., 0., 0., 0., 0., 0.), (12., 12., 0.06, 10.1, 15000., 1.0, 120., keblat.period+0.1, keblat.tpe+5., opt_lcpars[5]+0.05, opt_lcpars[6]+0.05, opt_lcpars[7]+0.1, 1., 1., 1., 1., 0.0001, 0.02))

bounds = ((0.1, 0.1, 0.001, 6.0, 10., 0., 118., opt_lcpars[3]*0.999, opt_lcpars[4]*0.99, opt_lcpars[5]-0.1*abs(opt_lcpars[5]), opt_lcpars[6]-0.05*abs(opt_lcpars[6]), opt_lcpars[7]*0.99, 0., 0., 0., 0., 0., 0.), (12., 12., 0.06, 10.1, 15000., 2.0, 120., opt_lcpars[3]*1.001, opt_lcpars[4]*1.01, opt_lcpars[5]+0.1*abs(opt_lcpars[5]), opt_lcpars[6]+0.05*abs(opt_lcpars[6]), opt_lcpars[7]*1.01, 1., 1., 1., 1., 0.0002, 0.02))

msum_trials = [opt_lcpars[0]] + [1.0, 1.5, 2.0]
mrat_trials = [opt_lcpars[2]**(1/0.8)] + [0.3, 0.5, 0.9]
allpars_bestchi2 = 1e25
opt_all_counter=0
if (keblat.z0 <= bounds[1][2]) and (keblat.z0 >= bounds[0][2]):
    z0 = keblat.z0 * 1.0
else:
    z0 = keblat.zsun

# for i_msum, i_mrat in list(itertools.product(msum_trials, mrat_trials)) + [(opt_allpars0[0]+opt_allpars0[1],
#                                                                             opt_allpars0[1]/opt_allpars0[0])]:
#     m1 = i_msum/(1.+i_mrat)
#     m2 = i_msum/(1.+1./i_mrat)
#     ### HERE ###
#     for i_age in get_age_trials(m1):
#         opt_allpars0[:6] = [m1, m2, z0, i_age, dist, ebv]
#         if m1<0.1 or m1>12. or m2<0.1 or m2>12.:
#             continue
#         print "Init allpars: ", opt_all_counter, opt_allpars0[:6]
#         res_huber = least_squares(lnlike_lmfit, opt_allpars0, bounds=bounds,
#                                   method='trf', loss='huber', f_scale=1., xtol=1e-8,
#                                   kwargs={'residual':True, 'qua': np.unique(keblat.quarter)})
#         allres = lnlike_lmfit(res_huber.x, qua=np.unique(keblat.quarter), polyorder=2, residual=True)
#         allpars_chi2 = np.sum(allres**2) / len(allres)
#         print "made it here"
#         if (opt_all_counter < 1) or (allpars_chi2 < allpars_bestchi2):
#             opt_allpars = res_huber.x*1.
#             allpars_bestchi2 = allpars_chi2
#         opt_all_counter+=1
#     if allpars_bestchi2 < 10:
#         break


opt_all_counter=0
for i_porder, i_b, i_ew, i_za, i_lc in list(itertools.product([1, 2], [True, False], [False, True], [False, True], [lc_constraints, None])):
    opt_allpars0 = opt_sedlc(fit_isopars, opt_lcpars, ebv_dist, ebv_arr, keblat.jd, keblat.phase, keblat.flux,
                        keblat.fluxerr, keblat.crowd, keblat.clip, mciso=mlpars, fit_ebv=True, set_upperb=upper_b,
                            init_porder=i_porder, init_varyb=i_b, init_varyew=i_ew, init_varyza=i_za, lc_constraints=i_lc)
    opt_allpars0 = np.asarray(opt_allpars0)
    allres = lnlike_lmfit(opt_allpars0, lc_constraints=i_lc, qua=np.unique(keblat.quarter), polyorder=2, residual=True)

    allpars_chi2 = np.sum(allres**2) / len(allres)

    if (opt_all_counter < 1) or (allpars_chi2 < allpars_bestchi2):
        opt_allpars = opt_allpars0*1.
        allpars_bestchi2 = allpars_chi2
    if allpars_bestchi2 < 10:
        break


print "Writing and making sedlc optimized plots..."
np.savetxt(prefix+'allpars.lmfit', opt_allpars)
plt.close('all')
if allpars_chi2>=1e10:
    print "The sedlc yielded no good fits. No sedlc_opt plots will be made. "
else:
    make_sedlc_plots(kic, opt_allpars, prefix, suffix='sedlc_opt', savefig=True, polyorder=2)


###################################################################################
##################################### END HERE ####################################
###################################################################################
#
# for i_b, i_rrat, i_ew in list(itertools.product(b_trials, rrat_trials, ew_trials)):
#     lcpars0 = np.array([rsum, rsum, i_rrat, period, tpe, i_ew[0], i_ew[1], i_b, i_rrat**(2.5),
#                         q1, q2, q3, q4])
#     upper_b = 2.*i_b if i_b==0.01 else 1.4
#
#     opt_lcpars0 = opt_lc(lcpars0, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, \
#                         keblat.clip, set_upperb=upper_b, prefix=prefix)
#
#     lcchi2 = np.sum(rez(opt_lcpars0, polyorder=2)**2)/np.sum(keblat.clip)
#     if (lcchi2 < bestlcchi2) or (lc_search_counts < 1):
#         bestlcchi2 = lcchi2*1.0
#         opt_lcpars = opt_lcpars0
#     lc_search_counts+=1
#
#     if bestlcchi2 < 10:
#         print "These init b, rrat, esinw, ecosw lcpars are: ", i_b, i_rrat, i_ew
#         break
#
#     # opt_lcpars0 = opt_lc(lcpars0, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, \
#     #                     keblat.clip, set_upperb=upper_b, prefix=prefix)
#     # lcchi2 = np.sum(rez(opt_lcpars0, polyorder=2)**2)/np.sum(keblat.clip)
#     # if lcchi2 < bestlcchi2:
#     #     bestlcchi2 = lcchi2*1.0
#     #     opt_lcpars = opt_lcpars0 * 1.0
#     #     make_lc_plots(kic, opt_lcpars0, prefix, polyorder=2, suffix='lc_opt')
#
# try:
#     make_lc_plots(kic, opt_lcpars, prefix, polyorder=2, suffix='lc_opt')
# except Exception, e:
#     print str(e)
#
# print "Saving lmfit lcpars..."
#
# np.savetxt(prefix+'lcpars.lmfit', opt_lcpars)
# opt_lcpars = np.append(opt_lcpars, np.log(np.median(abs(np.diff(keblat.flux)))))
# if np.isinf(k_lnprior(opt_lcpars)):
#     if (np.sqrt(opt_lcpars[5]**2+opt_lcpars[6]**2) > 1.):
#         opt_lcpars[5] = ew[0]
#         opt_lcpars[6] = ew[1]
#         opt_lcpars[3] = period
#         opt_lcpars[4] = tpe
#     print "Some stuff seems to be out of bounds... somehow... changing them to mean bound values"
#
# ndim = len(opt_lcpars)
# nwalkers = 64
# niter = 20000
# #opt_lcpars[-1] = np.log(np.median(abs(np.diff(keblat.flux))))
# p0_scale = np.ones(ndim)*1e-4
# p0_scale[3] = 1e-7
# p0_scale[4] = 1e-6
# p0_scale[[2, 8]] = 1e-5
# p0 = [opt_lcpars + p0_scale * opt_lcpars * np.random.randn(ndim) for ii in range(nwalkers)]
# p0 = np.array(p0)
#
# header = prefix+'lc_'#postml_'
# footer = str(nwalkers)+'x'+str(niter/1000)+'k'
# mcfile = header+footer+'.mcmc'
#
# outf = open(mcfile, "w")
# outf.close()
#
# start_time = time.time()
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, k_lnprob, threads=4)
#
# print("Running "+str(niter/1000)+"k MCMC chain")
#
# for res in sampler.sample(p0, iterations=niter, storechain=False):
#     if sampler.iterations % 10 == 0:
#         position = res[0]
#         outf = open(mcfile, "a")
#         for k in range(position.shape[0]):
#             blobs = np.array(res[3][k][1:-1].split(","), dtype=float)
#             outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(sampler.iterations,
#                        k, sampler.acceptance_fraction[k], str(res[1][k]),
#                        str(blobs[0]), str(blobs[1]), str(blobs[2]),
#                         " ".join([str(ii) for ii in position[k]])))
#         outf.close()
#     if sampler.iterations % 10000 == 0:
#         print "Time Elapsed since niter = ", sampler.iterations, time.time()-start_time
#
# print "Total Time Elapsed for MCMC Run = ", time.time()-start_time
#
# print "Tot. Acceptance Fraction = ", np.mean(sampler.acceptance_fraction)
# try:
#     print "Tot autocorr time = ", np.mean(sampler.acor)
# except:
#     print "Could not compute autocorr time..."
#
# isonames = ['msum', 'rsum', 'rrat', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'frat', 'q1', 'q2', 'q3', 'q4', 'lcerr']
# params, r1, temp1, logg1, mlpars = plot_mc(mcfile, header, footer, nwalkers, ndim, burnin=niter*3/4, plot=True, posteriors=True,
#                                     huber_truths=[None]*(len(isonames)+3), isonames=isonames, iso_extras=False)
#
# opt_lcpars = mlpars
# if opt_lcpars[-6]*2. > upper_b:
#     upper_b = opt_lcpars[-6]*2.
#
# m1, m2 = opt_lcpars[0]/(1.+opt_lcpars[2]), opt_lcpars[0]/(1.+1./opt_lcpars[2])
# mrat = m2/m1
# if m1>1.1:
#     m1 = 1.05
#     age = 9.7
# if m2>1.1:
#     m2 = m1*mrat if m1*mrat<1.0 else 0.95
# # if (m2 >= 1.1) and (m2 < 2):
# #     age = 8.5
# # if (m2 > 2):
# #     age = 7.1
# if keblat.magsobs[-1] > 12.:
#     dist = 1250.
# if keblat.magsobs[-1] > 14.:
#     dist = 1500.
#
# isopars0 = [m1, m2, keblat.z0, age, dist, ebv, h0, np.log(0.05)]
#
# lc_constraints = np.array([opt_lcpars[1]/opt_lcpars[0]**(1./3.), opt_lcpars[2], opt_lcpars[-5]])
# print "LC Constraints (rsum/msum^(1/3), rrat, frat) are: ", lc_constraints
# fit_isopars, opt_r1, opt_r2, opt_frat = opt_sed(isopars0, lc_constraints,
#                                                 ebv_dist, ebv_arr, fit_ebv=True,
#                                                 ret_lcpars=True, prefix=prefix)
#
# opt_isopars = keblat.ilnlike(fit_isopars, retpars=True)
# np.savetxt(prefix+'isopars.lmfit', opt_isopars)
#
# if np.isinf(ilnprob(opt_isopars, lc_constraints = lc_constraints)[0]):
#     print "Somehow isopars are out of bounds? You should double check this, DIana; using initial iso parameters instead"
#     opt_isopars = isopars0
#
# ndim = len(opt_isopars)
# nwalkers = 32
# niter = 20000
# p0_scale = np.ones(ndim)*1e-4
# if (opt_isopars[3] > 10.) or (opt_isopars[0] > 1.):
#     p0_scale[3] = 1e-5
# p0 = [opt_isopars + p0_scale * opt_isopars * np.random.randn(ndim) for ii in range(nwalkers)]
# p0 = np.array(p0)
# p0[:,6] = 119.
#
# header = prefix+'sed_'#postml_'
# footer = str(nwalkers)+'x'+str(niter/1000)+'k'
# mcfile = header+footer+'.mcmc'
#
# outf = open(mcfile, "w")
# outf.close()
#
# start_time = time.time()
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, ilnprob, args=(lc_constraints,), threads=4)
#
# print("Running "+str(niter/1000)+"k MCMC chain")
#
# for res in sampler.sample(p0, iterations=niter, storechain=False):
#     if sampler.iterations % 10 == 0:
#         position = res[0]
#         outf = open(mcfile, "a")
#         for k in range(position.shape[0]):
#             blobs = np.array(res[3][k][1:-1].split(","), dtype=float)
#             outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(sampler.iterations,
#                        k, sampler.acceptance_fraction[k], str(res[1][k]),
#                        str(blobs[0]), str(blobs[1]), str(blobs[2]),
#                         " ".join([str(ii) for ii in position[k]])))
#         outf.close()
#     if sampler.iterations % 10000 == 0:
#         print "Time Elapsed since niter = ", sampler.iterations, time.time()-start_time
#
# print "Total Time Elapsed for MCMC Run = ", time.time()-start_time
#
# print "Tot. Acceptance Fraction = ", np.mean(sampler.acceptance_fraction)
# try:
#     print "Tot autocorr time = ", np.mean(sampler.acor)
# except:
#     print "Could not compute autocorr time..."
#
# isonames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'isoerr']
# params, r1, temp1, logg1, mlpars = plot_mc(mcfile, header, footer, nwalkers, ndim, burnin=niter*3/4, plot=True, posteriors=True,
#                                     huber_truths=[None]*(len(isonames)) + list(lc_constraints), isonames=isonames, iso_extras=True)
#
# # now optimize all parameters
# opt_sedlc_user_controls = zip(np.array([1, 1, 2, 2]), [False, True, False, True])
# opt_lnlili = -1e25
# for ii in range(len(opt_sedlc_user_controls)):
#     opt_allpars0 = opt_sedlc(fit_isopars, opt_lcpars, ebv_dist, ebv_arr, keblat.jd, keblat.phase, keblat.flux,
#                         keblat.fluxerr, keblat.crowd, keblat.clip, mciso=mlpars, fit_ebv=True, set_upperb=upper_b,
#                             init_porder=opt_sedlc_user_controls[ii][0], init_varyb=opt_sedlc_user_controls[ii][1])
#
#     opt_lnlili0 = keblat.lnlike(opt_allpars0, lc_constraints=lc_constraints,
#                            qua=np.unique(keblat.quarter), polyorder=2)
#     if (ii==0) or (opt_lnlili0 > opt_lnlili):
#         opt_allpars = opt_allpars0*1.
#         opt_lnlili = opt_lnlili0*1.
#
# print "Writing and making sedlc optimized plots..."
# np.savetxt(prefix+'allpars.lmfit', opt_allpars)
# plt.close('all')
# if np.isinf(opt_lnlili):
#     print "The sedlc yielded no good fits. No sedlc_opt plots will be made. "
# else:
#     make_sedlc_plots(kic, opt_allpars, prefix, suffix='sedlc_opt', savefig=True, polyorder=2)

# ###################### pruning outliers ####################
# opt_allpars[-2] = 1e-4
# opt_allpars = np.append(opt_allpars[:-1], np.array([1e-3, np.mean(keblat.magsobs), opt_allpars[-1]]))
# ndim = len(opt_allpars)
# nwalkers = 64
# niter = 50000
#
# p0 = [opt_allpars + 1e-4 * opt_allpars * np.random.randn(ndim) for ii in range(nwalkers)]
# p0 = np.array(p0)
# p0[:,6] = 119.
# #p0[:,16] = 1e-4
# p0[:,12:16] = np.clip(p0[:,12:16], 0., 1.0)
#
# header = prefix+'sedlc_outliers_'#postml_'
# footer = str(nwalkers)+'x'+str(niter/1000)+'k'
# mcfile = header+footer+'.mcmc'
#
# outf = open(mcfile, "w")
# outf.write("""#{0} \n""".format(" ".join([str(ii) for ii in opt_allpars])))
# outf.close()
#
# start_time = time.time()
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, mix_lnprob, args=(2,), threads=4)
#
# print("Running "+str(niter/1000)+"k MCMC chain")
#
# for res in sampler.sample(p0, iterations=niter, storechain=False):
#     if sampler.iterations % 10 == 0:
#         position = res[0]
#         outf = open(mcfile, "a")
#         for k in range(position.shape[0]):
#             blobs = np.array(res[3][k][1:-1].split(","), dtype=float)
#             outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(sampler.iterations,
#                        k, sampler.acceptance_fraction[k], str(res[1][k]),
#                        str(blobs[0]), str(blobs[1]), str(blobs[2]),
#                         " ".join([str(ii) for ii in position[k]])))
#         outf.close()
#     if sampler.iterations % 10000 == 0:
#         print "Time Elapsed since niter = ", sampler.iterations, time.time()-start_time
#
# print "Total Time Elapsed for MCMC Run = ", time.time()-start_time
#
# print "Tot. Acceptance Fraction = ", np.mean(sampler.acceptance_fraction)
# try:
#     print "Tot autocorr time = ", np.mean(sampler.acor)
# except:
#     print "Could not compute autocorr time..."
#
#
# isonames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'q1', 'q2',
#             'q3', 'q4', 'lcerr', 'Pb', 'Yb', 'isoerr']
#
# params, r1, temp1, logg1, mlpars = plot_mc(mcfile, header, footer, nwalkers, ndim, burnin=35000, plot=True, posteriors=True,
#                                     huber_truths=[None]*(len(isonames)+3), isonames=isonames)
#
# plt.close('all')
#
# # ==============================================================================
# # end of code comment out
# # ==============================================================================
#
#
# ###################### pruning outliers ####################
# opt_allpars = np.loadtxt(prefix+'allpars.lmfit')
# if opt_allpars[-1] > 0.:
#     opt_allpars[-1] = -3.5
# opt_allpars[4] = np.log(opt_allpars[4])
# opt_allpars[-2] = 1e-4
# opt_allpars = np.append(opt_allpars[:-1], np.array([1e-3, np.mean(keblat.magsobs), opt_allpars[-1]]))
# ndim = len(opt_allpars)
# nwalkers = 64
# niter = 50000
#
# p0 = [opt_allpars + 1e-4 * opt_allpars * np.random.randn(ndim) for ii in range(nwalkers)]
# p0 = np.array(p0)
# p0[:,6] = 119.
# #p0[:,16] = 1e-4
# p0[:,12:16] = np.clip(p0[:,12:16], 0., 1.0)
#
# header = prefix+'sedlc_outliers2_'#postml_'
# footer = str(nwalkers)+'x'+str(niter/1000)+'k'
# mcfile = header+footer+'.mcmc'
#
# outf = open(mcfile, "w")
# outf.write("""#{0} \n""".format(" ".join([str(ii) for ii in opt_allpars])))
# outf.close()
#
# start_time = time.time()
# import emcee
# sampler = emcee.EnsembleSampler(nwalkers, ndim, mix_lnprob, args=(2,), threads=4)
#
# print("Running "+str(niter/1000)+"k MCMC chain")
#
# for res in sampler.sample(p0, iterations=niter, storechain=False):
#     if sampler.iterations % 10 == 0:
#         position = res[0]
#         outf = open(mcfile, "a")
#         for k in range(position.shape[0]):
#             blobs = np.array(res[3][k][1:-1].split(","), dtype=float)
#             outf.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(sampler.iterations,
#                        k, sampler.acceptance_fraction[k], str(res[1][k]),
#                        str(blobs[0]), str(blobs[1]), str(blobs[2]),
#                         " ".join([str(ii) for ii in position[k]])))
#         outf.close()
#     if sampler.iterations % 10000 == 0:
#         print "Time Elapsed since niter = ", sampler.iterations, time.time()-start_time
#
# print "Total Time Elapsed for MCMC Run = ", time.time()-start_time
#
# print "Tot. Acceptance Fraction = ", np.mean(sampler.acceptance_fraction)
# try:
#     print "Tot autocorr time = ", np.mean(sampler.acor)
# except:
#     print "Could not compute autocorr time..."
