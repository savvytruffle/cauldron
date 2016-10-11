import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import sys, os, itertools, glob
import lmfit
if lmfit.__version__[2] != '9':
    print "Version >= 0.9 of lmfit required..."
    sys.exit()
from lmfit import minimize, Parameters, report_fit
from matplotlib.ticker import MultipleLocator
from keblat import *
from scipy.optimize import fmin_l_bfgs_b
from helper_funcs import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

kic = int(sys.argv[1])
keblat = Keblat(preload=False)
keblat.loadiso2()
keblat.loadsed(sedfile='data/kepsedall_0216_full_test.dat')
keblat.loadvkeb(filename='data/kebproperties_0216_full.dat')
goodv = keblat.kiclookup(kic, target=keblat.vkeb[:, 0])
keblat.loadlc(kic, keblat.vkeb[goodv, [1, 2, 5, 6, 7]])
try:
    exclusive = excludelist[kic]
except:
    exclusive = []
magsobs, emagsobs, extinction, glat, z0 = keblat.getmags(kic)
ebv = extinction[0]
keblat.isoprep(magsobs, emagsobs, extinction, glat, z0, exclude=exclusive)
period, tpe = keblat.vkeb[goodv, 1], keblat.vkeb[goodv, 2]
ecosw, esinw = keblat.vkeb[goodv, -2], keblat.vkeb[goodv, -1]
frat = (keblat.vkeb[goodv, 4]/keblat.vkeb[goodv, 3])


kiclist, perlist, pdeplist, sdeplist, morphlist = np.loadtxt('data/kebproperties_0216_full.dat',
                                          usecols=(0, 1, 3, 4, 8), unpack=True, delimiter=';')
goodlist = (perlist>0)
goodlist_ind = np.where(kiclist[goodlist].astype(int) == kic)[0]
q1, q2, q3, q4 = 0.01, 0.01, 0.01, 0.01

if keblat.swidth < 0.:
    print "No secondary eclipses detected. Exiting."
    sys.exit()


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
        parnames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'isoerr', 'period', 'tpe',
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
    if len(guess) > 13:
        crowd_fits = keblat.broadcast_crowd(keblat.quarter, guess[-len(np.unique(keblat.quarter)):])
    else:
        crowd_fits = keblat.crowd
    lcmod, lcpol = keblat.lcfit(guess[:13], keblat.jd[keblat.clip],
                                keblat.quarter[keblat.clip], keblat.flux[keblat.clip],
                                keblat.fluxerr[keblat.clip], crowd_fits[keblat.clip],
                                polyorder=polyorder)
    if np.any(np.isinf(lcmod)):
        return np.ones(np.sum(keblat.clip)) * 1e10  # (1.-keblat.flux[keblat.clip])/keblat.dflux[keblat.clip]
    return (keblat.flux[keblat.clip] - lcmod * lcpol) / keblat.fluxerr[keblat.clip]

def estimate_rsum(rsum, period, eclipse_widths, msum=1.0):
    if msum is None:
        msum = rsum
    res = abs(rsum / compute_a(period, msum, unit_au=False) - eclipse_widths)
    return res

def opt_lc(lcpars0, jd, phase, flux, dflux, crowd, clip, set_upperb=2., fit_crowd=False,
           vary_msum=True):
    msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4 = lcpars0

    fit_params = Parameters()
    fit_params.add('esinw', value=esinw, min=-.999, max=0.999, vary=False)
    fit_params.add('ecosw', value=ecosw, min=-.999, max=0.999, vary=False)  # ecosw-0.05, max=ecosw+0.05, vary=False)
    fit_params.add('rsum', value=rsum, min=0.1, max=10000., vary=False)
    fit_params.add('rrat', value=rrat, min=1e-4, max=1., vary=False)
    fit_params.add('b', value=b, min=0., max=set_upperb, vary=False)
    fit_params.add('frat', value=frat, min=1e-6, max=1., vary=False)
    fit_params.add('msum', value=msum, min=0.2, max=24., vary=False)

    fit_params.add('period', value=period, min=period - 0.005, max=period + 0.005, vary=False)
    fit_params.add('tpe', value=tpe, min=tpe - 10., max=tpe + 10., vary=False)
    fit_params.add('q1', value=q1, min=0., max=1., vary=False)
    fit_params.add('q2', value=q2, min=0., max=1., vary=False)
    fit_params.add('q3', value=q3, min=0., max=1., vary=False)
    fit_params.add('q4', value=q4, min=0., max=1., vary=False)

    fit_params['rrat'].vary = True
    fit_params['rsum'].vary = True
    fit_params['b'].vary = True
    fit_params['frat'].vary = True
    fit_params['esinw'].vary = True
    fit_params['ecosw'].vary = True
    fit_params['tpe'].vary = True

    fit_kws = {'maxfev': 100 * (len(fit_params) + 1)}

    if fit_crowd:
        print "Fitting crowding parameters..."
        for ii in fit_params.keys():
            fit_params[ii].vary = True
        for ii in range(len(np.unique(keblat.quarter))):
            fit_params.add('cr' + str(np.unique(keblat.quarter)[ii]), value=keblat._crowdsap[ii], min=0.1, max=1.0)
        result0 = minimize(rez, fit_params, kws={'polyorder': 2}, iter_cb=MinimizeStopper(10), **fit_kws)
        report_fit(result0)

        for ii in result0.params.keys():
            result0.params[ii].vary = True
        niter = 0
        redchi2 = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=2)) ** 2) / np.sum(
            keblat.clip)
        guess = get_pars2vals(result0.params, partype='lc')
        while (redchi2 > 1.) and (niter < 5):
            result0 = minimize(rez, result0.params, kws={'polyorder': 2}, iter_cb=MinimizeStopper(10), **fit_kws)
            current_chi = np.sum(
                (rez(get_pars2vals(result0.params, partype='lc'), polyorder=2)) ** 2) / np.sum(
                keblat.clip)
            if current_chi < redchi2:
                redchi2 = current_chi * 1.0
                report_fit(result0)
                guess = get_pars2vals(result0.params, partype='lc')
            niter += 1

        return guess

    print "=========================================================================="
    print "==================== Starting LIGHTCURVE ONLY fit... ====================="
    print "=========================================================================="

    result0 = minimize(rez, fit_params, kws={'polyorder': 1}, iter_cb=MinimizeStopper(10), **fit_kws)
    report_fit(result0)

    fit_params = result0.params
    redchi2 = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=1)) ** 2) / np.sum(
        keblat.clip)

    fit_params['msum'].vary = vary_msum
    fit_params['tpe'].vary = True
    fit_params['period'].vary = True
    fit_params['b'].vary = True
    fit_params['frat'].vary = True
    fit_params['esinw'].vary = True
    fit_params['ecosw'].vary = True
    fit_params['rsum'].vary = True
    fit_params['rrat'].vary = True
    fit_params['q1'].vary = True
    fit_params['q2'].vary = True
    fit_params['q3'].vary = True
    fit_params['q4'].vary = True

    result0 = minimize(rez, fit_params, kws={'polyorder': 1}, iter_cb=MinimizeStopper(10), **fit_kws)
    #    current_redchi = np.sum((result0.residual)**2) / (len(result0.residual)-result0.nfev)
    current_redchi = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=1)) ** 2) / np.sum(
        keblat.clip)

    if current_redchi < redchi2:
        redchi2 = current_redchi * 1.
        fit_params = result0.params
        print "polyorder = 1: ", current_redchi, result0.redchi
        report_fit(result0)

    result0 = minimize(rez, fit_params, kws={'polyorder': 2}, iter_cb=MinimizeStopper(10), **fit_kws)
    current_redchi = np.sum((rez(get_pars2vals(result0.params, partype='lc'), polyorder=2)) ** 2) / np.sum(
        keblat.clip)
    if current_redchi < redchi2:
        redchi2 = current_redchi * 1.
        fit_params = result0.params
        print "polyorder = 2: ", current_redchi, result0.redchi
        report_fit(result0)

    guess = get_pars2vals(fit_params, partype='lc')
    return guess

def flatbottom(x, y, sep, swidth):
    check = (x < sep + swidth / 3.) * (x > sep - swidth / 3.)
    grady = np.gradient(y)
    grady_m = np.polyfit(x[check], grady[check], 1)[0]
    if abs(grady_m) < 0.1:
        return 0.01
    elif abs(grady_m) > 10.:
        return 0.4
    else:
        return 0.1

def guess_rrat(sdep, pdep):
    if (pdep > 0.2):
        val = sdep / pdep * 1.4
        if val > 1.:
            return 0.95
        elif val < 0.5:
            return 0.7
        return val
    else:
        return sdep / pdep

def ew_search_lmfit(ew_trials, pars0, argpars, fit_ecosw=True):
    def ew_search(ew, _pars0=None, _argpars=None):
        def tse_residuals(ew, *args):
            esinw, ecosw = ew
            e = np.sqrt(esinw ** 2 + ecosw ** 2)
            w = np.arctan2(esinw, ecosw)
            if e > 1:
                return 1e3
            period, tpe, tse0 = args[0], args[1], args[2]
            tse = tpe - sudarsky(np.pi / 2. - w, e, period) + sudarsky(-np.pi / 2. - w, e, period)
            return (tse % period - tse0 % period) / 0.01  # **2

        pars = np.array(_pars0).copy()
        pars[5], pars[6] = ew['esinw'].value, ew['ecosw'].value
        mod, _ = keblat.lcfit(pars, keblat.jd, keblat.quarter, keblat.flux, keblat.fluxerr, keblat.crowd,
                              polyorder=0)
        if np.any(np.isinf(mod)):
            return np.ones_like(keblat.jd + 1) * 1e10
        return np.append((keblat.flux - mod) / keblat.fluxerr, tse_residuals((pars[5], pars[6]), *_argpars))

    fit_ew = Parameters()
    fit_ew.add('esinw', value=pars0[5], min=-0.9, max=0.9, vary=True)
    fit_ew.add('ecosw', value=pars0[6], min=-0.9, max=0.9, vary=fit_ecosw)
    chisq = 1e18
    for ii in range(len(ew_trials)):
        fit_ew['esinw'].value = ew_trials[ii][0]
        fit_ew['ecosw'].value = ew_trials[ii][1]
        result = minimize(ew_search, fit_ew, kws={'_pars0': pars0, '_argpars': argpars})
        if result.redchi < chisq or ii == 0:
            chisq = result.redchi * 1.0
            ew_best = result.params['esinw'].value, result.params['ecosw'].value
            print "Better redchi: ", chisq, result.redchi, ew_best
    return ew_best


class Gloop(object):
    def __init__(self, kic, interactive=True):
        # self.ax = ax
        # self.fig = ax.figure
        self.kic = kic
        self.fitlc = False
        self.fitsed = False
        self.fitrv = False
        self.lcpars = None
        self.rvpars = None
        self.period = period
        self.tpe = tpe
        self.ecosw, self.esinw = ecosw, esinw
        self.sdepth = sdeplist[goodlist][goodlist_ind]
        self.pdepth = pdeplist[goodlist][goodlist_ind]
        self.swidth = keblat.swidth
        self.pwidth = keblat.pwidth
        self.sep = keblat.sep
        #keblat = Keblat(preload=False)

        self.tolerance=10
        check_dir_exists('kics/'+str(kic))
        self.prefix = 'kics/'+str(kic)+'/'
        keblat.start_errf(self.prefix+'lc.err')
        #self.norm = plt.Normalize(0)
        if interactive:
            self.x, self.xx = [], []
            self.y, self.yy = [], []
            self.fig, self.ax, self.ax2l, self.ax2r = self.setup_axes()
            #self.points = self.ax.scatter([], [], s=200, color='red', marker='o', picker=self.tolerance)
            #self.lines = self.ax.plot([], [], color='red', linestyle='-', lw=2, alpha=0.5, picker=self.tolerance)
            print("===================================================="
                  "============= Entering interactive mode ============"
                  "====================================================")
            fitcat = raw_input("What component would you like to fit? Enter 'lc', 'rv', or 'sed': ")
            if isinstance(fitcat, list):
                if 'lc' in fitcat:
                    self.fitlc = True
                elif 'rv' in fitcat:
                    self.fitrv = True
                elif 'sed' in fitcat:
                    self.fitsed = True
                else:
                    print("Certain aspects of your user input not recognized. Try again 'lc', 'rv', or 'sed': ")
            else:
                if fitcat == 'lc':
                    self.fitlc = True
                elif fitcat == 'rv':
                    self.fitlc = True
                elif fitcat == 'sed':
                    self.fitsed = True

            if self.fitlc:
                self.inter_lc()

    def connect_mpl(self):
        connect = self.fig.canvas.mpl_connect
        #connect('button_press_event', self.on_click)
        self.cid_key = connect('key_press_event', self.on_key)
        self.cid_pick = connect('pick_event', self.on_pick)

            #connect('motion_notify_event', self.on_hover)

    def disconnect_mpl(self):
        self.fig.canvas.mpl_disconnect(self.cid_key)
        self.fig.canvas.mpl_disconnect(self.cid_pick)

    @staticmethod
    def pause():
        programPause = raw_input("Press the <ENTER> key to continue...")

    def inter_lc(self):
        print("===================================================="
              "============= Light curve fitting / opt ============"
              "====================================================")
        guess_first = raw_input("Make initial fit based on auto. gen guesses? 'y' to continue, 'n' to make manual user inputs: ")
        while guess_first == 'n':
            self.connect_mpl()
            print("Estimate the period and time of primary eclipse by selecting 2 - 3 consecutive centers of primary eclipses.")
            self.pause()
            self.period = np.median(np.diff(self.x, axis=0))
            self.tpe = self.x[0][0]
            print("Period = {0} d, time of primary eclipse = {1}".format(self.period, self.tpe))
            #self.ax2l.plot((keblat.jd - self.tpe)%self.period / self.period, keblat.flux, 'k.')
            self.clear_xy()
            print("Estimate the primary eclipse duration and depth by selecting 'x' edges and 'y' edges")
            self.pause()
            self.pwidth = np.median(np.diff(self.x, axis=0)[::2]) / 2. / self.period
            self.pdepth = np.median(abs(np.diff(self.y, axis=0))[::2])
            print("pdep = {0}, pwidth = {1} d = {2} phase".format(self.pdepth, self.pwidth * self.period, self.pwidth))
            keblat.pwidth=self.pwidth
            keblat.updatephase(self.tpe, self.period)
            pe = (keblat.phase[keblat.clip] >= -1.2 * keblat.pwidth) * (keblat.phase[keblat.clip] <= 1.2 * keblat.pwidth)
            self.ax2l.plot(keblat.phase[keblat.clip][pe], keblat.flux[keblat.clip][pe], 'k.')
            self.ax2l.set_ylim((keblat.flux.min(), keblat.flux.max()))
            self.ax2l.set_xlim((keblat.phase[keblat.clip][pe].min(), keblat.phase[keblat.clip][pe].max()))
            self.clear_xy()
            print("Estimate the time of secondary eclipse by selecting the center of a secondary eclipse.")
            self.pause()
            self.sep = (self.x[0][0] - self.tpe) % self.period / self.period
            print("tse = {0} d, sep = {1} phase".format(self.x[0][0], self.sep))
            keblat.sep = self.sep
            self.clear_xy()
            print("Estimate the secondary eclipse duration and depth by selecting 'x' edges and 'y' edges")
            self.pause()
            self.swidth = np.median(np.diff(self.x, axis=0)[::2]) / 2. / self.period
            self.sdepth = np.median(abs(np.diff(self.y, axis=0))[::2])
            keblat.swidth=self.swidth
            print("sdep = {0}, swidth = {1} d = {2} phase".format(self.sdepth, self.swidth * self.period, self.swidth))
            self.clear_xy()
            se = (keblat.phase[keblat.clip] >= -1.2 * keblat.swidth + keblat.sep) * (keblat.phase[keblat.clip] <= 1.2 * keblat.swidth + keblat.sep)
            self.ax2r.plot(keblat.phase[keblat.clip][se], keblat.flux[keblat.clip][se], 'k.')
            self.ax2r.set_ylim((keblat.flux.min(), keblat.flux.max()))
            self.ax2r.set_xlim((keblat.phase[keblat.clip][se].min(), keblat.phase[keblat.clip][se].max()))

            self.ecosw = (self.sep * 2 - 1.) / np.pi/4.
            self.esinw = (self.pwidth / self.swidth - 1.) / (self.pwidth / self.swidth + 1.)
            print("esinw = {0}, ecosw = {1}".format(self.esinw, self.ecosw))
            self.clear_xy()
            self.disconnect_mpl()
            guess_first = raw_input("Continue? Enter 'n' to re-fit, any other key to continue")

        rsum = fmin_l_bfgs_b(estimate_rsum, 1.0,
                             args=(self.period, 2 * (keblat.pwidth + keblat.swidth)),
                             bounds=[(1e-3, 1e3)], approx_grad=True)[0][0]
        b = flatbottom(keblat.phase[keblat.clip], keblat.flux[keblat.clip], keblat.sep, keblat.swidth)
        if self.sdepth < 0.05 and self.pdepth < 0.05:
            b = 1.0
        rrat = guess_rrat(self.sdepth, self.pdepth)
        frat = rrat ** (2.5)

        ew_trials = [[self.esinw, self.ecosw], [-self.esinw, self.ecosw], [-0.521, self.ecosw],
                     [-0.332, self.ecosw], [-0.142, self.ecosw], [0.521, self.ecosw], [0.332, self.ecosw],
                     [0.142, self.ecosw]]
        lcpars0 = np.array([rsum, rsum, rrat, self.period, self.tpe, self.esinw, self.ecosw, 1.0, frat, q1, q2, q3, q4])
        ew = ew_search_lmfit(ew_trials, lcpars0, (self.period, self.tpe, self.tpe + keblat.sep * self.period),
                             fit_ecosw=False)

        b_trials = [0.01, 0.1, 0.4]
        rrat_trials = [0.3, 0.7, 0.95]

        b_trials = [b] + [float(jj) for jj in np.array(b_trials)[~np.in1d(b_trials, b)]]
        rrat_trials = [rrat] + [float(jj) for jj in np.array(rrat_trials)[~np.in1d(rrat_trials, rrat)]]
        lc_search_counts = 0
        bestlcchi2 = 1e25
        for i_b, i_rrat in list(itertools.product(b_trials, rrat_trials)):
            lcpars0 = np.array([rsum, rsum, i_rrat, self.period, self.tpe, ew[0], ew[1], i_b, i_rrat ** (2.5),
                                q1, q2, q3, q4])
            upper_b = 2. * i_b if i_b == 0.01 else 3.0

            opt_lcpars0 = opt_lc(lcpars0, keblat.jd, keblat.phase, keblat.flux, keblat.fluxerr, keblat.crowd, \
                                 keblat.clip, set_upperb=upper_b)

            lcchi2 = np.sum(rez(opt_lcpars0, polyorder=2) ** 2) / np.sum(keblat.clip)
            if (lcchi2 < bestlcchi2) or (lc_search_counts < 1):
                # print "Saving from this run:", lcchi2, bestlcchi2, lc_search_counts
                bestlcchi2 = lcchi2 * 1.0
                opt_lcpars = opt_lcpars0
            lc_search_counts += 1

            if (bestlcchi2 < 15) and opt_lcpars[2] <= 1.0:
                # print "These init b, rrat, esinw, ecosw lcpars are: ", i_b, i_rrat, ew
                break
        self.lcpars = opt_lcpars

        self.overplot_lcfit(self.lcpars, polyorder=2)
        #
        # try:
        #     make_lc_plots(kic, opt_lcpars, prefix, polyorder=2, suffix='lc_opt')
        # except Exception, e:
        #     print str(e)

    def overplot_lcfit(self, lcpars, polyorder=2):
        keblat.updatephase(lcpars[4], lcpars[3])
        lcmod, lcpol = keblat.lcfit(lcpars[:13], keblat.jd[keblat.clip], keblat.quarter[keblat.clip],
                                    keblat.flux[keblat.clip], keblat.fluxerr[keblat.clip],
                                    keblat.crowd[keblat.clip], polyorder=polyorder)
        lcres = keblat.flux[keblat.clip] - lcmod * lcpol
        if len(lcpars) == 14:
            lcerr = lcpars[-1]
            if lcerr<0: lcerr=np.exp(lcerr)
        else:
            lcerr=0.
        pe = (keblat.phase[keblat.clip] >= -1.2*keblat.pwidth) * (keblat.phase[keblat.clip] <= 1.2*keblat.pwidth)
        se = (keblat.phase[keblat.clip] >= -1.2*keblat.swidth+keblat.sep) * (keblat.phase[keblat.clip] <= 1.2*keblat.swidth+keblat.sep)

        self.ax.clear()
        self.ax2l.clear()
        self.ax2r.clear()

        self.ax.errorbar(keblat.jd, keblat.flux, keblat.fluxerr, fmt='k.', ecolor='gray')
        self.ax.plot(keblat.jd[keblat.clip], lcmod*lcpol, 'r.')
        divider = make_axes_locatable(self.ax)
        axb = divider.append_axes("bottom", size=1.0, pad=0., sharex=self.ax)
        axb.errorbar(keblat.jd[keblat.clip], lcres, keblat.dflux[keblat.clip], fmt='k.', ecolor='gray')

        divider = make_axes_locatable(self.ax2l)
        ax2lb = divider.append_axes("bottom", size=1.0, pad=0., sharex=self.ax2l)

        self.ax2l.errorbar(keblat.phase[keblat.clip][pe], keblat.flux[keblat.clip][pe] / lcpol[pe],
                     keblat.dflux[keblat.clip][pe], fmt='k.', ecolor='gray')
        self.ax2l.plot(keblat.phase[keblat.clip][pe], lcmod[pe], 'r.')
        ax2lb.errorbar(keblat.phase[keblat.clip][pe], lcres[pe], keblat.dflux[keblat.clip][pe], fmt='k.', ecolor='gray')

        divider = make_axes_locatable(self.ax2r)
        ax2rb = divider.append_axes("bottom", size=1.0, pad=0., sharex=self.ax2r)
        self.ax2r.errorbar(keblat.phase[keblat.clip][se], keblat.flux[keblat.clip][se] / lcpol[se],
                     keblat.dflux[keblat.clip][se], fmt='k.', ecolor='gray')
        self.ax2r.plot(keblat.phase[keblat.clip][se], lcmod[se], 'r.')
        ax2rb.errorbar(keblat.phase[keblat.clip][se], lcres[se], keblat.dflux[keblat.clip][se], fmt='k.', ecolor='gray')

        return

    def setup_axes(self):
        fig = plt.figure(figsize=(10,8))
        ax = plt.subplot2grid((2,2),(0,0),colspan=2)
        ax2a = plt.subplot2grid((2,2),(1,0))
        ax2b = plt.subplot2grid((2,2),(1,1))
        ax.errorbar(keblat.jd, keblat.flux, keblat.dflux, fmt='k.', ecolor='gray')
        ax.set_title("KIC {0} normalized SAP flux".format(self.kic))
        ax.set_xlim((keblat.jd[0], keblat.jd[-1]))
        print "Finished setting up plot"
        return fig, ax, ax2a, ax2b

    def on_key(self, event):
        if event.inaxes == self.ax:
            if event.key == 'x':
                print event.key, event.xdata
                self.x.append([event.xdata, event.xdata])
                #self.y.append([0,1])#list(self.ax.get_ylim()))
                #self.ax.plot([event.xdata, event.xdata], self.ax.get_ylim(), 'r-', lw=2, alpha=0.5, picker=5)
                self.ax.axvline(event.xdata, color='red', lw=2, alpha=0.5, picker=5)
            elif event.key == 'y':
                #self.x.append([0,1])
                self.y.append([event.ydata, event.ydata])
                self.ax.axhline(event.ydata, color='red', lw=2, alpha=0.5, picker=5)
            elif event.key == 'o':
                print "o"
            elif event.key == '?':
                print("Right click to select and overplot your input response. ? for help.")
            self.fig.canvas.draw()
        else:
            return

    def clear_xy(self):
        self.xx.append(self.x)
        self.yy.append(self.y)
        [self.ax.lines.pop() for j in range(len(self.x)+len(self.y))]
        self.fig.canvas.draw()
        self.x = []
        self.y = []


    def on_hover(self, event):
        #print "on hover ", event.x, event.y
        #print event.key
        return

    @staticmethod
    def get_obj_ind(base, targ):
        base_ids = np.array([id(base_i) for base_i in base])
        return np.where(base_ids == id(targ))[0]

    def on_pick(self, event):
        if event.inaxes == self.ax:
            #print("on pick", event.name, event.mouseevent, event.artist)
            if isinstance(event.artist, plt.Line2D):
                thisline = event.artist
                xdata = list(thisline.get_xdata())
                ydata = list(thisline.get_ydata())
                ind = event.ind
                # print('onpick1 line:', zip(np.take(xdata, ind), np.take(ydata, ind)))
                # 		if event.mouseevent.button == 3:
                obind = self.get_obj_ind(self.ax.lines, event.artist)
                #print('herf', ind, len(self.ax.lines), obind, id(event.artist))
                # print('onpick1 line:', ind, xdata, ydata, zip(np.take(xdata, ind), np.take(ydata, ind)))
                self.ax.lines[obind].remove()
                #print("onpick1 line:", list(xdata), list(ydata), self.x, self.y)
                if xdata in self.x: self.x.remove(xdata)
                if ydata in self.y: self.y.remove(ydata)
                #print('herf', len(self.ax.lines))
            # elif isinstance(event.artist, plt.Text):
            #     text = event.artist
            #     print('onpick1 text:', text.get_text())
            self.fig.canvas.draw()
        else:
            return


    def on_click(self, event):
        if event.inaxes == self.ax:
            print("on click", event.name, event.key)
            if event.button == 3:
                self.x.append([event.xdata, event.xdata])
                #self.y.append([0,1])#list(self.ax.get_ylim()))
                #self.ax.plot([event.xdata, event.xdata], self.ax.get_ylim(), 'r-', lw=2, alpha=0.5, picker=5)
                self.ax.axvline(event.xdata, color='red', lw=2, alpha=0.5, picker=5)
                # 		for ii in uin:
                # 			ax1.axvline(ii, color='red', lw=2, alpha=0.5, picker=5)
                print("event:", event, event.button, event.xdata, event.ydata)
            else:
                print("Nope")
            self.fig.canvas.draw()
        else:
            return
