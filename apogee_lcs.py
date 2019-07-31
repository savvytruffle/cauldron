import matplotlib
matplotlib.use('agg')
import sys, itertools, time, os
print(os.uname())
import emcee
from emcee.utils import MPIPool
#from helper_funcs import *
from eb_fitting import *

kic = int(sys.argv[1])
try:
    prefix = sys.argv[2]
except:
    prefix = '/astro/store/gradscratch/tmp/windemut/eb_fitting_mini/'

prefix = prefix + str(kic)+'/'
#kic = int(float(sys.argv[1]))
#period, tpe, esinw, ecosw, rsum, rrat, b, frat, q1, q2, q3, q4 = np.array(sys.argv[2:-2], dtype=float)
#nwalkers, niter = int(sys.argv[-2]), int(sys.argv[-1])

clobber_lc=False #overwrite LC only fits?
clobber_sed=True #overwrite SED only fits?
kiclist, perlist, pdeplist, sdeplist, morphlist = np.loadtxt('data/kebproperties_0216.dat',
                                          usecols=(0, 1, 3, 4, 8), unpack=True, delimiter=';')

#goodlist = (morphlist<0.6) & (pdeplist>0.1) & (sdeplist>0.01) & (perlist > 1.)
goodlist = (perlist>0)

excludelist = get_excludelist(fname='data/sed_flag_file_0328')

keblat = Keblat(preload=False)
#keblat.loadiso2()
keblat.loadiso2(isoname='isodata_jun4.dat')
keblat.loadsed(sedfile='data/kepsedall_0216.dat')
keblat.loadvkeb(filename='data/kebproperties_0216.dat')
goodlist_ind = np.where(kiclist[goodlist].astype(int) == kic)[0]
if len(goodlist_ind)>1:
    goodlist_ind=goodlist_ind[0]

goodv = keblat.kiclookup(kic, target=keblat.vkeb[:, 0])
keblat.loadlc(kic, keblat.vkeb[goodv, [1, 2, 5, 6, 7]], clip_tol=2.0, local_db=True)
keblat.morph = keblat.vkeb[goodv, 8]
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
if ecosw == 0: ecosw = 1e-5
if esinw == 0: esinw = 1e-5
frat = (keblat.vkeb[goodv, 4]/keblat.vkeb[goodv, 3])

ebv_arr, ebv_sig, ebv_dist_bounds, ebv_bounds = None, None, None, None#get3dmap(kic)

if keblat.swidth < 0.:
    print "No secondary eclipses detected. Exiting."
    sys.exit()

#if np.median(np.unique(keblat.crowd)) < 0.5:
#    print "Crowding > 0.5. Exiting."
#    sys.exit()

if (max(keblat.pwidth, keblat.swidth) > 0.04) and (keblat.pwidth+keblat.swidth > 0.091):
    clip_tol = 1.4
#elif ((keblat.pw3idth > 0.01) and (keblat.swidth > 0.01)) or (keblat.pwidth+keblat.swidth>:
#    clip_tol = 1.5
else:
    clip_tol = 1.7

print("Clip tolerance = {0}".format(clip_tol))
keblat.updatephase(tpe, period, clip_tol=clip_tol)


check_dir_exists(prefix)
keblat.start_errf(prefix+'lcfit.err')

# rvdata = np.loadtxt('data/{0}.rv'.format(kic), delimiter=';')

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
#opt_lcrvpars = opt_lcrv(keblat,msum=opt_rvpars[0], mrat=opt_rvpars[1],
#                         rsum=lcpars2[1], rrat=lcpars2[2], period=lcpars2[3],
#                         tpe=lcpars2[4], esinw=lcpars2[5], ecosw=lcpars2[6],
#                         b=lcpars2[7], frat=lcpars2[8], q1=lcpars2[-4],
#                         q2=lcpars2[-3], q3=lcpars2[-2], q4=lcpars2[-1],
#                         lcerr=0.0, k0=opt_rvpars[-2], rverr=0.)

q1, q2, q3, q4 = 0.01, 0.01, 0.01, 0.01
age, h0, dist = 9.2, 119., 850.
chunks = identify_gaps(keblat.cadnum, retbounds_inds=True)
chunks = np.delete(chunks, np.where(np.diff(chunks)<2)[0])
lcchi2_threshold = 5./np.nanmedian(np.array([np.nanmedian(abs(keblat.flux[chunks[ii]:chunks[ii+1]] -
                                                          np.nanmedian(keblat.flux[chunks[ii]:chunks[ii+1]])))
                                          for ii in range(len(chunks)-1)]))

keblat.pars['lcerr'] = 1e-5

print(blah)

################################### LC FITTING ########################################
if not os.path.isfile(prefix+'lcpars.lmfit2') or clobber_lc:

    # make initial guesses for rsum and f2/f1, assuming main sequence equal mass binary
    rsum = scipy.optimize.fmin_l_bfgs_b(estimate_rsum, 1.0,
                                        args=(period, 2*(keblat.pwidth+keblat.swidth)),
                                        bounds=[(1e-3, 1e3)], approx_grad=True)[0][0]

    # ew = scipy.optimize.fmin(tse_residuals, np.array([1e-3, ecosw]),
    #                          args=(period, tpe, tpe+keblat.sep*period))
    b = flatbottom(keblat.phase[keblat.clip], keblat.flux[keblat.clip], keblat.sep, keblat.swidth)
    rrat = guess_rrat(sdeplist[goodlist][goodlist_ind], pdeplist[goodlist][goodlist_ind])
    frat = rrat**(2.5)
    if rsum > 10:
        msum = 2.0
    else:
        msum = rsum
    
    ew_trials = [[esinw, ecosw], [-esinw, ecosw]]
    for jj in np.linspace(0.01, .5, 5):
        ew_trials = ew_trials + [[jj, ecosw], [-jj, ecosw]]
        #[[esinw, ecosw], [-esinw, ecosw], [-0.521, ecosw], [-0.332, ecosw], [-0.142, ecosw], [0.521, ecosw], [0.332, ecosw], [0.142, ecosw], [-.2]]
    lcpars0 = np.array([msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4])
    ew = ew_search_lmfit(ew_trials, keblat, lcpars0, (period, tpe, tpe+keblat.sep*period), fit_ecosw=False, polyorder=1)

    b_trials = [0.01, 0.1, 0.4, 0.8, 1.2]
    rrat_trials = [0.3, 0.7, 0.95]


    b_trials = [b] + [float(jj) for jj in np.array(b_trials)[~np.in1d(b_trials, b)]]
    rrat_trials = [rrat] + [float(jj) for jj in np.array(rrat_trials)[~np.in1d(rrat_trials, rrat)]]
    lc_search_counts=0
    bestlcchi2 = 1e25

    ###################################################################################
    ########################### LC ONLY OPTIMIZATION FIRST ############################
    ###################################################################################
    keblat.updatebounds('period', 'tpe', partol=0.1)
#    if pdeplist[goodlist][goodlist_ind]<0.1:
#        keblat.parbounds['rrat'] = [1e-6, 1.]
#        keblat.parbounds['frat'] = [1e-8, 1.]
    if sdeplist[goodlist][goodlist_ind] < 0.08 and pdeplist[goodlist][goodlist_ind] < 0.08:
        keblat.parbounds['rrat'] = [1e-6, 1]
        keblat.parbounds['frat'] = [1e-8, 1]

    if abs(ecosw) < 0.015:
        keblat.parbounds['ecosw'] = [-0.02, 0.02]
        keblat.parbounds['esinw'] = [-.9, .9]
    else:
        keblat.updatebounds('ecosw', partol=0.1)
        keblat.parbounds['esinw'] = [-np.sqrt(.9**2-ecosw**2), np.sqrt(.9**2-ecosw**2)]

    for i_b, i_rrat, i_ew in list(itertools.product(b_trials, rrat_trials, [ew, [esinw, ecosw]])):
#    for i_b, i_rrat, i_ew in list(itertools.product(b_trials, rrat_trials, [[esinw, ecosw]])):
        # lcpars0 = np.array([rsum, rsum, i_rrat, period, tpe, ew[0], ew[1], i_b, i_rrat**(2.5),
        #                     q1, q2, q3, q4])
        #upper_b = 2.*i_b if i_b==0.01 else 3.0
        #keblat.parbounds['b'][1] = upper_b
        opt_lcpars0 = opt_lc(keblat, msum=msum, rsum=rsum, rrat=i_rrat, period=period, tpe=tpe, esinw=i_ew[0],
                             ecosw=i_ew[1], b=i_b, frat=i_rrat**2.5, q1=q1, q2=q2, q3=q3, q4=q4)

        lcchi2 = np.sum(rez(opt_lcpars0, keblat, polyorder=2)**2)/(np.sum(keblat.clip) - len(opt_lcpars0) - 1)
        if (lcchi2 < bestlcchi2) or (lc_search_counts < 1):
            print "Saving from this run:", lcchi2, bestlcchi2, lc_search_counts
            bestlcchi2 = lcchi2*1.0
            opt_lcpars = opt_lcpars0.copy()
        lc_search_counts+=1

        if (bestlcchi2 <= 1.5) and opt_lcpars[2]<=1.0:
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
        keblat.plot_lc(opt_lcpars, prefix, polyorder=2, suffix='lc_opt2', savefig=True)
    except Exception, e:
        print str(e)


    if bestlcchi2 < lcchi2_threshold:
        print "Saving lmfit lcpars..."
        np.savetxt(prefix+'lcpars.lmfit2', opt_lcpars)
    else:
        print("Bestlcchi2 = {0}, exiting.".format(bestlcchi2))
        np.savetxt(prefix+'lcpars.lmfit2', opt_lcpars)
        sys.exit()
else:
    print "Loading lcpars lmfit"
    opt_lcpars = np.loadtxt(prefix+'lcpars.lmfit2')

bestlcchi2 = np.sum(rez(opt_lcpars[:13], keblat, polyorder=2)**2)/(np.sum(keblat.clip) - len(opt_lcpars) - 1)

msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4 = opt_lcpars[:13]
#opt_lcpars0 = opt_lc(keblat, msum=msum, rsum=rsum, rrat=rrat, period=period, 
#                     tpe=tpe, esinw=esinw, ecosw=ecosw, b=b, frat=frat, q1=q1, 
#                     q2=q2, q3=q3, q4=q4, vary_msum=False, fit_crowd=True)
opt_lcpars0 = opt_lc_crowd(keblat, msum=msum, rsum=rsum, rrat=rrat, period=period, 
                           tpe=tpe, esinw=esinw, ecosw=ecosw, b=b, frat=frat, 
                           q1=q1, q2=q2, q3=q3, q4=q4, vary_msum=False)
crowd_fits = keblat.broadcast_crowd(keblat.quarter, opt_lcpars0[13:])
keblat.crowd = crowd_fits
lcchi2 = np.sum(rez(opt_lcpars0[:13], keblat, polyorder=2)**2)/(np.sum(keblat.clip) - len(opt_lcpars0) - 1)
if (lcchi2 < bestlcchi2) and ((abs(opt_lcpars0[13:]-1)<1e-8).sum() < 1):
    print "Saving from this crowding fit run:", lcchi2, bestlcchi2
    bestlcchi2 = lcchi2*1.0
    opt_lcpars = opt_lcpars0[:13].copy()
    crowd = opt_lcpars0[13:]
    np.savetxt(prefix+'lcpars_wcrowd.lmfit2', opt_lcpars0)
else:
    crowd_fits = keblat.broadcast_crowd(keblat.quarter, keblat._crowdsap)
    keblat.crowd = crowd_fits
    crowd = keblat._crowdsap.ravel()

_, _ = keblat.lcfit(opt_lcpars[:13], keblat.jd, keblat.quarter, keblat.flux, keblat.dflux, keblat.crowd, polyorder=0)

keblat.updatephase(keblat.pars['tpe'], keblat.pars['period'], clip_tol=keblat.clip_tol)

opt_lcpars = np.append(opt_lcpars, np.log(np.median(abs(np.diff(keblat.flux)))))

####################### RV FITS ##########################
t, rv1, rv1err, rv2, rv2err = np.loadtxt('/astro/users/windemut/cauldron/data/{}.rv'.format(kic), usecols=(2,3,4,5,6), unpack=True)

esinw_array = np.sort(np.concatenate((np.linspace(0.1, 0.6, 4), -np.linspace(0.1, 0.6, 4), [0])))
for ii in range(len(esinw_array)):
    for jj in range(len(esinw_array)):
        if jj==0:
            lw=3
        else:
            lw=1
        _rvpars = rvpars.copy()
        _rvpars[4] = esinw_array[ii]
        _rvpars[5] = esinw_array[jj]
        rvmod = keblat.rvfit(_rvpars, rvt)
        print('C{}-'.format(ii), '{}'.format(1-jj/10.), '{}'.format(esinw_array[ii]))
        plt.plot(np.linspace(0,1,100), rvmod[0], 'C{}-'.format(ii), lw=lw, alpha=1-jj/10., label='{}'.format(esinw_array[ii]))

opt_rvpars = opt_rv(keblat, msum=m1+m2, mrat=m2/m1, period=opt_lcpars[3], tpe=opt_lcpars[4], esinw=opt_lcpars[5], ecosw=opt_lcpars[6], inc=keblat.pars['inc'], k0=k0, rverr=0)
m1, m2, k0 = keblat.rvprep(t, rv1*1e3, rv2*1e3, rv1err*1e3, rv2err*1e3)


def lnprob_lcrv(lcrvpars0, BF_constraint=None, retro=False):
    lcrvpars = lcrvpars0.copy()
    lp = keblat.lnprior_lcrv(lcrvpars)
    if np.isinf(lp):
        return -np.inf
    if BF_constraint is not None:
        lp += -0.5*((lcrvpars[9]-BF_constraint)/0.2)**2
    lcrvpars[-1] = np.exp(lcrvpars[-1])
    lcrvpars[-3] = np.exp(lcrvpars[-3])
    ll = keblat.lnlike_lcrv(lcrvpars, qua=np.unique(keblat.quarter), retro=retro)
    if (np.isnan(ll) or np.isinf(ll)):
        return -np.inf
    return lp + ll

nwalkers=128
niter=100000
clobber=False

BF_ratios = np.array([0.62, 0.77, 0.46, 0.39, 0.997, 0.65, 1./0.89])
BF_ratios_kics = np.array([5285607, 6864859, 6778289, 6449358, 4285087, 6131659, 6781535])
if ((1-keblat.pe_depth)<=0.1 and (1-keblat.se_depth)<=0.1) or (pdeplist[goodlist_ind]<=0.1 and sdeplist[goodlist_ind]<=0.1):
    BF_constraint=BF_ratios[BF_ratios_kics==kic][0]
else:
    BF_constraint=None

header = prefix+'lcrv_BF{}_'.format(0 if BF_constraint is None else 1)
footer = str(nwalkers)+'x'+str(niter/1000)+'k_final'
mcfile = header+footer+'.mcmc'

pars = opt_lcrvpars.copy()
#pars[-3] = np.log(opt_lcrvpars[-3]) 
#pars[-1] = np.log(np.nanmedian(keblat.rv1_err_obs))


ndim=len(pars)
p0_scale = np.ones(ndim)*1e-7
p0 = np.array([pars + p0_scale*pars*np.random.randn(ndim) for ii in range(nwalkers)])
for ii in range(ndim):
    p0[:,ii] = np.clip(p0[:,ii], keblat.parbounds[parnames_dict['lcrv'][ii]][0], keblat.parbounds[parnames_dict['lcrv'][ii]][1])
ll0 = np.zeros(nwalkers)
for ii in range(nwalkers):
    ll0[ii] = lnprob_lcrv(p0[ii,:], BF_constraint=BF_constraint, retro=retro)
if np.any(np.isinf(ll0)):
    print("Initial Gaussian ball of parameters yield -inf lnprob, check bounds")
    sys.exit()
if os.path.isfile(mcfile) and not clobber:
    print("File {0} already exists... do you want to clobber?".format(mcfile))
    sys.exit()
outf=open(mcfile, "w")
outf.close()

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_lcrv, threads=4, args=(BF_constraint,retro))
start_time = time.time()
print("Running {0}k MCMC chain".format(niter/1000))
for res in sampler.sample(p0, iterations=niter, storechain=False):
    if sampler.iterations % 10 == 0:
        position = res[0]
        outf = open(mcfile, "a")
        for k in range(position.shape[0]):
            outf.write("{0} {1} {2} {3} {4}\n".format(sampler.iterations,
                       k, sampler.acceptance_fraction[k], res[1][k],
                        " ".join([str(ii) for ii in position[k]])))
        outf.close()
    if sampler.iterations % 10000 == 0:
        print("Time elapsed since niter={0}:{1}".format(sampler.iterations, 
              time.time()-start_time))
print("Total time elapsed for MCMC run:{0}".format(time.time()-start_time))
print("Total acceptance fraction:{0}".format(np.mean(sampler.acceptance_fraction)))
try:
    print("Total autocorr time:{0}".format(np.mean(sampler.acor)))
except:
    print("Could not compute autocorr time...")


burnin=None
data=np.loadtxt(mcfile)
isonames = parnames_dict['lcrv']
blob_names=None
iwalker = np.arange(nwalkers)
afrac = np.empty((data.shape[0]/nwalkers, nwalkers))
logli = afrac*0.
params = np.empty((data.shape[0]/nwalkers, nwalkers, len(isonames)))
strays = []

for jj in iwalker:
    afrac[:,jj] = data[jj::nwalkers,2]
    logli[:,jj] = data[jj::nwalkers,3]
    if len(afrac[:,jj][(afrac[:,jj]<0.1)])>=0.66*len(afrac[:,jj]):
        strays.append(jj)
    for ii in range(len(isonames)):
        params[:, jj, ii] = data[jj::nwalkers, ii+4]

mostlike = np.where(logli == np.nanmax(logli))
mlpars = params[:,:,:][mostlike][0]
print "Max likelihood out of all samples: ", logli[:,:][mostlike]
for kk in range(len(isonames)):
    print("""{0} = {1}""".format(str(isonames[kk]), mlpars[kk]))
if burnin is None:
    print "Burn-in = when logli first crosses median value"
    burnin = np.where(np.nanmedian(logli, axis=1) >= np.nanmean(logli))[0][0] * 10
    
_, bad_walkers, walker_percentiles = get_stray_walkers(params, nwalkers, ndim, burnin=burnin)
strays = iwalker[bad_walkers>1]
print("{} bad/stray walkers = {}".format(len(strays), strays))

keep = iwalker[~np.in1d(iwalker, strays)]
if len(strays)>=0.33*nwalkers:
    keep = iwalker


fig = plt.figure(figsize=(16, 16))
for ii in range(len(isonames)):
    ax = fig.add_subplot(int(len(isonames)/2)+1, 2, ii+1)
    ax.plot(params[:, :, ii])
    ax.plot(np.nanmean(params[:,:,ii].T, axis=0), 'k-', lw=2, alpha=0.2)
    ax.plot([burnin, burnin], plt.ylim(), 'y-', lw=2.0)
    ax.set_xlabel('N/10 iteration')
    ax.set_ylabel(isonames[ii])
    divider = make_axes_locatable(ax)
    axhist = divider.append_axes("right", size=1.2, pad=0.1, sharey=ax)
    axhist.hist(params[:,:,ii], 100, histtype='step', alpha=0.6, normed=True,
                orientation='horizontal')
    axhist.hist(params[:,:,ii].ravel(), 100, histtype='step', color='k',
                normed=True, orientation='horizontal')
    plt.setp(axhist.get_yticklabels(), visible=False)
ax = fig.add_subplot(int(len(isonames)/2)+1, 2, len(isonames)+1)
ax.plot(logli)
ax.axvline(burnin, color='y', lw=2.0)
ax.set_xlabel('N/10 iteration')
ax.set_ylabel('logL')
plt.suptitle("KIC {} Parameter Trace".format(keblat.kic))
plt.savefig(header+footer+'_parameters.png')


bfpars = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
             zip(*np.nanpercentile(params[burnin:,keep,:].reshape((-1, ndim)),
                                [16, 50, 84], axis=0)))
try:
    blobpars = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
               zip(*np.nanpercentile(blobs[burnin:,keep,:].reshape((-1, len(blob_names))),
                                  [16, 50, 84], axis=0)))
except:
    print("No blobs.")
    
print("MCMC result: ")
print("Accep. Frac = ", np.mean(afrac[burnin:, keep]))
bffile = open(header+footer+'mcmc.pars', "w")
bffile.write("""{}\n""".format(" ".join([str(mp) for mp in mlpars])))
for kk in range(len(isonames)):
    print("""{0} = {1[0]} +{1[1]} -{1[2]}""".format(str(isonames[kk]),
          bfpars[kk]))
    bffile.write("""#{0} = {1[0]} +{1[1]} -{1[2]}\n""".format(str(isonames[kk]),
          bfpars[kk]))

changeofvar_names = ['m1', 'm2', 'r1', 'r2', 'inc', 'e']
params_changeofvar = np.zeros((params.shape[0], params.shape[1], len(changeofvar_names)))
params_changeofvar[:,:,0], params_changeofvar[:,:,1] = keblat.sumrat_to_12(params[:,:,0], params[:,:,1])
params_changeofvar[:,:,2], params_changeofvar[:,:,3] = keblat.sumrat_to_12(params[:,:,2], params[:,:,3])
params_changeofvar[:,:,4] = keblat.get_inc(params[:,:,8], params_changeofvar[:,:,2], keblat.get_a(params[:,:,4], params[:,:,0]))
params_changeofvar[:,:,5] = np.sqrt(params[:,:,6]**2 + params[:,:,7]**2)

bfpars_changeofvar = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
             zip(*np.percentile(params_changeofvar.reshape((-1, len(changeofvar_names))),
                                [16, 50, 84], axis=0)))
for kk in range(len(changeofvar_names)):
    print("""{0} = {1[0]} +{1[1]} -{1[2]}""".format(str(changeofvar_names[kk]),
          bfpars_changeofvar[kk]))
    bffile.write("""#{0} = {1[0]} +{1[1]} -{1[2]}\n""".format(str(changeofvar_names[kk]),
          bfpars_changeofvar[kk]))
bffile.close()


mlpars[[-1, -3]] = np.exp(mlpars[[-1, -3]])
keblat.plot_lcrv(mlpars, header+footer, '', savefig=True, thresh=12., retro=retro)
plt.savefig(header+footer+'.eps')
np.savetxt(header+footer+'.crowd', crowd)


import corner

plt.figure(figsize=(16, 16))
thin_by = np.clip((params.shape[0]-burnin)*params.shape[1]/50000, 1, 50000)
print("burned-in param matrix is {}; thinning by {}".format(params[burnin:, :, :].shape, thin_by))
post_inds = np.arange(len(isonames))
post_inds = np.delete(post_inds, np.where(np.nanstd(params[burnin::thin_by, :, :], axis=(0,1)) < 1e-12)[0])
try:
    corner.corner(params[burnin::thin_by, keep, :][:,:,post_inds].reshape((-1, len(post_inds))), 
                  labels=np.array(isonames)[post_inds], quantiles=[0.16, 0.5, 0.84],
                  show_titles=True, title_kwargs={"fontsize": 11})
    plt.suptitle("KIC {}".format(keblat.kic))
    plt.savefig(header+footer+'_posteriors.png')
except Exception, e:
    print(str(e), post_inds)

print "Making ACF plots now."
c=5
tau = np.zeros((ndim))
fig = plt.figure(figsize=(16, 16))
fig.subplots_adjust(hspace=0.0)
x = np.arange(params.shape[0])
for ii in range(ndim):
    ax = fig.add_subplot(ndim/3+1, 3, ii+1)            
    print("{}".format(isonames[ii]))
    tau[ii], mean_of_acfs, acf = get_acf_tau(params[:, keep, ii], c=c)
    ax.plot(acf, alpha=0.5)
    ax.plot(mean_of_acfs, 'k-', lw=1.5, alpha=0.8, label='mean(acfs)')
    ax.text(tau[ii]*1.1, 0.8, '{}'.format(int(tau[ii])))
    ax.plot(x, np.exp(-x/tau[ii]), '-', lw=1.5, alpha=0.8)
    ax.set_xlabel('N/10 iteration lag')
    ax.set_ylabel('{}'.format(isonames[ii]))
plt.suptitle("KIC {} ACF".format(keblat.kic))
plt.savefig(header+footer+'_ACF.png')


#################### CODE for KIC 6864859 heartbeat star #####################
mod,poly = keblat.jd*0., keblat.jd*0.
mod[keblat.clip],poly[keblat.clip] = keblat.lcfit(lcpars[:13], keblat.jd[keblat.clip], keblat.quarter[keblat.clip], 
                            keblat.flux[keblat.clip], keblat.dflux[keblat.clip], keblat.crowd[keblat.clip], 
                            polyorder=2)


for ii in range(9):
    plt.subplot(3, 3, ii+1)
    clip = (keblat.jd >= (keblat.tpe+keblat.period*ii - 2*keblat.pwidth*keblat.period)) * \
    (keblat.jd <= (keblat.tse +keblat.period*ii+ 2.*keblat.swidth*keblat.period)) * (keblat.quality<2)
#    mod,poly = keblat.lcfit(lcpars[:13], keblat.jd[clip], keblat.quarter[clip], 
#                            keblat.flux[clip], keblat.dflux[clip], keblat.crowd[clip], 
#                            polyorder=1, ooe=False)
    plt.plot(keblat.jd[clip], keblat.flux[clip], '.')
    plt.plot(keblat.jd[clip*keblat.clip], (mod*poly)[keblat.clip * clip], '-')
    #plt.ylim((0.999, 1.001))

_Ncad = int(np.ceil((keblat.tse + 2*keblat.swidth*keblat.period -
                     keblat.tpe - 2*keblat.pwidth*keblat.period) / (0.0204)))
_Neclipses = int(np.ceil((keblat.jd[-1]-keblat.jd[0])/period))+1
_tgrid = np.arange(keblat.tpe - 2*keblat.pwidth*keblat.period, 
                   keblat.tse + 2.1*keblat.swidth*keblat.period, 0.0204)
for ii in range(_Neclipses)[1:]:
    _tgrid = np.vstack((_tgrid,
                        np.arange(keblat.tpe+keblat.period*ii - 2*keblat.pwidth*keblat.period, 
                                  keblat.tse+keblat.period*ii + 2.1*keblat.swidth*keblat.period, 0.0204)))

_fgrid = _tgrid.copy()*np.nan
for ii in range(_Neclipses):
#    clip = (keblat.jd >= (keblat.tpe+keblat.period*ii - 1.2*keblat.pwidth*keblat.period)) * (keblat.jd <= (keblat.tse +keblat.period*ii + 1.2*keblat.swidth*keblat.period))

    clip = ((abs(keblat.jd-(keblat.tpe+keblat.period*ii))<2*keblat.pwidth*keblat.period) | (abs(keblat.jd-(keblat.tse +keblat.period*ii))<2*keblat.swidth*keblat.period)) * (keblat.quality<2)
    mod,poly = keblat.lcfit(lcpars[:13], keblat.jd[clip], keblat.quarter[clip], 
                            keblat.flux[clip], keblat.dflux[clip], keblat.crowd[clip], 
                            polyorder=1)
    linfit = np.poly1d(np.polyfit(keblat.jd[clip][mod>0.999], keblat.flux[clip][mod>0.999], 1))
#    _fgrid[ii,:] = np.interp(_tgrid[ii,:], keblat.jd[~clip * (keblat.quality<2)], keblat.flux[~clip * (keblat.quality<2)] / linfit(keblat.jd[~clip * (keblat.quality<2)]))
    
    _fgrid[ii,:] = np.interp(_tgrid[ii,:], keblat.jd[(keblat.quality<2)], keblat.flux[(keblat.quality<2)] / linfit(keblat.jd[(keblat.quality<2)]))

gap = np.where(np.diff(keblat.jd[(keblat.quality<2)])>7.*0.0204)[0]
bad = np.zeros(_tgrid.shape[0]*_tgrid.shape[1]).astype(bool)
for ii in range(len(gap)-1):
    bad = bad | ((_tgrid.ravel()>=keblat.jd[keblat.quality<2][gap[ii]]) * (_tgrid.ravel()<=keblat.jd[keblat.quality<2][gap[ii]+1]))

bad = bad.reshape((_tgrid.shape[0], _tgrid.shape[1]))

plt.figure(figsize=(8,5))
for ii in np.arange(_Neclipses)[~bad_eclipses]:
    good = ~(bad[ii,:] | np.append((np.diff(_fgrid[ii,:])==0), True))
    phase = ((_tgrid[ii,:][good])-keblat.tpe)%keblat.period/keblat.period
    phase[phase>0.8]-=1.
    plt.plot(phase, _fgrid[ii,:][good]+ii*0.0001, '.')
phase_peri = ((keblat.tpe - keblat.sudarsky(np.pi/2. - omega, e, period)) - keblat.tpe) % keblat.period/keblat.period
plt.axvline(phase_peri, linestyle='-', lw=2, alpha=0.7, color='0.1')

plt.text(0.001, 1.0035, '$\mathrm{PE}$', va='center', ha='center', rotation=35, fontsize=11)
plt.text(0.125, 1.0035, '$\mathrm{SE}$', va='center', ha='center', rotation=35, fontsize=11)
plt.text(phase_peri+0.001, 0.9997, '$\mathrm{periastron}$', fontsize=11)
plt.xlim((min(phase), max(phase)))
plt.ylim((0.9994, 1.004))
plt.xlabel('$\mathrm{Phase\ (P=40.8778\ d)}$')
plt.ylabel('$\mathrm{Offset Flux}$')
plt.title('$\mathrm{Eccentric\ Eclipsing\ System\ KIC\ 6864859}$')


#################### little snippet of code for cauldron on kic 6449 ####################
opt_lcrvpars = np.loadtxt(prefix_in+'lcrv_mass1.0_fix.lmfit')
keblat.plot_lcrv(opt_lcrvpars, prefix, savefig=False)
lcpars = [keblat.pars[zz] for zz in parnames_dict['lc']]
rvpars_case1 = [keblat.pars[zz] for zz in parnames_dict['rv']]
rv1_case1, rv2_case1 = keblat.rvfit(rvpars_case1, keblat.rv_t)

opt_lcrvpars = np.loadtxt(prefix_in+'lcrv_mass4.9_fix.lmfit')
keblat.plot_lcrv(opt_lcrvpars, prefix, savefig=False)
lcpars2 = [keblat.pars[zz] for zz in parnames_dict['lc']]
rvpars_case2 = [keblat.pars[zz] for zz in parnames_dict['rv']]
rv1_case2, rv2_case2 = keblat.rvfit(rvpars_case2, keblat.rv_t)


opt_lcrvpars = np.loadtxt(prefix_in+'lcrv_mass2.5_fix.lmfit')
keblat.plot_lcrv(opt_lcrvpars, prefix, savefig=False)
lcpars3 = [keblat.pars[zz] for zz in parnames_dict['lc']]
rvpars_case3 = [keblat.pars[zz] for zz in parnames_dict['rv']]
rv1_case3, rv2_case3 = keblat.rvfit(rvpars_case3, keblat.rv_t)

timestamps = np.array([2456557.73275,  2456559.72268, 2456584.63158, 2456585.63008, 
                       2456757.89224, 2456760.90501, 2456761.87212, 2456763.88043, 
                       2456784.82126, 2456786.79775, 2456787.80865, 2456815.78483, 
                       2456816.76558, 2456818.76389, 2456819.76152])
timestamps -= 2454833.

## load bfouts per visit ##
bfout = []
for ii in range(15):
    bfout.append(np.loadtxt('data/6449358BFOut_{}.txt'.format(ii+1)))
bfout = np.array(bfout)

bcv = np.loadtxt('data/6449358_bcv.txt')
fig=plt.figure()
for ii in range(15):                                                                                                                                                                 
    ax = fig.add_subplot(4,4,ii+1)
    ax.plot(bfout[ii, :, 0], bfout[ii, :, 1])
    ax.text(-244, -0.016, r'$\phi$='+str((timestamps[ii]-keblat.tpe)%keblat.period/keblat.period)[:5])
    ax.set_xlim((-250, 250))
    ax.axvline(rv2_case1[ii]*1e-3-bcv[ii,1], color='green', ls='--', 
                label=str(np.round(rvpars_case1[0], 1))+', '+str(np.round(rvpars_case1[1], 2)))
    ax.axvline(rv2_case3[ii]*1e-3-bcv[ii,1], color='orange', ls='--', 
                label=str(np.round(rvpars_case3[0], 1)) +', '+str(np.round(rvpars_case3[1], 2)))
    ax.axvline(rv2_case2[ii]*1e-3-bcv[ii,1], color='red', ls='--', 
                label=str(np.round(rvpars_case2[0], 1)) +', '+str(np.round(rvpars_case2[1], 2)))
fig.text(0.5, 0.04, 'Uncorrected Radial Velocity (km s$^{-1}$)', ha='center')
fig.text(0.04, 0.5, 'Broadening Function', va='center', rotation='vertical')
plt.legend(loc='upper left')

