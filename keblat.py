# -*- coding: utf-8 -*-

from glob import glob
import numpy as np
#from scipy import interpolate 
#from ext_func.occultquad import occultquad
#from ext_func.rsky import rsky
import matplotlib.pyplot as plt
import os.path
import time
from collections import OrderedDict
try:
    from helper_funcs import poly_lc_cwrapper, rsky, occultquad
except:
    print('Exception in from helper_funcs import poly_lc_cwrapper, rsky, occultquad')
    import datetime, os
    print('Timestamp: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
    print('Current dir: {0}'.format(os.getcwd()))
    print('helper_funcs.py exists? {0}'.format(os.path.isfile('helper_funcs.py')))
    print('helpers_linux.so exists? {0}'.format(os.path.isfile('helpers_linux.so')))
    # print(os.uname())
    try:
        from helper_funcs import poly_lc_cwrapper, rsky, occultquad
        print("worked here")
        # print(os.uname())
    except:
        time.sleep(30)
        # print(os.uname())
        print("just slept for 30")
        from helper_funcs import poly_lc_cwrapper, rsky, occultquad
        print("worked after sleep")

#from nbody import occultquad, rsky
from scipy.optimize import minimize as sp_minimize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import emcee

#############################################
################ constants ##################
#############################################
TWOPI = 2.0*np.pi
d2y = 365.242 #days in a year
r2au = 0.0046491 #solar radii to AU
parnames_dict = {'lc': ['msum', 'rsum', 'rrat', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'frat', 'q1',
                        'q2', 'q3', 'q4'],
                 'sed': ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'isoerr'],
                 'rv': ['msum', 'mrat', 'period', 'tpe', 'esinw', 'ecosw', 'inc', 'k0', 'rverr'],
                 'lcsed': ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'period', 'tpe','esinw',
                           'ecosw', 'b', 'q1', 'q2', 'q3', 'q4', 'lcerr', 'isoerr'],
                 'lcrv': ['msum', 'mrat', 'rsum', 'rrat', 'period', 'tpe', 'esinw', 'ecosw', 'b', 'frat',
                          'q1', 'q2', 'q3', 'q4', 'lcerr', 'k0', 'rverr']}
#############################################

class Keblat(object):
    """An EB fitting routine
    
    """
    def __init__(self, vkeb=False, preload=True):
        self.iso = None
        self.sed = None
        self.kic = None
        self.jd = None
        self.flux = None
        self.dflux = None
        self.vkeb = None
        self.fluxerr = None
        self.cadence = None
        self.rv1_obs = None
        self.rv2_obs = None
        self.rv1_err_obs = None
        self.rv2_err_obs = None
        self.rv_t = None
        self.exp = None
        self.clip_tol = 1.5
        self.res = []
        self.ldtype = 1 # quadratic limb darkening
        self.parnames = ['m1', 'm2', 'z0', 'age', 'dist', 'ebv', 'h0', 'period', 'tpe',
                         'esinw', 'ecosw', 'b', 'q1', 'q2', 'q3', 'q4', 'lcerr', 'isoerr',
                         'k0', 'rverr', 'msum', 'mrat', 'rsum', 'rrat', 'r1', 'r2', 'inc', 
                         'frat', 'e']

        self.pars = OrderedDict(zip(self.parnames, [None]*len(self.parnames)))

        self.parbounds = OrderedDict([('m1', [.1, 12.]), ('m2', [.1, 12.]), 
                          ('z0', [0.001, 0.06]), ('age', [6., 10.1]),
                          ('dist', [10., 15000.]), ('ebv', [0., 1.]),
                          ('h0', [119-1., 119+1.]), ('period', [0.05, 1000.]),
                          ('tpe', [0., 1e8]), ('esinw', [-.99, .99]),
                          ('ecosw', [-.99, .99]), ('b', [-7., 7.]), ('q1', [0., 1.]),
                          ('q2', [0., 1.]), ('q3', [0., 1.]), ('q4', [0., 1.]), 
                          ('lcerr', [0., 0.01]), ('isoerr', [0., 0.3]),
                                      ('k0', [-1e8, 1e8]), ('rverr', [0., 1e4]),
                                      ('msum', [0.01, 30.]), ('mrat', [0.0085, 10.]),
                                      ('rsum', [0.1, 1e6]), ('rrat', [1e-6, 1e3]),
                                      ('r1', [0.01, 1e6]), ('r2', [0.01, 1e6]),
                                      ('inc', [0., np.pi/2.]), ('frat', [1e-8, 2e2]),
                                    ('e', [0., .999])])
        if preload:
            self.loadiso2()
            self.loadsed()
            self.loadvkeb()
        if vkeb:
            self.loadvkeb()
        # solar values in cgs units
        self.tsun = 5778.
        self.lsun = 3.846e33
        self.rsun = 6.9598e10
        self.zsun = 0.01524
        self.msun = 1.989e33
        self.gsun = 6.6726e-8 * self.msun / (self.rsun**2)
        self.message='Initialized Properly'
        self.armstrongT1 = None
        self.armstrongdT1 = None
        self.armstrongT2 = None
        self.armstrongdT2 = None

    def updatebounds(self, *args, **kwargs):
        """Forces boundaries of specified arg parameters s.t. they are constrainted to 2% of parameter value

        Parameters
        ----------
        args : str
                for possible parameter names, see keblat.pars

        Example
        -------
        keblat.updatebounds('period', 'tpe', 'esinw', 'ecosw')

        """
        partol = kwargs.pop('partol', 0.02)
        for i in args:
            self.parbounds[i] = [self.pars[i]-abs(self.pars[i])*partol, self.pars[i]+abs(self.pars[i])*partol]
        if self.rv1_obs is not None and self.rv2_obs is not None:
            self.parbounds['k0'] = [min(np.nanmin(self.rv1_obs), np.nanmin(self.rv2_obs)),
                                     max(np.nanmax(self.rv1_obs), np.nanmax(self.rv2_obs))]
        return

    def updatepars(self, **kwargs):
        """Update free parameters to present values"""
        assert type(kwargs) is dict, "Must be dict"
        for key, val in kwargs.iteritems():
            if key in self.pars:
                self.pars[key] = val
            else:
                print("""{} is not a free parameter, sorry.""".format(key))
        #print "Updated pars"
        return
        # for ii in range(len(kwargs.keys())):
        #     if kwargs.keys()[ii] in self.pars:
        #         self.pars[kwargs.keys()[ii]] = kwargs.values()[ii]
        #     else:
        #         print("""{} is not a free parameter, sorry.""".format(kwargs.keys()[ii]))
        # print "made changes to updatepars"
        # return
        # #self.updatephase(self.tpe, self.period, self.clip_tol)

    def getpars(self, partype='allpars'):
        """Returns the values of the parameters for lc only, iso only, and iso+lc fits

        Parameters
        ----------
        partype : str (optional)
            'lc', 'iso', or 'allpars' (default).

        Returns
        -------
        x : array
            of length 13 (lc), 6 (iso), 18 (all)

        """
        if partype == 'allpars':
            return np.asarray(self.pars.values())
        elif partype == 'lc':
            lcpars = np.array([self.pars['m1'] + self.pars['m2'], self.pars['rsum'], self.pars['rrat'],
                    self.pars['period'], self.pars['tpe'], self.pars['esinw'], self.pars['ecosw'],
                    self.pars['b'], self.pars['frat'], self.pars['q1'], self.pars['q2'], self.pars['q3'],
                    self.pars['q4'], self.pars['lcerr']])
            return lcpars
        elif partype == 'sed':
            isopars = np.array([self.pars['m1'], self.pars['m2'], self.pars['z0'], self.pars['age'],
                                self.pars['dist'], self.pars['ebv'], self.pars['h0'], self.pars['isoerr']])
            return isopars
        elif partype == 'rv':
            inc = self.get_inc(self.pars['b'], self.pars['r1'],
                               self.get_a(self.pars['period'], self.pars['msum']))
            rvpars = np.array([self.pars['msum'], self.pars['mrat'], self.pars['period'], self.pars['tpe'],
                               self.pars['esinw'], self.pars['ecosw'], inc, self.pars['k0'], self.pars['rverr']])
            return rvpars
        elif partype == 'lcsed':
            return np.array([self.pars['m1'], self.pars['m2'], self.pars['z0'], self.pars['age'],
                             self.pars['dist'], self.pars['ebv'], self.pars['h0'], self.pars['period'],
                             self.pars['tpe'], self.pars['esinw'], self.pars['ecosw'], self.pars['b'],
                             self.pars['q1'], self.pars['q2'], self.pars['q3'], self.pars['q4'],
                             self.pars['lcerr'], self.pars['isoerr']])
        else:
            print "Not recognized parameter type, try again"
            return False

    def load_armstrong_temps(self, kic):
        """ Loads in temperature estimates from Armstrong et al and assigns them to the class

        Parameters
        ----------
        kic: int
            The KIC # of target EB

        Returns
        -------
        armstrong: array (4,)
            T1, dT1, T2, dT2 of EB components

        """
        armstrong = np.loadtxt('data/armstrong_keb_cat.csv', delimiter=',', usecols=(0,1,2,3,4))
        indx = self.kiclookup(kic, armstrong[:,0])
        self.armstrongT1 = armstrong[indx,1]
        self.armstrongdT1 = armstrong[indx,2]
        self.armstrongT2 = armstrong[indx,3]
        self.armstrongdT2 = armstrong[indx,4]
        return armstrong[indx, 1:5]

    def kiclookup(self, kic, target=None):
        """Function which returns the index of target ndarray that matches input KIC #
        
        Parameters
        ----------
        kic : int
            KIC # of target EB
        
        target : ndarray
            target array to search for matching input KIC #
            if None (default), search for match in SED matrix
            
        Returns
        -------
        sedidx : int or None
            index of array which corresponds to input KIC #
        
        Examples
        --------
        >>> keblat.kiclookup(9837578)
        38
        """
        if target is None:
            if self.sed is None:
                self.loadsed()
            target = self.sed[:, 0]
        sedidx = np.where(target.astype(int) == kic)[0]
        if len(sedidx)==0:
            print "No data matched to your input KIC number, double check it!"
            return None
        return sedidx[0]
        
    def _loadiso_defunct(self, isodir = 'data/'):
        """Function which loads Padova isochrone data
        
        Parameters
        ----------
        isodir : str
            directory which stores isochrone file
        
        Returns
        -------
        True : boolean
            if isodata loading is successful, stored in keblat.iso
        
        Examples
        --------
        >>> keblat.loadiso()
        Isodata.cat already exists; loaded.
        True
        """
        
        self.isodict = {'z': 0, 'logt': 1, 'mini': 2, 'mact': 3, 'logl': 4, 
                        'logte': 5, 'logg': 6, 'mbol': 7, 'mkep': 8, 'gmag': 9, 
                        'rmag': 10, 'imag': 11, 'zmag': 12, 'md51': 13, 
                        'Jmag': 14, 'Hmag': 15, 'Kmag': 16, 'Umag': 17, 
                        'Bmag': 18, 'Vmag': 19, 'w1': 20, 'w2': 21}
        fmt = ['%.4f', '%.2f', '%.8f', '%.4f', '%.4f', '%.4f', '%.4f', '%.3f', 
               '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', 
               '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f']
        isodatafile = isodir + 'isodata.dat'
        if os.path.isfile(isodatafile):
            iso = np.loadtxt(isodatafile, delimiter=',')
            self.iso = iso 
            print "Isodata.cat already exists; loaded."

        else:            
            sdssfiles = isodir + 'pad*.sdss.dat'
            noaofiles = isodir + 'pad*.noao.dat'
            wisefiles = isodir + 'pad*.wise.dat'
            
            sdss = glob(sdssfiles)
            noao = glob(noaofiles)
            wise = glob(wisefiles)
            
            sdss = np.sort(sdss)
            nsdss = len(sdss)
            noao = np.sort(noao)
            wise = np.sort(wise)
            
            sdss1 = np.loadtxt(sdss[0])
            noao1 = np.loadtxt(noao[0])
            wise1 = np.loadtxt(wise[0])
            iso = np.hstack((sdss1[:, :-2], noao1[:, [9,10,12]], 
                             wise1[:, 18:-4]))
            # things will be in mbol, mkep, g, r, i, z, 
            # d51, J, H, Ks, U, B, V, w1, w2
            for ii in range(1, nsdss):
                sdss1 = np.loadtxt(sdss[ii])
                noao1 = np.loadtxt(noao[ii])
                wise1 = np.loadtxt(wise[ii])
                iso1 = np.hstack((sdss1[:, :-2], 
                                  noao1[:, [9,10,12]], wise1[:, 18:-4]))
                iso = np.concatenate((iso, iso1))
            import operator
            self.iso = np.array(sorted(iso, key=operator.itemgetter(0, 1, 3)))
            np.savetxt(isodatafile, self.iso, delimiter=',', fmt=fmt)
            print "Isodata.dat created."

        self.zvals = np.unique(self.iso[:, self.isodict['z']])
        self.tvals = np.unique(self.iso[:, self.isodict['logt']])
        self.maxz, self.minz = np.max(self.zvals), np.min(self.zvals)
        self.maxt, self.mint = np.max(self.tvals), np.min(self.tvals)  
        return True
        
    def loadiso2(self, isodir = 'data/', ipnames=None):
        """Function which loads Padova isochrone data
        
        Parameters
        ----------
        isodir : str
            directory which stores isochrone file
        
        Returns
        -------
        True : boolean
            if isodata loading is successful, stored in keblat.iso
        
        Examples
        --------
        >>> keblat.loadiso()
        Isodata.cat already exists; loaded.
        True
        """
        
        self.isodict = {'z': 0, 'logt': 1, 'mini': 2, 'mact': 3, 'logl': 4, 
                        'logte': 5, 'logg': 6, 'mbol': 7, 'gmag': 8, 
                        'rmag': 9, 'imag': 10, 'zmag': 11, 'Umag': 12,  
                        'Bmag': 13, 'Vmag': 14, 'Jmag': 15, 'Hmag': 16, 
                        'Kmag': 17, 'w1': 18, 'w2': 19, 'mkep':20}
#        fmt = ['%.3f', '%.2f', '%.4f', '%.4f', '%.4f', '%.3f', '%.3f', '%.3f', 
#               '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', 
#               '%.3f', '%.3f', '%.3f', '%.3f',  '%.3f']
        isodatafile = isodir + 'isodata_final.dat'
        if os.path.isfile(isodatafile):
            iso = np.loadtxt(isodatafile, delimiter=',')
            self.iso = iso 
            print "Isodata.cat already exists; loaded."
        else:
            print "Isodata.cat does not exist."
            return False
        self.zvals = np.unique(self.iso[:, self.isodict['z']])
        self.tvals = np.unique(self.iso[:, self.isodict['logt']])
        self.maxz, self.minz = np.max(self.zvals), np.min(self.zvals)
        self.maxt, self.mint = np.max(self.tvals), np.min(self.tvals)

        if ipnames is None:
            self.ipname = np.array(['mkep', 'logl', 'logte', 'logg'])
        else:
            self.ipname = ipnames
        self.ipinds = np.array([self.isodict[ii] for ii in self.ipname]).astype(int)

        return True
        
    def loadsed(self, sedfile='data/kepsedall_0216.dat'):
        """Loads SED file 
        
        Parameters
        ----------
        sedfile : str
            name of SED file
        
        Returns
        -------
        True : boolean
            if file loads successfully, info stored in keblat.sed array
        False : boolean
            if file doesn't exist
        
        References
        ----------
        SED information from various sources: 
        RA, Dec, Glon, Glat, Teff, logg, Fe/H from KIC
        UBV (mag, error) from The Kepler Field UBV Survey (Everett+ 2012)
        griz (mag, error) from SDSS DR7 or DR9
        JHK (mag, error) from 2MASS
        w1w2 (mag, error) from WISE
        E(B-V), std E(B-V) from Schlegel dust maps (IRSA IPAC DUST)
        
        Examples
        --------
        >>> keblat.loadsed()
        True
        """

        if os.path.isfile(sedfile):
            self.sed = np.loadtxt(sedfile, delimiter=';')
        else:
            print "SED/obs file does not exist. Try again."
            return False
#        self.redconvert = np.array([4.107, 3.641, 2.682, 3.303, 2.285, 1.698, 
#            1.263, 0.189, 0.146, 0.723, 0.460, 0.310])

# distance grid from PANSTARRS 3D dust map (in pc)
        self.ebv_dist = np.array([63.09573445, 79.43282347, 100., 125.89254118,
            158.48931925, 199.5262315, 251.18864315, 316.22776602,
            398.10717055, 501.18723363, 630.95734448, 794.32823472, 1000.,
           1258.92541179, 1584.89319246, 1995.26231497, 2511.88643151,
           3162.27766017, 3981.07170553, 5011.87233627, 6309.5734448,
           7943.28234724, 10000., 12589.25411794, 15848.93192461,
          19952.62314969, 25118.8643151, 31622.77660168, 39810.71705535,
          50118.72336273, 63095.73444802])
        return True
    
    def getmags(self, kic, w2=True):
        """Fetches SED & extinction data corresponding to input KIC #
        
        Parameters
        ----------
        kic : int
            KIC #

        Returns
        -------
        magobs : array
            observed magnitudes in SED data file
        emagsobs: array
            associated uncertainties in observed magnitudes
        extinction: tuple
            E(B-V) and associated error from Schlegel dust maps
        glat : float
            galactic latitude of KIC
            
        Examples
        --------
        >>> magsobs, emagsobs, extinction, glat = keblat.getmags(9837578)

        """
        sedidx = self.kiclookup(kic)
        if np.equal(sedidx,None):
            print "No data matched to your input KIC number, double check it!"
            return
        magsobs = self.sed[sedidx, 8:-2:2].ravel()
        emagsobs = self.sed[sedidx, 9:-2:2].ravel()
        extinction = (self.sed[sedidx, -2], self.sed[sedidx, -1])
        glat = self.sed[sedidx, 4]
        zkic = 10**self.sed[sedidx, 7] * self.zsun
        if np.isnan(zkic):
            zkic = self.zsun
        return magsobs, emagsobs, extinction, glat, zkic
        
    def isoprep(self, magsobs, emagsobs, extinction, glat, zkic, exclude=[]):
        """Preps isochrone fitting for a given set of SED observations 
        and extinction parameters and gets rid of bad data
        
        Parameters
        ----------
        magsobs : array
            observed magnitudes
        emagsobs : array
            uncertainties of magnitudes
        extinction : tuple
            E(B-V) and std E(B-V)
        glat : float
            galactic longitude of target
        w2 : boolean
            True if want to include wise 2 band data in interpolation       
            
        Returns
        -------
        True : boolean
            if prep successful
            
        Examples
        --------
        >>> magsobs, emagsobs, extinction, glat = keblat.getmags(9837578)
        >>> keblat.isoprep(magsobs, emagsobs, extinction, glat)
        True
        """
        self.glat = glat
        (self.ebv, self.debv) = extinction
        if np.isfinite(zkic):
            self.z0 = zkic
        else:
            self.z0 = self.zsun
        fullmagnames = np.array(['Umag','Bmag','Vmag','gmag','rmag','imag',
                                 'zmag','w1', 'w2', 'Jmag','Hmag','Kmag'])
        fullmaglams = np.array([3733.9, 4308.9, 5516.6, 4716.7, 6165.1, 7475.9,
                            8922.9, 33200., 45700., 12300., 16400., 21600.])
        redconvert = np.array([4.107, 3.641, 2.682, 3.303, 2.285, 1.698, 1.263,
                               0.189, 0.146, 0.723, 0.460, 0.310])

        # according to Everett Howell 2012, typical limits for U, B, V
        self.exclude = exclude
        if magsobs[0] > 18.7:
            self.exclude += ['Umag']
        if magsobs[1] > 19.3:
            self.exclude += ['Bmag']
        if magsobs[2] > 19.1:
            self.exclude += ['Vmag']
        nodata = (magsobs == -99.999)

        # if np.any(magsobs[3:7] > abs(magsobs[0])+0.2):
        #     exclude = exclude + list(np.array(['gmag', 'rmag', 'imag', 'zmag'])[(magsobs[3:7] > abs(magsobs[0])+0.2)])
        suspect = (abs(magsobs-np.median(magsobs[~nodata])) > 4.) * (magsobs != -99.999)
        UBVset = np.array(['Umag','Bmag','Vmag'])
        grizset = np.array(['gmag','rmag','imag', 'zmag'])
        if np.in1d(UBVset, fullmagnames[suspect]).sum() > 0:
            suspect = suspect | ((abs(magsobs-np.median(magsobs[~nodata]))>2.) * (np.in1d(fullmagnames, UBVset)))
        if np.in1d(grizset, fullmagnames[suspect]).sum() > 0:
            suspect = suspect | ((abs(magsobs-np.median(magsobs[~nodata]))>2.) * (np.in1d(fullmagnames, grizset)))
        if np.sum(suspect)>0:
            self.exclude = self.exclude + list(fullmagnames[suspect])
        print "Excluding ", self.exclude

        self.exclude = np.unique(self.exclude)
#        for ii in range(len(exclude)):
#            nodata = nodata | (fullmagnames == exclude[ii])
        # artificially inflate bad data points
        outlier_mask = np.array([(fullmagnames == ii) for ii in self.exclude])
        if len(outlier_mask)>0:
            emagsobs[np.sum(outlier_mask, axis=0).astype(bool)] = 4.0

        obsinds = np.arange(len(fullmagnames))[~nodata]
        self.ipname = np.concatenate((fullmagnames[~nodata], 
                                np.array(['mkep', 'logl', 'logte', 'logg'])))
        self.ipinds = np.array([self.isodict[ii] for ii in self.ipname]).astype(int)
        
        self.maglams = fullmaglams[~nodata]
        self.magsobs = magsobs[obsinds]
        self.emagsobs = emagsobs[obsinds]

        if len(self.emagsobs)>0:
            self.emagsobs[self.emagsobs == -99.999] = np.clip(np.max(self.emagsobs), 0.025, 0.25)
        self.a_lam = redconvert[obsinds]

        return True
    
    def isoterpol(self, m0, z0, age):
        """Bi-linear interpolation given a mass, age, and metallicity 
        using Padova isochrones
        
        Parameters
        ----------
        m0 : float
            mass of star in solar mass
        z0 : float
            metallicity (i.e., z=0.017 for the Sun)
        age : float
            age in log10 (i.e., 6 for 1 Myr)
        
        Returns
        -------
        radius : float
            interpolated radius of star in solar radii
        fp : float array
            interpolated magnitudes, temp, and kepler flux
            
        Examples
        --------
        >>> keblat.isoterpol(1., 0.017, 9.4)
        """
        #age = np.log10(age)

        if (z0<=self.maxz) & (z0>=self.minz) & (age<=self.maxt) & (age>=self.mint):
            zbounds = np.digitize([z0], self.zvals)
            if zbounds == len(self.zvals):
                zbounds = np.concatenate((zbounds-2, zbounds-1))
            elif zbounds == 0:
                zbounds = np.concatenate((zbounds, zbounds+1))
            else:
                zbounds = np.concatenate((zbounds-1, zbounds))
            tbounds = np.digitize([age], self.tvals)
            if tbounds == len(self.tvals):
                tbounds = np.concatenate((tbounds-2, tbounds-1))
            elif tbounds == 0:
                tbounds = np.concatenate((tbounds, tbounds+1))
            else:
                tbounds = np.concatenate((tbounds-1, tbounds))
            
        else:
#            print "Error: z or age out of bounds! Note ",self.minz, 
#                "<z<", self.maxz, "and ", self.mint, "<logt<", self.maxt
            return (-np.inf, -np.inf)       

#        intergrid = np.empty((2, 2), dtype=object)
        intergrid = np.empty((2, 2, len(self.ipinds)))

        side = np.array(['left', 'right'])
        for ii in range(2):       
            goodindsz = np.searchsorted(self.iso[:, self.isodict['z']], 
                                        self.zvals[zbounds], side=side[ii])
            for jj in range(2):
                goodindst = np.searchsorted(self.iso[goodindsz[0]:goodindsz[1], 
                    self.isodict['logt']], self.tvals[tbounds], side=side[jj])
                mvals = np.unique(self.iso[goodindsz[0]:goodindsz[1], 
                            self.isodict['mact']][goodindst[0]:goodindst[1]])
                minm, maxm = np.min(mvals), np.max(mvals)
                if (m0<=maxm) & (m0>=minm):
                    mbounds = np.digitize([m0], mvals)
                    if mbounds == len(mvals):
                        mbounds = np.concatenate((mbounds-2, mbounds-1))
                    elif mbounds == 0:
                        mbounds = np.concatenate((mbounds, mbounds+1))
                    else:
                        mbounds = np.concatenate((mbounds-1, mbounds))
                    
                    goodindsm0 = np.searchsorted(self.iso[goodindsz[0]:goodindsz[1], self.isodict['mact']][goodindst[0]:goodindst[1]], mvals[mbounds])
                    goodindsm1 = np.searchsorted(self.iso[goodindsz[0]:goodindsz[1], self.isodict['mact']][goodindst[0]:goodindst[1]], mvals[mbounds], side='right')
                    try:
#                        intergrid[ii, jj] = interpolate.interp1d(self.iso[goodindsz[0]:goodindsz[1], self.isodict['mact']][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]], self.iso[goodindsz[0]:goodindsz[1]][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]][:, self.ipinds], axis=0, bounds_error=False)
                        intergrid[ii, jj, :] = np.array([np.interp(m0, self.iso[goodindsz[0]:goodindsz[1], self.isodict['mact']][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]], zz) for zz in self.iso[goodindsz[0]:goodindsz[1]][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]][:, self.ipinds].T])
                        #print "GOOD:", m0, age, z0, mvals[mbounds], self.tvals[tbounds], self.zvals[zbounds]
                        #print self.iso[goodindsz[0]:goodindsz[1], self.isodict['mini']][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]]
                        #print self.iso[goodindsz[0]:goodindsz[1]][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]][:, self.ipinds]
                    except Exception as e: 
                        print(e)
                        print "BAD:", m0, age, z0, mvals[mbounds], self.tvals[tbounds], self.zvals[zbounds]
                        #print goodindsz, goodindst, goodindsm0, goodindsm1
                        #print self.iso[goodindsz[0]:goodindsz[1], self.isodict['mact']][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]]
                        #print self.iso[goodindsz[0]:goodindsz[1]][goodindst[0]:goodindst[1]][goodindsm0[0]:goodindsm1[1]][:, self.ipinds]
                else:
                    #print "Error: mass out of bounds! Note 0.1<M/Msun<12"
                    return (-np.inf, -np.inf)

#        fq11 = intergrid[0,0](m0) #low z, low t
#        fq21 = intergrid[1,0](m0) #high z, low t
#        fq12 = intergrid[0,1](m0) #low z, high t
#        fq22 = intergrid[1,1](m0) #high z, high t
            
        zdiff = np.diff(self.zvals[zbounds])
        tdiff = np.diff(self.tvals[tbounds])
        fr1 = ((self.zvals[zbounds][1] - z0) / zdiff * intergrid[0,0,:]) + \
                ((z0 - self.zvals[zbounds][0]) / zdiff * intergrid[1,0,:])
        fr2 = ((self.zvals[zbounds][1] - z0) / zdiff * intergrid[0,1,:]) + \
                ((z0 - self.zvals[zbounds][0]) / zdiff * intergrid[1,1,:])
        fp = ((self.tvals[tbounds][1] - age) / tdiff * fr1) + \
                ((age - self.tvals[tbounds][0]) / tdiff * fr2)
        ipdict = dict(zip(self.ipname, list(np.arange(len(self.ipname)))))
        r1 = np.sqrt( 10**fp[ipdict['logl']] ) * (self.tsun/(10**fp[ipdict['logte']]))**2
        r2 = np.sqrt( m0 * self.gsun / (10**fp[ipdict['logg']]) )
        if abs(r1/r2 - 1.) > 0.05:
            print "Greater than 5% difference in radii: ", r1, r2, r1/r2
            print m0, z0, age
            return (-np.inf, -np.inf)
        radius = (r1+r2)/2.
        #lum = 10**self.fp[ipdict['logl']]
        return radius, fp
    
    def isofit(self, isopars, marginalize_distance=False):
        """Returns extincted, interpolated magnitudes corresponding to observed 
        wavelengths
        
        Parameters
        ----------
        isopars : float array
            parameters = (m1, m2, z0, age, dist, ebv, h0, isoerr)
            masses in solar mass, age in log10, dist & scaleheight in pc
        
        Returns
        -------
        magsmod : float array
            dust-extincted model mags 
        
        Examples
        --------
        >>> keblat.isofit([1., 1., 0.017, 9.4, 800., 0.032, 119., 0.04])
        """
        
        m1, m2, z0, age, dist, ebv, h0, isoerr = isopars
        #age = np.log10(age)
        self.r1, fp1 = self.isoterpol(m1, z0, age)
        self.r2, fp2 = self.isoterpol(m2, z0, age)
        if np.isinf(self.r1) or np.isinf(self.r2):
            #print "^Bad."
            return -np.inf

        mags1 = fp1[:len(self.magsobs)]
        mags2 = fp2[:len(self.magsobs)]
        self.temp1 = fp1[-2]
        self.temp2 = fp2[-2]
        self.logg1 = fp1[-1]
        self.logg2 = fp2[-1]
        self.frat = 10**((fp2[-4]-fp1[-4])/(-2.5))

        self.updatepars(m1=m1, m2=m2, z0=z0, age=age, dist=dist,
                        ebv=ebv, h0=h0, isoerr=isoerr, msum=m1+m2, mrat=m2/m1,
                        frat=self.frat, r1=self.r1, r2=self.r2,
                        rsum=self.r1+self.r2, rrat=self.r2/self.r1)

        absmagsmod = mags1 - 2.5 * np.log10(1. + 10**((mags1-mags2)/2.5))

        if marginalize_distance or dist<10. or dist>15000.:
            def magsmod_fn(x):
                magsmod = absmagsmod + 5. * np.log10(x / 10.) + \
                       self.a_lam * ebv * (1. - np.exp(-x * np.sin(self.glat * np.pi/180.) / h0))
                return np.sum(((magsmod-self.magsobs)/self.emagsobs)**2)
            res = sp_minimize(magsmod_fn, dist, method='L-BFGS-B', bounds=((10., 15000.),))
            dist = res.x
            self.message='marginalized distance'
        magsmod = absmagsmod + 5. * np.log10(dist / 10.) + \
            self.a_lam * ebv * (1. - np.exp(-dist * np.sin(self.glat * \
            np.pi/180.) / h0))
        return magsmod

    def loadlc(self, kic, properties, user_specified = None, 
               pdc = False, lc = True, clip_tol = 1.5, raw=False):
        """Loads Kepler SAP from database
        
        Parameters
        ----------
        kic : int
            KIC #
        properties : list
            EB period, time of PE, PE width, SE width, separation btw. eclipses
            note: this information should come from loadvkeb() function
        user_specified : 5 x Ndata ndarray
            Default is None. If NOT none, user_specified must be an ndarray 
            with the following structure: 
                    JD, FLUX, DFLUX, QUARTER, CROWDSAP
            where each of the five specified arrays have length Ndata
        pdc : boolean
            if True, returns Kepler's PDC data instead of SAP. Default False
        lc : boolean
            if True, returns long cadence data (else short cadence data). Default True
        clip_tol : float
            specifies tolerance around each eclipse to fit the data and model. Default = 1.5, ie includes
            eclipse and 1/2 eclipse durations before and after eclipse.
        raw : boolean
            if True, returns raw counts, else median divided. Default False

        Returns
        -------
        True : boolean
            if jd, phase, flux, fluxerr, crowd, clip loads successfully
        
        Examples
        --------
        >>> keblat.loadvkeb()
        >>> kic = 9837578
        >>> goodv = np.where(keblat.vkeb[:, 0] == kic)[0]
        >>> keblat.loadlc(kic, keblat.vkeb[goodv, [1, 2, 5, 6, 7]])
        True
        """
        self.kic = kic
        (self.period, self.tpe, self.pwidth, self.swidth, self.sep) = properties

        if self.tpe > 50000.:
            print "Your time of primary eclipse is > 50,000. Subtracting 54833 (Kepler BJD offset) from input value."
            self.tpe -= 54833.
        
        if user_specified is not None:
            self.jd = user_specified[0, :]
            self.flux, self.dflux = user_specified[1, :], user_specified[2, :]
            self.quarter, self.crowd = user_specified[3, :], user_specified[4, :]
            try:
                self.quality = user_specified[5, :]
            except:
                self.quality = self.jd*0.
            self.cadnum = None
        else:
            try:
                from loadlc_db import loadlc_db
                self.jd, self.flux, self.dflux, self.cadnum, self.quarter, self.quality, self._crowdsap = loadlc_db(kic, usepdc = pdc, lc = lc, raw = raw)
                self.crowd = self.broadcast_crowd(self.quarter, self._crowdsap)
            except:
                import kplr
                star = kplr.API().star(kic).get_light_curves(short_cadence=not lc)
                self.jd, self.flux, self.dflux = np.array([], dtype='float64'), np.array([], dtype='float64'), np.array([], dtype='float64')
                self.cadnum, self.quarter = np.array([], dtype='float64'), np.array([], dtype='float64')
                self.quality, self.crowd = np.array([], dtype='float64'), np.array([], dtype='float64')
                for ii in range(len(star)):
                    _s = star[ii].open(clobber=False)
                    _Npts = _s[1].data.shape[0]
                    self.jd = np.append(self.jd, _s[1].data['TIME'])
                    self.flux = np.append(self.flux, _s[1].data['SAP_FLUX'])
                    self.dflux = np.append(self.dflux, _s[1].data['SAP_FLUX_ERR'])
                    self.cadnum = np.append(self.cadnum, _s[1].data['CADENCENO'])
                    self.quarter = np.append(self.quarter, np.zeros(_Npts, dtype=int) + _s[0].header['QUARTER'])
                    self.quality = np.append(self.quality, _s[1].data['SAP_QUALITY'])
                    self.crowd = np.append(self.crowd, np.zeros(_Npts) + _s[1].header['CROWDSAP'])
                    naninds = np.isnan(self.flux)
                    self.jd = self.jd[~naninds]
                    self.dflux = self.dflux[~naninds]
                    self.cadnum = self.cadnum[~naninds]
                    self.quarter = self.quarter[~naninds]
                    self.quality = self.quality[~naninds]
                    self.crowd = self.crowd[~naninds]
                    self.flux = self.flux[~naninds]
        a = np.nanmedian(self.flux)
        if a>1:
            self.flux = self.flux/a
            self.dflux = abs(self.dflux/a)
#         self.phase = ((self.jd-self.tpe) % self.period)/self.period
#         self.phase[self.phase<-np.clip(self.pwidth*3., 0., 0.2)]+=1.
#         self.phase[self.phase>np.clip(self.sep+self.swidth*3., self.sep, 1.0)]-=1.
#         self.clip_tol = clip_tol
#         self.clip = (abs(self.phase)<self.clip_tol*self.pwidth) | \
#                     (abs(self.phase-self.sep)<self.clip_tol*self.swidth)
        self.updatephase(self.tpe, self.period, clip_tol=clip_tol)

        self.cadence = 0.0006811
        self.exp = 1

        self.fluxerr_tol = np.nanmedian(abs(np.diff(self.flux)))
        self.fluxerr = self.dflux.copy()
        self.updateErrors()
        if lc:
            self.cadence = 0.0204305556
            self.exp = 30
        print "LC data for KIC {0} loaded.".format(kic)
        return True

    def updateErrors(self, qtol=8, etol=10.):
        self.dflux = self.fluxerr.copy()
        self.dflux[(self.quality>qtol)] = etol*self.fluxerr_tol

    def updatephase(self, tpe, period, clip_tol=1.5):
        self.tpe = tpe
        self.period = period
        self.phase = ((self.jd-self.tpe) % self.period)/self.period
        self.phase[self.phase<-np.clip(self.pwidth*2., 0., 0.2)]+=1.
        self.phase[self.phase>1.-np.clip(self.pwidth*2., 0., 0.2)] -= 1.
        self.phase[self.phase>np.clip(self.sep+self.swidth*2., self.sep, 1.0)]-=1.
        self.clip_tol = clip_tol
        self.clip = ((abs(self.phase)<self.clip_tol*self.pwidth) | \
                     (abs(self.phase-1.0)<self.clip_tol*self.pwidth) | \
                    (abs(self.phase-self.sep)<self.clip_tol*self.swidth))
        #self.clip = self.clip * (self.quality < 0)
        return True
        
    def loadvkeb(self, 
                 filename='data/kebproperties_0216.dat',  
                 user_specified=None):
        """Loads information from Villanova Kepler EB database into keblat.vkeb, 
        including period, time of PE, PE width, SE width, sep, ecosw, esinw
        
        Parameters
        ----------
        filename : str
            name of villanova KEB properties file
        
        Returns
        -------
        True : boolean
            if keblat.vkeb loads successfully, has structure of:
            kic#, period, bjd0, pdepth, sdepth, pwidth, swidth, sep, ecosw, esinw
            
        Examples
        --------
        >>> keblat.loadvkeb()
        >>> print keblat.vkeb        
        """
        if user_specified is not None:
            self.vkeb = user_specified
        else:
            self.vkeb = np.loadtxt(filename,delimiter=";", 
                                   usecols=(0,1,2,3,4,5,6,7))
        # ^ kic#, period, bjd0, pdepth, sdepth, pwidth, swidth, sep
        ecosw = (self.vkeb[:, -1]*2. - 1.) * np.pi/4.
        esinw = ((self.vkeb[:, -3] / self.vkeb[:, -2]) - 1.) / ((self.vkeb[:, -3] / self.vkeb[:, -2]) + 1.)

#        switched = (self.vkeb[:, -3]<self.vkeb[:, -2])
#        esinw[switched] = ((self.vkeb[:, -2][switched] / self.vkeb[:, -3][switched]) - 1.) / \
#                          ((self.vkeb[:, -2][switched] / self.vkeb[:, -3][switched]) + 1.)


        self.vkeb[:, 2][self.vkeb[:,2]>50000] -= 54833.0
        self.vkeb = np.hstack((self.vkeb, ecosw[:, np.newaxis], esinw[:, np.newaxis]))
        return True        

    def start_errf(self, erfname):
        """Starts the error file
        Parameters
        ----------
        erfname : str
            name of error file
        """
        self.erfname = erfname
        errf = open(self.erfname, "w")
        errf.close()
        return True

    # Computes a template eclipse light curve with Mandel & Agol (2002) algorithm
    def _lctemplate_slow(self, lcpars, period, omega, e, a, inc, bgr, ldcoeffs,
                   rrat, tc, t0, cadence, exp, pe=True):
        """Computes a template Mandel & Agol (2002) eclipse lightcurve with
        correction for long-cadence binning (Kipping 2013)

        Parameters
        ----------
        period : float
            period of EB
        omega : float
            argument of periastron
        e : float
            eccentricity
        a : float
            semi-major axis in AU
        inc : float
            inclination; 90 = edge on
        bgr : float
            radius of star being eclipsed in solar radius
        ldcoeffs: float array
            limb darkening coefficients; tuple if quadratic,
            4 elements if quartic non-linear law
        rrat : float
            ratio of radii (eclipsing/eclipsed)
        tc : float
            time of center of eclipse (either PE or SE)
        t0 : float
            time of periastron passage
        cadence : float
            cadence (in days)
        exp : int
            number of Kepler exposures to bin (LC = 30, SC = 1)
        pe : boolean
            if True, template for PE

        Returns
        -------
        tmean : float array
            array of times for PE (or SE)
        resmean : float array
            array of flux values for PE (or SE)

        """

        if pe:
            half0 = (self.pwidth*period/2.)*self.clip_tol
        else:
            half0 = (self.swidth*period/2.)*self.clip_tol

#        half0 = period/2.

        t = np.linspace(-half0, +half0, int(2*half0/0.0006811) + 1)
        t += tc
        maf = rsky(e, period, t0, 1e-8, t)
        r = a*(1.-e**2) / (1.+e*np.cos(maf))
        z = r*np.sqrt(1.-np.sin(omega+maf)**2*np.sin(inc)**2)
        finite_z = np.isfinite(z)
        if np.sum(~finite_z)>0:
            print maf, e, a, inc, z
            if np.sum(~finite_z) < 0.05*len(z):
                z[~finite_z] = np.interp(t[~finite_z], t[finite_z], z[finite_z])
            else:
                exit
        if self.ldtype == 0:
            # print "non-linear case:", ldcoeffs
            res = self.occultnltheta(z/(bgr*r2au), rrat, ldcoeffs)
        else:
            # print "Quad case:", ldcoeffs
            res = occultquad(z/(bgr*r2au), ldcoeffs[0], ldcoeffs[1], rrat, len(z))

        bad = (res<-1e-4) | (res-1.>1e-4)
        if np.sum(bad)>0:
            cz = z[bad]/(bgr*r2au)
            interp_res = np.interp(t[bad], t[~bad], res[~bad])
            errfile = open(self.erfname, "a")
            errfile.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(min(z)/(bgr*r2au), max(z)/(bgr*r2au), ldcoeffs[0], ldcoeffs[1], rrat, " ".join([str(ii) for ii in cz]), " ".join([str(ii) for ii in res[bad]]), " ".join([str(ii) for ii in interp_res])))
            errfile.close()
            res[bad] = interp_res
        if cadence < 0.02:
            return t, res

        tt = t[:, np.newaxis] * np.ones(exp)
        rres = res[:, np.newaxis] * np.ones(exp)
        tt[:, 0] = t
        rres[:, 0] = res
        for ii in range(1, exp):
            tt[:-ii, ii] = t[ii:]
            rres[:-ii, ii] = res[ii:]
        tmean = np.mean(tt[:-exp+1, :], axis=1)
        resmean = np.mean(rres[:-exp+1, :], axis=1)
        return tmean, resmean

    def lctemplate(self, lcpars, period, omega, e, a, inc, bgr, ldcoeffs, 
                   rrat, tc, t0, cadence, exp, pe=True):
        """Computes a template Mandel & Agol (2002) eclipse lightcurve with
        correction for long-cadence binning (Kipping 2013)

        Parameters
        ----------
        period : float
            period of EB
        omega : float
            argument of periastron
        e : float
            eccentricity
        a : float
            semi-major axis in AU
        inc : float
            inclination; 90 = edge on
        bgr : float
            radius of star being eclipsed in solar radius
        ldcoeffs: float array
            limb darkening coefficients; tuple if quadratic, 
            4 elements if quartic non-linear law
        rrat : float
            ratio of radii (eclipsing/eclipsed)
        tc : float
            time of center of eclipse (either PE or SE)
        t0 : float
            time of periastron passage
        cadence : float
            cadence (in days)
        exp : int
            number of Kepler exposures to bin (LC = 30, SC = 1)
        pe : boolean
            if True, template for PE
            
        Returns
        -------
        tmean : float array
            array of times for PE (or SE)
        resmean : float array
            array of flux values for PE (or SE)
        
        """
        
        if pe:
            half0 = self.pwidth*period * 2.
        else:
            half0 = self.swidth*period * 2.

        #half0 = period/4.

        t = np.linspace(-half0, +half0, int(2*half0/0.0006811) + 1)
        t += tc
        maf = rsky(e, period, t0, 1e-8, t)
        r = a*(1.-e**2) / (1.+e*np.cos(maf))
        z = r*np.sqrt(1.-np.sin(omega+maf)**2*np.sin(inc)**2)
        if not np.all(np.isfinite(z)):
            badz = np.isfinite(z)
            print maf[~badz], e, a, inc, z[~badz], r[~badz]
            z[~badz] = np.interp(t[~badz], t[badz], z[badz])
            errfile = open(self.erfname, "a")
            errfile.write("inf z: {0} {1} {2} {3} {4} {5} {6}\n".format(e, a, inc, t[~badz], maf[~badz], r[~badz], z[~badz]))
            errfile.close()
        if self.ldtype == 0:
            # print "non-linear case:", ldcoeffs
            res = self.occultnltheta(z/(bgr*r2au), rrat, ldcoeffs)
        else:
            # print "Quad case:", ldcoeffs
            res = occultquad(z/(bgr*r2au), ldcoeffs[0], ldcoeffs[1], rrat, len(z))

        bad = (res<-1e-4) | (res-1.>1e-4)        
        if np.sum(bad)>0:
            cz = z[bad]/(bgr*r2au)
            interp_res = np.interp(t[bad], t[~bad], res[~bad])
            errfile = open(self.erfname, "a")
            errfile.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(min(z)/(bgr*r2au), max(z)/(bgr*r2au), ldcoeffs[0], ldcoeffs[1], rrat, " ".join([str(ii) for ii in cz]), " ".join([str(ii) for ii in res[bad]]), " ".join([str(ii) for ii in interp_res])))
            errfile.close()
            res[bad] = interp_res
        if cadence < 0.02:
            return t, res
        as_strided = np.lib.stride_tricks.as_strided
        tt = as_strided(t, (len(t)+1-self.exp, self.exp), (t.strides * 2))
        rres = as_strided(res, (len(res)+1-self.exp, self.exp), (res.strides * 2))
        if pe:
            self.pe_dur = tt[(rres<1.)][-1] - tt[(rres<1.)][0] if np.sum(rres<1.)>1 else 0
            self.pe_depth = np.min(rres)
        else:
            self.se_dur = tt[(rres<1.)][-1] - tt[(rres<1.)][0] if np.sum(rres<1.)>1 else 0
            self.se_depth = np.min(rres)
        return np.mean(tt, axis=1), np.mean(rres, axis=1)
        
    # LIGHT CURVE MODEL
    def _lcfit_slow(self, lcpars, jd, phase, flux, dflux, crowd, 
              polyorder=2):
        """Computes light curve model
        
        Parameters
        ----------
        lcpars : float array
            parameters for LC fitting: 
            msum, rsum, rratio, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4
        jd : float array
            time array
        phase : float array
            corresponding phase
        flux : float array
            observed flux
        dflux : float array
            flux error
        crowd : float array
            array of crowding values (additional flux)
        polyorder : int
            order of polynomial to detrend lightcurve
        
        Returns
        -------
        totmod : float array
            array of model fluxes
        totpol : float array
            array of polynomials for detrending
        """
        # r1, r2, frat derive from m1, m2, z0, t0, dist, E(B-V), scaleheight
        msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, \
            q1, q2, q3, q4 = lcpars       
#        self.updatepars(m1=m1, m2=m2, period=period, tpe=tpe, esinw=esinw, 
#                    ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4)
        # LD transformations (Kipping 2013)
        c1 = 2.*np.sqrt(q1)*q2
        c2 = np.sqrt(q1)*(1.-2.*q2)
        c3 = 2.*np.sqrt(q3)*q4
        c4 = np.sqrt(q3)*(1.-2.*q4)
        ldcoeffs1 = np.array([c1, c2])
        ldcoeffs2 = np.array([c3, c4])
            
#        if r2 > r1:
#            r1, r2 = r2, r1
#            m1, m2 = m2, m1
#            frat = 1./frat
        omega=np.arctan2(esinw,ecosw)
        e=np.sqrt(esinw**2+ecosw**2)

        # nip it at the bud.
        if (e>=1.):
            #print "e>=1", e
            return -np.inf, -np.inf
            


        r1 = rsum/(1.+rrat)
        r2 = rsum/(1.+1./rrat)
        a = ((period/d2y)**2 * (msum))**(1./3.)
        inc = np.arccos(b*r1/(a/r2au))
        
        if np.isnan(inc):
            #print "inc is nan", inc
            return -np.inf, -np.inf
        
        fpe = np.pi/2. - omega
        fse = -np.pi/2. - omega

        # transform time of center of PE to time of periastron (t0)
        # from Eq 9 of Sudarsky et al (2005)
        t0 = tpe - (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
            (e*np.sin(fpe)/(1.+e*np.cos(fpe)) - 2.*(1.-e**2)**(-0.5) * \
            np.arctan(np.sqrt(1.-e**2) * np.tan((fpe)/2.) / (1.+e)))
        tse = t0 + (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
            (e*np.sin(fse)/(1.+e*np.cos(fse)) - 2.*(1.-e**2)**(-0.5) * \
            np.arctan(np.sqrt(1.-e**2) * np.tan((fse)/2.) / (1.+e)))

        # if tse<tpe:
        #     tse+=period
            
        tempt1, tempres1 = self.lctemplate(lcpars, period, omega, e, a, inc, r1,
                                           ldcoeffs1, r2/r1, tpe, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=True)
        tempt2, tempres2 = self.lctemplate(lcpars, period, omega, e, a, inc, r2,
                                           ldcoeffs2, r1/r2, tse, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=False)

        tempt1 = tempt1 % period
        tempt2 = tempt2 % period
        tempres1 = (tempres1 - 1.)/(1. + frat) + 1.
        tempres2 = (tempres2 - 1.)/(1. + 1./frat) + 1.

        sorting1 = np.argsort(tempt1)
        sorting2 = np.argsort(tempt2)

        tempres1 = tempres1[sorting1]
        tempt1 = tempt1[sorting1]
        tempres2 = tempres2[sorting2]
        tempt2 = tempt2[sorting2]

        #not including crowdsap term.
        #tempres1 = (tempres1 + frat) / (1.+frat)
        #tempres2 = (tempres2 * frat + 1.) / (1. + frat)
        totpol, totmod = np.ones(len(jd)), np.ones(len(jd))

        if polyorder>0:
    
            # mask out continuum data
            #clip = ((abs(phase)<1.5*self.pwidth) | (abs(phase-self.sep) < 1.5*self.swidth))
            clip = (jd>0)
    
            chunk = np.array(np.where(np.diff(jd[clip]) > self.pwidth*period))[0]
            #put in dummy first and last element # placeholders
            
            chunk = np.append(chunk, len(jd[clip])-2).flatten()
            _, chunk3 = np.unique(np.searchsorted(jd[clip][chunk], jd), return_index=True)
            chunk=chunk3
            chunk[-1]+=1
    #        plt.plot(self.jd, self.flux, 'ro', self.jd[clip], self.flux[clip], 'go', self.jd[chunk[:-1]],self.flux[chunk[:-1]], 'mo')
    #        plt.show()
            for i in range(len(chunk)-1):
                #print i, chunk[i], chunk[i+1], self.jd[chunk[i]:chunk[i+1]]
                t = jd[chunk[i]:chunk[i+1]]
                f = flux[chunk[i]:chunk[i+1]]
                ef = dflux[chunk[i]:chunk[i+1]]
                crow = crowd[chunk[i]:chunk[i+1]]
                maf = rsky(e, period, t0, 1e-8, t)
                npts=len(maf)
                #use this version for full lightcurve treatment...
                r = a*(1.-e**2) / (1.+e*np.cos(maf))
                zcomp = np.sin(omega+maf) * np.sin(inc) 
                #z = r*np.sqrt(1.-zcomp**2)
                pe = ((r*zcomp>0.))# & (z <= 1.05*(r1+r2)*r2au))
                se = ((r*zcomp<0.))# & (z <= 1.05*(r1+r2)*r2au))
                model = np.ones(npts)
    #            sse = (((maf+omega) % (TWOPI))>np.pi)
    #            ppe = (((maf+omega) % (TWOPI))<=np.pi)
                
    #            plt.plot(t, f, 'ro', t[pe], f[pe], 'go', t[se], f[se], 'bo')
    #            plt.title(str(i))
    #            plt.show()
                if pe.any():
#                    shift = period * np.round((np.mean(t[pe]) - tpe)/period)
                    model[pe] = np.interp(t[pe]%period, tempt1, tempres1)
                    model[pe] = (model[pe] - 1.) * crow[pe] + 1.
    #                print "PE: mean(t[pe]), tpe, (mean(t[pe])-tpe)/period, round ver"
    #                print np.mean(t[pe]), tpe, (np.mean(t[pe]) - tpe)/period, np.round((np.mean(t) - tpe)/period), len(t[pe]), len(f[pe]), len(tempt1), len(tempres1), len(model[pe])
    #                plt.plot(t[pe]-shift, f[pe], 'ro', tempt1, (tempres1-1.)*crow[pe][0] + 1., 'bo', t[pe]-shift, model[pe], 'go')
    #                plt.title('pe')
    #                plt.show()
    #                plt.close('all')
                if se.any():
#                    shift = period * np.round((np.mean(t[se]) - tse)/period)
                    model[se] = np.interp(t[se]%period, tempt2, tempres2)
                    model[se] = (model[se] - 1.) * crow[se] + 1.
    
    #                print "SE"
    #                print np.mean(t[se]), tse, (np.mean(t[se]) - tse)/period, np.round((np.mean(t[se]) - tse)/period)
    #                plt.plot(t[se]-shift, f[se], 'ro', tempt2, (tempres2-1.)*crow[se][0] + 1., 'bo', t[se]-shift, model[se], 'go')
    #                plt.title('se')
    #                plt.show()
    #                plt.close('all')
    #            else:
    #                print "This data bundle does not belong to SE or PE"
                # marginalization (2nd order polynomial fit to residuals)

                bad = (model<1)
                tt = t[~bad]
                mmodel = model[~bad]
                ff = f[~bad]
                eef = ef[~bad]
                nnpts = len(ff)
                tnew = tt - np.mean(tt)
                #if len(t[~bad]) < 1:
                    #print "Npts ooe = ",len(t[~bad])
    
                # Bk = sum over i (D_i/M_i)(tdiff_i)^k / (sigma_i/M_i)^2
                # matrix 3 rows x npts columns since quadratic polynomial 
                # fit requires 3 coeffs
                if bad[0] or bad[-1]:
                    poly_remember = polyorder
                    polyorder=1
                #number of 'i' data or model points; polynomial order
                order_pow = np.arange(polyorder+1)
                t_pow = tnew[:,np.newaxis]**order_pow
                Bk = np.ones(shape=(polyorder+1,nnpts))*((ff/mmodel)/(eef/mmodel)**2)
                Bk*=t_pow.T
                #sum along 'i' (or along each row)
                Bksum = np.sum(Bk,axis=1)
                #Mjk = sum over i (tdiff_i)^j * (tdiff_i)^k / (sigma_i/M_i)^2
                #construct 3 rows x npts columns 
                Mj = np.ones(shape=(polyorder+1,nnpts))/(eef/mmodel)**2
                Mj*=t_pow.T
                #transform from 2D (j rows x i columns) to 3D (k x j x i)
                t_pow_3d = tnew[:,np.newaxis,np.newaxis]**order_pow
                Mjk = t_pow_3d.T * Mj[np.newaxis,:,:]
                #now sum along 'i' 
                Mjksum = np.sum(Mjk,axis=2)
                #do matrix inversion solver thing to get polynomial coeffs
                try:
                    Aj = np.linalg.lstsq(Mjksum,Bksum)[0]
                    pol = np.polyval(Aj[::-1],t-np.mean(t))
                except: 
                    pol = np.ones(npts)
                #Aj = np.dot(np.linalg.pinv(Mjksum), Bksum)
    #                plt.plot(t, f, 'ro', t, model*pol, 'go')
    #                plt.plot(t, pol, 'ms', tt, np.polyval(Aj[::-1],tnew), 'cs')
    #                plt.show()
                if bad[0] or bad[-1]:
                    polyorder = poly_remember
                totmod[chunk[i]:chunk[i+1]] = model
                totpol[chunk[i]:chunk[i+1]] = pol

        else:
            maf = rsky(e, period, t0, 1e-8, jd)
            r = a*(1.-e**2) / (1.+e*np.cos(maf))
            zcomp = np.sin(omega+maf) * np.sin(inc) 
            #z = r*np.sqrt(1.-zcomp**2)
            pe = ((r*zcomp>0.)) #& (z <= 1.05*(r1+r2)*r2au))
            se = ((r*zcomp<0.)) #& (z <= 1.05*(r1+r2)*r2au))
            tt = jd % period
            if pe.any():
                totmod[pe] = np.interp(tt[pe], tempt1, tempres1)
                totmod[pe] = (totmod[pe] - 1.) * crowd[pe] + 1.
            if se.any():
                totmod[se] = np.interp(tt[se], tempt2, tempres2)
                totmod[se] = (totmod[se] - 1.) * crowd[se] + 1.
        #     if np.sum(totmod[se]-1.) == 0.:
        #         return np.ones_like(totmod), totpol
        # if np.sum(totmod-1.) == 0.:
        #     return totmod, totmod
        return totmod, totpol

    def _lcfit_slow2(self, lcpars, jd, phase, flux, dflux, crowd, 
              polyorder=2):
        """Computes light curve model
        
        Parameters
        ----------
        lcpars : float array
            parameters for LC fitting: 
            msum, rsum, rratio, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4
        jd : float array
            time array
        phase : float array
            corresponding phase
        flux : float array
            observed flux
        dflux : float array
            flux error
        crowd : float array
            array of crowding values (additional flux)
        polyorder : int
            order of polynomial to detrend lightcurve
        
        Returns
        -------
        totmod : float array
            array of model fluxes
        totpol : float array
            array of polynomials for detrending
        """
        # r1, r2, frat derive from m1, m2, z0, t0, dist, E(B-V), scaleheight
        msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, \
            q1, q2, q3, q4 = lcpars       
#        self.updatepars(m1=m1, m2=m2, period=period, tpe=tpe, esinw=esinw, 
#                    ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4)
        # LD transformations (Kipping 2013)
        c1 = 2.*np.sqrt(q1)*q2
        c2 = np.sqrt(q1)*(1.-2.*q2)
        c3 = 2.*np.sqrt(q3)*q4
        c4 = np.sqrt(q3)*(1.-2.*q4)
        ldcoeffs1 = np.array([c1, c2])
        ldcoeffs2 = np.array([c3, c4])
            
#        if r2 > r1:
#            r1, r2 = r2, r1
#            m1, m2 = m2, m1
#            frat = 1./frat
        omega=np.arctan2(esinw,ecosw)
        e=np.sqrt(esinw**2+ecosw**2)

        # nip it at the bud.
        if (e>=1.):
            #print "e>=1", e
            return -np.inf, -np.inf
            


        r1 = rsum/(1.+rrat)
        r2 = rsum/(1.+1./rrat)
        a = ((period/d2y)**2 * (msum))**(1./3.)
        inc = np.arccos(b*r1/(a/r2au))
        
        if np.isnan(inc):
            #print "inc is nan", inc
            return -np.inf, -np.inf
        
        fpe = np.pi/2. - omega
        fse = -np.pi/2. - omega

        # transform time of center of PE to time of periastron (t0)
        # from Eq 9 of Sudarsky et al (2005)
        t0 = tpe - (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
            (e*np.sin(fpe)/(1.+e*np.cos(fpe)) - 2.*(1.-e**2)**(-0.5) * \
            np.arctan(np.sqrt(1.-e**2) * np.tan((fpe)/2.) / (1.+e)))
        tse = t0 + (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
            (e*np.sin(fse)/(1.+e*np.cos(fse)) - 2.*(1.-e**2)**(-0.5) * \
            np.arctan(np.sqrt(1.-e**2) * np.tan((fse)/2.) / (1.+e)))

        # if tse<tpe:
        #     tse+=period
            
        tempt1, tempres1 = self.lctemplate(lcpars, period, omega, e, a, inc, r1,
                                           ldcoeffs1, r2/r1, tpe, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=True)
        tempt2, tempres2 = self.lctemplate(lcpars, period, omega, e, a, inc, r2,
                                           ldcoeffs2, r1/r2, tse, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=False)

        tempt1 = tempt1 % period
        tempt2 = tempt2 % period
        tempres1 = (tempres1 - 1.)/(1. + frat) + 1.
        tempres2 = (tempres2 - 1.)/(1. + 1./frat) + 1.

        sorting1 = np.argsort(tempt1)
        sorting2 = np.argsort(tempt2)

        tempres1 = tempres1[sorting1]
        tempt1 = tempt1[sorting1]
        tempres2 = tempres2[sorting2]
        tempt2 = tempt2[sorting2]

        #not including crowdsap term.
        #tempres1 = (tempres1 + frat) / (1.+frat)
        #tempres2 = (tempres2 * frat + 1.) / (1. + frat)
        totpol, totmod = np.ones(len(jd)), np.ones(len(jd))

        if polyorder>0:
    
            # mask out continuum data
            #clip = ((abs(phase)<1.5*self.pwidth) | (abs(phase-self.sep) < 1.5*self.swidth))
            clip = (jd>0)
    
            chunk = np.array(np.where(np.diff(jd[clip]) > self.pwidth*period))[0]
            #put in dummy first and last element # placeholders
            
            chunk = np.append(chunk, len(jd[clip])-2).flatten()
            _, chunk3 = np.unique(np.searchsorted(jd[clip][chunk], jd), return_index=True)
            chunk=chunk3
            chunk[-1]+=1
    #        plt.plot(self.jd, self.flux, 'ro', self.jd[clip], self.flux[clip], 'go', self.jd[chunk[:-1]],self.flux[chunk[:-1]], 'mo')
    #        plt.show()
            for i in range(len(chunk)-1):
                #print i, chunk[i], chunk[i+1], self.jd[chunk[i]:chunk[i+1]]
                t = jd[chunk[i]:chunk[i+1]]
                f = flux[chunk[i]:chunk[i+1]]
                ef = dflux[chunk[i]:chunk[i+1]]
                crow = crowd[chunk[i]:chunk[i+1]]
                maf = rsky(e, period, t0, 1e-8, t)
                npts=len(maf)
                #use this version for full lightcurve treatment...
                r = a*(1.-e**2) / (1.+e*np.cos(maf))
                zcomp = np.sin(omega+maf) * np.sin(inc) 
                #z = r*np.sqrt(1.-zcomp**2)
                pe = ((r*zcomp>0.))# & (z <= 1.05*(r1+r2)*r2au))
                se = ((r*zcomp<0.))# & (z <= 1.05*(r1+r2)*r2au))
                model = np.ones(npts)
    #            sse = (((maf+omega) % (TWOPI))>np.pi)
    #            ppe = (((maf+omega) % (TWOPI))<=np.pi)
                
    #            plt.plot(t, f, 'ro', t[pe], f[pe], 'go', t[se], f[se], 'bo')
    #            plt.title(str(i))
    #            plt.show()
                if pe.any():
#                    shift = period * np.round((np.mean(t[pe]) - tpe)/period)
                    model[pe] = np.interp(t[pe]%period, tempt1, tempres1)
                    model[pe] = (model[pe] - 1.) * crow[pe] + 1.
    #                print "PE: mean(t[pe]), tpe, (mean(t[pe])-tpe)/period, round ver"
    #                print np.mean(t[pe]), tpe, (np.mean(t[pe]) - tpe)/period, np.round((np.mean(t) - tpe)/period), len(t[pe]), len(f[pe]), len(tempt1), len(tempres1), len(model[pe])
    #                plt.plot(t[pe]-shift, f[pe], 'ro', tempt1, (tempres1-1.)*crow[pe][0] + 1., 'bo', t[pe]-shift, model[pe], 'go')
    #                plt.title('pe')
    #                plt.show()
    #                plt.close('all')
                if se.any():
#                    shift = period * np.round((np.mean(t[se]) - tse)/period)
                    model[se] = np.interp(t[se]%period, tempt2, tempres2)
                    model[se] = (model[se] - 1.) * crow[se] + 1.
    
    #                print "SE"
    #                print np.mean(t[se]), tse, (np.mean(t[se]) - tse)/period, np.round((np.mean(t[se]) - tse)/period)
    #                plt.plot(t[se]-shift, f[se], 'ro', tempt2, (tempres2-1.)*crow[se][0] + 1., 'bo', t[se]-shift, model[se], 'go')
    #                plt.title('se')
    #                plt.show()
    #                plt.close('all')
    #            else:
    #                print "This data bundle does not belong to SE or PE"
                # marginalization (2nd order polynomial fit to residuals)

                bad = (model<1)
                tt = t[~bad]
                mmodel = model[~bad]
                ff = f[~bad]
                eef = ef[~bad]
                nnpts = len(ff)
                tnew = tt - np.mean(tt)
                #if len(t[~bad]) < 1:
                    #print "Npts ooe = ",len(t[~bad])
    
                # Bk = sum over i (D_i/M_i)(tdiff_i)^k / (sigma_i/M_i)^2
                # matrix 3 rows x npts columns since quadratic polynomial 
                # fit requires 3 coeffs
                if bad[0] or bad[-1]:
                    polyorder=1
                #number of 'i' data or model points; polynomial order
                order_pow = np.arange(polyorder+1)
                t_pow = tnew[:,np.newaxis]**order_pow
                Bk = np.ones(shape=(polyorder+1,nnpts))*((ff-mmodel)/(eef)**2)
                Bk*=t_pow.T
                #sum along 'i' (or along each row)
                Bksum = np.sum(Bk,axis=1)
                #Mjk = sum over i (tdiff_i)^j * (tdiff_i)^k / (sigma_i/M_i)^2
                #construct 3 rows x npts columns 
                Mj = np.ones(shape=(polyorder+1,nnpts))/(eef)**2
                Mj*=t_pow.T
                #transform from 2D (j rows x i columns) to 3D (k x j x i)
                t_pow_3d = tnew[:,np.newaxis,np.newaxis]**order_pow
                Mjk = t_pow_3d.T * Mj[np.newaxis,:,:]
                #now sum along 'i' 
                Mjksum = np.sum(Mjk,axis=2)
                #do matrix inversion solver thing to get polynomial coeffs
                try:
                    Aj = np.linalg.lstsq(Mjksum,Bksum)[0]
                    pol = np.polyval(Aj[::-1],t-np.mean(t))
                except: 
                    pol = np.ones(npts)
                #Aj = np.dot(np.linalg.pinv(Mjksum), Bksum)
    #                plt.plot(t, f, 'ro', t, model*pol, 'go')
    #                plt.plot(t, pol, 'ms', tt, np.polyval(Aj[::-1],tnew), 'cs')
    #                plt.show()
    
                totmod[chunk[i]:chunk[i+1]] = model
                totpol[chunk[i]:chunk[i+1]] = pol
        else:
            maf = rsky(e, period, t0, 1e-8, jd)
            r = a*(1.-e**2) / (1.+e*np.cos(maf))
            zcomp = np.sin(omega+maf) * np.sin(inc) 
            #z = r*np.sqrt(1.-zcomp**2)
            pe = ((r*zcomp>0.)) #& (z <= 1.05*(r1+r2)*r2au))
            se = ((r*zcomp<0.)) #& (z <= 1.05*(r1+r2)*r2au))
            tt = jd % period
            if pe.any():
                totmod[pe] = np.interp(tt[pe], tempt1, tempres1)
                totmod[pe] = (totmod[pe] - 1.) * crowd[pe] + 1.
            if se.any():
                totmod[se] = np.interp(tt[se], tempt2, tempres2)
                totmod[se] = (totmod[se] - 1.) * crowd[se] + 1.
        #     if np.sum(totmod[se]-1.) == 0.:
        #         return np.ones_like(totmod), totpol
        # if np.sum(totmod-1.) == 0.:
        #   return totmod, totmod
        return totmod, totpol

    def rvprep(self, t, rv1, rv2, drv1, drv2):
        """Stores observed radial velocity data points

        Parameters
        ----------
        t : float array or scalar
            times of observations
        rv1 : float array or scalar
            RV of primary in m/s
        rv2 : float array or scalar
            RV of secondary in m/s
        dr1 : float array of scalar
            RV err of primary in m/s
        dr2 : float array or scalar
            RV err of secondary in m/s

        Returns
        -------
        m1, m2, k0 : tuple
            guesses for the masses and systemic velocity of the binary from the RV semi-amplitudes
        """
        self.rv1_obs = rv1
        self.rv2_obs = rv2
        self.rv1_err_obs = drv1
        self.rv2_err_obs = drv2
        self.rv_t = t
        self.bad1 = np.isnan(self.rv1_obs)
        self.bad2 = np.isnan(self.rv2_obs)
        k1 = (np.nanmax(self.rv1_obs) - np.nanmin(self.rv1_obs))/2.
        k2 = (np.nanmax(self.rv2_obs) - np.nanmin(self.rv2_obs))/2.
        k0 = np.nanmedian(np.append(self.rv1_obs, self.rv2_obs))
        try:
            m1, m2 = self.rvpars_guess_mass(k1, k2, self.pars['period'],
                                        np.sqrt(self.pars['esinw']**2 + self.pars['ecosw']**2))
        except:
            print("Could not compute m1, m2 b/c esinw, ecosw not loaded to keblat.pars")
            m1, m2 = 1.0, 1.0
        return m1, m2, k0

    def _rvfit_old(self, rvpars, t):
        m1, m2, period, tpe, esinw, ecosw, inc, k0, rverr = rvpars

        a = ((period / d2y) ** 2 * (m1+m2)) ** (1. / 3.)
        e = np.sqrt(esinw**2+ecosw**2)
        omega = np.arctan2(esinw, ecosw)

        fpe = np.pi/2. - omega

        t0 = tpe - (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
            (e*np.sin(fpe)/(1.+e*np.cos(fpe)) - 2.*(1.-e**2)**(-0.5) * \
            np.arctan(np.sqrt(1.-e**2) * np.tan((fpe)/2.) / (1.+e)))


        maf = rsky(e, period, t0, 1e-8, t)

        # egative sign b/c observers' ref. frame is flipped
        vr2 = -np.sqrt(8.875985e12 * m1**2 / (m1+m2) / (a * (1-e**2))) * np.sin(inc) * \
              (np.cos(omega+maf) + e * np.cos(omega))
        omega+=np.pi # periapse of primary is 180 offset
        vr1 = -np.sqrt(8.875985e12 * m2**2 / (m1+m2) / (a * (1-e**2))) * np.sin(inc) * \
              (np.cos(omega+maf) + e * np.cos(omega))
        return vr1+k0, vr2+k0

    def rvfit(self, rvpars, t):
        """Computes the radial velocities of each binary component.

        Parameters
        ----------
        rvpars : float array or list
                msum, mrat, period, tpe, esinw, ecosw, inc, k0, rverr
        t : float array or scalar
                times of observations to compute RV

        Returns
        -------
        vr1, vr2 : tuple
                RVs in m/s, where vr1 and vr2 are of shape/type of input t
        """

        msum, mrat, period, tpe, esinw, ecosw, inc, k0, rverr = rvpars
        e = np.sqrt(esinw**2+ecosw**2)
        omega = np.arctan2(esinw, ecosw)

        fpe = np.pi/2. - omega

        m1, m2 = self.sumrat_to_12(msum, mrat)

        self.updatepars(msum=msum, mrat=mrat, m1=m1, m2=m2, period=period, tpe=tpe,
                        esinw=esinw, ecosw=ecosw, k0=k0, rverr=rverr, inc=inc)

        t0 = tpe - (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
            (e*np.sin(fpe)/(1.+e*np.cos(fpe)) - 2.*(1.-e**2)**(-0.5) * \
            np.arctan(np.sqrt(1.-e**2) * np.tan((fpe)/2.) / (1.+e)))

        maf = rsky(e, period, t0, 1e-8, t)
        amp = 29794.509 / np.sqrt(1-e**2) * (period/d2y)**(-1/3.) / (m1+m2)**(2/3.)
        vr2 = -amp * m1 * np.sin(inc) * \
              (np.cos(omega+maf) + e * np.cos(omega))
        omega+=np.pi # periapse of primary is 180 offset
        vr1 = -amp * m2 * np.sin(inc) * \
              (np.cos(omega+maf) + e * np.cos(omega))
        return vr1+k0, vr2+k0


    def lcfit(self, lcpars, jd, quarter, flux, dflux, crowd,
              polyorder=2, ooe=True):
        """Computes light curve model
        
        Parameters
        ----------
        lcpars : float array
            parameters for LC fitting: 
            msum, rsum, rratio, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4
        jd : float array
            time array
        quarter : float array
            corresponding kepler quarter for a given time
        flux : float array
            observed flux
        dflux : float array
            flux error
        crowd : float array
            array of crowding values (additional flux)
        polyorder : int
            order of polynomial to detrend lightcurve
        
        Returns
        -------
        totmod : float array
            array of model fluxes
        totpol : float array
            array of polynomials for detrending
        """
        # r1, r2, frat derive from m1, m2, z0, t0, dist, E(B-V), scaleheight
        msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, \
            q1, q2, q3, q4 = lcpars

        # LD transformations (Kipping 2013)
        c1 = 2.*np.sqrt(q1)*q2
        c2 = np.sqrt(q1)*(1.-2.*q2)
        c3 = 2.*np.sqrt(q3)*q4
        c4 = np.sqrt(q3)*(1.-2.*q4)
        ldcoeffs1 = np.array([c1, c2])
        ldcoeffs2 = np.array([c3, c4])
            
#        if r2 > r1:
#            r1, r2 = r2, r1
#            m1, m2 = m2, m1
#            frat = 1./frat
        omega=np.arctan2(esinw,ecosw)
        e=np.sqrt(esinw**2+ecosw**2)

        # nip it at the bud.
        if (e>=1.):
            #print "e>=1", e
            return -np.inf, -np.inf
            
        # r1 = rsum/(1.+rrat)
        # r2 = rsum/(1.+1./rrat)
        r1, r2 = self.sumrat_to_12(rsum, rrat)
        a = self.get_a(period, msum)
        inc = self.get_inc(b, r1, a)
        #inc = np.arccos(b*r1/(a/r2au))
        
        if np.isnan(inc):
            #print "inc is nan", inc
            return -np.inf, -np.inf

        self.updatepars(msum=msum, rsum=rsum, rrat=rrat, period=period, tpe=tpe,
                       esinw=esinw, ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4,
                       frat=frat, r1=r1, r2=r2, inc=inc, e=e)
        fpe = np.pi/2. - omega
        fse = -np.pi/2. - omega

        # transform time of center of PE to time of periastron (t0)
        # from Eq 9 of Sudarsky et al (2005)
        t0 = tpe - self.sudarsky(fpe, e, period)
        tse = t0 + self.sudarsky(fse, e, period)
        # t0 = tpe - (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
        #     (e*np.sin(fpe)/(1.+e*np.cos(fpe)) - 2.*(1.-e**2)**(-0.5) * \
        #     np.arctan(np.sqrt(1.-e**2) * np.tan((fpe)/2.) / (1.+e)))
        # tse = t0 + (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
        #     (e*np.sin(fse)/(1.+e*np.cos(fse)) - 2.*(1.-e**2)**(-0.5) * \
        #     np.arctan(np.sqrt(1.-e**2) * np.tan((fse)/2.) / (1.+e)))
        self.tpe = tpe
        self.tse = tse
        # if tse<tpe:
        #     tse+=period
            
        tempt1, tempres1 = self.lctemplate(lcpars, period, omega, e, a, inc, r1,
                                           ldcoeffs1, r2/r1, tpe, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=True)

        tempt2, tempres2 = self.lctemplate(lcpars, period, omega, e, a, inc, r2,
                                           ldcoeffs2, r1/r2, tse, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=False)

        tempt1 = tempt1 % period
        tempt2 = tempt2 % period
        tempres1 = (tempres1 - 1.)/(1. + frat) + 1.
        tempres2 = (tempres2 - 1.)/(1. + 1./frat) + 1.

        sorting1 = np.argsort(tempt1)
        sorting2 = np.argsort(tempt2)

        tempres1 = tempres1[sorting1]
        tempt1 = tempt1[sorting1]
        tempres2 = tempres2[sorting2]
        tempt2 = tempt2[sorting2]

        #not including crowdsap term.
        #tempres1 = (tempres1 + frat) / (1.+frat)
        #tempres2 = (tempres2 * frat + 1.) / (1. + frat)
        totmod, totpol = np.ones(len(jd)), np.ones(len(jd))

        maf = rsky(e, period, t0, 1e-8, jd)
        r = a*(1.-e**2) / (1.+e*np.cos(maf))
        zcomp = np.sin(omega+maf) * np.sin(inc) 
        pe = ((r*zcomp>0.)) #& (z <= 1.05*(r1+r2)*r2au))
        se = ((r*zcomp<0.)) #& (z <= 1.05*(r1+r2)*r2au))
        tt = jd % period
        if pe.any():
            totmod[pe] = np.interp(tt[pe], tempt1, tempres1)
            totmod[pe] = (totmod[pe] - 1.) * crowd[pe] + 1.
        if se.any():
            totmod[se] = np.interp(tt[se], tempt2, tempres2)
            totmod[se] = (totmod[se] - 1.) * crowd[se] + 1.

        if polyorder>0:
            if (self.sep-self.clip_tol*(self.pwidth+self.swidth) < self.pwidth):
                chunk = np.array(np.where(np.diff(jd) > np.median(np.diff(jd))*4.))[0]
            else:
                chunk = np.array(np.where(np.diff(jd) > self.pwidth*period))[0]
            #put in dummy first and last element # placeholders
            
            chunk = np.append(chunk, len(jd)-2).flatten()
            _, chunk3 = np.unique(np.searchsorted(jd[chunk], jd), return_index=True)
            chunk=chunk3
            chunk[-1]+=1
            chunk = np.unique(np.sort(np.append(chunk, np.where(np.diff(quarter)>0)[0]+1)))
            totpol = poly_lc_cwrapper(jd, flux, dflux, totmod, chunk, porder=polyorder, ooe=ooe)
#            phase = ((jd - tpe) % period) / period
#            sorting = np.argsort(phase)
#            nopoly = (totpol[sorting] == 1.)
#            if (np.sum(nopoly)>0) and (np.sum(nopoly)<len(totpol)*0.1):
#                _totpol = totpol[sorting]
#                tmp = np.interp(phase[sorting][nopoly], phase[sorting][~nopoly], flux[sorting][~nopoly]/totpol[sorting][~nopoly])
#                #print np.sum(nopoly), np.sum(~nopoly)
#                _totpol[nopoly] = flux[sorting][nopoly] / tmp
#                totpol[sorting] = _totpol
        return totmod, totpol

    def _lcfit(self, lcpars, jd, quarter, flux, dflux, crowd,
              polyorder=2, ooe=True, flares=None):
        """Computes light curve model
        
        Parameters
        ----------
        lcpars : float array
            parameters for LC fitting: 
            msum, rsum, rratio, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4
        jd : float array
            time array
        quarter : float array
            corresponding kepler quarter for a given time
        flux : float array
            observed flux
        dflux : float array
            flux error
        crowd : float array
            array of crowding values (additional flux)
        polyorder : int
            order of polynomial to detrend lightcurve
        
        Returns
        -------
        totmod : float array
            array of model fluxes
        totpol : float array
            array of polynomials for detrending
        """
        # r1, r2, frat derive from m1, m2, z0, t0, dist, E(B-V), scaleheight
        msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, \
            q1, q2, q3, q4 = lcpars

        # LD transformations (Kipping 2013)
        c1 = 2.*np.sqrt(q1)*q2
        c2 = np.sqrt(q1)*(1.-2.*q2)
        c3 = 2.*np.sqrt(q3)*q4
        c4 = np.sqrt(q3)*(1.-2.*q4)
        ldcoeffs1 = np.array([c1, c2])
        ldcoeffs2 = np.array([c3, c4])
            
#        if r2 > r1:
#            r1, r2 = r2, r1
#            m1, m2 = m2, m1
#            frat = 1./frat
        omega=np.arctan2(esinw,ecosw)
        e=np.sqrt(esinw**2+ecosw**2)

        # nip it at the bud.
        if (e>=1.):
            #print "e>=1", e
            return -np.inf, -np.inf
            
        # r1 = rsum/(1.+rrat)
        # r2 = rsum/(1.+1./rrat)
        r1, r2 = self.sumrat_to_12(rsum, rrat)
        a = self.get_a(period, msum)
        inc = self.get_inc(b, r1, a)
        #inc = np.arccos(b*r1/(a/r2au))
        
        if np.isnan(inc):
            #print "inc is nan", inc
            return -np.inf, -np.inf

        self.updatepars(msum=msum, rsum=rsum, rrat=rrat, period=period, tpe=tpe,
                       esinw=esinw, ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4,
                       frat=frat, r1=r1, r2=r2, inc=inc)
        fpe = np.pi/2. - omega
        fse = -np.pi/2. - omega

        # transform time of center of PE to time of periastron (t0)
        # from Eq 9 of Sudarsky et al (2005)
        t0 = tpe - self.sudarsky(fpe, e, period)
        tse = t0 + self.sudarsky(fse, e, period)
        # t0 = tpe - (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
        #     (e*np.sin(fpe)/(1.+e*np.cos(fpe)) - 2.*(1.-e**2)**(-0.5) * \
        #     np.arctan(np.sqrt(1.-e**2) * np.tan((fpe)/2.) / (1.+e)))
        # tse = t0 + (-np.sqrt(1.-e**2) * period / (2.*np.pi)) * \
        #     (e*np.sin(fse)/(1.+e*np.cos(fse)) - 2.*(1.-e**2)**(-0.5) * \
        #     np.arctan(np.sqrt(1.-e**2) * np.tan((fse)/2.) / (1.+e)))
        self.tpe = tpe
        self.tse = tse
        # if tse<tpe:
        #     tse+=period
            
        tempt1, tempres1 = self.lctemplate(lcpars, period, omega, e, a, inc, r1,
                                           ldcoeffs1, r2/r1, tpe, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=True)

        tempt2, tempres2 = self.lctemplate(lcpars, period, omega, e, a, inc, r2,
                                           ldcoeffs2, r1/r2, tse, t0,
                                           cadence = self.cadence,
                                           exp = self.exp, pe=False)

        tempt1 = tempt1 % period
        tempt2 = tempt2 % period
        tempres1 = (tempres1 - 1.)/(1. + frat) + 1.
        tempres2 = (tempres2 - 1.)/(1. + 1./frat) + 1.

        sorting1 = np.argsort(tempt1)
        sorting2 = np.argsort(tempt2)

        tempres1 = tempres1[sorting1]
        tempt1 = tempt1[sorting1]
        tempres2 = tempres2[sorting2]
        tempt2 = tempt2[sorting2]

        #not including crowdsap term.
        #tempres1 = (tempres1 + frat) / (1.+frat)
        #tempres2 = (tempres2 * frat + 1.) / (1. + frat)
        totmod, totpol = np.ones(len(jd)), np.ones(len(jd))

        maf = rsky(e, period, t0, 1e-8, jd)
        r = a*(1.-e**2) / (1.+e*np.cos(maf))
        zcomp = np.sin(omega+maf) * np.sin(inc) 
        pe = ((r*zcomp>0.)) #& (z <= 1.05*(r1+r2)*r2au))
        se = ((r*zcomp<0.)) #& (z <= 1.05*(r1+r2)*r2au))
        tt = jd % period

        if pe.any():
            totmod[pe] = np.interp(tt[pe], tempt1, tempres1)
            totmod[pe] = (totmod[pe] - 1.) * crowd[pe] + 1.
            
        if se.any():
            totmod[se] = np.interp(tt[se], tempt2, tempres2)
            totmod[se] = (totmod[se] - 1.) * crowd[se] + 1.
        
        if polyorder>0:
            mult=1.7
            clip = (abs((tt - tpe%period))<self.pe_dur%period*mult) | (abs((tt - tpe%period))>period-self.pe_dur%period*mult) | (abs((tt - tse%period))<self.se_dur%period*mult) | (abs((tt - tse%period))>period-self.se_dur%period*mult)
            chunk = np.where(abs(np.diff(clip.astype(int)))>0)[0]
            chunk = np.append(chunk, 0)
            chunk = np.unique(np.sort(np.append(chunk, len(tt))))
            self.chunk = chunk
            if flares is None:
                totpol = poly_lc_cwrapper(jd, flux, dflux, totmod, chunk, porder=polyorder, ooe=ooe)
            else:
                totpol = np.ones(len(jd))
                for ii in range(len(chunk)-1):
                    bad = (totmod[chunk[ii]:chunk[ii+1]] < 1) | flares[chunk[ii]:chunk[ii+1]]
                    tt = jd[chunk[ii]:chunk[ii+1]][~bad]
                    mmodel = totmod[chunk[ii]:chunk[ii+1]][~bad]
                    ff = flux[chunk[ii]:chunk[ii+1]][~bad]
                    eef = dflux[chunk[ii]:chunk[ii+1]][~bad]
                    tnew = tt - np.mean(tt)
                    nnpts = len(ff)
                    npts = len(bad)
                    if bad[0] or bad[-1]:
                        poly_remember = polyorder
                        polyorder=1
                    order_pow = np.arange(polyorder+1)
                    t_pow = tnew[:,np.newaxis]**order_pow
                    Bk = np.ones(shape=(polyorder+1,nnpts))*((ff/mmodel)/(eef/mmodel)**2)
                    Bk*=t_pow.T
                    #sum along 'i' (or along each row)
                    Bksum = np.sum(Bk,axis=1)
                    #Mjk = sum over i (tdiff_i)^j * (tdiff_i)^k / (sigma_i/M_i)^2
                    #construct 3 rows x npts columns 
                    Mj = np.ones(shape=(polyorder+1,nnpts))/(eef/mmodel)**2
                    Mj*=t_pow.T
                    #transform from 2D (j rows x i columns) to 3D (k x j x i)
                    t_pow_3d = tnew[:,np.newaxis,np.newaxis]**order_pow
                    Mjk = t_pow_3d.T * Mj[np.newaxis,:,:]
                    #now sum along 'i' 
                    Mjksum = np.sum(Mjk,axis=2)
                    #do matrix inversion solver thing to get polynomial coeffs
                    try:
                        Aj = np.linalg.lstsq(Mjksum,Bksum)[0]
                        pol = np.polyval(Aj[::-1],jd[chunk[ii]:chunk[ii+1]]-np.mean(jd[chunk[ii]:chunk[ii+1]]))
                    except: 
                        pol = np.ones(npts)
                    #Aj = np.dot(np.linalg.pinv(Mjksum), Bksum)
        #                plt.plot(t, f, 'ro', t, model*pol, 'go')
        #                plt.plot(t, pol, 'ms', tt, np.polyval(Aj[::-1],tnew), 'cs')
        #                plt.show()
                    if bad[0] or bad[-1]:
                        polyorder = poly_remember
                    totpol[chunk[ii]:chunk[ii+1]] = pol         
#            phase = ((jd - tpe) % period) / period
#            sorting = np.argsort(phase)
#            nopoly = (totpol[sorting] == 1.)
#            if (np.sum(nopoly)>0) and (np.sum(nopoly)<len(totpol)*0.1):
#                _totpol = totpol[sorting]
#                tmp = np.interp(phase[sorting][nopoly], phase[sorting][~nopoly], flux[sorting][~nopoly]/totpol[sorting][~nopoly])
#                #print np.sum(nopoly), np.sum(~nopoly)
#                _totpol[nopoly] = flux[sorting][nopoly] / tmp
#                totpol[sorting] = _totpol
        return totmod, totpol

    def _getvals_defunct2(self, fit_params, partype='allpars'):
        """Grabs the values of input

        Parameters
        ----------
        fit_params : float array or dict
                    if dict, updates keblat parameters and returns pars according to input keys
                    if array, returns array

        """
        if type(fit_params) is dict:
            self.updatepars(**fit_params)
            return self.getpars(partype=partype)
        elif type(fit_params) is np.ndarray:
            return fit_params
        else:
            self.updatepars(**fit_params.valuesdict())
            return self.getpars(partype=partype)

#    def _getvals_defunct(self, fit_params, ntotpars, lctype):
#        """Fetch values from parameters
#
#        Parameters
#        ----------
#        fit_params : dict or float array
#            can input lmfit's Parameter class or just an array of parameter vals
#        ntotpars : int
#            total number of parameters to be fit
#        lcpars : Boolean
#            if True, returns set of parameters for light-curve fitting only
#            if False, returns all parameters
#
#        Returns
#        -------
#        guess : float array
#            array of values for each parameter
#        """
#
#        guess=np.empty(ntotpars)
#        for jj in range(ntotpars):
#            try:
#                guess[jj] = fit_params[self.pars.keys()[jj]].value
#            except KeyError:
#                guess[jj] = self.pars.values()[jj]
#            except ValueError:
#                guess[jj] = fit_params[jj]
#            except IndexError:
#                guess[jj] = fit_params[jj]
#        if lcpars:
#            return np.array([guess[0], guess[1], self.r1+self.r2,
#                                           self.r2/self.r1, guess[7], guess[8],
#                                           guess[9], guess[10], guess[11],
#                                            self.frat, guess[12], guess[13],
#                                            guess[14], guess[15]])
#        return guess

    def ilnlike(self, fisopars, lc_constraints=[], ebv_dist=None, ebv_arr=None,
                residual=False, retpars=False):
        """Computes log likelihood of isochrone fit portion of KEBLAT
        
        Parameters
        ----------
        fisopars : dict or float array
            either from lmfit parameter class or just an array of vals
        residual: boolean
            True if want to return residual array
            False if return loglikelihood val
        
        Returns
        -------
        loglike : float
            returns -np.inf if model mags have invalid values
        isores : float array
            if residual = True
        
        """
        
        parnames = np.array(['m1', 'm2', 'z0', 'age', 'dist',
                             'ebv', 'h0', 'isoerr'])
        isopars=np.empty(8)
        for jj in range(8):
            try:
                isopars[jj] = fisopars[parnames[jj]].value
            except KeyError:
                isopars[jj] = 0.0
            except ValueError:
                isopars[jj] = fisopars[jj]
            except IndexError:
                isopars[jj] = fisopars[jj]
        if ebv_arr is not None:
            isopars[5] = np.interp(isopars[4], ebv_dist, ebv_arr)
        if retpars:
            return isopars
        # m1 = isopars[0] / (1. + isopars[1])
        # m2 = isopars[0] / (1. + 1./isopars[1])
        # isopars[0], isopars[1] = m1, m2
        # print m1, m2, isopars[0], isopars[1]
        isoerr = np.exp(isopars[-1])
        magsmod = self.isofit(isopars)
 #/ np.sqrt(self.emagsobs**2 + isoerr**2)
        if np.any(np.isinf(magsmod)):
            if residual:
                return np.ones(len(self.magsobs) + len(lc_constraints))*1e12
            return -np.inf, str((0,0,0))

        lc_inputs = np.array([(self.r1+self.r2)/(isopars[0]+isopars[1])**(1./3.), self.r2/self.r1, self.frat])
        lc_uncertainty = np.array([0.002, 0.002, 0.002])
        if np.any(np.isinf(lc_inputs)):
            if residual:
                return np.ones(len(self.magsobs) + len(lc_constraints))*1e12
            else:
            	return -np.inf, str((0,0,0))

        isores = np.concatenate(((magsmod - self.magsobs) / np.sqrt(self.emagsobs**2 + isoerr**2),
                                 (lc_inputs-lc_constraints)/(lc_uncertainty*lc_constraints))) #/ np.sqrt(self.emagsobs**2 + isoerr**2)
        for ii, dii, jj in zip([self.armstrongT1, self.armstrongT2],
                               [self.armstrongdT1, self.armstrongdT2],
                               [10**self.temp1, 10**self.temp2]):
            if ii is not None:
                if dii is None:
                    dii=0.05*ii
                isores = np.append(isores, (ii-jj)/dii)
        if residual:
            return isores

        chisq = np.sum(isores**2) + np.sum(np.log(TWOPI * \
                    (self.emagsobs**2 + isoerr**2)))
        chisq += ((isopars[5] - self.ebv)/(self.debv))**2
        #chisq += ((isopars[6] - 119.)/(15.))**2
        return -0.5*chisq, str((self.r1, self.r2, self.frat))
    
    def ilnprior(self, isopars):
        """Returns log-prior for isochrone fitting portion only
        
        Parameters
        ----------
        isopars : float array
            array of val for (m1, m2, z0, age, dist, ebv, h0, isoerr)        
        
        Returns
        lnp : float
            0 if within flat priors, -np.inf if outside
        """
        #m1, m2, z0, age, dist, ebv, h0, lnisoerr = isopars
        bounds = np.array([(.1, 12.), (.1, 12.), (0.001, 0.06), (6., 10.1),
                           (10., 15000.), (0.0, 1.0), (119-20., 119+20.),
                            (-8, 1.)])
        pcheck = np.all((np.array(isopars) >= bounds[:, 0]) & \
                        (np.array(isopars) <= bounds[:, 1]))
        if pcheck:
            return 0.0 + isopars[3]*np.log(10.) + np.log(np.log(10.))
        else:
            return -np.inf
    
    def ilnprob(self, isopars, lc_constraints=None):
        """Returns log probability (loglike + logprior)
        
        Parameters
        ----------
        isopars : float array
            array of val for (m1, m2, z0, age, dist, ebv, h0, isoerr)        
        
        Returns
        lnprob : float
            -np.inf if invalid values     
        """
        lp = self.ilnprior(isopars)
        if np.isinf(lp):
            return -np.inf, str((-np.inf, -np.inf, -np.inf))
        ll, blobs = self.ilnlike(isopars, lc_constraints=lc_constraints)
        if (np.isnan(ll) or np.isinf(ll)):
            return -np.inf, str((-np.inf, -np.inf, -np.inf))
        return lp + ll, blobs

    def lnlike(self, allpars, lc_constraints=None, qua=[1], polyorder=2,
               residual=False, ld4=False, clip=None, ooe=True):
        """Returns loglikelihood for both SED + LC fitting
        
        Parameters
        ----------
        allpars : float array
            array of parameter values
        qua : list
            list of quarters to include in LC analysis; default = [1]
        polyorder : int
            order of polynomial for detrending; default = 2
        residual : boolean
            True if want to return (weighted) residual array
        ld4 : boolean
            True if want quartic non-linear LD law; default = False (quadratic)
            
        Returns
        -------
        lili : float
            log likelihood value of model; returns -np.inf if invalid
        res : float array
            residuals of model and data if residuals = True
        """
        
        if ld4:
            self.ldtype=0
#            allpars = self.getvals(fitpars, 22)
            m1, m2, z0, age, dist, ebv, h0, period, tpe, esinw, \
                ecosw, b, q11, q12, q13, q14, q21, q22, q23, q24, \
                lcerr, isoerr  = allpars
            ldcoeffs = np.array([q11, q12, q13, q14, q21, q22, q23, q24])        
        else:
            self.ldtype=1
#            allpars = self.getvals(fitpars, 18)
            m1, m2, z0, age, dist, ebv, h0, period, tpe, esinw, \
                ecosw, b, q1, q2, q3, q4, lcerr, isoerr = allpars
            ldcoeffs = np.array([q1, q2, q3, q4])
            # m1 = msum / (1. + mrat)
            # m2 = msum / (1. + 1./mrat)
            self.updatepars(m1=m1, m2=m2, z0=z0, age=age, dist=dist, ebv=ebv,
                            h0=h0, period=period, tpe=tpe, esinw=esinw, 
                            ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4, 
                            lcerr=lcerr, isoerr=isoerr, msum=m1+m2, mrat=m2/m1)
        # do isochrone matching first
        isopars = [m1, m2, z0, age, dist, ebv, h0, isoerr]
        magsmod = self.isofit(isopars)
        if np.any(np.isinf(magsmod)):
            return -np.inf #/ np.sqrt(self.emagsobs**2 + isoerr**2)

        lc_inputs = np.array([(self.r1+self.r2)/(m1+m2)**(1./3.), self.r2/self.r1, self.frat])

        if np.any(np.isinf(lc_inputs)):
            return -np.inf
        if lc_constraints is None:
            lc_priors = np.array([])
        else:
            lc_priors = (lc_inputs-lc_constraints)/(0.03*lc_constraints)
        isores = np.concatenate(((magsmod - self.magsobs) / np.sqrt(self.emagsobs**2 + isoerr**2),
                                 lc_priors)) #/ np.sqrt(self.emagsobs**2 + isoerr**2)
        for ii, dii, jj in zip([self.armstrongT1, self.armstrongT2],
                               [self.armstrongdT1, self.armstrongdT2],
                               [10**self.temp1, 10**self.temp2]):
            if ii is not None:
                if dii is None:
                    dii=0.05*ii
                isores = np.append(isores, (ii-jj)/dii)
        chisq = np.sum(isores**2) + np.sum(np.log((self.emagsobs**2 + isoerr**2)))

        #now the light curve fitting part
        lcpars = np.concatenate((np.array([m1+m2, self.r1+self.r2, 
                                           self.r2/self.r1, period, 
                                           tpe, esinw, ecosw, b, self.frat]),
                                           ldcoeffs))

        # conds = [self.quarter == ii for ii in np.array(qua)]
        # conds = np.sum(np.array(conds), axis=0)
        # conds = np.array(conds, dtype=bool)
        #self.clip2 = (self.clip) #* conds
        if clip is None:
            clip = self.clip
        lcmod, lcpol = self.lcfit(lcpars, self.jd[clip],
                                  self.quarter[clip], self.flux[clip],
                                self.fluxerr[clip], self.crowd[clip],
                                polyorder=polyorder, ooe=ooe)
        lcres = (lcmod*lcpol - self.flux[clip]) / np.sqrt(self.fluxerr[clip]**2 + lcerr**2)

        if np.any(np.isinf(lcmod)):
            return -np.inf
        totres = np.concatenate((isores, lcres))
        self.chi = np.sum(totres**2)
        
        chisq += np.sum(lcres**2) + \
                    np.sum(np.log((self.fluxerr[clip]**2 + lcerr**2)))
        chisq += ((isopars[5] - self.ebv)/(self.debv))**2
        chisq += ((isopars[2] - self.z0)/(0.2 * np.log(10) * self.z0))**2
        #chisq += ((isopars[6] - 119.)/(15.))**2        

        if residual:
            # print lc_priors, totres.shape, isores.shape, self.magsobs.shape, self.jd[self.clip].shape, \
            #     self.clip.sum(), clip.sum(), lcres.shape, lcmod.shape, lcpol.shape
            return totres
        lili = -0.5 * chisq
        return lili

    @staticmethod
    def rvpars_guess_mass(k1, k2, P, e):
        mrat = k1/k2
        msum = ((k1+k2)/29794.509)**3 * (1-e**2)**(3./2.) * (P/d2y)
        return msum/(1.+mrat), msum/(1.+1./mrat)

    def lnlike_lcrv(self, lcrvpars, qua=[1], polyorder=2, residual=False, clip=None, ooe=True):
        self.ldtype = 1
        #            allpars = self.getvals(fitpars, 18)
        msum, mrat, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4, lcerr, k0, rverr = lcrvpars
        m1, m2 = self.sumrat_to_12(msum, mrat)
        self.updatepars(m1=m1, m2=m2, period=period, tpe=tpe, esinw=esinw,
                        ecosw=ecosw, b=b, q1=q1, q2=q2, q3=q3, q4=q4,
                        lcerr=lcerr, k0=k0, rverr=rverr, frat=frat, 
                        msum=msum, mrat=mrat)

        lcpars = [msum, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4]
        # conds = [self.quarter == ii for ii in np.array(qua)]
        # conds = np.sum(np.array(conds), axis=0)
        # conds = np.array(conds, dtype=bool)
        # self.clip2 = (self.clip) #* conds
        if clip is None:
            clip = self.clip
        lcmod, lcpol = self.lcfit(lcpars, self.jd[clip],
                                  self.quarter[clip], self.flux[clip],
                                  self.fluxerr[clip], self.crowd[clip],
                                  polyorder=polyorder, ooe=ooe)
        lcres = (lcmod * lcpol - self.flux[clip]) / np.sqrt(self.fluxerr[clip] ** 2 + lcerr ** 2)

        if np.any(np.isinf(lcmod)):
            return -np.inf

        rvpars = [msum, mrat, period, tpe, esinw, ecosw, self.pars['inc'], k0, rverr]
        rv1, rv2 = self.rvfit(rvpars, self.rv_t)

        rvres = np.concatenate(((rv1[~self.bad1] - self.rv1_obs[~self.bad1]) /
                                np.sqrt(self.rv1_err_obs[~self.bad1] ** 2 + rverr ** 2),
                                (rv2[~self.bad2] - self.rv2_obs[~self.bad2]) /
                                np.sqrt(self.rv2_err_obs[~self.bad2] ** 2 + rverr ** 2)))

        totres = np.concatenate((lcres, rvres))
        self.chi = np.sum(totres ** 2)
        if residual:
            return totres

        chisq = np.sum(lcres ** 2) + \
                 np.sum(np.log((self.fluxerr[clip] ** 2 + lcerr ** 2)))
        chisq += np.sum(rvres ** 2) + \
                 np.sum(np.log((self.rv1_err_obs[~self.bad1] ** 2 + rverr ** 2))) + \
                 np.sum(np.log((self.rv2_err_obs[~self.bad2] ** 2 + rverr ** 2)))

        lili = -0.5 * chisq
        return lili

    def lnprior_lcrv(self, lcrvpars):
        #lcerr, rverr in log space
        msum, mrat, rsum, rrat, period, tpe, esinw, ecosw, b, frat, q1, q2, q3, q4, lcerr, k0, rverr = lcrvpars
        e = np.sqrt(esinw**2 + ecosw**2)
        pars2check = np.array([msum, mrat, rsum, rrat, period, tpe, e, b, frat, 
                               q1, q2, q3, q4, lcerr, k0, rverr])
        bounds = np.array([self.parbounds[ii] for ii in ['msum', 'mrat', 'rsum', 
                           'rrat', 'period', 'tpe', 'e', 'b', 'frat', 
                               'q1', 'q2', 'q3', 'q4', 'lcerr', 'k0', 'rverr']])
        bounds[-1,:] = [-8, 12]
        bounds[-3,:] = [-12, -2]
        pcheck = np.all((pars2check >= bounds[:,0]) & (pars2check <= bounds[:,1]))
        if pcheck:
            return 0.0
        else:
            return -np.inf
    
    def lnprior(self, allpars):
        """Returns logprior of SED + LC model
        
        Parameters
        ----------
        allpars : float array
            full list of parameters:
            (m1, m2, z0, age, dist, ebv, h0, period, tpe, 
            esinw, ecosw, b, q1, q2, q3, q4, lcerr, isoerr)
            
        Returns
        -------
        res : float
            0 if within flat priors, -np.inf if without
        """
        m1, m2, z0, age, dist, ebv, h0, period, tpe, esinw, ecosw, \
            b, q1, q2, q3, q4, lcerr, isoerr = allpars
        e = np.sqrt(esinw**2 + ecosw**2)
        pars2check = np.array([m1, m2, z0, age, dist, ebv, h0, \
            period, tpe, e, b, q1, q2, q3, q4, lcerr, isoerr])
        bounds = np.array([(.1, 12.), (.1, 12.), (0.001, 0.06), (6., 10.1),
                           (10., 15000.), (0.0, 1.0), (119-20., 119+20.),
                            (5., 3000.), (0., 1e8), (0., 0.99), (0., 3.), (0.,1.),
                            (0.,1.), (0.,1.), (0.,1.), (0., 0.05), (0., 0.3)])

        pcheck = np.all((pars2check >= bounds[:,0]) & \
                        (pars2check <= bounds[:,1]))

        if pcheck:
            return 0.0 + age*np.log(10.) + np.log(np.log(10.))
        else:
            return -np.inf
    
    def lnprob(self, allpars, qua=[1]):
        """Returns logprob of SED + LC model
        
        Parameters
        ----------
        allpars : float array
            full list of parameters:
            (m1, m2, z0, age, dist, ebv, h0, period, tpe, 
            esinw, ecosw, b, q1, q2, q3, q4, lcerr, isoerr)
        qua : list
            list of integers for Kepler quarters to be included in analysis
            
        Returns
        -------
        lnprob, blob : tuple 
            lnprob is lnprior + lnlike
            blob is a string of chi^2, r1, r2, f2/f1 values
        """
        
        lp = self.lnprior(allpars)
        if np.isinf(lp):
            return -np.inf, str((-np.inf, -np.inf, -np.inf, -np.inf))
        ll = self.lnlike(allpars, qua=qua)
        if (np.isnan(ll) or np.isinf(ll)):
            return -np.inf, str((-np.inf, -np.inf, -np.inf, -np.inf))
        return lp + ll, str((self.chi, self.r1, self.r2, self.frat))

    def lnprob_lcrv(self, lcrvpars, qua=[1]):
        lp = self.lnprior_lcrv(lcrvpars)
        if np.isinf(lp):
            return -np.inf
        lcrvpars[-1] = np.exp(lcrvpars[-1])
        lcrvpars[-3] = np.exp(lcrvpars[-3])
        ll = self.lnlike_lcrv(lcrvpars, qua=qua)
        if (np.isnan(ll) or np.isinf(ll)):
            return -np.inf
        return lp + ll
        
    def run_emcee(self, pars, mcfile, p0_scale=None, nwalkers=64, niter=40000):
        assert self.rv_t is not None, "no rv data found"
        assert self.jd is not None, "no lc data found"
        self.ndim = len(pars)
        self.nwalkers=nwalkers
        self.niter=niter
        if p0_scale is None:
            p0_scale = np.ones(self.ndim)*1e-4
        self.p0 = [pars + p0_scale*pars*np.random.randn(self.ndim) for ii in range(self.nwalkers)]
        if os.path.isfile(mcfile):
            print("File {0} already exists... do you want to clobber?".format(mcfile))
            return
        outf=open(mcfile, "w")
        outf.close()
        start_time = time.time()
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, 
                                             self.lnprob_lcrv, threads=4,
                                             args=[np.unique(self.quarter)])
        print("Running {0}k MCMC chain".format(self.niter/1000))
        for res in self.sampler.sample(self.p0, iterations=self.niter, storechain=False):
            if self.sampler.iterations % 10 == 0:
                position = res[0]
                outf = open(mcfile, "a")
                for k in range(position.shape[0]):
                    outf.write("{0} {1} {2} {3} {4}\n".format(self.sampler.iterations,
                               k, self.sampler.acceptance_fraction[k], res[1][k],
                                " ".join([str(ii) for ii in position[k]])))
                outf.close()
            if self.sampler.iterations % 10000 == 0:
                print("Time elapsed since niter={0}:{1}".format(self.sampler.iterations, 
                      time.time-start_time))
        print("Total time elapsed for MCMC run:{0}".format(time.time()-start_time))
        print("Total acceptance fraction:{0}".format(np.mean(self.sampler.acceptance_fraction)))
        try:
            print("Total autocorr time:{0}".format(np.mean(self.sampler.acor)))
        except:
            print("Could not compute autocorr time...")
        return

    def run_emcee2(self, pars, mcfile, p0_scale=None, nwalkers=64, niter=40000):
        assert self.rv_t is not None, "no rv data found"
        assert self.jd is not None, "no lc data found"
        self.ndim = len(pars)
        self.nwalkers=nwalkers
        self.niter=niter
        if p0_scale is None:
            p0_scale = np.ones(self.ndim)*1e-4
        self.p0 = [pars + p0_scale*pars*np.random.randn(self.ndim) for ii in range(self.nwalkers)]
        start_time = time.time()
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, 
                                             self.lnprob_lcrv, threads=4,
                                             args=[np.unique(self.quarter)])
        self.sampler.run_mcmc(self.p0, self.niter)
        print("Running {0}k MCMC chain".format(self.niter/1000))
        
        print("Total time elapsed for MCMC run:{0}".format(time.time()-start_time))
        print("Total acceptance fraction:{0}".format(np.mean(self.sampler.acceptance_fraction)))
        try:
            print("Total autocorr time:{0}".format(np.mean(self.sampler.acor)))
        except:
            print("Could not compute autocorr time...")
        return
        
    def plot_sed(self, mlpars, prefix, suffix='sed', lc_constraints=None, 
                       savefig=True):
        isoerr = mlpars[-1]
        if isoerr < 0:
            isoerr = np.exp(isoerr)
        magsmod = self.isofit(mlpars)
        #magsmod, r, T, logg = isofit_single(mlpars[:5])
        if np.any(np.isinf(magsmod)):
            print "Input isopars give -inf magsmod: ", mlpars, magsmod
            return magsmod
        plt.figure()
        plt.subplot(311)
        plt.errorbar(self.maglams, self.magsobs, 
                     np.sqrt(self.emagsobs**2 + isoerr**2), fmt='k.')
        plt.plot(self.maglams, magsmod, 'r.')
        plt.ylabel('magnitude')
    
        plt.subplot(312)
        plt.errorbar(self.maglams, self.magsobs-magsmod, 
                     np.sqrt(self.emagsobs**2 + isoerr**2), fmt='k.')
        plt.xlabel('wavelength (angstrom)')
        plt.ylabel('data-model')
        plt.suptitle('KIC '+str(self.kic) +' SED (only)')
    
        plt.subplot(313)
        if lc_constraints is not None:
            plt.plot(lc_constraints, 'ko')
        lc_inputs = np.array([(self.r1 +self.r2)/(mlpars[0]+mlpars[1])**(1./3.), 
                              self.r2/self.r1, self.frat])
        plt.plot(lc_inputs, 'ro')
        plt.xlim((-1, 3))
        if savefig:
            plt.savefig(prefix + suffix+'.png')
        return magsmod

        
    def plot_lc(self, lcpars, prefix, suffix='lc', savefig=True, 
                      polyorder=2, ooe=True, clip_tol=1.5):
        self.updatephase(lcpars[4], lcpars[3], clip_tol=clip_tol)
    
        lcmod, lcpol = self.lcfit(lcpars[:13], self.jd[self.clip], self.quarter[self.clip],
                                    self.flux[self.clip], self.dflux[self.clip],
                                    self.crowd[self.clip], polyorder=polyorder, ooe=ooe)
    
        lcres = self.flux[self.clip] - lcmod*lcpol
    
        lcerr=0
        if len(lcpars)==14:
            lcerr=lcpars[-1]
            if lcerr<0:
                lcerr=np.exp(lcerr)
    
        pe = (self.phase[self.clip] >= -1.2*self.pwidth) * (self.phase[self.clip] <= 1.2*self.pwidth)
        se = (self.phase[self.clip] >= -1.2*self.swidth+self.sep) * (self.phase[self.clip] <= 1.2*self.swidth+self.sep)
    
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
    
        ax1.errorbar(self.phase[self.clip][pe], self.flux[self.clip][pe]/lcpol[pe],
                     self.dflux[self.clip][pe], fmt='k.', ecolor='0.9')
        ax1.plot(self.phase[self.clip][pe], lcmod[pe], 'r.')
        ax2.errorbar(self.phase[self.clip][se], self.flux[self.clip][se]/lcpol[se],
                     self.dflux[self.clip][se], fmt='k.', ecolor='0.9')
        ax2.plot(self.phase[self.clip][se], lcmod[se], 'r.')
        ax3.errorbar(self.phase[self.clip][pe], lcres[pe],
                     self.dflux[self.clip][pe], fmt='k.', ecolor='0.9')
        ax4.errorbar(self.phase[self.clip][se], lcres[se],
                     self.dflux[self.clip][se], fmt='k.', ecolor='0.9')
    
    
    
        ax1.set_xlim((-1.2*self.pwidth, 1.2*self.pwidth))
    
    
        ax2.set_xlim((-1.2*self.swidth+self.sep, 1.2*self.swidth+self.sep))
    
    
        fig.suptitle('KIC '+str(self.kic)+' LC (only)')
        plt.subplots_adjust(hspace=0)
        if savefig:
            plt.savefig(prefix + suffix+'.png')
        return True

    def plot_sedlc(self, allpars, prefix, suffix='sedlc', savefig=True, 
                         polyorder=2, ooe=True):
        if allpars[-1] < 0:
            allpars[-1] = np.exp(allpars[-1])
        if allpars[4] < 10:
            allpars[4] = np.exp(allpars[4])
        if allpars[16] < 0:
            allpars[16] = np.exp(allpars[16])
        residuals = self.lnlike(allpars, lc_constraints=None, qua=np.unique(self.quarter),
                                  polyorder=polyorder, residual=True)
        lcpars = self.getpars(partype='lc')[:13]
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(self.maglams, self.magsobs, self.emagsobs, fmt='k.')
        ax.plot(self.maglams, residuals[:len(self.maglams)] * \
                 np.sqrt(self.emagsobs**2 + self.pars['isoerr']**2) + self.magsobs, 'r.')
        for ii in range(len(self.maglams)):
            ax.text(self.maglams[ii], self.magsobs[ii], self.ipname[ii].replace('mag', ''))
        ax.set_ylabel('Magnitude')
    
        divider = make_axes_locatable(ax)
        ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
        ax2.errorbar(self.maglams, residuals[:len(self.maglams)] * \
                 np.sqrt(self.emagsobs**2 + self.pars['isoerr']**2),
                     np.sqrt(self.emagsobs**2 + self.pars['isoerr']**2), fmt='k.')
        ax2.set_xlabel('Wavelength (Angstrom)')
        ax2.set_ylabel('Data - Model')
        ax2.set_ylim((-0.3, 0.3))
        plt.setp(ax.get_xticklabels(), visible=False)
    
        plt.suptitle('KIC '+str(self.kic)+' SED (simultaneous)')
        if savefig:
            plt.savefig(prefix+suffix+'_SED.png')
    
        self.updatephase(lcpars[4], lcpars[3], clip_tol=self.clip_tol)
        lcmod, lcpol = self.lcfit(lcpars, self.jd[self.clip], self.quarter[self.clip],
        				self.flux[self.clip], self.dflux[self.clip],
        				self.crowd[self.clip], polyorder=polyorder, ooe=ooe)
    
        lcres = self.flux[self.clip] - lcmod*lcpol
    
        fig = plt.figure(figsize=(16, 16))
        ax = fig.add_subplot(121)
        ax.errorbar(self.phase[self.clip], self.flux[self.clip]/lcpol,
                     self.dflux[self.clip], fmt='k.', ecolor='0.9')
        ax.plot(self.phase[self.clip], lcmod, 'r.')
        ax.set_xlim((-1.2*self.pwidth, 1.2*self.pwidth))
        ax.set_ylabel('Kepler Flux')
    
        divider = make_axes_locatable(ax)
        axb = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
        axb.errorbar(self.phase[self.clip], lcres,
                     np.sqrt(self.dflux[self.clip]**2 + self.pars['lcerr']**2), fmt='k.', ecolor='0.9')
    
        axb.set_xlim((-1.2*self.pwidth, 1.2*self.pwidth))
        axb.set_ylabel('Data - Model')
        axb.set_xlabel('Phase (Primary Eclipse)')
        #axb.set_yticklabels(axb.yaxis.get_majorticklabels()[1:])
    
        ax2 = fig.add_subplot(122)
        ax2.errorbar(self.phase[self.clip], self.flux[self.clip]/lcpol,
                     self.dflux[self.clip], fmt='k.', ecolor='0.9')
        ax2.plot(self.phase[self.clip], lcmod, 'r.')
        ax2.set_xlim((-1.2*self.swidth+self.sep, 1.2*self.swidth+self.sep))
    
        divider2 = make_axes_locatable(ax2)
        ax2b = divider2.append_axes("bottom", size=2.0, pad=0, sharex=ax2)
        ax2b.errorbar(self.phase[self.clip], lcres,
                     np.sqrt(self.dflux[self.clip]**2 + self.pars['lcerr']**2), fmt='k.', ecolor='0.9')
    
        ax2b.set_xlim((-1.2*self.swidth+self.sep, 1.2*self.swidth+self.sep))
        ax2b.set_xlabel('Phase (Secondary Eclipse)')
    
        #ax2b.set_yticklabels(ax2b.yaxis.get_majorticklabels()[1:])
    
    
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.suptitle('KIC '+str(self.kic)+' LC (simultaneous)')
        if savefig:
            plt.savefig(prefix+suffix+'_LC.png')
        return True

    def plot_lcrv(self, allpars, prefix, suffix='lcrv', savefig=True, polyorder=2, ooe=True):
        residuals = self.lnlike_lcrv(allpars, qua=np.unique(self.quarter), polyorder=polyorder,
                                residual=True)
        lcpars = self.getpars(partype='lc')[:13]
        lcmod, lcpol = self.lcfit(lcpars, self.jd[self.clip], self.quarter[self.clip],
        				self.flux[self.clip], self.dflux[self.clip],
        				self.crowd[self.clip], polyorder=2, ooe=ooe)
    
        phase = ((self.jd[self.clip]-lcpars[4]) % lcpars[3])/lcpars[3]
        phase[phase<-np.clip(self.pwidth*3., 0., 0.2)]+=1.
        phase[phase>np.clip(self.sep+self.swidth*3., self.sep, 1.0)]-=1.
    
        lcres = self.flux[self.clip] - lcmod*lcpol
    
        fig = plt.figure(figsize=(16, 16))
        ax = fig.add_subplot(121)
        ax.plot(self.phase[self.clip], self.flux[self.clip]/lcpol, 'g.', alpha=0.4)
        ax.errorbar(phase, self.flux[self.clip]/lcpol,
                     self.dflux[self.clip], fmt='k.', ecolor='gray')
        ax.plot(phase, lcmod, 'r.')
        ax.set_xlim((-1.2*self.pwidth, 1.2*self.pwidth))
        ax.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))
        ax.set_ylabel('Kepler Flux')
    
        divider = make_axes_locatable(ax)
        axb = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
        axb.plot(self.phase[self.clip], lcres, 'g.', alpha=0.4)
        axb.errorbar(phase, lcres,
                     np.sqrt(self.dflux[self.clip]**2 + self.pars['lcerr']**2), fmt='k.', ecolor='gray')
    
        axb.set_xlim((-1.2*self.pwidth, 1.2*self.pwidth))
        axb.set_ylim((np.min(lcres), np.max(lcres)))
        axb.set_ylabel('Data - Model')
        axb.set_xlabel('Phase (Primary Eclipse)')
        #axb.set_yticklabels(axb.yaxis.get_majorticklabels()[1:])
    
        ax2 = fig.add_subplot(122)
        ax2.plot(self.phase[self.clip], self.flux[self.clip]/lcpol, 'g.', alpha=0.4)
        ax2.errorbar(phase, self.flux[self.clip]/lcpol,
                     self.dflux[self.clip], fmt='k.', ecolor='gray')
        ax2.plot(phase, lcmod, 'r.')
        ax2.set_xlim((-1.2*self.swidth+self.sep, 1.2*self.swidth+self.sep))
        ax2.set_ylim((np.min(lcmod)*0.98, np.max(lcmod)*1.02))
    
        divider2 = make_axes_locatable(ax2)
        ax2b = divider2.append_axes("bottom", size=2.0, pad=0, sharex=ax2)
        ax2b.plot(self.phase[self.clip], lcres, 'g.', alpha=0.4)
        ax2b.errorbar(phase, lcres,
                     np.sqrt(self.dflux[self.clip]**2 + self.pars['lcerr']**2), fmt='k.', ecolor='gray')
    
        ax2b.set_xlim((-1.2*self.swidth+self.sep, 1.2*self.swidth+self.sep))
        ax2b.set_ylim((np.min(lcres), np.max(lcres)))
        ax2b.set_xlabel('Phase (Secondary Eclipse)')
    
        #ax2b.set_yticklabels(ax2b.yaxis.get_majorticklabels()[1:])
    
    
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.suptitle('KIC '+str(self.kic)+' LC (simultaneous)')
        if savefig:
            plt.savefig(prefix+suffix+'_LC.png')
    
        # rvpars = self.getpars(partype='rv')
        rvpars = np.array([self.pars['msum'], self.pars['mrat'], self.pars['period'], self.pars['tpe'],
                           self.pars['esinw'], self.pars['ecosw'], self.pars['inc'], self.pars['k0'], self.pars['rverr']])
    
        rv_fit = self.rvfit(rvpars, self.rv_t)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rvphase = (self.rv_t - self.pars['tpe'])%self.pars['period']/self.pars['period']
        ax.errorbar(rvphase[~self.bad1], self.rv1_obs[~self.bad1], self.rv1_err_obs[~self.bad1], fmt='b*')
        ax.errorbar(rvphase[~self.bad2], self.rv2_obs[~self.bad2], self.rv2_err_obs[~self.bad2], fmt='r*')
        rvt = np.linspace(0, 1, 100)*self.pars['period']+self.pars['tpe']
        rvmod = self.rvfit(rvpars, rvt)
        ax.plot(np.linspace(0, 1, 100), rvmod[0], 'b-')
        ax.plot(np.linspace(0, 1, 100), rvmod[1], 'r-')
        ax.set_ylabel('RV [m/s]')
        #ax.set_yticklabels(ax.yaxis.get_majorticklabels()[:-1])
    
        divider = make_axes_locatable(ax)
        ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
        ax2.errorbar(rvphase[~self.bad1], (self.rv1_obs-rv_fit[0])[~self.bad1], np.sqrt(self.rv1_err_obs**2+rvpars[-1]**2)[~self.bad1], fmt='b.')
        ax2.errorbar(rvphase[~self.bad2], (self.rv2_obs-rv_fit[1])[~self.bad2], np.sqrt(self.rv2_err_obs**2+rvpars[-1]**2)[~self.bad2], fmt='r.')
    
        ax2.set_xlabel('Phase')
        ax2.set_ylabel('Data - Model')
        #ax2.set_yticklabels(ax2.yaxis.get_majorticklabels()[:-1])
        plt.setp(ax.get_xticklabels(), visible=False)
        #plt.setp(ax2.get_yticklabels()[-1], visible=False)
    
        plt.suptitle('KIC '+str(self.kic)+' RV (simultaneous)')
        if savefig:
            plt.savefig(prefix+suffix+'_RV.png')
    
        return True


    def plot_rv(self, rvpars, prefix, suffix='rv', savefig=True):  
        rv_fit = self.rvfit(rvpars, self.rv_t)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        rvphase = (self.rv_t - self.pars['tpe'])%self.pars['period']/self.pars['period']
        ax.errorbar(rvphase[~self.bad1], self.rv1_obs[~self.bad1], self.rv1_err_obs[~self.bad1], fmt='b*')
        ax.errorbar(rvphase[~self.bad2], self.rv2_obs[~self.bad2], self.rv2_err_obs[~self.bad2], fmt='r*')
        rvt = np.linspace(0, 1, 100)*self.pars['period']+self.pars['tpe']
        rvmod = self.rvfit(rvpars, rvt)
        ax.plot(np.linspace(0, 1, 100), rvmod[0], 'b-')
        ax.plot(np.linspace(0, 1, 100), rvmod[1], 'r-')
        ax.set_ylabel('RV [m/s]')
        #ax.set_yticklabels(ax.yaxis.get_majorticklabels()[:-1])
    
        divider = make_axes_locatable(ax)
        ax2 = divider.append_axes("bottom", size=2.0, pad=0, sharex=ax)
        ax2.errorbar(rvphase[~self.bad1], (self.rv1_obs-rv_fit[0])[~self.bad1], np.sqrt(self.rv1_err_obs**2+rvpars[-1]**2)[~self.bad1], fmt='b.')
        ax2.errorbar(rvphase[~self.bad2], (self.rv2_obs-rv_fit[1])[~self.bad2], np.sqrt(self.rv2_err_obs**2+rvpars[-1]**2)[~self.bad2], fmt='r.')
    
        ax2.set_xlabel('Phase')
        ax2.set_ylabel('Data - Model')
        #ax2.set_yticklabels(ax2.yaxis.get_majorticklabels()[:-1])
        plt.setp(ax.get_xticklabels(), visible=False)
        #plt.setp(ax2.get_yticklabels()[-1], visible=False)
    
        plt.suptitle('KIC '+str(self.kic)+' RV (simultaneous)')
        if savefig:
            plt.savefig(prefix+suffix+'_RV.png')
    
        return True
        

    @staticmethod
    def get_pars2vals(fisopars, partype='lc', crow=[]):
        try:
            parnames = parnames_dict[partype]
            if len(crow)>0:
                parnames += ['cr' + str(ii) for ii in crow]
        except:
            print("You entered: {0}. Partype options are 'lc', 'sed', 'rv', 'lcsed', 'lcrv'. Try again.".format(partype))
            return
    
        parvals = np.zeros(len(parnames))
        novalue = (len(fisopars) == len(parnames))
        #print fisopars, type(fisopars), parnames
        if isinstance(fisopars, lmfit.parameter.Parameters):
            for j in range(len(parnames)):
                parvals[j] = fisopars[parnames[j]].value
        elif isinstance(fisopars, dict):
            for j in range(len(parnames)):
                parvals[j] = fisopars[parnames[j]]
        else:
            for j in range(len(parnames)):
                parvals[j] = fisopars[j]*novalue
        return parvals
        
    @staticmethod
    def broadcast_crowd(qrt, cro):
        """Takes in crowding parameters for each quarter and broadcasts to all cadences

        Parameters
        ----------
        qrt: array of length equal to time array
            array that specifies the quarter # of each cadence
        cro: array of length N quarters
            array that specifies the crowing values for each unique quarter

        Returns
        -------
        cro_casted: array of length equal to time array
            output array with crowding values broadcasted to each cadence
        """
        cro = cro.ravel()
        good = np.where(np.diff(qrt)>0)[0] + 1
        good = np.insert(good, 0, 0)
        good = np.insert(good, len(good), len(qrt))
        cro_casted = qrt*1.0
        for ii in range(len(good)-1):
            cro_casted[good[ii]:good[ii+1]] = cro[ii]
        return cro_casted

        
    @staticmethod
    def get_P(a, msum):
        return np.sqrt(a**3 / msum) * d2y

    @staticmethod
    def get_a(period, msum):
        """Computes semi-major axis given period and total mass

        Parameters
        ----------
        period: Period of system in days
        msum: sum of masses in Msun

        Returns
        -------
        a: semi-major axis in AU
        """
        return ((period/d2y)**2 * (msum))**(1./3.)

    @staticmethod
    def sumrat_to_12(xsum, xrat):
        """Transforms sum and ratios of 2 quantities into individual components

        Parameters
        ----------
        xsum: x1 + x2 value
        xrat: x2/x1 value

        Returns
        -------
        x1, x2
        """
        x1 = xsum / (1+xrat)
        x2 = xsum / (1 + 1./xrat)
        return x1, x2

    @staticmethod
    def sudarsky(theta, e, period):
        """Computes time as a function of mean anomaly+omega given period, eccentricity

        Parameters
        ----------
        theta : float array or scalar
                value of omega + maf (arg. periapse + mean anomaly) in radians
        e : scalar
            eccentricity of binary
        period : scalar
            period of binary

        Returns
        -------
        t : float array or scalar
            time given theta, e, P
        """
        tt = (-np.sqrt(1. - e ** 2) * period / TWOPI) * \
             (e * np.sin(theta) / (1. + e * np.cos(theta)) - 2. * (1. - e ** 2) ** (-0.5) *
              np.arctan(np.sqrt(1. - e ** 2) * np.tan((theta) / 2.) / (1. + e)))
        return tt

    @staticmethod
    def get_inc(b, r1, a):
        """Transforms impact parameter into inclination

        Parameters
        ----------
        b: impact parameter
        r1: radius of primary star (in Rsun)
        a: semi-major axis in AU

        Returns
        -------
        inc: inclination of system"""
        return np.arccos(b*r1 / (a/r2au))

    @staticmethod
    def find_closest(base, target):
        """Returns indices of closest match between a base and target arrays
        
        Parameters
        ----------
        base : array
            base list of items to compare
        target : array
            list of target items to which to compare to the base array
            
        Returns
        -------
        indx : int
            index of matched item
        
        Examples
        --------
        >>> a = np.arange(10)
        >>> b = 4.3
        >>> print 
        """
        
        #A must be sorted
        idx = base.searchsorted(target)
        idx = np.clip(idx, 1, len(base)-1)
        left = base[idx-1]
        right = base[idx]
        idx -= target - left < right - target
        return idx


    @staticmethod
    def feh2z(feh):
        """Converts [Fe/H] value to Z metallicity fraction (zsun * 10**feh)

        Parameters
        ----------
        feh : array or scalar
            [Fe/H] value of object

        Returns
        -------
        z : array or scalar
            Z metallicity"""
        return 0.01524 * 10**feh

    @staticmethod
    def z2feh(z):
        """Converts Z metallicity value to [Fe/H] (log10(z/zsun))

        Parameters
        ----------
        z : array or scalar

        Returns
        -------
        feh : array or scalar
        """
        return np.log10(z/0.01524)
        
    @staticmethod
    def kipping_q2u(q1, q2):
        c1 = 2.*np.sqrt(q1)*q2
        c2 = np.sqrt(q1)*(1.-2.*q2)
        return c1, c2
    
    @staticmethod
    def kipping_u2q(u1, u2):
        q1 = (u1+u2)**2
        q2 = 0.5 * u1 / (u1+u2)
        return q1, q2