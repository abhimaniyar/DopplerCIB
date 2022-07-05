from headers_constants import *


class ProfHODMore15:
    # def __init__(self, U, MassFunc):
    def __init__(self, data_var_iv, cosmo_var_iv, gal_exp):
        self.dv = data_var_iv
        self.uni = cosmo_var_iv
        self.mh = self.uni.mass
        self.z = self.uni.z  # self.dv.z
        self.ell = self.dv.ell
        self.hmfmz = self.uni.hmf  # self.dv.hmf
        # self.biasmz = self.uni.bias_m_z
        self.biasmz = self.uni.interp_bias(self.z)
        self.unfw = self.uni.interp_nfw(self.ell, self.z)
        # self.dn_dz = self.dv.dn_dz()
        # self.Pk = self.uni.pkinterpz(self.z)
        self.gal_exp = gal_exp
        """HOD for CMASS
        from More Miyatake +15
        """
        # HOD params, More+15, average of the 3 columns in Table 1
        self.alphaInc = 0.51
        self.log10mInc = 13.84
        # self.log10mMin = 12.7  # 12.00  # 13.42
        """
        log10mMin determines somehow the minimum mass detectable by the survey.
        Of course, DESI survyes are going to detect galaxies much better than
        CMASS. So we need to go to lower lowest detectable mass than CMASS for
        these surveys. These values for DESI LRG and DESI ELG are determined
        such that total number of galaxies approximately detected (N_tot function below) by
        LRG in 0.4 < z < 1 are 10 times more than CMASS, ELG in 0.6 < z < 1.6
        are 30 times more than CMASS. For ExtDESI_ELG, we keep the number of
        galaxies detected as same as DESI ELG i.e. 30 times more than CMASS,
        but spread on abroader redshift range i.e. 0.6 < z < 3.2.
        """
        if self.gal_exp == 'CMASS':
            self.log10mMin = 13.42
        elif self.gal_exp == 'DESI_LRG':
            self.log10mMin = 12.5
        elif self.gal_exp == 'DESI_ELG':
            self.log10mMin = 12.1
        elif self.gal_exp == 'ExtDESI_ELG':
            self.log10mMin = 12.7
        self.sLog10m = np.sqrt(0.49)
        self.kappa = 1.10
        self.mMinHod = 10.**(13.42)
        self.m1 = 10.**(14.43)
        self.alpha_ = 1.09
        self.pOff = 0.357
        self.rOff = 2.3

        self.mMin = 0.
        self.mMax = np.inf

        # self.galexp = gal_exp
        # super(ProfHODMore15, self).__init__(U, MassFunc)

    def __str__(self):
        return "hodmore15"

# ############################################################################
    def dn_dz(self):
        if self.gal_exp == 'CMASS':
            # """
            # this is CMASS data
            data = np.loadtxt("data_files/dn_dz_cmass.txt")
            Z = data[:, 0]
            Dndz = data[:, 1]
        elif self.gal_exp == 'DESI_ELG' or self.gal_exp == 'DESI_LRG':
            # """
            # DESI_ELG and DESI_LRG data. Please note however, we are using the same HOD as that
            # of CMASS derived in More+15. Not sure how good this approx is
            # DESi data in a redshift bin is number of galaxies per square
            # degrees.
            # 1. This can be converted into ngal i.e. angular number density
            # which is number of galaxies per unit steradian. To go from per
            # square degree to per steradian => divide by (pi/180)**2
            # 2. dN/dz which is number of galaxies per redshift per steradian. So
            # dN/dz = ngal/\Delta z where \Delta z is redshift bin width i.e. 0.1
            # in this case. and ngal is calculated as in number 1 step.
            # data = np.loadtxt("data/dndz_DESI_ELG.txt")
            # data = np.loadtxt("data/dndz_DESI_LRG.txt")
            data = np.loadtxt('data_files/dndz_'+self.gal_exp+'.txt')
            Z1, Z2 = data[:, 0], data[:, 1]
            Z = (Z1+Z2)/2.
            deltaz = Z2 - Z1
            # converting from per squre degrees to per steradian
            ngal = data[:, 2]/(np.pi/180.)**2
            Dndz = ngal/deltaz
        else:
            data = np.loadtxt('data_files/dndz_DESI_ELG.txt')
            Z1, Z2 = data[:, 0], data[:, 1]
            Z = Z1+Z2
            deltaz = 2*(Z2 - Z1)
            # converting from per squre degrees to per steradian
            ngal = data[:, 2]/(np.pi/180.)**2
            # ngal *= 10
            Dndz = ngal/deltaz
            # """
        f = UnivariateSpline(Z, Dndz, k=1, s=0, ext=1)
        dndz = lambda z: f(z) * (z >= np.min(Z)) * (z <= np.max(Z))
        return dndz

    def fInc(self, m):
        """More+15 assume that for a given halo mass,
        a fixed fraction of CMASS galaxies are seen
        how physical is this?
        """
        # m = self.mh
        result = np.min([1., 1.+self.alphaInc*(np.log10(m) - self.log10mInc)])
        result = np.max([0., result])
        return result

    def Ncen(self, m):
        """number of central galaxies per halo, between 0 and 1
        """
        # m = self.mh
        result = np.log10(m) - self.log10mMin
        result /= self.sLog10m
        result = 0.5*(1.+special.erf(result))
        result *= self.fInc(m)
        return result

    def Nsat(self, m):
        """number of satellite galaxies per halo
        """
        # m = self.mh
        result = (m - self.kappa * self.mMinHod)
        if result > 0.:
            result /= self.m1
            result **= self.alpha_
            result *= self.Ncen(m)
        else:
            result = 0.
        return result

    def nbargal(self):
        m = self.mh
        hmf = self.hmfmz
        N_g = np.zeros(len(m))
        for mm in range(len(m)):
            N_g[mm] = self.Ncen(m[mm])+self.Nsat(m[mm])
        dm = np.log10(m[1] / m[0])
        integral = hmf*N_g[:, None]
        return intg.simps(integral, dx=dm, axis=0, even='avg')

    def N_tot(self):
        z = self.z
        # nbarg = self.nbargal()
        fnbar = interp1d(z, self.nbargal(), kind='linear',
                         bounds_error=False, fill_value=0.)
        if self.gal_exp == 'CMASS':
            zmin, zmax = 0.4, 0.7
        elif self.gal_exp == 'DESI_LRG':
            zmin, zmax = 0.4, 1.0
        elif self.gal_exp == 'DESI_ELG':
            zmin, zmax = 0.6, 1.6
        elif self.gal_exp == 'ExtDESI_ELG':
            zmin, zmax = 0.6, 3.2
        zint = np.linspace(zmin, zmax, 30)
        chi = self.uni.chi(zint)
        dchidz = self.uni.dchi_dz(zint)
        integral = fnbar(zint)*4*np.pi*chi**2*dchidz
        res = intg.simps(integral, x=zint, axis=0, even='avg')
        return res

    def Nbargal(self):
        Dndz = self.dn_dz()
        res = lambda z: intg.simps(Dndz(z), x=z, even='avg')
        return res

    def window_gal(self):
        dndz = self.dn_dz()
        nbarg = self.Nbargal()
        # this should be divided by dchi_dz() but as it is defined later, we
        # do the division later while calculating the total window func
        res = lambda z: dndz(z)/nbarg(z)
        return res

    def p1h_gal(self, ucen, unfw):
        m = self.mh
        z = self.z
        hmf = self.hmfmz
        # ucen = self.uni.self.nfw_u  # unfw
        # unfw = self.uni.self.nfw_u  # unfw
        num = np.zeros((len(unfw[0, :, 0]), len(m), len(z)))
        for mi in range(len(m)):
            num[:, mi, :] = 2*self.Ncen(m[mi])*ucen[mi, :, :]*self.Nsat(m[mi])*unfw[mi, :, :]
            num[:, mi, :] += (self.Nsat(m[mi])*unfw[mi, :, :])**2
        denom = self.nbargal()**2
        integral = hmf*num/denom
        # print (integral.shape)
        # print (np.shape(m))
        # dm = np.log10(m)
        dlog10m = np.log10(m[1] / m[0])
        res = intg.simps(integral, dx=dlog10m, axis=1, even='avg')
        return res

    def p2h_gal(self, ucen, unfw, power):
        m = self.mh
        z = self.z
        hmf = self.hmfmz
        # ucen = self.uni.self.nfw_u  # unfw
        # unfw = self.uni.self.nfw_u  # unfw
        biasmz = self.biasmz
        # power = self.uni.pkinterpz(self.z)
        
        num = np.zeros((len(unfw[0, :, 0]), len(m), len(z)))
        for mi in range(len(m)):
            num[:, mi, :] = self.Ncen(m[mi])*ucen[mi, :, :]+self.Nsat(m[mi])*unfw[mi, :, :]
        num *= hmf*biasmz
        denom = self.nbargal()
        integral = num/denom
        # dm = np.log10(m)
        dlog10m = np.log10(m[1] / m[0])
        res = (intg.simps(integral, dx=dlog10m, axis=1, even='avg'))**2
        res *= power
        return res

    def cl1h_gal(self):
        z = self.z
        ucen = self.unfw
        unfw = self.unfw
        dchidz = self.uni.dchi_dz(self.z)
        chi = self.uni.chi(self.z)

        window = self.window_gal()
        wind_gal = window(z)/dchidz
        geo = dchidz*wind_gal**2/chi**2
        oneh = self.p1h_gal(ucen, unfw)
        # print ("one halo gal Pk %s" % (oneh[-1, -10:]))
        integral = geo*oneh
        res = intg.simps(integral, x=z, axis=-1, even='avg')
        # print ("one halo gal cl %s" % (res[-10:]))
        return res

    def cl2h_gal(self):
        z = self.z
        ucen = self.unfw
        unfw = self.unfw
        power = self.uni.Pk_array(self.ell, self.z)
        dchidz = self.uni.dchi_dz(self.z)
        chi = self.uni.chi(self.z)

        window = self.window_gal()
        wind_gal = window(z)/dchidz
        geo = dchidz*wind_gal**2/chi**2
        twoh = self.p2h_gal(ucen, unfw, power)
        integral = geo*twoh
        res = intg.simps(integral, x=z, axis=-1, even='avg')
        return res

    def clshot_gal(self):
        if self.gal_exp == 'CMASS':
            # """
            # this is CMASS data
            data = np.loadtxt("data/dn_dz_cmass.txt")
            Z = data[:, 0]
            Dndz = data[:, 1]
            ngaltot = intg.simps(Dndz, x=Z, even='avg')
        elif self.gal_exp == 'DESI_ELG' or self.gal_exp == 'DESI_LRG':
            data = np.loadtxt('data/dndz_'+self.gal_exp+'.txt')
            # converting from per squre degrees to per steradian
            ngal = data[:, 2]/(np.pi/180.)**2
            ngaltot = np.sum(ngal)
        else:
            data = np.loadtxt('data/dndz_DESI_ELG.txt')
            # converting from per squre degrees to per steradian
            ngal = data[:, 2]/(np.pi/180.)**2
            ngaltot = np.sum(ngal)
            # ngaltot *= 10
            """
            # data = np.loadtxt("data/dndz_DESI_ELG.txt")
            data = np.loadtxt("data/dndz_DESI_LRG.txt")
            # add all the ngal contributions coming from different redshifts
            # and then take its inverse to get the shot noise in angular power
            # spectrum i.e. C_l^{SN} = 1./ngal where ngal is the angular number
            # density i.e. number of galaxies per steradian.
            # shot = self.Nbargal()
            # """
        shot = 1./ngaltot
        return shot

    def cl_galtot(self):
        oneh = self.cl1h_gal()
        twoh = self.cl2h_gal()
        shot = np.repeat(self.clshot_gal(), len(self.ell))
        tot = oneh+twoh+shot
        return tot

    def cl1h_vel(self):
        z = self.z
        ucen = self.unfw
        unfw = self.unfw
        dchidz = self.uni.dchi_dz(self.z)
        chi = self.uni.chi(self.z)

        window = self.window_gal()
        wind_gal = window(z)/dchidz
        geo = dchidz*wind_gal**2/chi**2
        oneh = self.p1h_gal(ucen, unfw)
        oneh *= self.uni.beta2(self.z)
        # print ("one halo vell Pk %s" % (oneh[-1, -10:]))
        integral = geo*oneh
        res = intg.simps(integral, x=z, axis=-1, even='avg')
        # print ("one halo vel cl %s" % (res[-10:]))
        return res

    def cl2h_vel(self):
        z = self.z
        ucen = self.unfw
        unfw = self.unfw
        power = self.uni.Pk_array(self.ell, self.z)
        dchidz = self.uni.dchi_dz(self.z)
        chi = self.uni.chi(self.z)

        window = self.window_gal()
        wind_gal = window(z)/dchidz
        geo = dchidz*wind_gal**2/chi**2
        twoh = self.p2h_gal(ucen, unfw, power)
        twoh *= self.uni.beta2(self.z)
        integral = geo*twoh
        res = intg.simps(integral, x=z, axis=-1, even='avg')
        return res

    def clshot_vel(self):
        """
        SHOT NOISE FOR VELOCITY POWER SPECTRUM? DO WE WEIGHT BETA2 AT EVERY
        REDSHIFT BIN WITH NUMBER DENSITY OF THE GALAXIES?
        """
        """
        For now, I am going to weigh the shot noise contribution coming from
        every redshift by <beta^2>. And then get the shot noise in angular
        space.
        """
        
        if self.gal_exp == 'CMASS':
            # """
            # this is CMASS data
            data = np.loadtxt("data/dn_dz_cmass.txt")
            Z = data[:, 0]
            Dndz = data[:, 1]
            beta2 = self.uni.beta2(Z)
            res = Dndz*beta2
            ngaltot = intg.simps(res, x=Z, even='avg')
        else:
            data = np.loadtxt('data/dndz_'+self.gal_exp+'.txt')
            # converting from per squre degrees to per steradian
            ngal = data[:, 2]/(np.pi/180.)**2
            Z = (data[:, 0]+data[:, 1])/2.
            beta2 = self.uni.beta2(Z)
            res = ngal*beta2
            ngaltot = np.sum(res)
            """
            # data = np.loadtxt("data/dndz_DESI_ELG.txt")
            data = np.loadtxt("data/dndz_DESI_LRG.txt")
            # add all the ngal contributions coming from different redshifts
            # and then take its inverse to get the shot noise in angular power
            # spectrum i.e. C_l^{SN} = 1./ngal where ngal is the angular number
            # density i.e. number of galaxies per steradian.
            # shot = self.Nbargal()
            # """
        shot = 1./ngaltot
        # print (shot)
        return shot

    def cl_veltot(self):
        oneh = self.cl1h_vel()
        twoh = self.cl2h_vel()
        # shot = np.repeat(self.clshot_vel(), len(self.ell))
        tot = oneh+twoh  # +shot
        return tot
