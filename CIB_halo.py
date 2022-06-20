from headers_constants import *


class Cib_halo:

    def __init__(self, data_var_iv, cosmo_var_iv):
        self.dv = data_var_iv
        self.uni = cosmo_var_iv
        self.z = self.uni.z  # self.dv.z
        self.z_c = self.dv.z_c
        self.mh = self.uni.mass
        # self.nu0min = self.dv.nu0min
        # self.nu0max = self.dv.nu0max
        # self.nu0 = np.linspace(self.nu0min, self.nu0max, 200)
        self.nu0 = self.dv.nu0
        # self.nucen = self.dv.nucen
        # self.snu_eff = self.dv.snu
        self.ell = self.dv.ell
        self.snu_unfilt = self.dv.unfiltered_snu(self.nu0, self.z)
        self.snu_filt = self.dv.snufilt(self.z)
        self.snu = self.snu_unfilt
        # print (self.snu_unfilt[4, 10:13])
        self.cosmo = cosmo
        # self.deltah = deltah
        self.Meffmax = self.dv.Meffmax
        self.etamax = self.dv.etamax
        self.sigmaMh = self.dv.sigmaMh
        self.tau = self.dv.tau
        self.hmfmz = self.uni.hmf  # self.dv.hmf
        # self.biasmz = self.uni.bias_m_z
        self.biasmz = self.uni.interp_bias(self.z)
        self.sig_z = np.array([max(self.z_c - r, 0.) for r in self.z])
        self.sigpow = self.sigmaMh - self.tau*self.sig_z
        self.unfw = self.uni.interp_nfw(self.ell, self.z)
        self.fc = self.dv.fc
        self.cc = self.dv.cc

    def sfr_mhdot(self, mhalo):
        """ SFR/Mhdot lognormal distribution wrt halomass """
        if hasattr(mhalo, "__len__"):
            a = np.zeros((len(mhalo), len(self.z)))
            for i in range(len(mhalo)):
                if mhalo[i] < self.Meffmax:
                    a[i, :] = self.etamax * np.exp(-(np.log(mhalo[i]) - np.log(self.Meffmax))**2 / (2 * self.sigmaMh**2))
                else:
                    a[i, :] = self.etamax * np.exp(-(np.log(mhalo[i]) - np.log(self.Meffmax))**2 / (2 * self.sigpow**2))
        else:
            if mhalo < self.Meffmax:
                a = self.etamax * np.exp(-(log(mhalo) - log(self.Meffmax))**2 / (2 * self.sigmaMh**2))
            else:
                a = self.etamax * np.exp(-(log(mhalo) - log(self.Meffmax))**2 / (2 * self.sigpow**2))
        return a

    def Mdot(self, mhalo):
        use_mean = True
        if use_mean:
            a = 46.1*(1 + 1.11*self.z) * \
                np.sqrt(self.cosmo.Om0 * (1 + self.z)**3 + self.cosmo.Ode0)
            b = (mhalo / 1.0e12)**1.1
            return np.outer(b, a)
        else:
            a = 25.3*(1 + 1.65*self.z) * \
                np.sqrt(self.cosmo.Om0*(1 + self.z)**3 + self.cosmo.Ode0)
            b = (mhalo / 1.0e12)**1.1
            return np.outer(b, a)

    def sfr(self, mhalo):
        sfrmhdot = self.sfr_mhdot(mhalo)
        mhdot = self.Mdot(mhalo)
        f_b = self.cosmo.Ob(self.z)/self.cosmo.Om(self.z)
        return mhdot * f_b * sfrmhdot

    def djc_dlnMh(self):
        fsub = 0.134
        """fraction of the mass of the halo that is in form of
        sub-halos. We have to take this into account while calculating the
        star formation rate of the central halos. It should be calulated by
        accounting for this fraction of the subhalo mass in the halo mass
        central halo mass in this case is (1-f_sub)*mh where mh is the total
        mass of the halo.
        for a given halo mass, f_sub is calculated by taking the first moment
        of the sub-halo mf and and integrating it over all the subhalo masses
        and dividing it by the total halo mass.
        """
        # a = np.zeros((len(self.snu_eff[:, 0]), len(self.mh), len(self.z)))
        snu = self.snu
        a = np.zeros((len(snu[:, 0]), len(self.mh), len(self.z)))
        rest = self.sfr(self.mh*(1-fsub))*(1 + self.z) *\
            self.cosmo.comoving_distance(self.z).value**2/KC
        # print (rest[50, 10:13])
        for f in range(len(snu[:, 0])):
            a[f, :, :] = rest*snu[f, :]
        return a

    def subhmf(self, mhalo, ms):
        # subhalo mass function from (https://arxiv.org/pdf/0909.1325.pdf)
        return 0.13*(ms/mhalo)**(-0.7)*np.exp(-9.9*(ms/mhalo)**2.5)*np.log(10)
    # np.log(10) added in the end as in the integration we are integrating with
    # respect to dlogm to the base 10.

    def msub(self, mhalo):
        """
        for a given halo mass mh, the subhalo masses would range from
        m_min to mh. For now, m_min has been taken as 10^5 solar masses
        """
        log10msub_min = 5
        if np.log10(mhalo) <= log10msub_min:
            raise ValueError, "halo mass %d should be greater than subhalo mass \
%d." % (np.log10(mhalo), log10msub_min)
        else:
            logmh = np.log10(mhalo)
            logmsub = np.arange(log10msub_min, logmh, 0.1)
            return 10**logmsub
        
    def djsub_dlnMh(self):
        """
        for subhalos, the SFR is calculated in two ways and the minimum of the
        two is assumed.
        """
        fsub = 0.134
        # a = np.zeros((len(self.snu_eff[:, 0]), len(self.mh), len(self.z)))
        snu = self.snu
        a = np.zeros((len(snu[:, 0]), len(self.mh), len(self.z)))
        # sfrmh = self.sfr(mh)
        for i in range(len(self.mh)):
            ms = self.msub(self.mh[i]*(1-fsub))
            dlnmsub = np.log10(ms[1] / ms[0])
            sfrI = self.sfr(ms)  # dim(len(ms), len(z))
            sfrII = self.sfr(self.mh[i]*(1-fsub))*ms[:, None]/(self.mh[i]*(1-fsub))
            # sfrII = sfrmh[i] * ms / mh[i]
            sfrsub = np.zeros((len(ms), len(self.z)))
            for j in range(len(ms)):
                sfrsub[j, :] = np.minimum(sfrI[j, :], sfrII[j, :])
            integral = self.subhmf(self.mh[i], ms)[:, None]*sfrsub / KC
            intgn = intg.simps(integral, dx=dlnmsub, axis=0)
            a[:, i, :] = snu*(1 + self.z)*intgn *\
                self.cosmo.comoving_distance(self.z).value**2
        return a

    def jbar(self):
        djdlogmh = self.djc_dlnMh()+self.djsub_dlnMh()
        # print (self.djc_dlnMh()[4, 50, 10:13])
        # print (self.djsub_dlnMh()[4, 50, 10:13])
        dm = np.log10(self.mh[1] / self.mh[0])
        integral = djdlogmh*self.hmfmz  # /(self.mh[None, :, None]*np.log(10))
        res = intg.simps(integral, dx=dm, axis=1, even='avg')
        # res = intg.simps(integral, x=np.log10(self.mh), axis=1, even='avg')
        # res = intg.simps(integral, x=self.mh, axis=1, even='avg')
        # print ("jbar = %s" % (res))
        return res

    def window_cib(self):
        a = 1/(1+self.z)
        return a*self.jbar()

    def onehalo_int(self):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        Cl_1h = np.zeros((nfreq, nfreq, nl))
        dj_cen, dj_sub = self.djc_dlnMh(), self.djsub_dlnMh()
        u = self.unfw
        geo = self.uni.dchi_dz(self.z)/(self.uni.chi(self.z)*(1+self.z))**2
        dm = np.log10(self.mh[1] / self.mh[0])
        fcxcc = self.fc*self.cc
        for i in range(nl):
            for f in range(nfreq):
                rest1 = (dj_cen[f, :, :]*dj_sub*u[:, i, :] + dj_cen *
                         dj_sub[f, :, :]*u[:, i, :] + dj_sub[f, :, :] *
                         dj_sub*u[:, i, :]**2)*self.hmfmz
                intg_mh = intg.simps(rest1, dx=dm, axis=1, even='avg')
                # intg_mh[:, :5] = 0  # cutting contribn from 0 to 0.5 redshift
                intg_z = intg.simps(intg_mh*geo, x=self.z, axis=-1, even='avg')
                Cl_1h[f, :, i] = fcxcc[f]*intg_z*fcxcc
        return Cl_1h

    def J_nu(self):  # , Meffmax, etamax, sigmaMh, alpha):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        Jnu = np.zeros((nfreq, len(self.z), nl))
        dj_cen, dj_sub = self.djc_dlnMh(), self.djsub_dlnMh()
        u = self.unfw
        dm = np.log10(self.mh[1] / self.mh[0])
        for i in range(nl):
            rest1 = (dj_cen + dj_sub*u[:, i, :])*self.biasmz*self.hmfmz
            intg_mh = intg.simps(rest1, dx=dm, axis=1, even='avg')
            Jnu[:, :, i] = intg_mh
        return Jnu

    def twohalo_int(self):  # , snu_eff, ell, unfw, power, fc, cc):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        Cl_2h = np.zeros((nfreq, nfreq, nl))
        Jv = self.J_nu()
        geo = self.uni.dchi_dz(self.z)/(self.uni.chi(self.z)*(1+self.z))**2
        power = self.uni.Pk_array(self.ell, self.z)
        pk_geo = power*geo  # multiplying by geometric factor here itself
        # because in the next step we are going to transpose power spectrum
        # and then we can't multiply dimension (z,ell) with (z,)
        pkt = np.transpose(pk_geo)  # power dimensions are (ell, z) and
        # here we need dimensions of (z, ell). o transpose
        for f in range(nfreq):
            rest1 = Jv*Jv[f, :, :]*pkt
            # rest1[:, :5, :] = 0  # cutting contribn from 0 to 0.5 redshift
            intg_z = intg.simps(rest1, x=self.z, axis=1, even='avg')
            fcxcc = self.fc*self.cc
            Cl_2h[f, :, :] = fcxcc[f]*intg_z*fcxcc[:, None]
        return Cl_2h

    def shot_cib(self):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        shot = np.zeros((nfreq, nfreq, nl))
        # 100, 143, 217, 353, 545, 857 Planck
        # values for 100 and 143 i.e. 3 and 7 are fake
        sar = self.dv.shotpl
        # """
        freq = np.array([100., 143., 217., 353., 545., 857.])
        sa = np.array([1.3*0.116689509208305475, 1.3*0.8714424869942087, 14.,
                       357., 2349., 7407.])
        # shot noise values for 100 and 143 GHz are taken to be approx 1.3
        # times what matt's model predicts. Predictions are given in
        # onehaloandshot.py code in workdrive/cib_hod directory. Factor of 1.3
        # is selected because best fit CIB shot noise from Planck+Herschel
        # data come out to be in the range of (1.15-1.79) times that predicted
        # by Matt's model
        sa[2:] = sar
        res = interp1d(freq, sa, kind='linear', bounds_error=False, fill_value="extrapolate")
        shotval = res(self.nu0)
        """
        freq = self.nu0
        shotval = sar
        # """
        for i in range(nfreq):
            # shot[i, i, :] = sa[i]
            shot[i, i, :] = shotval[i]
            # print (self.nu0[i], shotval[i])
        if max(self.nu0) > max(freq):
            print ("shot noise values for frequencies higher than 857 GHz extrapolated using the values for Planck")
        # shot = np.repeat(sa[2], len(self.ell))
        return shot

    def cl_cibtot(self):
        oneh = self.onehalo_int()  # snu_eff, ell, unfw, fc, cc)
        twoh = self.twohalo_int()  # snu_eff, ell, unfw, power, fc, cc)
        shot = self.shot_cib()  # ell)
        # print ("CIB x CIB")
        # print (oneh[4, 4, -2:], twoh[4, 4, -2:], shot[4, 4, -1:])
        tot = oneh+twoh+shot
        return tot

    def J_nu_iv(self):  # , Meffmax, etamax, sigmaMh, alpha):
        # integrated differential emissivity over all the masses
        dj_cen, dj_sub = self.djc_dlnMh(), self.djsub_dlnMh()
        intgral1 = dj_cen+dj_sub
        intgral1 *= self.hmfmz
        # dm = np.log10(self.mh[1] / self.mh[0])
        # return intg.simps(intgral1, dx=dm, axis=1, even='avg')
        return intg.simps(intgral1, x=np.log10(self.mh), axis=1, even='avg')

    def Iv(self):  # , Meffmax, etamax, sigmaMh, alpha):
        # jnu = self.J_nu_iv()
        # dchi_dz = c_light/(self.cosmo.H0*np.sqrt((self.cosmo.Om0)*(1+self.z)**3 + self.cosmo.Ode0)).value
        jnu = self.jbar()
        dchi_dz = self.uni.dchi_dz(self.z)
        intgral2 = dchi_dz*jnu/(1+self.z)
        # result = self.cc_cibmean*self.freq_Iv*intg.simps(intgral2, x=self.z,
        #                                                  axis=-1, even='avg')
        result = intg.simps(intgral2, x=self.z, axis=-1, even='avg')
        # print ('z: ', self.z)
        # print ('Inu: ', result)
        # print ("I_nu all gal Mat SED (Jy) = %s" % (result))
        # print (np.log(result))
        Ivint = interp1d(self.nu0, result, kind='linear', bounds_error=False,
                         fill_value="extrapolate")
        # result *= ghz*nW/w_jy  # nWm^2/sr
        return Ivint

    def alpha(self):
        """
        Doing integration by parts in the numerator. That gives first term as
        I_obs * nu0 * W(nu0) evaluated at nu0_min and nu0_max. If we assume
        that the filter at both nu0_min and nu0_max is zero, then this first
        term disappears. So I am solving only the second part of integration
        by parts.
        """
        if self.dv.name == 'Planck_only':
            """
            first = self.Iv()(self.nu0max)*self.nu0max*self.filt[self.nucen](self.nu0max)
            second = self.Iv()(self.nu0min)*self.nu0min*self.filt[self.nucen](self.nu0min)
            num1 = first-second
            numint = self.Iv()(self.nu0)
            numint *= self.filt[self.nucen](self.nu0) + self.nu0*self.filtgrad[self.nucen](self.nu0)
            num2 = intg.simps(numint, x=self.nu0, axis=0, even='avg')
            num = num1-num2
            """
            Inugrad = np.gradient(self.Iv()(self.nu0), self.nu0)
            numint = Inugrad*self.nu0*self.filt[self.nucen](self.nu0)
            num = intg.simps(numint, x=self.nu0, axis=0, even='avg')
            # """
            denomint = self.Iv()(self.nu0)*self.filt[self.nucen](self.nu0)
            denom = intg.simps(denomint, x=self.nu0, axis=0, even='avg')
            res = num/denom
        else:
            # print ("here")
            a = np.gradient(np.log(self.Iv()(self.nu0)), np.log(self.nu0), edge_order=2)
            # print (len(self.z))
            # print (a)
            # res *= self.nu0/self.Iv()(self.nu0)
            res = interp1d(self.nu0, a, kind='linear', bounds_error=False, fill_value="extrapolate")
        return res

    def Tdust(self, z):
        T0 = 24.4  # Planck CIB 2013 paper
        alpha = 0.36
        result = T0*(1.+z)**alpha
        return result


    def B_nu(self, z, nu):
        # nu in Hz and not GHz
        Td = self.Tdust(z)
        res = 2.*h_p*nu**3/(c_light*1.e3)**2
        x = h_p*nu/k_B/Td
        # print (min(x), max(x))
        res /= (np.exp(x) - 1)
        return res

    def mod_blackbody(self, z, nu):
        beta = 1.75
        Bnu = self.B_nu(z, nu)
        result = Bnu
        result *= nu**beta
        # result *= w_jy  # Watt to Jy
        return result

    def nu0_z(self, z):
        """
        for mod blackboy approximation, for the SED we have
        dlntheta/dlnnu = -gamma for nu=nu0.
        Here theta is the modified blackbody spectrum. In order to find nu0
        which isredshift dependent, we need to take a derivative and solve
        for this numerically. In the end it comes out in the form
        (x-(3+beta+gamma))e(x-(3+beta+gamma)) = -(3+beta+gamma)e(-(3+beta+gamma))
        The solution is x-(3+beta+gamma) = W(-(3+beta+gamma)e(-(3+beta+gamma)))
        here W is Lambert's W fnction which is implemented in scipy.
        x = hnu/KT
        """
        beta = 1.75
        gamma = 1.7
        y = -(3+beta+gamma)*np.exp(-(3+beta+gamma))
        xx = lambertw(y)
        x = xx + (3+beta+gamma)
        Td = self.Tdust(z)
        nu0z = np.real(x*k_B*Td/h_p)
        return nu0z

    def Theta(self, nu, z):
        # calculating only for nu < nu0
        # normalised SED such that theta(nu0) = 1
        num = self.mod_blackbody(z, nu)
        # nu0 = 2000*1.e9  # Hz => approximately taken from Fig.1 of 2010.16405
        nu0z = self.nu0_z(z)
        denom = self.mod_blackbody(z, nu0z)
        return num/denom

    def phiz(self, z):
        delta = 3.6
        return (1+z)**delta

    def sigmaM(self, m):
        log10Meff = 12.6
        # Mmin = 1e8  # randomly chosen
        sigmaLM2 = 0.5

        exp_coef = (np.log10(m) - log10Meff)**2/(2*sigmaLM2)
        norm = m/np.sqrt(2*np.pi*sigmaLM2)
        sigmam = norm*np.exp(-1.*exp_coef)
        return sigmam

    def LMrel(self, nu, m, z):
        # nu *= 1.e9  # frequency in Hz not GHz
        # best fit values of parameyers taken from 2010.16405
        L0 = 6.4e-8  # JyMpc^2/M_sun/Hz

        phi_z = self.phiz(z)
        # print ("phi(z) = %s" % (phi_z))
        sigmaM = self.sigmaM(m)
        # print ("sigmaM = %s" % (sigmaM))
        Theta = self.Theta(nu, z)
        # print ("Theta = %s" % (Theta))
        res = L0*phi_z*sigmaM*Theta
        # print ("L = %s" % (res))
        return res

    def jnu_modblack(self):
        # dn_dm = self.hmfmz/(self.mh*np.log(10))
        nz = len(self.z)
        Lnu = np.zeros((len(self.nu0), len(self.mh), nz))
        for iz in range(nz):
            for iv in range(len(self.nu0)):
                nurest = self.nu0[iv]*(1+self.z[iz])*1.e9
                Lnu[iv, :, iz] = self.LMrel(nurest, self.mh, self.z[iz])
        integrand = Lnu*self.hmfmz/(4.*np.pi)
        res = intg.simps(integrand, x=np.log10(self.mh), axis=1, even='avg')
        # print ("Jnu modblack at 100 GHz = %s" % (res[0, :]))
        return res
        
    def Iv_modblack(self):
        jnu = self.jnu_modblack()
        dchi_dz = self.uni.dchi_dz(self.z)
        intgral2 = dchi_dz*jnu/(1+self.z)
        # result = self.cc_cibmean*self.freq_Iv*intg.simps(intgral2, x=self.z,
        #                                                  axis=-1, even='avg')
        result = intg.simps(intgral2, x=self.z, axis=-1, even='avg')
        # print ("I_nu all gal mod blackbod (Jy) = %s" % (result))
        # print (np.log(result))
        Ivint = interp1d(self.nu0, result, kind='linear', bounds_error=False,
                         fill_value="extrapolate")
        return Ivint

    def alpha_modblack(self):
        a = np.gradient(np.log(self.Iv_modblack()(self.nu0)), np.log(self.nu0), edge_order=2)
        # print (a)
        # res *= self.nu0/self.Iv()(self.nu0)
        res = interp1d(self.nu0, a, kind='linear', bounds_error=False, fill_value="extrapolate")
        return res

    def deltaIdoppler(self, thetagal, nu, m, z):
        # nu has to be in Hz and not GHz
        # solid angle = A/chi^2=pi*(chi*theta)**2/chi**2=pi*theta**2
        thetagal /= 60.  # in degrees
        thetagal *= np.pi/180.  # radians
        Omega = np.pi*thetagal**2

        nurest = nu*(1+z)  # rest frame freq
        Lnurest = self.LMrel(nurest, m, z)
        chi = self.uni.chi(z)
        flux = Lnurest/(4*np.pi*chi**2*(1+z))
        # print ("flux for single gal at 100 GHz = %s" % (flux))
        Inu = flux/Omega
        # Inu2 = self.Iv_modblack()
        # print ("I_nu single gal (Jy) = %s" % (Inu))
        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            # alpha_cib = self.alpha()(self.nu0)
            # print ("alpha Mat = %s" % (alpha_cib))
            alpha_cib = self.alpha_modblack()(nu/1.e9)  # nu in GHz
            # print ("alpha mod blackbod = %s" % (alpha_cib))
        # print ("alpha_{CIB}: %s" % (alpha_cib))
        # print ("alpha = %s" % (alpha_cib))
        beta2 = self.uni.beta2(self.z)
        beta = np.sqrt(beta2)
        fbeta = interp1d(self.z, beta, kind='linear', bounds_error=False, fill_value="extrapolate")
        betaz = fbeta(z)
        # print ("beta = %s" % (betaz))

        res = (3.-alpha_cib)*Inu
        res *= betaz

        return res

    def deltaIdoppler_approx(self):
        """
        It is not really clear if the integration over redshift while
        calculating CIB intensity should be over the whole redshift or just
        a small shell around which you are taking the velocity component.
        """
        """
        (3-alpha)*beta*I_obs(M, z)
        Here I_obs is calculated for specific emissivity dj(M, nu, z)/dlogM
        rather than j(nu, M)
        i.e. I(M, nu) = int {dz dchi/dz a dj(M, nu, z)/dlogM}
        Then you can select the mass with maximum value of I for max effect.
        Finally, you should select beta at a redshift for typical CMASS or
        DESI galaxies.
        """
        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)
        # print ("alpha_{CIB}: %s" % (alpha_cib))
        beta2 = self.uni.beta2(self.z)
        beta = np.sqrt(beta2)
        # print (beta[60])  # z~1.0
        if len(self.z) == 210:
            ind = 60
            print ("beta = %s" % (beta[ind]))
        else:
            ind = 9
        """
        Inu = self.Iv()(self.nu0)
        # print (Inu)
        res = (3.-alpha_cib)*Inu
        # res *= Inu
        """

        deltaz = 0.05 # probably only need to multiply a*beta*dj/dlogm by delta chi where delta chi is around a thin shell around z
        deltachi = self.uni.chi(self.z[ind]+deltaz)-self.uni.chi(self.z[ind]-deltaz)

        dj_cen, dj_sub = self.djc_dlnMh(), self.djsub_dlnMh()
        djtot = dj_cen+dj_sub  # nu x mh x z
        # print (djtot.shape)
        djtot_T = np.transpose(djtot, (0, 2, 1))  # nu x z x mh
        # print (djtot_T.shape, np.log10(self.mh).shape)
        fdjtot = interp1d(np.log10(self.mh), djtot_T, kind='linear', bounds_error=False, fill_value="extrapolate")
        hmfmz_T = self.hmfmz.T  # z x mh
        fhmfmz = interp1d(np.log10(self.mh), hmfmz_T, kind='linear', bounds_error=False, fill_value="extrapolate")
        mh_cmass = 4.e13
        deltam = mh_cmass*0.1
        mhbin = np.linspace(mh_cmass-deltam, mh_cmass+deltam, 20)
        logmhbin = np.log10(mhbin)
        jnu_num_int = fdjtot(logmhbin)*fhmfmz(logmhbin)
        jnu_num = intg.simps(jnu_num_int, x=logmhbin, axis=-1, even='avg')

        jnu_denom_int = fhmfmz(logmhbin).T*self.uni.chi(self.z[ind])**2*deltachi
        jnu_denom = intg.simps(jnu_denom_int, x=logmhbin, axis=0, even='avg')

        # print (jnu_num.shape, jnu_denom.shape)
        jnu = jnu_num/jnu_denom  # nu x z
        # jnu = jnu.T  # nu x z
        """
        dchi_dz = self.uni.dchi_dz(self.z)
        intgral2 = dchi_dz*jnu/(1+self.z)
        Inu = intg.simps(intgral2, x=self.z, axis=-1, even='avg')
        res = (3.-alpha_cib)*Inu[:, ind]
        # """

        Inu = jnu/(1+self.z)*deltachi  # nu x z
        res = (3.-alpha_cib)*Inu[:, ind]
        # res = (3.-alpha_cib[:, None])*Inu[:, ind]
        # """

        res *= beta[ind]
        # print (res.shape)
        # fdjtot = interp1d(self.z, djtot, kind='linear', bounds_error=False, fill_value="extrapolate")
        return res

    def deltaIdoppler_all(self):
        """
        (3-alpha)*beta*I_obs
        """
        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)
        # print ("alpha_{CIB}: %s" % (alpha_cib))
        beta2 = self.uni.beta2(self.z)
        beta = np.sqrt(beta2)
        # print (beta[60])  # z~1.0
        Inu = self.Iv()(self.nu0)
        # print (Inu)
        res = (3.-alpha_cib)*Inu
        # res *= Inu
        if len(self.z) == 210:
            res *= beta[60]
            print ("beta = %s" % (beta[60]))
        else:
            res *= beta[9]
        return res

    def onehalo_int_Doppler(self):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        Cl_1h = np.zeros((nfreq, nfreq, nl))
        dj_cen, dj_sub = self.djc_dlnMh(), self.djsub_dlnMh()
        u = self.unfw

        geo = self.uni.dchi_dz(self.z)/(self.uni.chi(self.z)*(1+self.z))**2
        dm = np.log10(self.mh[1] / self.mh[0])
        fcxcc = self.fc*self.cc

        beta2 = self.uni.beta2(self.z)

        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)

        alpha_fac = (3-alpha_cib)  # **2

        for i in range(nl):
            for f in range(nfreq):
                rest1 = (dj_cen[f, :, :]*dj_sub*u[:, i, :] + dj_cen *
                         dj_sub[f, :, :]*u[:, i, :] + dj_sub[f, :, :] *
                         dj_sub*u[:, i, :]**2)*self.hmfmz
                intg_mh = intg.simps(rest1, dx=dm, axis=1, even='avg')
                # intg_mh[:, :5] = 0  # cutting contribn from 0 to 0.5 redshift
                res = intg_mh*geo*beta2
                intg_z = intg.simps(res, x=self.z, axis=-1, even='avg')

                Cl_1h[f, :, i] = alpha_fac[f]*fcxcc[f]*intg_z*fcxcc*alpha_fac
        return Cl_1h

    def twohalo_int_Doppler(self):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        Cl_2h = np.zeros((nfreq, nfreq, nl))
        Jv = self.J_nu()
        fcxcc = self.fc*self.cc
        geo = self.uni.dchi_dz(self.z)/(self.uni.chi(self.z)*(1+self.z))**2
        power = self.uni.Pk_array(self.ell, self.z)
        pk_geo = power*geo  # multiplying by geometric factor here itself
        # because in the next step we are going to transpose power spectrum
        # and then we can't multiply dimension (z,ell) with (z,)

        beta2 = self.uni.beta2(self.z)
        pk_geo *= beta2

        pkt = np.transpose(pk_geo)  # power dimensions are (ell, z) and
        # here we need dimensions of (z, ell). o transpose

        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)

        alpha_fac = (3-alpha_cib)  # **2

        for f in range(nfreq):
            rest1 = Jv*Jv[f, :, :]*pkt
            # rest1[:, :5, :] = 0  # cutting contribn from 0 to 0.5 redshift
            intg_z = intg.simps(rest1, x=self.z, axis=1, even='avg')
            Cl_2h[f, :, :] = alpha_fac[f]*fcxcc[f]*intg_z*fcxcc[:, None]*alpha_fac[:, None]
        return Cl_2h

    def cl_cibDoppler_tot(self):
        oneh = self.onehalo_int_Doppler()
        twoh = self.twohalo_int_Doppler()
        # print ("CIB x CIB Doppler")
        # print (oneh[4, 4, -2:], twoh[4, 4, -2:])
        tot = oneh+twoh
        return tot
