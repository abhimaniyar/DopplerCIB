from headers_constants import *
from CIB_halo import *
from Gal_halo import *


class CIBxgal(ProfHODMore15, Cib_halo):

    def __init__(self, data_var_iv, cosmo_var_iv, gal_exp, r_l):
        # halomassfunc.__init__(k, pk, z, cosmo, delta_h)
        # nfw_fourier.__init__(k, pk, z, cosmo, delta_h)
        ProfHODMore15.__init__(self, data_var_iv, cosmo_var_iv, gal_exp)
        Cib_halo.__init__(self, data_var_iv, cosmo_var_iv)
        self.r_l = r_l
        # self.gal = ProfHODMore15  # .__init__(self, data_var_iv, cosmo_var_iv)
        # self.cib = Cib_halo  # .__init__(self, data_var_iv, cosmo_var_iv)
        # self.dv = data_var_iv
        # self.uni = self.gal.uni  # cosmo_var_iv
        # self.snu_eff = cib.snu_unfilt
        # self.ucen = cib.unfw
        # self.unfw = cib.unfw
        # self.
        # self.pk = cib.

    """
    def ucib(self, mh, z):
        nfw_u = np.zeros((len(mh), len(self.k), len(z)))
        for r in range(len(z)):
            test_nfw = nfw_fourier(self.k, self.Pk[:, r], z[r], self.cosmo,
                                   self.delta_h)
            nfw_u[:, :, r] = test_nfw.nfwfourier_u(mh)
        return nfw_u
    """

    def cibterm(self, unfw):
        snu_eff = self.snu
        # unfw = self.unfw
        nfreq = len(snu_eff[:, 0])
        djcentral = self.djc_dlnMh()  #snu_eff)
        djsub = self.djsub_dlnMh()  # snu_eff)
        nl = len(unfw[0, :, 0])
        cibres = np.zeros((nfreq, nl, len(self.mh), len(self.z)))
        for i in range(nl):
            cibres[:, i, :, :] = djcentral+djsub*unfw[:, i, :]
        return cibres

    def galterm(self, ucen, unfw):
        # unfw = self.unfw
        # ucen = self.unfw
        result = np.zeros((len(unfw[0, :, 0]), len(self.mh), len(self.z)))
        for m in range(len(self.mh)):
            result[:, m, :] = self.Ncen(self.mh[m])*ucen[m, :, :]+self.Nsat(self.mh[m])*unfw[m, :, :]
        return result

    def cibgalcross_pk_1h(self, ucen, unfw):
        # snu_eff = self.snu_unfilt
        # unfw = self.unfw
        # ucen = self.unfw
        res1 = self.hmfmz*self.cibterm(unfw)*self.galterm(ucen, unfw)
        # dm = np.log10(self.mh[1] / self.mh[0])
        # intg_mh = intg.simps(res1, dx=dm, axis=2, even='avg')
        intg_mh = intg.simps(res1, x=np.log10(self.mh), axis=2, even='avg')
        return intg_mh/self.nbargal()/self.jbar()[:, None, :]

    def cibgalcross_pk_2h(self, ucen, unfw, pk):
        """
        Check the calculation first using two different integrals for both CIB
        and galaxies, then check with a single integral for both after
        multiplying. If both are same, then we can save one integral.
        """
        # pk = self.uni.Pk_array(self.ell, self.z)
        # print (self.hmfmz.shape, self.biasmz.shape)
        res1 = self.hmfmz*self.biasmz*self.cibterm(unfw)  # snu_eff, unfw)
        dm = np.log10(self.mh[1] / self.mh[0])
        intg_mh1 = intg.simps(res1, dx=dm, axis=2, even='avg')
        res2 = self.hmfmz*self.biasmz*self.galterm(ucen, unfw)  # ucen, unfw)
        intg_mh2 = intg.simps(res2, dx=dm, axis=1, even='avg')
        res3 = intg_mh1*intg_mh2*pk
        return res3/self.nbargal()/self.jbar()[:, None, :]

    def window_cibxgal(self):
        window = self.window_gal()
        w_gal = window(self.z)/self.uni.dchi_dz(self.z)
        # print (self.z, w_gal)
        w_cib = self.window_cib()
        # print (w_gal[10:13])
        # print (w_cib[4, 10:13])
        return w_cib*w_gal

    def cibgalcross_cell_1h(self):  # , snu_eff, ucen, unfw):
        # ucen = self.gal.unfw
        # unfw = self.gal.unfw
        ucen = self.unfw
        unfw = self.unfw
        w1w2 = self.window_cibxgal()
        # print (ucen[50, 10, 10:14])
        # print (unfw[50, 10, 10:14])
        # print (w1w2.shape)
        # print (w1w2[4, 10:14])
        geo = self.uni.dchi_dz(self.z)/self.uni.chi(self.z)**2
        fact = geo*w1w2
        res = fact[:, None, :]*self.cibgalcross_pk_1h(ucen, unfw)
        return intg.simps(res, x=self.z, axis=2, even='avg')

    def cibgalcross_cell_2h(self):  # , snu_eff, ucen, unfw, pk):
        # ucen = self.gal.unfw
        # unfw = self.gal.unfw
        ucen = self.unfw
        unfw = self.unfw
        power = self.uni.Pk_array(self.ell, self.z)
        w1w2 = self.window_cibxgal()
        geo = self.uni.dchi_dz(self.z)/self.uni.chi(self.z)**2
        fact = geo*w1w2
        res = fact[:, None, :]*self.cibgalcross_pk_2h(ucen, unfw, power)
        return intg.simps(res, x=self.z, axis=2, even='avg')

    def cibgalcross_cell_shot(self):  # , snu_eff, ucen, unfw, pk):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        shotcib = np.zeros((nfreq, nl))
        # 100, 143, 217, 353, 545, 857 Planck
        # values for 100 and 143 i.e. 3 and 7 are fake
        sar = self.dv.shotpl
        # """
        freq = np.array([100., 143., 217., 353., 545., 857.])
        sa = np.array([3., 7., 14., 357., 2349., 7407.])
        sa[2:] = sar
        res = interp1d(freq, sa, kind='linear', bounds_error=False, fill_value="extrapolate")
        shotval = res(self.nu0)
        """
        freq = self.nu0
        shotval = sar
        # """
        for i in range(nfreq):
            # shot[i, i, :] = sa[i]
            shotcib[i, :] = shotval[i]
            # print (self.nu0[i], shotval[i])
        if max(self.nu0) > max(freq):
            print ("shot noise values for frequencies higher than 857 GHz extrapolated using the values for Planck")

        r_l = self.r_l
        # shotcib = self.shot_cib()
        shotgal = np.repeat(self.clshot_gal(), len(self.ell))
        crossshot = r_l*np.sqrt(shotcib*shotgal)
        return crossshot

    def cibgalcross_cell_tot(self):  # , snu_eff, ucen, unfw, pk):
        oneh = self.cibgalcross_cell_1h()  # snu_eff, ucen, unfw)
        twoh = self.cibgalcross_cell_2h()  # snu_eff, ucen, unfw, pk)
        shot = self.cibgalcross_cell_shot()
        # print ("CIB x gal")
        # print (oneh[5, -2:], twoh[4, -2:], shot[4, -1:])
        tot = oneh+twoh+shot
        return tot

    def cibvelcross_cell_1h(self):  # , snu_eff, ucen, unfw):
        # ucen = self.gal.unfw
        # unfw = self.gal.unfw
        ucen = self.unfw
        unfw = self.unfw
        w1w2 = self.window_cibxgal()

        beta2 = self.uni.beta2(self.z)

        geo = self.uni.dchi_dz(self.z)/self.uni.chi(self.z)**2
        fact = geo*w1w2
        res = fact[:, None, :]*self.cibgalcross_pk_1h(ucen, unfw)
        res *= beta2

        result = intg.simps(res, x=self.z, axis=2, even='avg')
        # print ("calculating velcity 1 halo")
        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)
        # print ('alpha: ', alpha_cib)
        result *= (3.-alpha_cib)[:, None]

        return result

    def cibvelcross_cell_2h(self):  # , snu_eff, ucen, unfw, pk):
        # ucen = self.gal.unfw
        # unfw = self.gal.unfw
        ucen = self.unfw
        unfw = self.unfw
        power = self.uni.Pk_array(self.ell, self.z)
        w1w2 = self.window_cibxgal()

        beta2 = self.uni.beta2(self.z)

        geo = self.uni.dchi_dz(self.z)/self.uni.chi(self.z)**2
        fact = geo*w1w2
        res = fact[:, None, :]*self.cibgalcross_pk_2h(ucen, unfw, power)
        res *= beta2

        result = intg.simps(res, x=self.z, axis=2, even='avg')
        # print (result.shape)
        # print (result[:, 0])
        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)
        # print (alpha_cib)
        # print ((3.-alpha_cib))
        result *= (3.-alpha_cib)[:, None]

        return result

    def cibvelcross_cell_shot(self):
        nfreq = len(self.nu0)
        nl = len(self.ell)
        shotcib = np.zeros((nfreq, nl))
        # 100, 143, 217, 353, 545, 857 Planck
        # values for 100 and 143 i.e. 3 and 7 are fake
        sar = self.dv.shotpl
        freq = self.nu0
        # freq = np.array([100., 143., 217., 353., 545., 857.])
        # sa = np.array([3., 7., 14., 357., 2349., 7407.])
        # sa[2:] = sar
        # res = interp1d(freq, sa, kind='linear', bounds_error=False, fill_value="extrapolate")
        # shotval = res(self.nu0)
        shotval = sar

        if self.dv.name == 'Planck_only':
            alpha_cib = self.alpha()
        else:
            alpha_cib = self.alpha()(self.nu0)
        # print ("alpha_cib: ", alpha_cib)
        for i in range(nfreq):
            shotcib[i, :] = shotval[i]

        if max(self.nu0) > max(freq):
            print ("shot noise values for frequencies higher than 857 GHz extrapolated using the values for Planck")

        r_l = self.r_l
        # shotcib = self.shot_cib()
        shotvel = np.repeat(self.clshot_vel(), len(self.ell))
        """
        Because alpha > 3 near 217 GHz, 3-alpha ges negative giving a negative
        value of cross-shot noise.
        This is true for 1-halo and 2-halo termsa s well.
        """
        crossshot = r_l*np.sqrt(shotcib*shotvel)*(3-alpha_cib[:, None])
        # print ("shotcib:", shotcib)
        # print (shotvel)
        return crossshot

    def cibvelcross_cell_tot(self):  # , snu_eff, ucen, unfw, pk):
        oneh = self.cibvelcross_cell_1h()  # snu_eff, ucen, unfw)
        twoh = self.cibvelcross_cell_2h()  # snu_eff, ucen, unfw, pk)
        # shot = self.cibvelcross_cell_shot()
        tot = oneh+twoh  # +shot
        return tot

    def dvdz(self):
        return self.dchi_dz()*self.cosmo.comoving_distance(self.z).value**2

    def dndz2(self):
        a = self.nbargal(self.mh, self.hmfmz)
        return a*self.dvdz()

    def snr__g(self, n_nu, fsky, instcib, instgal):

        cross = self.cibgalcross_cell_tot()
        cibwhite = self.dv.whitenoise(n_nu)
        # print (cibwhite)
        """
        cibtot = instcib.cl_cibtot()+cibwhite
        galtot = instgal.cl_galtot()
        num = (2*self.ell+1)*fsky*cross[n_nu, :]**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = cross[n_nu, :]**2 + cibtot[n_nu, n_nu, :]*galtot
        snr2 = np.cumsum(num/denom)
        snr = np.sqrt(snr2)
        """
        fcross = UnivariateSpline(self.ell, cross[n_nu, :], k=1, s=0)
        # fcibw = UnivariateSpline(self.ell, cibwhite, k=1, s=0)
        cibonly = instcib.cl_cibtot()[n_nu, n_nu, :]
        cibtot = cibonly + cibwhite
        nl = len(self.ell)
        
        # for i in range(nl):
        #     if cibonly[i] < 1.e-4*cibwhite[i]:
        #         cibtot[i] = 1.e10  # some random big number

        fcib = UnivariateSpline(self.ell, cibtot, k=1, s=0)
        # print (cibtot[n_nu, n_nu, :])
        galtot = instgal.cl_galtot()
        fgal = UnivariateSpline(self.ell, galtot, k=1, s=0)
        # fcrossshot = UnivariateSpline(self.ell, crossshot, k=1, s=0)
        # el = np.linspace(50., 40000., 40000-50+1)
        # deltal = l[2]-l[1]
        # fcr = fcross(el)
        # fctot = fcibtot(el)
        # fw = fcibw(el)
        # fg = fgal(el)
        lcen = np.zeros(nl-1)
        fcrossbin = np.zeros(nl-1)
        # fwhitebin = np.zeros(nl-1)
        fcibtotbin = np.zeros(nl-1)
        fgaltotbin = np.zeros(nl-1)
        # fcrossshotbin = np.zeros(nl-1)
        for i in range(nl-1):
            l1 = self.ell[i]
            l2 = self.ell[i+1]
            lcen[i] = (l1+l2)/2.
            el = np.linspace(l1, l2, int(l2-l1+1))
            deltal = l2-l1
            fcrossbin[i] = np.sum(fcross(el))/deltal
            fcibtotbin[i] = np.sum(fcib(el))/deltal
            # fwhitebin[i] = np.sum(fcibw(el))/deltal
            fgaltotbin[i] = np.sum(fgal(el))/deltal
            # fcrossshotbin[i] = np.sum(fcrossshot(el))/deltal

        cibbintot = fcibtotbin  # + fwhitebin
        # crosstot = fcrossbin + fcrossshotbin
        # print (fcibtotbin)
        # print (fwhitebin)
        # print (cibbintot)
        # print (fgaltotbin)
        # print (fcrossbin)
        num = (2*lcen+1)*deltal*fsky*fcrossbin**2
        # for i in range(nl-1):
        #     if cibbintot[i] > 1.e20:
        #         num[i] = 0.
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = fcrossbin**2 + cibbintot*fgaltotbin
        # print (cibbintot)
        # print (fgaltotbin)
        # print (num)
        # print (denom)
        snr2 = np.cumsum(num/denom)
        snrbin = np.sqrt(snr2)

        # return snr, snrbin
        return snrbin

    def snr__v(self, n_nu, fsky, instcib, instgal):

        cross = self.cibvelcross_cell_tot()
        cibwhite = self.dv.whitenoise(n_nu)
        # print (cibwhite)
        """
        cibtot = instcib.cl_cibtot()+cibwhite
        veltot = instgal.cl_veltot()
        # print (veltot[:10])
        num = (2*self.ell+1)*fsky*cross[n_nu, :]**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = cross[n_nu, :]**2 + cibtot[n_nu, n_nu, :]*veltot
        snr2 = np.cumsum(num/denom)
        snr = np.sqrt(snr2)
        """

        fcross = UnivariateSpline(self.ell, cross[n_nu, :], k=1, s=0)
        fcibw = UnivariateSpline(self.ell, cibwhite, k=1, s=0)

        cibonly = instcib.cl_cibtot()[n_nu, n_nu, :]
        cibtot = cibonly + cibwhite

        fcib = UnivariateSpline(self.ell, cibtot, k=1, s=0)

        veltot = instgal.cl_veltot()
        fvel = UnivariateSpline(self.ell, veltot, k=1, s=0)
        # el = np.linspace(50., 40000., 40000-50+1)
        # deltal = l[2]-l[1]
        # fcr = fcross(el)
        # fctot = fcibtot(el)
        # fw = fcibw(el)
        # fg = fgal(el)
        nl = len(self.ell)
        lcen = np.zeros(nl-1)
        fcrossbin = np.zeros(nl-1)
        fwhitebin = np.zeros(nl-1)
        fcibtotbin = np.zeros(nl-1)
        fveltotbin = np.zeros(nl-1)
        for i in range(nl-1):
            l1 = self.ell[i]
            l2 = self.ell[i+1]
            lcen[i] = (l1+l2)/2.
            el = np.linspace(l1, l2, int(l2-l1+1))
            deltal = l2-l1
            fcrossbin[i] = np.sum(fcross(el))/deltal
            fcibtotbin[i] = np.sum(fcib(el))/deltal
            fwhitebin[i] = np.sum(fcibw(el))/deltal
            fveltotbin[i] = np.sum(fvel(el))/deltal

        cibbintot = fcibtotbin  # + fwhitebin
        # print (fcibtotbin)
        # print (fwhitebin)
        # print (cibbintot)
        # print (fgaltotbin)
        # print (fcrossbin)
        num = (2*lcen+1)*deltal*fsky*fcrossbin**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = fcrossbin**2 + cibbintot*fveltotbin
        snr2 = np.cumsum(num/denom)
        snrbin = np.sqrt(snr2)
        # return snr, snrbin
        return snrbin

    def snr_g(self, n_nu, fsky):
        """
        total cib and galaxy power spectra here are the ones calculated for the
        specs for galaxy and cib given for calculating their cross. So this
        might over-estimate SNR. For example, if galaxy is being calculated
        from redshift 1-2, then total CIB power will correspond to only that
        redshift. This is not really right, because the total CIB power
        spectrum entering the SNR formula is actually from the real CIB map
        where CIB comes from all redshifts i.e. 0-10. So that is why this
        calculation here is approximate and gives slightly higher SNR.
        This has been corrected in snr__g module below.
        """
        cross = self.cibgalcross_cell_tot()
        cibwhite = self.dv.whitenoise(n_nu)
        # print (cibwhite)
        cibtot = self.cl_cibtot()+cibwhite
        galtot = self.cl_galtot()
        num = (2*self.ell+1)*fsky*cross[n_nu, :]**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = cross[n_nu, :]**2 + cibtot[n_nu, n_nu, :]*galtot
        snr2 = np.cumsum(num/denom)
        snr = np.sqrt(snr2)

        fcross = UnivariateSpline(self.ell, cross[n_nu, :], k=1, s=0)
        fcibw = UnivariateSpline(self.ell, cibwhite, k=1, s=0)
        fcib = UnivariateSpline(self.ell, cibtot[n_nu, n_nu, :], k=1, s=0)
        fgal = UnivariateSpline(self.ell, galtot, k=1, s=0)
        # el = np.linspace(50., 40000., 40000-50+1)
        # deltal = l[2]-l[1]
        # fcr = fcross(el)
        # fctot = fcibtot(el)
        # fw = fcibw(el)
        # fg = fgal(el)
        nl = len(self.ell)
        lcen = np.zeros(nl-1)
        fcrossbin = np.zeros(nl-1)
        fwhitebin = np.zeros(nl-1)
        fcibtotbin = np.zeros(nl-1)
        fgaltotbin = np.zeros(nl-1)
        for i in range(nl-1):
            l1 = self.ell[i]
            l2 = self.ell[i+1]
            lcen[i] = (l1+l2)/2.
            el = np.linspace(l1, l2, int(l2-l1+1))
            deltal = l2-l1
            fcrossbin[i] = np.sum(fcross(el))/deltal
            fcibtotbin[i] = np.sum(fcib(el))/deltal
            fwhitebin[i] = np.sum(fcibw(el))/deltal
            fgaltotbin[i] = np.sum(fgal(el))/deltal

        cibbintot = fcibtotbin + fwhitebin
        # print (fcibtotbin)
        # print (fwhitebin)
        # print (cibbintot)
        # print (fgaltotbin)
        # print (fcrossbin)
        num = (2*lcen+1)*deltal*fsky*fcrossbin**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = fcrossbin**2 + cibbintot*fgaltotbin
        snr2 = np.cumsum(num/denom)
        snrbin = np.sqrt(snr2)
        return snr, snrbin

    def snr_v(self, n_nu, fsky):

        cross = self.cibvelcross_cell_tot()
        cibwhite = self.dv.whitenoise(n_nu)
        # print (cibwhite)
        cibtot = self.cl_cibtot()+cibwhite
        veltot = self.cl_veltot()
        num = (2*self.ell+1)*fsky*cross[n_nu, :]**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = cross[n_nu, :]**2 + cibtot[n_nu, n_nu, :]*veltot
        snr2 = np.cumsum(num/denom)
        snr = np.sqrt(snr2)

        fcross = UnivariateSpline(self.ell, cross[n_nu, :], k=1, s=0)
        fcibw = UnivariateSpline(self.ell, cibwhite, k=1, s=0)
        fcib = UnivariateSpline(self.ell, cibtot[n_nu, n_nu, :], k=1, s=0)
        fvel = UnivariateSpline(self.ell, veltot, k=1, s=0)
        # el = np.linspace(50., 40000., 40000-50+1)
        # deltal = l[2]-l[1]
        # fcr = fcross(el)
        # fctot = fcibtot(el)
        # fw = fcibw(el)
        # fg = fgal(el)
        nl = len(self.ell)
        lcen = np.zeros(nl-1)
        fcrossbin = np.zeros(nl-1)
        fwhitebin = np.zeros(nl-1)
        fcibtotbin = np.zeros(nl-1)
        fveltotbin = np.zeros(nl-1)
        for i in range(nl-1):
            l1 = self.ell[i]
            l2 = self.ell[i+1]
            lcen[i] = (l1+l2)/2.
            el = np.linspace(l1, l2, int(l2-l1+1))
            deltal = l2-l1
            fcrossbin[i] = np.sum(fcross(el))/deltal
            fcibtotbin[i] = np.sum(fcib(el))/deltal
            fwhitebin[i] = np.sum(fcibw(el))/deltal
            fveltotbin[i] = np.sum(fvel(el))/deltal

        cibbintot = fcibtotbin + fwhitebin
        # print (fcibtotbin)
        # print (fwhitebin)
        # print (cibbintot)
        # print (fgaltotbin)
        # print (fcrossbin)
        num = (2*lcen+1)*deltal*fsky*fcrossbin**2
        # print (np.shape(cross), np.shape(cibtot), np.shape(galtot))
        denom = fcrossbin**2 + cibbintot*fveltotbin
        snr2 = np.cumsum(num/denom)
        snrbin = np.sqrt(snr2)
        return snr, snrbin
