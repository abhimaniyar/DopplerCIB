from headers_constants import *
from input_var_cibmean import *
from Inu_cib import *
from cosmo_related import *
import time

time0 = time.time()

# color corrections for 100, 143, 217, 353, 545, 857 and 3000 GHz for Planck
# cc_pl = np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995])  # , 0.960])
cc_pl = np.ones(6)
fc_pl = np.ones(len(cc_pl))

"""
Calculating the observed CIB intensity for halos with a given mass at given
redshifts for different Planck frequencies. The SEDs used for the Planck
channels are bandpassed, so the observed intensities are calculated as
they would be observed with Planck frequency channels at 100, 143, 353, 545,
847 GHz as well as 3000 GHz for IRAS. Intensity is calculated in nW/m^2/sr.
"""
nuarray = np.array([100., 143., 217., 353., 545., 857.])  # , 1000., 1500., 2000., 2500., 3000.])
# nuarray = np.array([220., 280., 350., 410., 850.])  # in GHz
deltanu = 800
nucen = 100
"""
CAN SPPED UP THE COMPUTATION BY MAKING NUCEN AS AN ARRAY. IT JUST GIVES THE
CENTRAL FREQUENCY OF THE FILTER. SO WE CAN CALCULATE ALPHA FOR ALL THE
PLANCK FREQUENCIES TOGETHER IN NEED BE.
"""
nu0min = 50.  # nucen-deltanu/2.
nu0max = 3000.  # nucen+deltanu/2.
steps = 20  # 2000 nu0max-nu0min+1  # nu0max-nu0min+1  # 200
nu0 = np.linspace(nu0min, nu0max, int(steps))  # nuarray  # np.linspace(nu0min, nu0max, int(steps))

cib_exp = 'CCAT'  # CCAT  # Planck
gal_exp = 'CMASS'  # CMASS  # DESI_ELG  # DESI_LRG

if cib_exp == 'CCAT':
    ell = np.linspace(50, 40000, 200)
    fwhm = 0.5  # 4.8 Planck arcmin # 0.5 CCAT-Prime
    sensitivity = 1.2  # 13.5 Planck Jy/sr  # 1.2 CCAT-Prime
elif cib_exp == 'Planck':
    ell = np.linspace(50, 3000, 15)
    fwhm = 4.8  # 4.8 Planck arcmin # 0.5 CCAT-Prime
    sensitivity = 13.5  # 13.5 Planck Jy/sr  # 1.2 CCAT-Prime

strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nolens_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"
# strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"

cibres = "data_files/one_halo_bestfit_"+strfig
# clres = np.loadtxt('data/%s.txt' % (string))


custom = {'name': cib_exp,
          'do_cibmean': 1,
          'cc': cc_pl,
          'fc': fc_pl,
          'snuaddr': 'data_files/filtered_snu_planck.fits',
          'nu0min': nu0min, 'nu0max': nu0max,
          'nucen': str(int(nucen)),
          'nu0': nu0,
          'ell': ell,
          'fwhm': fwhm,
          'sensitivity': sensitivity,
          'cibpar_resfile': cibres}

exp = custom
# cc_cibmean = np.array([0.97125, 0.999638, 1.00573, 0.959529, 0.973914, 0.988669, 0.987954, 1., 1.])
# freq_iv = np.array([1875, 3000, 667, 353, 600, 857, 1200, 250, 242])
# freq_iv = np.array([100., 143., 217., 353., 545., 857., 3000.])
# snuaddr: 'data_files/filtered_snu_cib_15_new.fits'

# ell = np.linspace(150., 2000., 20)
redshifts = np.loadtxt('data_files/redshifts.txt')

zsource = 2.
z1 = np.linspace(min(redshifts), zsource, 10)
# z2 = np.linspace(min(redshifts), 1.5, 80)
# z3 = np.linspace(1.51, max(redshifts), 30)
# z11 = np.concatenate((z2, z3))
# zn = np.linspace(min(redshifts), 3., 130)
z = redshifts  # z1  # redshifts # zn # z11

logmass = np.arange(6, 15.005, 0.1)
mass = 10**logmass

do_powerspec = 0

driver_uni = cosmo_var_iv(mass, z, do_powerspec)
driver = data_var_iv(exp)  # , z)  # , ell)

cibmean = I_nu_cib(driver, driver_uni)
# I_nu = cibmean.Iv()(nu0)
# alpha = cibmean.alpha()


# freq = ['100', '143', '217', '353', '545', '857', '3000']
# Iv_cen = np.array([13.63, 12.61, 1.64, 0.46, 2.8, 6.6, 10.1, 0.08, 0.05])
# for i in range(len(I_nu)):
#     print "Intensity is %f nW/m^2/sr at %s GHz" % (I_nu[i], freq[i])
# print (I_nu*nu0*ghz*nW/w_jy)


def plot_Inu_freq(exp, mass, zs):
    # driver_uni = cosmo_var_iv(mass, z)
    nuarray = exp['nu0']
    driver = data_var_iv(exp)  # , z)  # , ell)
    nz = len(zs)

    fig = plt.figure(figsize=(10.5, 7))
    ax = fig.add_subplot(111)

    for i_z in range(nz):
        z1 = np.linspace(min(redshifts), zs[i_z], 10)
        z = z1  # z1  # redshifts # zn # z11
        driver_uni = cosmo_var_iv(mass, z, do_powerspec)
        cibmean = I_nu_cib(driver, driver_uni)
        I_nu = cibmean.Iv()(nuarray)
        # alpha = cibmean.alpha()
        col = plt.cm.rainbow(i_z/float(nz))
        ax.plot(nuarray, I_nu, c=col, label=r'$z_s = %s$' % (zs[i_z]))
        # ax.plot(np.log(nuarray), np.log(I_nu), c=col, label=r'$z_s = %s$' % (zs[i_z]))

        Inu_mod = mod_blackbody(ghz*nuarray, zs[i_z])
        ax.plot(nuarray, Inu_mod, c=col, ls='--')

        Inu_mod2 = Theta(ghz*nuarray, zs[i_z])
        ax.plot(nuarray, Inu_mod2, c=col, ls='-.')

    ax.legend(fontsize='18', loc='lower left', frameon=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\nu_{\rm obs}$ [GHz]', fontsize=24)
    ax.set_ylabel(r'$I_\nu$ [Jy]', fontsize=24)
    # ax.set_ylim((1e-9))  # , 4.e-6))
    # ax.set_xlim((5., 4.e3))
    # ax.set_title(r'$z_{\rm source} = %s$' % (zsource))
    ax.tick_params(axis='both', labelsize=20)
    plt.show()


zs = np.array([0.5, 1.5, 3., 5.])
# plot_Inu_freq(exp, mass, zs)


def plot_alpha_freq_z(exp, mass, zs):
    fig = plt.figure(figsize=(10.5, 7))
    ax = fig.add_subplot(111)
    if exp['name'] == 'Planck_only':
        freq = ['100', '143', '217', '353', '545', '857']  # , '3000']
        nuarray = np.array([100., 143., 217., 353., 545., 857.])  # , 3000.])
        nf = len(freq)
        nz = len(zs)
        alpha = np.zeros((nf, nz))
    
        # lines = ["-", "--", "-."]
        # cl = ["g", "b", "r", "k"]
        for i_z in range(nz):
            z1 = np.linspace(min(redshifts), zs[i_z], 10)
            z = z1  # z1  # redshifts # zn # z11
            driver_uni = cosmo_var_iv(mass, z, do_powerspec)
            for i_f in range(nf):
                Planck['nucen'] = freq[i_f]
                exp_n = Planck
                driver = data_var_iv(exp_n)  # , mass)
                cibmean = I_nu_cib(driver, driver_uni)
                alpha[i_f, i_z] = cibmean.alpha()
                # alpha[:, i_z] = cibmean.alpha()
    
            # ax.plot(nuarray, alpha[:, i_z], cl[i_z], label='z = %s' % (zs[i_z]))
            col = plt.cm.rainbow(i_z/float(nz))
            ax.plot(nuarray, alpha[:, i_z], c=col, label=r'$z_s = %s$' % (zs[i_z]))
    else:
        driver = data_var_iv(exp)  # , mass)
        nuarray = exp['nu0']
        # nuarray = np.array([100., 143., 217., 353., 545., 857.])  # , 1000., 1500., 2000., 2500., 3000.])
        # nuarray = np.array([220., 280., 350., 410., 850.])  # in GHz
        nz = len(zs)
        for i_z in range(nz):
            z1 = np.linspace(min(redshifts), zs[i_z], 50)
            """
            if i_z == 0:
                # z1 = np.linspace(0.44, zs[i_z], 15)
                # z1 = np.loadtxt('data/cmass_redshift.txt')
                z1 = np.linspace(zs[i_z]-0.25, zs[i_z]+0.25, 50)
            elif i_z == 1:
                # z1 = np.linspace(0.6, zs[i_z], 50)
                # data = np.loadtxt('data/dndz_DESI_ELG.txt')
                # Z1, Z2 = data[:, 0], data[:, 1]
                # reds = (Z1+Z2)/2.
                # z1 = reds[5:17]
                z1 = np.linspace(zs[i_z]-0.25, zs[i_z]+0.25, 50)
            # """
            z = z1  # z1  # redshifts # zn # z11
            driver_uni = cosmo_var_iv(mass, z, do_powerspec)
            cibmean = I_nu_cib(driver, driver_uni)
            # alpha = cibmean.alpha()
            alpha = cibmean.alpha()(nuarray)
            # print (alpha[50], alpha[50+43], alpha[50+117], alpha[495], alpha[807])
            # print ("alpha: ", alpha)
            col = plt.cm.rainbow(i_z/float(nz))
            ax.plot(nuarray, alpha, c=col, label=r'$z_s < %s$' % (zs[i_z]))

            # alpha_mod_black = cibmean.alpha_modblackbod(ghz*nuarray, zs[i_z])
            # print ("at z = %s alpha = %s" % (zs[i_z], alpha_mod_black))
            # ax.plot(nuarray, alpha_mod_black, c=col, ls='--')

            alpha_mod_black2 = cibmean.alpha_modblackbod_inu()(nuarray)
            ax.plot(nuarray, alpha_mod_black2, c=col, ls='-.')
            # print (alpha_mod_black[40], alpha_mod_black[100], alpha_mod_black[150])
        ax.plot([], [], c='k', ls='--', label="Modified blackbody")
    # ax.plot(nuarray, alpha, 'r')
    ax.legend(fontsize='18', frameon=False)  # , labelspacing=0.1)
    # ax.set_xscale('log')
    # ax.set_yscale('log', nonposy='mask')
    ax.set_xlabel(r'$\nu_{\rm obs}$ [GHz]', fontsize=24)
    ax.set_ylabel(r'$\alpha_{\nu_0}$', fontsize=24)
    # ax.set_ylim((1.0))  # , 4.e-6))
    ax.set_xlim((50.))  # , 4.e3))
    ax.tick_params(axis='both', labelsize=20)
    plt.show()
    # plt.savefig('output/Figures/alpha-nu-zs-nonfilteredInugradient-SEDtemplate_modblackbod.pdf', bbox_inches="tight")


# zs = np.array([0.69, 1.6])  # highest redshift of CMASS and DESI EL where you
# can still detect significant number of galaxies
zs = np.array([0.5, 1.5, 3., 5.])  # , 10.])
# plot_alpha_freq_z(exp, mass, zs)


def plot_Inu_deltaInu_freq(exp, mass, zs):
    # driver_uni = cosmo_var_iv(mass, z)
    nuarray = exp['nu0']
    driver = data_var_iv(exp)  # , z)  # , ell)
    nz = len(zs)

    fig = plt.figure(figsize=(10.5, 7))
    ax = fig.add_subplot(111)

    for i_z in range(nz):
        z1 = np.linspace(min(redshifts), zs[i_z], 10)
        z = z1  # z1  # redshifts # zn # z11
        driver_uni = cosmo_var_iv(mass, z, do_powerspec)
        cibmean = I_nu_cib(driver, driver_uni)
        I_nu = cibmean.Iv()(nuarray)
        deltaInu = cibmean.delta_Inu(zs[i_z])(nuarray)
        indneg = np.where((deltaInu < 0.))
        indpos = np.where((deltaInu > 0.))
        # print (indpos)
        # print (indneg)
        ipos1 = np.where(indpos[0] < np.min(indneg[0]))[0]
        # print (ipos1)
        ipos2 = np.where(indpos[0] > np.max(indneg[0]))[0]
        # alpha = cibmean.alpha()
        col = plt.cm.rainbow(i_z/float(nz))
        ax.plot(nuarray, I_nu, c=col, ls='-', label=r'$z_s = %s$' % (zs[i_z]))
        # ax.plot(nuarray[indneg], np.abs(deltaInu[indneg]), c=col, ls='-.')
        # ax.plot(nuarray[ipos1], deltaInu[ipos1], c=col, ls='--')  # , label=r'$z_s = %s$' % (zs[i_z]))
        # ax.plot(nuarray[ipos2], deltaInu[ipos2], c=col, ls='--')
        # ax.scatter(nuarray[indpos], deltaInu[indpos], c=col, label=r'$z_s = %s$' % (zs[i_z]))

        # Inu_mod = mod_blackbody(ghz*nuarray, zs[i_z])
        # ax.plot(nuarray, Inu_mod, c=col, ls='--')

        # Inu_mod2 = Theta(ghz*nuarray, zs[i_z])
        # ax.plot(nuarray, Inu_mod2, c=col, ls='-.')

    ax.legend(fontsize='18', frameon=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\nu_{\rm obs}$ [GHz]', fontsize=24)
    ax.set_ylabel(r'$I_\nu$ [Jy/sr]', fontsize=24)
    # ax.set_ylim((1e-9))  # , 4.e-6))
    # ax.set_xlim((5., 4.e3))
    # ax.set_title(r'$z_{\rm source} = %s$' % (zsource))
    ax.tick_params(axis='both', labelsize=20)
    plt.show()


def plot_Inu_deltaInu_freq_singlegal(exp, m_d, zs, thetagal):
    # driver_uni = cosmo_var_iv(mass, z)
    # solid angle = A/chi^2=pi*(chi*theta)**2/chi**2=pi*theta**2
    thetagal /= 60.  # in degrees
    thetagal *= np.pi/180.  # radians
    Omega = np.pi*thetagal**2

    nuarray = exp['nu0']
    driver = data_var_iv(exp)  # , z)  # , ell)
    nz = len(zs)
    nv = len(nuarray)

    fig = plt.figure(figsize=(10.5, 7))
    ax = fig.add_subplot(111)

    for i_z in range(nz):
        z1 = np.linspace(min(redshifts), zs[i_z], 10)
        z = z1  # z1  # redshifts # zn # z11
        # z = np.array([zs[i_z]])  # z1  # redshifts # zn # z11
        driver_uni = cosmo_var_iv(mass, z, do_powerspec)
        cibmean = I_nu_cib(driver, driver_uni)

        # I_nu1 = cibmean.Iv()(nuarray)
        """
        pool = Pool(ncpus=4)

        def f1(nuarr):
            # nvv = len(nuarr)
            # I_nuu = np.zeros(nvv)
            # deltaInuu = np.zeros(nvv)
            # for n_nu in range(nuarr):
                # I_nuu[n_nu], deltaInuu[n_nu] = cibmean.Iv_deltaIv_modblack_singlegal(thetagal, nuarr[n_nu]*1.e9, m_d, zs[i_z])
            I_nuu, deltaInuu = cibmean.Iv_deltaIv_modblack_singlegal(thetagal, nuarr*1.e9, m_d, zs[i_z])
            return I_nuu, deltaInuu

        Inuu = []
        deltainuu = []
        for i_nu, delta_inu in pool.map(f1, nuarray):
            Inuu.append(i_nu)
            deltainuu.append(delta_inu)
        # I_nu, deltaInu = np.array(pool.map(f1, nuarray))

        I_nu = np.asarray(Inuu)[:, 0]
        deltaInu = np.asarray(deltainuu)[:, 0]
        """
        I_nu = np.zeros(nv)
        deltaInu = np.zeros(nv)
        for n_nu in range(nv):
            I_nu[n_nu], deltaInu[n_nu] = cibmean.Iv_deltaIv_modblack_singlegal(thetagal, nuarray[n_nu]*1.e9, m_d, zs[i_z])
        # """

        # print (deltaInu[:, 0])
        indneg = np.where((deltaInu < 0.))
        indpos = np.where((deltaInu > 0.))
        # print (indpos)
        # print (indneg)
        # ipos1 = np.where(indpos[0] < np.min(indneg[0]))[0]
        # print (ipos1)
        # ipos2 = np.where(indpos[0] > np.max(indneg[0]))[0]
        # alpha = cibmean.alpha()
        Inuomega = I_nu*Omega
        deltaInuomega = deltaInu*Omega

        col = plt.cm.rainbow(i_z/float(nz))
        ax.plot(nuarray, Inuomega, c=col, ls='-', label=r'$z_s = %s$' % (zs[i_z]))
        ax.plot(nuarray[indneg], np.abs(deltaInuomega[indneg]), c=col, ls='-.')
        ax.plot(nuarray[indpos], deltaInuomega[indpos], c=col, ls='--')

        # ax.plot(nuarray, I_nu1, c=col, ls='--', label=r'$z_s = %s$' % (zs[i_z]))
        # ax.plot(nuarray, I_nu, c=col, ls='-')

    # ax.plot([], [], c='k', ls='-', label="Eq. A21")

    ax.legend(fontsize='18', frameon=False)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\nu_{\rm obs}$ [GHz]', fontsize=24)
    ax.set_ylabel(r'$I_\nu \: \Omega$ [Jy]', fontsize=24)
    # ax.set_ylim((1e-9))  # , 4.e-6))
    # ax.set_xlim((5., 4.e3))
    # ax.set_title(r'$z_{\rm source} = %s$' % (zsource))
    ax.tick_params(axis='both', labelsize=20)
    plt.show()
    # plt.savefig('output/Figures/Inu_deltaInu_Doppler_cib_singlegal.pdf', bbox_inches="tight")


zs = np.array([0.5, 1.5, 3., 5.])
# plot_Inu_deltaInu_freq(exp, mass, zs)
M_D = np.array([1.e13])
# z_D = 0.5
thetagal = 0.5  # cmass galaxies subtend ~0.5 arcmin angle (Fig.11 of 2009.05557)
plot_Inu_deltaInu_freq_singlegal(exp, M_D, zs, thetagal)

print(time.time()-time0)
