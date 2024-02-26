from headers_constants import *
from CIB_halo import *
from Gal_halo import *
from CIBxGal_halo import *
from input_var_cibmean import *
from cosmo_related import *
# from Inu_cib import *
from unit_conversions import factinterp
import time

t00 = time.time()

logmass = np.arange(7, 16.005, 0.1)
# logmass = np.arange(10, 16.005, 0.1)
mass = 10**logmass

t0 = time.time()

# ################## CIB term only #######################

nuarray = np.array([100., 143., 217., 353., 545., 857.])  # , 1500., 2500.])
deltanu = 800
nucen = 100
"""
CAN SPPED UP THE COMPUTATION BY MAKING NUCEN AS AN ARRAY. IT JUST GIVES THE
CENTRAL FREQUENCY OF THE FILTER. SO WE CAN CALCULATE ALPHA FOR ALL THE
PLANCK FREQUENCIES TOGETHER IN NEED BE.
"""
nu0min = 50.  # nucen-deltanu/2.
nu0max = 3000.  # nucen+deltanu/2.
steps = nu0max-nu0min+1  # nu0max-nu0min+1  # 200
# nu0 = nuarray  # nuarray  # np.linspace(nu0min, nu0max, 200)


# color corrections for 100, 143, 217, 353, 545, 857 and 3000 GHz for Planck
# cc_pl = np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995])  # , 0.960])

cib_exp = 'CCAT'  # CCAT  # Planck
gal_exp = 'DESI_ELG'  # CMASS  # DESI_ELG  # DESI_LRG  # ExtDESI_ELG

# specs for CCAT-Prime from https://cornell.app.box.com/s/iag34277e8j0eywu8qmc0a8x3z4ko1i3/file/801743165445
# CCAT-prime Collaboration, J.Low Temp.Phys., 199, 1089 (2020), Table 1.
# freq = [27, 39, 93, 145, 225, 280, 350, 410, 850] # in GHz
# 27-145 GHz correspond to Simon's Observatory. Rest CCAT
# fwhm = [91, 63, 30, 17, 11, 9, 0.58, 0.5, 0.23] # in arcmin
# sensitivity = [35, 21, 2.6, 3.3, 6.3, 16, 105, 372, 5.7e5]# in microK arcmin

"""
How to go from Jy^2/sr to uK-arcmin:
See procedure in unit_conversion.py
Calculate MJy_KCMB from factinterp in unit_conversion.py for a given frequency
Now let us say we have x Jy^2/sr
x Jy^2/sr = x*1.e-12/(MJy_KCMB^2) KCMB^2 where 1.e-12 is for Jy^2 to MJy^2
Now once you have it in KCMB^2 unit, take squre root i.e.
sqrt(x*1.e-12/(MJy_KCMB^2)) => KCMB unit
Now multiply by 1e6 to get in muKCMB.
Then divide by (np.pi/180.) to get in degrees and multiply by 60 to get in
arcmin. So final units are muKCMB-arcmin. Example:
for 220 GHz => factinterp(220.*1e9) = 481.70441667
If CCAT has white noise = sensitivity^2 = 4.2 Jy^2/sr that comes out to be
np.sqrt(4.2*1.e-12/481.70441667**2)*1e6*60./(np.pi/180.) = 14.63 muKCMB-arcmin
"""

if cib_exp == 'CCAT':
    # ell = np.linspace(50, 40000, 200)
    # 1908.10451v2 and Fiona's paper: 2010.16405
    nu0 = np.array([220., 280., 350., 410., 850.])  # in GHz
    fwhm = np.array([0.95, 0.75, 0.58, 0.5, 0.23])  # in arcmin
    fwhm_rad = fwhm*(np.pi/180.)/60.
    # fwhm = 0.5  # 4.8 Planck arcmin # 0.5 CCAT-Prime => previous
    sensitivity = np.array([15., 27., 105., 372., 5.7e5])  # in microK arcmin
    sensitivity_Jy2_sr = np.array([4.2, 11.8, 85.1, 468., 69483.])
    sensitivity_Jy_sr = np.sqrt(sensitivity_Jy2_sr)
    # sensitivity = 1.2  # 13.5 Planck Jy/sr  # 1.2 CCAT-Prime => previous
    # ell = np.linspace(100, 50000, 200)
    # ell = np.linspace(100, 50000, 400)
    ell = np.logspace(2, np.log10(5e4), num=400)
elif cib_exp == 'Planck':
    # ell = np.linspace(50, 3000, 15)
    nu0 = np.array([100., 143., 217., 353., 545., 857.])  # , 1000., 1500., 2000., 2500., 3000.])
    # fwhm obtained from https://wiki.cosmos.esa.int/planckpla/index.php/Effective_Beams
    fwhm = np.array([9.651, 7.248, 4.990, 4.818, 4.682, 4.325])  # in arcmin
    fwhm_rad = fwhm*(np.pi/180.)/60.
    # fwhm = 4.8  # 4.8 Planck arcmin # 0.5 CCAT-Prime
    # 1303.5067 for sensitivity (Tab:4 muK deg & kJy/sr deg)
    sensitivity_muKarcmin = np.array([108., 48., 60., 210., 1137., 29075.])
    sensitivity_Jy2_sr = np.array([58., 26.929, 72., 305, 369., 369.])
    sensitivity_Jy_sr = np.sqrt(sensitivity_Jy2_sr)
    # sensitivity = 13.5  # 13.5 Planck Jy/sr  # 1.2 CCAT-Prime
    ell = np.linspace(100, 5000, 99)

cc = np.ones(len(nu0))
fc = np.ones(len(nu0))

strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nolens_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"
# strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"

cibres = "data_files/one_halo_bestfit_"+strfig
# clres = np.loadtxt('data/%s.txt' % (string))

"""
Planck = {'name': 'Planck_only',
          'do_cibmean': 0,
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
"""

custom = {'name': cib_exp,
          'do_cibmean': 0,
          'cc': cc,
          'fc': fc,
          'snuaddr': 'data_files/filtered_snu_planck.fits',
          'nu0min': nu0min, 'nu0max': nu0max,
          'nucen': str(int(nucen)),
          'nu0': nu0,
          'ell': ell,
          'fwhm': fwhm,
          'sensitivity': sensitivity_Jy_sr,
          'cibpar_resfile': cibres}

exp = custom

redshifts = np.loadtxt('data_files/redshifts.txt')

if gal_exp == 'CMASS':
    reds = np.loadtxt('data/cmass_redshift.txt')
elif gal_exp == 'DESI_ELG' or gal_exp == 'DESI_LRG':
    data = np.loadtxt('data/dndz_'+gal_exp+'.txt')
    Z1, Z2 = data[:, 0], data[:, 1]
    reds = (Z1+Z2)/2.
else:
    data = np.loadtxt('data/dndz_DESI_ELG.txt')
    Z1, Z2 = data[:, 0], data[:, 1]
    reds = (Z1+Z2)/2.
    reds *= 2.  # just as a test to see if we expand the redshift
    # distribution of DESI ELG galaxies out to redshift 0 < z < 4 instead of
    # 0 < z < 2 which is currently in place. With the same value of
    # N_elg.

z = redshifts  # redshifts  # reds
do_powerspec = 1

driver_uni = cosmo_var_iv(mass, z, do_powerspec)
# driver_uni.plot_beta2()
driver = data_var_iv(exp)  # , z)  # , ell)

instcib = Cib_halo(driver, driver_uni)

"""
cl_cibtot = instcib.cl_cibtot()
cl_cibDopplertot = instcib.cl_cibDoppler_tot()


def plot_cib(nu1, nu2, cibinstance, ellinput):
    instcib = cibinstance
    cl1h = instcib.onehalo_int()
    cl2h = instcib.twohalo_int()
    clshot = instcib.shot_cib()
    cl_cibtot = cl1h+cl2h+clshot

    cl1h_Doppler = instcib.onehalo_int_Doppler()
    cl2h_Doppler = instcib.twohalo_int_Doppler()
    cltot_Doppler = cl1h_Doppler+cl2h_Doppler

    f, ax = plt.subplots(figsize=(10, 10))
    ax.plot(ellinput, cl1h[nu1, nu2, :], 'r--', label='1h-CIB')
    ax.plot(ellinput, cl2h[nu1, nu2, :], 'b--', label='2h-CIB')
    ax.plot(ellinput, clshot[nu1, nu2, :], 'g--', label='shot-CIB')
    ax.plot(ellinput, cl_cibtot[nu1, nu2, :], 'k--', label='tot-CIB')
    ax.plot(ellinput, cl1h_Doppler[nu1, nu2, :], 'r', label='1h-Doppler')
    ax.plot(ellinput, cl2h_Doppler[nu1, nu2, :], 'b', label='2h-Doppler')
    ax.plot(ellinput, cltot_Doppler[nu1, nu2, :], 'k', label='tot-Doppler')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper right', prop={'size': 14}, frameon=False)
    ax.set_ylabel(r'$C_\ell\: [Jy^2\: {sr}^{-1}]$', fontsize=14)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=14)


def plot_foregrounds(nu1, nu2, cibinstance, ellinput, elloutput):
    instcib = cibinstance
    clcib1h = interp1d(ellinput, instcib.onehalo_int()[nu1, nu2, :],
                       kind='linear', bounds_error=False, fill_value="extrapolate")
    clcib2h = interp1d(ellinput, instcib.twohalo_int()[nu1, nu2, :],
                       kind='linear', bounds_error=False, fill_value="extrapolate")
    clcibshot = interp1d(ellinput, instcib.shot_cib()[nu1, nu2, :],
                         kind='linear', bounds_error=False, fill_value="extrapolate")

    cl1h_Doppler = interp1d(ellinput, instcib.onehalo_int_Doppler()[nu1, nu2, :],
                            kind='linear', bounds_error=False, fill_value="extrapolate")
    cl2h_Doppler = interp1d(ellinput, instcib.twohalo_int_Doppler()[nu1, nu2, :],
                            kind='linear', bounds_error=False, fill_value="extrapolate")

    cmbdata = np.loadtxt('data/Manu/lensedCls.dat')  # lensed cmb
    dltt = interp1d(cmbdata[:, 0], cmbdata[:, 1], kind='linear',
                    bounds_error=False, fill_value="extrapolate")

    fg_loc = 'data/foreground_power_spectra/'

    tszres = np.loadtxt(fg_loc+'tsz/cltsz_planck_%sGHzx%sGHz_B1p41_Jy2.txt' % (str(int(custom['nu0'][nu1])), str(int(custom['nu0'][nu2]))))
    cltsz = interp1d(tszres[:, 0], tszres[:, -1], kind='linear',
                     bounds_error=False, fill_value="extrapolate")

    cibxtszres = np.loadtxt(fg_loc+'cibxtsz/clcibxtsz_planck_%sGHzx%sGHz_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_B1p41_Jy2.txt' % (str(int(custom['nu0'][nu1])), str(int(custom['nu0'][nu2]))))
    clcibxtsz = interp1d(cibxtszres[:, 0], cibxtszres[:, -1], kind='linear',
                         bounds_error=False, fill_value="extrapolate")

    MJy_KCMB2 = factinterp(custom['nu0'][nu1]*1.e9)*factinterp(custom['nu0'][nu2]*1.e9)

    dl = lambda l: l*(l+1.)/(2.*np.pi)
    el = elloutput
    f, ax = plt.subplots(figsize=(10, 10))

    clcibtot = clcib1h(el)+clcib2h(el)+clcibshot(el)
    cltot_Doppler = cl1h_Doppler(el)+cl2h_Doppler(el)

    ax.plot(el, clcibtot*dl(el)/MJy_KCMB2, 'r', label='CIB')
    ax.plot(el, cltsz(el)*dl(el)/MJy_KCMB2, 'b', label='tSZ')
    ax.plot(el, clcibxtsz(el)*dl(el)/MJy_KCMB2, 'g', label='CIB x tSZ')

    ax.plot(el, cltot_Doppler*dl(el)/MJy_KCMB2, 'r--', label='Doppler CIB')

    ax.plot(el, dltt(el), 'k', label='Lensed CMB')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper right', prop={'size': 14}, frameon=False)
    # ax.set_ylabel(r'$C_\ell\: [Jy^2\: {sr}^{-1}]$', fontsize=14)
    ax.set_ylabel(r'$D_\ell\: [\mu {\rm K}^2]$', fontsize=14)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=14)


# [100., 143., 217., 353., 545., 857.]
nu1, nu2 = 1, 1
# plot_cib(nu1, nu2, instcib, exp['ell'])
# ellout = np.linspace(2000, 5000., 3001)
# plot_foregrounds(nu1, nu2, instcib, exp['ell'], ellout)

# MJy_KCMB = factinterp(nu0*1.e9)
cl_100, cl_143, cl_217 = cl_cibtot[0, 0], cl_cibtot[1, 1], cl_cibtot[2, 2]
cl_100_muK2 = cl_100/factinterp(100*1.e9)**2
cl_143_muK2 = cl_143/factinterp(143*1.e9)**2
cl_217_muK2 = cl_217/factinterp(217*1.e9)**2

# cl_cibDoppler_muK2 = cl_cibDopplertot*MJy_KCMB**2
cl_dp_100, cl_dp_143, cl_dp_217 = cl_cibDopplertot[0, 0], cl_cibDopplertot[1, 1], cl_cibDopplertot[2, 2]
cldp_100_muK2 = cl_dp_100/factinterp(100*1.e9)**2
cldp_143_muK2 = cl_dp_143/factinterp(143*1.e9)**2
cldp_217_muK2 = cl_dp_217/factinterp(217*1.e9)**2

i1 = max(np.where(ell < 3000)[0])
i2 = min(np.where(ell >= 3000)[0])
l1 = ell[i1]
l2 = ell[i2]
d1 = np.abs(3000 - l1)
d2 = np.abs(3000 - l2)
if d1 < d2:
    ind3000 = i1
else:
    ind3000 = i2

print ("")
print "Doppler power spectrum in muK2 at 100, 143, and 217 GHz for Planck at l=%s is %s, %s, and %s:" % (ell[ind3000], cldp_100_muK2[ind3000], cldp_143_muK2[ind3000], cldp_217_muK2[ind3000])
print "CIB power spectrum in muK2 at 100, 143, and 217 GHz for Planck at l=%s is %s, %s, and %s:" % (ell[ind3000], cl_100_muK2[ind3000], cl_143_muK2[ind3000], cl_217_muK2[ind3000])

dl = ell*(ell+1)/(2*np.pi)

dl_100_muK2 = cl_100_muK2*dl
dl_143_muK2 = cl_143_muK2*dl
dl_217_muK2 = cl_217_muK2*dl

dldp_100_muK2 = cldp_100_muK2*dl
dldp_143_muK2 = cldp_143_muK2*dl
dldp_217_muK2 = cldp_217_muK2*dl

print ("")
print "Doppler power spectrum Dl in muK2 at 100, 143, and 217 GHz for Planck at l=%s is %s, %s, and %s" % (ell[ind3000], dldp_100_muK2[ind3000], dldp_143_muK2[ind3000], dldp_217_muK2[ind3000])
print "CIB power spectrum Dl in muK2 at 100, 143, and 217 GHz for Planck at l=%s is %s, %s, and %s" % (ell[ind3000], dl_100_muK2[ind3000], dl_143_muK2[ind3000], dl_217_muK2[ind3000])
# """
# djc = instcib.djc_dlnMh()
# print (instcib.alpha()(nu0))
# print (instcib.deltaIdoppler())

instgal = ProfHODMore15(driver, driver_uni, gal_exp)
# cl_galtot = instgal.cl_galtot()
# cl_veltot = instgal.cl_veltot()

"""
def plot_gal(galinstance, ellinput):
    instgal = galinstance
    cl1h = instgal.cl1h_gal()
    cl2h = instgal.cl2h_gal()
    shot = instgal.clshot_gal()
    clshot = np.repeat(shot, len(ellinput))
    cl_galtot = cl1h+cl2h+clshot

    f, ax = plt.subplots(figsize=(10, 10))
    ax.plot(ellinput, cl1h, 'r--', label='1h-gal')
    ax.plot(ellinput, cl2h, 'b--', label='2h-gal')
    ax.plot(ellinput, clshot, 'g--', label='shot')
    ax.plot(ellinput, cl_galtot, 'k--', label='tot-gal')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper right', prop={'size': 14}, frameon=False)
    ax.set_ylabel(r'$C_\ell$', fontsize=14)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=14)


plot_gal(instgal, exp['ell'])
# """

z = reds  # redshifts  # reds
# z = np.linspace(min(reds), max(reds), 50)
driver_uni = cosmo_var_iv(mass, z, do_powerspec)
# instgalz = ProfHODMore15(driver, driver_uni, gal_exp)
# print (instgalz.nbargal())
# print (instgalz.N_tot())
"""
r_l is the correlation coefficient between CIB and galaxies. It is being
used to calculate the cross shot nise between CIB and galaxies, and CIB and
velocities. Set r_l = 0 if you do not want cross-shot noise
"""
r_l = 1.0
# instcross = CIBxgal(instgal, instcib)
instcross = CIBxgal(driver, driver_uni, gal_exp, r_l)
# print (instcross.alpha()(nu0))
# print (instcross.deltaIdoppler())
# cl_crosstot = instcross.cibgalcross_cell_tot()

"""
def plot_cross(nu1, crossinstance, ellinput):
    inst = crossinstance
    cl1h = inst.cibgalcross_cell_1h()
    cl2h = inst.cibgalcross_cell_2h()
    clshot = inst.cibgalcross_cell_shot()
    # clshot = np.repeat(shot, len(ellinput))
    cl_crosstot = cl1h+cl2h+clshot

    f, ax = plt.subplots(figsize=(10, 10))
    ax.plot(ellinput, cl1h[nu1, :], 'r--', label='1h-cross')
    ax.plot(ellinput, cl2h[nu1, :], 'b--', label='2h-cross')
    ax.plot(ellinput, clshot[nu1, :], 'g--', label='shot')
    ax.plot(ellinput, cl_crosstot[nu1, :], 'k--', label='tot-cross')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper right', prop={'size': 14}, frameon=False)
    ax.set_ylabel(r'$C_\ell\: [Jy\: {sr}^{-1}]$', fontsize=14)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=14)


# [100., 143., 217., 353., 545., 857.]
nu1 = 5
plot_cross(nu1, instcross, exp['ell'])
# """

"""
def plot_all(nu1, crossinstance, galinstance, cibinstance, ellinput):
    instd = crossinstance
    cl1hd = instd.cibvelcross_cell_1h()
    cl2hd = instd.cibvelcross_cell_2h()
    cl_crosstotd = cl1hd+cl2hd

    instg = galinstance
    cl1hv = instg.cl1h_vel()
    cl2hv = instg.cl2h_vel()
    cl_crosstotv = cl1hv+cl2hv

    cl1hg = instg.cl1h_gal()
    cl2hg = instg.cl2h_gal()
    clshotg = instg.clshot_gal()
    cl_crosstotg = cl1hg+cl2hg+clshotg

    instc = cibinstance
    cl1hc = instc.onehalo_int()
    cl2hc = instc.twohalo_int()
    clshotc = instc.shot_cib()
    # clwhite = driver.whitenoise(nu1)
    cl_crosstotc = cl1hc+cl2hc+clshotc  # +clwhite

    cl1hdb = instc.onehalo_int_Doppler()
    cl2hdb = instc.twohalo_int_Doppler()
    cldbtot = cl1hdb+cl2hdb
    
    # lines = ['-', '--', '-.']
    # cols = ['r', 'b', 'g']

    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True, figsize=(10, 18))
    # for i in range(3):
        # ax.plot(ellinput, cl1h[nu1, :], c=cols[i], ls=lines[i], label='1h-cross')
        # ax.plot(ellinput, cl2h[nu1, :], 'b--', label='2h-cross')
        # ax.plot(ellinput, clshot[nu1, :], 'g--', label='shot')
    # cibfac = 1e-4
    # galvelfac = 1e4
    # ax1.plot(ellinput, cibfac*cl_crosstotc[nu1, nu1, :], 'r', label=r'$10^{-4} \times C_\ell^{\Delta \rm{I} \Delta \rm{I}}$')

    ax1.plot(ellinput, cl2hc[nu1, nu1, :], 'b--', label='2-halo')  #, r'$C_\ell^{\Delta \rm{I} \Delta \rm{I}}$')
    ax1.plot(ellinput, cl1hc[nu1, nu1, :], 'b-.', label='1-halo')
    # ax1.plot(ellinput, clwhite, 'b:', label='Detector noise')# , r'$C_\ell^{\Delta \rm{I^D} q_\gamma}$')
    ax1.plot(ellinput, cl_crosstotc[nu1, nu1, :], 'b-', label='total')  # , r'$C_\ell^{q_\gamma q_\gamma}$')
    ax1.set_ylim(ymin=1e2, ymax=1e6)
    ax1.set_ylabel(r'$C_\ell^{\Delta I \Delta I}\: [(\rm{Jy/sr})^{2}]$', fontsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend(loc='upper right', prop={'size': 18}, frameon=False, ncol=3, labelspacing=0.2)
    ax1.tick_params(axis='both', labelsize=18)

    ax2.plot(ellinput, cl2hd[nu1, :], 'b--') # , label=r'$C_\ell^{\Delta \rm{I^D} q_\gamma}$')
    ax2.plot(ellinput, cl1hd[nu1, :], 'b-.')
    ax2.plot(ellinput, cl_crosstotd[nu1, :], 'b-')
    ax2.set_ylim(ymin=1e-9)
    ax2.set_ylabel(r'$C_\ell^{\Delta I^{\rm DB} q_\gamma}\: [\rm{Jy}\: \rm{sr}^{-1}]$', fontsize=16)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.tick_params(axis='both', labelsize=18)

    ax3.plot(ellinput, cl2hv, 'b--') # , label=r'$C_\ell^{\Delta \rm{I^D} q_\gamma}$')
    ax3.plot(ellinput, cl1hv, 'b-.')
    ax3.plot(ellinput, cl_crosstotv, 'b-')
    ax3.set_ylim(ymin=1e-17, ymax=1.5e-12)
    ax3.set_ylabel(r'$C_\ell^{q_\gamma q_\gamma}$', fontsize=16)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.tick_params(axis='both', labelsize=18)

    ax4.plot(ellinput, cl1hdb[nu1, nu1, :], 'b--') # , label=r'$C_\ell^{\Delta \rm{I^D} q_\gamma}$')
    ax4.plot(ellinput, cl2hdb[nu1, nu1, :], 'b-.')
    ax4.plot(ellinput, cldbtot[nu1, nu1, :], 'b-')
    ax4.set_ylim(ymin=1e-4, ymax=1)
    ax4.set_xlim(xmin=1e2)
    ax4.set_ylabel(r'$C_\ell^{\Delta I^{\rm DB} \Delta I^{\rm DB}}\: [(\rm{Jy/sr})^{2}]$', fontsize=16)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.tick_params(axis='both', labelsize=18)

    ax5.plot(ellinput, cl2hg, 'b--') # , label=r'$C_\ell^{\Delta \rm{I^D} q_\gamma}$')
    ax5.plot(ellinput, cl1hg, 'b-.')
    ax5.plot(ellinput, cl_crosstotg, 'b-')
    ax5.set_ylim(ymin=1e-10)
    ax5.set_xlim(xmin=1e2)
    ax5.set_ylabel(r'$C_\ell^{\delta_g \delta_g}$', fontsize=16)
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    ax5.tick_params(axis='both', labelsize=18)

    # ax.legend(loc='upper right', prop={'size': 14}, frameon=False)
    # ax.set_ylabel(r'$C_\ell\: [Jy\: {sr}^{-1}]$', fontsize=14)
    ax5.set_xlabel('Multipole' r'$\;\ell$', fontsize=18)
    ax1.set_title(r'CCAT-Prime 850 GHz & DESI-ELG', fontsize=18)
    plt.savefig('output/Figures/allpowerspectra_CCAT850_DESIELG.pdf', bbox_inches="tight")

# [100., 143., 217., 353., 545., 857.] => Planck
# [220., 280., 350., 410., 850.] => CCAT Prime
nu1 = 4
# plot_all(nu1, instcross, instgal, instcib, exp['ell'])
# """

#"""

def plot_cibxvel(nu1, crossinstance, galinstance, cibinstance, lmin, lmax, deltal, fsky):
    instd = crossinstance
    cl1hd = instd.cibvelcross_cell_1h()
    cl2hd = instd.cibvelcross_cell_2h()
    cl_crosstotd = cl1hd+cl2hd
    # fclcrossd = UnivariateSpline(exp['ell'], cl_crosstotd[nu1, :], k=1, s=0)
    fclcross1hd = interp1d(exp['ell'], cl1hd[nu1, :], kind='linear', bounds_error=False, fill_value="extrapolate")
    fclcross2hd = interp1d(exp['ell'], cl2hd[nu1, :], kind='linear', bounds_error=False, fill_value="extrapolate")
    fclcrossd = interp1d(exp['ell'], cl_crosstotd[nu1, :], kind='linear', bounds_error=False, fill_value="extrapolate")

    clwhite = driver.whitenoise(nu1)
    clcib = cibinstance.cl_cibtot()[nu1, nu1, :]
    clcibtot = clcib+clwhite
    # fclcibtot = UnivariateSpline(exp['ell'], clcibtot, k=1, s=0)
    fclcibtot = interp1d(exp['ell'], clcibtot, kind='linear', bounds_error=False, fill_value="extrapolate")

    clvel = galinstance.cl_veltot()
    # fclvel = UnivariateSpline(exp['ell'], clvel, k=1, s=0)
    fclvel = interp1d(exp['ell'], clvel, kind='linear', bounds_error=False, fill_value="extrapolate")

    numnoise = cl_crosstotd[nu1, :]**2 + clcibtot*clvel
    denmnoise = (2.*exp['ell']+1)*fsky
    noise = np.sqrt(numnoise/denmnoise)
    fnoise = interp1d(exp['ell'], noise, kind='linear', bounds_error=False, fill_value="extrapolate")

    # deltal = 100.
    ellplt = np.linspace(lmin, lmax, int(lmax-lmin+1))
    ell2 = np.arange(lmin, lmax, deltal)
    nbin = len(ell2)
    cibtbin = np.zeros(nbin)
    velbin = np.zeros(nbin)
    crossdbin = np.zeros(nbin)
    noisebin2 = np.zeros(nbin)
    lcen = np.zeros(nbin)
    for i in range(nbin):
        l1 = ell2[i]
        if i == int(nbin-1):
            l2 = lmax
        else:
            l2 = ell2[i+1]
        lcen[i] = (l1+l2)/2.
        el = np.linspace(l1, l2, int(l2-l1+1))
        # deltal = l2-l1
        cibtbin[i] = np.sum(fclcibtot(el))/deltal
        velbin[i] = np.sum(fclvel(el))/deltal
        crossdbin[i] = np.sum(fclcrossd(el))/deltal
        noisebin2[i] = np.sum(fnoise(el))/deltal
    # print (noisebin)
    # print (fclcrossd(lcen))
    # lines = ['-', '--', '-.']
    # cols = ['r', 'b', 'g']
    num = crossdbin**2 + cibtbin*velbin
    denom = (2.*lcen+1)*deltal*fsky
    noisebin = np.sqrt(num/denom)

    # print (noisebin)
    # print (noisebin2)

    f, ax = plt.subplots(figsize=(8, 5))

    # ax.plot(ellinput, cl2hd[nu1, :], 'k--', label=r'2-halo')
    # ax.plot(ellinput, cl1hd[nu1, :], 'k-.', label=r'1-halo')
    # ax.plot(exp['ell'], cl_crosstotd[nu1, :], 'k-.')
    ax.plot(ellplt, fclcrossd(ellplt), c='tab:orange', ls='-', lw=3, label=r'total')
    ax.plot(ellplt, fclcross1hd(ellplt), c='tab:orange', ls='-.', lw=3, label=r'1-halo')
    ax.plot(ellplt, fclcross2hd(ellplt), c='tab:orange', ls='--', lw=3, label=r'2-halo')
    yerrorbin = noisebin/2
    ax.errorbar(lcen, fclcrossd(lcen), yerr=yerrorbin, fmt='o', c='b', markersize='3', capsize=3, capthick=1) # fmt='o', ':'
    # yerror = fnoise(ellplt)/2.
    # ax.fill_between(ellplt, fclcrossd(ellplt)-yerror, fclcrossd(ellplt)+yerror, alpha=0.7)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(ymin=5e-9, ymax=2e-7)
    ax.set_xlim(xmin=1e2)
    ax.legend(loc='upper right', prop={'size': 24}, frameon=False)
    ax.set_ylabel(r'$C_\ell^{\Delta I^{\rm DB} q_\gamma}\: [\rm{Jy}\: \rm{sr}^{-1}]$', fontsize=24)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=24)
    ax.set_title(r'CCAT-Prime 850 GHz x DESI-ELG', fontsize=24)
    ax.tick_params(axis='both', labelsize=24)
    # plt.savefig('output/Figures/DBCIBcrossmomentum_witherrorbars_CCAT850_DESIELG_1h2h.pdf', bbox_inches="tight")

# [100., 143., 217., 353., 545., 857.] => Planck
# [220., 280., 350., 410., 850.] => CCAT Prime
nu1 = 4
lmin = 100
lmax = 5e4
fsky = 0.4
deltal = 1200
plot_cibxvel(nu1, instcross, instgal, instcib, lmin, lmax, deltal, fsky)
# """

# zs = np.array([0.05, 1.5, 3., 5.])
# nu0_z = instcross.nu0_z(zs)
"""
MJy_KCMB = factinterp(nu0*1.e9)
# multiplying by 1e-6 because units are Jy/sr and we want MJy/sr
# then dividing by MJy_KCMB factor to get the units in KCMB
# then multiplying by 1e6 to get muKCMB
# cl_velcrosstot = instcross.cibvelcross_cell_tot()
# cl_velcrosstot_KCMB = cl_velcrosstot*1.e-6/MJy_KCMB[:, None]
# cl_velcrosstot_muKCMB = cl_velcrosstot_KCMB*1.e6
# ell[118] = 3000
# deltaI_Jy = instcross.deltaIdoppler()  # deltaIdoppler_all

# mass, redshift and freq value for Doppler calculation
M_D = 1.e13
z_D = 0.5
thetagal = 0.5  # cmass galaxies subtend ~0.5 arcmin angle (Fig.11 of 2009.05557)
deltaI_Jy = instcross.deltaIdoppler(thetagal, nu0*1.e9, M_D, z_D)
deltaI_KCMB = deltaI_Jy*1.e-6/MJy_KCMB  # [:, None]
deltaI_muKCMB = deltaI_KCMB*1.e6
print ("delta I in muKCMB for %s GHz = %s" % (nu0, deltaI_muKCMB))
print ("delta I in Jy for %s GHz = %s" % (nu0, deltaI_Jy))
# """

"""
# write the code such that it calculates SNR for all frequencies? Maybe
# make it loop over all the frequencies?
# CCAT: nu0 = np.array([220., 280., 350., 410., 850.])  # in GHz
# Planck: nu0 = np.array([100., 143., 217., 353., 545., 857.])
# n_nu = 4
fsky = 0.4
for n_nu in range(len(custom['nu0'])):
    # snrg1, snrg2 = instcross.snr_g(n_nu, fsky)
    # snrg11, snrg12 = instcross.snr__g(n_nu, fsky, instcib, instgal)
    snrg12 = instcross.snr__g(n_nu, fsky, instcib, instgal)
    print ("SNR for CIB x gal at %s GHz for %s x %s is" % (str(int(custom['nu0'][n_nu])), custom['name'], gal_exp), snrg12[-1])

print ("")

for n_nu in range(len(custom['nu0'])):
    # snrv1, snrv2 = instcross.snr_v(n_nu, fsky)
    # snrv11, snrv12 = instcross.snr__v(n_nu, fsky, instcib, instgal)
    snrv12 = instcross.snr__v(n_nu, fsky, instcib, instgal)
    print ("SNR for CIB x vel at %s GHz for %s x %s is" % (str(int(custom['nu0'][n_nu])), custom['name'], gal_exp), snrv12[-1])
# """

# driver_uni.plot_beta2()

"""
def plot_alpha_freq_z(exp, mass, zs):
    fig = plt.figure(figsize=(10.5, 7))
    ax = fig.add_subplot(111)
    driver = data_var_iv(exp)  # , mass)
    nuarray = exp['nu0']
    nz = len(zs)
    for i_z in range(nz):
        z1 = np.linspace(min(redshifts), zs[i_z], 50)
        z = z1  # z1  # redshifts # zn # z11
        driver_uni = cosmo_var_iv(mass, z, do_powerspec)
        # cibmean = I_nu_cib(driver, driver_uni)
        instcib = Cib_halo(driver, driver_uni)
        # alpha = cibmean.alpha()
        col = plt.cm.rainbow(i_z/float(nz))
        alpha_mod_black = instcib.alpha_modblack()(nuarray)
        # print (nuarray.shape, alpha_mod_black.shape)
        ax.plot(nuarray, alpha_mod_black, c=col, label=r'$z_s = %s$' % (zs[i_z]))

        print ("at z = %s alpha = %s" % (zs[i_z], alpha_mod_black))

        # alpha_mod_black2 = alpha_modblackbod_Theta(zs[i_z], ghz*nuarray)
        # ax.plot(nuarray, alpha_mod_black2, c=col, ls='-.')
        # print (alpha_mod_black[40], alpha_mod_black[100], alpha_mod_black[150])
    ax.plot([], [], c='k', ls='--', label="Modified blackbody")
    # ax.plot(nuarray, alpha, 'r')
    ax.legend(fontsize='18', loc='lower left', frameon=False)  # , labelspacing=0.1)
    # ax.set_xscale('log')
    # ax.set_yscale('log', nonposy='mask')
    ax.set_xlabel(r'$\nu_{\rm obs}$ [GHz]', fontsize=24)
    ax.set_ylabel(r'$\alpha_{\nu_0}$', fontsize=24)
    # ax.set_ylim((1.0))  # , 4.e-6))
    ax.set_xlim((50.))  # , 4.e3))
    ax.tick_params(axis='both', labelsize=20)
    plt.show()
    # plt.savefig('output/Figures/alpha-nu-zs-nonfilteredInugradient-general.pdf', bbox_inches="tight")


zs = np.array([0.5, 1.5, 3., 5.])  # , 10.])
# zs = np.array([0.5])  # , 1.5, 3., 5.])
plot_alpha_freq_z(exp, mass, zs)
# """

t1 = time.time()
print t1 - t0
