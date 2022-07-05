from headers_constants import *
from CIB_halo import *
from Gal_halo import *
from CIBxGal_halo import *
from input_var_cibmean import *
from cosmo_related import *
from unit_conversions import factinterp
import time

logmass = np.arange(7, 16.005, 0.1)
# logmass = np.arange(10, 16.005, 0.1)
mass = 10**logmass

t0 = time.time()

# ################## CIB term only #######################

cib_exp = 'Planck'  # CCAT  # Planck
gal_exp = 'CMASS'  # CMASS  # DESI_ELG  # DESI_LRG  # ExtDESI_ELG

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
    # 1908.10451v2
    nu0 = np.array([220., 280., 350., 410., 850.])  # in GHz
    nus = ['220', '280', '350', '410', '850']
    fwhm = np.array([0.95, 0.75, 0.58, 0.5, 0.23])  # in arcmin
    fwhm_rad = fwhm*(np.pi/180.)/60.
    # fwhm = 0.5  # 4.8 Planck arcmin # 0.5 CCAT-Prime => previous
    sensitivity = np.array([15., 27., 105., 372., 5.7e5])  # in microK arcmin
    sensitivity_Jy2_sr = np.array([4.2, 11.8, 85.1, 468., 69483.])
    sensitivity_Jy_sr = np.sqrt(sensitivity_Jy2_sr)
    # sensitivity = 1.2  # 13.5 Planck Jy/sr  # 1.2 CCAT-Prime => previous
    ell = np.linspace(100, 50000, 200)
elif cib_exp == 'Planck':
    nu0 = np.array([100., 143., 217., 353., 545., 857.])  # , 1000., 1500., 2000., 2500., 3000.])
    nus = ['100', '143', '217', '353', '545', '857']
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

# color corrections for 100, 143, 217, 353, 545, 857 and 3000 GHz for Planck
cc_pl = np.array([1.076, 1.017, 1.119, 1.097, 1.068, 0.995])  # , 0.960])

cc = cc_pl  #  np.ones(len(nu0))
fc = np.ones(len(nu0))

strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nolens_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"
# strfig = "allcomponents_lognormal_sigevol_1p5zcutoff_nospire_fcpl_onlyautoshotpar_no3000_gaussian600n857n1200_planck_spire_hmflog10.txt"

cibres = "data_files/one_halo_bestfit_"+strfig
# clres = np.loadtxt('data_files/%s.txt' % (string))


custom = {'name': cib_exp,
          'do_cibmean': 0,
          'cc': cc,
          'fc': fc,
          'snuaddr': 'data_files/filtered_snu_planck.fits',
          # 'nu0min': nu0min, 'nu0max': nu0max,
          # 'nucen': str(int(nucen)),
          'nu0': nu0,
          'ell': ell,
          'fwhm': fwhm,
          'sensitivity': sensitivity_Jy_sr,
          'cibpar_resfile': cibres,
          'nu_string': nus}

exp = custom

redshifts = np.loadtxt('data_files/redshifts.txt')

if gal_exp == 'CMASS':
    reds = np.loadtxt('data_files/cmass_redshift.txt')
elif gal_exp == 'DESI_ELG' or gal_exp == 'DESI_LRG':
    data = np.loadtxt('data_files/dndz_'+gal_exp+'.txt')
    Z1, Z2 = data[:, 0], data[:, 1]
    reds = (Z1+Z2)/2.
else:
    data = np.loadtxt('data_files/dndz_DESI_ELG.txt')
    Z1, Z2 = data[:, 0], data[:, 1]
    reds = (Z1+Z2)/2.
    reds *= 2.  # just as a test to see if we expand the redshift
    # distribution of DESI ELG galaxies out to redshift 0 < z < 4 instead of
    # 0 < z < 2 which is currently in place. With the same value of
    # N_elg.

z = redshifts  # redshifts  # reds
do_powerspec = 1  # calculates the halo mass function and halo bias required
# for the calculation of the halo model power spectra.

driver = data_var_iv(exp)  # , z)  # , ell)

# ################# CIB power spectra ########################

"""
driver_uni = cosmo_var_iv(mass, z, do_powerspec)

instcib = Cib_halo(driver, driver_uni)

cl_cibtot = instcib.cl_cibtot()


def plot_cib(nu1, nu2, cibinstance, ellinput):
    instcib = cibinstance
    cl1h = instcib.onehalo_int()
    cl2h = instcib.twohalo_int()
    clshot = instcib.shot_cib()
    cl_cibtot = cl1h+cl2h+clshot

    f, ax = plt.subplots(figsize=(10, 10))
    ax.plot(ellinput, cl1h[nu1, nu2, :], 'r--', label='1h-CIB')
    ax.plot(ellinput, cl2h[nu1, nu2, :], 'b--', label='2h-CIB')
    ax.plot(ellinput, clshot[nu1, nu2, :], 'g--', label='shot-CIB')
    ax.plot(ellinput, cl_cibtot[nu1, nu2, :], 'k--', label='tot-CIB')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='upper right', prop={'size': 14}, frameon=False)
    ax.set_ylabel(r'$C_\ell\: [Jy^2\: {sr}^{-1}]$', fontsize=14)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=14)


# [100., 143., 217., 353., 545., 857.]
nu1, nu2 = 1, 1
# plot_cib(nu1, nu2, instcib, exp['ell'])


# MJy_KCMB = factinterp(nu0*1.e9)
cl_100, cl_143, cl_217 = cl_cibtot[0, 0], cl_cibtot[1, 1], cl_cibtot[2, 2]
cl_100_muK2 = cl_100/factinterp(100*1.e9)**2
cl_143_muK2 = cl_143/factinterp(143*1.e9)**2
cl_217_muK2 = cl_217/factinterp(217*1.e9)**2

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
print "CIB power spectrum in muK2 at 100, 143, and 217 GHz for Planck at l=%s is %s, %s, and %s:" % (ell[ind3000], cl_100_muK2[ind3000], cl_143_muK2[ind3000], cl_217_muK2[ind3000])

dl = ell*(ell+1)/(2*np.pi)

dl_100_muK2 = cl_100_muK2*dl
dl_143_muK2 = cl_143_muK2*dl
dl_217_muK2 = cl_217_muK2*dl

print ("")
print "CIB power spectrum Dl in muK2 at 100, 143, and 217 GHz for Planck at l=%s is %s, %s, and %s" % (ell[ind3000], dl_100_muK2[ind3000], dl_143_muK2[ind3000], dl_217_muK2[ind3000])
# """

# ################# galaxy power spectra ########################

"""
instgal = ProfHODMore15(driver, driver_uni, gal_exp)
# cl_galtot = instgal.cl_galtot()
# cl_veltot = instgal.cl_veltot()


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


# ############ CIB x galaxies #####################################

z = reds  # redshifts  # reds
# z = np.linspace(min(reds), max(reds), 50)
driver_uni = cosmo_var_iv(mass, z, do_powerspec)

"""
r_l is the correlation coefficient between CIB and galaxies. It is being
used to calculate the cross shot nise between CIB and galaxies.
Set r_l = 0 if you do not want cross-shot noise or some other value if you
want the shot noise to be not 100% correlated.
"""
r_l = 1.0
# instcross = CIBxgal(instgal, instcib)
instcross = CIBxgal(driver, driver_uni, gal_exp, r_l)

# cl_crosstot = instcross.cibgalcross_cell_tot()

# """
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
    ax.legend(loc='upper right', prop={'size': 18}, frameon=False)
    ax.set_ylabel(r'$C^{\rm CIB x gal}_\ell\: [Jy\: {sr}^{-1}]$', fontsize=20)
    ax.set_xlabel('Multipole' r'$\;\ell$', fontsize=20)
    ax.set_title('%s CIB %s GHz x %s galaxy' % (exp['name'], exp['nu_string'][nu1], gal_exp), fontsize=20)
    plt.show()


# CCAT: nu0 = np.array([220., 280., 350., 410., 850.])  # in GHz
# Planck: nu0 = np.array([100., 143., 217., 353., 545., 857.])
nu1 = 3
plot_cross(nu1, instcross, exp['ell'])
# """

"""
# CCAT: nu0 = np.array([220., 280., 350., 410., 850.])  # in GHz
# Planck: nu0 = np.array([100., 143., 217., 353., 545., 857.])
# n_nu = 4
fsky = 0.4
for n_nu in range(len(custom['nu0'])):
    # snrg1, snrg2 = instcross.snr_g(n_nu, fsky)
    # snrg11, snrg12 = instcross.snr__g(n_nu, fsky, instcib, instgal)
    snrg12 = instcross.snr__g(n_nu, fsky, instcib, instgal)
    print ("SNR for CIB x gal at %s GHz for %s x %s is" % (str(int(custom['nu0'][n_nu])), custom['name'], gal_exp), snrg12[-1])

# """

t1 = time.time()
print t1 - t0
