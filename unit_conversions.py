# unit conversions
import numpy as np
from scipy.interpolate import interp1d
c = 3.e8  # m/s
hp = 6.63e-34  # SI
kB = 1.38e-23  # SI
Tcmb = 2.726  # K
Jansky = 1.e-26  # W/m^2/Hz
Jy = 1.e-26  # [W/m^2/Hz]


def blackbody(nu, T):
    '''blackbody function
    input: nu [Hz], T thermo temperature of the black body [K]
    output in SI: [W / Hz / m^2 / sr]
    '''
    x = hp*nu/(kB*T)
    result = 2.*hp*nu**3/c**2
    result /= np.exp(x) - 1.
    return result


def dBdT(nu, T):
    '''d(blackbody)/dT, such that
    dI = d(blackbody)/dT * dT
    input: nu [Hz], T thermo temperature of the black body [K]
    output in SI: [W / Hz / m^2 / sr / K]
    '''
    x = hp*nu/(kB*T)
    result = 2.*hp**2*nu**4
    result /= kB*T**2*c**2
    result *= np.exp(x) / (np.exp(x) - 1.)**2
    return result


def dlnBdlnT(nu, T):
    '''dlnBlackbody/dlnT
    input: nu [Hz], T thermo temperature of the black body [K]
    output is dimensionless
    '''
    x = hp*nu/(kB*T)
    return x * np.exp(x) / (np.exp(x) - 1.)


def dBdTrj(nu, T):
    '''d(blackbody)/dTrj, where Trj is the Rayleigh-Jeans brightness
    temperature,
    such that:
    dI = d(blackbody)/dTrj * dTrj
    input: nu [Hz], T thermo temperature [K]
    output in SI: [W / Hz / m^2 / sr / Krj]
    '''
    result = 2. * nu**2 * kB / c**2
    return result


def convertIntSITo(nu, kind="intSI"):
    '''kind: "intSI", "intJy/sr", "tempKcmb", "tempKrj"
    '''
    if kind == "intSI":
        result = 1.
    elif kind == "intJy/sr":
        result = 1. / Jy
    elif kind == "tempKcmb":
        result = 1. / dBdT(nu, Tcmb)
    elif kind == "tempKrj":
        result = 1. / dBdTrj(nu, Tcmb)
    return result


# Intensities [Jy/sr], as a function of freqs
# convert from muK to Jy/sr
nNu = 501
Nu = np.logspace(np.log10(0.1), np.log10(1.e4), nNu, 10.)*1.e9   # in Hz
factor = 1.e-6  # convert to K
factor /= convertIntSITo(Nu, kind="tempKcmb")  # convert to SI
factor *= convertIntSITo(Nu, kind="intJy/sr")  # convert to Jy/sr

factinterp = interp1d(Nu, factor, kind='linear', bounds_error=False,
                      fill_value=0.)
# print (factinterp(100*1.e9), factinterp(857*1.e9))
