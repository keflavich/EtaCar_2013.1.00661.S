import numpy as np
import radio_beam
from astropy import constants
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.stats import mad_std
from astroquery.splatalogue import Splatalogue
import regions

distance = 2.3*u.kpc

fn = "../FITS/EtaCar_CO2-1_mom0.fits"

fh = fits.open(fn)
mywcs = wcs.WCS(fh[0].header)
pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
data = fh[0].data * u.Unit(fh[0].header['BUNIT'])
#beam = radio_beam.Beam.from_fits_header(fh[0].header)
#ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam


# https://ned.ipac.caltech.edu/level5/Sept13/Bolatto/Bolatto4.html
# N(H2) is the output
co_xfactor = 2e20 * u.cm**-2/(u.K * u.km/u.s)
alpha_co = 4.3 * u.M_sun /(u.K * u.km/u.s * u.pc**2)

cutout_region = regions.io.read_ds9("../regions/etacar_cutout.reg")[0]

pixregion = cutout_region.to_pixel(mywcs)
pixmask = pixregion.to_mask()

jtok = u.brightness_temperature(radio_beam.Beam.from_fits_header(fh[0].header), 230*u.GHz, )
cutout = (pixmask.cutout(data)/(u.km/u.s)).to(u.K, equivalencies=jtok)*u.km/u.s

cutout_nh = co_xfactor * cutout

cutout_mass = (cutout_nh * (u.Da*pixscale**2*distance**2)).to(u.M_sun,
                                                              u.dimensionless_angles())
cutout_mass = (cutout * pixscale**2 * distance**2 * alpha_co).to(u.M_sun,
                                                                 u.dimensionless_angles())

mass = np.nansum(cutout_mass)
print("Mass inferred: {0:0.3g}".format(mass))
print("Peak integrated intensity: {0:0.3g}".format(np.nanmax(cutout)))
print("Mad-std integrated intensity: {0:0.3g}".format(mad_std(cutout, ignore_nan=True)))
print("Average integrated intensity: {0:0.3g}".format(np.nanmean(cutout)))
print("Std Dev integrated intensity: {0:0.3g}".format(np.nanstd(cutout)))
print("Area: {0:0.3g}".format((pixscale**2*np.isfinite(cutout).sum()).to(u.arcsec**2)))
pix_area_pc = (pixscale**2*distance**2).to(u.pc**2, u.dimensionless_angles())
area_pc = (np.isfinite(cutout).sum()*pix_area_pc)
print("Area: {0:0.3g}".format(area_pc))



# use Mangum eqn 90 (C18O) to approximate N(CO)
def J(T, nu=230.538*u.GHz):
    return (constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*u.Quantity(T, u.K)))-1)).to(u.K)

def Nco_c18o_10(Tex, flux):
    Tex = u.Quantity(Tex, u.K)
    flux = u.Quantity(flux, u.K*u.km/u.s)
    return (2.48e14*(Tex.value+0.88)*np.exp(5.27/Tex.value) * (np.exp(5.27/Tex.value)-1)**-1 * flux / (J(Tex)-J(2.7315*u.K)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2)

def Qrot(Tex, B0=54891.420*u.MHz):
    Tex = u.Quantity(Tex, u.K)
    return (constants.k_B*Tex)/(constants.h*B0) + 1/3.

def Nu_thin(flux,  tex, line_strength=1.1e-19*u.esu*u.cm, freq=230*u.GHz, Tbg=2.7315*u.K, fillingfactor=1.0, tau=None):
    """
    Derived from Eqn 30 + Eqn 80 of Mangum 2015
    """
    tex = u.Quantity(tex, u.K)
    assert flux.unit.is_equivalent(u.K*u.km/u.s)
    assert line_strength.unit.is_equivalent(u.esu*u.cm)
    hnu = constants.h * freq
    h = constants.h
    kbt = constants.k_B * tex
    term1 = (3*h/(8*np.pi**3 * line_strength**2))
    term3 = 1. / (np.exp(hnu/kbt) - 1)
    term4 = 1./(tex-Tbg)
    term5 = flux.to(u.K*u.km/u.s) / fillingfactor
    term6 = 1 if tau is None else tau/(1-np.exp(-tau))
    return (term1*term3*term4*term5*term6).to(u.cm**-2)

def Ntot(flux, tex, Ju, degeneracy_u, Eu, B0=54891.420*u.MHz, **kwargs):
    tex = u.Quantity(tex, u.K)
    #S = Ju/(2*Ju+1)
    # degeneracy_u = 2J + 1
    # Ju = S * degeneracy_u
    return Nu_thin(flux, tex, **kwargs) * Qrot(tex, B0) / Ju * np.exp(Eu/tex)

# molecular constants
# C18O 1-0
# the PASJ version is wrong, the arxiv version updated the constants
muC18O = 0.11079 * u.Debye
muC18O = 1.1079e-19 * u.esu * u.cm # Jeff had 1.098 (that was the wrong version)
B0C18O = 54891.420 * u.MHz # apparently Jeff had used the 12C16O version (again, old version)
EuC18O = 5.27 * u.K
nuC18O = 109.782182 * u.GHz

# calculate prefactor
prefactor_mangum = 4.79e13 # old, incorrect version
prefactor_mangum = 2.48e14
Ju = 1
prefactor = ((3 * constants.h) / (8 * np.pi**3 * muC18O**2 * Ju) * (constants.k_B / (constants.h * B0C18O))).to((u.km/u.s)**-1 * u.cm**-2 * u.K**-1)
prefactor_ = ((3 * 6.62607004e-27) / (8 * 3.14159**3 * 1.1079e-19**2 * 1) * (1.38064852e-16 / (6.62607004e-27 * 54891.420e6))) * 1e5
print("Prefactor for C18O eqn 90:  Jeff={0:0.3g}  me={1:0.3g}  ratio={2:0.4g}"
      .format(prefactor_mangum, prefactor, (prefactor/prefactor_mangum).value,))

# sanity check
print("Ntot full calc vs jeff: {0:0.4g} : {1:0.4g}   r={2}"
      .format(Ntot(1*u.K*u.km/u.s, 100*u.K, degeneracy_u=3, freq=nuC18O, Eu=EuC18O, line_strength=muC18O, B0=B0C18O, Ju=1), Nco_c18o_10(100*u.K, 1*u.K*u.km/u.s),
              Ntot(1*u.K*u.km/u.s, 100*u.K, degeneracy_u=3, freq=nuC18O, Eu=EuC18O, line_strength=muC18O, B0=B0C18O, Ju=1)/ Nco_c18o_10(100*u.K, 1*u.K*u.km/u.s))
             )

# CO 1-0
mu = 0.11011 * u.Debye
mu = 1.1011e-19 * u.esu*u.cm
B0 = 57635.968 * u.MHz
Eu = Splatalogue.query_lines(100*u.GHz, 235*u.GHz, chemical_name=' CO ')['E_U (K)'][-1] * u.K
Eu = 5.27 * u.K
nu = 230.538 * u.GHz / 2.
Ju = 1
degeneracy_u = 2*Ju+1


# CO 2-1
mu = 0.11011 * u.Debye
mu = 1.1011e-19 * u.esu*u.cm
B0 = 57635.968 * u.MHz
Eu = Splatalogue.query_lines(100*u.GHz, 235*u.GHz, chemical_name=' CO ')['E_U (K)'][-1] * u.K
Eu = 16.59608 * u.K
nu = 230.538 * u.GHz
Ju = 2
degeneracy_u = 2*Ju+1


prefactor = ((3 * constants.h) / (8 * np.pi**3 * mu**2 * Ju) * (constants.k_B / (constants.h * B0))).to((u.km/u.s)**-1 * u.cm**-2 * u.K**-1).value
def Nco_12co_21(Tex, flux, eu_k=Eu.value, nu=nu):
    Tex = u.Quantity(Tex, u.K)
    flux = u.Quantity(flux, u.K*u.km/u.s)
    return (prefactor*(Tex.value+0.92)*np.exp(eu_k/Tex.value) * (np.exp(constants.h*nu/(constants.k_B*Tex))-1)**-1 * flux / (J(Tex)-J(2.7315*u.K)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2)

print("Prefactor for 12CO 2-1 = {0:0.3g}" .format(prefactor))
print("Ntot full calc vs prefactor: {0:0.4g} : {1:0.4g}   r={2}"
      .format(Ntot(1*u.K*u.km/u.s, 100*u.K, degeneracy_u=2*Ju+1, freq=nu, Eu=Eu, line_strength=mu, B0=B0, Ju=Ju), Nco_12co_21(100*u.K, 1*u.K*u.km/u.s),
              Ntot(1*u.K*u.km/u.s, 100*u.K, degeneracy_u=2*Ju+1, freq=nu, Eu=Eu, line_strength=mu, B0=B0, Ju=Ju)/ Nco_12co_21(100*u.K, 1*u.K*u.km/u.s))
             )
print()


for Tex in (150, 1500, 2000):
    print("Tex = {0}".format(Tex))
    pk = np.nanmax(cutout)
    print("Peak Nco={0:0.3g}".format(Nco_12co_21(Tex, pk)))
    print("Peak Nh2=1e4*Nco={0:0.3g}".format(Nco_12co_21(Tex, pk)*1e4))
    print("Peak Nh2=1e4*Nco={0:0.3g}".format(Ntot(pk, Tex, degeneracy_u=degeneracy_u, freq=nu, Eu=Eu, line_strength=mu, Ju=Ju)*1e4))
    print("Peak Nh2=1e5*Nco={0:0.3g}".format(Nco_12co_21(Tex, pk)*1e5))

    Ncomap = Nco_12co_21(Tex, cutout).astype('float128')
    ok = np.isfinite(Ncomap)
    Ncototal = np.nansum((Ncomap[ok] * pix_area_pc).decompose())
    print('total co particles = {0}'.format(Ncototal))
    print('total co mass = {0:0.3g}'.format((Ncototal*(12+16)*u.Da).to(u.M_sun)))
    print('total h2 mass using X_CO=1e-5 = {0:0.3g}'.format(1e5*(Ncototal*(12+16)/2.*u.Da).to(u.M_sun)))
    print('total h2 mass using X_CO=1e-4 = {0:0.3g}'.format(1e4*(Ncototal*(12+16)/2.*u.Da).to(u.M_sun)))

    print()
