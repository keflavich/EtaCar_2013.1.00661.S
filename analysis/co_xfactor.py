import numpy as np
import radio_beam
from astropy import constants
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.stats import mad_std
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
Tex = 150
print("T_ex = {0}".format(Tex))
def J(T, nu=230.538*u.GHz):
    return (constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*u.Quantity(T, u.K)))-1)).to(u.K)
Nco = (4.79e13*(Tex+0.88)*np.exp(5.27/Tex) * (np.exp(5.27/Tex)-1)**-1 * (np.nanmax(cutout)) / (J(Tex)-J(2.73)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2)
print("Peak Nco={0:0.3g}".format(Nco))
print("Peak Nh2=1e4*Nco={0:0.3g}".format(Nco*1e4))
print("Peak Nh2=1e5*Nco={0:0.3g}".format(Nco*1e5))

Ncomap = (4.79e13*(Tex+0.88)*np.exp(5.27/Tex) * (np.exp(5.27/Tex)-1)**-1 * ((cutout)) / (J(Tex)-J(2.73)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2).astype('float128')
ok = np.isfinite(Ncomap)
Ncototal = np.nansum((Ncomap[ok] * pix_area_pc).decompose())
print('total co particles = {0}'.format(Ncototal))
print('total co mass = {0:0.2g}'.format((Ncototal*(12+16)*u.Da).to(u.M_sun)))
print('total h2 mass using X_CO=1e-5 = {0:0.2g}'.format(1e5*(Ncototal*(12+16)*u.Da).to(u.M_sun)))
print('total h2 mass using X_CO=1e-4 = {0:0.2g}'.format(1e4*(Ncototal*(12+16)*u.Da).to(u.M_sun)))


Tex = 1500
print("T_ex = {0}".format(Tex))
def J(T, nu=230.538*u.GHz):
    return (constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*u.Quantity(T, u.K)))-1)).to(u.K)
Nco = (4.79e13*(Tex+0.88)*np.exp(5.27/Tex) * (np.exp(5.27/Tex)-1)**-1 * (np.nanmax(cutout)) / (J(Tex)-J(2.73)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2)
print("Peak Nco={0:0.3g}".format(Nco))
print("Peak Nh2=1e4*Nco={0:0.3g}".format(Nco*1e4))
print("Peak Nh2=1e5*Nco={0:0.3g}".format(Nco*1e5))

Ncomap = (4.79e13*(Tex+0.88)*np.exp(5.27/Tex) * (np.exp(5.27/Tex)-1)**-1 * ((cutout)) / (J(Tex)-J(2.73)) / (u.km/u.s) * u.cm**-2).to(u.cm**-2).astype('float128')
ok = np.isfinite(Ncomap)
Ncototal = np.nansum((Ncomap[ok] * pix_area_pc).decompose())
print('total co particles = {0}'.format(Ncototal))
print('total co mass = {0:0.2g}'.format((Ncototal*(12+16)*u.Da).to(u.M_sun)))
print('total h2 mass using X_CO=1e-5 = {0:0.2g}'.format(1e5*(Ncototal*(12+16)*u.Da).to(u.M_sun)))
print('total h2 mass using X_CO=1e-4 = {0:0.2g}'.format(1e4*(Ncototal*(12+16)*u.Da).to(u.M_sun)))
