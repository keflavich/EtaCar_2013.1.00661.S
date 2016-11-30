from astropy import units as u
from astropy import wcs
from astropy.io import fits
import regions

distance = 2.3*u.kpc

fn = "../FITS/EtaCar_CO2-1_mom0.fits"

fh = fits.open(fn)
mywcs = wcs.WCS(fh[0].header)
pixscale = (mywcs.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
#beam = radio_beam.Beam.from_fits_header(fh[0].header)
#ppbeam = (beam.sr/(pixscale**2*u.deg**2)).decompose().value / u.beam


# https://ned.ipac.caltech.edu/level5/Sept13/Bolatto/Bolatto4.html
co_xfactor = 1e20 * u.cm**-2/(u.K * u.km/u.s)

cutout_region = regions.io.read_ds9("../regions/etacar_cutout.reg")[0]

pixregion = cutout_region.to_pixel(mywcs)
pixmask = pixregion.to_mask()

cutout = pixmask.cutout(fh[0].data) * u.K * u.km/u.s

cutout_nh2 = co_xfactor * cutout

cutout_mass = (cutout_nh2 * (2.8*u.Da*pixscale**2*distance**2)).to(u.M_sun,
                                                                   u.dimensionless_angles())
