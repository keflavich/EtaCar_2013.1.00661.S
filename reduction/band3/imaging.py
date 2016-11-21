
# continuum
for spw in '0123':
    myimagebase = 'EtaCar_band3_spw{0}_continuum_clean'.format(spw)
    tclean(vis='uid___A002_X9baf64_X5d_raw.ms.split.cal', field='Eta_Carinae',
           spw=spw, specmode='mfs', imsize=[384,384], cell='0.25arcsec',
           niter=5000,
           imagename=myimagebase)

    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band3_allspw_continuum_clean_uniform'.format(spw)
tclean(vis='uid___A002_X9baf64_X5d_raw.ms.split.cal', field='Eta_Carinae',
       spw='', specmode='mfs', imsize=[384,384], cell='0.15arcsec',
       robust=-2, weighting='briggs',
       niter=5000,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band3_allspw_continuum_clean'.format(spw)
tclean(vis='uid___A002_X9baf64_X5d_raw.ms.split.cal', field='Eta_Carinae',
       spw='', specmode='mfs', imsize=[384,384], cell='0.25arcsec',
       niter=5000,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band3_allspw_continuum_clean_taylor'.format(spw)
tclean(vis='uid___A002_X9baf64_X5d_raw.ms.split.cal', field='Eta_Carinae',
       spw='', specmode='mfs', imsize=[384,384], cell='0.25arcsec',
       deconvolver='mtmfs',
       nterms=2,
       niter=5000,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.tt0', pbimage=myimagebase+'.pb', outfile=myimagebase+'.tt0.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.tt0.pbcor', fitsimage=myimagebase+'.tt0.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.tt1', fitsimage=myimagebase+'.tt1.fits', dropdeg=True, overwrite=True)


for spw in '0123':
    myimagebase = 'EtaCar_band3_spw{0}_line_clean'.format(spw)
    tclean(vis='uid___A002_X9baf64_X5d_raw.ms.split.cal', field='Eta_Carinae',
           spw=spw, specmode='cube', imsize=[384,384], cell='0.25arcsec',
           outframe='LSRK',
           niter=5000,
           imagename=myimagebase)

    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band3_CO_clean_uniform'.format(spw)
tclean(vis='uid___A002_X9baf64_X5d_raw.ms.split.cal', field='Eta_Carinae',
       spw='0', specmode='cube', imsize=[384,384], cell='0.15arcsec',
       robust=-2, weighting='briggs',
       outframe='LSRK',
       restfreq='115.271208GHz',
       reffreq='115.271208GHz',
       niter=5000,
       start='-250km/s',
       nchan=400,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
