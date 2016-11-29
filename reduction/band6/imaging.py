inp_vis = 'uid___A002_X9d26c8_X8c5.ms.split.cal'
cell = '0.1arcsec'
imsize = [384,384]

# continuum
for spw in '0123':
    myimagebase = 'EtaCar_band6_spw{0}_continuum_clean'.format(spw)
    tclean(vis=inp_vis, field='Eta_Carinae',
           spw=spw, specmode='mfs', imsize=imsize, cell=cell,
           niter=5000,
           imagename=myimagebase)

    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band6_allspw_continuum_clean'.format(spw)
tclean(vis=inp_vis, field='Eta_Carinae',
       spw='', specmode='mfs', imsize=imsize, cell=cell,
       niter=5000,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band6_allspw_continuum_clean_uniform'.format(spw)
tclean(vis=inp_vis, field='Eta_Carinae',
       spw='', specmode='mfs', imsize=imsize, cell=cell,
       weighting='briggs', robust=-2,
       niter=5000,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)



myimagebase = 'EtaCar_band6_allspw_continuum_clean_taylor'.format(spw)
tclean(vis=inp_vis, field='Eta_Carinae',
       spw='', specmode='mfs', imsize=imsize, cell=cell,
       deconvolver='mtmfs',
       nterms=2,
       niter=5000,
       imagename=myimagebase)

#impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True)
#exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.image.tt0', fitsimage=myimagebase+'.image.tt0.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True)


for spw in '0123':
    myimagebase = 'EtaCar_band6_spw{0}_line_clean'.format(spw)
    tclean(vis=inp_vis, field='Eta_Carinae',
           spw=spw, specmode='cube', imsize=imsize, cell=cell,
           outframe='LSRK',
           niter=5000,
           imagename=myimagebase)

    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)

myimagebase = 'EtaCar_band6_CO2-1_clean_uniform'
tclean(vis=inp_vis, field='Eta_Carinae',
       spw='0', specmode='cube', imsize=imsize, cell=cell,
       robust=-2, weighting='briggs',
       outframe='LSRK',
       restfreq='230.538GHz',
       reffreq='230.538GHz',
       niter=5000,
       start='-250km/s',
       nchan=400,
       imagename=myimagebase)

impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True)
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True)
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True)
