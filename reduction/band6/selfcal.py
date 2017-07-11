field = 'Eta_Carinae'
inp_vis = 'uid___A002_X9d26c8_X8c5.ms.split.cal'
prev_vis = 'original.split.cal'
split(vis=inp_vis, outputvis=prev_vis, datacolumn='data', field=field)

cell = '0.1arcsec'
imsize = [384,384]

def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)

def makefits_cube(myimagebase):
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image', fitsimage=myimagebase+'.image.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image


for iternum, threshold, solint in [(0, '1 Jy', '30s'),
                                   (1, '1 Jy', 'int'),
                                   (2, '0.5 Jy', 'int'),
                                   (3, '0.25 Jy', 'int'),
                                   (4, '0.1 Jy', 'int'),
                                   (5, '0.05 Jy', 'int'),]:
    vis = 'selfcal{0}.ms'.format(iternum)
    phasetable = 'phase_{0}.cal'.format(iternum)
    success = split(vis=prev_vis, outputvis=vis)
    if not success:
        success = split(vis=prev_vis, outputvis=vis, datacolumn='data')
    if not success:
        raise ValueError("Failed to split")

    myimagebase = 'EtaCar_band6_allspw_continuum_clean_taylor_selfcal{0}'.format(iternum)
    tclean(vis=vis, field=field,
           spw='', specmode='mfs', imsize=imsize, cell=cell,
           deconvolver='mtmfs',
           nterms=2,
           niter=50000,
           threshold=threshold,
           weighting='briggs',
           robust=0,
           savemodel='modelcolumn',
           imagename=myimagebase)

    makefits(myimagebase)

    rmtables([phasetable])
    gaincal(vis=vis, caltable=phasetable, solint='int', gaintype='G',
            calmode='p')

    applycal(vis=vis, field=field, gaintable=[phasetable],
             interp="linear", applymode='calflag', calwt=False)

    prev_vis = vis

inp_vis = vis
for spw in '0123':
    myimagebase = 'EtaCar_band6_spw{0}_selfcal_line_clean'.format(spw)
    os.system('rm -rf ' + myimagebase + "*/")
    tclean(vis=inp_vis,
           field='Eta_Carinae',
           spw=spw,
           specmode='cube',
           imsize=imsize,
           cell=cell,
           outframe='LSRK',
           niter=50000,
           threshold='50 mJy',
           weighting='briggs',
           robust=0,
           imagename=myimagebase)

    makefits_cube(myimagebase)

myimagebase = 'EtaCar_band6_CO2-1_selfcal_clean_uniform'
tclean(vis=inp_vis, field='Eta_Carinae',
       spw='0', specmode='cube', imsize=imsize, cell=cell,
       robust=-2, weighting='briggs',
       outframe='LSRK',
       restfreq='230.538GHz',
       reffreq='230.538GHz',
       niter=50000,
       threshold='50 mJy',
       start='-250km/s',
       nchan=400,
       imagename=myimagebase)

makefits_cube(myimagebase)
