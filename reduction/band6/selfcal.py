field = 'Eta_Carinae'
prev_vis = inp_vis = 'uid___A002_X9d26c8_X8c5.ms.split.cal'
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


for iternum, threshold, solint in [(0, '1 Jy', '30s'),
                                   (1, '1 Jy', 'int'),
                                   (2, '0.5 Jy', 'int'),
                                   (3, '0.25 Jy', 'int'),
                                   (4, '0.1 Jy', 'int'),
                                   (5, '0.05 Jy', 'int'),]:
    vis = 'selfcal{0}.ms'.format(iternum)
    phasetable = 'phase_{0}.cal'.format(iternum)
    split(prev_vis, vis)

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
