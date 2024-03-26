import numpy as np
import os
import time

##
def get_convo2target(convo_file, ref_image):
    """
    Convolve input image to the CASA simulated image resolution.
    """
    # Get beam info from reference image
    hdr = imhead(ref_image, mode='summary')
    beam_major = hdr['restoringbeam']['major']
    beam_minor = hdr['restoringbeam']['minor']
    beam_PA = hdr['restoringbeam']['positionangle']
    ref_unit = hdr['unit']

    outfile_tmp = convo_file+'_convtmp'
    outfile_convo = convo_file+'_conv'

    # Convolution into the same beam as the reference
    os.system('rm -rf %s' % outfile_tmp)
    imsmooth(imagename=convo_file, outfile=outfile_tmp, kernel='gauss', 
             major=beam_major, minor=beam_minor, pa=beam_PA, targetres=True)
    imregrid(imagename=outfile_tmp, template=ref_image, output=outfile_convo, overwrite=True)
    imhead(outfile_convo, mode='put', hdkey='Bunit', hdvalue=ref_unit)
    os.system('rm -rf %s' % outfile_tmp)


##
def CASA2fits(CASAfile):
    """
    Export CASA image to FITS format.
    """
    exportfits(imagename=CASAfile, fitsimage=CASAfile + ".fits", overwrite=True)


##
def convert_JypB_JypP(imagename):
    """
    Convert image brightness unit from Jy/beam to Jy/pixel.
    """
    myimhead = imhead(imagename)

    # Get beam and pixel information
    bmaj = myimhead['restoringbeam']['major']['value']  # in Arcsec
    bmin = myimhead['restoringbeam']['minor']['value']  # in Arcsec
    pix = abs(myimhead['incr'][0]) * 206265.  # in Arcsec

    # Compute beam and pixel areas
    beam_area = (np.pi * bmaj * bmin) / (4 * np.log(2))
    pix_area = pix ** 2

    # Convert Jy/beam to Jy/pixel
    toJyPerPix = pix_area / beam_area
    SDEfficiency = 1.0
    fluxExpression = "(IM0 * {0:f} / {1:f})".format(toJyPerPix, SDEfficiency)

    # Generate scaled image
    scaled_name = imagename + '.Jyperpix'
    os.system('rm -rf ' + scaled_name)

    immath(imagename=imagename,
                outfile=scaled_name,
                mode='evalexpr',
                expr=fluxExpression)

    hdval = 'Jy/pixel'
    scaled_name = imhead(imagename=scaled_name,
                        mode='put',
                        hdkey='BUNIT',
                        hdvalue=hdval)

    return scaled_name
