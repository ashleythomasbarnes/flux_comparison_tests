
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.modeling import models, fitting
import numpy as np
from astropy.modeling import models, fitting
from astropy.io import fits
import astropy.wcs as wcs
import astropy.units as u
from astropy.convolution import convolve_fft
import numpy as np
from radio_beam import Beam
from astropy.modeling import models, fitting

import warnings 
warnings.filterwarnings('ignore')

def remove_padding(data):
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        return data
    
    # Crop the data array
    cropped_data = data[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    return cropped_data

def remove_padding_2arrays(data1, data2):
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data1, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data1, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        return data1
    
    # Crop the data array
    cropped_data = data2[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    return cropped_data


def fit_2d_gaussian(image, x_stddev=None, y_stddev=None):

    image = np.squeeze(image)
    image[np.isnan(image)] = 0

    # Get the center of the image
    shape_x, shape_y = image.shape
    center_x, center_y = np.array(image.shape) // 2
    
    # Create x, y indices grid for the sub-image
    x, y = np.array(np.mgrid[:shape_x, :shape_y], dtype=np.int32)
    
    # Initialize the Gaussian2D model
    g_init = models.Gaussian2D(amplitude=np.nanmax(image), 
                               x_mean=center_x, y_mean=center_y, 
                               x_stddev=x_stddev, y_stddev=y_stddev)
    
    # Fit the model to the sub-image
    fit_g = fitting.LevMarLSQFitter(calc_uncertainties=True)
    g = fit_g(g_init, x, y, image)

    # Calculate the uncertainty on the fit using the covariance matrix
    cov_matrix = fit_g.fit_info['param_cov']
    if cov_matrix is not None:
        g_err = np.sqrt(np.diag(cov_matrix))
    else:
        g_err = np.zeros_like(g.parameters)

    g_errl = models.Gaussian2D(amplitude=g.amplitude.value - np.abs(g_err[0]), 
                        x_mean=g.x_mean.value, 
                        y_mean=g.y_mean.value, 
                        x_stddev=g.x_stddev.value -  np.abs(g_err[3]), 
                        y_stddev=g.y_stddev.value -  np.abs(g_err[4]),
                        theta=g.theta.value)

    g_errh = models.Gaussian2D(amplitude=g.amplitude.value + np.abs(g_err[0]), 
                                x_mean=g.x_mean.value, 
                                y_mean=g.y_mean.value, 
                                x_stddev=g.x_stddev.value + np.abs(g_err[3]), 
                                y_stddev=g.y_stddev.value + np.abs(g_err[4]),
                                theta=g.theta.value)
    
    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g(x, y), dtype=np.float32)
    sum_tot = np.nansum(fitted_data)

    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g_errl(x, y), dtype=np.float32)
    sum_errl = np.nansum(fitted_data)

    # Calculate the sum of the fitted Gaussian
    fitted_data = np.array(g_errh(x, y), dtype=np.float32)
    sum_errh = np.nansum(fitted_data)

    # Compute the reduced chi-squared
    residuals = image - g(x, y)
    chi_squared = np.sum((residuals / image.std())**2)
    degrees_of_freedom = image.size - len(g.parameters)
    reduced_chi_squared = chi_squared / degrees_of_freedom

    return (g, g_err, (sum_tot*u.Jy, sum_errl*u.Jy, sum_errh*u.Jy), reduced_chi_squared)

def get_smooth(hdu, hdu_obs):
    
    # Create a WCS object from the input HDU header
    wcs_ = wcs.WCS(hdu.header)

    # Calculate the pixel scale in degrees
    pixscale = wcs.utils.proj_plane_pixel_area(wcs_.celestial) ** 0.5 * u.deg
    print(f"[INFO] Pixel scale: {pixscale.to('arcsec'):.2f} arcsec")

    # Define the initial and desired beams
    initial_beam = Beam(0*u.arcsec)
    desired_beam = Beam.from_fits_header(hdu_obs.header) 
    
    # Create the convolution kernel
    convolution_kernel = desired_beam.deconvolve(initial_beam).as_kernel(pixscale)

    # Convolve the image with the kernel to smooth it
    smoothed_data = convolve_fft(hdu.data, convolution_kernel, preserve_nan=True, allow_huge=True)
    output_hdu = fits.PrimaryHDU(np.array(smoothed_data, dtype=np.float32), hdu.header)

    return(output_hdu)