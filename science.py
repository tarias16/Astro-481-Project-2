#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: science.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
from astropy.stats import sigma_clip
from scipy.stats import mode
from astroscrappy import detect_cosmics
import os


def reduce_science_frame(
    science_filename,
    median_bias_filename,
    median_flat_filename,
    reduced_science_filename,
    save_loc,
    illumination_filename = None,
    median_dark_filename = None
):
    """This function must:

    - Accept a science frame filename as science_filename.
    - Accept a median bias frame filename as median_bias_filename (the one you created
      using create_median_bias).
    - Accept a median flat frame filename as median_flat_filename (the one you created
      using create_median_flat).
    - Accept a median dark frame filename as median_dark_filename (the one you created
      using create_median_dark).
    - Read all files.
    - Subtract the bias frame from the science frame.
    - Subtract the dark frame from the science frame. Remember to multiply the
      dark frame by the exposure time of the science frame. The exposure time can
      be found in the header of the FITS file.
    - Correct the science frame using the flat frame.
    - Optionally, remove cosmic rays.
    - Save the resulting reduced science frame to a FITS file with the filename
      reduced_science_filename.
    - Return the reduced science frame as a 2D numpy array.

    """
    cwd = os.getcwd()
    path = os.path.abspath(science_filename)
    os.chdir(path)

    images = os.listdir(path)

    v_images =[]
    b_images = []
    for i in range(0, len(images)):
      header = fits.getheader(images[i])
      if header['FILTER'] == 'B':
          b_images.append(images[i])
      elif header['FILTER'] == 'V':
         v_images.append(images[i])

    for i in range(0,len(b_images)):
      os.chdir(path)
      header = fits.getheader(b_images[i])
      # Gets the science fits data from science_filename
      science_frame = (fits.getdata(b_images[i])).astype('float32')

      # Gets the exposure time of the science image
      sci_exp = (fits.getheader(b_images[i]))['EXPTIME']

      os.chdir(cwd)
      # Reduces the science image by subtracting the provided median bias and dark file
      reduced_science = science_frame - fits.getdata(median_bias_filename)

      # Subtracts dark if dark is provided
      if median_dark_filename:
        reduced_science = reduced_science - (fits.getdata(median_dark_filename)*sci_exp)

      # Divides out the normalizes median_flat image
      reduced_science /= fits.getdata(f'{median_flat_filename}_B.fits')

      # Divides out illumination
      if illumination_filename:
        reduced_science /= fits.getdata(f'{illumination_filename}_B.fits')
      # Performs cosmic ray detection
      mask, output = detect_cosmics(reduced_science)

      # Sets up fits write to file

      os.chdir(save_loc)
      primary = fits.PrimaryHDU(data=output.data, header=header)
      # primary.header["AIRMASS"] = airmass
      hdul = fits.HDUList([primary])
      hdul.writeto(f'{reduced_science_filename}_B_{i+1}.fits', overwrite=True)
    
    for i in range(0,len(v_images)):
      os.chdir(path)
      header = fits.getheader(v_images[i])
      # Gets the science fits data from science_filename
      science_frame = (fits.getdata(v_images[i])).astype('float32')

      # Gets the exposure time of the science image
      sci_exp = (fits.getheader(v_images[i]))['EXPTIME']

      os.chdir(cwd)
      # Reduces the science image by subtracting the provided median bias and dark file
      reduced_science = science_frame - fits.getdata(median_bias_filename)

      # Subtracts dark if dark is provided
      if median_dark_filename:
        reduced_science = reduced_science - (fits.getdata(median_dark_filename)*sci_exp)

      # Divides out the normalizes median_flat image
      reduced_science /= fits.getdata(f'{median_flat_filename}_V.fits')

      # Divides out illumination
      if illumination_filename:
        reduced_science /= fits.getdata(f'{illumination_filename}_V.fits')
      # Performs cosmic ray detection
      mask, output = detect_cosmics(reduced_science)

      # Sets up fits write to file
      os.chdir(save_loc)
      primary = fits.PrimaryHDU(data=output.data, header=header)
      # primary.header["AIRMASS"] = airmass
      hdul = fits.HDUList([primary])
      hdul.writeto(f'{reduced_science_filename}_V_{i+1}.fits', overwrite=True)

    os.chdir(cwd)
    return output
