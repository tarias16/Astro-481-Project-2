#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: ptc.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
import numpy 
from astropy.io import fits
import os

def calculate_gain(files):
    """This function must:

    - Accept a list of files that you need to calculate the gain
      (two files should be enough, but what kind?).
    - Read the files and calculate the gain in e-/ADU.
    - Return the gain in e-/ADU.

    """
    cwd = os.getcwd()
    path = os.path.abspath(files)
    os.chdir(path)
    flats = os.listdir(path)

    # Defines list of flats for calculation
    flat_gain_list = []

    # Reads the flat files provided
    for i in range(0, 2):
        data = fits.getdata(flats[i])
        flat_gain_list.append((data).astype('float32'))

    os.chdir(cwd)

    # Takes the difference in signl between the flat images
    flat_diff = flat_gain_list[0] - flat_gain_list[1]

    # Calculates the variation in the difference between the flats
    flat_var = numpy.var(flat_diff)

    # Defines the average signal in each image
    S = 0.5 * numpy.mean(flat_gain_list[0] + flat_gain_list[1])

    # Uses signal and varience to calculate the gain
    gain = (2*S / flat_var).item()

    return gain


def calculate_readout_noise(files, gain):
    """This function must:

    - Accept a list of files that you need to calculate the readout noise
      (two files should be enough, but what kind?).
    - Accept the gain in e-/ADU as gain. This should be the one you calculated
      in calculate_gain.
    - Read the files and calculate the readout noise in e-.
    - Return the readout noise in e-.

    """

    cwd = os.getcwd()
    path = os.path.abspath(files)
    os.chdir(path)
    biases = os.listdir(path)


    # Defines list for bias fiels
    bias_rn_list = []

    # Reads provided bias files
    for i in range(0, 2):
        data = fits.getdata(biases[i])
        bias_rn_list.append((data).astype('float32'))

    os.chdir(cwd)
    
    # Finds difference in signal between images
    bias_diff = bias_rn_list[0] - bias_rn_list[1]

    # Calculates varience of bias_diff
    bias_diff_var = numpy.var(bias_diff)

    # Calculates the readout noise in ADU
    rn_adu = numpy.sqrt(bias_diff_var / 2)

    # Uses the provided gain from calculate_gain() to convert readout noise from ADU to e-
    readout_noise = (rn_adu * gain).item()

    return readout_noise
