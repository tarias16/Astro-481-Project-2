#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: reduction.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
from astropy.io import fits
import os
import shutil

def rename_images(folder):
    path = os.path.abspath(folder)
    files = os.listdir(folder)
   
    for i in range(0,len(files)):
        files[i] = os.path.join(path, files[i])
        file = fits.getheader(files[i])
        time = file['DATE-OBS']
        time = time.replace(':', '')
        time = time.replace('.', '')
        comm = file['COMMENT']
        comm= str(comm)
        comm = comm.replace('/', '')
        
        # os.rename(files[i], os.path.join(path, f'{file['IMAGETYP']}-{comm}-{file['FILTER']}-{time}.fits'))
        os.rename(files[i], os.path.join(path, f'{comm}-{file['FILTER']}-{time}.fits'))

    return

def sort_images(folder, location):
    path = os.path.abspath(folder)
    files = []
    for entry in os.listdir(folder):
         full_path = os.path.join(folder, entry)
         if os.path.isfile(full_path):
            files.append(entry)
    
    loc = os.path.abspath(location)

    for i in range(0,len(files)):
        files[i] = os.path.join(path,files[i])
        header = fits.getheader(files[i])

        if header['COMMENT'] == 'AC AND':
            shutil.move(files[i], loc)
    return
