import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
import astropy.units as u
import ccdproc as ccdp
import os
import pathlib
from ccdproc import ImageFileCollection
from astropy.visualization import hist
import itertools
from astropy.stats import sigma_clip, mad_std
import time
import sys


def create_directories():
    if not os.path.exists('Trimmed_Flat'):
        os.makedirs('Trimmed_Flat')
        print('created directory Trimmed_Flat')
    if not os.path.exists('Master_Files'):
        os.makedirs('Master_Files')
        print('created directory Master_Files')
    if not os.path.exists('Trimmed_Flat/subflatsmed'):
        os.makedirs('Trimmed_Flat/subflatsmed')
        print('created directory Trimmed_Flat/subflatsmed')
    if not os.path.exists('Trimmed_Flat/subflatssig'):
        os.makedirs('Trimmed_Flat/subflatssig')
        print('created directory Trimmed_Flat/subflatssig')


def trim_flat(refresh = '2'):
    flatcollection = ImageFileCollection('HD115709/flat_SII', ext=4)
    flag = 0
    tflatpathlist = []
    if refresh == '1':
        for ccdf, flatn in flatcollection.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            if flag == 0:
                print('all flats will be trimmed to :', ccdf.meta['trimsec'])
                flag = 1

            print('trimming', flatn)

            tflat = ccdp.trim_image(ccdf, fits_section=str(ccdf.meta['trimsec']))
            tflat.meta['imtype'] = ('trimmed flat', 'type of image')
            tflat.meta['taxis1'] = (2048, 'dimension1')
            tflat.meta['taxis2'] = (4096, 'dimension2')
            tflat.write('Trimmed_Flat/' + flatn[0:8] + '_trim.fits', overwrite=True)
            tflatpathlist.append('Trimmed_Flat/' + flatn[0:8] + '_trim.fits')
        print('created', len(tflatpathlist), 'trimmed flats')
    elif refresh == '2':
        try:
            tflatcollection = ImageFileCollection('Trimmed_Flat')
            tflatpathlist = tflatcollection.files_filtered(imtype='trimmed flat', include_path=True)
            print('found', len(tflatpathlist), 'trimmed flats')
        except:
            print('can\'t locate trimmed flats, create or check directory')
            sys.exit(0)
    return flatcollection, tflatpathlist


def sub_bias(refresh='2',bias='2'):
    tflatcollection = ImageFileCollection('Trimmed_Flat')
    if bias == '1':
        biaspath = 'Master_Files/mbias_median.fits'
        dest = 'Trimmed_Flat/subflatsmed/'
    elif bias == '2':
        biaspath = 'Master_Files/mbias.fits'
        dest = 'Trimmed_Flat/subflatssig/'
    if refresh == '1':
        subflatpathlist = []
        mbias = CCDData.read(biaspath, unit='adu')
        for ccdf, flatn in tflatcollection.ccds(imtype='trimmed flat', return_fname=True):
            subflat=ccdp.subtract_bias(ccdf,mbias,add_keyword='subbias')
            subflat.meta['imtype']=('subbias flat','bias subtracted flat')
            subflat.write(dest+flatn[0:8]+'_subbias.fits')
            subflatpathlist.append(dest+flatn[0:8]+'_subbias.fits')
    else:
        try:
            subflatcollection = ImageFileCollection(dest)
            subflatpathlist = subflatcollection.files_filtered(imtype='subbias flat', include_path=True)
            print('found', len(subflatpathlist), 'subflats')
        except:
            print('can\'t locate subflats, create or check directory')
            sys.exit()
    return tflatcollection, subflatpathlist


# create directories to save files
create_directories()

# trim flat files, no refresh returns existing path list
print('do you want to trim the flats? (1. Yes / 2. Read existing files)\n')
tfref = input()
flatcollection, tflatpathlist = trim_flat(tfref)


# subtract bias from flats
refresh=input('do you want to subtract bias from flats? (1. Yes / 2. Read existing files)\n')
bias=input('select which bias to use (1. Median / 2. Sigma clipped average): \n')
tflatcollection, subflatpathlist = sub_bias(refresh,bias)
