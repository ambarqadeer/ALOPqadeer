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
    if not os.path.exists('Trimmed_Bias'):
        os.makedirs('Trimmed_Bias')
        print('created directory Trimmed_Bias')
    if not os.path.exists('Master_Files'):
        os.makedirs('Master_Files')
        print('created directory Master_Files')


def trim_bias(refresh='2'):
    biascollection = ImageFileCollection('HD115709/bias', ext=4)
    flag = 0
    if refresh == '1':
        tbiaspathlist = []
        for ccdb, biasn in biascollection.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            if flag == 0:
                print('all biases will be trimmed to :', ccdb.meta['trimsec'])
                flag = 1

            print('trimming', biasn)

            tbias = ccdp.trim_image(ccdb, fits_section=str(ccdb.meta['trimsec']))
            tbias.meta['imtype'] = ('trimmed bias', 'type of image')
            tbias.meta['taxis1'] = (2048, 'dimension1')
            tbias.meta['taxis2'] = (4096, 'dimension2')
            tbias.write('Trimmed_Bias/' + biasn[0:8] + '_trim.fits', overwrite=True)
            tbiaspathlist.append('Trimmed_Bias/' + biasn[0:8] + '_trim.fits')
        print('created', len(tbiaspathlist), 'trimmed biases')
    else:
        try:
            tbiascollection = ImageFileCollection('Trimmed_Bias')
            tbiaspathlist = tbiascollection.files_filtered(imtype='trimmed bias', include_path=True)
            print('found', len(tbiaspathlist), 'trimmed bias')
        except:
            print('can\'t locate trimmed biases, create or check directory')
            sys.exit()
    return biascollection, tbiaspathlist


def bias_combine(refresh='2', method='2'):
    tbiascollection = ImageFileCollection('Trimmed_Bias')
    combtime = 0
    if refresh == '1':
        print('found', len(tbiascollection.values('file')), 'trimmed biases')
        start = time.time()
        if method == '1':
            combined_bias = ccdp.combine(tbiascollection.files_filtered(
                imtype='trimmed bias', include_path=True),
                method='median')
            combbiaspath = 'Master_Files/mbias_median.fits'
            combined_bias.meta['combined'] = 'median'
            combtime = time.time() - start
            print('combination took', combtime, 'seconds')
            combined_bias.write(combbiaspath, overwrite=True)
        elif method == '2':
            combined_bias = ccdp.combine(tbiascollection.files_filtered(imtype='trimmed bias', include_path=True),
                                         sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                         sigma_clip_func=np.nanmedian,sigma_clip_dev_func=mad_std)
            combbiaspath = 'Master_Files/mbias.fits'
            combined_bias.meta['combined'] = 'sigma_clip average'
            combtime = time.time() - start
            print('combination took', combtime, 'seconds')
            combined_bias.write(combbiaspath, overwrite=True)
    else:
        try:
            if method == '1':
                combined_bias = CCDData.read('Master_Files/mbias_median.fits', unit='adu')
                combbiaspath = 'Master_Files/mbias_median.fits'
            elif method == '2':
                combined_bias = CCDData.read('Master_Files/mbias.fits', unit='adu')
                combbiaspath = 'Master_Files/mbias.fits'
        except:
            print('can\'t locate master bias, create or check directory')
            sys.exit()

    return tbiascollection, combined_bias, combbiaspath, combtime


# create directories to save files
create_directories()

# trim bias files, no refresh returns existing path list
print('do you want to trim the biases again? (1. Yes / 2. Read existing files)\n')
tbref = input()
biascollection, tbiaspathlist = trim_bias(tbref)

# combine bias files to form master bias, no refresh returns sigma_clip average mbias
print('do you want to create master bias again? (1. Yes / 2. Read existing files)\n')
mbref = input()
print('select combination method? (1. median / 2. sigma clipped average)\n')
method = input()
tbiascollection, mbias, combbiaspath, combtime = bias_combine(mbref, method)

