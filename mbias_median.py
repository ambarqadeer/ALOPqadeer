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


def create_directories():
    if not os.path.exists('Trimmed_Bias'):
        os.makedirs('Trimmed_Bias')
        print('created directory Trimmed_Bias')
    if not os.path.exists('Master_Files'):
        os.makedirs('Master_Files')
        print('created directory Master_Files')


# pathlib.Path().absolute()
def trim_bias(refresh=False):
    biascollection = ImageFileCollection('HD115709/bias', ext=4)
    flag = 0
    if refresh == True:
        for ccdb, biasn in biascollection.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):

            print('trimming', biasn)

            ccdb.header['imtype'] = ('bias', 'type of image')
            if flag == 0:
                print('all biases will be trimmed to :', ccdb.meta['trimsec'])
                flag = 1
            tbias = ccdp.trim_image(ccdb, fits_section=str(ccdb.meta['trimsec']))
            tbias.meta['imtype'] = ('trimmed bias', 'type of image')
            tbias.meta['taxis1'] = (2048, 'dimension1')
            tbias.meta['taxis2'] = (4096, 'dimension2')
            tbias.write('Trimmed_Bias/' + biasn[0:8] + '_trim.fits', overwrite=True)
    return biascollection


def bias_combine(refresh=0):
    tbiascollection = ImageFileCollection('Trimmed_Bias')
    print('found', len(tbiascollection.values('file')), 'trimmed biases')
    start=time.time()
    combined_bias = ccdp.combine(tbiascollection.files_filtered(imtype='trimmed bias', include_path=True),
                                 method='median')
    print(time.time()-start)
    combined_bias.meta['combined'] = True
    combined_bias.write('Master_Files/mbias_median.fits', overwrite=True)
    return tbiascollection


create_directories()
biascollection = trim_bias(1)
tbiascollection = bias_combine(1)
