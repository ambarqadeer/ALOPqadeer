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


if not os.path.exists('Master_Files/Windowed'):
    os.makedirs('Master_Files/Windowed')
    print('folder \'Master_Files/Windowed\' created')
if not os.path.exists('Final_Science'):
    os.makedirs('Final_Science')
    print('folder \'Final_Science\' created')

#read science images and find RTDATSEC
sciencecollection=ImageFileCollection('HD115709/SII',ext=1)
for sci in sciencecollection.ccds(ccd_kwargs={'unit': 'adu'}):
    sciwin=sci.meta['RTDATSEC']
    break

cmpath='Master_Files/Windowed/'
cspath='Final_Science/'

#crop master biases and master flats
def crop_masters(path=cmpath):
    mastercollection=ImageFileCollection('Master_Files')
    for image, imname in mastercollection.ccds(imtype='trimmed bias',return_fname=True):
        trimage=ccdp.trim_image(image,fits_section=str(sciwin))
        trimage.meta['trimwind']=(str(sciwin),'readout window')
        trimage.meta['imtype'] = ('mbias', 'windowed master bias')
        trimage.write(path+imname,overwrite=True)
    for image, imname in mastercollection.ccds(imtype='subflat',return_fname=True):
        trimage=ccdp.trim_image(image,fits_section=str(sciwin))
        trimage.meta['trimwind']=(str(sciwin),'readout window')
        trimage.meta['imtype'] = ('mflat', 'windowed master flat')
        trimage.write(path+imname,overwrite=True)

tmastercollection=ImageFileCollection(cmpath)


'''check shape'''
#check shape
# print(sci.data.shape)
# for master in tmastercollection.ccds():
#     print(master.data.shape)


crop_masters()


for mbias, mbiasn in tmastercollection.ccds(imtype='mbias', combined='sigma_clip average',return_fname=True):
    print('using',mbiasn)
    for mflat, mflatn in tmastercollection.ccds(imtype='mflat', flatcom='sigma', return_fname=True):
        print('using', mflatn)
        start = time.time()
        for sci, scin in sciencecollection.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            print('correcting',scin)
            subsci=ccdp.subtract_bias(sci,mbias,add_keyword={'subsci':'sigmbias'})
            corsci=ccdp.flat_correct(subsci,mflat,norm_value=float(mflat.header['normmed']))
            corsci.write(cspath+scin,overwrite=True)

print (time.time()-start)