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


flatlist=ImageFileCollection('Trimmed_Flat')
print(flatlist.summary)

for file in flatlist.data():
    print(file.mean())

a=flatlist.summary['file'][0]
b=flatlist.summary['file'][4]

ccd1=CCDData.read('Trimmed_Flat/'+a,unit='adu')
ccd2=CCDData.read('Trimmed_Flat/'+b,unit='adu')
ccd=ccd1.divide(ccd2)
print(ccd.data.mean())
flatt=ccdp.trim_image(ccd,fits_section='[235:1564,1046:2509]')

import itertools
#
# a=[]
#
# for flatn,flat in enumerate(itertools.islice(flatlist.ccds(ccd_kwargs={'unit':'adu'}),2)):
#     flatt=ccdp.trim_image(flat,fits_section='[235:1564,1046:2509]')
#     a.append(flatt)
#
# #print(a.shape)
# ratio=a[0].data/a[1].data
# print(ratio.shape)
mask=ccdp.ccdmask(flatt,findbadcolumns=True)
# a=np.argwhere(mask==1)
# print(a)
# print(type(mask))
#
maski=np.multiply(mask,1)
hdu=fits.PrimaryHDU(maski)
hdul=fits.HDUList([hdu])
hdul.writeto('mask.fits',overwrite=True)



# plt.figure(figsize=(30,15))
# fig, axs=plt.subplots(1,2)
#
#
# axs[0].imshow(a[0].data,origin='lower')
# axs[1].imshow(mask,origin='lower')
# fig.show()