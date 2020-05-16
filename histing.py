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


a = CCDData.read('Master_Files/mbias_median.fits', unit='adu')
b = CCDData.read('Master_Files/mbias.fits', unit='adu')
median = np.ndarray.flatten(a.data)
sigclip = np.ndarray.flatten(b.data)


med=sigma_clip(median,sigma=3,cenfunc='median',stdfunc=mad_std,masked=False)
sig=sigma_clip(sigclip,sigma=3,cenfunc='median',stdfunc=mad_std,masked=False)


plt.figure(figsize=(10,20))
fig, axs=plt.subplots(2,1)


axs[0].hist(med,bins=20,label='median',histtype='stepfilled',log=False)
axs[1].hist(sig,bins=50,label='sigclip',histtype='stepfilled',log=False)
fig.show()

'''
hist(s, bins='freedman', label='average', alpha=0.5)
plt.show()
'''