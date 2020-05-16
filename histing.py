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

tbiaslist = sorted(glob.glob('Trimmed_Bias/r*.fits'))
print(len(tbiaslist))
a = CCDData.read(tbiaslist[0], unit='adu')
print(tbiaslist[0][-18:-10])
f = np.ndarray.flatten(a.data)
print(np.nanmean(f[1000000:-1000000]))
s=sigma_clip(f,cenfunc='mean',masked=False)


hist(s,bins='freedman')
plt.show()

'''
hist(s, bins='freedman', label='average', alpha=0.5)
plt.show()
'''