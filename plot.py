import sp
from astropy.io import ascii
import matplotlib.pyplot as plt
import time
import ccdproc as ccdp
from astropy.nddata import CCDData
import numpy as np
import glob
import itertools
from astropy.table import Table
import astropy.table as table
from astropy.convolution import convolve, Box1DKernel


import Logs
log1=Logs.parse_logfile('HD115709/log_files/run_log_20180428.int')
log2=Logs.parse_logfile('HD115709/log_files/run_log_20180430.int')

log1t=Table(list(log1['Observations'].values()),names=tuple(log1['Observations'].keys()))

log2t=Table(list(log2['Observations'].values()),names=tuple(log2['Observations'].keys()))


# logc=table.vstack([log1t,log2t])




fluxt=ascii.read('fluxes.csv')

# print(fluxt)
c=[]
for i, val in enumerate(fluxt['fname']):
    c.append(int(val[1:8]))

fluxt['Run']=c

final1=table.join(fluxt,log1t,keys='Run')
final2=table.join(fluxt,log2t,keys='Run')


final1['truefluxcalib']=final1['fluxcalib']/final1['Exptime']
final2['truefluxcalib']=final2['fluxcalib']/final2['Exptime']

final1['truefluxtarget']=final1['fluxtarget']/final1['Exptime']
final2['truefluxtarget']=final2['fluxtarget']/final2['Exptime']

final1['diff']=final1['truefluxcalib']-final1['truefluxtarget']
final2['diff']=final2['truefluxcalib']-final2['truefluxtarget']


# final1['smadiff']=convolve(list(final1['truefluxcalib']),kernel=Box1DKernel(10),boundary='extend')
# final2['smadiff']=convolve(list(final2['truefluxcalib']),kernel=Box1DKernel(10),boundary='extend')


ascii.write(final1,'final1.csv',overwrite=True,format='csv')
ascii.write(final2,'final2.csv',overwrite=True,format='csv')


# print(final)


fig, axs=plt.subplots(1,2,figsize=(40,10))
axs[0].scatter(final1['MJD'],final1['truefluxtarget']/1.15e6)
axs[1].scatter(final2['MJD'],final2['truefluxtarget']/1.15e6)
plt.show()
