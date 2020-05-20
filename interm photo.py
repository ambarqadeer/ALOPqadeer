import sp
from astropy.io import ascii
import time
import ccdproc as ccdp
from astropy.nddata import CCDData
import numpy as np
import glob
import itertools
from astropy.table import Table
from astropy import table
import Logs
from astropy.convolution import convolve, Box1DKernel
import matplotlib.pyplot as plt


science=ccdp.ImageFileCollection('Final_Science/median')

# sciencelist=glob.glob('Corrected_Science/r*fit')
#ascii.write(sc,'sc_summary.csv',format='csv',overwrite=True)
#found=Table(names=('fname','x','y','sh','ro'),dtype=('str','float32','float32','float32','float32'))
trackerf=[]
trackera=[]
start=time.time()
for imsd,imsn in science.data(return_fname=True):
    print('assessing',imsn)
    try:
        x,y,flux,sh,ro=sp.find(imsd,hmin=4000,fwhm=20,sharplim=([-2,2]),roundlim=([-2,2]))
        for i in range(len(x)):
            trackerf.append([imsn,x[i],y[i],sh[i],ro[i]])
    except Exception as e:
        print(imsn,'did not meet one of the criteria')
        print(e)
    try:
        flux,errap,sky,skerr=sp.aper(imsd,x,y,phpadu=2.9,apr=[5],skyrad=[15,20],flux=True)
        trackera.append([imsn,flux[0][0],flux[1][0],sky[0],sky[1]])
    except Exception as e:
        print(imsn,'returned an error:')
        print(e)

tf=np.reshape(trackerf,(-1,5))
ta=np.reshape(trackera,(-1,5))
print(time.time()-start)
foundf=Table(tf,names=('fname','x','y','sh','ro'),dtype=('str','float64','float64','float64','float64'))
fluxt=Table(ta,names=('fname','fluxcalib','fluxtarget','skycalib','skytarget')
            ,dtype=('str','float64','float64','float64','float64'))
ascii.write(foundf,'stars.csv',format='csv',overwrite=True)
ascii.write(fluxt,'fluxes.csv',format='csv',overwrite=True)

'''manually removed 2582'''


log1=Logs.parse_logfile('HD115709/log_files/run_log_20180428.int')
log2=Logs.parse_logfile('HD115709/log_files/run_log_20180430.int')

log1t=Table(list(log1['Observations'].values()),names=tuple(log1['Observations'].keys()))

log2t=Table(list(log2['Observations'].values()),names=tuple(log2['Observations'].keys()))


# logc=table.vstack([log1t,log2t])




# fluxt=ascii.read('fluxes.csv')

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


final1['relflux']=final1['diff']/np.nanmedian(final1['diff'])
final2['relflux']=final2['diff']/np.nanmedian(final2['diff'])




final1['smadiff']=convolve(list(final1['relflux']),kernel=Box1DKernel(10),boundary='extend')
final2['smadiff']=convolve(list(final2['relflux']),kernel=Box1DKernel(10),boundary='extend')


ascii.write(final1,'final1.csv',overwrite=True,format='csv')
ascii.write(final2,'final2.csv',overwrite=True,format='csv')


# print(final)


fig, axs=plt.subplots(1,2,figsize=(40,10))
axs[0].scatter(final1['MJD'],final1['smadiff'])
axs[1].scatter(final2['MJD'],final2['smadiff'])
plt.show()