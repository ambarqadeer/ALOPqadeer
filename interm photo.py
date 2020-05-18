import sp
from astropy.io import ascii
import time
import ccdproc as ccdp
from astropy.nddata import CCDData
import numpy as np
import glob
import itertools
from astropy.table import Table


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
        x,y,flux,sh,ro=sp.find(imsd,hmin=4000,fwhm=17,sharplim=([-2,2]),roundlim=([-2,2]))
        for i in range(len(x)):
            trackerf.append([imsn,x[i],y[i],sh[i],ro[i]])
    except Exception as e:
        print(imsn,'did not meet one of the criteria')
        print(e)
    try:
        flux,errap,sky,skerr=sp.aper(imsd,x,y,phpadu=2.9,apr=[10],skyrad=[15,20],flux=True)
        trackera.append([imsn,flux[0][0],flux[1][0],sky[0],sky[1]])
    except Exception as e:
        print(imsn,'returned an error:')
        print(e)

tf=np.reshape(trackerf,(-1,5))
ta=np.reshape(trackera,(-1,5))
print(time.time()-start)
foundf=Table(tf,names=('fname','x','y','sh','ro'))
founda=Table(ta,names=('fname','fluxcalib','fluxtarget','skycalib','skytarget'))
ascii.write(foundf,'stars.csv',format='csv',overwrite=True)
ascii.write(founda,'fluxes.csv',format='csv',overwrite=True)

'''manually removed 2582'''

