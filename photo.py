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
tracker=[]
start=time.time()
for imsd,imsn in science.data(return_fname=True):
    print('assessing',imsn)
    try:
        x,y,flux,sh,ro=sp.find(imsd,hmin=4000,fwhm=15,sharplim=([-2,2]),roundlim=([-2,2]))
        for i in range(len(x)):
            tracker.append([imsn,x[i],y[i],sh[i],ro[i]])
    except Exception as e:
        print(imsn,'did not meet one of the criteria')
        print(e)
t=np.reshape(tracker,(-1,5))
print(time.time()-start)
found=Table(t,names=('fname','x','y','sh','ro'))
ascii.write(found,'stars.csv',format='csv',overwrite=True)


'''
stars=ascii.read('stars.csv')
print(len(stars))
files=science.summary['file']
c=[]
for i in files:
    c.append(list(stars['fname']).count(i))
print(len(c))
counts=Table([files,c],names=('files','starcount'))
# ascii.write(counts,'starscounts.csv',format='csv',overwrite=True)
cmask=counts['starcount'] != 2
mfiles=counts['files'][cmask]
starmask=[]
for i in stars['fname']:
    starmask.append(i in mfiles)

exstars=stars[starmask]
ascii.write(exstars,'exstars.csv',format='csv',overwrite=True)

rex=ascii.read('exstars.csv')
uex= np.unique(rex['fname'])


tracker2=[]
start=time.time()
for u in uex:
    ccd=CCDData.read('Final_Science/median/'+u)
    try:
        print('asessing',u)
        x, y, flux, sh, ro = sp.find(ccd.data, hmin=4000, fwhm=17, sharplim=([-2, 2]), roundlim=([-2, 2]))
        for i in range(len(x)):
            tracker2.append([u, x[i], y[i], sh[i], ro[i]])
    except Exception as e:
        print(u, 'did not meet one of the criteria')
        print(e)
t2 = np.reshape(tracker2, (-1, 5))
print(time.time()-start)
found2=Table(t2,names=('fname','x','y','sh','ro'))
ascii.write(found2,'find2.csv',format='csv',overwrite=True)
'''




