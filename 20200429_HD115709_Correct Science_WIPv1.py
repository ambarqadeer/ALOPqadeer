# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 13:20:36 2020

@author: ambar
"""

'''
Imports
'''

import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
import astropy.units as u
import ccdproc
import os

##############################################################

'''
Function Definitions
'''

# function creating a master bias and saving it as 'tmaster_bias.fits '
def createmasterbias():
    #create biaslist
    biaslist=sorted(glob.glob('HD115709/bias/r*.fit'))
    print('Number of biases used =', len(biaslist))

    
    #open bias zero for cube dimensions
    hdub=fits.open(biaslist[0])
    bias0=hdub[4].data
    print(biaslist[0],'is open, shape:',bias0.shape)
    hdub.close()
    print(biaslist[0],'is closed')
    
    #create biascube with shape of science and len(biaslist)
    biascube=np.zeros((science0.shape[0],science0.shape[1],len(biaslist)),dtype=bias0.dtype)
    print('biascube created with shape :', biascube.shape)
    
    
    #crop and stack biases (biasc= counter, and biasn = name)
    for biasc, biasn in enumerate(biaslist):
        print('Open :',biasn)
        hdu=fits.open(biaslist[biasc])
        if hdu[4].header['CHIPNAME']=='A5382-1-7':
            bias=CCDData(hdu[4].data,unit=u.adu)
            tbias=ccdproc.trim_image(bias,fits_section=sciencewindow)
            biascube[:,:,biasc]=tbias.data
            hdu.close()
        else:
            hdu.close()
            print(biasn, 'returned wrong chipname for selected extension')
            exit
    
    #take median and write to disk
    mbias=np.nanmedian(biascube,2)
    ccdmbias=CCDData(mbias,unit=u.adu)
    ccdmbias.write('tmaster_bias.fits',overwrite=True)
    print('master bias shape :',ccdmbias.data.shape)


# function creating a master flat and saving it as 'tmaster_flat.fits '
# note that this is not normalized, the finction used to correct the science normalizes internally
def createmasterflat():
    
    #create flatlist
    flatlist=sorted(glob.glob('HD115709/flat_SII/r*.fit'))
    print('Number of flats used =', len(flatlist))
    
    #open flat0 for cube len
    hduf=fits.open(flatlist[0])
    flat0=hduf[4].data
    print(flatlist[0],'is open, shape:', flat0.shape)
    hduf.close()
    print(flatlist[0],'is closed')
    
    #create flatcube with shape of science and len(flatlist)
    flatcube=np.zeros((science0.shape[0],science0.shape[1],len(flatlist)),dtype=flat0.dtype)
    
    #convert trim and populate flatcube after bias correction
    for flatc, flatn in enumerate(flatlist):
        print('Open :', flatn)
        hdu=fits.open(flatlist[flatc])
        if hdu[4].header['CHIPNAME']=='A5382-1-7':
            ccdflat=CCDData(hdu[4].data,unit=u.adu)
            ccdtflat=ccdproc.trim_image(ccdflat,fits_section=sciencewindow)
            ccdcflat=ccdproc.subtract_bias(ccdtflat,ccdmbias_use)
            flatcube[:,:,flatc]=ccdcflat.data
            hdu.close()
        else:
            hdu.close()
            print(flatn, 'returned wrong chipname for selected extension')
            exit
    #take median and write to disk
    mflat=np.nanmedian(flatcube,2)
    ccdmflat=CCDData(mflat,unit=u.adu)
    ccdmflat.write('tmaster_flat.fits', overwrite=True) # write the fits to disk

#Corrects science images and saves them to the folder 'Corrected_Science'
#MAke sure the folder exists
def correctscience():
    for sciencec, sciencen in enumerate(sciencelist):
        hdu=fits.open(sciencelist[sciencec])
        print('open', sciencen)
        if hdu[1].header['CHIPNAME']=='A5382-1-7':
            ccdscience=CCDData(hdu[1].data,unit=u.adu)
            if not os.path.exists('Corrected_Science'):
                os.makedirs('Corrected_Science')
            cccdscience=ccdproc.flat_correct(ccdproc.subtract_bias(ccdscience,ccdmbias_use),ccdmflat_use,norm_value=np.nanmedian(ccdmflat_use))
            path='Corrected_Science/'+sciencen[-12:]
            cccdscience.write(path,overwrite=True)
            hdu.close()
        else:
            hdu.close()
            print(sciencen, 'returned wrong chipname for selected extension')
            exit
##############################################################

'''
Main
'''

#open 1 science frame to get shape
sciencelist=sorted(glob.glob('HD115709/SII/r*.fit'))
hdus=fits.open(sciencelist[0])
science0=hdus[1].data
print(sciencelist[0],'is open, shape:',science0.shape)
#sciencewindow=hdus[1].header['RTDATSEC']
sciencewindow='[288:1617,1046:2509]'
print(sciencewindow)
hdus.close()
print(sciencelist[0],'is closed')



#check for existing master bias and load into 'ccdmbias_use', if it doesn't exist, make
try:
    #get master bias data
    hdumb=fits.open('tmaster_bias.fits')
    ccdmbias_use=CCDData(hdumb[0].data,unit=u.adu)
    print('tmaster_bias.fits is open, shape:',ccdmbias_use.shape)
    hdumb.close()
    print('tmaster_bias.fits is closed')
    print('bias loaded in cdmbias_use')

except:
    createmasterbias()
    hdumb=fits.open('tmaster_bias.fits')
    ccdmbias_use=CCDData(hdumb[0].data,unit=u.adu)
    print('tmaster_bias.fits is open, shape:',ccdmbias_use.shape)
    hdumb.close()
    print('tmaster_bias.fits is closed')
    print('bias loaded in cdmbias_use')

#check for existing master flat and load into 'ccdmflat_use', if it doesn't exist, make
try:
    #get master flat data
    hdumf=fits.open('tmaster_flat.fits')
    ccdmflat_use=CCDData(hdumf[0].data,unit=u.adu)
    print('tmaster_flat.fits is open,shape:',ccdmflat_use.shape)
    hdumf.close()
    print('tmaster_flat.fits is closed')
    print('flat loaded in ccdmflat_use')
except:
    createmasterflat()
    hdumf=fits.open('tmaster_flat.fits')
    ccdmflat_use=CCDData(hdumf[0].data,unit=u.adu)
    print('tmaster_flat.fits is open,shape:',ccdmflat_use.shape)
    hdumf.close()
    print('tmaster_flat.fits is closed')
    print('flat loaded in ccdmflat_use')

print(np.nanmedian(ccdmflat_use))
correctscience()
