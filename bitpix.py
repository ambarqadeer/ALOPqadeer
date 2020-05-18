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


sci=CCDData.read('HD115709/SII/r1382286.fit',unit='adu')
print(sci.header['bitpix'])
