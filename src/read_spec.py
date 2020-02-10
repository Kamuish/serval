#! /usr/bin/env python
__author__ = 'Mathias Zechmeister'
__version__ = '2019-03-01'

import glob
import gzip
import os
import sys
import tarfile
import warnings
from collections import namedtuple

import time
import astropy.io.fits as pyfits

#import fitsio
import numpy as np

from src.utils import pause, stop
import src.brv_we14py as brv_we14py
import src.sunrise as sunrise
from src.utils import FitsClass 

import logging 
logger = logging.getLogger(__name__)
class nameddict(dict):
   """
   Examples
   --------
   >>> nameddict({'a':1, 'b':2})
   {'a': 1, 'b': 2}
   >>> x = nameddict(a=1, b=2)
   >>> x.a
   1
   >>> x['a']
   1
   >>> x.translate(3)
   ['a', 'b']

   """
   __getattr__ = dict.__getitem__

   def translate(self, x):
      return [name for name,f in list(self.items()) if (f & x) or f==x==0]

# bpmap flagging
flag = nameddict(
   ok=       0, # good pixel
   nan=      1, # nan flux in pixel of spectrum
   neg=      2, # significant negative flux in pixel of spectrum (f < -3*f_err < 0 && )
   sat=      4, # too high flux (saturated)
   atm=      8, # telluric at wavelength of spectrum
   sky=     16, # sky emission at wavelength of spectrum
   out=     32, # region outside the template
   clip=    64, # clipped value
   lowQ=   128, # low Q in stellar spectrum (mask continuum)
   badT=   256, # bad corresponding region in the template spectrum
)

# spectrum flagging
sflag = nameddict(
   nosci=    1, # nosci frame, e.g. calibration files
   iod=      2,
   config=   4, # bad configuration/setup (e.g. HARPS eggs, pre, post)
   dist=    16, # coordinates too much off
   daytime= 32, # not within nautical twilight
   lowSN=   64, # too low S/N
   hiSN=   128, # too high S/N
   led=    256, # LED on during observation (CARM_VIS)
   rvnan=  512
)


# bpmap flags
flag_cosm = flag.sat  # @ FEROS for now use same flag as sat
def_wlog = True
brvrefs = ['DRS', 'MH', 'WEhtml', 'WEidl', 'WE']

class Spectrum:
   """
   Compiles information from frame (filesname), instrument and target.
   and provides method to read the data.

   Examples
   --------
   >>> from read_spec import *
   >>> x = Spectrum('/media/5748244b-c6d4-49bb-97d4-1ae0a8ba6aed/data/HARPS/red_DRS/gj551/HARPS.2013-05-06T01:54:55.600.tar', fib='A',orders=-1,drs=1)
   >>> x = Spectrum('/home/zechmeister/data/harps/gj588/data/HARPS.2013-05-11T08:46:58.520_e2ds_A.fits')

   """
   brvref = 'MH' # can be overwritten by serval.py depending on user input
   def __init__(self, filename, inst='HARPS', pfits=True, orders=None, wlog=def_wlog, drs=False, fib=None, targ=None, verb=False):
      """
      pfits : fits reader
             True: pyfits (safe, e.g. to get the full header of the start refence spectrum)
             2: my own simple pure python fits reader, faster than pyfits
             It was developed for HARPS e2ds, but the advent ADP.*.tar was difficult.

      """
      self.filename = filename
      self.drsname = 'DRS'
      self.drs = drs
      self.flag = 0
      self.bflag = 0
      self.tmmean = 0
      self.f = None
      self.w = None
      self.e = None
      self.bpmap = None
      self.mod = None
      self.fib = fib
      self.fo = None
      self.wo = None
      self.eo = None
      self.bo = None
      self.drift = np.nan
      self.e_drift = np.nan
      self.utc = None
      self.ra = None
      self.de = None
      self.obs = type('specdata', (object,), dict(lat=None, lon=None))
      self.airmass = np.nan
      self.inst = inst
      self.scan = self.inst.scan
      self.data = self.inst.data

      if '.gz' in filename: 
         pfits=True

      self.ccf = type('ccf',(), dict(rvc=np.nan, err_rvc=np.nan, bis=np.nan, fwhm=np.nan, contrast=np.nan, mask=0, header=0))

      # scan fits header for times, modes, snr, etc.
      self.scan(self, filename, pfits=pfits)

      if verb:
         logger.info("scan %s:"%self.instname, self.timeid, self.header['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode)

      if inst.name != self.instname:
         pause('WARNING:', filename, 'from', self.inst, ', but mode is', inst)
      self.obj = self.header['OBJECT']

      ### Barycentric correction ###
      # start and end computed only for WE and WEhtml
      self.bjd_sart, self.berv_start = np.nan, np.nan
      self.bjd_end, self.berv_end = np.nan, np.nan
      if targ and targ.name == 'cal':
         self.bjd, self.berv = self.drsbjd, 0.
      elif targ and targ.ra:  # unique coordinates
         if self.brvref == 'MH':
            # fastest version
            #sys.path.append(os.environ['HOME']+'/programs/BarCor/')
            sys.path.insert(1, sys.path[0]+os.sep+'BarCor')            # bary now in src/BarCor
            import bary
            self.bjd, self.berv = bary.bary(self.dateobs, targ.ra, targ.de, inst.name, epoch=2000, exptime=self.exptime*2* self.tmmean, pma=targ.pmra, pmd=targ.pmde)
         if self.brvref in ('WEhtml', 'WEidl', 'WE'):
            # Wright & Eastman (2014) via online or idl request
            # cd /home/raid0/zechmeister/idl/exofast/bary
            # export ASTRO_DATA=/home/raid0/zechmeister/
            #idl -e 'print, bjd2utc(2457395.24563d, 4.58559072, 44.02195596), fo="(D)"'
            jd_utc = [self.mjd + 2400000.5 + self.exptime*self.tmmean/24./3600]
            jd_utcs = [self.mjd + 2400000.5, jd_utc[0], self.mjd + 2400000.5 + self.exptime/24./3600]
            ra = (targ.ra[0] + targ.ra[1]/60. + targ.ra[2]/3600.) * 15  # [deg]
            de = (targ.de[0] + np.copysign(targ.de[1]/60. + targ.de[2]/3600., targ.de[0]))       # [deg]
            obsname = inst.obsname #{'CARM_VIS':'ca', 'CARM_NIR':'ca', 'FEROS':'eso', 'HARPS':'eso', 'HARPN':'lapalma', 'HPF':'hpf'}[inst]
            obsloc = inst.obsloc if hasattr(inst, 'obsloc') else {}
            if self.brvref == 'WE':
               # pure python version
               #self.bjd, self.berv = brv_we14py.bjdbrv(jd_utc=jd_utc[0], ra=ra, dec=de, obsname=obsname, pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0])
               (_, self.bjd, _), (self.berv_start, self.berv, self.berv_end) = brv_we14py.bjdbrv(jd_utc=jd_utcs, ra=ra, dec=de, obsname=obsname, pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0], **obsloc)
            elif self.brvref == 'WEhtml':
               self.bjd = brv_we14html.utc2bjd(jd_utc=jd_utc, ra=ra, dec=de)
               #self.berv = brv_we14html.bvc(jd_utc=jd_utc, ra="%s+%s+%s"%targ.ra, dec="%s+%s+%s"%targ.de, obsname='ca', pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0], raunits='hours', deunits='degrees')[0]
               self.berv_start, self.berv, self.berv_end = brv_we14html.bvc(jd_utc=jd_utcs, ra="%s+%s+%s"%targ.ra, dec="%s+%s+%s"%targ.de, obsname='ca', pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0], raunits='hours', deunits='degrees')
            else:
               self.bjd, self.berv = bjdbrv(jd_utc=jd_utc[0], ra=ra, dec=de, obsname=obsname, pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0])

            self.berv /= 1000.   # m/s to km/s
         if self.brvref == 'DRS':
            self.bjd, self.berv = self.drsbjd, self.drsberv
            self.brvref = self.drsname
      else:
         self.bjd, self.berv = self.drsbjd, self.drsberv
         self.brvref = self.drsname
      if self.fib == 'B':
         self.berv = np.nan
         self.drsberv = np.nan

      self.header['HIERARCH SERVAL BREF'] = (self.brvref, 'Barycentric code')
      self.header['HIERARCH SERVAL BJD'] = (self.bjd, 'Barycentric Julian Day')
      self.header['HIERARCH SERVAL BERV'] = (self.berv, '[km/s] Barycentric correction')


      if self.utc:
         date = self.utc.year, self.utc.month, self.utc.day
         sunris = sunrise.sun(*date, lon=self.obs.lon, lat=self.obs.lat)
         sunset = sunrise.sun(*date, lon=self.obs.lon, lat=self.obs.lat, rise=False)
         ut = self.utc.hour + self.utc.minute/60. +  self.utc.second/3600.
         utend = ut + self.exptime/3600.
         dark = ((sunset < ut <= utend < sunris) or   # |****S__UT__R***|
                 (sunris < sunset < ut <= utend) or   # |__R*******S__UT|
                 (utend < sunris < sunset < ut)  or   # |T__R********S_U|
                 (ut < utend < sunris < sunset))      # |UT__R********S_|
         #if not dark:
            #self.flag |= sflag.daytime
      if orders is not None:
         self.read_data(orders=orders, wlog=wlog)

   def __get_item__(self, order):
      # spectrum needs to be dict like
      return self.get_data(order)

   def get_data(self, orders=np.s_[:], wlog=def_wlog, verb=False, **kwargs):
      """Returns only data."""
      o = orders
      if self.w:
         w, f, e, b = self.w[o], self.f[o], self.e[o], self.bpmap[o]
      else:
         w, f, e, b = self.data(self, orders=orders, **kwargs)
         w = np.log(w) if wlog else w.astype(np.float)
         f = f.astype(float)
         e = e.astype(float)
         self.bflag |= np.bitwise_or.reduce(b.ravel())
      return type('specdata', (object,),
               dict(w=w, f=f, e=e, bpmap=b, berv=self.berv, o=o))

   def read_data(self, verb=False, **kwargs):
      """Read only data, no header, store it as an attribute."""
      data = self.get_data(**kwargs)
      if isinstance(kwargs.get('orders'), int):
         self.wo, self.fo, self.eo, self.bo = data.w, data.f, data.e, data.bpmap
      else:
         self.w, self.f, self.e, self.bpmap = data.w, data.f, data.e, data.bpmap


class Inst:
   def __init__(self, inst):
      pass

tarmode = 5   # FIXME: temporary chaging the TARMODE
# 0  - extract physically (works for all, slow)
#      tarfits directory not cleaned
# 1  - treat tarfile as normal file which has offset, seek() and open()
#      but no size limit function, therefore not supported by pyfits
# 2  - as 1, but with try to overload build-in open (not working for both)
# 3  - tar.extractfile (use its build-in functions for read, offset and size, slower read!?) not supported by pyfits
# 5  - use tar.extractfile, works with pyfits, doesn't work with myfits (np.fromfile accepts only real files and files object, not tar)
# 6  - open as normal file, does not work

def file_from_tar(file_name, inst='HARPS', fib=None, **kwargs):
   """
   Returns a file-like object.
   Reading header and data should use the same fits reader

   >>> s = '/home/astro115/carmenes/data/HARPS/DRS/gj876/ADP.2014-12-07T19:00:30.953.tar'
   >>> tar = tarfile.open(s)
   >>> extr = tar.getmember('data/reduced/2014-12-06/HARPS.2014-12-07T00:35:29.894_e2ds_A.fits')

   """
   pat = {'HARPS': {'A': '_e2ds_A.fits', 'B': '_e2ds_B.fits'}[fib],
          'FEROS': {'A': '.1061.fits', 'B': '.1062.fits'}[fib]} [inst]

   keep_open = False
   tar =  tarfile.open(file_name)
   for member in tar.getmembers():
      if pat in member.name: 
         extr = member

   if tarmode in (1,2,3) and kwargs.get('pfits') == 2:
      # We could use tar.extractfile(extr) but this requires reopen the tar file
      # and does not support "with open()".
      # Instead we will use offset, seek() and open().
      # but both seems to be slower than physical extraction
      extr.mother = s
      s = extr
   elif tarmode in (0,5) and kwargs.get('pfits') == 2:
      # Open the tar file as a normal file and store the position and size of the fits file
      # as attributes.

      s = FitsClass(filepath = file_name, 
                  offset_data = extr.offset_data,
                  name = extr.name,
                  size=extr.size,
                  mode = 1
      )
      # actually creating those new properties

   elif tarmode == 5:
      # does not work with np.fromfile
      extracted_member  = tar.extractfile(extr)
      s = FitsClass( filepath = extracted_member, name = extr.name, mode = 2, **{'tar_file': tar})
      keep_open = True
   else:
      tar.extract(extr, path='tarfits')   # extract physically
      s = 'tarfits/'+extr.name
   
   if not keep_open:
      tar.close()
   return s


class imhead(dict):
   """
   Returns fitsheader as a dict.
   
   imhead runs faster than pyfits. It does less checks, but might be
   more error prone.

   Keyword arguments:
   filename -- filename of the fitsfile

   Returns
   -------
   imhead: dict of values and dict of comments
   EXTHDRSZ: size of the header [bytes]

   Examples
   --------
   >>> x = imhead(filename)
   >>> x['OBJECT']
   >>> x.comments['OBJECT']

   """
   def __init__(self, s, *args, **kwargs):
      hdr = {}
      self.comments = {}
      extpos = kwargs.get('extpos', 0)  # position of the extension with the fitsfile
      count = kwargs.get('count', -1) 

      #if not len(args):
         #args = ''   # read all
         #stop( '\n\nimhead: there must be are args; reading all keys not yet supported\n')
      args += ('NAXIS', 'BITPIX')
      NR = 0
      if isinstance(s, (tarfile.TarInfo, tarfile.ExFileObject)):
         if tarmode==1:
            fi = open(s.mother)
         if tarmode==3:
            fi = tarfile.open(s.mother).extractfile(s)
      elif isinstance(s, (FitsClass)):   # tarmode == 5,6  # FIXME: see if it's using the FITSClass; previously checked for tarfiles and file
         fi = s
      elif s.endswith('.gz'):   # normal file
         fi = gzip.open(s)
      else:   # normal file
         fi = open(s)

      if not extpos and hasattr(s, 'offset_data'):
         extpos = s.offset_data
      #pause()

      
      if isinstance(fi, FitsClass):  # FIXME: temporary solution to this problem; do the other file types work properly? 

         args = tuple((str.encode(i) for i in args))
         with fi as file:
            file.seek(extpos)
            for card in iter(lambda:file.read(80), ''):   # read in 80 byte blocks 
               # Header files have the header with 80 bytes blocks -> this reads the header line by line
               NR += 1
               if card.startswith(b'END'): 
                  break
               if card.startswith(args):
                  #key, val, comment = card.replace("= ","/ ",1).split("/ ")
                  key, val, comment = card.replace(b' /',b'= ',1).split(b'= ')   # does not parse HARPN.2012-09-03T05-29-56.622 "HIERARCH TNG OBS TARG NAME = 'M1'/ no comment"
                  # key, val, comment = card.replace('/','= ',1).split('= ')
                  hdr[key.strip().decode('utf-8')] = val.strip(b"' ").decode('utf-8') if b"'" in val else float(val) if b'.' in val else int(val)
                  #hdr[key.strip()] = val.strip("' ") if any(i in val for i in "'TF") else float(val) if '.' in val and not 'NaN.' in val else int(val) if val.strip(" -").isdigit() else val
                  self.comments[key.strip().decode('utf-8')] = comment.decode('utf-8')
                  count -= 1
                  if count==0: 
                     args = () # all found; do not check cards anymore; only read to end
               
            
      else:
         logger.fatal("I should implement this")
         raise NotImplementedError

      
      hsz = 2880 * ((NR-1)//36 + 1)

      dsz = 0     # EXTDATSZ
      self.NAXIS = hdr.get('NAXIS', 0)
      if self.NAXIS:
         self.BITPIX = hdr['BITPIX']
         self.NAXIS1 = hdr['NAXIS1']
         dsz = abs(self.BITPIX)   # data size
         dsz *= self.NAXIS1
         if self.NAXIS > 1:
            self.NAXIS2 = hdr['NAXIS2']
            dsz *= self.NAXIS2
         dsz = ((dsz//8-1)//2880 +1) * 2880

      self.EXTHDRSZ = hsz   # add hdr size
      self.EXTDATA = extpos + hsz
      self.EXTEND = extpos + hsz + dsz
      super(imhead, self).__init__(hdr)


class getext(dict):
   """
   Fast reading of single orders.

   Uses imhead to scan the extension headers and to get
   the data types and sizes.

   Examples
   --------
   >>> hdu = getext(filename)
   >>> hdu['SPEC']
   >>> hdu.getdata('SPEC', o=13)

   """
   def __init__(self, s):
      """
      Extracts from header for each extension the data position, type, and shape
      and saves the extension information as dictionary.

      """
      self.fileobj = s
      ext = {}
      filepos = 0
      fileend = -1
      i = 0
      if isinstance(s, tarfile.TarInfo) and s.tarmode:
         self.fileobj = s.mother
      if hasattr(s, 'offset_data'):
         filepos = s.offset_data
         fileend = s.offset_data + s.size 


      while True:
         exthdr = imhead(self.fileobj, 'EXTNAME', extpos=filepos)
         if exthdr.EXTEND <= filepos:
            print("no new records")
            break   # no new records
         filepos = exthdr.EXTEND
         ext[i]  = exthdr

         if filepos == fileend: 
            break   # end of filemember in tarfile
         i += 1
      super(getext, self).__init__(ext)

   def getdata(self, extname=None, order=np.s_[:]):
      if extname is None:  # search first extension with data
         extname = 0
         while not self[extname].NAXIS: extname += 1
      ext = self[extname]
      dtype = {-32: '>f4', -64: '>f8'}[ext.BITPIX]
      dsize = {-32: 4, -64: 8}[ext.BITPIX]
      was_open = True
      if isinstance(self.fileobj, (tarfile.ExFileObject, FitsClass)):
         funit = self.fileobj
      else:
         funit = self.fileobj
         print(self.fileobj, type(self.fileobj))

      with self.fileobj as funit:
         if order == np.s_[:]: # read all
            funit.seek(ext.EXTDATA)
            data = np.fromfile(funit, dtype=dtype, count=ext.NAXIS1*ext.NAXIS2).reshape((ext.NAXIS2,ext.NAXIS1))
         else:
            funit.seek(ext.EXTDATA+order*ext.NAXIS1*dsize)
            data = np.fromfile(funit, dtype=dtype, count=ext.NAXIS1)
      return data


def bary(obj,bjd, exptime):
   from subprocess import Popen, PIPE, STDOUT
   #bjd = 2450000.5+ obj['MJD-OBS']
   #bjd = obj['HIERARCH ESO OBS START']
   #obj = obj['OBJECT']
   #pause(bjd)
   YY, MM, DD, hh, mm, ss = bjd.replace('T', ' ').replace(':', ' ').replace('-', ' ').split(' ')
   # check alias names
   with open('./src/bary/star_alias.txt','r') as f:
      for line in f:
         if obj in line:   obj=line.split()[0]
   # prepare input file
   s = ["Select observatory (LASILLA / PARANAL):",
        "LASILLA",
        "number of observations:",
        "1",
        "coord. (0) or object(1)?",
        "1",
        "Objectname in list:",
        obj,
        "Epoch:",
        "2000.0",
        "R.A. (hh mm ss)",
        "00 00 00",
        "Dec. (deg mm ss)",
        "00 00 00",
        "Date (dd mm yyyy)",
        DD+" "+MM+" "+YY, #"13 08 2000",
        "Time (hh mm ss)",
        hh+" "+mm+" "+ss, #"07 35 49",
        "Exptime (sec)",
        str(int(exptime))]
   #pause(exptime)
   p = Popen(['./src/bary/eph'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
   eph = p.communicate(input="\n".join(s)+"\n")[0]
   if  'Star not in starlist' in eph: 
      logger.error("Star not starList")
      stop(obj, 'Star not in starlist')
   berv = float(eph.split("\n")[-3].lstrip(" DELTAV=").rstrip("m/s"))
   return berv

def airtovac(wave_air):
   """
   taken from idl astrolib
   ;+
   ; NAME:
   ;       AIRTOVAC
   ; PURPOSE:
   ;       Convert air wavelengths to vacuum wavelengths 
   ; EXPLANATION:
   ;       Wavelengths are corrected for the index of refraction of air under 
   ;       standard conditions.  Wavelength values below 2000 A will not be 
   ;       altered.  Uses relation of Ciddor (1996).
   ;
   ; CALLING SEQUENCE:
   ;       AIRTOVAC, WAVE_AIR, [ WAVE_VAC]
   ;
   ; INPUT/OUTPUT:
   ;       WAVE_AIR - Wavelength in Angstroms, scalar or vector
   ;               If this is the only parameter supplied, it will be updated on
   ;               output to contain double precision vacuum wavelength(s). 
   ; OPTIONAL OUTPUT:
   ;        WAVE_VAC - Vacuum wavelength in Angstroms, same number of elements as
   ;                 WAVE_AIR, double precision
   ;
   ; EXAMPLE:
   ;       If the air wavelength is  W = 6056.125 (a Krypton line), then 
   ;       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
   ;
   ; METHOD:
   ;	Formula from Ciddor 1996, Applied Optics 62, 958
   ;
   ; NOTES:
   ;       Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in air* are
   ;       not altered.
   ; REVISION HISTORY
   ;       Written W. Landsman                November 1991
   ;       Use Ciddor (1996) formula for better accuracy in the infrared 
   ;           Added optional output vector, W Landsman Mar 2011
   ;       Iterate for better precision W.L./D. Schlegel  Mar 2011
   ;-
   """

   wave_vac = wave_air * 1.0
   g = wave_vac > 2000     #Only modify above 2000 A

   if np.sum(g):
      for iter in [0, 1]:
         if isinstance(g, np.ndarray):
            sigma2 = (1e4/wave_vac[g])**2.     #Convert to wavenumber squared
            # Compute conversion factor
            fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
                               1.67917e-3 / (57.362 - sigma2)
            wave_vac[g] = wave_air[g] * fact              #Convert Wavelength
         else: # scalar version
            sigma2 = (1e4/wave_vac)**2.     #Convert to wavenumber squared
            # Compute conversion factor
            fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
                               1.67917e-3 / (57.362 - sigma2)
            wave_vac = wave_air * fact              #Convert Wavelength

   return wave_vac


if __name__ == "__main__":
   pass


