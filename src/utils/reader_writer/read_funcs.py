
from collections import namedtuple
import tarfile 
from src.utils import FitsClass

import glob
from src.read_spec import imhead


def read_template(filename):
   hdu = pyfits.open(filename)
   return hdu[2].data, hdu[1].data, hdu[0].header  # wave, flux

def read_harps_ccf(path_to_file):
    ccf = namedtuple('ccf', 'rvc err_rvc bis fwhm contrast mask')
    tar = None
    if ".tar" in path_to_file:
        tar = tarfile.open(path_to_file)
        extr = None
        for member in tar.getmembers():
            if 'A.fits' in member.name:
                if '_ccf_' in member.name and not extr: 
                    extr = member
                if '_bis_' in member.name: 
                    extr = member   # prefer bis 

        if not extr: 
            return ccf(0,0,0,0,0,0)
        s = FitsClass(path_to_file, extr.name, extr.offset_data, extr.size)
    else:
        s = glob.glob(s.replace("_e2ds","_bis_*")) + glob.glob(s.replace("_e2ds","_ccf_*"))
        if s: s = s[0]   # prefer bis
        else: 
            return ccf(*[np.nan]*6)
    HIERARCH = 'HIERARCH '

    hdr = imhead(s, HIERARCH+'ESO DRS CCF RVC', HIERARCH+'ESO DRS CCF CONTRAST', HIERARCH+'ESO DRS CCF FWHM', HIERARCH+'ESO DRS CCF MASK', HIERARCH+'ESO DRS DVRMS',HIERARCH+'ESO DRS BIS SPAN')
  

    if tar:
        tar.close()

    rvc = hdr[HIERARCH+'ESO DRS CCF RVC']   # [km/s]
    contrast = hdr.get(HIERARCH+'ESO DRS CCF CONTRAST', np.nan)
    fwhm = hdr[HIERARCH+'ESO DRS CCF FWHM']
    mask = hdr[HIERARCH+'ESO DRS CCF MASK']
    e_rvc = hdr.get(HIERARCH+'ESO DRS DVRMS', np.nan) / 1000.   # [km/s]
    bis = hdr.get(HIERARCH+'ESO DRS BIS SPAN', np.nan)
    return ccf(rvc, e_rvc, bis, fwhm, contrast, mask)
