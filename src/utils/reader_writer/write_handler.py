from .write_funcs import * 
import warnings
import astropy.io.fits as pyfits 


def write_handler(write_type, filename, **kwargs):
    """
        Allow to control the writting process;
        Either store the template, the res of the fits file;

    """
    mappings = {
        'template': write_template,
        'res': write_res,
        'fits' : write_fits
    }
    try:
        hddref = kwargs['hddref']
    except:
        hddref = None
    
    try:
        clobber = kwargs['clobber']
    except:
        clobber = False if write_type != 'fits' else True 

    if not kwargs['header'] and kwargs['hdrref']: 
        header = pyfits.getheader(hdrref) 
    else:
        header = kwargs['header']

    if write_type != 'fits':
 
        hdu = pyfits.PrimaryHDU(header=header)
        warnings.resetwarnings() # supress nasty overwrite warning http://pythonhosted.org/pyfits/users_guide/users_misc.html
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        hdu.writeto(filename, overwrite=clobber, output_verify='fix')
        warnings.resetwarnings()
        warnings.filterwarnings('always', category=UserWarning, append=True)
    mappings[write_type](filename = filename,
                        header_use = header,
                        clobber_value = clobber, 
                        **kwargs)