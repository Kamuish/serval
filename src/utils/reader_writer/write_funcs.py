import astropy.io.fits as pyfits 

def write_template(filename, header, clobber, **kwargs):

    flux = kwargs['flux']
    wave = kwargs['wave']
    if isinstance(flux, np.ndarray):
        pyfits.append(filename, flux)
        pyfits.append(filename, wave)
    else:
        # pad arrays with zero to common size
        maxpix = max(arr.size for arr in flux if isinstance(arr, np.ndarray))
        flux_new = np.zeros((len(flux), maxpix))
        wave_new = np.zeros((len(flux), maxpix))
        for o,arr in enumerate(flux):
            if isinstance(arr, np.ndarray): flux_new[o,:len(arr)] = arr
        for o,arr in enumerate(wave):
            if isinstance(arr, np.ndarray): wave_new[o,:len(arr)] = arr
        pyfits.append(filename, flux_new)
        pyfits.append(filename, wave_new)

    pyfits.setval(filename, 'EXTNAME', value='SPEC', ext=1)
    pyfits.setval(filename, 'EXTNAME', value='WAVE', ext=2)
    #fitsio.write(filename, flux)

def write_res(filename, header, **kwargs):

    datas = kwargs['data']
    extnames = kwargs['extnames']
    for i,extname in enumerate(extnames):
        data = datas[extname]
        if isinstance(data, np.ndarray):
        pyfits.append(filename, data)
        else:
        1/0

        pyfits.setval(filename, 'EXTNAME', value=extname, ext=i+1)
    #fitsio.write(filename, flux)

def write_fits(filename, header, clobber, **kwargs):

    data = kwargs['data']
    warnings.resetwarnings() # supress nasty overwrite warning http://pythonhosted.org/pyfits/users_guide/users_misc.html
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    pyfits.writeto(filename, data, header, clobber=clobber, output_verify='fix')
    warnings.resetwarnings()
    warnings.filterwarnings('always', category=UserWarning, append=True)
