import argparse
from .consts import lines
import numpy as np 

def arg2slice(arg):
   """Convert string argument to a slice."""
   # We want four cases for indexing: None, int, list of ints, slices.
   # Use [] as default, so 'in' can be used.
   if isinstance(arg, str):
      arg = eval('np.s_['+arg+']')
   return [arg] if isinstance(arg, int) else arg


def build_parser(description, epilog, default, pmin, pmax, brvrefs, insts, oset):

    parser = argparse.ArgumentParser(description=description, epilog=epilog, add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('obj', help='Tag, output directory and file prefix (e.g. Object name).')
    argopt('dir_or_inputlist', help='Directory name with reduced data fits/tar or a file listing the spectra (only suffixes .txt or .lis accepted).', nargs='?')
    argopt('-targ', help='Target name requested in simbad for coordinates, proper motion, parallax and absolute RV.')
    argopt('-targrade', help='Target coordinates: [ra|hh:mm:ss.sss de|de:mm:ss.sss].', nargs=2, default=[None,None])
    argopt('-targpm', help='Target proper motion: pmra [mas/yr] pmde [mas/yr].', nargs=2, type=float, default=[0.0,0.0])
    argopt('-targplx', help='Target parallax', type=float, default='nan')
    argopt('-targrv', help='[km/s] Target RV guess (for index measures) [float, "drsspt", "drsmed", "targ", None, "auto"]. None => no measure; targ => from simbad, hdr; auto => first from headers, second from simbad))', default={'CARM_NIR':None, 'else':'auto'})
    argopt('-atmmask', help='Telluric line mask ('' for no masking)'+default, default='auto', dest='atmfile')
    argopt('-atmwgt', help='Downweighting factor for coadding in telluric regions'+default, type=float, default=None)
    argopt('-atmspec', help='Telluric spectrum  (in fits format, e.g. lib/stdatmos_vis30a090rh0780p000t.fits) to correct spectra by simple division.'+default, type=str, default=None)
    argopt('-brvref', help='Barycentric RV code reference', choices=brvrefs, type=str, default='WE')
    argopt('-msklist', help='Ascii table with vacuum wavelengths to mask.', default='') # [flux and width]
    argopt('-mskwd', help='[km/s] Broadening width for msklist.', type=float, default=4.)
    argopt('-mskspec', help='Ascii 0-1 spectrum.'+default, default='')
    argopt('-ccf',  help='mode ccf [with files]', nargs='?', const='th_mask_1kms.dat', type=str)
    argopt('-ccfmode', help='type for ccf template', nargs='?', default='box',
                        choices=['box', 'binless', 'gauss', 'trapeze'])
    argopt('-coadd', help='coadd method'+default, default='post3',
                    choices=['post3'])
    argopt('-coset', help='index for order in coadding (default: oset)', type=arg2slice)
    argopt('-co_excl', help='orders to exclude in coadding (default: o_excl)', type=arg2slice)
    argopt('-ckappa', help='kappa sigma (or lower and upper) clip value in coadding. Zero values for no clipping'+default, nargs='+', type=float, default=(4.,4.))
    argopt('-deg',  help='degree for background polynomial', type=int, default=3)
    argopt('-distmax', help='[arcsec] Max distance telescope position from target coordinates.', nargs='?', type=float, const=30.)
    argopt('-driftref', help='reference file for drift mode', type=str)
    argopt('-fib',  help='fibre', choices=['', 'A', 'B', 'AB'], default='')
    argopt('-inst', help='instrument '+default, default='HARPS', choices=insts)
    argopt('-nset', '-iset', help='slice for file subset (e.g. 1:10, ::5)', default=':', type=arg2slice)
    argopt('-kapsig', help='kappa sigma clip value'+default, type=float, default=3.0)
    argopt('-last', help='use last template (-tpl <obj>/template.fits)', action='store_true')
    argopt('-look', help='slice of orders to view the fit [:]', nargs='?', default=[], const=':', type=arg2slice)
    argopt('-looki', help='list of indices to watch', nargs='*', choices=sorted(lines.keys()), default=[]) #, const=['Halpha'])
    argopt('-lookt', help='slice of orders to view the coadd fit [:]', nargs='?', default=[], const=':', type=arg2slice)
    argopt('-lookp', help='slice of orders to view the preRV fit [:]', nargs='?', default=[], const=':', type=arg2slice)
    argopt('-lookssr', help='slice of orders to view the ssr function [:]', nargs='?', default=[], const=':', type=arg2slice)
    argopt('-lookmlRV', help='chi2map and master', nargs='?', default=[], const=':', type=arg2slice)
    argopt('-lookmlCRX', help='chi2map and CRX fit ', nargs='?', default=[], const=':', type=arg2slice)
    argopt('-nclip', help='max. number of clipping iterations'+default, type=int, default=2)
    argopt('-niter', help='number of RV iterations'+default, type=int, default=2)
    argopt('-oset', help='index for order subset (e.g. 1:10, ::5)'+default, default=oset, type=arg2slice)
    argopt('-o_excl', help='Orders to exclude (e.g. 1,10,3)', default={"CARM_NIR":"0,2,12,13,16,17,18,19,20,21,22,23,24,25,26,27,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,49,51,53,54,55", "else":[]}, type=arg2slice)
    #argopt('-outmod', help='output the modelling results for each spectrum into a fits file',  choices=['ratio', 'HARPN', 'CARM_VIS', 'CARM_NIR', 'FEROS', 'FTS'])
    argopt('-ofac', help='oversampling factor in coadding'+default, default=1., type=float)
    argopt('-ofacauto', help='automatic knot spacing with BIC.', action='store_true')
    argopt('-outchi', help='output of the chi2 map', nargs='?', const='_chi2map.fits')
    argopt('-outfmt', help='output format of the fits file (default: None; const: fmod err res wave)', nargs='*', choices=['wave', 'waverest', 'err', 'fmod', 'res', 'spec', 'bpmap', 'ratio'], default=None)
    argopt('-outsuf', help='output suffix', default='_mod.fits')
    argopt('-pmin', help='Minimum pixel'+default, default=pmin, type=int)
    argopt('-pmax', help='Maximum pixel'+default, default=pmax, type=int)
    argopt('-pspline', help='pspline as coadd filter [smooth value]', nargs='?', const=0.0000001, dest='pspllam', type=float)
    argopt('-pmu', help='analog to GP mean. Default no GP penalty. Without the mean in each order. Otherwise this value.', nargs='?', const=True, type=float)
    argopt('-pe_mu', help='analog to GP mean deviation', default=5., type=float)
    argopt('-reana', help='flag reanalyse only', action='store_true')
    argopt('-rvwarn', help='[km/s] warning threshold in debug'+default, default=2., type=float)
    argopt('-safemode', help='does not pause or stop, optional level 1  2 (reana)', nargs='?', type=int, const=1, default=False)
    argopt('-skippre', help='Skip pre-RVs.', action='store_true')
    argopt('-skymsk', help='Sky emission line mask ('' for no masking)'+default, default='auto', dest='skyfile')
    argopt('-snmin', help='minimum S/N (considered as not bad and used in template building)'+default, default=10, type=float)
    argopt('-snmax', help='maximum S/N (considered as not bad and used in template building)'+default, default=400, type=float)
    argopt('-tfmt', help='output format of the template. nmap is a an estimate for the number of good data points for each knot. ddspec is the second derivative for cubic spline reconstruction. (default: spec sig wave)', nargs='*', choices=['spec', 'sig', 'wave', 'nmap', 'ddspec'], default=['spec', 'sig', 'wave'])
    argopt('-tpl',  help="template filename or directory, if None or integer a template is created by coadding, where highest S/N spectrum or the filenr is used as start tpl for the pre-RVs", nargs='?')
    argopt('-tplrv', help='[km/s] template RV. By default taken from the template header and set to 0 km/s for phoe tpl.[float, "tpl", "drsspt", "drsmed", "targ", None, "auto"]', default='auto')
    argopt('-tset',  help="slice for file subset in template creation", default=':', type=arg2slice)
    argopt('-verb', help='verbose', action='store_true')
    v_lo, v_hi, v_step = -5.5, 5.6, 0.1
    argopt('-vrange', help='[km/s] velocity grid around targrv (v_lo, v_hi, v_step)'+default, nargs='*', default=(v_lo, v_hi, v_step), type=float)
    argopt('-vtfix', help='fix RV in template creation', action='store_true')
    argopt('-wfix', help='fix wavelength solution', action='store_true')
    argopt('-debug', help='debug flag', nargs='?', default=0, const=1)
    argopt('-bp',   help='break points', nargs='*', type=int)
    argopt('-pdb',  help='debug post_mortem', action='store_true')
    argopt('-cprofile', help='profiling', action='store_true')
    # use add_help=false and re-add with more arguments
    argopt('-?', '-h', '-help', '--help',  help='show this help message and exit', action='help')
    #parser.__dict__['_option_string_actions']['-h'].__dict__['option_strings'] += ['-?', '-help']

    return parser