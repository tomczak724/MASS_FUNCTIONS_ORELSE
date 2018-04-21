
import os
import sys
import time
import pylab
import subprocess
from collections import namedtuple
#from threedhst import eazyPy

#   >>> python py_mk_MC_zphot_catalogs.py /home/atomczak/DATA/ORELSE/catalogs nep200_v0.0.5 125
#   >>> python py_mk_MC_zphot_catalogs.py 
#       /home/atomczak/DATA/ORELSE/catalogs
#       nep200_v0.0.5
#       125

def readcat(name, dtype=float, delimiter=None, comments='#'):
    '''
    Reads an ascii catalog into a named tuple. It assumes that the
    row immediately prior to the data has the names of each column.
    Knows how to handle gzipped files.

    dtype  ------>  data type for columns
    delimiter --->  delimiter values that separates data entries
    comments  --->  character that indicates comment lines
    ''' 
    ###  gunzip if necessary
    gzipped = 0
    if name[-3:] == '.gz':
        gzipped = 1
        subprocess.call('gunzip %s' % name, shell=1)
        name = name[:-3]

    dat = open(name, 'r')

    ###  IDENTIFYING HEADER PARAMETERS
    lines = dat.readlines()
    for i, line in enumerate(lines):
        if lines[i+1][0] != comments:
            header_params = line[1:].split(delimiter)
            break
    del dat, lines

    custom_catalog = namedtuple('custom_catalog', header_params)
    catalog = custom_catalog(*pylab.loadtxt(name, dtype=dtype, delimiter=delimiter, unpack=1))

    ###  re-gzip if necessary
    if gzipped:
        subprocess.call('gzip %s' % name, shell=1)

    return catalog

def tomczak_sampler(pdf, N=1):
    '''
    Draws random samples from the provided PDF.

    Parameters
    ----------

    pdf : array
        2d array of the PDF. The first and second columns
        must be x and P(x) respectively. P(x) is not required
        to be normalized beforehand.
    N : int
        Number of desired samplings.

    Returns
    -------
    sample : array
        Random samplings from the provided PDF.
    '''
    
    x, px = pdf[:,0], pdf[:,1]
    px = px / pylab.trapz(px, x)

    u = []
    for i in range(len(x)):
        xprime = x[pylab.arange(i)]
        pxprime = px[pylab.arange(i)]
        uprime = pylab.trapz(pxprime, xprime)
        u.append(uprime)
    u = pylab.array(u)

        
    sample = []
    urand = pylab.rand(N)

    for ui in urand:

        du = abs(u - ui)
        ind = pylab.where(du == min(du))[0][0]
        sample.append(x[ind])

    return pylab.array(sample)

def readEazyBinary(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
    tempfilt, coeffs, temp_sed, pz = readEazyBinary(MAIN_OUTPUT_FILE='photz', \
                                                OUTPUT_DIRECTORY='./OUTPUT', \
                                                CACHE_FILE = 'Same')

    Read Eazy BINARY_OUTPUTS files into structure data.
    
    If the BINARY_OUTPUTS files are not in './OUTPUT', provide either a relative or absolute path
    in the OUTPUT_DIRECTORY keyword.
    
    By default assumes that CACHE_FILE is MAIN_OUTPUT_FILE+'.tempfilt'.
    Specify the full filename if otherwise. 
    """
    
    #root='COSMOS/OUTPUT/cat3.4_default_lines_zp33sspNoU'
    
    root = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE
    
    ###### .tempfilt
    if CACHE_FILE == 'Same':
        CACHE_FILE = root+'.tempfilt'
    
    if os.path.exists(CACHE_FILE) is False:
        print ('File, %s, not found.' %(CACHE_FILE))
        return -1,-1,-1,-1
    
    f = open(CACHE_FILE,'rb')
    
    s = pylab.fromfile(file=f,dtype=pylab.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    tempfilt = pylab.fromfile(file=f,dtype=pylab.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
    lc = pylab.fromfile(file=f,dtype=pylab.double,count=NFILT)
    zgrid = pylab.fromfile(file=f,dtype=pylab.double,count=NZ)
    fnu = pylab.fromfile(file=f,dtype=pylab.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    efnu = pylab.fromfile(file=f,dtype=pylab.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    
    f.close()
    
    tempfilt  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}
    
    ###### .coeff
    f = open(root+'.coeff','rb')
    
    s = pylab.fromfile(file=f,dtype=pylab.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    coeffs = pylab.fromfile(file=f,dtype=pylab.double,count=NTEMP*NOBJ).reshape((NOBJ,NTEMP)).transpose()
    izbest = pylab.fromfile(file=f,dtype=pylab.int32,count=NOBJ)
    tnorm = pylab.fromfile(file=f,dtype=pylab.double,count=NTEMP)
    
    f.close()
    
    coeffs = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
              'coeffs':coeffs,'izbest':izbest,'tnorm':tnorm}
              
    ###### .temp_sed
    f = open(root+'.temp_sed','rb')
    s = pylab.fromfile(file=f,dtype=pylab.int32, count=3)
    NTEMP=s[0]
    NTEMPL=s[1]
    NZ=s[2]
    templam = pylab.fromfile(file=f,dtype=pylab.double,count=NTEMPL)
    temp_seds = pylab.fromfile(file=f,dtype=pylab.double,count=NTEMPL*NTEMP).reshape((NTEMP,NTEMPL)).transpose()
    da = pylab.fromfile(file=f,dtype=pylab.double,count=NZ)
    db = pylab.fromfile(file=f,dtype=pylab.double,count=NZ)
    
    f.close()
    
    temp_sed = {'NTEMP':NTEMP,'NTEMPL':NTEMPL,'NZ':NZ,\
              'templam':templam,'temp_seds':temp_seds,'da':da,'db':db}
              
    ###### .pz
    if os.path.exists(root+'.pz'):
        f = open(root+'.pz','rb')
        s = pylab.fromfile(file=f,dtype=pylab.int32, count=2)
        NZ=s[0]
        NOBJ=s[1]
        chi2fit = pylab.fromfile(file=f,dtype=pylab.double,count=NZ*NOBJ).reshape((NOBJ,NZ)).transpose()

        ### This will break if APPLY_PRIOR No
        s = pylab.fromfile(file=f,dtype=pylab.int32, count=1)
        
        if len(s) > 0:
            NK = s[0]
            kbins = pylab.fromfile(file=f,dtype=pylab.double,count=NK)
            priorzk = pylab.fromfile(file=f, dtype=pylab.double, count=NZ*NK).reshape((NK,NZ)).transpose()
            kidx = pylab.fromfile(file=f,dtype=pylab.int32,count=NOBJ)
            pz = {'NZ':NZ,'NOBJ':NOBJ,'NK':NK, 'chi2fit':chi2fit, 'kbins':kbins, 'priorzk':priorzk,'kidx':kidx}
        else:
            priorzk = pylab.ones((1,NZ))
            kidx = pylab.zeros(NOBJ)
            pz = {'NZ':NZ,'NOBJ':NOBJ,'NK':0, 'chi2fit':chi2fit, 'kbins':[0], 'priorzk':priorzk,'kidx':kidx}
            #pz = None
        
        f.close()
        
    else:
        pz = None
        
    ###### Done.    
    return tempfilt, coeffs, temp_sed, pz

def getEazyPz(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
    zgrid, pz = getEazyPz(idx, \
                      MAIN_OUTPUT_FILE='photz', \
                      OUTPUT_DIRECTORY='./OUTPUT', \
                      CACHE_FILE='Same')
                      
    Get Eazy p(z) for object #idx.
    """
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                                    OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                                    CACHE_FILE = CACHE_FILE)
    
    if pz is None:
        return None, None
    
    ###### Get p(z|m) from prior grid
    kidx = pz['kidx'][idx]
    #print kidx, pz['priorzk'].shape
    if (kidx > 0) & (kidx < pz['priorzk'].shape[1]):
        prior = pz['priorzk'][:,kidx]
    else:
        prior = pylab.ones(pz['NZ'])
        
    ###### Convert Chi2 to p(z)
    pzi = pylab.exp(-0.5*(pz['chi2fit'][:,idx]-min(pz['chi2fit'][:,idx])))*prior
    if pylab.sum(pzi) > 0:
        pzi/=pylab.trapz(pzi, tempfilt['zgrid'])
    
    ###### Done
    return tempfilt['zgrid'],pzi





data_dir = sys.argv[1]
field = sys.argv[2]
n_mc = int(sys.argv[3])


cat = readcat('%s/%s/%s.cat' % (data_dir, field, field))
ngal = len(cat.id)


row = '# '
for i_mc in range(n_mc):
	row += ' z%04i' % (i_mc + 1)
print row

for i_gal in range(ngal):

	pz = getEazyPz(i_gal, MAIN_OUTPUT_FILE=field, OUTPUT_DIRECTORY='%s/%s/EAZY' % (data_dir, field))
	samples = tomczak_sampler(pylab.array(zip(pz[0], pz[1])), N=n_mc)
	for samp in samples:
		row += ' %.4f' % samp
	print row












