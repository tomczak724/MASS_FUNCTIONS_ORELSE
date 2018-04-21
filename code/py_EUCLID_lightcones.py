
import os
import glob
import h5py
import numpy
import shutil
import subprocess
import progressbar
from collections import namedtuple



###########################################################
###  The following code was written to read in HDF5 files
###  from the Lagos models in the directory:
###    /Users/atomczak/DATA/EUCLID/EUCLID_1_DEEP
###########################################################

if False:

	class lightcone(object):

		def __init__(self):

			self.BCDM = numpy.array([])
			self.BoT = numpy.array([])
			self.DHaloID = numpy.array([])
			self.EW_tot_Halpha = numpy.array([])
			self.EW_tot_Hbeta = numpy.array([])
			self.EW_tot_Hgamma = numpy.array([])
			self.EW_tot_Lyalpha = numpy.array([])
			self.EW_tot_NII6583 = numpy.array([])
			self.EW_tot_OII3727 = numpy.array([])
			self.EW_tot_OIII4959 = numpy.array([])
			self.EW_tot_OIII5007 = numpy.array([])
			self.GalaxyID = numpy.array([])
			self.L_tot_Halpha = numpy.array([])
			self.L_tot_Hbeta = numpy.array([])
			self.L_tot_Hgamma = numpy.array([])
			self.L_tot_Lyalpha = numpy.array([])
			self.L_tot_NII6583 = numpy.array([])
			self.L_tot_OII3727 = numpy.array([])
			self.L_tot_OIII4959 = numpy.array([])
			self.L_tot_OIII5007 = numpy.array([])
			self.dec = numpy.array([])
			self.idrep = numpy.array([])
			self.is_central = numpy.array([])
			self.magBo_tot_ext = numpy.array([])
			self.magBr_tot_ext = numpy.array([])
			self.magDYo_tot_ext = numpy.array([])
			self.magDYr_tot_ext = numpy.array([])
			self.magDgo_tot_ext = numpy.array([])
			self.magDgr_tot_ext = numpy.array([])
			self.magDio_tot_ext = numpy.array([])
			self.magDir_tot_ext = numpy.array([])
			self.magDro_tot_ext = numpy.array([])
			self.magDrr_tot_ext = numpy.array([])
			self.magDzo_tot_ext = numpy.array([])
			self.magDzr_tot_ext = numpy.array([])
			self.magEHo_tot_ext = numpy.array([])
			self.magEHr_tot_ext = numpy.array([])
			self.magEJo_tot_ext = numpy.array([])
			self.magEJr_tot_ext = numpy.array([])
			self.magEKo_tot_ext = numpy.array([])
			self.magEKr_tot_ext = numpy.array([])
			self.magEYo_tot_ext = numpy.array([])
			self.magEYr_tot_ext = numpy.array([])
			self.magIo_tot_ext = numpy.array([])
			self.magIr_tot_ext = numpy.array([])
			self.magRo_tot_ext = numpy.array([])
			self.magRr_tot_ext = numpy.array([])
			self.magUo_tot_ext = numpy.array([])
			self.magUr_tot_ext = numpy.array([])
			self.magVo_tot_ext = numpy.array([])
			self.magVr_tot_ext = numpy.array([])
			self.maggPo_tot_ext = numpy.array([])
			self.maggPr_tot_ext = numpy.array([])
			self.maggSo_tot_ext = numpy.array([])
			self.maggSr_tot_ext = numpy.array([])
			self.magiPo_tot_ext = numpy.array([])
			self.magiPr_tot_ext = numpy.array([])
			self.magiSo_tot_ext = numpy.array([])
			self.magiSr_tot_ext = numpy.array([])
			self.magrPo_tot_ext = numpy.array([])
			self.magrPr_tot_ext = numpy.array([])
			self.magrSo_tot_ext = numpy.array([])
			self.magrSr_tot_ext = numpy.array([])
			self.maguSo_tot_ext = numpy.array([])
			self.maguSr_tot_ext = numpy.array([])
			self.magwPo_tot_ext = numpy.array([])
			self.magwPr_tot_ext = numpy.array([])
			self.magyPo_tot_ext = numpy.array([])
			self.magyPr_tot_ext = numpy.array([])
			self.magzPo_tot_ext = numpy.array([])
			self.magzPr_tot_ext = numpy.array([])
			self.magzSo_tot_ext = numpy.array([])
			self.magzSr_tot_ext = numpy.array([])
			self.mhhalo = numpy.array([])
			self.mstardot = numpy.array([])
			self.mstardot_average = numpy.array([])
			self.mstardot_burst = numpy.array([])
			self.mstars_bulge = numpy.array([])
			self.mstars_disk = numpy.array([])
			self.mstars_tot = numpy.array([])
			self.phi = numpy.array([])
			self.r = numpy.array([])
			self.ra = numpy.array([])
			self.rbulge = numpy.array([])
			self.rdisk = numpy.array([])
			self.theta = numpy.array([])
			self.vdisk = numpy.array([])
			self.vxgal = numpy.array([])
			self.vygal = numpy.array([])
			self.vzgal = numpy.array([])
			self.xgal = numpy.array([])
			self.ygal = numpy.array([])
			self.z_cos = numpy.array([])
			self.z_obs = numpy.array([])
			self.zgal = numpy.array([])

	Lagos = lightcone()

	###  Initializing progress bar
	widgets = ['  Running: ', progressbar.Percentage(), 
	           ' ', progressbar.Bar(marker=progressbar.RotatingMarker()), 
	           ' ', progressbar.ETA(), 
	           ' ', progressbar.FileTransferSpeed()]
	pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(files)+1).start()

	for i_file in range(len(files)):

		pbar.update(i_file+1)

		fopen = h5py.File(files[i_file], "r")
		try:
			dset = fopen.get('Data')
			Lagos.BCDM = numpy.append(Lagos.BCDM, dset["BCDM"].value)
			Lagos.BoT = numpy.append(Lagos.BoT, dset["BoT"].value)
			Lagos.DHaloID = numpy.append(Lagos.DHaloID, dset["DHaloID"].value)
			Lagos.EW_tot_Halpha = numpy.append(Lagos.EW_tot_Halpha, dset["EW_tot_Halpha"].value)
			Lagos.EW_tot_Hbeta = numpy.append(Lagos.EW_tot_Hbeta, dset["EW_tot_Hbeta"].value)
			Lagos.EW_tot_Hgamma = numpy.append(Lagos.EW_tot_Hgamma, dset["EW_tot_Hgamma"].value)
			Lagos.EW_tot_Lyalpha = numpy.append(Lagos.EW_tot_Lyalpha, dset["EW_tot_Lyalpha"].value)
			Lagos.EW_tot_NII6583 = numpy.append(Lagos.EW_tot_NII6583, dset["EW_tot_NII6583"].value)
			Lagos.EW_tot_OII3727 = numpy.append(Lagos.EW_tot_OII3727, dset["EW_tot_OII3727"].value)
			Lagos.EW_tot_OIII4959 = numpy.append(Lagos.EW_tot_OIII4959, dset["EW_tot_OIII4959"].value)
			Lagos.EW_tot_OIII5007 = numpy.append(Lagos.EW_tot_OIII5007, dset["EW_tot_OIII5007"].value)
			Lagos.GalaxyID = numpy.append(Lagos.GalaxyID, dset["GalaxyID"].value)
			Lagos.L_tot_Halpha = numpy.append(Lagos.L_tot_Halpha, dset["L_tot_Halpha"].value)
			Lagos.L_tot_Hbeta = numpy.append(Lagos.L_tot_Hbeta, dset["L_tot_Hbeta"].value)
			Lagos.L_tot_Hgamma = numpy.append(Lagos.L_tot_Hgamma, dset["L_tot_Hgamma"].value)
			Lagos.L_tot_Lyalpha = numpy.append(Lagos.L_tot_Lyalpha, dset["L_tot_Lyalpha"].value)
			Lagos.L_tot_NII6583 = numpy.append(Lagos.L_tot_NII6583, dset["L_tot_NII6583"].value)
			Lagos.L_tot_OII3727 = numpy.append(Lagos.L_tot_OII3727, dset["L_tot_OII3727"].value)
			Lagos.L_tot_OIII4959 = numpy.append(Lagos.L_tot_OIII4959, dset["L_tot_OIII4959"].value)
			Lagos.L_tot_OIII5007 = numpy.append(Lagos.L_tot_OIII5007, dset["L_tot_OIII5007"].value)
			Lagos.dec = numpy.append(Lagos.dec, dset["dec"].value)
			Lagos.idrep = numpy.append(Lagos.idrep, dset["idrep"].value)
			Lagos.is_central = numpy.append(Lagos.is_central, dset["is_central"].value)
			Lagos.magBo_tot_ext = numpy.append(Lagos.magBo_tot_ext, dset["magBo_tot_ext"].value)
			Lagos.magBr_tot_ext = numpy.append(Lagos.magBr_tot_ext, dset["magBr_tot_ext"].value)
			Lagos.magDYo_tot_ext = numpy.append(Lagos.magDYo_tot_ext, dset["magDYo_tot_ext"].value)
			Lagos.magDYr_tot_ext = numpy.append(Lagos.magDYr_tot_ext, dset["magDYr_tot_ext"].value)
			Lagos.magDgo_tot_ext = numpy.append(Lagos.magDgo_tot_ext, dset["magDgo_tot_ext"].value)
			Lagos.magDgr_tot_ext = numpy.append(Lagos.magDgr_tot_ext, dset["magDgr_tot_ext"].value)
			Lagos.magDio_tot_ext = numpy.append(Lagos.magDio_tot_ext, dset["magDio_tot_ext"].value)
			Lagos.magDir_tot_ext = numpy.append(Lagos.magDir_tot_ext, dset["magDir_tot_ext"].value)
			Lagos.magDro_tot_ext = numpy.append(Lagos.magDro_tot_ext, dset["magDro_tot_ext"].value)
			Lagos.magDrr_tot_ext = numpy.append(Lagos.magDrr_tot_ext, dset["magDrr_tot_ext"].value)
			Lagos.magDzo_tot_ext = numpy.append(Lagos.magDzo_tot_ext, dset["magDzo_tot_ext"].value)
			Lagos.magDzr_tot_ext = numpy.append(Lagos.magDzr_tot_ext, dset["magDzr_tot_ext"].value)
			Lagos.magEHo_tot_ext = numpy.append(Lagos.magEHo_tot_ext, dset["magEHo_tot_ext"].value)
			Lagos.magEHr_tot_ext = numpy.append(Lagos.magEHr_tot_ext, dset["magEHr_tot_ext"].value)
			Lagos.magEJo_tot_ext = numpy.append(Lagos.magEJo_tot_ext, dset["magEJo_tot_ext"].value)
			Lagos.magEJr_tot_ext = numpy.append(Lagos.magEJr_tot_ext, dset["magEJr_tot_ext"].value)
			Lagos.magEKo_tot_ext = numpy.append(Lagos.magEKo_tot_ext, dset["magEKo_tot_ext"].value)
			Lagos.magEKr_tot_ext = numpy.append(Lagos.magEKr_tot_ext, dset["magEKr_tot_ext"].value)
			Lagos.magEYo_tot_ext = numpy.append(Lagos.magEYo_tot_ext, dset["magEYo_tot_ext"].value)
			Lagos.magEYr_tot_ext = numpy.append(Lagos.magEYr_tot_ext, dset["magEYr_tot_ext"].value)
			Lagos.magIo_tot_ext = numpy.append(Lagos.magIo_tot_ext, dset["magIo_tot_ext"].value)
			Lagos.magIr_tot_ext = numpy.append(Lagos.magIr_tot_ext, dset["magIr_tot_ext"].value)
			Lagos.magRo_tot_ext = numpy.append(Lagos.magRo_tot_ext, dset["magRo_tot_ext"].value)
			Lagos.magRr_tot_ext = numpy.append(Lagos.magRr_tot_ext, dset["magRr_tot_ext"].value)
			Lagos.magUo_tot_ext = numpy.append(Lagos.magUo_tot_ext, dset["magUo_tot_ext"].value)
			Lagos.magUr_tot_ext = numpy.append(Lagos.magUr_tot_ext, dset["magUr_tot_ext"].value)
			Lagos.magVo_tot_ext = numpy.append(Lagos.magVo_tot_ext, dset["magVo_tot_ext"].value)
			Lagos.magVr_tot_ext = numpy.append(Lagos.magVr_tot_ext, dset["magVr_tot_ext"].value)
			Lagos.maggPo_tot_ext = numpy.append(Lagos.maggPo_tot_ext, dset["maggPo_tot_ext"].value)
			Lagos.maggPr_tot_ext = numpy.append(Lagos.maggPr_tot_ext, dset["maggPr_tot_ext"].value)
			Lagos.maggSo_tot_ext = numpy.append(Lagos.maggSo_tot_ext, dset["maggSo_tot_ext"].value)
			Lagos.maggSr_tot_ext = numpy.append(Lagos.maggSr_tot_ext, dset["maggSr_tot_ext"].value)
			Lagos.magiPo_tot_ext = numpy.append(Lagos.magiPo_tot_ext, dset["magiPo_tot_ext"].value)
			Lagos.magiPr_tot_ext = numpy.append(Lagos.magiPr_tot_ext, dset["magiPr_tot_ext"].value)
			Lagos.magiSo_tot_ext = numpy.append(Lagos.magiSo_tot_ext, dset["magiSo_tot_ext"].value)
			Lagos.magiSr_tot_ext = numpy.append(Lagos.magiSr_tot_ext, dset["magiSr_tot_ext"].value)
			Lagos.magrPo_tot_ext = numpy.append(Lagos.magrPo_tot_ext, dset["magrPo_tot_ext"].value)
			Lagos.magrPr_tot_ext = numpy.append(Lagos.magrPr_tot_ext, dset["magrPr_tot_ext"].value)
			Lagos.magrSo_tot_ext = numpy.append(Lagos.magrSo_tot_ext, dset["magrSo_tot_ext"].value)
			Lagos.magrSr_tot_ext = numpy.append(Lagos.magrSr_tot_ext, dset["magrSr_tot_ext"].value)
			Lagos.maguSo_tot_ext = numpy.append(Lagos.maguSo_tot_ext, dset["maguSo_tot_ext"].value)
			Lagos.maguSr_tot_ext = numpy.append(Lagos.maguSr_tot_ext, dset["maguSr_tot_ext"].value)
			Lagos.magwPo_tot_ext = numpy.append(Lagos.magwPo_tot_ext, dset["magwPo_tot_ext"].value)
			Lagos.magwPr_tot_ext = numpy.append(Lagos.magwPr_tot_ext, dset["magwPr_tot_ext"].value)
			Lagos.magyPo_tot_ext = numpy.append(Lagos.magyPo_tot_ext, dset["magyPo_tot_ext"].value)
			Lagos.magyPr_tot_ext = numpy.append(Lagos.magyPr_tot_ext, dset["magyPr_tot_ext"].value)
			Lagos.magzPo_tot_ext = numpy.append(Lagos.magzPo_tot_ext, dset["magzPo_tot_ext"].value)
			Lagos.magzPr_tot_ext = numpy.append(Lagos.magzPr_tot_ext, dset["magzPr_tot_ext"].value)
			Lagos.magzSo_tot_ext = numpy.append(Lagos.magzSo_tot_ext, dset["magzSo_tot_ext"].value)
			Lagos.magzSr_tot_ext = numpy.append(Lagos.magzSr_tot_ext, dset["magzSr_tot_ext"].value)
			Lagos.mhhalo = numpy.append(Lagos.mhhalo, dset["mhhalo"].value)
			Lagos.mstardot = numpy.append(Lagos.mstardot, dset["mstardot"].value)
			Lagos.mstardot_average = numpy.append(Lagos.mstardot_average, dset["mstardot_average"].value)
			Lagos.mstardot_burst = numpy.append(Lagos.mstardot_burst, dset["mstardot_burst"].value)
			Lagos.mstars_bulge = numpy.append(Lagos.mstars_bulge, dset["mstars_bulge"].value)
			Lagos.mstars_disk = numpy.append(Lagos.mstars_disk, dset["mstars_disk"].value)
			Lagos.mstars_tot = numpy.append(Lagos.mstars_tot, dset["mstars_tot"].value)
			Lagos.phi = numpy.append(Lagos.phi, dset["phi"].value)
			Lagos.r = numpy.append(Lagos.r, dset["r"].value)
			Lagos.ra = numpy.append(Lagos.ra, dset["ra"].value)
			Lagos.rbulge = numpy.append(Lagos.rbulge, dset["rbulge"].value)
			Lagos.rdisk = numpy.append(Lagos.rdisk, dset["rdisk"].value)
			Lagos.theta = numpy.append(Lagos.theta, dset["theta"].value)
			Lagos.vdisk = numpy.append(Lagos.vdisk, dset["vdisk"].value)
			Lagos.vxgal = numpy.append(Lagos.vxgal, dset["vxgal"].value)
			Lagos.vygal = numpy.append(Lagos.vygal, dset["vygal"].value)
			Lagos.vzgal = numpy.append(Lagos.vzgal, dset["vzgal"].value)
			Lagos.xgal = numpy.append(Lagos.xgal, dset["xgal"].value)
			Lagos.ygal = numpy.append(Lagos.ygal, dset["ygal"].value)
			Lagos.z_cos = numpy.append(Lagos.z_cos, dset["z_cos"].value)
			Lagos.z_obs = numpy.append(Lagos.z_obs, dset["z_obs"].value)
			Lagos.zgal = numpy.append(Lagos.zgal, dset["zgal"].value)
		except:
			pass








###########################################################
###  The following code was written to read in HDF5 files
###  from the Bower06 models in the directory:
###    /Users/atomczak/DATA/EUCLID/EUCLID_100sqdeg
###########################################################


def read_Bower06_catalog_data(name, dtype=float, delimiter=None, comments='#'):
    '''
    Reads a Bower06.* catalog into a named tuple. It assumes that the
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
            header_params = line[1:].split(delimiter)[1:]
            break
    dat.close()
    del dat, lines

    custom_catalog = namedtuple('custom_catalog', header_params)
    catalog = custom_catalog(*numpy.loadtxt(name, dtype=dtype, delimiter=delimiter, unpack=1))

    ###  re-gzip if necessary
    if gzipped:
        subprocess.call('gzip %s' % name, shell=1)

    return catalog


class Bower06_catalog(object):

	def __init__(self, filename):

		self.filename = filename

		###  reading the 12 header parameters
		dat = open(filename, 'r')
		lines = dat.readlines()

		self.omega0       = float(lines[0].split()[2])
		self.lambda0      = float(lines[1].split()[2])
		self.omegab       = float(lines[2].split()[2])
		self.h0           = float(lines[3].split()[2])
		self.zrange       = [float(lines[4].split()[2]), float(lines[4].split()[3])]
		self.view_radius  = float(lines[5].split()[2])
		self.field_centre = [float(lines[6].split()[2]), float(lines[6].split()[3])]
		self.n_sel_id     = int(lines[7].split()[2])
		self.sel_id       = str(lines[8].split()[2])
		self.low_lim      = float(lines[9].split()[2])
		self.upp_lim      = float(lines[10].split()[2])
		self.nprops       = int(lines[11].split()[2])

		dat.close()
		del lines

		self.data = read_Bower06_catalog_data(filename)




files =  glob.glob('/Users/atomczak/DATA/EUCLID/EUCLID_100sqdeg/Bower06*')

catalogs = [Bower06_catalog(f) for f in files]
















