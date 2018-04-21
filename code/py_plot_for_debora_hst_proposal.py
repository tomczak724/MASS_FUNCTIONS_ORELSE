
import os
import mypy
import time
import pylab
import pickle
import my_classes
import subprocess
from mypy import massfunc
from scipy import optimize
import matplotlib.patheffects as PathEffects
from astropy import units, constants, cosmology, modeling
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties


def fqu_arctan(lmassax, m0, beta):
	return pylab.arctan(beta*(lmassax - m0)) * (1. / pylab.pi) + 0.5

def fqu_exp(lmassax, a, b):
	return 1 - pylab.exp(-(10**lmassax / a)**b)

def fqu_peng2010_z1(lmassax, loverdens):
	'''
	Parameterization of quiescent fraction from Peng+2010
	for highest zbin 0.7 < z < 1
	'''
	p1, p2, p3, p4 = 10**1.9, 0.64, 10**10.89, 0.63
	return 1 - pylab.exp(-(10*loverdens / p1)**p2 - (10**lmassax / p3)**p4)


def line(x, m, b): return m * x + b

def smf_TOT_parameterized_leja2015(lmassax, zi):
	'''
	Returns the number densities at the input M* values
	based on the parameterized SMFs from Leja+2015. See
	their Section 2.
	'''
	a1, a2 = -0.39, -1.53
	lphistar1 = -2.64 + 0.07*zi - 0.28*zi**2
	lphistar2 = -3.11 - 0.18*zi - 0.03*zi**2
	lmstar = 10.72 - 0.13*zi + 0.11*zi**2
	return dschechter(lmassax, lmstar, a1, a2, 10**lphistar1, 10**lphistar2)

def schechter_mf(xaxis, alpha, mstar, phistar):
    """
    #  DESCRIPTION:
    #    Returns the values for a Schechter mass function
    #    from a given mass-axis and Schechter parameters.
    #
    #  INPUTS:
    #    xaxis = input mass value(s)
    #    alpha = Schechter parameters
    #    mstar = Schechter parameters
    #    phistar = Schechter parameters    
    """
    return pylab.log(10) * phistar * 10**((xaxis-mstar)*(1+alpha)) * pylab.exp(-10**(xaxis-mstar))
def dschechter(lmassax, lmstar, a1, a2, phistar1, phistar2):
	factor1 = pylab.log(10) * pylab.exp(-10**(lmassax - lmstar)) * 10**(lmassax - lmstar)
	factor2 = phistar1 * 10**(a1*(lmassax - lmstar)) + phistar2 * 10**(a2*(lmassax - lmstar))
	return factor1 * factor2

def dschechter2(lmassax, lmstar, a1, a2, phistar1, phistar2):
	mratio = 10**lmassax / 10**lmstar
	factor1 = pylab.exp(-(mratio))
	factor2 = phistar1 * mratio**a1 + phistar2 * mratio**a2
	return factor1 * factor2

def schechter_mf_log(xaxis, alpha, mstar, phistar):
    """
    #  DESCRIPTION:
    #    Returns the values for a Schechter mass function
    #    from a given mass-axis and Schechter parameters.
    #
    #  INPUTS:
    #    xaxis = input mass value(s)
    #    alpha = Schechter parameters
    #    mstar = Schechter parameters
    #    phistar = Schechter parameters    
    """
    return pylab.log10(pylab.log(10) * phistar * 10**((xaxis-mstar)*(1+alpha)) * pylab.exp(-10**(xaxis-mstar)))


def add_inset(subplot, rect=[0.4, 0.4, 0.2, 0.2]):
        '''
        This script creates an Axes instance within the
        coordinate frame of the provided subplot.
        '''

        ###  coordinates of the subplot in the figure window's coordinate frame
        box = subplot.get_position()
        xlo, xhi = box.x0, box.x1
        ylo, yhi = box.y0, box.y1
        dx, dy = (xhi - xlo), (yhi - ylo)

        ###  width/height of requested axes in the figure window's coordinate frame
        sub_dx = dx * rect[2]
        sub_dy = dy * rect[3]

        ###  position of requested axes in the figure window's coordinate frame
        sub_xlo = xlo + dx * rect[0]
        sub_ylo = ylo + dy * rect[1]

        inset = subplot.figure.add_axes([sub_xlo, sub_ylo, sub_dx, sub_dy])
        return inset


def illustris_merger_rates(lmass, mu, z=0.8):
	###  lmass = mass of the descendent galaxy
	###  mu    = mass ratio of merger (always <1)
	A0 = 10**-2.2287
	M0 = 2.e11
	eta = 2.4644
	alpha0 = 0.2241
	alpha1 = -1.1759
	beta0 = -1.2595
	beta1 = 0.0611
	gamma = -0.0477
	delta0 = 0.7668
	delta1 = -0.4695

	Az = A0 * (1 + z)**eta
	alphaz = alpha0 * (1 + z)**alpha1
	betaz = beta0 * (1 + z)**beta1
	deltaz = delta0 * (1 + z)**delta1

	factor1 = Az * (10**lmass / 10**10)**alphaz
	factor2 = (1 + (10**lmass / M0)**deltaz)
	factor3 = mu**(betaz + gamma * pylab.log10(10**lmass / 10**10))

	return factor1 * factor2 * factor3




fit_algorithm = modeling.fitting.LevMarLSQFitter()
class single_schechter_model(modeling.Fittable1DModel):
	lmstar = modeling.Parameter(default=10.8)
	alpha = modeling.Parameter(default=-1.)
	phistar = modeling.Parameter(default=10**-2.)

	@staticmethod
	def evaluate(x, lmstar, alpha, phistar):
		return pylab.log(10) * phistar * 10**((x-lmstar)*(1+alpha)) * pylab.exp(-10**(x-lmstar))
class double_schechter_model(modeling.Fittable1DModel):
	lmstar = modeling.Parameter(default=10.8)
	alpha1 = modeling.Parameter(default=-1.)
	alpha2 = modeling.Parameter(default=-2.54)
	phistar1 = modeling.Parameter(default=10**-2.)
	phistar2 = modeling.Parameter(default=10**-4.)

	@staticmethod
	def evaluate(x, lmstar, alpha1, alpha2, phistar1, phistar2):
		factor1 = pylab.log(10) * pylab.exp(-10**(x - lmstar)) * 10**(x - lmstar)
		factor2 = phistar1 * 10**(alpha1*(x - lmstar)) + phistar2 * 10**(alpha2*(x - lmstar))
		return factor1 * factor2


class single_schechter_model_log(modeling.Fittable1DModel):
	lmstar = modeling.Parameter(default=10.8)
	alpha = modeling.Parameter(default=-1.)
	phistar = modeling.Parameter(default=10**-2.)

	@staticmethod
	def evaluate(x, lmstar, alpha, phistar):
		return pylab.log10(pylab.log(10) * phistar * 10**((x-lmstar)*(1+alpha)) * pylab.exp(-10**(x-lmstar)))





def illustris_merger_rates(lmass, mu, z=0.8):
	###  lmass = mass of the descendent galaxy
	###  mu    = mass ratio of merger (always <1)
	A0 = 10**-2.2287
	M0 = 2.e11
	eta = 2.4644
	alpha0 = 0.2241
	alpha1 = -1.1759
	beta0 = -1.2595
	beta1 = 0.0611
	gamma = -0.0477
	delta0 = 0.7668
	delta1 = -0.4695

	Az = A0 * (1 + z)**eta
	alphaz = alpha0 * (1 + z)**alpha1
	betaz = beta0 * (1 + z)**beta1
	deltaz = delta0 * (1 + z)**delta1

	factor1 = Az * (10**lmass / 10**10)**alphaz
	factor2 = (1 + (10**lmass / M0)**deltaz)
	factor3 = mu**(betaz + gamma * pylab.log10(10**lmass / 10**10))

	return factor1 * factor2 * factor3





data_dir = '../data'

class massFunction_field:
	def __init__(self, name, area, lmassbars, zlos, zhis):
		self.name = name
		self.area = area    # unobscured survey area in arcmin**2
		self.lmassbars = lmassbars
		self.zlos = zlos
		self.zhis = zhis
		self.zbars = (zhis + zlos) / 2.
		self.dz = (zhis - zlos)
		self.number_counts = pylab.array([pylab.zeros(len(lmassbars)) for zi in range(len(zlos))])
		self.nlo_poisson = pylab.array([pylab.zeros(len(lmassbars)) for zi in range(len(zlos))])
		self.nhi_poisson = pylab.array([pylab.zeros(len(lmassbars)) for zi in range(len(zlos))])

class massFunction_voronoi_slices:
	def __init__(self, name, zclust, sigmaz, area_zspec_arcmin2, zlos, zhis):
		self.name = name
		self.zclust = zclust
		self.sigmaz = sigmaz
		self.area_zspec_arcmin2 = area_zspec_arcmin2

		self.zlos = zlos                  # lower zlimits for voronoi slices
		self.zhis = zhis                  # upper zlimits for voronoi slices
		self.zbars = (zhis + zlos) / 2.
		self.dz = (zhis - zlos) / 2.

		self.dm = 0.25
		self.lmassbins = pylab.arange(8.0-self.dm/2., 11.75+self.dm, self.dm)
		self.lmassbars = (self.lmassbins[1:] + self.lmassbins[:-1]) / 2.

		self.areas_full_arcmin2 = pylab.zeros(len(self.zbars))               # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_field_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_v_lt_0_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_n05v00_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_00v05_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_05v10_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_10v15_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_15v20_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice

		self.volumes_full_Mpc3 = pylab.zeros(len(self.zbars))               # volumes corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_field_Mpc3 = pylab.zeros(len(self.zbars))              # volumes corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_v_lt_0_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_n05v00_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_00v05_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_05v10_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_10v15_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_15v20_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice

		self.ngals_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.d_overdens = 0.5
		self.overdens_bins = pylab.array([-0.5, 0., 0.5, 1., 1.5, 2.])
		self.overdens_bars = (self.overdens_bins[1:] + self.overdens_bins[:-1]) / 2.

		self.ngals_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.areas_all_arcmin2 = [self.areas_n05v00_arcmin2, self.areas_00v05_arcmin2, self.areas_05v10_arcmin2, self.areas_10v15_arcmin2, self.areas_15v20_arcmin2]
		self.volumes_all_Mpc3 = [self.volumes_n05v00_Mpc3, self.volumes_00v05_Mpc3, self.volumes_05v10_Mpc3, self.volumes_10v15_Mpc3, self.volumes_15v20_Mpc3]
		self.ngals_all = [self.ngals_n05v00, self.ngals_00v05, self.ngals_05v10, self.ngals_10v15, self.ngals_15v20]
		self.elo_ngals_all = [self.elo_ngals_n05v00, self.elo_ngals_00v05, self.elo_ngals_05v10, self.elo_ngals_10v15, self.elo_ngals_15v20]
		self.ehi_ngals_all = [self.ehi_ngals_n05v00, self.ehi_ngals_00v05, self.ehi_ngals_05v10, self.ehi_ngals_10v15, self.ehi_ngals_15v20]
		self.phis_all = [self.phis_n05v00, self.phis_00v05, self.phis_05v10, self.phis_10v15, self.phis_15v20]



		##################
		###  SF galaxies
		##################

		self.ngals_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_all_sf = [self.ngals_n05v00_sf, self.ngals_00v05_sf, self.ngals_05v10_sf, self.ngals_10v15_sf, self.ngals_15v20_sf]
		self.elo_ngals_all_sf = [self.elo_ngals_n05v00_sf, self.elo_ngals_00v05_sf, self.elo_ngals_05v10_sf, self.elo_ngals_10v15_sf, self.elo_ngals_15v20_sf]
		self.ehi_ngals_all_sf = [self.ehi_ngals_n05v00_sf, self.ehi_ngals_00v05_sf, self.ehi_ngals_05v10_sf, self.ehi_ngals_10v15_sf, self.ehi_ngals_15v20_sf]
		self.phis_all_sf = [self.phis_n05v00_sf, self.phis_00v05_sf, self.phis_05v10_sf, self.phis_10v15_sf, self.phis_15v20_sf]


		##################
		###  QU galaxies
		##################

		self.ngals_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_all_qu = [self.ngals_n05v00_qu, self.ngals_00v05_qu, self.ngals_05v10_qu, self.ngals_10v15_qu, self.ngals_15v20_qu]
		self.elo_ngals_all_qu = [self.elo_ngals_n05v00_qu, self.elo_ngals_00v05_qu, self.elo_ngals_05v10_qu, self.elo_ngals_10v15_qu, self.elo_ngals_15v20_qu]
		self.ehi_ngals_all_qu = [self.ehi_ngals_n05v00_qu, self.ehi_ngals_00v05_qu, self.ehi_ngals_05v10_qu, self.ehi_ngals_10v15_qu, self.ehi_ngals_15v20_qu]
		self.phis_all_qu = [self.phis_n05v00_qu, self.phis_00v05_qu, self.phis_05v10_qu, self.phis_10v15_qu, self.phis_15v20_qu]

class mf_data:
	def __init__(self, name, zclust, sigmaz, voronoi_bin=''):
		self.name = name
		self.zclust = zclust
		self.sigmaz = sigmaz
		self.masslimits = mypy.readcat('../data/completeness/masslimits_%s_master.dat' % name)

		self.data = mypy.readcat('%s/MF_%s_LSS%s.dat' % (data_dir, name, voronoi_bin))
		#self.data_2filter = mypy.readcat('%s/MF_%s_LSS%s_2filter.dat' % (data_dir, name, voronoi_bin))
		self.mf_field = pickle.load(open('../data/MF_%s_field.pickle' % name, 'rb'))
		self.mf_voronoi_slices = pickle.load(open('../data/MF_%s_voronoi_slices.pickle' % name, 'rb'))























##############################################
###  Plot MFs in overdensity bins by summing 
###  across all Voronoi redshift slices.
##############################################

plot_fields = []
plot_fields.append(mf_data('N200',    0.691,  0.027))
plot_fields.append(mf_data('SC1324',  0.755,  0.033))
plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields.append(mf_data('N5281',   0.818,  0.029))
plot_fields.append(mf_data('SC1604',  0.910,  0.029))
plot_fields.append(mf_data('SC0910',  1.110,  0.035))
plot_fields.append(mf_data('SC0849',  1.261,  0.029))

cmap = pylab.cm.cool
voronoi_bins = ['_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']
voronoi_labels_merged = ['0.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']

voronoi_bins = ['_n05v00', '_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['-0.5 < log(1+$\delta_{\mathrm{gal}}$) < 0.0', '0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']

voronoi_bins = ['_n10v05', '_n05v00', '_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['-1.0 < log(1+$\delta_{\mathrm{gal}}$) < -0.5', '-0.5 < log(1+$\delta_{\mathrm{gal}}$) < 0.0', '0.0 < log(1+$\delta_{\mathrm{gal}}$) < 0.5', '0.5 < log(1+$\delta_{\mathrm{gal}}$) < 1.0', '1.0 < log(1+$\delta_{\mathrm{gal}}$) < 1.5', '1.5 < log(1+$\delta_{\mathrm{gal}}$) < 2.0']

###  Field Mass Functions - ZFOURGE 2014
###  /Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table1_TOT.dat
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qu = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]

###  GCLASS
vdB13 = mypy.readcat('../data/table3_vanderBurg+2013.dat')


if True:

	#outname = '../output/massFunctions_from_summed_voronoi_slices2.pdf'
	#pdf = PdfPages(outname)



	lmassbins = plot_fields[0].mf_voronoi_slices.lmassbins
	lmassbars = plot_fields[0].mf_voronoi_slices.lmassbars
	dm = lmassbins[1] - lmassbins[0]

	overdens_bins = plot_fields[0].mf_voronoi_slices.overdens_bins
	overdens_bars = plot_fields[0].mf_voronoi_slices.overdens_bars
	d_overdens = overdens_bins[1] - overdens_bins[0]

	###  master 3d array of number counts. Axes are LSSfield, overdensity bin, massbin, 
	master_number_counts = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_number_counts_sf = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_number_counts_qu = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_number_counts_field = pylab.zeros((len(plot_fields), len(lmassbars)))

	master_volumes = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_volumes_sf = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_volumes_qu = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_volumes_field = pylab.zeros((len(plot_fields), len(lmassbars)))

	master_nslices = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))   # total number of voronoi slices

	for fi in range(len(plot_fields)):
		f =  plot_fields[fi]

		###  have to skip by twos because of overlappipng voronoi zslices
		for zi in range(0, len(f.mf_voronoi_slices.zbars), 2):

			###  skipping if zslice is behind cluster
			#if f.zclust < f.mf_voronoi_slices.zlos[zi] - 0.2: continue

			###  skipping redshifts of choice
			#if f.mf_voronoi_slices.zlos[zi] < 1: continue

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim95)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_ssp)

			###  iterating over massbins, checking for mass completeness
			for mi in range(len(lmassbars)):

				if lmassbins[mi] > lmasslimit:
					master_number_counts_field[fi][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					master_volumes_field[fi][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_n10v05[zi][mi]
					master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_n10v05_Mpc3[zi]

					master_number_counts[fi][1][mi] += f.mf_voronoi_slices.ngals_n05v00[zi][mi]
					master_volumes[fi][1][mi] += f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi]

					master_number_counts[fi][2][mi] += f.mf_voronoi_slices.ngals_00v05[zi][mi]
					master_volumes[fi][2][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]

					master_number_counts[fi][3][mi] += f.mf_voronoi_slices.ngals_05v10[zi][mi]
					master_volumes[fi][3][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

					master_number_counts[fi][4][mi] += f.mf_voronoi_slices.ngals_10v15[zi][mi]
					master_volumes[fi][4][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

					master_number_counts[fi][5][mi] += f.mf_voronoi_slices.ngals_15v20[zi][mi]
					master_volumes[fi][5][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

					###  counting number of voronoi slices in mass-overdens bins
					if f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi] > 0: master_nslices[fi][0][mi] += 1
					if f.mf_voronoi_slices.volumes_00v05_Mpc3[zi] > 0: master_nslices[fi][1][mi] += 1
					if f.mf_voronoi_slices.volumes_05v10_Mpc3[zi] > 0: master_nslices[fi][2][mi] += 1
					if f.mf_voronoi_slices.volumes_10v15_Mpc3[zi] > 0: master_nslices[fi][3][mi] += 1
					if f.mf_voronoi_slices.volumes_15v20_Mpc3[zi] > 0: master_nslices[fi][4][mi] += 1

				if lmassbins[mi] > lmasslimit:
					master_number_counts_sf[fi][0][mi] += f.mf_voronoi_slices.ngals_n10v05_sf[zi][mi]
					master_volumes_sf[fi][0][mi] += f.mf_voronoi_slices.volumes_n10v05_Mpc3[zi]

					master_number_counts_sf[fi][1][mi] += f.mf_voronoi_slices.ngals_n05v00_sf[zi][mi]
					master_volumes_sf[fi][1][mi] += f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi]

					master_number_counts_sf[fi][2][mi] += f.mf_voronoi_slices.ngals_00v05_sf[zi][mi]
					master_volumes_sf[fi][2][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]

					master_number_counts_sf[fi][3][mi] += f.mf_voronoi_slices.ngals_05v10_sf[zi][mi]
					master_volumes_sf[fi][3][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

					master_number_counts_sf[fi][4][mi] += f.mf_voronoi_slices.ngals_10v15_sf[zi][mi]
					master_volumes_sf[fi][4][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

					master_number_counts_sf[fi][5][mi] += f.mf_voronoi_slices.ngals_15v20_sf[zi][mi]
					master_volumes_sf[fi][5][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

				if lmassbins[mi] > lmasslimit_ssp:
					master_number_counts_qu[fi][0][mi] += f.mf_voronoi_slices.ngals_n10v05_qu[zi][mi]
					master_volumes_qu[fi][0][mi] += f.mf_voronoi_slices.volumes_n10v05_Mpc3[zi]

					master_number_counts_qu[fi][1][mi] += f.mf_voronoi_slices.ngals_n05v00_qu[zi][mi]
					master_volumes_qu[fi][1][mi] += f.mf_voronoi_slices.volumes_n05v00_Mpc3[zi]

					master_number_counts_qu[fi][2][mi] += f.mf_voronoi_slices.ngals_00v05_qu[zi][mi]
					master_volumes_qu[fi][2][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]

					master_number_counts_qu[fi][3][mi] += f.mf_voronoi_slices.ngals_05v10_qu[zi][mi]
					master_volumes_qu[fi][3][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

					master_number_counts_qu[fi][4][mi] += f.mf_voronoi_slices.ngals_10v15_qu[zi][mi]
					master_volumes_qu[fi][4][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

					master_number_counts_qu[fi][5][mi] += f.mf_voronoi_slices.ngals_15v20_qu[zi][mi]
					master_volumes_qu[fi][5][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]




	###  plotting results

	fig = pylab.figure(figsize=(15.3, 6.8))
	sp2 = fig.add_subplot(122)
	sp1 = fig.add_subplot(121)

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number / Mpc$^{3}$ / dex$^{}$')

	sp1.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])
	sp2.axis([8.4, 11.9, 0.5, 1000])

	fig.subplots_adjust(left=0.09)

	ydatas = []
	fits_tot, covs_tot = [], []
	line_fits_tot, line_covs_tot = [], []
	mlim_hi = 11.55



	###  deciding which overdensity bins to plot
	o_range = range(1, len(overdens_bars))
	offie = pylab.linspace(-len(o_range)*0.01/2, len(o_range)*0.01/2, len(o_range))

	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		inds0 = pl.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)
		sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)





		###  fitting (d)schechter functions and plotting
		try:

			xdata, ydata = lmassbars[inds], (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0))[vi][inds]
			ydata_errs = 10**(elo+ehi)[inds]/2.
			ydata_errs = (elo+ehi)[inds]/2.


			###  fitting log-linear relation to ratio
			fit, cov = optimize.curve_fit(line, xdata, pylab.log10(y[inds] / ydatas[0][inds]), p0=[-8, 0.5], sigma=(elo+ehi)[inds]/2.)
			line_fits_tot.append(fit)
			line_covs_tot.append(cov.diagonal()**0.5)


			###  fitting schechter funtion to SMF
			#fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1., 11., 10**-3], sigma=ydata_errs, maxfev=10**5)
			#ymod = schechter_mf(xmod, *fit)

			###  fitting dschechter funtion to SMF
			p0 = [11, -1., -1., 10**-3, 10**-4]
			p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
			fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=p0, sigma=ydata_errs, maxfev=10**5)
			ymod = dschechter(xmod, *fit)
			a = sp1.plot(xmod + offie[vi_counter], ymod / dm, color=color_i, lw=1)
			fits_tot.append(fit)
			covs_tot.append(cov.diagonal()**0.5)


		except:
			pass


	sp1.legend(loc=3, numpoints=1, fontsize=14)

















#############################
###  Plotting toy model SMFs
#############################

simdata = pickle.load(open('../data/merger_simulation_scenario6_7m12.pickle', 'rb'))


inds_lmass = pylab.find((simdata.lmassax > 8.5) & (simdata.lmassax < 11.55))


###  figure for paper
cmap = pylab.cm.jet
#sp_cbar = add_inset(sp2, [0.6, 0.85, 0.33, 0.05])
#sp_cbar2 = add_inset(sp2, [0.6, 0.65, 0.33, 0.05])

sp2.axis([8.1, 11.8, 5*10**1, 3*10**5])
sp2.axis([8, 11.8, 5*10**1, 2*10**5])

#sp2.grid()
sp2.minorticks_on()
sp2.set_yscale('log')
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp2.set_ylabel('Number')

#sp2.plot(lmass, ngal0, color='k', lw=7, label='initial state', zorder=99)


#sp_cbar.set_title('$M_{*,\,\mathrm{ICL}} \, / \, M_{*,\,\mathrm{total}}$', fontsize=18, position=[.5, 1.25])
#sp_cbar.set_yticks([])
#imdatx, imdaty = pylab.meshgrid(pylab.linspace(simdata.f_icm.min(), simdata.f_icm.max(), int(simdata.f_icm.max()*100)+1), [0, 0])
#sp_cbar.imshow(imdatx, cmap=cmap)
#sp_cbar.set_aspect('auto')

#xticks = sp_cbar.get_xticks()

#sp_cbar.set_xticks([0, 6, 12, 18])
#sp_cbar.set_xticklabels(['0', '0.06', '0.12', '0.18'], fontsize=15)



###  f_icm to plot
f_icms = pylab.arange(0.02, 0.17, 0.02)
f_icms = pylab.arange(0.03, 0.17, 0.02)
f_mergeds = pylab.zeros(len(f_icms))

for i_f, f in enumerate(f_icms[:-1]):

	weights = 1. / abs(simdata.f_icm - f)**5
	smf = pylab.average(simdata.ngals_icm, axis=0, weights=weights)
	f_mergeds[i_f] = pylab.average(simdata.fraction_merged, weights=weights)

	c = f / simdata.f_icm.max()
	sp2.plot(simdata.lmassax[inds_lmass], smf[inds_lmass], marker='o', ms=4, mew=0.5, color=cmap(c))



	font = FontProperties()
	font.set_family('sans-serif')
	#if i_f == 0:
		#sp2.text(8.65-0.2, smf[inds_lmass][0]*2, 'fraction\nmerged', fontproperties=font,
		#	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)

	t = sp2.text(8.65-0.05, smf[inds_lmass][0], '%.1f%%' % (f_mergeds[i_f]*100), fontproperties=font,
		          verticalalignment='center', horizontalalignment='right', fontsize=15, fontweight='normal')


smf = simdata.ngals_icm[0]

c = 0.
sp2.plot(simdata.lmassax[inds_lmass], smf[inds_lmass],
	          marker='', lw=6, color=cmap(c))


font = FontProperties()
font.set_family('sans-serif')
sp2.text(8.65-0.2, smf[inds_lmass][0]*2.3, 'fraction\nmerged', fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)
sp2.text(8.65-0.05, smf[inds_lmass][0], '%i%%' % (0), fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='right', fontsize=15, fontweight='normal')











