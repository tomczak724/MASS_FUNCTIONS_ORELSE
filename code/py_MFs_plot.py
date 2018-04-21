
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




#################################
###  Defining functions etc.  ###
#################################

if True:

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
























#######################################################
###  Plot full MFs in 0.65<z<0.85 and 0.85<z<1.2 bins
###  Compare to EDisCS and GCLASS
#######################################################

plot_fields = []
plot_fields.append(mf_data('N200',    0.691,  0.027))
plot_fields.append(mf_data('SC1324',  0.755,  0.033))
plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields.append(mf_data('N5281',   0.818,  0.029))

plot_fields.append(mf_data('SC1604',  0.910,  0.029))
plot_fields.append(mf_data('SC0910',  1.110,  0.035))
plot_fields.append(mf_data('SC0849',  1.261,  0.029))

lmassbars = plot_fields[0].data.lmass
dm = lmassbars[1] - lmassbars[0]

if False:

	fig = pylab.figure(figsize=(11.9, 6.8))
	sp2 = fig.add_subplot(122)
	sp1 = fig.add_subplot(121)

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number  [dex$^{-1}$]')

	fig.subplots_adjust(wspace=0, left=0.09)

	t = sp1.text(0.55, 0.9, '0.65 < z < 0.85', transform=sp1.transAxes)
	t = sp2.text(0.55, 0.9, '0.85 < z < 1.3', transform=sp2.transAxes)

	sp1.axis([9.2, 11.95, 5, 5000])
	sp2.axis([9.2, 11.95, 5, 5000])


	ntot_zbin1 = pylab.zeros(len(lmassbars))
	ntot_zbin2 = pylab.zeros(len(lmassbars))


	for i in range(len(plot_fields)):

		f = plot_fields[i]

		if 0.65 < f.zclust < 0.85:
			ntot_zbin1 += f.data.Ntot
			print f.name, 1
		elif 0.85 < f.zclust < 1.3:
			ntot_zbin2 += f.data.Ntot
			print f.name, 2



	nhinlo_zbin1 = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in ntot_zbin1])
	nhinlo_zbin2 = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in ntot_zbin2])


	sp1.errorbar(lmassbars, ntot_zbin1/dm, xerr=dm/2., yerr=[nhinlo_zbin1[:,1]/dm, nhinlo_zbin1[:,0]/dm],
		         ls='', marker='o', ms=9, mfc='k', ecolor='k', zorder=10, label='ORELSE')

	sp2.errorbar(lmassbars, ntot_zbin2/dm, xerr=dm/2., yerr=[nhinlo_zbin2[:,1]/dm, nhinlo_zbin2[:,0]/dm],
		         ls='', marker='o', ms=9, mfc='k', ecolor='k', zorder=10)





	###  ploting results from literature

	v13 = mypy.readcat('../data/figure2_Vulcani+2013.dat')
	vdB13 = mypy.readcat('../data/table3_vanderBurg+2013.dat')


	###  Vulcani + 2013
	x, y, ylo, yhi = v13.lmass, v13.logphi, v13.logphi-v13.logphi_lo, v13.logphi_hi-v13.logphi
	y = 10**y

	massax_for_scale = pylab.arange(10.3, 11.51, 0.2)
	v13_y = pylab.interp(massax_for_scale, v13.lmass, y)
	t16_y = pylab.interp(massax_for_scale, lmassbars, ntot_zbin1)
	weights = 1. / pylab.interp(massax_for_scale, v13.lmass, ylo)

	scale_factor = pylab.average(t16_y/v13_y, weights=weights)

	c = '#ff6666'
	y = y * scale_factor / dm
	ylo, yhi = y * ylo, y * yhi
	sp1.errorbar(x+0.01, y, yerr=[ylo, yhi], ls='', marker='^', ms=8, mew=1.5, mfc='w',
		         mec=c, ecolor=c, elinewidth=1.5, zorder=1, label='EDisCS: 0.4 < z < 0.8')

	leg1 = sp1.legend(loc=3, numpoints=1, fontsize=15)




	###  van der Burg + 2013
	x, y, ylo, yhi = vdB13.lmass, vdB13.ntot, vdB13.nlo_tot, vdB13.nhi_tot

	massax_for_scale = pylab.arange(10.1, 11.51, 0.2)
	vdB13_y = pylab.interp(massax_for_scale, vdB13.lmass, y)
	t16_y = pylab.interp(massax_for_scale, lmassbars, ntot_zbin2)
	weights = pylab.interp(massax_for_scale, vdB13.lmass, y/ylo)

	scale_factor = pylab.average(t16_y/vdB13_y, weights=weights)

	c = '#6666ff'
	y = y * scale_factor / dm
	ylo, yhi = ylo * scale_factor / dm, yhi * scale_factor / dm
	sp2.errorbar(x+0.01, y, yerr=[ylo, yhi], ls='', marker='s', ms=6, mew=1.5, mfc='w',
		         mec=c, ecolor=c, elinewidth=1.5, zorder=1, label='GCLASS: 0.85 < z < 1.2')

	leg2 = sp2.legend(loc=3, numpoints=1, fontsize=15)

	pylab.savefig('../output/massfunctions_EDisCS+GCLASS.png')
	pylab.close()
	print 'wrote to: ../output/massfunctions_EDisCS+GCLASS.png'













################################################
###  Plot full MFs of each field independently
################################################

plot_fields = []
plot_fields.append(mf_data('N200',    0.691,  0.027))
plot_fields.append(mf_data('SC1324',  0.755,  0.033))
plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields.append(mf_data('N5281',   0.818,  0.029))
plot_fields.append(mf_data('SC1604',  0.910,  0.029))
plot_fields.append(mf_data('SC0910',  1.110,  0.035))
plot_fields.append(mf_data('SC0849',  1.261,  0.029))

lmassbars = plot_fields[0].data.lmass
dm = lmassbars[1] - lmassbars[0]

if False:

	fig = pylab.figure(figsize=(7., 6.8))
	sp1 = fig.add_subplot(111)

	for i in range(len(plot_fields)):

		f = plot_fields[i]

		sp1.minorticks_on()
		sp1.set_yscale('log')

		sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number  [dex$^{-1}$]')

		t = sp1.text(0.7, 0.9, f.name, fontweight='bold', transform=sp1.transAxes)

		sp1.axis([9.2, 11.85, 5, 2000])

		ntot = f.data.Ntot

		nhinlo = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in ntot])

		sp1.errorbar(lmassbars, ntot/dm, xerr=dm/2., yerr=[nhinlo[:,1]/dm, nhinlo[:,0]/dm],
			         ls='', marker='o', ms=9, mfc='k', ecolor='k', zorder=10, label=f.name)

		outname = '../output/MF_%s.png' % f.name
		fig.savefig(outname)
		print 'wrote to: %s' % outname
		sp1.clear()

	pylab.close()
	print ''














################################################
###  Plot field (ie non-LSS) MFs of each field
################################################

plot_fields = []
plot_fields.append(mf_data('N200',    0.691,  0.027))
plot_fields.append(mf_data('SC1324',  0.755,  0.033))
plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields.append(mf_data('N5281',   0.818,  0.029))
plot_fields.append(mf_data('SC1604',  0.910,  0.029))
#plot_fields.append(mf_data('CL1429',  0.920,  0.051))
plot_fields.append(mf_data('SC0910',  1.110,  0.035))
plot_fields.append(mf_data('SC0849',  1.261,  0.029))

###  Field Mass Functions - ZFOURGE 2014
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'

zfourge_massfunc_tot = mypy.readcat('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = mypy.readcat('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = mypy.readcat('%s/table1_QUI.dat' % zfourge_dir)

zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]

if False:

	outname = '../output/massFunctions_field.pdf'
	pdf = PdfPages(outname)

	fig = pylab.figure(figsize=(17., 13.5))
	sp3 = fig.add_subplot(2, 3, 3)
	sp6 = fig.add_subplot(2, 3, 6)
	sp2 = fig.add_subplot(2, 3, 2)
	sp5 = fig.add_subplot(2, 3, 5)
	sp1 = fig.add_subplot(2, 3, 1)
	sp4 = fig.add_subplot(2, 3, 4)
	sps = [sp1, sp2, sp3, sp4, sp5, sp6]

	fig.subplots_adjust(wspace=0, hspace=0, top=0.94, left=0.09, bottom=0.09)

	dschechter_fits = []
	dschechter_errs = []

	volumes_all = pylab.zeros((len(plot_fields[0].mf_field.zbars), len(plot_fields)))
	ngal_bins_all = pylab.zeros((len(plot_fields), len(plot_fields[0].mf_field.zbars), len(plot_fields[0].mf_field.lmassbars)))

	for fi in range(len(plot_fields)):

		mypy.progress_bar(fi, len(plot_fields))

		f = plot_fields[fi]
		title = '%s :   z$_{LSS}$ = %.3f' % (f.name, f.zclust)
		sp2.set_title(title)

		for zi in range(len(f.mf_field.zbars)):

			zlo, zhi = f.mf_field.zlos[zi], f.mf_field.zhis[zi]
			volume = massfunc.vcomoving_slice(f.mf_field.area, zlo, zhi)


			dm = f.mf_field.lmassbars[1] - f.mf_field.lmassbars[0]
			phi_bins = f.mf_field.number_counts[zi] / volume / dm
			ephi_lo = phi_bins * f.mf_field.nlo_poisson[zi] / f.mf_field.number_counts[zi]
			ephi_lo[pylab.isnan(ephi_lo)] = 0
			ephi_hi = phi_bins * f.mf_field.nhi_poisson[zi] / f.mf_field.number_counts[zi]
			ephi_hi[pylab.isnan(ephi_hi)] = (1. / volume / dm) * f.mf_field.nhi_poisson[zi][pylab.isnan(ephi_hi)] / 1.


			###  adding galaxies to total count if the LSS is not in this zbin
			zphot_lo = f.zclust - 1. * f.sigmaz * (1 + f.zclust)
			zphot_hi = f.zclust + 1. * f.sigmaz * (1 + f.zclust)
			if not (zlo < zphot_lo < zhi or zlo < zphot_hi < zhi):
				volumes_all[zi][fi] += volume
				ngal_bins_all[fi][zi] += f.mf_field.number_counts[zi]


			###  plotting it up!!!
			sp = sps[zi]
			sp.grid()
			sp.minorticks_on()
			axis_lims = [7.7, 11.8, -5.4, -0.9]
			sp.axis(axis_lims)
			sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
			sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
			sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)

			if (zlo < zphot_lo < zhi) or (zlo < zphot_hi < zhi):
				sp.fill_between(axis_lims[0:2], axis_lims[2], axis_lims[3], color='#ffe6e6', zorder=0)


			sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


			###  plotting zfourge first
			i_mf = 3*zi + 1
			inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
			sp.fill_between(zfourge_massfunc_tot[inds, 0],
				            zfourge_massfunc_tot[inds, i_mf] + zfourge_massfunc_tot[inds, i_mf+1], 
				            zfourge_massfunc_tot[inds, i_mf] - zfourge_massfunc_tot[inds, i_mf+2], 
				            color='#bfbfbf', zorder=1)
			sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')
	

			###  plotting orelse next
			y = pylab.log10(phi_bins)
			ylo = pylab.log10(phi_bins) - pylab.log10(phi_bins - ephi_lo)
			yhi = pylab.log10(phi_bins + ephi_hi) - pylab.log10(phi_bins)
			sp.errorbar(f.mf_field.lmassbars, y, xerr=dm/2., yerr=[ylo, yhi], ls='', marker='o', ms=7, mew=1.5, \
				        elinewidth=1.5, mfc='#cc00ff', ecolor='k', label='ORELSE today', zorder=2)

			if zi == 3:
				sp.legend(loc=3, fontsize=18, title=title)


			###  fitting double Schechter
			minds = pylab.find(f.mf_field.lmassbars >= 9.5)
			x = f.mf_field.lmassbars[minds]
			y = phi_bins[minds]
			yerr = (ephi_lo[minds] + ephi_hi[minds]) / 2.
			fit, cov = optimize.curve_fit(massfunc.dschechter_mf, x, y, p0=[10.78, -0.98, -2.54, -1.9, -4.29], sigma=yerr)

			xax = pylab.linspace(9.5, 11.6, 100)
			yax = massfunc.dschechter_mf(xax, *fit)
			sp.plot(xax, pylab.log10(yax), color='#cc00ff', lw=1.5)




		pdf.savefig()
		for sp in sps:
			sp.clear()

	print '\n'



	###  plotting all fields averaged together
	for zi in range(len(f.mf_field.zbars)):

		zlo, zhi = f.mf_field.zlos[zi], f.mf_field.zhis[zi]
		title = 'All fields minus LSS'



		ngal = ngal_bins_all.sum(axis=0)[zi]
		vol = volumes_all[zi].sum()
		nlo_poisson, nhi_poisson = [], []
		for n in ngal:
			nhi, nlo = mypy.massfunc.confidence_interval(n)
			nlo_poisson.append(nlo)
			nhi_poisson.append(nhi)
		nlo_poisson, nhi_poisson = pylab.array(nlo_poisson), pylab.array(nhi_poisson)

		phi = ngal * 1. / vol / dm

		ephi_lo = phi * nlo_poisson / ngal
		ephi_lo[pylab.isnan(ephi_lo)] = 0
		ephi_hi = phi * nhi_poisson / ngal
		ephi_hi[pylab.isnan(ephi_hi)] = (1. / vol / dm) * nhi_poisson[pylab.isnan(ephi_hi)] / 1.





		sp = sps[zi]
		sp.grid()
		sp.minorticks_on()
		axis_lims = [7.7, 11.8, -5.4, -0.9]
		sp.axis(axis_lims)
		sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
		sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)


		sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


		###  plotting zfourge first
		i_mf = 3*zi + 1
		inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
		sp.fill_between(zfourge_massfunc_tot[inds, 0],
			            zfourge_massfunc_tot[inds, i_mf] + zfourge_massfunc_tot[inds, i_mf+1], 
			            zfourge_massfunc_tot[inds, i_mf] - zfourge_massfunc_tot[inds, i_mf+2], 
		    	        color='#bfbfbf', zorder=1)
		sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')


		###  plotting orelse next
		y = pylab.log10(phi)
		ylo = pylab.log10(phi) - pylab.log10(phi - ephi_lo)
		yhi = pylab.log10(phi + ephi_hi) - pylab.log10(phi)
		sp.errorbar(f.mf_field.lmassbars, y, xerr=dm/2., yerr=[ylo, yhi], ls='', marker='o', ms=7, mew=1.5, \
			        elinewidth=1.5, mfc='k', ecolor='k', label='ORELSE today', zorder=2)

		if zi == 3:
			sp.legend(loc=3, fontsize=18, title=title)


		###  fitting double Schechter
		minds = pylab.find(f.mf_field.lmassbars >= 9.5)
		x = f.mf_field.lmassbars[minds]
		y = phi[minds]
		yerr = (ephi_lo[minds] + ephi_hi[minds]) / 2.
		fit, cov = optimize.curve_fit(massfunc.dschechter_mf, x, y, p0=[10.78, -0.98, -2.54, -1.9, -4.29], sigma=yerr)
		dschechter_fits.append(fit)
		dschechter_errs.append(cov.diagonal()**0.5)

		xax = pylab.linspace(9.5, 11.6, 100)
		yax = massfunc.dschechter_mf(xax, *fit)
		sp.plot(xax, pylab.log10(yax), color='r', lw=1.5)



	pdf.savefig()
	pdf.close()
	pylab.close()
	print 'wrote to: %s' % outname





	###  plotting schechter params vs redshift
	dschechter_fits = pylab.array(dschechter_fits)
	dschechter_errs = pylab.array(dschechter_errs)

	zfourge_params = mypy.readcat('/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table2_TOT.dat')
	inds = pylab.arange(len(f.mf_field.zbars))


	fig = pylab.figure(figsize=(23, 8.5))
	sp1 = fig.add_subplot(131)
	sp2 = fig.add_subplot(132)
	sp3 = fig.add_subplot(133)

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp3.minorticks_on()

	sp1.grid()
	sp2.grid()
	sp3.grid()

	fig.subplots_adjust(top=0.94, left=0.08)



	###  Mstar
	sp1.set_xlabel('redshfit')
	sp1.set_ylabel('log( M$^*$ / M$_{\odot}$ )')
	sp1.axis([-0.2, 2.3, 10.35, 11.4])

	y = dschechter_fits[:,0]
	yerr = dschechter_errs[:,0]
	sp1.errorbar(f.mf_field.zbars + 0.015, y, xerr=f.mf_field.dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
		         ecolor='r', elinewidth=2, label='ORELSE today')


	sp1.errorbar(f.mf_field.zbars - 0.015, zfourge_params.Mstar[inds], xerr=f.mf_field.dz/2., yerr=zfourge_params.e_Mstar[inds], \
		         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')






	###  alpha2
	sp2.set_xlabel('redshfit')
	sp2.set_ylabel('alpha 2')
	sp2.axis([-0.2, 2.3, -2.2, -1])

	y = dschechter_fits[:,3]
	yerr = dschechter_errs[:,3]
	sp2.errorbar(f.mf_field.zbars + 0.015, y, xerr=f.mf_field.dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
		         ecolor='r', elinewidth=2, label='ORELSE today')


	sp2.errorbar(f.mf_field.zbars - 0.015, zfourge_params.alpha2[inds], xerr=f.mf_field.dz/2., yerr=zfourge_params.e_alpha2[inds], \
		         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')







	###  log(phistar1 + phistar2)
	sp3.set_xlabel('redshfit')
	sp3.set_ylabel('log( $\Phi_1^*$ + $\Phi_2^*$ )')
	sp3.axis([-0.2, 2.3, -3.4, -2.1])

	y = pylab.log10(10**dschechter_fits[:,2] + 10**dschechter_fits[:,4])

	dphi1 = (10**dschechter_errs[:,2] - 1) * 10**dschechter_fits[:,2]
	dphi2 = (10**dschechter_errs[:,4] - 1) * 10**dschechter_fits[:,4]
	yerr = pylab.sqrt((dphi1**2 + dphi2**2) / pylab.log(10)**2 / (10**dschechter_fits[:,2] + 10**dschechter_fits[:,4])**2)

	sp3.errorbar(f.mf_field.zbars + 0.015, y, xerr=f.mf_field.dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
		         ecolor='r', elinewidth=2, label='ORELSE today')


	y = pylab.log10(10**zfourge_params.phistar1 + 10**zfourge_params.phistar2)

	dphi1 = (10**zfourge_params.e_phistar1 - 1) * 10**zfourge_params.phistar1
	dphi2 = (10**zfourge_params.e_phistar2 - 1) * 10**zfourge_params.phistar2
	yerr = pylab.sqrt((dphi1**2 + dphi2**2) / pylab.log(10)**2 / (10**zfourge_params.phistar1 + 10**zfourge_params.phistar2)**2)

	sp3.errorbar(f.mf_field.zbars - 0.015, y[inds], xerr=f.mf_field.dz/2., yerr=yerr[inds], \
		         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')

	sp3.legend(loc=3, fontsize=18)


	outname = '../output/schechter_params.png'
	fig.savefig(outname)
	pylab.close()
	print 'wrote to: %s\n' % outname


















#######################################################
###  Plot full MFs in 4 equal-size bins of Voronoi
###  overdensity between 0 < log(1+delta) < 2
#######################################################

cmap = pylab.cm.cool
voronoi_bins = ['_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['0.0 < log(1+$\delta$) < 0.5', '0.5 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 1.5', '1.5 < log(1+$\delta$) < 2.0']

lmassbars = plot_fields[0].data.lmass
dm = lmassbars[1] - lmassbars[0]

if False:

	##################################
	###  ALL galaxies in voronoi bins
	##################################
	nsims = 1*10**3
	print '\n\nALL galaxies\n\n'
	if True:
		fig = pylab.figure(figsize=(13.5, 7.))
		sp2 = fig.add_subplot(122)
		sp1 = fig.add_subplot(121)

		sp1.minorticks_on()
		sp2.minorticks_on()
		sp1.set_yscale('log')

		sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number  [dex$^{-1}$]')
		sp2.set_xlabel('log( M$^*$ / M$_{\odot}$ )')
		sp2.set_ylabel(r'$\alpha$', size=26)

		fig.subplots_adjust(left=0.09)

		sp1.grid()
		sp2.grid()

		sp1.axis([9.2, 11.85, 20, 5000])
		sp2.axis([10.5, 11.55, -1.55, -0.3])
		sp1.axis([9.2, 11.85, 5, 5000])
		sp2.axis([10.5, 11.55, -1.65, -0.2])


		###  numbers of galaxies in each mass-overdensity bin for each field
		ngals_zspec = pylab.zeros((len(voronoi_bins), 8, len(lmassbars)))
		ngals_zphot = pylab.zeros((len(voronoi_bins), 8, len(lmassbars)))

		###  survey areas for each fied and each vbin
		areas_vbin = pylab.zeros((8, len(voronoi_bins)))


		offie = [-0.015, -0.005, 0.005, 0.015]
		for vi in range(len(voronoi_bins)):

			plot_fields = []
			plot_fields.append(mf_data('N200',    0.691,  0.027, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC1324',  0.755,  0.033, voronoi_bins[vi]))
			plot_fields.append(mf_data('RCS0224', 0.772,  0.027, voronoi_bins[vi]))
			plot_fields.append(mf_data('RXJ1716', 0.813,  0.021, voronoi_bins[vi]))
			plot_fields.append(mf_data('N5281',   0.818,  0.029, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC1604',  0.910,  0.029, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC0910',  1.110,  0.035, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC0849',  1.261,  0.029, voronoi_bins[vi]))

			ntot = pylab.zeros(len(lmassbars))
			nlos = pylab.zeros((len(plot_fields), len(lmassbars)))
			nhis = pylab.zeros((len(plot_fields), len(lmassbars)))

			for fi in range(len(plot_fields)):
				f = plot_fields[fi]
				ntot += f.data.Ntot
				nlos[fi] += f.data.Nlo_tot
				nhis[fi] += f.data.Nhi_tot

				ngals_zspec[vi][fi] += f.data.Nzspec
				ngals_zphot[vi][fi] += f.data.Nzphot_corr

				#areas_vbin[fi][vi] += 



			nlos = (nlos**2).sum(axis=0)**0.5
			nhis = (nhis**2).sum(axis=0)**0.5

			color_i = cmap(vi/(len(voronoi_bins)-1.))

			sp1.errorbar(lmassbars+offie[vi], ntot/dm, xerr=dm/2., yerr=[nlos/dm, nhis/dm],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, ecolor=color_i, label=voronoi_labels[vi])



			###  fitting single-Schecther function
			xdata, ydata = lmassbars, ntot/dm
			ydata_errs = (nlos + nhis) / 2. / dm

			fit0, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1, 10.5, pylab.median(ydata)], sigma=ydata_errs)
			sp2.plot(fit0[1], fit0[0], ls='', marker='o', mew=3, ms=9, mec='k', mfc=color_i, zorder=10)

			xax = pylab.linspace(xdata.min(), xdata.max(), 100)
			yax = schechter_mf(xax, *fit0)
			sp1.plot(xax, yax, lw=1.5, color=color_i)

			'''

			schmodel = single_schechter_model(lmstar=10.8, alpha=-1., phistar=200.)
			schfit = fit_algorithm(schmodel, xdata, ydata, weights=1./ydata_errs**2)
			sp1.plot(xax, schfit(xax), lw=1.5, color=color_i)
			sp2.plot(schfit.lmstar.value, schfit.alpha.value, ls='', marker='o', mew=3, ms=9, mec='k', mfc=color_i, zorder=10)

			'''


			###  running MC simulations
			print '\nrunning %i MC simulations for %s' % (nsims, voronoi_labels[vi])
			xdata, ydata = lmassbars, ntot
			ydata_errs = (nlos + nhis) / 2.

			mc_ntot = pylab.poisson(lam=ydata, size=(nsims, len(ydata)))
			bestfits = []

			for mci in range(len(mc_ntot)):

				mypy.progress_bar(mci, nsims)

				ydata_mci = mc_ntot[mci]
				fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata_mci, p0=[-1, 10.5, pylab.median(ydata)], sigma=ydata_errs)
				bestfits.append(fit)

				#schmodel = single_schechter_model(lmstar=10.8, alpha=-1., phistar=200.)
				#schfit = fit_algorithm(schmodel, xdata, ydata_mci, weights=ydata_errs)
				#bestfits.append([schfit.alpha.value, schfit.lmstar.value, schfit.phistar.value])



			bestfits = pylab.array(bestfits)
			print ''

			### plotting results from MC sims
			#sp2.plot(bestfits[:,1], bestfits[:,0], 'ko', ms=1, mec=color_i, zorder=1)

			gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
			xlo, xhi = 10.5, 11.8
			ylo, yhi = -1.7, -0.1
			nbins = 100
			xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
			hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

			levels = hist2d.max() * gauss_levels[0:2]

			contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, colors=[color_i, color_i], linewidths=[3, 1])


		leg1 = sp1.legend(loc=3, numpoints=1, fontsize=14)




		###  values from GCLASS (van der Burg+2013)
		###  table 4: Total-cluster
		vdb_alpha = -0.46
		vdb_alpha_err = [[0.26], [0.08]]
		vdb_mstar = 10.72
		vdb_mstar_err = [[0.02], [0.09]]

		asdf = sp2.errorbar(vdb_mstar, vdb_alpha, xerr=vdb_mstar_err, yerr=vdb_alpha_err, 
			                ls='', marker='x', ms=10, mec='k', mew=2, mfc='k', ecolor='k', label='GCLASS: 0.85<z<1.2')


		###  values from EDisCS (Vulcani+2013)
		###  figure 2, table 5: EDisCS cluster regions 
		V_alpha = -1.03
		V_alpha_err = 0.08
		V_mstar = 11.15
		V_mstar_err = 0.07

		asdf = sp2.errorbar(V_mstar, V_alpha, xerr=V_mstar_err, yerr=V_alpha_err, 
			                ls='', marker='D', ms=8, mec='k', mew=2, mfc='gray', ecolor='k', label='EDisCS: 0.4<z<0.8')



		###  values from COSMOS-Xray (Giodini+2012)
		###  table 3: high-mass groups, 0.8<z<1.0
		g_alpha = -1.21
		g_alpha_err = [[0.08], [0.08]]
		g_mstar = 10.77
		g_mstar_err = [[0.12], [0.12]]

		asdf = sp2.errorbar(g_mstar, g_alpha, xerr=g_mstar_err, yerr=g_alpha_err, 
			                ls='', marker='s', ms=8, mec='k', mew=2, mfc='g', ecolor='k', label='COSMOS-Xray: 0.8<z<1.0')

		sp2.legend(loc=3, numpoints=1, fontsize=15)

		outname = '../output/MF_voronoiBins.png'
		fig.savefig(outname)
		print 'wrote to: %s\n' % outname
		pylab.close()





		############################################################
		###  plots of sepc completeness and fraction in vbins vs M*
		############################################################

		fig = pylab.figure(figsize=(8., 9.))
		sp1 = pylab.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
		sp2 = pylab.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
		fig.subplots_adjust(hspace=0, wspace=0)

		sp1.set_ylabel('spec. completeness')
		sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp2.set_ylabel('Fraction')

		sp1.grid()
		sp2.grid()
		sp1.minorticks_on()
		sp2.minorticks_on()

		sp1.axis([9.2, 11.95, -0.05, 1.05])
		sp2.axis([9.2, 11.95, -0.05, 1.05])


		for vi in range(len(voronoi_bins)):

			color_i = cmap(vi/(len(voronoi_bins)-1.))


			###  spec completeness
			nzspec = pylab.sum(ngals_zspec[vi], axis=0)
			nzphot = pylab.sum(ngals_zphot[vi], axis=0)

			sp1.plot(lmassbars, nzspec / (nzspec + nzphot),
				     ls='-', marker='o', ms=7, mew=1.5, lw=2, color=color_i, label=voronoi_labels[vi])


			###  fractions
			n_vbin = pylab.sum(ngals_zspec[vi], axis=0) + pylab.sum(ngals_zphot[vi], axis=0)
			n_tot = pylab.sum((pylab.sum(ngals_zspec, axis=1) + pylab.sum(ngals_zphot, axis=1)), axis=0)

			sp2.plot(lmassbars, n_vbin / n_tot,
				     ls='-', marker='s', ms=7, mew=1.5, lw=2, color=color_i, label=voronoi_labels[vi])


		sp1.legend(loc=2, fontsize=13)

		outname2 = '../output/spec_completeness_vbins.png'
		pylab.savefig(outname2)
		pylab.close()
		print '\nwrote to:   %s' % outname2












	#################################
	###  SF galaxies in voronoi bins
	#################################
	print '\n\nSF galaxies\n\n'
	if True:
		fig = pylab.figure(figsize=(13.5, 7.))
		sp2 = fig.add_subplot(122)
		sp1 = fig.add_subplot(121)

		sp1.minorticks_on()
		sp2.minorticks_on()
		sp1.set_yscale('log')

		sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number  [dex$^{-1}$]')
		sp2.set_xlabel('log( M$^*$ / M$_{\odot}$ )')
		sp2.set_ylabel(r'$\alpha$', size=26)

		fig.subplots_adjust(left=0.09)

		sp1.grid()
		sp2.grid()

		sp1.axis([9.2, 11.85, 20, 5000])
		sp2.axis([10.5, 11.55, -1.55, -0.3])
		sp1.axis([9.2, 11.85, 5, 5000])
		sp2.axis([10.5, 11.55, -1.65, -0.2])

		offie = [-0.015, -0.005, 0.005, 0.015]
		for vi in range(len(voronoi_bins)):

			plot_fields = []
			plot_fields.append(mf_data('N200',    0.691,  0.027, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC1324',  0.755,  0.033, voronoi_bins[vi]))
			plot_fields.append(mf_data('RCS0224', 0.772,  0.027, voronoi_bins[vi]))
			plot_fields.append(mf_data('RXJ1716', 0.813,  0.021, voronoi_bins[vi]))
			plot_fields.append(mf_data('N5281',   0.818,  0.029, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC1604',  0.910,  0.029, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC0910',  1.110,  0.035, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC0849',  1.261,  0.029, voronoi_bins[vi]))

			nsf = pylab.zeros(len(lmassbars))
			nlos = pylab.zeros((len(plot_fields), len(lmassbars)))
			nhis = pylab.zeros((len(plot_fields), len(lmassbars)))

			for i in range(len(plot_fields)):
				f = plot_fields[i]
				nsf += f.data.Nsf
				nlos[i] += f.data.Nlo_sf
				nhis[i] += f.data.Nhi_sf

			nlos = (nlos**2).sum(axis=0)**0.5
			nhis = (nhis**2).sum(axis=0)**0.5

			color_i = cmap(vi/(len(voronoi_bins)-1.))

			sp1.errorbar(lmassbars+offie[vi], nsf/dm, xerr=dm/2., yerr=[nlos/dm, nhis/dm],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, ecolor=color_i, label=voronoi_labels[vi])



			###  fitting single-Schecther function
			xdata, ydata = lmassbars, nsf/dm
			ydata_errs = (nlos + nhis) / 2. / dm
			fit0, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1, 10.5, pylab.median(ydata)], sigma=ydata_errs)
			sp2.plot(fit0[1], fit0[0], ls='', marker='o', mew=3, ms=9, mec='k', mfc=color_i, zorder=10)

			xax = pylab.linspace(xdata.min(), xdata.max(), 100)
			yax = schechter_mf(xax, *fit0)
			sp1.plot(xax, yax, lw=1.5, color=color_i)



			###  running MC simulations
			print '\nrunning %i MC simulations for %s' % (nsims, voronoi_labels[vi])
			xdata, ydata = lmassbars, nsf
			ydata_errs = (nlos + nhis) / 2.

			mc_nsf = pylab.poisson(lam=ydata, size=(nsims, len(ydata)))
			bestfits = []

			for mci in range(len(mc_nsf)):

				mypy.progress_bar(mci, nsims)

				ydata_mci = mc_nsf[mci]
				try:
					fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata_mci, p0=[-1, 10.5, pylab.median(ydata)], sigma=ydata_errs)
					bestfits.append(fit)
				except:
					pass

			bestfits = pylab.array(bestfits)
			print ''

			### plotting results from MC sims
			#sp2.plot(bestfits[:,1], bestfits[:,0], 'ko', ms=1, mec=color_i, zorder=1)

			gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
			xlo, xhi = 10.5, 11.8
			ylo, yhi = -1.7, -0.1
			nbins = 100
			xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
			hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

			levels = hist2d.max() * gauss_levels[0:2]

			contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, colors=[color_i, color_i], linewidths=[3, 1])


		leg1 = sp1.legend(loc=3, numpoints=1, fontsize=14)

	outname = '../output/MF_voronoiBins_sf.png'
	fig.savefig(outname)
	print 'wrote to: %s\n' % outname
	pylab.close()












	#################################
	###  QU galaxies in voronoi bins
	#################################
	lmassbars = plot_fields[0].data.lmass[2:]
	print '\n\nQU galaxies\n\n'
	if True:
		fig = pylab.figure(figsize=(13.5, 7.))
		sp2 = fig.add_subplot(122)
		sp1 = fig.add_subplot(121)

		sp1.minorticks_on()
		sp2.minorticks_on()
		sp1.set_yscale('log')

		sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number  [dex$^{-1}$]')
		sp2.set_xlabel('log( M$^*$ / M$_{\odot}$ )')
		sp2.set_ylabel(r'$\alpha$', size=26)

		fig.subplots_adjust(left=0.09)

		sp1.grid()
		sp2.grid()

		sp1.axis([9.2, 11.85, 20, 5000])
		sp2.axis([10.5, 11.55, -1.55, -0.3])
		sp1.axis([9.2, 11.85, 5, 5000])
		sp2.axis([10.5, 11.55, -1.65, -0.2])

		offie = [-0.015, -0.005, 0.005, 0.015]
		for vi in range(len(voronoi_bins)):

			plot_fields = []
			plot_fields.append(mf_data('N200',    0.691,  0.027, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC1324',  0.755,  0.033, voronoi_bins[vi]))
			plot_fields.append(mf_data('RCS0224', 0.772,  0.027, voronoi_bins[vi]))
			plot_fields.append(mf_data('RXJ1716', 0.813,  0.021, voronoi_bins[vi]))
			plot_fields.append(mf_data('N5281',   0.818,  0.029, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC1604',  0.910,  0.029, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC0910',  1.110,  0.035, voronoi_bins[vi]))
			plot_fields.append(mf_data('SC0849',  1.261,  0.029, voronoi_bins[vi]))

			nqu = pylab.zeros(len(lmassbars))
			nlos = pylab.zeros((len(plot_fields), len(lmassbars)))
			nhis = pylab.zeros((len(plot_fields), len(lmassbars)))

			for i in range(len(plot_fields)):
				f = plot_fields[i]
				nqu += f.data.Nqu[2:]
				nlos[i] += f.data.Nlo_qu[2:]
				nhis[i] += f.data.Nhi_qu[2:]

			nlos = (nlos**2).sum(axis=0)**0.5
			nhis = (nhis**2).sum(axis=0)**0.5

			color_i = cmap(vi/(len(voronoi_bins)-1.))

			sp1.errorbar(lmassbars+offie[vi], nqu/dm, xerr=dm/2., yerr=[nlos/dm, nhis/dm],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, ecolor=color_i, label=voronoi_labels[vi])



			###  fitting single-Schecther function
			xdata, ydata = lmassbars, nqu/dm
			ydata_errs = (nlos + nhis) / 2. / dm
			fit0, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1, 10.5, pylab.median(ydata)], sigma=ydata_errs)
			sp2.plot(fit0[1], fit0[0], ls='', marker='o', mew=3, ms=9, mec='k', mfc=color_i, zorder=10)

			xax = pylab.linspace(xdata.min(), xdata.max(), 100)
			yax = schechter_mf(xax, *fit0)
			sp1.plot(xax, yax, lw=1.5, color=color_i)



			###  running MC simulations
			print '\nrunning %i MC simulations for %s' % (nsims, voronoi_labels[vi])
			xdata, ydata = lmassbars, nqu
			ydata_errs = (nlos + nhis) / 2.

			mc_nqu = pylab.poisson(lam=ydata, size=(nsims, len(ydata)))
			bestfits = []

			for mci in range(len(mc_nqu)):

				mypy.progress_bar(mci, nsims)

				ydata_mci = mc_nqu[mci]
				try:
					fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata_mci, p0=[-1, 10.5, pylab.median(ydata)], sigma=ydata_errs)
					bestfits.append(fit)
				except:
					pass

			bestfits = pylab.array(bestfits)
			print ''

			### plotting results from MC sims
			#sp2.plot(bestfits[:,1], bestfits[:,0], 'ko', ms=1, mec=color_i, zorder=1)

			gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
			xlo, xhi = 10.5, 11.8
			ylo, yhi = -1.7, -0.1
			nbins = 100
			xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
			hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

			levels = hist2d.max() * gauss_levels[0:2]

			contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, colors=[color_i, color_i], linewidths=[3, 1])

		leg1 = sp1.legend(loc=3, numpoints=1, fontsize=14)

	outname = '../output/MF_voronoiBins_qu.png'
	fig.savefig(outname)
	print 'wrote to: %s\n' % outname
	pylab.close()






















########################################################
###  Plot Field MFs by summing up MFs in each Voronoi
###  slice below an overdensity threshold.
########################################################

plot_fields = []
plot_fields.append(mf_data('N200',    0.691,  0.027))
plot_fields.append(mf_data('SC1324',  0.755,  0.033))
plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields.append(mf_data('N5281',   0.818,  0.029))
plot_fields.append(mf_data('SC1604',  0.910,  0.029))
plot_fields.append(mf_data('SC0910',  1.110,  0.035))
plot_fields.append(mf_data('SC0849',  1.261,  0.029))

###  Field Mass Functions - ZFOURGE 2014
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]

zfourge_params_tot = mypy.readcat('%s/table2_TOT.dat' % zfourge_dir)

if False:

	outname = '../output/massFunctions_from_summed_voronoi_slices.pdf'
	pdf = PdfPages(outname)

	fig = pylab.figure(figsize=(17., 13.5))
	sp3 = fig.add_subplot(2, 3, 3)
	sp6 = fig.add_subplot(2, 3, 6)
	sp2 = fig.add_subplot(2, 3, 2)
	sp5 = fig.add_subplot(2, 3, 5)
	sp1 = fig.add_subplot(2, 3, 1)
	sp4 = fig.add_subplot(2, 3, 4)
	sps = [sp1, sp2, sp3, sp4, sp5, sp6]

	fig.subplots_adjust(wspace=0, hspace=0, top=0.94, left=0.09, bottom=0.09)

	lmassbins = plot_fields[0].mf_voronoi_slices.lmassbins
	lmassbars = plot_fields[0].mf_voronoi_slices.lmassbars
	dm = lmassbins[1] - lmassbins[0]

	zbins = pylab.array([[0.20, 0.50],
                         [0.50, 0.75],
                         [0.75, 1.00],
                         [1.00, 1.25],
                         [1.25, 1.50],
                         [1.50, 2.00]])
	zbars = (zbins[:,1] + zbins[:,0]) / 2.
	dz = (zbins[:,1] - zbins[:,0]) / 2.


	###  master 3d array of number counts. Axes are zbin, massbin, LSSfield
	master_number_counts = pylab.zeros((len(plot_fields), len(zbins), len(lmassbars)))
	master_volumes = pylab.zeros((len(plot_fields), len(zbins), len(lmassbars)))


	for fi in range(len(plot_fields)):

		mypy.progress_bar(fi, len(plot_fields))

		f = plot_fields[fi]
		title = '%s :   z$_{LSS}$ = %.3f' % (f.name, f.zclust)
		sp2.set_title(title)

		for zi in range(len(zbars)):

			zlo, zhi = zbins[zi]

			sp = sps[zi]
			sp.grid()
			sp.minorticks_on()
			axis_lims = [7.7, 11.8, -5.4, -0.9]
			sp.axis(axis_lims)
			sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
			sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
			sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)

			###  plotting red bg for zbins that prob include LSS galaxies
			zphot_lo = f.zclust - 1. * f.sigmaz * (1 + f.zclust)
			zphot_hi = f.zclust + 1. * f.sigmaz * (1 + f.zclust)
			if (zlo < zphot_lo < zhi) or (zlo < zphot_hi < zhi):
				sp.fill_between(axis_lims[0:2], axis_lims[2], axis_lims[3], color='#ffe6e6', zorder=0)


			###  plotting zfourge first
			i_mf = 3*zi + 1
			inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
			sp.fill_between(zfourge_massfunc_tot[inds, 0],
				            zfourge_massfunc_tot[inds, i_mf] + zfourge_massfunc_tot[inds, i_mf+1], 
				            zfourge_massfunc_tot[inds, i_mf] - zfourge_massfunc_tot[inds, i_mf+2], 
				            color='#bfbfbf', zorder=1)
			sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')
	


			###  ORELSE next: running through Voronoi slices
			zlo_vbins, zhi_vbins = 9, 0
			for vi in range(0, len(f.mf_voronoi_slices.zbars), 2):
				###  skipping bins with possible LSS contamination
				if (zlo < zphot_lo < zhi) or (zlo < zphot_hi < zhi): continue

				###  checking if Voronoi slice is within the zbin
				if f.mf_voronoi_slices.zlos[vi] >= zlo and f.mf_voronoi_slices.zhis[vi] <= zhi:
					if f.mf_voronoi_slices.zlos[vi] < zlo_vbins: zlo_vbins = f.mf_voronoi_slices.zlos[vi]
					if f.mf_voronoi_slices.zhis[vi] > zhi_vbins: zhi_vbins = f.mf_voronoi_slices.zhis[vi]
					###  iterating over massbins, checking for mass completeness
					for mi in range(len(lmassbars)):
						lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[vi], f.masslimits.z, f.masslimits.masslim_empirical)
						if lmassbins[mi] > lmasslimit:
							master_number_counts[fi][zi][mi] += f.mf_voronoi_slices.ngals_field[vi][mi]
							master_volumes[fi][zi][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[vi]
							#master_number_counts[fi][zi][mi] += f.mf_voronoi_slices.ngals_00v05[vi][mi]
							#master_volumes[fi][zi][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[vi]

			###  plotting zwindow text
			if master_number_counts[fi][zi].max() > 0:
				sp.text(0.5, 0.85, '%.2f < z < %.2f\n%.2f < z < %.2f' % (zlo, zhi, zlo_vbins, zhi_vbins), transform=sp.transAxes, fontweight='bold')
			else:
				sp.text(0.5, 0.85, '%.2f < z < %.2f\n' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')

			elo_poisson, ehi_poisson = mypy.massfunc.CI(master_number_counts[fi][zi], sigma=1)
			
			phi = master_number_counts[fi][zi] / master_volumes[fi][zi] / dm
			elo_phi = phi * (elo_poisson / master_number_counts[fi][zi])
			ehi_phi = phi * (ehi_poisson / master_number_counts[fi][zi])

			lphi = pylab.log10(phi)
			elo_lphi = pylab.log10(phi / (phi - elo_phi))
			ehi_lphi = pylab.log10((phi + ehi_phi) / phi)

			inds_plot = pylab.find(master_number_counts[fi][zi] > 0)
			sp.errorbar(lmassbars[inds_plot], lphi[inds_plot], xerr=dm/2., yerr=[elo_lphi[inds_plot], ehi_lphi[inds_plot]],
				        ls='', marker='o', ms=8, mfc='r', ecolor='k')


			if zi == 3:
				sp.legend(loc=3, fontsize=18, title=title)


		pdf.savefig()
		for sp in sps:
			sp.clear()

	print '\n'







	'''
	fig = pylab.figure(figsize=(17., 13.5))
	sp3 = fig.add_subplot(2, 3, 3)
	sp6 = fig.add_subplot(2, 3, 6)
	sp2 = fig.add_subplot(2, 3, 2)
	sp5 = fig.add_subplot(2, 3, 5)
	sp1 = fig.add_subplot(2, 3, 1)
	sp4 = fig.add_subplot(2, 3, 4)
	sps = [sp1, sp2, sp3, sp4, sp5, sp6]

	fig.subplots_adjust(wspace=0, hspace=0, top=0.94, left=0.09, bottom=0.09)
	'''

	###  plotting volume-average of all fields
	single_schechter_models = []
	double_schechter_models = []
	dschechter_fits = pylab.zeros((len(zbins), 5))
	dschechter_errs = pylab.zeros((len(zbins), 5))
	summed_number_counts = pylab.sum(master_number_counts, axis=0)
	summed_volumes = pylab.sum(master_volumes, axis=0)
	for zi in range(len(f.mf_field.zbars)-1):

		zlo, zhi = f.mf_field.zlos[zi], f.mf_field.zhis[zi]
		title = 'All fields minus LSS'

		ngal = summed_number_counts[zi]
		vol = summed_volumes[zi]
		elo_poisson, ehi_poisson = mypy.massfunc.CI(ngal, sigma=1)
		inds_plot = pylab.find((ngal > 0))# & (lmassbars > 9.49))

		phi = ngal * 1. / vol / dm
		elo_phi = phi * (elo_poisson / ngal)
		ehi_phi = phi * (ehi_poisson / ngal)

		lphi = pylab.log10(phi)
		elo_lphi = pylab.log10(phi / (phi - elo_phi))
		ehi_lphi = pylab.log10((phi + ehi_phi) / phi)

		###  scaling errors by volume per massbin
		#elo_lphi += pylab.log10((vol.max() / vol)**0.5)
		#ehi_lphi += pylab.log10((vol.max() / vol)**0.5)



		sp = sps[zi]
		sp.grid()
		sp.minorticks_on()
		axis_lims = [7.7, 11.8, -5.4, -0.9]
		sp.axis(axis_lims)
		sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
		sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)


		sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


		###  plotting zfourge first
		i_mf = 3*zi + 1
		inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
		sp.fill_between(zfourge_massfunc_tot[inds, 0],
			            zfourge_massfunc_tot[inds, i_mf] + zfourge_massfunc_tot[inds, i_mf+1], 
			            zfourge_massfunc_tot[inds, i_mf] - zfourge_massfunc_tot[inds, i_mf+2], 
		    	        color='#bfbfbf', zorder=1)
		sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')

		###  plotting ZFOURGE 2x schechter
		p = [zfourge_params_tot.Mstar[zi], zfourge_params_tot.alpha1[zi], zfourge_params_tot.alpha2[zi], 10**zfourge_params_tot.phistar1[zi], 10**zfourge_params_tot.phistar2[zi]]
		xmod = pylab.linspace(7.95, 11.8, 1000)
		sp.plot(xmod, pylab.log10(dschechter(xmod, *p)), color='gray')


		###  plotting orelse next
		sp.errorbar(lmassbars[inds_plot], lphi[inds_plot], xerr=dm/2., yerr=[elo_lphi[inds_plot], ehi_lphi[inds_plot]],
			        ls='', marker='o', ms=7, mew=1.5, elinewidth=1.5, mfc='k', ecolor='k', label='ORELSE today', zorder=2)

		if zi == 3:
			sp.legend(loc=3, fontsize=18, title=title)


		###  fitting double Schechter
		model_single = single_schechter_model()
		model_double = double_schechter_model()
		try:
			x = lmassbars[inds_plot]
			y = phi[inds_plot]
			yerr = ((elo_lphi[inds_plot] + ehi_lphi[inds_plot]) / 2.)
			fit_single = fit_algorithm(model_single, x, y, weights=yerr)
			fit_double = fit_algorithm(model_double, x, y, weights=yerr)
			single_schechter_models.append(fit_single)
			double_schechter_models.append(fit_double)

			xax = pylab.linspace(x.min(), x.max(), 100)
			yax = fit_double(xax)
			sp.plot(xax, pylab.log10(yax), color='r', lw=1.5)
		except: pass






	pdf.savefig()
	pdf.close()
	pylab.close()
	print 'wrote to: %s' % outname






















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


	data_table_tot = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
	data_table_sf  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
	data_table_qu  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99



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

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

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

	sp1.grid()
	sp2.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
	sp2.set_ylabel('Ratio')

	sp1.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])
	sp2.axis([8.4, 11.9, 0.5, 1000])

	fig.subplots_adjust(left=0.09)

	ydatas = []
	fits_tot, covs_tot = [], []
	line_fits_tot, line_covs_tot = [], []
	mlim_hi = 11.55



	###  plotting ZFOURGE
	lmassbars_zf = zfourge_massfunc_tot[:,0]

	area_zf = 316.   # arcmin**2
	volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
	volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
	volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
	volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

	lmf_050z075 = zfourge_massfunc_tot[:,4]
	lmf_075z100 = zfourge_massfunc_tot[:,7]
	lmf_100z125 = zfourge_massfunc_tot[:,10]

	lmf_err_050z075 = (zfourge_massfunc_tot[:,5] + zfourge_massfunc_tot[:,6]) / 2.
	lmf_err_075z100 = (zfourge_massfunc_tot[:,8] + zfourge_massfunc_tot[:,9]) / 2.
	lmf_err_100z125 = (zfourge_massfunc_tot[:,11] + zfourge_massfunc_tot[:,12]) / 2.

	lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
	lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
	lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

	inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))


	#sp1.fill_between(lmassbars_zf[inds_zf], 10**(lmf_050z125[inds_zf]-lmf_err_050z125[inds_zf]), 10**(lmf_050z125[inds_zf]+lmf_err_050z125[inds_zf]), color='#cccccc', label='ZFOURGE1')



	###  deciding which overdensity bins to plot
	o_range = range(1, len(overdens_bars))
	offie = pylab.linspace(-len(o_range)*0.01/2, len(o_range)*0.01/2, len(o_range))

	chi2s_single = []
	chi2s_double = []

	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)
		sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  adding to data_table
		data_table_tot[:,0] = lmassbars
		data_table_tot[inds,3*vi_counter+1] = y[inds]
		data_table_tot[inds,3*vi_counter+2] = elo[inds]
		data_table_tot[inds,3*vi_counter+3] = ehi[inds]



		###  plotting ratio
		sp2.errorbar(x[inds] + offie[vi_counter], y[inds] / ydatas[0][inds], xerr=dm/2., yerr=[elo[inds] / ydatas[0][inds], ehi[inds] / ydatas[0][inds]],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i)



		###  fitting (d)schechter functions and plotting
		try:

			xdata, ydata = lmassbars[inds], (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0))[vi][inds]
			ydata_errs = 10**(elo+ehi)[inds]/2.
			ydata_errs = (elo+ehi)[inds]/2.


			###  fitting log-linear relation to ratio
			fit, cov = optimize.curve_fit(line, xdata, pylab.log10(y[inds] / ydatas[0][inds]), p0=[-8, 0.5], sigma=(elo+ehi)[inds]/2.)
			line_fits_tot.append(fit)
			line_covs_tot.append(cov.diagonal()**0.5)
			if vi > 0:
				#a = sp2.plot(xmod + offie[vi_counter], 10**line(xmod, *fit), color=color_i, lw=1)
				sp2.axvline(0, color=color_i, lw=2, label='slope = %.2f' % fit[0])


			###  fitting schechter funtion to SMF
			fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1., 11., 10**-3], sigma=ydata_errs, maxfev=10**5)
			#ymod = schechter_mf(xmod, *fit)

			chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
			chi2s_single.append(chi2.sum())


			###  fitting dschechter funtion to SMF
			p0 = [11, -1., -1., 10**-3, 10**-4]
			p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
			fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=p0, sigma=ydata_errs, maxfev=10**5)
			ymod = dschechter(xmod, *fit)
			a = sp1.plot(xmod + offie[vi_counter], ymod / dm, color=color_i, lw=1)
			fits_tot.append(fit)
			covs_tot.append(cov.diagonal()**0.5)

			chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
			chi2s_double.append(chi2.sum())


		except:
			pass



	fits_tot = pylab.array(fits_tot)
	covs_tot = pylab.array(covs_tot)

	chi2s_single = pylab.array(chi2s_single)
	chi2s_double = pylab.array(chi2s_double)

	###  print TeX table of dSchechter params
	if False:
		fits_tot[:,3] *= 10**3
		fits_tot[:,4] *= 10**3
		covs_tot[:,3] *= 10**3
		covs_tot[:,4] *= 10**3
		order = pylab.array([0, 1, 3, 2, 4])
		for vi_counter, vi in enumerate(o_range):
			s = '%s & ' % voronoi_labels[vi]
			for f, c in zip(fits_tot[vi_counter][order], covs_tot[vi_counter][order]):
				s += '$%5.2f\pm%4.2f$ & ' % (f, c)
			print s[:-2] + ' \\\\'


	sp1.axvline(0, color='#cccccc', lw=7, label='ZFOURGE')
	sp1.legend(loc=3, numpoints=1, fontsize=14)

	#t = sp2.text(0.05, 0.9, 'Total', fontweight='bold', fontsize=22, color='k', transform=sp2.transAxes)
	#sp2.legend(loc=2, title='Total', fontsize=16)




	###  writing data_table to ascii file
	outname = '../data/massFunctions_TOT_Tomczak+2017'
	for f in plot_fields: outname += '_' + f.name
	outer = open(outname + '.dat', 'w')

	outer.write('#    Galaxy stellar mass functions for all galaxies\n')
	outer.write('#\n')
	outer.write('#    Column  Units            Description\n')
	outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
	outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
	outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
	outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
	outer.write('#\n')
	outer.write('#    Overdensity bins\n')
	outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
	outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
	outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
	outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
	outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
	outer.write('#\n')
	outer.write('# lmass')
	for vi_counter, vi in enumerate(o_range):
		if len(voronoi_bins[vi]) == 6:
			outer.write('   phi%s' % (voronoi_bins[vi]))
			outer.write('   elo%s' % (voronoi_bins[vi]))
			outer.write('   ehi%s' % (voronoi_bins[vi]))
		if len(voronoi_bins[vi]) == 7:
			outer.write('  phi%s' % (voronoi_bins[vi]))
			outer.write('  elo%s' % (voronoi_bins[vi]))
			outer.write('  ehi%s' % (voronoi_bins[vi]))
	outer.write('\n')

	for row in data_table_tot:
		outer.write('%6.2f ' % row[0])
		for item in row[1:]:
			outer.write('  %.4e' % item)
		outer.write('\n')

	outer.close()











	##############################
	###  plotting for SF galaxies
	##############################

	fig = pylab.figure(figsize=(15.5, 12.7))
	sp2 = fig.add_subplot(222)
	sp1 = fig.add_subplot(221)
	sp4 = fig.add_subplot(224)
	sp3 = fig.add_subplot(223)

	#sp1.grid()
	#sp2.grid()
	#sp3.grid()
	#sp4.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()
	sp3.minorticks_on()
	sp4.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')
	sp3.set_yscale('log')
	sp4.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
	sp2.set_ylabel('Ratio')
	sp3.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
	sp4.set_ylabel('Ratio')

	sp1.axis([8.4, 11.9, 8*10**-6, 4*10**-1])
	sp2.axis([8.4, 11.9, 0.5, 2000])
	sp3.axis([8.4, 11.9, 8*10**-6, 4*10**-1])
	sp4.axis([8.4, 11.9, 0.5, 2000])

	fig.subplots_adjust(left=0.09, hspace=0, bottom=0.08)




	###  plotting ZFOURGE
	lmassbars_zf = zfourge_massfunc_sf[:,0]

	area_zf = 316.   # arcmin**2
	volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
	volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
	volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
	volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

	lmf_050z075 = zfourge_massfunc_sf[:,4]
	lmf_075z100 = zfourge_massfunc_sf[:,7]
	lmf_100z125 = zfourge_massfunc_sf[:,10]

	lmf_err_050z075 = (zfourge_massfunc_sf[:,5] + zfourge_massfunc_sf[:,6]) / 2.
	lmf_err_075z100 = (zfourge_massfunc_sf[:,8] + zfourge_massfunc_sf[:,9]) / 2.
	lmf_err_100z125 = (zfourge_massfunc_sf[:,11] + zfourge_massfunc_sf[:,12]) / 2.

	lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
	lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
	lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

	inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))


	sp1.fill_between(lmassbars_zf[inds_zf], 10**(lmf_050z125[inds_zf]-lmf_err_050z125[inds_zf]), 10**(lmf_050z125[inds_zf]+lmf_err_050z125[inds_zf]), color='#cccccc', label='ZFOURGE')
	sp1.axvline(0, color='#cccccc', lw=7, label='ZFOURGE')





	chi2s_single = []
	chi2s_double = []


	ydatas = []
	fits_sf, covs_sf = [], []
	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0) / dm)[vi], master_number_counts_sf.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)

		n_error = n * 1.
		y_error = y * 1.
		ninds = pylab.find(n_error == 0)
		n_error[ninds] = 1.
		y_error[ninds] = 1. / master_volumes_sf.sum(axis=0)[vi][ninds] / dm

		elo = y_error * (ylo / n_error)
		ehi = y_error * (yhi / n_error)
		sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  adding to data_table
		data_table_sf[:,0] = lmassbars
		data_table_sf[inds,3*vi_counter+1] = y[inds]
		data_table_sf[inds,3*vi_counter+2] = elo[inds]
		data_table_sf[inds,3*vi_counter+3] = ehi[inds]





		###  plotting ratio
		sp2.errorbar(x[inds] + offie[vi_counter], y[inds] / ydatas[0][inds], xerr=dm/2., yerr=[elo[inds] / ydatas[0][inds], ehi[inds] / ydatas[0][inds]],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])



		###  fitting (d)schechter functions and plotting

		xdata, ydata = lmassbars[inds], (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0))[vi][inds]
		ydata_errs = (elo+ehi)[inds]/2.

		fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)
		ymod = schechter_mf(xmod, *fit)

		chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
		chi2s_single.append(chi2.sum())



		fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[11, -1., -1., 10**-3, 10**-4], sigma=ydata_errs, maxfev=10**4)
		#ymod = dschechter(xmod, *fit)

		chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
		chi2s_double.append(chi2.sum())

		
		a = sp1.plot(xmod, ymod / dm, color=color_i, lw=1)
		fits_sf.append(fit)
		covs_sf.append(cov.diagonal()**0.5)

	fits_sf = pylab.array(fits_sf)
	covs_sf = pylab.array(covs_sf)

	chi2s_single = pylab.array(chi2s_single)
	chi2s_double = pylab.array(chi2s_double)

	#sp1.legend(loc=3, numpoints=1, fontsize=16)

	t = sp2.text(0.05, 0.9, 'Star-Forming', fontweight='bold', fontsize=22, color='b', transform=sp2.transAxes)


	###  writing data_table to ascii file
	outname = '../data/massFunctions_SF_Tomczak+2017'
	for f in plot_fields: outname += '_' + f.name
	outer = open(outname + '.dat', 'w')

	outer.write('#    Galaxy stellar mass functions for star-forming galaxies\n')
	outer.write('#\n')
	outer.write('#    Column  Units            Description\n')
	outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
	outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
	outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
	outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
	outer.write('#\n')
	outer.write('#    Overdensity bins\n')
	outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
	outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
	outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
	outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
	outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
	outer.write('#\n')
	outer.write('# lmass')
	for vi_counter, vi in enumerate(o_range):
		if len(voronoi_bins[vi]) == 6:
			outer.write('   phi%s' % (voronoi_bins[vi]))
			outer.write('   elo%s' % (voronoi_bins[vi]))
			outer.write('   ehi%s' % (voronoi_bins[vi]))
		if len(voronoi_bins[vi]) == 7:
			outer.write('  phi%s' % (voronoi_bins[vi]))
			outer.write('  elo%s' % (voronoi_bins[vi]))
			outer.write('  ehi%s' % (voronoi_bins[vi]))
	outer.write('\n')

	for row in data_table_sf:
		outer.write('%6.2f ' % row[0])
		for item in row[1:]:
			outer.write('  %.4e' % item)
		outer.write('\n')

	outer.close()















	##############################
	###  plotting for QU galaxies
	##############################




	###  plotting ZFOURGE
	lmassbars_zf = zfourge_massfunc_qu[:,0]

	area_zf = 316.   # arcmin**2
	volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
	volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
	volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
	volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

	lmf_050z075 = zfourge_massfunc_qu[:,4]
	lmf_075z100 = zfourge_massfunc_qu[:,7]
	lmf_100z125 = zfourge_massfunc_qu[:,10]

	lmf_err_050z075 = (zfourge_massfunc_qu[:,5] + zfourge_massfunc_qu[:,6]) / 2.
	lmf_err_075z100 = (zfourge_massfunc_qu[:,8] + zfourge_massfunc_qu[:,9]) / 2.
	lmf_err_100z125 = (zfourge_massfunc_qu[:,11] + zfourge_massfunc_qu[:,12]) / 2.

	lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
	lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
	lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

	inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))


	sp3.fill_between(lmassbars_zf[inds_zf], 10**(lmf_050z125[inds_zf]-lmf_err_050z125[inds_zf]), 10**(lmf_050z125[inds_zf]+lmf_err_050z125[inds_zf]), color='#cccccc', label='ZFOURGE')
	sp3.axvline(0, color='#cccccc', lw=7, label='ZFOURGE')




	chi2s_single = []
	chi2s_double = []


	qu_inds = []   # indices where qu SMFs are mass-complete for each of the overdensity bins

	ydatas = []
	fits_qu, covs_qu = [], []
	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		#inds = pylab.find((master_volumes_qu.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		inds = pylab.find((master_number_counts_qu.sum(axis=0)[vi] > 0) & 
			              (master_volumes_qu.sum(axis=0)[vi] > 0) & 
			              (lmassbars < mlim_hi))
		qu_inds.append(inds)
		x, y, n = lmassbars, (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0) / dm)[vi], master_number_counts_qu.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)
		sp3.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  adding to data_table
		data_table_qu[:,0] = lmassbars
		data_table_qu[inds,3*vi_counter+1] = y[inds]
		data_table_qu[inds,3*vi_counter+2] = elo[inds]
		data_table_qu[inds,3*vi_counter+3] = ehi[inds]





		###  plotting ratio
		sp4.errorbar(x[inds] + offie[vi_counter], y[inds] / ydatas[0][inds], xerr=dm/2., yerr=[elo[inds] / ydatas[0][inds], ehi[inds] / ydatas[0][inds]],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])



		###  fitting (d)schechter functions and plotting

		#if vi_counter in [0, 1]:
		if 0:
			xdata, ydata = lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds]
			ydata_errs = (elo+ehi)[inds]/2.

			p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
			fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[10.83, -0.5, -1.55, 10**-2.81, 10**-4.40], sigma=pylab.log10(ydata_errs))
			ymod = dschechter(xmod, *fit)
			aa = sp3.plot(xmod, ymod / dm, color=color_i, lw=1)[0]

		else:
			fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=pylab.log10((elo+ehi)[inds]/2.), maxfev=10**6)
			fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)
			ymod = schechter_mf(xmod, *fit)
			a = sp3.plot(xmod, ymod / dm, color=color_i, lw=1)


		fits_qu.append(fit)
		covs_qu.append(cov.diagonal()**0.5)


		###  calculating chi2s
		xdata, ydata = lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds]
		ydata_errs = (elo+ehi)[inds]/2.

		fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)

		chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
		chi2s_single.append(chi2.sum())


		p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
		fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[10.83, -0.5, -1.55, 10**-2.81, 10**-4.40], sigma=pylab.log10(ydata_errs), maxfev=10**5)

		chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
		chi2s_double.append(chi2.sum())




	fits_qu = pylab.array(fits_qu)
	covs_qu = pylab.array(covs_qu)

	chi2s_single = pylab.array(chi2s_single)
	chi2s_double = pylab.array(chi2s_double)



	###  print TeX table of dSchechter params
	if True:
		print '\n'
		fits_sf[:,2] *= 10**3
		covs_sf[:,2] *= 10**3
		order = pylab.array([1, 0, 2])
		for vi in range(len(o_range)):
			s = '%s & ' % voronoi_labels[vi]
			for f, c in zip(fits_sf[vi][order], covs_sf[vi][order]):
				s += '$%5.2f\pm%4.2f$ & ' % (f, c)
			print s[:-2] + ' \\\\'

		print '\n'
		fits_qu[:,2] *= 10**3
		covs_qu[:,2] *= 10**3
		order = pylab.array([1, 0, 2])
		for vi in range(len(o_range)):
			s = '%s & ' % voronoi_labels[vi]
			for f, c in zip(fits_qu[vi][order], covs_qu[vi][order]):
				s += '$%5.2f\pm%4.2f$ & ' % (f, c)
			print s[:-2] + ' \\\\'




	sp1.legend(loc=3, numpoints=1, fontsize=14)

	t = sp4.text(0.05, 0.9, 'Quiescent', fontweight='bold', fontsize=22, color='r', transform=sp4.transAxes)


	###  writing data_table to ascii file
	outname = '../data/massFunctions_QU_Tomczak+2017'
	for f in plot_fields: outname += '_' + f.name
	outer = open(outname + '.dat', 'w')

	outer.write('#    Galaxy stellar mass functions for quiescent galaxies\n')
	outer.write('#\n')
	outer.write('#    Column  Units            Description\n')
	outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
	outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
	outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
	outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
	outer.write('#\n')
	outer.write('#    Overdensity bins\n')
	outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
	outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
	outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
	outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
	outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
	outer.write('#\n')
	outer.write('# lmass')
	for vi_counter, vi in enumerate(o_range):
		if len(voronoi_bins[vi]) == 6:
			outer.write('   phi%s' % (voronoi_bins[vi]))
			outer.write('   elo%s' % (voronoi_bins[vi]))
			outer.write('   ehi%s' % (voronoi_bins[vi]))
		if len(voronoi_bins[vi]) == 7:
			outer.write('  phi%s' % (voronoi_bins[vi]))
			outer.write('  elo%s' % (voronoi_bins[vi]))
			outer.write('  ehi%s' % (voronoi_bins[vi]))
	outer.write('\n')

	for row in data_table_qu:
		outer.write('%6.2f ' % row[0])
		for item in row[1:]:
			outer.write('  %.4e' % item)
		outer.write('\n')

	outer.close()

	pylab.close()
	pylab.close()























###  test




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

mlim_hi = 11.51

if True:





    lmassbins = plot_fields[0].mf_voronoi_slices.lmassbins
    lmassbars = plot_fields[0].mf_voronoi_slices.lmassbars
    dm = lmassbins[1] - lmassbins[0]

    overdens_bins = plot_fields[0].mf_voronoi_slices.overdens_bins
    overdens_bars = plot_fields[0].mf_voronoi_slices.overdens_bars
    d_overdens = overdens_bins[1] - overdens_bins[0]


    data_table_tot = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
    data_table_sf  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
    data_table_qu  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99



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

            lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
            lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

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





    ##############################
    ###  plotting for SF galaxies
    ##############################

    ###  deciding which overdensity bins to plot
    o_range = range(1, len(overdens_bars))
    offie = pylab.linspace(-len(o_range)*0.01/2, len(o_range)*0.01/2, len(o_range))




    fig = pylab.figure(figsize=(15.5, 12.7))
    sp1 = fig.add_subplot(111)

    #sp1.grid()
    #sp2.grid()
    #sp3.grid()
    #sp4.grid()
    sp1.minorticks_on()
    sp1.set_yscale('log')

    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')

    sp1.axis([8.4, 11.9, 8*10**-6, 4*10**-1])






    chi2s_single = []
    chi2s_double = []


    ydatas = []
    fits_sf, covs_sf = [], []
    for vi_counter, vi in enumerate(o_range):
    #for vi in [0]:

        color_i = cmap(vi_counter/(len(o_range)-1.))

        inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
        inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
        x, y, n = lmassbars, (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0) / dm)[vi], master_number_counts_sf.sum(axis=0)[vi]
        ydatas.append(y)
        ylo, yhi = mypy.massfunc.CI(n)

        n_error = n * 1.
        y_error = y * 1.
        ninds = pylab.find(n_error == 0)
        n_error[ninds] = 1.
        y_error[ninds] = 1. / master_volumes_sf.sum(axis=0)[vi][ninds] / dm

        elo = y_error * (ylo / n_error)
        ehi = y_error * (yhi / n_error)
        xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)





        ###  fitting (d)schechter functions and plotting

        xdata, ydata = lmassbars[inds], (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0))[vi][inds]
        ydata_errs = (elo+ehi)[inds]/2.

        fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)
        ymod = schechter_mf(xmod, *fit)

        chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
        chi2s_single.append(chi2.sum())



        fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[11, -1., -1., 10**-3, 10**-4], sigma=ydata_errs, maxfev=10**4)
        #ymod = dschechter(xmod, *fit)

        chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
        chi2s_double.append(chi2.sum())

        
        a = sp1.plot(xmod, ymod / dm, color=color_i, lw=1)
        fits_sf.append(fit)
        covs_sf.append(cov.diagonal()**0.5)

    fits_sf = pylab.array(fits_sf)
    covs_sf = pylab.array(covs_sf)

    chi2s_single = pylab.array(chi2s_single)
    chi2s_double = pylab.array(chi2s_double)

    #sp1.legend(loc=3, numpoints=1, fontsize=16)
















    ##############################
    ###  plotting for QU galaxies
    ##############################




    chi2s_single = []
    chi2s_double = []


    qu_inds = []   # indices where qu SMFs are mass-complete for each of the overdensity bins

    ydatas = []
    fits_qu, covs_qu = [], []
    for vi_counter, vi in enumerate(o_range):
    #for vi in [0]:

        color_i = cmap(vi_counter/(len(o_range)-1.))

        #inds = pylab.find((master_volumes_qu.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
        inds = pylab.find((master_number_counts_qu.sum(axis=0)[vi] > 0) & 
                          (master_volumes_qu.sum(axis=0)[vi] > 0) & 
                          (lmassbars < mlim_hi))
        qu_inds.append(inds)
        x, y, n = lmassbars, (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0) / dm)[vi], master_number_counts_qu.sum(axis=0)[vi]
        ydatas.append(y)
        ylo, yhi = mypy.massfunc.CI(n)
        elo = y * (ylo / n)
        ehi = y * (yhi / n)
        xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



        ###  adding to data_table
        data_table_qu[:,0] = lmassbars
        data_table_qu[inds,3*vi_counter+1] = y[inds]
        data_table_qu[inds,3*vi_counter+2] = elo[inds]
        data_table_qu[inds,3*vi_counter+3] = ehi[inds]






        ###  fitting (d)schechter functions and plotting

        #if vi_counter in [0, 1]:
        if 0:
            xdata, ydata = lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds]
            ydata_errs = (elo+ehi)[inds]/2.

            p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
            fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[10.83, -0.5, -1.55, 10**-2.81, 10**-4.40], sigma=pylab.log10(ydata_errs))
            ymod = dschechter(xmod, *fit)
            aa = sp1.plot(xmod, ymod / dm, color=color_i, lw=1)[0]

        else:
            fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=pylab.log10((elo+ehi)[inds]/2.), maxfev=10**6)
            fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)
            ymod = schechter_mf(xmod, *fit)
            a = sp1.plot(xmod, ymod / dm, color=color_i, lw=1)


        fits_qu.append(fit)
        covs_qu.append(cov.diagonal()**0.5)


        ###  calculating chi2s
        xdata, ydata = lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds]
        ydata_errs = (elo+ehi)[inds]/2.

        fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)

        chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
        chi2s_single.append(chi2.sum())


        p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
        fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[10.83, -0.5, -1.55, 10**-2.81, 10**-4.40], sigma=pylab.log10(ydata_errs), maxfev=10**5)

        chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
        chi2s_double.append(chi2.sum())




    fits_qu = pylab.array(fits_qu)
    covs_qu = pylab.array(covs_qu)

    chi2s_single = pylab.array(chi2s_single)
    chi2s_double = pylab.array(chi2s_double)



    sp1.legend(loc=3, numpoints=1, fontsize=14)

































##################################################
###  plotting SF/Qui fractions:  4 overdens bins
##################################################

if False:

	if False:
		fig = pylab.figure()
		sp = fig.add_subplot(111)
		sp.grid()

		f_qu_matrix = pylab.zeros((len(overdens_bars), len(lmassbars)))
		flo_qu_matrix = pylab.zeros((len(overdens_bars), len(lmassbars)))
		fhi_qu_matrix = pylab.zeros((len(overdens_bars), len(lmassbars)))

		arctan_fits = []
		arctan_covs = []

		offie = [-0.015, +0.005, -0.005, 0.015]
		for vi in range(len(overdens_bars)):

			sp.minorticks_on()

			sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
			sp.set_ylabel('Quiescent fraction')
			sp.axis([8.4, 11.8, -0.1, 1.1])

			color_i = cmap(vi/(len(voronoi_bins)-1.))

			n_sf = master_number_counts_sf.sum(axis=0)[vi]
			n_qu = master_number_counts_qu.sum(axis=0)[vi]
			f_qu = n_qu / (n_sf + n_qu)

			nlo_sf, nhi_sf = mypy.massfunc.CI(n_sf)
			nlo_qu, nhi_qu = mypy.massfunc.CI(n_qu)

			err_sf = (nlo_sf + nhi_sf) / 2.
			err_qu = (nlo_qu + nhi_qu) / 2.

			flo_qu = pylab.sqrt(err_qu**2 * (1. / (n_sf + n_qu) - n_qu / (n_sf + n_qu)**2)**2 + err_sf**2 * (n_qu / (n_sf + n_qu)**2)**2)
			fhi_qu = pylab.sqrt(err_qu**2 * (1. / (n_sf + n_qu) - n_qu / (n_sf + n_qu)**2)**2 + err_sf**2 * (n_qu / (n_sf + n_qu)**2)**2)

			###  cropping errors in fraction at 0% and 100%
			finds = pylab.find(f_qu - flo_qu < 0.)
			flo_qu[finds] = f_qu[finds]
			finds = pylab.find(f_qu + fhi_qu > 1.)
			fhi_qu[finds] = 1. - f_qu[finds]


			inds = pylab.find(((n_sf + n_qu) > 0) & (lmassbars < 11.55))

			sp.errorbar(lmassbars[qu_inds[vi]] + offie[vi], f_qu[qu_inds[vi]], xerr=dm/2., yerr=[flo_qu[qu_inds[vi]], fhi_qu[qu_inds[vi]]],
					    ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])


			###  adding to matix
			f_qu_matrix[vi] += f_qu
			flo_qu_matrix[vi] += flo_qu
			fhi_qu_matrix[vi] += fhi_qu



			###  fitting arctangent
			fit, cov = optimize.curve_fit(fqu_arctan, lmassbars[qu_inds[vi]], f_qu_matrix[vi][qu_inds[vi]], p0=[10., 1.], 
				                          sigma=(flo_qu_matrix[vi][qu_inds[vi]] + fhi_qu_matrix[vi][qu_inds[vi]])/2., maxfev=10**5)

			arctan_fits.append(fit)
			arctan_covs.append(cov.diagonal()**0.5)

			xmod = pylab.linspace(9.1, 11.6, 1000)
			ymod = fqu_arctan(xmod, *fit)
			sp.plot(xmod, ymod, color=color_i, lw=1)




		arctan_fits = pylab.array(arctan_fits)
		arctan_covs = pylab.array(arctan_covs)



		a = sp.plot(vdB13.lmass, vdB13.nqu / vdB13.ntot, color='r', alpha=0.5, lw=3, label='GCLASS')[0]
		sp.legend(loc=4, numpoints=1, fontsize=14)






	#################################################
	###  plotting SF/Qui fractions:  2 overdens bins
	#################################################

	if False:
		fig = pylab.figure()
		sp = fig.add_subplot(111)
		sp.grid()

		f_qu_matrix = pylab.zeros((2, len(lmassbars)))
		flo_qu_matrix = pylab.zeros((2, len(lmassbars)))
		fhi_qu_matrix = pylab.zeros((2, len(lmassbars)))

		arctan_fits2 = []
		arctan_covs2 = []

		offie = [+0.005, -0.005]
		for vi in [0, 2]:

			sp.minorticks_on()

			sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
			sp.set_ylabel('Quiescent fraction')
			sp.axis([8.4, 11.8, -0.1, 1.1])

			color_i = pylab.cm.cool((vi*2+1)/6.)

			n_sf = master_number_counts_sf.sum(axis=0)[vi] + master_number_counts_sf.sum(axis=0)[vi+1]
			n_qu = master_number_counts_qu.sum(axis=0)[vi] + master_number_counts_qu.sum(axis=0)[vi+1]
			f_qu = n_qu / (n_sf + n_qu)

			nlo_sf, nhi_sf = mypy.massfunc.CI(n_sf)
			nlo_qu, nhi_qu = mypy.massfunc.CI(n_qu)

			err_sf = (nlo_sf + nhi_sf) / 2.
			err_qu = (nlo_qu + nhi_qu) / 2.

			flo_qu = pylab.sqrt(err_qu**2 * (1. / (n_sf + n_qu) - n_qu / (n_sf + n_qu)**2)**2 + err_sf**2 * (n_qu / (n_sf + n_qu)**2)**2)
			fhi_qu = pylab.sqrt(err_qu**2 * (1. / (n_sf + n_qu) - n_qu / (n_sf + n_qu)**2)**2 + err_sf**2 * (n_qu / (n_sf + n_qu)**2)**2)

			###  cropping errors in fraction at 0% and 100%
			finds = pylab.find(f_qu - flo_qu < 0.)
			flo_qu[finds] = f_qu[finds]
			finds = pylab.find(f_qu + fhi_qu > 1.)
			fhi_qu[finds] = 1. - f_qu[finds]


			inds = pylab.find(((n_sf + n_qu) > 0) & (lmassbars < 11.55))

			sp.errorbar(lmassbars[inds] + offie[vi/2], f_qu[inds], xerr=dm/2., yerr=[flo_qu[inds], fhi_qu[inds]],
					    ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels_merged[vi/2])


			###  adding to matix
			f_qu_matrix[vi/2] += f_qu
			flo_qu_matrix[vi/2] += flo_qu
			fhi_qu_matrix[vi/2] += fhi_qu


			###  fitting arctangent
			fit, cov = optimize.curve_fit(fqu_arctan, lmassbars[qu_inds[vi+1]], f_qu_matrix[vi/2][qu_inds[vi+1]], p0=[10., 1.], 
				                          sigma=(flo_qu_matrix[vi/2][qu_inds[vi+1]] + fhi_qu_matrix[vi/2][qu_inds[vi+1]])/2.)

			arctan_fits2.append(fit)
			arctan_covs2.append(cov.diagonal()**0.5)

			xmod = pylab.linspace(9.1, 11.6, 1000)
			ymod = fqu_arctan(xmod, *fit)
			sp.plot(xmod, ymod, color=color_i, lw=1)




		arctan_fits2 = pylab.array(arctan_fits2)
		arctan_covs2 = pylab.array(arctan_covs2)

		a = sp.plot(vdB13.lmass, vdB13.nqu / vdB13.ntot, color='r', alpha=0.5, lw=3, label='GCLASS')[0]
		sp.legend(loc=4, numpoints=1, fontsize=14)






	#################################
	###  plotting SF/Qui fractions
	###  4 overdens bins
	###  multipanel
	#################################

	fig = pylab.figure(figsize=(14., 7.75))
	sp2 = pylab.subplot2grid((6, 10), (0, 3), rowspan=3, colspan=3)
	sp1 = pylab.subplot2grid((6, 10), (0, 0), rowspan=3, colspan=3)
	sp4 = pylab.subplot2grid((6, 10), (3, 3), rowspan=3, colspan=3)
	sp3 = pylab.subplot2grid((6, 10), (3, 0), rowspan=3, colspan=3)
	sps = [sp1, sp2, sp3, sp4]

	sp5 = pylab.subplot2grid((6, 10), (2, 7), rowspan=3, colspan=3)

	sp5.minorticks_on()
	sp5.set_xlabel(r'$m_0$')
	sp5.set_ylabel(r'$\beta$', rotation=1)

	#sp5.set_title(r'$f_{QU} = \pi^{-1} \times arctan( \beta (m - m_0) ) + 0.5$')
	sp5.axis([9.3, 11.2, 0, 2.5])

	for i in range(len(sps)):
		sp = sps[i]
		sp.grid()
		sp.minorticks_on()
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp.set_ylabel('Quiescent fraction')
		sp.axis([8.85, 11.9, -0.1, 1.1])
		sp.text(0.05, 0.9, voronoi_labels[i], transform=sp.transAxes, fontsize=14, fontweight='bold')

	fig.subplots_adjust(wspace=0, hspace=0, left=0.08)



	f_qu_matrix = pylab.zeros((len(overdens_bars), len(lmassbars)))
	flo_qu_matrix = pylab.zeros((len(overdens_bars), len(lmassbars)))
	fhi_qu_matrix = pylab.zeros((len(overdens_bars), len(lmassbars)))

	fqu_fits = []
	fqu_covs = []

	offie = [-0.015, +0.005, -0.005, 0.015]
	for vi in range(len(overdens_bars)):

		sp = sps[vi]

		color_i = cmap(vi/(len(voronoi_bins)-1.))

		n_sf = master_number_counts_sf.sum(axis=0)[vi]
		n_qu = master_number_counts_qu.sum(axis=0)[vi]
		nd_sf = master_number_counts_sf.sum(axis=0)[vi] / master_volumes_sf.sum(axis=0)[vi]
		nd_qu = master_number_counts_qu.sum(axis=0)[vi] / master_volumes_qu.sum(axis=0)[vi]

		f_qu = nd_qu / (nd_sf + nd_qu)

		nlo_sf, nhi_sf = mypy.massfunc.CI(n_sf)
		nlo_qu, nhi_qu = mypy.massfunc.CI(n_qu)

		err_sf = (nlo_sf + nhi_sf) / 2.
		err_qu = (nlo_qu + nhi_qu) / 2.

		flo_qu = pylab.sqrt(err_qu**2 * (1. / (n_sf + n_qu) - n_qu / (n_sf + n_qu)**2)**2 + err_sf**2 * (n_qu / (n_sf + n_qu)**2)**2)
		fhi_qu = pylab.sqrt(err_qu**2 * (1. / (n_sf + n_qu) - n_qu / (n_sf + n_qu)**2)**2 + err_sf**2 * (n_qu / (n_sf + n_qu)**2)**2)

		###  cropping errors in fraction at 0% and 100%
		finds = pylab.find(f_qu - flo_qu < 0.)
		flo_qu[finds] = f_qu[finds]
		finds = pylab.find(f_qu + fhi_qu > 1.)
		fhi_qu[finds] = 1. - f_qu[finds]


		inds = pylab.find(((n_sf + n_qu) > 0) & (lmassbars < 11.55))

		sp.errorbar(lmassbars[qu_inds[vi]] + offie[vi], f_qu[qu_inds[vi]], xerr=dm/2., yerr=[flo_qu[qu_inds[vi]], fhi_qu[qu_inds[vi]]],
				    ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])


		###  adding to matix
		f_qu_matrix[vi] += f_qu
		flo_qu_matrix[vi] += flo_qu
		fhi_qu_matrix[vi] += fhi_qu



		###  fitting function
		#fit, cov = optimize.curve_fit(fqu_arctan, lmassbars[qu_inds[vi]], f_qu_matrix[vi][qu_inds[vi]], p0=[10., 1.], 
		#	                          sigma=(flo_qu_matrix[vi][qu_inds[vi]] + fhi_qu_matrix[vi][qu_inds[vi]])/2.)
		fit, cov = optimize.curve_fit(fqu_exp, lmassbars[qu_inds[vi]], f_qu_matrix[vi][qu_inds[vi]], p0=[10**11., 1.], 
			                          sigma=(flo_qu_matrix[vi][qu_inds[vi]] + fhi_qu_matrix[vi][qu_inds[vi]])/2.)

		fqu_fits.append(fit)
		fqu_covs.append(cov.diagonal()**0.5)

		###  plotting best-fit
		xmod = pylab.linspace(lmassbars[qu_inds[vi]][0]-dm/2., lmassbars[qu_inds[vi]][-1]+dm/2., 1000)
		ymod = fqu_exp(xmod, *fit)
		#sp.plot(xmod, ymod, color=color_i, lw=1.5)

		###  plotting fit parameters
		sp5.errorbar(fit[0], fit[1], xerr=cov[0][0]**0.5, yerr=cov[1][1]**0.5, 
				     ls='', mew=1.5, marker='s', ms=10, mfc=color_i, mec='k', elinewidth=1.5, ecolor='k', label=voronoi_labels[vi])


		###  plotting Peng+2010 prediction
		p10 = fqu_peng2010_z1(xmod, overdens_bars[vi])
		sp.plot(xmod, p10, color='gray', ls='-', lw=1.5)


	fqu_fits = pylab.array(fqu_fits)
	fqu_covs = pylab.array(fqu_covs)























##############################################
###  Plot CI contours for Schecther params
###  of the SF/QU mass functions
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
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]

###  GCLASS
vdB13 = mypy.readcat('../data/table3_vanderBurg+2013.dat')

master_bestfits_sf = [[] for i in range(5)]
master_bestfits_qu = [[] for i in range(5)]

if True:

	lmassbins = plot_fields[0].mf_voronoi_slices.lmassbins
	lmassbars = plot_fields[0].mf_voronoi_slices.lmassbars
	dm = lmassbins[1] - lmassbins[0]

	overdens_bins = plot_fields[0].mf_voronoi_slices.overdens_bins
	overdens_bars = plot_fields[0].mf_voronoi_slices.overdens_bars
	d_overdens = overdens_bins[1] - overdens_bins[0]


	###  master 3d array of number counts. Axes are LSSfield, overdensity bin, massbin, 
	master_number_counts_sf = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_number_counts_qu = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))

	master_volumes_sf = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))
	master_volumes_qu = pylab.zeros((len(plot_fields), len(overdens_bars), len(lmassbars)))

	for fi in range(len(plot_fields)):
		f =  plot_fields[fi]

		###  have to skip by twos because of overlappipng voronoi zslices
		for zi in range(0, len(f.mf_voronoi_slices.zbars), 2):

			###  skipping if zslice is behind cluster
			#if f.zclust < f.mf_voronoi_slices.zlos[zi] - 0.2: continue

			###  skipping redshifts of choice
			#if f.mf_voronoi_slices.zlos[zi] < 0.85: continue

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

			###  iterating over massbins, checking for mass completeness
			for mi in range(len(lmassbars)):

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

	#sp1.grid()
	#sp2.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()

	sp1.set_xlabel('log( $M^*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( $M^*$ / M$_{\odot}$ )')
	sp1.set_ylabel(r'$\alpha$', rotation=0, fontsize=26)
	sp2.set_ylabel(r'$\alpha$', rotation=0, fontsize=26)

	#sp1.axis([10.4, 11.9, -2, 0])
	#sp2.axis([10.4, 11.9, -2, 0])
	sp1.axis([10.35, 11.55, -1.75, 0.25])
	sp2.axis([10.35, 11.55, -1.75, 0.25])

	fig.subplots_adjust(left=0.09)



	ydatas = []
	mlim_hi = 11.55

	nsims0 = [4000, 5000, 35000, 50000, 60000]
	
	nres0 = 100   # resolution for plotting 2d contours
	rfactors = [1., 1., 1., 1., 0.5]

	###  deciding which overdensity bins to plot
	o_range = range(1, len(overdens_bars))
	offie = pylab.linspace(-len(o_range)*0.01/2, len(o_range)*0.01/2, len(o_range))

	for vi_counter, vi in enumerate(o_range):

		color_i = cmap(vi_counter/(len(o_range)-1.))
		nsims = nsims0[vi_counter]

		print '\n\n'
		print 'Star-forming galaxies'
		print '%s' % (voronoi_labels[vi])
		print 'Nsims = %i' % nsims

		inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi) & (lmassbars > 9.9))

		x, y, n = lmassbars, (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0) / dm)[vi], master_number_counts_sf.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)

		n_error = n * 1.
		y_error = y * 1.
		ninds = pylab.find(n_error == 0)
		n_error[ninds] = 1.
		y_error[ninds] = 1. / master_volumes_sf.sum(axis=0)[vi][ninds] / dm

		elo = y_error * (ylo / n_error)
		ehi = y_error * (yhi / n_error)


		#pl.errorbar(x[inds], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]], ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i)



		###  resampling from Poisson distribution
		mc_nsf = pylab.poisson(lam=n, size=(nsims, len(n)))
		mci = 0
		for mci in range(nsims):

			mypy.progress_bar(mci, nsims)

			try:
				y_mci = (mc_nsf[mci] / master_volumes_sf.sum(axis=0) / dm)[vi]

				###  adding extra 10% error for non-Poisson sources
				#y_mci = y_mci + y_mci * 0.1 * pylab.randn(len(y_mci))
				#yinds = pylab.find(y_mci < 0)
				#y_mci[yinds] = 0.

				fit_mci, cov_mci = optimize.curve_fit(schechter_mf, lmassbars[inds], y_mci[inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.))
				master_bestfits_sf[vi_counter].append(fit_mci)
			except:
				pass


		bestfits = pylab.array(master_bestfits_sf[vi_counter])


		###  plotting
		gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
		xlo, xhi = 10.2, 11.5
		ylo, yhi = -1.7, -0.1
		nbins = int(nres0 * rfactors[vi_counter])
		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp1.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp1.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])


	print ''







	'''

	###  quick plotting
	for ii in range(len(master_bestfits_sf)):

		color_i = cmap(ii/(len(o_range)-1.))
		bestfits = pylab.array(master_bestfits_sf[ii])

		gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
		xlo, xhi = 10.2, 11.5
		ylo, yhi = -1.7, -0.1

		nbins = int(nres0 * rfactors[ii])

		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp1.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp1.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])

	
	'''







	##############################
	###  plotting for QU galaxies
	##############################

	qu_inds = []   # indices where qu SMFs are mass-complete for each of the overdensity bins

	ydatas = []

	nsims0 = [8000, 5000, 25000, 40000, 40000]
	
	nres0 = 100   # resolution for plotting 2d contours
	rfactors = [1., 1., 1., 0.75, 0.5]

	for vi_counter, vi in enumerate(o_range):

		color_i = cmap(vi_counter/(len(o_range)-1.))
		nsims = nsims0[vi_counter]

		#inds = pylab.find((master_volumes_qu.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		inds = pylab.find((master_number_counts_qu.sum(axis=0)[vi] > 0) & 
			              (master_volumes_qu.sum(axis=0)[vi] > 0) & 
			              (lmassbars < mlim_hi))
		inds = pylab.find((master_number_counts_qu.sum(axis=0)[vi] > 0) & 
			              (master_volumes_qu.sum(axis=0)[vi] > 0) & 
			              (lmassbars < mlim_hi) &
			              (lmassbars > 9.9))

		print '\n\n'
		print 'Quiescent galaxies'
		print '%s' % (voronoi_labels[vi])
		print 'Nsims = %i' % nsims

		qu_inds.append(inds)
		x, y, n = lmassbars, (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0) / dm)[vi], master_number_counts_qu.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)


		#pl.errorbar(x[inds], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]], ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i)



		###  resampling from Poisson distribution
		mc_nqu = pylab.poisson(lam=n, size=(nsims, len(n)))
		mci = 0
		for mci in range(nsims):

			mypy.progress_bar(mci, nsims)

			try:
				y_mci = (mc_nqu[mci] / master_volumes_qu.sum(axis=0) / dm)[vi]

				###  adding extra 10% error for non-Poisson sources
				#y_mci = y_mci + y_mci * 0.1 * pylab.randn(len(y_mci))
				#yinds = pylab.find(y_mci < 0)
				#y_mci[yinds] = 0.

				fit_mci, cov_mci = optimize.curve_fit(schechter_mf, lmassbars[inds], y_mci[inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.))
				master_bestfits_qu[vi_counter].append(fit_mci)
			except:
				pass


		bestfits = pylab.array(master_bestfits_qu[vi_counter])

		###  plotting
		gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
		xlo, xhi = 10.5, 11.8
		ylo, yhi = -1.1, 0.5
		nbins = int(nres0 * rfactors[vi_counter])
		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])


	print ''





	'''

	###  quick plotting
	for ii in range(len(master_bestfits_qu)):

		color_i = cmap(ii/(len(o_range)-1.))
		bestfits = pylab.array(master_bestfits_qu[ii])

		gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
		xlo, xhi = 10.5, 11.8
		ylo, yhi = -1.1, 0.5

		nbins = int(nres0 * rfactors[ii])

		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(bestfits[:,1], bestfits[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])

	'''



	t1 = sp1.text(0.95, 0.95, 'Star-Forming', fontweight='bold', horizontalalignment='right', verticalalignment='top', fontsize=22, color='b', transform=sp1.transAxes)
	t2 = sp2.text(0.95, 0.95, 'Quiescent', fontweight='bold', horizontalalignment='right', verticalalignment='top', fontsize=22, color='r', transform=sp2.transAxes)










	###  plotting 4-panels
	###      top 2: schechter params of SF/QU SMFs
	###      bot 2: schechter params of SF/QU SMFs at >10^10 Msol


	fig = pylab.figure(figsize=(13.375, 11.2))
	sp2 = fig.add_subplot(222)
	sp1 = fig.add_subplot(221)
	sp4 = fig.add_subplot(224)
	sp3 = fig.add_subplot(223)

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp3.minorticks_on()
	sp4.minorticks_on()

	sp3.set_xlabel('log( ${\it M}^*$ / M$_{\odot}$ )')
	sp4.set_xlabel('log( ${\it M}^*$ / M$_{\odot}$ )')
	sp1.set_ylabel(r'$\alpha$', rotation=0, fontsize=26)
	sp2.set_ylabel(r'$\alpha$', rotation=0, fontsize=26)
	sp3.set_ylabel(r'$\alpha$', rotation=0, fontsize=26)
	sp4.set_ylabel(r'$\alpha$', rotation=0, fontsize=26)

	sp1.axis([10.35, 11.55, -1.75, 0.25])
	sp2.axis([10.35, 11.55, -1.75, 0.25])
	sp3.axis([10.35, 11.55, -1.75, 0.25])
	sp4.axis([10.35, 11.55, -1.75, 0.25])

	fig.subplots_adjust(left=0.07, bottom=0.07, top=0.99, right=0.99, hspace=0)


	for vi_counter, vi in enumerate(o_range):

		color_i = cmap(vi_counter/(len(o_range)-1.))
		gauss_levels = pylab.array([0.6065, 0.1353, 0.0111])   # 1,2,3 sigma levels of a Gaussian with peak=1
		print '%i / %i' % (vi_counter+1, len(o_range))


		params_sf = pylab.loadtxt('../data/MCsims_schechter_params_sf%s.txt' % voronoi_bins[vi])
		params_qu = pylab.loadtxt('../data/MCsims_schechter_params_qu%s.txt' % voronoi_bins[vi])
		params_sf_gt_10 = pylab.loadtxt('../data/MCsims_schechter_params_sf%s_mass_gt_10.txt' % voronoi_bins[vi])
		params_qu_gt_10 = pylab.loadtxt('../data/MCsims_schechter_params_qu%s_mass_gt_10.txt' % voronoi_bins[vi])




		###  SF galaxies

		xlo, xhi = 10.2, 11.5
		ylo, yhi = -1.7, -0.1

		nbins = int(nres0 * rfactors[vi_counter])
		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(params_sf[:,1], params_sf[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp1.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp1.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])





		###  QU galaxies

		xlo, xhi = 10.5, 11.8
		ylo, yhi = -1.1, 0.5

		nbins = int(nres0 * rfactors[vi_counter])
		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(params_qu[:,1], params_qu[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp2.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])






		###  SF galaxies >10**10

		xlo, xhi = 10.2, 11.5
		ylo, yhi = -1.7, 0.25

		nbins = int(nres0 * rfactors[vi_counter])
		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(params_sf_gt_10[:,1], params_sf_gt_10[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp3.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp3.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])





		###  QU galaxies >10*10

		xlo, xhi = 10.5, 11.8
		ylo, yhi = -1.1, 0.5

		nbins = int(nres0 * rfactors[vi_counter])
		xgrid, ygrid = pylab.meshgrid(pylab.linspace(xlo, xhi, nbins), pylab.linspace(ylo, yhi, nbins))
		hist2d, xedges, yedges = pylab.histogram2d(params_qu_gt_10[:,1], params_qu_gt_10[:,0], bins=(nbins, nbins), range=([xlo, xhi], [ylo, yhi]))

		levels = hist2d.max() * gauss_levels[0:2]
		contours = sp4.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-','--'])
		contours = sp4.contour(xgrid, ygrid, hist2d.T, levels, 
			                     colors=[color_i, color_i], 
			                     linewidths=[3, 2], 
			                     linestyles=['-',':'])



	t1 = sp1.text(0.95, 0.95, 'Star-Forming', fontweight='bold', horizontalalignment='right', verticalalignment='top', fontsize=22, color='b', transform=sp1.transAxes)
	t2 = sp2.text(0.95, 0.95, 'Quiescent', fontweight='bold', horizontalalignment='right', verticalalignment='top', fontsize=22, color='r', transform=sp2.transAxes)

	t3a = sp3.text(0.95, 0.95, 'Star-Forming', fontweight='bold', horizontalalignment='right', verticalalignment='top', fontsize=22, color='b', transform=sp3.transAxes)
	t4a = sp4.text(0.95, 0.95, 'Quiescent', fontweight='bold', horizontalalignment='right', verticalalignment='top', fontsize=22, color='r', transform=sp4.transAxes)
	t3b = sp3.text(0.95, 0.88, 'log(M$_*$)>10', horizontalalignment='right', verticalalignment='top', fontsize=18, transform=sp3.transAxes)
	t4b = sp4.text(0.95, 0.88, 'log(M$_*$)>10', horizontalalignment='right', verticalalignment='top', fontsize=18, transform=sp4.transAxes)





























######################################
###  Fitting toy models to TOTAL SMFs
######################################

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
voronoi_labels_merged = ['0.0 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 2.0']

#simdata = pickle.load(open('../data/merger_simulation_scenario4.pickle', 'rb'))
simdata = pickle.load(open('../data/merger_simulation_scenario6_7m12.pickle', 'rb'))
#simdata = pickle.load(open('../data/merger_simulation_scenario6_9m12.pickle', 'rb'))
#simdata = pickle.load(open('../merger_simulation_scenario6_test.pickle', 'rb'))
cmap = pylab.cm.cool


###  Field Mass Functions - ZFOURGE 2014
###  /Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table1_TOT.dat
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
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
			#if f.mf_voronoi_slices.zlos[zi] < 0.85: continue

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

			###  iterating over massbins, checking for mass completeness
			for mi in range(len(lmassbars)):

				if lmassbins[mi] > lmasslimit:
					master_number_counts_field[fi][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					master_volumes_field[fi][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_00v05[zi][mi]
					master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]
					#master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					#master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][1][mi] += f.mf_voronoi_slices.ngals_05v10[zi][mi]
					master_volumes[fi][1][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

					master_number_counts[fi][2][mi] += f.mf_voronoi_slices.ngals_10v15[zi][mi]
					master_volumes[fi][2][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

					master_number_counts[fi][3][mi] += f.mf_voronoi_slices.ngals_15v20[zi][mi]
					master_volumes[fi][3][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

					###  counting number of voronoi slices in mass-overdens bins
					if f.mf_voronoi_slices.volumes_00v05_Mpc3[zi] > 0: master_nslices[fi][0][mi] += 1
					if f.mf_voronoi_slices.volumes_05v10_Mpc3[zi] > 0: master_nslices[fi][1][mi] += 1
					if f.mf_voronoi_slices.volumes_10v15_Mpc3[zi] > 0: master_nslices[fi][2][mi] += 1
					if f.mf_voronoi_slices.volumes_15v20_Mpc3[zi] > 0: master_nslices[fi][3][mi] += 1






	###  plotting results

	#fig = pylab.figure(figsize=(8.25, 15.))
	fig = pylab.figure(figsize=(17., 9.))

	sp3 = pylab.subplot2grid((3, 6), (0, 4), rowspan=2, colspan=2)
	sp2 = pylab.subplot2grid((3, 6), (0, 2), rowspan=2, colspan=2)
	sp1 = pylab.subplot2grid((3, 6), (0, 0), rowspan=2, colspan=2)
	sp6 = pylab.subplot2grid((3, 6), (2, 4), rowspan=1, colspan=2)
	sp5 = pylab.subplot2grid((3, 6), (2, 2), rowspan=1, colspan=2)
	sp4 = pylab.subplot2grid((3, 6), (2, 0), rowspan=1, colspan=2)

	sps123 = [sp1, sp2, sp3]
	sps456 = [sp4, sp5, sp6]

	for sp in sps123 + sps456:
		#sp.grid()
		sp.minorticks_on()
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp.set_ylabel('log( Residual )')
		sp.axis([8.4, 11.9, -0.5, 0.5])

	for sp in sps123:
		sp.set_yscale('log')
		sp.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
		sp.axis([8.4, 11.9, 5*10**-4, 4*10**-1])

	for sp in sps456:
		sp.axhline(0, color='k', ls='--')

	fig.subplots_adjust(hspace=0, wspace=0, left=0.09, bottom=0.1, top=0.95)

	ydatas = []
	fits_tot, covs_tot = [], []
	line_fits_tot, line_covs_tot = [], []
	mlim_hi = 11.55
	colors = ['b', 'c', 'm', 'r']
	for vi in [1, 2, 3]:

		color_i = cmap(vi/(len(voronoi_bins)-1.))
		color_i = cmap((vi+1)/4.)

		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)

		###  plotting data
		sps123[vi-1].errorbar(x[inds], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				           ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi], zorder=2)
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  fitting toy models
		inds_model = pylab.find((simdata.lmassax >= x[inds].min()) & (simdata.lmassax <= x[inds].max()))
		scale_factors = pylab.zeros(len(simdata.fraction_merged))
		chi2reds = pylab.zeros(len(simdata.fraction_merged))

		log_errors = (pylab.log10(y/(y-elo)) + pylab.log10((y+ehi)/y)) / 2.

		weights = 1. / ((elo + ehi)/2.)**2
		weights = 1. / log_errors**2


		bestfit_models = pylab.zeros((len(simdata.fraction_merged), len(inds)))


		for i_model in range(len(simdata.fraction_merged)):

			c = cmap(i_model * 1. / (len(simdata.fraction_merged)-1))
			
			numer = pylab.sum(y[inds]  * weights[inds])
			denom = pylab.sum(simdata.ngals_icm[i_model][inds_model] * weights[inds])

			scale = numer / denom
			scale_factors[i_model] = scale
			bestfit_models[i_model] = simdata.ngals_icm[i_model][inds_model] * scale

			chi2reds[i_model] = pylab.sum((pylab.log10(y[inds]) - pylab.log10(simdata.ngals_icm[i_model][inds_model] * scale))**2 / log_errors[inds]**2) / len(inds)
			'''
			scale = pylab.average(y[inds] / simdata.ngals[i_model][inds_model], weights=weights[inds])
			scale_factors[i_model] = scale
			bestfit_models[i_model] = simdata.ngals[i_model][inds_model] * scale

			chi2reds[i_model] = pylab.sum((pylab.log10(y[inds]) - pylab.log10(simdata.ngals[i_model][inds_model] * scale))**2 / log_errors[inds]**2) / len(inds)
			'''

			#sps123[vi-1].plot(simdata.lmassax[inds_model], simdata.ngals[i_model][inds_model] * scale, lw=1.5, color=c)


		Pchi = pylab.exp(-0.5 * chi2reds)
		Pchi = 1. / chi2reds**2

		#Pchi_norm = Pchi / pylab.trapz(Pchi, simdata.fraction_merged)
		#f = pylab.trapz(Pchi_norm * simdata.fraction_merged, simdata.fraction_merged)

		Pchi_norm = Pchi / pylab.trapz(Pchi, simdata.f_icm)
		f = pylab.trapz(Pchi_norm * simdata.f_icm, simdata.f_icm)

		Pchi_norm_scale = Pchi / pylab.trapz(Pchi, scale_factors)
		scale_weighted = pylab.trapz(Pchi_norm_scale * scale_factors, scale_factors)

		print pylab.average(simdata.f_icm, weights=Pchi)

		### cumulative probabilities
		cPchi = pylab.array([pylab.trapz(Pchi_norm[:i], (simdata.f_icm)[:i]) for i in range(len(simdata.fraction_merged))])
		p16 = pylab.interp(0.16, cPchi, simdata.f_icm) * 100
		p50 = pylab.interp(0.50, cPchi, simdata.f_icm) * 100
		p84 = pylab.interp(0.84, cPchi, simdata.f_icm) * 100

		bestfit = pylab.average(bestfit_models, axis=0, weights=Pchi)
		sps123[vi-1].plot(simdata.lmassax[inds_model], bestfit, lw=5, color='k', zorder=1)
		sps123[vi-1].plot(simdata.lmassax[inds_model], bestfit, lw=2, color='gray', zorder=1)


		###  plotting t=0 toy model SMF
		#sps123[vi-1].plot(simdata.lmassax[inds_model], simdata.ngals_icm[0][inds_model] * scale_weighted, lw=1, color='gray', zorder=1)


		###  calculating reduced chi2 of bestfit model
		chi2red_best = pylab.sum((pylab.log10(bestfit) - pylab.log10(y[inds]))**2 / log_errors[inds]**2) / len(inds)
		print '%s,  f_icm=%.3f,  min(chi**2)=%.3f' % (voronoi_labels[vi], f, chi2red_best)



		###  adding text in panel
		#sps123[vi-1].text(0.05, 0.06, '%s\nFraction Mergerd: $%2i^{+%i}_{\;-%i}$%%' % (voronoi_labels[vi], f*100, p84-p50, p50-p16), transform=sps123[vi-1].transAxes)
		sps123[vi-1].set_title(voronoi_labels[vi], y=1.015, fontsize=19)
		sps123[vi-1].text(0.17, 0.06, '$M_{*,\mathrm{ICL}} \, / \, M_{*,\mathrm{total}}$ = %2i%%' % (round(f*100)), transform=sps123[vi-1].transAxes, fontsize=16)

		sps123[vi-1].plot([0.07, 0.14], [0.14, 0.14], lw=5, color='k', transform=sps123[vi-1].transAxes, zorder=1)
		sps123[vi-1].plot([0.07, 0.14], [0.14, 0.14], lw=2, color='gray', transform=sps123[vi-1].transAxes, zorder=1)
		#sps123[vi-1].text(0.17, 0.14, 'Toy model ($\chi^2_{\mathrm{red.}}$ = %.2f)' % chi2red_best, transform=sps123[vi-1].transAxes, verticalalignment='center', fontsize=16)
		sps123[vi-1].text(0.17, 0.14, 'Model ($\chi^2_{\mathrm{red.}}$ = %.2f)' % chi2red_best, transform=sps123[vi-1].transAxes, verticalalignment='center', fontsize=16)


		###  plotting residual
		sps456[vi-1].errorbar(x[inds], pylab.log10(y[inds]) - pylab.log10(bestfit), xerr=dm/2., yerr=[pylab.log10(y / (y-elo))[inds], pylab.log10((y+ehi) / y)[inds]],
				           ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi], zorder=2)








#####################
###  Plotting figure
#####################

inds_lmass = pylab.find((simdata.lmassax > 8.5) & (simdata.lmassax < 11.55))


###  figure for paper
cmap = pylab.cm.jet
fig_paper = pylab.figure(figsize=(8.9, 7.5))
sp_paper = fig_paper.add_subplot(111)
sp_cbar = add_inset(sp_paper, [0.6, 0.85, 0.33, 0.05])
#sp_cbar2 = add_inset(sp_paper, [0.6, 0.65, 0.33, 0.05])

sp_paper.axis([8.1, 11.8, 5*10**1, 3*10**5])

sp_paper.grid()
sp_paper.minorticks_on()
sp_paper.set_yscale('log')
sp_paper.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp_paper.set_ylabel('Number')

#sp_paper.plot(lmass, ngal0, color='k', lw=7, label='initial state', zorder=99)


sp_cbar.set_title('$M_{*,\,\mathrm{ICL}} \, / \, M_{*,\,\mathrm{total}}$', fontsize=18, position=[.5, 1.25])
sp_cbar.set_yticks([])
imdatx, imdaty = pylab.meshgrid(pylab.linspace(simdata.f_icm.min(), simdata.f_icm.max(), int(simdata.f_icm.max()*100)+1), [0, 0])
sp_cbar.imshow(imdatx, cmap=cmap)
sp_cbar.set_aspect('auto')

xticks = sp_cbar.get_xticks()

sp_cbar.set_xticks([0, 6, 12, 18])
sp_cbar.set_xticklabels(['0', '0.06', '0.12', '0.18'], fontsize=15)



###  f_icm to plot
f_icms = pylab.arange(0.02, 0.17, 0.02)
f_icms = pylab.arange(0.03, 0.17, 0.02)
f_mergeds = pylab.zeros(len(f_icms))

for i_f, f in enumerate(f_icms):

	weights = 1. / abs(simdata.f_icm - f)**5
	smf = pylab.average(simdata.ngals_icm, axis=0, weights=weights)
	f_mergeds[i_f] = pylab.average(simdata.fraction_merged, weights=weights)

	c = f / simdata.f_icm.max()
	sp_paper.plot(simdata.lmassax[inds_lmass], smf[inds_lmass], marker='o', ms=4, mew=0.5, color=cmap(c))



	font = FontProperties()
	font.set_family('sans-serif')
	#if i_f == 0:
		#sp_paper.text(8.65-0.2, smf[inds_lmass][0]*2, 'fraction\nmerged', fontproperties=font,
		#	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)

	t = sp_paper.text(8.65-0.05, smf[inds_lmass][0], '%.1f%%' % (f_mergeds[i_f]*100), fontproperties=font,
		          verticalalignment='center', horizontalalignment='right', fontsize=15, fontweight='normal')


smf = simdata.ngals_icm[0]

c = 0.
sp_paper.plot(simdata.lmassax[inds_lmass], smf[inds_lmass],
	          marker='', lw=6, color=cmap(c))


font = FontProperties()
font.set_family('sans-serif')
sp_paper.text(8.65-0.2, smf[inds_lmass][0]*2.3, 'fraction\nmerged', fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)
sp_paper.text(8.65-0.05, smf[inds_lmass][0], '%i%%' % (simdata.fraction_merged[0]*100), fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='right', fontsize=15, fontweight='normal')






































######################################
###  Fitting toy models to TOTAL SMFs
###  modified for 1st referee report
######################################

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
voronoi_labels_merged = ['0.0 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 2.0']

#simdata = pickle.load(open('../data/merger_simulation_scenario4.pickle', 'rb'))
simdata = pickle.load(open('../data/merger_simulation_scenario6_7m12.pickle', 'rb'))
#simdata = pickle.load(open('../data/merger_simulation_scenario6_9m12.pickle', 'rb'))
#simdata = pickle.load(open('../merger_simulation_scenario6_test.pickle', 'rb'))
cmap = pylab.cm.cool


###  Field Mass Functions - ZFOURGE 2014
###  /Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table1_TOT.dat
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
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
			#if f.mf_voronoi_slices.zlos[zi] < 0.85: continue

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

			###  iterating over massbins, checking for mass completeness
			for mi in range(len(lmassbars)):

				if lmassbins[mi] > lmasslimit:
					master_number_counts_field[fi][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					master_volumes_field[fi][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_00v05[zi][mi]
					master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]
					#master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					#master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][1][mi] += f.mf_voronoi_slices.ngals_05v10[zi][mi]
					master_volumes[fi][1][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

					master_number_counts[fi][2][mi] += f.mf_voronoi_slices.ngals_10v15[zi][mi]
					master_volumes[fi][2][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

					master_number_counts[fi][3][mi] += f.mf_voronoi_slices.ngals_15v20[zi][mi]
					master_volumes[fi][3][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

					###  counting number of voronoi slices in mass-overdens bins
					if f.mf_voronoi_slices.volumes_00v05_Mpc3[zi] > 0: master_nslices[fi][0][mi] += 1
					if f.mf_voronoi_slices.volumes_05v10_Mpc3[zi] > 0: master_nslices[fi][1][mi] += 1
					if f.mf_voronoi_slices.volumes_10v15_Mpc3[zi] > 0: master_nslices[fi][2][mi] += 1
					if f.mf_voronoi_slices.volumes_15v20_Mpc3[zi] > 0: master_nslices[fi][3][mi] += 1






	###  plotting results

	#fig = pylab.figure(figsize=(8.25, 15.))
	fig = pylab.figure(figsize=(17., 9.))

	sp3 = pylab.subplot2grid((3, 6), (0, 4), rowspan=2, colspan=2)
	sp2 = pylab.subplot2grid((3, 6), (0, 2), rowspan=2, colspan=2)
	sp1 = pylab.subplot2grid((3, 6), (0, 0), rowspan=2, colspan=2)
	sp6 = pylab.subplot2grid((3, 6), (2, 4), rowspan=1, colspan=2)
	sp5 = pylab.subplot2grid((3, 6), (2, 2), rowspan=1, colspan=2)
	sp4 = pylab.subplot2grid((3, 6), (2, 0), rowspan=1, colspan=2)

	sps123 = [sp1, sp2, sp3]
	sps456 = [sp4, sp5, sp6]

	for sp in sps123 + sps456:
		#sp.grid()
		sp.minorticks_on()
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp.set_ylabel('log( Residual )')
		sp.axis([8.4, 11.9, -0.5, 0.5])

	for sp in sps123:
		sp.set_yscale('log')
		sp.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
		sp.axis([8.4, 11.9, 5*10**-4, 1.5*10**0])

	for sp in sps456:
		sp.axhline(0, color='k', ls='--')

	fig.subplots_adjust(hspace=0, wspace=0, left=0.09, bottom=0.1, top=0.95)

	ydatas = []
	fits_tot, covs_tot = [], []
	line_fits_tot, line_covs_tot = [], []
	mlim_hi = 11.55
	colors = ['b', 'c', 'm', 'r']
	for vi in [1, 2, 3]:

		color_i = cmap(vi/(len(voronoi_bins)-1.))
		color_i = cmap((vi+1)/4.)

		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)

		###  plotting data
		sps123[vi-1].errorbar(x[inds], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				           ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi], zorder=2)
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  fitting toy models
		inds_model = pylab.find((simdata.lmassax >= x[inds].min()) & (simdata.lmassax <= x[inds].max()))
		scale_factors = pylab.zeros(len(simdata.fraction_merged))
		chi2reds = pylab.zeros(len(simdata.fraction_merged))

		log_errors = (pylab.log10(y/(y-elo)) + pylab.log10((y+ehi)/y)) / 2.

		weights = 1. / ((elo + ehi)/2.)**2
		weights = 1. / log_errors**2


		bestfit_models = pylab.zeros((len(simdata.fraction_merged), len(inds)))


		for i_model in range(len(simdata.fraction_merged)):

			c = cmap(i_model * 1. / (len(simdata.fraction_merged)-1))
			
			numer = pylab.sum(y[inds]  * weights[inds])
			denom = pylab.sum(simdata.ngals_icm[i_model][inds_model] * weights[inds])

			scale = numer / denom
			scale_factors[i_model] = scale
			bestfit_models[i_model] = simdata.ngals_icm[i_model][inds_model] * scale

			chi2reds[i_model] = pylab.sum((pylab.log10(y[inds]) - pylab.log10(simdata.ngals_icm[i_model][inds_model] * scale))**2 / log_errors[inds]**2) / len(inds)
			'''
			scale = pylab.average(y[inds] / simdata.ngals[i_model][inds_model], weights=weights[inds])
			scale_factors[i_model] = scale
			bestfit_models[i_model] = simdata.ngals[i_model][inds_model] * scale

			chi2reds[i_model] = pylab.sum((pylab.log10(y[inds]) - pylab.log10(simdata.ngals[i_model][inds_model] * scale))**2 / log_errors[inds]**2) / len(inds)
			'''

			#sps123[vi-1].plot(simdata.lmassax[inds_model], simdata.ngals[i_model][inds_model] * scale, lw=1.5, color=c)


		Pchi = pylab.exp(-0.5 * chi2reds)
		Pchi = 1. / chi2reds**2

		#Pchi_norm = Pchi / pylab.trapz(Pchi, simdata.fraction_merged)
		#f = pylab.trapz(Pchi_norm * simdata.fraction_merged, simdata.fraction_merged)

		Pchi_norm = Pchi / pylab.trapz(Pchi, simdata.f_icm)
		f = pylab.trapz(Pchi_norm * simdata.f_icm, simdata.f_icm)

		Pchi_norm_scale = Pchi / pylab.trapz(Pchi, scale_factors)
		scale_weighted = pylab.trapz(Pchi_norm_scale * scale_factors, scale_factors)

		print pylab.average(simdata.f_icm, weights=Pchi)

		### cumulative probabilities
		cPchi = pylab.array([pylab.trapz(Pchi_norm[:i], (simdata.f_icm)[:i]) for i in range(len(simdata.fraction_merged))])
		p16 = pylab.interp(0.16, cPchi, simdata.f_icm) * 100
		p50 = pylab.interp(0.50, cPchi, simdata.f_icm) * 100
		p84 = pylab.interp(0.84, cPchi, simdata.f_icm) * 100

		bestfit = pylab.average(bestfit_models, axis=0, weights=Pchi)
		sps123[vi-1].plot(simdata.lmassax[inds_model], bestfit, lw=5, color='k', zorder=1)
		sps123[vi-1].plot(simdata.lmassax[inds_model], bestfit, lw=2, color='gray', zorder=1)


		###  plotting t=0 toy model SMF
		if vi == 2:
			xt = simdata.lmassax[inds_model][0] + 0.2
			yt = simdata.ngals_icm[0][inds_model][0] * scale_weighted / 1.041
			t = sps123[vi-1].text(xt, yt, 'initial state', color='gray', rotation=-28)
		sps123[vi-1].plot(simdata.lmassax[inds_model], simdata.ngals_icm[0][inds_model] * scale_weighted, lw=1.5, ls=':', color='gray', zorder=1)
		sps123[vi-1].plot(simdata.lmassax[inds_model], simdata.ngals_icm[0][inds_model] * scale_weighted, lw=1.5, ls='--', color='gray', zorder=1)


		###  calculating reduced chi2 of bestfit model
		chi2red_best = pylab.sum((pylab.log10(bestfit) - pylab.log10(y[inds]))**2 / log_errors[inds]**2) / len(inds)
		print '%s,  f_icm=%.3f,  min(chi**2)=%.3f' % (voronoi_labels[vi], f, chi2red_best)



		###  adding text in panel
		#sps123[vi-1].text(0.05, 0.06, '%s\nFraction Mergerd: $%2i^{+%i}_{\;-%i}$%%' % (voronoi_labels[vi], f*100, p84-p50, p50-p16), transform=sps123[vi-1].transAxes)
		sps123[vi-1].set_title(voronoi_labels[vi], y=1.015, fontsize=19)
		sps123[vi-1].text(0.17, 0.06, '$M_{*,\mathrm{ICL}} \, / \, M_{*,\mathrm{total}}$ = %2i%%' % (round(f*100)), transform=sps123[vi-1].transAxes, fontsize=16)

		sps123[vi-1].plot([0.07, 0.14], [0.14, 0.14], lw=5, color='k', transform=sps123[vi-1].transAxes, zorder=1)
		sps123[vi-1].plot([0.07, 0.14], [0.14, 0.14], lw=2, color='gray', transform=sps123[vi-1].transAxes, zorder=1)
		#sps123[vi-1].text(0.17, 0.14, 'Toy model ($\chi^2_{\mathrm{red.}}$ = %.2f)' % chi2red_best, transform=sps123[vi-1].transAxes, verticalalignment='center', fontsize=16)
		sps123[vi-1].text(0.17, 0.14, 'Model ($\chi^2_{\mathrm{red.}}$ = %.2f)' % chi2red_best, transform=sps123[vi-1].transAxes, verticalalignment='center', fontsize=16)


		###  plotting residual
		sps456[vi-1].errorbar(x[inds], pylab.log10(y[inds]) - pylab.log10(bestfit), xerr=dm/2., yerr=[pylab.log10(y / (y-elo))[inds], pylab.log10((y+ehi) / y)[inds]],
				           ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi], zorder=2)



























######################################
###  Fitting toy models to TOTAL SMFs
######################################

plot_fields = []
plot_fields.append(mf_data('N200',    0.691,  0.027))
plot_fields.append(mf_data('SC1324',  0.755,  0.033))
plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
plot_fields.append(mf_data('N5281',   0.818,  0.029))
plot_fields.append(mf_data('SC1604',  0.910,  0.029))
plot_fields.append(mf_data('SC0910',  1.110,  0.035))
plot_fields.append(mf_data('SC0849',  1.261,  0.029))

voronoi_bins = ['_00v05', '_05v10', '_10v15', '_15v20']
voronoi_labels = ['0.0 < log(1+$\delta$) < 0.5', '0.5 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 1.5', '1.5 < log(1+$\delta$) < 2.0']
voronoi_labels_merged = ['0.0 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 2.0']

#simdata = pickle.load(open('../data/merger_simulation_scenario4.pickle', 'rb'))
simdata = pickle.load(open('../data/merger_simulation_scenario6_7m12.pickle', 'rb'))
#simdata = pickle.load(open('../data/merger_simulation_scenario6_9m12.pickle', 'rb'))
#simdata = pickle.load(open('../merger_simulation_scenario6_test.pickle', 'rb'))
cmap = pylab.cm.cool


###  Field Mass Functions - ZFOURGE 2014
###  /Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table1_TOT.dat
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
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
			#if f.mf_voronoi_slices.zlos[zi] < 0.85: continue

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

			###  iterating over massbins, checking for mass completeness
			for mi in range(len(lmassbars)):

				if lmassbins[mi] > lmasslimit:
					master_number_counts_field[fi][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					master_volumes_field[fi][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_00v05[zi][mi]
					master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_00v05_Mpc3[zi]
					#master_number_counts[fi][0][mi] += f.mf_voronoi_slices.ngals_field[zi][mi]
					#master_volumes[fi][0][mi] += f.mf_voronoi_slices.volumes_field_Mpc3[zi]

					master_number_counts[fi][1][mi] += f.mf_voronoi_slices.ngals_05v10[zi][mi]
					master_volumes[fi][1][mi] += f.mf_voronoi_slices.volumes_05v10_Mpc3[zi]

					master_number_counts[fi][2][mi] += f.mf_voronoi_slices.ngals_10v15[zi][mi]
					master_volumes[fi][2][mi] += f.mf_voronoi_slices.volumes_10v15_Mpc3[zi]

					master_number_counts[fi][3][mi] += f.mf_voronoi_slices.ngals_15v20[zi][mi]
					master_volumes[fi][3][mi] += f.mf_voronoi_slices.volumes_15v20_Mpc3[zi]

					###  counting number of voronoi slices in mass-overdens bins
					if f.mf_voronoi_slices.volumes_00v05_Mpc3[zi] > 0: master_nslices[fi][0][mi] += 1
					if f.mf_voronoi_slices.volumes_05v10_Mpc3[zi] > 0: master_nslices[fi][1][mi] += 1
					if f.mf_voronoi_slices.volumes_10v15_Mpc3[zi] > 0: master_nslices[fi][2][mi] += 1
					if f.mf_voronoi_slices.volumes_15v20_Mpc3[zi] > 0: master_nslices[fi][3][mi] += 1











	###  plotting results

	#fig = pylab.figure(figsize=(8.25, 15.))
	fig = pylab.figure(figsize=(17., 9.))

	sp3 = pylab.subplot2grid((3, 6), (0, 4), rowspan=2, colspan=2)
	sp2 = pylab.subplot2grid((3, 6), (0, 2), rowspan=2, colspan=2)
	sp1 = pylab.subplot2grid((3, 6), (0, 0), rowspan=2, colspan=2)
	sp6 = pylab.subplot2grid((3, 6), (2, 4), rowspan=1, colspan=2)
	sp5 = pylab.subplot2grid((3, 6), (2, 2), rowspan=1, colspan=2)
	sp4 = pylab.subplot2grid((3, 6), (2, 0), rowspan=1, colspan=2)

	sps123 = [sp1, sp2, sp3]
	sps456 = [sp4, sp5, sp6]


	for sp in sps123 + sps456:
		sp.grid()
		sp.minorticks_on()
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp.set_ylabel('log( Residual )')
		sp.axis([8.4, 11.9, -0.5, 0.5])

	for sp in sps123:
		sp.set_yscale('log')
		sp.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
		sp.axis([8.4, 11.9, 5*10**-4, 4*10**-1])

	for sp in sps456:
		sp.axhline(0, color='k', ls='--')

	fig.subplots_adjust(hspace=0, wspace=0, left=0.09, bottom=0.1)



	###  plot for merger counts
	fig_mergers = pylab.figure(figsize=(16.6, 6.2))
	sp3_mergers = fig_mergers.add_subplot(133)
	sp2_mergers = fig_mergers.add_subplot(132)
	sp1_mergers = fig_mergers.add_subplot(131)
	sps_mergers = [sp1_mergers, sp2_mergers, sp3_mergers]
	ymaxes = []
	mlo = 8.9

	fig_mergers.subplots_adjust(left=0.07, top=0.93, wspace=0)

	for sp in sps_mergers:
		sp.grid()
		sp.minorticks_on()
		sp.axis([9.3, 11.9, -0.05, 2.1])
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp.set_ylabel(r'$\langle$ $N_{mergers}$ $\rangle$')

	sp1_mergers.plot([0, 1], [0, 1], lw=2, color='k', marker='s', mfc='w', ms=9, ls='-', label='Major mergers\n1:4 < $\mu_*$ < 1:1')
	sp1_mergers.plot([0, 1], [0, 1], lw=2, color='k', marker='o', mfc='w', ms=6, ls='--', label='Minor mergers\n1:10 < $\mu_*$ < 1:4')

	sp2_mergers.plot([0, 1], [0, 1], lw=3, color='gray', ls='-',  label='Illustris major mergers')
	sp2_mergers.plot([0, 1], [0, 1], lw=3, color='gray', ls='--', label='Illustris minor mergers')

	sp1_mergers.legend(loc=2, fontsize=17)
	sp2_mergers.legend(loc=2, fontsize=17)



	###  arrays to store merger numbers vs. stellar mass for different environment bins
	n_majors_all = []
	n_minors_all = []
	n_vminors_all = []

	ydatas = []
	fits_tot, covs_tot = [], []
	line_fits_tot, line_covs_tot = [], []
	mlim_hi = 11.55
	colors = ['b', 'c', 'm', 'r']
	for vi in [1, 2, 3]:
		sps_mergers[vi-1].set_title(voronoi_labels[vi])

		color_i = pylab.cm.cool(vi/3.)
		color_i = cmap((vi+1)/4.)

		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)

		###  plotting data
		sps123[vi-1].errorbar(x[inds], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				           ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi], zorder=2)
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  fitting toy models
		inds_model = pylab.find((simdata.lmassax >= x[inds].min()) & (simdata.lmassax <= x[inds].max()))
		scale_factors = pylab.zeros(len(simdata.fraction_merged))
		chi2reds = pylab.zeros(len(simdata.fraction_merged))

		log_errors = (pylab.log10(y/(y-elo)) + pylab.log10((y+ehi)/y)) / 2.

		weights = 1. / ((elo + ehi)/2.)**2
		weights = 1. / log_errors**2

		bestfit_models = pylab.zeros((len(simdata.fraction_merged), len(inds)))


		for i_model in range(len(simdata.fraction_merged)):

			c = cmap(i_model * 1. / (len(simdata.fraction_merged)-1))
			
			numer = pylab.sum(y[inds]  * weights[inds])
			denom = pylab.sum(simdata.ngals_icm[i_model][inds_model] * weights[inds])

			scale = numer / denom
			scale_factors[i_model] = scale
			bestfit_models[i_model] = simdata.ngals_icm[i_model][inds_model] * scale

			chi2reds[i_model] = pylab.sum((pylab.log10(y[inds]) - pylab.log10(simdata.ngals_icm[i_model][inds_model] * scale))**2 / log_errors[inds]**2) / len(inds)
			'''
			scale = pylab.average(y[inds] / simdata.ngals[i_model][inds_model], weights=weights[inds])
			scale_factors[i_model] = scale
			bestfit_models[i_model] = simdata.ngals[i_model][inds_model] * scale

			chi2reds[i_model] = pylab.sum((pylab.log10(y[inds]) - pylab.log10(simdata.ngals[i_model][inds_model] * scale))**2 / log_errors[inds]**2) / len(inds)
			'''

			#sps123[vi-1].plot(simdata.lmassax[inds_model], simdata.ngals[i_model][inds_model] * scale, lw=1.5, color=c)


		Pchi = pylab.exp(-0.5 * chi2reds)
		Pchi = 1. / chi2reds**2

		#Pchi_norm = Pchi / pylab.trapz(Pchi, simdata.fraction_merged)
		#f = pylab.trapz(Pchi_norm * simdata.fraction_merged, simdata.fraction_merged)

		Pchi_norm = Pchi / pylab.trapz(Pchi, simdata.f_icm)
		f_icl = pylab.trapz(Pchi_norm * simdata.f_icm, simdata.f_icm)
		f_merged = pylab.average(simdata.fraction_merged, axis=0, weights=Pchi)

		### cumulative probabilities
		cPchi = pylab.array([pylab.trapz(Pchi_norm[:i], (simdata.f_icm)[:i]) for i in range(len(simdata.fraction_merged))])
		p16 = pylab.interp(0.16, cPchi, simdata.f_icm) * 100
		p50 = pylab.interp(0.50, cPchi, simdata.f_icm) * 100
		p84 = pylab.interp(0.84, cPchi, simdata.f_icm) * 100

		bestfit = pylab.average(bestfit_models, axis=0, weights=Pchi)
		sps123[vi-1].plot(simdata.lmassax[inds_model], bestfit, lw=5, color='k', zorder=1)
		sps123[vi-1].plot(simdata.lmassax[inds_model], bestfit, lw=2.5, color='gray', zorder=1)
		print '%s   f_icl=%.1f%%   f_merged=%.2f%%' % (voronoi_labels[vi], f_icl*100, f_merged*100)

		#sps123[vi-1].text(0.05, 0.06, '%s\nFraction Mergerd: $%2i^{+%i}_{\;-%i}$%%' % (voronoi_labels[vi], f*100, p84-p50, p50-p16), transform=sps123[vi-1].transAxes)
		sps123[vi-1].text(0.05, 0.06, '%s\n$M_{*,\mathrm{ICL}} \, / \, M_{*,\mathrm{total}}$ = %2i%%' % (voronoi_labels[vi], round(f_icl*100)), transform=sps123[vi-1].transAxes)


		###  plotting residual
		sps456[vi-1].errorbar(x[inds], pylab.log10(y[inds]) - pylab.log10(bestfit), xerr=dm/2., yerr=[pylab.log10(y / (y-elo))[inds], pylab.log10((y+ehi) / y)[inds]],
				           ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi], zorder=2)





		###  plotting mean number of major/minor mergers vs. lmass
		sp_merger = sps_mergers[vi-1]

		n_majors_marginalized = pylab.average(simdata.n_major_mergers, axis=0, weights=Pchi)
		n_minors_marginalized = pylab.average(simdata.n_minor_mergers, axis=0, weights=Pchi)
		n_vminors_marginalized = pylab.average(simdata.n_vminor_mergers, axis=0, weights=Pchi)

		inds = pylab.find((simdata.lmassax > mlo) & (simdata.lmassax < 11.55))
		n_majors_all.append(pylab.array(zip(simdata.lmassax[inds], n_majors_marginalized[inds])))
		n_minors_all.append(pylab.array(zip(simdata.lmassax[inds], n_minors_marginalized[inds])))
		n_vminors_all.append(pylab.array(zip(simdata.lmassax[inds], n_vminors_marginalized[inds])))

		sp_merger.plot(simdata.lmassax[inds], n_majors_marginalized[inds], color=color_i, lw=2, marker='s', ms=9, ls='-')
		sp_merger.plot(simdata.lmassax[inds], n_minors_marginalized[inds], color=color_i, lw=2, marker='o', ms=6, ls='--')
		#sp_merger.plot(simdata.lmassax[inds], n_vminors_marginalized[inds], color=color_i, lw=2, marker='^', ms=6, ls=':')

		ymaxes.append(max(n_majors_marginalized[inds].tolist() + n_minors_marginalized[inds].tolist()))


	for sp in sps_mergers:
		sp.axis([simdata.lmassax[inds].min()-0.4, 11.9, -0.05, max(ymaxes) * 1.3])






###############################################
###  plotting number of mergers from Illustris
###############################################

z_axis = pylab.linspace(0.8, 5, 100)
t_axis = cosmo.lookback_time(z_axis).value
mu_axis = pylab.linspace(0.1, 1, 100)

dN_dmu_dt = pylab.zeros((len(mu_axis), len(t_axis)))


n_vminor_mergers_illustris = pylab.zeros(len(lmassbars))
n_minor_mergers_illustris = pylab.zeros(len(lmassbars))
n_major_mergers_illustris = pylab.zeros(len(lmassbars))
for i_lmass in range(len(lmassbars)):

	lm = lmassbars[i_lmass]
	n_vminor, n_minor, n_major = 0., 0., 0.
	for i_mu in range(len(mu_axis)-1):
		dmu = mu_axis[i_mu+1] - mu_axis[i_mu]
		for i_t in range(len(t_axis)-1):
			dt = t_axis[i_t+1] - t_axis[i_t]

			l1 = illustris_merger_rates(lmass=lm, mu=mu_axis[i_mu], z=z_axis[i_t])
			l2 = illustris_merger_rates(lmass=lm, mu=mu_axis[i_mu], z=z_axis[i_t+1])
			l3 = illustris_merger_rates(lmass=lm, mu=mu_axis[i_mu+1], z=z_axis[i_t])
			l4 = illustris_merger_rates(lmass=lm, mu=mu_axis[i_mu+1], z=z_axis[i_t+1])

			if 0.01 < (mu_axis[i_mu+1] + mu_axis[i_mu]) / 2. < 0.1:
				n_vminor += pylab.average([l1, l2, l3, l4]) * dmu * dt
			elif 0.1 < (mu_axis[i_mu+1] + mu_axis[i_mu]) / 2. < 0.25:
				n_minor += pylab.average([l1, l2, l3, l4]) * dmu * dt
			else:
				n_major += pylab.average([l1, l2, l3, l4]) * dmu * dt

	n_vminor_mergers_illustris[i_lmass] = n_vminor
	n_minor_mergers_illustris[i_lmass] = n_minor
	n_major_mergers_illustris[i_lmass] = n_major

ymaxes.append(n_major_mergers_illustris.max())




for sp in sps_mergers:
	sp.set_ylim(-0.05, max(ymaxes) * 1.2)
	sp.set_ylim(-0.05, 4.1)
	sp.plot(lmassbars, n_major_mergers_illustris, color='gray', lw=3, ls='-')
	sp.plot(lmassbars, n_minor_mergers_illustris, color='gray', lw=3, ls='--')









#############################################################
###  Plotting Nmergers figure: N_illustris, N_toy versus M*
#############################################################

fig_mergers2 = pylab.figure(figsize=(16.6, 6.2))
sp3_mergers = fig_mergers2.add_subplot(133)
sp2_mergers = fig_mergers2.add_subplot(132)
sp1_mergers = fig_mergers2.add_subplot(131)
sps_mergers = [sp1_mergers, sp2_mergers, sp3_mergers]
ymaxes = []
mlo = 8.9

fig_mergers2.subplots_adjust(left=0.07, top=0.93, wspace=0)

for sp in sps_mergers:
	#sp.grid()
	sp.minorticks_on()
	sp.axis([8.7, 11.9, -0.05, 4.1])
	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp.set_ylabel(r'$\langle$ $N_{mergers}$ $\rangle$')
	#sp.set_ylabel('Fractional merger excess')
	#sp.axhline(1, color='k')

sp1_mergers.plot([0, 1], [0, 1], lw=2, color='k', marker='o', mfc='w', ms=6, ls='--', label='Minor mergers\n1:10 < $\mu_*$ < 1:4')
sp1_mergers.plot([0, 1], [0, 1], lw=2, color='k', marker='s', mfc='w', ms=9, ls='-', label='Major mergers\n1:4 < $\mu_*$ < 1:1')

sp2_mergers.plot([0, 1], [0, 1], lw=3, color='gray', ls='--', label='Minor mergers')
sp2_mergers.plot([0, 1], [0, 1], lw=3, color='gray', ls='-',  label='Major mergers')

leg1 = sp1_mergers.legend(loc=2, fontsize=17)
leg2 = sp2_mergers.legend(loc=2, fontsize=17, title='Field mergers at 0.8 < z < 5')
leg2.get_title().set_fontsize(17)

sp1_mergers.set_title('0.5 < log(1+$\delta_{gal}$) < 1.0', fontsize=18, y=1.015)
sp2_mergers.set_title('1.0 < log(1+$\delta_{gal}$) < 1.5', fontsize=18, y=1.015)
sp3_mergers.set_title('1.5 < log(1+$\delta_{gal}$) < 2.0', fontsize=18, y=1.015)


for i, vi in enumerate([1, 2, 3]):

	###  plotting mean number of major/minor mergers vs. lmass
	sp_merger = sps_mergers[i]
	color_i = cmap((vi+1)/4.)

	nlss_majors = n_majors_all[i][:,1]
	nlss_minors = n_minors_all[i][:,1]
	nfield_majors = pylab.interp(n_majors_all[i][:,0], lmassbars, n_major_mergers_illustris)
	nfield_minors = pylab.interp(n_minors_all[i][:,0], lmassbars, n_minor_mergers_illustris)
	ntot_majors = nlss_majors + nfield_majors
	ntot_minors = nlss_minors + nfield_minors

	sp_merger.plot(n_majors_all[i][:,0], nlss_majors, color=color_i, lw=2, marker='s', ms=9, ls='-')
	sp_merger.plot(n_minors_all[i][:,0], nlss_minors, color=color_i, lw=2, marker='o', ms=6, ls='--')

	sp_merger.plot(n_majors_all[i][:,0], nfield_majors, color='gray', lw=3, ls='-')
	sp_merger.plot(n_minors_all[i][:,0], nfield_minors, color='gray', lw=3, ls='--')












#####################################################################################
###  Plotting revised Nmergers figure: (N_illustris + N_toy) / N_illustris versus M*
#####################################################################################

fig_mergers2 = pylab.figure(figsize=(16.6, 6.2))
sp3_mergers = fig_mergers2.add_subplot(133)
sp2_mergers = fig_mergers2.add_subplot(132)
sp1_mergers = fig_mergers2.add_subplot(131)
sps_mergers = [sp1_mergers, sp2_mergers, sp3_mergers]
ymaxes = []
mlo = 8.9

fig_mergers2.subplots_adjust(left=0.07, top=0.93, wspace=0)

for sp in sps_mergers:
	#sp.grid()
	sp.minorticks_on()
	sp.axis([8.7, 11.9, 0.9, 2.35])
	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	#sp.set_ylabel(r'$\langle$ $N_{mergers}$ $\rangle$')
	sp.set_ylabel('Fractional merger excess')
	sp.axhline(1, color='k')

sp1_mergers.plot([0, 1], [0, 1], lw=2, color='k', marker='o', mfc='w', ms=6, ls='--', label='Minor mergers\n1:10 < $\mu_*$ < 1:4')
sp1_mergers.plot([0, 1], [0, 1], lw=2, color='k', marker='s', mfc='w', ms=9, ls='-', label='Major mergers\n1:4 < $\mu_*$ < 1:1')

#sp2_mergers.plot([0, 1], [0, 1], lw=3, color='gray', ls='--', label='Illustris minor mergers')
#sp2_mergers.plot([0, 1], [0, 1], lw=3, color='gray', ls='-',  label='Illustris major mergers')

sp1_mergers.legend(loc=1, fontsize=17)
#sp2_mergers.legend(loc=2, fontsize=17)

sp1_mergers.set_title('0.5 < log(1+$\delta_{gal}$) < 1.0', fontsize=18, y=1.015)
sp2_mergers.set_title('1.0 < log(1+$\delta_{gal}$) < 1.5', fontsize=18, y=1.015)
sp3_mergers.set_title('1.5 < log(1+$\delta_{gal}$) < 2.0', fontsize=18, y=1.015)




for i, vi in enumerate([1, 2, 3]):


	###  plotting mean number of major/minor mergers vs. lmass
	sp_merger = sps_mergers[i]
	color_i = cmap((vi+1)/4.)

	nlss_majors = n_majors_all[i][:,1]
	nlss_minors = n_minors_all[i][:,1]
	nfield_majors = pylab.interp(n_majors_all[i][:,0], lmassbars, n_major_mergers_illustris)
	nfield_minors = pylab.interp(n_minors_all[i][:,0], lmassbars, n_minor_mergers_illustris)
	ntot_majors = nlss_majors + nfield_majors
	ntot_minors = nlss_minors + nfield_minors

	sp_merger.plot(n_majors_all[i][:,0], ntot_majors / nfield_majors, color=color_i, lw=2, marker='s', ms=9, ls='-')
	sp_merger.plot(n_minors_all[i][:,0], ntot_minors / nfield_minors, color=color_i, lw=2, marker='o', ms=6, ls='--')















#############################################
###  Plotting toy model example SMFs figure
#############################################

inds_lmass = pylab.find((simdata.lmassax > 8.5) & (simdata.lmassax < 11.55))


###  figure for paper
cmap = pylab.cm.jet
fig_paper = pylab.figure(figsize=(8.9, 7.5))
sp_paper = fig_paper.add_subplot(111)
sp_cbar = add_inset(sp_paper, [0.6, 0.85, 0.33, 0.05])
#sp_cbar2 = add_inset(sp_paper, [0.6, 0.65, 0.33, 0.05])

sp_paper.axis([8.1, 11.8, 5*10**1, 3*10**5])

sp_paper.grid()
sp_paper.minorticks_on()
sp_paper.set_yscale('log')
sp_paper.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp_paper.set_ylabel('Number')

#sp_paper.plot(lmass, ngal0, color='k', lw=7, label='initial state', zorder=99)


sp_cbar.set_title('$M_{*,\,\mathrm{ICL}} \, / \, M_{*,\,\mathrm{total}}$', fontsize=18, position=[.5, 1.25])
sp_cbar.set_yticks([])
imdatx, imdaty = pylab.meshgrid(pylab.linspace(simdata.f_icm.min(), simdata.f_icm.max(), int(simdata.f_icm.max()*100)+1), [0, 0])
sp_cbar.imshow(imdatx, cmap=cmap)
sp_cbar.set_aspect('auto')

xticks = sp_cbar.get_xticks()

sp_cbar.set_xticks([0, 6, 12, 18])
sp_cbar.set_xticklabels(['0', '0.06', '0.12', '0.18'], fontsize=15)



###  f_icm to plot
f_icms = pylab.arange(0.02, 0.17, 0.02)
f_icms = pylab.arange(0.03, 0.17, 0.02)
f_mergeds = pylab.zeros(len(f_icms))

for i_f, f in enumerate(f_icms):

	weights = 1. / abs(simdata.f_icm - f)**5
	smf = pylab.average(simdata.ngals_icm, axis=0, weights=weights)
	f_mergeds[i_f] = pylab.average(simdata.fraction_merged, weights=weights)

	c = f / simdata.f_icm.max()
	sp_paper.plot(simdata.lmassax[inds_lmass], smf[inds_lmass], marker='o', ms=4, mew=0.5, color=cmap(c))



	font = FontProperties()
	font.set_family('sans-serif')
	#if i_f == 0:
		#sp_paper.text(8.65-0.2, smf[inds_lmass][0]*2, 'fraction\nmerged', fontproperties=font,
		#	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)

	t = sp_paper.text(8.65-0.05, smf[inds_lmass][0], '%.1f%%' % (f_mergeds[i_f]*100), fontproperties=font,
		          verticalalignment='center', horizontalalignment='right', fontsize=15, fontweight='normal')


smf = simdata.ngals_icm[0]

c = 0.
sp_paper.plot(simdata.lmassax[inds_lmass], smf[inds_lmass],
	          marker='', lw=6, color=cmap(c))


font = FontProperties()
font.set_family('sans-serif')
sp_paper.text(8.65-0.2, smf[inds_lmass][0]*2.3, 'fraction\nmerged', fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)
sp_paper.text(8.65-0.05, smf[inds_lmass][0], '%i%%' % (simdata.fraction_merged[0]*100), fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='right', fontsize=15, fontweight='normal')



































###################################
###  Plot SF/Qui MFs for the LSSs
###################################

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

###  Field Mass Functions - ZFOURGE 2014
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]

zfourge_massfunc_tot = mypy.readcat('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = mypy.readcat('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = mypy.readcat('%s/table1_QUI.dat' % zfourge_dir)

if False:

	outname = '../output/MF_sfqui.pdf'
	pdf = PdfPages(outname)

	fig = pylab.figure(figsize=(12., 9.))
	sp3 = pylab.subplot2grid((3, 2), (0, 1), rowspan=2, colspan=1)
	sp4 = pylab.subplot2grid((3, 2), (2, 1), rowspan=1, colspan=1)
	sp1 = pylab.subplot2grid((3, 2), (0, 0), rowspan=2, colspan=1)
	sp2 = pylab.subplot2grid((3, 2), (2, 0), rowspan=1, colspan=1)
	fig.subplots_adjust(hspace=0, wspace=0, left=0.1)

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp3.minorticks_on()
	sp4.minorticks_on()
	sp1.set_yscale('log')
	sp3.set_yscale('log')

	sp1.grid()
	sp2.grid()
	sp3.grid()
	sp4.grid()

	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number  [dex$^{-1}$]')
	sp2.set_ylabel('QU Fraction')

	sp1.axis([9.2, 11.95, 8, 20000])
	sp2.axis([9.2, 11.95, -0.05, 1.05])
	sp3.axis([9.2, 11.95, 8, 20000])
	sp4.axis([9.2, 11.95, -0.05, 1.05])


	ntot_fields = []
	nsf_fields = []
	nqu_fields = []
	ntot_fields_spec = []
	nsf_fields_spec = []
	nqu_fields_spec = []
	for i in range(len(plot_fields)):
		f = plot_fields[i]
		ntot_fields.append(f.data.Ntot)
		nsf_fields.append(f.data.Nsf)
		nqu_fields.append(f.data.Nqu)
		ntot_fields_spec.append(f.data.Nzspec)
		nsf_fields_spec.append(f.data.Nzspec_sf)
		nqu_fields_spec.append(f.data.Nzspec_qu)
	ntot_fields = pylab.array(ntot_fields)
	nsf_fields = pylab.array(nsf_fields)
	nqu_fields = pylab.array(nqu_fields)
	ntot_fields_spec = pylab.array(ntot_fields_spec)
	nsf_fields_spec = pylab.array(nsf_fields_spec)
	nqu_fields_spec = pylab.array(nqu_fields_spec)




	###  First page for all fields combined
	dm = f.data.lmass[1] - f.data.lmass[0]
	sp1.plot(f.data.lmass, pylab.sum(ntot_fields, axis=0)/dm, color='k', lw=2, marker='o', ms=5, mew=1., label='All')
	sp2.plot(f.data.lmass, pylab.sum(nqu_fields, axis=0) / (pylab.sum(nsf_fields, axis=0) + pylab.sum(nqu_fields, axis=0)), 
		     color='r', lw=2, marker='o', ms=5, mew=1.)
	sp1.plot(f.data.lmass, pylab.sum(nsf_fields, axis=0)/dm, color='b', lw=2, marker='o', ms=5, mew=1., label='SF')
	sp1.plot(f.data.lmass, pylab.sum(nqu_fields, axis=0)/dm, color='r', lw=2, marker='o', ms=5, mew=1., label='QU')

	sp1.legend(loc=3, title='spec+phot', numpoints=1, fontsize=18)



	sp3.plot(f.data.lmass, pylab.sum(ntot_fields_spec, axis=0)/dm, color='k', lw=2, marker='o', ms=5, mew=1., label='All')
	sp4.plot(f.data.lmass, pylab.sum(nqu_fields_spec, axis=0) / (pylab.sum(nsf_fields_spec, axis=0) + pylab.sum(nqu_fields_spec, axis=0)), 
		     color='r', lw=2, marker='o', ms=5, mew=1.)
	sp3.plot(f.data.lmass, pylab.sum(nsf_fields_spec, axis=0)/dm, color='b', lw=2, marker='o', ms=5, mew=1., label='SF')
	sp3.plot(f.data.lmass, pylab.sum(nqu_fields_spec, axis=0)/dm, color='r', lw=2, marker='o', ms=5, mew=1., label='QU')

	sp3.legend(loc=3, title='spec', numpoints=1, fontsize=18)


	pdf.savefig()
	sp1.clear()
	sp2.clear()
	sp3.clear()
	sp4.clear()



	for i in range(len(plot_fields)):

		f = plot_fields[i]
		dm = f.data.lmass[1] - f.data.lmass[0]

		sp1.minorticks_on()
		sp2.minorticks_on()
		sp3.minorticks_on()
		sp4.minorticks_on()
		sp1.set_yscale('log')
		sp3.set_yscale('log')

		sp1.grid()
		sp2.grid()
		sp3.grid()
		sp4.grid()

		sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp1.set_ylabel('Number  [dex$^{-1}$]')
		sp2.set_ylabel('QU Fraction')

		sp1.axis([9.2, 11.95, 8, 20000])
		sp2.axis([9.2, 11.95, -0.05, 1.05])
		sp3.axis([9.2, 11.95, 8, 20000])
		sp4.axis([9.2, 11.95, -0.05, 1.05])





		###  All galaxies
		sp1.plot(f.data.lmass, pylab.sum(ntot_fields, axis=0)/dm, color='k', lw=1)
		sp2.plot(f.data.lmass, pylab.sum(nqu_fields, axis=0) / (pylab.sum(nsf_fields, axis=0) + pylab.sum(nqu_fields, axis=0)), color='r', lw=1)

		ntot = ntot_fields[i]
		elo, ehi = mypy.massfunc.CI(ntot)
		sp1.errorbar(f.data.lmass, ntot/dm, xerr=0/2., yerr=[elo/dm, ehi/dm],
			         ls='', marker='o', ms=9, mfc='k', ecolor='k', zorder=10, label='All')


		sp3.plot(f.data.lmass, pylab.sum(ntot_fields_spec, axis=0)/dm, color='k', lw=1)
		sp4.plot(f.data.lmass, pylab.sum(nqu_fields_spec, axis=0) / (pylab.sum(nsf_fields_spec, axis=0) + pylab.sum(nqu_fields_spec, axis=0)), color='r', lw=1)

		ntot_spec = ntot_fields_spec[i]
		elo, ehi = mypy.massfunc.CI(ntot_spec)
		sp3.errorbar(f.data.lmass, ntot_spec/dm, xerr=0/2., yerr=[elo/dm, ehi/dm],
			         ls='', marker='o', ms=9, mfc='k', ecolor='k', zorder=10, label='All')





		###  SF galaxies
		sp1.plot(f.data.lmass, pylab.sum(nsf_fields, axis=0)/dm, color='b', lw=1)

		nsf = nsf_fields[i]
		elo, ehi = mypy.massfunc.CI(nsf)
		sp1.errorbar(f.data.lmass, nsf/dm, xerr=0/2., yerr=[elo/dm, ehi/dm],
			         ls='', marker='o', ms=9, mfc='b', ecolor='b', zorder=10, label='SF')


		sp3.plot(f.data.lmass, pylab.sum(nsf_fields_spec, axis=0)/dm, color='b', lw=1)

		nsf_spec = nsf_fields_spec[i]
		elo, ehi = mypy.massfunc.CI(nsf_spec)
		sp3.errorbar(f.data.lmass, nsf_spec/dm, xerr=0/2., yerr=[elo/dm, ehi/dm],
			         ls='', marker='o', ms=9, mfc='b', ecolor='b', zorder=10, label='SF')




		###  QU galaxies
		sp1.plot(f.data.lmass, pylab.sum(nqu_fields, axis=0)/dm, color='r', lw=1)

		nqu = nqu_fields[i]
		elo, ehi = mypy.massfunc.CI(nqu)
		sp1.errorbar(f.data.lmass, nqu/dm, xerr=0/2., yerr=[elo/dm, ehi/dm],
			         ls='', marker='o', ms=9, mfc='r', ecolor='r', zorder=10, label='QU')

		sp2.errorbar(f.data.lmass, nqu / (nsf + nqu),
			         ls='', marker='o', ms=9, mfc='r', ecolor='r', zorder=10)


		sp3.plot(f.data.lmass, pylab.sum(nqu_fields_spec, axis=0)/dm, color='r', lw=1)

		nqu_spec = nqu_fields_spec[i]
		elo, ehi = mypy.massfunc.CI(nqu_spec)
		sp3.errorbar(f.data.lmass, nqu_spec/dm, xerr=0/2., yerr=[elo/dm, ehi/dm],
			         ls='', marker='o', ms=9, mfc='r', ecolor='r', zorder=10, label='QU')

		sp4.errorbar(f.data.lmass, nqu_spec / (nsf_spec + nqu_spec),
			         ls='', marker='o', ms=9, mfc='r', ecolor='r', zorder=10)




		###  mass completeness limits
		mlim_tot = pylab.interp(f.zclust, f.masslimits.z, f.masslimits.masslim_empirical)
		mlim_ssp = pylab.interp(f.zclust, f.masslimits.z, f.masslimits.masslim_model_ssp)
		mlim_exp = pylab.interp(f.zclust, f.masslimits.z, f.masslimits.masslim_model_exp)
		sp2.fill_between([0, mlim_tot], -1, 2, color='#cccccc')
		sp4.fill_between([0, mlim_tot], -1, 2, color='#cccccc')
		sp2.fill_between([mlim_tot, mlim_ssp], -1, 2, color='#ffe6e6')
		sp4.fill_between([mlim_tot, mlim_ssp], -1, 2, color='#ffe6e6')
		sp2.fill_between([mlim_tot, mlim_exp], -1, 2, color='#ffb3b3')
		sp4.fill_between([mlim_tot, mlim_exp], -1, 2, color='#ffb3b3')



		sp1.legend(loc=1, title='%s\nspec+phot' % f.name, numpoints=1, fontsize=13)
		sp3.legend(loc=1, title='%s\nspec' % f.name, numpoints=1, fontsize=13)

		pdf.savefig()
		sp1.clear()
		sp2.clear()
		sp3.clear()
		sp4.clear()

	pdf.close()
	pylab.close()
	print '\nwrote to:   %s' % outname










	###  plots of QU fraction and sepc completeness vs M*

	outname2 = '../output/spec_completeness.pdf'
	pdf = PdfPages(outname2)

	fig = pylab.figure(figsize=(8., 9.))
	sp1 = pylab.subplot2grid((2, 1), (0, 0), rowspan=1, colspan=1)
	sp2 = pylab.subplot2grid((2, 1), (1, 0), rowspan=1, colspan=1)
	fig.subplots_adjust(hspace=0, wspace=0)

	sp1.text(0.05, 0.9, 'All fields', fontsize=20, transform=sp1.transAxes)
	sp1.set_ylabel('spec. completeness')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_ylabel('QU Fraction')

	sp1.grid()
	sp2.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()

	sp1.axis([9.2, 11.95, -0.05, 1.05])
	sp2.axis([9.2, 11.95, -0.05, 1.05])



	###  spec completeness
	sp1.plot(f.data.lmass, pylab.sum(nsf_fields_spec, axis=0) / pylab.sum(nsf_fields, axis=0),
		     ls='-', marker='s', ms=7, mew=1.5, lw=1.5, color='b', label='SF')
	sp1.plot(f.data.lmass, pylab.sum(nqu_fields_spec, axis=0) / pylab.sum(nqu_fields, axis=0),
		     ls='-', marker='s', ms=7, mew=1.5, lw=1.5, color='r', label='QU')





	###  QU fractions
	###  ZFOURGE
	zfinds = pylab.find(zfourge_massfunc_qui.lphi_z3 > -98)
	sp2.fill_between(zfourge_massfunc_qui.lmass[zfinds],
		             10**(zfourge_massfunc_qui.lphi_z3[zfinds] - zfourge_massfunc_qui.elo_z3[zfinds]) / 10**(zfourge_massfunc_tot.lphi_z3[zfinds]),
		             10**(zfourge_massfunc_qui.lphi_z3[zfinds] + zfourge_massfunc_qui.ehi_z3[zfinds]) / 10**(zfourge_massfunc_tot.lphi_z3[zfinds]),
		             color='#ffcccc')
	sp2.plot(zfourge_massfunc_qui.lmass[zfinds], 
		     10**(zfourge_massfunc_qui.lphi_z3[zfinds]) / 10**(zfourge_massfunc_tot.lphi_z3[zfinds]),
		     color='r', lw=2, label='ZFOURGE:\n0.75<z<1.00')


	###  ORELSE
	sp2.plot(f.data.lmass, pylab.sum(nqu_fields, axis=0) / (pylab.sum(nsf_fields, axis=0) + pylab.sum(nqu_fields, axis=0)), 
		     mfc='r', color='k', marker='o', ms=9, lw=2, mew=2.5, label='ORELSE')


	sp2.legend(loc=4, fontsize=17)
	pdf.savefig()
	sp1.clear()
	sp2.clear()




	###  looping over fields
	for fi in range(len(nsf_fields)):

		sp1.text(0.05, 0.9, plot_fields[fi].name, fontsize=20, transform=sp1.transAxes)
		sp1.set_ylabel('spec. completeness')
		sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
		sp2.set_ylabel('QU Fraction')

		sp1.grid()
		sp2.grid()
		sp1.minorticks_on()
		sp2.minorticks_on()

		sp1.axis([9.2, 11.95, -0.05, 1.05])
		sp2.axis([9.2, 11.95, -0.05, 1.05])



		###  spec completeness
		sp1.plot(f.data.lmass, nsf_fields_spec[fi] / nsf_fields[fi],
			     ls='-', marker='s', ms=7, mew=1.5, lw=1.5, color='b', label='SF')
		sp1.plot(f.data.lmass, nqu_fields_spec[fi] / nqu_fields[fi],
			     ls='-', marker='s', ms=7, mew=1.5, lw=1.5, color='r', label='QU')





		###  QU fractions
		###  ZFOURGE
		zfinds = pylab.find(zfourge_massfunc_qui.lphi_z3 > -98)
		sp2.fill_between(zfourge_massfunc_qui.lmass[zfinds],
			             10**(zfourge_massfunc_qui.lphi_z3[zfinds] - zfourge_massfunc_qui.elo_z3[zfinds]) / 10**(zfourge_massfunc_tot.lphi_z3[zfinds]),
			             10**(zfourge_massfunc_qui.lphi_z3[zfinds] + zfourge_massfunc_qui.ehi_z3[zfinds]) / 10**(zfourge_massfunc_tot.lphi_z3[zfinds]),
			             color='#ffcccc')
		sp2.plot(zfourge_massfunc_qui.lmass[zfinds], 
			     10**(zfourge_massfunc_qui.lphi_z3[zfinds]) / 10**(zfourge_massfunc_tot.lphi_z3[zfinds]),
			     color='r', lw=2, label='ZFOURGE:\n0.75<z<1.00')


		###  ORELSE
		sp2.plot(f.data.lmass, nqu_fields[fi] / (nsf_fields[fi] + nqu_fields[fi]), 
			     mfc='r', color='k', marker='o', ms=9, lw=2, mew=2.5, label='ORELSE')


		sp2.legend(loc=4, fontsize=17)
		pdf.savefig()
		sp1.clear()
		sp2.clear()

	pdf.close()
	pylab.close()
	print '\nwrote to:   %s' % outname2

























































##############################################
###  Testing redshift variation in SMFs by
###  splitting ORESLE into two zbins.
##############################################

plot_fields = []
#plot_fields.append(mf_data('N200',    0.691,  0.027))
#plot_fields.append(mf_data('SC1324',  0.755,  0.033))
#plot_fields.append(mf_data('RCS0224', 0.772,  0.027))
#plot_fields.append(mf_data('RXJ1716', 0.813,  0.021))
#plot_fields.append(mf_data('N5281',   0.818,  0.029))
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

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

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

	fig = pylab.figure(figsize=(20, 6.8))
	sp3 = fig.add_subplot(133)
	sp2 = fig.add_subplot(132)
	sp1 = fig.add_subplot(131)

	sp1.grid()
	sp2.grid()
	sp3.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()
	sp3.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')
	sp3.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
	sp2.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
	sp3.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')

	sp1.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])
	sp2.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])
	sp3.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])

	sp1.text(0.04, 0.04, 'Total\nz > 0.9', fontweight='bold', color='k', transform=sp1.transAxes, zorder=100)
	sp2.text(0.04, 0.04, 'Star-forming\nz > 0.9', fontweight='bold', color='b', transform=sp2.transAxes, zorder=100)
	sp3.text(0.04, 0.04, 'Quiescent\nz > 0.9', fontweight='bold', color='r', transform=sp3.transAxes, zorder=100)

	fig.subplots_adjust(left=0.09, wspace=0)

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

		inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)
		#sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
		#		     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])

		sp1.fill_between(x[inds], y[inds]-elo[inds], y[inds]+ehi[inds],
				         color=color_i, label=voronoi_labels[vi], zorder=99-vi_counter)



		###  plotting ratio
		#sp2.errorbar(x[inds] + offie[vi_counter], y[inds] / ydatas[0][inds], xerr=dm/2., yerr=[elo[inds] / ydatas[0][inds], ehi[inds] / ydatas[0][inds]],
		#		         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i)


	sp1.legend(loc=3, numpoints=1, fontsize=14)

	#t = sp2.text(0.05, 0.9, 'Total', fontweight='bold', fontsize=22, color='k', transform=sp2.transAxes)
	#sp2.legend(loc=2, title='Total', fontsize=16)










	##############################
	###  plotting for SF galaxies
	##############################



	ydatas = []
	fits_sf, covs_sf = [], []
	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		inds0 = pl.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0) / dm)[vi], master_number_counts_sf.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)

		n_error = n * 1.
		y_error = y * 1.
		ninds = pylab.find(n_error == 0)
		n_error[ninds] = 1.
		y_error[ninds] = 1. / master_volumes_sf.sum(axis=0)[vi][ninds] / dm

		elo = y_error * (ylo / n_error)
		ehi = y_error * (yhi / n_error)
		#sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
		#		     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		#xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)

		sp2.fill_between(x[inds], y[inds]-elo[inds], y[inds]+ehi[inds],
				         color=color_i, label=voronoi_labels[vi], zorder=99-vi_counter)






	##############################
	###  plotting for QU galaxies
	##############################


	qu_inds = []   # indices where qu SMFs are mass-complete for each of the overdensity bins

	ydatas = []
	fits_qu, covs_qu = [], []
	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		#inds = pylab.find((master_volumes_qu.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		inds = pylab.find((master_number_counts_qu.sum(axis=0)[vi] > 0) & 
			              (master_volumes_qu.sum(axis=0)[vi] > 0) & 
			              (lmassbars < mlim_hi))
		qu_inds.append(inds)
		x, y, n = lmassbars, (master_number_counts_qu.sum(axis=0) / master_volumes_qu.sum(axis=0) / dm)[vi], master_number_counts_qu.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)
		#sp3.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
		#		     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		#xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)

		sp3.fill_between(x[inds], y[inds]-elo[inds], y[inds]+ehi[inds],
				         color=color_i, label=voronoi_labels[vi], zorder=99-vi_counter)


































#############################
###  Plot of SF-SMF for Lori
#############################

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


	data_table_tot = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
	data_table_sf  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
	data_table_qu  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99



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

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

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








	##############################
	###  plotting for SF galaxies
	##############################

	fig = pylab.figure(figsize=(8.2, 7.4))
	sp1 = fig.add_subplot(111)

	#sp1.grid()
	#sp2.grid()
	#sp3.grid()
	#sp4.grid()
	sp1.minorticks_on()
	sp1.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')

	sp1.axis([8.4, 11.9, 8*10**-6, 4*10**-1])

	fig.subplots_adjust(left=0.14, bottom=0.11)




	###  plotting ZFOURGE
	lmassbars_zf = zfourge_massfunc_sf[:,0]

	area_zf = 316.   # arcmin**2
	volume_zf_050z075 = mypy.massfunc.vcomoving_slice(area_zf, 0.50, 0.75)
	volume_zf_075z100 = mypy.massfunc.vcomoving_slice(area_zf, 0.75, 1.00)
	volume_zf_100z125 = mypy.massfunc.vcomoving_slice(area_zf, 1.00, 1.25)
	volumes_zf_050z125 = (volume_zf_050z075, volume_zf_075z100, volume_zf_100z125)

	lmf_050z075 = zfourge_massfunc_sf[:,4]
	lmf_075z100 = zfourge_massfunc_sf[:,7]
	lmf_100z125 = zfourge_massfunc_sf[:,10]

	lmf_err_050z075 = (zfourge_massfunc_sf[:,5] + zfourge_massfunc_sf[:,6]) / 2.
	lmf_err_075z100 = (zfourge_massfunc_sf[:,8] + zfourge_massfunc_sf[:,9]) / 2.
	lmf_err_100z125 = (zfourge_massfunc_sf[:,11] + zfourge_massfunc_sf[:,12]) / 2.

	lmf_050z125 = pylab.log10(pylab.average([10**lmf_050z075, 10**lmf_075z100, 10**lmf_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0))
	lmf_err_050z125 = pylab.average([lmf_err_050z075, lmf_err_075z100, lmf_err_100z125], weights=[volume_zf_050z075, volume_zf_075z100, volume_zf_100z125], axis=0)
	lmf_err_050z125 *= (pylab.average(volumes_zf_050z125) / pylab.sum(volumes_zf_050z125))**0.5

	inds_zf = pylab.find((lmf_050z075 > -98) & (lmf_075z100 > -98) & (lmf_100z125 > -98))


	#sp1.fill_between(lmassbars_zf[inds_zf], 10**(lmf_050z125[inds_zf]-lmf_err_050z125[inds_zf]), 10**(lmf_050z125[inds_zf]+lmf_err_050z125[inds_zf]), color='#cccccc', label='ZFOURGE')
	#sp1.axvline(0, color='#cccccc', lw=7, label='ZFOURGE')






	ydatas = []
	fits_sf, covs_sf = [], []
	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		inds0 = pl.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes_sf.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0) / dm)[vi], master_number_counts_sf.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)

		n_error = n * 1.
		y_error = y * 1.
		ninds = pylab.find(n_error == 0)
		n_error[ninds] = 1.
		y_error[ninds] = 1. / master_volumes_sf.sum(axis=0)[vi][ninds] / dm

		elo = y_error * (ylo / n_error)
		ehi = y_error * (yhi / n_error)
		sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  adding to data_table
		data_table_sf[:,0] = lmassbars
		data_table_sf[inds,3*vi_counter+1] = y[inds]
		data_table_sf[inds,3*vi_counter+2] = elo[inds]
		data_table_sf[inds,3*vi_counter+3] = ehi[inds]





		###  fitting (d)schechter functions and plotting

		fit, cov = optimize.curve_fit(schechter_mf, lmassbars[inds], (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0))[vi][inds], p0=[-1, 11, 10**-3], sigma=((elo+ehi)[inds]/2.), maxfev=10**6)
		ymod = schechter_mf(xmod, *fit)

		#xdata, ydata = lmassbars[inds], (master_number_counts_sf.sum(axis=0) / master_volumes_sf.sum(axis=0))[vi][inds]
		#ydata_errs = 10**(elo+ehi)[inds]/2.
		#fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=[11, -1., -1., 10**-3, 10**-4], sigma=ydata_errs)
		#ymod = dschechter(xmod, *fit)
		
		a = sp1.plot(xmod, ymod / dm, color=color_i, lw=1)
		fits_sf.append(fit)
		covs_sf.append(cov.diagonal()**0.5)

	fits_sf = pylab.array(fits_sf)
	covs_sf = pylab.array(covs_sf)

	sp1.legend(loc=3, numpoints=1, fontsize=16)

	#t = sp2.text(0.05, 0.9, 'Star-Forming', fontweight='bold', fontsize=22, color='b', transform=sp2.transAxes)
	t = sp1.text(0.97, 0.97, 'Star-Forming\nGalaxies', fontweight='bold', fontsize=22, color='b', 
		         transform=sp1.transAxes, horizontalalignment='right', verticalalignment='top', multialignment='center')
















































########################################
###  Plots for Keck Science Meeting  ###
########################################

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


	data_table_tot = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
	data_table_sf  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
	data_table_qu  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99



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

			lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
			lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

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

	#sp1.grid()
	#sp2.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$')
	sp2.set_ylabel('Ratio')

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

	chi2s_single = []
	chi2s_double = []

	for vi_counter, vi in enumerate(o_range):
	#for vi in [0]:

		color_i = cmap(vi_counter/(len(o_range)-1.))

		inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
		inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
		x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
		ydatas.append(y)
		ylo, yhi = mypy.massfunc.CI(n)
		elo = y * (ylo / n)
		ehi = y * (yhi / n)
		#sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
		#		     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
				     ls='-', lw=1, color=color_i, mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
		xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



		###  adding to data_table
		data_table_tot[:,0] = lmassbars
		data_table_tot[inds,3*vi_counter+1] = y[inds]
		data_table_tot[inds,3*vi_counter+2] = elo[inds]
		data_table_tot[inds,3*vi_counter+3] = ehi[inds]



		###  plotting ratio
		sp2.errorbar(x[inds] + offie[vi_counter], y[inds] / ydatas[0][inds], xerr=dm/2., yerr=[elo[inds] / ydatas[0][inds], ehi[inds] / ydatas[0][inds]],
				         ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i)



		###  fitting (d)schechter functions and plotting
		try:

			xdata, ydata = lmassbars[inds], (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0))[vi][inds]
			ydata_errs = 10**(elo+ehi)[inds]/2.
			ydata_errs = (elo+ehi)[inds]/2.


			###  fitting log-linear relation to ratio
			fit, cov = optimize.curve_fit(line, xdata, pylab.log10(y[inds] / ydatas[0][inds]), p0=[-8, 0.5], sigma=(elo+ehi)[inds]/2.)
			line_fits_tot.append(fit)
			line_covs_tot.append(cov.diagonal()**0.5)
			if vi > 0:
				#a = sp2.plot(xmod + offie[vi_counter], 10**line(xmod, *fit), color=color_i, lw=1)
				sp2.axvline(0, color=color_i, lw=2, label='slope = %.2f' % fit[0])


			###  fitting schechter funtion to SMF
			fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1., 11., 10**-3], sigma=ydata_errs, maxfev=10**5)
			#ymod = schechter_mf(xmod, *fit)

			chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
			chi2s_single.append(chi2.sum())


			###  fitting dschechter funtion to SMF
			p0 = [11, -1., -1., 10**-3, 10**-4]
			p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
			fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=p0, sigma=ydata_errs, maxfev=10**5)
			ymod = dschechter(xmod, *fit)
			#a = sp1.plot(xmod + offie[vi_counter], ymod / dm, color=color_i, lw=1)
			fits_tot.append(fit)
			covs_tot.append(cov.diagonal()**0.5)

			chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
			chi2s_double.append(chi2.sum())


		except:
			pass



	fits_tot = pylab.array(fits_tot)
	covs_tot = pylab.array(covs_tot)

	chi2s_single = pylab.array(chi2s_single)
	chi2s_double = pylab.array(chi2s_double)

	###  print TeX table of dSchechter params
	if False:
		fits_tot[:,3] *= 10**3
		fits_tot[:,4] *= 10**3
		covs_tot[:,3] *= 10**3
		covs_tot[:,4] *= 10**3
		order = pylab.array([0, 1, 3, 2, 4])
		for vi_counter, vi in enumerate(o_range):
			s = '%s & ' % voronoi_labels[vi]
			for f, c in zip(fits_tot[vi_counter][order], covs_tot[vi_counter][order]):
				s += '$%5.2f\pm%4.2f$ & ' % (f, c)
			print s[:-2] + ' \\\\'


	#sp1.axvline(0, color='#cccccc', lw=7, label='ZFOURGE')
	sp1.legend(loc=3, numpoints=1, fontsize=14)

	#t = sp2.text(0.05, 0.9, 'Total', fontweight='bold', fontsize=22, color='k', transform=sp2.transAxes)
	#sp2.legend(loc=2, title='Total', fontsize=16)




	###  writing data_table to ascii file
	outname = '../data/massFunctions_TOT_Tomczak+2017'
	for f in plot_fields: outname += '_' + f.name
	outer = open(outname + '.dat', 'w')

	outer.write('#    Galaxy stellar mass functions for all galaxies\n')
	outer.write('#\n')
	outer.write('#    Column  Units            Description\n')
	outer.write('#    lmass   Msol             Log stellar mass at center of mass-bin\n')
	outer.write('#    phi     cMpc^-3 dex^-1   Number density of galaxies, normalized by comoving volume\n')
	outer.write('#    elo     cMpc^-3 dex^-1   Lower 1 sigma error in phi\n')
	outer.write('#    ehi     cMpc^-3 dex^-1   Upper 1 sigma error in phi\n')
	outer.write('#\n')
	outer.write('#    Overdensity bins\n')
	outer.write('#    _n05d00    -0.5 < log(1+delta) < 0.0\n')
	outer.write('#    _00d05      0.0 < log(1+delta) < 0.5\n')
	outer.write('#    _05d10      0.5 < log(1+delta) < 1.0\n')
	outer.write('#    _10d15      1.0 < log(1+delta) < 1.5\n')
	outer.write('#    _15d20      1.5 < log(1+delta) < 2.0\n')
	outer.write('#\n')
	outer.write('# lmass')
	for vi_counter, vi in enumerate(o_range):
		if len(voronoi_bins[vi]) == 6:
			outer.write('   phi%s' % (voronoi_bins[vi]))
			outer.write('   elo%s' % (voronoi_bins[vi]))
			outer.write('   ehi%s' % (voronoi_bins[vi]))
		if len(voronoi_bins[vi]) == 7:
			outer.write('  phi%s' % (voronoi_bins[vi]))
			outer.write('  elo%s' % (voronoi_bins[vi]))
			outer.write('  ehi%s' % (voronoi_bins[vi]))
	outer.write('\n')

	for row in data_table_tot:
		outer.write('%6.2f ' % row[0])
		for item in row[1:]:
			outer.write('  %.4e' % item)
		outer.write('\n')

	outer.close()




























































######################
###  Plot for Roy  ###
######################

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


    data_table_tot = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
    data_table_sf  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99
    data_table_qu  = pylab.zeros((len(lmassbars), 1 + 3*(len(overdens_bars)-1))) - 99



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

            lmasslimit = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_empirical)
            lmasslimit_ssp = pylab.interp(f.mf_voronoi_slices.zbars[zi], f.masslimits.z, f.masslimits.masslim_model_ssp)

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

    fig = pylab.figure(figsize=(7.5, 6.8))
    sp1 = fig.add_subplot(111)

    sp1.minorticks_on()
    sp1.set_yscale('log')

    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=20)
    sp1.set_ylabel('Number / cMpc$^{3}$ / dex$^{}$', size=20)

    sp1.axis([8.4, 11.9, 2.8*10**-5, 4*10**-1])

    fig.subplots_adjust(left=0.14)

    ydatas = []
    fits_tot, covs_tot = [], []
    line_fits_tot, line_covs_tot = [], []
    mlim_hi = 11.55



    ###  deciding which overdensity bins to plot
    o_range = range(1, len(overdens_bars))
    offie = pylab.linspace(-len(o_range)*0.01/2, len(o_range)*0.01/2, len(o_range))

    chi2s_single = []
    chi2s_double = []

    for vi_counter, vi in enumerate(o_range):
    #for vi in [0]:

        color_i = cmap(vi_counter/(len(o_range)-1.))

        inds0 = pylab.find((lmassbars >= 9.75) & (lmassbars < mlim_hi))
        inds = pylab.find((master_volumes.sum(axis=0)[vi] > 0) & (lmassbars < mlim_hi))
        x, y, n = lmassbars, (master_number_counts.sum(axis=0) / master_volumes.sum(axis=0) / dm)[vi], master_number_counts.sum(axis=0)[vi]
        ydatas.append(y)
        ylo, yhi = mypy.massfunc.CI(n)
        elo = y * (ylo / n)
        ehi = y * (yhi / n)
        sp1.errorbar(x[inds] + offie[vi_counter], y[inds], xerr=dm/2., yerr=[elo[inds], ehi[inds]],
                     ls='', mew=1.5, marker='o', ms=9, mfc='none', mec=color_i, elinewidth=1.5, ecolor=color_i, label=voronoi_labels[vi])
        xmod = pylab.linspace(lmassbars[inds][0] - dm/2., lmassbars[inds][-1] + dm/2., 1000)



        ###  adding to data_table
        data_table_tot[:,0] = lmassbars
        data_table_tot[inds,3*vi_counter+1] = y[inds]
        data_table_tot[inds,3*vi_counter+2] = elo[inds]
        data_table_tot[inds,3*vi_counter+3] = ehi[inds]






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
            fit, cov = optimize.curve_fit(schechter_mf, xdata, ydata, p0=[-1., 11., 10**-3], sigma=ydata_errs, maxfev=10**5)
            #ymod = schechter_mf(xmod, *fit)

            chi2 = (ydata - schechter_mf(xdata, *fit))**2 / ydata_errs**2
            chi2s_single.append(chi2.sum())


            ###  fitting dschechter funtion to SMF
            p0 = [11, -1., -1., 10**-3, 10**-4]
            p0 = [10.64, 0.12, -1.49, 5.7e-4, 2.6e-4]
            fit, cov = optimize.curve_fit(dschechter, xdata, ydata, p0=p0, sigma=ydata_errs, maxfev=10**5)
            ymod = dschechter(xmod, *fit)
            a = sp1.plot(xmod + offie[vi_counter], ymod / dm, color=color_i, lw=1)
            fits_tot.append(fit)
            covs_tot.append(cov.diagonal()**0.5)

            chi2 = (ydata - dschechter(xdata, *fit))**2 / ydata_errs**2
            chi2s_double.append(chi2.sum())


        except:
            pass



    fits_tot = pylab.array(fits_tot)
    covs_tot = pylab.array(covs_tot)

    chi2s_single = pylab.array(chi2s_single)
    chi2s_double = pylab.array(chi2s_double)

    sp1.legend(loc=3, numpoints=1, fontsize=14)
















