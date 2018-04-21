
import time
import mypy
import numpy
import pickle
import my_classes
from scipy import optimize
from matplotlib import pyplot
from matplotlib.font_manager import FontProperties

print '\nGenerating cosmology...'
t0 = time.time()
from astropy import cosmology
cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)
age_of_universe = cosmo.age(0).value
z_master = numpy.linspace(0., 10., 30000)
t_master = cosmo.lookback_time(z_master).value
tf = time.time()
print 'done!   %i seconds\n' % (tf-t0)


def line(x, m, b): return m * x + b

def dschechter(lmassax, lmstar, a1, a2, phistar1, phistar2):
	factor1 = numpy.log(10) * numpy.exp(-10**(lmassax - lmstar)) * 10**(lmassax - lmstar)
	factor2 = phistar1 * 10**(a1*(lmassax - lmstar)) + phistar2 * 10**(a2*(lmassax - lmstar))
	return factor1 * factor2


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









#####################################################################
###  Reading photometric catalog to obtain reference for galaxy ages
#####################################################################

print '\nReading XLSS005 catalogs for reference ages...'
t0 = time.time()
data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
version = 'xlss005_v0.0.3'
cat_xlss005 = mypy.readcat('%s/%s/%s.cat.gz' % (data_dir, version, version))
fout_xlss005 = mypy.readcat('%s/%s/%s_ZFparams.fout.gz' % (data_dir, version, version))

inds_ages = pylab.find((cat_xlss005.use == 1) & 
	                   (fout_xlss005.z > 0.5) & 
	                   (fout_xlss005.z < 1.5) & 
	                   -numpy.isnan(fout_xlss005.lmass))

lmass_xlss005 = fout_xlss005.lmass[inds_ages]
lage_xlss005 = fout_xlss005.lage[inds_ages]

del cat_xlss005, fout_xlss005
tf = time.time()
print 'done!   %i seconds\n' % (tf-t0)








def fquiescent_leja2015(z, lmass):
	'''
	Returns the Quiescent fraction at the given redshift
	and stellar mass making use of the parameterizations
	from Leja+2015. In brief, this function calculates
	the QUIESCENT and STAR-FORMING stellar mass functions
	at z, lmass and returns the ratio: QU / (QU + SF)
	'''
	#logphi1_tot  = -2.64 + 0.07 * z - 0.28 * z**2
	#logphi2_tot  = -3.11 - 0.18 * z - 0.03 * z**2
	#logMstar_tot = 10.72 - 0.13 * z + 0.11 * z**2
	#alpha1_tot = -0.39
	#alpha2_tot = -1.53

	logphi1_sf  = -2.88 + 0.11 * z - 0.31 * z**2
	logphi2_sf  = -3.48 + 0.07 * z - 0.11 * z**2
	logMstar_sf = 10.67 - 0.02 * z + 0.10 * z**2
	alpha1_sf = -0.97
	alpha2_sf = -1.58

	logphi1_qu  = -2.51 - 0.33 * z - 0.07 * z**2
	logphi2_qu  = -3.54 - 2.31 * z + 0.27 * z**2
	logMstar_qu = 10.70
	alpha1_qu = -0.10
	alpha2_qu = -1.69

	smf_sf = dschechter(lmass, logMstar_sf, alpha1_sf, alpha2_sf, 10**logphi1_sf, 10**logphi2_sf)
	smf_qu = dschechter(lmass, logMstar_qu, alpha1_qu, alpha2_qu, 10**logphi1_qu, 10**logphi2_qu)

	return smf_qu / (smf_qu + smf_sf)


def sfrmass_relation_tomczak2016(z, lmass):
	'''
	Returns the SFR at the corresponding
	log(M*) and z for star-forming galaxies
	from Tomczak et al. (2016)
	'''
	s0    = 0.448 + 1.220 * z - 0.174 * z**2
	logM0 = 9.458 + 0.865 * z - 0.132 * z**2
	gamma = 1.091

	logSFR = s0 - numpy.log10(1 + (10**lmass / 10**logM0)**-gamma)
	return 10**logSFR



def f_loss(t):
	'''
	Returns the fraction of stellar mass lost
	due to passive evolution for a stellar 
	population after a time interval t (in yr).
	From Equation 16 of Moster+2012
	'''
	return 0.05 * numpy.log((t + 3.e5) / 3.e5)





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
	factor3 = mu**(betaz + gamma * numpy.log10(10**lmass / 10**10))

	return factor1 * factor2 * factor3









QUENCHING_TIMESCALE_DATA = mypy.readcat('../data/figure8_Fillingham+2017.dat')

class galaxy:
	'''
	Description:
	    This class is designed to encapsulate all of the 
	    relevant information of a single simulated galaxy.

	Attributes:
	    z_array                   # Array of redshifts at each time-step
	    t_array                   # Array of lookback-times at each time-step in Gyr
	    lmass0                    # Initial seed stellar mass at z=z0
	    lage0                     # Log age of seed stellar mass at z=z0
	    lmass_form                # Stellar mass at formation inferred from lage0
	    mass_in_situ              # Array of stellar mass formed in situ at each time-step
	    mass_ex_situ              # Array of stellar mass accreted ex situ at each time-step
	    lage_ex_situ              # Array of 
	    lmass_total               # Array of total stellar mass at each time-step
	    quenched                  # Boolean for if galaxy is star-forming or quenched
	    sfr                       # Current star-formation rate in Msol/yr
	    n_vminor_mergers          # Total number of very minor mergers experienced
	    n_minor_mergers           # Total number of minor mergers experienced
	    n_major_mergers           # Total number of major mergers experienced
	    mass_from_vminor_mergers  # Total mass accreted from very minor mergers
	    mass_from_minor_mergers   # Total mass accreted from minor mergers
	    mass_from_major_mergers   # Total mass accreted from major mergers
	    merged_in_this_timestep   # A flag that indicates if a galaxy was involved in a merger
	    quenching_timescale       # Quenching time-scale defined at galaxy seed mass in Gyr
	    quenching_countdown       # Randomly generated time interval after which this galaxy will be quenched in Gyr
	'''
	def __init__(self, z0, lmass0, lage0=None, quenched=False):

		self.z_array = numpy.array([z0])
		self.t_array = cosmo.lookback_time([z0]).value
		self.lmass0 = lmass0
		self.lmass_total = numpy.array([lmass0])
		self.lage0 = lage0

		self.lmass_form = numpy.log10(10**lmass0 / (1 - f_loss(10**lage0)))

		self.mass_in_situ = numpy.array([0.])
		self.mass_ex_situ = numpy.array([0.])
		self.lage_ex_situ = numpy.array([0.])

		self.quenched = numpy.array([quenched])

		if quenched:
			self.sfr = numpy.array([0.])
		else:
			self.sfr = numpy.array([sfrmass_relation_tomczak2016(z0, lmass0)])

		self.n_vminor_mergers = 0
		self.n_minor_mergers = 0
		self.n_major_mergers = 0
		self.mass_from_vminor_mergers = 0.
		self.mass_from_minor_mergers = 0.
		self.mass_from_major_mergers = 0.

		self.merged_in_this_timestep = False

		self.quenching_timescale = numpy.interp(lmass0, QUENCHING_TIMESCALE_DATA.lmass, QUENCHING_TIMESCALE_DATA.qenching_timescale)
		self.quenching_countdown = -1
		while self.quenching_countdown < 0:
			self.quenching_countdown = numpy.random.randn() * 1. + self.quenching_timescale



	def get_lmass(self):
		'''
		Returns the current total stellar mass:

		total  =  (seed mass) + (ex situ mass) + (in situ mass) - (mass losses)
		'''

		age_of_seed_mass = 10**g.lage0 + (self.t_array[0] - self.t_array[-1])*10**9

		mass_summation = 10**self.lmass_form * (1 - f_loss(age_of_seed_mass))
		mass_summation += self.mass_ex_situ.sum()

		###  adding past bursts of SF with mass-loss
		for i in range(1, len(self.mass_in_situ)):
			dt = self.t_array[i-1] - self.t_array[-1]
			mass_summation += self.mass_in_situ[i] * (1 - f_loss(dt*10**9))
			mass_summation += self.mass_ex_situ[i] * (1 - f_loss(dt*10**9 + 10**lage_ex_situ[i]))

		return numpy.log10(mass_summation)



class mock_catalog:
	'''
	Description:
	    This class is designed to encapsulate all of the 
	    information of a toy model simulation from start 
	    to finish and every time-step inbetween.

	Attributes:
	    lmassbins          # Stellar mass bin boundaries for SMF
	    lmassbars          # Stellar mass bin centers for SMF

	    z_start   = 5.     # Redshift at which the simulation "begins"
	    z_end     = 0.8    # Redshift at which the simulation "ends"
	    timestep  = 100.   # Duration of each time-step in Myr

	    t_start            # Lookback-time at z_start in Gyr
	    t_end              # Lookback-time at z_end in Gyr
	    t_array            # Array of lookback-times at each time-step in Gyr
	    z_array            # Array of redshifts at each time-step
	    n_timesteps        # Total number of time-steps in simulation

	    galaxies           # Array containing simulated galaxies
	    ngal_initial       # Total number of galaxies from simulation start
	    mass_running       # Total amount of stellar mass at each time-step
	    mass_in_mergers    # Total amount of stellar mass from merging at each time-step
	    mass_in_ICL        # Total amount of stellar mass in the ICL at each time-step
	    fquiescent         # Array contining fquiescent in each mass-bin at each time-step
	    nstep_mergers      # Number of mergers to occur for each time-step

	    SMFs               # Array contaning the numbers of gals per massbin for each time-step
	'''
	def __init__(self, lmassbins, lmassbars, z_start=5., z_end=0.8, timestep=100.):

		self.z_start  = z_start
		self.z_end    = z_end
		self.timestep = timestep

		self.t_start = cosmo.lookback_time(z_start).value
		self.t_end   = cosmo.lookback_time(z_end).value

		self.t_array = numpy.arange(self.t_start, self.t_end, -timestep/10.**3)
		self.z_array = numpy.interp(self.t_array, t_master, z_master)
		self.n_timesteps = len(self.t_array)



		############################################
		###  Obtaining merger rates from Illustris
		###
		###  nstep_mergers contains the RELATIVE
		###  number of mergers for each time-step
		###  in the simulation, integrated from
		###  Table 1 of Rodriguez-Gomez+2015
		############################################

		mu_axis = numpy.linspace(0.1, 1, 100)
		self.nstep_mergers = numpy.array([])

		for i_t in range(len(self.t_array)):

			z = self.z_array[i_t]
			n_mergers = 0

			for i_mu in range(len(mu_axis)-1):
				dmu = mu_axis[i_mu+1] - mu_axis[i_mu]
				for i_lmass in range(len(lmassbars)-1):
					dlmass = lmassbars[i_lmass+1] - lmassbars[i_lmass]

					corner1 = illustris_merger_rates(lmass=lmassbars[i_lmass+1], mu=mu_axis[i_mu],   z=z)
					corner2 = illustris_merger_rates(lmass=lmassbars[i_lmass],   mu=mu_axis[i_mu],   z=z)
					corner3 = illustris_merger_rates(lmass=lmassbars[i_lmass+1], mu=mu_axis[i_mu+1], z=z)
					corner4 = illustris_merger_rates(lmass=lmassbars[i_lmass],   mu=mu_axis[i_mu+1], z=z)

					n_mergers += numpy.average([corner1, corner2, corner3, corner4]) * dmu * (timestep/10.**3)

			self.nstep_mergers = numpy.append(self.nstep_mergers, n_mergers)

		###  normalizing
		self.nstep_mergers /= self.nstep_mergers.sum()



		self.galaxies = []
		self.ngal_initial = 0
		self.mass_running    = numpy.zeros(self.n_timesteps)
		self.mass_in_mergers = numpy.zeros(self.n_timesteps)
		self.mass_in_ICL     = numpy.zeros(self.n_timesteps)

		self.fquiescent = []

		self.lmassbins = lmassbins
		self.lmassbars = lmassbars
		self.SMFs = []

	def get_galaxy_lmasses(self):
		'''
		Returns an array of the masses of the current galaxies list
		'''
		return numpy.array([g.get_lmass() for g in self.galaxies])

	def get_f_ICL(self):
		'''
		Returns an array of the fraction of total stellar mass in the ICL at each time-step
		'''
		return self.mass_in_ICL / self.mass_running

	def generate_SMF(self):
		'''
		Generates the current SMF and appends it to the SMFs attribute
		'''
		lmasses = numpy.array(self.get_galaxy_lmasses())
		digi = numpy.digitize(lmasses, self.lmassbins)
		bincount = numpy.bincount(digi, minlength=len(self.lmassbins)+1)
		self.SMFs.append(bincount[1:-1])

	def generate_fquiescent(self):
		'''
		Generates the current fquiescent vs. stellar mass
		'''

		self.fquiescent.append(numpy.zeros(len(self.lmassbars)))

		lmasses = numpy.array(self.get_galaxy_lmasses())
		digi = numpy.digitize(lmasses, self.lmassbins)

		quenched_array = numpy.array([g.quenched[-1] for g in self.galaxies])
		for i in range(1, len(self.lmassbins)):
			inds_lmassbin = numpy.where(digi == i)[0]
			inds_lmassbin_quenched = numpy.where(quenched_array[inds_lmassbin])[0]

			if len(inds_lmassbin) > 0:
				self.fquiescent[-1][i-1] = len(inds_lmassbin_quenched) * 1. / len(inds_lmassbin)











######################
###  SMF of the field
######################

params_n05v00 = [10.77, -0.14, -1.52, 0.00021, 0.00014]
params_00v05  = [10.87, -0.59, -1.60, 0.00045, 0.00012]

dm = 0.25
lmassbins = numpy.arange(7.-dm/2, 11.55+dm/2., dm)
lmassbars = (lmassbins[:-1] + lmassbins[1:]) / 2.
phi_n05v00 = dschechter(lmassbars, *params_n05v00)
phi_00v05  = dschechter(lmassbars, *params_00v05)
phi_n05v05 = (phi_n05v00 + phi_00v05) / 2.
phi_n05v05_normed = (1 * phi_n05v05 / phi_n05v05.min()).astype(int)



####################################
###  Generating simulation object(s)
####################################

print 'Initializing mock catalog...'
t0 = time.time()

simulation = mock_catalog(lmassbins, lmassbars)
age_at_z_start = age_of_universe - simulation.t_start

###  generating simulated galaxies
for i in range(len(phi_n05v05_normed)):
	ni = phi_n05v05_normed[i]
	mi = lmassbars[i]
	for j in range(ni):

		###  generating quasi-randomized stellar mass
		m_rand = numpy.random.rand() * dm + (mi - dm/2.)

		###  quenching galaxies in order to match the quiescent
		###  fraction vs. stellar mass at z_start
		q = False
		if numpy.random.rand() < fquiescent_leja2015(simulation.z_start, m_rand):
			q = True

		###  grabbing a random age, not to exceed the 
		###  age of the universe at z_start
		inds = numpy.where((lmass_xlss005 > m_rand-dm/2.) &
			               (lmass_xlss005 < m_rand+dm/2.) &
			               (lage_xlss005 < numpy.log10(age_at_z_start * 10**9)))[0]
		lage = lage_xlss005[inds][numpy.random.randint(0, len(inds))]


		simulation.galaxies.append(galaxy(simulation.z_start, m_rand, lage0=lage, quenched=q))
		simulation.ngal_initial += 1

###  randomizing order of galaxies list
numpy.random.shuffle(simulation.galaxies)


###  storing initial fquiescent vs. stellar mass
simulation.generate_fquiescent()


###  storing initial SMF of simulation
simulation.generate_SMF()
simulation.mass_running[0] += (10**simulation.get_galaxy_lmasses()).sum()


###  setting the number of mergers that should occur at each time-step.
###    (1) Define such that at the end the number of galaxies
###        is 1% of the number at the simulation start.
###    (2) Distribute in a way to match the merger-rates
###        found in Illustris as a function of redshift
simulation.nstep_mergers *= (simulation.ngal_initial * 0.995)
simulation.nstep_mergers = simulation.nstep_mergers.astype(int)


tf = time.time()
print 'done!   %i seconds\n' % (tf-t0)


















#########################
###  Running simulation
#########################

print 'Running merger simulation...'
t0 = time.time()

for i_timestep in range(1, simulation.n_timesteps):



	mypy.progress_bar(i_timestep, simulation.n_timesteps)




	###########################################
	###  updating properties for all galaxies
	###########################################

	for i_galaxy in range(len(simulation.galaxies)):

		g = simulation.galaxies[i_galaxy]

		###  reset merger flag
		g.merged_in_this_timestep = False

		###  update lmass, z, and t arrays
		g.z_array = numpy.append(g.z_array, simulation.z_array[i_timestep])
		g.t_array = numpy.append(g.t_array, simulation.t_array[i_timestep])
		g.quenched = numpy.append(g.quenched, g.quenched[-1])


		#############################################
		###  grow galaxy via in situ star-formation
		#############################################

		mass_growth = g.sfr[-1] * (simulation.timestep * 10**6)
		g.mass_in_situ = numpy.append(g.mass_in_situ, mass_growth)


		######################################
		###  adding element to ex situ arrays
		######################################

		g.mass_ex_situ = numpy.append(g.mass_ex_situ, 0.)
		g.lage_ex_situ = numpy.append(g.lage_ex_situ, 0.)





		###  update SFR array if not already quenched
		if not g.quenched[-1]:
			g.sfr = numpy.append(g.sfr, sfrmass_relation_tomczak2016(simulation.z_array[i_timestep], g.get_lmass()))
		else:
			g.sfr = numpy.append(g.sfr, 0.)





	
	######################################################
	###  Quenching galaxies   -   Matching fquiescent(z)
	######################################################

	###  quenching galaxies at random to match the
	###  fquiescent vs. stellar mass from Leja+2015
	
	fquiescent_from_leja2015 = fquiescent_leja2015(simulation.z_array[i_timestep], simulation.lmassbars)
	d_fquiescent = fquiescent_from_leja2015 - simulation.fquiescent[-1]
	digi_lmass = numpy.digitize(simulation.get_galaxy_lmasses(), simulation.lmassbins)

	for i_lmassbin in range(len(simulation.lmassbars)):

		###  if fquiescent in this bin is already above
		###  the prediction from Leja+2015 then pass
		if d_fquiescent[i_lmassbin] <= 0: continue


		###  calculate the necessary number of galaxies to quench
		###  in this mass-bin to match Leja+2015
		n2quench = int(simulation.SMFs[-1][i_lmassbin] * d_fquiescent[i_lmassbin] + 0.5)


		###  identify galaxies in this mass-bin
		inds_lmassbin = numpy.where(digi_lmass == i_lmassbin+1)[0]


		###  quench n2quench galaxies that are not already quenched
		while n2quench > 0:
			for i_galaxy in inds_lmassbin:
				if n2quench == 0:
					break
				if not simulation.galaxies[i_galaxy].quenched[-1]:
					simulation.galaxies[i_galaxy].quenched[-1] = True
					n2quench -= 1
	




	########################################################
	###  Quenching galaxies   -   Quenching countdown timer
	########################################################

	dt = simulation.t_array[0] - simulation.t_array[i_timestep]
	for g in simulation.galaxies:

		if dt > g.quenching_countdown:
			g.quenched[-1] = True










	######################################################################
	###  iterating through necessary number of mergers for this time-step
	######################################################################

	for i_merger in range(simulation.nstep_mergers[i_timestep-1]):


		###  identifying an "acceptable" pair of galaxies to merge.
		###  need to enforce minor mergers occur at 3x the rate
		###  of major mergers.
		ready_to_merge = False
		while not ready_to_merge:

			###  grabbing random galaxies
			ij_rand = numpy.random.randint(0, len(simulation.galaxies), 2)

			###  ... but not the same galaxy twice
			while ij_rand[0] == ij_rand[1]:
				ij_rand = numpy.random.randint(0, len(simulation.galaxies), 2)

			###  ... and not if either galaxy was already involved in a merger in this time-step
			while simulation.galaxies[ij_rand[0]].merged_in_this_timestep or simulation.galaxies[ij_rand[1]].merged_in_this_timestep:
				ij_rand = numpy.random.randint(0, len(simulation.galaxies), 2)



			###  checking merger probability
			lmasses = numpy.array([simulation.galaxies[ij_rand[0]].get_lmass(), 
				                   simulation.galaxies[ij_rand[1]].get_lmass()])
			mratio = 10**(lmasses.min() - lmasses.max())

			###  definition of major merger 1:4
			if mratio > 0.25 and numpy.random.rand() < 1./3:
				ready_to_merge = True
			elif mratio < 0.25:
				ready_to_merge = True


		###  now we have selected our galaxies to merge
		###  merging galaxy masses, with 30% mass-loss to ICL
		i_more_massive = ij_rand[lmasses == lmasses.max()][0]
		i_less_massive = ij_rand[lmasses == lmasses.min()][0]

		gal0 = simulation.galaxies[i_more_massive]
		gal1 = simulation.galaxies[i_less_massive]

		#simulation.galaxies[i_more_massive].lmass[-1] = numpy.log10(10**lmasses[0] + 
		#	                                                        10**lmasses[1] -
		#	                                                        0.3 * 10**lmasses.min())

		gal0.mass_ex_situ[-1] = 10**gal1.lmass_form * 0.3
		gal0.lage_ex_situ[-1] = gal1.lage0
		gal0.mass_in_situ += gal1.mass_in_situ * 0.3


		###  now that this galaxy was involved in a merger, it quenches !!!
		#gal0.quenched[-1] = True


		gal0.merged_in_this_timestep = True


		###  incrementing merger counter of galaxy
		if 0.01 <= mratio < 0.1:   gal0.n_vminor_mergers += 1
		elif 0.1 <= mratio < 0.25: gal0.n_minor_mergers += 1
		elif mratio >= 0.25:       gal0.n_major_mergers += 1


		###  adding mass to ICL from merger
		simulation.mass_in_ICL[i_timestep:] += 0.3 * 10**gal1.get_lmass()

		###  popping merged galaxy
		pop = simulation.galaxies.pop(i_less_massive)


	###  done with merging



	#################################
	###  Updating lmass_total arrays
	#################################
	for i_galaxy in range(len(simulation.galaxies)):
		g = simulation.galaxies[i_galaxy]
		g.lmass_total = numpy.append(g.lmass_total, g.get_lmass())





	###  recording total running stellar mass
	simulation.mass_running[i_timestep] += (10**simulation.get_galaxy_lmasses()).sum()
	simulation.mass_running[i_timestep] += simulation.mass_in_ICL[i_timestep]


	###  storing SMF and fquiescent
	simulation.generate_SMF()
	simulation.generate_fquiescent()




tf = time.time()
dt = (tf-t0)
hours = int(dt / 3600)
minutes = int((dt % 3600) / 60)
seconds = int((dt % 3600) % 60)
print '\n\ndone!   %02ih %02im %02is\n' % (hours, minutes, seconds)




















def plot_diagnostic(g):
	'''
	Generats a diagnostic plot of a galaxy object
	'''

	fig = pyplot.figure(figsize=(13., 8.))
	sp1 = fig.add_subplot(211)
	sp2 = fig.add_subplot(212)

	fig.subplots_adjust(hspace=0)

	sp1.grid()
	sp2.grid()

	sp1.minorticks_on()
	sp2.minorticks_on()

	sp1.set_ylabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('lookback time [Gyr]')
	sp2.set_ylabel('log( SFR  [M$_{\odot}$ / yr] )')

	dt = abs(g.t_array[0] - g.t_array[1])
	sp1.set_xlim(g.t_array.min() - 10*dt, g.t_array.max() + 10*dt)
	sp2.set_xlim(g.t_array.min() - 10*dt, g.t_array.max() + 10*dt)


	###  plotting masses
	sp1.plot(g.t_array, g.lmass_total, ls='-', marker='o', color='b', label='Total')

	inds = numpy.where(g.mass_in_situ > 0)[0]
	sp1.plot(g.t_array[inds], numpy.log10(g.mass_in_situ[inds]),
		     ls='-', marker='o', color='c', label='in situ')


	###  plotting SFRs
	inds = numpy.where(g.sfr > 0)[0]
	sp2.plot(g.t_array[inds], numpy.log10(g.sfr[inds]),
		     ls='-', marker='o', color='b', label='in situ')


	###  plotting vlines at mergers
	sp1. axvline(0, color='gray', lw=2, label='merger')
	for i in range(len(g.t_array)):
		if g.mass_ex_situ[i] > 0:
			sp1.axvline(g.t_array[i], color='gray', lw=2)


	#   plotting vline at quenching
	for i in range(len(g.t_array)):
		if g.quenched[i]:
			sp1.axvline(g.t_array[i], color='r', lw=2, label='quench')
			sp2.axvline(g.t_array[i], color='r', lw=2, label='quench')
			break

	return fig, sp1, sp2



