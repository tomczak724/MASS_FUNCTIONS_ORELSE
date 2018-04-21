
import sys
import time
import mypy
import numpy
import pylab
import pickle
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



F_MASS_LOSS = 0.3  # Fraction of stellar mass to loss per merger event


def line(x, m, b): return m * x + b

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

inds_ages = numpy.where((cat_xlss005.use == 1) & 
	                    (fout_xlss005.z > 0.5) & 
	                    (fout_xlss005.z < 1.5) & 
	                    -numpy.isnan(fout_xlss005.lmass))[0]

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











class stellar_mass_particle(object):
	'''
	Description:
	    This class effectively represents a constituent
	    unit of stellar mass for the galaxy class.

	Attributes:
	    t_form       # Lookback time at formation in Gyr
	    mass_form    # Stellar mass of particle at formation
	'''
	def __init__(self, t_form, mass_form):
		self.t_form = t_form
		self.mass_form = mass_form

	def get_mass(self, t):
		'''
		Returns the stellar mass at lookback time t (in Gyr)
		accounting for mass loss du to passive evolution.
		'''
		dt = self.t_form - t
		if dt < 0:
			raise(IOError("Cannot return stellar mass at time before formation"))

		return self.mass_form * (1 - f_loss(dt*10**9))





QUENCHING_TIMESCALE_DATA = mypy.readcat('../data/figure8_Fillingham+2017.dat')

class galaxy(object):
	'''
	Description:
	    This class is designed to encapsulate all of the 
	    relevant information of a single simulated galaxy.

	Attributes:
	    z_init                    # Redshift of initialization
	    lmass_init                # Seed stellar mass at z_init
	    lage_init                 # Log age of seed stellar mass at z_init
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
	def __init__(self, z_init, lmass_init, lage_init, quenched=False):

		self.z_init = z_init
		
		self.mass_form = 10**lmass_init / (1 - f_loss(10**lage_init))
		self.t_form = cosmo.lookback_time(z_init).value + 10**(lage_init-9)
		self.z_form = numpy.interp(self.t_form, t_master, z_master)

		self.seed_particle = stellar_mass_particle(self.t_form, self.mass_form)

		self.in_situ_particles = []
		self.ex_situ_particles = []

		self.quenched = numpy.array([quenched])

		if quenched:
			self.sfr = numpy.array([0.])
		else:
			self.sfr = numpy.array([sfrmass_relation_tomczak2016(z_init, lmass_init)])

		self.n_vminor_mergers = 0
		self.n_minor_mergers = 0
		self.n_major_mergers = 0
		self.mass_from_vminor_mergers = 0.
		self.mass_from_minor_mergers = 0.
		self.mass_from_major_mergers = 0.

		self.merged_in_this_timestep = False

		self.quenching_timescale = numpy.interp(lmass_init, QUENCHING_TIMESCALE_DATA.lmass, QUENCHING_TIMESCALE_DATA.qenching_timescale)
		self.quenching_countdown = -1
		while self.quenching_countdown < 0:
			self.quenching_countdown = numpy.random.randn() * 1. + self.quenching_timescale

	def get_lmass_total(self, t):
		'''
		Returns the total stellar mass at lookback time t (in Gyr)

		total  =  (seed mass) + (ex situ mass) + (in situ mass) - (mass losses)
		'''

		dt = self.seed_particle.t_form - t
		if dt < 0:
			raise(IOError("Cannot return stellar mass at time before formation"))


		###  Start with mass of seed particle
		mass_summation = self.seed_particle.mass_form * (1 - f_loss(dt*10**9))


		###  Next add mass particles from in situ SF
		for particle in self.in_situ_particles:
			mass_summation += particle.get_mass(t)


		###  Next add mass particles from ex situ mergers
		for particle in self.ex_situ_particles:
			mass_summation += particle.get_mass(t)


		return numpy.log10(mass_summation)




class mock_catalog(object):
	'''
	Description:
	    This class is designed to encapsulate all of the 
	    information of a toy model simulation from start 
	    to finish and every time-step inbetween.

	Attributes:
	    lmassbins          # Stellar mass bin boundaries for SMF
	    lmassbars          # Stellar mass bin centers for SMF

	    z_start   = 5.     # Redshift at which the simulation "begins"
	    z_final   = 0.8    # Redshift at which the simulation "ends"
	    timestep  = 100.   # Duration of each time-step in Myr

	    t_start            # Lookback-time at z_start in Gyr
	    t_final            # Lookback-time at z_final in Gyr
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

	    fquiescent_ORELSE_interp  # Quiescent fraction interpolated between the measured ORELSE high-denstiy and 0% fqu at z=z_start

	    SMFs               # Array contaning the numbers of gals per massbin for each time-step
	'''
	def __init__(self, lmassbins, lmassbars, z_start=5., z_final=0.8, timestep=100.):

		self.z_start  = z_start
		self.z_final  = z_final
		self.timestep = timestep

		self.t_start = cosmo.lookback_time(z_start).value
		self.t_final   = cosmo.lookback_time(z_final).value

		self.t_array = numpy.arange(self.t_start, self.t_final, -timestep/10.**3)
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
		self.ngal_timesteps = numpy.array([])
		self.mass_running    = numpy.zeros(self.n_timesteps)
		self.mass_in_mergers = numpy.zeros(self.n_timesteps)
		self.mass_in_ICL     = numpy.zeros(self.n_timesteps)

		self.fquiescent = []
		self.fquiescent_ORELSE_interp = [numpy.zeros(len(lmassbars))]

		self.lmassbins = lmassbins
		self.lmassbars = lmassbars
		self.SMFs = []

	def get_galaxy_lmasses(self, t):
		'''
		Returns an array of the masses of the current
		galaxies list at lookback time t (in Gyr)
		'''
		return numpy.array([g.get_lmass_total(t) for g in self.galaxies])

	def get_f_ICL(self):
		'''
		Returns an array of the fraction of total stellar mass in the ICL at each time-step
		'''
		return self.mass_in_ICL / self.mass_running

	def generate_SMF(self, t):
		'''
		Generates the current SMF at lookback time t (in Gyr)
		and appends it to the SMFs attribute
		'''
		lmasses = numpy.array(self.get_galaxy_lmasses(t))
		digi = numpy.digitize(lmasses, self.lmassbins)
		bincount = numpy.bincount(digi, minlength=len(self.lmassbins)+1)
		self.SMFs.append(bincount[1:-1])

	def generate_fquiescent(self, t):
		'''
		Generates the current fquiescent vs. stellar mass at lookback time t (in Gyr)
		'''

		self.fquiescent.append(numpy.zeros(len(self.lmassbars)))

		lmasses = numpy.array(self.get_galaxy_lmasses(t))
		digi = numpy.digitize(lmasses, self.lmassbins)

		quenched_array = numpy.array([g.quenched[-1] for g in self.galaxies])
		for i in range(1, len(self.lmassbins)):
			inds_lmassbin = numpy.where(digi == i)[0]
			inds_lmassbin_quenched = numpy.where(quenched_array[inds_lmassbin])[0]

			if len(inds_lmassbin) > 0:
				self.fquiescent[-1][i-1] = len(inds_lmassbin_quenched) * 1. / len(inds_lmassbin)




























#########################################################
###  Looping over fmerge (faction of galaxies to merge)
#########################################################


fmerge = 0.5





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
		#if numpy.random.rand() < fquiescent_leja2015(simulation.z_start, m_rand):
		#	q = True

		###  grabbing a random age, not to exceed the 
		###  age of the universe at z_start
		inds = numpy.where((lmass_xlss005 > m_rand-dm/2.) &
			               (lmass_xlss005 < m_rand+dm/2.) &
			               (lage_xlss005 < numpy.log10(age_at_z_start * 10**9)))[0]
		lage = lage_xlss005[inds][numpy.random.randint(0, len(inds))]


		simulation.galaxies.append(galaxy(simulation.z_start, m_rand, lage_init=lage, quenched=q))
		simulation.ngal_initial += 1

simulation.ngal_timesteps = numpy.append(simulation.ngal_timesteps, len(simulation.galaxies))

###  randomizing order of galaxies list
numpy.random.shuffle(simulation.galaxies)


###  storing initial fquiescent vs. stellar mass
simulation.generate_fquiescent(simulation.t_start)


###  storing initial SMF of simulation
simulation.generate_SMF(simulation.t_start)
simulation.mass_running[0] += (10**simulation.get_galaxy_lmasses(simulation.t_start)).sum()


###  setting the number of mergers that should occur at each time-step.
###    (1) Define such that at the end the number of galaxies
###        is 1% of the number at the simulation start.
###    (2) Distribute in a way to match the merger-rates
###        found in Illustris as a function of redshift

#simulation.nstep_mergers *= (simulation.ngal_initial * 0.995)
simulation.nstep_mergers *= (simulation.ngal_initial * fmerge)

simulation.nstep_mergers = simulation.nstep_mergers.astype(int)


tf = time.time()
dt = (tf-t0)
hours = int(dt / 3600)
minutes = int((dt % 3600) / 60)
seconds = int((dt % 3600) % 60)
print '\n\ndone!   %02ih %02im %02is\n' % (hours, minutes, seconds)









'''
###  Quenching all galaxies to test the scenario w/o star-formation
for g in simulation.galaxies:
	g.quenched[-1] = True
'''






#########################
###  Running simulation
#########################

print 'Running merger simulation...'
t0 = time.time()

for i_timestep in range(1, simulation.n_timesteps):

	t_step = simulation.t_array[i_timestep]


	mypy.progress_bar(i_timestep, simulation.n_timesteps)




	###########################################
	###  updating properties for all galaxies
	###########################################

	for i_galaxy in range(len(simulation.galaxies)):

		g = simulation.galaxies[i_galaxy]

		###  reset merger flag
		g.merged_in_this_timestep = False





		#############################################
		###  grow galaxy via in situ star-formation
		#############################################

		###  update SFR array if not already quenched
		g.quenched = numpy.append(g.quenched, g.quenched[-1])
		if not g.quenched[-1]:
			mass_growth = g.sfr[-1] * (simulation.timestep * 10**6)
			g.in_situ_particles.append(stellar_mass_particle(t_step, mass_growth))
			g.sfr = numpy.append(g.sfr, sfrmass_relation_tomczak2016(simulation.z_array[i_timestep], g.get_lmass_total(t_step)))
		else:
			g.sfr = numpy.append(g.sfr, 0.)





	'''
	#########################################
	###  Quenching galaxies:  Prescription 1
	###    Matching fquiescent of Leja+2015
	#########################################

	###  quenching galaxies at random to match the
	###  fquiescent vs. stellar mass from Leja+2015
	
	fquiescent_from_leja2015 = fquiescent_leja2015(simulation.z_array[i_timestep], simulation.lmassbars)
	d_fquiescent = fquiescent_from_leja2015 - simulation.fquiescent[-1]
	digi_lmass = numpy.digitize(simulation.get_galaxy_lmasses(t_step), simulation.lmassbins)

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
	'''



	'''
	###################################################
	###  Quenching galaxies:  Prescription 2
	###    Quenching countdown timer, Fillingham+2017
	###################################################

	dt = simulation.t_start - t_step
	for g in simulation.galaxies:

		if dt > g.quenching_countdown:
			g.quenched[-1] = True
	'''



	###################################################
	###  Quenching galaxies:  Prescription 3
	###    Matching interpolated fquiescent
	###    between ORELSE at z~0.8 and assuming
	###    0% fquiescent at z_start
	###################################################

	###          alpha  lmstar phistar
	params_sf = [-0.82, 10.83, 8.57*10**-3]
	params_qu = [-0.52, 11.04, 24.2*10**-3]

	phi_ORELSE_sf = schechter_mf(simulation.lmassbars, *params_sf)
	phi_ORELSE_qu = schechter_mf(simulation.lmassbars, *params_qu)
	fquiescent_ORELSE_15v20 = phi_ORELSE_qu / (phi_ORELSE_qu + phi_ORELSE_sf)


	dt_total = simulation.t_start - cosmo.lookback_time(0.8).value
	dt_tstep = simulation.t_start - t_step


	###  interpolating fquiescent_ORELSE between t_final and t_now
	fq_tstep = numpy.average([fquiescent_ORELSE_15v20, 
		                      numpy.zeros(len(simulation.lmassbars))],
		                      axis=0, weights=[dt_tstep / dt_total, 1 - dt_tstep / dt_total])

	simulation.fquiescent_ORELSE_interp.append(fq_tstep)


	###  calculating numbers of galaxies to quench to match fquiescent_ORELSE
	d_fquiescent = fq_tstep - simulation.fquiescent[-1]
	digi_lmass = numpy.digitize(simulation.get_galaxy_lmasses(t_step), simulation.lmassbins)

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
			lmasses = numpy.array([simulation.galaxies[ij_rand[0]].get_lmass_total(t_step), 
				                   simulation.galaxies[ij_rand[1]].get_lmass_total(t_step)])
			mratio = 10**(lmasses.min() - lmasses.max())

			###  definition of major merger 1:4
			if mratio > 0.25 and numpy.random.rand() < 1./3:
				ready_to_merge = True
			elif mratio < 0.25:
				ready_to_merge = True


		###  now we have selected our galaxies to merge
		i_more_massive = ij_rand[lmasses == lmasses.max()][0]
		i_less_massive = ij_rand[lmasses == lmasses.min()][0]

		gal0 = simulation.galaxies[i_more_massive]
		gal1 = simulation.galaxies[i_less_massive]



		###  stripping F_MASS_LOSS of less massive galaxy's particles to ICL
		simulation.mass_in_ICL[i_timestep:] += F_MASS_LOSS * gal1.seed_particle.get_mass(t_step)
		gal1.seed_particle.mass_form *= (1 - F_MASS_LOSS)

		for particle in gal1.in_situ_particles:
			simulation.mass_in_ICL[i_timestep:] += F_MASS_LOSS * particle.get_mass(t_step)
			particle.mass_form *= (1 - F_MASS_LOSS)

		for particle in gal1.ex_situ_particles:
			simulation.mass_in_ICL[i_timestep:] += F_MASS_LOSS * particle.get_mass(t_step)
			particle.mass_form *= (1 - F_MASS_LOSS)




		###  adding mass particles of less massive galaxy to more massive galaxy
		gal0.ex_situ_particles.append(gal1.seed_particle)
		gal0.ex_situ_particles += gal1.in_situ_particles
		gal0.ex_situ_particles += gal1.ex_situ_particles


		###  updating merger flag
		gal0.merged_in_this_timestep = True


		###  incrementing merger counter of galaxy
		if 0.01 <= mratio < 0.1:   gal0.n_vminor_mergers += 1
		elif 0.1 <= mratio < 0.25: gal0.n_minor_mergers += 1
		elif mratio >= 0.25:       gal0.n_major_mergers += 1


		###  popping merged galaxy
		pop = simulation.galaxies.pop(i_less_massive)


	###  done with merging
	simulation.ngal_timesteps = numpy.append(simulation.ngal_timesteps, len(simulation.galaxies))




	###  recording total running stellar mass
	simulation.mass_running[i_timestep] += (10**simulation.get_galaxy_lmasses(t_step)).sum()
	simulation.mass_running[i_timestep] += simulation.mass_in_ICL[i_timestep]


	###  storing SMF and fquiescent
	simulation.generate_SMF(t_step)
	simulation.generate_fquiescent(t_step)




tf = time.time()
dt = (tf-t0)
hours = int(dt / 3600)
minutes = int((dt % 3600) / 60)
seconds = int((dt % 3600) % 60)
print '\n\ndone!   %02ih %02im %02is\n' % (hours, minutes, seconds)









################################################################
###  Generats a diagnostic plot of an individual galaxy object
################################################################
def plot_diagnostic(g):
	'''
	Generats a diagnostic plot of an individual galaxy object
	'''

	fig = pyplot.figure(figsize=(13., 8.))
	sp1 = fig.add_subplot(211)
	sp2 = fig.add_subplot(212)

	fig.subplots_adjust(hspace=0)

	#sp1.grid()
	#sp2.grid()

	sp1.minorticks_on()
	sp2.minorticks_on()

	sp1.set_ylabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('time [Gyr]')
	sp2.set_ylabel('log( SFR  [M$_{\odot}$ / yr] )')

	dt = abs(simulation.t_array[0] - simulation.t_array[1])
	sp1.set_xlim((age_of_universe-simulation.t_array).min() - 3*dt, (age_of_universe-simulation.t_array).max() + 3*dt)
	sp2.set_xlim((age_of_universe-simulation.t_array).min() - 3*dt, (age_of_universe-simulation.t_array).max() + 3*dt)


	###  plotting masses
	sp1.step(age_of_universe-simulation.t_array, numpy.log10(g.mass_array), 
		     ls='-', lw=3, marker='', color='b', label='Total', where='mid')


	###  plotting SFRs
	inds = numpy.where(g.sfr > 0)[0]
	sp2.step(age_of_universe-simulation.t_array[inds], numpy.log10(g.sfr[inds]),
		     ls='-', lw=3, marker='', color='b', label='in situ', where='mid')


	###  plotting vlines at mergers
	#sp1. axvline(0, color='gray', lw=2, label='merger')
	#for i in range(len(g.t_array)):
	#	if g.mass_ex_situ[i] > 0:
	#		sp1.axvline(g.t_array[i], color='gray', lw=2)


	#   plotting vline at quenching
	for i in range(len(simulation.t_array)):
		if g.quenched[i]:
			sp1.axvline(age_of_universe-simulation.t_array[i], color='r', lw=2, label='quench')
			sp2.axvline(age_of_universe-simulation.t_array[i], color='r', lw=2, label='quench')

			sp1.plot(age_of_universe-simulation.t_array[i], numpy.log10(g.mass_array[i]),
				     marker='x', markersize=12, mew=3, color='r')
			break

	#for ti in g.t_vminor_mergers: sp1.axvline(age_of_universe-ti, lw=2, color='k', ls=':')
	#for ti in g.t_minor_mergers: sp1.axvline(age_of_universe-ti, lw=2, color='k', ls=':')
	#for ti in g.t_major_mergers: sp1.axvline(age_of_universe-ti, lw=2, color='k', ls='-')


	for ti in g.t_vminor_mergers + g.t_minor_mergers:
		mi = numpy.interp(ti, simulation.t_array[::-1], g.mass_array[::-1])
		sp1.plot(age_of_universe-ti, numpy.log10(mi), marker='s', mfc='none', mec='g', mew=3, markersize=9)

	for ti in g.t_major_mergers:
		mi = numpy.interp(ti, simulation.t_array[::-1], g.mass_array[::-1])
		sp1.plot(age_of_universe-ti, numpy.log10(mi), marker='o', mfc='none', mec='k', mew=3, markersize=12)




	return fig, sp1, sp2




























######################################
###  plotting best-fit models to data
######################################

if False:
	smfdata = mypy.readcat('../data/massFunctions_TOT_Tomczak+2017.dat')
	cmap = pyplot.cm.cool


	fig = pyplot.figure(figsize=(16.4, 6.6375))

	sp3 = fig.add_subplot(133)
	sp2 = fig.add_subplot(132)
	sp1 = fig.add_subplot(131)

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp3.minorticks_on()

	sp1.set_yscale('log')
	sp2.set_yscale('log')
	sp3.set_yscale('log')

	sp1.axis([8.4, 11.7, 8e-4, 3e-1])
	sp2.axis([8.4, 11.7, 8e-4, 3e-1])
	sp3.axis([8.4, 11.7, 8e-4, 3e-1])

	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('$\Phi$ / Mpc$^3$ / dex')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')

	fig.subplots_adjust(wspace=0, left=0.07)


	def plot_smf_with_fit(sp, phi, elo, ehi, color_i):


		inds = phi > 0

		x_data = smfdata.lmass[inds]
		y_data = phi[inds]
		elo_data = elo[inds]
		ehi_data = ehi[inds]

		sp.errorbar(x_data, y_data, xerr=0.25/2, yerr=[elo_data, ehi_data], 
			ls='', marker='o', mew=2, ms=11, mfc='none', mec=color_i, ecolor=color_i, zorder=10)




		###  fitting toy models
		scale_factors = pylab.zeros(simulation.n_timesteps)
		chi2s = pylab.zeros(simulation.n_timesteps)

		log_errors = (pylab.log10(y_data/(y_data-elo_data)) + pylab.log10((y_data+ehi_data)/y_data)) / 2.

		weights = 1. / ((elo_data + ehi_data)/2.)**2
		weights = 1. / log_errors**2


		bestfit_models = numpy.zeros((simulation.n_timesteps, len(x_data)))


		for i_model in range(simulation.n_timesteps):

			y_model = simulation.SMFs[i_model]
			y_model = numpy.interp(x_data, simulation.lmassbars, y_model)

			numer = pylab.sum(y_data  * weights)
			denom = pylab.sum(y_model * weights)

			scale = numer / denom
			scale_factors[i_model] = scale
			bestfit_models[i_model] = y_model * scale

			chi2s[i_model] = pylab.sum((pylab.log10(y_data) - pylab.log10(y_model * scale))**2 / log_errors**2)


		Pchi = 1. / chi2s**2
		Pchi = pylab.exp(-0.5 * chi2s)



		###  plotting best-fit model SMF
		i_best = pylab.find(Pchi == Pchi.max())[0]
		sp.plot(x_data, bestfit_models[i_best], color='r', lw=1.5)


		###  plotting marginalized best-fit SMF
		best_smf = pylab.average(bestfit_models, axis=0, weights=Pchi)

		sp.plot(x_data, best_smf, color='k', lw=5, zorder=1)
		sp.plot(x_data, best_smf, color='gray', lw=2, zorder=1)


		###  plotting stats
		best_f_icl = pylab.average(simulation.get_f_ICL(), axis=0, weights=Pchi)


		###  plotting info text
		t = 'f$_{ICL}$ = %.1f%%' % (best_f_icl*100)
		t += '\nN$_{initial}$ = %i' % simulation.ngal_initial
		t += '\nN$_{best-fit}$ = %i' % pylab.average(simulation.ngal_timesteps, axis=0, weights=Pchi)
		t += ' (%.1f%%)' % (100.*pylab.average(simulation.ngal_timesteps, axis=0, weights=Pchi)/simulation.ngal_initial)

		sp.text(0.03, 0.03, t,
			    horizontalalignment='left', verticalalignment='bottom', transform=sp.transAxes)





	plot_smf_with_fit(sp1, smfdata.phi_05v10, smfdata.elo_05v10, smfdata.ehi_05v10, cmap(2./4.))
	plot_smf_with_fit(sp2, smfdata.phi_10v15, smfdata.elo_10v15, smfdata.ehi_10v15, cmap(3./4.))
	plot_smf_with_fit(sp3, smfdata.phi_15v20, smfdata.elo_15v20, smfdata.ehi_15v20, cmap(4./4.))












########################################################
###  Plotting SFR-M* and fquiescent goals of simulation
########################################################

smfs_orelse_tot = mypy.readcat('../data/massFunctions_TOT_Tomczak+2017.dat')
smfs_orelse_qu  = mypy.readcat('../data/massFunctions_QU_Tomczak+2017.dat')
smfs_orelse_sf  = mypy.readcat('../data/massFunctions_SF_Tomczak+2017.dat')

inds_orelse = numpy.where(smfs_orelse_qu.phi_15v20 > 0)

if True:

    #fig = pyplot.figure(figsize=(16.3, 8.))
    fig = pyplot.figure(figsize=(16.3, 7.7))

    sp1 = fig.add_subplot(121)
    sp2 = fig.add_subplot(122)

    fig.subplots_adjust(left=0.08)

    sp1.minorticks_on()
    sp2.minorticks_on()

    sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
    sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')

    sp1.set_ylabel('log( SFR [M$_{\odot}$ / yr] )')
    sp2.set_ylabel('Quiescent Fraction')

    sp1.axis([simulation.lmassbars.min()-0.3,
              simulation.lmassbars.max()+0.3,
              -3., 3.])
    sp2.axis([simulation.lmassbars.min()-0.3,
              simulation.lmassbars.max()+0.3,
              -0.03, 1.03])

    cmap = pyplot.cm.Spectral

    #for i in range(simulation.n_timesteps):
    for i in [0, simulation.n_timesteps-1]:

        color_i = cmap(i / (simulation.n_timesteps - 1.))

        if i < 4:
            sfrs_i = sfrmass_relation_tomczak2016(simulation.z_array[4], simulation.lmassbars)
        else:
            sfrs_i = sfrmass_relation_tomczak2016(simulation.z_array[i], simulation.lmassbars)

        lw = 1.5
        lw = 3
        sp1.plot(simulation.lmassbars, numpy.log10(sfrs_i), color=color_i, lw=lw)
        sp2.plot(simulation.lmassbars, simulation.fquiescent_ORELSE_interp[i], color=color_i, lw=lw, ls='-')


    #dx, dy = 40, 5
    #sp_cbar = add_inset(sp1, rect=[0.5, 0.08, dx/100., dy/100.])
    #sp_cbar.set_title('redshift')
    #sp_cbar.set_yticks([])
    
    #xx, yy = numpy.meshgrid(numpy.linspace(simulation.z_start, simulation.z_final, dx), range(dy))
    #sp_cbar.imshow(xx, interpolation='nearest', cmap=cmap, extent=[xx[0][-1], xx[0][0], yy[0][0], yy[-1][0]])
    #sp_cbar.set_aspect('auto')
    #sp_cbar.set_xticks([5, 4, 3, 2, 1])


    ###  plotting ORELSE fqu
    #xdata = smfs_orelse_qu.lmass[inds_orelse]
    #ydata = smfs_orelse_qu.phi_15v20[inds_orelse] / smfs_orelse_tot.phi_15v20[inds_orelse]
    #ydata = smfs_orelse_qu.phi_15v20[inds_orelse] / (smfs_orelse_qu.phi_15v20[inds_orelse] + smfs_orelse_sf.phi_15v20[inds_orelse])
    #sp2.errorbar(xdata, ydata, xerr=0.25/2, yerr=0, ls='', marker='o', ms=11, mfc='none', mew=2, mec='k', ecolor='k')


    ###  adding white vlines to mimic dashed lines. This is for indicating the mass-completeness limits
    #for lmass_i in numpy.arange(9.5+0.25/2, 5, -0.06):
    #    sp2.axvline(lmass_i, color='w', lw=3)












######################################
###  Plotting final SMF of simulation
######################################
if True:

	t = 'Toy Model:'
	t += '\nz$_{start}$ = %.1f  -->  z$_{final}$ = %.1f' % (simulation.z_start, simulation.z_final)
	t += '\ntime-step = %i Myr (%i steps)' % (simulation.timestep, simulation.n_timesteps)
	t += '\nN$_{start}$ = %i' % simulation.ngal_initial
	t += '\nN$_{final}$ = %i (%.1f%% merged)' % (simulation.ngal_timesteps[-1], 100.*(1-simulation.ngal_timesteps[-1]/simulation.ngal_initial))


	fig = pyplot.figure(figsize=(14., 11.3))
	sp1 = fig.add_subplot(111)

	fig.subplots_adjust(left=0.08, bottom=0.07)

	sp1.minorticks_on()
	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number')

	sp1.set_yscale('log')
	sp1.grid()

	sp1.set_xlim([simulation.lmassbars.min()-0.3,
		          simulation.lmassbars.max()+0.3])


	sp1.plot(simulation.lmassbars, simulation.SMFs[0], color='gray', ls='--', lw=4, label='z$_{start}$')
	sp1.plot(simulation.lmassbars, simulation.SMFs[-1], color='k', lw=4, label='z$_{final}$')


	#sp1.text(0.03, 0.03, t, transform=sp1.transAxes, bbox={'facecolor':'w'})
	sp1.legend(loc=3, title=t, ncol=2)


	fmerged = 100. * (1-simulation.ngal_timesteps[-1] / simulation.ngal_initial)
	fig.savefig('../output/simulations/smf_start2final_fmerged%02i.pdf' % fmerged)
	pickle.dump(simulation, open('../output/simulations/simulation_merged%02i.pickle' % fmerged, 'wb'))
	pyplot.close()




















