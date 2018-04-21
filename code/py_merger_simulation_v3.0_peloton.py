
import os
import sys
import time
import numpy
import pickle
import subprocess
from scipy import optimize
from matplotlib import pyplot
from collections import namedtuple
from matplotlib.font_manager import FontProperties

######################################
###  Setting pyplot rc parameters  ###
######################################
if True:
	pyplot.ioff()
	pyplot.rcParams[u'agg.path.chunksize'] = 0
	pyplot.rcParams[u'animation.avconv_args'] = []
	pyplot.rcParams[u'animation.avconv_path'] = u'avconv'
	pyplot.rcParams[u'animation.bitrate'] = -1
	pyplot.rcParams[u'animation.codec'] = u'mpeg4'
	pyplot.rcParams[u'animation.convert_args'] = []
	pyplot.rcParams[u'animation.convert_path'] = u'convert'
	pyplot.rcParams[u'animation.ffmpeg_args'] = []
	pyplot.rcParams[u'animation.ffmpeg_path'] = u'ffmpeg'
	pyplot.rcParams[u'animation.frame_format'] = u'png'
	pyplot.rcParams[u'animation.mencoder_args'] = []
	pyplot.rcParams[u'animation.mencoder_path'] = u'mencoder'
	pyplot.rcParams[u'animation.writer'] = u'ffmpeg'
	pyplot.rcParams[u'axes.axisbelow'] = False
	pyplot.rcParams[u'axes.color_cycle'] = [u'b', u'g', u'r', u'c', u'm', u'y', u'k']
	pyplot.rcParams[u'axes.edgecolor'] = u'k'
	pyplot.rcParams[u'axes.facecolor'] = u'w'
	pyplot.rcParams[u'axes.formatter.limits'] = [-7, 7]
	pyplot.rcParams[u'axes.formatter.use_locale'] = False
	pyplot.rcParams[u'axes.formatter.use_mathtext'] = False
	#pyplot.rcParams[u'axes.formatter.useoffset'] = True
	pyplot.rcParams[u'axes.grid'] = False
	#pyplot.rcParams[u'axes.grid.which'] = u'major'
	pyplot.rcParams[u'axes.hold'] = True
	pyplot.rcParams[u'axes.labelcolor'] = u'k'
	pyplot.rcParams[u'axes.labelsize'] = u'medium'
	pyplot.rcParams[u'axes.labelweight'] = u'normal'
	pyplot.rcParams[u'axes.linewidth'] = 2.5
	pyplot.rcParams[u'axes.titlesize'] = u'large'
	#pyplot.rcParams[u'axes.titleweight'] = u'normal'
	pyplot.rcParams[u'axes.unicode_minus'] = True
	pyplot.rcParams[u'axes.xmargin'] = 0.0
	pyplot.rcParams[u'axes.ymargin'] = 0.0
	pyplot.rcParams[u'axes3d.grid'] = True
	pyplot.rcParams[u'backend'] = u'TkAgg'
	pyplot.rcParams[u'backend.qt4'] = u'PyQt4'
	#pyplot.rcParams[u'backend.qt5'] = u'PyQt5'
	pyplot.rcParams[u'backend_fallback'] = True
	pyplot.rcParams[u'contour.negative_linestyle'] = u'dashed'
	pyplot.rcParams[u'docstring.hardcopy'] = False
	pyplot.rcParams[u'examples.directory'] = u''
	pyplot.rcParams[u'figure.autolayout'] = False
	pyplot.rcParams[u'figure.dpi'] = 80.0
	pyplot.rcParams[u'figure.edgecolor'] = u'w'
	pyplot.rcParams[u'figure.facecolor'] = u'0.75'
	pyplot.rcParams[u'figure.figsize'] = [8.0, 6.0]
	pyplot.rcParams[u'figure.frameon'] = True
	pyplot.rcParams[u'figure.max_open_warning'] = 20
	pyplot.rcParams[u'figure.subplot.bottom'] = 0.12
	pyplot.rcParams[u'figure.subplot.hspace'] = 0.2
	pyplot.rcParams[u'figure.subplot.left'] = 0.15
	pyplot.rcParams[u'figure.subplot.right'] = 0.97
	pyplot.rcParams[u'figure.subplot.top'] = 0.97
	pyplot.rcParams[u'figure.subplot.wspace'] = 0.2
	pyplot.rcParams[u'font.cursive'] = [u'Apple Chancery', u'Textile', u'Zapf Chancery', u'Sand', u'cursive']
	pyplot.rcParams[u'font.family'] = [u'serif']
	pyplot.rcParams[u'font.fantasy'] = [u'Comic Sans MS', u'Chicago', u'Charcoal', u'ImpactWestern', u'fantasy']
	pyplot.rcParams[u'font.monospace'] = [u'Bitstream Vera Sans Mono', u'DejaVu Sans Mono', u'Andale Mono', u'Nimbus Mono L', u'Courier New', u'Courier', u'Fixed', u'Terminal', u'monospace']
	pyplot.rcParams[u'font.sans-serif'] = [u'Bitstream Vera Sans', u'DejaVu Sans', u'Lucida Grande', u'Verdana', u'Geneva', u'Lucid', u'Arial', u'Helvetica', u'Avant Garde', u'sans-serif']
	pyplot.rcParams[u'font.serif'] = [u'Bitstream Vera Serif', u'DejaVu Serif', u'New Century Schoolbook', u'Century Schoolbook L', u'Utopia', u'ITC Bookman', u'Bookman', u'Nimbus Roman No9 L', u'Times New Roman', u'Times', u'Palatino', u'Charter', u'serif']
	pyplot.rcParams[u'font.size'] = 18.0
	pyplot.rcParams[u'font.stretch'] = u'normal'
	pyplot.rcParams[u'font.style'] = u'normal'
	pyplot.rcParams[u'font.variant'] = u'normal'
	pyplot.rcParams[u'font.weight'] = u'normal'
	pyplot.rcParams[u'grid.alpha'] = 1.0
	pyplot.rcParams[u'grid.color'] = u'k'
	pyplot.rcParams[u'grid.linestyle'] = u':'
	pyplot.rcParams[u'grid.linewidth'] = 0.5
	pyplot.rcParams[u'image.aspect'] = u'equal'
	pyplot.rcParams[u'image.cmap'] = u'BuGn'
	pyplot.rcParams[u'image.interpolation'] = u'bilinear'
	pyplot.rcParams[u'image.lut'] = 256
	pyplot.rcParams[u'image.origin'] = u'lower'
	pyplot.rcParams[u'image.resample'] = False
	pyplot.rcParams[u'interactive'] = True
	pyplot.rcParams[u'keymap.all_axes'] = [u'a']
	pyplot.rcParams[u'keymap.back'] = [u'left', u'c', u'backspace']
	pyplot.rcParams[u'keymap.forward'] = [u'right', u'v']
	pyplot.rcParams[u'keymap.fullscreen'] = [u'f', u'ctrl+f']
	pyplot.rcParams[u'keymap.grid'] = [u'g']
	pyplot.rcParams[u'keymap.home'] = [u'h', u'r', u'home']
	pyplot.rcParams[u'keymap.pan'] = [u'p']
	pyplot.rcParams[u'keymap.quit'] = [u'ctrl+w', u'cmd+w']
	pyplot.rcParams[u'keymap.save'] = [u's', u'ctrl+s']
	pyplot.rcParams[u'keymap.xscale'] = [u'k', u'L']
	pyplot.rcParams[u'keymap.yscale'] = [u'l']
	pyplot.rcParams[u'keymap.zoom'] = [u'o']
	pyplot.rcParams[u'legend.borderaxespad'] = 0.5
	pyplot.rcParams[u'legend.borderpad'] = 0.4
	pyplot.rcParams[u'legend.columnspacing'] = 2.0
	pyplot.rcParams[u'legend.fancybox'] = False
	pyplot.rcParams[u'legend.fontsize'] = u'large'
	#pyplot.rcParams[u'legend.framealpha'] = 1.0
	pyplot.rcParams[u'legend.frameon'] = True
	pyplot.rcParams[u'legend.handleheight'] = 0.7
	pyplot.rcParams[u'legend.handlelength'] = 2.5
	pyplot.rcParams[u'legend.handletextpad'] = 0.8
	pyplot.rcParams[u'legend.isaxes'] = True
	pyplot.rcParams[u'legend.labelspacing'] = 0.5
	pyplot.rcParams[u'legend.loc'] = u'upper right'
	pyplot.rcParams[u'legend.markerscale'] = 1.0
	pyplot.rcParams[u'legend.numpoints'] = 2
	pyplot.rcParams[u'legend.scatterpoints'] = 1
	pyplot.rcParams[u'legend.shadow'] = True
	pyplot.rcParams[u'lines.antialiased'] = True
	pyplot.rcParams[u'lines.color'] = u'b'
	pyplot.rcParams[u'lines.dash_capstyle'] = u'butt'
	pyplot.rcParams[u'lines.dash_joinstyle'] = u'round'
	pyplot.rcParams[u'lines.linestyle'] = u'-'
	pyplot.rcParams[u'lines.linewidth'] = 2.0
	pyplot.rcParams[u'lines.marker'] = u'None'
	pyplot.rcParams[u'lines.markeredgewidth'] = 1.5
	pyplot.rcParams[u'lines.markersize'] = 6.0
	pyplot.rcParams[u'lines.solid_capstyle'] = u'projecting'
	pyplot.rcParams[u'lines.solid_joinstyle'] = u'round'
	pyplot.rcParams[u'mathtext.bf'] = u'serif:bold'
	pyplot.rcParams[u'mathtext.cal'] = u'cursive'
	pyplot.rcParams[u'mathtext.default'] = u'regular'
	pyplot.rcParams[u'mathtext.fallback_to_cm'] = True
	pyplot.rcParams[u'mathtext.fontset'] = u'cm'
	pyplot.rcParams[u'mathtext.it'] = u'serif:italic'
	pyplot.rcParams[u'mathtext.rm'] = u'serif'
	pyplot.rcParams[u'mathtext.sf'] = u'sans\\-serif'
	pyplot.rcParams[u'mathtext.tt'] = u'monospace'
	#pyplot.rcParams[u'nbagg.transparent'] = True
	pyplot.rcParams[u'patch.antialiased'] = True
	pyplot.rcParams[u'patch.edgecolor'] = u'k'
	pyplot.rcParams[u'patch.facecolor'] = u'b'
	pyplot.rcParams[u'patch.linewidth'] = 1.0
	pyplot.rcParams[u'path.effects'] = []
	pyplot.rcParams[u'path.simplify'] = True
	pyplot.rcParams[u'path.simplify_threshold'] = 0.1111111111111111
	pyplot.rcParams[u'path.sketch'] = None
	pyplot.rcParams[u'path.snap'] = True
	pyplot.rcParams[u'pdf.compression'] = 6
	pyplot.rcParams[u'pdf.fonttype'] = 42
	pyplot.rcParams[u'pdf.inheritcolor'] = False
	pyplot.rcParams[u'pdf.use14corefonts'] = False
	pyplot.rcParams[u'pgf.debug'] = False
	pyplot.rcParams[u'pgf.preamble'] = []
	pyplot.rcParams[u'pgf.rcfonts'] = True
	pyplot.rcParams[u'pgf.texsystem'] = u'xelatex'
	pyplot.rcParams[u'plugins.directory'] = u'.matplotlib_plugins'
	pyplot.rcParams[u'polaraxes.grid'] = True
	pyplot.rcParams[u'ps.distiller.res'] = 6000
	pyplot.rcParams[u'ps.fonttype'] = 42
	pyplot.rcParams[u'ps.papersize'] = u'letter'
	pyplot.rcParams[u'ps.useafm'] = False
	pyplot.rcParams[u'ps.usedistiller'] = False
	pyplot.rcParams[u'savefig.bbox'] = None
	pyplot.rcParams[u'savefig.directory'] = u'~'
	pyplot.rcParams[u'savefig.dpi'] = 100.0
	pyplot.rcParams[u'savefig.edgecolor'] = u'w'
	pyplot.rcParams[u'savefig.facecolor'] = u'w'
	pyplot.rcParams[u'savefig.format'] = u'png'
	pyplot.rcParams[u'savefig.frameon'] = True
	pyplot.rcParams[u'savefig.jpeg_quality'] = 95
	pyplot.rcParams[u'savefig.orientation'] = u'portrait'
	pyplot.rcParams[u'savefig.pad_inches'] = 0.1
	#pyplot.rcParams[u'savefig.transparent'] = False
	pyplot.rcParams[u'svg.fonttype'] = u'path'
	pyplot.rcParams[u'svg.image_inline'] = True
	pyplot.rcParams[u'svg.image_noscale'] = False
	pyplot.rcParams[u'text.antialiased'] = True
	pyplot.rcParams[u'text.color'] = u'k'
	#pyplot.rcParams[u'text.dvipnghack'] = None
	pyplot.rcParams[u'text.hinting'] = True
	pyplot.rcParams[u'text.hinting_factor'] = 8
	pyplot.rcParams[u'text.latex.preamble'] = []
	pyplot.rcParams[u'text.latex.preview'] = False
	pyplot.rcParams[u'text.latex.unicode'] = False
	pyplot.rcParams[u'text.usetex'] = False
	pyplot.rcParams[u'timezone'] = u'UTC'
	pyplot.rcParams[u'tk.window_focus'] = False
	pyplot.rcParams[u'toolbar'] = u'toolbar2'
	pyplot.rcParams[u'verbose.fileo'] = u'sys.stdout'
	pyplot.rcParams[u'verbose.level'] = u'silent'
	pyplot.rcParams[u'webagg.open_in_browser'] = True
	pyplot.rcParams[u'webagg.port'] = 8988
	pyplot.rcParams[u'webagg.port_retries'] = 50
	pyplot.rcParams[u'xtick.color'] = u'k'
	pyplot.rcParams[u'xtick.direction'] = u'in'
	pyplot.rcParams[u'xtick.labelsize'] = u'medium'
	pyplot.rcParams[u'xtick.major.pad'] = 4.0
	pyplot.rcParams[u'xtick.major.size'] = 7.0
	pyplot.rcParams[u'xtick.major.width'] = 1.5
	pyplot.rcParams[u'xtick.minor.pad'] = 4.0
	pyplot.rcParams[u'xtick.minor.size'] = 4.0
	pyplot.rcParams[u'xtick.minor.width'] = 1.5
	pyplot.rcParams[u'ytick.color'] = u'k'
	pyplot.rcParams[u'ytick.direction'] = u'in'
	pyplot.rcParams[u'ytick.labelsize'] = u'medium'
	pyplot.rcParams[u'ytick.major.pad'] = 4.0
	pyplot.rcParams[u'ytick.major.size'] = 7.0
	pyplot.rcParams[u'ytick.major.width'] = 1.5
	pyplot.rcParams[u'ytick.minor.pad'] = 4.0
	pyplot.rcParams[u'ytick.minor.size'] = 4.0
	pyplot.rcParams[u'ytick.minor.width'] = 1.5






sim_number = int(sys.argv[1])     # id number for this particular simulation run
fmerge     = float(sys.argv[2])   # fraction of galaxies to merger between z_start and z_final

print '\nsim_number = %i'   % sim_number
print   'fmerge     = %.2f' % fmerge
print ''


#####################################
###  Defining functions and such  ###
#####################################
if True:

	cosmo_data = numpy.loadtxt('/home/atomczak/DATA/cosmology_FlatLambdaCDM_H0_70_Om0_03.dat')
	z_master = cosmo_data[:,0]
	t_master = cosmo_data[:,1]
	age_of_universe = 13.4617




	F_MASS_LOSS = 0.3  # Fraction of stellar mass to loss per merger event


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
	    catalog = custom_catalog(*numpy.loadtxt(name, dtype=dtype, delimiter=delimiter, unpack=1))

	    ###  re-gzip if necessary
	    if gzipped:
	        subprocess.call('gzip %s' % name, shell=1)

	    return catalog

	def line(x, m, b): return m * x + b

	def schechter_mf(xaxis, alpha, mstar, phistar):
	    '''
	    #  DESCRIPTION:
	    #    Returns the values for a Schechter mass function
	    #    from a given mass-axis and Schechter parameters.
	    #
	    #  INPUTS:
	    #    xaxis = input mass value(s)
	    #    alpha = Schechter parameters
	    #    mstar = Schechter parameters
	    #    phistar = Schechter parameters    
	    '''
	    return numpy.log(10) * phistar * 10**((xaxis-mstar)*(1+alpha)) * numpy.exp(-10**(xaxis-mstar))

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
	data_dir = '/home/atomczak/DATA/ORELSE/catalogs'
	version = 'xlss005_v0.0.3'
	cat_xlss005 = readcat('%s/%s/%s.cat' % (data_dir, version, version))
	fout_xlss005 = readcat('%s/%s/%s_ZFparams.fout' % (data_dir, version, version))

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





	QUENCHING_TIMESCALE_DATA = readcat('/home/atomczak/PROJECTS/MASSFUNCTIONS_ORELSE/data/figure8_Fillingham+2017.dat')

	class galaxy(object):
		'''
		Description:
		    This class is designed to encapsulate all of the 
		    relevant information of a single simulated galaxy.

		Attributes:
		    z_init                    # Redshift of initialization
		    lmass_init                # Seed stellar mass at z_init
		    lage_init                 # Log age of seed stellar mass at z_init
		    mass_form                 # Stellar mass at formation inferred from lage0
		    mass_array                # Array of the running stellar mass at each time-step
		    quenched                  # Boolean for if galaxy is star-forming or quenched
		    sfr                       # Current star-formation rate in Msol/yr
		    n_vminor_mergers          # Total number of very minor mergers experienced
		    n_minor_mergers           # Total number of minor mergers experienced
		    n_major_mergers           # Total number of major mergers experienced
		    t_vminor_mergers          # Lookback time(s) at which very minor mergers occurred in Gyr
		    t_minor_mergers           # Lookback time(s) at which minor mergers occurred in Gyr
		    t_major_mergers           # Lookback time(s) at which major mergers occurred in Gyr
		    merged_in_this_timestep   # A flag that indicates if a galaxy was involved in a merger
		    quenching_timescale       # Quenching time-scale defined at galaxy seed mass in Gyr
		    quenching_countdown       # Randomly generated time interval after which this galaxy will quench SF in Gyr
		'''
		def __init__(self, id, z_init, lmass_init, lage_init, quenched=False):

			self.id = id
			self.z_init = z_init
			
			self.mass_form = 10**lmass_init / (1 - f_loss(10**lage_init))
			self.mass_array = numpy.array([10**lmass_init])
			self.t_form = numpy.interp(z_init, z_master, t_master) + 10**(lage_init-9)
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
			self.n_minor_mergers  = 0
			self.n_major_mergers  = 0

			self.t_vminor_mergers = []
			self.t_minor_mergers  = []
			self.t_major_mergers  = []

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




	def get_galaxy_lmasses(galaxies_list, t):
		'''
		Returns an array of the masses of the provided
		galaxies list at lookback time t (in Gyr)
		'''
		return numpy.array([g.get_lmass_total(t) for g in galaxies_list])



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

		    ngal_initial       # Total number of galaxies from simulation start
		    ngal_timesteps     # Total number of galaxies at each time-step
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

			self.t_start = numpy.interp(self.z_start, z_master, t_master)
			self.t_final = numpy.interp(self.z_final, z_master, t_master)

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


			###  enforcing constant merger rate
			#self.nstep_mergers = numpy.ones(self.n_timesteps)

			###  normalizing
			self.nstep_mergers /= self.nstep_mergers.sum()


			self.ngal_initial = 0
			self.ngal_timesteps = numpy.array([0])
			self.mass_running    = numpy.zeros(self.n_timesteps)
			self.mass_in_mergers = numpy.zeros(self.n_timesteps)
			self.mass_in_ICL     = numpy.zeros(self.n_timesteps)

			self.fquiescent = []
			self.fquiescent_ORELSE_interp = [numpy.zeros(len(lmassbars))]

			self.lmassbins = lmassbins
			self.lmassbars = lmassbars
			self.SMFs = []


		def get_f_ICL(self):
			'''
			Returns an array of the fraction of total stellar mass in the ICL at each time-step
			'''
			return self.mass_in_ICL / self.mass_running

		def generate_SMF(self, galaxies_list, t):
			'''
			Generates the current SMF at lookback time t (in Gyr)
			and appends it to the SMFs attribute
			'''
			lmasses = numpy.array(get_galaxy_lmasses(galaxies_list, t))
			digi = numpy.digitize(lmasses, self.lmassbins)
			bincount = numpy.bincount(digi, minlength=len(self.lmassbins)+1)
			self.SMFs.append(bincount[1:-1])

		def generate_fquiescent(self, galaxies_list, t):
			'''
			Generates the current fquiescent vs. stellar mass at lookback time t (in Gyr)
			'''

			self.fquiescent.append(numpy.zeros(len(self.lmassbars)))

			lmasses = numpy.array(get_galaxy_lmasses(galaxies_list, t))
			digi = numpy.digitize(lmasses, self.lmassbins)

			quenched_array = numpy.array([g.quenched[-1] for g in galaxies_list])
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
lmassbins = numpy.arange(6.-dm/2, 11.55+dm/2., dm)
lmassbars = (lmassbins[:-1] + lmassbins[1:]) / 2.
phi_n05v00 = dschechter(lmassbars, *params_n05v00)
phi_00v05  = dschechter(lmassbars, *params_00v05)
phi_n05v05 = (phi_n05v00 + phi_00v05) / 2.
phi_n05v05_normed = (1 * phi_n05v05 / phi_n05v05.min()).astype(int)



####################################
###  Generating simulation object(s)
####################################

print 'Initializing mock catalog with N0=%i galaxies...' % phi_n05v05_normed.sum()
t0 = time.time()

galaxies_list = []
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


		galaxies_list.append(galaxy(simulation.ngal_initial, simulation.z_start, m_rand, lage_init=lage, quenched=q))
		simulation.ngal_initial += 1
		simulation.ngal_timesteps[0] += 1


###  randomizing order of galaxies list
numpy.random.shuffle(galaxies_list)


###  storing initial fquiescent vs. stellar mass
simulation.generate_fquiescent(galaxies_list, simulation.t_start)


###  storing initial SMF of simulation
simulation.generate_SMF(galaxies_list, simulation.t_start)
simulation.mass_running[0] += (10**get_galaxy_lmasses(galaxies_list, simulation.t_start)).sum()


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
print 'done!   %02ih %02im %02is\n' % (hours, minutes, seconds)









'''
###  Quenching all galaxies to test the scenario w/o star-formation
for g in galaxies_list:
	g.quenched[-1] = True
'''






#########################
###  Running simulation
#########################

print 'Running merger simulation...'
t0 = time.time()

for i_timestep in range(1, simulation.n_timesteps):

	t_step = simulation.t_array[i_timestep]





	###########################################
	###  updating properties for all galaxies
	###########################################

	for i_galaxy in range(len(galaxies_list)):

		g = galaxies_list[i_galaxy]

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

		g.mass_array = numpy.append(g.mass_array, 10**g.get_lmass_total(t_step))



	'''
	#########################################
	###  Quenching galaxies:  Prescription 1
	###    Matching fquiescent of Leja+2015
	#########################################

	###  quenching galaxies at random to match the
	###  fquiescent vs. stellar mass from Leja+2015
	
	fquiescent_from_leja2015 = fquiescent_leja2015(simulation.z_array[i_timestep], simulation.lmassbars)
	d_fquiescent = fquiescent_from_leja2015 - simulation.fquiescent[-1]
	digi_lmass = numpy.digitize(get_galaxy_lmasses(galaxies_list, t_step), simulation.lmassbins)

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
				if not galaxies_list[i_galaxy].quenched[-1]:
					galaxies_list[i_galaxy].quenched[-1] = True
					n2quench -= 1
	'''



	'''
	###################################################
	###  Quenching galaxies:  Prescription 2
	###    Quenching countdown timer, Fillingham+2017
	###################################################

	dt = simulation.t_start - t_step
	for g in galaxies_list:

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


	###  calculating for low overdensity bin
	#params_sf = [10.76, -0.21, -1.52, 0.17*10**-3, 0.17*10**-3]
	#params_qu = [10.91, -0.71, -2.78, 0.31*10**-3, 0.0002*10**-3]
	#phi_ORELSE_sf = dschechter(simulation.lmassbars, *params_sf)
	#phi_ORELSE_qu = dschechter(simulation.lmassbars, *params_qu)

	params_sf = [-1.48, 11.06, 0.13*10**-3]
	params_qu = [-0.78, 10.93, 0.29*10**-3]

	phi_ORELSE_sf = schechter_mf(simulation.lmassbars, *params_sf)
	phi_ORELSE_qu = schechter_mf(simulation.lmassbars, *params_qu)
	fquiescent_ORELSE_00v05 = phi_ORELSE_qu / (phi_ORELSE_qu + phi_ORELSE_sf)



	dt_total = simulation.t_start - numpy.interp(0.8, z_master, t_master)
	dt_tstep = simulation.t_start - t_step


	###  interpolating fquiescent_ORELSE between t_final and t_now
	fq_tstep = numpy.average([fquiescent_ORELSE_00v05, 
		                      numpy.zeros(len(simulation.lmassbars))],
		                      axis=0, weights=[dt_tstep / dt_total, 1 - dt_tstep / dt_total])

	simulation.fquiescent_ORELSE_interp.append(fq_tstep)


	###  calculating numbers of galaxies to quench to match fquiescent_ORELSE
	d_fquiescent = fq_tstep - simulation.fquiescent[-1]
	digi_lmass = numpy.digitize(get_galaxy_lmasses(galaxies_list, t_step), simulation.lmassbins)

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
				if not galaxies_list[i_galaxy].quenched[-1]:
					galaxies_list[i_galaxy].quenched[-1] = True
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
			ij_rand = numpy.random.randint(0, len(galaxies_list), 2)

			###  ... but not the same galaxy twice
			while ij_rand[0] == ij_rand[1]:
				ij_rand = numpy.random.randint(0, len(galaxies_list), 2)

			###  ... and not if either galaxy was already involved in a merger in this time-step
			while galaxies_list[ij_rand[0]].merged_in_this_timestep or galaxies_list[ij_rand[1]].merged_in_this_timestep:
				ij_rand = numpy.random.randint(0, len(galaxies_list), 2)



			###  checking merger probability
			lmasses = numpy.array([galaxies_list[ij_rand[0]].get_lmass_total(t_step), 
				                   galaxies_list[ij_rand[1]].get_lmass_total(t_step)])
			mratio = 10**(lmasses.min() - lmasses.max())

			###  definition of major merger 1:4
			if mratio > 0.25 and numpy.random.rand() < 1./3:
				ready_to_merge = True
			elif mratio < 0.25:
				ready_to_merge = True


		###  now we have selected our galaxies to merge
		i_more_massive = ij_rand[lmasses == lmasses.max()][0]
		i_less_massive = ij_rand[lmasses == lmasses.min()][0]

		gal0 = galaxies_list[i_more_massive]
		gal1 = galaxies_list[i_less_massive]



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
		if 0.01 <= mratio < 0.1:
			gal0.n_vminor_mergers += 1
			gal0.t_vminor_mergers.append(t_step)
		elif 0.1 <= mratio < 0.25:
			gal0.n_minor_mergers += 1
			gal0.t_minor_mergers.append(t_step)
		elif mratio >= 0.25:
			gal0.n_major_mergers += 1
			gal0.t_major_mergers.append(t_step)


		###  popping merged galaxy
		pop = galaxies_list.pop(i_less_massive)


	###  done with merging
	simulation.ngal_timesteps = numpy.append(simulation.ngal_timesteps, len(galaxies_list))




	###  recording total running stellar mass
	simulation.mass_running[i_timestep] += (10**get_galaxy_lmasses(galaxies_list, t_step)).sum()
	simulation.mass_running[i_timestep] += simulation.mass_in_ICL[i_timestep]


	###  storing SMF and fquiescent
	simulation.generate_SMF(galaxies_list, t_step)
	simulation.generate_fquiescent(galaxies_list, t_step)



tf = time.time()
dt = (tf-t0)
hours = int(dt / 3600)
minutes = int((dt % 3600) / 60)
seconds = int((dt % 3600) % 60)
print 'done!   %02ih %02im %02is\n' % (hours, minutes, seconds)





###################################################################################
###  Saving pickles of mock_catalog and each galaxy in the final galaxies list  ###
###################################################################################

print 'saving pickles...'
t0 = time.time()

outdir = '/home/atomczak/PROJECTS/MASSFUNCTIONS_ORELSE/output/simulation_fmerge%02i_fquiescentField' % (100*fmerge)

try:
	os.mkdir(outdir)
	#os.mkdir(outdir + '/galaxies')
except:
	pass



###  save pickle of mock catalog
pickle.dump(simulation, open('%s/simulation_merged%02i_iter%02i.pickle' % (outdir, 100*fmerge, sim_number), 'wb'))



'''
###  save pickles of 1 random galaxy from each massbin
digi_lmass = numpy.digitize(get_galaxy_lmasses(galaxies_list, simulation.t_final), simulation.lmassbins)
for i_lmassbin in range(1, len(simulation.lmassbins)):

	inds_galaxies = numpy.where(digi_lmass == i_lmassbin)[0]

	###  make sure there is at least one galaxy in this massbin, if not then pass to next massbin
	if len(inds_galaxies) > 0:
		g_rand = numpy.random.choice(galaxies_list)
		pickle.dump(g_rand, open('%s/galaxies/galaxy%07i_iter%02i.pickle' % (outdir, g_rand.id, sim_number), 'wb'))
'''



tf = time.time()
dt = (tf-t0)
hours = int(dt / 3600)
minutes = int((dt % 3600) / 60)
seconds = int((dt % 3600) % 60)
print 'done!   %02ih %02im %02is\n' % (hours, minutes, seconds)













