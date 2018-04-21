
import mypy
import pylab
import pickle
import my_classes
from scipy import optimize
from matplotlib.font_manager import FontProperties


def line(x, m, b): return m * x + b

def dschechter(lmassax, lmstar, a1, a2, phistar1, phistar2):
	factor1 = pylab.log(10) * pylab.exp(-10**(lmassax - lmstar)) * 10**(lmassax - lmstar)
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



table11_leja2015 = mypy.readcat('../data/table11_Leja+2015.dat')
def get_merge_prob(lmass):
	return pylab.interp(lmass, table11_leja2015.lmass, 
		                       table11_leja2015.f_destroyed_Gyr / table11_leja2015.f_destroyed_Gyr.max())

class galaxy:
	def __init__(self, lmass):
		self.lmass = lmass
		self.lmass0 = lmass
		self.n_vminor_mergers = 0
		self.n_minor_mergers = 0
		self.n_major_mergers = 0
		self.mass_from_vminor_mergers = 0.
		self.mass_from_minor_mergers = 0.
		self.mass_from_major_mergers = 0.

class mock_catalog_class:
	def __init__(self):
		self.galaxies = []
		self.n_initial = 0.
		self.mass_initial = 0.
		self.mass_in_mergers = 0.
		self.mass_in_icm = 0.
		self.n_initial_gt_9 = 0
		self.n_merged_gt_9 = 0

class mock_catalog_class2:
	def __init__(self, n_initial):
		self.lmass = pylab.zeros(n_initial)
		self.n_initial = n_initial
		self.n_minor_mergers = 0.
		self.n_major_mergers = 0.


###  Total SMF of overdensity bin 0<log(1+d)<0.5
lmass = pylab.array([8.75, 9., 9.25, 9.5, 9.75, 10., 10.25, 10.5, 10.75, 11., 11.25, 11.5, 11.75])
dm = lmass[1] - lmass[0]
phi = pylab.array([1.9662e-02, 1.5473e-02, 1.1043e-02, 7.7730e-03, 6.2030e-03, 4.7015e-03,
	               3.6567e-03, 3.2601e-03, 2.3632e-03, 1.4071e-03, 6.0336e-04, 1.2318e-04, 2.9229e-05])


lmass = pylab.array([7., 7.25, 7.5, 7.75, 8., 8.25, 8.5, 8.75, 9., 9.25, 9.5,
                     9.75, 10., 10.25, 10.5, 10.75, 11., 11.25, 11.5, 11.75])
phi_model = pylab.array([  3.69646470e-02,   2.78002222e-02,   2.09050888e-02,
                     1.57164690e-02,   1.18109411e-02,   8.86993140e-03,
                     6.65381628e-03,   4.98256058e-03,   3.72138301e-03,
                     2.77036826e-03,   2.05721273e-03,   1.53215646e-03,
                     1.16293233e-03,   9.24215390e-04,   7.72929473e-04,
                     6.18448557e-04,   3.69364077e-04,   1.08863332e-04,
                     8.19019897e-06]) / dm

phi = pylab.array([  3.69646470e-02/dm,   2.78002222e-02/dm,   2.09050888e-02/dm,
                     1.57164690e-02/dm,   1.18109411e-02/dm,   8.86993140e-03/dm,
                     6.65381628e-03/dm,   1.9662e-02, 1.5473e-02, 1.1043e-02, 
                     7.7730e-03, 6.2030e-03, 4.7015e-03, 3.6567e-03, 3.2601e-03, 
                     2.3632e-03, 1.4071e-03, 6.0336e-04, 1.2318e-04, 2.9229e-05])



###  dschechter parameters for -0.5<overdens<0 and 0<overdens<0.5 bins
p0 = [10.77, -0.14, -1.52, 0.00021, 0.00014]
p1 = [10.87, -0.59, -1.60, 0.00045, 0.00012]

lmass = pylab.array([7., 7.25, 7.5, 7.75, 8., 8.25, 8.5, 8.75, 9., 9.25, 9.5,
                     9.75, 10., 10.25, 10.5, 10.75, 11., 11.25, 11.5, 11.75])

phi0 = dschechter(lmass, *p0)
phi1 = dschechter(lmass, *p1)
phi = 10**(pylab.log10((phi0*phi1)**0.5))






##################################
###  scenario 1: "stupid" merging
###   - completely random
###   - no mass-loss
##################################

if False:

	###  generating simulated galaxies
	ngal0 = (phi / (phi.min() / 10) * 4).astype(int)
	galaxies = []
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			galaxies.append(galaxy(m_rand))



	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp2 = fig.add_subplot(122)
	sp1 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp1.grid()
	sp2.grid()
	sp1.minorticks_on()
	sp2.minorticks_on()
	sp1.set_yscale('log')
	sp2.set_yscale('log')
	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp1.set_ylabel('Number')
	sp2.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp1.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp1.transAxes)
	sp1.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	ax = sp1.axis()
	sp1.axis([ax[0], ax[1], ax[2]/3., ax[3]])

	sp2.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)





	###  executing random galaxy mergers
	chunk_size = ngal0.sum() / 10


	###  meging 80% of all galaxies in chunks
	fraction_merged = []
	line_fits, line_covs = [], []
	for i_bulk_merge in range(9):

		for i_merge in range(chunk_size):


			###  grabbing random galaxies
			ij_rand = pylab.randint(0, len(galaxies), 2)
			while ij_rand[0] == ij_rand[1]:
				ij_rand = pylab.randint(0, len(galaxies), 2)


			###  merging galaxy masses
			galaxies[ij_rand[1]].lmass = pl.log10(10**galaxies[ij_rand[1]].lmass + 10**galaxies[ij_rand[0]].lmass)


			###  popping merged galaxy
			pop = galaxies.pop(ij_rand[0])


		###  plotting evolved SMF
		ms = pl.array([g.lmass for g in galaxies])
		h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
		sp1.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()), zorder=98-i_bulk_merge)

		###  plotting ratio to inital state
		sp2.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()))


		###  fitting lines to ratio SMFs
		fraction_merged.append((1 - len(galaxies) * 1. / ngal0.sum()))
		fit, cov = optimize.curve_fit(line, lmass, pylab.log10(h[0] * 1. / ngal0))
		line_fits.append(fit)
		line_covs.append(cov.diagonal()**0.5)

		#sp2.plot(lmass, 10**line(lmass, *fit), lw=1, color=cmap(len(galaxies) * 1. / ngal0.sum()))

	sp1.legend(loc=3, fontsize=15, ncol=2, title='merging with no mass-loss')

	fraction_merged = pylab.array(fraction_merged)
	line_fits = pylab.array(line_fits)
	line_covs = pylab.array(line_covs)

	ax1, ax2 = sp1.axis(), sp2.axis()

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_SMF_1.png')
	pylab.close()



	###  plotting % merged vs. slope
	fig = pylab.figure()
	sp3 = fig.add_subplot(111)

	sp3.grid()
	sp3.minorticks_on()
	sp3.set_xlabel('slope')
	sp3.set_ylabel('fraction merged')
	sp3.axis([0, 1, -0.03, 1.03])

	sp3.plot(line_fits[:,0], fraction_merged, 'k', lw=3)

	sp3.text(0.06, 0.06, 'merging with\nno mass-loss', transform=sp1.transAxes)

	measured_slopes = pylab.array([0.39, 0.63, 0.69])
	for i, slope in enumerate(measured_slopes):

		color_i = pylab.cm.cool((i+1)/3.)

		f_interp = pylab.interp(slope, line_fits[:,0], fraction_merged)
		sp3.plot([slope, slope], [-10, f_interp], color=color_i, marker='o', ms=7)
		sp3.plot([0, slope], [f_interp, f_interp], color=color_i, ls='--')

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_slopes_1.png')
	pylab.close()















#############################################
###  scenario 2: 
###   - completely random
###   - 30% mass-loss in less-massive galaxy
#############################################


if False:


	###  generating simulated galaxies
	ngal0 = (phi / (phi.min() / 10) * 4).astype(int)
	galaxies = []
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			galaxies.append(galaxy(m_rand))







	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp5 = fig.add_subplot(122)
	sp4 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp4.grid()
	sp5.grid()
	sp4.minorticks_on()
	sp5.minorticks_on()
	sp4.set_yscale('log')
	sp5.set_yscale('log')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp5.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_ylabel('Number')
	sp5.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp4.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp4.transAxes)
	sp4.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	#sp4.axis(ax1)
	#sp5.axis(ax2)

	sp5.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)





	###  executing random galaxy mergers
	chunk_size = ngal0.sum() / 10


	###  meging 80% of all galaxies in chunks
	fraction_merged = []
	line_fits, line_covs = [], []
	for i_bulk_merge in range(9):

		for i_merge in range(chunk_size):


			###  grabbing random galaxies
			ij_rand = pylab.randint(0, len(galaxies), 2)
			while ij_rand[0] == ij_rand[1]:
				ij_rand = pylab.randint(0, len(galaxies), 2)


			###  merging galaxy masses
			mlower = min([galaxies[ij_rand[1]].lmass, galaxies[ij_rand[0]].lmass])
			galaxies[ij_rand[1]].lmass = pylab.log10(10**galaxies[ij_rand[1]].lmass + 
				                                     10**galaxies[ij_rand[0]].lmass -
				                                     0.3 * 10**mlower)


			###  popping merged galaxy
			pop = galaxies.pop(ij_rand[0])


		###  plotting evolved SMF
		ms = pl.array([g.lmass for g in galaxies])
		h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
		sp4.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()), zorder=98-i_bulk_merge)

		###  plotting ratio to inital state
		sp5.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()))


		###  fitting lines to ratio SMFs
		fraction_merged.append((1 - len(galaxies) * 1. / ngal0.sum()))
		fit, cov = optimize.curve_fit(line, lmass, pylab.log10(h[0] * 1. / ngal0))
		line_fits.append(fit)
		line_covs.append(cov.diagonal()**0.5)

		#sp2.plot(lmass, 10**line(lmass, *fit), lw=1, color=cmap(len(galaxies) * 1. / ngal0.sum()))

	sp4.legend(loc=3, fontsize=15, ncol=2, title='merging with 30% mass-loss')

	fraction_merged = pylab.array(fraction_merged)
	line_fits = pylab.array(line_fits)
	line_covs = pylab.array(line_covs)

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_SMF_2.png')
	pylab.close()







	###  plotting % merged vs. slope
	fig = pylab.figure()
	sp6 = fig.add_subplot(111)

	sp6.grid()
	sp6.minorticks_on()
	sp6.set_xlabel('slope')
	sp6.set_ylabel('fraction merged')
	sp6.axis([0, 1, -0.03, 1.03])

	sp6.plot(line_fits[:,0], fraction_merged, 'k', lw=3)

	sp6.text(0.06, 0.06, 'merging with\n30% mass-loss', transform=sp1.transAxes)


	measured_slopes = pylab.array([0.39, 0.63, 0.69])
	for i, slope in enumerate(measured_slopes):

		color_i = pylab.cm.cool((i+1)/3.)

		f_interp = pylab.interp(slope, line_fits[:,0], fraction_merged)
		sp6.plot([slope, slope], [-10, f_interp], color=color_i, marker='o', ms=7)
		sp6.plot([0, slope], [f_interp, f_interp], color=color_i, ls='--')

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_slopes_2.png')
	pylab.close()












#############################################
###  scenario 3:
###   - mass-dependent merger probability
###   - 30% mass-loss in less-massive galaxy
#############################################


if True:


	###  generating simulated galaxies
	ngal0 = (phi / (phi.min() / 10) * 4).astype(int)
	galaxies = []
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			galaxies.append(galaxy(m_rand))







	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp5 = fig.add_subplot(122)
	sp4 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp4.grid()
	sp5.grid()
	sp4.minorticks_on()
	sp5.minorticks_on()
	sp4.set_yscale('log')
	sp5.set_yscale('log')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp5.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_ylabel('Number')
	sp5.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp4.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp4.transAxes)
	sp4.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	sp4.axis(ax1)
	sp5.axis(ax2)

	sp5.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)





	###  meging galaxies
	fraction_merged = []
	line_fits, line_covs = [], []
	while len(galaxies) > 0.09 * ngal0.sum():


		###  grabbing random galaxies
		ij_rand = pylab.randint(0, len(galaxies), 2)
		while ij_rand[0] == ij_rand[1]:
			ij_rand = pylab.randint(0, len(galaxies), 2)


		###  checking merger probability
		masses = pylab.array([galaxies[ij_rand[1]].lmass, galaxies[ij_rand[0]].lmass])
		mlower = masses.min()
		i_mlower = pylab.find(masses == mlower)[0]
		if pylab.rand() > galaxies[ij_rand[i_mlower]].merge_prob:
			continue



		###  merging galaxy masses
		galaxies[ij_rand[1]].lmass = pylab.log10(10**galaxies[ij_rand[1]].lmass + 
			                                     10**galaxies[ij_rand[0]].lmass -
			                                     0.3 * 10**mlower)
		galaxies[ij_rand[1]].merge_prob = get_merge_prob(galaxies[ij_rand[1]].lmass)


		###  popping merged galaxy
		pop = galaxies.pop(ij_rand[0])




		if len(galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 10):

			###  plotting evolved SMF
			ms = pl.array([g.lmass for g in galaxies])
			h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
			sp4.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()), zorder=98-i_bulk_merge)

			###  plotting ratio to inital state
			sp5.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()))


			###  fitting lines to ratio SMFs
			fraction_merged.append((1 - len(galaxies) * 1. / ngal0.sum()))
			fit, cov = optimize.curve_fit(line, lmass, pylab.log10(h[0] * 1. / ngal0))
			line_fits.append(fit)
			line_covs.append(cov.diagonal()**0.5)

			#sp2.plot(lmass, 10**line(lmass, *fit), lw=1, color=cmap(len(galaxies) * 1. / ngal0.sum()))

	sp4.legend(loc=3, fontsize=15, ncol=2, title='Leja+2015 rates + 30% mass-loss')

	fraction_merged = pylab.array(fraction_merged)
	line_fits = pylab.array(line_fits)
	line_covs = pylab.array(line_covs)

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_SMF_3.png')
	pylab.close()







	###  plotting % merged vs. slope
	fig = pylab.figure()
	sp6 = fig.add_subplot(111)

	sp6.grid()
	sp6.minorticks_on()
	sp6.set_xlabel('slope')
	sp6.set_ylabel('fraction merged')
	sp6.axis([0, 1, -0.03, 1.03])

	sp6.plot(line_fits[:,0], fraction_merged, 'k', lw=3)

	sp6.text(0.06, 0.06, 'scenario 3', transform=sp1.transAxes)


	measured_slopes = pylab.array([0.39, 0.63, 0.69])
	for i, slope in enumerate(measured_slopes):

		color_i = pylab.cm.cool((i+1)/3.)

		f_interp = pylab.interp(slope, line_fits[:,0], fraction_merged)
		sp6.plot([slope, slope], [-10, f_interp], color=color_i, marker='o', ms=7)
		sp6.plot([0, slope], [f_interp, f_interp], color=color_i, ls='--')

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_slopes_3.png')
	pylab.close()












###############################################
###  scenario 4:
###   - Lotz+2011 rate: minor merger 3x likely
###   - 30% mass-loss in less-massive galaxy
###############################################


if True:


	###  generating simulated galaxies
	ngal0 = (phi / (phi.min() / 10) * 4).astype(int)
	galaxies = []
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			galaxies.append(galaxy(m_rand))



	simdata = my_classes.merger_simulation(lmass)



	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp5 = fig.add_subplot(122)
	sp4 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp4.grid()
	sp5.grid()
	sp4.minorticks_on()
	sp5.minorticks_on()
	sp4.set_yscale('log')
	sp5.set_yscale('log')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp5.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_ylabel('Number')
	sp5.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp4.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp4.transAxes)
	sp4.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	#sp4.axis(ax1)
	#sp5.axis(ax2)

	sp5.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)



	###  meging galaxies
	fraction_merged = []
	line_fits, line_covs = [], []
	while len(galaxies) > 0.09 * ngal0.sum():


		###  grabbing random galaxies
		ij_rand = pylab.randint(0, len(galaxies), 2)
		while ij_rand[0] == ij_rand[1]:
			ij_rand = pylab.randint(0, len(galaxies), 2)



		###  checking merger probability
		masses = pylab.array([galaxies[ij_rand[1]].lmass, galaxies[ij_rand[0]].lmass])
		mratio = 10**(masses.min() - masses.max())
		mlower = masses.min()
		i_mlower = pylab.find(masses == mlower)[0]

		###  definition of major merger 1:4
		if mratio > 0.25:
			p = pylab.rand()
			if p > 1./3:
				continue



		###  merging galaxy masses
		galaxies[ij_rand[1]].lmass = pylab.log10(10**galaxies[ij_rand[1]].lmass + 
			                                     10**galaxies[ij_rand[0]].lmass -
			                                     0.3 * 10**mlower)
		galaxies[ij_rand[1]].merge_prob = get_merge_prob(galaxies[ij_rand[1]].lmass)


		###  popping merged galaxy
		pop = galaxies.pop(ij_rand[0])




		#if len(galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 10):
		if len(galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 20):
		#if len(galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 100):

			###  plotting evolved SMF
			ms = pl.array([g.lmass for g in galaxies])
			h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
			sp4.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()))

			###  saving simulation data
			simdata.fraction_merged = pylab.append(simdata.fraction_merged, 1 - len(galaxies) * 1. / ngal0.sum())
			simdata.ngals.append(h[0])

			###  plotting ratio to inital state
			sp5.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()))


			###  fitting lines to ratio SMFs
			fraction_merged.append((1 - len(galaxies) * 1. / ngal0.sum()))
			fit, cov = optimize.curve_fit(line, lmass, pylab.log10(h[0] * 1. / ngal0))
			line_fits.append(fit)
			line_covs.append(cov.diagonal()**0.5)

			#sp2.plot(lmass, 10**line(lmass, *fit), lw=1, color=cmap(len(galaxies) * 1. / ngal0.sum()))

	sp4.legend(loc=3, fontsize=15, ncol=2, title='Lotz+2011 rates + 30% mass-loss')

	fraction_merged = pylab.array(fraction_merged)
	line_fits = pylab.array(line_fits)
	line_covs = pylab.array(line_covs)

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_SMF_4.png')
	pylab.close()







	###  plotting % merged vs. slope
	fig = pylab.figure()
	sp6 = fig.add_subplot(111)

	sp6.grid()
	sp6.minorticks_on()
	sp6.set_xlabel('slope')
	sp6.set_ylabel('fraction merged')
	sp6.axis([0, 1, -0.03, 1.03])

	sp6.plot(line_fits[:,0], fraction_merged, 'k', lw=3)

	sp6.text(0.06, 0.06, 'scenario 3', transform=sp6.transAxes)


	measured_slopes = pylab.array([0.39, 0.63, 0.69])
	for i, slope in enumerate(measured_slopes):

		color_i = pylab.cm.cool((i+1)/3.)

		f_interp = pylab.interp(slope, line_fits[:,0], fraction_merged)
		sp6.plot([slope, slope], [-10, f_interp], color=color_i, marker='o', ms=7)
		sp6.plot([0, slope], [f_interp, f_interp], color=color_i, ls='--')

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_slopes_4.png')
	pylab.close()


	with open('../data/merger_simulation_scenario4.pickle', 'wb') as outer:
		pickle.dump(simdata, outer)





	###  figure for paper
	cmap = pylab.cm.jet
	fig_paper = pylab.figure(figsize=(9.5, 7.5))
	sp_paper = fig_paper.add_subplot(111)
	sp_cbar = add_inset(sp_paper, [0.6, 0.85, 0.33, 0.05])

	sp_paper.grid()
	sp_paper.minorticks_on()
	sp_paper.set_yscale('log')
	sp_paper.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp_paper.set_ylabel('Number')

	sp_paper.plot(lmass, ngal0, color='k', lw=7, label='initial state', zorder=99)


	sp_cbar.set_title('Fraction Merged', fontsize=18)
	sp_cbar.set_yticks([])
	imdatx, imdaty = pylab.meshgrid(pylab.linspace(1, 0, 101), [0, 0])
	sp_cbar.imshow(imdatx, cmap=cmap)
	sp_cbar.set_aspect('auto')

	for i in range(len(simdata.fraction_merged)):
		sp_paper.plot(lmass, simdata.ngals[i], marker='o', ms=4, mew=0.5, color=cmap(1-simdata.fraction_merged[i]))












###############################################
###  scenario 4b:
###   - Same as scenario 4 but now tracking 
###     ex-situ mass growth
###############################################


if True:


	###  generating simulated galaxies
	mock_catalog = mock_catalog_class()
	ngal0 = (phi / (phi.min() / 10) * 4).astype(int)
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			mock_catalog.galaxies.append(galaxy(m_rand))
			mock_catalog.mass_initial += 10**m_rand


	###  adding initial number counts
	mock_catalog.n_initial = ngal0.sum()
	for g in mock_catalog.galaxies:
		if g.lmass > 9.:
			mock_catalog.n_initial_gt_9 += 1




	simdata = my_classes.merger_simulation(lmass)



	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp5 = fig.add_subplot(122)
	sp4 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp4.grid()
	sp5.grid()
	sp4.minorticks_on()
	sp5.minorticks_on()
	sp4.set_yscale('log')
	sp5.set_yscale('log')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp5.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_ylabel('Number')
	sp5.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp4.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp4.transAxes)
	sp4.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	#sp4.axis(ax1)
	#sp5.axis(ax2)

	sp5.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)



	###  meging galaxies
	hit = 0
	fraction_merged = []
	line_fits, line_covs = [], []

	while mock_catalog.n_merged_gt_9 * 1. / mock_catalog.n_initial_gt_9 < 0.9:

		mypy.progress_bar(int(mock_catalog.n_merged_gt_9 * 100 / mock_catalog.n_initial_gt_9), 100)

		###  grabbing random galaxies
		ij_rand = pylab.randint(0, len(mock_catalog.galaxies), 2)
		while ij_rand[0] == ij_rand[1]:
			ij_rand = pylab.randint(0, len(mock_catalog.galaxies), 2)



		###  checking merger probability
		masses = pylab.array([mock_catalog.galaxies[ij_rand[1]].lmass, mock_catalog.galaxies[ij_rand[0]].lmass])
		mratio = 10**(masses.min() - masses.max())
		mlower, mhigher = masses.min(), masses.max()
		i_mlower = pylab.find(masses == mlower)[0]
		i_mhigher = pylab.find(masses == mhigher)[0]



		###  definition of major merger 1:4
		###  enforcing minor mergers 3x more frequent than major
		if mratio > 0.25:
			p = pylab.rand()
			if p > 1./3:
				continue



		###  merging galaxy masses
		mock_catalog.galaxies[ij_rand[i_mhigher]].lmass = pylab.log10(10**mock_catalog.galaxies[ij_rand[1]].lmass + 
			                                             10**mock_catalog.galaxies[ij_rand[0]].lmass -
			                                             0.3 * 10**mlower)
		mock_catalog.galaxies[ij_rand[i_mhigher]].merge_prob = get_merge_prob(mock_catalog.galaxies[ij_rand[i_mhigher]].lmass)

		if mratio > 0.25:
			mock_catalog.galaxies[ij_rand[i_mhigher]].n_major_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].mass_from_major_mergers += 10**mock_catalog.galaxies[ij_rand[i_mlower]].lmass * 0.7
		else:
			mock_catalog.galaxies[ij_rand[i_mhigher]].n_minor_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].mass_from_minor_mergers += 10**mock_catalog.galaxies[ij_rand[i_mlower]].lmass * 0.7



		###  If seed mass of less massive galaxy is >10**9 increment merger count
		if mock_catalog.galaxies[ij_rand[i_mlower]].lmass0 >= 9.:
			mock_catalog.n_merged_gt_9 += 1


		###  popping merged galaxy
		pop = mock_catalog.galaxies.pop(ij_rand[i_mlower])




		###  saving snapshot: poor method
		if int(100 - len(mock_catalog.galaxies) * 100. / ngal0.sum()) == 77 and not hit:

			hit = 1
			snap_lmasses = []
			snap_n_minor = []
			snap_n_major = []
			snap_mass_minor = []
			snap_mass_major = []

			for g in mock_catalog.galaxies:
				snap_lmasses.append(g.lmass)
				snap_n_minor.append(g.n_minor_mergers)
				snap_n_major.append(g.n_major_mergers)
				snap_mass_minor.append(g.mass_from_minor_mergers)
				snap_mass_major.append(g.mass_from_major_mergers)

			snap_lmasses = pylab.array(snap_lmasses)
			snap_n_minor = pylab.array(snap_n_minor)
			snap_n_major = pylab.array(snap_n_major)
			snap_mass_minor = pylab.array(snap_mass_minor)
			snap_mass_major = pylab.array(snap_mass_major)







		###  plotting intermittently

		#if len(mock_catalog.galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 10):
		#if len(mock_catalog.galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 20):
		#if len(mock_catalog.galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 100):

		if mock_catalog.n_merged_gt_9 in range(mock_catalog.n_initial_gt_9 / 10, mock_catalog.n_initial_gt_9, mock_catalog.n_initial_gt_9 / 20):

			fmerged = mock_catalog.n_merged_gt_9 * 1. / mock_catalog.n_initial_gt_9

			###  plotting evolved SMF
			ms = pylab.array([g.lmass for g in mock_catalog.galaxies])
			h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
			sp4.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(fmerged), label='%i%% merged' % round(fmerged * 100))

			###  saving simulation data
			simdata.fraction_merged = pylab.append(simdata.fraction_merged, 1 - len(mock_catalog.galaxies) * 1. / ngal0.sum())
			simdata.ngals.append(h[0])

			###  plotting ratio to inital state
			sp5.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(fmerged), label='%i%% merged' % round(100 - len(mock_catalog.galaxies) * 100. / ngal0.sum()))


			###  fitting lines to ratio SMFs
			fraction_merged.append((1 - len(mock_catalog.galaxies) * 1. / ngal0.sum()))
			fit, cov = optimize.curve_fit(line, lmass, pylab.log10(h[0] * 1. / ngal0))
			line_fits.append(fit)
			line_covs.append(cov.diagonal()**0.5)

			#sp2.plot(lmass, 10**line(lmass, *fit), lw=1, color=cmap(len(mock_catalog.galaxies) * 1. / ngal0.sum()))

	#sp4.legend(loc=3, fontsize=15, ncol=2, title='Lotz+2011 rates + 30% mass-loss')

	fraction_merged = pylab.array(fraction_merged)
	line_fits = pylab.array(line_fits)
	line_covs = pylab.array(line_covs)
























	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-15/simulation_SMF_4.png')
	pylab.close()







	###  plotting % merged vs. slope
	fig = pylab.figure()
	sp6 = fig.add_subplot(111)

	sp6.grid()
	sp6.minorticks_on()
	sp6.set_xlabel('slope')
	sp6.set_ylabel('fraction merged')
	sp6.axis([0, 1, -0.03, 1.03])

	sp6.plot(line_fits[:,0], fraction_merged, 'k', lw=3)

	sp6.text(0.06, 0.06, 'scenario 3', transform=sp6.transAxes)


	measured_slopes = pylab.array([0.39, 0.63, 0.69])
	for i, slope in enumerate(measured_slopes):

		color_i = pylab.cm.cool((i+1)/3.)

		f_interp = pylab.interp(slope, line_fits[:,0], fraction_merged)
		sp6.plot([slope, slope], [-10, f_interp], color=color_i, marker='o', ms=7)
		sp6.plot([0, slope], [f_interp, f_interp], color=color_i, ls='--')

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-15/simulation_slopes_4.png')
	pylab.close()


	with open('../data/merger_simulation_scenario4.pickle', 'wb') as outer:
		pickle.dump(simdata, outer)





	###  figure for paper
	cmap = pylab.cm.jet
	fig_paper = pylab.figure(figsize=(9.5, 7.5))
	sp_paper = fig_paper.add_subplot(111)
	sp_cbar = add_inset(sp_paper, [0.6, 0.85, 0.33, 0.05])

	sp_paper.grid()
	sp_paper.minorticks_on()
	sp_paper.set_yscale('log')
	sp_paper.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp_paper.set_ylabel('Number')

	sp_paper.plot(lmass, ngal0, color='k', lw=7, label='initial state', zorder=99)


	sp_cbar.set_title('Fraction Merged', fontsize=18)
	sp_cbar.set_yticks([])
	imdatx, imdaty = pylab.meshgrid(pylab.linspace(1, 0, 101), [0, 0])
	sp_cbar.imshow(imdatx, cmap=cmap)
	sp_cbar.set_aspect('auto')

	for i in range(len(simdata.fraction_merged)):
		sp_paper.plot(lmass, simdata.ngals[i], marker='o', ms=4, mew=0.5, color=cmap(1-simdata.fraction_merged[i]))

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-15/simulation_fmerged_4.png')








###############################################
###  scenario 5:
###   - Inverted Lotz+2011 rate: major merger 3x likely
###   - 30% mass-loss in less-massive galaxy
###############################################


if True:


	###  generating simulated galaxies
	ngal0 = (phi / (phi.min() / 10) * 4).astype(int)
	galaxies = []
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			galaxies.append(galaxy(m_rand))







	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp5 = fig.add_subplot(122)
	sp4 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp4.grid()
	sp5.grid()
	sp4.minorticks_on()
	sp5.minorticks_on()
	sp4.set_yscale('log')
	sp5.set_yscale('log')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp5.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_ylabel('Number')
	sp5.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp4.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp4.transAxes)
	sp4.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	sp4.axis(ax1)
	sp5.axis(ax2)

	sp5.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)





	###  meging galaxies
	fraction_merged = []
	line_fits, line_covs = [], []
	while len(galaxies) > 0.09 * ngal0.sum():


		###  grabbing random galaxies
		ij_rand = pylab.randint(0, len(galaxies), 2)
		while ij_rand[0] == ij_rand[1]:
			ij_rand = pylab.randint(0, len(galaxies), 2)



		###  checking merger probability
		masses = pylab.array([galaxies[ij_rand[1]].lmass, galaxies[ij_rand[0]].lmass])
		mratio = 10**(masses.min() - masses.max())
		mlower = masses.min()
		i_mlower = pylab.find(masses == mlower)[0]

		###  definition of major merger 1:4
		if mratio < 0.25:
			p = pylab.rand()
			if p > 1./3:
				continue



		###  merging galaxy masses
		galaxies[ij_rand[1]].lmass = pylab.log10(10**galaxies[ij_rand[1]].lmass + 
			                                     10**galaxies[ij_rand[0]].lmass -
			                                     0.3 * 10**mlower)
		galaxies[ij_rand[1]].merge_prob = get_merge_prob(galaxies[ij_rand[1]].lmass)


		###  popping merged galaxy
		pop = galaxies.pop(ij_rand[0])




		if len(galaxies) in range(ngal0.sum() / 10, ngal0.sum() * 9 / 10, ngal0.sum() / 10):

			###  plotting evolved SMF
			ms = pl.array([g.lmass for g in galaxies])
			h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
			sp4.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()), zorder=98-i_bulk_merge)

			###  plotting ratio to inital state
			sp5.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(len(galaxies) * 1. / ngal0.sum()), label='%i%% merged' % round(100 - len(galaxies) * 100. / ngal0.sum()))


			###  fitting lines to ratio SMFs
			fraction_merged.append((1 - len(galaxies) * 1. / ngal0.sum()))
			fit, cov = optimize.curve_fit(line, lmass, pylab.log10(h[0] * 1. / ngal0))
			line_fits.append(fit)
			line_covs.append(cov.diagonal()**0.5)

			#sp2.plot(lmass, 10**line(lmass, *fit), lw=1, color=cmap(len(galaxies) * 1. / ngal0.sum()))

	sp4.legend(loc=3, fontsize=15, ncol=2, title='Lotz+2011 inverted + 30% mass-loss')

	fraction_merged = pylab.array(fraction_merged)
	line_fits = pylab.array(line_fits)
	line_covs = pylab.array(line_covs)

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_SMF_5.png')
	pylab.close()







	###  plotting % merged vs. slope
	fig = pylab.figure()
	sp6 = fig.add_subplot(111)

	sp6.grid()
	sp6.minorticks_on()
	sp6.set_xlabel('slope')
	sp6.set_ylabel('fraction merged')
	sp6.axis([0, 1, -0.03, 1.03])

	sp6.plot(line_fits[:,0], fraction_merged, 'k', lw=3)

	sp6.text(0.06, 0.06, 'scenario 3', transform=sp1.transAxes)


	measured_slopes = pylab.array([0.39, 0.63, 0.69])
	for i, slope in enumerate(measured_slopes):

		color_i = pylab.cm.cool((i+1)/3.)

		f_interp = pylab.interp(slope, line_fits[:,0], fraction_merged)
		sp6.plot([slope, slope], [-10, f_interp], color=color_i, marker='o', ms=7)
		sp6.plot([0, slope], [f_interp, f_interp], color=color_i, ls='--')

	pylab.savefig('/Users/atomczak/Dropbox/show-and-tell/2017-02-02/simulation_slopes_5.png')
	pylab.close()












##################################################
###  scenario 6:
###   - Progress will be checked by the growth 
###     of the ICL fraction. That is the fraction
###     of the total stellar mass that is in the
###     ICL (as contributed from mergers)
##################################################


if True:

	###  fraction of stellar mass of smaller galaxy the becomes unbound per merger
	merger_stripped_fraction = 0.3


	###  generating simulated galaxies
	mock_catalog = mock_catalog_class()
	ngal0 = (phi / phi.min() * 4).astype(int)
	print '\ngenerating %i simulated galaxies...' % ngal0.sum()
	for i in range(len(ngal0)):
		ni = ngal0[i]
		mi = lmass[i]
		for j in range(ni):
			m_rand = pylab.rand() * dm + (mi - dm/2.)
			mock_catalog.galaxies.append(galaxy(m_rand))
			mock_catalog.mass_initial += 10**m_rand

	mock_catalog.n_initial = len(mock_catalog.galaxies)
	print 'done!\n'





	simdata = my_classes.merger_simulation(lmass)
	simdata.n_initial = ngal0.sum()
	simdata.mass_initial = mock_catalog.mass_initial

	###  adding to simulation pickle
	simdata.fraction_merged = pylab.append(simdata.fraction_merged, 0.)
	simdata.mass_in_icm = pylab.append(simdata.mass_in_icm, 0.)
	simdata.f_icm = pylab.append(simdata.f_icm, 0.)
	simdata.ngals_icm.append(ngal0 * 1.)




	###  initiallizing figure
	fig = pylab.figure(figsize=(15.3, 6.8))
	sp5 = fig.add_subplot(122)
	sp4 = fig.add_subplot(121)
	fig.subplots_adjust(left=0.09)

	sp4.grid()
	sp5.grid()
	sp4.minorticks_on()
	sp5.minorticks_on()
	sp4.set_yscale('log')
	sp5.set_yscale('log')
	sp4.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp5.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
	sp4.set_ylabel('Number')
	sp5.set_ylabel('Ratio')
	cmap = pylab.cm.jet

	sp4.text(0.65, 0.92, 'N$_0$ = %i' % ngal0.sum(), transform=sp4.transAxes)
	sp4.plot(lmass, ngal0, color='k', lw=5, label='initial state', zorder=99)
	#sp4.axis(ax1)
	#sp5.axis(ax2)

	sp5.plot(lmass, ngal0 / ngal0, color='k', lw=5, label='initial state', zorder=99)



	###  tracking icm growth for plotting
	percent_hits = pylab.arange(0.02, 0.19, 0.02).astype(str)
	percent_hits = pylab.arange(0.005, 0.19, 0.001).astype(str)
	hits = pylab.zeros(len(percent_hits))
	f_icm_progress = []


	f_icm_max = float(percent_hits[-1])
	while mock_catalog.mass_in_icm / mock_catalog.mass_initial < f_icm_max:

		f_icm = mock_catalog.mass_in_icm / mock_catalog.mass_initial
		f_icm_progress.append(round(f_icm, 2))
		#mypy.progress_bar(int(f_icm * 100), 100)


		###  grabbing random galaxies
		ij_rand = pylab.randint(0, len(mock_catalog.galaxies), 2)
		while ij_rand[0] == ij_rand[1]:
			ij_rand = pylab.randint(0, len(mock_catalog.galaxies), 2)



		###  checking merger probability
		masses = pylab.array([mock_catalog.galaxies[ij_rand[1]].lmass, mock_catalog.galaxies[ij_rand[0]].lmass])
		mratio = 10**(masses.min() - masses.max())
		mlower, mhigher = masses.min(), masses.max()
		i_mlower = pylab.find(masses == mlower)[0]
		i_mhigher = pylab.find(masses == mhigher)[0]



		###  definition of major merger 1:4
		###  enforcing minor mergers 3x more frequent than major
		if mratio > 0.25:
			p = pylab.rand()
			if p > 1./3:
				continue



		###  merging galaxy masses
		mock_catalog.galaxies[ij_rand[i_mhigher]].lmass = pylab.log10(10**mock_catalog.galaxies[ij_rand[1]].lmass + 
			                                                          10**mock_catalog.galaxies[ij_rand[0]].lmass -
			                                                          merger_stripped_fraction * 10**mlower)

		mock_catalog.mass_in_icm += merger_stripped_fraction * 10**mlower

		###  Track number of major/minor mergers for the more massive galaxy
		if mratio >= 0.25:
			mock_catalog.n_major_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].n_major_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].mass_from_major_mergers += 10**mock_catalog.galaxies[ij_rand[i_mlower]].lmass * 0.7
		elif 0.1 <= mratio < 0.25:
			mock_catalog.n_minor_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].n_minor_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].mass_from_minor_mergers += 10**mock_catalog.galaxies[ij_rand[i_mlower]].lmass * 0.7
		elif 0.01 <= mratio < 0.1:
			mock_catalog.n_vminor_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].n_vminor_mergers += 1
			mock_catalog.galaxies[ij_rand[i_mhigher]].mass_from_vminor_mergers += 10**mock_catalog.galaxies[ij_rand[i_mlower]].lmass * 0.7

		
		###  counting mergers at lmass > 10**9
		if min([mock_catalog.galaxies[ij_rand[0]].lmass0, mock_catalog.galaxies[ij_rand[1]].lmass0]) >= 9.: mock_catalog.n_merged_gt_9 += 1


		###  popping merged galaxy
		pop = mock_catalog.galaxies.pop(ij_rand[i_mlower])





		###  plotting if in percent_hits
		if str(round(f_icm, 3)) in percent_hits:
			i_hit = percent_hits.tolist().index(str(round(f_icm, 3)))
			if hits[i_hit] == 0:
				print 'hit ', round(f_icm, 3)
				hits[i_hit] = 1

				###  plotting
				c = f_icm / f_icm_max
				ms = pylab.array([g.lmass for g in mock_catalog.galaxies])
				h = pylab.histogram(ms, range=(lmass[0]-dm/2., lmass[-1]+dm/2.), bins=len(lmass))
				sp4.plot(lmass, h[0], 'r', marker='o', ms=5, mew=0.5, color=cmap(c), label='%i%% in icm' % round(f_icm * 100))


				###  plotting ratio to inital state
				sp5.plot(lmass, h[0] * 1. / ngal0, 'r', marker='o', ms=5, mew=0.5, color=cmap(c), label='%i%% merged' % round(100 - len(mock_catalog.galaxies) * 100. / ngal0.sum()))


				###  adding to simulation pickle
				simdata.fraction_merged = pylab.append(simdata.fraction_merged, 1 - len(mock_catalog.galaxies) * 1. / mock_catalog.n_initial)
				simdata.mass_in_icm = pylab.append(simdata.mass_in_icm, mock_catalog.mass_in_icm)
				simdata.f_icm = pylab.append(simdata.f_icm, mock_catalog.mass_in_icm / mock_catalog.mass_initial)
				simdata.ngals_icm.append(h[0])



				###  adding mean numbers of major/minor mergers in bins of lmass
				lmasses = pylab.array([g.lmass for g in mock_catalog.galaxies])
				lmasses0 = pylab.array([g.lmass0 for g in mock_catalog.galaxies])
				n_vmin = pylab.array([g.n_vminor_mergers for g in mock_catalog.galaxies])
				n_min = pylab.array([g.n_minor_mergers for g in mock_catalog.galaxies])
				n_maj = pylab.array([g.n_major_mergers for g in mock_catalog.galaxies])
				mass_min = pylab.array([g.mass_from_minor_mergers for g in mock_catalog.galaxies])
				mass_maj = pylab.array([g.mass_from_major_mergers for g in mock_catalog.galaxies])

				digi = pylab.digitize(lmasses, lmassbins_merger_sim)
				n_vminors, n_minors, n_majors  = pylab.zeros(len(lmassbars_merger_sim)), pylab.zeros(len(lmassbars_merger_sim)), pylab.zeros(len(lmassbars_merger_sim))
				for di in range(1, len(lmassbins_merger_sim)):
					sinds = pylab.find(digi == di)
					n_vminors[di-1] = pylab.mean(n_vmin[sinds])
					n_minors[di-1] = pylab.mean(n_min[sinds])
					n_majors[di-1] = pylab.mean(n_maj[sinds])
				simdata.n_vminor_mergers.append(n_vminors)
				simdata.n_minor_mergers.append(n_minors)
				simdata.n_major_mergers.append(n_majors)



	with open('../merger_simulation_scenario6_test.pickle', 'wb') as outer:
		pickle.dump(simdata, outer)















#####################
###  Plotting figure
#####################


simdata = pickle.load(open('../merger_simulation_scenario6_test.pickle', 'rb'))

inds_lmass = pylab.find((simdata.lmassax > 8.5) & (simdata.lmassax < 11.55))


###  figure for paper
cmap = pylab.cm.jet
fig_paper = pylab.figure(figsize=(9.5, 7.5))
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
sp_paper.text(8.65-0.2, smf[inds_lmass][0]*2, 'fraction\nmerged', fontproperties=font,
	          verticalalignment='bottom', horizontalalignment='center', fontsize=22)
sp_paper.text(8.65-0.05, smf[inds_lmass][0], '%.1f%%' % (simdata.fraction_merged[0]*100), fontproperties=font,
	          verticalalignment='center', horizontalalignment='right', fontsize=15, fontweight='normal')

























