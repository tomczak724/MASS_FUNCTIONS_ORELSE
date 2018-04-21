
import mypy
import pylab
from mypy import massfunc
from scipy import optimize





zbins = pylab.array([[0.20, 0.50],
                     [0.50, 0.75],
                     [0.75, 1.00],
                     [1.00, 1.25],
                     [1.25, 1.50],
                     [1.50, 2.00]])
zbars = (zbins[:,1] + zbins[:,0]) / 2.
dz = (zbins[:,1] - zbins[:,0])




dat = pylab.loadtxt('../output/smf_7fields.dat')
dm = dat[1,0] - dat[0,0]



fig = pylab.figure(figsize=(16.1, 13.))
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

for zi in range(len(zbins)):

	zlo, zhi = zbins[zi]

	sp = sps[zi]
	sp.grid()
	sp.minorticks_on()
	axis_lims = [8.7, 11.9, -5.4, -0.9]
	sp.axis(axis_lims)
	#sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
	sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)


	sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


	###  plotting orelse
	x = dat[:,0]
	y = pylab.log10(dat[:,2*zi+1])
	yerr = pylab.log10(dat[:,2*zi+1]) - pylab.log10(dat[:,2*zi+1] - dat[:,2*zi+2])
	sp.errorbar(x, y, xerr=dm/2., yerr=yerr, ls='', marker='o', ms=7, mew=1.5, \
		        elinewidth=1.5, mfc='k', ecolor='k', label='ORELSE', zorder=2)


	###  fitting double Schechter
	x = dat[:,0]
	y = dat[:,2*zi+1]
	yerr = dat[:,2*zi+2]
	fit, cov = optimize.curve_fit(massfunc.dschechter_mf, x, y, p0=[10.78, -0.98, -2.54, -1.9, -4.29], sigma=yerr)
	dschechter_fits.append(fit)
	dschechter_errs.append(cov.diagonal()**0.5)

	xmodel = pylab.linspace(9, 11.7, 1000)
	ymodel = pylab.log10(massfunc.dschechter_mf(xmodel, *fit))
	sp.plot(xmodel, ymodel, color='r', label='2x Schechter fit', zorder=1)

	if zi == 3:
		title = 'All fields minus LSS'
		leg = sp.legend(loc=3, fontsize=15, title=title)
		leg.get_title().set_fontsize(15)

fig.savefig('../output/schechter_fits.png')
pylab.close()












###  plotting schechter params vs. redshift

dschechter_fits = pylab.array(dschechter_fits)
dschechter_errs = pylab.array(dschechter_errs)

zfourge_params = mypy.readcat('/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table2_TOT.dat')
inds = pylab.arange(len(zbins))


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
sp1.errorbar(zbars + 0.015, y, xerr=dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
	         ecolor='r', elinewidth=2, label='ORELSE today')


sp1.errorbar(zbars - 0.015, zfourge_params.Mstar[inds], xerr=dz/2., yerr=zfourge_params.e_Mstar[inds], \
	         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')






###  alpha2
sp2.set_xlabel('redshfit')
sp2.set_ylabel('alpha 2')
sp2.axis([-0.2, 2.3, -2.2, -1])

y = dschechter_fits[:,3]
yerr = dschechter_errs[:,3]
sp2.errorbar(zbars + 0.015, y, xerr=dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
	         ecolor='r', elinewidth=2, label='ORELSE today')


sp2.errorbar(zbars - 0.015, zfourge_params.alpha2[inds], xerr=dz/2., yerr=zfourge_params.e_alpha2[inds], \
	         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')







###  log(phistar1 + phistar2)
sp3.set_xlabel('redshfit')
sp3.set_ylabel('log( $\Phi_1^*$ + $\Phi_2^*$ )')
sp3.axis([-0.2, 2.3, -3.4, -2.1])

y = pylab.log10(10**dschechter_fits[:,2] + 10**dschechter_fits[:,4])

dphi1 = (10**dschechter_errs[:,2] - 1) * 10**dschechter_fits[:,2]
dphi2 = (10**dschechter_errs[:,4] - 1) * 10**dschechter_fits[:,4]
yerr = pylab.sqrt((dphi1**2 + dphi2**2) / pylab.log(10)**2 / (10**dschechter_fits[:,2] + 10**dschechter_fits[:,4])**2)

sp3.errorbar(zbars + 0.015, y, xerr=dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
	         ecolor='r', elinewidth=2, label='ORELSE today')


y = pylab.log10(10**zfourge_params.phistar1 + 10**zfourge_params.phistar2)

dphi1 = (10**zfourge_params.e_phistar1 - 1) * 10**zfourge_params.phistar1
dphi2 = (10**zfourge_params.e_phistar2 - 1) * 10**zfourge_params.phistar2
yerr = pylab.sqrt((dphi1**2 + dphi2**2) / pylab.log(10)**2 / (10**zfourge_params.phistar1 + 10**zfourge_params.phistar2)**2)

sp3.errorbar(zbars - 0.015, y[inds], xerr=dz/2., yerr=yerr[inds], \
	         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')

sp3.legend(loc=3, fontsize=18)


fig.savefig('../output/schechter_params.png')
pylab.close()
















