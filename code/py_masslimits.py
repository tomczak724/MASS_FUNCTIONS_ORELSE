
import mypy
import pylab
import ezgal
from matplotlib.backends.backend_pdf import PdfPages





data_dir = '../data/completeness'

class masslim_data:
	def __init__(self, name, name2, zclust, sigmaz, filterinfo=['', -1], maglim=0):
		self.name = name
		self.name2 = name2
		self.zclust = zclust
		self.sigmaz = sigmaz
		self.filterinfo = filterinfo
		self.maglim = maglim
		self.data = mypy.readcat('%s/masslimits_%s_%s.dat' % (data_dir, name, filterinfo[0]))





#########################################
###  mass limits from maximally old SSP
#########################################

ezgal_dir = '/Users/atomczak/pydir/ezgal'
bc03_dir = '/Users/atomczak/DATA/TEMPLATES/bc03/models'
model_ssp = ezgal.model('%s/bc03_ssp_z_0.02_chab.model' % bc03_dir)
model_exp = ezgal.model('%s/bc03_exp_1.0_z_0.02_chab.model' % bc03_dir)

model_ssp.add_filter('%s/data/filters/sloan_r' % ezgal_dir)    # index = 0
model_ssp.add_filter('%s/data/filters/sloan_i' % ezgal_dir)    # index = 1
model_ssp.add_filter('%s/data/filters/sloan_z' % ezgal_dir)    # index = 2
model_ssp.add_filter('%s/data/filters/ukidss_y' % ezgal_dir)   # index = 3

model_exp.add_filter('%s/data/filters/sloan_r' % ezgal_dir)    # index = 0
model_exp.add_filter('%s/data/filters/sloan_i' % ezgal_dir)    # index = 1
model_exp.add_filter('%s/data/filters/sloan_z' % ezgal_dir)    # index = 2
model_exp.add_filter('%s/data/filters/ukidss_y' % ezgal_dir)   # index = 3

zform = 5.






masslim_fields = []
masslim_fields.append(masslim_data('nep200_v0.0.6',       'N200',    0.691,  0.027, ['I', 1], 24.67))
masslim_fields.append(masslim_data('rxj1221+4918_v0.0.2', 'RXJ1221', 0.700,  0.023, ['I', 1], 24.35))
masslim_fields.append(masslim_data('sc1324_v0.0.3',       'SC1324',  0.755,  0.033, ['I', 1], 24.31))
masslim_fields.append(masslim_data('rcs0224-0002_v0.0.3', 'RCS0224', 0.772,  0.027, ['I', 1], 25.77))

masslim_fields.append(masslim_data('cl1350+6007_v0.0.1',  'CL1350',  0.804,  0.035, ['R', 0], 24.29))
#masslim_fields.append(masslim_data('cl1350+6007_v0.0.1',  'CL1350',  0.804,  0.035, ['R', 0], 25.07))
#masslim_fields.append(masslim_data('cl1350+6007_v0.0.1',  'CL1350',  0.804,  0.035, ['I', 1], 23.53))
#masslim_fields.append(masslim_data('cl1350+6007_v0.0.1',  'CL1350',  0.804,  0.035, ['Z', 2], 22.86))

masslim_fields.append(masslim_data('rxj1716+6708_v0.0.8', 'RXJ1716', 0.813,  0.021, ['I', 1], 25.43))
masslim_fields.append(masslim_data('nep5281_v0.0.3',      'N5281',   0.818,  0.029, ['Y', 3], 23.40))
masslim_fields.append(masslim_data('sg0023+0423_v0.2.0',  'SG0023',  0.845,  0.025, ['I', 1], 25.19))
masslim_fields.append(masslim_data('sc1604_v0.0.4',       'SC1604',  0.910,  0.029, ['I', 1], 25.07))
masslim_fields.append(masslim_data('cl1429+4241_v0.0.3',  'CL1429',  0.920,  0.037, ['Y', 3], 23.19))
masslim_fields.append(masslim_data('cl1137+3007_v0.0.1',  'CL1137',  0.959,  0.032, ['Z', 2], 24.70))
masslim_fields.append(masslim_data('cl0910+5422_v0.0.4',  'SC0910',  1.110,  0.035, ['I', 1], 25.76))
masslim_fields.append(masslim_data('rxj1053+5735_v0.0.2', 'RXJ1053', 1.140,  0.031, ['Z', 2], 25.12))
masslim_fields.append(masslim_data('sc0849+4452_v0.0.3',  'SC0849',  1.261,  0.029, ['Z', 2], 25.08))
masslim_fields.append(masslim_data('xlss005_v0.0.3',      'XLSS005', 1.000,  0.024, ['I', 1], 25.82))




outname = '../output/massLimits.pdf'
pdf = PdfPages(outname)

fig = pylab.figure(figsize=(8.2, 6.2))
sp1 = fig.add_subplot(111)

for i in range(len(masslim_fields)):

	f = masslim_fields[i]
	zax = f.data.z

	sp1.grid()
	sp1.minorticks_on()
	sp1.axis([0.52, 1.45, 8.1, 11.9])
	sp1.set_xlabel('redshift', size=22)
	sp1.set_ylabel('log( M$_{\mathrm{limit}}$ / M$_{\odot}$ )', size=22)

	sp1.plot(f.data.z, f.data.masslim95, ls='-', marker='o', color='#666666', ms=8, label='80% completeness')
	sp1.axvline(f.zclust, color='k', ls='--', label='LSS redshift')

	sp1.fill_between([masslim_fields[0].zclust, masslim_fields[-1].zclust], 0, 100, color='#e6e6e6')



	###  calculating and plott masslimit from maximally old SSP
	masses_ssp = model_ssp.get_masses(zform, zs=zax, nfilters=1)
	masses_exp = model_exp.get_masses(zform, zs=zax, nfilters=1)
	mags_ssp = model_ssp.get_apparent_mags(zform, zs=zax)
	mags_exp = model_exp.get_apparent_mags(zform, zs=zax)
	masslim_ssp = pylab.log10(10**((mags_ssp[:,f.filterinfo[1]] - f.maglim) / 2.5) * masses_ssp)
	masslim_exp = pylab.log10(10**((mags_exp[:,f.filterinfo[1]] - f.maglim) / 2.5) * masses_exp)

	sp1.plot(zax, masslim_ssp, color='r', label='SSP: zform=%.1f' % zform)
	sp1.plot(zax, masslim_exp, color='r', ls='--', label=r'SFH exp. $\tau$=1Gyr')


	sp1.legend(loc=2, numpoints = 1, fontsize=16, title='%s: %s-band' % (f.name, f.filterinfo[0]))





	###################################################
	###  writing master mass completeness limits file
	###################################################
	outer = open('%s/masslimits_%s_master.dat' % (data_dir, f.name2), 'w')
	outer.write('#  z  masslim_empirical  masslim_model_ssp  masslim_model_exp\n')
	for mi in range(len(zax)):
		outer.write('  %.3f' % zax[mi])
		outer.write('  %.3f' % f.data.masslim95[mi])
		outer.write('  %.3f' % masslim_ssp[mi])
		outer.write('  %.3f\n' % masslim_exp[mi])
	outer.close()



	pdf.savefig()
	sp1.clear()

pdf.close()
pylab.close()
print '\nwrote to: %s' % outname



















