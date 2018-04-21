
import os
import sys
import mypy
import pylab
import subprocess
from mypy import massfunc
from scipy import optimize
from astropy import wcs, cosmology
from threedhst import eazyPy
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages

font = FontProperties()
font.set_family('sans-serif')

cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

def gunzip_read_gzip(path2file, readcat=0, readzout=0, dtype=float, delimiter=None, comments='#'):

	if readcat == readzout:
		raise IOError('Need to specify either readcat or readzout!')

	print '  reading: %s' % os.path.basename(path2file)
	subprocess.call('gunzip %s' % path2file, shell=1)
	if readcat:
		outer = mypy.readcat(path2file[:-3], dtype=dtype, delimiter=delimiter, comments=comments)
	elif readzout:
		outer = mypy.readzout(path2file[:-3])
	subprocess.call('gzip %s' % path2file[:-3], shell=1)

	return outer

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

	inset = subplot.figure.add_axes([sub_xlo, sub_ylo, sub_dx, sub_dy], transform=subplot.transAxes)
	return inset





data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'

class field:
	def __init__(self, name, version, zclust, sigmaz):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz

		self.cat = gunzip_read_gzip('%s/%s/%s.cat.gz' % (data_dir, version, version), readcat=1)
		self.zout = gunzip_read_gzip('%s/%s/%s.zout.gz' % (data_dir, version, version), readzout=1)
		#self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (data_dir, version, version), readcat=1)
		#self.restframe = gunzip_read_gzip('%s/%s/%s.restframe.gz' % (data_dir, version, version), readcat=1)
		#self.zgrid_pz, self.pzs = eazyPy_tomczak.getEazyPz_all(MAIN_OUTPUT_FILE=version, OUTPUT_DIRECTORY='%s/%s/EAZY' % (data_dir, version))
		print ''


fields = []
fields.append(field('RXJ1757', 'nep200_v0.0.7',       0.691,  0.027))
fields.append(field('SC1324',  'sc1324_v0.0.4',       0.755,  0.033))
fields.append(field('RCS0224', 'rcs0224-0002_v0.0.4', 0.772,  0.027))
fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.8', 0.813,  0.021))
fields.append(field('RXJ1821', 'nep5281_v0.0.4',      0.818,  0.029))
fields.append(field('SC1604',  'sc1604_v0.0.6',       0.910,  0.029))
fields.append(field('SC0910',  'cl0910+5422_v0.0.5',  1.110,  0.035))
fields.append(field('SC0849',  'sc0849+4452_v0.0.4',  1.261,  0.029))
print ''













fig = pylab.figure(figsize=(16.5, 12.7))
fig.subplots_adjust(hspace=0, wspace=0, left=0.08, bottom=0.08)

axis1 = [-0.04, 1.65, -0.04, 1.65]
axis2 = [-0.04, 1.65, -0.28, 0.28]

sp_master          = pylab.subplot2grid((10, 13), (2, 0), rowspan=6, colspan=6)
sp_master_resid    = pylab.subplot2grid((10, 13), (8, 0), rowspan=2, colspan=6)
sp_master_colorbar = add_inset(sp_master, [0.15, 1.16, 0.7, 0.06])

sp_master.minorticks_on()
sp_master_resid.minorticks_on()
sp_master_colorbar.yaxis.set_visible(0)
sp_master_colorbar.set_title('Number')

sp_master.set_ylabel('photometric redshift')
sp_master_resid.set_xlabel('spectroscopic redshift')
sp_master_resid.set_ylabel('$\Delta z$ / (1 + z$_{spec}$)')

sp3 =   pylab.subplot2grid((10, 13), (0, 11), rowspan=2, colspan=2)
sp2 =   pylab.subplot2grid((10, 13), (0, 9),  rowspan=2, colspan=2)
sp1 =   pylab.subplot2grid((10, 13), (0, 7),  rowspan=2, colspan=2)

sp6 =   pylab.subplot2grid((10, 13), (2, 11), rowspan=2, colspan=2)
sp5 =   pylab.subplot2grid((10, 13), (2, 9),  rowspan=2, colspan=2)
sp4 =   pylab.subplot2grid((10, 13), (2, 7),  rowspan=2, colspan=2)

sp9 =   pylab.subplot2grid((10, 13), (4, 11), rowspan=2, colspan=2)
sp8 =   pylab.subplot2grid((10, 13), (4, 9),  rowspan=2, colspan=2)
sp7 =   pylab.subplot2grid((10, 13), (4, 7),  rowspan=2, colspan=2)

sp12 =  pylab.subplot2grid((10, 13), (6, 11), rowspan=2, colspan=2)
sp11 =  pylab.subplot2grid((10, 13), (6, 9),  rowspan=2, colspan=2)
sp10 =  pylab.subplot2grid((10, 13), (6, 7),  rowspan=2, colspan=2)

sp15 =  pylab.subplot2grid((10, 13), (8, 11), rowspan=2, colspan=2)
sp14 =  pylab.subplot2grid((10, 13), (8, 9),  rowspan=2, colspan=2)
sp13 =  pylab.subplot2grid((10, 13), (8, 7),  rowspan=2, colspan=2)

sps = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10, sp11, sp12, sp13, sp14, sp15]
sps = [sp4, sp5, sp7, sp8, sp10, sp11, sp13, sp14, sp1, sp2, sp3, sp6, sp9, sp12, sp15]
for sp in sps:
	sp.axis(axis1)
	sp.tick_params(size=4.5)
	sp.xaxis.set_ticklabels([])
	sp.yaxis.set_ticklabels([])
	sp_master.axis(axis1)
	sp_master_resid.axis(axis2)



zspec_master = []
zphot_master = []
for i in range(len(fields)):

	f = fields[i]
	sp = sps[i]

	inds = pylab.find(f.cat.z_spec > 0)
	zspec = f.cat.z_spec[inds]
	zphot = f.zout.z_peak[inds]
	zspec_master += zspec.tolist()
	zphot_master += zphot.tolist()




zspec_master = pylab.array(zspec_master)
zphot_master = pylab.array(zphot_master)
dz1pluszspec_master = (zphot_master - zspec_master) / (1 + zspec_master)

nmad = mypy.nmad(dz1pluszspec_master)
foutlier = len(pylab.find(abs(dz1pluszspec_master) >= 0.15)) * 1. / len(dz1pluszspec_master)





###  major panel: all fields

nbins = 120
clo, chi = 0, 1.05   ###  colorbar limits
cmap = pylab.cm.spectral
hist2d, xedges, yedges = pylab.histogram2d(zspec_master, zphot_master, bins=(nbins, nbins), range=([axis1[0], axis1[1]], [axis1[2], axis1[3]]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

nlo, nhi = hist2d.min(), hist2d.max()
hist2d[pylab.where(hist2d == 0)] = pylab.nan
hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)])

zphotzspec_map = sp_master.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=cmap)
zphotzspec_map.set_clim(0, chi*zphotzspec_map.get_clim()[1])

sp_master.set_aspect('auto')

#  something more simple
#sp_master.plot(zspec_master, zphot_master, 'ko', ms=1)

label = '%s\nN$_{total}$ = %i\n$\sigma_{NMAD}$ = %.3f\nf$_{outlier}$ = %.3f' % ('All Fields', len(zspec_master), nmad, foutlier)
t = sp_master.text(0.03, 0.75, label, transform=sp_master.transAxes, fontsize=24, color='r',
	               path_effects=[PathEffects.withStroke(linewidth=1.5, foreground='k')])






###  colorbar
cmap_image = pylab.meshgrid(pylab.linspace(0., 1., 101), range(2))[0]
c = sp_master_colorbar.imshow(cmap_image, cmap=cmap)
c.set_clim(clo, chi)
sp_master_colorbar.set_aspect('auto')

xticks = sp_master_colorbar.get_xticks()
sp_master_colorbar.xaxis.set_ticklabels((xticks / 100. * nhi).astype(int))
sp_master_colorbar.xaxis.set_tick_params(size=4.5)





###  residual plot
hist2d, xedges, yedges = pylab.histogram2d(zspec_master, dz1pluszspec_master, bins=(nbins, nbins/3), range=([axis2[0], axis2[1]], [axis2[2], axis2[3]]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

nlo, nhi = hist2d.min(), hist2d.max()
hist2d[pylab.where(hist2d == 0)] = pylab.nan
hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)])

zphotzspec_resid_map = sp_master_resid.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=cmap)
zphotzspec_resid_map.set_clim(0, chi*zphotzspec_resid_map.get_clim()[1])

sp_master_resid.set_aspect('auto')

#l1.set_visible(0), l2.set_visible(0), l3.set_visible(0), l4.set_visible(0), l5.set_visible(0), l6.set_visible(0)

l1 = sp_master_resid.axhline(0, lw=4, color='#b3b3b3')
l2 = sp_master_resid.axhline(0, lw=2, color='r')

l3 = sp_master_resid.axhline(0.15, lw=4, color='#b3b3b3')
l4 = sp_master_resid.axhline(0.15, lw=2, color='r', ls='--')

l5 = sp_master_resid.axhline(-0.15, lw=4, color='#b3b3b3')
l6 = sp_master_resid.axhline(-0.15, lw=2, color='r', ls='--')





###  plotting individual fields in subpanels
for i in range(len(fields)):

	f = fields[i]
	inds = pylab.find(f.cat.z_spec > 0)
	zspec = f.cat.z_spec[inds]
	zphot = f.zout.z_peak[inds]
	dz1pluszspec = (zphot - zspec) / (1 + zspec)

	nmad = mypy.nmad(dz1pluszspec)
	foutlier = len(pylab.find(abs(dz1pluszspec) >= 0.15)) * 1. / len(dz1pluszspec)

	sp = sps[i]

	#sp.plot(zspec, zphot, ls ='', marker='o', mec='gray', ms=1)
	sp.plot(zspec, zphot, ls ='', marker='o', ms=2, mew=0.001, mfc='gray', mec='gray')

	label = '%s\n%i\n%.3f\n%.3f' % (f.name, len(inds), nmad, foutlier)
	t = sp.text(0.05, 0.61, label, transform=sp.transAxes, fontsize=14, color='k',
		        path_effects=[PathEffects.withStroke(linewidth=1., foreground='r')])

















'''
OLD SUBPLOT GRIDS

fig = pylab.figure(figsize=(20, 11))
fig.subplots_adjust(hspace=0, wspace=0, left=0.08, bottom=0.08)
sp_master       = pylab.subplot2grid((8, 15), (0, 0), rowspan=6, colspan=6)
sp_master_resid = pylab.subplot2grid((8, 15), (6, 0), rowspan=2, colspan=6)
sp_master.minorticks_on()
sp_master_resid.minorticks_on()
sp4 =   pylab.subplot2grid((8, 15), (0, 13), rowspan=2, colspan=2)
sp3 =   pylab.subplot2grid((8, 15), (0, 11), rowspan=2, colspan=2)
sp2 =   pylab.subplot2grid((8, 15), (0, 9),  rowspan=2, colspan=2)
sp1 =   pylab.subplot2grid((8, 15), (0, 7),  rowspan=2, colspan=2)
sp8 =   pylab.subplot2grid((8, 15), (2, 13), rowspan=2, colspan=2)
sp7 =   pylab.subplot2grid((8, 15), (2, 11), rowspan=2, colspan=2)
sp6 =   pylab.subplot2grid((8, 15), (2, 9),  rowspan=2, colspan=2)
sp5 =   pylab.subplot2grid((8, 15), (2, 7),  rowspan=2, colspan=2)
sp12 =  pylab.subplot2grid((8, 15), (4, 13), rowspan=2, colspan=2)
sp11 =  pylab.subplot2grid((8, 15), (4, 11), rowspan=2, colspan=2)
sp10 =  pylab.subplot2grid((8, 15), (4, 9),  rowspan=2, colspan=2)
sp9 =   pylab.subplot2grid((8, 15), (4, 7),  rowspan=2, colspan=2)
sp16 =  pylab.subplot2grid((8, 15), (6, 13), rowspan=2, colspan=2)
sp15 =  pylab.subplot2grid((8, 15), (6, 11), rowspan=2, colspan=2)
sp14 =  pylab.subplot2grid((8, 15), (6, 9),  rowspan=2, colspan=2)
sp13 =  pylab.subplot2grid((8, 15), (6, 7),  rowspan=2, colspan=2)


'''







#fig = pylab.figure(figsize=(16.5, 12.7))
fig = pylab.figure(figsize=(10, 8))
fig.subplots_adjust(hspace=0, wspace=0, left=0.08, bottom=0.08)

axis1 = [-0.04, 1.65, -0.04, 1.65]
axis2 = [-0.04, 1.65, -0.28, 0.28]

sp_master          = pylab.subplot2grid(((8, 10)), (1, 0), rowspan=5, colspan=5, aspect=1)
sp_master_resid    = pylab.subplot2grid(((8, 10)), (6, 0), rowspan=2, colspan=5)
sp_master_colorbar = add_inset(sp_master, [0.15, 1.16, 0.7, 0.06])

sp_master.minorticks_on()
sp_master_resid.minorticks_on()
sp_master_colorbar.yaxis.set_visible(0)
sp_master_colorbar.set_title('Number')

sp_master.set_ylabel('photometric redshift')
sp_master_resid.set_xlabel('spectroscopic redshift')
sp_master_resid.set_ylabel('$\Delta z$ / (1 + z$_{spec}$)')

sp2 =  pylab.subplot2grid(((8, 10)), (0, 8),  rowspan=2, colspan=2)
sp1 =  pylab.subplot2grid(((8, 10)), (0, 6),  rowspan=2, colspan=2)

sp4 =  pylab.subplot2grid(((8, 10)), (2, 8),  rowspan=2, colspan=2)
sp3 =  pylab.subplot2grid(((8, 10)), (2, 6),  rowspan=2, colspan=2)

sp6 =  pylab.subplot2grid(((8, 10)), (4, 8),  rowspan=2, colspan=2)
sp5 =  pylab.subplot2grid(((8, 10)), (4, 6),  rowspan=2, colspan=2)

sp8 =  pylab.subplot2grid(((8, 10)), (6, 8),  rowspan=2, colspan=2)
sp7 =  pylab.subplot2grid(((8, 10)), (6, 6),  rowspan=2, colspan=2)















#fig = pylab.figure(figsize=(16.5, 12.7))
#fig = pylab.figure(figsize=(10, 8))
fig = pylab.figure(figsize=(11.42, 9.75))
fig.subplots_adjust(hspace=0, wspace=0, left=0.09, bottom=0.08, top=0.99, right=0.99)

axis1 = [-0.04, 1.7, -0.04, 1.7]
axis2 = [-0.04, 1.7, -0.28, 0.28]

sp_master          = pylab.subplot2grid(((12, 14)), (2, 0), rowspan=7, colspan=7, aspect=1)
sp_master_resid    = pylab.subplot2grid(((12, 14)), (9, 0), rowspan=3, colspan=7)
sp_master_colorbar = add_inset(sp_master, [0.15, 1.16, 0.7, 0.06])

sp_master.minorticks_on()
sp_master.xaxis.set_ticklabels([''])
#sp_master.xaxis.set_ticklabels(['', '0', '0.5', '1.0', '1.5'])
#sp_master.yaxis.set_ticklabels(['', '0', '0.5', '1.0', '1.5'])
sp_master_resid.minorticks_on()
sp_master_colorbar.yaxis.set_visible(0)
sp_master_colorbar.set_title('Number')

sp_master.set_ylabel('photometric redshift')
sp_master_resid.set_xlabel('spectroscopic redshift')
sp_master_resid.set_ylabel('$\Delta z$ / (1 + z$_{spec}$)')

sp2 =  pylab.subplot2grid(((12, 14)), (0, 11),  rowspan=3, colspan=3)
sp1 =  pylab.subplot2grid(((12, 14)), (0, 8),  rowspan=3, colspan=3)

sp4 =  pylab.subplot2grid(((12, 14)), (3, 11),  rowspan=3, colspan=3)
sp3 =  pylab.subplot2grid(((12, 14)), (3, 8),  rowspan=3, colspan=3)

sp6 =  pylab.subplot2grid(((12, 14)), (6, 11),  rowspan=3, colspan=3)
sp5 =  pylab.subplot2grid(((12, 14)), (6, 8),  rowspan=3, colspan=3)

sp8 =  pylab.subplot2grid(((12, 14)), (9, 11),  rowspan=3, colspan=3)
sp7 =  pylab.subplot2grid(((12, 14)), (9, 8),  rowspan=3, colspan=3)


sps = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8]
for sp in sps:
	sp.axis(axis1)
	sp.tick_params(size=4.5)
	sp.set_xlabel('z$_{spec}$', size=15)
	sp.set_ylabel('z$_{phot}$', size=15)
	sp.xaxis.set_ticklabels(['', '0', '0.5', '1.0', '1.5'], size=13)
	sp.yaxis.set_ticklabels(['', '0', '0.5', '1.0', '1.5'], size=13)
	sp_master.axis(axis1)
	sp_master_resid.axis(axis2)





zspec_master = []
zphot_master = []
nmads = pylab.zeros(len(fields))
sigzs = pylab.zeros(len(fields))
for i in range(len(fields)):

	f = fields[i]
	sp = sps[i]

	inds = pylab.find(f.cat.z_spec > 0)
	zspec = f.cat.z_spec[inds]
	zphot = f.zout.z_peak[inds]
	zspec_master += zspec.tolist()
	zphot_master += zphot.tolist()

	dz1pluszspec = (zphot - zspec) / (1 + zspec)
	nmad = mypy.nmad(dz1pluszspec)
	foutlier = len(pylab.find(abs(dz1pluszspec) >= 0.15)) * 1. / len(dz1pluszspec)

	###  measuring scatter and f_outlier
	h = pylab.histogram(dz1pluszspec, range=(-0.5, 0.5), bins=65)
	x, y = (h[1][:-1]+h[1][1:])/2., h[0]
	weights = y**0.5

	### replacing weights of 0 -> 1
	winds = pylab.find(weights == 0)
	weights[winds] = 1.

	fit, cov = optimize.curve_fit(mypy.gauss1d, x, y, sigma=1./weights)

	nmads[i] = abs(nmad)
	sigzs[i] = abs(fit[2])



zspec_master = pylab.array(zspec_master)
zphot_master = pylab.array(zphot_master)
dz1pluszspec_master = (zphot_master - zspec_master) / (1 + zspec_master)
zinds = pylab.find((zspec_master > 0.5) & (zspec_master < 1.3))

###  measuring scatter and f_outlier
h = pylab.histogram(dz1pluszspec_master, range=(-0.5, 0.5), bins=65)
x, y = (h[1][:-1]+h[1][1:])/2., h[0]
weights = y**0.5

### replacing weights of 0 -> 1
winds = pylab.find(weights == 0)
weights[winds] = 1.

fit, cov = optimize.curve_fit(mypy.gauss1d, x, y, sigma=1./weights)



nmad = mypy.nmad(dz1pluszspec_master)
foutlier = len(pylab.find(abs(dz1pluszspec_master) >= 0.15)) * 1. / len(dz1pluszspec_master)





###  major panel: all fields

nbins = 120
clo, chi = 0, 1.05   ###  colorbar limits
cmap = pylab.cm.spectral
hist2d, xedges, yedges = pylab.histogram2d(zspec_master, zphot_master, bins=(nbins, nbins), range=([axis1[0], axis1[1]], [axis1[2], axis1[3]]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

nlo, nhi = hist2d.min(), hist2d.max()
hist2d[pylab.where(hist2d == 0)] = pylab.nan
hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)])

zphotzspec_map = sp_master.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=cmap)
#zphotzspec_map = sp_master.imshow(hist2d.T[::-1], extent=extent, interpolation='nearest', cmap=cmap)
zphotzspec_map.set_clim(0, chi*zphotzspec_map.get_clim()[1])

sp_master.set_aspect('auto')

#  something more simple
#sp_master.plot(zspec_master, zphot_master, 'ko', ms=1)

label = '%s\nN$_{total}$ = %i\n$\sigma_{NMAD}$ = %.3f\nf$_{outlier}$ = %.3f' % ('All Fields', len(zspec_master), nmad, foutlier)
label = '%s\nN$_{total}$ = %i\n$\sigma_{\Delta z_{spec}}$ = %.3f\nf$_{outlier}$ = %.3f' % ('All Fields', len(zspec_master), abs(fit[2]), foutlier)
#t = sp_master.text(0.03, 0.97, label, transform=sp_master.transAxes, fontsize=24, color='r',
#	               horizontalalignment='left', verticalalignment='top',
#	               path_effects=[PathEffects.withStroke(linewidth=1.5, foreground='k')])
t0 = sp_master.text(0.03, 0.97, label, transform=sp_master.transAxes, fontsize=22, color='k',
	               horizontalalignment='left', verticalalignment='top', fontproperties=font,
	               path_effects=[PathEffects.withStroke(linewidth=2, foreground='r')])






###  colorbar
cmap_image = pylab.meshgrid(pylab.linspace(0., 1., 101), range(2))[0]
c = sp_master_colorbar.imshow(cmap_image, cmap=cmap)
c.set_clim(clo, chi)
sp_master_colorbar.set_aspect('auto')

xticks = sp_master_colorbar.get_xticks()
sp_master_colorbar.xaxis.set_ticklabels((xticks / 100. * nhi).astype(int))
sp_master_colorbar.xaxis.set_tick_params(size=4.5)





###  residual plot
hist2d, xedges, yedges = pylab.histogram2d(zspec_master, dz1pluszspec_master, bins=(nbins, nbins/3), range=([axis2[0], axis2[1]], [axis2[2], axis2[3]]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

nlo, nhi = hist2d.min(), hist2d.max()
hist2d[pylab.where(hist2d == 0)] = pylab.nan
hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)])

zphotzspec_resid_map = sp_master_resid.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=cmap)
zphotzspec_resid_map.set_clim(0, chi*zphotzspec_resid_map.get_clim()[1])

sp_master_resid.set_aspect('auto')

#l1.set_visible(0), l2.set_visible(0), l3.set_visible(0), l4.set_visible(0), l5.set_visible(0), l6.set_visible(0)

l1 = sp_master_resid.axhline(0, lw=4, color='#b3b3b3')
l2 = sp_master_resid.axhline(0, lw=2, color='r')

l3 = sp_master_resid.axhline(0.15, lw=4, color='#b3b3b3')
l4 = sp_master_resid.axhline(0.15, lw=2, color='r', ls='--')

l5 = sp_master_resid.axhline(-0.15, lw=4, color='#b3b3b3')
l6 = sp_master_resid.axhline(-0.15, lw=2, color='r', ls='--')





###  plotting individual fields in subpanels
for i in range(len(fields)):

	f = fields[i]
	inds = pylab.find(f.cat.z_spec > 0)
	zspec = f.cat.z_spec[inds]
	zphot = f.zout.z_peak[inds]
	dz1pluszspec = (zphot - zspec) / (1 + zspec)

	nmad = mypy.nmad(dz1pluszspec)
	foutlier = len(pylab.find(abs(dz1pluszspec) >= 0.15)) * 1. / len(dz1pluszspec)

	sp = sps[i]

	#sp.plot(zspec, zphot, ls ='', marker='o', mec='gray', ms=1)
	sp.plot(zspec, zphot, ls ='', marker='o', ms=2, mew=0.001, mfc='gray', mec='gray')

	label = '%s\n%i\n%.3f\n%.3f' % (f.name, len(inds), sigzs[i], foutlier)
	t = sp.text(0.05, 0.97, label, transform=sp.transAxes, fontsize=15, color='k',
	            horizontalalignment='left', verticalalignment='top', fontproperties=font,
		        path_effects=[PathEffects.withStroke(linewidth=1.5, foreground='r')])




























'''
###  calculating stats for each field separately

outname = '../output/asdf_zphotzspec.pdf'
pdf = PdfPages(outname)


nmads = pylab.zeros(len(fields))
sigzs = pylab.zeros(len(fields))
for i in range(len(fields)):

	f = fields[i]
	inds = pylab.find(f.cat.z_spec > 0)
	zspec = f.cat.z_spec[inds]
	zphot = f.zout.z_peak[inds]
	dz1pluszspec = (zphot - zspec) / (1 + zspec)
	zinds = pylab.find((zspec > 0.5) & (zspec < 1.2))


	nmad = mypy.nmad(dz1pluszspec)
	foutlier = len(pylab.find(abs(dz1pluszspec) >= 0.15)) * 1. / len(dz1pluszspec)


	###  measuring scatter and f_outlier
	h = pylab.histogram(dz1pluszspec, range=(-0.5, 0.5), bins=65)
	x, y = (h[1][:-1]+h[1][1:])/2., h[0]
	weights = y**0.5

	### replacing weights of 0 -> 1
	winds = pylab.find(weights == 0)
	weights[winds] = 1.

	fit, cov = optimize.curve_fit(mypy.gauss1d, x, y, sigma=1./weights)
	#print '%.1f' % (fit[2] * 100)
	print '%.1f' % (nmad * 100)

	nmads[i] = abs(nmad)
	sigzs[i] = abs(fit[2])



	pylab.hist(dz1pluszspec, histtype='step', lw=2, color='k', range=(-0.4, 0.4), bins=int(len(dz1pluszspec)**0.5*2))

	xmod = pylab.linspace(-0.4, 0.4, 1000)
	ymod = mypy.gauss1d(xmod, *fit)
	pylab.plot(xmod, ymod, 'r')

	pylab.grid()
	pylab.xlim(-0.4, 0.4)

	pdf.savefig()
	pylab.close()

pdf.close()
print '\n!open %s' % outname
'''











