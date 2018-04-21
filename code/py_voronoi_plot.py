
import mypy
import pylab
import matplotlib
from astropy import wcs
from scipy import signal, ndimage
from astropy.io import fits
from matplotlib import collections
import matplotlib.patches as patches
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.patheffects as PathEffects
from matplotlib.font_manager import FontProperties

font = FontProperties()
font.set_family('sans-serif')


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



voronoi_dir = '/Users/atomczak/DATA/ORELSE/Voronoi'
voronoi_dir = '/Users/atomczak/GoogleDrive/ORELSE/Voronoi_Maps/Number_Density_and_Overdensity_Maps'
voronoi_dir = '/Volumes/SASQUATCH/Backups.backupdb/ukulele/2017-04-08-000910/Macintosh HD/Users/atomczak/DATA/ORELSE/Voronoi'
data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
version = 'rxj1716+6708_v0.0.8'


rxj1716_voronoi_maps = fits.open('%s/RXJ1716.mastermedVoronoi.overdens.100iterations.fits' % voronoi_dir)
ind = 18   #  index of voronoi map centered on zclust
voronoi_map = rxj1716_voronoi_maps[ind].data
voronoi_wcs = wcs.WCS(rxj1716_voronoi_maps[18].header)
z1, z2 = rxj1716_voronoi_maps[ind].header['Z1'], rxj1716_voronoi_maps[ind].header['Z2']

xgrid = pylab.arange(0., voronoi_map.shape[0], 1)
ygrid = pylab.arange(0., voronoi_map.shape[1], 1)
ragrid, decgrid = voronoi_wcs.wcs_pix2world(xgrid, ygrid, 1)

kernel = mypy.gauss2d(3, 3, 1, 1, 1, 0.5, 0.5)
kernel /= kernel.sum()
voronoi_map = signal.convolve2d(voronoi_map, kernel, mode='same')


cat = mypy.readcat('%s/%s/%s.cat.gz' % (data_dir, version, version))
fout = mypy.readcat('%s/%s/%s.fout.gz' % (data_dir, version, version))
zout = mypy.readzout('%s/%s/%s.zout.gz' % (data_dir, version, version))











###  coordinates and size of large box (size in arcmin)
ra0, dec0, box_size_amin, cutout_size_amin = 259.1919, 67.147616, 14., 2.5
x0, y0 = voronoi_wcs.wcs_world2pix([ra0], [dec0], 1)
x0, y0 = int(x0 + 0.5), int(y0 + 0.5)

###  pixel scale of voronoi map
xs, ys = voronoi_wcs.wcs_pix2world([0, 0], [0, 1], 1)
px_scale = mypy.radec_sep(xs[0], ys[0], xs[1], ys[1])
box_size_px = int(box_size_amin * 60 / px_scale)






fig = pylab.figure(figsize=(18.1125, 8.35))

sp2 = pylab.subplot2grid((5, 12), (1, 3), rowspan=3, colspan=3, aspect=1)
sp1 = pylab.subplot2grid((5, 12), (1, 0), rowspan=3, colspan=3, aspect=1, title='Voronoi Tesselation')
sp3 = pylab.subplot2grid((5, 12), (0, 7), rowspan=5, colspan=5, aspect=1)

fig.subplots_adjust(wspace=0.1, left=0.02)
sp1.xaxis.set_visible(0)
sp1.yaxis.set_visible(0)
sp2.xaxis.set_visible(0)
sp2.yaxis.set_visible(0)

sp3.set_xlabel('R.A.$_{\mathrm{J2000}}$ [$\degree$]')
sp3.set_ylabel('Dec.$_{\mathrm{J2000}}$ [$\degree$]')

sp2.set_title('Overdensity Map')

###  plotting large voronoi map
m1 = voronoi_map[y0-box_size_px/2:y0+box_size_px/2, x0-box_size_px/2:x0+box_size_px/2]
vmin, vmax = m1.min(), m1.max()
d = sp3.imshow(m1, cmap=pylab.cm.rainbow, vmin=vmin, vmax=vmax, interpolation='nearest')

m1_side = m1.shape[0]
sp3.axis([-0.5, -0.5+m1_side, -0.5, -0.5+m1_side])


info =  'RXJ1716\n'
info += '%.3f < z < %.3f' % (z1, z2)
t0 = sp3.text(0.6, 0.88, info, transform=sp3.transAxes, fontweight='bold', color='w', fontproperties=font, fontsize=21,
	    	 path_effects=[PathEffects.withStroke(linewidth=4., foreground='k')])





###  plotting 5' scale bar
x1 = 5
x2 = x1 + (5 * 60) / px_scale

sp3.plot([x1, x2], [x1, x1], color='k', lw=7)
sp3.plot([x1, x2], [x1, x1], color='w', lw=3)

a = sp3.text((x1+x2)/2., x1+1, "5'", fontsize=22, fontproperties=font, color='w', 
	         horizontalalignment='center', verticalalignment='bottom',
	    	 path_effects=[PathEffects.withStroke(linewidth=4., foreground='k')])






###  plotting smaller voronoi cutout
dx = int(box_size_px * cutout_size_amin / box_size_amin + 0.5)
m2 = voronoi_map[y0-dx/2:y0+dx/2, x0-dx/2:x0+dx/2]
sp2.imshow(m2, cmap=pylab.cm.rainbow, vmin=vmin, vmax=vmax, interpolation='nearest')

m2_side = m2.shape[0]
sp2.axis([-0.5, -0.5+m2_side, -0.5, -0.5+m2_side])



square = patches.Rectangle((box_size_px/2-dx/2-1, box_size_px/2-dx/2-1), dx, dx, \
	                       fill=0, lw=3, color='k')
sp3.add_patch(square)










###  galaxies
mags = 25 - 2.5 * pylab.log10(cat.fluxauto_i)
inds = pylab.find((cat.ra > ra0-cutout_size_amin) & (cat.ra < ra0+cutout_size_amin) & 
                  (cat.dec > dec0-cutout_size_amin) & (cat.dec < dec0+cutout_size_amin) & 
                  (cat.use == 1) & (cat.z_spec > z1) & (cat.z_spec < z2) & (mags > 18) & (mags < 24.5) & (zout.odds > 0.8))
#                  (cat.use == 1) & (fout.z > z1) & (fout.z < z2) & (mags > 18) & (mags < 24.5) & (zout.odds > 0.8))


x_gals, y_gals = voronoi_wcs.wcs_world2pix(cat.ra[inds], cat.dec[inds], 1)
#x_gals /= pylab.cos(dec0 * pylab.pi / 180.)
x_gals -= (x0 - dx/2 + 1)
y_gals -= (y0 - dx/2 + 1)

sp1.plot(x_gals, y_gals, 'ko', ms=3)
sp2.plot(x_gals, y_gals, 'ko', ms=3)


###  plotting voronoi grid
points = pylab.array([[x_gals[i], y_gals[i]] for i in range(len(x_gals))])
vor = Voronoi(points)
sp1.axis([-0.5, -0.5+m2_side, -0.5, -0.5+m2_side])

for reg in vor.regions:

	if len(reg) > 0:
		reg.append(reg[0])
		inds = pylab.find(pylab.array(reg) > -1)
		for i in range(len(inds) - 1):

			xy1, xy2 = vor.vertices[reg[inds[i]]], vor.vertices[reg[inds[i+1]]]
			sp1.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], lw=1, color='gray')



transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(sp3.transData.transform([box_size_px/2 - dx/2-1, box_size_px/2 + dx/2]))
coord2 = transFigure.transform(sp2.transData.transform([-0.5+m2_side, -0.5+m2_side]))
coord3 = transFigure.transform(sp3.transData.transform([box_size_px/2 - dx/2-1, box_size_px/2 - dx/2-1]))
coord4 = transFigure.transform(sp2.transData.transform([-0.5+m2_side, -0.5]))

line1 = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                               transform=fig.transFigure, color='k', lw=3)
line2 = matplotlib.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),
                               transform=fig.transFigure, color='k', lw=3)
fig.lines = [line1, line2]




###  setting axis values for larger plot

ra_labels0 = pylab.arange(round(ragrid.min(), 1), round(ragrid.max(), 1), 0.1)
x_labels0 = pylab.interp(ra_labels0, ragrid[::-1], xgrid[::-1])
x_labels0 -= (x0 - box_size_px/2)


dec_labels0 = pylab.arange(round(decgrid.min(), 1), round(decgrid.max(), 1), 0.05)
y_labels0 = pylab.interp(dec_labels0, decgrid, ygrid)
y_labels0 -= (y0 - box_size_px/2)


pylab.xticks(x_labels0, ra_labels0)
pylab.yticks(y_labels0, dec_labels0)
sp3.axis([-0.5, -0.5+m1_side, -0.5, -0.5+m1_side])






###  colorbar

sp4 = add_inset(sp2, rect=[0.1, -0.27, 0.8, 0.08])
sp4.yaxis.set_visible(0)
sp4.tick_params(length=4)
sp4.set_xlabel('log(1 + $\delta_{\mathrm{gal}}$)')

n = 100
cbar = pylab.meshgrid(pylab.linspace(vmin, vmax, n), [0, 0])[0]
sp4.imshow(cbar, cmap=pylab.cm.rainbow)
sp4.set_aspect('auto')

dgrid = pylab.arange(round(vmin, 0), round(vmax, 0), 0.6)
x_labels0 = pylab.interp(dgrid, pylab.linspace(vmin, vmax, n), pylab.arange(n))
pylab.xticks(x_labels0, dgrid)








































###########################
###  repeating for sc1604
###########################

voronoi_dir = '/Users/atomczak/DATA/ORELSE/Voronoi'
voronoi_dir = '/Users/atomczak/GoogleDrive/ORELSE/Voronoi_Maps/Number_Density_and_Overdensity_Maps'
voronoi_dir = '/Volumes/SASQUATCH/Backups.backupdb/ukulele/2017-04-08-000910/Macintosh HD/Users/atomczak/DATA/ORELSE/Voronoi'
data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
version = 'sc1604_v0.0.6'


sc1604_voronoi_maps = fits.open('%s/SC1604.mastermedVoronoi.overdens.100iterations.fits' % voronoi_dir)
ind = 49   #  index of voronoi map centered on zclust
voronoi_map = sc1604_voronoi_maps[ind].data
voronoi_wcs = wcs.WCS(sc1604_voronoi_maps[49].header)
z1, z2 = sc1604_voronoi_maps[ind].header['Z1'], sc1604_voronoi_maps[ind].header['Z2']

xgrid = pylab.arange(0., voronoi_map.shape[0], 1)
ygrid = pylab.arange(0., voronoi_map.shape[0], 1)
ragrid, decgrid = voronoi_wcs.wcs_pix2world(xgrid, ygrid, 1)

kernel = mypy.gauss2d(3, 3, 1, 1, 1, 0.5, 0.5)
kernel /= kernel.sum()
voronoi_map = signal.convolve2d(voronoi_map, kernel, mode='same')


cat = mypy.readcat('%s/%s/%s.cat.gz' % (data_dir, version, version))
fout = mypy.readcat('%s/%s/%s.fout.gz' % (data_dir, version, version))
zout = mypy.readzout('%s/%s/%s.zout.gz' % (data_dir, version, version))















###  coordinates and size of large box (size in arcmin)
ra0, dec0, box_size_amin, cutout_size_amin = 241.11382, 43.323246, 14., 2.5
ra0, dec0, box_size_amin, cutout_size_amin = 241.14, 43.345, 14., 2.5
x0, y0 = voronoi_wcs.wcs_world2pix([ra0], [dec0], 1)
x0, y0 = int(x0 + 0.5), int(y0 + 0.5)

###  pixel scale of voronoi map
xs, ys = voronoi_wcs.wcs_pix2world([0, 0], [0, 1], 1)
px_scale = mypy.radec_sep(xs[0], ys[0], xs[1], ys[1])
box_size_px = int(box_size_amin * 60 / px_scale)






fig = pylab.figure(figsize=(18.1125, 8.35))

sp2 = pylab.subplot2grid((5, 12), (1, 3), rowspan=3, colspan=3, aspect=1)
sp1 = pylab.subplot2grid((5, 12), (1, 0), rowspan=3, colspan=3, aspect=1, title='Voronoi Tesselation')
sp3 = pylab.subplot2grid((5, 12), (0, 7), rowspan=5, colspan=5, aspect=1)

fig.subplots_adjust(wspace=0.1, left=0.02)
sp1.xaxis.set_visible(0)
sp1.yaxis.set_visible(0)
sp2.xaxis.set_visible(0)
sp2.yaxis.set_visible(0)

sp3.set_xlabel('R.A.$_{\mathrm{J2000}}$ [$\degree$]')
sp3.set_ylabel('Dec.$_{\mathrm{J2000}}$ [$\degree$]')

sp2.set_title('Overdensity Map')

###  plotting large voronoi map
m1 = voronoi_map[y0-box_size_px/2:y0+box_size_px/2, x0-box_size_px/2:x0+box_size_px/2]
vmin, vmax = m1.min(), m1.max()
d = sp3.imshow(m1, cmap=pylab.cm.rainbow, vmin=vmin, vmax=vmax, interpolation='nearest')

m1_side = m1.shape[0]
sp3.axis([-0.5, -0.5+m1_side, -0.5, -0.5+m1_side])


info =  'Cl1604\nCluster D + Filament\n'
info += '%.3f < z < %.3f' % (z1, z2)
t0 = sp3.text(0.55, 0.97, info, transform=sp3.transAxes, fontweight='bold', color='w', fontproperties=font, fontsize=21,
	    	  verticalalignment='top', path_effects=[PathEffects.withStroke(linewidth=4., foreground='k')])





###  plotting 5' scale bar
x1 = 5
x2 = x1 + (5 * 60) / px_scale

sp3.plot([x1, x2], [x1, x1], color='k', lw=7)
sp3.plot([x1, x2], [x1, x1], color='w', lw=3)

a = sp3.text((x1+x2)/2., x1+1, "5'", fontsize=22, fontproperties=font, color='w', 
	         horizontalalignment='center', verticalalignment='bottom',
	    	 path_effects=[PathEffects.withStroke(linewidth=4., foreground='k')])






###  plotting smaller voronoi cutout
dx = int(box_size_px * cutout_size_amin / box_size_amin + 0.5)
m2 = voronoi_map[y0-dx/2:y0+dx/2, x0-dx/2:x0+dx/2]
sp2.imshow(m2, cmap=pylab.cm.rainbow, vmin=vmin, vmax=vmax, interpolation='nearest')

m2_side = m2.shape[0]
sp2.axis([-0.5, -0.5+m2_side, -0.5, -0.5+m2_side])



square = patches.Rectangle((box_size_px/2-dx/2-1, box_size_px/2-dx/2-1), dx, dx, \
	                       fill=0, lw=3, color='k')
sp3.add_patch(square)










###  galaxies
mags = 25 - 2.5 * pylab.log10(cat.fluxauto_i)
inds = pylab.find((cat.ra > ra0-cutout_size_amin) & (cat.ra < ra0+cutout_size_amin) & 
                  (cat.dec > dec0-cutout_size_amin) & (cat.dec < dec0+cutout_size_amin) & 
                  (cat.use == 1) & (cat.z_spec > z1) & (cat.z_spec < z2) & (mags > 18) & (mags < 24.5) & (zout.odds > 0.8))
#                  (cat.use == 1) & (fout.z > z1) & (fout.z < z2) & (mags > 18) & (mags < 24.5) & (zout.odds > 0.8))


x_gals, y_gals = voronoi_wcs.wcs_world2pix(cat.ra[inds], cat.dec[inds], 1)
#x_gals /= pylab.cos(dec0 * pylab.pi / 180.)
x_gals -= (x0 - dx/2 + 1)
y_gals -= (y0 - dx/2 + 1)

sp1.plot(x_gals, y_gals, 'ko', ms=3)
sp2.plot(x_gals, y_gals, 'ko', ms=3)


###  plotting voronoi grid
points = pylab.array([[x_gals[i], y_gals[i]] for i in range(len(x_gals))])
vor = Voronoi(points)
sp1.axis([-0.5, -0.5+m2_side, -0.5, -0.5+m2_side])

for reg in vor.regions:

	if len(reg) > 0:
		reg.append(reg[0])
		inds = pylab.find(pylab.array(reg) > -1)
		for i in range(len(inds) - 1):

			xy1, xy2 = vor.vertices[reg[inds[i]]], vor.vertices[reg[inds[i+1]]]
			sp1.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], lw=1, color='gray')



transFigure = fig.transFigure.inverted()
coord1 = transFigure.transform(sp3.transData.transform([box_size_px/2 - dx/2-1, box_size_px/2 + dx/2.25]))
coord2 = transFigure.transform(sp2.transData.transform([-0.5+m2_side, -0.5+m2_side]))
coord3 = transFigure.transform(sp3.transData.transform([box_size_px/2 - dx/2-1, box_size_px/2 - dx/2-1]))
coord4 = transFigure.transform(sp2.transData.transform([-0.5+m2_side, -0.5]))

line1 = matplotlib.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),
                               transform=fig.transFigure, color='k', lw=3)
line2 = matplotlib.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),
                               transform=fig.transFigure, color='k', lw=3)
fig.lines = [line1, line2]




###  setting axis values for larger plot

ra_labels0 = pylab.arange(round(ragrid.min(), 1), round(ragrid.max(), 1), 0.1)
x_labels0 = pylab.interp(ra_labels0, ragrid[::-1], xgrid[::-1])
x_labels0 -= (x0 - box_size_px/2)


dec_labels0 = pylab.arange(round(decgrid.min(), 1), round(decgrid.max(), 1), 0.05)
y_labels0 = pylab.interp(dec_labels0, decgrid, ygrid)
y_labels0 -= (y0 - box_size_px/2)


pylab.xticks(x_labels0, ra_labels0)
pylab.yticks(y_labels0, dec_labels0)
sp3.axis([-0.5, -0.5+m1_side, -0.5, -0.5+m1_side])






###  colorbar

sp4 = add_inset(sp2, rect=[0.1, -0.27, 0.8, 0.08])
sp4.yaxis.set_visible(0)
sp4.tick_params(length=4)
sp4.set_xlabel('log(1 + $\delta_{\mathrm{gal}}$)')

n = 100
cbar = pylab.meshgrid(pylab.linspace(vmin, vmax, n), [0, 0])[0]
sp4.imshow(cbar, cmap=pylab.cm.rainbow)
sp4.set_aspect('auto')

dgrid = pylab.arange(round(vmin, 0), round(vmax, 0), 0.6)
x_labels0 = pylab.interp(dgrid, pylab.linspace(vmin, vmax, n), pylab.arange(n))
pylab.xticks(x_labels0, dgrid)












