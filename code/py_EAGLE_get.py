
import mypy
import pylab
import numpy as np
import eagleSqlTools as sql
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import wcs, cosmology, constants, units

cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

sim = ['RefL0100N1504', 100.]
connection = sql.connect("atomczak", password="inTiHGVg")


###  storing basic meta data of EAGLE snapshots
class snapshot:
	def __init__(self, id, expansionFactor, lookbackTime, redshift):
		self.id = id
		self.expansionFactor = expansionFactor
		self.lookbackTime = lookbackTime
		self.redshift = redshift
		self.z = redshift

snapshots_info = [[0, 0.05, 13.62, 20.00],
                 [1, 0.06, 13.53, 15.13],
                 [2, 0.09, 13.32, 9.99],
                 [3, 0.10, 13.25, 8.99],
                 [4, 0.11, 13.16, 8.07],
                 [5, 0.12, 13.04, 7.05],
                 [6, 0.14, 12.86, 5.97],
                 [7, 0.15, 12.75, 5.49],
                 [8, 0.17, 12.63, 5.04],
                 [9, 0.18, 12.46, 4.49],
                 [10, 0.20, 12.25, 3.98],
                 [11, 0.22, 12.01, 3.53],
                 [12, 0.25, 11.66, 3.02],
                 [13, 0.29, 11.16, 2.48],
                 [14, 0.31, 10.86, 2.24],
                 [15, 0.33, 10.53, 2.01],
                 [16, 0.37, 10.05, 1.74],
                 [17, 0.40, 9.49, 1.49],
                 [18, 0.44, 8.86, 1.26],
                 [19, 0.50, 7.93, 1.00],
                 [20, 0.54, 7.37, 0.87],
                 [21, 0.58, 6.71, 0.74],
                 [22, 0.62, 6.01, 0.62],
                 [23, 0.67, 5.19, 0.50],
                 [24, 0.73, 4.16, 0.37],
                 [25, 0.79, 3.23, 0.27],
                 [26, 0.85, 2.29, 0.18],
                 [27, 0.91, 1.34, 0.10],
                 [28, 1.00, 0.00, 0.001]]
snapshots = [snapshot(*info) for info in snapshots_info]






###  plotting galaxy positions (x vs. y) for each snapshot
pdf = PdfPages('../output/eagle_simulations/sky_positions.pdf')

fig = pylab.figure(figsize=(12., 11.))
sp = fig.add_subplot(111)

fig.subplots_adjust(left=0.09, bottom=0.09, top=0.95, right=0.95)

for si in range(len(snapshots)):

	mypy.progress_bar(si, len(snapshots))
	proper2comoving = (cosmo.kpc_proper_per_arcmin(snapshots[si].redshift) / cosmo.kpc_comoving_per_arcmin(snapshots[si].redshift)).value

	###  grabbing galaxy data
	snap = snapshots[si]
	myQuery = "SELECT \
	            SH.GalaxyID as GalaxyID, \
	            SH.LastProgID as LastProgID, \
	            SH.TopLeafID as TopLeafID, \
	            SH.DescendantID as DescendantID, \
	            SH.CentreOfMass_x as x, \
	            SH.CentreOfMass_y as y, \
	            AP.Mass_Star as stellar_mass, \
	            AP.SFR as sfr \
	                   FROM \
	            %s_Aperture as AP, \
	            %s_SubHalo as SH \
	                   WHERE \
	            SH.GalaxyID = AP.GalaxyID and \
	            AP.Mass_Star > 1e8 and \
	            AP.ApertureSize = 30 and \
	            SH.SnapNum = %i " % (sim[0], sim[0], snap.id)

	galaxyData = sql.execute_query(connection, myQuery)
	pylab.savetxt('../data/eagle_simulations/snap%02i.dat' % snap.id, galaxyData, 
		          fmt='%10i  %10i  %10i  %10i  %7.3f  %7.3f  %.3e  %.1e', 
		          header=' GalaxyID  LastProgID  TopLeafID  DescendantID  x_cMpc  y_cMpc  stellarMass  sfr')
	galaxyData = galaxyData.reshape(-1)


	###  grabbing group data
	myQuery = "SELECT \
				  FOF.GroupID as id, \
				  FOF.GroupCentreOfPotential_x as x, \
				  FOF.GroupCentreOfPotential_y as y, \
				  FOF.GroupMass as Mgroup, \
				  FOF.Group_M_Crit200 as M200, \
				  FOF.Group_R_Crit200 as R200 \
	           FROM \
				  %s_FOF as FOF \
	           WHERE \
				  FOF.GroupMass > 1e14 and \
				  FOF.SnapNum = %i " % (sim[0], snap.id)

	groupData = sql.execute_query(connection, myQuery)
	groupData = groupData.reshape(-1)




	###  plotting
	sp.minorticks_on()
	sp.axis([0, 100, 0, 100])
	sp.set_xlabel('x position [cMpc]')
	sp.set_ylabel('y position [cMpc]')

	sp.set_title('EAGLE simulation :    z = %.2f ,    log(M) > 9' % snap.redshift)


	###  plotting galaxies at >10**9 Msol
	pinds = pylab.find(galaxyData['stellar_mass'] > 10**9)
	sp.plot(galaxyData['x'][pinds], galaxyData['y'][pinds], 'ko', ms=1)

	#nbins = 200
	#hist2d, xedges, yedges = pylab.histogram2d(galaxyData['x'], galaxyData['y'], bins=(nbins, nbins), range=([0, 100], [0, 100]))
	#extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	#asdf = sp.imshow(hist2d.T, extent=extent, interpolation='nearest', cmap=pylab.cm.gray_r)
	#asdf.set_clim(0, 50)



	###  plotting groups/clusters at >10**14 Msol
	for gi in range(len(groupData)):
		if 10**14 < groupData['M200'][gi] < 10**14.5:
			circ = pylab.Circle((groupData['x'][gi], groupData['y'][gi]), groupData['R200'][gi] / proper2comoving / 1000, color='#00ccff')
			sp.add_patch(circ)
		if 10**14.5 <= groupData['M200'][gi]:
			circ = pylab.Circle((groupData['x'][gi], groupData['y'][gi]), groupData['R200'][gi] / proper2comoving / 1000, color='#ff8080')
			sp.add_patch(circ)

	sp.plot(-10, -10, 'ko', ms=22, mew=1, mfc='#ff8080', label='log(M$_{200}$) > 14.5')
	sp.plot(-10, -10, 'ko', ms=22, mew=1, mfc='#00ccff', label='log(M$_{200}$) > 14')


	sp.legend(loc=4, numpoints=1, fontsize=20)

	pdf.savefig()
	sp.clear()


pdf.close()
pylab.close()
























'''

myQuery = "SELECT \
			  FOF.GroupID as id, \
			  FOF.GroupCentreOfPotential_x as x, \
			  FOF.GroupCentreOfPotential_y as y \
           FROM \
			  %s_FOF as FOF \
           WHERE \
			  FOF.Group_M_Crit200 > 1e14 and \
			  FOF.SnapNum = 21 "%(sim[0])

groupData = sql.execute_query(con, myQuery)




myQuery = "SELECT \
            AP.GalaxyID as id, \
            SH.CentreOfMass_x as x, \
            SH.CentreOfMass_y as y, \
            SH.Redshift as redshift, \
			AP.Mass_Star as mass \
                   FROM \
			%s_SubHalo as SH, \
			%s_Aperture as AP \
                   WHERE \
			SH.GalaxyID = AP.GalaxyID and \
			AP.ApertureSize =  30 and \
			AP.Mass_Star > 1e8 and \
			SH.SnapNum = 21 "%(sim[0], sim[0])

galaxyData = sql.execute_query(con, myQuery)


'''









