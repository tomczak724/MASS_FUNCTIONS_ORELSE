
#  This script generates catalogs tabulating values measured
#  from Brian's Voronoi maps of galaxies in the photometric
#  catalogs. These new catalogs are line-matched.

import os
import sys
import math
import mypy
import glob
import time
import numpy
import pylab
import pickle
import subprocess
from astropy import wcs
from scipy import optimize
from astropy.io import fits
from matplotlib.path import Path
import shapely.geometry as geometry
import matplotlib.patheffects as PathEffects
from scipy.spatial import ConvexHull, Delaunay
from matplotlib.font_manager import FontProperties
from shapely.ops import cascaded_union, polygonize
from astropy import wcs, cosmology, constants, units
from matplotlib.backends.backend_pdf import PdfPages







#################################################
###  Defining functions and reading catalogs  ###
#################################################

if True:

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

    def uvj_select(uv, vj, z, uvj_slope=0.88, uvj_intercept=0.59):
        """
        #  INPUT:
        #    uv  --->  Restframe (U-V) color (in AB system)  [array-like]
        #    vj  --->  Restframe (V-J) color (in AB system)  [array-like]
        #    z  ---->  Galaxy redshifts  [array-like]
        #
        #  OUTPUT: quiescent=0, star-forming=1
        #
        #  EXAMPLE:
        #            [ SF, SF, Qui, SF,  ....  Qui, Qui]
        #     array( [ 1,  1,  0,   1,   ....  0,   0 ] )
        #
        #  based on Whitaker+2011
        """
        if type(uv)==list or type(uv)==float or type(uv)==int or type(uv)==long or type(uv)==numpy.float64: uv = numpy.array(uv)
        if type(vj)==list or type(vj)==float or type(vj)==int or type(vj)==long or type(vj)==numpy.float64: vj = numpy.array(vj)

        floor = pylab.zeros(len(z))
        wall = pylab.zeros(len(z))
        slope = pylab.zeros(len(z))
        intercept = pylab.zeros(len(z))

        floor[ pylab.where( (0.0<=z) & (z<1.5) ) ] = 1.3
        floor[ pylab.where( (1.5<=z) & (z<2.0) ) ] = 1.3
        floor[ pylab.where( (2.0<=z) ) ] = 1.2

        wall[ pylab.where( (0.0<=z) & (z<1.5) ) ] = 1.6
        wall[ pylab.where( (1.5<=z) & (z<2.0) ) ] = 1.5
        wall[ pylab.where( (2.0<=z) ) ] = 1.4
        wall += 10

        slope[ pylab.where( z<0.5 ) ] = uvj_slope
        slope[ pylab.where( 0.5<=z ) ] = uvj_slope

        intercept[ pylab.where( z<0.5 ) ] = uvj_intercept + 0.1
        intercept[ pylab.where( 0.5<=z ) ] = uvj_intercept

        outer = pylab.zeros(len(z))
        outer[ pylab.where( (uv<slope*vj+intercept) | (uv<floor) | (vj>wall) )[0] ] = 1
        return outer.astype(int)




    cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

    zgrid = pylab.arange(0.4, 1.4, 0.00001)
    vgrid = ((1+zgrid)**2 - 1) / ((1+zgrid)**2 + 1) * constants.c.to(units.km / units.s).value
    def vrecess(zi):
    	###  return the cosmo recessional velocity at input z
    	return pylab.interp(zi, zgrid, vgrid)
    def zrecess(vi):
    	###  return the cosmo redshift at input v [km/s]
    	return pylab.interp(vi, vgrid, zgrid)

    def PolyArea(x, y):
    	#  calculate the area of an arbitrary ploygon with given vertices
        return 0.5 * abs(pylab.dot(x, pylab.roll(y, 1)) - pylab.dot(y, pylab.roll(x, 1)))



    def checkPoint(ref_ra, ref_dec, check_ra, check_dec):
    	'''
    	This function test the ra, dec from one catalogs in the convex hull from the reference catalog
    	ConvexHull reference: 
    	http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html#scipy.spatial.ConvexHull
    	example: http://stackoverflow.com/questions/31404658/check-if-points-lies-inside-a-convex-hull

    	input:
    	  ref_ra, ref_dec ---- ra and dec from reference catalog
    	  check_ra, check_dec ---- ra dec from to be checked catalog

    	output:
    	  return a array:
    	  if in, return 1
    	  if not, return 0
    	'''
    	points = pylab.zeros((len(ref_ra), 2))
    	points[:,0] = ref_ra
    	points[:,1] = ref_dec
    	hull = ConvexHull(points)
    	hull_path = Path(points[hull.vertices])

    	check = pylab.zeros(len(check_ra)) - 1
    	for i in range(len(check_ra)):
    		check[i] = hull_path.contains_point((check_ra[i], check_dec[i]))
    	return check

    def alpha_shape(points,alpha):
        """
        Compute the alpha shape (concave hull) of a set
        of points.
     
        @param points: Iterable container of points.
        @param alpha: alpha value to influence the
            gooeyness of the border. Smaller numbers
            don't fall inward as much as larger numbers.
            Too large, and you lose everything!
        """
        
        if len(points) < 4:
            # When you have a triangle, there is no sense
            # in computing an alpha shape.
            return geometry.MultiPoint(list(points)).convex_hull
     
        def add_edge(edges, edge_points, coords, i, j):
            """
            Add a line between the i-th and j-th points,
            if not in the list already
            """
            if (i, j) in edges or (j, i) in edges:
                # already added
                return
            edges.add( (i, j) )
            edge_points.append(coords[ [i, j] ])
     
        #coords = pylab.array([point.coords[0] for point in points])
        coords = points
     
        tri = Delaunay(coords)
        edges = set()
        edge_points = []
        # loop over triangles:
        # ia, ib, ic = indices of corner points of the
        # triangle
        for ia, ib, ic in tri.vertices:
            pa = coords[ia]
            pb = coords[ib]
            pc = coords[ic]
     
            # Lengths of sides of triangle
            a = pylab.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
            b = pylab.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
            c = pylab.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
     
            # Semiperimeter of triangle
            s = (a + b + c)/2.0
     
            # Area of triangle by Heron's formula
            area = pylab.sqrt(s*(s-a)*(s-b)*(s-c))
            circum_r = a*b*c/(4.0*area)
     
            # Here's the radius filter.
            #print circum_r
            if circum_r < 1.0/alpha:
                add_edge(edges, edge_points, coords, ia, ib)
                add_edge(edges, edge_points, coords, ib, ic)
                add_edge(edges, edge_points, coords, ic, ia)
     
        m = geometry.MultiLineString(edge_points)
        triangles = list(polygonize(m))
        return cascaded_union(triangles), edge_points

    def CHullRandomPoint(poly,npnts):
        ra_LB,dec_LB,ra_UB,dec_UB=poly.bounds
        rand_arr = pylab.random.random(2*npnts)
        ra_rand = (ra_UB-ra_LB)*rand_arr[:npnts] + ra_LB
        dec_rand = (dec_UB-dec_LB)*rand_arr[npnts:] + dec_LB
        testinpoly=pylab.zeros(npnts,dtype='bool')
        for iCR in range(0,npnts):
            testinpoly[iCR]=poly.contains(geometry.Point(ra_rand[iCR],dec_rand[iCR]))
        while pylab.sum(testinpoly)<npnts:
            repinds=True-testinpoly
            rand_arr = pylab.random.random(2*(npnts-pylab.sum(testinpoly)))
            ra_rand[repinds] = (ra_UB-ra_LB)*rand_arr[:npnts-pylab.sum(testinpoly)] + ra_LB
            dec_rand[repinds] = (dec_UB-dec_LB)*rand_arr[npnts-pylab.sum(testinpoly):] + dec_LB
            tmptest=testinpoly[repinds]
            for iCR in range(0,npnts-pylab.sum(testinpoly)):
                tmptest[iCR]=poly.contains(geometry.Point(ra_rand[repinds][iCR],dec_rand[repinds][iCR]))
            testinpoly[repinds]=tmptest
        return ra_rand,dec_rand

    def CheckPoints(poly,ras,decs):
        g=[poly.contains(geometry.Point(ras[ip],decs[ip])) for ip in range(0,len(ras))]
        return pylab.arange(len(ras))[pylab.array(g)==True]






    class voronoi_catalog:
    	def __init__(self, version, ngals, nzbins):
    		self.version = version            # photometric catalog
    		self.zlos = pylab.zeros(nzbins)
    		self.zhis = pylab.zeros(nzbins)
    		self.zbars = pylab.zeros(nzbins)
    		self.ntot = pylab.zeros(nzbins)   # Median number of galaxies in bin from all iterat
    		self.zspec_fraction = pylab.zeros(nzbins)   # Median fraction of galaxies in bin with spec-zs
    		self.dens_matrix = pylab.zeros((ngals, nzbins)) - 99
    		self.overdens_matrix = pylab.zeros((ngals, nzbins)) - 99
    		self.bunit_dens = ''
    		self.bunit_overdens = ''

    class substructure:
    	def __init__(self, name, ra, dec, zmean, vsigma05mpc, err_vsigma05mpc, Ngal05mpc, vsigma1mpc, err_vsigma1mpc, Ngal1mpc, lmvir, err_lmvir):
    		'''
    		All information comes from ORELSE.dispersions
    		name                    name of substructure
    		ra                      luminosity weighted central RA
    		dec                     luminosity weighted central Dec
    		zmean                   mean redshift of substructure
    		vsigma05mpc      km/s   velocity dispersion for galaxies at <0.5 Mpc
    		err_vsigma05mpc  km/s   error on previous
    		Ngal05mpc               number of zspec galaxies at <0.5 Mpc
    		vsigma1mpc       km/s   velocity dispersion for galaxies at <1.0 Mpc
    		err_vsigma1mpc   km/s   error on previous
    		Ngal1mpc                number of zspec galaxies at <1.0 Mpc
    		lmvir            lmsol  log of M_vir (see Lemaux+2012)
    		err_lmvir        lmsol  error on previous
    		'''
    		self.name = name
    		self.ra = ra
    		self.dec = dec
    		self.zmean = zmean
    		self.vsigma05mpc = vsigma05mpc
    		self.vsigma1mpc = vsigma1mpc
    		self.err_vsigma05mpc = err_vsigma05mpc
    		self.err_vsigma1mpc = err_vsigma1mpc
    		self.Ngal05mpc = Ngal05mpc
    		self.Ngal1mpc = Ngal1mpc
    		self.lmvir = lmvir
    		self.err_lmvir = err_lmvir

    		self.vmean = vrecess(zmean)
    		self.vlo_3sig = self.vmean - 3*vsigma1mpc
    		self.vhi_3sig = self.vmean + 3*vsigma1mpc
    		self.zlo_3sig = zrecess(self.vlo_3sig)
    		self.zhi_3sig = zrecess(self.vhi_3sig)

    github_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
    class field:
    	def __init__(self, name, version, zclust, sigmaz, alpha=50, chi2red_thresh=10, uvj_slope=0.88, uvj_intercept=0.59):
    		self.name = name          # e.g. "NEP 200"
    		self.version = version    # e.g. "nep200_v0.0.4"
    		self.zclust = zclust      # cluster redshift
    		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
    		self.zspec_lo = 0         # lower redshift bound for specz
    		self.zspec_hi = 0         # upper redshift bound for specz
    		self.substructures = []
    		self.alpha = alpha

    		self.cat = gunzip_read_gzip('%s/%s/%s.cat.gz' % (github_dir, version, version), readcat=1)
    		self.zout = gunzip_read_gzip('%s/%s/%s.zout.gz' % (github_dir, version, version), readzout=1)
    		self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (github_dir, version, version), readcat=1)
    		self.restframe = gunzip_read_gzip('%s/%s/%s.restframe.gz' % (github_dir, version, version), readcat=1)
    		#self.rfcolors = gunzip_read_gzip('%s/%s/%s.restframe_colors.gz' % (github_dir, version, version), readcat=1)

    		self.masslimits = mypy.readcat('../data/completeness/masslimits_%s_master.dat' % name)


    		###  SETTING OBJECTS IDENTIFIED AS SECURE STARS FROM SPECTROSCOPY TO use=0
    		self.crossmatch = gunzip_read_gzip('%s/%s/%s.crossmatch.gz' % (github_dir, version, version), readcat=1, dtype=str)
    		self.star_inds = pylab.find(self.crossmatch.Q == '-1')
    		for i_star in self.star_inds:
    			id_phot_arr = self.crossmatch.id_phot[i_star].split(',')
    			for id_phot in id_phot_arr:
    				if id_phot != '-1':
    					self.cat.use[int(id_phot)-1] *= 0


    		print '  reading: %s_voronoi.pickle' % name
    		#self.voronoi = pickle.load(open('../data/%s_voronoi.pickle' % name, 'rb'))
    		self.voronoi = pickle.load(open('/Users/atomczak/Dropbox/ORELSE/DATA/VORONOI/old_prior_2_bugfix/%s_voronoi.pickle' % name, 'rb'))
    		#self.overdens_max = []
    		#for vi in range(len(self.voronoi.overdens_matrix)):
    		#	self.overdens_max.append(self.voronoi.overdens_matrix[vi].max())
    		#self.overdens_max = pylab.array(self.overdens_max)




    		###  assigning Voroini overdensity values to galaxies
    		self.overdensities = pylab.zeros(len(self.cat.id))

    		for i_gal in range(len(self.cat.id)):
    			overdensities_gal = []
    			for i_vmap in range(len(self.voronoi.zlos)):
    				zlo, zhi = self.voronoi.zlos[i_vmap], self.voronoi.zhis[i_vmap]
    				if zlo < self.fout.z[i_gal] < zhi:
    					overdensities_gal.append(self.voronoi.overdens_matrix[i_gal][i_vmap])

    				if len(overdensities_gal) > 0:
    					#self.overdensities[i_gal] = pylab.mean(overdensities_gal)
    					self.overdensities[i_gal] = max(overdensities_gal)
    				else:
    					self.overdensities[i_gal] = pylab.nan





    		###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
    		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
    		cinds = pylab.find((chi2red > chi2red_thresh) & (self.cat.z_spec < 0))
    		self.cat.use[cinds] = 0.


    		###  UVJ classification
    		self.UV_color = -2.5 * pylab.log10(self.restframe.restflux_U / self.restframe.restflux_V)
    		self.VJ_color = -2.5 * pylab.log10(self.restframe.restflux_V / self.restframe.restflux_J)
    		#self.uvj_class = mypy.uvj_select(self.UV_color, self.VJ_color, self.fout.z)
    		self.uvj_class = uvj_select(self.UV_color, self.VJ_color, self.fout.z, uvj_slope=uvj_slope, uvj_intercept=uvj_intercept)

    		###  UVJ classification based on EAZY in 2-filter mode
    		#self.UV_color_2filter = self.rfcolors.U_V_color
    		#self.VJ_color_2filter = self.rfcolors.V_J_color
    		#self.uvj_class_2filter = mypy.uvj_select(self.UV_color_2filter, self.VJ_color_2filter, self.fout.z)

    		###  NUVrJ classification
    		#self.NUVr_color = -2.5 * pylab.log10(self.restframe.restflux_NUV / self.restframe.restflux_r)
    		#self.rJ_color = -2.5 * pylab.log10(self.restframe.restflux_r / self.restframe.restflux_J)
    		#self.nuvrj_class = mypy.nuvrj_select(self.NUVr_color, self.rJ_color, self.fout.z)



    		xyrd1 = self.cat.x[0], self.cat.y[0], self.cat.ra[0], self.cat.dec[0]
    		xyrd2 = self.cat.x[1], self.cat.y[1], self.cat.ra[1], self.cat.dec[1]
    		d_arcsec = mypy.radec_sep(xyrd1[2], xyrd1[3], xyrd2[2], xyrd2[3])
    		d_pixel = ((xyrd1[0]-xyrd2[0])**2 + (xyrd1[1] - xyrd2[1])**2)**0.5 
    		self.px_scale = d_arcsec / d_pixel

    		### getting z band magnitude
    		try: self.zmagnitude = 25 - 2.5 * pylab.log10(self.cat.fluxauto_z)
    		except: pass
    		if name != 'SC1324':
    			try: self.zmagnitude = 25 - 2.5 * pylab.log10(self.cat.fluxauto_Zplus)
    			except: pass

    		###  setting spatial flag based on Convex Hull method
    		#zspecinds = pylab.find(self.cat.z_spec > 0)
    		#self.spatial_flag = checkPoint(self.cat.ra[zspecinds], self.cat.dec[zspecinds], self.cat.ra, self.cat.dec)
    		#self.inds_spatial = pylab.find(self.spatial_flag == 1)

    		###  setting spatial flag based on Concave Hull method (creates a "tighter" spatial selection than Cnvex Hull)
    		zspecinds = pylab.find(self.cat.z_spec > 0)
    		points = pylab.array(zip(self.cat.x[zspecinds], self.cat.y[zspecinds]))
    		self.concave_hull, self.edge_points = alpha_shape(points, alpha)
    		self.inds_spatial = CheckPoints(self.concave_hull.buffer(10), self.cat.x, self.cat.y)
    		self.area_pix2 = self.concave_hull.buffer(10).area
    		self.area_arcmin2 = self.area_pix2 * (self.px_scale/60.)**2
    		print ''


    fields = []
    fields.append(field('N200',    'nep200_v0.0.7',       0.691,  0.027, alpha=1./600, chi2red_thresh=7, uvj_intercept=0.49))
    fields.append(field('SC1324',  'sc1324_v0.0.4',       0.755,  0.033, alpha=1./600, chi2red_thresh=10))
    fields.append(field('RCS0224', 'rcs0224-0002_v0.0.4', 0.772,  0.027, alpha=1./500, chi2red_thresh=10))
    fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.8', 0.813,  0.021, alpha=1./500, chi2red_thresh=8))
    fields.append(field('N5281',   'nep5281_v0.0.4',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10))
    #fields.append(field('SG0023',  'sg0023+0423_v0.1.9',  0.845,  0.025, alpha=, chi2red_thresh=14))
    fields.append(field('SC1604',  'sc1604_v0.0.6',       0.910,  0.029, alpha=1./500, chi2red_thresh=10))
    #fields.append(field('CL1429',  'cl1429+4241_v0.0.1',  0.920,  0.051, alpha=1./500, chi2red_thresh=10))
    fields.append(field('SC0910',  'cl0910+5422_v0.0.5',  1.110,  0.035, alpha=1./500, chi2red_thresh=10))
    fields.append(field('SC0849',  'sc0849+4452_v0.0.4',  1.261,  0.029, alpha=1./600, chi2red_thresh=10, uvj_intercept=0.54))
    print ''



    for i in [12203, 13820]: fields[0].cat.use[i-1] *= 0
    for i in [4386, 10388, 10907]: fields[1].cat.use[i-1] *= 0
    for i in [17351, 26627, 38501, 43517, 51926, 55037, 55089, 61118]: fields[2].cat.use[i-1] *= 0
    for i in []: fields[3].cat.use[i-1] *= 0
    for i in [11331, 17375]: fields[4].cat.use[i-1] *= 0
    for i in [105231, 170586]: fields[5].cat.use[i-1] *= 0
    for i in [51114, 57254]: fields[6].cat.use[i-1] *= 0
    for i in [9419, 15619, 18725, 22467, 25342, 28122]: fields[7].cat.use[i-1] *= 0



    ii = 0
    fields[ii].substructures.append(substructure('RXJ1757', 269.33196, 66.525991, 0.6931, 541.1, 139.2, 17, 862.3, 107.9, 34, 14.832, 0.250))

    ii += 1
    fields[ii].substructures.append(substructure('Cluster-A*', 201.20097, 30.1920, 0.7556, 1019.2, 142.0, 25, 873.4, 110.8, 43, 14.833, 0.254))
    fields[ii].substructures.append(substructure('Cluster-B', 201.08830, 30.2156, 0.6979, 897.0, 394.7, 11, 677.1, 143.6, 13, 14.516, 0.424))
    fields[ii].substructures.append(substructure('Group-C', 201.00770, 30.4164, 0.7382, 143.8, 40.7, 6, 205.9, 90.1, 8, 12.955, 0.875))
    fields[ii].substructures.append(substructure('Cluster-I', 201.20096, 30.9680, 0.6957, 646.1, 113.0, 16, 891.5, 128.9, 27, 14.875, 0.291))

    ###  NOTE: The 3rd substructure is a serendipitous
    ###        group/cluster behind the main LSS
    ii += 1
    fields[ii].substructures.append(substructure('RCS0224-A*+', 36.15714, -0.0949, 0.7780, 713.3, 179.3, 14, 825.4, 193.2, 34, 14.754, 0.468))
    fields[ii].substructures.append(substructure('RCS0224-B', 36.14123, -0.0394, 0.7781, 815.5, 138.5, 24, 710.7, 58.8, 52, 14.559, 0.165))
    #fields[ii].substructures.append(substructure('RCS0224-S**', 36.32021, -0.0928, 0.8454, 204.2, 107.7, 7, 437.7, 115.8, 15, 13.911, 0.529))

    ii += 1
    fields[ii].substructures.append(substructure('RXJ1716-A*', 259.20162, 67.1392, 0.8116, 1150.1, 161.8, 25, 1150.2, 113.4, 62, 15.178, 0.197))
    fields[ii].substructures.append(substructure('RXJ1716-B+', 259.25229, 67.1533, 0.8127, 767.1, 112.0, 20, 693.5, 61.0, 40, 14.519, 0.176))

    ii += 1
    fields[ii].substructures.append(substructure('RXJ1821', 275.38451, 68.465768, 0.8168, 1146.1, 124.8, 27, 1119.6, 99.6, 52, 15.142, 0.178))

    ii += 1
    fields[ii].substructures.append(substructure('Cluster-A', 241.09311, 43.0821, 0.8984, 576.9, 120.8, 23, 722.4, 134.5, 35, 14.551, 0.372))
    fields[ii].substructures.append(substructure('Cluster-B', 241.10796, 43.2397, 0.8648, 792.9, 80.1, 28, 818.4, 74.2, 49, 14.722, 0.269))
    fields[ii].substructures.append(substructure('Group-C*', 241.03142, 43.2679, 0.9344, 1079.6, 289.6, 12, 453.5, 39.6, 32, 13.935, 0.039))
    fields[ii].substructures.append(substructure('Cluster-D+', 241.14094, 43.3539, 0.9227, 675.6, 180.1, 40, 688.2, 88.1, 70, 14.481, 0.256))
    fields[ii].substructures.append(substructure('Group-F++', 241.20104, 43.3684, 0.9331, 619.9, 135.0, 14, 541.9, 110.0, 20, 14.168, 0.406))
    fields[ii].substructures.append(substructure('Group-G', 240.92745, 43.4030, 0.9019, 398.5, 85.4, 7, 539.3, 124.0, 18, 14.169, 0.460))
    fields[ii].substructures.append(substructure('Group-H', 240.89890, 43.3669, 0.8528, 283.2, 71.7, 9, 287.0, 68.3, 10, 13.359, 0.074))
    fields[ii].substructures.append(substructure('Group-I', 240.79746, 43.3915, 0.9024, 163.0, 65.1, 5, 333.0, 129.4, 7, 13.541, 0.777))

    ###  CL 1429
    #ii += 1


    ii += 1
    fields[ii].substructures.append(substructure('0910-A', 137.51190, 54.3103, 1.1024, 506.8, 235.9, 10, 893.0, 285.8, 18, 14.777, 0.640))
    fields[ii].substructures.append(substructure('0910-B', 137.68463, 54.3736, 1.1016, 906.3, 192.4, 11, 795.5, 138.1, 22, 14.627, 0.347))

    ii += 1
    fields[ii].substructures.append(substructure('SC0849A+', 132.23463, 44.761780, 1.2622, 609.1,  178.1, 8,  708.0, 186.9, 13, 14.437, 0.344))
    fields[ii].substructures.append(substructure('SC0849B+', 132.29977, 44.865903, 1.2636, 433.3,  389.7, 6,  286.3, 64.9,  9,  13.257, 0.295))
    fields[ii].substructures.append(substructure('SC0849C+', 132.24443, 44.866012, 1.2631, 1006.8, 310.8, 17, 766.2, 97.6,  25, 14.540, 0.166))	
    fields[ii].substructures.append(substructure('SC0849D',  132.15110, 44.899400, 1.2700, 930.9,  147.7, 15, 788.7, 136.1, 24, 14.576, 0.225))
    fields[ii].substructures.append(substructure('SC0849E+', 132.27496, 44.959253, 1.2600, 366.0,  205.9, 10, 414.1, 107.8, 14, 13.739, 0.339))

    ###  Note: I'm removing sg0023 because we feel that it isn't belong in this 
    ###        study due, in large part, to the fact that it's a lower-mass system.
    #fields[].substructures.append(substructure('0023-A*', 6.02560, 4.3590, 0.8396, 507.0, 125.5, 10, 412.8, 119.2, 14, 13.836, 0.578))
    #fields[].substructures.append(substructure('0023-B1**', 5.97570, 4.3884, 0.8290, 106.2, 51.4, 12, 176.3, 29.6, 23, 12.730, 0.336))
    #fields[].substructures.append(substructure('0023-B2+', 5.96970, 4.3820, 0.8453, 231.3, 53.8, 18, 277.8, 41.0, 38, 13.319, 0.295))
    #fields[].substructures.append(substructure('0023-C', 5.92470, 4.3807, 0.8466, 543.8, 58.7, 22, 385.3, 54.3, 45, 13.744, 0.282))
    #fields[].substructures.append(substructure('0023-M++', 5.96740, 4.3199, 0.8472, 487.3, 84.8, 7, 418.8, 68.9, 14, 13.853, 0.330))



    print '\nreading expanded parameter space FAST catalogs...'
    fields[0].fout = mypy.readcat('%s/nep200_v0.0.7/nep200_v0.0.7_ZFparams.fout.gz' % github_dir)
    fields[1].fout = mypy.readcat('%s/sc1324_v0.0.4/sc1324_v0.0.4_ZFparams.fout.gz' % github_dir)
    fields[2].fout = mypy.readcat('%s/rcs0224-0002_v0.0.4/rcs0224-0002_v0.0.4_ZFparams.fout.gz' % github_dir)
    fields[3].fout = mypy.readcat('%s/rxj1716+6708_v0.0.8/rxj1716+6708_v0.0.8_ZFparams.fout.gz' % github_dir)
    fields[4].fout = mypy.readcat('%s/nep5281_v0.0.4/nep5281_v0.0.4_ZFparams.fout.gz' % github_dir)
    fields[5].fout = mypy.readcat('%s/sc1604_v0.0.6/sc1604_v0.0.6_ZFparams.fout.gz' % github_dir)
    #fields[6].fout = mypy.readcat('%s/cl1429+4241_v0.0.1/cl1429+4241_v0.0.1_ZFparams.fout.gz' % github_dir)
    fields[6].fout = mypy.readcat('%s/cl0910+5422_v0.0.5/cl0910+5422_v0.0.5_ZFparams.fout.gz' % github_dir)
    fields[7].fout = mypy.readcat('%s/sc0849+4452_v0.0.4/sc0849+4452_v0.0.4_ZFparams.fout.gz' % github_dir)
    print 'done!\n'






















fig = pylab.figure(figsize=(14, 6.8))
sp = fig.add_subplot(111, aspect='auto')

sp.minorticks_on()
sp.set_xlabel('redshift')
sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp.set_ylabel('log( 1 + $\delta_{\mathrm{gal}}$ )')





zlims = [0.55, 1.3]

dm = 0.25
dv = 0.5
vlims = pylab.arange(-0.5, 2.+dv, dv)
mlims = pylab.arange(8.75-dm/2., 11.75+dm/2., dm)

#for zi in zlims:
#	sp.plot([zi, zi], [min(vlims), max(vlims)], color='r', lw=2, zorder=2)


for vi in vlims:
	sp.plot([min(mlims), max(mlims)], [vi, vi], color='w', lw=4, zorder=2)
	sp.plot([min(mlims), max(mlims)], [vi, vi], color='r', lw=2, zorder=3)

for mi in mlims:
	sp.plot([mi, mi], [min(vlims), max(vlims)], color='w', lw=4, zorder=2)
	sp.plot([mi, mi], [min(vlims), max(vlims)], color='r', lw=2, zorder=3)


sp.axis([min(mlims)-dm*0.75, max(mlims)+dm*0.75, min(vlims)-dv*0.75, max(vlims)+dv*0.75])



###  version 1:  small black points

zspec_master = []
z_master = []
v_master = []
m_master = []

for i_f in range(len(fields)):

    f = fields[i_f]

    inds = pylab.find((f.fout.z[f.inds_spatial] > min(zlims)) & 
                      (f.fout.z[f.inds_spatial] < max(zlims)) &
                     -(pylab.isnan(f.fout.lmass[f.inds_spatial])) &
                      (f.overdensities[f.inds_spatial] < 2.2) &
                      (f.fout.lmass[f.inds_spatial] < 11.75) &
                      (f.fout.lmass[f.inds_spatial] > pylab.interp(f.fout.z[f.inds_spatial], f.masslimits.z, f.masslimits.masslim_empirical)) &
                      (f.cat.use[f.inds_spatial] == 1))

    z_master += f.fout.z[f.inds_spatial][inds].tolist()
    v_master += f.overdensities[f.inds_spatial][inds].tolist()
    m_master += f.fout.lmass[f.inds_spatial][inds].tolist()


    inds_spec = pylab.find((f.fout.z[f.inds_spatial] > min(zlims)) & 
                      (f.fout.z[f.inds_spatial] < max(zlims)) &
                      (f.cat.z_spec[f.inds_spatial] > 0) &
                     -(pylab.isnan(f.fout.lmass[f.inds_spatial])) &
                      (f.overdensities[f.inds_spatial] < 2.2) &
                      (f.fout.lmass[f.inds_spatial] < 11.75) &
                      (f.fout.lmass[f.inds_spatial] > pylab.interp(f.fout.z[f.inds_spatial], f.masslimits.z, f.masslimits.masslim_empirical)) &
                      (f.cat.use[f.inds_spatial] == 1))

    zspec_master += f.cat.z_spec[f.inds_spatial][inds_spec].tolist()




	#sp.plot(f.fout.lmass[f.inds_spatial][inds], f.overdensities[f.inds_spatial][inds], 
	#	    ls='', marker='o', ms=1, zorder=1)


zspec_master = pylab.array(zspec_master)
z_master = pylab.array(z_master)
v_master = pylab.array(v_master)
m_master = pylab.array(m_master)




cmap = pylab.cm.gray_r
nbins = 5
hist2d, xedges, yedges = pylab.histogram2d(m_master, v_master,
	                                       bins=(len(mlims)*nbins+2*nbins, len(vlims)*nbins+2*nbins), 
	                                       range=([min(mlims)-dm, max(mlims)+dm], [min(vlims)-dv, max(vlims)+dv]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

nlo, nhi = hist2d.min(), hist2d.max()
hist2d[pylab.where(hist2d == 0)] = pylab.nan
hist2d[pylab.where(hist2d > 0)] = pylab.log10(hist2d[pylab.where(hist2d > 0)])
nhi = hist2d[pylab.where(-pylab.isnan(hist2d))].max()

#colormap = sp.imshow(hist2d.T[::-1,:], extent=extent, vmin=-0.3, vmax=nhi+0.3, interpolation='nearest', cmap=cmap)
colormap = sp.imshow(hist2d.T, extent=extent, vmin=-0.3, vmax=nhi+0.3, interpolation='nearest', cmap=cmap)


"""
font = FontProperties()
font.set_family('sans-serif')

for vlo, vhi in zip(vlims[:-1], vlims[1:]):
	for mlo, mhi in zip(mlims[:-1], mlims[1:]):
		n = len(pylab.find((m_master > mlo) & (m_master < mhi) & 
			               (v_master > vlo) & (v_master < vhi)))

		t = sp.text(mlo+dm*0.05, vhi-dv*0.25, '%i' % n, fontsize=15, color='k', fontweight='bold', fontproperties=font,
	    	        path_effects=[PathEffects.withStroke(linewidth=4, foreground='w')])
"""




'''
for ix in range(len(hist2d)):
	for iy in range(len(hist2d[ix])):
		if -pylab.isnan(hist2d[ix][iy]) & (hist2d[ix][iy] < pylab.log10(3)):
			sinds = pylab.find((m_master > xedges[ix]) & (m_master < xedges[ix+1]) & 
				               (v_master > yedges[iy]) & (v_master < yedges[iy+1]))
			sp.plot(m_master[sinds], v_master[sinds], 'ko', ms=1)
'''



sp.set_aspect('auto')








