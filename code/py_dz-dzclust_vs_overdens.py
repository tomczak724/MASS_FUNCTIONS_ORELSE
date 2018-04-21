
import os
import mypy
import time
import pylab
import pickle
import subprocess
from matplotlib.path import Path
import shapely.geometry as geometry
from scipy.spatial import ConvexHull, Delaunay
from shapely.ops import cascaded_union, polygonize
import matplotlib.patheffects as PathEffects
from astropy import units, constants, cosmology
from matplotlib.backends.backend_pdf import PdfPages

t0 = time.time()
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


data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'





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
		name                   name of substructure
		ra                     luminosity weighted central RA
		dec                    luminosity weighted central Dec
		zmean                  mean redshift of substructure
		vsigma05mpc      km/s  velocity dispersion for galaxies at <0.5 Mpc
		err_vsigma05mpc  km/s  error on previous
		Ngal05mpc              number of zspec galaxies at <0.5 Mpc
		vsigma1mpc       km/s  velocity dispersion for galaxies at <1.0 Mpc
		err_vsigma1mpc   km/s  error on previous
		Ngal1mpc               number of zspec galaxies at <1.0 Mpc
		lmvir           lmsol  log of M_vir (see Lemaux+2012)
		err_lmvir       lmsol  error on previous
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

class field:
	def __init__(self, name, version, zclust, sigmaz, alpha=50):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz
		self.substructures = []

		self.cat = gunzip_read_gzip('%s/%s/%s.cat.gz' % (data_dir, version, version), readcat=1)
		self.zout = gunzip_read_gzip('%s/%s/%s.zout.gz' % (data_dir, version, version), readzout=1)
		self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (data_dir, version, version), readcat=1)
		self.restframe = gunzip_read_gzip('%s/%s/%s.restframe.gz' % (data_dir, version, version), readcat=1)

		print '  reading: %s_voronoi.pickle' % name
		self.voronoi = pickle.load(open('../data/%s_voronoi.pickle' % name, 'rb'))


		###  UVJ classification
		self.UV_color = -2.5 * pylab.log10(self.restframe.restflux_U / self.restframe.restflux_V)
		self.VJ_color = -2.5 * pylab.log10(self.restframe.restflux_V / self.restframe.restflux_J)
		self.uvj_class = mypy.uvj_select(self.UV_color, self.VJ_color, self.fout.z)


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
fields.append(field('N200',    'nep200_v0.0.5',       0.691,  0.027, alpha=1./600))
fields.append(field('SC1324',  'sc1324_v0.0.2',       0.755,  0.033, alpha=1./600))
fields.append(field('RCS0224', 'rcs0224-0002_v0.0.2', 0.772,  0.027, alpha=1./500))
fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.7', 0.813,  0.021, alpha=1./500))
fields.append(field('N5281',   'nep5281_v0.0.2',      0.818,  0.029, alpha=1./1000))
#fields.append(field('SG0023',  'sg0023+0423_v0.1.9',  0.845,  0.025, alpha=))
fields.append(field('SC1604',  'sc1604_v0.0.3',       0.910,  0.029, alpha=1./500))
fields.append(field('SC0910',  'cl0910+5422_v0.0.3',  1.110,  0.035, alpha=1./500))
fields.append(field('SC0849',  'sc0849+4452_v0.0.2',  1.261,  0.029, alpha=1./600))
print ''


fields[0].substructures.append(substructure('RXJ1757', 269.33196, 66.525991, 0.6931, 541.1, 139.2, 17, 862.3, 107.9, 34, 14.832, 0.250))

fields[1].substructures.append(substructure('Cluster-A*', 201.20097, 30.1920, 0.7556, 1019.2, 142.0, 25, 873.4, 110.8, 43, 14.833, 0.254))
fields[1].substructures.append(substructure('Cluster-B', 201.08830, 30.2156, 0.6979, 897.0, 394.7, 11, 677.1, 143.6, 13, 14.516, 0.424))
fields[1].substructures.append(substructure('Group-C', 201.00770, 30.4164, 0.7382, 143.8, 40.7, 6, 205.9, 90.1, 8, 12.955, 0.875))
fields[1].substructures.append(substructure('Cluster-I', 201.20096, 30.9680, 0.6957, 646.1, 113.0, 16, 891.5, 128.9, 27, 14.875, 0.291))

fields[2].substructures.append(substructure('RCS0224-A*+', 36.15714, -0.0949, 0.7780, 713.3, 179.3, 14, 825.4, 193.2, 34, 14.754, 0.468))
fields[2].substructures.append(substructure('RCS0224-B', 36.14123, -0.0394, 0.7781, 815.5, 138.5, 24, 710.7, 58.8, 52, 14.559, 0.165))
fields[2].substructures.append(substructure('RCS0224-S**', 36.32021, -0.0928, 0.8454, 204.2, 107.7, 7, 437.7, 115.8, 15, 13.911, 0.529))

fields[3].substructures.append(substructure('RXJ1716-A*', 259.20162, 67.1392, 0.8116, 1150.1, 161.8, 25, 1150.2, 113.4, 62, 15.178, 0.197))
fields[3].substructures.append(substructure('RXJ1716-B+', 259.25229, 67.1533, 0.8127, 767.1, 112.0, 20, 693.5, 61.0, 40, 14.519, 0.176))

fields[4].substructures.append(substructure('RXJ1821', 275.38451, 68.465768, 0.8168, 1146.1, 124.8, 27, 1119.6, 99.6, 52, 15.142, 0.178))

###  Note: substructures B, C, G, and H are excluded because they
###        are nearby in projection but are significantly separated
###        along the LoS. Included them would contaminate the selection. 
fields[5].substructures.append(substructure('Cluster-A', 241.09311, 43.0821, 0.8984, 576.9, 120.8, 23, 722.4, 134.5, 35, 14.551, 0.372))
fields[5].substructures.append(substructure('Cluster-B', 241.10796, 43.2397, 0.8648, 792.9, 80.1, 28, 818.4, 74.2, 49, 14.722, 0.269))
fields[5].substructures.append(substructure('Group-C*', 241.03142, 43.2679, 0.9344, 1079.6, 289.6, 12, 453.5, 39.6, 32, 13.935, 0.039))
fields[5].substructures.append(substructure('Cluster-D+', 241.14094, 43.3539, 0.9227, 675.6, 180.1, 40, 688.2, 88.1, 70, 14.481, 0.256))
fields[5].substructures.append(substructure('Group-F++', 241.20104, 43.3684, 0.9331, 619.9, 135.0, 14, 541.9, 110.0, 20, 14.168, 0.406))
fields[5].substructures.append(substructure('Group-G', 240.92745, 43.4030, 0.9019, 398.5, 85.4, 7, 539.3, 124.0, 18, 14.169, 0.460))
fields[5].substructures.append(substructure('Group-H', 240.89890, 43.3669, 0.8528, 283.2, 71.7, 9, 287.0, 68.3, 10, 13.359, 0.074))
fields[5].substructures.append(substructure('Group-I', 240.79746, 43.3915, 0.9024, 163.0, 65.1, 5, 333.0, 129.4, 7, 13.541, 0.777))

fields[6].substructures.append(substructure('0910-A', 137.51190, 54.3103, 1.1024, 506.8, 235.9, 10, 893.0, 285.8, 18, 14.777, 0.640))
fields[6].substructures.append(substructure('0910-B', 137.68463, 54.3736, 1.1016, 906.3, 192.4, 11, 795.5, 138.1, 22, 14.627, 0.347))

fields[7].substructures.append(substructure('SC0849A+', 132.23463, 44.761780, 1.2622, 609.1,  178.1, 8,  708.0, 186.9, 13, 14.437, 0.344))
fields[7].substructures.append(substructure('SC0849B+', 132.29977, 44.865903, 1.2636, 433.3,  389.7, 6,  286.3, 64.9,  9,  13.257, 0.295))
fields[7].substructures.append(substructure('SC0849C+', 132.24443, 44.866012, 1.2631, 1006.8, 310.8, 17, 766.2, 97.6,  25, 14.540, 0.166))	
fields[7].substructures.append(substructure('SC0849D',  132.15110, 44.899400, 1.2700, 930.9,  147.7, 15, 788.7, 136.1, 24, 14.576, 0.225))
fields[7].substructures.append(substructure('SC0849E+', 132.27496, 44.959253, 1.2600, 366.0,  205.9, 10, 414.1, 107.8, 14, 13.739, 0.339))


###  Note: I'm removing sg0023 because we feel that it isn't belong in this 
###        study due, in large part, to the fact that it's a lower-mass system.
#fields[].substructures.append(substructure('0023-A*', 6.02560, 4.3590, 0.8396, 507.0, 125.5, 10, 412.8, 119.2, 14, 13.836, 0.578))
#fields[].substructures.append(substructure('0023-B1**', 5.97570, 4.3884, 0.8290, 106.2, 51.4, 12, 176.3, 29.6, 23, 12.730, 0.336))
#fields[].substructures.append(substructure('0023-B2+', 5.96970, 4.3820, 0.8453, 231.3, 53.8, 18, 277.8, 41.0, 38, 13.319, 0.295))
#fields[].substructures.append(substructure('0023-C', 5.92470, 4.3807, 0.8466, 543.8, 58.7, 22, 385.3, 54.3, 45, 13.744, 0.282))
#fields[].substructures.append(substructure('0023-M++', 5.96740, 4.3199, 0.8472, 487.3, 84.8, 7, 418.8, 68.9, 14, 13.853, 0.330))
















###################
###  ALL GALAXIES
###################

zmag_thresh = 23.5

dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.5+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

overdens_bins = pylab.array([0., 0.5, 1., 1.5, 2.0])
overdens_bars = (overdens_bins[1:] + overdens_bins[:-1]) / 2.
voronoi_labels = ['0.0 < log(1+$\delta$) < 0.5', '0.5 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 1.5', '1.5 < log(1+$\delta$) < 2.0']




###   1.5 sigma_zphot
###   (zphot-zlss)/(1+zphot)  vs  (zspec-zlss)/(1+zspec)

n_truepos_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_trueneg_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_falsepos_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_falseneg_lmass = pylab.zeros((len(fields), len(lmassbars)))

n_truepos_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_trueneg_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_falsepos_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_falseneg_overdens = pylab.zeros((len(fields), len(overdens_bars)))

n_truepos_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_trueneg_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_falsepos_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_falseneg_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))

for i in range(len(fields)):

    f = fields[i]

    dzphot_thresh = 1.5 * f.sigmaz

    ###  checking if galaxy is in spatial selection regions
    zinds = pylab.find((f.cat.z_spec[f.inds_spatial] > 0) & 
                       (f.zmagnitude[f.inds_spatial] <= zmag_thresh) &
                       (f.fout.lmass[f.inds_spatial] > lmassbins[0]))

    zinds = f.inds_spatial[zinds]

    ra_spec = f.cat.ra[zinds]
    dec_spec = f.cat.dec[zinds]


    ###  Identifying which LSS substructure a galaxy belongs to.
    ###  This is simply defined as a galaxy needs to be <=0.5Mpc
    ###  of any substructure center. If there are two such
    ###  substructures, then the nearest is chosen.

    zspec = []
    zphot = []
    zlss = []
    lmass = []
    vor = []
    truefalse_array = []

    for i_gal in range(len(zinds)):

        ra0, dec0 = ra_spec[i_gal], dec_spec[i_gal]
        zspec_i = f.cat.z_spec[zinds][i_gal]
        zphot_i = f.zout.z_peak[zinds][i_gal]

        ###  First, a galaxy is deemed a member of the overall
        ###  LSS if it lies within 3vsigma of at least 1 substructure
        lss_member_status = False
        substructure_member_status = []
        rproj_list = []
        for ss in f.substructures:
            rproj_arcsec = mypy.radec_sep(ra0, dec0, ss.ra, ss.dec)
            asecperkpc = cosmo.arcsec_per_kpc_proper(ss.zmean)
            rproj_kpc = rproj_arcsec / asecperkpc.value
            rproj_list.append(rproj_kpc)

            if ss.zlo_3sig < zspec_i < ss.zhi_3sig:
                lss_member_status = True
                substructure_member_status.append(True)
            else:
                substructure_member_status.append(False)


        ###  if the galaxy is >2Mpc away from all substructures skip it
        #if min(rproj_list) >= 2000: continue


        ###  If it is an LSS member, identify which substructure
        ###  it is a member of. If only 1 substructure matches
        ###  then this is easy...
        rproj_list = pylab.array(rproj_list)
        substructure_member_status = pylab.array(substructure_member_status)
        ssinds = pylab.find(substructure_member_status)

        if len(ssinds) == 1:
            zlo, zhi = f.substructures[ssinds[0]].zlo_3sig, f.substructures[ssinds[0]].zhi_3sig
            zlss_i = f.substructures[ssinds[0]].zmean

        elif len(ssinds) > 1:
            sub_rproj_list = rproj_list[ssinds]
            rind = ssinds[sub_rproj_list.tolist().index(min(sub_rproj_list))]
            zlo, zhi = f.substructures[rind].zlo_3sig, f.substructures[rind].zhi_3sig
            zlss_i = f.substructures[rind].zmean

        elif len(ssinds) == 0:
            rind = rproj_list.tolist().index(min(rproj_list))
            zlo, zhi = f.substructures[rind].zlo_3sig, f.substructures[rind].zhi_3sig
            zlss_i = f.substructures[rind].zmean


        ###  Finding voronoi zslice that is "centered" on the LSS
        zdif = abs(f.voronoi.zbars - zlss_i)
        vi = pylab.find(zdif == zdif.min())[0]
        vor.append(f.voronoi.overdens_matrix[zinds][i_gal][vi])


        ### now that we've grabbed the galaxy, determine if it is:
        ###   0 = true +
        ###   1 = true -
        ###   2 = false +
        ###   3 = false -
        zspec.append(zspec_i)
        zphot.append(zphot_i)
        zlss.append(zlss_i)
        lmass.append(f.fout.lmass[zinds][i_gal])
        dzphot_i = (zphot_i - zlss_i) / (1 + zphot_i)
        if zlo < zspec_i < zhi and abs(dzphot_i) <= dzphot_thresh:
            truefalse_array.append(0)
        if (zspec_i < zlo and abs(dzphot_i) > dzphot_thresh) or (zspec_i > zhi and abs(dzphot_i) > dzphot_thresh):
            truefalse_array.append(1)
        if (zspec_i < zlo and abs(dzphot_i) <= dzphot_thresh) or (zspec_i > zhi and abs(dzphot_i) <= dzphot_thresh):
            truefalse_array.append(2)
        if zlo < zspec_i < zhi and abs(dzphot_i) > dzphot_thresh:
            truefalse_array.append(3)

    truefalse_array = pylab.array(truefalse_array)


    #pylab.plot(lmass, vor, 'ko', ms=1)
    h = pylab.hist(vor, histtype='step', lw=2, bins=1000, cumulative=1, normed=1, label=f.name)
    #h = pylab.hist(vor, histtype='step', lw=2, bins=4, range=(0,2), label=f.name)



    inds_truepos = pylab.find(truefalse_array == 0)
    inds_trueneg = pylab.find(truefalse_array == 1)
    inds_falsepos = pylab.find(truefalse_array == 2)
    inds_falseneg = pylab.find(truefalse_array == 3)

    zspec = pylab.array(zspec)
    zphot = pylab.array(zphot)
    lmass = pylab.array(lmass)
    vor = pylab.array(vor)
    zlss = pylab.array(zlss)

    dzspec = (zspec - zlss) / (1 + zspec)
    dzphot = (zphot - zlss) / (1 + zphot)



    ###  binning by mass
    digi_truepos_lmass = pylab.digitize(lmass[inds_truepos], lmassbins)
    digi_trueneg_lmass = pylab.digitize(lmass[inds_trueneg], lmassbins)
    digi_falsepos_lmass = pylab.digitize(lmass[inds_falsepos], lmassbins)
    digi_falseneg_lmass = pylab.digitize(lmass[inds_falseneg], lmassbins)

    bincount_truepos_lmass = pylab.bincount(digi_truepos_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_trueneg_lmass = pylab.bincount(digi_trueneg_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_falsepos_lmass = pylab.bincount(digi_falsepos_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_falseneg_lmass = pylab.bincount(digi_falseneg_lmass, minlength=len(lmassbins)+1)[1:-1]

    n_truepos_lmass[i] += bincount_truepos_lmass
    n_trueneg_lmass[i] += bincount_trueneg_lmass
    n_falsepos_lmass[i] += bincount_falsepos_lmass
    n_falseneg_lmass[i] += bincount_falseneg_lmass



    ###  binning by overdensity
    digi_truepos_overdens = pylab.digitize(vor[inds_truepos], overdens_bins)
    digi_trueneg_overdens = pylab.digitize(vor[inds_trueneg], overdens_bins)
    digi_falsepos_overdens = pylab.digitize(vor[inds_falsepos], overdens_bins)
    digi_falseneg_overdens = pylab.digitize(vor[inds_falseneg], overdens_bins)

    bincount_truepos_overdens = pylab.bincount(digi_truepos_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_trueneg_overdens = pylab.bincount(digi_trueneg_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_falsepos_overdens = pylab.bincount(digi_falsepos_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_falseneg_overdens = pylab.bincount(digi_falseneg_overdens, minlength=len(overdens_bins)+1)[1:-1]

    n_truepos_overdens[i] += bincount_truepos_overdens
    n_trueneg_overdens[i] += bincount_trueneg_overdens
    n_falsepos_overdens[i] += bincount_falsepos_overdens
    n_falseneg_overdens[i] += bincount_falseneg_overdens


    for overdensi in range(1, len(overdens_bins)):

        inds1 = pylab.find(digi_truepos_overdens == overdensi)
        inds2 = pylab.find(digi_trueneg_overdens == overdensi)
        inds3 = pylab.find(digi_falsepos_overdens == overdensi)
        inds4 = pylab.find(digi_falseneg_overdens == overdensi)

        n_truepos_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_truepos_lmass[inds1], minlength=len(lmassbins)+1)[1:-1]
        n_trueneg_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_trueneg_lmass[inds2], minlength=len(lmassbins)+1)[1:-1]
        n_falsepos_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_falsepos_lmass[inds3], minlength=len(lmassbins)+1)[1:-1]
        n_falseneg_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_falseneg_lmass[inds4], minlength=len(lmassbins)+1)[1:-1]


pylab.legend(loc=2)
print 'done with +/- 1.5 sigma_zphot'
















###  frac_falsepos/frac_falseneg vs lmass
fig = pylab.figure(figsize=(16., 8.8))
sp3 = pylab.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
sp4 = pylab.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)
sp1 = pylab.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
sp2 = pylab.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1)

fig.subplots_adjust(hspace=0, wspace=0)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp4.minorticks_on()
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('fraction [%]')
sp2.set_ylabel('correction factor')
sp4.set_xlabel('log( 1 + $\delta_{gal}$ )')

sp1.grid()
sp2.grid()
sp3.grid()
sp4.grid()
sp1.axis([9.2, 11.8, -3, 85])
sp2.axis([9.2, 11.8, 0.25, 1.35])
sp3.axis([-0.1, 2.1, -3, 85])
sp4.axis([-0.1, 2.1, 0.25, 1.35])





###  plotting false +/- fractions as a function of lmass
nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass.sum(axis=0)])

frac_falsepos_all = n_falsepos_lmass.sum(axis=0) * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
frac_falseneg_all = n_falseneg_lmass.sum(axis=0) * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

c1, c2 = 'k', 'w'

sp1.errorbar(lmassbars+0.006, frac_falsepos_all, xerr=dm/2., yerr=[elo_falsepos_all, ehi_falsepos_all],
             ls='', marker='o', mew=2, ms=9, color=c1, mfc=c1, ecolor='k', capsize=0, label='false +')

sp1.errorbar(lmassbars-0.006, frac_falseneg_all, xerr=dm/2., yerr=[elo_falseneg_all, ehi_falseneg_all],
             ls='', marker='s', mew=2, ms=9, color=c2, mfc=c2, ecolor='k', capsize=0, label='false  -')

sp1.legend(loc=1, numpoints=1)


corr_factor = 1 - frac_falsepos_all/100. + frac_falseneg_all/100.
elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

sp2.errorbar(lmassbars, corr_factor, yerr=[elo_corrfactor, ehi_corrfactor],
             ls='', marker='o', ms=10, mew=2, color='k', mfc='#666666', ecolor='#666666', capsize=0)





###  ... as a function of overdensity: log(1 + delta)
nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_overdens.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_overdens.sum(axis=0)])

frac_falsepos_all = n_falsepos_overdens.sum(axis=0) * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
frac_falseneg_all = n_falseneg_overdens.sum(axis=0) * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

sp3.errorbar(overdens_bars+0.006, frac_falsepos_all, xerr=0.5/2., yerr=[elo_falsepos_all, ehi_falsepos_all],
             ls='', marker='o', mew=2, ms=9, color=c1, mfc=c1, ecolor='k', capsize=0, label='false +')

sp3.errorbar(overdens_bars-0.006, frac_falseneg_all, xerr=0.5/2., yerr=[elo_falseneg_all, ehi_falseneg_all],
             ls='', marker='s', mew=2, ms=9, color=c2, mfc=c2, ecolor='k', capsize=0, label='false  -')

sp3.legend(loc=1, numpoints=1)


corr_factor = 1 - frac_falsepos_all/100. + frac_falseneg_all/100.
elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

sp4.errorbar(overdens_bars, corr_factor, yerr=[elo_corrfactor, ehi_corrfactor],
             ls='', marker='o', ms=10, mew=2, color='k', mfc='#666666', ecolor='#666666', capsize=0)


outname = '../output/falsefractions.png'
fig.savefig(outname)
print 'wrote to: %s\n' % outname


###  printing correction factors
print '\ncorr_factors'
c_string = ''
for i in range(len(corr_factor)):
    c_string += '%.3f, ' % corr_factor[i]
print '[ ' + c_string[:-2] + '], \n\n'
























###  plotting fractions vs lmass for each overdens bin

pdf = PdfPages('../output/falsefractions_lmass_overdens.pdf')





fig = pylab.figure(figsize=(16.5, 8.4))
sp1 = pylab.subplot2grid((4, 9), (0, 0), rowspan=2, colspan=4)
sp2 = pylab.subplot2grid((4, 9), (2, 0), rowspan=2, colspan=4)
sp3 = pylab.subplot2grid((4, 9), (0, 5), rowspan=4, colspan=4)

fig.subplots_adjust(hspace=0, wspace=0, bottom=0.08, left=0.06)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('fraction')
sp2.set_ylabel('fraction')
sp3.set_ylabel('correction factor')

#sp1.grid()
#sp2.grid()
#sp3.grid()
sp1.axis([9.2, 11.8, -0.07, 0.9])
sp2.axis([9.2, 11.8, -0.07, 0.9])
sp3.axis([9.2, 11.8, 0, 2])

frac_falsepos_all = n_falsepos_lmass.sum(axis=0) * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
frac_falseneg_all = n_falseneg_lmass.sum(axis=0) * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass.sum(axis=0)])

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

sp1.fill_between(lmassbars, frac_falsepos_all-elo_falsepos_all, frac_falsepos_all+ehi_falsepos_all, color='#a6a6a6', zorder=1)
sp2.fill_between(lmassbars, frac_falseneg_all-elo_falseneg_all, frac_falseneg_all+ehi_falseneg_all, color='#a6a6a6', zorder=1)


t = sp1.text(0.8, 0.9, 'false +', transform=sp1.transAxes, fontweight='bold')
t = sp2.text(0.8, 0.9, 'false -', transform=sp2.transAxes, fontweight='bold')

corr_factors = pylab.zeros((len(overdens_bars), len(lmassbars)))
elo_corrfactors = pylab.zeros((len(overdens_bars), len(lmassbars)))
ehi_corrfactors = pylab.zeros((len(overdens_bars), len(lmassbars)))
offie = [-0.015, -0.005, 0.005, 0.015]
for j in range(4):

    frac_falsepos_all = n_falsepos_lmass_overdens.sum(axis=0)[j] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    frac_falseneg_all = n_falseneg_lmass_overdens.sum(axis=0)[j] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])

    nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass_overdens.sum(axis=0)[j]])
    nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass_overdens.sum(axis=0)[j]])

    ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])


    color_i = pylab.cm.cool(j/3.)

    sp1.errorbar(lmassbars+offie[j], frac_falsepos_all, xerr=dm/200., yerr=[elo_falsepos_all, ehi_falsepos_all],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])

    sp2.errorbar(lmassbars+offie[j], frac_falseneg_all, xerr=dm/200., yerr=[elo_falseneg_all, ehi_falseneg_all],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])


    corr_factor = 1 - frac_falsepos_all/1. + frac_falseneg_all/1.
    elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    corr_factors[j] += corr_factor
    elo_corrfactors[j] += elo_corrfactor
    ehi_corrfactors[j] += ehi_corrfactor

    sp3.errorbar(lmassbars+offie[j], corr_factor, xerr=dm/200., yerr=[elo_corrfactor, ehi_corrfactor],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])

sp3.legend(loc=3, fontsize=14, numpoints=1)





pdf.savefig()
sp1.clear()
sp2.clear()
sp3.clear()

pdf.close()
pylab.close()

###  printing correction factors
print '\ncorr_factors (voronoi bins)'
for i in range(len(corr_factors)):
    c_string = ''
    for j in range(len(corr_factors[i])):
        c_string += '%.3f, ' % corr_factors[i][j]
    print '[ ' + c_string[:-2] + '], '
print '\nelo_corrfactos (voronoi bins)'
for i in range(len(elo_corrfactors)):
    clo_string = ''
    for j in range(len(elo_corrfactors[i])):
        clo_string += '%.3f, ' % elo_corrfactors[i][j]
    print '[ ' + clo_string[:-2] + '], '
print '\nehi_corrfactos (voronoi bins)'
for i in range(len(ehi_corrfactors)):
    chi_string = ''
    for j in range(len(elo_corrfactors[i])):
        chi_string += '%.3f, ' % ehi_corrfactors[i][j]
    print '[ ' + chi_string[:-2] + '], '


print '\nwrote to: ../output/falsefractions_lmass_overdens.pdf'



print 'done with false +/-\n'




























































##############################
##############################
##############################
###  STAR-FORMING GALAXIES
##############################
##############################
##############################




zmag_thresh = 23.5

dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.5+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

overdens_bins = pylab.array([0., 0.5, 1., 1.5, 2.0])
overdens_bars = (overdens_bins[1:] + overdens_bins[:-1]) / 2.
voronoi_labels = ['0.0 < log(1+$\delta$) < 0.5', '0.5 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 1.5', '1.5 < log(1+$\delta$) < 2.0']




###   1.5 sigma_zphot
###   (zphot-zlss)/(1+zphot)  vs  (zspec-zlss)/(1+zspec)

n_truepos_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_trueneg_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_falsepos_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_falseneg_lmass = pylab.zeros((len(fields), len(lmassbars)))

n_truepos_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_trueneg_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_falsepos_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_falseneg_overdens = pylab.zeros((len(fields), len(overdens_bars)))

n_truepos_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_trueneg_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_falsepos_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_falseneg_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))

for i in range(len(fields)):

    f = fields[i]

    dzphot_thresh = 1.5 * f.sigmaz

    ###  checking if galaxy is in spatial selection regions
    zinds = pylab.find((f.cat.z_spec[f.inds_spatial] > 0) & 
                       (f.zmagnitude[f.inds_spatial] <= zmag_thresh) &
                       (f.fout.lmass[f.inds_spatial] > lmassbins[0]) &
                       (f.uvj_class[f.inds_spatial] == 1))

    zinds = f.inds_spatial[zinds]

    ra_spec = f.cat.ra[zinds]
    dec_spec = f.cat.dec[zinds]


    ###  Identifying which LSS substructure a galaxy belongs to.
    ###  This is simply defined as a galaxy needs to be <=0.5Mpc
    ###  of any substructure center. If there are two such
    ###  substructures, then the nearest is chosen.

    zspec = []
    zphot = []
    zlss = []
    lmass = []
    vor = []
    truefalse_array = []

    for i_gal in range(len(zinds)):

        ra0, dec0 = ra_spec[i_gal], dec_spec[i_gal]
        zspec_i = f.cat.z_spec[zinds][i_gal]
        zphot_i = f.zout.z_peak[zinds][i_gal]

        ###  First, a galaxy is deemed a member of the overall
        ###  LSS if it lies within 3vsigma of at least 1 substructure
        lss_member_status = False
        substructure_member_status = []
        rproj_list = []
        for ss in f.substructures:
            rproj_arcsec = mypy.radec_sep(ra0, dec0, ss.ra, ss.dec)
            asecperkpc = cosmo.arcsec_per_kpc_proper(ss.zmean)
            rproj_kpc = rproj_arcsec / asecperkpc.value
            rproj_list.append(rproj_kpc)

            if ss.zlo_3sig < zspec_i < ss.zhi_3sig:
                lss_member_status = True
                substructure_member_status.append(True)
            else:
                substructure_member_status.append(False)


        ###  if the galaxy is >2Mpc away from all substructures skip it
        #if min(rproj_list) >= 2000: continue


        ###  If it is an LSS member, identify which substructure
        ###  it is a member of. If only 1 substructure matches
        ###  then this is easy...
        rproj_list = pylab.array(rproj_list)
        substructure_member_status = pylab.array(substructure_member_status)
        ssinds = pylab.find(substructure_member_status)

        if len(ssinds) == 1:
            zlo, zhi = f.substructures[ssinds[0]].zlo_3sig, f.substructures[ssinds[0]].zhi_3sig
            zlss_i = f.substructures[ssinds[0]].zmean

        elif len(ssinds) > 1:
            sub_rproj_list = rproj_list[ssinds]
            rind = ssinds[sub_rproj_list.tolist().index(min(sub_rproj_list))]
            zlo, zhi = f.substructures[rind].zlo_3sig, f.substructures[rind].zhi_3sig
            zlss_i = f.substructures[rind].zmean

        elif len(ssinds) == 0:
            rind = rproj_list.tolist().index(min(rproj_list))
            zlo, zhi = f.substructures[rind].zlo_3sig, f.substructures[rind].zhi_3sig
            zlss_i = f.substructures[rind].zmean


        ###  Finding voronoi zslice that is "centered" on the LSS
        zdif = abs(f.voronoi.zbars - zlss_i)
        vi = pylab.find(zdif == zdif.min())[0]
        vor.append(f.voronoi.overdens_matrix[zinds][i_gal][vi])


        ### now that we've grabbed the galaxy, determine if it is:
        ###   0 = true +
        ###   1 = true -
        ###   2 = false +
        ###   3 = false -
        zspec.append(zspec_i)
        zphot.append(zphot_i)
        zlss.append(zlss_i)
        lmass.append(f.fout.lmass[zinds][i_gal])
        dzphot_i = (zphot_i - zlss_i) / (1 + zphot_i)
        if zlo < zspec_i < zhi and abs(dzphot_i) <= dzphot_thresh:
            truefalse_array.append(0)
        if (zspec_i < zlo and abs(dzphot_i) > dzphot_thresh) or (zspec_i > zhi and abs(dzphot_i) > dzphot_thresh):
            truefalse_array.append(1)
        if (zspec_i < zlo and abs(dzphot_i) <= dzphot_thresh) or (zspec_i > zhi and abs(dzphot_i) <= dzphot_thresh):
            truefalse_array.append(2)
        if zlo < zspec_i < zhi and abs(dzphot_i) > dzphot_thresh:
            truefalse_array.append(3)

    truefalse_array = pylab.array(truefalse_array)


    #pylab.plot(lmass, vor, 'ko', ms=1)
    #h = pylab.hist(vor, histtype='step', lw=2, bins=1000, cumulative=1, normed=1, label=f.name)
    #h = pylab.hist(vor, histtype='step', lw=2, bins=4, range=(0,2), label=f.name)



    inds_truepos = pylab.find(truefalse_array == 0)
    inds_trueneg = pylab.find(truefalse_array == 1)
    inds_falsepos = pylab.find(truefalse_array == 2)
    inds_falseneg = pylab.find(truefalse_array == 3)

    zspec = pylab.array(zspec)
    zphot = pylab.array(zphot)
    lmass = pylab.array(lmass)
    vor = pylab.array(vor)
    zlss = pylab.array(zlss)

    dzspec = (zspec - zlss) / (1 + zspec)
    dzphot = (zphot - zlss) / (1 + zphot)



    ###  binning by mass
    digi_truepos_lmass = pylab.digitize(lmass[inds_truepos], lmassbins)
    digi_trueneg_lmass = pylab.digitize(lmass[inds_trueneg], lmassbins)
    digi_falsepos_lmass = pylab.digitize(lmass[inds_falsepos], lmassbins)
    digi_falseneg_lmass = pylab.digitize(lmass[inds_falseneg], lmassbins)

    bincount_truepos_lmass = pylab.bincount(digi_truepos_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_trueneg_lmass = pylab.bincount(digi_trueneg_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_falsepos_lmass = pylab.bincount(digi_falsepos_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_falseneg_lmass = pylab.bincount(digi_falseneg_lmass, minlength=len(lmassbins)+1)[1:-1]

    n_truepos_lmass[i] += bincount_truepos_lmass
    n_trueneg_lmass[i] += bincount_trueneg_lmass
    n_falsepos_lmass[i] += bincount_falsepos_lmass
    n_falseneg_lmass[i] += bincount_falseneg_lmass



    ###  binning by overdensity
    digi_truepos_overdens = pylab.digitize(vor[inds_truepos], overdens_bins)
    digi_trueneg_overdens = pylab.digitize(vor[inds_trueneg], overdens_bins)
    digi_falsepos_overdens = pylab.digitize(vor[inds_falsepos], overdens_bins)
    digi_falseneg_overdens = pylab.digitize(vor[inds_falseneg], overdens_bins)

    bincount_truepos_overdens = pylab.bincount(digi_truepos_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_trueneg_overdens = pylab.bincount(digi_trueneg_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_falsepos_overdens = pylab.bincount(digi_falsepos_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_falseneg_overdens = pylab.bincount(digi_falseneg_overdens, minlength=len(overdens_bins)+1)[1:-1]

    n_truepos_overdens[i] += bincount_truepos_overdens
    n_trueneg_overdens[i] += bincount_trueneg_overdens
    n_falsepos_overdens[i] += bincount_falsepos_overdens
    n_falseneg_overdens[i] += bincount_falseneg_overdens


    for overdensi in range(1, len(overdens_bins)):

        inds1 = pylab.find(digi_truepos_overdens == overdensi)
        inds2 = pylab.find(digi_trueneg_overdens == overdensi)
        inds3 = pylab.find(digi_falsepos_overdens == overdensi)
        inds4 = pylab.find(digi_falseneg_overdens == overdensi)

        n_truepos_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_truepos_lmass[inds1], minlength=len(lmassbins)+1)[1:-1]
        n_trueneg_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_trueneg_lmass[inds2], minlength=len(lmassbins)+1)[1:-1]
        n_falsepos_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_falsepos_lmass[inds3], minlength=len(lmassbins)+1)[1:-1]
        n_falseneg_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_falseneg_lmass[inds4], minlength=len(lmassbins)+1)[1:-1]


#pylab.legend(loc=2)
print 'done with +/- 1.5 sigma_zphot'











###  frac_falsepos/frac_falseneg vs lmass
fig = pylab.figure(figsize=(16., 8.8))
sp3 = pylab.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
sp4 = pylab.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)
sp1 = pylab.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
sp2 = pylab.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1)

fig.subplots_adjust(hspace=0, wspace=0)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp4.minorticks_on()
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('fraction [%]')
sp2.set_ylabel('correction factor')
sp4.set_xlabel('log( 1 + $\delta_{gal}$ )')

sp1.grid()
sp2.grid()
sp3.grid()
sp4.grid()
sp1.axis([9.2, 11.8, -3, 85])
sp2.axis([9.2, 11.8, 0.25, 1.35])
sp3.axis([-0.1, 2.1, -3, 85])
sp4.axis([-0.1, 2.1, 0.25, 1.35])





###  plotting false +/- fractions as a function of lmass
nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass.sum(axis=0)])

frac_falsepos_all = n_falsepos_lmass.sum(axis=0) * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
frac_falseneg_all = n_falseneg_lmass.sum(axis=0) * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

c1, c2 = 'k', 'w'

sp1.errorbar(lmassbars+0.006, frac_falsepos_all, xerr=dm/2., yerr=[elo_falsepos_all, ehi_falsepos_all],
             ls='', marker='o', mew=2, ms=9, color=c1, mfc=c1, ecolor='k', capsize=0, label='false +')

sp1.errorbar(lmassbars-0.006, frac_falseneg_all, xerr=dm/2., yerr=[elo_falseneg_all, ehi_falseneg_all],
             ls='', marker='s', mew=2, ms=9, color=c2, mfc=c2, ecolor='k', capsize=0, label='false  -')

sp1.legend(loc=1, numpoints=1)


corr_factor = 1 - frac_falsepos_all/100. + frac_falseneg_all/100.
elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

sp2.errorbar(lmassbars, corr_factor, yerr=[elo_corrfactor, ehi_corrfactor],
             ls='', marker='o', ms=10, mew=2, color='k', mfc='#666666', ecolor='#666666', capsize=0)


###  printing correction factors
print '\ncorr_factors'
c_string = ''
for i in range(len(corr_factor)):
    c_string += '%.3f, ' % corr_factor[i]
print '[ ' + c_string[:-2] + '], \n\n'





###  ... as a function of overdensity: log(1 + delta)
nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_overdens.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_overdens.sum(axis=0)])

frac_falsepos_all = n_falsepos_overdens.sum(axis=0) * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
frac_falseneg_all = n_falseneg_overdens.sum(axis=0) * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

sp3.errorbar(overdens_bars+0.006, frac_falsepos_all, xerr=0.5/2., yerr=[elo_falsepos_all, ehi_falsepos_all],
             ls='', marker='o', mew=2, ms=9, color=c1, mfc=c1, ecolor='k', capsize=0, label='false +')

sp3.errorbar(overdens_bars-0.006, frac_falseneg_all, xerr=0.5/2., yerr=[elo_falseneg_all, ehi_falseneg_all],
             ls='', marker='s', mew=2, ms=9, color=c2, mfc=c2, ecolor='k', capsize=0, label='false  -')

sp3.legend(loc=1, numpoints=1)


corr_factor = 1 - frac_falsepos_all/100. + frac_falseneg_all/100.
elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

sp4.errorbar(overdens_bars, corr_factor, yerr=[elo_corrfactor, ehi_corrfactor],
             ls='', marker='o', ms=10, mew=2, color='k', mfc='#666666', ecolor='#666666', capsize=0)


outname = '../output/falsefractions_sf.png'
fig.savefig(outname)
print 'wrote to: %s\n' % outname















###  plotting fractions vs lmass for each overdens bin

pdf = PdfPages('../output/falsefractions_lmass_overdens_sf.pdf')





fig = pylab.figure(figsize=(16.5, 8.4))
sp1 = pylab.subplot2grid((4, 9), (0, 0), rowspan=2, colspan=4)
sp2 = pylab.subplot2grid((4, 9), (2, 0), rowspan=2, colspan=4)
sp3 = pylab.subplot2grid((4, 9), (0, 5), rowspan=4, colspan=4)

fig.subplots_adjust(hspace=0, wspace=0, bottom=0.08, left=0.06)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('fraction')
sp2.set_ylabel('fraction')
sp3.set_ylabel('correction factor')

#sp1.grid()
#sp2.grid()
#sp3.grid()
sp1.axis([9.2, 11.8, -0.07, 0.9])
sp2.axis([9.2, 11.8, -0.07, 0.9])
sp3.axis([9.2, 11.8, 0, 2])

frac_falsepos_all = n_falsepos_lmass.sum(axis=0) * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
frac_falseneg_all = n_falseneg_lmass.sum(axis=0) * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass.sum(axis=0)])

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

sp1.fill_between(lmassbars, frac_falsepos_all-elo_falsepos_all, frac_falsepos_all+ehi_falsepos_all, color='#a6a6a6', zorder=1)
sp2.fill_between(lmassbars, frac_falseneg_all-elo_falseneg_all, frac_falseneg_all+ehi_falseneg_all, color='#a6a6a6', zorder=1)


t = sp1.text(0.8, 0.9, 'false +', transform=sp1.transAxes, fontweight='bold')
t = sp2.text(0.8, 0.9, 'false -', transform=sp2.transAxes, fontweight='bold')

corr_factors = pylab.zeros((len(overdens_bars), len(lmassbars)))
elo_corrfactors = pylab.zeros((len(overdens_bars), len(lmassbars)))
ehi_corrfactors = pylab.zeros((len(overdens_bars), len(lmassbars)))
offie = [-0.015, -0.005, 0.005, 0.015]
for j in range(4):

    frac_falsepos_all = n_falsepos_lmass_overdens.sum(axis=0)[j] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    frac_falseneg_all = n_falseneg_lmass_overdens.sum(axis=0)[j] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])

    nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass_overdens.sum(axis=0)[j]])
    nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass_overdens.sum(axis=0)[j]])

    ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])


    color_i = pylab.cm.cool(j/3.)

    sp1.errorbar(lmassbars+offie[j], frac_falsepos_all, xerr=dm/200., yerr=[elo_falsepos_all, ehi_falsepos_all],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])

    sp2.errorbar(lmassbars+offie[j], frac_falseneg_all, xerr=dm/200., yerr=[elo_falseneg_all, ehi_falseneg_all],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])


    corr_factor = 1 - frac_falsepos_all/1. + frac_falseneg_all/1.
    elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    corr_factors[j] += corr_factor
    elo_corrfactors[j] += elo_corrfactor
    ehi_corrfactors[j] += ehi_corrfactor

    sp3.errorbar(lmassbars+offie[j], corr_factor, xerr=dm/200., yerr=[elo_corrfactor, ehi_corrfactor],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])

sp3.legend(loc=3, fontsize=14, numpoints=1)





pdf.savefig()
sp1.clear()
sp2.clear()
sp3.clear()

pdf.close()
pylab.close()

###  printing correction factors
print '\ncorr_factors (voronoi bins)'
for i in range(len(corr_factors)):
    c_string = ''
    for j in range(len(corr_factors[i])):
        c_string += '%.3f, ' % corr_factors[i][j]
    print '[ ' + c_string[:-2] + '], '
print '\nelo_corrfactos (voronoi bins)'
for i in range(len(elo_corrfactors)):
    clo_string = ''
    for j in range(len(elo_corrfactors[i])):
        clo_string += '%.3f, ' % elo_corrfactors[i][j]
    print '[ ' + clo_string[:-2] + '], '
print '\nehi_corrfactos (voronoi bins)'
for i in range(len(ehi_corrfactors)):
    chi_string = ''
    for j in range(len(elo_corrfactors[i])):
        chi_string += '%.3f, ' % ehi_corrfactors[i][j]
    print '[ ' + chi_string[:-2] + '], '


print '\nwrote to: ../output/falsefractions_lmass_overdens_sf.pdf'



print 'done with false +/-\n'




























































##############################
##############################
##############################
###  QUIESCENT GALAXIES
##############################
##############################
##############################




zmag_thresh = 23.5

dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.5+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

overdens_bins = pylab.array([0., 0.5, 1., 1.5, 2.0])
overdens_bars = (overdens_bins[1:] + overdens_bins[:-1]) / 2.
voronoi_labels = ['0.0 < log(1+$\delta$) < 0.5', '0.5 < log(1+$\delta$) < 1.0', '1.0 < log(1+$\delta$) < 1.5', '1.5 < log(1+$\delta$) < 2.0']




###   1.5 sigma_zphot
###   (zphot-zlss)/(1+zphot)  vs  (zspec-zlss)/(1+zspec)

n_truepos_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_trueneg_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_falsepos_lmass = pylab.zeros((len(fields), len(lmassbars)))
n_falseneg_lmass = pylab.zeros((len(fields), len(lmassbars)))

n_truepos_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_trueneg_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_falsepos_overdens = pylab.zeros((len(fields), len(overdens_bars)))
n_falseneg_overdens = pylab.zeros((len(fields), len(overdens_bars)))

n_truepos_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_trueneg_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_falsepos_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))
n_falseneg_lmass_overdens = pylab.zeros((len(fields), len(overdens_bars), len(lmassbars)))

for i in range(len(fields)):

    f = fields[i]

    dzphot_thresh = 1.5 * f.sigmaz

    ###  checking if galaxy is in spatial selection regions
    zinds = pylab.find((f.cat.z_spec[f.inds_spatial] > 0) & 
                       (f.zmagnitude[f.inds_spatial] <= zmag_thresh) &
                       (f.fout.lmass[f.inds_spatial] > lmassbins[0]) &
                       (f.uvj_class[f.inds_spatial] == 0))

    zinds = f.inds_spatial[zinds]

    ra_spec = f.cat.ra[zinds]
    dec_spec = f.cat.dec[zinds]


    ###  Identifying which LSS substructure a galaxy belongs to.
    ###  This is simply defined as a galaxy needs to be <=0.5Mpc
    ###  of any substructure center. If there are two such
    ###  substructures, then the nearest is chosen.

    zspec = []
    zphot = []
    zlss = []
    lmass = []
    vor = []
    truefalse_array = []

    for i_gal in range(len(zinds)):

        ra0, dec0 = ra_spec[i_gal], dec_spec[i_gal]
        zspec_i = f.cat.z_spec[zinds][i_gal]
        zphot_i = f.zout.z_peak[zinds][i_gal]

        ###  First, a galaxy is deemed a member of the overall
        ###  LSS if it lies within 3vsigma of at least 1 substructure
        lss_member_status = False
        substructure_member_status = []
        rproj_list = []
        for ss in f.substructures:
            rproj_arcsec = mypy.radec_sep(ra0, dec0, ss.ra, ss.dec)
            asecperkpc = cosmo.arcsec_per_kpc_proper(ss.zmean)
            rproj_kpc = rproj_arcsec / asecperkpc.value
            rproj_list.append(rproj_kpc)

            if ss.zlo_3sig < zspec_i < ss.zhi_3sig:
                lss_member_status = True
                substructure_member_status.append(True)
            else:
                substructure_member_status.append(False)


        ###  if the galaxy is >2Mpc away from all substructures skip it
        #if min(rproj_list) >= 2000: continue


        ###  If it is an LSS member, identify which substructure
        ###  it is a member of. If only 1 substructure matches
        ###  then this is easy...
        rproj_list = pylab.array(rproj_list)
        substructure_member_status = pylab.array(substructure_member_status)
        ssinds = pylab.find(substructure_member_status)

        if len(ssinds) == 1:
            zlo, zhi = f.substructures[ssinds[0]].zlo_3sig, f.substructures[ssinds[0]].zhi_3sig
            zlss_i = f.substructures[ssinds[0]].zmean

        elif len(ssinds) > 1:
            sub_rproj_list = rproj_list[ssinds]
            rind = ssinds[sub_rproj_list.tolist().index(min(sub_rproj_list))]
            zlo, zhi = f.substructures[rind].zlo_3sig, f.substructures[rind].zhi_3sig
            zlss_i = f.substructures[rind].zmean

        elif len(ssinds) == 0:
            rind = rproj_list.tolist().index(min(rproj_list))
            zlo, zhi = f.substructures[rind].zlo_3sig, f.substructures[rind].zhi_3sig
            zlss_i = f.substructures[rind].zmean


        ###  Finding voronoi zslice that is "centered" on the LSS
        zdif = abs(f.voronoi.zbars - zlss_i)
        vi = pylab.find(zdif == zdif.min())[0]
        vor.append(f.voronoi.overdens_matrix[zinds][i_gal][vi])


        ### now that we've grabbed the galaxy, determine if it is:
        ###   0 = true +
        ###   1 = true -
        ###   2 = false +
        ###   3 = false -
        zspec.append(zspec_i)
        zphot.append(zphot_i)
        zlss.append(zlss_i)
        lmass.append(f.fout.lmass[zinds][i_gal])
        dzphot_i = (zphot_i - zlss_i) / (1 + zphot_i)
        if zlo < zspec_i < zhi and abs(dzphot_i) <= dzphot_thresh:
            truefalse_array.append(0)
        if (zspec_i < zlo and abs(dzphot_i) > dzphot_thresh) or (zspec_i > zhi and abs(dzphot_i) > dzphot_thresh):
            truefalse_array.append(1)
        if (zspec_i < zlo and abs(dzphot_i) <= dzphot_thresh) or (zspec_i > zhi and abs(dzphot_i) <= dzphot_thresh):
            truefalse_array.append(2)
        if zlo < zspec_i < zhi and abs(dzphot_i) > dzphot_thresh:
            truefalse_array.append(3)

    truefalse_array = pylab.array(truefalse_array)


    #pylab.plot(lmass, vor, 'ko', ms=1)
    #h = pylab.hist(vor, histtype='step', lw=2, bins=1000, cumulative=1, normed=1, label=f.name)
    #h = pylab.hist(vor, histtype='step', lw=2, bins=4, range=(0,2), label=f.name)



    inds_truepos = pylab.find(truefalse_array == 0)
    inds_trueneg = pylab.find(truefalse_array == 1)
    inds_falsepos = pylab.find(truefalse_array == 2)
    inds_falseneg = pylab.find(truefalse_array == 3)

    zspec = pylab.array(zspec)
    zphot = pylab.array(zphot)
    lmass = pylab.array(lmass)
    vor = pylab.array(vor)
    zlss = pylab.array(zlss)

    dzspec = (zspec - zlss) / (1 + zspec)
    dzphot = (zphot - zlss) / (1 + zphot)



    ###  binning by mass
    digi_truepos_lmass = pylab.digitize(lmass[inds_truepos], lmassbins)
    digi_trueneg_lmass = pylab.digitize(lmass[inds_trueneg], lmassbins)
    digi_falsepos_lmass = pylab.digitize(lmass[inds_falsepos], lmassbins)
    digi_falseneg_lmass = pylab.digitize(lmass[inds_falseneg], lmassbins)

    bincount_truepos_lmass = pylab.bincount(digi_truepos_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_trueneg_lmass = pylab.bincount(digi_trueneg_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_falsepos_lmass = pylab.bincount(digi_falsepos_lmass, minlength=len(lmassbins)+1)[1:-1]
    bincount_falseneg_lmass = pylab.bincount(digi_falseneg_lmass, minlength=len(lmassbins)+1)[1:-1]

    n_truepos_lmass[i] += bincount_truepos_lmass
    n_trueneg_lmass[i] += bincount_trueneg_lmass
    n_falsepos_lmass[i] += bincount_falsepos_lmass
    n_falseneg_lmass[i] += bincount_falseneg_lmass



    ###  binning by overdensity
    digi_truepos_overdens = pylab.digitize(vor[inds_truepos], overdens_bins)
    digi_trueneg_overdens = pylab.digitize(vor[inds_trueneg], overdens_bins)
    digi_falsepos_overdens = pylab.digitize(vor[inds_falsepos], overdens_bins)
    digi_falseneg_overdens = pylab.digitize(vor[inds_falseneg], overdens_bins)

    bincount_truepos_overdens = pylab.bincount(digi_truepos_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_trueneg_overdens = pylab.bincount(digi_trueneg_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_falsepos_overdens = pylab.bincount(digi_falsepos_overdens, minlength=len(overdens_bins)+1)[1:-1]
    bincount_falseneg_overdens = pylab.bincount(digi_falseneg_overdens, minlength=len(overdens_bins)+1)[1:-1]

    n_truepos_overdens[i] += bincount_truepos_overdens
    n_trueneg_overdens[i] += bincount_trueneg_overdens
    n_falsepos_overdens[i] += bincount_falsepos_overdens
    n_falseneg_overdens[i] += bincount_falseneg_overdens


    for overdensi in range(1, len(overdens_bins)):

        inds1 = pylab.find(digi_truepos_overdens == overdensi)
        inds2 = pylab.find(digi_trueneg_overdens == overdensi)
        inds3 = pylab.find(digi_falsepos_overdens == overdensi)
        inds4 = pylab.find(digi_falseneg_overdens == overdensi)

        n_truepos_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_truepos_lmass[inds1], minlength=len(lmassbins)+1)[1:-1]
        n_trueneg_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_trueneg_lmass[inds2], minlength=len(lmassbins)+1)[1:-1]
        n_falsepos_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_falsepos_lmass[inds3], minlength=len(lmassbins)+1)[1:-1]
        n_falseneg_lmass_overdens[i][overdensi-1] += pylab.bincount(digi_falseneg_lmass[inds4], minlength=len(lmassbins)+1)[1:-1]


#pylab.legend(loc=2)
print 'done with +/- 1.5 sigma_zphot'











###  frac_falsepos/frac_falseneg vs lmass
fig = pylab.figure(figsize=(16., 8.8))
sp3 = pylab.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
sp4 = pylab.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)
sp1 = pylab.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
sp2 = pylab.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1)

fig.subplots_adjust(hspace=0, wspace=0)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp4.minorticks_on()
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('fraction [%]')
sp2.set_ylabel('correction factor')
sp4.set_xlabel('log( 1 + $\delta_{gal}$ )')

sp1.grid()
sp2.grid()
sp3.grid()
sp4.grid()
sp1.axis([9.2, 11.8, -3, 85])
sp2.axis([9.2, 11.8, 0.25, 1.35])
sp3.axis([-0.1, 2.1, -3, 85])
sp4.axis([-0.1, 2.1, 0.25, 1.35])





###  plotting false +/- fractions as a function of lmass
nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass.sum(axis=0)])

frac_falsepos_all = n_falsepos_lmass.sum(axis=0) * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
frac_falseneg_all = n_falseneg_lmass.sum(axis=0) * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 100. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

c1, c2 = 'k', 'w'

sp1.errorbar(lmassbars+0.006, frac_falsepos_all, xerr=dm/2., yerr=[elo_falsepos_all, ehi_falsepos_all],
             ls='', marker='o', mew=2, ms=9, color=c1, mfc=c1, ecolor='k', capsize=0, label='false +')

sp1.errorbar(lmassbars-0.006, frac_falseneg_all, xerr=dm/2., yerr=[elo_falseneg_all, ehi_falseneg_all],
             ls='', marker='s', mew=2, ms=9, color=c2, mfc=c2, ecolor='k', capsize=0, label='false  -')

sp1.legend(loc=1, numpoints=1)


corr_factor = 1 - frac_falsepos_all/100. + frac_falseneg_all/100.
elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

sp2.errorbar(lmassbars, corr_factor, yerr=[elo_corrfactor, ehi_corrfactor],
             ls='', marker='o', ms=10, mew=2, color='k', mfc='#666666', ecolor='#666666', capsize=0)


###  printing correction factors
print '\ncorr_factors'
c_string = ''
for i in range(len(corr_factor)):
    c_string += '%.3f, ' % corr_factor[i]
print '[ ' + c_string[:-2] + '], \n\n'





###  ... as a function of overdensity: log(1 + delta)
nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_overdens.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_overdens.sum(axis=0)])

frac_falsepos_all = n_falsepos_overdens.sum(axis=0) * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
frac_falseneg_all = n_falseneg_overdens.sum(axis=0) * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 100. / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

sp3.errorbar(overdens_bars+0.006, frac_falsepos_all, xerr=0.5/2., yerr=[elo_falsepos_all, ehi_falsepos_all],
             ls='', marker='o', mew=2, ms=9, color=c1, mfc=c1, ecolor='k', capsize=0, label='false +')

sp3.errorbar(overdens_bars-0.006, frac_falseneg_all, xerr=0.5/2., yerr=[elo_falseneg_all, ehi_falseneg_all],
             ls='', marker='s', mew=2, ms=9, color=c2, mfc=c2, ecolor='k', capsize=0, label='false  -')

sp3.legend(loc=1, numpoints=1)


corr_factor = 1 - frac_falsepos_all/100. + frac_falseneg_all/100.
elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))
ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_overdens.sum(axis=0) + n_falsepos_overdens.sum(axis=0))

sp4.errorbar(overdens_bars, corr_factor, yerr=[elo_corrfactor, ehi_corrfactor],
             ls='', marker='o', ms=10, mew=2, color='k', mfc='#666666', ecolor='#666666', capsize=0)


outname = '../output/falsefractions_qu.png'
fig.savefig(outname)
print 'wrote to: %s\n' % outname















###  plotting fractions vs lmass for each overdens bin

pdf = PdfPages('../output/falsefractions_lmass_overdens_qu.pdf')





fig = pylab.figure(figsize=(16.5, 8.4))
sp1 = pylab.subplot2grid((4, 9), (0, 0), rowspan=2, colspan=4)
sp2 = pylab.subplot2grid((4, 9), (2, 0), rowspan=2, colspan=4)
sp3 = pylab.subplot2grid((4, 9), (0, 5), rowspan=4, colspan=4)

fig.subplots_adjust(hspace=0, wspace=0, bottom=0.08, left=0.06)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()
sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp2.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp3.set_xlabel('log( M$_*$ / M$_{\odot}$ )')
sp1.set_ylabel('fraction')
sp2.set_ylabel('fraction')
sp3.set_ylabel('correction factor')

#sp1.grid()
#sp2.grid()
#sp3.grid()
sp1.axis([9.2, 11.8, -0.07, 0.9])
sp2.axis([9.2, 11.8, -0.07, 0.9])
sp3.axis([9.2, 11.8, 0, 2])

frac_falsepos_all = n_falsepos_lmass.sum(axis=0) * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
frac_falseneg_all = n_falseneg_lmass.sum(axis=0) * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass.sum(axis=0)])
nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass.sum(axis=0)])

ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))
elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 1. / (n_truepos_lmass.sum(axis=0) + n_falsepos_lmass.sum(axis=0))

sp1.fill_between(lmassbars, frac_falsepos_all-elo_falsepos_all, frac_falsepos_all+ehi_falsepos_all, color='#a6a6a6', zorder=1)
sp2.fill_between(lmassbars, frac_falseneg_all-elo_falseneg_all, frac_falseneg_all+ehi_falseneg_all, color='#a6a6a6', zorder=1)


t = sp1.text(0.8, 0.9, 'false +', transform=sp1.transAxes, fontweight='bold')
t = sp2.text(0.8, 0.9, 'false -', transform=sp2.transAxes, fontweight='bold')

corr_factors = pylab.zeros((len(overdens_bars), len(lmassbars)))
elo_corrfactors = pylab.zeros((len(overdens_bars), len(lmassbars)))
ehi_corrfactors = pylab.zeros((len(overdens_bars), len(lmassbars)))
offie = [-0.015, -0.005, 0.005, 0.015]
for j in range(4):

    frac_falsepos_all = n_falsepos_lmass_overdens.sum(axis=0)[j] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    frac_falseneg_all = n_falseneg_lmass_overdens.sum(axis=0)[j] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])

    nhinlo_falsepos_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falsepos_lmass_overdens.sum(axis=0)[j]])
    nhinlo_falseneg_tot = pylab.array([mypy.massfunc.confidence_interval(ni) for ni in n_falseneg_lmass_overdens.sum(axis=0)[j]])

    ehi_falsepos_all = nhinlo_falsepos_tot[:,0] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    elo_falsepos_all = nhinlo_falsepos_tot[:,1] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    ehi_falseneg_all = nhinlo_falseneg_tot[:,0] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    elo_falseneg_all = nhinlo_falseneg_tot[:,1] * 1. / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])


    color_i = pylab.cm.cool(j/3.)

    sp1.errorbar(lmassbars+offie[j], frac_falsepos_all, xerr=dm/200., yerr=[elo_falsepos_all, ehi_falsepos_all],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])

    sp2.errorbar(lmassbars+offie[j], frac_falseneg_all, xerr=dm/200., yerr=[elo_falseneg_all, ehi_falseneg_all],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])


    corr_factor = 1 - frac_falsepos_all/1. + frac_falseneg_all/1.
    elo_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,1]**2 + nhinlo_falseneg_tot[:,1]**2) / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    ehi_corrfactor = pylab.sqrt(nhinlo_falsepos_tot[:,0]**2 + nhinlo_falseneg_tot[:,0]**2) / (n_truepos_lmass_overdens.sum(axis=0)[j] + n_falsepos_lmass_overdens.sum(axis=0)[j])
    corr_factors[j] += corr_factor
    elo_corrfactors[j] += elo_corrfactor
    ehi_corrfactors[j] += ehi_corrfactor

    sp3.errorbar(lmassbars+offie[j], corr_factor, xerr=dm/200., yerr=[elo_corrfactor, ehi_corrfactor],
                 ls='-', marker='o', lw=1.5, mew=1.5, elinewidth=1.5, ms=9, color=color_i, mec='k', mfc=color_i, ecolor=color_i, capsize=0, label=voronoi_labels[j])

sp3.legend(loc=3, fontsize=14, numpoints=1)





pdf.savefig()
sp1.clear()
sp2.clear()
sp3.clear()

pdf.close()
pylab.close()

###  printing correction factors
print '\ncorr_factors (voronoi bins)'
for i in range(len(corr_factors)):
    c_string = ''
    for j in range(len(corr_factors[i])):
        c_string += '%.3f, ' % corr_factors[i][j]
    print '[ ' + c_string[:-2] + '], '
print '\nelo_corrfactos (voronoi bins)'
for i in range(len(elo_corrfactors)):
    clo_string = ''
    for j in range(len(elo_corrfactors[i])):
        clo_string += '%.3f, ' % elo_corrfactors[i][j]
    print '[ ' + clo_string[:-2] + '], '
print '\nehi_corrfactos (voronoi bins)'
for i in range(len(ehi_corrfactors)):
    chi_string = ''
    for j in range(len(elo_corrfactors[i])):
        chi_string += '%.3f, ' % ehi_corrfactors[i][j]
    print '[ ' + chi_string[:-2] + '], '


print '\nwrote to: ../output/falsefractions_lmass_overdens_qu.pdf'



print 'done with false +/-\n'



































