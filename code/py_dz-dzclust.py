
import os
import mypy
import time
import pylab
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

data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'



def selectionFunction_sg0023(cat):

	image = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/Cl_0023+0423/catalogs/detection_RIsup.fits'
	imwcs = wcs.WCS(image)

	###  center and radius(arcmin) of circle for selection region
	ra_center, dec_center, radius = 5.9679299, 4.371602, 8.

	inds = pylab.find(mypy.radec_sep(cat.ra, cat.dec, ra_center, dec_center) < (radius*60))


	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[6.0002043,4.443438,13.101802],
						 [6.0442164,4.4577949,14.926456],
						 [6.0607516,4.4383802,16.134604],
						 [6.0596828,4.4229565,11.242492],
						 [6.0636816,4.4075321,18.773037],
						 [6.0260797,4.4758797,8.532085],
						 [6.0111424,4.4806675,7.0114467],
						 [6.0004691,4.3984959,10.130726],
						 [6.012204,4.3894537,8.1548529],
						 [6.0130024,4.3602015,12.254177],
						 [6.0260685,4.3312143,12.638204],
						 [5.9529946,4.3602027,9.4850712],
						 [5.9439261,4.3796153,6.3277603],
						 [5.945526,4.3926459,7.6725706],
						 [5.9031293,4.2700502,18.44122],
						 [5.8625822,4.365514,11.637917],
						 [5.8505758,4.39636,11.531429],
						 [5.866039,4.4415702,11.457626],
						 [5.8908441,4.4487533,9.8495449],
						 [6.0970164,4.3796046,9.194863],
						 [6.0930157,4.3788075,4.9400946],
						 [5.9988648,4.3016978,7.724253]])

	area_arcsec2 = pylab.pi * (radius*60)**2
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2

def selectionFunction_sc1604(cat):

	image='/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/Cl_1604+4320/catalogs/CL1604_master_Rc.fits'
	imwcs = wcs.WCS(image)

	###  vertices of lines foriming polygon selection region
	ra1_a, dec1_a, ra2_a, dec2_a = 241.18999,43.319742,241.13575,43.032008
	ra1_b, dec1_b, ra2_b, dec2_b = 241.04834,43.042383,241.13575,43.032008
	ra1_c, dec1_c, ra2_c, dec2_c = 241.04834,43.042383,240.91837,43.295081
	ra1_d, dec1_d, ra2_d, dec2_d = 240.75818,43.365714,240.91837,43.295081
	ra1_e, dec1_e, ra2_e, dec2_e = 240.75818,43.365714,240.79529,43.455913
	ra1_f, dec1_f, ra2_f, dec2_f = 241.20824,43.440948,240.79529,43.455913
	ra1_g, dec1_g, ra2_g, dec2_g = 241.20824,43.440948,241.24855,43.325684
	ra1_h, dec1_h, ra2_h, dec2_h = 241.18999,43.319742,241.24855,43.325684



	(x1_a, y1_a), (x2_a, y2_a) = imwcs.wcs_world2pix(pylab.array(zip([ra1_a, ra2_a], [dec1_a, dec2_a])), 1)
	(x1_b, y1_b), (x2_b, y2_b) = imwcs.wcs_world2pix(pylab.array(zip([ra1_b, ra2_b], [dec1_b, dec2_b])), 1)
	(x1_c, y1_c), (x2_c, y2_c) = imwcs.wcs_world2pix(pylab.array(zip([ra1_c, ra2_c], [dec1_c, dec2_c])), 1)
	(x1_d, y1_d), (x2_d, y2_d) = imwcs.wcs_world2pix(pylab.array(zip([ra1_d, ra2_d], [dec1_d, dec2_d])), 1)
	(x1_e, y1_e), (x2_e, y2_e) = imwcs.wcs_world2pix(pylab.array(zip([ra1_e, ra2_e], [dec1_e, dec2_e])), 1)
	(x1_f, y1_f), (x2_f, y2_f) = imwcs.wcs_world2pix(pylab.array(zip([ra1_f, ra2_f], [dec1_f, dec2_f])), 1)
	(x1_g, y1_g), (x2_g, y2_g) = imwcs.wcs_world2pix(pylab.array(zip([ra1_g, ra2_g], [dec1_g, dec2_g])), 1)
	(x1_h, y1_h), (x2_h, y2_h) = imwcs.wcs_world2pix(pylab.array(zip([ra1_h, ra2_h], [dec1_h, dec2_h])), 1)
	xs = pylab.array([x1_a, x2_a, x1_b, x2_b, x1_c, x2_c, x1_d, x2_d, x1_e, x2_e, x1_f, x2_f, x1_g, x2_g, x1_h, x2_h])
	ys = pylab.array([y1_a, y2_a, y1_b, y2_b, y1_c, y2_c, y1_d, y2_d, y1_e, y2_e, y1_f, y2_f, y1_g, y2_g, y1_h, y2_h])

	slope_a = (y1_a - y2_a) / (x1_a - x2_a)
	slope_b = (y1_b - y2_b) / (x1_b - x2_b)
	slope_c = (y1_c - y2_c) / (x1_c - x2_c)
	slope_d = (y1_d - y2_d) / (x1_d - x2_d)
	slope_e = (y1_e - y2_e) / (x1_e - x2_e)
	slope_f = (y1_f - y2_f) / (x1_f - x2_f)
	slope_g = (y1_g - y2_g) / (x1_g - x2_g)
	slope_h = (y1_h - y2_h) / (x1_h - x2_h)

	norm_a = y1_a - slope_a * x1_a
	norm_b = y1_b - slope_b * x1_b
	norm_c = y1_c - slope_c * x1_c
	norm_d = y1_d - slope_d * x1_d
	norm_e = y1_e - slope_e * x1_e
	norm_f = y1_f - slope_f * x1_f
	norm_g = y1_g - slope_g * x1_g
	norm_h = y1_h - slope_h * x1_h

	###  selecting galaxies within a user-defind polygon enclosing the spectroscopic coverage
	inds = pylab.find((cat.y > slope_b * cat.x + norm_b) &
					  (cat.y < slope_e * cat.x + norm_e) &
					  (cat.y < slope_f * cat.x + norm_f) &
					  (cat.y < slope_g * cat.x + norm_g) &
					  ((cat.y > slope_a * cat.x + norm_a) | (cat.y > slope_h * cat.x + norm_h)) &
					  ((cat.y > slope_c * cat.x + norm_c) | (cat.y > slope_d * cat.x + norm_d)))



	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[241.05287,43.441308,25.845097],
						 [241.06541,43.41401,32.433931],
						 [241.14686,43.421204,18.177141],
						 [241.16931,43.423087,20.559791],
						 [241.14311,43.316167,33.185926],
						 [241.14619,43.265732,16.692676],
						 [241.01331,43.335491,18.023439],
						 [240.9647,43.396499,14.2458],
						 [240.90351,43.438506,17.323747],
						 [240.86569,43.357668,14.619093],
						 [241.0764,43.323762,13.507682],
						 [241.0665,43.321485,12.243407],
						 [241.16817,43.348006,10.962367],
						 [241.15879,43.351425,10.710095],
						 [240.99444,43.401458,11.073026],
						 [241.13639,43.386322,9.8980624],
						 [241.0962,43.266503,12.931588],
						 [241.05871,43.255123,11.746842],
						 [241.12848,43.252088,14.79412],
						 [241.13003,43.233886,20.697477],
						 [240.96029,43.262636,10.254993],
						 [240.96757,43.267574,8.1846382],
						 [240.88354,43.322437,13.506917],
						 [241.14769,43.195158,20.80547],
						 [241.1222,43.190997,13.509939],
						 [241.0936,43.219441,8.9580593],
						 [241.09204,43.206169,9.7404603],
						 [241.02857,43.204634,11.522894],
						 [241.11335,43.159146,13.839286],
						 [241.09931,43.099992,11.866212],
						 [241.0396,43.082537,15.065785],
						 [240.80567,43.366973,13.366928]])

	area_px2 = PolyArea(xs, ys)
	area_arcsec2 = area_px2 * px_scale**2
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2

def selectionFunction_sc1324(cat):

	print '\n##############################################################################'
	print '##  WARNING:  DO NOT USE THIS SPATIAL SELECTION FOR SPECTROSCOPIC ANALYSIS  ##'
	print '##############################################################################'

	image = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/SC_1324/catalogs/detection_i.fits'
	imwcs = wcs.WCS(image)

	###  center and radius(arcmin) of circle for selection region
	ra_center1, dec_center1, radius1 = 201.19367, 30.378647, 7.5
	ra_center2, dec_center2, radius2 = 201.19789, 30.797195, 7.5

	inds = pylab.find((mypy.radec_sep(cat.ra, cat.dec, ra_center1, dec_center1) < (radius1*60)) |
					  (mypy.radec_sep(cat.ra, cat.dec, ra_center2, dec_center2) < (radius2*60)))


	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[201.25003,30.730522,47.512735],
						 [201.22212,30.748129,13.458029],
						 [201.27612,30.772911,22.274966],
						 [201.10165,30.812886,23.675712],
						 [201.25291,30.909185,17.122665],
						 [201.15439,30.900116,12.10503],
						 [201.25132,30.829188,9.7452902],
						 [201.13392,30.846507,9.2278058],
						 [201.13487,30.824375,7.6672608],
						 [201.16874,30.769457,6.0614826],
						 [201.26223,30.47085,36.090217],
						 [201.25449,30.462586,12.670795],
						 [201.22354,30.383661,33.454697],
						 [201.20035,30.387129,27.711625],
						 [201.10457,30.344156,28.993961],
						 [201.13763,30.345778,11.809207],
						 [201.30168,30.31483,20.60916],
						 [201.2918,30.32177,8.7791076]])


	area_arcsec2 = pylab.pi * ((radius1*60)**2 + (radius2*60)**2)
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2

def selectionFunction_rxj1716(cat):

	image = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/RX_1716+6708/catalogs/detection_Riz.fits'
	imwcs = wcs.WCS(image)

	###  center and radius(arcmin) of circle for selection region
	ra_center, dec_center, radius = 259.22971, 67.143098, 7.5

	inds = pylab.find(mypy.radec_sep(cat.ra, cat.dec, ra_center, dec_center) < (radius*60))


	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[259.27315,67.137511,25.866736],
						 [259.30226,67.138466,15.023128],
						 [259.09091,67.166876,13.157666],
						 [259.10066,67.218895,13.300506],
						 [259.43251,67.196812,12.558967],
						 [259.30401,67.162694,11.047812],
						 [259.42721,67.163867,10.458859],
						 [259.35307,67.153612,8.4644722],
						 [259.18329,67.193086,8.7346083],
						 [259.15519,67.086793,10.356007],
						 [259.09615,67.113576,11.004811],
						 [259.05713,67.110635,10.120949],
						 [259.03645,67.098984,10.50559],
						 [259.33435,67.234392,10.655293],
						 [259.29348,67.246694,10.491968],
						 [259.31441,67.25702,10.165386],
						 [259.14271,67.101326,5.997919],
						 [259.13191,67.100998,7.3613226],
						 [259.13364,67.081616,9.624462],
						 [259.10055,67.243124,12.539747],
						 [259.04234,67.207214,9.2868871],
						 [259.15655,67.224092,8.8955654]])

	area_arcsec2 = pylab.pi * (radius*60)**2
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2

def selectionFunction_nep5281(cat):

	image = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/NEP5281/catalogs/detection_iz.fits'
	imwcs = wcs.WCS(image)

	###  center and radius(arcmin) of circle for selection region
	ra_center, dec_center, radius = 275.36968, 68.44678, 7.5

	inds = pylab.find(mypy.radec_sep(cat.ra, cat.dec, ra_center, dec_center) < (radius*60))


	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[275.27978,68.520525,46.056259],
						 [275.53037,68.517381,40.509531],
						 [275.50134,68.477453,24.628801],
						 [275.26677,68.477929,23.442467],
						 [275.268,68.46827,16.965804],
						 [275.59468,68.47603,26.884972],
						 [275.60796,68.486988,28.267962],
						 [275.42581,68.432277,25.402459],
						 [275.50188,68.374711,30.164172],
						 [275.47205,68.365075,27.880421],
						 [275.42809,68.39759,25.101274],
						 [275.63699,68.413618,36.739731],
						 [275.30528,68.383986,17.026886],
						 [275.27421,68.405926,17.166376],
						 [275.49268,68.43092,18.095503],
						 [275.38766,68.489806,16.166629],
						 [275.23104,68.445854,13.20395],
						 [275.32431,68.409895,15.990675],
						 [275.33392,68.352819,21.671266]])

	area_arcsec2 = pylab.pi * (radius*60)**2
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2

def selectionFunction_nep200(cat):

	image = '/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/NEP200/catalogs/detection_ri.fits'
	imwcs = wcs.WCS(image)

	###  center and radius(arcmin) of circle for selection region
	ra_center, dec_center, radius = 269.34685, 66.517873, 7.5

	inds = pylab.find(mypy.radec_sep(cat.ra, cat.dec, ra_center, dec_center) < (radius*60))


	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[269.46246,66.429833,70.075997],
						 [269.43695,66.446738,40.164039],
						 [269.46396,66.517388,42.679144],
						 [269.15691,66.559085,29.881235],
						 [269.13197,66.599941,20.804753],
						 [269.42969,66.613853,24.66493],
						 [269.26007,66.457846,15.792927],
						 [269.25119,66.452952,15.792927],
						 [269.32334,66.592537,13.907813],
						 [269.53687,66.57243,17.460344],
						 [269.29668,66.505863,16.241715],
						 [269.26984,66.532962,14.637499],
						 [269.22505,66.561377,11.433118],
						 [269.55164,66.462634,14.754796]])

	area_arcsec2 = pylab.pi * (radius*60)**2
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2

def selectionFunction_cl0910(cat):

	image='/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/Cl_0910+5422/images/sci/wregister/detection_image/my_detection_RcI+Z+.fits'
	imwcs = wcs.WCS(image)

	###  vertices of lines foriming polygon selection region
	ra1_a, dec1_a, ra2_a, dec2_a = 137.46145,54.333665,137.71293,54.468663
	ra1_b, dec1_b, ra2_b, dec2_b = 137.89374,54.400527,137.71354,54.469018
	ra1_c, dec1_c, ra2_c, dec2_c = 137.89313,54.400172,137.6829,54.293655
	ra1_d, dec1_d, ra2_d, dec2_d = 137.35995,54.202896,137.6829,54.293655
	ra1_e, dec1_e, ra2_e, dec2_e = 137.35995,54.202896,137.33306,54.292819
	ra1_f, dec1_f, ra2_f, dec2_f = 137.46145,54.333665,137.33306,54.292819

	(x1_a, y1_a), (x2_a, y2_a) = imwcs.wcs_world2pix(pylab.array(zip([ra1_a, ra2_a], [dec1_a, dec2_a])), 1)
	(x1_b, y1_b), (x2_b, y2_b) = imwcs.wcs_world2pix(pylab.array(zip([ra1_b, ra2_b], [dec1_b, dec2_b])), 1)
	(x1_c, y1_c), (x2_c, y2_c) = imwcs.wcs_world2pix(pylab.array(zip([ra1_c, ra2_c], [dec1_c, dec2_c])), 1)
	(x1_d, y1_d), (x2_d, y2_d) = imwcs.wcs_world2pix(pylab.array(zip([ra1_d, ra2_d], [dec1_d, dec2_d])), 1)
	(x1_e, y1_e), (x2_e, y2_e) = imwcs.wcs_world2pix(pylab.array(zip([ra1_e, ra2_e], [dec1_e, dec2_e])), 1)
	(x1_f, y1_f), (x2_f, y2_f) = imwcs.wcs_world2pix(pylab.array(zip([ra1_f, ra2_f], [dec1_f, dec2_f])), 1)
	xs = pylab.array([x1_a, x2_a, x1_b, x2_b, x1_c, x2_c, x1_d, x2_d, x1_e, x2_e, x1_f, x2_f])
	ys = pylab.array([y1_a, y2_a, y1_b, y2_b, y1_c, y2_c, y1_d, y2_d, y1_e, y2_e, y1_f, y2_f])

	slope_a = (y1_a - y2_a) / (x1_a - x2_a)
	slope_b = (y1_b - y2_b) / (x1_b - x2_b)
	slope_c = (y1_c - y2_c) / (x1_c - x2_c)
	slope_d = (y1_d - y2_d) / (x1_d - x2_d)
	slope_e = (y1_e - y2_e) / (x1_e - x2_e)
	slope_f = (y1_f - y2_f) / (x1_f - x2_f)

	norm_a = y1_a - slope_a * x1_a
	norm_b = y1_b - slope_b * x1_b
	norm_c = y1_c - slope_c * x1_c
	norm_d = y1_d - slope_d * x1_d
	norm_e = y1_e - slope_e * x1_e
	norm_f = y1_f - slope_f * x1_f

	###  selecting galaxies within a user-defind polygon enclosing the spectroscopic coverage
	inds = pylab.find(((cat.y < slope_a * cat.x + norm_a) | (cat.y < slope_f * cat.x + norm_f)) & 
					  (cat.y < slope_b * cat.x + norm_b) & 
					  (cat.y > slope_c * cat.x + norm_c) & 
					  (cat.y > slope_d * cat.x + norm_d) & 
					  (cat.y > slope_e * cat.x + norm_e) & 
					  (cat.use == 1))


	###  calculating the area of the user-defind polygon
	ra12, dec12 = imwcs.wcs_pix2world([0, 0], [0, 1], 1)
	px_scale = mypy.radec_sep(ra12[0], dec12[0], ra12[1], dec12[1])

	###  ra_center, dec_center, radius(arcsec)
	holes = pylab.array([[137.72017,54.439911,11.491746],
						 [137.81728,54.398832,6.7431164],
						 [137.71937,54.409576,9.3121025],
						 [137.76736,54.369088,12.255985],
						 [137.7616,54.351868,9.1073549],
						 [137.68401,54.371568,6.48],
						 [137.67523,54.374569,11.123555],
						 [137.65235,54.358788,9.5279771],
						 [137.64415,54.358564,9.6752827],
						 [137.62866,54.394229,19.453046],
						 [137.65501,54.374567,5.6216075],
						 [137.58682,54.325645,13.965681],
						 [137.5307,54.353927,8.7263942],
						 [137.58229,54.305084,8.101748],
						 [137.5724,54.299409,16.05717],
						 [137.61487,54.295104,6.9433894],
						 [137.59756,54.28876,6.0919288],
						 [137.44301,54.249783,10.550997],
						 [137.46675,54.292829,7.7258551],
						 [137.54163,54.2766,7.7839873],
						 [137.48524,54.330752,7.1368481],
						 [137.52003,54.250906,8.5717222],
						 [137.50938,54.250447,7.4277835],
						 [137.47348,54.241281,11.223961],
						 [137.41227,54.237608,12.044997],
						 [137.34195,54.281106,11.348039],
						 [137.86589,54.391751,6.305658]])

	area_px2 = PolyArea(xs, ys)
	area_arcsec2 = area_px2 * px_scale**2
	for hole in holes:
		area_arcsec2 -= pylab.pi * hole[2]**2
	area_arcmin2 = area_arcsec2 / 60.**2

	return inds, area_arcmin2






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
		#self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (data_dir, version, version), readcat=1)

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



















dzphot_thresh = 0.05
zmag_thresh = 23.5















###   1 sigma_zphot
###   (zphot-zlss)/(1+zphot)  vs  (zspec-zlss)/(1+zspec)

fig = pylab.figure(figsize=(13., 13.))

sp2 = fig.add_subplot(3, 3, 2)
sp1 = fig.add_subplot(3, 3, 1)

sp5 = fig.add_subplot(3, 3, 6)
sp4 = fig.add_subplot(3, 3, 5)
sp3 = fig.add_subplot(3, 3, 4)

sp8 = fig.add_subplot(3, 3, 9)
sp7 = fig.add_subplot(3, 3, 8)
sp6 = fig.add_subplot(3, 3, 7)

sps = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8]

fig.subplots_adjust(wspace=0, hspace=0)

rproj_search = 1.   # projected distance from substructure centers to search for galaxies
n_falsepos, n_truepos, n_falseneg, n_trueneg = [], [], [], []

for i in range(len(fields)):

	f = fields[i]

	dzphot_thresh = 1. * f.sigmaz

	zinds = pylab.find((f.cat.z_spec > 0) & (f.zmagnitude <= zmag_thresh))

	ra_spec = f.cat.ra[zinds]
	dec_spec = f.cat.dec[zinds]

	###  Identifying which LSS substructure a galaxy belongs to.
	###  This is simply defined as a galaxy needs to be <=0.5Mpc
	###  of any substructure center. If there are two such
	###  substructures, then the nearest is chosen.

	zspec = []
	zphot = []
	zlss = []
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
		if min(rproj_list) >= 1000: continue


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


		### now that we've grabbed the galaxy, determine if it is:
		###   0 = true +
		###   1 = true -
		###   2 = false +
		###   3 = false -
		zspec.append(zspec_i)
		zphot.append(zphot_i)
		zlss.append(zlss_i)
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






	zspec = pylab.array(zspec)
	zphot = pylab.array(zphot)
	zlss = pylab.array(zlss)

	dzspec = (zspec - zlss) / (1 + zspec)
	dzphot = (zphot - zlss) / (1 + zphot)



	sp = sps[i]

	sp.minorticks_on()
	sp.axis([-0.15, 0.15, -0.15, 0.15])
	sp.xaxis.set_ticklabels(['', '-0.1', '', '0.0', '', '0.1', ''])
	sp.yaxis.set_ticklabels(['', '-0.10', '-0.05', '0.00', '0.05', '0.10', ''])
	sp.set_xlabel('$(z_{spec} - z_{LSS}) / (1 + z_{spec})$')
	sp.set_ylabel('$(z_{phot} - z_{LSS}) / (1 + z_{phot})$')


	inds_truepos = pylab.find(truefalse_array == 0)
	inds_trueneg = pylab.find(truefalse_array == 1)
	inds_falsepos = pylab.find(truefalse_array == 2)
	inds_falseneg = pylab.find(truefalse_array == 3)


	sp.plot(dzspec[inds_truepos], dzphot[inds_truepos], ls='', marker='o', ms=3, mec='k', mfc='k', mew=1)
	sp.plot(dzspec[inds_trueneg], dzphot[inds_trueneg], ls='', marker='^', ms=4, mec='gray', mfc='gray', mew=1)
	sp.plot(dzspec[inds_falsepos], dzphot[inds_falsepos], ls='', marker='x', ms=6, mec='r', mfc='r', mew=1)
	sp.plot(dzspec[inds_falseneg], dzphot[inds_falseneg], ls='', marker='v', ms=4, mec='orange', mfc='orange', mew=1)

	sp.axhline(dzphot_thresh, ls='--', color='gray')
	sp.axhline(-dzphot_thresh, ls='--', color='gray')
	sp.axhline(0, color='#666666', lw=1., ls=':')




	n_truepos.append(len(inds_truepos))
	n_trueneg.append(len(inds_trueneg))
	n_falsepos.append(len(inds_falsepos))
	n_falseneg.append(len(inds_falseneg))

	frac_falsepos = n_falsepos[-1] * 1. / (n_truepos[-1] + n_falseneg[-1])
	frac_falseneg = n_falseneg[-1] * 1. / (n_truepos[-1] + n_falseneg[-1])




	label_false_stats = 'false + : %.2f\nfalse  - : %.2f' % (frac_falsepos, frac_falseneg)
	label_false_stats = 'false + : %2i%%\nfalse  - : %2i%%' % (frac_falsepos*100, frac_falseneg*100)
	t = sp.text(0.95, 0.05, label_false_stats, transform=sp.transAxes, fontsize=14, color='k', horizontalalignment='right',
				path_effects=[PathEffects.withStroke(linewidth=1., foreground='r')])



	label = '%s\nN$_{tot}$ = %i' % (f.name, len(zlss))
	t = sp.text(0.05, 0.85, label, transform=sp.transAxes, fontsize=14, color='k',
				path_effects=[PathEffects.withStroke(linewidth=1., foreground='#8080ff')])



	if i == 1:
		sp.plot(10, 10, ls='', marker='o', ms=7, mec='k', mfc='k', mew=1, label='true positive')
		sp.plot(10, 10, ls='', marker='^', ms=7, mec='gray', mfc='gray', mew=1, label='true negative')
		sp.plot(10, 10, ls='', marker='x', ms=7, mec='r', mfc='r', mew=1, label='false positive')
		sp.plot(10, 10, ls='', marker='v', ms=7, mec='orange', mfc='orange', mew=1, label='false negative')
		leg = sp.legend(loc=(1.12, 0.2), fontsize=16, numpoints=1)
		leg.set_title('z$_{mag}$ < %.1f\n$\Delta z = \pm 1\sigma_{zphot}$\nR$_{proj}$ < %.2f Mpc' % (zmag_thresh, rproj_search))


fig.savefig('../output/dzphot_dzspec_zlss_1.0sigma_zphot.png')
pylab.close()
print 'done with +/- 1 sigma_zphot'


###  frac_falsepos vs frac_falseneg (R<1Mpc)
fig = pylab.figure(figsize=(7.6, 7.4))
sp = fig.add_subplot(111, aspect=1)

sp.minorticks_on()
sp.set_xlabel('false positive fraction [%]')
sp.set_ylabel('false negative fraction [%]')

sp.grid()
sp.axis([0, 100, 0, 100])
sp.plot([0, 100], [0, 100], color='#cccccc', lw=1)

sym = ['o', 's', 'D', 'p', '^', '<', '>', 'v']
col = ['b', 'g', 'r', 'orange', 'm', 'gray', '#6600ff', '#33cccc']

for i in range(len(fields)):

	n_truepos_i = n_truepos[i]
	n_trueneg_i = n_trueneg[i]
	n_falsepos_i = n_falsepos[i]
	n_falseneg_i = n_falseneg[i]

	nhi_falsepos_i, nlo_falsepos_i = mypy.massfunc.confidence_interval(n_falsepos_i)
	nhi_falseneg_i, nlo_falseneg_i = mypy.massfunc.confidence_interval(n_falseneg_i)

	f_falsepos = n_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	f_falseneg = n_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)

	flo_falsepos = nlo_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	flo_falseneg = nlo_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)
	fhi_falsepos = nhi_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	fhi_falseneg = nhi_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)

	sp.errorbar(f_falsepos, f_falseneg, xerr=[[flo_falsepos], [fhi_falsepos]], yerr=[[flo_falseneg], [fhi_falseneg]], \
				ls='', marker=sym[i], mfc='w', mec=col[i], ecolor=col[i], elinewidth=1.5, capsize=0, ms=9, mew=2, label=fields[i].name)

n_truepos_tot = pylab.sum(n_truepos)
n_trueneg_tot = pylab.sum(n_trueneg)
n_falsepos_tot = pylab.sum(n_falsepos)
n_falseneg_tot = pylab.sum(n_falseneg)

nhi_falsepos_tot, nlo_falsepos_tot = mypy.massfunc.confidence_interval(n_falsepos_tot)
nhi_falseneg_tot, nlo_falseneg_tot = mypy.massfunc.confidence_interval(n_falseneg_tot)

f_falsepos = n_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
f_falseneg = n_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)

flo_falsepos = nlo_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
flo_falseneg = nlo_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)
fhi_falsepos = nhi_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
fhi_falseneg = nhi_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)

sp.errorbar(f_falsepos, f_falseneg, xerr=[[flo_falsepos], [fhi_falsepos]], yerr=[[flo_falseneg], [fhi_falseneg]], \
			ls='', marker='o', mfc='k', mec='k', ecolor='k', elinewidth=1.5, capsize=0, ms=9, mew=2, label='combined')

leg = sp.legend(title='z$_{mag}$ < %.1f\nR$_{proj}$ = 1 Mpc\n$\Delta z = \pm 1\sigma_{zphot}$' % (zmag_thresh), \
				loc=1, numpoints=1, fontsize=15, framealpha=0., shadow=0)
fig.savefig('../output/falsePosNeg_fractions_1.0sigma_zphot.png')
pylab.close()
print 'done with false +/-\n'






















###   2 sigma_zphot
###   (zphot-zlss)/(1+zphot)  vs  (zspec-zlss)/(1+zspec)

fig = pylab.figure(figsize=(13., 13.))

sp2 = fig.add_subplot(3, 3, 2)
sp1 = fig.add_subplot(3, 3, 1)

sp5 = fig.add_subplot(3, 3, 6)
sp4 = fig.add_subplot(3, 3, 5)
sp3 = fig.add_subplot(3, 3, 4)

sp8 = fig.add_subplot(3, 3, 9)
sp7 = fig.add_subplot(3, 3, 8)
sp6 = fig.add_subplot(3, 3, 7)

sps = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8]

fig.subplots_adjust(wspace=0, hspace=0)

rproj_search = 1.   # projected distance from substructure centers to search for galaxies
n_falsepos, n_truepos, n_falseneg, n_trueneg = [], [], [], []

for i in range(len(fields)):

	f = fields[i]

	dzphot_thresh = 2 * f.sigmaz

	zinds = pylab.find((f.cat.z_spec > 0) & (f.zmagnitude <= zmag_thresh))

	ra_spec = f.cat.ra[zinds]
	dec_spec = f.cat.dec[zinds]

	###  Identifying which LSS substructure a galaxy belongs to.
	###  This is simply defined as a galaxy needs to be <=0.5Mpc
	###  of any substructure center. If there are two such
	###  substructures, then the nearest is chosen.

	zspec = []
	zphot = []
	zlss = []
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
		if min(rproj_list) >= 1000: continue


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


		### now that we've grabbed the galaxy, determine if it is:
		###   0 = true +
		###   1 = true -
		###   2 = false +
		###   3 = false -
		zspec.append(zspec_i)
		zphot.append(zphot_i)
		zlss.append(zlss_i)
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





	zspec = pylab.array(zspec)
	zphot = pylab.array(zphot)
	zlss = pylab.array(zlss)

	dzspec = (zspec - zlss) / (1 + zspec)
	dzphot = (zphot - zlss) / (1 + zphot)



	sp = sps[i]

	sp.minorticks_on()
	sp.axis([-0.15, 0.15, -0.15, 0.15])
	sp.xaxis.set_ticklabels(['', '-0.1', '', '0.0', '', '0.1', ''])
	sp.yaxis.set_ticklabels(['', '-0.10', '-0.05', '0.00', '0.05', '0.10', ''])
	sp.set_xlabel('$(z_{spec} - z_{LSS}) / (1 + z_{spec})$')
	sp.set_ylabel('$(z_{phot} - z_{LSS}) / (1 + z_{phot})$')


	inds_truepos = pylab.find(truefalse_array == 0)
	inds_trueneg = pylab.find(truefalse_array == 1)
	inds_falsepos = pylab.find(truefalse_array == 2)
	inds_falseneg = pylab.find(truefalse_array == 3)


	sp.plot(dzspec[inds_truepos], dzphot[inds_truepos], ls='', marker='o', ms=3, mec='k', mfc='k', mew=1)
	sp.plot(dzspec[inds_trueneg], dzphot[inds_trueneg], ls='', marker='^', ms=4, mec='gray', mfc='gray', mew=1)
	sp.plot(dzspec[inds_falsepos], dzphot[inds_falsepos], ls='', marker='x', ms=6, mec='r', mfc='r', mew=1)
	sp.plot(dzspec[inds_falseneg], dzphot[inds_falseneg], ls='', marker='v', ms=4, mec='orange', mfc='orange', mew=1)

	sp.axhline(dzphot_thresh, ls='--', color='gray')
	sp.axhline(-dzphot_thresh, ls='--', color='gray')
	sp.axhline(0, color='#666666', lw=1., ls=':')




	n_truepos.append(len(inds_truepos))
	n_trueneg.append(len(inds_trueneg))
	n_falsepos.append(len(inds_falsepos))
	n_falseneg.append(len(inds_falseneg))

	frac_falsepos = n_falsepos[-1] * 1. / (n_truepos[-1] + n_falseneg[-1])
	frac_falseneg = n_falseneg[-1] * 1. / (n_truepos[-1] + n_falseneg[-1])




	label_false_stats = 'false + : %.2f\nfalse  - : %.2f' % (frac_falsepos, frac_falseneg)
	label_false_stats = 'false + : %2i%%\nfalse  - : %2i%%' % (frac_falsepos*100, frac_falseneg*100)
	t = sp.text(0.95, 0.05, label_false_stats, transform=sp.transAxes, fontsize=14, color='k', horizontalalignment='right',
				path_effects=[PathEffects.withStroke(linewidth=1., foreground='r')])



	label = '%s\nN$_{tot}$ = %i' % (f.name, len(zlss))
	t = sp.text(0.05, 0.85, label, transform=sp.transAxes, fontsize=14, color='k',
				path_effects=[PathEffects.withStroke(linewidth=1., foreground='#8080ff')])



	if i == 1:
		sp.plot(10, 10, ls='', marker='o', ms=7, mec='k', mfc='k', mew=1, label='true positive')
		sp.plot(10, 10, ls='', marker='^', ms=7, mec='gray', mfc='gray', mew=1, label='true negative')
		sp.plot(10, 10, ls='', marker='x', ms=7, mec='r', mfc='r', mew=1, label='false positive')
		sp.plot(10, 10, ls='', marker='v', ms=7, mec='orange', mfc='orange', mew=1, label='false negative')
		leg = sp.legend(loc=(1.12, 0.2), fontsize=16, numpoints=1)
		leg.set_title('z$_{mag}$ < %.1f\n$\Delta z = \pm 2\sigma_{zphot}$\nR$_{proj}$ < %.2f Mpc' % (zmag_thresh, rproj_search))


fig.savefig('../output/dzphot_dzspec_zlss_2.0sigma_zphot.png')
pylab.close()
print 'done with +/- 2 sigma_zphot'



###  frac_falsepos vs frac_falseneg (R<1Mpc)
fig = pylab.figure(figsize=(7.6, 7.4))
sp = fig.add_subplot(111, aspect=1)

sp.minorticks_on()
sp.set_xlabel('false positive fraction [%]')
sp.set_ylabel('false negative fraction [%]')

sp.grid()
sp.axis([0, 100, 0, 100])
sp.plot([0, 100], [0, 100], color='#cccccc', lw=1)

sym = ['o', 's', 'D', 'p', '^', '<', '>', 'v']
col = ['b', 'g', 'r', 'orange', 'm', 'gray', '#6600ff', '#33cccc']

for i in range(len(fields)):

	n_truepos_i = n_truepos[i]
	n_trueneg_i = n_trueneg[i]
	n_falsepos_i = n_falsepos[i]
	n_falseneg_i = n_falseneg[i]

	nhi_falsepos_i, nlo_falsepos_i = mypy.massfunc.confidence_interval(n_falsepos_i)
	nhi_falseneg_i, nlo_falseneg_i = mypy.massfunc.confidence_interval(n_falseneg_i)

	f_falsepos = n_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	f_falseneg = n_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)

	flo_falsepos = nlo_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	flo_falseneg = nlo_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)
	fhi_falsepos = nhi_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	fhi_falseneg = nhi_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)

	sp.errorbar(f_falsepos, f_falseneg, xerr=[[flo_falsepos], [fhi_falsepos]], yerr=[[flo_falseneg], [fhi_falseneg]], \
				ls='', marker=sym[i], mfc='w', mec=col[i], ecolor=col[i], elinewidth=1.5, capsize=0, ms=9, mew=2, label=fields[i].name)

n_truepos_tot = pylab.sum(n_truepos)
n_trueneg_tot = pylab.sum(n_trueneg)
n_falsepos_tot = pylab.sum(n_falsepos)
n_falseneg_tot = pylab.sum(n_falseneg)

nhi_falsepos_tot, nlo_falsepos_tot = mypy.massfunc.confidence_interval(n_falsepos_tot)
nhi_falseneg_tot, nlo_falseneg_tot = mypy.massfunc.confidence_interval(n_falseneg_tot)

f_falsepos = n_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
f_falseneg = n_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)

flo_falsepos = nlo_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
flo_falseneg = nlo_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)
fhi_falsepos = nhi_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
fhi_falseneg = nhi_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)

sp.errorbar(f_falsepos, f_falseneg, xerr=[[flo_falsepos], [fhi_falsepos]], yerr=[[flo_falseneg], [fhi_falseneg]], \
			ls='', marker='o', mfc='k', mec='k', ecolor='k', elinewidth=1.5, capsize=0, ms=9, mew=2, label='combined')


leg = sp.legend(title='z$_{mag}$ < %.1f\nR$_{proj}$ = 1 Mpc\n$\Delta z = \pm 2\sigma_{zphot}$' % (zmag_thresh), \
				loc=1, numpoints=1, fontsize=15, framealpha=0., shadow=0)
fig.savefig('../output/falsePosNeg_fractions_2.0sigma_zphot.png')
pylab.close()
print 'done with false +/-\n'






















###   1.5 sigma_zphot
###   (zphot-zlss)/(1+zphot)  vs  (zspec-zlss)/(1+zspec)

fig = pylab.figure(figsize=(13., 13.))

sp2 = fig.add_subplot(3, 3, 2)
sp1 = fig.add_subplot(3, 3, 1)

sp5 = fig.add_subplot(3, 3, 6)
sp4 = fig.add_subplot(3, 3, 5)
sp3 = fig.add_subplot(3, 3, 4)

sp8 = fig.add_subplot(3, 3, 9)
sp7 = fig.add_subplot(3, 3, 8)
sp6 = fig.add_subplot(3, 3, 7)

sps = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8]

fig.subplots_adjust(wspace=0, hspace=0)

rproj_search = 1.   # projected distance from substructure centers to search for galaxies
n_falsepos, n_truepos, n_falseneg, n_trueneg = [], [], [], []

for i in range(len(fields)):

	f = fields[i]

	dzphot_thresh = 1.5 * f.sigmaz

	zinds = pylab.find((f.cat.z_spec > 0) & (f.zmagnitude <= zmag_thresh))

	ra_spec = f.cat.ra[zinds]
	dec_spec = f.cat.dec[zinds]

	###  Identifying which LSS substructure a galaxy belongs to.
	###  This is simply defined as a galaxy needs to be <=0.5Mpc
	###  of any substructure center. If there are two such
	###  substructures, then the nearest is chosen.

	zspec = []
	zphot = []
	zlss = []
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
		if min(rproj_list) >= 1000: continue


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


		### now that we've grabbed the galaxy, determine if it is:
		###   0 = true +
		###   1 = true -
		###   2 = false +
		###   3 = false -
		zspec.append(zspec_i)
		zphot.append(zphot_i)
		zlss.append(zlss_i)
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




	zspec = pylab.array(zspec)
	zphot = pylab.array(zphot)
	zlss = pylab.array(zlss)

	dzspec = (zspec - zlss) / (1 + zspec)
	dzphot = (zphot - zlss) / (1 + zphot)



	sp = sps[i]

	sp.minorticks_on()
	sp.axis([-0.15, 0.15, -0.15, 0.15])
	sp.xaxis.set_ticklabels(['', '-0.1', '', '0.0', '', '0.1', ''])
	sp.yaxis.set_ticklabels(['', '-0.10', '-0.05', '0.00', '0.05', '0.10', ''])
	sp.set_xlabel('$(z_{spec} - z_{LSS}) / (1 + z_{spec})$')
	sp.set_ylabel('$(z_{phot} - z_{LSS}) / (1 + z_{phot})$')


	inds_truepos = pylab.find(truefalse_array == 0)
	inds_trueneg = pylab.find(truefalse_array == 1)
	inds_falsepos = pylab.find(truefalse_array == 2)
	inds_falseneg = pylab.find(truefalse_array == 3)


	sp.plot(dzspec[inds_truepos], dzphot[inds_truepos], ls='', marker='o', ms=3, mec='k', mfc='k', mew=1)
	sp.plot(dzspec[inds_trueneg], dzphot[inds_trueneg], ls='', marker='^', ms=4, mec='gray', mfc='gray', mew=1)
	sp.plot(dzspec[inds_falsepos], dzphot[inds_falsepos], ls='', marker='x', ms=6, mec='r', mfc='r', mew=1)
	sp.plot(dzspec[inds_falseneg], dzphot[inds_falseneg], ls='', marker='v', ms=4, mec='orange', mfc='orange', mew=1)

	sp.axhline(dzphot_thresh, ls='--', color='gray')
	sp.axhline(-dzphot_thresh, ls='--', color='gray')
	sp.axhline(0, color='#666666', lw=1., ls=':')




	n_truepos.append(len(inds_truepos))
	n_trueneg.append(len(inds_trueneg))
	n_falsepos.append(len(inds_falsepos))
	n_falseneg.append(len(inds_falseneg))

	frac_falsepos = n_falsepos[-1] * 1. / (n_truepos[-1] + n_falseneg[-1])
	frac_falseneg = n_falseneg[-1] * 1. / (n_truepos[-1] + n_falseneg[-1])




	label_false_stats = 'false + : %.2f\nfalse  - : %.2f' % (frac_falsepos, frac_falseneg)
	label_false_stats = 'false + : %2i%%\nfalse  - : %2i%%' % (frac_falsepos*100, frac_falseneg*100)
	t = sp.text(0.95, 0.05, label_false_stats, transform=sp.transAxes, fontsize=14, color='k', horizontalalignment='right',
				path_effects=[PathEffects.withStroke(linewidth=1., foreground='r')])



	label = '%s\nN$_{tot}$ = %i' % (f.name, len(zlss))
	t = sp.text(0.05, 0.85, label, transform=sp.transAxes, fontsize=14, color='k',
				path_effects=[PathEffects.withStroke(linewidth=1., foreground='#8080ff')])



	if i == 1:
		sp.plot(10, 10, ls='', marker='o', ms=7, mec='k', mfc='k', mew=1, label='true positive')
		sp.plot(10, 10, ls='', marker='^', ms=7, mec='gray', mfc='gray', mew=1, label='true negative')
		sp.plot(10, 10, ls='', marker='x', ms=7, mec='r', mfc='r', mew=1, label='false positive')
		sp.plot(10, 10, ls='', marker='v', ms=7, mec='orange', mfc='orange', mew=1, label='false negative')
		leg = sp.legend(loc=(1.12, 0.2), fontsize=16, numpoints=1)
		leg.set_title('z$_{mag}$ < %.1f\n$\Delta z = \pm 1.5\sigma_{zphot}$\nR$_{proj}$ < %.2f Mpc' % (zmag_thresh, rproj_search))


fig.savefig('../output/dzphot_dzspec_zlss_1.5sigma_zphot.png')
pylab.close()
print 'done with +/- 1.5 sigma_zphot'



###  frac_falsepos vs frac_falseneg (R<1Mpc)
fig = pylab.figure(figsize=(7.6, 7.4))
sp = fig.add_subplot(111, aspect=1)

sp.minorticks_on()
sp.set_xlabel('false positive fraction [%]')
sp.set_ylabel('false negative fraction [%]')

sp.grid()
sp.axis([0, 100, 0, 100])
sp.plot([0, 100], [0, 100], color='#cccccc', lw=1)

sym = ['o', 's', 'D', 'p', '^', '<', '>', 'v']
col = ['b', 'g', 'r', 'orange', 'm', 'gray', '#6600ff', '#33cccc']

for i in range(len(fields)):

	n_truepos_i = n_truepos[i]
	n_trueneg_i = n_trueneg[i]
	n_falsepos_i = n_falsepos[i]
	n_falseneg_i = n_falseneg[i]

	nhi_falsepos_i, nlo_falsepos_i = mypy.massfunc.confidence_interval(n_falsepos_i)
	nhi_falseneg_i, nlo_falseneg_i = mypy.massfunc.confidence_interval(n_falseneg_i)

	f_falsepos = n_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	f_falseneg = n_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)

	flo_falsepos = nlo_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	flo_falseneg = nlo_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)
	fhi_falsepos = nhi_falsepos_i * 100. / (n_truepos_i + n_falseneg_i)
	fhi_falseneg = nhi_falseneg_i * 100. / (n_truepos_i + n_falseneg_i)

	sp.errorbar(f_falsepos, f_falseneg, xerr=[[flo_falsepos], [fhi_falsepos]], yerr=[[flo_falseneg], [fhi_falseneg]], \
				ls='', marker=sym[i], mfc='w', mec=col[i], ecolor=col[i], elinewidth=1.5, capsize=0, ms=9, mew=2, label=fields[i].name)

n_truepos_tot = pylab.sum(n_truepos)
n_trueneg_tot = pylab.sum(n_trueneg)
n_falsepos_tot = pylab.sum(n_falsepos)
n_falseneg_tot = pylab.sum(n_falseneg)

nhi_falsepos_tot, nlo_falsepos_tot = mypy.massfunc.confidence_interval(n_falsepos_tot)
nhi_falseneg_tot, nlo_falseneg_tot = mypy.massfunc.confidence_interval(n_falseneg_tot)

f_falsepos = n_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
f_falseneg = n_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)

flo_falsepos = nlo_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
flo_falseneg = nlo_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)
fhi_falsepos = nhi_falsepos_tot * 100. / (n_truepos_tot + n_falseneg_tot)
fhi_falseneg = nhi_falseneg_tot * 100. / (n_truepos_tot + n_falseneg_tot)

sp.errorbar(f_falsepos, f_falseneg, xerr=[[flo_falsepos], [fhi_falsepos]], yerr=[[flo_falseneg], [fhi_falseneg]], \
			ls='', marker='o', mfc='k', mec='k', ecolor='k', elinewidth=1.5, capsize=0, ms=9, mew=2, label='combined')


leg = sp.legend(title='z$_{mag}$ < %.1f\nR$_{proj}$ = 1 Mpc\n$\Delta z = \pm 1.5\sigma_{zphot}$' % (zmag_thresh), \
				loc=1, numpoints=1, fontsize=15, framealpha=0., shadow=0)
fig.savefig('../output/falsePosNeg_fractions_1.5sigma_zphot.png')
pylab.close()
print 'done with false +/-\n'






















tf = time.time()
print '\n%.1f mins elapsed\n' % ((tf-t0)/60.)










