def checkPoint(ref_ra, ref_dec, check_ra, check_dec):
    """
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
        
        """
    from scipy.spatial import ConvexHull
    from matplotlib.path import Path
    points = np.zeros((len(ref_ra), 2))
    points[:,0] = ref_ra
    points[:,1] = ref_dec
    hull = ConvexHull(points)
    hull_path = Path(points[hull.vertices])
    
    check = np.zeros(len(check_ra)) - 1
    for i in range(len(check_ra)):
        check[i] = hull_path.contains_point((check_ra[i], check_dec[i]))
    return check