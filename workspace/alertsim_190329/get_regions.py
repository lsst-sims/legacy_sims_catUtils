from lsst.sims.utils import ModifiedJulianDate
from lsst.sims.utils import htmModule as htm

from lsst.sims.utils import galacticFromEquatorial
from lsst.sims.utils import equatorialFromGalactic

import palpy as palpy

import numpy as np

rng = np.random.RandomState(11723)

tx_dict = htm.getAllTrixels(6)
tx_list = []
for htmid in tx_dict:
    if htm.levelFromHtmid(htmid) == 6:
        tx_list.append(tx_dict[htmid])

######### Region 1
#1) ~600 deg2 along ecliptic plane, say |latitude| < 30 and a 100 deg stretch
#    in longitude, but with |galactic latitude| > 30 deg

gal_n_ra, gal_n_dec = equatorialFromGalactic(0.0, 90.0)
gal_s_ra, gal_s_dec = equatorialFromGalactic(0.0, -90.0)

# must be fully inside one of these two half spaces to have
# |galactic latitude| > 30
gal_n_30_hs = htm.halfSpaceFromRaDec(gal_n_ra, gal_n_dec, 60.0)
gal_s_30_hs = htm.halfSpaceFromRaDec(gal_s_ra, gal_s_dec, 60.0)

ecl_n_ra, ecl_n_dec = np.degrees(palpy.ecleq(0.0, 0.5*np.pi, 59580.0))
ecl_s_ra, ecl_s_dec = np.degrees(palpy.ecleq(0.0, -0.5*np.pi, 59580.0))

# must be inside both of these half spaces to have |ecliptic latitude|<30
ecl_n_30_hs = htm.halfSpaceFromRaDec(ecl_n_ra, ecl_n_dec, 120.0)
ecl_s_30_hs = htm.halfSpaceFromRaDec(ecl_s_ra, ecl_s_dec, 120.0)

ecl_e_ra, ecl_e_dec = np.degrees(palpy.ecleq(-0.2*np.pi, 0.0, 59580.0))

gal_lon, gal_lat = galacticFromEquatorial(ecl_e_ra, ecl_e_dec)
print('gal_lon %e gal_lat %e' % (gal_lon, gal_lat))
print('ra %e dec %e' % (ecl_e_ra, ecl_e_dec))

# half space bounding the 10 degree stretch of longitude
# must be inside this half space
ecl_e_hs = htm.halfSpaceFromRaDec(ecl_e_ra, ecl_e_dec, 12.0)

region_1_trixels = []
for tx in tx_list:
    lat_30_cut = False
    if gal_n_30_hs.contains_trixel(tx) == 'full':
        lat_30_cut = True
    if gal_s_30_hs.contains_trixel(tx) == 'full':
        lat_30_cut = True
    if not lat_30_cut:
        continue

    if (ecl_n_30_hs.contains_trixel(tx) == 'outside' or
        ecl_s_30_hs.contains_trixel(tx) == 'outside'):
        continue

    if ecl_e_hs.contains_trixel(tx) == 'outside':
        continue

    region_1_trixels.append(tx)

print('n_trixels %d' % len(region_1_trixels))

with open('data/region_1_trixels.txt', 'w') as out_file:
    out_file.write('# htmid ra dec ecliptic_lon ecliptic_lat glon glat\n')
    for tx in region_1_trixels:
        ra, dec = tx.get_center()
        lon, lat = np.degrees(palpy.eqecl(np.radians(ra), np.radians(dec),
                                          59580.0))

        glon, glat = galacticFromEquatorial(ra,dec)
        out_file.write('%d %e %e %e %e %e %e\n' %
                       (tx.htmid, ra, dec, lon, lat, glon,glat))


exit()
# test area
n_pts = 3000000
pts = rng.normal(0.0, 1.0, (3, n_pts))
norms = np.sqrt(np.sum(pts**2, axis=0))
assert len(norms) == n_pts
pts /= norms
pts = pts.transpose()

norms = np.sqrt(np.sum(pts**2, axis=1))
assert len(norms) == n_pts
assert np.abs(norms.min()-1.0)<1.0e-5
assert np.abs(norms.max()-1.0)<1.0e-5

contains = np.zeros(len(pts), dtype=bool)
for tx in region_1_trixels:
    outside = np.where(np.logical_not(contains))
    #print('n_out %e' % len(outside[0]))
    contains[outside] = tx.contains_pt(pts[outside])

valid  = np.where(contains)
frac = len(valid[0])/n_pts

area = frac*4.0*np.pi*(180.0/np.pi)**2
print('%e -- %e sq degrees' % (n_pts, area))
