import os
import numpy as np
import json
from lsst.utils import getPackageDir
from lsst.sims.utils import halfSpaceFromRaDec
from lsst.sims.utils import halfSpaceFromPoints
from lsst.sims.utils import cartesianFromSpherical, sphericalFromCartesian
from lsst.sims.catUtils.baseCatalogModels import GalaxyTileObj

__all__ = ["FatboyTiles"]


class FatboyTiles(object):
    """
    A class to store the fatboy galaxy tiles as a series of Half Spaces
    """

    def __init__(self):
        data_dir = os.path.join(getPackageDir('sims_catUtils'), 'data')
        data_file = os.path.join(data_dir, 'tile_data.txt')
        dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                          ('box', str, 500)])
        tile_data = np.genfromtxt(data_file, dtype=dtype, delimiter=';')

        self._tile_id = tile_data['id']
        self._tile_ra = {}
        self._tile_dec = {}
        for ii, rr, dd, in zip(tile_data['id'], tile_data['ra'], tile_data['dec']):
            self._tile_ra[ii] = rr
            self._tile_dec[ii] = dd

        self._tile_half_spaces = {}
        tol = 1.0e-10
        for tile_id, box in zip(tile_data['id'], tile_data['box']):
            local_hs_list = []
            box_corners = json.loads(box)
            ra_range = [c[0] for c in box_corners]
            ra_min = min(ra_range)
            ra_max = max(ra_range)
            dec_range = [c[1] for c in box_corners]
            dec_min = min(dec_range)
            dec_max = max(dec_range)
            #print(tile_id,dec_min,dec_max)
            for i_c1 in range(len(box_corners)):
                c1 = box_corners[i_c1]
                pt1 = cartesianFromSpherical(np.degrees(c1[0]), np.degrees(c1[1]))
                for i_c2 in range(i_c1+1, len(box_corners), 1):
                    hs = None
                    c2 = box_corners[i_c2]
                    pt2 = cartesianFromSpherical(np.degrees(c2[0]), np.degrees(c2[1]))
                    if np.abs(1.0-np.dot(pt1, pt2))<tol:
                        continue

                    dra = np.abs(c1[0]-c2[0])
                    ddec = np.abs(c1[1]-c2[1])

                    if dra<tol and ddec>tol:
                        # The RAs of the two corners is identical, but the Decs are
                        # different; this Half Space is defined by a Great Circle
                        if np.abs(c1[0]-ra_min)<tol:
                            inner_pt = (ra_min+0.001, dec_min+0.001)
                        else:
                            inner_pt = (ra_max-0.001, dec_min+0.001)
                        hs = halfSpaceFromPoints(c1, c2, inner_pt)
                    elif ddec<tol and dra>tol:
                        # The Decs of the two corners is identical, bu the RAs are
                        # different; this Half Space is defined by a line of constant
                        # Dec and should be centered at one of the poles
                        if np.abs(c1[1]-dec_min)<tol:
                            hs = halfSpaceFromRaDec(0.0, 90.0, 90.0-dec_min)
                        else:
                            hs = halfSpaceFromRaDec(0.0, -90.0, 90.0+dec_max)
                    else:
                        continue

                    if hs is None:
                        raise RuntimeError("Somehow Half Space == None")
                    local_hs_list.append(hs)
            self._tile_half_spaces[tile_id] = local_hs_list

    def tile_ra(self, tile_idx):
        return self._tile_ra[tile_idx]

    def tile_dec(self, tile_idx):
        return self._tile_dec[tile_idx]

    def find_all_tiles(self, ra, dec, radius):
        """
        ra, dec, radius are all in degrees

        returns a numpy array of tile IDs that intersect the circle
        """
        valid_id = []
        radius_rad = np.radians(radius)
        center_pt = cartesianFromSpherical(np.radians(ra), np.radians(dec))
        for tile_id in self._tile_half_spaces:
            hs_list = self._tile_half_spaces[tile_id]
            is_contained = True
            for hs in hs_list:
                if not hs.intersects_circle(center_pt, radius_rad):
                    is_contained = False
                    break
            if is_contained:
                valid_id.append(tile_id)
        return np.array(valid_id)
