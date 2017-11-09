from lsst.sims.utils import trixelFromHtmid
from lsst.sims.utils import levelFromHtmid
from lsst.sims.utils import angularSeparation

import numpy as np
rng = np.random.RandomState(182312)
n_pointings = 3*10**7
ra_list = rng.random_sample(n_pointings)*360.0
dec_list = rng.random_sample(n_pointings)*180.0-90.0

import time
t_start = time.time()
start_trixels = range(8,16)
n_bits_added = 12

pre_trixel = {}
for t0 in start_trixels:
    trix0 = trixelFromHtmid(t0)
    pre_trixel[t0] = trix0

ct = 0
for t0 in start_trixels:
    t0 = t0 << n_bits_added
    for dt in range(2**n_bits_added):
        htmid = t0 + dt
        ct += 1
        if htmid in pre_trixel:
            continue

        parent_id = htmid >> 2
        
        while parent_id not in pre_trixel:
            for n_right in range(2,n_bits_added,2):
                if htmid >> n_right in pre_trixel:
                    break
            to_gen = htmid >> n_right
            if to_gen in pre_trixel:
                trix0 = pre_trixel[to_gen]
            else:
                trix0= trixelFromHtmid(to_gen)
                pre_trixel[to_gen] = trix0

            pre_trixel[to_gen<<2] = trix0.get_child(0)
            pre_trixel[(to_gen<<2)+1] = trix0.get_child(1)
            pre_trixel[(to_gen<<2)+2] = trix0.get_child(2)
            pre_trixel[(to_gen<<2)+3] = trix0.get_child(3)

        trix0 = pre_trixel[parent_id]
        pre_trixel[(parent_id<<2)] = trix0.get_child(0)
        pre_trixel[(parent_id<<2)+1] = trix0.get_child(1)
        pre_trixel[(parent_id<<2)+2] = trix0.get_child(2)
        pre_trixel[(parent_id<<2)+3] = trix0.get_child(3)
        
        elapsed = time.time()-t_start
        #print('pregen %d len %d elapsed %e per %e' %
        #      (ct,len(pre_trixel),elapsed,elapsed/ct))
        

level = 7
level_seven_trixels = []
for htmid in pre_trixel:
    if levelFromHtmid(htmid) == level:
        level_seven_trixels.append(pre_trixel[htmid])

print('%d level %d trixels' % (len(level_seven_trixels), level))
print(8*(4**(level-1)))

t_start = time.time()

with open('trixel_mapping.txt', 'w') as out_file:
    for i_trixel in range(len(level_seven_trixels)):
        base_trixel = level_seven_trixels[i_trixel]
        neighbors = []
        for i_trixel_2 in range(len(level_seven_trixels)):
            if i_trixel_2 == i_trixel:
                continue
            is_neighbor = False
            for corner1 in base_trixel._corners:
                if is_neighbor:
                    break
                for corner2 in level_seven_trixels[i_trixel_2]._corners:
                    if np.abs(corner1-corner2).max()<0.001:
                        is_neighbor = True
                        break

            if is_neighbor:
                neighbors.append(level_seven_trixels[i_trixel_2].htmid)

        elapsed =time.time()-t_start
        print('%d took %e per %e-- %d' %
        (base_trixel.htmid, elapsed, elapsed/(i_trixel+1), len(neighbors)))
        out_file.write('%d ' % base_trixel.htmid)
        for hh in neighbors:
            out_file.write('%d ' % hh)
        out_file.write('\n')

exit()

ct = 0
for t0 in start_trixels:
    t0 = t0<<n_bits_added
    for dt in range(2**n_bits_added):
        htmid = t0+dt
        assert levelFromHtmid(htmid) == 7
        trixel = trixelFromHtmid(htmid)
        #trixel2 = pre_trixel[htmid+1]
        #np.testing.assert_array_equal(trixel._corners, trixel2._corners)
        #trixel = pre_trixel[htmid]
        ra0, dec0 = trixel.get_center()
        valid = np.where(angularSeparation(ra0, dec0, ra_list, dec_list)<1.5)

        ct += 1
        elapsed = time.time()-t_start
        print('%d took %e; per %e' % (ct, elapsed, elapsed/ct))
