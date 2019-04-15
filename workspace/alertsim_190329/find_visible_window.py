import h5py
import numpy as np

m5_single = {}
m5_single['u'] = 23.57
m5_single['g'] = 24.65
m5_single['r'] = 24.21
m5_single['i'] = 23.79
m5_single['z'] = 23.21
m5_single['y'] = 22.31

t_min = 100000.0
t_max = -t_min

abs_mag_grid = np.arange(-20.5,-18.1,0.1)

vis_log = open('data/visible_sn_tags.txt' ,'w')

with open('data/invisible_sn_tags.txt', 'w') as invisible_log:
    with h5py.File('data/sne_interp_models.h5', 'r') as in_file:
        param_mins = in_file['param_mins'].value
        m0 = param_mins[3]
        t_grid = in_file['t_grid'].value
        for kk in in_file.keys():
            if kk == 'param_mins':
                continue
            elif kk == 'd_params':
                continue
            elif kk == 'param_names':
                continue
            elif kk == 't_grid':
                continue

            is_visible = False
            mag_grid = in_file[kk].value
            for i_bp, bp in enumerate('ugrizy'):
                for abs_mag in abs_mag_grid:
                    d_mag = abs_mag-m0
                    visible = np.where(d_mag+mag_grid[i_bp]<=m5_single[bp])
                    if len(visible[0])==0:
                        continue
                    is_visible = True
                    t_v = t_grid[visible]
                    t0 = t_v.min()
                    t1 = t_v.max()
                    if t0<t_min:
                        t_min = t0
                    if t1>t_max:
                        t_max = t1
                        print(t_max,kk)

            if not is_visible:
                invisible_log.write('%s\n' % kk)
            else:
                vis_log.write('%s\n' % kk)

print('%e <= t <= %e' % (t_min, t_max))
vis_log.close()
