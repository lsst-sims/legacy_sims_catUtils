import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import h5py
import numpy as np

def make_2d_histogram(xx, yy, dx, dy):
    """
    returns indices and counts of unique points on the map
    """
    i_color1 = np.round(xx/dx).astype(int)
    i_color2 = np.round(yy/dy).astype(int)
    dex_reverse = np.array([i_color1, i_color2])
    dex_arr = dex_reverse.transpose()
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    dex_raw = np.ascontiguousarray(dex_arr).view(np.dtype((np.void, dex_arr.dtype.itemsize*dex_arr.shape[1])))
    _, unique_rows, unique_counts = np.unique(dex_raw, return_index=True, return_counts=True)

    return unique_rows, unique_counts


def plot_color(xx, yy, dx, dy):
    dexes, cts = make_2d_histogram(xx, yy, dx, dy)
    sorted_dex = np.argsort(cts)
    dexes = dexes[sorted_dex]
    cts = cts[sorted_dex]
    plt.scatter(xx[dexes], yy[dexes], c=cts, s=5,
                cmap=plt.cm.gist_ncar, edgecolor='')

    plt.colorbar()


def plot_color_mesh(xx, yy, dx, dy, vmin=None, vmax=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex %e %d %e %d' %
                           (xx.min(), i_x_arr.min(), yy.min(), i_y_arr.min()))

    x_mesh=np.arange(xx.min(),xx.max()+0.1,dx)
    y_mesh=np.arange(yy.min(),yy.max()+0.1,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    z_mesh = np.zeros(shape=x_mesh.shape, dtype=int)
    ct_1000b = 0

    for dex, ct in zip(dex_list, ct_list):
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]
        z_mesh[iy][ix] += ct

    z_mesh = np.ma.masked_where(z_mesh==0,z_mesh)
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax,
    cmap=plt.cm.gist_ncar)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=20)

    ct_1000 = 0
    big_min = np.round((2.8-xx.min())/dx).astype(int)
    big_max = x_mesh.shape[0]


def plot_color_mesh_set_color(xx, yy, color, dx, dy, vmin=None, vmax=None,
                              color_label=None):
    i_x_arr = np.round((xx-xx.min())/dx).astype(int)
    i_y_arr = np.round((yy-yy.min())/dy).astype(int)
    new_x = i_x_arr*dx
    new_y = i_y_arr*dy
    dex_list, ct_list = make_2d_histogram(new_x, new_y, dx, dy)

    if i_x_arr.min()<0 or i_y_arr.min()<0:
        raise RuntimeError('negative dex %e %d %e %d' %
                           (xx.min(), i_x_arr.min(), yy.min(), i_y_arr.min()))

    x_mesh=np.arange(xx.min(),xx.max()+0.1,dx)
    y_mesh=np.arange(yy.min(),yy.max()+0.1,dy)
    x_mesh,y_mesh = np.meshgrid(x_mesh,y_mesh,indexing='xy')
    z_mesh = np.ones(shape=x_mesh.shape, dtype=float)*(-999.0)

    for dex in dex_list:
        ix = i_x_arr[dex]
        iy = i_y_arr[dex]

        xmin = xx.min()+ix*dx-0.5*dx
        xmax = xx.min()+ix*dx+0.5*dx
        ymin = yy.min()+iy*dy-0.5*dy
        ymax = yy.min()+iy*dy+0.5*dy

        valid_pts = np.where(np.logical_and(xx>=xmin,
                             np.logical_and(xx<=xmax,
                             np.logical_and(yy>=ymin,yy<=ymax))))

        color_valid = color[valid_pts]

        if len(color_valid)>0:
            val = np.median(color_valid)
            z_mesh[iy][ix] = val

    z_mesh = np.ma.masked_where(z_mesh<-990.0,z_mesh)
    plt.pcolormesh(x_mesh,y_mesh,z_mesh, vmin=vmin, vmax=vmax)
                   #norm=matplotlib.colors.LogNorm(vmin=1.0,
                   #                               vmax=1.2e6))
    plt.colorbar(label=color_label)

if __name__ == "__main__":

    with h5py.File('sne_interp_test_data.h5', 'r') as in_file:
        plt.figure(figsize=(20,20))
        truth = in_file['mag_true'].value
        truth = np.where(np.isnan(truth), 30.0, truth)
        interp = in_file['mag_interp'].value
        interp = np.where(np.isnan(interp), 30.0, interp)
        for i_bp, bp in enumerate('ugrizy'):
            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(truth[i_bp],interp[i_bp],
                            0.05, 0.05)
            plt.axhline(25.0,color='r',linestyle='--')
            plt.axvline(25.0,color='r',linestyle='--')
            plt.title(bp,fontsize=20)
            plt.xlabel('truth', fontsize=20)
            plt.ylabel('interp', fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_validation.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(in_file['z'].value,interp[i_bp],
                            0.01, 0.05)
            plt.axhline(25.0, color='r', linestyle='--')
            plt.xlabel('redshift', fontsize=20)
            plt.ylabel('mag', fontsize=20)
            plt.title(bp, fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_mag_v_z.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(interp[i_bp], in_file['t'].value,
                            0.05, 0.05)
            plt.axvline(25.0, color='r', linestyle='--')
            plt.xlabel('mag', fontsize=20)
            plt.ylabel('t', fontsize=20)
            plt.title(bp, fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_mag_v_t.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(in_file['z'].value,
                            interp[i_bp]-truth[i_bp],
                            0.01, 0.05)
            plt.axhline(0.0, color='r', linestyle='--')
            plt.xlabel('redshift', fontsize=20)
            plt.ylabel('Dmag', fontsize=20)
            plt.title(bp, fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_Dmag_v_z.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            vanished = np.where(np.logical_and(truth[i_bp]<25.0,
                                               interp[i_bp]>25.0))

            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(in_file['z'].value[vanished],
                            interp[i_bp][vanished]-truth[i_bp][vanished],
                            0.01,0.05)
            plt.xlabel('redshift', fontsize=20)
            plt.ylabel('Dmag', fontsize=20)
            plt.title(bp,fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_Dmag_v_z_vanished.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            vanished = np.where(np.logical_and(truth[i_bp]<25.0,
                                               interp[i_bp]>25.0))

            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(in_file['c0'].value[vanished],
                            interp[i_bp][vanished]-truth[i_bp][vanished],
                            0.01,0.05)
            plt.xlabel('c0', fontsize=20)
            plt.ylabel('Dmag', fontsize=20)
            plt.title(bp,fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_Dmag_v_c0_vanished.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            vanished = np.where(np.logical_and(truth[i_bp]<25.0,
                                               interp[i_bp]>25.0))

            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(in_file['x1'].value[vanished],
                            interp[i_bp][vanished]-truth[i_bp][vanished],
                            0.01,0.05)
            plt.xlabel('x1', fontsize=20)
            plt.ylabel('Dmag', fontsize=20)
            plt.title(bp,fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_Dmag_v_x1_vanished.png')
        plt.close()

        plt.figure(figsize=(20,20))
        for i_bp, bp in enumerate('ugrizy'):
            vanished = np.where(np.logical_and(truth[i_bp]<25.0,
                                               interp[i_bp]>25.0))

            plt.subplot(3,2,i_bp+1)
            plot_color_mesh(in_file['t'].value[vanished],
                            interp[i_bp][vanished]-truth[i_bp][vanished],
                            0.01,0.05)
            plt.xlabel('t', fontsize=20)
            plt.ylabel('Dmag', fontsize=20)
            plt.title(bp,fontsize=20)
        plt.tight_layout()
        plt.savefig('sne_interp_Dmag_v_t_vanished.png')
        plt.close()
