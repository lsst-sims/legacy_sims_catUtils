import numpy as np
from MulticoreTSNE import MulticoreTSNE as TSNE

if __name__ == "__main__":

    n_features = 24
    feature_files = ['eb_features_180108.txt', 'rrly_features_180108.txt']
    dtype_list = []
    for ii in range(n_features):
        dtype_list.append(('f%d' % ii,float))
    dtype_list.append(('n_g', int))
    dtype_list.append(('n_i', int))
    
    dtype = np.dtype(dtype_list)
    data = None
    data_labels = None
    for i_file, file_name in enumerate(feature_files):
        local_data = np.genfromtxt(file_name, dtype=dtype)
        if data is None:
            data = local_data
            data_labels = np.array([i_file]*len(data))
        else:
            data = np.append(data, local_data, axis=0)
            data_labels = np.append(data_labels, np.array([i_file]*len(local_data)))

    raw_features = np.array([data['f%d' % ii] for ii in range(n_features)])
    assert len(data_labels) == len(raw_features[0])
    raw_features[n_features-3] = np.abs(raw_features[n_features-3])
    valid = np.where(np.logical_not(np.isnan(raw_features[n_features-2])))
    features = np.zeros((n_features, len(valid[0])), dtype=float)
    data_labels = data_labels[valid]
    for ii in range(n_features):
        features[ii] = raw_features[ii][valid]

    assert len(data_labels) == len(features[0])
    mean_features = np.array([np.mean(features[ii]) for ii in range(n_features)])
    
    covar = np.array([[np.mean((features[ii]-mean_features[ii])*(features[jj]-mean_features[jj]))
                       for ii in range(n_features)] for jj in range(n_features)])
    #print(covar.shape)
    #print(covar)
    #print(features.shape)
    e_val, e_vec = np.linalg.eig(covar)
    #print(e_val)
    #print((e_vec[:,2]**2).sum())
    e_vec_t = e_vec.transpose()
    samples = features.transpose()
    tsne_features = np.zeros((len(samples),n_features//2), dtype=float)
    sorted_dex = np.argsort(-1.0*np.abs(e_val))
    for i_s in range(len(samples)):
        for i_f in range(n_features//2):
            e_v_dex = sorted_dex[i_f]
            tsne_features[i_s][i_f] = np.dot(samples[i_s], e_vec_t[e_v_dex])
            tsne_features[i_s][i_f] /= np.sqrt(np.abs(e_val[e_v_dex]))
            assert not np.isnan(tsne_features[i_s][i_f])
            assert not np.isinf(tsne_features[i_s][i_f])
    print(tsne_features.shape)
    
    tsne_model = TSNE(n_jobs=10)
    tsne_result = tsne_model.fit_transform(tsne_features, perplexity=200.0)
    print(tsne_result.shape)
    
    for i_file, file_name in enumerate(feature_files):
        out_name =file_name.split('.')[0] + '_tsne_features.txt'
        with open(out_name, 'w') as out_file:
            valid = np.where(data_labels == i_file)
            tsne_valid = tsne_result[valid]
            for i_obj in range(len(tsne_valid)):
                out_file.write('%e %e\n' % (tsne_valid[i_obj][0], tsne_valid[i_obj][1]))
