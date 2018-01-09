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
    for file_name in feature_files:
        local_data = np.genfromtxt(file_name, dtype=dtype)
        if data is None:
            data = local_data
        else:
            data = np.append(data, local_data, axis=0)

    raw_features = np.array([data['f%d' % ii] for ii in range(n_features)])
    raw_features[n_features-3] = np.abs(raw_features[n_features-3])
    valid = np.where(np.logical_not(np.isnan(raw_features[n_features-2])))
    features = np.zeros((n_features, len(valid[0])), dtype=float)
    for ii in range(n_features):
        features[ii] = raw_features[ii][valid]

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
    tsne_features = np.zeros((len(samples),n_features), dtype=float)
    for i_s in range(len(samples)):
        for i_f in range(n_features):
            tsne_features[i_s][i_f] = np.dot(samples[i_s], e_vec_t[i_f])
    print(tsne_features.shape)
    
    tsne_model = TSNE(n_jobs=10)
    tsne_result = tsne_model.fit_transform(tsne_features)
    print(tsne_result.shape)
    
