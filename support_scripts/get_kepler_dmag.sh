#!/bin/bash

# these variables indicate from where to download the
# light curves
data_url=https://lsst-web.ncsa.illinois.edu/
data_dir=sim-data/parametrized_light_curves/
data_file=kplr_dmag_171204.txt

curl -O ${data_url}${data_dir}${data_file}

# this is the directory where the light curve file
# will be stored
destination_dir=${SIMS_DATA_DIR}/catUtilsData

# create destination_dir if it does not exist
if [ ! -e $destination_dir ]
then
    mkdir $destination_dir
fi

mv ${data_file} ${destination_dir}/${data_file}

echo "\nThe md5 checksum for "${destination_dir}/${data_file}" should be"
echo "9033d02f931551ccab837af1821a68f9"
