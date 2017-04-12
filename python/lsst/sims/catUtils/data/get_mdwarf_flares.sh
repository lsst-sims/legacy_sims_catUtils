#!/bin/bash

data_url=https://lsst-web.ncsa.illinois.edu/
data_dir=sim-data/mdwarf_flare_light_curves/
data_file=mdwarf_flare_light_curves_170412.npz

curl -O ${data_url}${data_dir}${data_file}

echo "\nThe md5 checksum for "${data_file}" should be"
echo "5d62eba4ad496ab700f11142ced348f7"
