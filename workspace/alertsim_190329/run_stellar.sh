out_dir=/astro/store/pogo4/danielsf/alertsim_190329/baseline/

python -W'ignore' simulate_stellar_alerts.py \
--circular_fov --htmid_file data/region_4_trixels.txt \
--out_name=${out_dir}region_4_stars.h5 \
--n_threads 12

echo "completed stars region 4"

python -W'ignore' simulate_stellar_alerts.py \
--circular_fov --htmid_file data/region_1_trixels.txt \
--out_name=${out_dir}region_1_stars.h5 \
--n_threads 12

echo "completed stars region 1"

python -W'ignore' simulate_stellar_alerts.py \
--circular_fov --htmid_file data/region_2_trixels.txt \
--out_name=${out_dir}region_2_stars.h5 \
--n_threads 12

echo "completed stars region 2"

python -W'ignore' simulate_stellar_alerts.py \
--circular_fov --htmid_file data/region_3_trixels.txt \
--out_name=${out_dir}region_3_stars.h5 \
--n_threads 12

echo "completed stars region 3"
