out_dir=/astro/store/pogo4/danielsf/alertsim_190329/baseline/

python -W'ignore' simulate_agn_alerts.py \
--circular_fov --htmid_file data/region_4_trixels.txt \
--out_name=${out_dir}region_4_agn.h5 \
--n_threads 12

echo "completed agn region 4"

python -W'ignore' simulate_agn_alerts.py \
--circular_fov --htmid_file data/region_1_trixels.txt \
--out_name=${out_dir}region_1_agn.h5 \
--n_threads 12

echo "completed agn region 1"

python -W'ignore' simulate_agn_alerts.py \
--circular_fov --htmid_file data/region_2_trixels.txt \
--out_name=${out_dir}region_2_agn.h5 \
--n_threads 12

echo "completed agn region 2"

python -W'ignore' simulate_agn_alerts.py \
--circular_fov --htmid_file data/region_3_trixels.txt \
--out_name=${out_dir}region_3_agn.h5 \
--n_threads 12

echo "completed agn region 3"
