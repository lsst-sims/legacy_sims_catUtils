out_dir=/astro/store/pogo4/danielsf/alertsim_190329/

python -W'ignore' simulate_sne_alerts.py \
--circular_fov --htmid_file data/region_4_trixels.txt \
--out_name=${out_dir}region_4_sne.pickle --fast_t0

echo "completed sne region 4"

python -W'ignore' simulate_sne_alerts.py \
--circular_fov --htmid_file data/region_1_trixels.txt \
--out_name=${out_dir}region_1_sne.pickle --fast_t0

echo "completed sne region 1"

python -W'ignore' simulate_sne_alerts.py \
--circular_fov --htmid_file data/region_2_trixels.txt \
--out_name=${out_dir}region_2_sne.pickle --fast_t0

echo "completed sne region 2"

python -W'ignore' simulate_sne_alerts.py \
--circular_fov --htmid_file data/region_3_trixels.txt \
--out_name=${out_dir}region_3_sne.pickle --fast_t0

echo "completed sne region 3"
