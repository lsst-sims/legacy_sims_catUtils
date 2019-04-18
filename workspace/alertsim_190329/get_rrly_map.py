import os
rrly_dir = os.path.join(os.environ['SIMS_SED_LIBRARY_DIR'], 'rrly_lc')
assert os.path.isdir(rrly_dir)

used_id = set()
with open('data/rrly_lc_map.txt', 'w') as out_file:
    for sub_dir in ['RRab', 'RRc']:
        full_dir = os.path.join(rrly_dir, sub_dir)
        contents = os.listdir(full_dir)
        for name in contents:
            if not name.endswith('per.txt'):
                continue
            map_name = os.path.join(sub_dir, name)
            lc_id = int(name.replace('_per.txt', ''))
            assert lc_id not in used_id
            used_id.add(lc_id)
            out_file.write('%d;%s\n' % (lc_id, map_name))
