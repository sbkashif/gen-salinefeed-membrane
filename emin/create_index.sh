cp ../init-files/system-resnr.gro .

gmx make_ndx -f system-resnr.gro <<EOF
2 | 3 | 4 | 5 | 6
name 25 Polyamide
q
EOF

p ~/scripts/find-freeze-indices/find_freeze_indices.py -gro_file system-resnr.gro -index_file index.ndx -polymer_resname "MPD|TMC"
