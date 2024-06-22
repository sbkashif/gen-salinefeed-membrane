

gmx grompp -f em.mdp -c  system-resnr.gro  -p  ../topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top -po em-padqi-mdout.mdp -o em-padqi.tpr -n index.ndx

export OMP_NUM_THREADS=24
gmx mdrun -s em-padqi.tpr -v -deffnm em-padqi -ntomp 24 -ntmpi 1

gmx trjconv -f em-padqi.gro -s em-padqi.tpr -o em-padqi-whole.gro -pbc whole <<EOF
0
EOF

