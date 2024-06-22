#Place the polymer in the ceneter of the box such that there is 0.3 nm space both at the top and bottom

#Extract the target X,Y,Z coordinates
old_x=$(sed -n  '$p' ../system-ionized.gro | tr -s ' ' | cut -d ' ' -f 2)
old_y=$(sed -n  '$p' ../system-ionized.gro | tr -s ' ' | cut -d ' ' -f 3)
old_z=$(sed -n  '$p' ../system-ionized.gro | tr -s ' ' | cut -d ' ' -f 4)

#XY dimensions same as equilibrated sturucture
x=$(awk "BEGIN {print $old_x}") 
y=$(awk "BEGIN {print $old_y}")
 
#Adding 0.3 nm on blank space on each side for graphene sheet
z=$(awk "BEGIN {print $old_z+2}") 

echo "X:$x"
echo "Y:$y"
echo "New Z:$z"



gmx editconf -f ../system-ionized.gro -o system-adjusted.gro -box $x $y $z -c


cat system-adjusted.gro ../gra-setup/walls-box.gro > system.gro

n_hydrated_polymer_atoms=$(sed -n 2p system-adjusted.gro) 

#echo $n_atoms

n_graphene_atoms=$(sed -n 2p ../gra-setup/walls-box.gro)

start_d=$(($n_hydrated_polymer_atoms + 3))
end_d=$(($n_hydrated_polymer_atoms + 5))

total=$(($n_hydrated_polymer_atoms+$n_graphene_atoms))
#echo $start_d
#echo $end_d

#delete overlapping text in the new .gro file
sed -i -e "${start_d},${end_d}d;" system.gro

#change number at the top of the .gro file to reflect additional atoms
sed -i "2s/.*/$total/"  system.gro

z=$(awk "BEGIN {print $old_z+2.0+20.0}") 


gmx editconf -f system.gro -o system-vacuum.gro -box $x $y $z -c



gmx editconf -f system-vacuum.gro -o system-resnr.gro -resnr 1



python ~/scripts/gen-salinefeed-membrane/modify_topology.py -top_file ../topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top -o ../topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top -gro system-resnr.gro -target_resname TOP

python ~/scripts/gen-salinefeed-membrane/modify_topology.py -top_file ../topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top -o ../topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top -gro system-resnr.gro -target_resname BTM





