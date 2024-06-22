#First command was run from Jupyter Notebook which have the output file cut.gro
#Jupyter notebook makes it easier to do some exploratory analysis on box sizes


#If you need to run it from here, uncomment the next line
# python ~/scripts/gen-salinefeed-membrane/obtain_polymer_domain.py \
#         -gro_file prod-padqi-npt-55-ns-copy-1.gro \
#         -dzp 0.6 \
#         -dzn 0.6 \
#         -o cut.gro \
#         -solvent_resnames "SOL" \
#         -polymer_resnames "MPD TMC" \
#         -n_atoms_water 4


gmx editconf    -f cut.gro \
                -o cut-renumbered.gro \
                -resnr 1

eval $(awk '/^[^;]/ { print "var"NR"="$1}' polymer_info.txt )

#var is read from polymer_info.txt

#2.8 is the desired permeate thickness, keeping consistency with Malmir et al.
perm_shift=$(awk "BEGIN {print 2.8-$var2; exit}")
echo "Premeate:$perm_shift"

#Adjusting the system to get the target permeate thickness
#pbc condition is to ensure that water molecules do not break at bondary
gmx editconf    -f cut-renumbered.gro \
                -o permeate-adjusted.gro \
                -translate 0 0 $perm_shift \
                -pbc no

#4.5 is the desired feed thickness, keeping consistency with Malmir et al.
feed_shift=$(awk "BEGIN {print $var4+$perm_shift+4.5; exit}")
echo "Feed:$feed_shift"

x=$(sed -n  '$p' cut-renumbered.gro | tr -s ' ' | cut -d ' ' -f 2)
y=$(sed -n  '$p' cut-renumbered.gro | tr -s ' ' | cut -d ' ' -f 3)
z=$(sed -n  '$p' cut-renumbered.gro | tr -s ' ' | cut -d ' ' -f 4)

#Extending the box length in the permeate side
#Making sure that the syystem is not centered
gmx editconf    -f permeate-adjusted.gro \
                -o feed-adjusted.gro \
                -box $x $y $feed_shift \
                -c no \
                -pbc no

feed_strip_water_count=$(($var8))
echo "Feed strip water count:$feed_strip_water_count"

permeate_strip_water_count=$(($var10))
echo "Permeate strip water count:$permeate_strip_water_count"

polymer_water_count=$(($var12))
echo "Polymer water count:$polymer_water_count"

#Determined from running determin_global_parameters.ipynb
total_water_mol=$((20938))

ions_each_type=$((231))

n_water_mol_add=$(awk "BEGIN {print $total_water_mol-$polymer_water_count-$feed_strip_water_count-$permeate_strip_water_count; exit}")
echo "Number of water molecules to add:$n_water_mol_add"

gmx solvate -cp feed-adjusted.gro -cs topology_solvated_TIP4P2005/tip4p.gro -o resolvated.gro -p topology_solvated_TIP4P2005/padqi-hydrolyzed.top -maxsol $n_water_mol_add

#cp padqi-hydrolyzed.top{,.bak}

gmx editconf -f resolvated.gro -o resolvated.gro -resnr 1

#Updating the topology file to account for updated number of water molecules

python ~/scripts/gen-salinefeed-membrane/modify_topology.py \
                -top_file topology_solvated_TIP4P2005/padqi-hydrolyzed.top \
                -o topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
                -gro resolvated.gro \
                -target_resname SOL

#Now adding ions in the gro file
#gmx solvate does not allow adding ions in specific z-range
#This is what this script takes care of

#First add NA ions, and then CL ions

python ~/scripts/domain-ionization/add-ions.py  -gro_file resolvated.gro \
                                                -o system-box-na.gro \
                                                -ion NA \
                                                -nn $ions_each_type \
                                                -nw 4

gmx editconf    -f system-box-na.gro \
                -o system-box-na-resnr.gro \
                -resnr 1

python ~/scripts/domain-ionization/add-ions.py \
            -gro_file system-box-na-resnr.gro \
            -o system-box-cl.gro \
            -ion CL \
            -nn $ions_each_type \
            -nw 4

gmx editconf    -f system-box-cl.gro \
                -o system-ionized.gro \
                -resnr 1

#Now update the topology file to account for the ions, and deleted water molecules

python ~/scripts/gen-salinefeed-membrane/modify_topology.py \
            -top_file topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
            -o topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
            -gro system-ionized.gro \
            -target_resname NA

python ~/scripts/gen-salinefeed-membrane/modify_topology.py \
            -top_file topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
            -o topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
            -gro system-ionized.gro \
            -target_resname CL

python ~/scripts/gen-salinefeed-membrane/modify_topology.py \
            -top_file topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
            -o topology_solvated_TIP4P2005/padqi-hydrolyzed-system.top \
            -gro system-ionized.gro \
            -target_resname SOL


