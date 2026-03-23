#!/bin/bash

equi_steps="1 2 3 4 5 6"
for i in $equi_steps
do

# now center molecule in whole xtc
gmx trjconv -s step4.${i}_equilibration.tpr -f step4.${i}_equilibration.xtc -o step4.${i}_equilibration_centered.xtc -pbc mol -n index.ndx -center -ur compact <<EOF
0
2
EOF

# get a pdb file of first frame
gmx trjconv -s step4.${i}_equilibration.tpr -f step4.${i}_equilibration.xtc -o step4.${i}_equilibration.pdb -pbc mol -dump 0 -n index.ndx -center -ur compact <<EOF
0
2
EOF

# RMSD of full prot aligned to first frame of first equilibration step
gmx rms -s step4.1_equilibration.pdb -f step4.${i}_equilibration_centered.xtc -o rmsd_step4.${i}_equilibration_prot.xvg <<EOF
1
1
EOF

done