#!/bin/bash

# How to calculate mmPBSA affinity values with the gromacs tool:

# Valdés-Tresanco, M.S., Valdés-Tresanco, M.E., Valiente, P.A. and Moreno E. 
# gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS. 
# Journal of Chemical Theory and Computation, 
# 2021 17 (10), 6281-6291 https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645

conda activate mmpbsa

# make index with prot and pocket-Calcium as one group
# for DC: find out which is the right CAL (which residue index) 
gmx make_ndx -f md_500ns.tpr -o mm.ndx <<EOF
1|ri [132/133] 
q
EOF

# for MRC there is only 1 calcium
gmx make_ndx -f md_500ns.tpr -o mm.ndx <<EOF
1|14
q
EOF



# calculation: choose new group with protein & calcium for receptor (here 18) and Ligand group for Ligand (here 13 or 14)
# CAREFUL: INDEX GROUP NUMBERS (MAN-DC=18 14)
# the mmpbsa.in input file can be found in the scripts folder
state=crystal
gmx_MMPBSA -O -i mmpbsa.in -cs md_500ns.tpr -ci mm.ndx -cg 18 13 -ct ${state}.xtc -cp topol.top -o mmpbsa_results.dat

