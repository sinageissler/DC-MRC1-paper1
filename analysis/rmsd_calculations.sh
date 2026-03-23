#!/bin/bash

############ RMSD calculation for pocket calcium and ligand (FIGURE S3) #############################
# these calculations were done on the concatenated trajectories of all three replica of each receptor-ligand combination

# Man - DC simulation 
# make an index for pocket
gmx make_ndx -f md_500ns.pdb -o rmsd.ndx <<EOF
ri 94-113
name 18 pocket
ri 132
name 19 pocketCAL
q
EOF

# center traj on pocket -- important to use -pbc res (not mol!!)
gmx trjconv -f combined_centered.xtc -s md_500ns.pdb -o pocket_centered.xtc -n rmsd.ndx -center -pbc res -ur compact <<EOF
18
0
EOF

# calculate RMSD cal to pocket
gmx rms -f pocket_centered.xtc -s md_500ns.pdb -o stats/cal_center_rmsd.xvg -n rmsd.ndx <<EOF
4
19
EOF

# calculate RMSD Lig to pocket
gmx rms -f pocket_centered.xtc -s md_500ns.pdb -o stats/rmsd_lig.xvg -n rmsd.ndx <<EOF
4
14
EOF

# DC Sims - bigger ligands (different indices)
# make an index for pocket
gmx make_ndx -f md_500ns.pdb -o rmsd.ndx <<EOF
ri 94-113
name 18 pocket
ri 133
name 19 pocketCAL
q
EOF

# center traj on pocket -- important to use -pbc res (not mol!!)
gmx trjconv -f combined_centered.xtc -s md_500ns.pdb -o pocket_centered.xtc -n rmsd.ndx -center -pbc res -ur compact <<EOF
18
0
EOF

# calculate RMSD cal to pocket
gmx rms -f pocket_centered.xtc -s md_500ns.pdb -o stats/cal_center_rmsd.xvg -n rmsd.ndx <<EOF
4
19
EOF

# calculate RMSD Lig to pocket
gmx rms -f pocket_centered.xtc -s md_500ns.pdb -o stats/rmsd_lig.xvg -n rmsd.ndx <<EOF
4
13
EOF


# MRC Sims (different indices)
# make an index for pocket
gmx make_ndx -f md_500ns.tpr -o rmsd.ndx <<EOF
ri 86-109
name 18 pocket
ri 124
name 19 pocketCAL
q
EOF

# center traj on pocket -- important to use -pbc res (not mol!!)
gmx trjconv -f combined_centered.xtc -s md_500ns.tpr -o pocket_centered.xtc -n rmsd.ndx -center -pbc res -ur compact <<EOF
18
0
EOF

# calculate RMSD cal to pocket
gmx rms -f pocket_centered.xtc -s md_500ns.pdb -o stats/cal_center_rmsd.xvg -n rmsd.ndx <<EOF
4
19
EOF

# calculate RMSD Lig to pocket
gmx rms -f pocket_centered.xtc -s md_500ns.pdb -o stats/rmsd_lig.xvg -n rmsd.ndx <<EOF
4
13
EOF


############# RMSD of protein (for figure 3) #######################
# these calculations were done on each replica individually

# RMSD of full prot
gmx rms -s md_500ns.pdb -f md_500ns_centered.xtc -o rmsd_prot.xvg <<EOF
1
1
EOF

# RMSD of pocket
# make an index for the pocket (ri 94-113 for DC, ri86-109 for MRC)
gmx make_ndx -f md_500ns.tpr -o rmsd.ndx <<EOF
ri 86-109
name 18 pocket
q
EOF

# calculate RMSD
gmx rms -f md_500ns_centered.xtc -s md_500ns.pdb -o rmsd_pocket.xvg -n pocket.ndx <<EOF
1
18
EOF

# RMSD of full prot aligned to first frame of first equilibration step (equilibration figure)
gmx rms -s step4.1_equilibration.pdb -f md_500ns_centered.xtc -o rmsd_prod_prot.xvg <<EOF
1
1
EOF


############## RMSF calculation on simulations without ligand (FIGURE 3b) ###################

gmx rmsf -f md_500ns_centered.xtc -o rmsf.xvg -res -s md_500ns.tpr<<EOF
1
EOF