#!/bin/bash

# How to do generate a state specific traj (binning together frames with the same binding state)
# this takes as input frame numbers determined by the jupyter notebook of each simulation


state=stateB # name of the state
mkdir ${state}
first_frame="59400" 
frame_nums="59600
59900
60000
60800
60900
61000
61100
61300
61600
61700
62800
63100
63500
64000"

# get first frame
gmx trjconv -s md_500ns.tpr -f bound.xtc -o ${state}.xtc -pbc mol -dump ${first_frame} -n index.ndx -center -ur compact <<EOF
0
2
EOF

# run loop over frame numbers
for i in $frame_nums
do

# get pdb file of everything but only 1 frame and molecule centered
gmx trjconv -s md_500ns.tpr -f bound.xtc -o ${state}/${i}.xtc -pbc mol -dump ${i} -n index.ndx -center -ur compact <<EOF
0
2
EOF

# concatenate all frames to one trajectory
gmx trjcat -f ${state}.xtc ${state}/${i}.xtc -o ${state}.xtc -cat

done