#!/bin/bash

# calculate the covariance matrix (FIGURE S2)
# for only the C-alpha atoms
gmx covar -f md_500ns_centered.xtc -s md_500ns.tpr -o eigenval_ca.xvg -v eigenvec_ca.trr<<EOF
1
3
EOF

# analyse covariance matrix to rmsf
# protein
gmx anaeig -v eigenvec_ca.trr -s md_500ns.tpr -rmsf eigrmsf_ca.xvg<<EOF
3
EOF