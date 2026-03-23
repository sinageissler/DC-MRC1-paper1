
# These are all the scripts that I used to produce and analyze the data of this paper:
A comparative investigation of the mannose binding interface in DC-SIGN and MRC1 carbohydrate recognition domains with all-atom molecular dynamics simulations, Sina Geissler and Sophie Sacquin-Mora, Laboratoire de Biochimie Théorique, Université Paris-Cité, CNRS


## Simulation scripts
- run scripts to run equilibration and production
- .mdp files

## Analysis Scripts

You can find several different notebooks for the calculations contributing to this paper (and the figures):
- Notebooks for calculations that I did separately for each simulation (to separate frames into binding states)
- A Notebook for calculations that I did for all simulations at once ('all_trajs_analysis')
- A Notebook for 'summary figures' of results obtained by calculations not done with python (Figures 5 and 7)

You can also find bash scripts for:
- RMSD and RMSF calculations with Gromacs
- generating trajectories for each binding state
- mmPBSA calculations