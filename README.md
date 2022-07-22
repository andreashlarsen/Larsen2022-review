# Scripts for Larsen2022-review
Scripts for review article:     
Larsen (2022). Molecular dynamics simulations of curved lipid membranes. Int. J. Mol. Sci. 28: 8098.     
https://www.mdpi.com/1422-0067/23/15/8098    

## Prerequisites 
All scripts  require GROMACS, and were tested in GROMACS 2021.4 on a GPU-supported workstation.   
The scripts can be run with minor modifications, e.g. change of paths etc.   
All simulateions in the examples were done in Martini2.2, but most of the methods are transfereable to be used with other force fields (more info in the paper).        

### EnCurv 
The EnCurv setup requires that EnCurv is installed with PLUMED, and PLUMED (with EnCurv) pathced before installing GROMACS. Tested with PLUMED version 2.7.2.    

## Notes

### Buckling notes
There are two scripts for membrane buckling, as described in the review. One uses pressure for buckling the membrane. The other scripts scales all coordinates along one axis for step-wise buckling of the membrane. 

### Use/cite
Please use the scripts for your research. I would appreciate acknowledgement by citing the review paper, and, of course, the relevant original papers. 
