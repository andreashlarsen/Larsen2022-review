# Larsen2022-review
scripts for review:
Larsen(2022). Molecular dynamics simulations of curved lipid membranes.   

## prerequisites 
all scripts  require GROMACS, and were tested in GROMACS 2021.4 on a gpu-supported workstation.   
The scripts can be run with minor modifications, e.g. change of paths etc.   
all simulateions in the examples were done in Martini, but most of the methods are transfereable to be used with other force fields    

### EnCurv
the EnCurv setup require encurv installed with PLUMED, and PLUMED pathced before installing GROMACS. 

### Buckling
There are two scripts for membrane buckling, as descripbed in the review. One uses pressure, and the other scripts scales one axis in order to do step-wise buckling of the membrane. 

