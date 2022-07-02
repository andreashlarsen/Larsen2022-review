# compile ts2cg files
#cd TS2CG1.1-master
#bash compile.sh
#cd ../

# go to working dir
dir=5
mkdir -p $dir
cd $dir

# run PCG
martini=../TS2CG1.1-master/Tutorials/files/Martini3.LIB
x=20
y=10
z=20
cat<<EOF > input.str

[Shape Data]
ShapeType 1D_PBC_Fourier
Box $x $y $z
WallRange 0 1 0 1
Density 3 3
Thickness 2.8
Mode 2 1 0 # amplitude,periods per box,??
End

[Lipids List]
Domain 0
POPC 1 1 0.75
End
EOF
../TS2CG1.1-master/PCG -str input.str -Bondlength 0.2 -LLIB $martini -defout system -function analytical_shape -Wall -WallBName WL 

# vizualize
cat<<EOF > system.pml
load system.gro
show spheres
color red, resname Wall
extract wall, resname Wall
set orthoscopic, on
show cell
EOF
pymol system.pml

# dry minimization
#echo 3 | gmx trjconv -f system.gro -o wall.gro -quiet
#echo 2 | gmx genrestr -f wall.gro -o posre.itp -quiet
cat<<EOF > posre.itp
; position restraints for Wall of System

[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
EOF
ffdir=../TS2CG1.1-master/Tutorials/files/itp/martini3
cat<<EOF > topol.add
#include "$ffdir/martini_v3.0.4.itp"
#include "$ffdir/martini_v3.0_phospholipids.itp"
#include "$ffdir/Wall.itp"
#ifdef POSRES
#include "posre.itp"
#endif
#include "$ffdir/martini_v3.0_solvents.itp"
EOF
cat topol.add system.top > topol.top
rm topol.add

gmx grompp -f ../mdp_files/min.mdp -c system.gro -r system.gro -p topol.top -o min_dry.tpr -quiet
CPU=4
GPU=1
pinoffset=8
gmx mdrun -deffnm min_dry -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

# vizualize
cat<<EOF > min_dry.pml
load min_dry.gro
show spheres
color red, resname Wall
extract wall, resname Wall
set orthoscopic, on
show cell
EOF
pymol min_dry.pml

# solvation
python2 ../insane.py -f min_dry.gro -sol W -o min_solv.gro -p water.top -salt 0 -x $x -y $y -z $z
tail -1 water.top > tmp
cat topol.top tmp > topol_solv.top
#python2 insane.py -SOL W -water.gro -p water.top -salt 0 -x $x -y $y -z $z # make water box with insane
#../TS2CG1.1-master/SOL -in min.gro -tem water.gro -unsize 3 -ion 0 0 -o solvated.gro -Rcutoff 0.32

# minimization of solvated system
gmx grompp -f ../mdp_files/min.mdp -c min_solv.gro -r min_solv.gro -p topol_solv.top -o min.tpr -quiet
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

# vizualize
cat<<EOF > min.pml
load min.gro
show spheres
color red, resname Wall
extract wall, resname Wall
extract SOL, resname W
set orthoscopic, on
show cell
EOF
pymol min.pml

# equilibration
gmx grompp -f ../mdp_files/eq.mdp -c min.gro -r min.gro -p topol_solv.top -o eq.tpr -quiet
gmx mdrun -deffnm eq -v -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

# monitor temperature
echo Temperature | gmx energy -f eq.edr -o temperature.xvg -quiet
xmgrace temperature.xvg

# vizualize
cat<<EOF > eq.pml
load eq.gro
show spheres
color red, resname Wall
extract wall, resname Wall
extract SOL, resname W
set orthoscopic, on
show cell
EOF
pymol eq.pml

# production
gmx grompp -f md.mdp -c eq.gro -r eq.gro -p topol_solv.top -o md.tpr -quiet

## clean up
rm \#*
cd ..

