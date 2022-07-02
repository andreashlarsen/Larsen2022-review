
## (de)select part of script
DOWNLOAD=0
PREPARE=0
NPT=0
RESCALE=1
NVT=1

## parallelisation
pinoffset=16
CPU=8
GPU=0

######################################################################################################
### DOWNLOAD #########################################################################################
######################################################################################################

if [ $DOWNLOAD -eq 1 ]
then
  ## download (and change) martini sctipts and itp files
  wget www.cgmartini.nl/images/tools/insane/insane.py
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.2.itp
  #mv martini_v2.2.itp martini_v2.2_repulsive_Qa.itp
  #manually add QaC and QaS beads to martini_v2.2_repulsive_Qa.itp
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp
  wget http://cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp
  wget http://cgmartini.nl/images/parameters/lipids/PC/POPC/martini_v2.0_POPC_02.itp
  mv martini_v2.0_POPC_02.itp martini_v2.0_POPCQ_02.itp
  sed -i -e 's/Qa /QaC/g' martini_v2.0_POPCQ_02.itp
  sed -i -e 's/POPC /POPCQ/g' martini_v2.0_POPCQ_02.itp
  wget http://cgmartini.nl/images/parameters/lipid/PS/POPS/martini_v2.0_POPS_02.itp
  mv martini_v2.0_POPS_02.itp martini_v2.0_POPSQ_02.itp
  sed -i -e 's/Qa /QaS/g' martini_v2.0_PQPSr_02.itp
  sed -i -e 's/POPS /POPSQ/g' martini_v2.0_POPSQ_02.itp
fi

######################################################################################################
### PREPARE ##########################################################################################
######################################################################################################

if [ $PREPARE -eq 1 ]
then

## python script for translate coordinates in gro file
cat << EOF > translate_gro_file.py
import numpy as np

filename_in = input('enter name of input .gro file:  ')
filename_out = input('enter name of output .gro file:  ')
trans_x = float(input('translate x with :  '))
trans_y = float(input('translate y with:  '))
trans_z = float(input('translate z with:  '))
vel = int(input('gro file include velocities? (0 or 1):  '))

print('\ninput file: %s' % filename_in)
print('translate x: %2.2f' % trans_x)
print('translate y: %2.2f' % trans_y)
print('translate z: %2.2f' % trans_z)
print('output file: %s' % filename_out)

n_vel = 0
if vel:
    n_vel = 3

f_in = open(filename_in,'r')
f_out = open(filename_out,'w')

line = f_in.readline()
while line:
    line_split = line.split()
    n = len(line_split)
    if n == 3 or ( n == 9 and vel == 0):
        x,y,z = float(line_split[0]),float(line_split[1]),float(line_split[2])
        x_trans,y_trans,z_trans = x+trans_x,y+trans_y,z+trans_z
        new_line = '%10.5f%10.5f%10.5f\n' % (x_trans,y_trans,z_trans)
    elif n in [3,5+n_vel,6+n_vel,7+n_vel]:
        idx = n-3-n_vel
        x,y,z = float(line_split[idx]),float(line_split[idx+1]),float(line_split[idx+2])
        x_trans,y_trans,z_trans = x+trans_x,y+trans_y,z+trans_z
        if vel:
            new_line = '%s%8.3f%8.3f%8.3f%s' %  (line[0:20],x_trans,y_trans,z_trans,line[44:])
        else:
            new_line = '%s%8.3f%8.3f%8.3f\n' %  (line[0:20],x_trans,y_trans,z_trans)
    else:
        new_line = line
    f_out.write(new_line)
    line = f_in.readline()
f_in.close()
f_out.close()
EOF

## input file for translating y-coordinate
cat << EOF > input_translate
bilayer2.gro
bilayer2_translate.gro
0.0
15.0
0.0
0
EOF

## translate y-coordinates in bilayer2.gro file
python3 translate_gro_file.py < input_translate

## vizualize pymol
cat << EOF > bilayer.gro
load bilayer2_translate.gro
load bilayer1.gro
remove resname W
remove resname WF
remove resname NA\+
show spheres
color blue, resname POPC
color red, resname POPS
extract PC1, (resname POPC and bilayer1)
extract PS1, (resname POPS and bilayer1)
extract PC2, (resname POPC and bilayer2)
extract PS2, (resname POPS and bilayer2)
EOF


## make bilayer with insane
python2 insane.py -o bilayer1.gro -p topol1.top -x 10 -y 15 -z 15 -l POPC -u POPS -sol W:90 -sol WF:10 -salt 0.05 
python2 insane.py -o bilayer2.gro -p topol2.top -x 10 -y 15 -z 15 -l POPS -u POPC -sol W:90 -sol WF:10 -salt 0.05 

head -2 bilayer1.gro >tmp
N1=$(tail -1 tmp)
head -2 bilayer2.gro >tmp
N2=$(tail -1 tmp)
N=$((N1+N2))
echo N1 = $N1 and N2 = $N2 and N = N1 + N2 = $N

## merge bilayer files
cat << EOF > tmp1
INSANE! Membrane UpperLeaflet>POPS=1 LowerLeaflet>POPC=1, vice versa in other half. 
$N
EOF
tail -n +3 bilayer1.gro > tmp
head -n -2 tmp > tmp2
tail -n +3 bilayer2_translate.gro > tmp3
cat tmp1 tmp2 tmp3 > bilayer_merge.gro

## change topology file
cat << EOF > tmp1
#include "martini_v2.2_repulsive_Qa.itp"
#include "martini_v2.0_ions.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "martini_v2.0_POPCQ_02.itp"
#include "martini_v2.0_POPSQ_02.itp"
EOF
tail -n +2 topol1.top > tmp2
tail -n +9 topol2.top > tmp3
cat tmp1 tmp2 tmp3 > topol.top
rm tmp tmp1 tmp2 tmp3 
sed -i -e 's/POPC /POPCQ/g' topol.top
sed -i -e 's/POPS /POPSQ/g' topol.top

## minimization
mkdir -p mdp_files
cat << EOF > mdp_files/min.mdp
integrator               = steep
nsteps                   = 20000
nstxout                  = 0
nstfout                  = 0
nstlog                   = 100 
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
EOF
gmx grompp -f mdp_files/min.mdp -c bilayer_merge.gro -p topol.top -o min.tpr -quiet -maxwarn 1 # non-matching names: CNO - CN0
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

echo Potential | gmx energy -f min.edr -o potential.xvg -quiet
#xmgrace potential.xvg

# vizualization script for pymol
cat << EOF > min.pml
remove all
load min.gro
remove resname W
remove resname WF
remove resname ION
extract PR, resname POPR
color red, PR
set orthoscopic, on
show spheres
show cell
EOF

## make index file with LIP and SOL
cat << EOF > index_input
2 | 3
name 7 LIP
4 | 5 | 6
name 8 SOL
q
EOF
gmx make_ndx -f min.gro -quiet < index_input
rm index_input

fi 

######################################################################################################
### NPT ##############################################################################################
######################################################################################################

if [ $NPT -eq 1 ]
then

# equilibrate 10 ns, NPT
cat << EOF > mdp_files/npt.mdp
integrator               = md
dt                       = 0.033 
nsteps                   = 303030303 ; 10 us
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 909000 ; ca every 30 ns 
nstenergy                = 909000
nstxout-compressed       = 909000
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL
tau_t                    = 1.0  1.0 
ref_t                    = 320 320 
Pcoupl                   = berendsen
Pcoupltype               = anisotropic
tau_p                    = 12.0
compressibility          = 0.0  0.0  3e-4  0  0  0  
ref_p                    = 0.0  0.0  1.0   0  0  0
gen_vel                  = yes
gen_temp                 = 320
gen_seed                 = 473529
constraints              = none 
constraint_algorithm     = Lincs
EOF
gmx grompp -f mdp_files/npt.mdp -c min.gro -p topol.top -n index.ndx -o npt.tpr -quiet
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 30000 # short run to get .gro fil
gmx mdrun -v -deffnm npt -cpi npt.cpt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset # continue

# vizualization scripts for pymol
cat << EOF > npt.pml
remove all
load npt.gro
load_traj npt.xtc
remove resname W
remove resname WF
remove resname ION
extract PC, resname POPCQ
color blue, PC
extract PS, resname POPSQ
color red, PS
color grey50, npt
set orthoscopic, on
show spheres
show cell, npt
EOF

fi

######################################################################################################
### NVT ##############################################################################################
######################################################################################################

if [ $NVT -eq 1 ]
then

## equilibrate with lower leaflet POPC's (named POPR) contrained 

echo -----------------------------------------------------------
echo NVT run, with posres
echo -----------------------------------------------------------

cat << EOF > mdp_files/nvt.mdp
define                   = -DPOSRES_POPR
refcoord_scaling         = com
integrator               = md
dt                       = 0.033 
nsteps                   = 90909091 ; 3 us
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 303000 ; every 10 ns 
nstenergy                = 303000
nstxout-compressed       = 303000
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL
tau_t                    = 1.0  1.0 
ref_t                    = 320 320 
Pcoupl                   = no
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
EOF
gmx grompp -f mdp_files/nvt.mdp -c npt_$n_ite.gro -r npt_$n_ite.gro -p topol.top -n index.ndx -o nvt.tpr -quiet
gmx mdrun -v -deffnm nvt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 10000 

## vizualization scripts for pymol
cat << EOF > script_nvt.pml
remove all
load nvt.gro
load_traj nvt.xtc
remove resname W
remove resname WF
remove resname ION
extract PR, resname POPR
color red, PR
set orthoscopic, on
show spheres
show cell
EOF

fi 

## clean up
rm \#*
rm step*.pdb
