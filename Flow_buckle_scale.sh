
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
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp
  wget http://cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp
  wget http://cgmartini.nl/images/parameters/lipids/PC/POPC/martini_v2.0_POPC_02.itp
  mv martini_v2.0_POPC_02.itp martini_v2.0_POPR.itp
  sed -i -e "s/POPC/POPR/g" martini_v2.0_POPR.itp
fi

######################################################################################################
### PREPARE ##########################################################################################
######################################################################################################

if [ $PREPARE -eq 1 ]
then

## make bilayer with insane
python2 insane.py -o bilayer.gro -p topol.top -x 10 -y 30 -z 15 -l POPC:60 -l POPS:20 -l CHOL:20 -sol W:90 -sol WF:10 -salt 0.05 

## change topology file
sed -i -e '0,/POPC /{s/POPC /POPR /}' topol.top
cat << EOF > tmp1
#include "martini_v2.2.itp"
#include "martini_v2.0_ions.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "martini_v2.0_POPR.itp"

#ifdef POSRES_POPR
#include "posres_popr.itp"
#endif
EOF
tail -n +2 topol.top > tmp2
cat tmp1 tmp2 > topol.top
rm tmp1 tmp2

# position restraint on PO4's of POPCs in lower leaflet (named POPR)
cat << EOF > posres_popr.itp
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   2    1         20         20         20
EOF

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
gmx grompp -f mdp_files/min.mdp -c bilayer.gro -p topol.top -o min.tpr -quiet -maxwarn 1 # non-matching names: CNO - CN0
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
2 | 3 | 4 | 5
name 9 LIP
6 | 7 | 8
name 10 SOL
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

## python script for rescale x
cat << EOF > scale_gro_file.py
import numpy as np

filename_in = input('enter name of input .gro file:  ')
filename_out = input('enter name of output .gro file:  ')
scale_x = float(input('enter scale factor x:  '))
scale_y = float(input('enter scale factor y:  '))
scale_z = float(input('enter scale factor z:  '))
vel = int(input('gro file include velocities? (0 or 1):  '))

print('\ninput file: %s' % filename_in)
print('scaling factor x: %2.2f' % scale_x)
print('scaling factor y: %2.2f' % scale_y)
print('scaling factor z: %2.2f' % scale_z)
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
    if n in [3,5+n_vel,6+n_vel,7+n_vel]:
        idx = n-3-n_vel
        x,y,z = float(line_split[idx]),float(line_split[idx+1]),float(line_split[idx+2])
        x_scale,y_scale,z_scale = x*scale_x,y*scale_y,z*scale_z
        if n == 3:
            new_line = '%10.5f%10.5f%10.5f\n' % (x_scale,y_scale,z_scale)
        else:
            if vel:
                new_line = '%s%8.3f%8.3f%8.3f%s' %  (line[0:20],x_scale,y_scale,z_scale,line[44:])
            else:
                new_line = '%s%8.3f%8.3f%8.3f\n' %  (line[0:20],x_scale,y_scale,z_scale)
    else:
        new_line = line
    f_out.write(new_line)
    line = f_in.readline()
f_in.close()
f_out.close()
EOF

# equilibrate 10 ns, NPT
cat << EOF > mdp_files/npt.mdp
integrator               = md
dt                       = 0.033 
nsteps                   = 3030303 ; 100 ns
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 30300 ; ca every 1 ns 
nstenergy                = 30000
nstxout-compressed       = 30000
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
gmx grompp -f mdp_files/npt.mdp -c min.gro -p topol.top -n index.ndx -o npt_0.tpr -quiet
gmx mdrun -v -deffnm npt_0 -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 10000

fi

######################################################################################################
### RESCALE ##########################################################################################
######################################################################################################

n_ite=10 # needed by NVT, therefore outside if statement

if [ $RESCALE -eq 1 ]
then

vel=1 # velocities in gro file? (1: yes, 0= no)
scale_final=0.7
cat << EOF > calc_scale.py
print(${scale_final}**(1/$n_ite))
EOF
scale_y=$(python3 calc_scale.py)
rm calc_scale.py
echo scale_final = $scale_final
echo scale_y = $scale_y
           

## rescale y-coordinate iteratively, with equilibration
for i in $(seq 1 $n_ite)
do
  echo -------------------------
  echo iteration $i
  echo -------------------------
  
  i_prev=$((i-1))
  
  ## input file for rescaling y-coordinate
  cat << EOF > input_scale
npt_$i_prev.gro
npt_${i_prev}_scale.gro
1.0
$scale_y
1.0
$vel
EOF
 
  # rescale
  python3 scale_gro_file.py < input_scale
  
  # equilibrate
  gmx grompp -f mdp_files/npt.mdp -c npt_${i_prev}_scale.gro -p topol.top -n index.ndx -o npt_$i.tpr -quiet
  gmx mdrun -v -deffnm npt_$i -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 10000
  
  # vizualization scripts for pymol
  cat << EOF > script_npt_$i.pml
remove all
load npt_$i.gro
load_traj npt_$i.xtc
remove resname W
remove resname WF
remove resname ION
extract PR, resname POPR
color red, PR
set orthoscopic, on
show spheres
show cell
EOF

done

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
