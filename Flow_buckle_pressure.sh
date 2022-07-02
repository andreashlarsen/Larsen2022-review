## follow protocol in Boyd2017

## (de)select part of script
DOWNLOAD=0

## parallelisation
pinoffset=16
CPU=8
GPU=0

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

## loop over different latteral pressures (in y-direction)
#for p_y in 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 
for p_y in 3.5 
do

## equilibration, no position restraint

echo -----------------------------------------------------------
echo NPT equilibration, no posres, pressure in y-dir is $p_y bar
echo -----------------------------------------------------------

cat << EOF > mdp_files/eq.mdp
integrator               = md
dt                       = 0.033 
nsteps                   = 90909091 ; 3 us 
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 303000 ; ca every 10 ns (600 frames)
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
Pcoupl                   = berendsen
Pcoupltype               = anisotropic
tau_p                    = 12.0
compressibility          = 0.0  3e-4  3e-4  0  0  0  
ref_p                    = 0.0  $p_y  1.0   0  0  0
gen_vel                  = yes
gen_temp                 = 320
gen_seed                 = 473529
constraints              = none 
constraint_algorithm     = Lincs
continuation             = no
EOF
gmx grompp -f mdp_files/eq.mdp -c min.gro -p topol.top -n index.ndx -o eq_p$p_y.tpr -quiet
gmx mdrun -v -deffnm eq_p$p_y -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

## equilibrate with lower leaflet POPC's (named POPR) contrained 

echo -----------------------------------------------------------
echo NVT run, with posres
echo -----------------------------------------------------------

cat << EOF > mdp_files/eq_res.mdp
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
EOF
gmx grompp -f mdp_files/eq_res.mdp -c eq_p$p_y.gro -r eq_p$p_y.gro -p topol.top -n index.ndx -o eq_res_p$p_y.tpr -quiet
gmx mdrun -v -deffnm eq_res_p$p_y -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset 

## vizualization scripts for pymol

cat << EOF > script_p$p_y.pml
remove all
load eq_p$p_y.gro
load_traj eq_p$p_y.xtc
remove resname W
remove resname WF
remove resname ION
extract PR, resname POPR
color red, PR
set orthoscopic, on
show spheres
show cell
EOF

cat << EOF > script_res_p$p_y.pml
remove all
load eq_res_p$p_y.gro
load_traj eq_res_p$p_y.xtc
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

## clean up
rm \#*
rm step*.pdb
