
source ~/.bashrc

pinoffset=0
CPU=4
GPU=0

# make bilayer with insane
/usr/bin/python2.7 insane.py -o bilayer.gro -p topol.top -x 30 -y 7 -z 15 -l POPC -sol W:9 -sol WF:1 -salt 0.1 -neutral

# change topology file
sed -i -e 's/#include "martini.itp"//g' topol.top
cat << EOF > topol.add
#include "$ffdir/martini_v2.2.itp"
#include "$ffdir/martini_v2.0_ions.itp"
#include "$ffdir/martini_v2.0_lipids_all_201506.itp"
EOF
cat topol.add topol.top > tmp
rm topol.add
mv tmp topol.top

# minimization
gmx grompp -f min.mdp -c bilayer.gro -p topol.top -o min.tpr -quiet
gmx mdrun -deffnm min -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# equilibration
cat << EOF > index.input
2
name 13 LIP
3 | 4 | 7
name 14 SOL
q
EOF
gmx make_ndx -f bilayer.gro -quiet < index.input
rm index.input
gmx grompp -f eq.mdp -c min.gro -r min.gro -p topol.top -o eq.tpr -n index.ndx -quiet
gmx mdrun -deffnm eq -v -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# encurv
gmx grompp -f md.mdp -p -c eq.gro -n -quiet -o md.tpr
gmx mdrun -v -deffnm md -plumed plumed3.dat -cpi md.cpt -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -quiet

# vizualize
echo 2 0 | gmx trjconv -f md.xtc -s md.tpr -o md_pbc.xtc -center -pbc whole -dt 200 -quiet
cat << EOF > md.pml
load md.gro
load_traj md_pbc.xtc
remove resname W resname WF
show spheres
EOF
pymol md.pml
rm md.mpl

## clean up
rm \#*

