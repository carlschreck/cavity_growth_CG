#!/bin/bash

##SBATCH -p shared
#SBATCH -p debug
#SBATCH -n 1  
##SBATCH -t 48:00:00
#SBATCH -t 00:30:00
#SBATCH -J LRE_VELS

# Set directories
rundir=/global/homes/c/cschreck/production/cavity_growth_CG
outdir=/global/project/projectdirs/m2282/cfs/cavity_growth_CG_divtime

# geometric parameters
ar=1.01
Lx=$1
Ly=$2
AP=$3

# cell parameters
celltype=$4
divtype=$5
P0=$6
att=$7

# run parameters
rate0=$8
desync=0.4
steps=$9
seed=-${10}
distrem=4.5
skip=${11}

# output files
if [ ${celltype} -eq 1 ]
then
cell="ellipse"
elif [ ${celltype} -eq 2 ]
then
cell="budding"
elif [ ${celltype} -eq 3 ]
then
cell="disk"
fi
if [ ${P0%.*} -eq -1 ]
then 
suffix=movie_${cell}_ar${ar}_div${divtype}_seed$((-$seed))_Lx${Lx}_Ly${Ly}_a${AP}_att${att}.dat
else
suffix=movie_${cell}_ar${ar}_div${divtype}_seed$((-$seed))_Lx${Lx}_Ly${Ly}_a${AP}_att${att}_P${P0}.dat
fi
prodfile=prod_$suffix
pressfile=press_$suffix
phifile=phi_$suffix
flowfile=flow_$suffix
treefile=clones_$suffix
agefile=age_$suffix
divfile=div_$suffix

cd $outdir

# run program
time $rundir/divide_cavity_physbox_combo.o <<EOF
  $ar
  $Lx
  $Ly
  $AP
  $celltype
  $divtype
  $P0
  $att
  $rate0
  $desync
  $steps
  $seed
  $distrem
  $skip
  $prodfile
  $pressfile
  $phifile
  $flowfile
  $treefile
  $agefile
  $divfile
EOF
