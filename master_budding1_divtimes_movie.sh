#!/bin/bash

dir=/global/homes/c/cschreck/production/cavity_growth_CG

Lx=6.0
Ly=16.0

P0=-1
att=0.0

cell=2
div=1

#rate=0.002
#steps=20000
rate=0.002
steps=4000
skip=4
seed=1

a=2.0
sbatch run_slurm_divtime_movie.sh $Lx $Ly $a $cell $div $P0 $att $rate $steps $seed $skip
