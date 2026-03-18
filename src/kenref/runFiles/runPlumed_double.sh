#!/bin/bash
# Updated to use GROMACS-4-PLUMED installation
source /smithlab/opt/gromacs-4-plumed/2025/2025.3/debug/AVX_512/bin/GMXRC
source /home/amr/git/plumed2/sourceme.sh

# Verify PLUMED is loaded
if [ -z "$PLUMED_KERNEL" ]; then
    echo "ERROR: PLUMED_KERNEL not set. Please check PLUMED installation."
    exit 1
fi
echo "PLUMED_KERNEL=$PLUMED_KERNEL"

#When you run multi simulations, you need to have the file in a relative folder
mpirun -n 2 gmx_mpi mdrun -s topol.tpr -multidir repl_01 repl_02 -nsteps 50 -plumed ../../runFiles/test_double.dat # -deffnm $fnm &>> $outname