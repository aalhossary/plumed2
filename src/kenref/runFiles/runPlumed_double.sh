source /smithlab/opt/gromacs-dev/2025/2025.3/debug/AVX_512/bin/GMXRC
source /home/amr/git/plumed2/sourceme.sh
#When you run multi simulations, you need to have the file in a relative folder
mpirun -n 2 gmx_mpi mdrun -s topol.tpr -multidir repl_01 repl_02 -nsteps 50 -plumed ../../runFiles/test_double.dat # -deffnm $fnm &>> $outname