source /smithlab/opt/gromacs-dev/2025/2025.3/debug/AVX_512/bin/GMXRC
source /home/amr/git/plumed2/sourceme.sh
#When you run multi simulations, you need to have the file in a relative folder
gmx_mpi mdrun -s /home/amr/CLionProjects/KEnRef/res/run-output/repl_01/topol.tpr -nsteps 50 -plumed ../runFiles/test.dat