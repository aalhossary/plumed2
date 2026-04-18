#!/bin/bash
# Updated to use GROMACS-4-PLUMED installation
source /smithlab/opt/gromacs-4-plumed/2025/2025.3/relwithdebinfo-asan/AVX_512/bin/GMXRC
#source /smithlab/opt/gromacs-4-plumed/2025/2025.3/debug/AVX_512/bin/GMXRC

#source /home/amr/git/plumed2/sourceme.sh
#export PLUMED_KERNEL=/smithlab/opt/plumed-dev/master/asan/AVX_512/lib/libplumedKernel.so
export PLUMED_KERNEL=/smithlab/opt/plumed-dev/master/relwithdebinfo-asan/AVX_512/lib/libplumedKernel.so
#export PLUMED_KERNEL=/smithlab/opt/plumed-dev/master/debug/AVX_512/lib/libplumedKernel.so

# Verify PLUMED is loaded
if [ -z "$PLUMED_KERNEL" ]; then
    echo "ERROR: PLUMED_KERNEL not set. Please check PLUMED installation."
    exit 1
fi
echo "PLUMED_KERNEL=$PLUMED_KERNEL"

#export ASAN_OPTIONS="detect_leaks=0:log_path=./asan.log:abort_on_error=1"
#export LD_PRELOAD=/smithlab/opt/llvm/20.1.1/lib/clang/20/lib/x86_64-unknown-linux-gnu/libclang_rt.asan.so

#When you run multi simulations, you need to have the file in a relative folder
#gmx_mpi mdrun -s /home/amr/CLionProjects/KEnRef/res/run-output/repl_01/topol.tpr -nsteps 50 -plumed ../runFiles/ubiquitin_1e8_.25_999_A.dat

#valgrind --tool=memcheck \
#    --error-exitcode=1 \
#    --leak-check=no \
#    --track-origins=yes \
#    --num-callers=30 \
#    --log-file=valgrind.%p.log \
#    gmx_mpi mdrun -cpi -pin on -ntomp 1 -cpt 30 \
#        -plumed ../runFiles/ubiquitin_1e8_.25_999_A.dat \
#        -deffnm test \
#        -s /smithlab/home/aalhossary/ubiquitin-plateaus-plumed/10nsstart+fitting/alef.tpr

export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1:log_path=asan.log
#export LD_PRELOAD=/smithlab/opt/llvm/20.1.1/lib/clang/20/lib/x86_64-unknown-linux-gnu/libclang_rt.asan.so
#PLUMED_KERNEL=/smithlab/opt/plumed-dev/master/asan/AVX_512/lib/libplumedKernel.so \
export PLUMED_LOAD_NODEEPBIND=1
gmx_mpi mdrun -cpi -pin on -ntomp 20 -cpt 30 \
    -plumed ../runFiles/ubiquitin_1e8_.25_999_A.dat \
    -deffnm test \
    -s /smithlab/home/aalhossary/ubiquitin-plateaus-plumed/10nsstart+fitting/alef.tpr


##cd /home/amr/git/plumed2/src/kenref/run-output
#PLUMED_KERNEL=/smithlab/opt/plumed-dev/master/relwithdebinfo/AVX_512/lib/libplumedKernel.so \
#rr record gmx_mpi mdrun -cpi -pin on -ntomp 1 -cpt 30 \
#    -plumed ../runFiles/ubiquitin_1e8_.25_999_A.dat \
#    -deffnm test \
#    -s /smithlab/home/aalhossary/ubiquitin-plateaus-plumed/10nsstart+fitting/alef.tpr
