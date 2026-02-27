cmake_minimum_required(VERSION 3.10)
project(PLUMED_IDE_HELPER)
set(CMAKE_CXX_STANDARD 17)

### comment the next line to make it float. Uncomment to make it double.
add_compile_options(-DKENREF_DOUBLE)

# 1. Define PLUMED include paths
# (These paths mimic what 'configure' found. Update paths if they differ on your system)
include_directories(
        src
#        src/lib
#        src/core
#        src/bias
#        src/cltools
#        src/wrapper
#        src/kenref
        # Add external paths found by configure (copy from Makefile.conf if needed)
        /smithlab/opt/fftw/3.3.10/include
        /smithlab/opt/openmpi/5.0.7_clang20.1.1_cuda12.4.1/include
)

# 2. Add your KEnRef include paths manually
# (Copy the content of KENREF_INCLUDE from your Makefile.conf, removing the -I flags)
include_directories(
        /smithlab/opt/kenref/debug/AVX_512/include/core
        /smithlab/opt/eigen/3.4.0/include/eigen3
)

# 3. Add source files
# Recursively find all C++ files so CLion indexes them
file(GLOB_RECURSE SOURCES "src/*/*.cpp" "src/*/*.h")

# 4. Create a fake target
add_executable(plumed_dummy ${SOURCES}
        src/kenref/config.h)

add_executable(nothing donothing.cpp src/kenref/config.h)

# 5. Define Macros used by PLUMED (Optional but helps)
add_definitions(-D__PLUMED_HAS_GROMACS -D__PLUMED_HAS_MPI -DKENREF_DOUBLE)