#I am not sure this first line is important
#export CC=mpicc; export CXX=mpicxx; export FC=mpif90
cd /home/amr/git/plumed2/ || exit 4
# Updated to use new kenref CORE path structure: /smithlab/opt/kenref/${build_type}/${ACCEL}
export PKG_CONFIG_PATH=/smithlab/opt/kenref/debug/AVX_512/lib/pkgconfig:/smithlab/opt/eigen/3.4.0/share/pkgconfig:$PKG_CONFIG_PATH
autoreconf --force
./configure CXXFLAGS="-stdlib=libc++ -O0 -fPIC -Wall -pedantic -std=c++17 -g" LDFLAGS="-L/smithlab/opt/llvm/20.1.1/lib/x86_64-unknown-linux-gnu"   LIBS="-Wl,--push-state,--no-as-needed -lc++ -lc++abi -lunwind -Wl,--pop-state -lpthread -ldl"   --prefix=/smithlab/opt/plumed-dev/2.10-clang20.1.1   --enable-modules=reset:+kenref  --enable-multi_sim --enable-debug   CXX=mpicxx CC=mpicc FC=mpif90
make clean
make -j18
#sudo make install
