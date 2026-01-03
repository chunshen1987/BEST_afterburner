#!/usr/bin/env bash

[ -d "smash" ] && rm -fr smash

machine="$(uname -s)"
case "${machine}" in
    Linux*)     number_of_cores=`nproc --all`;;
    Darwin*)    number_of_cores=`sysctl -n hw.ncpu`;;
    *)          number_of_cores=1;;
esac
number_of_cores_to_compile=$(( ${number_of_cores} > 4 ? 4 : ${number_of_cores} ))

# 1) Download the SMASH code
git clone -b SMASH-3.3 --depth=1 https://github.com/smash-transport/smash.git

# 2) Get the right version of Pythia right here (newer SMASH may need a different version)
(
    wget https://pythia.org/download/pythia83/pythia8316.tar.bz2
    tar xf pythia8316.tar.bz2 && rm pythia8316.tar.bz2
    cd pythia8316
    ./configure --cxx-common='-std=c++11 -O3 -fPIC'

    echo "Compiling PYTHIA using ${number_of_cores_to_compile} cores."
    make -j${number_of_cores_to_compile}
)

# 3) Compile SMASH
(
    cd smash
    rm cmake/FindGSL.cmake
    mkdir -p build && cd build
    cmake .. -DPythia_CONFIG_EXECUTABLE=../../pythia8316/bin/pythia8-config -DUSE_ROOT=OFF

    echo "Compiling SMASH using ${number_of_cores_to_compile} cores."
    make -j${number_of_cores} smash smash_shared
)

# 4) Get smash analysis scripts
#git clone --depth=1 https://github.com/smash-transport/smash-analysis.git
