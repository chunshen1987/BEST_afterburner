#!/usr/bin/env bash

# 1) Download the SMASH code
git clone --depth=1 https://github.com/smash-transport/smash.git

# 2) Get the right version of Pythia right here (newer SMASH may need a different version)
(
  wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8303.tgz
  tar xf pythia8303.tgz && rm pythia8303.tgz
  cd pythia8303
  ./configure --cxx-common='-std=c++11 -march=native -mfpmath=sse -O3 -fPIC'
  make
)

# 3) Compile SMASH
(
cd smash && mkdir build && cd build
cmake .. -DPythia_CONFIG_EXECUTABLE=../../pythia8303/bin/pythia8-config -DUSE_ROOT=OFF

machine="$(uname -s)"
case "${machine}" in
    Linux*)     number_of_cores=`nproc --all`;;
    Darwin*)    number_of_cores=`sysctl -n hw.ncpu`;;
    *)          number_of_cores=1;;
esac
number_of_cores_to_compile=$(( ${number_of_cores} > 20 ? 20 : ${number_of_cores} ))
echo "Compiling SMASH using ${number_of_cores_to_compile} cores."
make -j${number_of_cores} smash smash_shared
)

# 4) Get smash analysis scripts
git clone --depth=100 https://github.com/smash-transport/smash-analysis.git
