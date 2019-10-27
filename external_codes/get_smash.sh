#!/usr/bin/env bash

# 1) Download the SMASH code
git clone --depth=1 https://github.com/smash-transport/smash.git

# 2) Compile SMASH
cd smash && mkdir build && cd build
cmake .. # -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
number_of_cores=`nproc --all`
number_of_cores_to_compile=$(( ${number_of_cores} > 20 ? 20 : ${number_of_cores} ))
echo "Compiling SMASH using ${number_of_cores_to_compile} cores."
make -j${number_of_cores} smash smash_shared
