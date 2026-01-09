#!/usr/bin/env bash

machine="$(uname -s)"
case "${machine}" in
    Linux*)     number_of_cores=`nproc --all`;;
    Darwin*)    number_of_cores=`sysctl -n hw.ncpu`;;
    *)          number_of_cores=1;;
esac
number_of_cores_to_compile=$(( ${number_of_cores} > 4 ? 4 : ${number_of_cores} ))

(
    cd external_codes
    ./get_iSS.sh
)

mkdir -p build
(
    cd build
    rm -fr *
    cmake .. -DPythia_CONFIG_EXECUTABLE=../external_codes/pythia8315/bin/pythia8-config -DiSS=ON
    make -j${number_of_cores}
)
