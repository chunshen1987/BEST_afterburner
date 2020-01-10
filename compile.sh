#!/usr/bin/env bash

mkdir -p build
(
    cd build
    rm -fr *
    cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
    make
)
