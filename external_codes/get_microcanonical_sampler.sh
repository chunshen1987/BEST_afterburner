#!/usr/bin/env bash

[ -d "microcanonical_cooper_frye" ] && rm -fr microcanonical_cooper_frye

git clone --depth=1 https://github.com/chunshen1987/microcanonical_cooper_frye

cd microcanonical_cooper_frye && mkdir build && cd build
cmake .. && make -j4
