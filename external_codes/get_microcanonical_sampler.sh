#!/usr/bin/env bash

git clone --depth=1 https://github.com/doliinychenko/microcanonical_cooper_frye.git
cd microcanonical_cooper_frye && mkdir build && cd build
cmake .. && make


