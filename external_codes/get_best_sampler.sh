#!/usr/bin/env bash

# 1) Download the sampler code
git clone --depth=1 https://github.com/steinhor/best_sampler.git

# 2) Compile
cd best_sampler && mkdir build && cd build
cmake .. && make -j
