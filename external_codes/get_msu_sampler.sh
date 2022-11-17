#!/usr/bin/env bash

[ -d "msu_sampler" ] && rm -fr msu_sampler

# 1) Download the sampler code
git clone --depth=1 https://github.com/chunshen1987/msu_sampler

# 2) Compile
cd best_sampler && mkdir build && cd build
cmake .. && make -j4
