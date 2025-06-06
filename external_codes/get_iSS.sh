#!/usr/bin/env bash

[ -d "iSS" ] && rm -fr iSS

# download the code package
git clone --depth=1 https://github.com/chunshen1987/iSS -b dev iSS
