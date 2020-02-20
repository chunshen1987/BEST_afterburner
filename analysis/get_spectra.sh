#!/bin/bash

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

python_scripts_folder=${SCRIPTPATH}/../external_codes/smash-analysis/test/energy_scan
spectra_output_folder=smash_analyzed
mkdir -p ${spectra_output_folder}

SMASH_results_folder=./smash_output

python ${python_scripts_folder}/mult_and_spectra.py \
       --output_files ${spectra_output_folder}/yspectra.txt \
                      ${spectra_output_folder}/mtspectra.txt \
                      ${spectra_output_folder}/ptspectra.txt \
                      ${spectra_output_folder}/v2spectra.txt \
                      ${spectra_output_folder}/meanmt0_midrapidity.txt \
                      ${spectra_output_folder}/meanpt_midrapidity.txt \
                      ${spectra_output_folder}/midrapidity_yield.txt \
                      ${spectra_output_folder}/total_multiplicity.txt \
       --input_files  ${SMASH_results_folder}/*/particles_binary.bin*
#       --parallel true
