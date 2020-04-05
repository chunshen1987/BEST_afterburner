#!/bin/bash

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
analysis_scripts_folder=${SCRIPTPATH}/../external_codes/smash-analysis/test/energy_scan
smash_particles_path=${SCRIPTPATH}/../external_codes/smash/input/particles.txt
smash_decaymodes_path=${SCRIPTPATH}/../external_codes/smash/input/decaymodes.txt

for sampler in MSU iSS Microcanonical
do
  simulation_output_folder=sampler_${sampler}_SMASH_output
  mkdir -p ${simulation_output_folder}
  ./sampler_and_smash -c ${SCRIPTPATH}/config_samplers_comparison.yaml \
                      -o "Particles: ${smash_particles_path}" \
                      -o "DecayModes: ${smash_decaymodes_path}" \
                      -o "Sampler: {Type: ${sampler}}" \
                      -o "Output_Directory: ${simulation_output_folder}" > ${simulation_output_folder}/out.txt &
done
wait

for sampler in MSU iSS Microcanonical
do
  simulation_output_folder=sampler_${sampler}_SMASH_output
  spectra_output_folder=sampler_${sampler}_spectra
  mkdir -p ${spectra_output_folder}
  python ${analysis_scripts_folder}/mult_and_spectra.py \
         --output_files ${spectra_output_folder}/yspectra.txt \
                        ${spectra_output_folder}/mtspectra.txt \
                        ${spectra_output_folder}/ptspectra.txt \
                        ${spectra_output_folder}/v2spectra.txt \
                        ${spectra_output_folder}/meanmt0_midrapidity.txt \
                        ${spectra_output_folder}/meanpt_midrapidity.txt \
                        ${spectra_output_folder}/midrapidity_yield.txt \
                        ${spectra_output_folder}/total_multiplicity.txt \
         --input_files  ${simulation_output_folder}/*/particles_binary.bin*
         --parallel true
done
