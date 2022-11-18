# Wrapper for samplers + SMASH

The wrapper allows to run samplers and [SMASH](https://smash-transport.github.io) hadronic transport in a consistent manner under one configuration.

To download this repository run

```
git clone https://doliinychenko@bitbucket.org/bestcollaboration/best_afterburner.git
```

## Installing

1. Download the 3rd party software

```
cd external_codes
./get_msu_sampler.sh
./get_iSS.sh
./get_microcanonical_sampler.sh
./get_smash.sh
cd ..
```

  Do this one by one and make sure that it is actually downloaded and installed.
  If compilation went wrong, take a look at the cmake messages.
  Likely gsl, boost or Pythia is not found. Install these SMASH prequisites, see smash/Readme.md for details.

2. Compile the wrapper

```
mkdir build && cd build
cmake .. -DPythia_CONFIG_EXECUTABLE=../external_codes/pythia8307/bin/pythia8-config
make
```

  Troubleshooting: pay attention at cmake error messages.

## Configuring and running

  If you just run ./sampler_and_smash you will get an error: "Config file ../config.yaml not found".
  The wrapper need a configuration file, which by default is ../config.yaml.

  The configuration for ALL codes (SMASH, and all samplers) should be in the config.yaml file.
  An example of configuration file can be found at
  analysis/config_samplers_comparison.yaml

  Sampler can be run with the command-line options like

```
./sampler_and_smash -h                // This print help
./sampler_and_smash -c myconfig.yaml  // This runs sampler with configuration file myconfig.yaml
```

## Analysis of SMASH results

  A considerable amount of analysis tools is already developed by the SMASH team.
  For an example of analysis using these tools look ata analysis/get_spectra.sh
