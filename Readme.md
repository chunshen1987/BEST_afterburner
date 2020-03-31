# BEST Afterburner

This is a wrapper, that combines different particle samplers and a hadronic transport package [SMASH](https://smash-transport.github.io).


## Install

* Get and compile the [SMASH transport package](https://smash-transport.github.io)
```
    cd external_codes
    ./get_smash.sh
```

If compilation went wrong, take a look at the cmake messages. Likely gsl, boost or Pythia is not found. Install these SMASH prequisites, see smash/Readme.md for details.

* Get and compile the BEST particle sampler

```
    cd external_codes
    ./get_best_sampler.sh
```

* Get and compile the microcanonical particle sampler

```
    cd external_codes
    ./get_microcanonical_sampler.sh
```


* Compile the wrapper for particle sampler and SMASH

```
    mkdir build
    cd build
    cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
    make
```

The commands above are compiled in a bash script `compile.sh`. One can run this script to build the wrapper package.

If one wants to run the iSS code package to sample particles and feed into SMASH, one can run the provided bash script, `./compile_with_iSS.sh`.

* Configuring and running

The configuration for ALL codes (SMASH, and all samplers) is in the config.yaml file.
For explanations of parameters see the help of specific codes.
