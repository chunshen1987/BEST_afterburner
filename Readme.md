# BEST Afterburner

This is a code repostiory for a wrapper between individual particle sampler and hadronic transport package [SMASH](https://smash-transport.github.io). 


## Install

1. Get and compile the [SMASH transport package](https://smash-transport.github.io)

```
    cd external_codes
    ./get_smash.sh
```

If compilation went wrong, take a look at the cmake messages. Likely gsl, boost or Pythia is not found.
Install these SMASH prequisites, see smash/Readme.md for details.

2. Get and compile the BEST particle sampler

```
    cd external_codes
    ./get_best_sampler.sh
```

3. Compile the wrapper for particle sampler and SMASH

```
    mkdir build
    cd build
    cmake .. -DPythia_CONFIG_EXECUTABLE=${PYTHIA8DIR}/bin/pythia8-config
    make
```

