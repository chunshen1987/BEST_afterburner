To install the code:


1) Get and compile SMASH

  cd external_codes
  ./get_smash

If compilation went wrong, look and cmake messages. Likely gsl, boost or Pythia is not found.
Install these SMASH prequisites, see smash/Readme.md for details.

2) Get and compile the sampler

  ???

3) Compile the wrapper

mkdir build
cd build
cmake ..
make
