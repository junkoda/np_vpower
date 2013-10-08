np_vpower
=========

`np_vpower` calculates velocity power spectrum of cosmological
large-scale structure using the Nearest-particle method. The code
reads N-body particles or haloes from a simulation in a periodic box.

## Getting Started

1. The code uses the FFTW3 library. Set the `FFTW_DIR` variable in the
`Makefile`.

2. The code also uses the `BOOST` library to read command line
options, which is optional. Set the `BOOST_DIR` variable if necessary,
or you can disable using comandline options by setting
`USE_BOOST_PROGRAM_OPTIONS = 0`.

3. Make

       $ make

4. This package also includes a code to calculate momentum power spectrum.

       $ make momentum_power


## Authors
This code is written by Jun Koda.

1. Koda et al. 2013, *Are peculiar velocity surveys competitive as a
cosmological probe?*, in preparation.

Please reference this paper (when available) if you use this code for
scientific works.

## Licence
This code is distributed under the GPLv3 license.


