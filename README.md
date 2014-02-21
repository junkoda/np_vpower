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

4. Run. This code reads Gadget binary file

       $ ./np_vpower --gadget-file snp > vpower.txt

where `snp` is the name of the snapshot file. Modify `read_fof_text`
function in `halo_file.cpp` to read your own files of dark matter
particles or haloes.

## Command-line Options
You can change the parameters by setting command-line options

       $ ./np_vpower [options] <filename>

* --nc=128   : Set the number of grids for the Fourier Transform.
* --dk=0.01  : Set dk to change output k binning
* --kmax=1.0 : Set kmax to fix the value of maximum k in the output

The code calculates power spectra in real space by default. Set

* --reddshift-space

for power spectra in redshift space.

The code also calculates momentum power spectrum:

* --momentum-power

Simulation parameters are read from Gadget file, but if you read FoF
halo data from an ascii file, you need to give simulation boxsize,

* --boxsize=1000 : Box size in [/h Mpc]

You also need cosmological parameters for the redshift-space
distortions for FoF data,

* --omegam=0.273

which assumes flat LambdaCDM (omega_lambda= 1 - omegam).

You can select haloes in a mass range, log10(M) in [logMmin, logMmax], using,

* --m=7.5e11       : Particle mass
* --logMmin=12.0   : log10[Mmin/(1/h Solar mass)]
* --logMmax=13.0

See a brief description of options by running the code without arguments,

       $ ./np_vpower


## Authors
This code is written by Jun Koda.

1. Koda et al. 2013, *Are peculiar velocity surveys competitive as a
cosmological probe?*, MNRAS submitted,
[arXiv:1312.1022](http://arxiv.org/abs/1312.1022).

Please reference this paper if you use this code for scientific works.

## Licence
This code is distributed under the GPLv3 license.
