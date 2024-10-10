# pNMRsim
An NMR simulation tool optimised for dense dipolar coupled networks

For pre-compiled binaries (strongly recommended), do refer to the [this webpage](https://www.durham.ac.uk/departments/academic/chemistry/about-us/solid-state-nmr/research-profile/equipment-and-software/pnmrsim/).
The links in this page are sadly broken, but a zip file containing all the files can be downloaded from [this page](https://www.dur.ac.uk/paul.hodgkinson/pNMRsim).

This repository is the bleeding edge source code, currently a maintenance version of V18.05.13.

## Introduction

pNMRsim is a [SIMPSON](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson)-like simulation program which is suited to the simulation of extended systems of spin-1/2 nuclei - the p refers to the use of 'periodic' geometries to introduce additional symmetry and block diagonalisation. pNMRsim tries hard to only recalculate quantities when strictly necessary, and historically at least will generally be faster than the equivalent SIMPSON simulation, but this added complexity does increase the chances of bugs surfacing in odd situations.  In short, use at your own risk and always check that initial results make sense before relying on the output.

## Requirements

pNMRsim is dependent on the following components:

A moderatately up-to-date C++ compiler (supporting C++11).  The code is only tested only GCC compilers under Linux, and some tweaks may be required for other compilers / run-times.

The [`libcmatrix` NMR library](https://github.com/dch0ph/libcmatrix).  pNMRsim is in effect a "user interface" to these underlying library routines.  The library is updated as required by pNMRsim and so it is advisable to compile an up-to-date of the library when compiling pNMRsim. 

The `Minuit` library (in its C++ incarnation) is required for optimisation problems.  It is not required, however, for data fitting which uses `libcmatrix` routines.  The "version 2" form of the C++ Minuit library is recommended as more up to date, but it must be a stand-alone version rather than one integrated into ROOT (see the README for `libcmatrix`).

The "classic" version of the Boost Spirit library is used for parsing. You will need the standalone version of the final incarnation (1.8.5), not the version that is currently supplied with Boost. A slightly patched version is available at the pNMRsim web page.
 
[SIMPLOT](https://inano.au.dk/about/research-centers-and-projects/nmr/software/simpson) is useful for visualising pNMRsim output (although a variety of output formats are supported).


## Installation

It is most efficient to compile `libcmatrix` and pNMRsim one after the other since the same configuration options are involved.

The options to configure are:
```
--with-atlas, --with-acml, --with-mkl, --with-openblas  enable use of one of the ATLAS, ACML, MKL or OpenBLAS libraries 
--with-MPI  enable the use of MPI for parallel computations
--with-minuit  enable the use of Minuit
--with-sse  enable the use of SSE3 instructions for complex arithmetic 
```

See the installation instructions for `libcmatrix` for more details.


### pNMRsim specific options

`--enable-static`  Use static linkage (default off). This is recommended when compiling for Windows. It will generate stand-alone (but rather large) .exe files.

`--disable-CXX11`  Disable C++11 features (default enabled). Up-to-date compilers should be happy with the minimal use of C++11 features. Note that the flags used to enable C++11 features may need tweaking depending on compiler (uses `-std=c++11`)


1.  Ensure `libcmatrix` is compiled and (optionally) installed.

2.  Configure pNMRsim as above, making sure that `configure` can find the `libcmatrix` header files and library (`CPPFLAGS` and `LDFLAGS` respectively - again the `libcmatrix` INSTALL instructions discuss this).  You will need to add the path to the Spirit files by adding `-I`<Spirit directory (containing boost/spirit)>.

Hence a minimal configure run will typically be
`CPPFLAGS="-I../spirit-1.8.5-miniboost -I../libcmatrix/include" LDFLAGS="-L../libcmatrix/lib" ./configure`

3.  `make` will compile pNMRsim and pNMRproc.

These is currently an issue with configure failing to find `std::isnan` on Windows. If compile errors involving `isnan` are raised, try adding `-DISNAN_IN_STD` to `CPPFLAGS`. In general, compilation for Linux is strongly recommended.


### Other files
The `test` directory contains a varied and rather chaotic collection of input files.  Many have been hacked around for testing purposes so check that the input file actually does what suggests it does...
