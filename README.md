# ImpactX

[![CI Status](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ECP-WarpX/impactx/actions/workflows/ubuntu.yml)  
[![Supported Platforms](https://img.shields.io/badge/platforms-linux%20|%20osx%20|%20win-blue)](https://impactx.readthedocs.io/en/latest/install/users.html)
[![GitHub commits since last release](https://img.shields.io/github/commits-since/ECP-WarpX/impactx/latest/development.svg)](https://github.com/ECP-WarpX/impactx/compare/development)
[![Development Status](https://img.shields.io/badge/development%20status-pre--alpha-orange.svg)]()
[![Language: C++14](https://img.shields.io/badge/language-C%2B%2B14-orange.svg)](https://isocpp.org/)
[![Language: Python](https://img.shields.io/badge/language-Python-orange.svg)](https://python.org/)  
[![License ImpactX](https://img.shields.io/badge/license-BSD--3--Clause--LBNL-blue.svg)](https://spdx.org/licenses/BSD-3-Clause-LBNL.html)

ImpactX: the next generation of the [IMPACT-Z](https://github.com/impact-lbl/IMPACT-Z) code

## Documentation

TODO: In order to learn how to install and run the code, please see the online documentation:
https://impactx.readthedocs.io

* AMReX Doxygen: https://amrex-codes.github.io/amrex/doxygen/
* ImpactX Doxygen: ... TODO
* WarpX Doxygen: ... TODO

## Contributing

[![AMReX](https://img.shields.io/static/v1?label="runs%20on"&message="AMReX"&color="blueviolet")](https://amrex-codes.github.io/)

Our workflow is described in [CONTRIBUTING.rst](CONTRIBUTING.rst).

## Developer Environment

Please prepare you local development environment as follows.
Pick *one* of the methods below:

### Perlmutter (NERSC)

### Cori (NERSC)

```bash
```

### Homebrew (macOS)

```bash
brew update
brew install adios2      # for openPMD
brew install ccache
brew install cmake
brew install fftw
brew install git
brew install hdf5-mpi    # for openPMD
brew install libomp      # for OpenMP
brew install pkg-config  # for fftw
brew install open-mpi
```

### Apt (Debian/Ubuntu)

```bash
sudo apt update
sudo apt install build-essential ccache cmake g++ git libfftw3-mpi-dev libfftw3-dev libhdf5-openmpi-dev libopenmpi-dev pkg-config python3 python3-matplotlib python3-numpy python3-scipy
```

### Spack (Linux)

```bash
spack env create impactx-dev
spack env activate impactx-dev
spack add adios2        # for openPMD
spack add ccache
spack add cmake
spack add fftw
spack add hdf5          # for openPMD
spack add mpi
spack add pkgconfig     # for fftw
# optional:
# spack add cuda
spack install
```

(in new terminals, re-activate the environment with `spack env activate impactx-dev` again)

### Conda (Windows)

```bash
conda create -n impactx-dev -c conda-forge adios2 ccache cmake compilers git hdf5 fftw matplotlib ninja
conda activate impactx-dev

# compile with -DImpactX_MPI=OFF
```

## Get the Source Code

Before you start, you will need a copy of the ImpactX source code:

```bash
git clone git@github.com:ECP-WarpX/impactx.git
cd impactx
```

## Compile

```bash
# find dependencies & configure
cmake -S . -B build

# compile
cmake --build build -j 4
```

That's all!
ImpactX binaries are now in `build/bin/`.
Most people execute these binaries directly or copy them out.

You can inspect and modify build options after running `cmake -S . -B` build with either

```bash
ccmake build
```

or by adding arguments with `-D<OPTION>=<VALUE>` to the first CMake call, e.g.:

```bash
cmake -S . -B build -DImpactX_COMPUTE=CUDA -DImpactX_MPI=OFF
```

## Run

An executable ImpactX binary with the current compile-time options encoded in its file name will be created in `build/bin/`.

Additionally, a symbolic link named `impactx` can be found in that directory, which points to the last built ImpactX executable.

## Acknowledgements

TODO: LDRD et al.

## License

TODO
