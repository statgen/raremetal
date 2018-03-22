Raremetal: A tool for rare variants meta-analysis

(c) 2012-2018 Shuang Feng, Dajiang Liu, Sai Chen, Gonçalo Abecasis

For more information, please see the [RAREMETAL wiki](http://genome.sph.umich.edu/wiki/RAREMETAL).

## Installation
1. Clone this repository
2. Open a terminal and change to the directory containing this README file.
3. Enter the following sequence of commands to download dependencies and build the code (with tests run automatically during build):
 
```bash
cget install -f requirements.txt
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DBUILD_TESTS=1 ..
make
```

To run unit tests, ensure that you have run cmake with `-DBUILD_TESTS=1` and built at least once, 
then run `make` or `ctest --verbose`. 

If you encounter problems while building, see the [wiki instructions](https://genome.sph.umich.edu/wiki/RAREMETAL_DOWNLOAD_%26_BUILD) 
and [FAQ](https://genome.sph.umich.edu/wiki/RAREMETAL_FAQ) for more information.

### Requirements
Raremetal requires the following libraries to build. These should already be installed on most scientific 
computing clusters.

- libRMath (The R math library. If not installed, use `apt-get install r-mathlib` on ubuntu or `brew install r` 
    on Mac OS)
- zlib (`apt-get install zlib1g-dev` on Ubuntu)

On Mac OS, installing R as a standalone package is not sufficient to use libRMath in other executables. 
Consider `brew install r` or [r-devel](https://r.research.att.com/) in that case.   

## Changes
Newest version: v.4.14.1 released at 07/10/2017
View change log at: http://genome.sph.umich.edu/wiki/RAREMETAL_Change_Log

### RAREMETAL2
We are also working on a beta version that provides new features and options- RAREMETAL2. 

This code is still in the experimental stages, but a full release is planned in the future. See:
https://github.com/statgen/Raremetal2
