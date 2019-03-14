Raremetal: A tool for rare variants meta-analysis

(c) 2012-2019 Shuang Feng, Dajiang Liu, Sai Chen, Gonçalo Abecasis

For more information, please see the [RAREMETAL wiki](http://genome.sph.umich.edu/wiki/RAREMETAL).

## Installation
1. Clone this repository
2. Open a terminal and change to the directory containing this README file.
3. Enter the following sequence of commands to download dependencies and build the code:
 
```bash
$ cget install -f requirements.txt
$ mkdir build && cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DBUILD_TESTS=1 ..
$ make
```

To run unit tests, ensure that you have run cmake with `-DBUILD_TESTS=1` and built at least once, 
then run `make` or `ctest --verbose`. These unit tests were written on Mac OS and you may encounter small differences 
due to system precision on other platforms.

If you encounter problems while building, see the [wiki instructions](https://genome.sph.umich.edu/wiki/RAREMETAL_DOWNLOAD_%26_BUILD) 
and [FAQ](https://genome.sph.umich.edu/wiki/RAREMETAL_FAQ) for more information.

### Requirements
Raremetal requires the following libraries to build. These should already be installed on most scientific 
computing clusters.

- libRMath (The R math library. If not installed, use `apt-get install r-mathlib` on ubuntu or `brew install r` 
    on Mac OS)
- zlib (`apt-get install zlib1g-dev` on Ubuntu)

### Troubleshooting
On Mac OS, installing R as a standalone package is not sufficient to use libRMath in other executables. 
Consider `brew install r` or [r-devel](https://r.research.att.com/) in that case.   

On Mac OS, some `cget` dependencies (such as xz) may fail to install. In certain cases, this can be fixed by running an 
OS-version-specific command similar to the below:

`$ open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg`

## Changes
Newest version: v4.15.0 released March 14, 2019.
View change log at: http://genome.sph.umich.edu/wiki/RAREMETAL_Change_Log

### RAREMETAL2
We are also working on a beta version that provides new features and options- RAREMETAL2. 

This code is still in the experimental stages, but a full release is planned in the future. See:
https://github.com/statgen/Raremetal2
