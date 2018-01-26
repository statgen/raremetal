Raremetal: A tool for rare variants meta-analysis

(c) 2012-2017 Shuang Feng, Dajiang Liu, Sai Chen, Gonçalo Abecasis

For more information, please see the [RAREMETAL wiki](http://genome.sph.umich.edu/wiki/RAREMETAL).


**There is a known issue with the newest versions (>=4.14.0) of RAREMETAL, which may give incorrect p-values for some
 burden test calculation methods. We are working on a fix, but for now we recommend using an older version (<= 4.13.9).**


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

If you encounter problems while building, see the [wiki FAQ](https://genome.sph.umich.edu/wiki/RAREMETAL_FAQ) for 
more information.

## Changes
Newest version: v.4.14.1 released at 07/10/2017
View change log at: http://genome.sph.umich.edu/wiki/RAREMETAL_Change_Log

### RAREMETAL2
We are also working on a beta version that provides new features and options- RAREMETAL2. 

This code is still in the experimental stages, but a full release is planned in the future. See:
https://github.com/traxexx/Raremetal2
