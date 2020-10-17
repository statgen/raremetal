RAREMETAL: A tool for rare variants meta-analysis

(c) 2012-2018 Shuang Feng, Dajiang Liu, Sai Chen, Gonçalo Abecasis

For more information, please see the [RAREMETAL wiki](http://genome.sph.umich.edu/wiki/RAREMETAL).

## Installation

### Requirements

RAREMETAL requires the following libraries to build. These should already be installed on most scientific 
computing clusters.

#### Linux

On Linux, the following packages are needed:

- zlib
- liblzma
- cmake
- gcc/g++
- gfortran
- python3
- cget

As an example in Ubuntu specifically, you can install with the following commands:

```bash
sudo apt update
sudo apt install build-essential cmake gfortran locales zlib1g-dev liblzma-dev
sudo apt install python3 python3-setuptools python3-pip python3-virtualenv
sudo python3 -m pip install cget
```

#### MacOS

Use Homebrew to install the necessary packages:

```bash
brew install cmake zlib xz python3
brew cask install gfortran
sudo python3 -m pip install cget
```

You may need to manually install the C++ headers on MacOS 10.4 (Mojave) by running the following command: 

```bash
sudo installer -pkg /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg -target /
```

For Mac OS 10.5 (Catalina), C++ headers may be missing and you may receive excessive compiler warnings. The following settings may help:

```bash
# Specify path to MacOS C++ SDK header files
export CPATH=$(xcrun --show-sdk-path)/usr/include

# Disable compiler warnings introduced in Catalina
export CFLAGS="${CFLAGS} -Wno-nullability-completeness -Wno-expansion-to-defined -Wno-undef"
```

#### Windows

Windows is currently unsupported. You may have some luck with Windows Subsystem for Linux, or a VM.

### Installing

1. Clone this repository: `git clone https://github.com/statgen/raremetal.git`
2. Open a terminal and change to the directory containing this README file.
3. Enter the following sequence of commands to download dependencies and build the code:

    ```bash
    cget install -f requirements.txt
    mkdir build && cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DBUILD_TESTS=1 ..
    make
    ```

If you encounter problems while building, see the [wiki instructions](https://genome.sph.umich.edu/wiki/RAREMETAL_DOWNLOAD_%26_BUILD) 
and [FAQ](https://genome.sph.umich.edu/wiki/RAREMETAL_FAQ) for more information.

### Running tests

To run unit tests, ensure that you have run cmake with `-DBUILD_TESTS=1` and built at least once, then run: 

```bash
# To test RAREMETAL, run the following:
build/raremetal/tests/bin/tests-raremetal

# To test RAREMETALWORKER, run the following:
cd build
make test
```

Both sets of tests are automatically run by Travis CI as well, see `.travis.yml` for details.

## Changes

Newest version: v.4.14.1 released at 07/10/2017

View change log at: http://genome.sph.umich.edu/wiki/RAREMETAL_Change_Log

### RAREMETAL2

We are also working on a beta version that provides new features and options- RAREMETAL2. 

This code is still in the experimental stages, but a full release is planned in the future. See:
https://github.com/statgen/Raremetal2
