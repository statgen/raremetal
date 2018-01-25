#######
# A Dockerfile that demonstrates how to build from scratch in a clean environment.
# This is provided as an aid to development, and is not recommended as a means to deploy the binary for production use.
#######

# CSG clusters may run an older LTS version
FROM ubuntu:16.04
LABEL maintainer="University of Michigan Center for Statistical Genetics"
WORKDIR /code
COPY . /code

# TODO: Vim and tabix can be safely removed once shared directories are mounted
#  (for now they are provided as workflow tools to help verify whether the binary works as expected)
RUN apt-get update && \
    apt-get -y install apt-utils cmake build-essential gfortran \
        gsl-bin libgsl-dev zlib1g-dev \
        vim tabix

RUN apt-get -y install python python-pip && \
    pip install cget

# Build the binary once when container is first created.
RUN cget install -f requirements.txt
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake ..  && \
    make

# Eventually, we will expect to pass arguments. For now start a shell and keep it running,
#  so user can connect to the shell inside the container for their work
CMD /bin/bash