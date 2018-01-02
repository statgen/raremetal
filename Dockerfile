# Our clusters run Ubuntu and so shall we (LTS version may vary)
FROM ubuntu:16.04
LABEL maintainer="University of Michigan Center for Statistical Genetics"
WORKDIR /code
COPY . /code

# TODO: Vim and tabix can be safely removed as dependencies once shared directories are mounted
#  (will move analysis workflow out of the container)
RUN apt-get update && \
    apt-get -y install apt-utils make gcc g++ gfortran \
        gsl-bin libgsl-dev zlib1g-dev \
        vim tabix

# Build the binary once when container is first created.
# FIXME: TEMPORARY: Build dependent libraries first (remove this line + librmath + libstatgen once cmake is integrated)
RUN cd libStatGen/ && make clean && make && cd ../libRareMetal/ && make clean && make && cd ..
RUN make clean && make

# TODO: In future, mount VOLUMEs (eg for mirroring development changes inside the container or capturing results)

# Eventually, we will expect to pass arguments. For now start a shell and keep it running,
#  so user can connect to the shell inside the container for their work
CMD /bin/bash