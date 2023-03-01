FROM ubuntu:18.04

ARG BUILD_DATE
ARG GIT_SHA
ARG RAREMETAL_VERSION

LABEL org.label-schema.name="RAREMETAL"
LABEL org.label-schema.description="Meta-analysis of rare genetic variant aggregation tests"
LABEL org.label-schema.vendor="University of Michigan, Center for Statistical Genetics"
LABEL org.label-schema.url="https://github.com/statgen/RAREMETAL"
LABEL org.label-schema.vcs-url="https://github.com/statgen/RAREMETAL"
LABEL org.label-schema.schema-version="1.0"

# Install required packages for swiss to install. Many of swiss' dependencies
# require compiling C/C++ code.
RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  gfortran \
  python3 \
  python3-pip \
  zlib1g-dev \
  liblzma-dev \
  locales \
  && rm -rf /var/lib/apt/lists/* \
  && locale-gen en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Install necessary python packages (backports.lzma needed for cget to extract .xz archives)
RUN pip3 install cget

# Create a group and user to execute as, then drop root
ARG UID
ARG GID
RUN \
  if [ -n "$GID" ]; then \
    addgroup --gid $GID raremetal; \
  else \
    addgroup raremetal; \
  fi && \
  if [ -n "$UID" ]; then \
    adduser --gecos "User for running raremetal as non-root" --shell /bin/bash --disabled-password --uid $UID --ingroup raremetal raremetal; \
  else \
    adduser --gecos "User for running raremetal as non-root" --shell /bin/bash --disabled-password --ingroup raremetal raremetal; \
  fi

WORKDIR /home/raremetal
USER raremetal

# Copy source
COPY --chown=raremetal:raremetal . /home/raremetal/

# Compile & run tests
# We want the docker build to fail if the tests fail
RUN cget install -f requirements.txt
RUN mkdir build && \
  cd build && \
  cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DBUILD_TESTS=1 .. && \
  make && \
  cd .. && \
  build/raremetal/tests/bin/tests-raremetal

ENTRYPOINT ["build/bin/raremetal"]

# Frequently changing metadata here to avoid cache misses
LABEL org.label-schema.version=$RAREMETAL_VERSION
LABEL org.label-schema.vcs-ref=$GIT_SHA
LABEL org.label-schema.build-date=$BUILD_DATE
