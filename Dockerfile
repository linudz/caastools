## Emacs, make this -*- mode: sh; -*-
FROM ubuntu:jammy

# Set session as noninteractive and install required debian packages
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -qq -y --no-install-recommends \ 
        software-properties-common \
        dirmngr \
        ed \
        gpg-agent \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates

# Set installation as noninteractive and install required python/r dependencies
# is it a good idea to explicit python and r deffinitions? also rerconverge
# explore the idea of adding radian for an actually useful R console
RUN apt-get update \
        && DEBIAN_FRONTEND=noninteractive \
        apt-get install -y --no-install-recommends \
            python3 \
            python3-pip \
            python3-setuptools

# JIC
RUN pip3 install --upgrade pip

# Set CAASTools WD and build directory structure
WORKDIR /ct
RUN mkdir -p ./requirements ./modules ./scripts
ADD requirements/requirements.txt ./requirements/
ADD modules/ ./modules/
ADD scripts/ ./scripts/
# Make ct executable and add to $PATH
ADD ct .
RUN chmod +x ./ct
ENV PATH=/ct:$PATH

# Installing Python libraries (Discovery/Bootstrap/Resample)
RUN pip3 install -r requirements/requirements.txt

# NOTE: No R installation is included, use biocontainers/rerconverge package

CMD ["bash"]