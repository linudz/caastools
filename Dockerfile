## Emacs, make this -*- mode: sh; -*-
FROM ubuntu:jammy

## Set a default user. Available via runtime flag `--user docker` 
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (for rstudio or linked volumes to work properly). 
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

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

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Otherwise timedatectl will get called which leads to 'no systemd' inside Docker
ENV TZ UTC

# Set installation as noninteractive and install required python/r dependencies
# is it a good idea to explicit python and r deffinitions? also rerconverge
RUN apt-get update \
        && DEBIAN_FRONTEND=noninteractive \
        apt-get install -y --no-install-recommends \
            python3 \
            python3-pip \
            python3-setuptools \
            r-base \
            r-base-dev \
            r-recommended \
            r-cran-devtools \

        && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
        && rm -rf /var/lib/apt/lists/*

# JIC
RUN pip3 install --upgrade pip

# Set CAASTools WD and build directory structure
WORKDIR /ct
RUN mkdir -p ./requirements ./modules ./scripts
ADD requirements/requirements.txt ./requirements/
ADD requirements/requirements.r ./requirements/
ADD modules/ ./modules/
ADD scripts/ ./scripts/
# Make ct executable and add to $PATH
ADD ct .
RUN chmod +x ./ct
ENV PATH=/ct:$PATH

# Installing Python libraries (Discovery/Bootstrap/Resample)
RUN pip3 install -r requirements/requirements.txt

# Installing R libraries (BM Resampling)
RUN Rscript requirements/requirements.r

CMD ["bash"]