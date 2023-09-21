# Deriving the R-base 4.3.1 for repro
# Maybe use a more mainstream container
FROM thomaschln/r-devtools:latest

# Set session as noninteractive and install/upgrade common debian packages
# Set Python==3.9 for repro
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -qq -y --no-install-recommends \ 
        build-essential \
        libpq-dev \
        python3.9 \
        python3-pip \
        python3-setuptools \
        python3-dev \
        vim

RUN pip3 install --upgrade pip

WORKDIR /app
ADD requirements/requirements.txt .
ADD requirements/requirements.r .

# Installing Python libraries (Discovery/Bootstrap/Resample)
RUN pip3 install -r requirements.txt

# Installing R libraries (BM Resampling)
RUN Rscript requirements.r