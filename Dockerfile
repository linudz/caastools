# Deriving the R-base 4.3.1 for repro
# Maybe use a more mainstream container
# maaybe build something that already has devtools and fuck this
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
ADD requirements/requirements.txt ./requirements/
ADD requirements/requirements.r ./requirements/
ADD modules/ ./modules/
ADD scripts/ ./scripts/
ADD ct .
RUN chmod +x ./ct


# Installing Python libraries (Discovery/Bootstrap/Resample)
RUN pip3 install -r requirements/requirements.txt

# Installing R libraries (BM Resampling)
RUN Rscript requirements/requirements.r
# think about not automatizing this