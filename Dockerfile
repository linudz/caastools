# Deriving the latest RERConverge docker module
FROM wem26/rerconverge:latest

# Set session as noninteractive and install/upgrade common debian packages
# Set Python==3.9 for repro
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libpq-dev python3.9 python3-pip python3-setuptools python3-dev vim
RUN pip3 install --upgrade pip

# installing python libraries (Discovery/Bootstrap/Resample)
RUN pip3 install -r requirements.txt

# installing r libraries (already installed from module)
# RUN Rscript requirements.r





