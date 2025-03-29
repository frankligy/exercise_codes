# syntax=docker/dockerfile:1

FROM debian:buster-slim
USER root


RUN apt-get update \
    && apt-get install -y wget bzip2 parallel git perl build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# install the main python environment including ncbi orffinder
RUN wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh -O ~/anaconda.sh \
    && /bin/bash ~/anaconda.sh -b -p /opt/conda \
    && rm ~/anaconda.sh

ENV PATH /opt/conda/bin:$PATH

WORKDIR /
RUN conda create -n python_env python=3.8
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN /opt/conda/envs/python_env/bin/pip install pandas numpy pyteomics dash tqdm matplotlib spectrum_utils lxml

# copy the file
COPY ./static /static
COPY ./assets /assets

# copy the code
COPY init-shiny.sh /
COPY web.py /

# modify permission
RUN chmod -R 777 /static
RUN chmod -R 777 /assets
RUN chmod 777 /init-shiny.sh
RUN chmod 777 /web.py

# web
EXPOSE 3838
CMD /opt/conda/envs/python_env/bin/python3.8 /web.py 