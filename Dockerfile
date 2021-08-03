ARG BASEIMAGE=ontresearch/base-workflow-image:v0.1.1

FROM ubuntu:20.04 as clair
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        git \
        wget \
        bzip2 \
        make \
        g++ \
        libboost-graph-dev && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/bin

# install anaconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz -P /opt/models && \
    tar -zxvf /opt/models/clair3_models.tar.gz -C /opt/models && \
    rm /opt/models/clair3_models.tar.gz && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair3 python=3.6.10 -y

ENV PATH /opt/conda/envs/clair3/bin:$PATH
ENV CONDA_DEFAULT_ENV clair3

RUN /bin/bash -c "source activate clair3" && \
    conda install -c conda-forge pypy3.6 -y && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install mpmath==1.2.1 && \
    pip install tensorflow-cpu==2.2.0 && \
    pip install tensorflow-addons==0.11.2 tables==3.6.1 && \
    conda install -c anaconda pigz==2.4 -y && \
    conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y && \
    conda install -c conda-forge -c bioconda samtools=1.10 -y && \
    conda install -c conda-forge -c bioconda whatshap=1.0 -y && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip

RUN git clone https://github.com/HKU-BAL/Clair3.git .

RUN cd /opt/bin/preprocess/realign && \
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    echo "source activate clair3" > ~/.bashrc

FROM $BASEIMAGE as base
ARG ENVFILE=environment.yaml

COPY --from=clair /opt/models/ont/ /opt/models/ont/
COPY --from=clair /opt/bin/run_clair3.sh /opt/bin/
COPY --from=clair /opt/bin/clair3.py /opt/bin/
COPY --from=clair /opt/bin/scripts/clair3.sh /opt/bin/scripts/
COPY --from=clair /opt/bin/scripts/clair3_ont_quick_demo.sh /opt/bin/scripts/
COPY --from=clair /opt/ /opt/
COPY $ENVFILE $HOME/environment.yaml

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive


RUN \
    . $CONDA_DIR/etc/profile.d/mamba.sh \
    && micromamba activate \
    && micromamba install -n base --file $HOME/environment.yaml \
    && micromamba clean --all --yes \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME \
    && rm -rf $CONDA_DIR/conda-meta \
    && rm -rf $CONDA_DIR/include \
    && rm -rf $CONDA_DIR/lib/python3.*/site-packages/pip \
    && find $CONDA_DIR -name '__pycache__' -type d -exec rm -rf '{}' '+'

USER $WF_UID
WORKDIR $HOME


