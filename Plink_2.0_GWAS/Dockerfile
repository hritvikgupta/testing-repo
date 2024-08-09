FROM continuumio/miniconda3
WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends libz-dev g++ gcc git wget tar unzip make \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir plink

RUN wget -P plink https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip

WORKDIR plink

RUN unzip plink2_linux_x86_64_20240105.zip

RUN rm -rf plink2_linux_x86_64_20240105.zip

WORKDIR /app

RUN mv plink/plink2 /usr/bin

RUN git clone https://github.com/dummy-Informatics-and-Genomics/NEAT-Plots.git

RUN mv NEAT-Plots/manhattan-plot/ /app/

RUN git clone https://github.com/RitchieLab/biofilter.git

RUN mv biofilter/trunk/* /app/

RUN git clone https://github.com/RitchieLab/LOKI.git

RUN mv LOKI/trunk/* /app/

USER $CONDA_USER

RUN conda install -y -n base -c conda-forge scikit-learn dominate wget libtiff conda-build scipy pandas seaborn matplotlib numpy apsw sqlite && \
    conda clean --all --yes

WORKDIR /app/manhattan-plot/

RUN conda develop .

WORKDIR /app

RUN /opt/conda/bin/python setup.py install
