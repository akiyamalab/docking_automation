FROM condaforge/miniforge3:24.11.3-2
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y wget openssh-client

#### Anaconda ####
RUN mamba install -y -c conda-forge openbabel=3.1.1
RUN mamba install -y -c conda-forge fpocket

#### AutoDock Vina ####
RUN pip install vina
# ENV ADVINA_HOME /advina
# ENV ADVINA_VER  1.1.2
# ENV PATH $ADVINA_HOME/$ADVINA_VER/bin:$PATH
# RUN wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz
# RUN tar -zxvf autodock_vina_1_1_2_linux_x86.tgz \
#     && rm autodock_vina_1_1_2_linux_x86.tgz
# RUN mkdir -p $ADVINA_HOME
# RUN mv autodock_vina_1_1_2_linux_x86 $ADVINA_HOME/$ADVINA_VER

#### AutoDock FR (for prepare_ligand, prepare_receptor) ####
ENV ADFR_HOME /adfr
ENV ADFR_VER  1.0
ENV PATH $ADFR_HOME/$ADFR_VER/bin:$PATH
RUN wget https://ccsb.scripps.edu/adfr/download/1028/ADFRsuite_Linux-x86_64_1.0_install \
    && chmod +x ADFRsuite_Linux-x86_64_1.0_install \
    && ./ADFRsuite_Linux-x86_64_1.0_install --mode silent --prefix / \
    && rm ADFRsuite_Linux-x86_64_1.0_install
RUN tar -zxvf ADFRsuite_x86_64Linux_1.0.tar.gz \
    && rm ADFRsuite_x86_64Linux_1.0.tar.gz
RUN cd /ADFRsuite_x86_64Linux_1.0 \
    && ./install.sh -d $ADFR_HOME/$ADFR_VER \
    && rm -r /ADFRsuite_x86_64Linux_1.0

#### AlphaCutter
RUN apt update && apt install -y gcc
RUN apt install -y dssp
RUN mkdir -p /AlphaCutter
RUN wget https://raw.githubusercontent.com/johnnytam100/AlphaCutter/c1b8a96024f71509ac700ef232f632c13e40a886/AlphaCutter_v111.py -O /AlphaCutter/AlphaCutter_v111.py
# conda should be prioritized compared to adfr
ENV PATH /opt/conda/bin:/opt/conda/condabin:$PATH 

