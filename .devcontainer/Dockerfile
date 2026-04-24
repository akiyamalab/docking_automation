FROM condaforge/miniforge3:24.11.3-2
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y wget openssh-client

#### Anaconda ####
RUN mamba install -y -c conda-forge openbabel=3.1.1
# fpocket: conda-forgeパッケージが利用可能な場合はそちらを優先、
# 環境によってはソースビルドで対応（後段のRUN参照）
RUN mamba install -y -c conda-forge fpocket || true

#### fpocket source build (conda失敗時のフォールバック) ####
# CFLAGS: 新しい GCC (14+) は -Wincompatible-pointer-types を既定でエラー扱いするため、
# fpocket 同梱の qhull が通らない。警告に戻すことでビルドを通す。
RUN if ! command -v fpocket >/dev/null 2>&1; then \
      apt-get install -y cmake gcc g++ git \
      && git clone https://github.com/Discngine/fpocket.git /tmp/fpocket \
      && cd /tmp/fpocket \
      && make CFLAGS="-Wno-error=incompatible-pointer-types -Wno-error=implicit-function-declaration" -j$(nproc) \
      && cp /tmp/fpocket/bin/fpocket /usr/local/bin/fpocket \
      && cp /tmp/fpocket/bin/dpocket /usr/local/bin/dpocket \
      && cp /tmp/fpocket/bin/tpocket /usr/local/bin/tpocket \
      && cp /tmp/fpocket/bin/mdpocket /usr/local/bin/mdpocket \
      && rm -rf /tmp/fpocket; \
    fi

#### AutoDock Vina ####
RUN pip install vina

#### Dimorphite-DL (protonation at physiological pH) ####
RUN pip install dimorphite-dl==2.0.2
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
# ccsb.scripps.edu は TLS 証明書の SAN が一致しないため --no-check-certificate で回避
RUN wget --no-check-certificate https://ccsb.scripps.edu/adfr/download/1028/ADFRsuite_Linux-x86_64_1.0_install \
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

#### Uni-Dock (GPU docking backend, Phase 3+) ####
# nvcc (CUDA compiler) is not available in the runtime-only CUDA image,
# so we use the pre-built binary from GitHub Releases (CUDA 12.0 compatible with 12.x).
ENV CUDA_HOME=/usr/local/cuda-12.4
ENV LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
RUN pip install unidock-tools
RUN wget -q "https://github.com/dptech-corp/Uni-Dock/releases/download/1.1.0/unidock-1.1.0-cuda120-linux-x86_64" \
    -O /usr/local/bin/unidock \
    && chmod +x /usr/local/bin/unidock

# Remove cuda-compat so nvidia-container-runtime auto-mounts host libcuda (portable across driver versions)
RUN apt-get remove -y cuda-compat-12-4 2>/dev/null || true

#### Uni-Dock 2 (v0.6.x, separate conda env) ####
# - base env とは numpy/MDAnalysis ABI が衝突するため独立 env にする
# - quetz.dp.tech が唯一の配布経路 (pip 未対応)
# - openbabel/rdkit の pip wheel が system libstdc++ を引き込むと msys 拡張が
#   CXXABI_1.3.15 not found で落ちるため、libstdcxx-ng を --force-reinstall する
RUN mamba create -n unidock2 python=3.10 -c conda-forge -y \
    && mamba install -n unidock2 unidock2 cuda-version=12.0 \
       -c http://quetz.dp.tech:8088/get/baymax -c conda-forge -y \
    && mamba install -n unidock2 -c conda-forge openbabel -y \
    && mamba install -n unidock2 -c conda-forge libstdcxx-ng libgcc-ng --force-reinstall -y \
    && /opt/conda/envs/unidock2/bin/pip install vina h5py

