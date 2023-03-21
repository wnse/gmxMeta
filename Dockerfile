FROM centos:latest
MAINTAINER yangkai07@gmail.com

## SET WORKING DIRECTORY
WORKDIR /data

COPY . /gmxMeta

RUN cd /etc/yum.repos.d/ \
&& sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-* \
&& sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-* \
&& yum install wget -y \
&& yum install which -y \
&& yum install libXrender-0.9.10-7.el8.i686 -y \
&& yum install cairo-devel -y \
&& wget -O /etc/yum.repos.d/CentOS-Base.repo https://mirrors.aliyun.com/repo/Centos-vault-8.5.2111.repo \
&& yum makecache

## TIMEZONE
RUN cp /usr/share/zoneinfo/Asia/Shanghai /etc/localtime
RUN echo Asia/Shanghai > /etc/timezone

RUN cd /data \
&& wget -O Miniconda3-latest-Linux-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
&& sh Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
&& rm -rf Miniconda3-py38_22.11.1-1-Linux-x86_64.sh \
&& ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
        && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
&& /opt/conda/bin/conda clean -afy \
&& source /opt/conda/bin/activate \
&& conda install -y mamba -c conda-forge \
&& mamba install -y -c anaconda networkx numpy pandas openpyxl \
&& mamba install -y -c conda-forge matplotlib \
&& mamba create -y -n humann -c conda-forge -c bioconda -c biobakery metaphlan humann fastp \
&& mamba create -y -n rgi -c conda-forge -c bioconda -c defaults rgi \
&& mamba create -y -n megahit -c bioconda megahit

#CMD ["/bin/bash"]
ENV PATH /opt/conda/bin:$PATH
ENV CONDA_DEFAULT_ENV base

CMD ["source /opt/conda/bin/activate"]
