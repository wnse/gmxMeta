

# 流程

0. fastp
1. humann
2. humann_split_stratified_table
3. megahit
4. rgi
    rgi main -n 8 -i contig/ma-mix-2159yk.contigs.fa -o rgi_test_out/ma-mix-2159yk -a DIAMOND --low_quality
5. vfdb
    diamond blastx --db /data/database/VFDB/VFDB_setB_pro.dmnd --query ma-mix-2159yk.contigs.fa.temp.contigToORF.fsa --out ma-mix-2159yk.contigs.fa.temp.contigToORF.vfdb.txt --outfmt 6 --evalue 1e5 --threads 8 
6. bowtie2
    bowtie2-build ma-mix-2159yk.contigs.fa.temp.contigToORF.fsa orf
    bowtie2 -p 8 -x orf -1 ../ma-mix-2159yk_FKDL190725375-1a-9/ma-mix-2159yk_FKDL190725375-1a-9_1.fq.gz -2 ../ma-mix-2159yk_FKDL190725375-1a-9/ma-mix-2159yk_FKDL190725375-1a-9_2.fq.gz -S raw2orf.sam
     samtools view -bS -@ 4 raw2orf.sam |samtools sort -@ 4 -o raw2orf.bam -
     samtools index raw2orf.bam
     samtools idxstats raw2orf.bam > raw2orf.idxstats.txt


## script

1. 执行流程
2. 处理结果
    1. 物和丰度
    2. 物种多样性
    3. 代谢通路丰度
    4. 抗性基因丰度
    5. 毒力因子丰度
3. 生成格式化文件
    1. 配置文件
        1. 有益菌
        2. 有害菌
        3. 条件致病菌
        4. 肠型
        5. 多样性
        6. 抗性基因抗生素种类
        7. 毒力因子
        8. 疾病风险
        9. 代谢水平


## docker
```
FROM centos:latest
MAINTAINER yangkai07@gmail.com

## SET WORKING DIRECTORY
WORKDIR /data

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
&& mamba create -y -n humann -c conda-forge -c bioconda -c biobakery metaphlan humann fastp \
&& mamba create -y -n rgi -c conda-forge -c bioconda -c defaults rgi \
&& mamba create -y -n megahit -c bioconda megahit

CMD ["/bin/bash"]

```


# 建立conda环境

FROM centos:latest
MAINTAINER yangkai07@gmail.com

## SET WORKING DIRECTORY
WORKDIR /metagenome

# COPY . /metagenome
ENV PATH /opt/conda/bin:$PATH

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

RUN cd /cnv_ana \
&& sh Miniconda3-py38_22.11.1-1-Linux-x86_64.sh -b -p /opt/conda \
&& rm -rf Miniconda3-py38_22.11.1-1-Linux-x86_64.sh \
&& ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
        && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
&& /opt/conda/bin/conda clean -afy \
&& source /opt/conda/bin/activate \
&& conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ \
&& conda env create -f env_config.yml


ENV PATH /opt/conda/envs/pgt/bin:$PATH
ENV CONDA_DEFAULT_ENV pgt

# SHELL ["conda", "run", "-n", "pgt", "/bin/bash", "-c"]

RUN R -e "install.packages('BiocManager',repos='http://cran.us.r-project.org')" \
&& R -e "BiocManager::install('DNAcopy')" \
&& R -e "install.packages(c('ggplot2','cowplot','hash','gridExtra'),repos='http://cran.us.r-project.org')" \
&& rm -rf /tmp/downloaded_packages/

#CMD ["/bin/bash"]
#CMD ["source /opt/conda/bin/activate"]
CMD ["/opt/conda/envs/pgt/bin/gunicorn", "-c", "./gunicorn.conf.py", "wsgi:app"]