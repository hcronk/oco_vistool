FROM continuumio/anaconda3:5.3.0

RUN conda install -y -c conda-forge numpy=1.16.3
RUN conda install -y -c conda-forge numpy-base=1.16.3
RUN conda install -y -c conda-forge gdal=2.3.3
RUN conda install -y -c conda-forge cartopy=0.16.0
