FROM continuumio/miniconda3:4.6.14

RUN conda install -y -c conda-forge h5py=2.9.0 \
    pandas=0.25.0 cartopy=0.17.0 gdal=3.0.1
