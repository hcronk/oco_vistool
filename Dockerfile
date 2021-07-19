FROM continuumio/miniconda3:4.9.2

RUN conda install -y -c conda-forge netcdf4=1.5.7 h5py=2.10.0 \
    pandas=1.1.4 cartopy=0.18.0 owslib=0.24.1 pillow=7.2.0