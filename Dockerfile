FROM continuumio/miniconda3:4.9.2

RUN conda install -y -c conda-forge satpy=0.29.0 netcdf4=1.5.6 h5py=3.2.1 matplotlib=3.3.4 \
    pandas=1.1.4 cartopy=0.19.0 owslib=0.24.1 pillow=8.2.0 boto3=1.17.101
    
COPY . .

CMD ["python", "oco_vistool.py"]