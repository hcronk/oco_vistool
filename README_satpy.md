# setup for GOES-ABI data

The satpy package is used for displaying GOES-ABI data, along with cartopy 
for the matplotlib integration.


### satpy setup notes

The bare minimum install would be (tested with python 3.6, conda 4.6.14, June 2019)

`$ conda create -n <envname> python=3.6`

`$ conda activate <envname>`

`$ conda install -c conda-forge satpy`

`$ conda install -c conda-forge cartopy`

`$ conda install future`

`$ conda install -c conda-forge gdal`

### environment notes

satpy uses the dask library for fast multithreaded processing of the large ABI 
image arrays. While this gives satpy a fast throughput for image processing, 
for most typical vistool applications we are using a small sub region of the 
full disk or CONUS ABI images. Thus, it is usually better to limit the 
amount of CPU that dask will consume. The two important shell environment 
variables are `DASK_NUM_WORKERS` (to control the number of dask worker 
processes that are started) and `OMP_NUM_THREADS` to control the number of 
CPU threads that are consumed by each dask worker.

For example, by setting the following in the `.bash_profile` or similar:

`export OMP_NUM_THREADS=4`
`export DASK_NUM_WORKERS=2`

this would limit satpy to use 4x2 = 8 CPU on a multiprocessor system.

With limited testing, for small image creation (1 x 1 degree in lat/lon), 
there is very limited benefit by increasing the number of workers or 
threads. In these cases it might make sense to reduce the number of 
workers and threads to 1. It is advisable to test a few combinations 
on your system.

See the [satpy FAQ](https://satpy.readthedocs.io/en/latest/faq.html) for more information.


The satpy library also relies on pytroll for image processing, which uses 
an environment variable to control the internal "chunk" size for image 
processing. This value is very important for efficiency when making small 
image creation, since chunks that are outside the area of interest will be 
skipped. The default chunk size is 4096 if the environment variable is not 
set. For small image processing, a value on the order of 1024 or 512 may 
be optimal:

`export PYTROLL_CHUNK_SIZE=512`
