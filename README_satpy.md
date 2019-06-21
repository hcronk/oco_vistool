# setup for GOES-ABI data

The satpy package is used for displaying GOES-ABI data, along with cartopy 
for the matplotlib integration.


### satpy setup notes

The bare minimum install would be (tested with python 3.6)

`$ conda create -n satpy36 python=3.6`

`$ conda activate satpy36`

`$ conda install -c conda-forge satpy`

`$ conda install -c conda-forge cartopy`

TBD: this will need a revisit to make sure there is a coherent install that 
includes all dependencies for the vistool (the gdal install in particular).

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
See the [satpy FAQ](https://satpy.readthedocs.io/en/latest/faq.html) for more information.