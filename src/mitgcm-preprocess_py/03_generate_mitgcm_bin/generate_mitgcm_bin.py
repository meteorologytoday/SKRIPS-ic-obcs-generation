import numpy as np
import xarray as xr


def writeBinary(filename, data, dtype='>f4'):

    data = da_output.to_numpy().astype(dtype)

    if not data.flags.c_contiguous:
        print("Warning: Not row major. Coverting it now...")
        data = np.array(data, order='C', dtype=dtype)

    with open(filename, "wb") as f:
        f.write(data.tobytes())



def genMITgcmBinaries(init_cond_filename, open_bnd_filenames, varnames, output_dir, open_bnd_thickness=1, dtype='>f4'):

    ds_init_cond = xr.open_dataset(init_cond_filename)
    ds_open_bnd  = xr.open_mfdataset(open_bnd_filenames)

    # check if they share the same spatial size

    for coordname in ["lat", "lon", "z"]:

        coord_init_cond = ds_init_cond.coords[coordname].to_numpy()
        coord_open_bnd  = ds_open_bnd.coords[coordname].to_numpy()

        if len(coord_init_cond) != len(coord_open_bnd):
            raise Exception("Error: initial condition file and open boundary files do not have the same length in coordinate %s" % (coordname,))
        elif np.amax(np.abs(coord_init_cond - coord_open_bnd)) !=0 :
            raise Exception("Error: initial condition file and open boundary files do not have the same coordinate value of %s" % (coordname,))

