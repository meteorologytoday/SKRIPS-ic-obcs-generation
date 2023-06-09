import numpy as np
import xarray as xr
from pathlib import Path
import pandas as pd
import argparse


def writeBinary(filename, data, dtype='>f4'):

    data = data.astype(dtype)

    if not data.flags.c_contiguous:
        print("Warning: Not row major. Coverting it now...")
        data = np.array(data, order='C', dtype=dtype)

    with open(filename, "wb") as f:
        f.write(data.tobytes())



def genMITgcmInitCondAndOpenBnd(input_filenames, output_filenames, varname, output_dir, open_bnd_thickness=1, dtype='>f4', output_fmt='netcdf,binary'):
    
    output_fmts = output_fmt.replace(" ", "").split(",")

    # check filename mapping
    needed_keys = ['init_cond', 'open_bnd']
    for v in [input_filenames, output_filenames]:
        for needed_key in needed_keys:
            if not (needed_key in v):
                raise Exception("Error: missing key `%s` in either `input_filenames` or `output_filenames`." % (needed_key, ))

   


    # Create output dir if it does not exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)


    ds_init_cond = xr.open_dataset(input_filenames['init_cond'])

    print(input_filenames['open_bnd'])
    ds_open_bnd  = xr.open_mfdataset(input_filenames['open_bnd'], concat_dim="time", combine="nested")

    # check if they share the same spatial size
    for coordname in ["lat", "lon", "z"]:

        coord_init_cond = ds_init_cond.coords[coordname].to_numpy()
        coord_open_bnd  = ds_open_bnd.coords[coordname].to_numpy()

        if len(coord_init_cond) != len(coord_open_bnd):
            raise Exception("Error: initial condition file and open boundary files do not have the same length in coordinate %s" % (coordname,))
        elif np.amax(np.abs(coord_init_cond - coord_open_bnd)) !=0 :
            raise Exception("Error: initial condition file and open boundary files do not have the same coordinate value of %s" % (coordname,))


    
    # crop boundaries
    # Notice here they are DataArray not DataSet
    ds_open_bnd = dict(
        north = ds_open_bnd[varname][:, :, -1,  :],
        south = ds_open_bnd[varname][:, :,  0,  :],
        east  = ds_open_bnd[varname][:, :,  :, -1],
        west  = ds_open_bnd[varname][:, :,  :,  0],
    )
    


    
    if "binary" in output_fmts:

        print("Output binary...")
        
        filename = "%s.bin" % (output_filenames['init_cond'], )
        print("Output file: ", filename)
        writeBinary(filename, ds_init_cond[varname].to_numpy())
        for k in ds_open_bnd.keys():
        
            filename = "%s_%s.bin" % (output_filenames['open_bnd'], k,)
            print("Output file: ", filename)
            writeBinary(filename, ds_open_bnd[k].to_numpy())



    if "netcdf" in output_fmts:
    
        print("Output netcdf...")
        
        filename = "%s.nc" % (output_filenames['init_cond'], )
        print("Output file: ", filename)
        ds_init_cond.to_netcdf(filename)
        
        for k in ds_open_bnd.keys():
        
            filename = "%s_%s.nc" % (output_filenames['open_bnd'], k,)
            print("Output file: ", filename)
            ds_open_bnd[k].to_netcdf(filename, unlimited_dims="time")



if __name__ == "__main__":

  
    parser = argparse.ArgumentParser(
                    prog = 'generate_mitgcm_bin.py',
                    description = 'Generate initial binary file and open boundary conditions for mitgcm.',
    )

    parser.add_argument('--beg-date', type=str, required=True)
    parser.add_argument('--end-date', type=str, required=True)
    parser.add_argument('--input-dir', type=str, required=True)
    parser.add_argument('--output-dir', type=str, required=True)
    parser.add_argument('--varnames', type=str, nargs='+', default=['water_u', 'water_v', 'water_temp', 'salinity'], choices=['water_u', 'water_v', 'water_temp', 'salinity',])

    args = parser.parse_args()
    print(args)

 
    dt_fmt = "%Y-%m-%d_%H"
    dts = pd.date_range(args.beg_date, args.end_date, freq="D", inclusive="both")
    
    for varname in args.varnames:
        
        input_filenames = dict(
            init_cond = "%s/hycom_%s_%s.nc" % (args.input_dir, dts[0].strftime(dt_fmt), varname),
            open_bnd = [ "%s/hycom_%s_%s.nc" % (args.input_dir, dt.strftime(dt_fmt), varname) for dt in dts ]
        )

        output_filenames = dict(
            init_cond = "%s/init_cond_%s_%s" % (args.output_dir, varname, dts[0].strftime(dt_fmt)),
            open_bnd = "%s/open_bnd_%s_%s_%s" % (args.output_dir, varname, dts[0].strftime(dt_fmt), dts[-1].strftime(dt_fmt)),
        )


        genMITgcmInitCondAndOpenBnd(
            input_filenames,
            output_filenames,
            varname,
            args.output_dir,
            output_fmt="netcdf,binary"
        )

    
