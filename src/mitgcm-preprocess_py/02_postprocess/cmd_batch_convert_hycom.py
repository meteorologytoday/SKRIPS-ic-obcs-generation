import xarray as xr
import numpy as np

import argparse
import convert_grid

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input', type=str, required=True, help='Input hycom data file')
parser.add_argument('--output',  type=str, required=True, help='Output filename.')
parser.add_argument('--varname', type=str, required=True, help='Output filename.')
parser.add_argument('--grid-type', type=str, required=True, help='Output filename.')
parser.add_argument('--grid-dir', type=str, required=True, help='MITgcm grid folder.')
parser.add_argument('--iter-max', type=int, default=50, help='Max iteration when doing data filling.')
parser.add_argument('--check-rng', type=float, nargs=2, default=[-np.inf, np.inf], help='Check range.')
args = parser.parse_args()

print(args)

ds_hycom = xr.open_dataset(args.input)

ZC1 = - ds_hycom.coords["depth"].to_numpy()
YC1 =   ds_hycom.coords["lat"].to_numpy()
XC1 =   ds_hycom.coords["lon"].to_numpy()

interpolated_data, grid2 = convert_grid.convertGrid(ds_hycom[args.varname][0, :, :, :].to_numpy(), args.grid_type, XC1, YC1, ZC1, grid2_dir=args.grid_dir, fill_value=-999, iter_max=args.iter_max, check_rng = args.check_rng)

da_output =xr.DataArray(
    data = np.expand_dims(interpolated_data, 0),
    dims = ["time", "z", "lat", "lon"],
    coords = dict(
        lon=(["lon"], grid2['X']),
        lat=(["lat"], grid2['Y']),
        z=(["z"], grid2['Z']),
        time=ds_hycom.coords["time"],
    ),
).rename(args.varname) 

print("Output: ", args.output)
da_output.to_netcdf(args.output)



"""

print("See if i can output an okay binary file")

dtype = '>f4'
data = da_output.to_numpy().astype(dtype)

if not data.flags.c_contiguous:
    print("Warning: Not row major. Coverting it now...")
    data = np.array(data, order='C', dtype=dtype)


print("Output binary file: ", args.output2)
with open(args.output2, "wb") as f:
    f.write(data.tobytes())


"""
