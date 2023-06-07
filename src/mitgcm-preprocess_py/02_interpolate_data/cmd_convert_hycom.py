import xarray as xr
import numpy as np

import argparse
import convert_grid

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input', type=str, required=True, help='Input hycom data file')
parser.add_argument('--output', type=str, required=True, help='Output filename.')
parser.add_argument('--grid-type', type=str, required=True, help='Output filename.')
parser.add_argument('--grid-dir', type=str, required=True, help='MITgcm grid folder.')
args = parser.parse_args()

print(args)

ds_hycom = xr.open_dataset(args.input)

ZC1 = - ds_hycom.coords["depth"].to_numpy()
YC1 =   ds_hycom.coords["lat"].to_numpy()
XC1 =   ds_hycom.coords["lon"].to_numpy()
interpolated_data, grid2 = convert_grid.convertGrid(ds_hycom[args.varname][0, :, :, :].to_numpy(), args.grid_type, XC1, YC1, ZC1, grid2_dir=args.grid_dir, fill_value=0.0)

ds_output =xr.DataArray(
    data = np.expand_dims(interpolated_data, 0)
    dims = ["time", "z", "lat", "lon"],
    coords = dict(
        lon=(["lon"], grid2.X),
        lat=(["lat"], grid2.Y),
        z=(["z"], grid2.Z),
        time=ds_hycom.coords["time"],
    ),
) 

ds_output.to_netcdf(args.output)









