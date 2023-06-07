import xarray as xr
import numpy as np

import argparse


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--input', type=str, required=True, help='Input hycom data file')
parser.add_argument('--output', type=str, required=True, help='Output filename.')
parser.add_argument('--grid-dir', type=str, required=True, help='MITgcm grid folder.')
args = parser.parse_args()

print(args)

z =






