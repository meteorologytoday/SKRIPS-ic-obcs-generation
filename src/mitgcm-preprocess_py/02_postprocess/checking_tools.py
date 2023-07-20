import numpy as np
import scipy.interpolate
from scipy.interpolate import RegularGridInterpolator
import functools
import MITgcmutils
import argparse

def checkingInitialData(init_file, grid_dir, dtype='>f4', extra_missing_value=0):

    ocn_mask = MITgcmutils.rdmds('%s/hFacC' % (grid_dir,)) != 0
    
    data = np.fromfile(init_file, dtype=dtype)

    if len(data) != ocn_mask.size:
        raise Exception("Error: init file %s has %d numbers while the mask has %d." % (init_file, len(data), ocn_mask.size))

    data = data.reshape(ocn_mask.shape)

    ocn_cnt = np.sum(ocn_mask)
    valid_pts = ocn_mask & np.isfinite(data)

    finite_data_cnt = np.sum(valid_pts)

    if finite_data_cnt == ocn_cnt:
        
        print("All data on ocean grid points are finite.")

    else:
        raise Exception("Warning: %d out of %d ocean grid points are NaN." % (finite_data_cnt, ocn_cnt))


    if not np.isnan(extra_missing_value):
        print("Extra missing value received: %f" % (extra_missing_value,))
        extra_missing_pts = valid_pts & (data == extra_missing_value)
        extra_missing_cnt = np.sum(extra_missing_pts)
        
        if extra_missing_cnt == 0:
            
            print("No missing value found on all ocean grid points.")

        else:
            raise Exception("Warning: %d out of %d ocean grid points have missing value %f." % (extra_missing_cnt, ocn_cnt, extra_missing_value))

       
        

if __name__ == "__main__":


    parser = argparse.ArgumentParser(
                        prog = 'batch_convert_hycom.py',
                        description = 'Convert hycom data to mitgcm grid in batch.',
    )

    parser.add_argument('--init-file', type=str, required=True)
    parser.add_argument('--grid-dir', type=str, required=True)
    parser.add_argument('--extra-missing-value', type=float, default=np.nan)

    args = parser.parse_args()
    print(args)

    print("Checking file %s using grid dir %s" % (args.init_file, args.grid_dir, ))
    checkingInitialData(args.init_file, args.grid_dir, extra_missing_value=args.extra_missing_value)
