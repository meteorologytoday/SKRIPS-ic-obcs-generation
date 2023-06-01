import xarray as xr
import pandas as pd
from datetime import datetime
import hycom_share
import json
import numpy as np

def findfirst(arr):
    return np.argmax(arr)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.int64):
            return int(obj)
        else:
            return json.JSONEncoder.default(self, obj)

hycom_dataset_list = {

    "GLBv0.08" : [
        "expt_93.0",
#        "expt_92.9",
#        "expt_57.7",
#        "expt_92.8",
#        "expt_57.2",
#        "expt_56.3",
#        "expt_53.X",
    ],

# Warning: GLBa0.08 use a completely different
# set of time format and coordinate names. This is terrible.    
#    "GLBa0.08" : [
#        "expt_91.2",
#        "expt_91.1",
#        "expt_91.0",
#        "expt_90.9",
#        "expt_90.8",
#        "expt_90.6",
#   ],

}

def scanHycomInfo(dataset_list=None):

    if dataset_list is None:
        dataset_list = hycom_dataset_list

    scanned_info = {}
    
    for dataset_name, subset_names in dataset_list.items():
        
        scanned_info[dataset_name] = {}
        
        for subset_name in subset_names:
                
            url = "%s/%s/%s" % (
                hycom_share.incomplete_OpenDAP_URL,
                dataset_name,
                subset_name,
            )

            print("Scanning %s/%s ..." % (dataset_name, subset_name, ))

            info = {}            
            with xr.open_dataset(url, decode_times=False) as ds:

                hycom_times = ds.coords["time"].to_numpy()
                beg_dt = hycom_share.hycomTime2Datetime(hycom_times[0])
                end_dt = hycom_share.hycomTime2Datetime(hycom_times[-1])


                lat = ds.coords["lat"].to_numpy()
                
                # Detect if longitude does not start from zero
                lon = ds.coords["lon"].to_numpy() % 360
                dlon = lon[1:] - lon[:-1]

                lon_discont = np.abs(dlon) > 90
                num_of_discont = np.sum(lon_discont)
                if num_of_discont >= 2: 
                    raise Exception("More than one discontinuity in longitude this dataset.")
 
                elif num_of_discont == 0:
                    lon_beg = findfirst(lon_discont)
                
                else:
                    lon_beg = 0

                scanned_info[dataset_name][subset_name] = {
                    'time_rng' : (beg_dt.strftime("%Y-%m-%d"), end_dt.strftime("%Y-%m-%d")),
                    'time' : hycom_times,
                    'lat' : lat,
                    'lon' : lon,
                    'lon_beg' : lon_beg,
                }

            print("Detected time for %s/%s : %s ~ %s" % (dataset_name, subset_name, *scanned_info[dataset_name][subset_name]['time_rng']) )

        
             
    return scanned_info 

if __name__ == "__main__" : 

    print("Run as an independent program. Scan hycom info now.")
    result = scanHycomInfo()

    output_filenname = "hycom_info.json"
    print("Output info to file: %s" % (output_filenname,))
    with open(output_filenname, 'w') as f:
        json.dump(result, f, indent=4, cls=NumpyEncoder)


    output_filenname = "hycom_info.npz"
    print("Output info to file: %s" % (output_filenname,))
    np.savez(output_filenname, **result)


    print(dict(np.load(output_filenname, allow_pickle=True)))




