import pandas as pd
from datetime import datetime

incomplete_OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC'
 #GLBv0.08/expt_93.0';


hycom_beg_dt = pd.Timestamp(datetime.strptime('2000-01-01', "%Y-%m-%d"))

def hycomTime2Datetime(hycom_time: int):
    return hycom_beg_dt + pd.Timedelta(hycom_time, unit='h')
    
def datetime2hycomTime(dt):
    return int( (dt - hycom_beg_dt) / pd.Timedelta(1, unit='h') )
    

