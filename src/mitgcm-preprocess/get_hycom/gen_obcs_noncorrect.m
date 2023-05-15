clearvars -except input_json_file tool_root run_steps;

clear all
close all
clc

fpath = '/project/rus043_data/generate_ar_wave/mitgcm-preprocess/save_ic_obcs/';
hycompath1 = '/project/rus043_data/generate_ar_wave/mitgcm-preprocess/save_hycom/';

obcs = 'south';
hycompath = hycompath1;

% FMT.mat contains two strings: fmt = 'real*4' and Ieee = 'b'
eval(['load ' fpath 'FMT.mat']);
eval(['load ' hycompath '0021_00_' obcs '.mat']);
lon_hycom = D.Longitude
lat_hycom = D.Latitude
z_hycom = D.Depth;
clear D;

nx = length(lon_hycom);
ny = length(lat_hycom);
nz = length(z_hycom);
nt = 8; % number of days
t0 = 21; % index of starting day

%% writing into bin files
wrslice([hycompath 'lon_hycom_' obcs '.bin'],lon_hycom,1,fmt,Ieee);
wrslice([hycompath 'lat_hycom_' obcs '.bin'],lat_hycom,1,fmt,Ieee);
wrslice([hycompath 'depth_hycom_' obcs '.bin'],z_hycom,1,fmt,Ieee);

totaln = 0
for i = 0:nt-1
  for iHour = 0:0
    k = i+t0;
    t = num2str(k);    m = length(t);    n = 4 - m;    p = blanks(n);    p(:) = '0';
    h = num2str(iHour*3);    m = length(h);    n = 2 - m;    ph = blanks(n);    ph(:) = '0';
    frame = [num2str(p) t '_' num2str(ph) num2str(h) '_' obcs '.mat'];
    eval( ['load ' hycompath  frame]);
    time_hycom = D.Date;
    wrslice([hycompath 'time_hycom_' obcs '.bin'],time_hycom,totaln+1,fmt,Ieee);
    
    ctl = {'surf_el', 'water_u', 'water_v', 'water_temp', 'salinity'};
    for np = 1:length(ctl)
        param = ctl{np};
        eval(['tmp = ' 'D.' param ';']);
        tmp(isnan(tmp)) = 0;
        size(tmp);
        wrslice([hycompath 'hycom_' param '_' obcs '.bin'],tmp,totaln+1,fmt,Ieee);
        clear tmp param;
    end
    totaln = totaln + 1;
    clear D;
  end
end
