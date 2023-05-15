clc
clear all;
close all;

run ~/matlab_bin/pathdef.m

% My model paths.
fpath = '/home/t2hsu/models/skrips_testcase/case02_NEPacific/mitgcm-preprocess/save_ic_obcs/';
hycompath = '/home/t2hsu/models/skrips_testcase/case02_NEPacific/mitgcm-preprocess/save_hycom/';

eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);
% xc = XC(:,1); yc = YC(1,:); zc = RC;
% xf = XG(:,1); yf = YG(1,:);
% nxc = length(xc); nyc = length(yc);nzc = length(zc);

% HYCOM outputs
eval(['load ' hycompath '0021_2018012100.mat'])

hycom_lon = D.Longitude;
hycom_lat = D.Latitude;
hycom_z = D.Depth;
hycom_time = D.Date;
hycomu = D.water_u;
hycomv = D.water_v;
hycomt = D.water_temp;
hycoms = D.salinity;

date = ['ar_' hycom_time]

% Extract HYCOM fields at this time
[Tm,xc,yc,zc] = exinthycomic(fpath,'T',hycomt,hycom_time,hycom_lon,hycom_lat,hycom_z,xc,yc,xf,yf,zc);
[Sm,xc,yc,zc] = exinthycomic(fpath,'S',hycoms,hycom_time,hycom_lon,hycom_lat,hycom_z,xc,yc,xf,yf,zc);
[Um,xc,yc,zc] = exinthycomic(fpath,'U',hycomu,hycom_time,hycom_lon,hycom_lat,hycom_z,xc,yc,xf,yf,zc);
[Vm,xc,yc,zc] = exinthycomic(fpath,'V',hycomv,hycom_time,hycom_lon,hycom_lat,hycom_z,xc,yc,xf,yf,zc);

T = reshape(Tm,nxc,nyc,nzc);
S = reshape(Sm,nxc,nyc,nzc);
U = reshape(Um,nxc,nyc,nzc);
V = reshape(Vm,nxc,nyc,nzc);

% Save outputs
wrslice([fpath 'hycom_T_' date '.bin'],T,1,fmt,Ieee);
wrslice([fpath 'hycom_S_' date '.bin'],S,1,fmt,Ieee);
wrslice([fpath 'hycom_U_' date '.bin'],U,1,fmt,Ieee);
wrslice([fpath 'hycom_V_' date '.bin'],V,1,fmt,Ieee);

% Check if data are available over all model grid
check3dmask(fpath,['hycom_T_' date '.bin'],1,'T');
check3dmask(fpath,['hycom_S_' date '.bin'],1,'S');
check3dmask(fpath,['hycom_U_' date '.bin'],1,'U');
check3dmask(fpath,['hycom_V_' date '.bin'],1,'V');

