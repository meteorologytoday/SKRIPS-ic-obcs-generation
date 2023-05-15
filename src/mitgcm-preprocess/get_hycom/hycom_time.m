function [hycom_time_num] = hycom_time(date,hr)
%% used to get the correct HYCOM timestr from HYCOM mat file
matlab_datenum = datenum(date,'yyyymmdd');
hour_num = str2num(hr);

% for GOFS3.1
hycom_time_num = (matlab_datenum - 728872)*24 - 38736 + hour_num;


return;
