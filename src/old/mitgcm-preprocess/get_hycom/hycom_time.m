function [hycom_time_num] = hycom_time(date, hr, fmt)
    %% used to get the correct HYCOM timestr from HYCOM mat file
    matlab_datenum = datenum(date, fmt);
    hour_num = str2num(hr);
    
    % for GOFS3.1
    % As a test: datestr(728872 - 38736/24, 'yyyy-mm-dd') == '1991-03-01'
    hycom_time_num = (matlab_datenum - 728872)*24 - 38736 + hour_num;
return;
