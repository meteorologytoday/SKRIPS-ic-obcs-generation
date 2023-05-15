function [matlab_time_num] = hycom_revert_time(hycom_time)
    %% used to get the correct HYCOM timestr from HYCOM mat file
    %% Inverse function of hycom_time.m
    %% Notice that the output string contains hours
    
    matlab_time_num = (hycom_time + 38736)/24.0 + 728872;    
    matlab_time_num = datestr(matlab_time_num,'yyyymmddhh');

return;
