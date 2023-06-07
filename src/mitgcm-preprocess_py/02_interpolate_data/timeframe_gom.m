function [time_num] = timeframe_gom(a)
%% used to get the correct HYCOM timestr from HYCOM mat file
k = num2str(a);
yr = str2num(k(:,1:4));
mn = str2num(k(:,5:6));
dy = str2num(k(:,7:8));
hr = zeros(size(k,1),1);
mi = zeros(size(k,1),1);
sc = zeros(size(k,1),1);
time_num = datenum(yr,mn,dy,hr,mi,sc);
return;
