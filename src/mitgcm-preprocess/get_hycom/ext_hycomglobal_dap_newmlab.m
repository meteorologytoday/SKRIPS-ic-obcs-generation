function [] = ext_hycomglobal_dap_newmlab(opath,frame,str_date,hour,z,y,x)

D = struct('Date', [], 'Depth', [],'Latitude', [], 'Longitude', [],'ssh',[],'u',[],'v',[],'temperature',[],'salinity',[]);

t = num2str(frame);
p = blanks(4 - length(t));p(:) = '0';
frame = [num2str(p) t];
disp (frame)

prec='double';

xmin = x(1);
xmax = x(end);
ymin = y(1);
ymax = y(end);
zmin = z(1);
zmax = z(end);

deltax = xmax - xmin + 1 ;
deltay = ymax - ymin + 1 ;
deltaz = zmax - zmin + 1 ;
%% will need to update for an array for later
deltat = 1;

Z = zmax+1;

%% 2008-09-18 → 2009-05-06
%% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6
%% time = [0:230];
%% saving from Jan 1, 2009
%% time = [105:1:230];nt = 126;

%% 2009-05-07 → 2011-01-02
%% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8
%% time = [0:1:605];nt = 606;

%% 2011-01-03 → 2013-08-20
%% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9
%% time = [0:1:960];nt = 961;

%% 2013-08-21 → 2014-04-04
%% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0
%% time = [0:1:226];nt = 227;

%% 2014-04-05 → 2016-04-17
%% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1
%% time = [0:1:743];nt = 744;

%% 2016-04-18 → Present
%% http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2
%% time = [0:1:148];nt = 149;
ed = datenum(date);
indtime = timeframe_gom(str_date);
if ((datenum(2008,9,18) >= indtime) || (indtime <= datenum(2009,5,6)))
    OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.6';
elseif ((datenum(2009,5,7) >= indtime) || (indtime <= datenum(2011,1,2)))
    OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.8';
elseif ((datenum(2011,1,3) >= indtime) || (indtime <= datenum(2013,8,20)))
    OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_90.9';
elseif ((datenum(2013,8,21) >= indtime) || (indtime <= datenum(2014,4,4)))
    OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.0';
elseif ((datenum(2014,4,5) >= indtime) || (indtime <= datenum(2016,4,17)))
    OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.1';
elseif ((datenum(2016,4,18) >= indtime) || (indtime <= ed))
    OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2';
end

ncid = netcdf.open(OpenDAP_URL,'NOWRITE');
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);

for varid =1:numvars
    varlist(varid) = {netcdf.inqVar(ncid,varid-1)};
end


%%  First get the hycom date and time
param = 'Date';
ind=find(ismember(varlist,param));

Temp_Date = netcdf.getVar(ncid,ind-1);
time = find(Temp_Date==str2num(str_date))
hycom_date = Temp_Date(time)
hycom_time = timeframe_gom(hycom_date)
D.Date = hycom_time;

%% to make array indicies correct for
%% multidimensional arrays that start
%% at 0 for opendap
time = time-1;

param = 'Depth';
ind = find(ismember(varlist,param));

Temp_Depth = netcdf.getVar(ncid,ind-1,prec);
hycom_depth = -1*Temp_Depth;
%% HYCOM GOES ONLY UPTO 5500 m, so adding one more layer with same field as
hycom_depth(end+1) = -6500;
D.Depth = hycom_depth;

param = 'Latitude';
ind=find(ismember(varlist,param));

hycom_lat = netcdf.getVar(ncid,ind-1,[xmin,ymin],[deltax,deltay],prec);
D.Latitude = hycom_lat(1,:);
D.Latitude = permute(D.Latitude,[2,1]);

param = 'Longitude';
ind=find(ismember(varlist,param));

Temp_Lon = netcdf.getVar(ncid,ind-1,[xmin,ymin],[deltax,deltay],prec);
hycom_lon = Temp_Lon(:,1)
if hycom_lon < 0; hycom_lon = hycom_lon + 360; end
if hycom_lon > 360; hycom_lon = hycom_lon - 360; end
D.Longitude = hycom_lon;
D.Longitude = permute(D.Longitude,[2,1]);

params = {'ssh', 'u', 'v', 'temperature', 'salinity'};
for np = 1:length(params)
    param = params{np};
    ind=find(ismember(varlist,param));
    if (strcmp(param,'ssh') == 1)
        tmp = netcdf.getVar(ncid,ind-1,[xmin,ymin,time],[deltax,deltay,deltat],prec);
    else
        tmp = netcdf.getVar(ncid,ind-1,[xmin,ymin,zmin,time],[deltax,deltay,deltaz,deltat],prec);
        %% HYCOM GOES ONLY UPTO 5500 m, so adding one more layer with same field as
        tmp(:,:,Z+1) = tmp(:,:,Z);
    end
    
    %% miss_value = 1.2677 * 1.0e+30
    miss_value = 1.e30;fill = 0;
    tmp(tmp > miss_value)  = fill;
    eval(['D.' param ' = tmp;']);
    clear tmp param
end

%% saving orginal HYCOM  fields into mat file
varList = 'D';

%% for OB slices we dont append date
%% savename = [opath frame '_' num2str(hycom_date) '.mat'];
savename = [opath frame '.mat'];
eval(['save ',savename,' ',varList]);
netcdf.close(ncid);

