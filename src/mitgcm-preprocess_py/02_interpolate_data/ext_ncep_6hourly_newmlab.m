clear all
close all
clc

dpath1 = '/project_shared/DATA_DOWNLOADS/DATA/NCEP_NCAR-Reanalysis/Hourly_6hours/';
opath = '/home/rus043/kaust_project/codeNew/regional-coupled-ocean-atmosphere-model/coupler/script-preprocess/ic_obcs_20180101_new/';
fpath = opath;

eval ([ 'load ' fpath 'FMT.mat']);

ctl = {'air.2m', 'shum.2m', 'prate.sfc', 'dlwrf.sfc', 'ulwrf.sfc','dswrf.sfc', 'uswrf.sfc', 'nlwrs.sfc','nswrs.sfc','uwnd.10m','vwnd.10m',};
ctl_i = {'air', 'shum', 'prate', 'dlwrf', 'ulwrf','dswrf', 'uswrf', 'nlwrs','nswrs','uwnd','vwnd',};

prec='double';

np = length(ctl);
for jp = 1:np
    p = [ctl{jp}];
    param = [ctl_i{jp}];
    dpath = [dpath1 p '/'];
    dpath
    
    for it = [2016:2017]
        file = [p '.gauss.' num2str(it) '.nc']
        
        ncid = netcdf.open([dpath file],'NOWRITE');
        [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
        
        for varid =1:numvars
            varlist(varid) = {netcdf.inqVar(ncid,varid-1)};
        end
        
        ind=find(ismember(varlist,param));
        tmp = netcdf.getVar(ncid,ind-1,prec);
        
        [varname vartype vardimIDs varatts] = netcdf.inqVar(ncid,ind-1);
        
        for varatid =1:varatts
            attlist(varatid) = {netcdf.inqAttName(ncid,ind-1,varatid-1)};
        end
        
        miss = netcdf.getAtt(ncid,ind-1,'missing_value',prec)
        if miss < 0
            tmp(tmp <= miss) = NaN;
        else
            tmp(tmp >= miss) = NaN;
        end
        
        if isempty(find(ismember(attlist,'scale_factor')))
            scale = 1.0
        else
            scale = netcdf.getAtt(ncid,ind-1,'scale_factor',prec)
        end
        
        if isempty(find(ismember(attlist,'add_offset')))
            offset = 0.0
        else
            offset = netcdf.getAtt(ncid,ind-1,'add_offset',prec)
        end
        
        tmp = offset + (tmp.*scale);
        
        %% we need to flip as NCEPlat is flipped
        tmp = flipdim(tmp,2);
        size(tmp)
        wrslice([opath 'ncep_' param '_6h_' num2str(it) ],tmp,1,fmt,Ieee);
        netcdf.close(ncid);
    end
    clear p param tmp dpath
end
