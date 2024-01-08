clc
clear all
close all

run ~/matlab_bin/pathdef.m

%% My model paths.
fpath = '/home/t2hsu/models/skrips_testcase/case02_NEPacific/mitgcm-preprocess/save_ic_obcs/';
hycompath1 = '/home/t2hsu/models/skrips_testcase/case02_NEPacific/mitgcm-preprocess/save_hycom/';

%% Grid of the nested model
eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);
% xc = XC(:,1); yc = YC(1,:); zc = RC;
% xf = XG(:,1); yf = YG(1,:);
% nxc = length(xc); nyc = length(yc);nzc = length(zc);

nPx = 1;
nPy = 1;

%% from 20100102 to 20100103 every day
nt1 = 8*1;

%% Position of the open boundaries (as defined in the nested model)
%% North.
obcsN = 'o';
nsn1 = 0;
nsn2 = 800;

% Corner islands
if strcmp(obcsN,'o') == 1
    OB_Jnorth = [zeros(nsn1,1); ones(nsn2,1); zeros(nxc - nsn1 - nsn2,1)]';
end

%% South.
obcsS = 'o';
nss1 = 0;
nss2 = 800;
if strcmp(obcsS,'o') == 1
    OB_Jsouth = [zeros(nss1,1); ones(nss2,1); zeros(nxc - nss1 - nss2,1)]';
end

%% West.
obcsW = 'o';
nsw1 = 0;
nsw2 = 448;
if strcmp(obcsW,'o') == 1
    OB_Iwest = [zeros(nsw1,1); ones(nsw2,1); zeros(nyc - nsw1 - nsw2,1)]';
end

%% East.
obcsE = 'o';
nse1 = 0;
nse2 = 448;
if strcmp(obcsE,'o') == 1
    OB_Ieast = [zeros(nse1,1); ones(nse2,1); zeros(nyc - nse1 - nse2,1)]';
end

%% Initialization of the arrays
if strcmp(obcsN,'o') == 1
    TobcsN = zeros(nxc,nPy,nzc,nt1);
    SobcsN = zeros(nxc,nPy,nzc,nt1);
    UobcsN = zeros(nxc,nPy,nzc,nt1);
    VobcsN = zeros(nxc,nPy,nzc,nt1);
end

if strcmp(obcsE,'o') == 1
    TobcsE = zeros(nyc,nPx,nzc,nt1);
    SobcsE = zeros(nyc,nPx,nzc,nt1);
    UobcsE = zeros(nyc,nPx,nzc,nt1);
    VobcsE = zeros(nyc,nPx,nzc,nt1);
end

if strcmp(obcsS,'o') == 1
    TobcsS = zeros(nxc,nPy,nzc,nt1);
    SobcsS = zeros(nxc,nPy,nzc,nt1);
    UobcsS = zeros(nxc,nPy,nzc,nt1);
    VobcsS = zeros(nxc,nPy,nzc,nt1);
end

if strcmp(obcsW,'o') == 1
    TobcsW = zeros(nyc,nPx,nzc,nt1);
    SobcsW = zeros(nyc,nPx,nzc,nt1);
    UobcsW = zeros(nyc,nPx,nzc,nt1);
    VobcsW = zeros(nyc,nPx,nzc,nt1);
end

if strcmp(obcsN,'o') == 1
    disp('NORTHERN BOUNDARY IS PROCESSING');
    %% loading northern boundary
    obcs = 'north';
    hycompath = hycompath1;
    
    %% 2009-2016 aevery 7 day frames 439 (nrec)
    nx1 =  813; ny1 = 13; nz1 = 41;
    
    HY_N_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
    HY_N_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
    HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
    time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
    
    %% NORTH BOUNDARY SLICE
    %% Loop on the time steps
    indt = 0;
    report('Time = xxxxx')
    for it = 1:nt1
        report('\b\b\b\b\b%5i',it)
        indt = indt+1;
        indx = [1:nxc];
        indy = [nyc-4:nyc];
        
        OB_Jnorth(OB_Jnorth ~= 0 ) = find(indy == nyc);
        
        xc1 = xc(indx); nxc1 = length(xc1);
        yc1 = yc(indy); nyc1 = length(yc1);
        xf1 = xf(indx); nxf1 = length(xf1);
        yf1 = yf(indy); nyf1 = length(yf1);
        zc1 = zc;  nzc1 = length(zc1);
        
        HY_N_U = rdslice([hycompath 'hycom_water_u_' obcs '.bin'],[nx1 ny1 nz1 ],it,fmt,Ieee);
        HY_N_V = rdslice([hycompath 'hycom_water_v_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_N_T = rdslice([hycompath 'hycom_water_temp_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_N_S = rdslice([hycompath 'hycom_salinity_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
       
        % interpolation 
        [Tm,xc1,yc1,zc1] = exinthycomobcs(fpath,'T',HY_N_T,HY_N_LON,HY_N_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Sm,xc1,yc1,zc1] = exinthycomobcs(fpath,'S',HY_N_S,HY_N_LON,HY_N_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Um,xc1,yc1,zc1] = exinthycomobcs(fpath,'U',HY_N_U,HY_N_LON,HY_N_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Vm,xc1,yc1,zc1] = exinthycomobcs(fpath,'V',HY_N_V,HY_N_LON,HY_N_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        
        Tin = reshape(Tm,nxc1,nyc1,nzc1);
        Sin = reshape(Sm,nxc1,nyc1,nzc1);
        Uin = reshape(Um,nxc1,nyc1,nzc1);
        Vin = reshape(Vm,nxc1,nyc1,nzc1);
        clear Tm Sm Um Vm
        
        if strcmp(obcsN,'o') == 1
            inorth = find(OB_Jnorth~=0);
            for ii = inorth
                for iproc = 1:nPy
                    TobcsN(ii,iproc,:,indt) = Tin(ii,OB_Jnorth(ii),:);
                    SobcsN(ii,iproc,:,indt) = Sin(ii,OB_Jnorth(ii),:);
                    UobcsN(ii,iproc,:,indt) = Uin(ii,OB_Jnorth(ii),:);
                    VobcsN(ii,iproc,:,indt) = Vin(ii,OB_Jnorth(ii),:);
                end
            end
        end
        clear Tin Sin Uin Vin
    end
    
    clear  HY_N_LON HY_N_LAT HY_DEPTH time_hycom HY_N_U HY_N_V HY_N_T HY_N_S hycompath
    clear xc1 yc1 xf1 yf1 zc1 indx indy nxc1 nxf1 nyc1 nyf1 nzc1
    disp('NORTHERN BOUNDARY IS DONE');
end

if strcmp(obcsS,'o') == 1
    disp('SOUTHERN BOUNDARY IS PROCESSING');
    %% loading southern boundary
    obcs = 'south';
    hycompath = hycompath1;
    
    %% 2009-2016 aevery 7 day frames 439 (nrec)
    nx1 =  813; ny1 = 13; nz1 = 41;
    
    HY_S_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
    HY_S_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
    HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
    time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
    
    %% SOUTH BOUNDARY SLICE
    %% Loop on the time steps
    indt = 0;
    report('Time = xxxxx')
    for it = 1:nt1
        report('\b\b\b\b\b%5i',it)
        indt = indt+1;
        indx = [1:nxc];
        indy = [1:5];
        
        disp('EASTERN BOUNDARY LOOP...');
        OB_Jsouth(OB_Jsouth ~= 0 ) = find(indy == 1);
        
        xc1 = xc(indx); nxc1 = length(xc1);
        yc1 = yc(indy); nyc1 = length(yc1);
        xf1 = xf(indx); nxf1 = length(xf1);
        yf1 = yf(indy); nyf1 = length(yf1);
        zc1 = zc;  nzc1 = length(zc1);
        
        HY_S_U = rdslice([hycompath 'hycom_water_u_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_S_V = rdslice([hycompath 'hycom_water_v_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_S_T = rdslice([hycompath 'hycom_water_temp_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_S_S = rdslice([hycompath 'hycom_salinity_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        
        [Tm,xc1,yc1,zc1] = exinthycomobcs(fpath,'T',HY_S_T,HY_S_LON,HY_S_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Sm,xc1,yc1,zc1] = exinthycomobcs(fpath,'S',HY_S_S,HY_S_LON,HY_S_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Um,xc1,yc1,zc1] = exinthycomobcs(fpath,'U',HY_S_U,HY_S_LON,HY_S_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Vm,xc1,yc1,zc1] = exinthycomobcs(fpath,'V',HY_S_V,HY_S_LON,HY_S_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        
        nzc = length(zc1);
        
        Tin = reshape(Tm,nxc1,nyc1,nzc1);
        Sin = reshape(Sm,nxc1,nyc1,nzc1);
        Uin = reshape(Um,nxc1,nyc1,nzc1);
        Vin = reshape(Vm,nxc1,nyc1,nzc1);
        clear Tm Sm Um Vm
        
        if strcmp(obcsS,'o') == 1
            isouth = find(OB_Jsouth~=0);
            for ii = isouth
                for iproc = 1:nPy
                    TobcsS(ii,iproc,:,indt) = Tin(ii,OB_Jsouth(ii)  ,:);
                    SobcsS(ii,iproc,:,indt) = Sin(ii,OB_Jsouth(ii)  ,:);
                    UobcsS(ii,iproc,:,indt) = Uin(ii,OB_Jsouth(ii)  ,:);
                    VobcsS(ii,iproc,:,indt) = Vin(ii,OB_Jsouth(ii)+1,:);
                end
            end
        end
        clear Tin Sin Uin Vin
    end
    clear  HY_S_LON HY_S_LAT HY_DEPTH time_hycom HY_S_U HY_S_V HY_S_T HY_S_S hycompath
    clear xc1 yc1 xf1 yf1 zc1 indx indy nxc1 nxf1 nyc1 nyf1 nzc1
    disp('SOUTHERN BOUNDARY IS DONE');
end

if strcmp(obcsW,'o') == 1
    disp('WESTERN BOUNDARY IS PROCESSING');
    %% loading western boundary
    obcs = 'west';
    hycompath = hycompath1;
    
    %% 2009-2016 aevery 7 day frames 439 (nrec)
    nx1 =  13; ny1 = 637; nz1 = 41;
    
    HY_W_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
    HY_W_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
    HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
    time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
    
    %% WEST BOUNDARY SLICE
    %% Loop on the time steps
    indt = 0;
    report('Time = xxxxx')
    for it = 1:nt1
        report('\b\b\b\b\b%5i',it)
        indt = indt+1;
        indx = [1:5];
        indy = [1:nyc];
        
        OB_Iwest(OB_Iwest ~= 0 ) = find(indx == 1);
        
        xc1 = xc(indx); nxc1 = length(xc1);
        yc1 = yc(indy); nyc1 = length(yc1);
        xf1 = xf(indx); nxf1 = length(xf1);
        yf1 = yf(indy); nyf1 = length(yf1);
        zc1 = zc;  nzc1 = length(zc1);
        
        HY_W_U = rdslice([hycompath 'hycom_water_u_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_W_V = rdslice([hycompath 'hycom_water_v_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_W_T = rdslice([hycompath 'hycom_water_temp_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_W_S = rdslice([hycompath 'hycom_salinity_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        
        [Tm,xc1,yc1,zc1] = exinthycomobcs(fpath,'T',HY_W_T,HY_W_LON,HY_W_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Sm,xc1,yc1,zc1] = exinthycomobcs(fpath,'S',HY_W_S,HY_W_LON,HY_W_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Um,xc1,yc1,zc1] = exinthycomobcs(fpath,'U',HY_W_U,HY_W_LON,HY_W_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Vm,xc1,yc1,zc1] = exinthycomobcs(fpath,'V',HY_W_V,HY_W_LON,HY_W_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        
        nzc = length(zc1);
        
        Tin = reshape(Tm,nxc1,nyc1,nzc1);
        Sin = reshape(Sm,nxc1,nyc1,nzc1);
        Uin = reshape(Um,nxc1,nyc1,nzc1);
        Vin = reshape(Vm,nxc1,nyc1,nzc1);
        clear Tm Sm Um Vm
        
        if strcmp(obcsW,'o') == 1
            jwest = find(OB_Iwest~=0);
            for jj = jwest
                for iproc = 1:nPx
                    TobcsW(jj,iproc,:,indt) = Tin(OB_Iwest(jj)  ,jj,:);
                    SobcsW(jj,iproc,:,indt) = Sin(OB_Iwest(jj)  ,jj,:);
                    UobcsW(jj,iproc,:,indt) = Uin(OB_Iwest(jj)+1,jj,:);
                    VobcsW(jj,iproc,:,indt) = Vin(OB_Iwest(jj)  ,jj,:);
                end
            end
        end
        clear Tin Sin Uin Vin
    end
    clear  HY_W_LON HY_W_LAT HY_DEPTH time_hycom HY_W_U HY_W_V HY_W_T HY_W_S hycompath
    clear xc1 yc1 xf1 yf1 zc1 indx indy nxc1 nxf1 nyc1 nyf1 nzc1
    disp('WESTERN BOUNDARY IS DONE');
end

if strcmp(obcsE,'o') == 1
    disp('EASTERN BOUNDARY IS PROCESSING');
    %% loading eastern boundary
    obcs = 'east';
    hycompath = hycompath1;
    
    %% 2009-2016 aevery 7 day frames 439 (nrec)
    nx1 =  13; ny1 = 637; nz1 = 41;
    
    HY_E_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
    HY_E_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
    HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
    time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
    
    %% EAST BOUNDARY SLICE
    %% Loop on the time steps
    indt = 0;
    report('Time = xxxxx')
    for it = 1:nt1
        report('\b\b\b\b\b%5i',it)
        disp('EASTERN BOUNDARY LOOP...');
        indt = indt+1;
        indx = [nxc-4:nxc];
        indy = [1:nyc];
        
        disp('EASTERN BOUNDARY LOOP 0.5...');
        OB_Ieast(OB_Ieast ~= 0 ) = find(indx == nxc);
        
        disp('EASTERN BOUNDARY LOOP 0.8...');
        xc1 = xc(indx); nxc1 = length(xc1);
        yc1 = yc(indy); nyc1 = length(yc1);
        xf1 = xf(indx); nxf1 = length(xf1);
        yf1 = yf(indy); nyf1 = length(yf1);
        zc1 = zc;  nzc1 = length(zc1);
        
        disp('EASTERN BOUNDARY LOOP 1...');
        HY_E_U = rdslice([hycompath 'hycom_water_u_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_E_V = rdslice([hycompath 'hycom_water_v_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_E_T = rdslice([hycompath 'hycom_water_temp_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        HY_E_S = rdslice([hycompath 'hycom_salinity_' obcs '.bin'],[nx1 ny1 nz1],it,fmt,Ieee);
        
        disp('EASTERN BOUNDARY LOOP 2...');
        [Tm,xc1,yc1,zc1] = exinthycomobcs(fpath,'T',HY_E_T,HY_E_LON,HY_E_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Sm,xc1,yc1,zc1] = exinthycomobcs(fpath,'S',HY_E_S,HY_E_LON,HY_E_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Um,xc1,yc1,zc1] = exinthycomobcs(fpath,'U',HY_E_U,HY_E_LON,HY_E_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        [Vm,xc1,yc1,zc1] = exinthycomobcs(fpath,'V',HY_E_V,HY_E_LON,HY_E_LAT,HY_DEPTH,xc1,yc1,xf1,yf1,indx,indy,zc1);
        
        nzc = length(zc1);
        
        Tin = reshape(Tm,nxc1,nyc1,nzc1);
        Sin = reshape(Sm,nxc1,nyc1,nzc1);
        Uin = reshape(Um,nxc1,nyc1,nzc1);
        Vin = reshape(Vm,nxc1,nyc1,nzc1);
        clear Tm Sm Um Vm
        
        disp('');
        disp('EASTERN BOUNDARY LOOP 3...');
        if strcmp(obcsE,'o') == 1
            jeast = find(OB_Ieast~=0);
            for jj = jeast
                for iproc = 1:nPx
                    TobcsE(jj,iproc,:,indt) = Tin(OB_Ieast(jj)  ,jj,:);
                    SobcsE(jj,iproc,:,indt) = Sin(OB_Ieast(jj)  ,jj,:);
                    UobcsE(jj,iproc,:,indt) = Uin(OB_Ieast(jj)  ,jj,:);
                    VobcsE(jj,iproc,:,indt) = Vin(OB_Ieast(jj)  ,jj,:);
                end
            end
        end
        clear Tin Sin Uin Vin
    end
    clear  HY_E_LON HY_E_LAT HY_DEPTH time_hycom HY_E_U HY_E_V HY_E_T HY_E_S hycompath
    clear xc1 yc1 xf1 yf1 zc1 indx indy nxc1 nxf1 nyc1 nyf1 nzc1
    disp('EASTERN BOUNDARY IS DONE');
end


%% Save outputs

if strcmp(obcsN,'o') == 1
    wrslice([fpath 'Rs_TobcsN.bin' ],TobcsN,1,fmt,Ieee);
    wrslice([fpath 'Rs_SobcsN.bin' ],SobcsN,1,fmt,Ieee);
    wrslice([fpath 'Rs_UobcsN.bin' ],UobcsN,1,fmt,Ieee);
    wrslice([fpath 'Rs_VobcsN.bin' ],VobcsN,1,fmt,Ieee);
end

if strcmp(obcsS,'o') == 1
    wrslice([fpath 'Rs_TobcsS.bin' ],TobcsS,1,fmt,Ieee);
    wrslice([fpath 'Rs_SobcsS.bin' ],SobcsS,1,fmt,Ieee);
    wrslice([fpath 'Rs_UobcsS.bin' ],UobcsS,1,fmt,Ieee);
    wrslice([fpath 'Rs_VobcsS.bin' ],VobcsS,1,fmt,Ieee);
end

if strcmp(obcsE,'o') == 1
    wrslice([fpath 'Rs_TobcsE.bin' ],TobcsE,1,fmt,Ieee);
    wrslice([fpath 'Rs_SobcsE.bin' ],SobcsE,1,fmt,Ieee);
    wrslice([fpath 'Rs_UobcsE.bin' ],UobcsE,1,fmt,Ieee);
    wrslice([fpath 'Rs_VobcsE.bin' ],VobcsE,1,fmt,Ieee);
end

if strcmp(obcsW,'o') == 1
    wrslice([fpath 'Rs_TobcsW.bin' ],TobcsW,1,fmt,Ieee);
    wrslice([fpath 'Rs_SobcsW.bin' ],SobcsW,1,fmt,Ieee);
    wrslice([fpath 'Rs_UobcsW.bin' ],UobcsW,1,fmt,Ieee);
    wrslice([fpath 'Rs_VobcsW.bin' ],VobcsW,1,fmt,Ieee);
end

mean(TobcsE(:))
mean(TobcsS(:))
mean(SobcsE(:))
mean(SobcsS(:))

disp('BOUNDARY BIN FILES (UNCORRECTED) ARE GENERATED');
