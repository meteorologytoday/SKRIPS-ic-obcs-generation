%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO CALCULATE VOLUME TRANSPORT AND HEAT TRANSPORT ACCROSS OPEN BOUNDARIES
%% FROM GLOBAL MODEL (HYCOM MODEL USED)
%% THIS SCRIPT ONLY CONSIDERS NORTH/SOUTH BOUNDARY SLICES
%% MERIDIAONAL TRANSPORT
clc
clear all
close all

%% Initialization
fpath = '/project/rus043_data/ar_2018_more_vertical/script-pre-ar2018/ic_obcs_20180123_sponge13/';
hycompath1 = '/home/rus043/HYCOM/2018_ar';
deltaX = 0.08;
deltaY = 0.08;

%% Grid of the nested model
eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);
% xc = XC(:,1); yc = YC(1,:); zc = RC;
% xf = XG(:,1); yf = YG(1,:);
% nxc = length(xc); nyc = length(yc);nzc = length(zc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READING NORTHERN BOUNDARY
%% loading northern boundary
obcs = 'north';

hycompath = [hycompath1 '_' obcs '/'];

%% 20180101 to 20180301 every day frames 28 (nt1)
nx1 =  445; ny1 = 13; nz1 = 41; nt1 = 61;
HY_N_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
HY_N_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
HY_N_V = rdslice([hycompath 'hycom_v_' obcs '.bin'],[nx1 ny1 nz1 nt1],1,fmt,Ieee);
HY_N_T = rdslice([hycompath 'hycom_temperature_' obcs '.bin'],[nx1 ny1 nz1 nt1],1,fmt,Ieee);
date = ['ccs_' datestr(time_hycom(1),7) datestr(time_hycom(1),28)];

%% HYCOM GRID INFO
x = HY_N_LON; a = length(HY_N_LON);
y = HY_N_LAT; b = length(HY_N_LAT);

dz = diff(-HY_DEPTH);
z = cumsum(dz)-dz/2; c = length(HY_DEPTH);
%% ADDING AN ADDITIONAL BOTTOM LAYER,AS HYCOM STARTS WITH '0'
dz = [dz; dz(end)];

%% HYCOM solution outputs.
hycomtfile = HY_N_T;
hycomvfile = HY_N_V;

%% Lats and longs where to compute HYCOM transports.
%% this is with respect to MITGCM model boundaries
%% dlat = dlon = 0.01 degrees for GOM GRID
%% (NORTH, SOUTH  Open Boundaries)
%% NORTH BOUND
%% tablat_N = [round(yc(end))];
tablat_N = [(yc(end))+0.5*deltaY];

%% ###################################################################
%% GETTING THE AREAS/VOLUMES
%% ###################################################################
%% getting x increments
DX = diff(x);
%% getting y increments
DY = diff(y);
%% converting into meters (1 deg = 111111 m)
dx = [DX; DX(end)].*111111;
dy = [DY; DY(end)].*111111;
dV  = zeros(a,b,c);
dVT = zeros(a,b,c);
dA  = zeros(a,b,c);

%% computing total VOLUME of hycom grid
for i = 1:a
    for j = 1:b
        for k = 1:c
            %% dx(long) is corrected for each dy (lat) by cos(lat) factor
            dV(i,j,k) = dx(i)*cos(y(j)*pi/180)*dy(j)*dz(k);
        end
    end
end

%%  area SLICES are computed along Longitudes
for i = 1:a
    for j = 1:b
        for k = 1:c
            %% dx(long) is corrected for each dy (lat) by cos(lat) factor
            dA(i,j,k) = dx(i)*cos(y(j)*pi/180)*dz(k);
        end
    end
end

%% SURFACE AREA
for i = 1:a
    for j = 1:b
        %% dx(long) is corrected for each dy (lat) by cos(lat) factor
        AS(i,j) = dx(i)*cos(y(j)*pi/180)*dy(j);
    end
end

%% TO GET LAND MASK OF HYCOM (vel file)
%% missing values in HYCOM are assigned to '0' during extraction
ssh = squeeze(hycomvfile(:,:,:,1));
for k = 1:c
    [xx,yy] = find(squeeze(ssh(:,:,k)) == 0);
    nn = length(xx);
    for ii = 1:nn
        dA(xx(ii),yy(ii),k) = 0;
    end
end
%% duiplicating for volume tranport
dVT = dV;

%% TO GET LAND MASK OF HYCOM (temp file)
%% missing values in HYCOM are assigned to '0' during extraction
ssh = squeeze(hycomtfile(:,:,:,1));
for k = 1:c
    [xx,yy] = find(squeeze(ssh(:,:,k)) == 0);
    nn = length(xx);
    for ii = 1:nn
        dVT(xx(ii),yy(ii),k)=0;
    end
end

%% getting MITGCM model lats to fix north boundary
%% NORTHERN BOUNDARY >> tablat_N
phi1 = tablat_N;
% lam1 = round(xc(1));
% lam2 = round(xc(end));
lam1 = (xc(1)-0.5*deltaX);
lam2 = (xc(end)+0.5*deltaY);

%% identifying North extremes
i1 = find(x-lam1 > 0.0001);
i2 = find(x-lam2 > 0.0001);
j1 = find(y-phi1 < 0.0001);
j1 = j1(end);

%% plotting to check the HYCOM raw analysis
figure
orient tall
subplot(1,2,1);
dummy = squeeze(hycomvfile(i1:i2,j1,:,1))';
dummy(dummy == 0) = NaN;
contourf(x(i1:i2),HY_DEPTH',dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridional Velocity at ' num2str(round(tablat_N)) '^oN'],'FontSize',14);
xlabel('Longitude (deg) ','FontSize',14); ylabel('Depth (m)','FontSize',14)

subplot(1,2,2)
dummy = squeeze(hycomtfile(i1:i2,j1,:,1))';
dummy(dummy == 0) = NaN;
contourf(x(i1:i2),HY_DEPTH',dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridional Temp at ' num2str(round(tablat_N)) '^oN'],'FontSize',14);
eval(['print -depsc ' fpath 'HYCOM_VEL_TEMP_MER_' date '_' num2str(round(tablat_N)) 'N']);
xlabel('Longitude (deg) ','FontSize',14); ylabel('Depth (m)','FontSize',14)

for l = 1:nt1
    %% Set T as temp
    T = zeros(a,b,c);
    ssh = squeeze(hycomtfile(:,:,:,l));
    %% MAKE VOLUME WEIGHTED Averaged
    ssh = ssh.*dVT;
    for i = i1:i2
        for k = 1:c
            %% getting adjacent y
            j = j1-1;
            if((dVT(i,j,k)+dVT(i,j+1,k))>0)
                T(i,j+1,k) = (ssh(i,j,k)+ssh(i,j+1,k))/(dVT(i,j,k)+dVT(i,j+1,k));
            end
            if(ssh(i,j,k)==0 & dVT(i,j+1,k)>0)
                T(i,j+1,k) =  (ssh(i,j+1,k))/(dVT(i,j+1,k));
            end
            if(ssh(i,j+1,k)==0 & dVT(i,j,k)>0)
                T(i,j+1,k) =  (ssh(i,j,k))/(dVT(i,j,k));
            end
        end
    end
    %% Set V. as Velocity
    V = squeeze(hycomvfile(:,:,:,l));
    %% Compute transports.
    %% TEMP TRANPOSRT
    VT_N(l,1) = 1024*4000*sum(sum(T(i1:i2,j1,:).*V(i1:i2,j1,:).*dA(i1:i2,j1,:)));
    %% Volmue Transport
    Vtrans_N(l,1) = sum(sum(V(i1:i2,j1,:).*dA(i1:i2,j1,:)));
end

%% Plotting the transports
figure
orient tall
subplot(2,1,1);
set(gca,'FontSize',14);
plot([1:nt1],VT_N(:,1)*1e-15,'b','LineWidth',2);
grid on;
title(['Meridional Heat Transport (x 10^{15} W) - Mean = ' ...
    num2str(mean(VT_N(:,1))*1e-15)],'FontSize',14);

subplot(2,1,2);
set(gca,'FontSize',14);
plot([1:nt1],Vtrans_N(:,1)*1e-6,'b','LineWidth',2);
grid on;
title(['Meridional Transport (Sv) - Mean = ' ...
    num2str(mean(Vtrans_N(:,1))*1e-6)],'FontSize',14);
suptitle(['HYCOM - ' num2str(round(tablat_N)) '^oN']);
eval(['print -depsc ' fpath 'HYCOM_HF_MT_' date '_' num2str(round(tablat_N)) 'N']);

clear DX DY dx dy dV dVT dA AS a b c x y dz z hycomtfile hycomufile hycomvfile
clear  time_hycom HY_DEPTH HY_LON HY_LAT HY_N_LON HY_N_LAT HY_N_V HY_N_T HY_N_ETA hycompath

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% READING SOUTHERN BOUNDARY SLICE

%% loading southern boundary
obcs = 'south';
hycompath = [hycompath1 '_' obcs '/'];

%% 20180101 to 20180301 every day frames 28 (nt1)
nx1 =  445; ny1 = 13; nz1 = 41; nt1 = 61;
HY_S_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
HY_S_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
HY_S_V = rdslice([hycompath 'hycom_v_' obcs '.bin'],[nx1 ny1 nz1 nt1],1,fmt,Ieee);
HY_S_T = rdslice([hycompath 'hycom_temperature_' obcs '.bin'],[nx1 ny1 nz1 nt1],1,fmt,Ieee);
date = ['ccs_' datestr(time_hycom(1),7) datestr(time_hycom(1),28)];

%% HYCOM GRID INFO
x = HY_S_LON;   a = length(HY_S_LON);
y = HY_S_LAT;   b = length(HY_S_LAT);

dz = diff(-HY_DEPTH);
z = cumsum(dz)-dz/2;    c = length(HY_DEPTH);
%% ADDING AN ADDITIONAL BOTTOM LAYER,AS HYCOM STARTS WITH '0'
dz = [dz; dz(end)];

%% HYCOM solution outputs.
hycomtfile = HY_S_T;
hycomvfile = HY_S_V;

%% Lats and longs where to compute HYCOM transports.
%% this is with respect to MITGCM model boundaries
%% dlat = dlon = 0.01 degrees for GOM GRID
%% (NORTH, SOUTH  Open Boundaries)
%% south BOUND
%% tablat_S = [round(yc(1))];
tablat_S = [(yc(1))-0.5*deltaY];

%% ###################################################################
%% GETTING THE AREAS/VOLUMES
%% ###################################################################
%% getting x increments
DX = diff(x);
%% getting y increments
DY = diff(y);
%% converting into meters (1 deg = 111111 m)
dx = [DX; DX(end)].*111111;
dy = [DY; DY(end)].*111111;
dV  = zeros(a,b,c);
dVT = zeros(a,b,c);
dA  = zeros(a,b,c);
%% computing total VOLUME of hycom grid
for i = 1:a
    for j = 1:b
        for k = 1:c
            %% dx(long) is corrected for each dy (lat) by cos(lat) factor
            dV(i,j,k) = dx(i)*cos(y(j)*pi/180)*dy(j)*dz(k);
        end
    end
end
%%  area SLICES are computed along Longitudes
for i = 1:a
    for j = 1:b
        for k = 1:c
            %% dx(long) is corrected for each dy (lat) by cos(lat) factor
            dA(i,j,k) = dx(i)*cos(y(j)*pi/180)*dz(k);
        end
    end
end
%% SURFACE AREA
for i = 1:a
    for j = 1:b
        %% dx(long) is corrected for each dy (lat) by cos(lat) factor
        AS(i,j) = dx(i)*cos(y(j)*pi/180)*dy(j);
    end
end

%% TO GET LAND MASK OF HYCOM (vel file)
%% missing values in HYCOM are assigned to '0' during extraction
ssh = squeeze(hycomvfile(:,:,:,1));
for k = 1:c
    [xx,yy] = find(squeeze(ssh(:,:,k)) == 0);
    nn = length(xx);
    for ii = 1:nn
        dA(xx(ii),yy(ii),k) = 0;
    end
end
%% duiplicating for volume tranport
dVT = dV;

%% TO GET LAND MASK OF HYCOM (temp file)
%% missing values in HYCOM are assigned to '0' during extraction
ssh = squeeze(hycomtfile(:,:,:,1));
for k = 1:c
    [xx,yy] = find(squeeze(ssh(:,:,k)) == 0);
    nn = length(xx);
    for ii = 1:nn
        dVT(xx(ii),yy(ii),k)=0;
    end
end

%% getting MITGCM model lats to fix north boundary
%% NORTHERN BOUNDARY >> tablat_N
phi1 = tablat_S;
%% dlon = 0.1 dlat = 0.1
%% lam1 = round(xc(1));
lam1 = (xc(1)-0.5*deltaX);
%% dlon = 0.1, dlat = 0.1
%% lam2 = round(xc(end));
lam2 = (xc(end)+0.5*deltaX);
%% identifying North extremes
i1 = find(x-lam1 > 0.0001);
i2 = find(x-lam2 > 0.0001);
j1 = find(y-phi1 < 0.0001);
j1 = j1(end);

%% plotting to check the HYCOM raw analysis
figure
orient tall
subplot(1,2,1);
dummy = squeeze(hycomvfile(i1:i2,j1,:,1))';
dummy(dummy == 0) = NaN;
contourf(x(i1:i2),HY_DEPTH',dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridional Velocity at ' num2str(round(tablat_S)) '^oN'],'FontSize',14);
xlabel('Longitude (deg) ','FontSize',14); ylabel('Depth (m)','FontSize',14)

subplot(1,2,2)
dummy = squeeze(hycomtfile(i1:i2,j1,:,1))';
dummy(dummy == 0) = NaN;
contourf(x(i1:i2),HY_DEPTH',dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridional Temp at ' num2str(round(tablat_S)) '^oN'],'FontSize',14);
eval(['print -depsc ' fpath 'HYCOM_VEL_TEMP_MER_' date '_' num2str(round(tablat_N)) 'N']);
xlabel('Longitude (deg) ','FontSize',14); ylabel('Depth (m)','FontSize',14)

for l = 1:nt1
    %% Set T as temp
    T = zeros(a,b,c);
    ssh = squeeze(hycomtfile(:,:,:,l));
    %% MAKE VOLUME WEIGHTED Averaged
    ssh = ssh.*dVT;
    for i = i1:i2
        for k = 1:c
            %% getting adjacent y
            j = j1-1;
            if((dVT(i,j,k)+dVT(i,j+1,k))>0)
                T(i,j+1,k) = (ssh(i,j,k)+ssh(i,j+1,k))/(dVT(i,j,k)+dVT(i,j+1,k));
            end
            if(ssh(i,j,k)==0 & dVT(i,j+1,k)>0)
                T(i,j+1,k) =  (ssh(i,j+1,k))/(dVT(i,j+1,k));
            end
            if(ssh(i,j+1,k)==0 & dVT(i,j,k)>0)
                T(i,j+1,k) =  (ssh(i,j,k))/(dVT(i,j,k));
            end
        end
    end
    %% Set V. as Velocity
    V = squeeze(hycomvfile(:,:,:,l));
    %% Compute transports.
    %% TEMP TRANPOSRT
    VT_S(l,1) = 1024*4000*sum(sum(T(i1:i2,j1,:).*V(i1:i2,j1,:).*dA(i1:i2,j1,:)));
    %% Volmue Transport
    Vtrans_S(l,1) = sum(sum(V(i1:i2,j1,:).*dA(i1:i2,j1,:)));
end

%% Plotting the transports
figure
orient tall
subplot(2,1,1);
set(gca,'FontSize',14);
plot([1:nt1],VT_S(:,1)*1e-15,'*b','LineWidth',2);
grid on;
title(['Meridional Heat Transport (x 10^{15} W) - Mean = ' ...
    num2str(mean(VT_S(:,1))*1e-15)],'FontSize',14);

subplot(2,1,2);
set(gca,'FontSize',14);
plot([1:nt1],Vtrans_S(:,1)*1e-6,'*b','LineWidth',2);
grid on;
title(['Meridional Transport (Sv) - Mean = ' ...
    num2str(mean(Vtrans_S(:,1))*1e-6)],'FontSize',14);
suptitle(['HYCOM - ' num2str(round(tablat_S)) '^oN']);
eval(['print -depsc ' fpath 'HYCOM_HF_MT_' date '_' num2str(round(tablat_S)) 'N']);

clear DX DY dx dy dV dVT dA AS a b c x y dz z hycomtfile hycomufile hycomvfile
clear  time_hycom HY_LON HY_LAT HY_DEPTH HY_S_LON HY_S_LAT HY_S_V HY_S_T HY_S_ETA

%% saving for volume adjusment in MITgcm later
vtot_MER = [VT_N VT_S Vtrans_N Vtrans_S];
eval(['save ' fpath 'HYCOM_HF_MT_' date ' vtot_MER -ascii']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
