%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TO CALCULATE VOLUME TRANSPORT AND HEAT TRANSPORT ACCROSS OPEN BOUNDARIES
%% FROM GLOBAL MODEL (HYCOM MODEL USED)
%% THIS SCRIPT ONLY CONSIDERS EAST/WEST BOUNDARY SLICES
%% ZONAL TRANSPORT
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
%% READING WESTERN BOUNDARY
%% loading western boundary
obcs = 'west';

hycompath = [hycompath1 '_' obcs '/'];

%% 20180101 to 20180302 every day frames 61 (nt1)
nx1 =  13; ny1 = 270; nz1 = 41; nt1 = 61;
HY_W_LON = rdslice([hycompath 'lon_hycom_' obcs '.bin'],[nx1 1],1,fmt,Ieee);
HY_W_LAT = rdslice([hycompath 'lat_hycom_' obcs '.bin'],[ny1 1],1,fmt,Ieee);
HY_DEPTH = rdslice([hycompath 'depth_hycom_' obcs '.bin'],[nz1 1],1,fmt,Ieee);
time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
HY_W_U = rdslice([hycompath 'hycom_u_' obcs '.bin'],[nx1 ny1 nz1 nt1],1,fmt,Ieee);
HY_W_T = rdslice([hycompath 'hycom_temperature_' obcs '.bin'],[nx1 ny1 nz1 nt1],1,fmt,Ieee);
date = ['ccs_' datestr(time_hycom(1),7) datestr(time_hycom(1),28)];

%% HYCOM GRID INFO
x = HY_W_LON;   a = length(HY_W_LON);
y = HY_W_LAT;   b = length(HY_W_LAT);

dz = diff(-HY_DEPTH);
z = cumsum(dz)-dz/2;    c = length(HY_DEPTH);

%% ADDING AN ADDITIONAL BOTTOM LAYER,AS HYCOM STARTS WITH '0'
dz = [dz; dz(end)];
%% HYCOM solution outputs.
hycomtfile = HY_W_T;
hycomufile = HY_W_U;

%% Lats and longs where to compute HYCOM transports.
%% this is with respect to MITGCM model boundaries
%% dlat = dlon = 0.01 degrees for GOM GRID
%% (EAST, WEST Open Boundaries)
%% WEST BOUND
%% tablon_W = [round(xc(1))];
tablon_W = [(xc(1))+0.5*deltaX];

%% ###################################################################
%% GETTING THE AREAS/VOLUMES
%% ###################################################################
%% getting dy along HYCOM latitude
%% getting x increments
DX = diff(x);
%% getting y increments
DY = diff(y);
%% converting into meters (1 deg = 111111 m)
dx = [DX; DX(end)].*111111;
%% converting into meters (1 deg = 111111 m)
dy = [DY; DY(end)].*ones(b,1)*111111;

dV  = zeros(a,b,c);
dVT = zeros(a,b,c);
dA  = zeros(a,b,c);

%% COMPUTING SLICE AREAS ALONG LATITUDE
for i = 1:a
    for j = 1:b
        for k = 1:c
            dA(i,j,k) = dy(j)*dz(k);
        end
    end
end

%% computing total VOLUME of hycom grid
for i = 1:a
    for j = 1:b
        for k = 1:c
            %% dx(long) is corrected for each dy (lat) by cos(lat) factor
            dV(i,j,k) = dx(i)*cos(y(j)*pi/180)*dy(j)*dz(k);
        end
    end
end

%% MASKING HYCOM VALUES
%% missing values in HYCOM are assigned to '0' during extraction
ssh = squeeze(hycomufile(:,:,:,1));
for k = 1:c
    [xx,yy] = find(ssh(:,:,k) == 0);
    nn = length(xx);
    for ii = 1:nn
        dA(xx(ii),yy(ii),k) = 0;
    end
end

%% duiplicating for volume tranport
dVT = dV;

%% ###################################################################
%% GETTING THE model bounds
%% ###################################################################
phi1 = tablon_W;
%% lom1 = round(yc(1));
%% lom2 = round(yc(end));
lom1 = (yc(1)-0.5*deltaY);
lom2 = (yc(end)+0.5*deltaY);

%% identifying east extreme
j1 = find(y - lom1 > 0.0001);
j2 = find(y - lom2 > 0.0001);
i1 = find(x - phi1 < 0.0001);
i1 = i1(end);

%% plotting to check the HYCOM raw analysis
figure
orient tall
subplot(1,2,1);
dummy = squeeze(hycomufile(i1,j1:j2,:,1))';
dummy(dummy == 0) = NaN;
contourf(y(j1:j2),HY_DEPTH',dummy)
colorbar
set(gca,'FontSize',14);
title(['Zonal Velocity at ' num2str(round(tablon_W)) '^oE'],'FontSize',14);
xlabel('Latitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

subplot(1,2,2)
dummy = squeeze(hycomtfile(i1,j1:j2,:,1))';
dummy(dummy == 0) = NaN;
contourf(y(j1:j2),HY_DEPTH',dummy)
colorbar
set(gca,'FontSize',14);
title(['Zonal Temp at ' num2str(round(tablon_W)) '^oE'],'FontSize',14);
eval(['print -depsc ' fpath 'HYCOM_VEL_TEMP_ZON_' date '_' num2str(round(tablon_W)) 'E']);
xlabel('Latitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l = 1:nt1
    %% Set T.for temp
    %% MAKE VOLUME WEIGHTED Averaged
    T = zeros(a,b,c);
    ssh = squeeze(hycomtfile(:,:,:,l));
    ssh = ssh.*dVT;
    for j = j1:j2
        for k = 1:c
            i = i1-1;
            if((dVT(i,j,k)+dVT(i+1,j,k))>0)
                T(i+1,j,k) = (ssh(i,j,k)+ssh(i+1,j,k))/(dVT(i,j,k)+dVT(i+1,j,k));
            end
            if(ssh(i,j,k)==0 & dVT(i+1,j,k)>0)
                T(i+1,j,k) =  (ssh(i+1,j,k))/(dVT(i+1,j,k));
            end
            if(ssh(i+1,j,k)==0 & dVT(i,j,k)>0)
                T(i+1,j,k) =  (ssh(i,j,k))/(dVT(i,j,k));
            end
        end
    end
    %% Set U.for velocity
    U = squeeze(hycomufile(:,:,:,l));
    %% Compute transports.
    %% TEMP TRSNPORT
    UT_W(l,1) = 1024*4000*sum(sum(T(i1,j1:j2,:).*U(i1,j1:j2,:).*dA(i1,j1:j2,:)));
    %% VOLUME TRANSPORT
    Utrans_W(l,1) = sum(sum(U(i1,j1:j2,:).*dA(i1,j1:j2,:)));
end

%% Plotting the Zonal Trnaports
figure
orient tall
subplot(2,1,1);
set(gca,'FontSize',14);
plot([1:nt1],UT_W(:,1)*1e-15,'b','LineWidth',2);
grid on;
title(['Zonal Heat Transport (x 10^{15} W) - Mean = ' ...
    num2str(mean(UT_W(:,1))*1e-15)],'FontSize',14);
subplot(2,1,2);
grid on;
set(gca,'FontSize',14);
plot([1:nt1],Utrans_W(:,1)*1e-6,'b','LineWidth',2);
grid on;
title(['Zonal Transport (Sv) - Mean = ' ...
    num2str(mean(Utrans_W(:,1))*1e-6)],'FontSize',14);
suptitle(['HYCOM - ' num2str(round(tablon_W)) '^oE']);
eval(['print -depsc ' fpath 'HYCOM_HF_ZT_' date '_' num2str(round(tablon_W)) 'E']);

clear DX DY dx dy dV dVT dA AS a b c x y dz z hycomtfile hycomufile hycomvfile
clear  HY_TIME HY_DEPTH HY_W_LON HY_W_LAT HY_W_U HY_W_T HY_W_ETA hycompath

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% saving for volume adjusment in MITgcm later
utot_ZON = [UT_W Utrans_W];
eval(['save ' fpath 'HYCOM_HF_ZT_' date ' utot_ZON -ascii']);
