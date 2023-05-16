%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correct the velocity fields: for MITGCM MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc

run ~/matlab_bin/pathdef.m
fpath = '/project/rus043_data/ar_2018_more_vertical/script-pre-ar2018/ic_obcs_20180123_sponge13/';
hycompath1 = '/home/rus043/HYCOM/2018_ar';

%% Grid of the nested model
eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);

%% READING SOUTHHERN BOUNDARY SLICE
%% 2017 Feb aevery 1 day frames 28 (nt1)
nt1 = 61;
obcs = 'south';
hycompath = [hycompath1 '_' obcs '/'];

time_hycom = rdslice([hycompath 'time_hycom_' obcs '.bin'],[nt1 1],1,fmt,Ieee);
date = ['ccs_' datestr(time_hycom(1),7) datestr(time_hycom(1),28)];

nPx = 1;
nPy = 1;

%% density of water
rho = 1024;
%% CP constant for TEMP tranport
cp  = 4000;

%% Lats and longs where to compute HYCOM transports.
%% this is with respect to MITGCM model boundaries
%% dlat = dlon = 0.01 degrees for GOM GRID
%% (NORTH, SOUTH  Open Boundaries)
%% NORTH BOUND
tablat_N = [(yc(end))];
%% south BOUND
tablat_S = [(yc(1))];
%% EAST BOUND
tablon_E = [(xc(end))];
%% WEST BOUND
tablon_W = [(xc(1))];

%% Grid - Mask
msk  = rdmds([fpath 'MASK/hFacC']);
msku = rdmds([fpath 'MASK/hFacW']);
mskv = rdmds([fpath 'MASK/hFacS']);

msk(msk~=0) = 1;
msku(msku~=0) = 1;
mskv(mskv~=0) = 1;

dxv = rdmds([fpath 'MASK/DXG']);
dxc = rdmds([fpath 'MASK/DXG']);
dyu = rdmds([fpath 'MASK/DYG']);
dyc = rdmds([fpath 'MASK/DYG']);

%% this mask file is diff from pmask.bin
%% this is obtained by spininig the model(actual model mask, which is
%% created by the MITGCM model
%% (initially run the model to crash the run and get this file, this is writen in STDOUT.0001)

drf = squeeze(rdmds([fpath 'MASK/DRF']));
drc = squeeze(rdmds([fpath 'MASK/DRC']));

vs  = zeros(nxc,nzc,nt1);
vn  = zeros(nxc,nzc,nt1);
uw  = zeros(nyc,nzc,nt1);

hts = zeros(nxc,nzc,nt1);
htn = zeros(nxc,nzc,nt1);
hte = zeros(nyc,nzc,nt1);

%% Position of the open boundaries (as defined in the nested model)
%% North.
obcsN = 'o';
nsn1 = 13;
nsn2 = 406;
if strcmp(obcsN,'o') == 1
    OB_Jnorth = [nyc*zeros(nsn1,1); nyc*ones(nsn2,1); nyc*zeros(nxc - nsn1 - nsn2,1)]';
end

%% South.
obcsS = 'o';
nss1 = 13;
nss2 = 406;
if strcmp(obcsS,'o') == 1
    OB_Jsouth = [zeros(nss1,1); ones(nss2,1); zeros(nxc - nss1 - nss2,1)]';
end

%% West.
obcsW = 'o';
nsw1 = 13;
nsw2 = 230;
if strcmp(obcsW,'o') == 1
    OB_Iwest = [zeros(nsw1,1); ones(nsw2,1); zeros(nyc - nsw1 - nsw2,1)]';
end
%% East.
obcsE = 'n';

%% HYCOM fluxes %%% RAW ANALYSIS
%% Meridional transport.
eval(['load ' fpath 'HYCOM_HF_MT_' date]);
vhycomtotm = eval(['HYCOM_HF_MT_' date]);
%% writen in different way for each slice (refer hycomtransp_NS.m)
hthycomn = vhycomtotm(:,1);
hthycoms = vhycomtotm(:,2);
vthycomn = vhycomtotm(:,3);
vthycoms = vhycomtotm(:,4);

%% Zonal transport.
eval(['load ' fpath 'HYCOM_HF_ZT_' date]);
uhycomtotz = eval(['HYCOM_HF_ZT_' date]);
%% writen in different way for each slice (refer hycomtransp_EW.m)
hthycomw = uhycomtotz(:,1);
uthycomw = uhycomtotz(:,2);

%% Balance transport: Take out 1/? of the total transport from each
%% of the three boundaries.
totransp =  vthycomn - vthycoms - uthycomw;

%% computing hycom alone tranport
totransp_hycom = vthycomn - vthycoms - uthycomw;
vthycoms1 = vthycoms;
vthycomn1 = vthycomn;
uthycomw1 = uthycomw;

hycom_s = mean(vthycoms1);
hycom_n = mean(vthycomn1);
hycom_w = mean(uthycomw1);

%% balancing transport
abstotransp = abs(vthycomn) + abs(vthycoms) + abs(uthycomw);
cn = abs(vthycomn./abstotransp);
cs = abs(vthycoms./abstotransp);
cw = abs(uthycomw./abstotransp);

vthycomn = vthycomn - cn.*totransp;
vthycoms = vthycoms + cs.*totransp;
uthycomw = uthycomw + cw.*totransp;

%% Open boundaries files
%% read interpolated fields at open boundaries

%% North
if strcmp(obcsN,'o') == 1
    obnv = rdda([fpath 'Rs_VobcsN_' date '_nPy' num2str(nPy) '.bin'],[nxc nPy nzc nt1],1,fmt,Ieee);
    obnt = rdda([fpath 'Rs_TobcsN_' date '_nPy' num2str(nPy) '.bin'],[nxc nPy nzc nt1],1,fmt,Ieee);
end

%% South
if strcmp(obcsS,'o') == 1
    obsv = rdda([fpath 'Rs_VobcsS_' date '_nPy' num2str(nPy) '.bin'],[nxc nPy nzc nt1],1,fmt,Ieee);
    obst = rdda([fpath 'Rs_TobcsS_' date '_nPy' num2str(nPy) '.bin'],[nxc nPy nzc nt1],1,fmt,Ieee);
end

%% West
if strcmp(obcsW,'o') == 1
    obwu = rdda([fpath 'Rs_UobcsW_' date '_nPx' num2str(nPx) '.bin'],[nyc nPx nzc nt1],1,fmt,Ieee);
    obwt = rdda([fpath 'Rs_TobcsW_' date '_nPx' num2str(nPx) '.bin'],[nyc nPx nzc nt1],1,fmt,Ieee);
end

%% Compute the model transport at the boundaries
%% HERE HYCOM IS INTERPOLATED TO MODEL USUNGI INTERP2 (BILINEAR
%% INTERPOLATION), OR INTERP2GAUSS>> AND THE TRANSPORT IS CALCULATED WITH
%% RESPECT TO MITGCM MODEL GRID

%% Southern transport
if strcmp(obcsS,'o') == 1
    surfs = 0;
    for ii = 2:nxc-1
        if OB_Jsouth(ii) ~= 0
            ij = OB_Jsouth(ii);
            for ik = 1:nzc
                surfs = surfs + dxv(ii,ij+1)*drf(ik)*mskv(ii,ij+1,ik);
                vs(ii,ik,:)  = obsv(ii,1,ik,:)*dxv(ii,ij+1)*drf(ik)*mskv(ii,ij+1,ik);
                hts(ii,ik,:) = rho*cp*obst(ii,1,ik,:).*obsv(ii,1,ik,:)*dxc(ii,ij)*drf(ik)*mskv(ii,ij+1,ik);
            end
        end
    end
    vstot  = squeeze(sum(sum(vs)));
    htstot = squeeze(sum(sum(hts)));
end

%% NORTHERN transport
if strcmp(obcsN,'o') == 1
    surfn = 0;
    for ii = 2:nxc-1
        if OB_Jnorth(ii) ~= 0
            ij = OB_Jnorth(ii);
            for ik = 1:nzc
                surfn = surfn + dxv(ii,ij)*drf(ik)*mskv(ii,ij,ik);
                vn(ii,ik,:)  = obnv(ii,1,ik,:)*dxv(ii,ij)*drf(ik)*mskv(ii,ij,ik);
                htn(ii,ik,:) = rho*cp*obnv(ii,1,ik,:).*obnt(ii,1,ik,:)*dxc(ii,ij)*drf(ik)*mskv(ii,ij,ik);
            end
        end
    end
    vntot  = squeeze(sum(sum(vn)));
    htntot = squeeze(sum(sum(htn)));
end

%% Western transport
if strcmp(obcsW,'o') == 1
    surfw = 0;
    for ij = 2:nyc-1
        if OB_Iwest(ij) ~= 0
            ii = OB_Iwest(ij);
            for ik = 1:nzc
                surfw = surfw + dyu(ii+1,ij)*drf(ik)*msku(ii+1,ij,ik);
                uw(ij,ik,:)  = obwu(ij,1,ik,:)*dyu(ii+1,ij)*drf(ik)*msku(ii+1,ij,ik);
                htw(ij,ik,:) = rho*cp*obwu(ij,1,ik,:).*obwt(ij,1,ik,:)*dyc(ii,ij)*drf(ik)*msku(ii+1,ij,ik);
            end
        end
    end
    uwtot  = squeeze(sum(sum(uw)));
    htwtot = squeeze(sum(sum(htw)));
end

%% ADDING THE DIFFERENCE BETWEEN THE HYCOM RAW AND INTERPOLATED MITGCM
%% TRANSPORTS TO MITGCM TRANSPORT, NORMALIZED BY THE AREA SLICE

%% --> North : vntotcor = vntot+cor = vthycomn;
%% Compute correction
if strcmp(obcsN,'o') == 1
    vcorn = (vthycomn - vntot)/surfn;
    obnvcor = obnv;
    for it = 1:nt1
        for ii = 2:nxc-1
            if OB_Jnorth(ii) ~= 0
                ij = OB_Jnorth(ii);
                for ik = 1:nzc
                    obnvcor(ii,:,ik,it) = (obnvcor(ii,:,ik,it) + vcorn(it))*mskv(ii,ij,ik);
                end
            end
        end
    end
    
    %% Compute new transport
    for ii = 2:nxc-1
        if OB_Jnorth(ii) ~= 0
            ij = OB_Jnorth(ii);
            for ik = 1:nzc
                vncor(ii,ik,:)  = obnvcor(ii,1,ik,:)*dxv(ii,ij)*drf(ik)*mskv(ii,ij,ik);
                htncor(ii,ik,:) = rho*cp*obnvcor(ii,1,ik,:).*obnt(ii,1,ik,:)*dxc(ii,ij)*drf(ik)*mskv(ii,ij,ik);
            end
        end
    end
    vncortot  = squeeze(sum(sum(vncor)));
    htncortot = squeeze(sum(sum(htncor)));
end

%% --> South : vstotcor = vstot+cor = vthycoms;
%% Compute correction
if strcmp(obcsS,'o') == 1
    vcors = (vthycoms - vstot)/surfs;
    obsvcor = obsv;
    
    for it = 1:nt1
        for ii = 2:nxc-1
            if OB_Jsouth(ii) ~= 0
                ij = OB_Jsouth(ii);
                for ik = 1:nzc
                    obsvcor(ii,:,ik,it) = (obsvcor(ii,:,ik,it) + vcors(it))*mskv(ii,ij+1,ik);
                end
            end
        end
    end
    
    %% Compute new transport
    for ii = 2:nxc-1
        if OB_Jsouth(ii) ~= 0
            ij = OB_Jsouth(ii);
            for ik = 1:nzc
                vscor(ii,ik,:)  = obsvcor(ii,1,ik,:)*dxv(ii,ij+1)*drf(ik)*mskv(ii,ij+1,ik);
                htscor(ii,ik,:) = rho*cp*obst(ii,1,ik,:).*obsvcor(ii,1,ik,:)*dxc(ii,ij)*drf(ik)*mskv(ii,ij+1,ik);
            end
        end
    end
    vscortot  = squeeze(sum(sum(vscor)));
    htscortot = squeeze(sum(sum(htscor)));
end

%% --> West : uwtotcor = uwtot+cor = uthycomw;
%% Compute correction
if strcmp(obcsW,'o') == 1
    ucorw = (uthycomw - uwtot)/surfw;
    obwucor = obwu;
    for it = 1:nt1
        for ij = 2:nyc-1
            if OB_Iwest(ij) ~= 0
                ii = OB_Iwest(ij);
                for ik = 1:nzc
                    obwucor(ij,:,ik,it) = (obwucor(ij,:,ik,it) + ucorw(it))*msku(ii+1,ij,ik);
                end
            end
        end
    end
    
    %% Compute new transport
    for ij = 2:nyc-1
        if OB_Iwest(ij) ~= 0
            ii = OB_Iwest(ij);
            for ik = 1:nzc
                uwcor(ij,ik,:)  = obwucor(ij,1,ik,:)*dyu(ii+1,ij)*drf(ik)*msku(ii+1,ij,ik);
                htwcor(ij,ik,:) = rho*cp*obwt(ij,1,ik,:).*obwucor(ij,1,ik,:)*dyc(ii,ij)*drf(ik)*msku(ii+1,ij,ik);
            end
        end
    end
    uwcortot  = squeeze(sum(sum(uwcor)));
    htwcortot = squeeze(sum(sum(htwcor)));
end

%% PLOTTING MITGCM MODEL INTERPOLATED TRANSPORT
figure
subplot(1,2,1)
dummy = squeeze(obnvcor(:,1,:,1))';
dummy(dummy == 0) = NaN;
contourf(xc,zc,dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridiaonal Velocity at ' num2str(round(tablat_N)) '^oN'],'FontSize',14);
xlabel('Longitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

subplot(1,2,2)
dummy = squeeze(obnt(:,1,:,1))';
dummy(dummy == 0) = NaN;
contourf(xc,zc,dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridiaonal Temp at ' num2str(round(tablat_N)) '^oN'],'FontSize',14);
eval(['print -depsc ' fpath 'MITGCM_INTER_VEL_TEMP_MER_' date '_' num2str(round(tablat_N)) 'N']);
xlabel('Longitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

figure
subplot(1,2,1)
dummy = squeeze(obsvcor(:,1,:,1))';
dummy(dummy == 0) = NaN;
contourf(xc,zc,dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridiaonal Velocity at ' num2str(round(tablat_S)) '^oN'],'FontSize',14);
xlabel('Longitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

subplot(1,2,2)
dummy = squeeze(obst(:,1,:,1))';
dummy(dummy == 0) = NaN;
contourf(xc,zc,dummy)
colorbar
set(gca,'FontSize',14);
title(['Meridiaonal Temp at ' num2str(round(tablat_S)) '^oN'],'FontSize',14);
eval(['print -depsc ' fpath 'MITGCM_INTER_VEL_TEMP_MER_' date '_' num2str(round(tablat_S)) 'N']);
xlabel('Longitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

figure
subplot(1,2,1)
dummy  = squeeze(obwucor(:,1,:,1))';
dummy(dummy == 0) = NaN;
contourf(yc,zc,dummy)
colorbar
set(gca,'FontSize',14);
title(['Zonal Velocity at ' num2str(round(tablon_W)) '^oE'],'FontSize',14);
xlabel('Latitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

subplot(1,2,2)
dummy = squeeze(obwt(:,1,:,1))';
dummy(dummy == 0) = NaN;
contourf(yc,zc,dummy)
colorbar
set(gca,'FontSize',14);
title(['Zonal Temp at ' num2str(round(tablon_W)) '^oE'],'FontSize',14);
eval(['print -depsc ' fpath 'MITGCM_INTER_VEL_TEMP_ZON_' date '_' num2str(round(tablon_W)) 'E']);
xlabel('Latitude (deg) ','FontSize',14)
ylabel('Depth (m)','FontSize',14)

%% All
%% BALACING THE CORRECTED TRANSPORT SO AS TO MAKE THE SUMMATION IS ZERO
%% VOLUME TRANSPORT IS CONSERVED
totransp_uncor    =  vntot - vstot - uwtot;
totransp_cor = vncortot - vscortot - uwcortot;

%% Plots
mehts = mean([htscortot htstot hthycoms]);
mehtn = mean([htncortot htntot hthycomn]);
mehtw = mean([htwcortot htwtot hthycomw]);

mevs = mean([vscortot vstot vthycoms]);
mevn = mean([vncortot vntot vthycomn]);
meuw = mean([uwcortot uwtot uthycomw]);

mevs2 = num2str(round(mevs'*100*1e-6)/100);
mevn2 = num2str(round(mevn'*100*1e-6)/100);
meuw2 = num2str(round(meuw'*100*1e-6)/100);

mehts2 = num2str(round(mehts'*100*1e-15)/100);
mehtn2 = num2str(round(mehtn'*100*1e-15)/100);
mehtw2 = num2str(round(mehtw'*100*1e-15)/100);

%% South
figure
subplot(2,1,1); hold on;
set(gca,'FontSize',14)
h1 = plot([1:nt1],vscortot*1e-6,'b',[1:nt1],vstot*1e-6,'r',[1:nt1],vthycoms*1e-6,'g--','LineWidth',2); grid on;
title(['Meridional Transport (Sv) - ' num2str(round(tablat_S)) '^oN - means = ' mevs2(1,:) ' ' mevs2(2,:) ' ' mevs2(3,:)])
legend(h1,'Corrected','Non-Corrected','HYCOM');
subplot(2,1,2); hold on; grid on;
set(gca,'FontSize',14)
plot([1:nt1],htscortot*1e-15,'b',[1:nt1],htstot*1e-15,'r',[1:nt1],hthycoms*1e-15,'g','LineWidth',2); grid on;
title(['Meridional Heat Transport (10^{15} W) - ' num2str(round(tablat_S)) '^oN - means = ' mehts2(1,:) ' ' mehts2(2,:) ' ' mehts2(3,:)])
legend('Corrected','Non-Corrected','HYCOM');
eval(['print -depsc ' fpath 'COR_VEL_TEMP_MER_' date '_' num2str(round(tablat_S)) 'N']);

%% North
figure
subplot(2,1,1); hold on;
set(gca,'FontSize',14)
h1 = plot([1:nt1],vncortot*1e-6,'b',[1:nt1],vntot*1e-6,'r',[1:nt1],vthycomn*1e-6,'g--','LineWidth',2); grid on;
title(['Meridional Transport (Sv) - ' num2str(round(tablat_N)) '^oN - means = ' mevn2(1,:) ' ' mevn2(2,:) ' ' mevn2(3,:)])
legend(h1,'Corrected','Non-Corrected','HYCOM');
subplot(2,1,2); hold on; grid on;
set(gca,'FontSize',14)
plot([1:nt1],htncortot*1e-15,'b',[1:nt1],htntot*1e-15,'r',[1:nt1],hthycomn*1e-15,'g','LineWidth',2); grid on;
title(['Meridional Heat Transport (10^{15} W) - ' num2str(round(tablat_N)) '^oN - means = ' mehtn2(1,:) ' ' mehtn2(2,:) ' ' mehtn2(3,:)])
legend('Corrected','Non-Corrected','HYCOM');
eval(['print -depsc ' fpath 'COR_VEL_TEMP_MER_' date '_' num2str(round(tablat_N)) 'N']);

%% West
figure
subplot(2,1,1); hold on;
set(gca,'FontSize',14)
plot([1:nt1],uwcortot*1e-6,'b',[1:nt1],uwtot*1e-6,'r',[1:nt1],uthycomw*1e-6,'g--','LineWidth',2);grid on;
title(['Zonal Transport (Sv) - ' num2str(round(tablon_W)) '^oE - means = ' meuw2(1,:) ' ' meuw2(2,:) ' ' meuw2(3,:)])
legend('Corrected','Non-Corrected','HYCOM');
subplot(2,1,2); hold on; grid on;
set(gca,'FontSize',14)
plot([1:nt1],htwcortot*1e-15,'b',[1:nt1],htwtot*1e-15,'r',[1:nt1],hthycomw*1e-15,'g--','LineWidth',2);grid on;
title(['Zonal Heat Transport (10^{15} W) - ' num2str(round(tablon_W)) '^oE - means = ' mehtw2(1,:) ' ' mehtw2(2,:) ' ' mehtw2(3,:)])
legend('Corrected','Non-Corrected','HYCOM');
eval(['print -depsc ' fpath 'COR_VEL_TEMP_ZON_' date '_' num2str(round(tablon_W)) 'E']);

figure
hold on;
set(gca,'FontSize',14);
plot([1:nt1],vscortot*1e-6,'g--',[1:nt1],vncortot*1e-6,'r--',[1:nt1],uwcortot*1e-6,'b--',[1:nt1],totransp_cor*1e-6,'k','LineWidth',2)
title(['Transport (Sv) means: S: ' mevs2(1,:) ',N: ' mevn2(1,:) ',W: ' meuw2(1,:)])
legend('South','North','West','Total'); grid on;
eval(['print -depsc ' fpath 'COR_VEL_TOTAL_' date '_' num2str(round(tablat_S)) '_' num2str(round(tablat_N)) '-' num2str(round(tablon_E)) '-' num2str(round(tablon_W))]);

figure
subplot(2,1,1);
hold on; grid on;
set(gca,'FontSize',14); plot([1:nt1],vstot*1e-6,'g--',[1:nt1],vntot*1e-6,'r--',[1:nt1],uwtot*1e-6,'b--',[1:nt1],totransp_uncor*1e-6,'k','LineWidth',2)
title(['SV Before Correction, Means: S: ' num2str(1e-6*mevs(2)) ',N: ' num2str(1e-6*mevn(2)) ', W: ' num2str(1e-6*meuw(2))])
legend('South','North','West','Total');
subplot(2,1,2);
hold on; grid on;
set(gca,'FontSize',14); plot([1:nt1],vscortot*1e-6,'g--',[1:nt1],vncortot*1e-6,'r--',[1:nt1],uwcortot*1e-6,'b--',[1:nt1],totransp_cor*1e-6,'k','LineWidth',2)
title(['Sv After Correction, Means: S: ' num2str(1e-6*mevs(1)) ',N: ' num2str(1e-6*mevn(1)) ',W: ' num2str(1e-6*meuw(1))])
legend('South','North','West','Total');
eval(['print -depsc ' fpath 'obcall_transp_' date ]);

figure
set(gca,'FontSize',14); plot([1:nt1],vthycoms1*1e-6,'g--',[1:nt1],vthycomn1*1e-6,'r--',[1:nt1],uthycomw1*1e-6,'b--',[1:nt1],totransp_hycom*1e-6,'k','LineWidth',2)
grid on
title(['SV HYCOM, Means: S: ' num2str(1e-6*hycom_s) ',N: ' num2str(1e-6*hycom_n) ', W: ' num2str(1e-6*hycom_w)])
legend('South','North','West','Total');
eval(['print -depsc ' fpath 'HYCOM_transp_' date ]);

%% Write output
%% North
if strcmp(obcsN,'o') == 1
    wrslice([fpath 'Rs_VobcsN_' date '_nPy' num2str(nPy) '_c1.bin' ],obnvcor,1,fmt,Ieee);
end

%% East
if strcmp(obcsE,'o') == 1
    wrslice([fpath 'Rs_UobcsE_' date '_nPx' num2str(nPy) '_c1.bin' ],obeucor,1,fmt,Ieee);
end

%% West
if strcmp(obcsW,'o') == 1
    wrslice([fpath 'Rs_UobcsW_' date '_nPx' num2str(nPy) '_c1.bin' ],obwucor,1,fmt,Ieee);
end

%% South
if strcmp(obcsS,'o') == 1
    wrslice([fpath 'Rs_VobcsS_' date '_nPy' num2str(nPy) '_c1.bin' ],obsvcor,1,fmt,Ieee);
end
