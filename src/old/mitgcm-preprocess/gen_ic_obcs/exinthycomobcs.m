%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract HYCOM fields for a given domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tm,xc1,yc1,zc1] = exinthycomobcs(fpath,var,data,lon_hycom,lat_hycom,z_hycom,xc1,yc1,xf1,yf1,IX,IY,zc1)

warning('off')

eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);
% xc = XC(:,1); yc = YC(1,:); zc = RC;
% nxc = length(xc); nyc = length(yc);nzc = length(zc);

if strcmp(var,'T') == 1
    msk = rdmds([fpath 'MASK/hFacC']);
elseif strcmp(var,'S') == 1
    msk = rdmds([fpath 'MASK/hFacC']);
elseif strcmp(var,'U') == 1
    msk = rdmds([fpath 'MASK/hFacW']);
elseif strcmp(var,'V') == 1
    msk = rdmds([fpath 'MASK/hFacS']);
end

msk(msk~=0) = 1;

%% getting mask only for slice
msk = msk(IX,IY,:);

%% getting slice dimensions
LX = length(xc1);
LY = length(yc1);
LZ = length(zc1);

TL = squeeze(data);
xl = lon_hycom;
yl = lat_hycom;
zl = z_hycom;
[nxl,nyl,nzl] = size(TL);

%% Fill-in missing-data as much as possible horizontally then vertically
%% Layer 1.
TL(:,:,1) = xyexpand(TL(:,:,1),25);
%% Layer 2.
TL(:,:,2) = xyexpand( TL(:,:,2),25);
%% Layers 3 to nz.
for k = 3:nzl,
    TTL = xyexpand( TL(:,:,k),25);
    TTLe = 2*TL(:,:,k-1) - TL(:,:,k-2);
    TTL(isnan(TTL)) = TTLe(isnan(TTL));
    TL(:,:,k) = TTL;
end
clear TTL TTLe;

%% Interpolate data vertically.
for j = 1:nyl
    TLz(:,j,:) = interp1(zl,sq(TL(:,j,:))',zc1)';
end

%% Set nans to zeros.
TLz(isnan(TLz)) = 0;

%% Interpolate horizontaly.
Tm = zeros(LX,LY,LZ);
report('XY-Interpolating k = xxx');
for k = 1:LZ,
    report('\b\b\b%3i',k)
    if strcmp(var,'T') == 1 | strcmp(var,'S') == 1
        Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yc1,xc1,'linear');
    elseif strcmp(var,'U') == 1
        Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yc1,xf1,'linear');
    elseif strcmp(var,'V') == 1
        Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yf1,xc1,'linear');
    end
end

%% Mask.
Tm(isnan(Tm)) = 0;
Tm = Tm.*msk;
