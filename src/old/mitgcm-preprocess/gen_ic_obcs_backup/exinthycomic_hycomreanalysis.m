%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract HYCOM fields for a given domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tm,xc,yc,zc] = exinthycomic_hycomreanalysis(fpath,var,data,nrec,lon_hycom,lat_hycom,z_hycom,xc,yc,xf,yf,zc)

warning('off')

eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);
xc = XC(:,1); yc = YC(1,:); zc = RC;
xf = XG(:,1); yf = YG(1,:);
nxc = length(xc); nyc = length(yc);nzc = length(zc);

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

TL = squeeze(data);
xl = lon_hycom;
yl = lat_hycom;
zl = z_hycom;
[nxl,nyl,nzl] = size(TL);

%% Fill-in missing-data as much as possible horizontally then vertically
%% Layer 1.
TL(:,:,1) = xyexpand(TL(:,:,1),100);
%% Layer 2.
TL(:,:,2) = xyexpand( TL(:,:,2),100);
%% Layers 3 to nz.
for k = 3:nzl
    TTL = xyexpand( TL(:,:,k),100);
    TTLe = 2*TL(:,:,k-1) - TL(:,:,k-2);
    TTL(isnan(TTL)) = TTLe(isnan(TTL));
    TL(:,:,k) = TTL;
end
clear TTL TTLe;

%% Interpolate data vertically.
for j = 1:nyl
    TLz(:,j,:) = interp1(zl,sq(TL(:,j,:))',zc)';
end

%% Set nans to zeros.
TLz(isnan(TLz)) = 0;

%% Interpolate horizontaly.
Tm = zeros(nxc,nyc,nzc);
report('XY-Interpolating k = xxx');
for k = 1:nzc,
    report('\b\b\b%3i',k)
    if strcmp(var,'T') == 1 | strcmp(var,'S') == 1
        Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yc,xc,'linear');
    elseif strcmp(var,'U') == 1
        Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yc,xf,'linear');
    elseif strcmp(var,'V') == 1
        Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yf,xc,'linear');
    end
end

%% Mask.
Tm(isnan(Tm)) = 0;
Tm = Tm.*msk;

%% fixing zeros in U/V to small value
if strcmp(var,'U') || strcmp(var,'V')
    val = 1*10^-5;
    tmp = Tm;
    tmp(msk == 0) = nan;
    tmp(tmp == 0) = val;
    tmp(isnan(tmp)) = 0;
    Tm = tmp;
    clear tmp
end

