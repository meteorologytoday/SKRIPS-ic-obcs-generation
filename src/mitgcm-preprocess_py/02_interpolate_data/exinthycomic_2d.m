%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract HYCOM fields for a given domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tm,xc,yc] = exinthycomic_2d(fpath,var,data,nrec,lon_hycom,lat_hycom,xc,yc,xf,yf,fmt,Ieee);
warning('off')

eval(['load ' fpath 'FMT.mat']);
eval(['load ' fpath 'grid.mat']);
xc = XC(:,1); yc = YC(1,:);
xf = XG(:,1); yf = YG(1,:);
nxc = length(xc); nyc = length(yc);

msk  = rdmds([fpath 'MASK/hFacC']);
msk(msk~=0) = 1;
msk = squeeze(msk(:,:,1));

TL = squeeze(data);
xl = lon_hycom;
yl = lat_hycom;
[nxl,nyl] = size(TL);

%% interp2guass needs a dummy mask
maskinter2guass = ones(nxc,nyc);

%% Interpolate horizontaly.
T = zeros(nxc,nyc);
report('XY-Interpolating k = xxx');
if (strcmp(var,'SSH') == 1)
    TTL = xyexpand(TL(:,:),100);
    TTL(isnan(TTL) == 1) = 0;
    T = interp2(yl',xl', squeeze(TTL),yc,xc,'linear');
    T(isnan(T)) = 0;
end
report('\n')
clear TL TTL;

%% Mask and save outputs.
kdep = zeros(nxc,nyc);
report('Masking k = xxx');
Ti = T.*msk(:,:);
Tm = reshape(Ti,[nxc*nyc 1]);
Tm(isnan(Tm)) = 0;

