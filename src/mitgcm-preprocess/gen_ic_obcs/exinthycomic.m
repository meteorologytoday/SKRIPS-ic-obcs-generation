%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract HYCOM fields for a given domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tm,xc,yc,zc] = exinthycomic(varname, data, lon_hycom, lat_hycom, z_hycom, xc, yc, xf, yf, zc, mask_dir)

    %warning('off')

    %eval(['load ' fpath 'FMT.mat']);
    %eval(['load ' fpath 'grid.mat']);
    % xc = XC(:,1); yc = YC(1,:); zc = RC;
    % xf = XG(:,1); yf = YG(1,:);
    % nxc = length(xc); nyc = length(yc);nzc = length(zc);

    if strcmp(varname,'T') == 1
        msk = rdmds([mask_dir '/hFacC']);
    elseif strcmp(varname,'S') == 1
        msk = rdmds([mask_dir '/hFacC']);
    elseif strcmp(varname,'U') == 1
        msk = rdmds([mask_dir '/hFacW']);
    elseif strcmp(varname,'V') == 1
        msk = rdmds([mask_dir '/hFacS']);
    end

    msk(msk~=0) = 1;

    nxc = length(xc);
    nyc = length(yc);
    nzc = length(zc);

    TL = squeeze(data);
    xl = lon_hycom;
    yl = lat_hycom;
    zl = z_hycom;
    [nxl,nyl,nzl] = size(TL);

    %% Fill-in missing-data as much as possible horizontally then vertically
    %% Layer 1.
    TL(:,:,1) = xyexpand(TL(:,:,1),50);
    %% Layer 2.
    TL(:,:,2) = xyexpand( TL(:,:,2),50);
    %% Layers 3 to nz.
    for k = 3:nzl
        %report('filling missing data: %3i',k)
        fprintf('Fill the missing data of the %d-th layer.\n', k);
        TTL = xyexpand( TL(:,:,k),50);
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
    for k = 1:nzc
        fprintf('Interpolating horizontally the missing data of the %d-th layer.\n', k);
        if strcmp(varname,'T') == 1 || strcmp(varname,'S') == 1
            Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yc,xc,'linear');
        elseif strcmp(varname,'U') == 1
            Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yc,xf,'linear');
        elseif strcmp(varname,'V') == 1
            Tm(:,:,k) = interp2(yl',xl', squeeze(TLz(:,:,k)),yf,xc,'linear');
        end
    end

    %% Mask.
    Tm(isnan(Tm)) = 0;
    Tm = Tm.*msk;

end

