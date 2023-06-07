%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract HYCOM fields for a given domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tm,xc,yc,zc] = hycom2modelgrid(varname, data, lon_hycom, lat_hycom, z_hycom, xc, yc, xf, yf, zc, mask_dir, idx_lon, idx_lat)
%% data               : binary HYCOM data file
%% xc, yc, xf, yf, zc : target grid
 

    nxc = length(xc);
    nyc = length(yc);
    nzc = length(zc);
  
    fprintf('(nxc, nyc, nzc) = (%d, %d, %d)\n', nxc, nyc, nzc);
 
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

    % Slice the msk according to idx_lon and idx_lat
    if idx_lon == 0
        idx_lon = 1:nxc;
    end

    if idx_lat == 0
        idx_lat = 1:nyc;
    end

    msk = msk(idx_lon, idx_lat,:);
  
    % Compute the filled HYCOM data first (TL)
    TL = squeeze(data);
    xl = lon_hycom;
    yl = lat_hycom;
    zl = z_hycom;
    [nxl,nyl,nzl] = size(TL);

    % Fill-in HYCOM missing-data as much as possible horizontally
    % BEFORE interpolating to target grid
    iterations = 50;
    for k = 1:nzl
        fprintf('Fill the missing data of the %d-th layer.\n', k);

        if k <= 2    
            TL(:,:,k) = xyexpand(TL(:,:,1), iterations);

        else
            
            TTL = xyexpand( TL(:,:,k), iterations);
            
            % The east-most grid point will be extrapolated if missing
            TTLe = 2*TL(:,:,k-1) - TL(:,:,k-2);
            TTL(isnan(TTL)) = TTLe(isnan(TTL));
            
            TL(:,:,k) = TTL;
        end

    end

    % Now, interpolate data vertically to the target grid
    for j = 1:nyl
        TLz(:,j,:) = interp1(zl, sq(TL(:,j,:))', zc)';
    end

    %% Set nans to zeros.
    TLz(isnan(TLz)) = 0;

    %% Interpolate horizontaly.
    Tm = zeros(nxc, nyc, nzc);
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

