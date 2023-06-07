function [status, Tmsk, Pmsk] = check3dmask(varname, datafile, mask_dir, nrec, gd, fmt, Ieee)

    %% check3dmask(datafile)

    %% Load 3D data in "datafile", calculate a missing
    %% value mask and compare to pmask.bin.
    %% They should always be the same!!!

    %% Created  11/11/99 by adcroft@mit.edu
    %% Modified          by
    %% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

    %% gd : the datastructure loaded from grid.mat

    %% Read model grid and masks.

    if strcmp(varname,'T') == 1
        Pmsk = rdmds([mask_dir '/hFacC']);
    elseif strcmp(varname,'S') == 1
        Pmsk = rdmds([mask_dir '/hFacC']);
    elseif strcmp(varname,'U') == 1
        Pmsk = rdmds([mask_dir '/hFacW']);
    elseif strcmp(varname,'V') == 1
        Pmsk = rdmds([mask_dir '/hFacS']);
    end
    Pmsk(Pmsk~=0) = 1;
    sum_Pmsk = sum(Pmsk(:));

    nxc = gd.nxc;
    nyc = gd.nyc;
    nzc = gd.nzc;

    %% Read data.
    T = rdslice(datafile,[nxc nyc nzc],nrec,fmt,Ieee);
    Tmsk = ones([nxc nyc nzc]);
    Tmsk(T==0) = 0;
    
    sum_Tmsk = sum(Tmsk(:));

    % Pmsk is the mask from mitgcm
    % Tmsk is the inferred mask from `datafile`


    %% Find mismatchs.
    mismatch = find(Pmsk~=Tmsk);
    
    fprintf('Sum of mitgcm mask : %d\n', sum_Pmsk);
    fprintf('Sum of datafile inferred mask : %d\n', sum_Tmsk);


    %% Stop if any mismatch.
    if isempty(mismatch)
        fprintf('check3dmask(%i th-record of "%s") passed  OK\n',nrec,datafile);
        status = 0;
        %report('check3dmask(%i th-record of "%s") passed  OK\n',nrec,datafile)
    else

        fprintf('I found %i points where %i th-record of %s mis-matched %s\n',...
            prod(size(mismatch)),nrec, datafile, [mask_dir '/pmask.bin']);

        %report('I found %i points where %i th-record of %s mis-matched %s\n',...
        %    prod(size(mismatch)),nrec,[fpath datafile],[fpath 'pmask.bin']);

        fprintf('check3dmask FAILED!!\n');

        status = 1;
    end

end

