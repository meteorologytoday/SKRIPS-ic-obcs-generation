function status = check3dmask(fpath,datafile,nrec,var)

    %% check3dmask(datafile)

    %% Load 3D data in "datafile", calculate a missing
    %% value mask and compare to pmask.bin.
    %% They should always be the same!!!

    %% Created  11/11/99 by adcroft@mit.edu
    %% Modified          by
    %% Maintained by adcroft@mit.edu, abiastoch@ucsd.edu

    %% Read model grid and masks.
    eval(['load ' fpath 'FMT.mat']);
    eval(['load ' fpath 'grid.mat']);
    % xc = XC(:,1); yc = YC(1,:); zc = RC;
    % nxc = length(xc); nyc = length(yc);nzc = length(zc);

    if strcmp(var,'T') == 1
        Pmsk = rdmds([fpath 'MASK/hFacC']);
    elseif strcmp(var,'S') == 1
        Pmsk = rdmds([fpath 'MASK/hFacC']);
    elseif strcmp(var,'U') == 1
        Pmsk = rdmds([fpath 'MASK/hFacW']);
    elseif strcmp(var,'V') == 1
        Pmsk = rdmds([fpath 'MASK/hFacS']);
    end
    Pmsk(Pmsk~=0) = 1;
    sum(Pmsk(:))

    %% Read data.
    T = rdslice([fpath datafile],[nxc nyc nzc],nrec,fmt,Ieee);
    Tmsk = ones([nxc nyc nzc]);
    Tmsk(T==0) = 0;
    sum(Tmsk(:))

    %% Finf mismatchs.
    mismatch = find(Pmsk~=Tmsk);

    %% Stop if any mismatch.
    if isempty(mismatch)
        fprintf('check3dmask(%i th-record of "%s") passed  OK\n',nrec,datafile);
        status = 0;
        %report('check3dmask(%i th-record of "%s") passed  OK\n',nrec,datafile)
    else

        save('Tmsk.mat','Tmsk');
        save('Pmsk.mat','Pmsk');
        fprintf('I found %i points where %i th-record of %s mis-matched %s\n',...
            prod(size(mismatch)),nrec,[fpath datafile],[fpath 'pmask.bin']);

        %report('I found %i points where %i th-record of %s mis-matched %s\n',...
        %    prod(size(mismatch)),nrec,[fpath datafile],[fpath 'pmask.bin']);

        fprintf('check3dmask FAILED!!\n');

        status = 1;
    end

end

