function getHycomData(start_date, end_date, region, opath, OpenDAP_URL, f0)

        
    zl = [0:1:39];
    if ( strcmp(region, 'all') == 1 ) 
        xl = [2244:1:3056];
        yl = [1720:1:2356];

    elseif (strcmp(region, 'north') == 1)
        xl = [2244:1:3056];
        yl = [2344:1:2356];

    elseif (strcmp(region, 'east') == 1)
        xl = [3044:1:3056];
        yl = [1720:1:2356];

    elseif (strcmp(region, 'south') == 1)
        xl = [2244:1:3056];
        yl = [1720:1:1732];

    elseif (strcmp(region, 'west') == 1)
        xl = [2244:1:2256];
        yl = [1720:1:2356];

    else
        error(['Unrecognized region: ', region]);
    end

    display(['Region: ', region])
    display(['Start date: ', start_date])
    display(['End date:   ', end_date])

    %% TO EXTRACT HYCOM DATA FROM HYCOM DATA SERVER USING OPENDAP TOOLS %%
    %% THIS SCRIPT IS MODIFIED FOR HYCOM_GLOBAL VERSION
    %% THIS IS FOR CCS

    % xl from 0 to 360, 0.08 degree interval, x = 2250 for 180 degree
    % yl from -80 to -40 and 40 to 80, 0.08 degree interval
    % yl from -40 to 40, 0.04 degree interval
    % y = 1500 for 0 degree

    hours = {'00','03','06','09','12','15','18','21'};
    date_inc = 1;

    format = 'yyyymmdd';
    time = [datenum(start_date,format):date_inc:datenum(end_date(1:8),format)];
    nt = length(time);

    for i = 1:nt
      for ih = 1:1



        f1 = f0+i-1
        tl = datestr(time(i),format);

        frame = f1;
        str_date = tl;
        hour = hours{ih};
        z = zl;
        y = yl;
        x = xl;

        D = struct('Date', [], 'Depth', [],'Latitude', [], 'Longitude', [],'ssh',[],'u',[],'v',[],'temperature',[],'salinity',[]);
        
        t = num2str(frame);
        p = blanks(4 - length(t));p(:) = '0';
        frame = [num2str(p) t];
        disp (frame)
        
        prec='double';
        
        xmin = x(1);
        xmax = x(end);
        ymin = y(1);
        ymax = y(end);
        zmin = z(1);
        zmax = z(end);
        
        deltax = xmax - xmin + 1 ;
        deltay = ymax - ymin + 1 ;
        deltaz = zmax - zmin + 1 ;
        %% will need to update for an array for later
        deltat = 1;
        
        Z = zmax+1;

        
        ncid = netcdf.open(OpenDAP_URL,'NOWRITE')
        [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid)
        
        for varid =1:numvars
            varlist(varid) = {netcdf.inqVar(ncid,varid-1)};
        end
        
        
        %%  First get the hycom date and time
        param = 'time';
        ind=find(ismember(varlist,param))
        
        Temp_Date = netcdf.getVar(ncid,ind-1)
        str_date_convert = hycom_time(str_date,hour)
        stime = find(Temp_Date==str_date_convert)
        if isempty(stime)
          continue;
        end 

        hycom_date = Temp_Date(stime)
        hycom_date_save = hycom_revert_time(hycom_date)
        D.Date = hycom_date_save;
        
        %% to make array indicies correct for
        %% multidimensional arrays that start
        %% at 0 for opendap
        stime = stime-1;
        
        param = 'depth';
        ind = find(ismember(varlist,param));
        
        Temp_Depth = netcdf.getVar(ncid,ind-1,prec);
        hycom_depth = -1*Temp_Depth;
        %% HYCOM GOES ONLY UPTO 5500 m, so adding one more layer with same field as
        hycom_depth(end+1) = -6500;
        D.Depth = hycom_depth;
        
        param = 'lat';
        ind=find(ismember(varlist,param));
        
        hycom_lat = netcdf.getVar(ncid,ind-1,ymin,deltay,prec);
        hycom_lat
        D.Latitude = hycom_lat;
        D.Latitude = permute(D.Latitude,[2,1]);
        
        param = 'lon';
        ind=find(ismember(varlist,param));
        
        Temp_Lon = netcdf.getVar(ncid,ind-1,xmin,deltax,prec);
        hycom_lon = Temp_Lon;
        hycom_lon
        if hycom_lon < 0; hycom_lon = hycom_lon + 360; end
        if hycom_lon > 360; hycom_lon = hycom_lon - 360; end
        D.Longitude = hycom_lon;
        D.Longitude = permute(D.Longitude,[2,1]);
        
        params = {'surf_el', 'water_u', 'water_v', 'water_temp', 'salinity'};
        for np = 1:length(params)
            param = params{np};
            ind=find(ismember(varlist,param));
            if (strcmp(param,'surf_el') == 1)
                tmp = netcdf.getVar(ncid,ind-1,[xmin,ymin,stime],[deltax,deltay,deltat],prec);
            else
                tmp = netcdf.getVar(ncid,ind-1,[xmin,ymin,zmin,stime],[deltax,deltay,deltaz,deltat],prec);
                %% HYCOM GOES ONLY UPTO 5500 m, so adding one more layer with same field as
                tmp(:,:,Z+1) = tmp(:,:,Z);
            end

            tmpScaleName = netcdf.inqAttName(ncid,ind-1,6)
            tmpScale = netcdf.getAtt(ncid,ind-1,tmpScaleName)
            tmpOffsetName = netcdf.inqAttName(ncid,ind-1,7)
            tmpOffset = netcdf.getAtt(ncid,ind-1,tmpOffsetName)
            tmpScale = double(tmpScale);
            tmpOffset = double(tmpOffset);
            tmp = tmp*tmpScale+tmpOffset;

            miss_value = netcdf.getAtt(ncid,ind-1,'missing_value');
            miss_value = double(miss_value)*tmpScale+tmpOffset;
            miss_e = abs(1*tmpScale);
            field1 = tmp < miss_value + miss_e;
            field2 = tmp > miss_value - miss_e;
            bad_field = logical(floor((field1+field2+0.1)/2));
            bad = 0;
            tmp(bad_field) = bad;
            
            eval(['D.' param ' = tmp;']);
            clear tmp param
        end
        
        %% saving orginal HYCOM  fields into mat file
        varList = 'D';
        
        if (strcmp(region, 'all') == 1)
            savename = [opath '/' frame '_' num2str(hycom_date_save) '.mat'];
        else
            savename = [opath '/' frame '_' hour '_' region '.mat'];
        end

        fprintf('Output: %s', savename)

        eval(['save ',savename,' ',varList]);
        netcdf.close(ncid);

      end
    end

    return
end
