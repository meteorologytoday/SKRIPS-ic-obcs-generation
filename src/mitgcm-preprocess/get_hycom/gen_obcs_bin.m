function gen_obcs_bin(obcs, input_dir, output_dir, start_date, end_date, fmt, Ieee)

    time = [datenum(start_date, 'yyyy-mm-dd'):1:datenum(end_date, 'yyyy-mm-dd')];
    
    % FMT.mat contains two strings: fmt = 'real*4' and Ieee = 'b'
    % eval(['load ' fpath 'FMT.mat']);

    % Load the first timestep to get spatial coordinates
        

    start_time_label = datestr(time(1), 'yyyy-mm-dd_hh');
    example_file = sprintf('%s/hycom_%s_%s.mat', input_dir, obcs, start_time_label);
    D = load( sprintf('%s/hycom_%s_%s.mat', input_dir, obcs, start_time_label) );
    D = D.D;

    lon_hycom = D.Longitude;
    lat_hycom = D.Latitude;
    z_hycom = D.Depth;
    clear D;

    nx = length(lon_hycom);
    ny = length(lat_hycom);
    nz = length(z_hycom);
    nt = length(time);

    % Output lat lon depth as bin files
    fprintf('Output spatial coordinates\n');
    wrslice(sprintf('%s/lon_%s.bin', output_dir, obcs),   lon_hycom, 1, fmt, Ieee);
    wrslice(sprintf('%s/lat_%s.bin', output_dir, obcs),   lat_hycom, 1, fmt, Ieee);
    wrslice(sprintf('%s/depth_%s.bin', output_dir, obcs), z_hycom, 1, fmt, Ieee);

    % currently I only do hour = 00
    hours = {'00'};
    
            
    params = {'surf_el', 'water_u', 'water_v', 'water_temp', 'salinity'};
    % Loop through the rest of the dates
    totaln = 0
    for i = 1:nt
        for ih = 1:length(hours) 

            %fprintf('%s', hours{ih}); 
            time_label = sprintf('%s_%s', datestr(time(i), 'yyyy-mm-dd'), hours{ih});
            
            mat_filename = sprintf('%s/hycom_%s_%s.mat', input_dir, obcs, time_label);
            fprintf('Loading mat file: %s\n', mat_filename);
            D = load(mat_filename);
            D = D.D;
 
            time_hycom = D.Date;

            % implicitly append the new datetime string
            wrslice(sprintf('%s/time_%s.bin', output_dir, obcs), time_hycom, totaln+1, fmt, Ieee);
             
            for np = 1:length(params)
                param = params{np};
                eval(['tmp = ' 'D.' param ';']);
                tmp(isnan(tmp)) = 0;
                
                output_file = sprintf('%s/%s_%s.bin', output_dir, param, obcs);
                fprintf('Output file: %s\n', output_file);
                wrslice(output_file, tmp, totaln+1, fmt, Ieee);
                clear tmp param;
            end
            totaln = totaln + 1;
            clear D;

        end
    end
end
