function status = gen_ic_hycom(grid_dir, mask_dir, hycom_data_dir, output_dir, target_date)
%% `grid_dir`       : should contain the file grid.mat
%% `mask_dir`       : should contain mask files generated by mitgcm
%% `hycom_data_dir` : should contain the initial condition mat file
%% `output_dir`     : the output directory
%% `target_date`    : should be a formatted string `yyyy-mm-dd_hh`

    status = 1;

    run('pathdef.m');
    FMT_json = read_json('../FMT.json');
    fmt = FMT_json.fmt;
    Ieee = FMT_json.Ieee;
    
    grid_file = sprintf('%s/grid.mat', grid_dir);
    fprintf('Loading grid file: %s\n', grid_file);
    gd = load(grid_file);

    % HYCOM outputs
    ic_file = sprintf('%s/hycom_%s.mat', hycom_data_dir, target_date);
    fprintf('Loading initial condition data file: %s\n', ic_file);

    if exist(ic_file) ~= 2
        fprintf('Error reading file: %s\n', ic_file);
        fprintf('Not doing anything.\n');

        return
    end
    D = load(ic_file);
    D = D.D;

    hycom_lon = D.Longitude;
    hycom_lat = D.Latitude;
    hycom_z = D.Depth;
    hycom_time = D.Date;
    hycomu = D.water_u;
    hycomv = D.water_v;
    hycomt = D.water_temp;
    hycoms = D.salinity;

    % Extract HYCOM fields at this time
    varnames = {'T', 'S', 'U', 'V'};
    vars = {hycomt, hycoms, hycomu, hycomv};
    maskchk = [];    
    interpolated_vars = {};
    
    for i = 1:length(varnames)
        varname = varnames{i};
        var = vars{i};

        fprintf('Processing variable: %s\n', varname);
        
        output_bin_file = sprintf('%s/hycom_%s_%s.bin', output_dir, varname, target_date);
        if exist(output_bin_file) ~= 0
            fprintf('Target file %s already exists. Skip.\n', output_bin_file);
            continue;
        end


        [interpolated_var, xc,yc,zc] = hycom2modelgrid(varname, var, hycom_lon, hycom_lat, hycom_z, gd.xc, gd.yc, gd.xf, gd.yf, gd.zc, mask_dir, 0, 0);

        %interpolated_var = reshape(interpolated_var, gd.nxc, gd.nyc, gd.nzc);
        interpolated_vars{i} = interpolated_var;
    

        output_bin_file = sprintf('%s/hycom_%s_%s.bin', output_dir, varname, target_date);
        fprintf('Writing variable %s to file %s\n', varname, output_bin_file);
        wrslice(output_bin_file, interpolated_var, 1, fmt, Ieee);
    
        fprintf('Check if we still have missing data...\n');
        [status, Tmsk, Pmsk] = check3dmask(varname, output_bin_file, mask_dir, 1, gd, fmt, Ieee);

        if (status ~= 0)
            fprintf('Warning: There are mismatches of the variable %s\n', varname);

            var_err_chk.Tmsk = Tmsk;
            var_err_chk.Pmsk = Pmsk;
            output_err_chk_file = sprintf('%s/ch3ck3dmask_err_%s.mat', output_dir, varname);
            fprintf('Warning: Saving the data for error checking to file %s\n', output_err_chk_file);
            save(output_err_chk_file, 'var_err_chk');
            
        else
            fprintf('No missing data found.\n');
        end
    end

    status = 0;
    
end
