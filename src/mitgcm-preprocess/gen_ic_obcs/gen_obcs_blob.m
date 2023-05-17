function gen_obcs_blob(bnd, input_dir, output_dir, start_date, end_date)
%% `bnd`          : side of the boundary. Accepting 'north', 'south', 'west', 'east'.
%% `input_dir`    : should contain the file grid.mat
%% `output_dir`   : the output directory
%% `start_date`   : should be a formatted string `yyyy-mm-dd_hh`
%% `end_date`     : should be a formatted string `yyyy-mm-dd_hh`

    time = [datenum(start_date, format):1:datenum(end_date, format)];
    nt = length(time);
    varnames = {'T', 'S', 'U', 'V'};
    
    for i = 1:length(varnames)

        varname = varnames{i};
        output_filename = sprintf('%s/obcs_%s_%s.bin', output_dir, varname, bnd);

        fprintf('Generating file: %s\n', output_filename);
        delete(output_filename);
        system(sprintf('touch %s', output_filename));
        for t = 1:length(time)
            target_date = datestr(time(t), 'yyyy-mm-dd_hh');
            obcs_filename = sprintf('%s/obcs_%s_%s_%s.bin', input_dir, varname, bnd, target_date);
            fprintf('Dumping file: %s\n', obcs_filename);
            system(sprintf('cat %s >> %s', obcs_filename, output_filename));
        end

    end

end
