clear; clc;

run('init.m');

gridgen_nml_file = input_json.gridgen_nml_file;
nml_grid_init = read_namelist(sprintf('%s/%s', input_json.workspace, gridgen_nml_file), 'GRID_INIT');

data_dir  = nml_grid_init.data_dir;

start_date = input_json.hycom.start_date;
end_date   = input_json.hycom.end_date;

interp_data_dir = sprintf('%s/%s', input_json.workspace, input_json.interp_data_dir);
ic_dir = sprintf('%s/initial_conditions', interp_data_dir);
obc_dir = sprintf('%s/open_boundary_conditions', interp_data_dir);

final_output_data_dir = sprintf('%s/%s', input_json.workspace, input_json.output_dir);


fprintf('Start date: %s\n', start_date);
fprintf('End   date: %s\n', end_date);



fprintf('Making output dir: %s\n', final_output_data_dir);
mkdir(final_output_data_dir);


% Make open boundary condition files
bnds = {'north', 'south', 'west', 'east'};
for i = 1:length(bnds)
    bnd = bnds{i};
    fprintf('Making %s boundary condition file.\n', bnd);
    gen_obcs_blob(bnd, obc_dir, final_output_data_dir, start_date, end_date);
end

% Copy necessary files...
copyfile(sprintf('%s/bathymetry_ar_50v.bin', data_dir), final_output_data_dir); 

varnames = {'T', 'S', 'U', 'V'};
for i = 1:length(varnames)
    varname = varnames{i};
    copyfile(sprintf('%s/hycom_%s_%s_00.bin', ic_dir, varname, start_date), final_output_data_dir); 
end



