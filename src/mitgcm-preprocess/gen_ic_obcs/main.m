clearvars -except input_json_file tool_root run_steps;

input_json = read_json(input_json_file);
start_date = input_json.hycom.start_date;
end_date   = input_json.hycom.end_date;

hycom_data_dir = sprintf('%s/%s', input_json.workspace, input_json.hycom.data_dir);
mask_dir = sprintf('%s/mitgcm_mask', input_json.workspace);

gridgen_nml_file = sprintf('%s/%s', input_json.workspace, input_json.gridgen_nml_file);
nml_grid_init = read_namelist(gridgen_nml_file, 'GRID_INIT');

grid_dir = nml_grid_init.data_dir;

final_output_dir = sprintf('%s/%s', input_json.workspace, input_json.output_dir);
fprintf('Making output dir: %s\n', final_output_dir);
mkdir(final_output_dir);



% Make initial condition file
target_date_ic = sprintf('%s_00', input_json.hycom.start_date);
fprintf('Making initial condition date: %s\n', target_date_ic);
%gen_ic_hycom(grid_dir, mask_dir, hycom_data_dir, final_output_dir, target_date_ic);


% Make boundary condition files
time = [datenum(start_date, format):1:datenum(end_date, format)];
nt = length(time);
bnds = {'north', 'south', 'east', 'west'};
thickness = 1;
corner_island_flag=0;

fprintf('Making boundary conditions... \n');
fprintf('thickness = %d\n', thickness);
fprintf('corner_island_flag = %d\n', corner_island_flag);
for t = 1:length(time)
    target_date = datestr(time(t), 'yyyy-mm-dd_hh');

    for i = 1:length(bnds)
        bnd = bnds{i};

        fprintf('Making boundary condition bnd=%s, date=%s\n', bnd, target_date);
        gen_obcs_hycom(bnd, thickness, corner_island_flag, grid_dir, mask_dir, hycom_data_dir, final_output_dir, t, target_date);

    end
end


    



