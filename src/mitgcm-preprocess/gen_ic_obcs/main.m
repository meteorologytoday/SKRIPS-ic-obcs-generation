clearvars -except input_json_file tool_root run_steps;

input_json = read_json(input_json_file);
gridgen_nml_file = sprintf('%s/%s', input_json.workspace, input_json.gridgen_nml_file);
nml_grid_init = read_namelist(gridgen_nml_file, 'GRID_INIT');

final_output_dir = input_json.output_dir;

grid_dir = nml_grid_init.data_dir;
hycom_data_dir = sprintf('%s/%s', input_json.workspace, input_json.hycom.data_dir);
mask_dir = sprintf('%s/mitgcm_mask', input_json.workspace);


fprintf('Making output dir: %s\n', final_output_dir);
mkdir(final_output_dir);

target_date_ic = sprintf('%s_00', input_json.hycom.start_date);
gen_ic_hycom(grid_dir, mask_dir, hycom_data_dir, input_json.output_dir, target_date_ic);


