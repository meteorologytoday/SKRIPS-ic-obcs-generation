clearvars -except input_json_file tool_root run_steps;
close all
clc

input_json = read_json(input_json_file);
FMT_json = read_json('../FMT.json');

start_date = input_json.hycom.start_date;
end_date   = input_json.hycom.end_date;

input_dir  = sprintf('%s/%s', input_json.workspace, input_json.hycom.data_dir);
output_dir = sprintf('%s/%s', input_json.workspace, input_json.hycom.ic_obcs_dir);

disp('Converting mat files into binary files.');

fprintf('start_date : %s\n', start_date);
fprintf('end_date   : %s\n', end_date);


