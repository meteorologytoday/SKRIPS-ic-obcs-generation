input_json_file = '/cw3e/mead/projects/csg102/t2hsu/AR_projects/project01/produce_ic_obcs_008/input.json'

input_json = read_json(input_json_file);
tool_root = input_json.tool_root;

%addpath([tool_root '/src']);
%addpath('/cw3e/mead/projects/csg102/t2hsu/scripps_kaust_model/MITgcm_c67m/utils/matlab/');


