input_json_file = '/cw3e/mead/projects/csg102/t2hsu/AR_projects/project01/case04_ISOLATED_CYCLONE/produice_ic_obcs_case04/input.json'
    
tool_root = '/cw3e/mead/projects/csg102/t2hsu/SKRIPS-case-generation'

addpath([tool_root '/src/matlab_utils']);


input_json = read_json(input_json_file);
tool_root = input_json.tool_root;
%addpath('/cw3e/mead/projects/csg102/t2hsu/scripps_kaust_model/MITgcm_c67m/utils/matlab/');


