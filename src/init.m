input_json_file = '';  %;'/home/t2hsu/projects/AR_airsea_coupling/case01/01_produce_case/input.json';

input_json = read_json(input_json_file);

tool_root = input_json.tool_root;

addpath([tool_root '/src']);
addpath([tool_root '/src/mitgcm-preprocess/gen_ic_obcs']);


