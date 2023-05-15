tool_root = '/home/t2hsu/projects/SKRIPS-case-generation';
input_json_file = '/home/t2hsu/projects/SKRIPS-case-generation/sample/input.json';

addpath([tool_root '/src']);

%run([ tool_root '/src/step01_ww3.m' ])
%run([ tool_root '/src/mitgcm-preprocess/gen_mesh_no_mask/genMesh.m' ])
run([ tool_root '/src/mitgcm-preprocess/get_hycom/main.m' ])

%run([ tool_root '/src/mitgcm-preprocess/get_hycom/gen_obcs_noncorrect.m'])