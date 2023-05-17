tool_root = '/home/t2hsu/projects/SKRIPS-case-generation';
input_json_file = '/home/t2hsu/projects/SKRIPS-case-generation/sample/input.json';

addpath([tool_root '/src']);

input_json = read_json(input_json_file);

run_steps = {'step06', 'step07'};
for i = 1:length(run_steps)

    run_step = run_steps{i};
    display(['Running step: ' run_step]);
    if strcmp(run_step, 'step05') == 1
        
        system(sprintf('cp -r %s/src/mitgcm-preprocess/set-mask %s', tool_root, input_json.workspace))
        
        disp('Now you need to compile mitgcm.')
        disp('Now you need to copy bathmetry binary file into the run folder.')
        disp('Now run mitgcm')

        disp('I assume everything is done.')
    elseif strcmp(run_step, 'step06') == 1
        
        % copy mesh

    elseif strcmp(run_step, 'step07') == 1
        
        % Generate initial condition and boundary conditions with interpolation
        run([tool_root '/src/mitgcm-preprocess/gen_ic_obcs/main.m'])
        
    else
        error('Unknown step: %s', run_step);
    end
end




