
run('init.m');

run_steps = {'step01', 'step02', 'step03', 'step04'};

%run_steps = { 'step02', 'step03', 'step04'};


for i = 1:length(run_steps)

    run_step = run_steps{i};
    display(['Running step: ' run_step]);
    if strcmp(run_step, 'step01') == 1
        % Use wavewatch3 (ww3) utility to generated needed bathymetry information
        % and intermediate files. Most importantly, `depth_mitgcm.mat`.
        run([ tool_root '/src/step01_ww3.m' ])

    elseif strcmp(run_step, 'step02') == 1
        % Use the ww3-generated intermediate files to produce
        % bathymetry binary file,  and grid.mat 
        run([ tool_root '/src/mitgcm-preprocess/gen_mesh_no_mask/genMesh.m' ])

    elseif strcmp(run_step, 'step03') == 1
        % Download Hycom data: initial condition (1 data), boundary
        % conditions (multiple data)
        run([ tool_root '/src/mitgcm-preprocess/get_hycom/main.m' ])

    elseif strcmp(run_step, 'step04') == 1
        run([ tool_root '/src/mitgcm-preprocess/get_hycom/main2.m'])
        %run([ tool_root '/src/mitgcm-preprocess/get_hycom/gen_obcs_noncorrect.m'])
    else
        error('Unknown step: %s', run_step);
    end
end




