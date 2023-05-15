clearvars -except input_json_file tool_root run_steps;
close all
clc

input_json = read_json(input_json_file);
opath = sprintf('%s/%s', input_json.workspace, input_json.hycom.data_dir);
OpenDAP_URL = 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0';

start_date = input_json.hycom.start_date;
end_date   = input_json.hycom.end_date;


fprintf('Output directory: %s\n', opath);
fprintf('start_date : %s\n', start_date);
fprintf('end_date   : %s\n', end_date);

mkdir(opath);

% Boundaries for all dates
regions = {'north', 'east', 'south', 'west'};
for i = 1:length(regions)
    region = regions{i};
    display(['Grabbing boundary data of region: ' region]);
    getHycomData(start_date, end_date, region, opath, OpenDAP_URL);
end 

% Initial condition for the first date
disp('Grabbing initial condition data');
getHycomData(start_date, start_date, 'all', opath, OpenDAP_URL);


disp('Download data finished.');



