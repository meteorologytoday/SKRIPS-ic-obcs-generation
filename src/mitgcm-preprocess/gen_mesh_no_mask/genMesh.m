clearvars -except input_json_file tool_root;

input_json = read_json(input_json_file);
gridgen_nml_file = sprintf('%s/%s', input_json.workspace, input_json.gridgen_nml_file)

nml_grid_init = read_namelist(gridgen_nml_file, 'GRID_INIT');

data_dir = nml_grid_init.data_dir
nml_rect_nml = read_namelist([ data_dir '/' nml_grid_init.fname '.meta' ], 'RECT_NML');



nxc = nml_rect_nml.nx;
nyc = nml_rect_nml.ny;
nzc = 50;

deltaX = nml_rect_nml.sx;
deltaY = nml_rect_nml.sy;

% These are the edges
lat_lo = nml_rect_nml.y0 - deltaY / 2;
lat_hi = lat_lo + deltaY * nyc;

lon_lo = nml_rect_nml.x0 - deltaX / 2;
lon_hi = lon_lo + deltaX * nxc;

deg2m = 4*(27794+0.3683+3.8010e-5-3.0093e-10); % ??????? = 1.1118e5

dlat = zeros(nxc,1) + deltaY;
dlon = zeros(nyc,1) + deltaX;

DXF = zeros(nxc,nyc)+deltaX*deg2m;
DYF = zeros(nxc,nyc)+deltaY*deg2m;
zA  = zeros(nxc,nyc)+deltaX*deg2m*deltaY*deg2m;
xc = lon_lo+0.5*deltaX:deltaX:lon_hi-0.5*deltaX;xc = xc';
xf = lon_lo:deltaX:lon_hi-deltaX;xf = xf';
yc = lat_lo+0.5*deltaY:deltaY:lat_hi-0.5*deltaY;
yf = lat_lo:deltaY:lat_hi-deltaY;

dz =-1*[4.,5.,5.,5.,6.,6.,7.,8.,8.,9.,...
10.,11.,12.,13.,14.,15.,17.,18.,20.,21.,...
23.,25.,28.,30.,32.,35.,38.,42.,45.,49.,...
53.,58.,63.,68.,74.,81.,88.,95.,103.,112.,...
122.,150.,150.,150.,200.,200.,200.,250.,250.,300.];
zf = zeros(1,50);
zc = zeros(1,50);

zf(1) = -4;
zc(1) = -2;
for i = 2:50
  zf(i) = zf(i-1)+dz(i);
  zc(i) = zf(i-1)+0.5*dz(i);
end

% % ROMS mesh streching
% % Song and Haidvogel(1994)
% % https://www.myroms.org/wiki/Vertical_S-coordinate
% sigma_f = 0:1/nzc:1;
% sigma_c = 0.5*1/nzc:1/nzc:1-0.5*1/nzc;
% theta_b = 0.0;
% theta_s = 8.0;
% c_sigma_f = (1-theta_b)*sinh(theta_s*sigma_f)/sinh(theta_s)...
%           +theta_b*(tanh(theta_s*(sigma_f+0.5))/2.0/tanh(0.5*theta_s)-0.5);
% c_sigma_c = (1-theta_b)*sinh(theta_s*sigma_c)/sinh(theta_s)...
%           +theta_b*(tanh(theta_s*(sigma_c+0.5))/2.0/tanh(0.5*theta_s)-0.5);
% 
% zf = c_sigma_f*(-6400);
% zc = c_sigma_c*(-6400);
% dz = diff(zf);

clear c_sigma_* delta* theta_* sigma_* deg2m

for i = 1:nzc
  [i, dz(i), sum(dz(1:i))]
end

% plot dz and other parameters
output_filename = sprintf('%s/grid.mat', data_dir);
fprintf("Going to output: %s \n", output_filename);
fileID = fopen([ data_dir '/delZ.txt' ],'w');
fprintf(fileID,'%f, \n',dz);
fclose(fileID);
save(output_filename);

% plot bathymetry
depth_mitgcm_filename = sprintf('%s/depth_mitgcm.mat', data_dir);
fprintf('Going to load: %s \n', depth_mitgcm_filename);
load(depth_mitgcm_filename);

output_filename = sprintf('%s/bathymetry_ar_50v.bin', data_dir);
fprintf('Going to output: %s\n', output_filename);
fileID = fopen(output_filename,'w','b');
fwrite(fileID, depth_mitgcm', 'real*4');
fclose(fileID);
