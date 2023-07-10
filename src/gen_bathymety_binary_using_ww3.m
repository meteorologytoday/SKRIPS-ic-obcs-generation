
run('init.m')


clearvars -except input_json_file tool_root;

clc;
close all;

set(groot,'DefaultFigureColormap',jet)

fprintf('Going to read JSON file: %s\n', input_json_file);

input_json = read_json(input_json_file);

gridgen_nml_file = input_json.gridgen_nml_file;
gridgen_nml_file_fullpath = sprintf('%s/%s', input_json.workspace, gridgen_nml_file);

nml_grid_init = read_namelist(gridgen_nml_file_fullpath, 'GRID_INIT');
nml_outgrid   = read_namelist(gridgen_nml_file_fullpath, 'OUTGRID');

bin_dir  = nml_grid_init.bin_dir;
ref_dir  = nml_grid_init.ref_dir;
data_dir  = nml_grid_init.data_dir;

dx = nml_outgrid.dx;
dy = nml_outgrid.dy;

lon_west  = nml_outgrid.lon_west;
lon_east  = nml_outgrid.lon_east;

lat_south = nml_outgrid.lat_south;
lat_north = nml_outgrid.lat_north;

fname = nml_grid_init.fname;
fnameb = nml_grid_init.fnameb;

edge_lon_west = lon_west - dx / 2; 
edge_lon_east = lon_east + dx / 2; 
edge_lat_south = lat_south - dy / 2; 
edge_lat_north = lat_north + dy / 2; 

% LAKE_TOL is the threshold count of the connected water surface grids.
% If the count of the connected grids is lower than LAKE_TOL, it will be
% categorized as lake rather than the ocean.
if isfield(input_json, 'lake_grids_threshold')
    LAKE_TOL = input_json.lake_grids_threshold;
else
    LAKE_TOL = 10000;
end



fprintf('fname        : %s\n', fname);
fprintf('fnameb       : %s\n', fnameb);
fprintf('bin_dir      : %s\n', bin_dir);
fprintf('ref_dir      : %s\n', ref_dir);
fprintf('data_dir     : %s\n', data_dir);
fprintf('dlon : %.2f\n', dx);
fprintf('dlat : %.2f\n', dy);

fprintf('lon_west  : %.2f\n', lon_west);
fprintf('lon_east  : %.2f\n', lon_east);
fprintf('lat_south : %.2f\n', lat_south);
fprintf('lat_north : %.2f\n', lat_north);

fprintf('edge_lon_west  : %.2f\n', edge_lon_west);
fprintf('edge_lon_east  : %.2f\n', edge_lon_east);
fprintf('edge_lat_south : %.2f\n', edge_lat_south);
fprintf('edge_lat_north : %.2f\n', edge_lat_north);

fprintf('Lake grids threshold: %d\n', LAKE_TOL);


mkdir(data_dir);

% Make softlinks to have global file
basegrid_dir = sprintf('%s/src/ww3_gridgen/data', tool_root);
fprintf('Making soft links of basegrid files from %s to data_dir\n', basegrid_dir);
system(sprintf('ln -s %s/%s.meta %s/', basegrid_dir, input_json.basegrid, data_dir));
system(sprintf('ln -s %s/%s.mask %s/', basegrid_dir, input_json.basegrid, data_dir));

lon1d = (lon_west:dx:lon_east);
lat1d = (lat_south:dy:lat_north);
[lon,lat] = meshgrid(lon1d,lat1d);

fprintf('Number of grids in (lon, lat) = (%d, %d)\n', length(lon1d), length(lat1d) );


addpath(bin_dir,'-END');
load([ref_dir,'/coastal_bound_inter.mat']);
bound

figure(1);clf;
for i = 1:1000
  plot(bound(i).x,bound(i).y,'.','MarkerSize',0.5);
  hold all;
end;

CUT_OFF = 0.0; % Cut-off depth to distinguish between dry & wet cells
LIM_BATHY = 0.4; % Base bathymetry cells needing to be wet for the target cell to be considered wet.
DRY_VAL = 999999; % Depth value set for fry cells

ref_grid = 'gebco360'; % Name of the file without the '.nc' extension
xvar = 'lon'; % Name of the variable defining longitudes in file
yvar = 'lat'; % Name of the variable defining latitudes in file
zvar = 'elevation'; % Name of the variable defining depths in file

depth = generate_grid('rect',lon,lat,ref_dir,ref_grid,LIM_BATHY,CUT_OFF,DRY_VAL,xvar,yvar,zvar);
figure(1);clf;
d=depth;d(d==DRY_VAL)=nan; pcolor(lon,lat,d); shading flat; colorbar

m1 = ones(size(depth));
m1(depth == DRY_VAL) = 0;

figure(1);clf;
pcolor(lon,lat,m1);shading flat;caxis([0 3]);colorbar

% compute the boundaries
lon_start = min(min(lon))-dx;
lon_end = max(max(lon))+dx;
lat_start = min(min(lat))-dy;
lat_end = max(max(lat))+dy;
coord = [lat_start lon_start lat_end lon_end];
MIN_DIST = 10; % minimum distance between edge of polygon and boundary

[b,N] = compute_boundary(coord,bound);

figure(1);clf;
for i = 1:N
  plot(b(i).x,b(i).y);
  hold on;
end;

% split up boundary polygons
SPLIT_LIM = 2;
b_split = split_boundary(b,SPLIT_LIM,MIN_DIST);
Nb = length(b_split);
figure(1);clf;
for i = 1:Nb
  plot(b_split(i).x,b_split(i).y);
  hold on;
end;

% clean up the initial mask
LIM_VAL = 0.5;
OFFSET = max([dx dy]);
m2 = clean_mask(lon,lat,m1,b_split,LIM_VAL,OFFSET);

figure(1);clf;
pcolor(lon,lat,m2);shading flat;caxis([0 3]);colorbar

% remove artificially generated lakes
% LAKE_TOL is a count of grids that is labeled water
% The function `remove_lake` will first group connected
% water pixels with labels. Then, it removes the water
% body that is too small, measured by LAKE_TOL. Please
% read the `remove_lake.m` file for detailed usage.
%LAKE_TOL = 10000;
LAKE_TOL = 10000;
IS_GLOBAL = 0;
[m4, mask_map] = remove_lake(m2,LAKE_TOL,IS_GLOBAL);

figure(1);clf;
pcolor(lon,lat,m4);shading flat;caxis([0 3]);colorbar

figure(1);clf;
pcolor(lon,lat,mask_map);shading flat;caxis([-1 6]);colorbar

% generate obstruction grids
OBSTR_OFFSET = 1;
[sx1,sy1] = create_obstr(lon,lat,b,m4,OBSTR_OFFSET,OBSTR_OFFSET);

sx1(find(m4==0))=NaN;
sy1(find(m4==0))=NaN;
figure(1);clf;
pcolor(lon,lat,sx1);shading flat;caxis([0 1]); colorbar
figure(2);clf;
pcolor(lon,lat,sy1);shading flat;caxis([0 1]); colorbar

% save files
depth_scale = 1000;
obstr_scale = 100;
d = round((depth)*depth_scale);
write_ww3file([data_dir,'/',fname,'.bot'],d);
write_ww3file([data_dir,'/',fname,'.mask_nobound'],m4);
d1 = round((sx1)*obstr_scale);
d2 = round((sy1)*obstr_scale);
write_ww3obstr([data_dir,'/',fname,'.obst'],d1,d2);
write_ww3meta([data_dir,'/',fname], gridgen_nml_file_fullpath, 'RECT', lon,lat,1/depth_scale,1/obstr_scale,1.0);

depth_mitgcm = d/1000;
save([ data_dir, '/', 'depth_mitgcm_before_masked.mat'], "depth_mitgcm");

depth_mitgcm = m4.*d/1000;
save([ data_dir, '/', 'depth_mitgcm.mat'], "depth_mitgcm");

output_filename = sprintf('%s/bathymetry.bin', data_dir);
fprintf('Going to output: %s\n', output_filename);
fileID = fopen(output_filename,'w','b');
fwrite(fileID, depth_mitgcm', 'real*4');
fclose(fileID);

% create_grid('/home/rus043/test_ww3_grid/gridgen.ZA-7M.nml');
%create_boundary(gridgen_nml_file_fullpath);
