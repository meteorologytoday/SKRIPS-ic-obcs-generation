clear all;
close all;

addpath('./mitgcm-preprocess/gen_ic_obcs');

varname = 'S';
cbrng = [33, 35];
cmap = 'redblue';
input_dir = '/cw3e/mead/projects/csg102/t2hsu/AR_projects/project01/case01/caseInput/mitGCM_bin';
%/home/t2hsu/projects/SKRIPS-case-generation/sample/output_final';
input_dir = '/home/t2hsu/verify';

date_str='2017-02-01_00';
obcs_skips = 0
%ic_filename = sprintf('%s/initial_conditions/hycom_%s_%s.bin', input_dir, varname, date_str);
%bnd_east_filename = sprintf('%s/open_boundary_conditions/Rs_%s_obcs_east_%s.bin', input_dir, varname, date_str);
%bnd_west_filename = sprintf('%s/open_boundary_conditions/Rs_%s_obcs_west_%s.bin', input_dir, varname, date_str);
%bnd_north_filename = sprintf('%s/open_boundary_conditions/Rs_%s_obcs_north_%s.bin', input_dir, varname, date_str);
%bnd_south_filename = sprintf('%s/open_boundary_conditions/Rs_%s_obcs_south_%s.bin', input_dir, varname, date_str);

ic_filename = sprintf('%s/hycom_%s_%s.bin', input_dir, varname, date_str);
bnd_east_filename = sprintf('%s/obcs_%s_east.bin', input_dir, varname);
bnd_west_filename = sprintf('%s/obcs_%s_west.bin', input_dir, varname);
bnd_north_filename = sprintf('%s/obcs_%s_north.bin', input_dir, varname);
bnd_south_filename = sprintf('%s/obcs_%s_south.bin', input_dir, varname);


fprintf('ic_filename = %s\n', ic_filename);
fprintf('bnd_east_filename = %s\n', bnd_east_filename);
fprintf('bnd_west_filename = %s\n', bnd_west_filename);
fprintf('bnd_north_filename = %s\n', bnd_north_filename);
fprintf('bnd_south_filename = %s\n', bnd_south_filename);

magnify = 1000;
nx = 800;
ny = 400;
nz = 50;

fmt = 'real*4';
Ieee = 'b';

ic = rdslice(ic_filename, [nx, ny, nz], 1, fmt, Ieee);

rec = 1 + obcs_skips
bnd_east  = rdslice(bnd_east_filename,  [ny, nz], rec, fmt, Ieee);
bnd_west  = rdslice(bnd_west_filename,  [ny, nz], rec, fmt, Ieee);
bnd_north = rdslice(bnd_north_filename, [nx, nz], rec, fmt, Ieee);
bnd_south = rdslice(bnd_south_filename, [nx, nz], rec, fmt, Ieee);

figure;
ax = subplot(1, 1, 1);
s = pcolor(ax, ic(:, :, 1)');
s.EdgeColor = 'none';
colorbar(ax);
caxis(cbrng);
[~, fname, fext] = fileparts(ic_filename);
title(ax, sprintf('Init Cond: %s%s', fname, fext), 'Interpreter', 'none');
xlabel(ax, 'lon');
ylabel(ax, 'lat');

figure; clf;
ax = subplot(1, 1, 1);
s = pcolor(ax, flipud(bnd_west'));
s.EdgeColor = 'none';
colorbar(ax);
caxis(cbrng);
[~, fname, fext] = fileparts(bnd_west_filename);
title(ax, sprintf('West Bnd: %s%s', fname, fext), 'Interpreter', 'none');
xlabel(ax, 'lat');
ylabel(ax, 'z');
hold off;


figure(3); clf;
ax = subplot(1, 1, 1);
s = pcolor(ax, flipud(bnd_east'));
s.EdgeColor = 'none';
colorbar(ax);
caxis(cbrng);
[~, fname, fext] = fileparts(bnd_east_filename);
title(ax, sprintf('East Bnd: %s%s', fname, fext), 'Interpreter', 'none');
xlabel(ax, 'lat');
ylabel(ax, 'z');

figure(4); clf;
ax = subplot(1, 1, 1);
s = pcolor(ax, flipud(bnd_south'));
s.EdgeColor = 'none';
colorbar(ax);
caxis(cbrng);
[~, fname, fext] = fileparts(bnd_south_filename);
title(ax, sprintf('South Bnd: %s%s', fname, fext), 'Interpreter', 'none');
xlabel(ax, 'lon');
ylabel(ax, 'z');


figure(5); clf;
ax = subplot(1, 1, 1);
s = pcolor(ax, flipud(bnd_north'));
s.EdgeColor = 'none';
colorbar(ax);
caxis(cbrng);
[~, fname, fext] = fileparts(bnd_north_filename);
title(ax, sprintf('North Bnd: %s%s', fname, fext), 'Interpreter', 'none');
xlabel(ax, 'lon');
ylabel(ax, 'z');

