clear all;
close all;

addpath('./mitgcm-preprocess/gen_ic_obcs');

varname = 'T';
cbrng = [31, 35];
cmap = 'redblue';

filename = sprintf('test.bin')
fprintf('Plotting filename = %s\n', filename)

nx = 800;
ny = 400;
nz = 50;

fmt = 'real*4';
Ieee = 'b';

data = rdslice(filename, [nx, ny, nz], 1, fmt, Ieee);

figure;
ax = subplot(1, 1, 1);
s = pcolor(ax, data(:, :, 1)');
s.EdgeColor = 'none';
colorbar(ax);
caxis(cbrng);
[~, fname, fext] = fileparts(filename);
title(ax, sprintf('File: %s%s', fname, fext), 'Interpreter', 'none');
xlabel(ax, 'lon');
ylabel(ax, 'lat');

