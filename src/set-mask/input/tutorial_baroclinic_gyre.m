% plot the figures using matlab

run '~/matlab_bin/pathdef.m'

% temperature
depth=rdmds('Depth');
maskInC=rdmds('maskInC');
maskInS=rdmds('maskInS');
maskInW=rdmds('maskInW');
DYC=rdmds('DYC');
DYG=rdmds('DYG');
DXC=rdmds('DXC');
DXG=rdmds('DXG');
DRC=rdmds('DRC');
DRF=rdmds('DRF');
PHrefC=rdmds('PHrefC');
PHrefF=rdmds('PHrefF');
RAC=rdmds('RAC');
RAS=rdmds('RAS');
RAW=rdmds('RAW');
RAZ=rdmds('RAZ');
RC=rdmds('RC');
RF=rdmds('RF');
RhoRef=rdmds('RhoRef');
XC=rdmds('XC');
XG=rdmds('XG');
YC=rdmds('YC');
YG=rdmds('YG');

% plot depth contour
figure()
surf(XC,YC,maskInC)
view([0,90])
title('maskInC')
xlabel('X')
ylabel('Y')
colorbar

% plot depth contour
figure()
surf(XC,YC,maskInS)
view([0,90])
title('maskInS')
xlabel('X')
ylabel('Y')
colorbar

% plot depth contour
figure()
surf(XC,YC,maskInW)
view([0,90])
title('maskInW')
xlabel('X')
ylabel('Y')
colorbar
