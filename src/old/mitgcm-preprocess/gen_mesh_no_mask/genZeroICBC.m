clear all;
load('grid_as.mat');

iniU = zeros(nxc,nyc,nzc);
iniV = zeros(nxc,nyc,nzc);
iniT = 15+zeros(nxc,nyc,nzc);
iniS = 35+zeros(nxc,nyc,nzc);

mean(iniT(:))

fileID = fopen('iniU.bin','w','b');
fwrite(fileID,iniU,'real*4');
fclose(fileID);
fileID = fopen('iniV.bin','w','b');
fwrite(fileID,iniV,'real*4');
fclose(fileID);
fileID = fopen('iniT.bin','w','b');
fwrite(fileID,iniT,'real*4');
fclose(fileID);
fileID = fopen('iniS.bin','w','b');
fwrite(fileID,iniS,'real*4');
fclose(fileID);

nt = 100
obcU_n = zeros(nxc,nyc,nt);
obcU_s = zeros(nxc,nyc,nt);
obcU_e = zeros(nxc,nyc,nt);
obcU_w = zeros(nxc,nyc,nt);

obcV_n = zeros(nxc,nyc,nt);
obcV_s = zeros(nxc,nyc,nt);
obcV_e = zeros(nxc,nyc,nt);
obcV_w = zeros(nxc,nyc,nt);

obcT_n = 15+zeros(nxc,nyc,nt);
obcT_s = 15+zeros(nxc,nyc,nt);
obcT_e = 15+zeros(nxc,nyc,nt);
obcT_w = 15+zeros(nxc,nyc,nt);

obcS_n = 35+zeros(nxc,nyc,nt);
obcS_s = 35+zeros(nxc,nyc,nt);
obcS_e = 35+zeros(nxc,nyc,nt);
obcS_w = 35+zeros(nxc,nyc,nt);

fileID = fopen('obcU_n.bin','w','b');
fwrite(fileID,obcU_n,'real*4');
fclose(fileID);
fileID = fopen('obcU_s.bin','w','b');
fwrite(fileID,obcU_s,'real*4');
fclose(fileID);
fileID = fopen('obcU_e.bin','w','b');
fwrite(fileID,obcU_e,'real*4');
fclose(fileID);
fileID = fopen('obcU_w.bin','w','b');
fwrite(fileID,obcU_w,'real*4');
fclose(fileID);

fileID = fopen('obcV_n.bin','w','b');
fwrite(fileID,obcV_n,'real*4');
fclose(fileID);
fileID = fopen('obcV_s.bin','w','b');
fwrite(fileID,obcV_s,'real*4');
fclose(fileID);
fileID = fopen('obcV_e.bin','w','b');
fwrite(fileID,obcV_e,'real*4');
fclose(fileID);
fileID = fopen('obcV_w.bin','w','b');
fwrite(fileID,obcV_w,'real*4');
fclose(fileID);

fileID = fopen('obcT_n.bin','w','b');
fwrite(fileID,obcT_n,'real*4');
fclose(fileID);
fileID = fopen('obcT_s.bin','w','b');
fwrite(fileID,obcT_s,'real*4');
fclose(fileID);
fileID = fopen('obcT_e.bin','w','b');
fwrite(fileID,obcT_e,'real*4');
fclose(fileID);
fileID = fopen('obcT_w.bin','w','b');
fwrite(fileID,obcT_w,'real*4');
fclose(fileID);

fileID = fopen('obcS_n.bin','w','b');
fwrite(fileID,obcS_n,'real*4');
fclose(fileID);
fileID = fopen('obcS_s.bin','w','b');
fwrite(fileID,obcS_s,'real*4');
fclose(fileID);
fileID = fopen('obcS_e.bin','w','b');
fwrite(fileID,obcS_e,'real*4');
fclose(fileID);
fileID = fopen('obcS_w.bin','w','b');
fwrite(fileID,obcS_w,'real*4');
fclose(fileID);
