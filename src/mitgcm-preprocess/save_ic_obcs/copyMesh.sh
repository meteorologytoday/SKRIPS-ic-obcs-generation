#!/bin/csh -f

rm -rf grid.mat
rm -rf MASK

cp ../gen_mesh_no_mask/grid.mat .

mkdir MASK
cd MASK

cp ../../set-mask/run/Depth* .
cp ../../set-mask/run/DR* .
cp ../../set-mask/run/DX* .
cp ../../set-mask/run/DY* .
cp ../../set-mask/run/hFac* .
cp ../../set-mask/run/RA* .
cp ../../set-mask/run/RC* .
cp ../../set-mask/run/RF* .
cp ../../set-mask/run/X* .
cp ../../set-mask/run/Y* .

cd ..
