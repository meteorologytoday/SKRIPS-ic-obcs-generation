#!/bin/bash

mitgcm_rundir=$1
output_dir=$2

if [ "$mitgcm_rundir" = "" ] ; then
    echo "Error: Missing the first argument"
fi

if [ "$output_dir" = "" ] ; then
    echo "Error: Missing the second argument"
fi


echo "Mitgcm rundir : $mitgcm_rundir"
echo "Output dir    : $output_dir"


echo "Removing output dir $output_dir"
echo "Will sleep for 10 seconds before continuing."
sleep 10


echo "Creating output dir $output_dir"
mkdir -p $output_dir
cd $output_dir

echo "Copying files..."

set -x
cp $mitgcm_rundir/Depth* .
cp $mitgcm_rundir/DR* .
cp $mitgcm_rundir/DX* .
cp $mitgcm_rundir/DY* .
cp $mitgcm_rundir/hFac* .
cp $mitgcm_rundir/RA* .
cp $mitgcm_rundir/RC* .
cp $mitgcm_rundir/RF* .
cp $mitgcm_rundir/X* .
cp $mitgcm_rundir/Y* .


echo "Finished."
