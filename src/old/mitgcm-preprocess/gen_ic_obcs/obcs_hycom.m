%% This is a script that executes a series of commands for
%% extracting values from the global model outputs,
%% interpolating, and adjusting (for conservation).

%% Created  03/01/02 by vthierry@ucsd.edu
clear all;
close all
clc

run ~/matlab_bin/pathdef.m

%% Extracts the values form the global model
%% ALL FOUR SLICES (NSEW) INTERPOLATED TO MIT GRID(VERT/HOR)
gen_obcs_hycom

%% Calulates the mass transport and the heat transport
%% at the northern and southern boundaries
%% in the global model ...

%% FOR NORTH AND SOUTH SLICES
hycomtransp_NS

%% FOR EAST AND WEST SLICES
hycomtransp_EW

%% ... and accordingly corrects the
%% interpolated fields on the higher resolution grid.
cor_obcs_hycom_transp

