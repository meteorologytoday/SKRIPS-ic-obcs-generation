function report(varargin)
%% This is a simple device for writing to
%% stdout (the terminal) with out getting
%% annoying formatting problems like "ans ="

%% Takes the same arguments as SPRINTF()

fprintf(1,varargin{:})
