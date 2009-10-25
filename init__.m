function[piecewise_interpolation] = init__()
% init__ -- Initialization file for piecewise_interpolation package
%
% [nodes] = init__()

piecewise_interpolation = recurse_files(pwd);
piecewise_interpolation.grid_tools = matlab_import('grid_tools');
piecewise_interpolation.debug = matlab_import('debug');

pwd_addpath('classes');
