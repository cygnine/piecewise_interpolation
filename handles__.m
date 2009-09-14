function[hs,pathadditions,name] = handles__()
% handles__ -- constructs path tree for piecewise_interpolation module
%
% [hs,pathadditions,name] = handles__()
%
%     Returns directory pointers for common module in hs. pathadditions is a
%     cell array with a string in each element indicated paths to add to the
%     global path structure. name is the name of this module.

name = 'piecewise_interpolation';

% This is by default
hs.base = fileparts(mfilename('fullpath'));

pathadditions = cell(0);
%pathadditions{end+1} = fullfile(hs.base,'classes');
