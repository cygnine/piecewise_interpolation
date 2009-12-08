% Shows example for usage of convolution

clear;
close all;

global packages;
repnodes = packages.piecewise_interpolation.grid_tools.replicate_local_nodes;
% convolution = packages.piecewise_interpolation.convolution;

f_exact = @(x) cos(x); % The input function
Nq = 10; % Number of quadrature points per cell

% This is the convolution kernel:
g = @(x) (1-cos(x)).*log(2*(1-cos(x))) + 3/2*cos(x) - 1;

% Construct a piecewise polynomial object:
x = linspace(0, 2*pi, 100).';

local_nodes = [-2/3; -1/3; 1/3; 2/3];  % -1 ---> left-hand side, +1 ---> right-hand side
% Generate a grid with local nodes inside each cell:
xf = repnodes(local_nodes, x);
xf = xf(:);
yf = f_exact(xf);

% Use the (xf,yf) data to construct a pwpoly on cells defined by x
temp.cell_boundaries = x;
temp.x = xf;
temp.y = yf;
temp.N = 4;  % Used four points per cell to generate data

f = PiecewisePolynomial(temp);  % The constructor

temp = f.opoly_opt;

Gf = convolution(f.modal_coefficients, g, xf, 'cells', f.cell_boundaries, ...
     'alpha', temp.alpha, 'beta', temp.beta, 'Nq', Nq);
