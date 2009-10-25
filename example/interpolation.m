% Example files for piecewise_interpolation

% This module can be used for piecewise polynomial interpolation of any order.
% The stencil is chosen by default, but if you dig through the implementation of
% the functions used in this example, and the helpstring of
% finite_difference.difference_stencil, you can figure out how to be very picky
% about what stencil is used for interpolation. Note that the ENO
% stencil-choosing rubric is essentially a special case of this picky stencil
% choosing.

clear;
close all;
global packages;
interpolate = packages.piecewise_interpolation.poly_interpolation;

% Let's do a pretty easy function:
f = @(x) atan(2000*x);

% (x,y) are the data points
x = linspace(-0.1,0.1,250);
y = f(x);

% Where do we want to interpolate?
z = linspace(-0.1,0.1,5e4);

% Let's look at the interpolant, the default is cubic
plot(x,y,'b.', z, interpolate(x,y,z), 'r-');
xlabel('x');
ylabel('Cubic interpolation');
axis([-0.01 0.01, -pi/2, pi/2]);

figure();
% You can give optional input 'k': polynomial order of interpolation
plot(x,y,'b.', z, interpolate(x,y,z, 'k',10), 'r-')
xlabel('x');
ylabel('10th order interpolation');
axis([-0.01 0.01, -pi/2, pi/2]);

% Note that this is not an advertisement for using 10th-order interpolation.
% This example is smooth, so nothing horrible happens, but when you start
% interpolating with high-order approximations when either your data is
% irregular or the function does not seem to be smooth, high-order interpolation
% starts giving horrible oscillations -- this is exactly what eno is supposed to
% help remedy.
