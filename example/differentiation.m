% Example files for piecewise_interpolation

% This module can be used for piecewise polynomial (finite-difference)
% differentiation of any order.  The stencil is chosen by default, but if you
% dig through the implementation of the functions used in this example, and the
% helpstring of finite_difference.difference_stencil, you can figure out how to
% be very picky about what stencil is used for interpolation. Note that the ENO
% stencil-choosing rubric is essentially a special case of this picky stencil
% choosing.

clear;
close all;
global handles;
differentiate = handles.piecewise_interpolation.poly_differentiation;

% Let's do a pretty easy function:
f = @(x) atan(200*x);
df = @(x) 200./(1+(200*x).^2);

% (x,y) are the data points
x = linspace(-0.1,0.1,250);
y = f(x);
dy = df(x);

% Where do we want to differentiate?
z = linspace(-0.1,0.1,5e4);

% Let's look at the derivative, the default is cubic differentiation
plot(x,dy,'b.', z, differentiate(x,y,z), 'r-');
xlabel('x');
ylabel('Cubic interpolation - derivative');
temp = axis; axis([-0.01 0.01, temp(3:4)]);

figure();
% You can give optional input 'k': polynomial order of interpolation
plot(x,dy,'b.', z, differentiate(x,y,z, 'k',10), 'r-')
xlabel('x');
ylabel('10th order interpolation - derivative');
temp = axis; axis([-0.01 0.01, temp(3:4)]);
