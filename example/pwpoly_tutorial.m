% Shows examples of usage for PiecewisePolynomial class.

% The construction of a PiecewisePolynomial object requires at least 4 things:
%   - A vector with monotonically increasing values delineating cell boundaries
%   - A number N indicating the number of degrees of freedom per cell (i.e.
%     denotes a polynomial of order N-1 on each cell)
%   - A vector x with the nodal locations of input data
%   - A vector y of same length of x with evaluations of the function at x
%
% The construction of the PiecewisePolynomial object does the following:
%   - Creates cells based on the cell_boundaries input, and allocates N degrees
%     of freedom to each cell. The number of cells is the property K of the
%     object.
%   - Bins the nodal locations x into the cells. Goes through each cell and uses
%     the (x,y) data in that cell to construct an interpolating polynomial.
%   - If fewer than N (x,y) data points lie in a cell, an interpolating
%     polynomial of maximal order given the cell data is constructed (but a total
%     of N degrees of freedom are still allocated).
%   - If more than N (x,y) data points lie in a cell, the interpolant does not
%     exist (unless the data comes exactly from a polynomial of degree N-1 or
%     less, which is not checked). In this case the constructor throws an error.
%
% If you wish to inspect any of the code for PiecewisePolynomial objects, it's
% in piecewise_interpolation/classes/@PiecewisePolynomial.

clear;
close all;

global packages;
repnodes = packages.piecewise_interpolation.grid_tools.replicate_local_nodes;

f_exact = @(x) sin(pi*x);
df_exact = @(x) pi*cos(pi*x);

% Construct a piecewise polynomial object:
x = sort(randn([100 1]));

local_nodes = [-0.5; 0; 0.5];  % -1 ---> left-hand side, +1 ---> right-hand side
% Generate a grid with local nodes inside each cell:
xf = repnodes(local_nodes, x);
xf = xf(:);
yf = f_exact(xf);

% Use the (xf,yf) data to construct a pwpoly on cells defined by x
temp.cell_boundaries = x;
temp.x = xf;
temp.y = yf;
temp.N = 3;  % Used three points per cell to generate data

f = PiecewisePolynomial(temp);  % The constructor

% You can take derivatives and integrals:
df = diff(f);  % Generates another PiecewisePolynomial
If = int(f);  % Computes the definite integral over the support of the function (exact)

% You can check that If is close to the exact answer -1/pi*(cos(max(x)) - cos(min(x)))

% You can add together and multiply PiecewisePolynomial's as well:
df_plus_f = df + f;
df_times_f = df*f;

% Even better, you can plot PiecewisePolynomial objects directly:
figure;
plot(df_plus_f, 'b.-'); hold on; plot(xf, f_exact(xf) + df_exact(xf), 'r');
title('Plot of f + df');
figure;
plot(df_times_f, 'b.-'); hold on; plot(xf, f_exact(xf).*df_exact(xf), 'r');
title('Plot of f * df')

% You can see that the class method 'plot' for PiecewisePolynomial uses it's
% Infinite Wisdom (c) to choose nice points for plotting. However, another nice
% thing about PiecewisePolynomial objects is that you can evaluate them exactly
% as you would a function: f(some_points) works just fine. So you can evaluate
% the function wherever you wish, or you can force it to plot at locations that
% you care about:
figure;
f1 = f_exact(xf).*df_exact(xf);
f2 = df_times_f(xf);
semilogy(xf, abs(f1-f2));
title('Plot of pointwise error of f * df');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PiecewisePolynomial objects have been shoehorned into some ENO codes as well.
% There's a new m-file, eno.eno_reconstruction that spits out a
% PiecewisePolynomial object that is the reconstruction from the input data
% (x,y). It forms this by using x as the cell_boundaries, and uses ENO to
% evaluate the function inside the cells.

eno = packages.eno.eno_reconstruction.handle;

exact_u = @(x) x>0;  % Yay step function

u = eno(x, exact_u(x), 'k', 4);  % Force piecewise quartic (for no good reason)

% Now u is a piecewise polynomial. We can do whatever we want with it. How close
% is u' to a delta function?
du = diff(u);
xff = linspace(min(x), max(x), 1e4);
figure; 
subplot(2,1,1);
plot(xff, u(xff));
title('Plot of step function');
subplot(2,1,2);
plot(xff, du(xff));
title('Plot of derivative of step function');

% If you *really* want the coefficients for the monomials, the method
% monomial_coefficients will extract them for you. It computes them on the fly,
% and if you experiment, you'll probably see that matlab will scream at you
% about ill-conditioning, because that's exactly what happens when you try to do
% stuff with monomials. The method returns an N x K matrix, where each column
% has the expansion coefficients for the first N monomials.

% u.monomial_coefficients
