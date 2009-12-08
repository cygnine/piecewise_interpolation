function[h] = convolution(piece_poly, g, x, varargin)
% convolution -- Convolution of a piecewise polynomial
%
% [h] = convolution(piece_poly, x, {cells, alpha=0, beta=0, Nq=2*(N+1)})
%
%     Computes the linear convolution of apiecewise polynomial function with a
%     given function and evaluates it at the points x. The convolution kernel is
%     g. By "linear", we mean that the value of g outside the interval of
%     approximation is what is used to evaluate the kernel at locations outside
%     the interval. If you want it to behave like true linear truncated linear
%     convolution, force the input g to be zero outside the interval of
%     interest. If g is periodic, the periodic convolution will be computed.  No
%     points of x may coincide with the knots defining the piecewise function;
%     an error is thrown in this case. g should be vectorized in it's single
%     input.
%
%     The input piece_poly is an N x K matrix of coefficients. piece_poly(:,k)
%     are the N modal coefficients for the first N (local) Jacobi(alpha,beta)
%     polynomials on cell k. The cell boundaries are given by the (K+1)-length
%     vector cells, which is mandatory. Nq is the number of quadrature nodes
%     used on each interval. 
%
%     This function assumes that interval of approximation is [0,2*pi].

[N,K] = size(piece_poly);

global packages;
inputs = {'cells', 'alpha', 'beta', 'Nq'};
defaults = {[], 0, 0, 2*(N+1)};
opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});

opoly = packages.speclab.orthopoly1d;
jac = packages.speclab.orthopoly1d.jacobi;
gq = jac.quad.gauss_quadrature.handle;
evalpoly = opoly.eval_polynomial_standard.handle;

if isempty(opt.cells)
  error('You must define the cells over which the piecewise polynomial is defined');
end

Nq = opt.Nq;
x_size = size(x);
x = x(:);

M = length(x);

[garbage, bin] = histc(x, opt.cells);
bin(x==opt.cells(end)) = K;

% Generate global common vertices
jopt.alpha = 0;
jopt.beta = 0;
[r,w] = gq(Nq,jopt);

cell_scale = diff(opt.cells.')/2;
cell_shift = (opt.cells(2:end).' + opt.cells(1:(end-1)).')/2;
vertices = repmat(r, [1, K])*spdiags(cell_scale(:), 0, K,K);
vertices = vertices + repmat(cell_shift, [Nq 1]);

% Now quadrature can be done on each cell:
H = repmat(x, [1, Nq*K]) - repmat(vertices(:).', [M 1]);
H = g(H);
%H = (1-cos(H)).*log(2*(1-cos(H))) + 3/2*cos(H) - 1;
tempw = repmat(w,[1,K])*spdiags(cell_scale.',0,K,K);
H = H*spdiags(tempw(:),0,K*Nq, K*Nq);

% Evaluate the function at the nodes:
jopt.alpha = opt.alpha;
jopt.beta = opt.beta;
[recurrence_a, recurrence_b] = jac.coefficients.recurrence(N,jopt);
standard_polys = evalpoly(r,recurrence_a, recurrence_b,0:(N-1));
f_vertices = standard_polys*piece_poly;

% The convolution:
h = H*f_vertices(:);
