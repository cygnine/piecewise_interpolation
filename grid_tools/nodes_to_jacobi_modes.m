function[modes] = nodes_to_jacobi_modes(cell_boundaries,N,x,y,varargin)
% nodes_to_jacobi_modes -- converts nodal evaluations into modal coefficients
%
% [modes] = nodes_to_modes(cell_boundaries,N,x,y,{alpha=0,beta=0})
%     Creates local cellwise modal coefficient data via interpolation. The
%     input cell_boundaries is a nondecreasing vector of values delineating the
%     cell boundaries. The number of cells K is inferred from this. The total
%     number of degrees of freedom on each cell is N (i.e. an (N-1)st order
%     polynomial). The data is specified by the tuple (x,y). The modal
%     coefficients correspond to L2-normalized Jacobi polynomials.

global packages;
opt = packages.labtools.input_schema({'alpha', 'beta'}, {0,0}, [], varargin{:});
jac = packages.speclab.orthopoly1d.jacobi;
eval_jac = packages.speclab.orthopoly1d.eval_polynomial_standard;
compute_scaleshift = packages.piecewise_interpolation.grid_tools.compute_scaleshift.handle;

[recurrence_a,recurrence_b] = jac.coefficients.recurrence(N+1,opt);

x = x(:);
y = y(:);
[bin_sums, bin_id] = histc(x, cell_boundaries);

assert(length(x)==length(y), 'Nodal data x and y must have the same length');
assert(max(bin_sums)<=N,  'Too many data points lie inside a cell: cannot uniquely interpolate');
K = length(cell_boundaries(:))-1;

[jacobians, cell_shift] = compute_scaleshift(cell_boundaries);

% erase previous data
modes = zeros([N, K]);

r = (x - cell_shift(bin_id))./jacobians(bin_id);
polys = eval_jac(r, recurrence_a, recurrence_b, 0:(N-1));
for q = 1:K;
  flags = (bin_id==q);

  Nr = bin_sums(q);
  if Nr>0
    modes(1:Nr,q) = inv(polys(flags,1:Nr))*y(flags);
  end
end
