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

global handles;
opt = handles.common.input_schema({'alpha', 'beta'}, {0,0}, [], varargin{:});
jac = handles.speclab.orthopoly1d.jacobi;
eval_jac = handles.speclab.orthopoly1d.eval_polynomial;
sss = handles.speclab.common.standard_scaleshift_1d.handle;
compute_scaleshift = handles.piecewise_interpolation.grid_tools.compute_scaleshift.handle;

[recurrence_a,recurrence_b] = jac.coefficients.recurrence(N+1,opt);

[bin_sums, bin_id] = histc(x, cell_boundaries);
x = x(:);
y = y(:);

assert(length(x)==length(y), 'Nodal data x and y must have the same length');
assert(max(bin_sums)<=N,  'Too many data points lie inside a cell: cannot uniquely interpolate');
K = length(cell_boundaries(:))-1;

[jacobians, cell_shift] = compute_scaleshift(cell_boundaries);

% erase previous data
modes = zeros([N, K]);

for q = 1:K;
  flags = (bin_id==q);
  r = sss(x(flags), 'scale', jacobians(q), 'shift', cell_shift(q));

  Nr = length(r);
  if Nr>0
    polys = eval_jac(r, recurrence_a, recurrence_b, 0:(Nr-1));
    modes(1:Nr,q) = inv(polys)*y(flags);
  end
end
