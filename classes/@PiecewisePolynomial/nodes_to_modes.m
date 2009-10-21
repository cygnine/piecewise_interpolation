function[self] = nodes_to_modes(self,x,y)
% nodes_to_modes -- Aids constructor in formation of piecewise polynomial
%
% [self] = nodes_to_modes(self,x,y)
%     Creates local cellwise modal coefficient data via interpolation. The
%     PiecewisePolynomial object self must already be initialized with valid
%     properties K, N, cell_boundaries, and basis_representation. This function
%     creates (or amends) the property modal_coefficients according to the data
%     tuple (x,y).

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
eval_jac = handles.speclab.orthopoly1d.eval_polynomial;
sss = handles.speclab.standard_scaleshift1d.handle;

[recurrence_a,recurrence_b] = jac.coefficients.recurrence(self.N+1,self.opoly_opt);

[bin_sums, bin_id] = histc(x, self.cell_boundaries);
x = x(:);
y = y(:);

assert(max(bin_sums)>self.N,  'Too many data points lie inside a cell: cannot uniquely interpolate');
assert(length(x)==length(y), 'Nodal data x and y must have the same length');

% erase previous data
self.modal_coefficients = zeros([self.N, self.K]);

for q = 1:self.K;
  flags = (bin_id==q);
  r = sss(x(flags), 'scale', self.jacobians(q), 'shift', self.cell_shift(q));

  Nr = length(r);
  polys = eval_jac(r, recurrence_a, recurrence_b, 0:(Nr-1));
  self.modal_coefficients(1:Nr,q) = inv(polys)*y(flags);
end

self = rmfield(self,{'x','y'});
