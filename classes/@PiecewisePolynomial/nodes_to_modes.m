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
eval_jac = handles.speclab.orthopoly1d.eval_polynomial_standard.handle;
sss = handles.speclab.standard_scaleshift1d.handle;

[recurrence_a,recurrence_b] = jac.coefficients.recurrence(self.N+1,self.opoly_opt);

[bin_sums, bin_id] = histc(x, self.cell_boundaries);
x = x(:);
y = y(:);

assert(max(bin_sums)>self.N,  'Too many data points lie inside a cell: cannot uniquely interpolate');
assert(length(x)==length(y), 'Nodal data x and y must have the same length');

% erase previous data
self.modal_coefficients = zeros([self.N, self.K]);

r = (x - self.cell_shifts(bin_id))./self.jacobians(bin_id);
polys = eval_jac(r, recurrence_a, recurrence_b, 0:N);
for q = 1:self.K;
  flags = (bin_id==q);
  self.modal_coefficients(1:bin_sums(q),q) = inv(polys(flags,1:bin_sums(q)))*y(flags);
end
