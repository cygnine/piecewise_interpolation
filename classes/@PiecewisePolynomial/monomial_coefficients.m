function[m] = monomial_coefficients(self)
% monomial_coefficients -- compute global monomial coefficients
%
% m = monomial_coefficients(self)
% m = self.monomial_coefficients
%
%     Uses quadrature to compute the (global) monomial coefficients on each
%     cell.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
pwtools = handles.piecewise_interpolation.grid_tools;

[r,w] = jac.quad.gauss_quadrature(self.N,self.opoly_opt);
polys = jac.eval.eval_jacobi_poly(r,0:(self.N-1), self.opoly_opt);
vinv = polys'*spdiags(w,0,self.N,self.N);

x = pwtools.replicate_local_nodes(r, self.cell_boundaries);

connection_matrix = zeros([self.N*self.K, self.N]);
m = zeros([self.N, self.K]);

for q = 1:self.N;
  temp = x.^(q-1);
  temp = vinv*temp;
  connection_matrix(:,q) = temp(:);
end
for q = 1:self.K
  r1 = (q-1)*self.N + 1;
  r2 = q*self.N;
  connection_matrix(r1:r2,:) = inv(connection_matrix(r1:r2,:));
  m(:,q) = connection_matrix(r1:r2,:)*self.modal_coefficients(:,q);
end
