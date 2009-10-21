function[y] = int(self)
% int -- Area under the curve of piecewise polynomial
% 
% y = int(self)
% self.int
%
%     Uses a quadrature rule exact for the order of the polynomial on each cell
%     to compute the exact integral.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
pwtools = handles.piecewise_interpolation.grid_tools;

switch self.basis_representation
case 'jacobi'
  [r,w] = jac.quad.gauss_quadrature(self.N, self.opoly_opt);

  x = pwtools.replicate_local_nodes(r,self.cell_boundaries);
  z = self.evaluate(x);
  y = (w'*z)*self.jacobians;
otherwise
  error('Not yet implemented')
end
