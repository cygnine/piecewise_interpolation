function[y] = diff(self)
% diff -- Differentiates a piecewise polynomial
%
% y = diff(self)
%
%     Returns a PiecewisePolynomial y that is the derivative of self. The number
%     of degrees of freedom per cell, y.N, is one less than self.N.

global packages;
mode_diff = packages.speclab.orthopoly1d.jacobi.operators.stiffness_operator;

temp = self.modal_coefficients;
switch self.basis_representation
case 'jacobi'
  S = packages.speclab.orthopoly1d.jacobi.operators.stiffness_matrix(self.N, self.opoly_opt);
  %temp(:,q) = mode_diff(temp(:,q),self.opoly_opt);
  temp = S*temp*spdiags(1./self.jacobians,0,self.K,self.K);;
  %temp(:,q) = temp(:,q)/self.jacobians(q);
otherwise
  error('Not yet implemented')
end

y = PiecewisePolynomial('N', self.N-1, 'K', self.K, 'cell_boundaries',...
    self.cell_boundaries, 'basis_representation', self.basis_representation, ...
    'jacobi_alpha', self.opoly_opt.alpha, 'jacobi_beta', self.opoly_opt.beta, ...
    'modal_coefficients', temp(1:(self.N-1),:));
