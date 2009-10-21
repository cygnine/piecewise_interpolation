function[y] = diff(self)
% diff -- Differentiates a piecewise polynomial
%
% y = diff(self)
%
%     Returns a PiecewisePolynomial y that is the derivative of self. The number
%     of degrees of freedom per cell, y.N, is one less than self.N.

global handles;
mode_diff = handles.speclab.orthopoly1d.jacobi.operators.stiffness_operator;

temp = self.modal_coefficients;
for q = 1:self.K;
  switch self.basis_representation
  case 'jacobi'
    temp(:,q) = mode_diff(temp(:,q),self.opoly_opt);
    temp(:,q) = temp(:,q)/self.jacobians(q);
  otherwise
    error('Not yet implemented')
  end
end

y = PiecewisePolynomial('N', self.N-1, 'K', self.K, 'cell_boundaries',...
    self.cell_boundaries, 'basis_representation', self.basis_representation, ...
    'jacobi_alpha', self.opoly_opt.alpha, 'jacobi_beta', self.opoly_opt.beta, ...
    'modal_coefficients', temp(1:(self.N-1),:));
