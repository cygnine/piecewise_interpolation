function[y] = int(self)
% int -- Area under the curve of piecewise polynomial
% 
% y = int(self)
% self.int
%
%     Uses a quadrature rule exact for the order of the polynomial on each cell
%     to compute the exact integral.

switch self.basis_representation
case 'jacobi'
  alpha = self.opoly_opt.alpha;
  beta = self.opoly_opt.beta;

  % Value of P_0^(alpha,beta), to normalized the modal coefficients
  p0 = 1/sqrt(2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2));

  % Just need to add up 0th order modal coefficients
  y = p0*self.modal_coefficients(1,:)*2*self.jacobians;
otherwise
  error('Not yet implemented')
end
