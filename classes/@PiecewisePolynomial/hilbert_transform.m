function[h] = hilbert_transform(self,x)
% hilbert_transform -- evaluates the Hilbert transform
%
% h = hilbert_transform(self, x)
%
%     Using high-order quadrature, this function evaluate the Hilbert transform
%     of the PiecewisePolynomial object with the kernel g using a high-order
%     quadrature. 

global packages;

switch self.basis_representation
case 'jacobi'
  h = packages.hilbert_transform.piecewise_polynomial_transform(...
        self.modal_coefficients, ...
        x, ...
        'cells', self.cell_boundaries, ...
        'alpha', self.opoly_opt.alpha, ...
        'beta', self.opoly_opt.beta);
otherwise
  error('Not yet implemented');
end
