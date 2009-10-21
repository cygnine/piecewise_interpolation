function[z] = plus(self,other)
% plus -- sums two piecewise polynomials
%
% z = plus(x,y)
% z = x.plus(y)
%
%     Sums together two piecewise polynomials, or a piecewise polynomial and a
%     scalar double.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;

switch self.basis_representation
case 'jacobi'
  if isa(other,'PiecewisePolynomial')
    if not(self.conforming(other))
      error('Cannot add together non-conforming polynomials');
    end

    % If they're conforming, addition is easy
    if self.N>=other.N
      z = self;
      z.modal_coefficients(1:other.N,:) = z.modal_coefficients(1:other.N,:) + ...
                                          other.modal_coefficients;
    else
      z = other;
      z.modal_coefficients(1:self.N,:) = z.modal_coefficients(1:self.N,:) + ...
                                          self.modal_coefficients;
    end
  elseif isa(other, 'double')
    val = jac.eval.eval_jacobi_poly(0,0,self.opoly_opt);
    self.modal_coefficients(1,:) = self.modal_coefficients(1,:) + other/val;
    z = self;
  else
    error('Cannot add objects')
  end
otherwise
  error('Not yet implemented')
end
