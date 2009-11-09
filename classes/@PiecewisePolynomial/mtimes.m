function[z] = mtimes(self,other)
% plus -- multiplies two piecewise polynomials
%
% z = mtimes(x,y)
%
%     Multiplies together two piecewise polynomials, or a piecewise polynomial and a
%     scalar double. The output z inherits the representation from the first
%     input.

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
pwtools = packages.piecewise_interpolation.grid_tools;

if isa(self,'double')
  z = other;
  z.modal_coefficients = self*(z.modal_coefficients);
  return
end

switch self.basis_representation
case 'jacobi'
  if isa(other,'PiecewisePolynomial')
    if not(self.conforming(other))
      error('Cannot multiply together non-conforming polynomials');
    end

    % If they're conforming, can multiply via quadrature
    temp.N = self.N + other.N - 1;
    temp.jacobi_alpha = self.opoly_opt.alpha;
    temp.jacobi_beta = self.opoly_opt.beta;
    temp.basis_representation = 'jacobi';
    temp.K = self.K;
    temp.cell_boundaries = self.cell_boundaries;

    [r,w] = jac.quad.gauss_quadrature(temp.N, self.opoly_opt);
    x = pwtools.replicate_local_nodes(r, self.cell_boundaries);
    polys = jac.eval.eval_jacobi_poly(r,0:(temp.N-1), self.opoly_opt);

    y = self.evaluate(x).*other.evaluate(x);
    temp.modal_coefficients = polys'*spdiags(w,0,temp.N,temp.N)*y;

    z = PiecewisePolynomial(temp);
  elseif isa(other, 'double')
    self.modal_coefficients = (self.modal_coefficients)*other;
    z = self;
  else
    error('Cannot add objects')
  end
otherwise
  error('Not yet implemented')
end
