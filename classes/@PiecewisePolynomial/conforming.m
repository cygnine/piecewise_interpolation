function[tf] = conforming(self,other)
% conforming -- tests if two PiecewisePolynomials are conforming
%
% tf = conforming(self,other)
%
%     Tests to see if two PiecewisePolynomial objects have the same domain and
%     cell boundaries. If they don't, simple operations like addition and
%     multiplication throw errors.

tol = 1e-12;

if self.K ~= other.K
  tf = false;
  return;
elseif max(abs((self.cell_boundaries - other.cell_boundaries)))>tol
  tf = false;
  return;
else
  tf = true;
end
