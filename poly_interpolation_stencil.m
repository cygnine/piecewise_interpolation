function[fz] = poly_interpolation_stencil(x,y,z,stencil)
% [FZ] = POLY_INTERPOLATION_STENCIL(X,Y,Z,STENCIL)
%
%     Performs piecewise polynomial interpolation on the grid points (X,Y) using
%     the interpolation stencil STENCIL. The interpolant is evaluated at the
%     points Z.

global handles;
newton = handles.speclab.NewtonPolynomials;

% Compute x values
XInput = x(stencil);

% Use stencil to compute interpolants
if size(stencil,2)==1
  dd = reshape(y,[1,n]);
else
  dd = newton.divided_difference(XInput.',y(stencil.'));
end

% Determine indicators for locations of nodes
% Temporarily redefine x as the bin separators to include all real numbers
x= [-inf; x(2:(end-1)); inf];
[temp,bin] = histc(z,x);

fz = zeros(size(z));

% For loops are bad, but I can't figure out how to vectorize this
% Compute matrix of locations to interpolate to
for q = 1:size(stencil,1)
  flags = bin==q;
  if any(flags)
    fz(flags) = newton.newton_evaluate(XInput(q,:),dd(:,q),z(flags));
  end
end
