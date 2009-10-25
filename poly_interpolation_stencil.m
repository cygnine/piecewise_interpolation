function[fz] = poly_interpolation_stencil(x,y,z,stencil,varargin)
% poly_interpolation_stencil -- piecewise polynomial interpolation of any order
%
% [fz] = poly_interpolation_stencil(x,y,z,stencil,{stencil_periodicity=false,interval=false})
%
%     Performs piecewise polynomial interpolation on the grid points (x,y) using
%     the interpolation stencil stencil. The interpolant is evaluated at the
%     points z.
%
%     The optional input stencil_periodicity determines whether or not a
%     periodic interpolation rubric is chosen. If it's false, no periodicity is
%     assumed. If periodicity is sought, input the value given as output from
%     finite_difference/difference_stencil.m. If periodicity is given, then the
%     second optional input interval must be a 2-vector denoted the bounding
%     interval of periodicity.

global packages;
newton = packages.speclab.newton_polynomials;

opt = packages.labtools.input_schema({'stencil_periodicity','interval'}, ...
         {false,false}, [], varargin{:});

% Compute x values
if opt.stencil_periodicity ~= false
  xmax = opt.interval(2); xmin = opt.interval(1);
  % Compute x values
  XInput = x(stencil);
  inds = opt.stencil_periodicity==1;
  % For indices that wrap down to 1:
  XInput(inds) = xmax + (XInput(inds) - xmin);

  inds = opt.stencil_periodicity==-1;
  % For indices that wrap up to n:
  XInput(inds) = xmin - (xmax - XInput(inds));

  % For histogram finding
  x = [x;xmax + x(1) - xmin]
else
  XInput = x(stencil);

  % For histogram finding
  x= [-inf; x(2:(end-1)); inf];
end

% Use stencil to compute interpolants
if size(stencil,2)==1
  dd = reshape(y,[1,n]);
else
  dd = newton.divided_difference(XInput.',y(stencil.'));
end

% Determine indicators for locations of nodes
% Temporarily redefine x as the bin separators to include all real numbers
[temp,bin] = histc(z,x);

fz = zeros(size(z));

% For loops are bad, but I can't figure out how to vectorize this
% Compute matrix of locations to interpolate to
for q = 1:size(stencil,1)
  flags = bin==q;
  if any(flags)
    fz(flags) = newton.newton_evaluate(XInput(q,:).',dd(:,q),z(flags));
  end
end
