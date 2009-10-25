function[fz] = poly_interpolation(x,y,z,varargin)
% poly_interpolation -- piecewise polynomial interpolation of any order
%
% [fz] = poly_interpolation(x,y,z,{k=3,bias=true,interval=false})
%
%     Performs piecewise polynomial interpolation on the grid points (x,y). The
%     interpolation scheme is identical for all sub-intevals (except the
%     endpoints). Order k (i.e. nearest k+1 points) interpolation is used. The
%     interpolation is evaluated at the locations z. The optional input bias
%     determines which direction to bias in the case of odd k. The default is
%     true = bias to the left. Set it to false to bias to the right.
%
%     The optional input interval denotes whether the interpolation stencil is
%     periodic. If a 2-vector is given, periodicity is assumed over that
%     bounding interval.

global packages;
fd = packages.finite_difference;
pw = packages.piecewise_interpolation;

opt = packages.labtools.input_schema({'k','bias','interval'},{3,true,false},[],varargin{:});

% Force column vector
x = x(:);
y = y(:);

n = length(x);

if (opt.bias==false) && (mod(opt.k,2)==1)
  r = 1;
else
  r = 0;
end

if opt.interval ~= false
  periodic=true;
else
  periodic=false;
end

[stencil,stencil_periodicity] = fd.difference_stencil(n,opt.k,'r',r*ones([n,1]),'periodic',periodic);

fz = pw.poly_interpolation_stencil(x,y,z,stencil,...
        'stencil_periodicity', stencil_periodicity, 'interval', opt.interval);
