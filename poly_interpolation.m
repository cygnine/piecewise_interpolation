function[fz] = poly_interpolation(x,y,z,varargin)
% [FZ] = POLY_INTERPOLATION(X,Y,Z,{K=3,BIAS=true,INTERVAL=false})
%
%     Performs piecewise polynomial interpolation on the grid points (X,Y). The
%     interpolation scheme is identical for all sub-intevals (except the
%     endpoints). Order K (i.e. nearest K+1 points) interpolation is used. The
%     interpolation is evaluated at the locations Z. The optional input BIAS
%     determines which direction to bias in the case of odd K. The default is
%     true = bias to the left. Set it to false to bias to the right.
%
%     The optional input INTERVAL denotes whether the interpolation stencil is
%     periodic. If a 2-vector is given, periodicity is assumed over that
%     bounding interval.

global handles;
fd = handles.FiniteDifference;
pw = handles.PiecewiseInterpolation;

opt = handles.common.InputSchema({'k','bias','interval'},{3,true,false},[],varargin{:});

% Force column vector
x = x(:);
y = y(:);

n = length(x);

if (opt.bias==false) && (mod(opt.k,2)==1)
  r = 1;
else
  r = 0;
end

if interval ~= false
  periodic=true;
else
  periodic=false;
end

[stencil,stencil_periodicity] = fd.difference_stencil(n,opt.k,'r',r*ones([n,1]),'periodic',periodic);

fz = pw.poly_interpolation_stencil(x,y,z,stencil...
        'stencil_periodicity', stencil_periodicity, 'interval', interval);
