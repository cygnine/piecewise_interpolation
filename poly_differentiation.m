function[fz] = poly_differentiation(x,y,z,varargin)
% [FZ] = POLY_DIFFERENTIATION(X,Y,Z,{K:3,BIAS=TRUE)
%
%     Performs piecewise polynomial interpolation on the grid points (X,Y). The
%     interpolation scheme is identical for all sub-intevals (except the
%     endpoints). Order K (i.e. nearest K+1 points) interpolation is used. The
%     interpolation is differentiated and evaluated at the locations Z. The
%     optional input BIAS determines which direction to bias in the case of odd
%     K. The default is true = bias to the left. Set it to false to bias to the
%     right.

global handles;
fd = handles.FiniteDifference;
pw = handles.PiecewiseInterpolation;

opt = handles.common.InputSchema({'k','bias'},{3,true},[],varargin{:});

% Force column vector
x = x(:);
y = y(:);

n = length(x);

if (opt.bias==false) && (mod(opt.k,2)==1)
  r = 1;
else
  r = 0;
end

stencil = fd.difference_stencil(n,opt.k,'r',r*ones([n,1]));

fz = pw.poly_differentiation_stencil(x,y,z,stencil);
