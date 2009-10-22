function[varargout] = plot(self,varargin)
% plot -- Plots the piecewise polynomial
%
% plot(self,{plotstyle})
% handle = plot(self,{plotstyle})
%
%     Plots the piecewise polynomial over its domain. The optional input
%     plotstyle is a common Matlab plotting string.

global handles;
pwtools = handles.piecewise_interpolation.grid_tools;
gq = handles.speclab.orthopoly1d.jacobi.quad.gauss_quadrature;

[r,w] = gq(self.N, self.opoly_opt);

x = pwtools.replicate_local_nodes(r, self.cell_boundaries);

if isempty(varargin)
  varargout{1} = plot(x(:), self.evaluate(x(:)));
else
  varargout{1} = plot(x(:), self.evaluate(x(:)), varargin{1});
end
