function[varargout] = plot(self)
% plot -- Plots the piecewise polynomial
%
% plot(self)
% handle = plot(self)
%
%     Plots the piecewise polynomial over its domain.

global handles;
pwtools = handles.piecewise_interpolation.grid_tools;
gq = handles.speclab.orthopoly1d.jacobi.quad.gauss_quadrature;

[r,w] = gq(self.N, self.opoly_opt);

x = pwtools.replicate_local_nodes(r, self.cell_boundaries);

varargout{1} = plot(x(:), self.evaluate(x(:)));
