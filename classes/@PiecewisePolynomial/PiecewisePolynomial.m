classdef PiecewisePolynomial

  properties
    K = 1; % The number of disjoint elements
    N = 3; % The polynomial order on each element
    global_interval = [0,0];
    local_intervals = zeros([K,2]);
  end
  
  methods

  function self = PiecewisePolynomial()
  % PiecewisePolynomial -- A piecewise polynomial functional representation
  %
  % obj = PiecewisePolynomial({K=1, N=3, global_interval=[-1,1], ...
  %                            local_interval=[-1;1], nodal_data = [0;0;0;0]})
  %
  %     Constructs a PiecewisePolynomial functional object with K elements, and
  %     N degree of freedom on each element. The internal data structure storage
  %     calls for a nodal N-point Chebyshev-Gauss representation on each
  %     subinterval. 
  %
  %     A this stage, N-adaptivity is not supported; if you wish to have
  %     different order polynomials on elements, just define N as the maximal
  %     order, and upsample the lower-order elements accordingly.

    global handles; 
    inputs = {'K','N','global_interval', 'local_interval','nodal_data'};
    defaults = {0,3,[-1,1],[-1;1],[0;0;0;0]};
    speclab = handles.speclab;

    self = handles.common.input_schema(inputs, defaults, [], varargin{:});

    % Jacobians: dr/dx, where r is the standard interval [1;1], x is the local one
    self.jacobians = 2./diff(self.local_interval,2);

    % Store Chebyshev vandermonde, inverse vandermonde, and differentiation
    % matrix.
    self.matrices.vandermonde = 
  end
end
