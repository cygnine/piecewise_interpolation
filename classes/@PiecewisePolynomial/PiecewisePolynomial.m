classdef PiecewisePolynomial
% PiecewisePolynomial -- A piecewise polynomial functional representation
%
% obj = PiecewisePolynomial({K=1, N=3, cell_boundaries=[-1,1], ...
%                            standard_interval=[-1,1], x=[], y=[]... 
%                            basis_representation='jacobi', jacobi_alpha=0,...
%                            jacobi_beta=0, modal_coefficients=[]})
%
%     Constructs a PiecewisePolynomial functional object with K elements, and
%     N degree of freedom on each element. The default internal data structure
%     storage calls for a nodal N-point Legendre-Gauss representation on each
%     subinterval. 
%
%     The input cell_boundaries determines the piecewise nature; it has length
%     (K+1), is non-decreasing, and defines the boundaries of each cell.
%
%     The optional input basis_representation determines the internal
%     representation of the basis functions. The default 'jacobi' means that
%     the polynomials are stored as modal coefficients of L2-normalized Jacobi
%     polynomials of order (alpha,beta) (default (0,0)). Other possibilities
%     not yet supported are:
%       - 'newton'  Newton polynomial representation
%       - 'monomial'  Monomial basis representation
%
%     A this stage, N-adaptivity is not supported; if you wish to have
%     different order polynomials on elements, just define N as the maximal
%     order, and upsample the lower-order elements accordingly.


properties
  K = 1; % The number of disjoint elements
  N = 4; % The number of degrees of freedom on each element
  global_interval = [-1,-1];
  standard_interval = [-1,1];
  cell_boundaries = [-1,1];
  opoly_opt = struct('alpha', 0, 'beta', 0);
  jacobians = 1;
  cell_shifts = 0;

  modal_coefficients = zeros([4 1]);
  basis_representation = 'jacobi';
end
  
methods

  function self = PiecewisePolynomial(varargin)
  % PiecewisePolynomial -- Class constructor

    persistent input_schema repgrid find_modes scaleshift
    if isempty(input_schema)
      from labtools import input_schema
      from piecewise_interpolation.grid_tools import replicate_local_nodes as repgrid
      from piecewise_interpolation.grid_tools import nodes_to_jacobi_modes as find_modes
      from piecewise_interpolation.grid_tools import compute_scaleshift as scaleshift
    end

    inputs = {'K','N','cell_boundaries', 'standard_interval','x', 'y', ...
              'basis_representation', 'jacobi_alpha', 'jacobi_beta', ...
              'modal_coefficients'};
    defaults = {1,3,[-1,1],[-1,1],[],[],'jacobi',0,0, []};

    %jac = packages.speclab.orthopoly1d.jacobi;
    %repgrid = packages.piecewise_interpolation.grid_tools.replicate_local_nodes.handle;
    %find_modes = packages.piecewise_interpolation.grid_tools.nodes_to_jacobi_modes;
    %scaleshift = packages.piecewise_interpolation.grid_tools.compute_scaleshift;

    opt = input_schema(inputs, defaults, [], varargin{:});
    [self.K, self.N, self.cell_boundaries, self.standard_interval, ...
      self.basis_representation] = deal(opt.K, opt.N, opt.cell_boundaries,...
      opt.standard_interval, opt.basis_representation);
    self.cell_boundaries = self.cell_boundaries(:);
    self.K = length(self.cell_boundaries) - 1;
    self.global_interval = [min(self.cell_boundaries), max(self.cell_boundaries)];
    self.opoly_opt.alpha = opt.jacobi_alpha;
    self.opoly_opt.beta = opt.jacobi_beta;

    % Allocation:
    self.modal_coefficients = zeros([self.N, self.K]);

    % Jacobians: dr/dx, where r is the standard interval [1;1], x is the local one
    [self.jacobians, self.cell_shifts] = scaleshift(self.cell_boundaries);

    switch self.basis_representation
    case 'jacobi'
      if ~isempty(opt.modal_coefficients)
        assert(prod(size(opt.modal_coefficients))==self.N*self.K, 'There must be exactly N*K modal coefficients');
        self.modal_coefficients = reshape(opt.modal_coefficients, [self.N, self.K]);
      elseif ~isempty(opt.x) & ~isempty(opt.y)
        self.modal_coefficients = find_modes(self.cell_boundaries, self.N, ...
           opt.x, opt.y, self.opoly_opt);
      end
    case 'newton'
      error('This basis type is not yet supported')
    case 'monomial'
      error('This basis type is not yet supported')
    otherwise
      error('Unrecognized basis type')
    end
  end

  varargout = subsref(self,varargin);
  y = evaluate(self,x);
  varargout = plot(self,varargin);
  y = diff(self);
  y = int(self);
  z = plus(self,other);
  z = mtimes(self,other);
  tf = conforming(self,other);
  m = monomial_coefficients(self);
  h = hilbert_transform(self,x,varargin);
end
end
