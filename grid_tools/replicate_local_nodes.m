function[global_nodes] = replicate_local_nodes(local_nodes, cell_vertices)
% replicate_local_nodes -- creates global set of nodes from local node template
%
% global_nodes = replicate_local_nodes(local_nodes, cell_vertices)
%
%     Given a length-(K+1) vector cell_vertices containing non-decreasing
%     values, and a length-N vector local_nodes, this function returns an N x K
%     matrix where column k corresponds to the vector of nodes local_nodes
%     mapped from [-1,1] -----> [cell_vertices(k), cell_vertices(k+1)].

global packages;
scaleshift = packages.piecewise_interpolation.grid_tools.compute_scaleshift.handle;

local_nodes = local_nodes(:);
N = length(local_nodes);

K = length(cell_vertices(:)) - 1;

global_nodes = zeros([N K]);

[cell_scale, cell_shift] = scaleshift(cell_vertices);

global_nodes = repmat(local_nodes, [1,K])*spdiags(cell_scale,0,K,K);
global_nodes = global_nodes + repmat(cell_shift', [N, 1]);
