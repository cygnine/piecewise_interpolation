function[scales,shifts] = compute_scaleshift(cell_boundaries)
% compute_scaleshift -- computes Jacobians and scaling factors
%
% [scales,shifts] = compute_scaleshift(cell_boundaries) 
%
%     Given the length-(K+1) vector cell_boundaries, returns the K affine scale
%     and shift parameters such that [-1,1]*scale(q) + shift(q) ------->
%     [cell_boundaries(q), cell_boundaries(q+1)].

cell_boundaries = cell_boundaries(:);
K = length(cell_boundaries) - 1;

scales = zeros([K 1]);
shifts = zeros([K 1]);

scales = diff(cell_boundaries)/2;
%shifts = mean([cell_boundaries(1:K), cell_boundaries(2:end)],2);
shifts = (cell_boundaries(1:K) + cell_boundaries(2:end))/2;

assert(all(scales)>=0, 'The input cell_boundaries must be a monotonically increasing vector');
