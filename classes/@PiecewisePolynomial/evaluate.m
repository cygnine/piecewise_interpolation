function[y] = evaluate(self,x)
% evaulate -- Evaluates piecewise polynomial function
%
% [y] = evaluate(x)
%       Evaluates the piecewise polynomial at the global locations x. If any
%       location x is outside the global_interval, the polynomial is
%       extrapolated.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
eval_jac = handles.speclab.orthopoly1d.eval_polynomial;
sss = handles.speclab.common.standard_scaleshift_1d.handle;

xsize = size(x);
x = x(:);
y = 0*x;
M = length(x);

cell_temp = [-Inf; self.cell_boundaries(2:(end-1)); Inf];
[bin_sums,bin] = histc(x, cell_temp);

local_coordinates = [];
cell_end_indices = zeros([self.K 1]);
cumulative_bin_sum = 0;
flags = false([M 1]);

switch self.basis_representation
case 'jacobi'
  [recurrence_a,recurrence_b] = jac.coefficients.recurrence(self.N+1,self.opoly_opt);

  for q = 1:self.K;
    flags = (bin==q);
    r = sss(x(flags),'scale',self.jacobians(q), 'shift', self.cell_shifts(q));

    %% Update indices so we know positions in local_coordinates
    %cumulative_bin_sum = cumulative_bin_sum + bin_sums(q);
    %cell_end_indices(q) = cumulative_bin_sum;
%
%    % Calculate local coordinates
%    if q==1
%      local_coordinates(1:cumulative_bin_sum) = r;
%    else
%      local_coordinates(cell_end_indices(q-1)+1:cumulative_bin_sum) = r;
%    end

    polys = eval_jac(r, recurrence_a, recurrence_b,0:(self.N-1));
    y(flags) = polys*self.modal_coefficients(:,q);
  end
end

y = reshape(y, xsize);
