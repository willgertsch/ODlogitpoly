% worstparam_constraints
% constraints for paramters for worst case D opt criterion
% make sure the coef are within user specified bounds and round powers to
% the closest match in the set of powers
% x: current solution
% b_up, b_low: vectors of upper and lower bounds for coef
% pset: set of powers
function x = worstparam_constraints(x, b_up, b_low, pset)

    % keep the coef within bounds
    x(1:3) = max(min(x', b_up), ...
        b_low)';
    
    % round powers to closest element of pset


end