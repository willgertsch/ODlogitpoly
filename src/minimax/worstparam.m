% given a range of coef and a set of powers, find worst case param and
% fractional powers using a metaheuristic for a given design
% xi: a design with n design points. First n/2 entries of array are design
% points and the rest are weights
% b_low, b_up: lower and upper bounds for beta coef
% pset: set of possible fractional powers
% algo: algorithm options
function [b0, b1, b2, p1, p2] = worstparam(xi, b_low, b_up, pset, algo)

    % parse algorithm options
    method = [algo(1), algo(4:end)];
    maxFE = algo{2};
    swarm = algo{3};
    
    % function initialize swarm; uniformly gen points in coef interval
    init_low = min([b_low, pset]); 
    init_up = max([b_up, pset]);
    init = @(n, d) [unifrnd(init_low, init_up, n, pts),...
        unifrnd(0, 1, n, pts)];
    
    % set upper and lower bounds for coef and powers
    pmax = max(pset);
    pmin = min(pset);
    upperbounds = [b_up, max(pset), pmax];
    lowerbounds = [b_low, min(pset), pmin];
    
    % call PlatEMO to find optimal values for parameters
    [Dec, Obj, ~] = platemo('objFcn', @(x,d) worstparam_obj(xi, x), ...
        'lower', lowerbounds, 'upper', upperbounds, ...
        'algorithm', method, ...
        'decFcn', @(x,d) worstparam_constraints(x, b_up, b_low, pset), ...
        'maxFE', maxFE, ...
        'N', swarm, ...
        'initFcn', init);

end
