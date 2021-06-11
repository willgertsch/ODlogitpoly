% repair function to make weights sum to 1
% TODO: simplex design space
function vars = sum_constraints(vars, upper, lower)
    

    numpts = length(vars)/2;
    x = vars(1:numpts);
    w = vars(numpts+1:end);
    
    % scale weights
    S = sum(w);
    vars(numpts+1:end) = abs(w./S);
    
    % keep the bounding for the x's
    % this code is copied from UserProblem.m in PlatEMO
    vars(1:numpts) = max(min(x', repmat(upper, numpts, 1)), ...
        repmat(lower, numpts, 1))';
    
return
end