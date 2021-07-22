% compute objective function for worst parameter finding
% D-optimality
% xi: vector of design points(1:n/2) and weights (n/2+1:end)
% beta: parameter values
% p:
function obj = worstparam_obj(xi, beta, p)
    
    % calc number of design points
    numpts = length(xi)/2;
    
    % design points and weights
    x = xi(1:numpts);
    w = xi(numpts+1:end);
    
    % compute eta
    lbeta = length(beta);
    
    % no checking here unlike in logistic.m
    % assuming a fractional poly is used
    eta = fracpoly(x, 2, beta, p);
    p1 = p(1);
    p2 = p(2);
    
    % weight function
    sigma = exp(eta)./(1 + exp(eta)).^2;
    
    % compute info matrix
    % technically, this is conditional on the polynomial powers
    M = 0;
    for i = 1:numpts
        
        m12 = H_j(x(i), 2, [1, p]);
        m13 = H_j(x(i), 3, [1, p]);
        m23 = H_j(x(i), 2, [1, p]) * H_j(x(i), 3, [1, p]);
        
        M_i = w(i) * sigma(i) * ...
            [1 m12 m13; ...
            m12 H_j(x(i), 2, [1, p])^2 m23; ...
            m13 m23 H_j(x(i), 3, [1, p])^2];
        
        M = M + M_i;
        
    end
    
    % D-optimality maximizes the log determinant
    % for a worst case, we minimize the log determinant
    obj = log(det(M));
    
end