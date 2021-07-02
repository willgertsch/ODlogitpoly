% more general logistic objective
function obj = logistic(vars, beta, opt, p)

    % Figure out how many design points
    numpts = length(vars)/2;
    
    % separate out design points and weights
    x = vars(1:numpts);
    w = vars(numpts+1:end);
    
    % compute eta
    lbeta = length(beta);
    if ~isempty(p) && lbeta == 3 % degree 2 frac poly
        eta = fracpoly(x, 2, beta, p);
        p1 = p(1);
        p2 = p(2);
    elseif lbeta == 2
       eta = beta(1) + beta(2) * x;
    elseif lbeta == 3
        eta = beta(1) + beta(2) * x + beta(3) * x.^2;
    else
        disp("Polnomial not supported")
        return
    end
    
    % weight function
    % need to update for higher degree poly
    sigma = exp(eta)./(1 + exp(eta)).^2;
    
    % Design matrix
    M = 0;
    for i = 1:numpts
        if lbeta == 3 && ~isempty(p)
            % take advantage of symmetry
            m12 = x(i)^p1;
            m13 = x(i)^p2;
            m23 = x(i)^(p1+p2);
            M_i = w(i)*sigma(i) * ...
                [1 m12 m13;...
                m12 x(i)^(2*p1) m23;...
                m13 m23 x(i)^(2*p2)];
        elseif lbeta == 3
            M_i = w(i)*sigma(i) * [1 x(i) x(i)^2; x(i) x(i)^2 x(i)^3;...
                x(i)^2 x(i)^3 x(i)^4];
        elseif lbeta == 2 
            M_i = [w(i)*sigma(i) w(i)*sigma(i)*x(i); ...
            w(i)*sigma(i)*x(i) w(i)*sigma(i)*x(i)^2];
        end
        
        M = M + M_i;
    end
    
    
    % switch between optimality objectives
    if (opt == 'D')
        obj = -log(det(M));
        
    elseif (opt == 'A')
        
        obj = trace(inv(M));
        
    elseif (opt == 'E') % experiment with E-opt
        obj = -min(eig(M));
    else % default at D optimality
        obj = -log(det(M));
    end

end