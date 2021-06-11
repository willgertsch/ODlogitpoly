% more general logistic objective
function obj = logistic(vars, beta, opt)

    % Figure out how many design points
    numpts = length(vars)/2;
    
    % separate out design points and weights
    x = vars(1:numpts);
    w = vars(numpts+1:end);
    
    % compute eta
    lbeta = length(beta);
    if lbeta == 2
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
        if lbeta == 3
            M_i = w(i)*sigma(i) * [1 x(i) x(i)^2; x(i) x(i)^2 x(i)^3;...
                x(i)^2 x(i)^3 x(i)^4];
        else 
            M_i = [w(i)*sigma(i) w(i)*sigma(i)*x(i); ...
            w(i)*sigma(i)*x(i) w(i)*sigma(i)*x(i)^2];
        end
        
        M = M + M_i;
    end
    
    
    % switch between optimality objectives
    if (opt == 'D')
        obj = -log(det(M));
        
    elseif (opt == 'A')
        
        obj = trace(inv(M + eye(2)));
        
    elseif (opt == 'E') % experiment with E-opt
        obj = -min(eig(M));
    else % default at D optimality
        obj = -log(det(M));
    end

end