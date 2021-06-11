% sensitivity function for two parameter logistic model
function y = ch_logistic(z, vars, beta, opt)
    
    numpts = length(vars)/2;
    x = vars(1:numpts);
    w = vars(numpts+1:end);
    
%     sigma = exp(beta(1) + beta(2)*xvars)./(1 + exp(beta(1) + beta(2)*xvars)).^2;
%     g = exp(beta(1) + beta(2)*x)/(1+exp(beta(1) + beta(2)*x))^2;

     % compute eta
    lbeta = length(beta);
    if lbeta == 2
       eta = beta(1) + beta(2) * x;
       etaz = beta(1) + beta(2) * z;
    elseif lbeta == 3
        eta = beta(1) + beta(2) * x + beta(3) * x.^2;
        etaz = beta(1) + beta(2) * z + beta(3) * z^2;
    else
        disp("Polynomial not supported")
        return
    end
    
    % weight functions
    % need to update for higher degree poly
    sigma = exp(eta)./(1 + exp(eta)).^2;
    g = exp(etaz)/(1+exp(etaz))^2;
    
    % Design matrix
    M = 0;
    for i = 1:numpts
        if lbeta == 3
            M_i = w(i)*sigma(i) .* [1 x(i) x(i)^2; x(i) x(i)^2 x(i)^3;...
                x(i)^2 x(i)^3 x(i)^4];
        else 
            M_i = [w(i)*sigma(i) w(i)*sigma(i)*x(i); ...
            w(i)*sigma(i)*x(i) w(i)*sigma(i)*x(i)^2];
        end
        M = M + M_i;
    end
    
    %detM = det(M);
    
    Minv = inv(M);
    
    % switch optimality
    if (opt == 'D')
        if lbeta == 3
            y = g * [1 z z^2] * Minv * [1;z;z^2] - 3;
        else
            y = g * [1 z] * Minv * [1;z] - 2;
        end
        
    elseif (opt == 'A')
        if lbeta == 3
            y = g * [1 z z^2] * inv(M + 1e-10 * eye(2))^2 * [1;z;z^2] - trace(Minv);
        else
           y = g * [1 z] * inv(M + 1e-10 * eye(2))^2 * [1;z] - trace(Minv);
        end
        
    elseif (opt == 'E')
        
        [V, D] = eig(M); 
        [minlambda, index] = min(diag(D));
        v = V(:, index);
        E = v * v';
        
        if lbeta == 3
            y = g * [1 z z^2] * E * [1;z;z^2] - minlambda;
        else
            y = g * [1 z] * E * [1;z] - minlambda;
        end
        
    else
        if lbeta == 3
            y = g * [1 z z^2] * Minv * [1;z;z^2] - 2;
        else
            y = g * [1 z] * Minv * [1;z] - 2;
        end
    end
    
    return
end