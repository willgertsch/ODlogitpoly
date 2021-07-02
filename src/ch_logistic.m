% sensitivity function for two parameter logistic model
function y = ch_logistic(z, vars, beta, opt, p)
    
    numpts = length(vars)/2;
    x = vars(1:numpts);
    w = vars(numpts+1:end);

     % compute eta
    lbeta = length(beta);
    if lbeta == 3 && ~isempty(p)
        eta = fracpoly(x, 2, beta, p);
        etaz = fracpoly(z, 2, beta, p);
        p1 = p(1);
        p2 = p(2);
    elseif lbeta == 2
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
            M_i = w(i)*sigma(i) .* [1 x(i) x(i)^2; x(i) x(i)^2 x(i)^3;...
                x(i)^2 x(i)^3 x(i)^4];
        else 
            M_i = [w(i)*sigma(i) w(i)*sigma(i)*x(i); ...
            w(i)*sigma(i)*x(i) w(i)*sigma(i)*x(i)^2];
        end
        M = M + M_i;
    end
    
    %detM = det(M);
    
    %Minv = inv(M);
    
    % switch optimality
    if (opt == 'D')
        if lbeta == 3 && ~isempty(p)
            y = g * [1 z^p1 z^p2] * pinv(M) * [1;z^p1;z^p2] - 3;
        elseif lbeta == 3
            y = g * [1 z z^2] * (M \ [1;z;z^2]) - 3;
        else
            y = g * [1 z] * (M \ [1;z]) - 2;
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