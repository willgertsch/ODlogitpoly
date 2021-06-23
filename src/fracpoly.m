%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fractional Polynomials
% Reference:
% Royston P., Altman G.,
% "Regression using Fractional Polynomials of Continuous Covariates:
% Parsimonious Parameteric Modeling", Appl. Statist. 1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute fractional polynomial at x
% x: independent variable
% m: degree of polynomial
% beta: coefficients
% p: powers
function eta = fracpoly(x, m, beta, p)
    
    % user doesn't specify power for oth term
    % add power for internal use
    p = [1, p];
    
    % intercept term
    eta = beta(1);
    
    % x terms
    % loop will run twice for m = 3 and once for m = 2
    for j = 2:m+1
        
        eta = eta + beta(j) * H_j(x, j, p);
    end
    
    return

end
   