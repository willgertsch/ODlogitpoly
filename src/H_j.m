% helper function for adding log terms to fractional polynomial
function h = H_j(x, j, p)
    
    % base case
    if j == 1
        h = 1;
        return
    end
    
    if p(j) ~= p(j-1)
        h = BT(x, p(j));
    elseif p(j) == p(j-1)
        % recursive step
        h = H_j(x, j-1, p) .* log(x);
    end
    
    return
end