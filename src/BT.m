% Box-Tidwell transformation for fractional polynomials
function xp = BT(x, p)
   
    if p ~= 0
        xp = x.^p;
    elseif p == 0
        xp = log(x);
    end
    
    return
    
end