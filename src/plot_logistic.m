function plot_logistic(xi, beta, opt, low, up, p)
    
    x = linspace(low, up, 1000);
    y = arrayfun(@(x) ch_logistic(x, xi, beta, opt, p), x);
    plot(x,y)
    xlim([low, up])
    yline(0, '--')
    n = length(xi)/2;
    for i = 1:n
        xline(xi(i), "--r", sprintf('x=%0.3f\n w=%0.3f', xi(i), xi(n+i)),...
            'LabelVerticalAlignment', 'bottom',...
            'LabelHorizontalAlignment', 'center');
    end
    if length(beta) == 2
        title(sprintf('%s = (%0.2f, %0.2f)', '\beta', beta(1), beta(2)))
    elseif length(beta) == 3
        title(sprintf('%s = (%0.2f, %0.2f, %0.2f)', '\beta', beta(1), beta(2), beta(3)))
    end
    
    return
end